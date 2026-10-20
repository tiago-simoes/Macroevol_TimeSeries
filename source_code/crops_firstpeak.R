# ============================================================================
#  crops_firstpeak.R  --  CROPS alternative to fixed-penalty (BIC) PELT
#  Source AFTER changepoint_pipeline.R and mc_uncertainty.R.
#
#  CROPS (Changepoints for a Range Of PenaltieS; Haynes et al. 2017, JCGS)
#  runs PELT for every penalty in [pen_min, pen_max] and returns ALL optimal
#  segmentations. Instead of fixing the penalty a priori (BIC), the number of
#  shifts is chosen from the "elbow" of the unpenalised cost vs number of
#  changepoints curve. The elbow is picked automatically as the point of
#  maximum perpendicular distance to the chord joining the curve's endpoints
#  (kneedle rule), on cost and n_cpts each rescaled to [0,1].
#  The first-peak rule is then applied to the elbow segmentation:
#    onset  = first upward shift (start of the rising phase that ends in the first peak; series start
#             if the series begins at peak level)
#    offset = first downward shift after the first peak (series end if no decline).
#  Optional n_fix (fixed model order, see model_order_overrides.csv) replaces the elbow by the solution with n_fix shifts.
#  Monte Carlo: replicates are CONDITIONED on the model order of the point estimate (n_fix =
#  observed elbow, or the fixed value), so the 95% MC range measures the timing uncertainty of the shifts
#  in the point-estimate segmentation; each replicate's own kneedle elbow is still recorded (n_elbow) to
#  report model-order stability.
# ============================================================================

# first-peak onset/offset from an arbitrary set of changepoints
peak_from_cpts <- function(v, age, cp) {
  bnds   <- c(0, sort(cp), length(v))
  starts <- utils::head(bnds, -1) + 1; ends <- bnds[-1]
  seg_mean <- mapply(function(s, e) mean(v[s:e]), starts, ends)
  is_peak <- function(i) {
    (if (i > 1) seg_mean[i] >= seg_mean[i - 1] else TRUE) &&
    (if (i < length(seg_mean)) seg_mean[i] >= seg_mean[i + 1] else TRUE)
  }
  pk <- which(vapply(seq_along(seg_mean), is_peak, logical(1)))[1]
  s <- starts[pk]; e <- ends[pk]
  # ONSET = FIRST UPWARD SHIFT. Every segment before the first local-max segment is
  # part of a strictly rising staircase (s1 < s2 < ... < s_pk), so the radiation starts at the first
  # changepoint whenever pk > 1.
  onset  <- if (pk == 1) age[s] else (age[ends[1]] + age[starts[2]]) / 2
  offset <- if (pk == length(seg_mean)) age[e] else (age[e] + age[starts[pk + 1]]) / 2
  peak_i <- s - 1 + which.max(v[s:e])
  list(onset = round(onset, 1), offset = round(offset, 1), peak_seg = pk,
       peak_age = round(age[peak_i], 1), seg_means = round(seg_mean, 3),
       starts = starts, ends = ends,
       shift_ages = round((age[ends[-length(ends)]] + age[starts[-1]]) / 2, 1))
}

seg_cost <- function(z, cp) {                     # unpenalised Normal-mean cost (RSS)
  bnds <- c(0, sort(cp), length(z))
  sum(mapply(function(s, e) sum((z[s:e] - mean(z[s:e]))^2),
             utils::head(bnds, -1) + 1, bnds[-1]))
}

# Exact least-squares segmentation with exactly n changepoints (segment length >= minseglen), by dynamic
# programming on the same Normal-mean (RSS) cost that PELT/CROPS minimise. Used when the model order is
# fixed (n_fix): CROPS only returns segmentations that are optimal for some penalty, so a given n can be
# missing from a replicate's CROPS path; this returns the optimal n-shift segmentation in every case.
exact_seg <- function(z, n, minseglen = 3) {
  N <- length(z); if (n <= 0) return(integer(0)); K <- n + 1
  if (N < K * minseglen) stop("series too short for n_fix")
  cs <- c(0, cumsum(z)); cs2 <- c(0, cumsum(z^2))
  cost <- function(a, b) (cs2[b + 1] - cs2[a]) - (cs[b + 1] - cs[a])^2 / (b - a + 1)
  F <- matrix(Inf, K, N); B <- matrix(NA_integer_, K, N)
  t1 <- minseglen:N; F[1, t1] <- cost(rep(1, length(t1)), t1)
  for (k in 2:K) for (t in (k * minseglen):N) {
    s <- ((k - 1) * minseglen):(t - minseglen)            # end of segment k-1
    val <- F[k - 1, s] + cost(s + 1, rep(t, length(s)))
    j <- which.min(val); F[k, t] <- val[j]; B[k, t] <- s[j]
  }
  cp <- integer(n); t <- N
  for (k in K:2) { t <- B[k, t]; cp[k - 1] <- t }
  cp
}

crops_first_peak <- function(grid, pen.range = c(0.5, 200), minseglen = 3, n_fix = NA) {
  # n_fix: fixed model order (override, or the point-estimate order inside the Monte Carlo) -> exact
  #        optimal segmentation with n_fix changepoints (exact_seg) instead of the kneedle elbow.
  v <- grid$value; age <- grid$age
  z <- (v - mean(v)) / stats::sd(v)
  utils::capture.output(                            # silence CROPS progress printing
    fit <- cpt.mean(z, method = "PELT", penalty = "CROPS",
                    pen.value = pen.range, minseglen = minseglen))
  full <- cpts.full(fit); pens <- pen.value.full(fit)
  sols <- lapply(seq_len(nrow(full)), function(i) {
    cp <- full[i, ]; cp <- sort(cp[!is.na(cp) & cp > 0]); cp
  })
  tab <- do.call(rbind, lapply(seq_along(sols), function(i) {
    cp <- sols[[i]]; fp <- peak_from_cpts(v, age, cp)
    data.frame(penalty = round(pens[i], 3), n_cpts = length(cp),
               cost = round(seg_cost(z, cp), 4),
               shift_ages_Ma = paste(fp$shift_ages, collapse = ";"),
               onset = fp$onset, offset = fp$offset, peak_age = fp$peak_age,
               stringsAsFactors = FALSE)
  }))
  # add the zero-changepoint solution if CROPS' upper penalty didn't reach it
  if (!any(tab$n_cpts == 0)) {
    fp0 <- peak_from_cpts(v, age, integer(0))
    tab <- rbind(tab, data.frame(penalty = NA, n_cpts = 0, cost = round(seg_cost(z, integer(0)), 4),
                                 shift_ages_Ma = "", onset = fp0$onset, offset = fp0$offset,
                                 peak_age = fp0$peak_age))
    sols <- c(sols, list(integer(0)))
  }
  # one row per n_cpts (keep lowest cost), ordered by n_cpts
  o   <- order(tab$n_cpts, tab$cost); tab <- tab[o, ]; sols <- sols[o]
  keep <- !duplicated(tab$n_cpts); tab <- tab[keep, ]; sols <- sols[keep]
  # kneedle elbow on (n_cpts, cost) rescaled to [0,1]
  x <- tab$n_cpts; y <- tab$cost
  if (length(x) >= 3) {
    xs <- (x - min(x)) / diff(range(x)); ys <- (y - min(y)) / diff(range(y))
    d  <- abs((ys[length(ys)] - ys[1]) * xs - (xs[length(xs)] - xs[1]) * ys +
              xs[length(xs)] * ys[1] - ys[length(ys)] * xs[1]) /
          sqrt((ys[length(ys)] - ys[1])^2 + (xs[length(xs)] - xs[1])^2)
    el <- which.max(d)
  } else el <- which.min(abs(x - 1))
  tab$elbow <- seq_len(nrow(tab)) == el
  sel <- el; cp_use <- sols[[el]]
  if (!is.na(n_fix)) {                                     # fixed model order: exact optimal n_fix-shift segmentation
    cp_use <- exact_seg(z, n_fix, minseglen)
    on_path <- which(tab$n_cpts == n_fix); sel <- if (length(on_path)) on_path[1] else NA
  }
  tab$selected <- seq_len(nrow(tab)) %in% sel
  fp <- peak_from_cpts(v, age, cp_use)
  c(fp, list(table = tab, n_cpts = length(cp_use), n_elbow = length(sols[[el]]), age = age))
}

mc_crops_first_peak <- function(df, n_rep = 1000, pen.range = c(0.5, 200),
                                minseglen = 3, step = 1, age_sd = 0,
                                nonneg = FALSE, seed = 1, n_fix = NA) {
  set.seed(seed)
  if (all(c("Err_lo", "Err_hi", "Err_type") %in% names(df))) {
    k <- ifelse(df$Err_type == "95CI", 1.96, 1)
    se_lo <- (df$Value - df$Err_lo) / k; se_hi <- (df$Err_hi - df$Value) / k
  } else se_lo <- se_hi <- df$Value_err
  se_lo[!is.finite(se_lo) | se_lo < 0] <- 0; se_hi[!is.finite(se_hi) | se_hi < 0] <- 0
  reps <- t(vapply(seq_len(n_rep), function(i) {
    z <- stats::rnorm(nrow(df))
    v <- df$Value + ifelse(z > 0, z * se_hi, z * se_lo); if (nonneg) v <- pmax(v, 0)
    a <- sort(pmax(df$Age_Ma + stats::rnorm(nrow(df), 0, age_sd), 0), decreasing = TRUE)
    g <- resample_grid(a, v[order(-df$Age_Ma)], step = step)
    fp <- tryCatch(crops_first_peak(g, pen.range, minseglen, n_fix), error = function(e) NULL)
    if (is.null(fp)) rep(NA, 5) else c(fp$onset, fp$offset, fp$peak_age, fp$n_cpts, fp$n_elbow)
  }, numeric(5)))
  # n_cpts = number of shifts used in the replicate; n_elbow = kneedle elbow of that replicate (model-order stability)
  colnames(reps) <- c("onset", "offset", "peak_age", "n_cpts", "n_elbow")
  reps <- as.data.frame(reps); reps$mc_mode <- if (is.na(n_fix)) "free" else "conditioned_exact"
  reps
}
