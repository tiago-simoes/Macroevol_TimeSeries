# ============================================================================
#  seg_firstpeak.R -- segmented (continuous piecewise-linear) regression for
#  clade-origin time series: BIC model order at 5-Myr resolution, first-peak
#  onset/offset (hybrid rule by curve type), model CIs and Monte-Carlo uncertainty.
#
#  Base R only (CRAN package `segmented` not required). Breakpoints are estimated
#  with Muggeo's (2003, Stat. Med. 22:3055) iterative linearisation -- the algorithm
#  of R `segmented` -- started from an exhaustive grid (k <= 2), a greedy forward
#  search (k-1 solution + one knot) and random restarts (cf. bootstrap restarting,
#  Wood 2001; Muggeo 2008); the minimum-RSS valid solution is kept.
#  Breakpoint SE (delta method at convergence, as segmented::confint.segmented):
#      Var(psi_j) = [Var(d_j) + (d_j/g_j)^2 Var(g_j) - 2 (d_j/g_j) Cov(d_j,g_j)] / g_j^2
#  where g_j = slope change and d_j = coefficient of -I(t > psi_j).
#
#  Time runs forward inside: t = (oldest age - Age_Ma), Myr since series start.
#  Model: y = b0 + b1 t + sum_j g_j (t - psi_j)_+ ; slope of segment j = b1 + sum_{i<j} g_i
#
#  Dense traced curves are averaged into bins of w = min(BIN_W, span / MIN_BINS) Myr (5 Myr, or
#  finer for short series so each keeps >= 20 bins; binned only if median spacing < 0.8 w) for
#  model-order selection (BIC, smallest k within DBIC_TIE of the minimum), and their SEs are
#  scaled to that effective sample size (x sqrt(n / n_bins)); breakpoints themselves are
#  estimated on the full-resolution data.
# ============================================================================

SEG_PAR <- list(BIN_W = 5, MIN_BINS = 20, BIN_GAP_FRAC = 0.8, KMAX = 6, MINSEG = 3, DBIC_TIE = 2,
                DELTA_MIN = 0.10,   # hump rule: a segment is rising/declining if its change >= 10% of the fitted range (and slope significant)
                RHO_MIN = 0.25,     # cumulative rule: slope change >= 25% of the larger adjacent slope (and significant)
                DD_MAX = 0.20, NET_MIN = 0.70)   # cumulative curve: drawdown <= 20% of range and net rise >= 70% of range

hinge <- function(t, psi) {
  if (!length(psi)) return(matrix(0, length(t), 0))
  vapply(psi, function(p) pmax(t - p, 0), numeric(length(t)))
}
vterm <- function(t, psi) vapply(psi, function(p) -as.numeric(t > p), numeric(length(t)))

valid_psi <- function(t, psi, minseg) {
  if (!length(psi)) return(TRUE)
  if (any(!is.finite(psi)) || is.unsorted(psi)) return(FALSE)
  if (any(psi <= min(t)) || any(psi >= max(t))) return(FALSE)
  cnt <- tabulate(findInterval(t, psi, left.open = TRUE) + 1L, nbins = length(psi) + 1L)
  all(cnt >= minseg)
}

rss_psi <- function(t, y, psi) {
  f <- .lm.fit(cbind(1, t, hinge(t, psi)), y)
  if (f$rank < 2 + length(psi)) Inf else sum(f$residuals^2)
}

# Muggeo iteration from one start; returns list(psi, rss) or NULL
muggeo <- function(t, y, psi0, minseg = 3, maxit = 60, tol = 1e-6) {
  psi <- sort(psi0); k <- length(psi)
  if (!k) return(list(psi = numeric(0), rss = rss_psi(t, y, psi)))
  if (!valid_psi(t, psi, minseg)) return(NULL)
  rss <- rss_psi(t, y, psi); span <- diff(range(t))
  for (it in seq_len(maxit)) {
    cf <- lm.fit(cbind(1, t, hinge(t, psi), vterm(t, psi)), y)$coefficients
    if (any(is.na(cf))) break
    gam <- cf[2 + seq_len(k)]; del <- cf[2 + k + seq_len(k)]
    stp <- ifelse(abs(gam) > 1e-12, del / gam, 0)
    h <- 1; moved <- FALSE
    for (j in 1:12) {                                  # step halving: valid knots, non-increasing RSS
      pn <- psi + h * stp
      if (valid_psi(t, pn, minseg)) { rn <- rss_psi(t, y, pn); if (rn <= rss * (1 + 1e-10)) { moved <- TRUE; break } }
      h <- h / 2
    }
    if (!moved) break
    conv <- max(abs(pn - psi)) < tol * span
    psi <- pn; rss <- rn
    if (conv) break
  }
  list(psi = psi, rss = rss)
}

cand_knots <- function(t, minseg, m = Inf) {
  u <- sort(unique(t)); mid <- (u[-1] + u[-length(u)]) / 2
  mid <- mid[vapply(mid, function(p) valid_psi(t, p, minseg), logical(1))]
  if (length(mid) > m) mid <- mid[unique(round(seq(1, length(mid), length.out = m)))]
  mid
}
random_psi <- function(t, k, minseg, tries = 200) {
  for (i in seq_len(tries)) { p <- sort(stats::runif(k, min(t), max(t))); if (valid_psi(t, p, minseg)) return(p) }
  NULL
}

# best k-breakpoint fit (minimum RSS over starts)
fit_k <- function(t, y, k, minseg = 3, n_random = 20, starts = list(), prev = NULL, grid = TRUE, refine = grid,
                  hops = if (refine) 10 else 0, m_cand = 150, n_refine = 5) {
  if (k == 0) return(list(psi = numeric(0), rss = rss_psi(t, y, numeric(0))))
  S <- starts
  S[[length(S) + 1]] <- stats::quantile(t, seq_len(k) / (k + 1), names = FALSE)
  if (grid && k == 1) {
    ck <- cand_knots(t, minseg); r <- vapply(ck, function(p) rss_psi(t, y, p), numeric(1))
    if (length(r)) S[[length(S) + 1]] <- ck[which.min(r)]
  }
  if (grid && k == 2) {
    ck <- cand_knots(t, minseg, 60); best <- Inf; bp <- NULL
    if (length(ck) >= 2) for (a in seq_along(ck)[-length(ck)]) for (b in (a + 1):length(ck)) {
      p <- ck[c(a, b)]; if (!valid_psi(t, p, minseg)) next
      r <- rss_psi(t, y, p); if (r < best) { best <- r; bp <- p } }
    if (!is.null(bp)) S[[length(S) + 1]] <- bp
  }
  if (grid && k >= 3 && length(t) <= 150) {               # (near-)exhaustive search on a candidate grid with C(m,k) <= 1e5
    m <- k + 1; while (choose(m + 1, k) <= 1e5) m <- m + 1
    ck <- cand_knots(t, minseg, m)
    if (length(ck) >= k) {
      cb <- utils::combn(length(ck), k); r <- rep(Inf, ncol(cb))
      for (c in seq_len(ncol(cb))) { p <- ck[cb[, c]]; if (valid_psi(t, p, minseg)) r[c] <- rss_psi(t, y, p) }
      for (c in head(order(r), 10)) if (is.finite(r[c])) S[[length(S) + 1]] <- ck[cb[, c]]
    }
  }
  if (!is.null(prev) && length(prev) == k - 1) {
    ck <- cand_knots(t, minseg, 80)
    r <- vapply(ck, function(p) { q <- sort(c(prev, p)); if (valid_psi(t, q, minseg)) rss_psi(t, y, q) else Inf }, numeric(1))
    if (length(r) && is.finite(min(r))) S[[length(S) + 1]] <- sort(c(prev, ck[which.min(r)]))
  }
  for (i in seq_len(n_random)) { p <- random_psi(t, k, minseg); if (!is.null(p)) S[[length(S) + 1]] <- p }
  sols <- list()
  for (s in S) {
    if (length(s) != k) next
    f <- tryCatch(muggeo(t, y, s, minseg), error = function(e) NULL)
    if (!is.null(f) && is.finite(f$rss)) sols[[length(sols) + 1]] <- f
  }
  if (!length(sols)) return(NULL)
  rs <- vapply(sols, `[[`, 0, "rss"); sols <- sols[order(rs)]
  if (refine && k >= 2) {                                 # refine the 5 best distinct local optima
    key <- vapply(sols, function(z) paste(round(z$psi, 2), collapse = ";"), "")
    sols <- sols[!duplicated(key)][seq_len(min(n_refine, sum(!duplicated(key))))]
    sols <- lapply(sols, function(z) coord_refine(t, y, z, minseg, m_cand))
    sols <- sols[order(vapply(sols, `[[`, 0, "rss"))]
  }
  best <- sols[[1]]
  if (hops > 0 && k >= 1) {                               # basin hopping: jitter the best knots, re-optimise
    span <- diff(range(t))
    for (h in seq_len(hops)) {
      p <- sort(best$psi + stats::rnorm(k, 0, 0.08 * span))
      if (!valid_psi(t, p, minseg)) next
      g <- tryCatch(muggeo(t, y, p, minseg), error = function(e) NULL)
      if (is.null(g)) next
      if (k >= 2) g <- coord_refine(t, y, g, minseg, m_cand)
      if (g$rss < best$rss * (1 - 1e-9)) best <- g
    }
  }
  best
}

# coordinate-descent refinement: move one knot at a time over the candidate grid (others fixed),
# until no single move lowers the RSS; then polish with Muggeo. Guards against local optima.
coord_refine <- function(t, y, f, minseg, m = 150, max_sweep = 10) {
  ck <- cand_knots(t, minseg, m); psi <- f$psi; rss <- f$rss; k <- length(psi)
  for (sw in seq_len(max_sweep)) {
    improved <- FALSE
    for (j in seq_len(k)) {
      lo <- if (j > 1) psi[j - 1] else -Inf; hi <- if (j < k) psi[j + 1] else Inf
      cj <- ck[ck > lo & ck < hi]
      if (!length(cj)) next
      r <- vapply(cj, function(p) { q <- psi; q[j] <- p; if (valid_psi(t, q, minseg)) rss_psi(t, y, q) else Inf }, numeric(1))
      if (min(r) < rss * (1 - 1e-9)) { psi[j] <- cj[which.min(r)]; rss <- min(r); improved <- TRUE }
    }
    if (!improved) break
  }
  g <- tryCatch(muggeo(t, y, psi, minseg), error = function(e) NULL)
  if (!is.null(g) && g$rss <= rss) g else list(psi = psi, rss = rss)
}

bic_pwl <- function(rss, n, k) { ll <- -n / 2 * (log(2 * pi) + log(max(rss, 1e-300) / n) + 1); -2 * ll + (2 * k + 3) * log(n) }

# model-order selection by BIC over 0..Kmax (on the series given: binned for dense curves)
select_k <- function(t, y, Kmax = 6, minseg = 3, n_random = 150, hops = 10, n_refine = 150, grid = TRUE, refine = TRUE) {
  n <- length(t); Kmax <- max(0, min(Kmax, floor(n / minseg) - 1))
  fits <- list(); prev <- NULL
  for (k in 0:Kmax) {
    f <- fit_k(t, y, k, minseg, n_random, prev = prev, hops = hops, n_refine = n_refine, grid = grid, refine = refine)
    if (is.null(f)) break
    fits[[k + 1]] <- f; prev <- f$psi
  }
  K <- length(fits) - 1                                   # backward pass: drop one knot from the (k+1) solution
  if (K >= 2) for (k in (K - 1):1) {
    up <- fits[[k + 2]]$psi
    for (j in seq_along(up)) {
      g <- tryCatch(muggeo(t, y, up[-j], minseg), error = function(e) NULL)
      if (!is.null(g) && k >= 2 && refine) g <- coord_refine(t, y, g, minseg)
      if (!is.null(g) && g$rss < fits[[k + 1]]$rss * (1 - 1e-9)) fits[[k + 1]] <- g
    }
  }
  tab <- data.frame(k = seq_along(fits) - 1, rss = vapply(fits, `[[`, 0, "rss"))
  tab$BIC <- mapply(function(r, k) bic_pwl(r, n, k), tab$rss, tab$k)
  tab$dBIC <- tab$BIC - min(tab$BIC)
  # parsimony: the smallest k within DBIC_TIE of the BIC minimum (Kass & Raftery 1995: dBIC <= 2 = no real support)
  list(table = tab, k = min(tab$k[tab$dBIC <= SEG_PAR$DBIC_TIE]), k_min = tab$k[which.min(tab$BIC)], fits = fits, Kmax = Kmax)
}

# 5-Myr binning of dense curves (time forward); returns list(t, y, binned)
# bin width = BIN_W Myr, or finer for short series so that each has >= MIN_BINS bins;
# a series is binned only if it is dense relative to that width (median spacing < 0.8 w)
bin_width <- function(t) min(SEG_PAR$BIN_W, diff(range(t)) / SEG_PAR$MIN_BINS)
bin_series <- function(t, y, w = bin_width(t)) {
  gap <- stats::median(abs(diff(sort(t))))
  if (!is.finite(gap) || gap >= SEG_PAR$BIN_GAP_FRAC * w) return(list(t = t, y = y, binned = FALSE, w = NA))
  b <- floor((t - min(t)) / w + 1e-9)
  list(t = as.numeric(tapply(t, b, mean)), y = as.numeric(tapply(y, b, mean)), binned = TRUE, w = w)
}

# inference for a fitted psi: slopes and SEs, slope changes, breakpoint SEs; se_mult/df_eff scale
# the SEs to the effective sample size (dense curves: sqrt(n / n_bins))
pwl_inference <- function(t, y, psi, se_mult = 1, n_eff = length(t)) {
  k <- length(psi); n <- length(t)
  X <- cbind(1, t, hinge(t, psi)); f <- lm.fit(X, y); cf <- f$coefficients
  s2 <- sum(f$residuals^2) / max(n - (2 + 2 * k), 1)
  V <- s2 * ginv(crossprod(X))
  L <- cbind(0, 1, matrix(0, k + 1, k)); if (k) for (j in 2:(k + 1)) L[j, 2 + seq_len(j - 1)] <- 1
  slope <- as.numeric(L %*% cf); slope_se <- sqrt(pmax(diag(L %*% V %*% t(L)), 0)) * se_mult
  gam <- if (k) as.numeric(cf[2 + seq_len(k)]) else numeric(0)
  gam_se <- if (k) sqrt(pmax(diag(V)[2 + seq_len(k)], 0)) * se_mult else numeric(0)
  psi_se <- rep(NA_real_, k)
  if (k) {
    Xa <- cbind(X, vterm(t, psi)); fa <- lm.fit(Xa, y); ca <- fa$coefficients
    if (!any(is.na(ca))) {
      Va <- sum(fa$residuals^2) / max(n - ncol(Xa), 1) * ginv(crossprod(Xa))
      ig <- 2 + seq_len(k); id <- 2 + k + seq_len(k); g <- ca[ig]; d <- ca[id]
      psi_se <- sqrt(pmax((diag(Va)[id] + (d / g)^2 * diag(Va)[ig] - 2 * (d / g) * diag(Va[id, ig, drop = FALSE])) / g^2, 0)) * se_mult
    }
  }
  df <- max(n_eff - (2 + 2 * k), 1)
  list(coef = cf, fitted = as.numeric(X %*% cf), slope = slope, slope_se = slope_se, gam = gam, gam_se = gam_se,
       psi = psi, psi_se = psi_se, tcrit = stats::qt(0.975, df), df = df, rss = sum(f$residuals^2))
}
ginv <- function(X) { s <- svd(X); d <- ifelse(s$d > max(s$d) * 1e-12, 1 / s$d, 0); s$v %*% (d * t(s$u)) }

# ---------------------------------------------------------------------------
# Onset / offset rules (segments j = 1..k+1, time forward)
#  hump (rates, disparity, hump-shaped diversity; PELT-style definitions):
#    segment class R / F / D = rising / flat / declining: slope significant (|slope| > tcrit SE)
#    AND |change over the segment| >= DELTA_MIN x fitted range.
#    first peak = end of the rising/plateau phase = start of the first D segment (series end if none);
#    onset  = start of the first R segment before it (series start if none, or if segment 1 is R);
#    offset = start of that first D segment (a plateau after the rise belongs to the burst).
#  cumulative (monotonic diversity curves, e.g. lineages through time; slope-change definitions):
#    breakpoint class A / Dc = acceleration / deceleration: slope change significant AND
#    >= RHO_MIN x the larger of the two adjacent slopes.
#    offset = first deceleration (series end if none); onset = first acceleration before it (series start if none).
# ---------------------------------------------------------------------------
first_peak_pwl <- function(t, y, psi, curve = c("hump", "cumulative"), t0, inf,
                           delta_min = SEG_PAR$DELTA_MIN, rho_min = SEG_PAR$RHO_MIN) {
  curve <- match.arg(curve)
  k <- length(psi); lo <- min(t); hi <- max(t)
  b <- c(lo, psi, hi); dur <- diff(b)
  R <- diff(range(inf$fitted)); if (!is.finite(R) || R <= 0) R <- 1
  sig <- abs(inf$slope) > inf$tcrit * inf$slope_se
  rel <- inf$slope * dur / R
  cls <- ifelse(sig & rel >= delta_min, "R", ifelse(sig & rel <= -delta_min, "D", "F"))
  bpc <- rep("", k); dec_end <- NA; n_steps <- 0
  if (curve == "hump") {
    jD <- which(cls == "D")[1]
    jP <- if (is.na(jD)) k + 1 else jD - 1
    rs <- if (jP >= 1) which(cls[seq_len(jP)] == "R") else integer(0)
    on_t <- if (length(rs)) b[rs[1]] else lo
    off_t <- if (is.na(jD)) hi else b[jD]
    n_steps <- if (length(rs)) sum(diff(c(-1, rs)) > 1) else 0          # separate rising phases before the peak
    if (!is.na(jD)) { jj <- jD; while (jj < k + 1 && cls[jj + 1] == "D") jj <- jj + 1; dec_end <- b[jj + 1] }
    w <- which(t >= on_t - 1e-9 & t <= off_t + 1e-9)
    pk_t <- if (length(w)) t[w][which.max(inf$fitted[w])] else on_t
    if (!is.na(jD) && jD == 1) pk_t <- lo
    main_t <- if (length(rs)) b[rs[which.max(inf$slope[rs])]] else NA      # start of the steepest rising segment
    later <- numeric(0)
  } else {
    if (k) {
      sg <- abs(inf$gam) > inf$tcrit * inf$gam_se
      big <- abs(inf$gam) >= rho_min * pmax(abs(inf$slope[-1]), abs(inf$slope[-(k + 1)]))
      bpc <- ifelse(sg & big & inf$gam > 0, "A", ifelse(sg & big & inf$gam < 0, "Dc", "0"))
    }
    jd <- which(bpc == "Dc")[1]
    ja <- which(bpc[seq_len(if (is.na(jd)) k else jd - 1)] == "A")
    on_t <- if (length(ja)) psi[ja[1]] else lo
    off_t <- if (is.na(jd)) hi else psi[jd]
    n_steps <- length(ja)
    w <- which(t >= on_t - 1e-9 & t <= off_t + 1e-9)                    # steepest segment inside the burst
    segs <- which(b[-length(b)] >= on_t - 1e-9 & b[-1] <= off_t + 1e-9)
    js <- if (length(segs)) segs[which.max(inf$slope[segs])] else 1
    pk_t <- (b[js] + b[js + 1]) / 2
    main_t <- b[js]
    later <- if (k && !is.na(jd) && jd < k) psi[which(bpc == "A" & seq_len(k) > jd)] else numeric(0)
  }
  list(onset = t0 - on_t, offset = t0 - off_t, peak_age = t0 - pk_t,
       decline_end = if (is.na(dec_end)) NA_real_ else t0 - dec_end,
       main_rise = if (is.na(main_t)) NA_real_ else t0 - main_t, later_acc = t0 - later,
       seg_class = cls, bp_class = bpc, n_steps = n_steps,
       left_cens = abs(on_t - lo) < 1e-9, right_cens = abs(off_t - hi) < 1e-9,
       start_decline = curve == "hump" && abs(off_t - lo) < 1e-9,
       onset_bp = if (abs(on_t - lo) < 1e-9) NA_integer_ else which.min(abs(psi - on_t)),
       offset_bp = if (abs(off_t - hi) < 1e-9 || abs(off_t - lo) < 1e-9) NA_integer_ else which.min(abs(psi - off_t)))
}

curve_type <- function(age, value, parameter, dd_max = SEG_PAR$DD_MAX, net_min = SEG_PAR$NET_MIN) {
  o <- order(-age); v <- value[o]; R <- diff(range(v))
  dd <- max(cummax(v) - v) / R; net <- (v[length(v)] - v[1]) / R
  list(type = if (parameter == "taxonomic_diversity" && dd <= dd_max && net >= net_min) "cumulative" else "hump",
       drawdown = dd, net = net)
}

# ---------------------------------------------------------------------------
# full analysis of one series (point estimate)
#   df: Age_Ma, Value (any order); parameter: standardised parameter name
# ---------------------------------------------------------------------------
seg_analyse <- function(age, value, curve, k_fix = NA, Kmax = SEG_PAR$KMAX, n_random = 50, starts = list()) {
  o <- order(-age); age <- age[o]; value <- value[o]
  t0 <- max(age); t <- t0 - age; n <- length(t)
  bs <- bin_series(t, value)
  minseg_full <- if (bs$binned) max(SEG_PAR$MINSEG, ceiling(bs$w / stats::median(diff(t)))) else SEG_PAR$MINSEG
  sel <- select_k(bs$t, bs$y, Kmax, SEG_PAR$MINSEG)
  k <- if (is.na(k_fix)) sel$k else k_fix
  # breakpoints on the full-resolution data: greedy sequence 1..k + binned solution as start
  prev <- NULL; f <- list(psi = numeric(0), rss = rss_psi(t, value, numeric(0)))
  if (k > 0) for (kk in seq_len(k)) {
    st <- starts
    if (length(sel$fits) >= kk + 1 && length(sel$fits[[kk + 1]]$psi) == kk) st[[length(st) + 1]] <- sel$fits[[kk + 1]]$psi
    f <- fit_k(t, value, kk, minseg_full, if (kk == k) n_random else 10, starts = st, prev = prev,
               refine = kk == k, hops = if (kk == k) 20 else 0)
    if (is.null(f)) stop("no valid ", kk, "-breakpoint fit")
    prev <- f$psi
  }
  n_eff <- length(bs$t); se_mult <- sqrt(n / n_eff)
  inf <- pwl_inference(t, value, f$psi, se_mult, n_eff)
  fp <- first_peak_pwl(t, value, f$psi, curve, t0, inf)
  list(age = age, value = value, t = t, t0 = t0, k = k, k_bic = sel$k, sel = sel, binned = bs$binned,
       bins = bs, n_eff = n_eff, se_mult = se_mult, minseg_full = minseg_full,
       psi = f$psi, psi_age = t0 - f$psi, inf = inf, fp = fp, curve = curve)
}

# ---------------------------------------------------------------------------
# Monte Carlo: perturb values (split-normal from error bars) and ages (Normal, sd = age_sd),
# refit with the point-estimate k (conditioned), same rule and curve type; the replicate's own
# BIC-best k (same binned resolution) is recorded for the first n_bic replicates (model-order stability).
# ---------------------------------------------------------------------------
MC_VERSION <- "kbic-refined-100"
mc_seg <- function(df, pe, n_rep = 1000, age_sd = 0, nonneg = FALSE, n_bic = 100, n_random = 4, seed = 1) {
  set.seed(seed)
  if (all(c("Err_lo", "Err_hi", "Err_type") %in% names(df))) {
    kk <- ifelse(df$Err_type == "95CI", 1.96, 1)
    se_lo <- (df$Value - df$Err_lo) / kk; se_hi <- (df$Err_hi - df$Value) / kk
  } else se_lo <- se_hi <- df$Value_err
  se_lo[!is.finite(se_lo) | se_lo < 0] <- 0; se_hi[!is.finite(se_hi) | se_hi < 0] <- 0
  k <- pe$k
  out <- matrix(NA_real_, n_rep, 7 + 2 * max(k, 1),
                dimnames = list(NULL, c("onset", "offset", "peak_age", "decline_end", "left_cens", "right_cens", "k_bic",
                                        paste0("psi", seq_len(max(k, 1))), paste0("psi_se", seq_len(max(k, 1))))))
  for (i in seq_len(n_rep)) {
    z <- stats::rnorm(nrow(df))
    v <- df$Value + ifelse(z > 0, z * se_hi, z * se_lo); if (nonneg) v <- pmax(v, 0)
    a <- pmax(df$Age_Ma + stats::rnorm(nrow(df), 0, age_sd), 0)
    o <- order(-a); a <- a[o]; v <- v[o]
    t0 <- max(a); t <- t0 - a
    st <- list(); p0 <- sort(t0 - pe$psi_age); if (length(p0) == k && k > 0) st[[1]] <- p0
    f <- tryCatch(fit_k(t, v, k, pe$minseg_full, n_random, starts = st, grid = FALSE), error = function(e) NULL)
    if (is.null(f)) next
    inf <- tryCatch(pwl_inference(t, v, f$psi, pe$se_mult, pe$n_eff), error = function(e) NULL)
    if (is.null(inf)) next
    fp <- tryCatch(first_peak_pwl(t, v, f$psi, pe$curve, t0, inf), error = function(e) NULL)
    if (is.null(fp)) next
    kb <- NA
    if (i <= n_bic) {
      bs <- bin_series(t, v)
      sk <- tryCatch(select_k(bs$t, bs$y, SEG_PAR$KMAX, SEG_PAR$MINSEG, n_random = 10, hops = 2, n_refine = 5, grid = FALSE, refine = TRUE), error = function(e) NULL)
      if (!is.null(sk)) kb <- sk$k
    }
    pa <- if (k) t0 - f$psi else NA; ps <- if (k) inf$psi_se else NA
    out[i, ] <- c(fp$onset, fp$offset, fp$peak_age, fp$decline_end, fp$left_cens, fp$right_cens, kb,
                  rep_len(pa, max(k, 1)), rep_len(ps, max(k, 1)))
  }
  out <- as.data.frame(out); out$mc_mode <- "conditioned_k"; out
}
