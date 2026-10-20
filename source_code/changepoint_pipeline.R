# ============================================================================
#  changepoint_pipeline.R
#  Detection of shifts in macroevolutionary time series with PELT
#  (Pruned Exact Linear Time), via the `changepoint` package.
#
#  Input : tidy CSV of digitized points, one row per point, columns:
#            Clade, Parameter, Age_Ma, Value, Value_err, Source
#  Output: changepoint ages (onset / offset), Monte-Carlo uncertainty,
#          and a diagnostic plot per series.
#
#  Requires: changepoint (>= 2.2; tested with 2.3), optionally changepoint.np for
#            the nonparametric PELT option.
# ============================================================================

## ---- dependencies ---------------------------------------------------------
# install.packages("changepoint")          # CRAN
# install.packages("changepoint.np")        # optional, nonparametric PELT
suppressPackageStartupMessages({
  library(changepoint)
})

# ---------------------------------------------------------------------------
# 1. Read + validate a digitized series
# ---------------------------------------------------------------------------
read_series <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("Age_Ma", "Value")
  if (!all(need %in% names(df)))
    stop("CSV must contain at least columns: ", paste(need, collapse = ", "))
  if (!"Value_err" %in% names(df)) df$Value_err <- 0        # no uncertainty given
  df <- df[stats::complete.cases(df$Age_Ma, df$Value), ]
  df <- df[order(-df$Age_Ma), ]                             # old -> young
  df
}

# ---------------------------------------------------------------------------
# 2. Resample to a regular time grid.
#    PELT operates on the *sequence* of values and treats observations as
#    evenly spaced. Digitized points are irregularly spaced in Ma, so we must
#    interpolate onto a constant-step grid BEFORE running PELT, then map the
#    returned changepoint indices back to ages. `step` is in Myr.
#    The returned vectors run time-forward (i.e. Age_Ma DECREASING).
# ---------------------------------------------------------------------------
resample_grid <- function(age, value, step = 1) {
  o     <- order(age)                                       # ascending age
  a     <- age[o]; v <- value[o]
  grid  <- seq(min(a), max(a), by = step)                   # ascending age
  gv    <- stats::approx(x = a, y = v, xout = grid)$y       # linear interp
  list(age = rev(grid), value = rev(gv))                    # -> time-forward
}

# ---------------------------------------------------------------------------
# 3. Core PELT call.
#    mode = "level"   -> cpt.mean    : shifts in the MEAN of the parameter
#                                       (use for disparity / rate levels)
#    mode = "meanvar" -> cpt.meanvar : shifts in mean AND variance
#                                       (regime shifts, e.g. a radiation)
#    mode = "slope"   -> cpt.mean on the first difference of the series
#                                       (use for cumulative diversity /
#                                        lineage-count curves, where the
#                                        biologically meaningful shift is a
#                                        change in the *rate of accumulation*,
#                                        i.e. the slope)
#
#    penalty: "MBIC" (default, conservative), "BIC"/"SIC", "AIC",
#             "Hannan-Quinn", "Manual" (supply pen.value), or "CROPS"
#             (supply pen.value = c(min,max) to explore a range).
#    minseglen: minimum number of grid steps between changepoints.
# ---------------------------------------------------------------------------
run_pelt <- function(grid, mode = "level",
                     penalty = "MBIC", pen.value = 0, minseglen = 3) {

  sig     <- grid$value
  cpt_age <- grid$age
  if (mode == "slope") {                                    # analyse the slope
    sig     <- diff(grid$value)
    cpt_age <- grid$age[-1]
  }

  fit <- switch(
    mode,
    "level"   = cpt.mean(sig, method = "PELT", penalty = penalty,
                         pen.value = pen.value, minseglen = minseglen),
    "slope"   = cpt.mean(sig, method = "PELT", penalty = penalty,
                         pen.value = pen.value, minseglen = minseglen),
    "meanvar" = cpt.meanvar(sig, test.stat = "Normal", method = "PELT",
                            penalty = penalty, pen.value = pen.value,
                            minseglen = minseglen),
    stop("mode must be one of 'level', 'meanvar', 'slope'")
  )

  idx  <- cpts(fit)                                         # segment-end indices
  ages <- sort(cpt_age[idx], decreasing = TRUE)             # oldest -> youngest
  list(fit = fit, ages = ages, cpt_age = cpt_age, signal = sig)
}

# ---------------------------------------------------------------------------
# 3b. FIRST-PEAK onset / offset  (the standard shift-timing output).
#
#    Onset and offset are the RISING and FALLING edges of the SAME first peak
#    -- not two independent events. PELT (cpt.mean) segments the series into
#    constant-mean levels; we take the first segment that is a local maximum in
#    segment-mean space (>= the segment before AND after it). Its older edge is
#    the onset (rates increase into the peak) and its younger edge is the offset
#    (rates fall back off the peak). If the series already starts at its peak
#    level (no lower segment before it), onset = the oldest age (clade origin).
#
#    Penalty: BIC is the standard choice used throughout. Under BIC an initial
#    "shoulder" of moderately elevated rate is absorbed into the first peak
#    (rather than split off), which is the intended behaviour here. The signal
#    is standardized first so the penalty scale is comparable across series.
#
#    Returns onset, offset, peak age/value, and the segment means.
# ---------------------------------------------------------------------------
first_peak <- function(grid, penalty = "BIC", pen.value = 0, minseglen = 3) {
  v <- grid$value; age <- grid$age                         # age is time-forward (decreasing)
  z <- (v - mean(v)) / stats::sd(v)                        # standardize
  fit <- cpt.mean(z, method = "PELT", penalty = penalty,
                  pen.value = pen.value, minseglen = minseglen)
  bnds   <- c(0, cpts(fit), length(z))
  starts <- utils::head(bnds, -1) + 1
  ends   <- bnds[-1]
  seg_mean <- mapply(function(s, e) mean(v[s:e]), starts, ends)

  is_peak <- function(i) {
    left  <- if (i > 1)               seg_mean[i] >= seg_mean[i - 1] else TRUE
    right <- if (i < length(seg_mean)) seg_mean[i] >= seg_mean[i + 1] else TRUE
    left && right
  }
  pk <- which(vapply(seq_along(seg_mean), is_peak, logical(1)))[1]
  s <- starts[pk]; e <- ends[pk]

  onset  <- if (pk == 1)               age[s] else (age[ends[pk - 1]] + age[s]) / 2
  offset <- if (pk == length(seg_mean)) age[e] else (age[e] + age[starts[pk + 1]]) / 2
  peak_i <- s - 1 + which.max(v[s:e])

  list(onset = round(onset, 1), offset = round(offset, 1),
       peak_age = round(age[peak_i], 1), peak_val = round(v[peak_i], 3),
       seg_means = round(seg_mean, 3), starts = starts, ends = ends, age = age, fit = fit)
}

# ---------------------------------------------------------------------------
# 4. CROPS penalty diagnostic.
#    Runs PELT across a range of penalties and reports how many changepoints
#    (and where) are returned at each. Use this to choose a defensible
#    penalty / check the stability of a solution (an "elbow" in the number of
#    changepoints vs penalty is the usual target). Returns a data.frame.
# ---------------------------------------------------------------------------
crops_scan <- function(grid, mode = "level",
                       pen.range = c(0.1, 30), minseglen = 3) {
  sig     <- if (mode == "slope") diff(grid$value) else grid$value
  cpt_age <- if (mode == "slope") grid$age[-1] else grid$age
  fit <- if (mode == "meanvar")
           cpt.meanvar(sig, test.stat = "Normal", method = "PELT",
                       penalty = "CROPS", pen.value = pen.range,
                       minseglen = minseglen)
         else
           cpt.mean(sig, method = "PELT", penalty = "CROPS",
                    pen.value = pen.range, minseglen = minseglen)
  full <- cpts.full(fit); pens <- pen.value.full(fit)
  out <- lapply(seq_len(nrow(full)), function(i) {
    ix <- full[i, ]; ix <- ix[!is.na(ix) & ix > 0]
    data.frame(penalty  = round(pens[i], 3),
               n_cpts   = length(ix),
               ages_Ma  = paste(round(sort(cpt_age[ix], decreasing = TRUE), 1),
                                collapse = ";"),
               stringsAsFactors = FALSE)
  })
  unique(do.call(rbind, out))
}

# ---------------------------------------------------------------------------
# 5. Analyse ONE digitized line -> first-peak onset & offset.
#    With a single line there is no across-line spread, so the reported
#    _Max / _Med / _Min collapse to the same value (per your convention:
#    the min/med/max only carry meaning when several lines are combined).
# ---------------------------------------------------------------------------
analyse_series <- function(path, penalty = "BIC", pen.value = 0,
                           minseglen = 3, step = 1, plot_file = NULL) {
  df <- read_series(path)
  g  <- resample_grid(df$Age_Ma, df$Value, step = step)
  fp <- first_peak(g, penalty = penalty, pen.value = pen.value, minseglen = minseglen)

  if (!is.null(plot_file)) {
    grDevices::png(plot_file, width = 1000, height = 620, res = 130)
    plot(df$Age_Ma, df$Value, pch = 16, col = "grey60", xlim = rev(range(df$Age_Ma)),
         xlab = "Age (Ma)", ylab = df$Parameter[1],
         main = paste(df$Clade[1], "-", df$Parameter[1], "(first peak)"))
    graphics::rect(fp$offset, graphics::par("usr")[3], fp$onset, graphics::par("usr")[4],
                   col = grDevices::adjustcolor("gold", 0.18), border = NA)
    graphics::lines(g$age, g$value, col = "#c0392b")
    graphics::abline(v = c(fp$onset, fp$offset), col = c("forestgreen", "purple"),
                     lty = 2, lwd = 2)
    grDevices::dev.off()
  }

  data.frame(
    Clade      = df$Clade[1]     %||% NA,
    Parameter  = df$Parameter[1] %||% NA,
    Penalty    = penalty,
    Peak_Age_Ma = fp$peak_age, Peak_Value = fp$peak_val,
    Onset_Max = fp$onset, Onset_Med = fp$onset, Onset_Min = fp$onset,
    Offset_Max = fp$offset, Offset_Med = fp$offset, Offset_Min = fp$offset,
    Source     = if ("Source" %in% names(df)) df$Source[1] else NA,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 6. Combine SEVERAL lines of the same figure (e.g. the six root-age curves of
#    Goswami et al. 2022 Fig.3A) into one onset/offset estimate WITH
#    uncertainty. Each line contributes one first-peak onset and one offset;
#    the across-line min / median / max become the _Min / _Med / _Max. This is
#    where the uncertainty columns acquire meaning.
#
#    `paths` : character vector of per-line digitized CSVs.
#    Lines whose first-peak cannot be resolved (too few points, or no clear
#    local-max segment) should be dropped by the caller; `min_points` guards
#    the most common failure (heavily occluded lines).
# ---------------------------------------------------------------------------
combine_lines <- function(paths, penalty = "BIC", minseglen = 3, step = 1,
                          min_points = 60, max_width_frac = 0.5,
                          clade = NULL, parameter = NULL) {
  per <- lapply(paths, function(p) {
    df <- read_series(p)
    span <- diff(range(df$Age_Ma))
    if (nrow(df) < min_points)                              # too occluded to read
      return(data.frame(line = basename(p), n = nrow(df), reliable = FALSE,
                        onset = NA, offset = NA, peak = NA,
                        note = "too few points", stringsAsFactors = FALSE))
    g  <- resample_grid(df$Age_Ma, df$Value, step = step)
    fp <- first_peak(g, penalty = penalty, minseglen = minseglen)
    width <- fp$onset - fp$offset
    ok    <- is.finite(width) && width <= max_width_frac * span   # a real peak, not a broad plateau
    data.frame(line = basename(p), n = nrow(df), reliable = ok,
               onset = fp$onset, offset = fp$offset, peak = fp$peak_age,
               note = if (ok) "" else "no resolvable peak (broad)",
               stringsAsFactors = FALSE)
  })
  per <- do.call(rbind, per)
  rel <- per[per$reliable & !is.na(per$onset), ]
  q <- function(x, p) round(as.numeric(stats::quantile(x, p, names = FALSE)), 1)
  combined <- data.frame(
    Clade      = clade     %||% NA,
    Parameter  = parameter %||% NA,
    N_lines    = nrow(rel),
    Lines_used = paste(rel$line, collapse = ";"),
    Onset_Max = max(rel$onset), Onset_Med = q(rel$onset, .5), Onset_Min = min(rel$onset),
    Offset_Max = max(rel$offset), Offset_Med = q(rel$offset, .5), Offset_Min = min(rel$offset),
    stringsAsFactors = FALSE
  )
  list(per_line = per, combined = combined)
}

# tiny null-coalescing helper
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
