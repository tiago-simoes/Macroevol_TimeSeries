# ============================================================================
#  mc_uncertainty.R  --  Monte-Carlo uncertainty for first-peak onset/offset
#  Companion to changepoint_pipeline.R (source that first).
#
#  Use when a parameter is reported as ONE series with per-point error bars
#  (e.g. Wright 2017 SI Figs S6/S8/S10), so there is no across-line spread to
#  build _Min/_Med/_Max from (cf. combine_lines()). Each replicate:
#    * Value  ~ split-normal(Value, lower/upper SE)  -- SE from the bar ends
#               (95% CI bars -> half-width/1.96; 1-SE bars -> half-width)
#    * Age    ~ Normal(Age, age_sd)                  -- digitisation/calibration
#               error of the x axis (Ma), order preserved
#  then resample_grid() + first_peak() exactly as in the point estimate.
#  _Min/_Med/_Max = 2.5% / 50% / 97.5% quantiles of the replicate ages.
# ============================================================================

mc_first_peak <- function(df, n_rep = 1000, penalty = "BIC", minseglen = 3,
                          step = 1, age_sd = 0, nonneg = FALSE, seed = 1) {
  set.seed(seed)
  # per-point lower / upper SE (in value units)
  if (all(c("Err_lo", "Err_hi", "Err_type") %in% names(df))) {
    k <- ifelse(df$Err_type == "95CI", 1.96, 1)
    se_lo <- (df$Value - df$Err_lo) / k
    se_hi <- (df$Err_hi - df$Value) / k
  } else {
    se_lo <- se_hi <- df$Value_err
  }
  se_lo[!is.finite(se_lo) | se_lo < 0] <- 0
  se_hi[!is.finite(se_hi) | se_hi < 0] <- 0

  reps <- t(vapply(seq_len(n_rep), function(i) {
    z   <- stats::rnorm(nrow(df))
    v   <- df$Value + ifelse(z > 0, z * se_hi, z * se_lo)
    if (nonneg) v <- pmax(v, 0)
    a   <- sort(df$Age_Ma + stats::rnorm(nrow(df), 0, age_sd), decreasing = TRUE)
    g   <- resample_grid(a, v[order(-df$Age_Ma)], step = step)
    fp  <- tryCatch(first_peak(g, penalty = penalty, minseglen = minseglen),
                    error = function(e) NULL)
    if (is.null(fp)) c(NA, NA, NA) else c(fp$onset, fp$offset, fp$peak_age)
  }, numeric(3)))
  colnames(reps) <- c("onset", "offset", "peak_age")
  as.data.frame(reps)
}

summarise_mc <- function(reps, probs = c(.025, .5, .975)) {
  q <- function(x) round(as.numeric(stats::quantile(x, probs, na.rm = TRUE)), 1)
  o <- q(reps$onset); f <- q(reps$offset); p <- q(reps$peak_age)
  data.frame(Onset_Max = o[3], Onset_Med = o[2], Onset_Min = o[1],
             Offset_Max = f[3], Offset_Med = f[2], Offset_Min = f[1],
             PeakAge_Max = p[3], PeakAge_Med = p[2], PeakAge_Min = p[1],
             MC_n_ok = sum(stats::complete.cases(reps)))
}
