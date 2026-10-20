# ============================================================================
#  run_key_segmented.R -- segmented (continuous piecewise-linear) regression alternative to PELT+CROPS:
#  breakpoints, first-peak onset/offset and uncertainty for every clade x parameter row of a key.
#  Reads the same series as source_code/run_key_crops.R and writes only inside results_segmented(SEG)/.
#
#  Run from the repository root:
#    Rscript source_code/run_key_segmented.R <key.csv> [table_name]
#  Clade origins:
#    Rscript source_code/run_key_segmented.R "ImageDigitization_pipeline/input_keys/key_origin_invert&vert.csv" Table_CladeOrigin_segmented
#  Mass extinctions (automatic when the key has an `extinction_event` column; run the clade origins first):
#    Rscript source_code/run_key_segmented.R "ImageDigitization_pipeline/input_keys/key_massext_invert&vert.csv" Table_MassExtinction_segmented
#  Inputs come from the ImageDigitization_pipeline submodule (DIGITISATION_REPO=<folder> points elsewhere).
#  Environment: N_REP (Monte-Carlo replicates, default 1000); REUSE_PE=1 / REUSE_MC=1 (reuse cached fits / replicates
#               that are newer than the data and were made with the same settings); ONLY=<series id regex> (run a
#               subset; the master table is then not rewritten); OUT_SUB=<folder> (output sub-folder).
#  Outputs (<sub> = clade_origin or mass_extinction):
#    results_segmented(SEG)/<sub>/results/<table>.csv     master table (all columns)
#    results_segmented(SEG)/<sub>/results/clade_tables/    <series>_SEG_{models,segments,MCreps}.csv, _SEG_pe.rds
#    results_segmented(SEG)/<sub>/figures/<series>_SEG.png BIC vs k; data, fitted segments, onset/offset with MC ranges
#  Mass-extinction mode: series = outputs/data_tables/mass_extinction/<series id>__<event slug>.csv; extra columns
#    Extinction_event and Extinction_age_Ma; clades in input_keys/massext_copy_from_origin.csv copied from the
#    clade-origin table; series with < 4 points after the extinction not analysed.
#  Optional overrides in data_digitized/digitisation/ (mass extinctions: *_massext.csv with an Extinction_event column):
#    model_order_overrides_segmented.csv (Clade, Parameter, k, Reason) and
#    curve_type_overrides.csv (Clade, Parameter, Curve_type = hump|cumulative, Reason);
#    Parameter may carry the data_type slug, e.g. taxonomic_diversity_BDD.
# ============================================================================
if (!isTRUE(l10n_info()[["UTF-8"]]))
  for (lc in c("C.UTF-8", "en_US.UTF-8", "English_United States.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)))) break
SD <- "results_segmented(SEG)"; DATA_DIR <- "data_digitized"   # analysis settings (optional overrides)
DIGI <- Sys.getenv("DIGITISATION_REPO", "ImageDigitization_pipeline")   # submodule: keys, digitised series, copy list
if (!dir.exists(file.path(DIGI, "outputs", "data_tables")))
  stop("ImageDigitization_pipeline not found in '", DIGI, "': run `git submodule update --init` or set DIGITISATION_REPO")
source(file.path("source_code", "seg_firstpeak.R"))
args <- commandArgs(trailingOnly = TRUE)
key_file   <- if (length(args) >= 1) args[1] else file.path(DIGI, "input_keys", "key_origin_invert&vert.csv")
tag        <- if (length(args) >= 2) args[2] else "Table_CladeOrigin_segmented"
N_REP <- as.integer(Sys.getenv("N_REP", "1000")); ONLY <- Sys.getenv("ONLY", "")
STAMP <- paste(unlist(SEG_PAR), collapse = "/")            # settings stamp stored with the MC replicates

rd <- function(f) {
  x <- readLines(f, warn = FALSE, encoding = "bytes")
  x <- if (all(validUTF8(x))) { Encoding(x) <- "UTF-8"; x } else iconv(x, "latin1", "UTF-8")   # UTF-8 files with accents: re-mark
  x <- sub("^﻿", "", x, useBytes = TRUE)
  utils::read.csv(text = x, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}
read_series <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"Value_err" %in% names(df)) df$Value_err <- 0
  df <- df[stats::complete.cases(df$Age_Ma, df$Value), ]
  df[order(-df$Age_Ma), ]
}
key <- rd(key_file); key <- key[!is.na(key$clade) & nzchar(trimws(key$clade)), ]
for (cc in c("group", "legend/line", "data_file", "data_type")) if (!cc %in% names(key)) key[[cc]] <- ""
key[is.na(key)] <- ""
MASS <- "extinction_event" %in% names(key)                     # mass-extinction key
SUB <- Sys.getenv("OUT_SUB", if (MASS) "mass_extinction" else "clade_origin")
RES <- file.path(SD, SUB, "results"); CT <- file.path(RES, "clade_tables"); FIG <- file.path(SD, SUB, "figures")
for (d in c(CT, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
SERIES_DIR <- file.path(DIGI, "outputs", "data_tables", if (MASS) "mass_extinction" else "clade_origin")
if (MASS) {
  key$extinction_event <- trimws(key$extinction_event)
  key$ext_age <- as.numeric(key[[grep("^extinction_age", names(key), value = TRUE)[1]]])
  cpy <- rd(file.path(DIGI, "input_keys", "massext_copy_from_origin.csv"))
  origin_tab <- rd(file.path(SD, "clade_origin", "results", "Table_CladeOrigin_segmented.csv"))
} else key$extinction_event <- ""
slug <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", trimws(x)))
std_param <- function(x) { x <- tolower(trimws(x))
  if (grepl("rate", x)) "phenotypic_rate" else if (grepl("dispar", x)) "phenotypic_disparity"
  else if (grepl("divers", x)) "taxonomic_diversity" else gsub(" +", "_", x) }
sfx <- if (MASS) "_massext" else ""
ovr_k  <- if (file.exists(f <- file.path(DATA_DIR, "digitisation", paste0("model_order_overrides_segmented", sfx, ".csv")))) rd(f) else NULL
ovr_ct <- if (file.exists(f <- file.path(DATA_DIR, "digitisation", paste0("curve_type_overrides", sfx, ".csv")))) rd(f) else NULL
evok <- function(o, ev) if (MASS && !is.null(o) && "Extinction_event" %in% names(o)) trimws(o$Extinction_event) == ev else TRUE

key$pr_std <- vapply(key$parameter, std_param, "")
key$dt_slug <- slug(key$data_type)
sid0 <- paste(key$clade, key$pr_std, key$dt_slug, key$extinction_event)   # within an event
dup <- duplicated(sid0) | duplicated(sid0, fromLast = TRUE)
nodt <- dup & !nzchar(key$dt_slug) & !(tolower(trimws(key$reference)) %in% c("na", ""))
key$dt_slug[nodt] <- slug(key$reference[nodt]); key$data_type[nodt] <- key$reference[nodt]

r1 <- function(x) round(x, 1)
fmtq <- function(x) { q <- r1(as.numeric(stats::quantile(x, c(.975, .5, .025), na.rm = TRUE))); if (all(is.na(x))) rep(NA, 3) else q }
out <- list()
for (i in seq_len(nrow(key))) {
  cl <- key$clade[i]; pr <- key$pr_std[i]; dt <- key$dt_slug[i]
  sid <- paste0(cl, "_", pr, if (nzchar(dt)) paste0("_", dt)); pr_dt <- sub(paste0("^", cl, "_"), "", sid)
  ev <- key$extinction_event[i]; sidf <- if (MASS) paste0(sid, "__", slug(ev)) else sid   # file / figure id
  if (nzchar(ONLY) && !grepl(ONLY, sidf)) next
  note <- character(0)
  fig <- key[["figure/table"]][i]; dfile <- trimws(key$data_file[i])
  base <- data.frame(Group = key$group[i], Clade = cl, Parameter = pr, Data_type = key$data_type[i],
                     Reference = key$reference[i], Figure = fig, Legend_line = key[["legend/line"]][i],
                     Key_data_file = dfile, stringsAsFactors = FALSE)
  if (MASS) base <- cbind(base[, 1:2], Extinction_event = ev, Extinction_age_Ma = key$ext_age[i], base[, -(1:2)])
  f <- file.path(SERIES_DIR, paste0(sidf, ".csv"))
  if (MASS && cl %in% cpy$Clade && !(tolower(trimws(key$reference[i])) %in% c("na", ""))) {   # copied from the clade-origin table
    o <- origin_tab[origin_tab$Clade == cl & origin_tab$Parameter == pr & origin_tab$Reference == key$reference[i], , drop = FALSE]
    if (nrow(o) != 1) stop("clade-origin row not found once: ", cl, " ", pr, " ", key$reference[i])
    o <- cbind(o[, 1:2], Extinction_event = ev, Extinction_age_Ma = key$ext_age[i], o[, -(1:2)])
    o$Flags <- paste(c(sprintf("Copied from the clade-origin analysis (results_segmented(SEG)/clade_origin): %s", cpy$Reason[cpy$Clade == cl][1]),
                       if (!is.na(o$Flags) && nzchar(o$Flags)) o$Flags), collapse = " | ")
    out[[i]] <- o; message(sprintf("%-24s %-38s copied from clade origin", cl, pr_dt)); next
  }
  if (tolower(trimws(key$reference[i])) %in% c("na", "") || !file.exists(f)) {
    msg <- if (tolower(trimws(key$reference[i])) %in% c("na", "")) "No published source for this parameter (key: na)" else paste("Data file not found:", f)
    out[[i]] <- cbind(base, Method = NA, Flags = paste(c(msg, note), collapse = " | ")); next
  }
  df <- read_series(f)
  if (MASS && nrow(df) < 4) {
    out[[i]] <- cbind(base, Data_source = df$Source[1], N_points = nrow(df), Method = NA,
                      Flags = paste(c(sprintf("Too few data points after the extinction (n = %d, %s Ma): not analysed", nrow(df),
                                              paste(round(df$Age_Ma, 1), collapse = ", ")), note), collapse = " | "))
    message(sprintf("%-24s %-38s %s: only %d points after the extinction - skipped", cl, pr_dt, ev, nrow(df))); next
  }
  oldest <- max(df$Age_Ma); youngest <- min(df$Age_Ma); span <- oldest - youngest
  ct <- curve_type(df$Age_Ma, df$Value, pr); curve <- ct$type; ct_note <- NULL
  oc <- if (!is.null(ovr_ct)) ovr_ct[ovr_ct$Clade == cl & ovr_ct$Parameter == pr_dt & evok(ovr_ct, ev), , drop = FALSE] else NULL
  if (!is.null(oc) && nrow(oc)) { curve <- trimws(oc$Curve_type[1]); ct_note <- sprintf("Curve type set to '%s' (override): %s", curve, oc$Reason[1]) }
  ok <- if (!is.null(ovr_k)) ovr_k[ovr_k$Clade == cl & ovr_k$Parameter == pr_dt & evok(ovr_k, ev), , drop = FALSE] else NULL
  k_fix <- if (!is.null(ok) && nrow(ok)) as.integer(ok$k[1]) else NA
  pe_file <- file.path(CT, paste0(sidf, "_SEG_pe.rds"))           # cached point estimate (REUSE_PE=1)
  pe <- NULL
  if (Sys.getenv("REUSE_PE") == "1" && file.exists(pe_file) && file.mtime(pe_file) > file.mtime(f)) {
    pe <- readRDS(pe_file); if (!identical(pe$stamp, paste(STAMP, curve, k_fix))) pe <- NULL }
  if (is.null(pe)) {
    set.seed(42)
    pe <- seg_analyse(df$Age_Ma, df$Value, curve, k_fix); pe$stamp <- paste(STAMP, curve, k_fix)
    saveRDS(pe, pe_file)
  }
  fp <- pe$fp; inf <- pe$inf; k <- pe$k
  # ---- tables: model selection + segments
  sel <- pe$sel$table
  sel$breakpoints_Ma <- vapply(pe$sel$fits, function(z) paste(r1(pe$t0 - z$psi), collapse = ";"), "")
  sel$selected <- sel$k == pe$k_bic
  utils::write.csv(sel, file.path(CT, paste0(sidf, "_SEG_models.csv")), row.names = FALSE)
  b <- c(pe$t0, pe$psi_age, youngest)
  seg <- data.frame(segment = seq_len(k + 1), start_Ma = r1(b[-length(b)]), end_Ma = r1(b[-1]),
                    slope_per_Myr = signif(inf$slope, 4), slope_SE = signif(inf$slope_se, 3),
                    change_frac_of_range = round(inf$slope * diff(-b) / diff(range(inf$fitted)), 3),
                    class = fp$seg_class,
                    breakpoint_after_Ma = c(r1(pe$psi_age), NA), breakpoint_SE = c(signif(inf$psi_se, 3), NA),
                    slope_change_class = c(fp$bp_class, rep("", 1))[seq_len(k + 1)])
  utils::write.csv(seg, file.path(CT, paste0(sidf, "_SEG_segments.csv")), row.names = FALSE)
  # ---- Monte Carlo (conditioned on k)
  mc_file <- file.path(CT, paste0(sidf, "_SEG_MCreps.csv"))
  reuse <- Sys.getenv("REUSE_MC") == "1" && file.exists(mc_file) && file.mtime(mc_file) > file.mtime(f)
  MSTAMP <- paste(STAMP, MC_VERSION)
  if (reuse) { reps <- utils::read.csv(mc_file); reuse <- nrow(reps) == N_REP && "stamp" %in% names(reps) && reps$stamp[1] == MSTAMP &&
                 reps$k_pe[1] == k && reps$curve[1] == curve }
  if (!reuse) {
    age_sd <- if ("Age_cal_sd" %in% names(df)) df$Age_cal_sd[1] else 0
    reps <- if (N_REP > 0) mc_seg(df, pe, N_REP, age_sd, nonneg = pr == "taxonomic_diversity") else
      data.frame(onset = NA, offset = NA, peak_age = NA, decline_end = NA, left_cens = NA, right_cens = NA, k_bic = NA, mc_mode = "none")
    reps$stamp <- MSTAMP; reps$k_pe <- k; reps$curve <- curve
    utils::write.csv(reps, mc_file, row.names = FALSE)
  }
  qo <- fmtq(reps$onset); qf <- fmtq(reps$offset); qp <- fmtq(reps$peak_age)
  kb <- table(reps$k_bic); kb_txt <- if (length(kb)) paste(sprintf("%s:%.0f%%", names(kb), 100 * kb / sum(kb)), collapse = ";") else ""
  lc_pct <- round(100 * mean(reps$left_cens == 1, na.rm = TRUE)); rc_pct <- round(100 * mean(reps$right_cens == 1, na.rm = TRUE))
  # ---- model-based 95% CI of the onset/offset breakpoints (delta method, effective-n scaled)
  ci <- function(j) if (is.na(j)) c(NA, NA) else r1(pe$psi_age[j] + c(1, -1) * inf$tcrit * inf$psi_se[j])
  oci <- ci(fp$onset_bp); fci <- ci(fp$offset_bp)
  bp_ci <- if (k) paste(sprintf("%.1f [%.1f-%.1f]", pe$psi_age, pe$psi_age + inf$tcrit * inf$psi_se, pe$psi_age - inf$tcrit * inf$psi_se), collapse = "; ") else ""
  # robustness of the point estimate to data error: PE outside the MC interquartile range
  tolr <- max(0.5, 0.01 * span)
  nonrob <- function(pe_v, x) { q <- stats::quantile(x, c(.25, .5, .75), na.rm = TRUE)
    if (N_REP > 0 && all(is.finite(q)) && (pe_v > q[3] + tolr || pe_v < q[1] - tolr)) q[2] else NA }
  nr_on <- nonrob(fp$onset, reps$onset); nr_off <- nonrob(fp$offset, reps$offset)
  et <- if ("Err_type" %in% names(df)) df$Err_type[1] else "unknown"
  no_rise <- curve == "hump" && !any(fp$seg_class == "R") && !fp$start_decline
  flag <- c(
    if (curve == "cumulative") sprintf("Cumulative (monotonic) diversity curve: slope-change rule (onset = first acceleration, offset = first deceleration); drawdown %.0f%% of range", 100 * ct$drawdown),
    ct_note,
    if (k == 0) "No breakpoint (BIC k = 0): single linear trend",
    if (curve == "hump" && fp$start_decline) "Series starts in decline: first peak at or before the oldest data point (onset = offset = series start)",
    if (no_rise) "No significant rising segment: onset/offset = series limits",
    if (curve == "cumulative" && fp$left_cens && !is.na(fp$offset_bp)) "Curve is steepest at its start: onset = oldest data point (left-censored)",
    if (!fp$start_decline && !no_rise && fp$left_cens && curve == "hump") "Onset = oldest data point (left-censored: series starts rising or at peak level)",
    if (fp$right_cens && !no_rise) "Offset = youngest data point (right-censored: no decline after the first peak)",
    if (curve == "cumulative" && length(fp$later_acc)) sprintf("Later acceleration(s) after the first burst at %s Ma", paste(r1(fp$later_acc), collapse = ", ")),
    if (curve == "hump" && !is.na(fp$main_rise) && abs(fp$main_rise - fp$onset) > 0.05)
      sprintf("Rise steepens at %.1f Ma (start of the steepest rising segment)", fp$main_rise),
    if (fp$n_steps > 1) sprintf("Stepwise %s: %d separate %s before the first peak; onset = the first",
                                if (curve == "hump") "rise" else "acceleration", fp$n_steps, if (curve == "hump") "rising phases" else "accelerations"),
    if (!is.na(k_fix)) sprintf("Number of breakpoints fixed at %d (BIC k = %d): %s", k, pe$k_bic, ok$Reason[1]),
    if (N_REP > 0 && is.finite(qo[1] - qo[3]) && (qo[1] - qo[3]) > 0.5 * span) "Onset 95% MC range > 50% of series span: weakly constrained",
    if (N_REP > 0 && length(kb) && max(kb) / sum(kb) < 0.6) "BIC model order unstable across MC replicates (<60% of 100 replicates agree; MC conditioned on the point-estimate k)",
    if (!is.na(nr_on)) sprintf("Onset point estimate not robust to data error: outside the MC interquartile range (MC median %.1f Ma)", nr_on),
    if (!is.na(nr_off)) sprintf("Offset point estimate not robust to data error: outside the MC interquartile range (MC median %.1f Ma)", nr_off),
    if (N_REP > 0 && lc_pct > 5 && lc_pct < 95) sprintf("MC: onset at the series start in %d%% of replicates", lc_pct),
    if (N_REP > 0 && rc_pct > 5 && rc_pct < 95) sprintf("MC: offset at the series end in %d%% of replicates", rc_pct),
    if (et == "none") "No error bars in source: MC reflects age jitter only; use model CI / breakpoint SE for timing uncertainty",
    if (et == "Poisson") "Value error assumed Poisson (sqrt(n)) for counts without published error bars",
    if (MASS && nrow(df) < 8) sprintf("Few data points after the extinction (n = %d): breakpoints poorly constrained", nrow(df)),
    note)
  meth <- sprintf("Segmented regression (continuous piecewise-linear, Muggeo); k = %d by BIC%s%s", k,
                  if (pe$binned) sprintf(" on %.1f-Myr bins (n = %d)", pe$bins$w, pe$n_eff) else sprintf(" on data points (n = %d)", pe$n_eff),
                  if (!is.na(k_fix)) sprintf(" [fixed; BIC k = %d]", pe$k_bic) else "")
  origin_bp <- pe$psi_age[pe$psi_age >= fp$offset - 0.01]
  out[[i]] <- cbind(base, data.frame(
    Data_source = df$Source[1], Error_type = et, N_points = nrow(df),
    Series_start_Ma = r1(oldest), Series_end_Ma = r1(youngest), Method = meth,
    Curve_type = curve, Rule = if (curve == "hump") "first rise -> first decline (PELT-style)" else "first acceleration -> first deceleration",
    BIC_k = k, Kmax = pe$sel$Kmax, Breakpoints_Ma = paste(r1(pe$psi_age), collapse = ";"), Breakpoints_CI95 = bp_ci,
    Origin_breakpoints_Ma = paste(r1(origin_bp), collapse = ";"),
    Segment_class = paste(fp$seg_class, collapse = ""), Slope_change_class = paste(fp$bp_class, collapse = ","),
    Segment_slopes = paste(signif(inf$slope, 3), collapse = ";"),
    Onset_PE = r1(fp$onset), Onset_CI95_older = oci[1], Onset_CI95_younger = oci[2],
    Offset_PE = r1(fp$offset), Offset_CI95_older = fci[1], Offset_CI95_younger = fci[2],
    Peak_Age_PE = r1(fp$peak_age), Main_rise_PE = r1(fp$main_rise), Decline_end_PE = r1(fp$decline_end),
    BIC_resolution = if (pe$binned) sprintf("%.1f-Myr bins (n = %d)", pe$bins$w, pe$n_eff) else sprintf("data points (n = %d)", pe$n_eff),
    Onset_Max = qo[1], Onset_Med = qo[2], Onset_Min = qo[3], Offset_Max = qf[1], Offset_Med = qf[2], Offset_Min = qf[3],
    PeakAge_Max = qp[1], PeakAge_Med = qp[2], PeakAge_Min = qp[3],
    MC_n_ok = sum(!is.na(reps$onset)), MC_k_bic_freq = kb_txt, MC_left_cens_pct = lc_pct, MC_right_cens_pct = rc_pct,
    Flags = paste(flag, collapse = " | "), stringsAsFactors = FALSE))

  # ---- figure
  png(file.path(FIG, paste0(sidf, "_SEG.png")), width = 1500, height = 620, res = 130)
  layout(matrix(1:2, 1), widths = c(1, 2.1)); par(mar = c(4.2, 4.2, 3, 1))
  plot(sel$k, sel$dBIC, type = "b", pch = 16, xlab = "Number of breakpoints (k)", ylab = expression(Delta * "BIC"),
       main = if (pe$binned) sprintf("BIC (%.1f-Myr bins, n = %d)", pe$bins$w, pe$n_eff) else sprintf("BIC (n = %d)", pe$n_eff))
  points(k, sel$dBIC[sel$k == k], pch = 21, cex = 2.4, col = "red", lwd = 2)
  lo <- df$Err_lo; hi <- df$Err_hi; if (is.null(lo)) { lo <- df$Value - df$Value_err; hi <- df$Value + df$Value_err }
  plot(df$Age_Ma, df$Value, pch = 16, cex = if (nrow(df) > 100) .4 else .9, col = "grey25",
       xlim = rev(range(c(df$Age_Ma, qo, qf), na.rm = TRUE)), ylim = range(c(lo, hi, inf$fitted), na.rm = TRUE),
       xlab = "Age (Ma)", ylab = pr_dt, main = sprintf("%s - %s%s: segmented, k = %d (%s)", cl, pr_dt,
                                                       if (MASS) paste0(" after ", ev, " (", key$ext_age[i], " Ma)") else "", k, curve))
  usr <- par("usr")
  if (N_REP > 0) { rect(qo[1], usr[3], qo[3], usr[4], col = adjustcolor("forestgreen", .18), border = NA)
                   rect(qf[1], usr[3], qf[3], usr[4], col = adjustcolor("purple", .18), border = NA) }
  if (nrow(df) <= 100) { ea <- which(hi - lo > 0); if (length(ea)) arrows(df$Age_Ma[ea], lo[ea], df$Age_Ma[ea], hi[ea], angle = 90, code = 3, length = .02, col = "grey55") }
  else polygon(c(df$Age_Ma, rev(df$Age_Ma)), c(lo, rev(hi)), col = adjustcolor("grey50", .25), border = NA)
  if (pe$binned) points(pe$t0 - pe$bins$t, pe$bins$y, pch = 0, cex = .9, col = "#e67e22")
  cols <- c(R = "#c0392b", F = "steelblue", D = "#2c3e50")
  xs <- c(pe$t0, pe$psi_age, youngest); tt <- pe$t0 - xs
  ys <- as.numeric(cbind(1, tt, hinge(tt, pe$psi)) %*% inf$coef)
  for (j in seq_len(k + 1)) segments(xs[j], ys[j], xs[j + 1], ys[j + 1], col = cols[fp$seg_class[j]], lwd = 3)
  if (k) { abline(v = pe$psi_age, col = "grey30", lty = 3)
    yb <- usr[3] + 0.04 * diff(usr[3:4])
    segments(pe$psi_age + inf$tcrit * inf$psi_se, yb, pe$psi_age - inf$tcrit * inf$psi_se, yb, col = "grey30", lwd = 2) }
  abline(v = c(fp$onset, fp$offset), col = c("forestgreen", "purple"), lty = 2, lwd = 2)
  legend("topright", bty = "n", cex = .72,
    legend = c(sprintf("onset %.1f Ma [95%% MC %.1f-%.1f]", fp$onset, qo[1], qo[3]),
               sprintf("offset %.1f Ma [95%% MC %.1f-%.1f]", fp$offset, qf[1], qf[3]),
               "fit: rising / flat / declining segment", "breakpoints (bar = model 95% CI)",
               if (pe$binned) sprintf("%.1f-Myr bin means (used for BIC)", pe$bins$w)),
    col = c("forestgreen", "purple", "#c0392b", "grey30", if (pe$binned) "#e67e22"),
    lty = c(2, 2, 1, 3, if (pe$binned) NA), pch = c(NA, NA, NA, NA, if (pe$binned) 0), lwd = 2)
  dev.off()
  message(sprintf("%-24s %-38s %s %s k=%d onset %6.1f [%6.1f-%6.1f]  offset %6.1f [%6.1f-%6.1f] %s",
                  cl, pr_dt, ev, substr(curve, 1, 3), k, fp$onset, qo[1], qo[3], fp$offset, qf[1], qf[3], paste(fp$seg_class, collapse = "")))
}
if (!nzchar(ONLY)) {
  cols <- unique(unlist(lapply(out, names)))
  res <- do.call(rbind, lapply(out[!vapply(out, is.null, TRUE)], function(d) { for (cc in setdiff(cols, names(d))) d[[cc]] <- NA; d[, cols] }))
  master <- file.path(RES, paste0(tag, ".csv"))
  utils::write.csv(res, master, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
