# ============================================================================
#  run_key_crops.R -- PELT + CROPS shift points (kneedle-elbow model order), first-peak onset/offset and
#  Monte-Carlo uncertainty for every clade x parameter row of a key.
#
#  Run from the repository root:
#    Rscript source_code/run_key_crops.R <key.csv> [table_name] [extra_series.csv]
#  Clade origins:
#    Rscript source_code/run_key_crops.R "ImageDigitization_pipeline/input_keys/key_origin_invert&vert.csv" Table_CladeOrigin
#  Mass extinctions (automatic when the key has an `extinction_event` column; run the clade origins first):
#    Rscript source_code/run_key_crops.R "ImageDigitization_pipeline/input_keys/key_massext_invert&vert.csv" Table_MassExtinction
#  Inputs come from the ImageDigitization_pipeline submodule (DIGITISATION_REPO=<folder> points elsewhere): keys in
#  input_keys/, digitised series in outputs/data_tables/clade_origin/ and mass_extinction/.
#  Outputs (<sub> = clade_origin or mass_extinction; OUT_SUB=<folder> overrides it):
#    results_PELT+CROP/<sub>/results/<table_name>.csv     master table (all columns)
#    results_PELT+CROP/<sub>/results/clade_tables/         per-series CROPS solutions and MC replicates
#    results_PELT+CROP/<sub>/figures/<series>_CROPS.png    CROPS elbow; data, segment means, onset/offset with MC ranges
#  Environment: REUSE_MC=1 reuses saved MC replicates that are newer than the data file.
#
#  key.csv      : group, clade, parameter, reference, figure/table [, legend/line, note, data_file, data_type].
#                 Rows with reference "na" are reported as "no data source". Several rows for one clade x parameter
#                 are kept apart by data_type (series <clade>_<parameter>_<data_type slug>.csv; slug = non-alphanumerics
#                 -> "_", e.g. continuous_data_GMM) or, when data_type is empty, by their reference.
#  Mass-extinction keys add extinction_event and extinction_age(Ma). Series =
#                 outputs/data_tables/mass_extinction/<series id>__<event slug>.csv (cropped to the time after the
#                 extinction); the tables get Extinction_event and Extinction_age_Ma after Clade; clades listed in
#                 input_keys/massext_copy_from_origin.csv are copied from the clade-origin table; series with < 4
#                 points after the extinction are not analysed.
#  extra_series : optional; same columns as the key, `parameter` = exact file suffix of the extra series.
#  Model-order overrides (optional): data_digitized/digitisation/model_order_overrides.csv (Clade, Parameter, N_cpts,
#                 Reason; mass extinctions: model_order_overrides_massext.csv, with an Extinction_event column) fix the
#                 number of shifts of the point estimate and of every MC replicate (reported in Method and Flags).
#  Series       : outputs/data_tables/clade_origin/<series id>.csv with columns Clade, Parameter, Age_Ma, Value,
#                 Value_err, Err_lo, Err_hi, Err_type, Source, Age_cal_sd.
#  Rule         : onset = first upward shift; offset = first decline after the first peak (later peaks ignored).
# ============================================================================
if (!isTRUE(l10n_info()[["UTF-8"]]))           # make sure accented names (e.g. "simões") survive in a C/latin locale
  for (lc in c("C.UTF-8", "en_US.UTF-8", "English_United States.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)))) break
suppressPackageStartupMessages(source(file.path("source_code", "changepoint_pipeline.R")))
source(file.path("source_code", "mc_uncertainty.R")); source(file.path("source_code", "crops_firstpeak.R"))
args <- commandArgs(trailingOnly = TRUE)
DIGI <- Sys.getenv("DIGITISATION_REPO", "ImageDigitization_pipeline")   # submodule: keys, digitised series, copy list
if (!dir.exists(file.path(DIGI, "outputs", "data_tables")))
  stop("ImageDigitization_pipeline not found in '", DIGI, "': run `git submodule update --init` or set DIGITISATION_REPO")
key_file   <- if (length(args) >= 1) args[1] else file.path(DIGI, "input_keys", "key_origin_invert&vert.csv")
tag        <- if (length(args) >= 2) args[2] else paste0("Table_", tools::file_path_sans_ext(basename(key_file)))
extra_file <- if (length(args) >= 3) args[3] else NA
N_REP <- 1000; PEN_RANGE <- c(0.5, 200); MINSEG <- 3; STEP <- 1
DATA_DIR <- "data_digitized"                   # analysis settings (data_digitized/digitisation/model_order_overrides.csv)

rd <- function(f) {                       # UTF-8, or Latin-1 (Excel exports, e.g. "simões"), locale-independent
  x <- readLines(f, warn = FALSE, encoding = "bytes")
  x <- if (all(validUTF8(x))) { Encoding(x) <- "UTF-8"; x } else iconv(x, "latin1", "UTF-8")   # UTF-8 files with accents: re-mark
  x <- sub("^\ufeff", "", x, useBytes = TRUE)   # drop a UTF-8 byte-order mark if present
  utils::read.csv(text = x, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}
key <- rd(key_file); key <- key[!is.na(key$clade) & nzchar(trimws(key$clade)), ]
for (cc in c("group", "legend/line", "data_file", "data_type")) if (!cc %in% names(key)) key[[cc]] <- ""
key[is.na(key)] <- ""
MASS <- "extinction_event" %in% names(key)                     # mass-extinction key
SUB <- Sys.getenv("OUT_SUB", if (MASS) "mass_extinction" else "clade_origin")
RES_DIR <- file.path("results_PELT+CROP", SUB, "results"); FIG_DIR <- file.path("results_PELT+CROP", SUB, "figures")
dir.create(file.path(RES_DIR, "clade_tables"), recursive = TRUE, showWarnings = FALSE); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
CT <- file.path(RES_DIR, "clade_tables")   # clade-specific tables
SERIES_DIR <- file.path(DIGI, "outputs", "data_tables", if (MASS) "mass_extinction" else "clade_origin")
if (MASS) {
  key$extinction_event <- trimws(key$extinction_event)
  key$ext_age <- as.numeric(key[[grep("^extinction_age", names(key), value = TRUE)[1]]])
  cpy <- rd(file.path(DIGI, "input_keys", "massext_copy_from_origin.csv"))
  origin_tab <- rd(file.path("results_PELT+CROP", "clade_origin", "results", "Table_CladeOrigin.csv"))
} else key$extinction_event <- ""
key$extra <- FALSE
if (!is.na(extra_file) && file.exists(extra_file)) {
  ex <- rd(extra_file); ex$extra <- TRUE
  for (cc in setdiff(names(key), names(ex))) ex[[cc]] <- ""
  key <- rbind(key, ex[, names(key)])
}
slug <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", trimws(x)))
ovr_file <- file.path(DATA_DIR, "digitisation", if (MASS) "model_order_overrides_massext.csv" else "model_order_overrides.csv")
ovr <- if (file.exists(ovr_file)) rd(ovr_file) else NULL

std_param <- function(x, extra = FALSE) {
  if (extra) return(trimws(x)); x <- tolower(trimws(x))
  if (grepl("rate", x)) "phenotypic_rate" else if (grepl("dispar", x)) "phenotypic_disparity"
  else if (grepl("divers", x)) "taxonomic_diversity" else gsub(" +", "_", x)
}

key$pr_std <- mapply(function(p, e) std_param(p, e), key$parameter, key$extra)
key$dt_slug <- slug(key$data_type)
sid0 <- paste(key$clade, key$pr_std, key$dt_slug, key$extinction_event)   # only rows whose series id would collide (within an event)
dup <- duplicated(sid0) | duplicated(sid0, fromLast = TRUE)
nodt <- dup & !nzchar(key$dt_slug) & !(tolower(trimws(key$reference)) %in% c("na", ""))
key$dt_slug[nodt] <- slug(key$reference[nodt])                  # same clade x parameter, no data_type -> use reference
key$data_type[nodt] <- key$reference[nodt]
out <- list()
for (i in seq_len(nrow(key))) {
  cl <- key$clade[i]; pr <- std_param(key$parameter[i], key$extra[i])
  dt <- key$dt_slug[i]; sid <- paste0(cl, "_", pr, if (nzchar(dt)) paste0("_", dt))   # series id
  pr_dt <- sub(paste0("^", cl, "_"), "", sid)
  ev <- key$extinction_event[i]; sidf <- if (MASS) paste0(sid, "__", slug(ev)) else sid   # file / figure id
  note <- character(0)
  fig <- key[["figure/table"]][i]; dfile <- trimws(key$data_file[i])
  if (nzchar(dfile) && !grepl("raw data|converted data", fig, ignore.case = TRUE))
    note <- c(note, paste0("Key data_file column: '", dfile, "'"))
  base <- data.frame(Group = key$group[i], Clade = cl, Parameter = pr, Data_type = key$data_type[i],
                     Reference = key$reference[i], Figure = fig, Legend_line = key[["legend/line"]][i],
                     Key_data_file = dfile, stringsAsFactors = FALSE)
  if (MASS) base <- cbind(base[, 1:2], Extinction_event = ev, Extinction_age_Ma = key$ext_age[i], base[, -(1:2)])
  f <- file.path(SERIES_DIR, paste0(sidf, ".csv"))
  if (MASS && cl %in% cpy$Clade && !(tolower(trimws(key$reference[i])) %in% c("na", ""))) {   # copied from the clade-origin table
    o <- origin_tab[origin_tab$Clade == cl & origin_tab$Parameter == pr & origin_tab$Reference == key$reference[i], , drop = FALSE]
    if (nrow(o) != 1) stop("clade-origin row not found once: ", cl, " ", pr, " ", key$reference[i])
    o <- cbind(o[, 1:2], Extinction_event = ev, Extinction_age_Ma = key$ext_age[i], o[, -(1:2)])
    o$Flags <- paste(c(sprintf("Copied from the clade-origin analysis (results_PELT+CROP/clade_origin): %s", cpy$Reason[cpy$Clade == cl][1]),
                       if (!is.na(o$Flags) && nzchar(o$Flags)) o$Flags), collapse = " | ")
    out[[i]] <- o; message(sprintf("%-24s %-42s copied from clade origin", cl, pr_dt)); next
  }
  if (tolower(trimws(key$reference[i])) %in% c("na", "") || !file.exists(f)) {
    msg <- if (tolower(trimws(key$reference[i])) %in% c("na", "")) "No published source for this parameter (key: na)"
           else paste("Data file not found:", f)
    out[[i]] <- cbind(base, Method = NA, Flags = paste(c(msg, note), collapse = " | "))
    next
  }
  df <- read_series(f)
  if (MASS && nrow(df) < 4) {
    out[[i]] <- cbind(base, Data_source = df$Source[1], N_points = nrow(df), Method = NA,
                      Flags = paste(c(sprintf("Too few data points after the extinction (n = %d, %s Ma): not analysed", nrow(df),
                                              paste(round(df$Age_Ma, 1), collapse = ", ")), note), collapse = " | "))
    message(sprintf("%-24s %-42s %s: only %d points after the extinction - skipped", cl, pr_dt, ev, nrow(df))); next
  }
  g <- resample_grid(df$Age_Ma, df$Value, step = STEP)
  o <- if (!is.null(ovr)) ovr[ovr$Clade == cl & ovr$Parameter == pr_dt &
                              (if (MASS && "Extinction_event" %in% names(ovr)) trimws(ovr$Extinction_event) == ev else TRUE), , drop = FALSE] else NULL
  n_fix <- if (!is.null(o) && nrow(o)) as.integer(o$N_cpts[1]) else NA
  cr <- crops_first_peak(g, PEN_RANGE, MINSEG, n_fix)
  utils::write.csv(cr$table, file.path(CT, paste0(sidf, "_CROPS_solutions.csv")), row.names = FALSE)
  age_sd <- if ("Age_cal_sd" %in% names(df)) df$Age_cal_sd[1] else 0
  mc_file <- file.path(CT, paste0(sidf, "_CROPS_MCreps.csv"))
  reuse <- Sys.getenv("REUSE_MC") == "1" && file.exists(mc_file) && file.mtime(mc_file) > file.mtime(f)
  if (reuse) { reps <- utils::read.csv(mc_file); reuse <- "mc_mode" %in% names(reps) && reps$mc_mode[1] == "conditioned_exact" && nrow(reps) == N_REP }
  if (!reuse)                                 # (REUSE_MC=1 resumes an interrupted run from saved replicates)
    reps <- mc_crops_first_peak(df, N_REP, PEN_RANGE, MINSEG, STEP, age_sd,
                                nonneg = pr == "taxonomic_diversity", n_fix = cr$n_cpts)   # MC conditioned on the PE model order
  utils::write.csv(reps, file.path(CT, paste0(sidf, "_CROPS_MCreps.csv")), row.names = FALSE)
  s <- summarise_mc(reps)
  ncp <- table(reps$n_elbow); ncp_txt <- paste(sprintf("%s:%.0f%%", names(ncp), 100 * ncp / sum(ncp)), collapse = ";")
  # ORIGIN focus: keep only shifts bounding the first peak (ignore later peaks)
  sh_origin <- cr$shift_ages[cr$shift_ages >= cr$offset - 0.01]
  # temporal resolution: onset/offset can only be placed between two sampled bins
  brk <- function(t) { old <- df$Age_Ma[df$Age_Ma >= t - 1e-6]; yng <- df$Age_Ma[df$Age_Ma <= t + 1e-6]
    c(if (length(old)) round(min(old), 1) else NA, if (length(yng)) round(max(yng), 1) else NA) }
  oldest <- max(df$Age_Ma); youngest <- min(df$Age_Ma)
  ob <- if (abs(cr$onset - oldest) < 1) rep(round(oldest, 1), 2) else brk(cr$onset)
  fb <- if (abs(cr$offset - youngest) < 1) rep(round(youngest, 1), 2) else brk(cr$offset)
  et <- if ("Err_type" %in% names(df)) df$Err_type[1] else "unknown"
  flag <- c(
    if (cr$n_cpts == 0) "No shift detected (CROPS elbow = 0 changepoints): no first peak; onset/offset = series limits",
    if (cr$n_cpts > 0 && abs(cr$onset - oldest) < 1) "Onset = oldest data point (left-censored: series starts at peak level)",
    if (cr$n_cpts > 0 && abs(cr$offset - youngest) < 1) "Offset = youngest data point (right-censored: no decline after peak)",
    if (cr$n_cpts >= 1 && abs(cr$offset - youngest) < 1 && abs(cr$onset - oldest) >= 1)
      "No first peak: monotonic rise to a plateau; onset = first upward shift, offset not interpretable",
    if (cr$peak_seg > 2) sprintf("Stepwise rise: first peak reached after %d upward shifts; onset = first upward shift (later steps in Origin_shift_ages_Ma)", cr$peak_seg - 1),
    if (!is.na(n_fix)) sprintf("Model order fixed at %d shifts (kneedle elbow = %d): %s", cr$n_cpts, cr$n_elbow, o$Reason[1]),
    if ((s$Onset_Max - s$Onset_Min) > 0.5 * (oldest - youngest)) "Onset 95% MC range > 50% of series span: weakly constrained",
    if (max(table(reps$n_elbow)) / sum(!is.na(reps$n_elbow)) < 0.6) "CROPS elbow unstable across MC replicates (<60% of replicates' own elbows agree on n shifts; MC conditioned on the point-estimate model order)",
    if (et == "none") "No error bars in source: MC reflects digitisation error only; use bin brackets for timing uncertainty",
    if (et == "Poisson") "Value error assumed Poisson (sqrt(n)) for counts without published error bars",
    if (MASS && nrow(df) < 8) sprintf("Few data points after the extinction (n = %d): shifts poorly constrained", nrow(df)),
    note)
  out[[i]] <- cbind(base, data.frame(
    Data_source = df$Source[1], Error_type = et, N_points = nrow(df),
    Series_start_Ma = round(oldest, 1), Series_end_Ma = round(youngest, 1),
    Method = if (is.na(n_fix)) "PELT+CROPS (elbow)" else sprintf("PELT+CROPS (fixed %d shifts; elbow %d)", cr$n_cpts, cr$n_elbow),
    Elbow_n_cpts = cr$n_cpts,
    Origin_shift_ages_Ma = paste(sh_origin, collapse = ";"),
    Onset_PE = cr$onset, Onset_bin_older = ob[1], Onset_bin_younger = ob[2],
    Offset_PE = cr$offset, Offset_bin_older = fb[1], Offset_bin_younger = fb[2],
    Peak_Age_PE = cr$peak_age, s, MC_n_cpts_freq = ncp_txt,
    Flags = paste(flag, collapse = " | "), stringsAsFactors = FALSE))

  png(file.path(FIG_DIR, paste0(sidf, "_CROPS.png")), width = 1500, height = 620, res = 130)
  layout(matrix(1:2, 1), widths = c(1, 2.1)); par(mar = c(4.2, 4.2, 3, 1))
  tb <- cr$table
  plot(tb$n_cpts, tb$cost, type = "b", pch = 16, xlab = "Number of changepoints",
       ylab = "Unpenalised cost (standardised RSS)", main = "CROPS elbow")
  points(tb$n_cpts[tb$selected], tb$cost[tb$selected], pch = 21, cex = 2.4, col = "red", lwd = 2)
  if (!is.na(n_fix)) { points(tb$n_cpts[tb$elbow], tb$cost[tb$elbow], pch = 24, cex = 1.8, col = "grey40", lwd = 2)
    legend("topright", bty = "n", cex = .7, pch = c(21, 24), col = c("red", "grey40"), legend = c("used (fixed)", "kneedle elbow")) }
  lo <- df$Err_lo; hi <- df$Err_hi
  plot(df$Age_Ma, df$Value, pch = 16, cex = if (nrow(df) > 100) .4 else 1,
       xlim = rev(range(c(df$Age_Ma, s$Onset_Max, s$Offset_Min), na.rm = TRUE)),
       ylim = range(c(lo, hi), na.rm = TRUE), xlab = "Age (Ma)", ylab = pr_dt,
       main = sprintf("%s - %s%s: CROPS elbow (%d shifts)", cl, pr_dt, if (MASS) paste0(" after ", ev, " (", key$ext_age[i], " Ma)") else "", cr$n_cpts))
  usr <- par("usr")
  rect(s$Onset_Max, usr[3], s$Onset_Min, usr[4], col = adjustcolor("forestgreen", .18), border = NA)
  rect(s$Offset_Max, usr[3], s$Offset_Min, usr[4], col = adjustcolor("purple", .18), border = NA)
  if (nrow(df) <= 100) arrows(df$Age_Ma, lo, df$Age_Ma, hi, angle = 90, code = 3, length = .03, col = "grey40")
  else polygon(c(df$Age_Ma, rev(df$Age_Ma)), c(lo, rev(hi)), col = adjustcolor("grey50", .25), border = NA)
  lines(g$age, g$value, col = "#c0392b")
  lines(g$age, rep(cr$seg_means, cr$ends - cr$starts + 1), col = "steelblue", lwd = 2, type = "s")
  abline(v = sh_origin, col = "grey20", lty = 3)
  abline(v = c(cr$onset, cr$offset), col = c("forestgreen", "purple"), lty = 2, lwd = 2)
  legend("topright", bty = "n", cex = .75,
    legend = c(sprintf("onset %.1f Ma [95%% MC %.1f-%.1f]", cr$onset, s$Onset_Max, s$Onset_Min),
               sprintf("offset %.1f Ma [95%% MC %.1f-%.1f]", cr$offset, s$Offset_Max, s$Offset_Min),
               "CROPS-elbow segment means", "origin shift points"),
    col = c("forestgreen", "purple", "steelblue", "grey20"), lty = c(2, 2, 1, 3), lwd = 2)
  dev.off()
  message(sprintf("%-24s %-42s %s onset %6.1f [%6.1f-%6.1f]  offset %6.1f [%6.1f-%6.1f]  n_cpts %d",
                  cl, pr_dt, ev, cr$onset, s$Onset_Max, s$Onset_Min, cr$offset, s$Offset_Max, s$Offset_Min, cr$n_cpts))
}
# bind rows with differing columns (na rows carry fewer columns)
cols <- unique(unlist(lapply(out, names)))
res <- do.call(rbind, lapply(out, function(d) { for (cc in setdiff(cols, names(d))) d[[cc]] <- NA; d[, cols] }))
utils::write.csv(res, file.path(RES_DIR, paste0(tag, ".csv")), row.names = FALSE, na = "", fileEncoding = "UTF-8")
