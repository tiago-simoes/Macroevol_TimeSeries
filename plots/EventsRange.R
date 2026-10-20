# EventsRange.R -- range plots of the onset-to-offset intervals (Fig. 1: clade origins; Fig. 2: mass extinctions).
# Reads plots/Data_Origin_Vert.csv, Data_Origin_Invert.csv, Data_MassExct_Vert.csv and Data_MassExct_Invert.csv
# (onset/offset ranges of the Summary tables in results_combined/) and writes plots/Fig1_Origins.pdf and
# plots/Fig2_MassExct.pdf. Run from the repository root (Rscript plots/EventsRange.R) or from plots/.
# Requires ggplot2, dplyr, deeptime, ggpattern and cowplot.
if (!file.exists("Data_Origin_Vert.csv") && file.exists(file.path("plots", "Data_Origin_Vert.csv"))) setwd("plots")
if (!interactive()) pdf(NULL)   # Rscript: no Rplots.pdf from the plots printed on screen

library(ggplot2)
library(dplyr)
library(deeptime)
library(ggpattern)
library(cowplot)   # composite figures (end of script)

########## Helpers for the clade-origin and mass-extinction plots ############
# Bands inside each clade row, top to bottom: clade origin, Div, Div2, Disp, Disp2, Rate, Rate2.
# Second-series columns (Div2_, Disp2_, Rate2_ + Onset_Max/Onset_Min/Offset_Max/Offset_Min) are optional:
# a clade with values there gets an extra band, and its bands are narrowed so the row keeps the same height
# (span 0.8 around the clade's y). Clades without a second series keep the original 0.2-high bands.
# Every box fills its band exactly (gap = 0). gap > 0 lowers the top of a second-series band by that fraction
# of the band height, leaving a white line between two same-coloured bands (the box is then shorter).
# origin = FALSE (mass-extinction plots): no clade-origin band. If a clade label has several table rows (all drawn
# in the same clade row), a second-series band is used in every row of that label if any of its rows has values
# there, so all its boxes share the same bands.
add_bands <- function(d, span = 0.8, gap = 0, origin = TRUE) {
  lev   <- c(if (origin) "Origin", "Div", "Div2", "Disp", "Disp2", "Rate", "Rate2")
  extra <- c("Div2", "Disp2", "Rate2")
  rng   <- c("_Onset_Max", "_Onset_Min", "_Offset_Max", "_Offset_Min")
  for (p in extra) for (s in rng) if (is.null(d[[paste0(p, s)]])) d[[paste0(p, s)]] <- NA_real_
  num <- grep("_(Onset|Offset)_(Max|Min)$|_Age$", names(d))      # empty columns are read as logical
  d[num] <- lapply(d[num], as.numeric)
  used <- matrix(sapply(lev, function(l)
    if (l %in% extra) ave(rowSums(!is.na(d[paste0(l, rng)])), d$Clade, FUN = max) > 0 else rep(TRUE, nrow(d))),
    nrow = nrow(d), dimnames = list(NULL, lev))
  h   <- span / rowSums(used)                                   # band height in each clade row
  pos <- matrix(t(apply(used, 1, cumsum)), nrow = nrow(d))      # band number, counted from the top
  top <- as.numeric(d$Clade) + span / 2
  for (j in seq_along(lev)) {
    d[[paste0(lev[j], "_ymax")]] <- ifelse(used[, j], top - (pos[, j] - 1 + gap * (lev[j] %in% extra)) * h, NA_real_)
    d[[paste0(lev[j], "_ymin")]] <- ifelse(used[, j], top - pos[, j] * h, NA_real_)
  }
  d
}

# Parameter p: solid box from Onset_Min to Offset_Max, hatched boxes for onset and offset uncertainty
range_layers <- function(p, fill) {
  v <- function(s) paste0(p, s)
  hatch <- function(from, to)
    geom_rect_pattern(aes(xmin = .data[[v(from)]], xmax = .data[[v(to)]],
                          ymin = .data[[v("_ymin")]], ymax = .data[[v("_ymax")]]),
                      pattern = "stripe", pattern_density = 0.8, pattern_angle = 90,
                      pattern_color = "white", pattern_fill = fill, alpha = 0.01,
                      pattern_spacing = 0.01, pattern_linetype = 1, pattern_size = 0.01, na.rm = TRUE)
  list(geom_rect(aes(xmin = .data[[v("_Onset_Min")]], xmax = .data[[v("_Offset_Max")]],
                     ymin = .data[[v("_ymin")]], ymax = .data[[v("_ymax")]]),
                 fill = fill, alpha = 0.7, na.rm = TRUE),
       hatch("_Onset_Max", "_Onset_Min"),
       hatch("_Offset_Max", "_Offset_Min"))
}

# Superscripts for clades with two series of one parameter: a small number right of the bar range
# (to the right of its youngest bound) of each of the two series; their meaning (method / data type) goes into the
# figure caption. Keys "<Clade label>|<band>": band <P> = first series, <P>2 = second series (columns <P>_ and <P>2_
# of the plot tables). The same number is used in every figure. A clade with a
# second series that is missing here stops the script, so that no pair is left unlabelled.
SUP <- c("Dinosauria (early)|Div"  = "1", "Dinosauria (early)|Div2"  = "2",   # Occur, PyRate | Phylo, FBDD
         "Aves (crown)|Rate"       = "3", "Aves (crown)|Rate2"       = "4",   # GMM (Skull), PCM(BT) | LMM (WB), PCM(BT)
         "Lepidosauria|Rate"       = "5", "Lepidosauria|Rate2"       = "6",   # GMM, PCM(BT) | Discrete (WB), Clock
         "Mammalia (early)|Div"    = "7", "Mammalia (early)|Div2"    = "8")   # Phylo, LTT | Phylo, BDD
sup_labels <- function(d, xlim, size = 2.0, nudge = 0.006) {
  dx  <- nudge * abs(diff(xlim))                 # gap between the bar and the number (0.6% of the age axis)
  rng <- c("_Onset_Max", "_Onset_Min", "_Offset_Max", "_Offset_Min")
  lab <- NULL
  for (p in c("Div", "Disp", "Rate")) {
    two <- which(!is.na(d[[paste0(p, "2_ymin")]]))          # rows of clades with a second series of p
    for (i in two) for (q in c(p, paste0(p, "2"))) {
      v <- unlist(d[i, paste0(q, rng)])
      if (all(is.na(v))) next
      key <- paste(as.character(d$Clade[i]), q, sep = "|")
      if (!key %in% names(SUP)) stop("no superscript for ", key, ": add it to SUP")
      lab <- rbind(lab, data.frame(x = min(v, na.rm = TRUE) - dx,                  # right of the youngest bound
                                   y = d[[paste0(q, "_ymin")]][i] + 0.6 * (d[[paste0(q, "_ymax")]][i] - d[[paste0(q, "_ymin")]][i]),
                                   label = SUP[[key]]))
    }
  }
  if (is.null(lab)) return(NULL)
  geom_text(data = lab, aes(x = x, y = y, label = label), hjust = 0, vjust = 0.5, size = size, inherit.aes = FALSE)
}

########## Origins Dynamic: vertebrates ############

data_orig <- read.csv("Data_Origin_Vert.csv")
# Set the desired order for Clade levels

# bottom to top; one row per clade (second series of a parameter = extra band in the row)
clade_order <- c("Placentalia", "Mammalia (early)", "Lepidosauria", "Aves (crown)",
                 "Aves (stem)", "Theropoda (non-avian)", "Dinosauria (early)",
                 "Archosauromorpha (stem)", "Tetrapoda (stem)", "Acanthomorpha")

stopifnot("Clade labels in the table missing from clade_order" = all(data_orig$Clade %in% clade_order))
data_orig$Clade <- factor(data_orig$Clade, levels = clade_order, ordered = TRUE)
data_orig <- add_bands(data_orig)   # y-positions of the bands in each clade row

# Generate y-positions for the dotted lines between each clade
y_positions <- seq(1.5, length(clade_order) - 0.5, by = 1)

### Plots

P1a <- ggplot(data_orig) +
## Geoscale
  coord_geo(xlim = c(400, 0), ylim = c(0.5, length(clade_order) + 0.5), expand = FALSE, clip = "on",
            dat = list("eras", "periods"), abbrv = list(FALSE, FALSE), 
            skip = c("Quaternary"),
            pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(4, 3), neg = FALSE) +
  
  geom_rect(aes(xmin = 400, xmax = 358.9, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 298.9, xmax = 252, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 200, xmax = 145, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 66, xmax = 23.03, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 2.58, xmax = 0, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +

  geom_hline(yintercept = y_positions, linetype = "dotted", color = "black") +
 
## Clade origin (box) and median age (diamond)
  geom_rect(aes(xmin = Max_Age, xmax = Min_Age, ymin = Origin_ymin, ymax = Origin_ymax),
            fill = "purple3", alpha = 0.7) +
  geom_point(aes(x = Median_Age, y = (Origin_ymin + Origin_ymax) / 2), color = "black", shape = 18, size = 2) +

## Onset to offset of each parameter (solid box) with hatched onset and offset uncertainty.
## Div2 / Disp2 / Rate2 = second series of a parameter: own band below the first, only where present
  range_layers("Div",   "cyan2") +
  range_layers("Div2",  "cyan2") +
  range_layers("Disp",  "orange2") +
  range_layers("Disp2", "orange2") +
  range_layers("Rate",  "tomato") +
  range_layers("Rate2", "tomato") +
  sup_labels(data_orig, c(400, 0)) +   # superscripts of the two-series pairs

  # Customize axis labels and theme
  labs(title = "Dynamics after origin of major clades: vertebrates",
       x = "Ma", 
       y = NULL, 
       color = "Event Type") +
  theme_minimal() +
  
  # Set custom y-axis and x-axis scales with the reversed order
  scale_y_continuous(breaks = 1:length(clade_order), labels = clade_order) +
  scale_x_continuous(limits = c(400, 0), breaks = seq(400, 0, -25), minor_breaks = seq(400, 0, -10), trans = "reverse") +
  
  # Keep only axis lines and labels
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),   # short ticks joining the age labels to the axis
    axis.ticks.length.x = unit(4, "pt"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

P1a
#ggsave("Vert_Origin.pdf", width = 10, height = 6)



########## Origins Dynamics: invertebrates ############

data_orig <- read.csv("Data_Origin_Invert.csv")

# Set the desired order for Clade levels
clade_order <- c("Echinoidea (Irregularia)", "Crinoidea (Eucladida)", "Echinodermata",
                 "Trilobita", "Pancrustacea (non-Hexapoda)", "Bivalvia")

stopifnot("Clade labels in the table missing from clade_order" = all(data_orig$Clade %in% clade_order))
data_orig$Clade <- factor(data_orig$Clade, levels = clade_order, ordered = TRUE)
data_orig <- add_bands(data_orig)   # y-positions of the bands in each clade row

# Generate y-positions for the dotted lines between each clade
y_positions <- seq(1.5, length(clade_order) - 0.5, by = 1)

### Plots

P1b <- ggplot(data_orig) +
  ## Geoscale
  coord_geo(xlim = c(600, 145), ylim = c(0.5, length(clade_order) + 0.5), expand = FALSE, clip = "on",
            dat = list("eras", "periods"), abbrv = list(FALSE, FALSE), 
            skip = c("Quaternary"),
            pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(4, 3), neg = FALSE) +
  
  geom_rect(aes(xmin = 538.8, xmax = 486.8, ymin = -Inf, ymax = Inf), #Cambrian
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 443.1, xmax = 419.6, ymin = -Inf, ymax = Inf), #Silurian
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 358.9, xmax = 298.9, ymin = -Inf, ymax = Inf), #Carb
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 251.9, xmax = 201.4, ymin = -Inf, ymax = Inf), #Triassic
            fill = "grey70", alpha = 0.01) +

  geom_hline(yintercept = y_positions, linetype = "dotted", color = "black") +
  
## Clade origin (box) and median age (diamond)
  geom_rect(aes(xmin = Max_Age, xmax = Min_Age, ymin = Origin_ymin, ymax = Origin_ymax),
            fill = "purple3", alpha = 0.7) +
  geom_point(aes(x = Median_Age, y = (Origin_ymin + Origin_ymax) / 2), color = "black", shape = 18, size = 2) +

## Onset to offset of each parameter (solid box) with hatched onset and offset uncertainty.
## Div2 / Disp2 / Rate2 = second series of a parameter: own band below the first, only where present
  range_layers("Div",   "cyan2") +
  range_layers("Div2",  "cyan2") +
  range_layers("Disp",  "orange2") +
  range_layers("Disp2", "orange2") +
  range_layers("Rate",  "tomato") +
  range_layers("Rate2", "tomato") +
  sup_labels(data_orig, c(600, 145)) +   # superscripts of the two-series pairs

  # Customize axis labels and theme
  labs(title = "Dynamics after origin of major clades: invertebrates",
       x = "Ma", 
       y = NULL, 
       color = "Event Type") +
  theme_minimal() +
  
  # Set custom y-axis and x-axis scales with the reversed order
  scale_y_continuous(breaks = 1:length(clade_order), labels = clade_order) +
  scale_x_continuous(limits = c(600, 145), breaks = seq(600, 145, -25), minor_breaks = seq(400, 145, -10), trans = "reverse") +
  
  # Keep only axis lines and labels
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),   # short ticks joining the age labels to the axis
    axis.ticks.length.x = unit(4, "pt"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

P1b
#ggsave("Invert_Origin.pdf", width = 10, height = 6)



##############################################################################



########## Mass Exct Dynamics: vertebrates ############
# Data_MassExct_Vert.csv / _Invert.csv: ranges of results_combined/Summary_MassExtinction.csv; one row per clade x extinction (Exct_Age = extinction age of the key).
# Same helpers as the clade-origin plots, without the clade-origin band: bands inside each clade row, top to bottom,
# Div, Div2, Disp, Disp2, Rate, Rate2 (second series only where present). span = 0.6 gives the three usual
# 0.2-high bands, centred on the clade row. A clade analysed after two extinctions has two table rows with the same
# label: both are drawn in the same clade row, so their bars may overlap.

data_mass <- read.csv("Data_MassExct_Vert.csv")

# Set the desired order for Clade levels (bottom to top)
clade_order <- c("Placentalia", "Mammalia (early)", "Aves (crown)", "Dinosauria (early)",
                 "Archosauromorpha (stem)", "Tetrapoda (stem)", "Acanthomorpha")

stopifnot("Clade labels in the table missing from clade_order" = all(data_mass$Clade %in% clade_order))
data_mass$Clade <- factor(data_mass$Clade, levels = clade_order, ordered = TRUE)
data_mass <- add_bands(data_mass, span = 0.6, origin = FALSE)   # y-positions of the bands in each clade row

# Mass extinctions (dashed red lines); every extinction age in the table must have a line
ext_lines <- c("end-Devonian (Famennian)" = 359, "end-Guadalupian" = 259.51, "end-Permian" = 251.9,
               "end-Triassic" = 201.4, "K-Pg" = 66)
stopifnot("Extinction age in the table without a dashed line" = all(data_mass$Exct_Age %in% ext_lines))

# Generate y-positions for the dotted lines between each clade
y_positions <- seq(1.5, length(clade_order) - 0.5, by = 1)

### Plots

P2a <- ggplot(data_mass) +
  ## Geoscale
  coord_geo(xlim = c(400, 0), ylim = c(0.5, length(clade_order) + 0.5), expand = FALSE, clip = "on",
            dat = list("eras", "periods"), abbrv = list(FALSE, FALSE), 
            skip = c("Quaternary"),
            pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(4, 3), neg = FALSE) +
  
  geom_rect(aes(xmin = 400, xmax = 358.9, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 298.9, xmax = 252, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 200, xmax = 145, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 66, xmax = 23.03, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 2.58, xmax = 0, ymin = -Inf, ymax = Inf), 
            fill = "grey70", alpha = 0.01) +
  geom_vline(xintercept = ext_lines, linetype = "dashed", color = "red", linewidth = 0.7) +
  
  geom_hline(yintercept = y_positions, linetype = "dotted", color = "black") +
  
## Onset to offset of each parameter (solid box) with hatched onset and offset uncertainty.
## Div2 / Disp2 / Rate2 = second series of a parameter: own band below the first, only where present
  range_layers("Div",   "cyan2") +
  range_layers("Div2",  "cyan2") +
  range_layers("Disp",  "orange2") +
  range_layers("Disp2", "orange2") +
  range_layers("Rate",  "tomato") +
  range_layers("Rate2", "tomato") +
  sup_labels(data_mass, c(400, 0)) +   # superscripts of the two-series pairs
  
  # Customize axis labels and theme
  labs(title = "Dynamics after mass extinctions: vertebrates",
       x = "Ma", 
       y = NULL, 
       color = "Event Type") +
  theme_minimal() +
  
  # Set custom y-axis and x-axis scales with the reversed order
  # (oob_keep: a box running past the axis limits is cut at the plot edge instead of being dropped)
  scale_y_continuous(breaks = 1:length(clade_order), labels = clade_order) +
  scale_x_continuous(limits = c(400, 0), breaks = seq(400, 0, -25), minor_breaks = seq(400, 0, -10), trans = "reverse",
                     oob = scales::oob_keep) +
  
  # Keep only axis lines and labels
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),   # short ticks joining the age labels to the axis
    axis.ticks.length.x = unit(4, "pt"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

P2a
#ggsave("Vert_MassExct.pdf", width = 10, height = 5)




########## Mass Exct Dynamics: invertebrates ############

data_mass <- read.csv("Data_MassExct_Invert.csv")

# Set the desired order for Clade levels (bottom to top; one row per clade, also for clades analysed after
# two extinctions)
clade_order <- c("Echinoidea (Irregularia)", "Crinoidea (Diplobathrida)", "Crinoidea (Eucladida)",
                 "Brachiopoda (Athyridida)", "Brachiopoda (Spiriferinida)", "Brachiopoda (Rhynchonellida)",
                 "Brachiopoda (Terebratulida)", "Trilobita", "Eurypterida", "Ammonoidea")

stopifnot("Clade labels in the table missing from clade_order" = all(data_mass$Clade %in% clade_order))
data_mass$Clade <- factor(data_mass$Clade, levels = clade_order, ordered = TRUE)
data_mass <- add_bands(data_mass, span = 0.6, origin = FALSE)   # y-positions of the bands in each clade row

# Mass extinctions (dashed red lines); every extinction age in the table must have a line
ext_lines <- c("end-Ordovician" = 443, "end-Devonian (Frasnian)" = 374, "end-Devonian (Famennian)" = 359,
               "end-Guadalupian" = 259.51, "end-Permian" = 251.9, "end-Triassic" = 201.4)   # end-Triassic: Echinoidea (Irregularia)
stopifnot("Extinction age in the table without a dashed line" = all(data_mass$Exct_Age %in% ext_lines))

# Generate y-positions for the dotted lines between each clade
y_positions <- seq(1.5, length(clade_order) - 0.5, by = 1)

### Plots

P2b <- ggplot(data_mass) +
  ## Geoscale
  coord_geo(xlim = c(475, 66), ylim = c(0.5, length(clade_order) + 0.5), expand = FALSE, clip = "on",
            dat = list("eras", "periods"), abbrv = list(FALSE, FALSE), 
            skip = c("Quaternary"),
            pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(4, 3), neg = FALSE) +
  
  geom_rect(aes(xmin = 538.8, xmax = 486.8, ymin = -Inf, ymax = Inf), #Cambrian
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 443.1, xmax = 419.6, ymin = -Inf, ymax = Inf), #Silurian
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 358.9, xmax = 298.9, ymin = -Inf, ymax = Inf), #Carb
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 251.9, xmax = 201.4, ymin = -Inf, ymax = Inf), #Triassic
            fill = "grey70", alpha = 0.01) +
  geom_rect(aes(xmin = 145, xmax = 66, ymin = -Inf, ymax = Inf), #Cretaceous
            fill = "grey70", alpha = 0.01) +
  geom_vline(xintercept = ext_lines, linetype = "dashed", color = "red", linewidth = 0.7) +
  
  geom_hline(yintercept = y_positions, linetype = "dotted", color = "black") +
  
## Onset to offset of each parameter (solid box) with hatched onset and offset uncertainty.
## Div2 / Disp2 / Rate2 = second series of a parameter: own band below the first, only where present
  range_layers("Div",   "cyan2") +
  range_layers("Div2",  "cyan2") +
  range_layers("Disp",  "orange2") +
  range_layers("Disp2", "orange2") +
  range_layers("Rate",  "tomato") +
  range_layers("Rate2", "tomato") +
  sup_labels(data_mass, c(475, 66)) +   # superscripts of the two-series pairs
  
  # Customize axis labels and theme
  labs(title = "Dynamics after mass extinctions: invertebrates",
       x = "Ma", 
       y = NULL, 
       color = "Event Type") +
  theme_minimal() +
  
  # Set custom y-axis and x-axis scales with the reversed order
  # (x-axis to 66 Ma so that the Brachiopoda offsets are shown; oob_keep: a box running past the
  #  axis limits, e.g. the Terebratulida no-offset disparity bar to 5.8 Ma, is cut at the plot edge instead of
  #  being dropped)
  scale_y_continuous(breaks = 1:length(clade_order), labels = clade_order) +
  scale_x_continuous(limits = c(475, 66), breaks = seq(475, 66, -25), minor_breaks = seq(475, 66, -10), trans = "reverse",
                     oob = scales::oob_keep) +
  
  # Keep only axis lines and labels
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),   # short ticks joining the age labels to the axis
    axis.ticks.length.x = unit(4, "pt"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

P2b
#ggsave("Invert_MassExct.pdf", width = 10, height = 6)






########################################################################################
########## Composite figures: Fig. 1 clade origins, Fig. 2 mass extinctions ############


row_in   <- 0.45   # height of one clade row (in)
extra_in <- 1.10   # title + geological timescale + age labels of one panel (in)
legend_in <- 0.60  # legend strip (in)

panel <- function(p, title) {
  p + labs(title = title) +
    theme(axis.title.x = element_text(size = 12, hjust = 1.045, margin = margin(t = -14)),   # "Ma" at the axis end
          plot.margin = margin(t = 5, r = 30, b = 5, l = 15))
}
panel_height <- function(p) nlevels(p$data$Clade) * row_in + extra_in

# Legend strip drawn as a small plot, in inches (x from 0 to the figure width): keys listed column by column
# (nrow keys per column), columns col_in wide and centred in the strip; each key = a 0.3-in symbol and its label.
# type = "box" (coloured box), "point" (diamond) or "dashed" (red dashed line)
legend_panel <- function(keys, nrow = 2, col_in = 3.3, width_in = 10, text_size = 10) {
  k <- seq_len(nrow(keys)) - 1
  ncol <- max(k %/% nrow) + 1
  keys$x <- (width_in - ncol * col_in) / 2 + (k %/% nrow) * col_in      # left end of the symbol
  keys$y <- -(k %% nrow)
  ggplot(keys) +
    geom_rect(data = keys[keys$type == "box", ],
              aes(xmin = x, xmax = x + 0.3, ymin = y - 0.12, ymax = y + 0.12, fill = colour), alpha = 0.7) +
    geom_point(data = keys[keys$type == "point", ], aes(x = x + 0.15, y = y), colour = "black", shape = 18, size = 3) +
    geom_segment(data = keys[keys$type == "dashed", ], aes(x = x + 0.15, xend = x + 0.15, y = y - 0.4, yend = y + 0.4),
                 colour = "red", linetype = "dashed", linewidth = 0.7) +
    geom_text(aes(x = x + 0.42, y = y, label = label), hjust = 0, size = text_size / .pt) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, width_in), ylim = c(-nrow + 0.5, 0.5), expand = FALSE) +
    theme_void()
}

### Fig. 1: clade origins
legend_origin <- legend_panel(data.frame(
  label  = c("Clade Origin (Median age)", "Clade Origin (95% HPD)", "Phenotypic Disparity Expansion",
             "Phenotypic Rate Elevated", "Diversification Peak"),
  type   = c("point", "box", "box", "box", "box"),
  colour = c("black", "purple3", "orange2", "tomato", "cyan2")))

h1 <- c(panel_height(P1b), panel_height(P1a), legend_in)
Fig1 <- plot_grid(panel(P1b, "Invertebrates"), panel(P1a, "Vertebrates"), legend_origin,
                  ncol = 1, rel_heights = h1, labels = c("a", "b", ""), label_size = 12,
                  label_fontface = "bold", align = "v", axis = "lr")
Fig1
ggsave("Fig1_Origins.pdf", Fig1, width = 10, height = sum(h1))

### Fig. 2: mass extinctions
legend_massext <- legend_panel(data.frame(
  label  = c("Mass Extinctions", "Diversification Peak", "Phenotypic Disparity Expansion",
             "Phenotypic Rate Elevated"),
  type   = c("dashed", "box", "box", "box"),
  colour = c("red", "cyan2", "orange2", "tomato")))

h2 <- c(panel_height(P2b), panel_height(P2a), legend_in)
Fig2 <- plot_grid(panel(P2b, "Invertebrates"), panel(P2a, "Vertebrates"), legend_massext,
                  ncol = 1, rel_heights = h2, labels = c("a", "b", ""), label_size = 12, 
                  label_fontface = "bold", align = "v", axis = "lr")
Fig2
ggsave("Fig2_MassExct.pdf", Fig2, width = 10, height = sum(h2))
