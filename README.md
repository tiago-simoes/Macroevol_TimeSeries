# Macroevolutionary time series: change-point analysis

Code and data to reproduce the change-point analysis of macroevolutionary time series (taxonomic diversity,
phenotypic disparity and phenotypic evolutionary rates): the timing of the first shift after the origin of each
clade and after each mass extinction, estimated with two methods (PELT+CROPS and segmented regression), the final
result tables and the figures.

The time series were digitised from published figures, or converted from raw data files, with the
[ImageDigitization_pipeline](https://github.com/tiago-simoes/ImageDigitization_pipeline) repository. It is included
here as a git submodule (`ImageDigitization_pipeline/`), pinned to the commit whose keys and data tables were
analysed; the scripts read the keys and the series directly from it.

## Repository layout

```
ImageDigitization_pipeline/  submodule (digitisation repository), used here for:
  input_keys/                keys of the two analyses (one row per clade x parameter [x extinction]) and the clades
                             whose clade-origin results are used for their extinction (massext_copy_from_origin.csv)
  outputs/data_tables/       clade_origin/<Clade>_<Parameter>[_<data_type>].csv (series of the clade-origin key) and
                             mass_extinction/<series>__<event>.csv (series cropped to the time after the extinction)
data_digitized/digitisation/ model-order override of the PELT+CROPS runs (model_order_overrides.csv)
source_code/                 R and Python scripts (below)
results_PELT+CROP/           PELT+CROPS results: clade_origin/ and mass_extinction/, each with results/ (master
                             table) and figures/ (one figure per series)
results_segmented(SEG)/      segmented-regression results, same layout; clade_origin/results/ also holds the
                             PELT vs SEG comparison
results_combined/            final tables (Table_*.csv, Summary_*.csv) and the mass-extinction Decision Table and
                             exceptions
plots/                       plot tables and the plotting script
```

Folder names contain `+`, `(`, `)` and file names contain `&`: quote them in shell commands.

### Scripts (`source_code/`)

| Script | Role |
|---|---|
| `changepoint_pipeline.R` | series reading, regridding to a 1-Myr grid, PELT helpers (`changepoint` package) |
| `crops_firstpeak.R` | CROPS penalty scan, kneedle elbow, first-peak rule, Monte Carlo conditioned on the model order |
| `mc_uncertainty.R` | Monte-Carlo replicate generation and quantile summaries |
| `run_key_crops.R` | PELT+CROPS runner: every row of a key -> master table, per-series tables and figures |
| `seg_firstpeak.R` | segmented-regression engine: breakpoints, BIC model order, onset/offset rules, Monte Carlo |
| `run_key_segmented.R` | segmented-regression runner (same inputs and outputs as `run_key_crops.R`) |
| `compare_methods.py` | PELT vs SEG comparison of the clade-origin results (table + figure) |
| `combine_results_massext.py` | final mass-extinction tables from both method tables and the Decision Table |

## Requirements

* R (tested with 4.3.3) with `changepoint` (tested 2.3) for PELT+CROPS; the segmented regression uses base R only.
* For the figures: R packages `ggplot2` (3.4.4), `dplyr` (1.1.4), `deeptime` (1.0.1), `ggpattern` (1.0.1),
  `cowplot` (1.1.3).
* Python 3 with `pandas`, `numpy` and `matplotlib`.

## Getting the code and data

```
git clone --recursive https://github.com/tiago-simoes/Macroevol_TimeSeries.git
```

(after a plain `git clone`, run `git submodule update --init`). The submodule brings the whole digitisation
repository, about 330 MB, most of it the source PDFs. To use another copy of it, set `DIGITISATION_REPO=<folder>`.

## Reproducing the results

Run from the repository root:

```
# 1. Clade origins
Rscript source_code/run_key_crops.R "ImageDigitization_pipeline/input_keys/key_origin_invert&vert.csv" Table_CladeOrigin
Rscript source_code/run_key_segmented.R "ImageDigitization_pipeline/input_keys/key_origin_invert&vert.csv" Table_CladeOrigin_segmented
python3 source_code/compare_methods.py

# 2. Mass extinctions (after step 1: clades in massext_copy_from_origin.csv copy their clade-origin rows)
Rscript source_code/run_key_crops.R "ImageDigitization_pipeline/input_keys/key_massext_invert&vert.csv" Table_MassExtinction
Rscript source_code/run_key_segmented.R "ImageDigitization_pipeline/input_keys/key_massext_invert&vert.csv" Table_MassExtinction_segmented
python3 source_code/combine_results_massext.py --summary

# 3. Figures
Rscript plots/EventsRange.R
```

The method runs write `results_PELT+CROP/<analysis>/results/Table_*.csv` and
`results_segmented(SEG)/<analysis>/results/Table_*_segmented.csv`, one figure per series in `figures/`, and per-series
tables (CROPS solutions, SEG models and segments, the 1000 Monte-Carlo replicates, cached fits) in
`results/clade_tables/` (not tracked). Each series uses fixed random seeds, so a re-run reproduces the tables exactly.
A full segmented-regression run takes about 1-1.5 h per analysis on one core; `ONLY=<series id regex>` runs a subset
in parallel processes, and `REUSE_PE=1 REUSE_MC=1` then assembles the master table from the cached fits.
`REUSE_MC=1` does the same for PELT+CROPS.

## Methods

### PELT+CROPS (`run_key_crops.R`)

* Each series is averaged onto a regular 1-Myr grid. PELT (Pruned Exact Linear Time, changes in mean) is run for every
  penalty in [0.5, 200] (CROPS; minimum segment length 3); the number of shifts is the elbow of the cost vs number
  of changepoints curve (kneedle rule). `model_order_overrides.csv` fixes the number of shifts of one series
  (Dinosauria_early PyRate diversity: 3 instead of the elbow's 2).
* First-peak rule: onset = first upward shift (start of the rising phase that ends in the first peak; series start if
  the series begins at peak level); offset = first decline after the first peak (series end if none). Later peaks
  are ignored.
* Monte Carlo (1000 replicates): values drawn from split-normal distributions built from the published error bars
  (95% intervals -> half-width / 1.96; SE bars as they are; Poisson error for counts without error bars), ages
  jittered by the calibration SD of the digitisation (`Age_cal_sd`). Each replicate is refitted with the model order
  of the point estimate; `Onset_Max/Med/Min` etc. are the 97.5/50/2.5% quantiles (Max = older bound). Each
  replicate's own elbow is recorded (`MC_n_cpts_freq`).

### Segmented regression (`run_key_segmented.R`)

* Continuous piecewise-linear regression of the value on time, y = b0 + b1 t + sum_j g_j (t - psi_j)+, with
  breakpoints estimated by Muggeo's (2003) iterative linearisation (the algorithm of the R package `segmented`,
  re-implemented in base R), started from an exhaustive grid of knots (k <= 2), a greedy forward search, a backward
  pass and random restarts; the minimum-RSS solution is kept, with at least 3 points per segment.
* Number of breakpoints: BIC over k = 0..6, keeping the smallest k within 2 BIC units of the minimum. Densely traced
  curves are binned for the BIC (bin width min(5 Myr, span / 20)); breakpoints are estimated on the full data.
* Curve types. *Hump* (rates, disparity, hump-shaped diversity): segments are rising or declining if their slope is
  significant and changes the value by >= 10% of the fitted range; onset = start of the first rising segment,
  offset = start of the first declining segment. *Cumulative* (monotonic diversity curves such as lineages through
  time, detected automatically): onset = first significant acceleration, offset = first deceleration.
* Uncertainty: delta-method 95% CI of each breakpoint, and the same Monte Carlo as for PELT+CROPS (refitted with the
  point-estimate k; the first 100 replicates also re-select k by BIC, `MC_k_bic_freq`).
* The optimiser was checked against the Python package `piecewise-regression` (Pilgrim 2021) on 8 series
  (`results_segmented(SEG)/clade_origin/results/validation_vs_piecewise_regression.csv`).

References: Killick et al. 2012 (JASA 107:1590-1598, PELT); Haynes et al. 2017 (JCGS 26:134-143, CROPS);
Muggeo 2003 (Stat Med 22:3055-3071); Kass & Raftery 1995 (JASA 90:773-795).

### Mass extinctions

Each series is cropped to the time after the extinction (`Age_Ma` <= extinction age of the key; node-based raw rates
re-binned in 5-Myr bins anchored at the extinction age), so onset and offset are the first shifts after the event,
with the same rules. Series with fewer than 4 points after the extinction are not analysed. The five clades in
`ImageDigitization_pipeline/input_keys/massext_copy_from_origin.csv` use their full clade-origin series and final
clade-origin results.

### Final tables (`results_combined/`)

* Onset and offset of each series are taken from the method chosen per series: by default the segmented-regression
  onset (earlier shifts, closer to the origin or extinction) and the PELT offset (tail of the decline). Exceptions
  replace a method value by the other method or by a manual range where neither method captures the first peak; the
  point estimate of a manual range is the method estimate when it lies inside the range, otherwise the midpoint.
  Series without a decline ("no offset") carry the end of the series as offset, so that the plots draw the range to
  the end of the series.
* `Table_CladeOrigin.csv` and `Summary_CladeOrigin.csv` are the final clade-origin results, assembled from the two
  method tables with these rules and then curated by hand. For every bound, `Table_CladeOrigin.csv` gives the
  decision (`Decision_onset`, `Decision_offset`), the method used (`Onset_method`, `Offset_method`), how the point
  estimate was obtained (`*_PE_basis`), the reason for each exception (`Exception_note`) and both methods' original
  values and flags. `Summary_CladeOrigin.csv` adds the data type and inference method of each source study and the
  clade-origin ages and references, and the time from the clade origin to the onset
  (`MinAge_Onset_after_CladeOrigin_Myr` = `Clade_origin_Min` - `Onset_Max`; `MaxAge_...` = `Clade_origin_Max` -
  `Onset_Min`).
* `Table_MassExtinction.csv` and `Summary_MassExtinction.csv` are built by `combine_results_massext.py` from the
  method tables, `Decision_Table_MassExtinction.csv` and `Decision_exceptions_MassExtinction.csv`; they add the time
  from the extinction to the onset and offset (`*_after_extinction_Myr`; negative where a reused clade-origin onset
  predates the extinction). An older 95% bound older than the extinction age (age jitter of the first points after
  the event) is capped at the extinction age.
* For Dinosauria_early and Theropoda_non-avian phenotypic rates, the final offsets (215.8 [218.8-212.8] and 208.5 Ma)
  are the PELT offsets of the series before the node-age cutoffs (`Series_start_Ma` 263.1 and 245.7 Ma); the method
  tables in this repository are re-runs on the cut series (82 and 69 points), which give 215.8 [218.8-213.8] and
  209.5 Ma. The onsets of both series are manual ranges.
* Ages are in Ma; `_Max` = older bound, `_Min` = younger bound of the 95% range.

Abbreviations in `Data_type`, `Inference_Method` and `Clade_origin_method`: BD, birth-death model; BDD, birth-death
diffusion model; BT, multivariate rate shifts with reversible-jump MCMC in BayesTraits; CalBayes, time-calibrated
Bayesian inference with relaxed clocks (node or tip dating, molecular and/or morphological data); Clock, clock rates
from Bayesian phylodynamics; Discrete, discrete characters; FBDD, fossilized birth-death diffusion model; GMM,
geometric morphometrics; LMM, linear morphometrics; LTT, lineages through time from a phylogenetic tree; MLR, maximum
local richness; MPD, mean pairwise disparity; Occur, fossil occurrences; PCM, phylogenetic comparative methods;
Phylo, time-calibrated phylogenetic data; PyRate, time-variable birth-death models in PyRate; RT, range-through
counts from raw occurrence data (no subsampling); RT (CONOP), RT with quantitative biostratigraphic correlation
(CONOP); SFBD, skyline fossilized birth-death model; SOR, sum of ranges; SOV, sum of variances; TESS-CoMET,
time-dependent birth-death model; WB, whole body.

### Figures (`plots/`)

The plot tables `Data_Origin_*.csv` and `Data_MassExct_*.csv` hold the onset and offset ranges of the two Summary
tables, one row per clade (a second series of one parameter goes in `Div2_` / `Rate2_` columns). `EventsRange.R`
draws `Fig1_Origins.pdf` (clade origins) and `Fig2_MassExct.pdf` (mass extinctions; both written to `plots/`, not
tracked): solid bars from the younger onset bound to the older offset bound, hatched onset and offset ranges.
