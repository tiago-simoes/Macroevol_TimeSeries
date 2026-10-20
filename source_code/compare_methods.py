"""Side-by-side comparison of first-peak onset/offset: PELT+CROPS (results_PELT+CROP/clade_origin/results/) vs segmented regression (results_segmented(SEG)/clade_origin/results/).
Usage (from the repository root; the three paths below are the defaults):
  python3 source_code/compare_methods.py "results_PELT+CROP/clade_origin/results/Table_CladeOrigin.csv" \
          "results_segmented(SEG)/clade_origin/results/Table_CladeOrigin_segmented.csv" "results_segmented(SEG)/clade_origin/results/Comparison_PELT_vs_segmented.csv"
(Mass extinctions: the comparison is part of results_combined/Table_MassExtinction.csv, source_code/combine_results_massext.py.)
Differences are segmented minus PELT, in Myr (positive = segmented estimate is OLDER).
'Overlap' = the two 95% Monte-Carlo ranges intersect."""
import sys, pandas as pd, numpy as np
pelt_f = sys.argv[1] if len(sys.argv) > 1 else "results_PELT+CROP/clade_origin/results/Table_CladeOrigin.csv"
seg_f = sys.argv[2] if len(sys.argv) > 2 else "results_segmented(SEG)/clade_origin/results/Table_CladeOrigin_segmented.csv"
out_f = sys.argv[3] if len(sys.argv) > 3 else "results_segmented(SEG)/clade_origin/results/Comparison_PELT_vs_segmented.csv"
rd = lambda f: pd.read_csv(f, keep_default_na=False, na_values=[""])
P, S = rd(pelt_f), rd(seg_f)
keys = ["Clade", "Parameter", "Data_type", "Reference"]
for d in (P, S):
    for k in keys: d[k] = d[k].fillna("").astype(str)
M = S.merge(P, on=keys, how="left", suffixes=("_SEG", "_PELT"))
f = lambda x: "na" if pd.isna(x) else f"{x:.1f}"
def rng(lo, hi): return "na" if pd.isna(lo) else f"{f(hi)}-{f(lo)}"
def overlap(a_min, a_max, b_min, b_max):
    if any(pd.isna(v) for v in (a_min, a_max, b_min, b_max)): return "na"
    return "yes" if max(a_min, b_min) <= min(a_max, b_max) + 1e-9 else "no"
def short_pelt(fl):
    fl = "" if pd.isna(fl) else fl; n = []
    if "left-censored" in fl: n.append("onset = series start")
    if "right-censored" in fl: n.append("offset = series end")
    if "Stepwise rise" in fl: n.append("stepwise rise")
    if "Model order fixed" in fl: n.append("k fixed")
    if "elbow unstable" in fl: n.append("elbow unstable")
    return "; ".join(n)
def short_seg(fl):
    fl = "" if pd.isna(fl) else fl; n = []
    if "Series starts in decline" in fl: n.append("starts in decline")
    elif "onset = oldest data point" in fl.lower(): n.append("onset = series start")
    if "offset = youngest data point" in fl.lower(): n.append("offset = series end")
    if "No breakpoint" in fl: n.append("k = 0")
    if "No significant rising" in fl: n.append("no significant rise")
    if "Cumulative" in fl: n.append("cumulative curve (slope-change rule)")
    if "steepens" in fl: n.append("rise steepens later (see Main_rise)")
    if "weakly constrained" in fl: n.append("onset weakly constrained")
    if "model order unstable" in fl: n.append("k unstable in MC")
    if "Onset point estimate not robust" in fl: n.append("onset PE not robust")
    if "Offset point estimate not robust" in fl: n.append("offset PE not robust")
    return "; ".join(n)
rows = []
for _, r in M.iterrows():
    base = dict(Group=r.get("Group_SEG", ""), Clade=r.Clade, Parameter=r.Parameter, Data_type=r.Data_type,
                Source=", ".join(x for x in (r.Reference, str(r.get("Figure_SEG", "") or "")) if x and x != "nan"))
    if pd.isna(r.get("Onset_PE_SEG")):
        rows.append(dict(base, Agreement="na (no data)")); continue
    don = r.Onset_PE_SEG - r.Onset_PE_PELT if not pd.isna(r.get("Onset_PE_PELT")) else np.nan
    dof = r.Offset_PE_SEG - r.Offset_PE_PELT if not pd.isna(r.get("Offset_PE_PELT")) else np.nan
    oo = overlap(r.Onset_Min_SEG, r.Onset_Max_SEG, r.Onset_Min_PELT, r.Onset_Max_PELT)
    fo = overlap(r.Offset_Min_SEG, r.Offset_Max_SEG, r.Offset_Min_PELT, r.Offset_Max_PELT)
    span = r.Series_start_Ma_SEG - r.Series_end_Ma_SEG
    tol = max(5.0, 0.05 * span)                      # PE difference regarded as small
    on_ok = oo == "yes" or (not pd.isna(don) and abs(don) <= tol)
    of_ok = fo == "yes" or (not pd.isna(dof) and abs(dof) <= tol)
    agr = "consistent" if on_ok and of_ok else "onset differs" if of_ok else "offset differs" if on_ok else "onset and offset differ"
    rows.append(dict(base,
        PELT_Onset=f(r.Onset_PE_PELT), PELT_Onset_95MC=rng(r.Onset_Min_PELT, r.Onset_Max_PELT),
        SEG_Onset=f(r.Onset_PE_SEG), SEG_Onset_95MC=rng(r.Onset_Min_SEG, r.Onset_Max_SEG), SEG_Onset_MC_median=f(r.Onset_Med_SEG),
        SEG_Onset_modelCI=rng(r.Onset_CI95_younger, r.Onset_CI95_older), dOnset_Myr=f(don), Onset_MC_overlap=oo,
        PELT_Offset=f(r.Offset_PE_PELT), PELT_Offset_95MC=rng(r.Offset_Min_PELT, r.Offset_Max_PELT),
        SEG_Offset=f(r.Offset_PE_SEG), SEG_Offset_95MC=rng(r.Offset_Min_SEG, r.Offset_Max_SEG), SEG_Offset_MC_median=f(r.Offset_Med_SEG),
        SEG_Offset_modelCI=rng(r.Offset_CI95_younger, r.Offset_CI95_older), dOffset_Myr=f(dof), Offset_MC_overlap=fo,
        SEG_k=int(r.BIC_k), SEG_curve=r.Curve_type, SEG_Main_rise=f(r.get("Main_rise_PE")),
        PELT_n_shifts=r.get("Elbow_n_cpts"), PELT_notes=short_pelt(r.get("Flags_PELT")), SEG_notes=short_seg(r.get("Flags_SEG")),
        Agreement=agr))
out = pd.DataFrame(rows)
out.to_csv(out_f, index=False)
print(out["Agreement"].value_counts().to_string())
print(out[["Clade", "Parameter", "Data_type", "PELT_Onset", "SEG_Onset", "PELT_Offset", "SEG_Offset", "Agreement"]].to_string(index=False))

# figure: PELT vs segmented point estimates with 95% MC ranges
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    d = M[~M.Onset_PE_SEG.isna() & ~M.Onset_PE_PELT.isna()]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.6))
    for a, w, ttl in ((ax[0], "Onset", "Onset (Ma)"), (ax[1], "Offset", "Offset (Ma)")):
        x, y = d[f"{w}_PE_PELT"], d[f"{w}_PE_SEG"]
        c0 = lambda v: np.clip(v.fillna(0).values, 0, None)
        a.errorbar(x, y, xerr=[c0(x - d[f"{w}_Min_PELT"]), c0(d[f"{w}_Max_PELT"] - x)],
                   yerr=[c0(y - d[f"{w}_Min_SEG"]), c0(d[f"{w}_Max_SEG"] - y)], fmt="o", ms=4, lw=.8,
                   color="#2c6e9b", ecolor="#9bbad1", capsize=0)
        lim = [0, max(x.max(), y.max()) * 1.03]
        a.plot(lim, lim, color="grey", lw=1, ls="--"); a.set_xlim(lim[::-1]); a.set_ylim(lim[::-1])
        a.set_xlabel(f"PELT+CROPS {ttl}"); a.set_ylabel(f"Segmented regression {ttl}"); a.set_title(ttl.split(" ")[0])
        for _, r in d.iterrows():
            if abs(r[f"{w}_PE_SEG"] - r[f"{w}_PE_PELT"]) > max(10, 0.1 * (r.Series_start_Ma_SEG - r.Series_end_Ma_SEG)):
                a.annotate(f"{r.Clade} {r.Parameter.split('_')[-1]}{(' ' + r.Data_type) if r.Data_type else ''}"[:38],
                           (r[f"{w}_PE_PELT"], r[f"{w}_PE_SEG"]), fontsize=6, alpha=.8, xytext=(3, 3), textcoords="offset points")
    fig.suptitle("First-peak onset/offset: PELT+CROPS vs segmented regression (95% MC ranges; dashed = 1:1)", fontsize=10)
    fig.tight_layout(); fig.savefig(out_f.replace(".csv", ".png"), dpi=150)
except Exception as e:
    print("figure skipped:", e)
