"""Mass extinctions: combine the PELT+CROPS and segmented-regression (SEG) results into the final tables
results_combined/Table_MassExtinction.csv and Summary_MassExtinction.csv, following the Decision Table.

Run from the repository root, after both method runs:
    python3 source_code/combine_results_massext.py              # Table_MassExtinction.csv
    python3 source_code/combine_results_massext.py --summary    # also Summary_MassExtinction.csv

Inputs : results_PELT+CROP/mass_extinction/results/Table_MassExtinction.csv
         results_segmented(SEG)/mass_extinction/results/Table_MassExtinction_segmented.csv
         results_combined/Decision_Table_MassExtinction.csv       best model for ONSET and OFFSET per series
         results_combined/Decision_exceptions_MassExtinction.csv  exceptions (rows marked * in the Decision Table)
         results_combined/Table_CladeOrigin.csv, Summary_CladeOrigin.csv   final clade-origin results (copied clades)
         ImageDigitization_pipeline/input_keys/massext_copy_from_origin.csv (submodule; DIGITISATION_REPO=<folder> points elsewhere)
Outputs: results_combined/Table_MassExtinction.csv (all columns, incl. both methods' values)
         results_combined/Summary_MassExtinction.csv (with --summary)
         If the Decision Table or the exceptions file is missing, a pre-filled Decision Table (defaults below) and an
         empty exceptions file are written; existing files are inputs and are never overwritten.

Series key: Clade + Extinction_event + Parameter + Data_type.
Decision labels: SEG, PELT (with * = see exception), no offset, no peak, na, too few points, clade origin, neither* .
Defaults of a pre-filled Decision Table: ONSET = SEG, OFFSET = PELT (general rule of the clade-origin analysis);
"clade origin" for the clades in massext_copy_from_origin.csv (their FINAL clade-origin values, full time series, are
copied; their onsets can be older than the extinction, giving a negative time after it); "na" (no published source)
and "too few points" (< 4 data points after the extinction): no values. "clade origin" can also be set for other
series (same trajectory as in the clade-origin analysis).
Rules:
 1. Values copied as estimated: point estimate (PE) + 95% Monte-Carlo range (Max = older, Min = younger).
 2. SEG onset + PELT offset with the PELT offset OLDER than the SEG onset -> SEG for both (noted).
 2b. A 95% MC older bound older than the extinction age (age jitter of the first post-extinction points) is capped at
     the extinction age (noted in Combination_notes and the PE basis).
 3. "no offset" (and an exception "none" for the offset): no decline after the first peak; the offset is set to the
    END of the time series (PE = Max = Min = Series_end_Ma), so that the plots draw the range to the end of the series.
    "na", "no peak", "too few points": no values.
 4. Exceptions override 1: Action = PELT / SEG / range (Older_Ma-Younger_Ma; PE = PE_Ma if given, else the first
    method in PE_check whose PE lies inside the range, else the midpoint) / none.
"""
import os, re, sys
import numpy as np, pandas as pd

PELT_F = os.path.join("results_PELT+CROP", "mass_extinction", "results", "Table_MassExtinction.csv")
SEG_F = os.path.join("results_segmented(SEG)", "mass_extinction", "results", "Table_MassExtinction_segmented.csv")
OUT = "results_combined"
DEC_F = os.path.join(OUT, "Decision_Table_MassExtinction.csv"); EXC_F = os.path.join(OUT, "Decision_exceptions_MassExtinction.csv")
ORIG_F = os.path.join(OUT, "Table_CladeOrigin.csv")
COPY_F = os.path.join(os.environ.get("DIGITISATION_REPO", "ImageDigitization_pipeline"), "input_keys", "massext_copy_from_origin.csv")

rd = lambda f: pd.read_csv(f, keep_default_na=False, na_values=[""], encoding="utf-8")
P, S, O = rd(PELT_F), rd(SEG_F), rd(ORIG_F)
COPY = pd.read_csv(COPY_F).Clade.tolist()
KEY = ["Clade", "Extinction_event", "Parameter", "Data_type"]
for d in (P, S, O):
    for c in [c for c in KEY + ["Reference"] if c in d.columns]: d[c] = d[c].fillna("").astype(str)
assert len(P) == len(S) and (P[KEY].values == S[KEY].values).all(), "PELT and SEG tables do not list the same series"
V = ["Onset_PE", "Onset_Max", "Onset_Min", "Offset_PE", "Offset_Max", "Offset_Min"]

def origin_row(r):
    o = O[(O.Clade == r.Clade) & (O.Parameter == r.Parameter) & (O.Reference == r.Reference)]
    if len(o) != 1: sys.exit(f"final clade-origin row not found once: {r.Clade} {r.Parameter} {r.Reference}")
    return o.iloc[0]

def status(r):
    if r.Clade in COPY and str(r.Reference).lower() != "na": return "clade origin"
    if pd.isna(r.get("Method")) or str(r.get("Method")) == "":
        fl = str(r.get("Flags", ""))
        return "too few points" if "Too few data points" in fl else "na"
    return ""

def agreement(p, s):
    """same classification as results_segmented(SEG)/R/compare_methods.py"""
    if any(pd.isna(x) for x in (p.Onset_PE, s.Onset_PE)): return "na (no data)"
    def ov(a0, a1, b0, b1):
        return not any(pd.isna(v) for v in (a0, a1, b0, b1)) and max(a0, b0) <= min(a1, b1) + 1e-9
    span = s.Series_start_Ma - s.Series_end_Ma; tol = max(5.0, 0.05 * span)
    on_ok = ov(s.Onset_Min, s.Onset_Max, p.Onset_Min, p.Onset_Max) or abs(s.Onset_PE - p.Onset_PE) <= tol
    of_ok = ov(s.Offset_Min, s.Offset_Max, p.Offset_Min, p.Offset_Max) or abs(s.Offset_PE - p.Offset_PE) <= tol
    return "consistent" if on_ok and of_ok else "onset differs" if of_ok else "offset differs" if on_ok else "onset and offset differ"

# ---------------------------------------------------------------- pre-filled Decision Table (only if missing)
if not os.path.exists(DEC_F):
    rows = []
    for i, p in P.iterrows():
        s = S.loc[i]; st = status(p)
        if st == "clade origin":
            o = origin_row(p); agr, on, off = o.Decision_Agreement, "clade origin", "clade origin"
            note = f"Final clade-origin decision: onset {o.Decision_onset}, offset {o.Decision_offset}"
        elif st in ("na", "too few points"):
            agr, on, off, note = ("na (no data)" if st == "na" else "too few points"), st, st, ""
        else:
            agr, on, off, note = agreement(p, s), "SEG", "PELT", "default (clade-origin general rule)"
        rows.append({"Group": p.Group, "Clade": p.Clade, "Extinction_event": p.Extinction_event, "Extinction_age_Ma": p.Extinction_age_Ma,
                     "Parameter": p.Parameter, "Data_type": p.Data_type, "Agreement": agr,
                     "ONSET_Best model": on, "OFFSET_Best model": off, "Prefill_note": note})
    pd.DataFrame(rows).to_csv(DEC_F, index=False, encoding="utf-8")
    print("pre-filled", DEC_F)
if not os.path.exists(EXC_F):
    pd.DataFrame(columns=["Clade", "Extinction_event", "Parameter", "Data_type", "Bound", "Action", "Older_Ma", "Younger_Ma",
                          "PE_Ma", "PE_check", "Note"]).to_csv(EXC_F, index=False, encoding="utf-8")
    print("empty template", EXC_F)

D = pd.read_csv(DEC_F, keep_default_na=False, encoding="utf-8-sig", dtype=str)
E = pd.read_csv(EXC_F, keep_default_na=False, encoding="utf-8-sig", dtype=str)
for d in (D, E):
    for c in d.columns: d[c] = d[c].astype(str).str.strip()

def vals(T, i, bound):
    r = T.loc[i]; b = bound.capitalize()
    return dict(PE=r[f"{b}_PE"], Max=r[f"{b}_Max"], Min=r[f"{b}_Min"])
def pick(method, i, bound):
    return dict(vals(P if method == "PELT" else S, i, bound), method=method, basis=f"{method} PE + 95% MC range")
def pi_range(ex, i, bound):
    older, younger = float(ex.Older_Ma), float(ex.Younger_Ma)
    if ex.PE_Ma:
        return dict(PE=float(ex.PE_Ma), Max=older, Min=younger, method="Manual range", basis="Manual range; PE = manual estimate")
    for m in [x for x in ex.PE_check.split(";") if x]:
        pe = vals(P if m == "PELT" else S, i, bound)["PE"]
        if pd.notna(pe) and younger - 1e-9 <= pe <= older + 1e-9:
            return dict(PE=pe, Max=older, Min=younger, method="Manual range", basis=f"Manual range; PE = {m} PE (inside range)")
    return dict(PE=(older + younger) / 2, Max=older, Min=younger, method="Manual range", basis="Manual range; PE = range midpoint (no method PE inside range)")

NONE = dict(PE=np.nan, Max=np.nan, Min=np.nan)
def series_end(r, label):
    """no offset: the offset is placed at the end (youngest point) of the analysed series"""
    e = float(r.Series_end_Ma)
    return dict(PE=e, Max=e, Min=e, method=label, basis=f"No offset (no decline): offset set to the end of the time series ({e:g} Ma; for plotting)")
rows = []
for i, r in P.iterrows():
    s = S.loc[i]
    dd = D[(D.Clade == r.Clade) & (D.Extinction_event == r.Extinction_event) & (D.Parameter == r.Parameter) & (D.Data_type == r.Data_type)]
    if len(dd) != 1: sys.exit(f"series not found once in the Decision Table: {r.Clade} | {r.Extinction_event} | {r.Parameter} | {r.Data_type}")
    dd = dd.iloc[0]; on_dec, off_dec = dd["ONSET_Best model"], dd["OFFSET_Best model"]
    ex = E[(E.Clade == r.Clade) & (E.Extinction_event == r.Extinction_event) & (E.Parameter == r.Parameter) & (E.Data_type == r.Data_type)]
    notes, ex_notes, out = [], [], {}
    rec = {c: r.get(c, "") for c in ["Group", "Clade", "Extinction_event", "Extinction_age_Ma", "Parameter", "Data_type", "Reference", "Figure",
                                     "Legend_line", "Key_data_file", "Data_source", "Error_type", "N_points", "Series_start_Ma", "Series_end_Ma"]}
    if on_dec == "clade origin" or off_dec == "clade origin":
        o = origin_row(r)
        rec.update(Decision_Agreement=o.Decision_Agreement, Decision_onset=f"clade origin ({o.Decision_onset})",
                   Decision_offset=f"clade origin ({o.Decision_offset})")
        for c in ["Onset_method", "Onset_PE", "Onset_Max", "Onset_Min", "Onset_PE_basis", "Offset_method", "Offset_PE", "Offset_Max",
                  "Offset_Min", "Offset_PE_basis", "Exception_note"]: rec[c] = o[c]
        cn = ("Copied from the FINAL clade-origin results (results_combined/Table_CladeOrigin.csv): " +
              ("the clade originates right after this extinction" if r.Clade in COPY else
               "same trajectory as in the clade-origin analysis (Decision Table)"))
        if pd.notna(o.Onset_PE) and float(o.Onset_PE) > float(r.Extinction_age_Ma) + 0.05:
            cn += f" | Note: the clade-origin onset ({float(o.Onset_PE):.1f} Ma) is older than the extinction age ({float(r.Extinction_age_Ma):g} Ma)"
        rec["Combination_notes"] = cn + (" | " + str(o.Combination_notes) if pd.notna(o.Combination_notes) and str(o.Combination_notes) else "")
        for c in V: rec[f"PELT_{c}"] = o[f"PELT_{c}"]
        for c in V: rec[f"SEG_{c}"] = o[f"SEG_{c}"]
        rec["PELT_n_shifts"] = r.get("Elbow_n_cpts", ""); rec["SEG_k"] = s.get("BIC_k", ""); rec["SEG_curve"] = s.get("Curve_type", "")
        rec["PELT_Flags"] = o.PELT_Flags; rec["SEG_Flags"] = o.SEG_Flags
        rows.append(rec); continue
    if on_dec in ("na", "no peak", "too few points"):
        out["onset"] = dict(NONE, method=on_dec, basis=""); out["offset"] = dict(NONE, method=on_dec, basis="")
        notes.append({"na": "No published source (key: na)", "no peak": "No clear peak: treated as na",
                      "too few points": "Too few data points after the extinction (< 4): not analysed"}[on_dec])
    else:
        for bound, dec in (("onset", on_dec), ("offset", off_dec)):
            e = ex[ex.Bound == bound]
            if len(e) > 1: sys.exit(f"duplicate exception {r.Clade} {r.Extinction_event} {r.Parameter} {bound}")
            if len(e):
                e = e.iloc[0]; ex_notes.append(f"{bound}: {e.Note}")
                if e.Action in ("PELT", "SEG"): out[bound] = dict(pick(e.Action, i, bound), basis=f"{e.Action} PE + 95% MC range (Exception)")
                elif e.Action == "range": out[bound] = pi_range(e, i, bound)
                elif e.Action == "none": out[bound] = dict(NONE, method="none", basis="") if bound == "onset" else series_end(r, "none")
                else: sys.exit(f"unknown exception action {e.Action}")
            elif dec.rstrip("*") in ("SEG", "PELT"):
                if dec.endswith("*"): sys.exit(f"{r.Clade} {r.Extinction_event} {r.Parameter} {bound} marked {dec} but no exception given")
                out[bound] = pick(dec, i, bound)
            elif dec == "no offset" and bound == "offset":
                out[bound] = series_end(r, "no offset")
                notes.append("No offset: no decline after the first peak (plateau or monotonic increase); offset set to the end of the time series")
            else:
                sys.exit(f"{r.Clade} {r.Extinction_event} {r.Parameter} {bound}: decision '{dec}' needs an exception entry")
        on, off = out["onset"], out["offset"]
        if (on["method"], off["method"]) == ("SEG", "PELT") and "exception" not in on["basis"] + off["basis"] \
                and pd.notna(off["PE"]) and pd.notna(on["PE"]) and off["PE"] > on["PE"]:
            out["offset"] = dict(pick("SEG", i, "offset"), basis="SEG PE + 95% MC range (PELT offset older than SEG onset)")
            notes.append(f"PELT offset ({off['PE']:.1f} Ma) older than SEG onset ({on['PE']:.1f} Ma): SEG used for both")
        # the series start at the extinction: a 95% MC bound older than the extinction age only reflects the age jitter
        # of the first post-extinction points -> capped at the extinction age (as onsets are capped at the clade origin)
        ea = float(r.Extinction_age_Ma)
        for b in ("onset", "offset"):
            v = out[b]
            if v["method"] in ("SEG", "PELT", "Manual range") and pd.notna(v["Max"]) and v["Max"] > ea + 0.05:
                notes.append(f"{b.capitalize()} 95% MC older bound {v['Max']:.1f} Ma capped at the extinction age ({ea:g} Ma): age jitter of the first post-extinction points")
                out[b] = dict(v, Max=ea, basis=v["basis"] + "; older bound capped at the extinction age")
        on, off = out["onset"], out["offset"]
        if pd.notna(on["PE"]) and pd.notna(off["PE"]) and off["PE"] > on["PE"] + 1e-9: notes.append("CHECK: offset PE older than onset PE")
        if pd.notna(on["Min"]) and pd.notna(off["Max"]) and on["Min"] < off["Max"] - 1e-9:
            notes.append(f"Onset range younger bound ({on['Min']:.1f}) is younger than the offset range older bound ({off['Max']:.1f}): ranges overlap")
        for b, v in out.items():
            if v["method"] in ("SEG", "PELT") and pd.notna(v["PE"]) and not (v["Min"] - 0.05 <= v["PE"] <= v["Max"] + 0.05):
                notes.append(f"{b.capitalize()} PE lies outside its 95% MC range (method flags: PE not robust)")
        # review hints (no effect on the values)
        pf, sf = str(r.get("Flags", "")), str(s.get("Flags", ""))
        if "right-censored" in pf and "right-censored" in sf: notes.append("Hint: both methods put the offset at the series end (no decline after the first peak) -> consider 'no offset'")
        if "left-censored" in pf and ("left-censored" in sf or "starts in decline" in sf.lower()):
            notes.append("Hint: both methods put the onset at the first data point after the extinction")
        if "Few data points after the extinction" in pf: notes.append("Hint: few data points after the extinction")
    rec.update(Decision_Agreement=dd.Agreement, Decision_onset=on_dec, Decision_offset=off_dec)
    for b in ("onset", "offset"):
        v = out[b]; B = b.capitalize()
        rec.update({f"{B}_method": v["method"], f"{B}_PE": v["PE"], f"{B}_Max": v["Max"], f"{B}_Min": v["Min"], f"{B}_PE_basis": v["basis"]})
    rec["Exception_note"] = " | ".join(ex_notes); rec["Combination_notes"] = " | ".join(notes)
    for c in V: rec[f"PELT_{c}"] = r[c]
    for c in V: rec[f"SEG_{c}"] = s[c]
    rec["PELT_n_shifts"] = r.get("Elbow_n_cpts", ""); rec["SEG_k"] = s.get("BIC_k", ""); rec["SEG_curve"] = s.get("Curve_type", "")
    rec["PELT_Flags"] = r.get("Flags", ""); rec["SEG_Flags"] = s.get("Flags", "")
    rows.append(rec)

T = pd.DataFrame(rows)
for c in [c for c in T.columns if re.search(r"_(PE|Max|Min)$", c)]: T[c] = pd.to_numeric(T[c], errors="coerce").round(1)
# time from the extinction to the onset / offset (Myr)
T.insert(T.columns.get_loc("Onset_PE_basis") + 1, "Onset_after_extinction_Myr", (T.Extinction_age_Ma - T.Onset_PE).round(1))
T.insert(T.columns.get_loc("Offset_PE_basis") + 1, "Offset_after_extinction_Myr", (T.Extinction_age_Ma - T.Offset_PE).round(1))
# no offset: the offset value is only the end of the series (for plotting), not a dated decline -> no time to offset
T.loc[T.Offset_method.isin(["no offset", "none"]), "Offset_after_extinction_Myr"] = np.nan
T.to_csv(os.path.join(OUT, "Table_MassExtinction.csv"), index=False, na_rep="", encoding="utf-8")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 100); pd.set_option("display.max_colwidth", 30)
print(T[["Clade", "Extinction_age_Ma", "Parameter", "Data_type", "Decision_onset", "Onset_PE", "Onset_Max", "Onset_Min",
         "Decision_offset", "Offset_PE", "Offset_Max", "Offset_Min"]].to_string(index=False))

if "--summary" in sys.argv:
    # Layout and wording follow results_combined/Summary_CladeOrigin.csv; Extinction_event and Extinction_age_Ma are
    # added after Clade. Data_type / Inference_Method are copied from the clade-origin Summary for series with the same
    # Clade + Parameter + source reference; the sources new to this analysis are listed in NEW_META.
    import unicodedata
    # numbers as in Summary_CladeOrigin.csv: one decimal, no trailing ".0"
    f = lambda x: "na" if pd.isna(x) else (lambda v: f"{0.0 if abs(v) < 0.05 else v:.1f}".rstrip("0").rstrip("."))(round(float(x), 1))
    norm = lambda x: "".join(c for c in unicodedata.normalize("NFKD", str(x)) if not unicodedata.combining(c)).strip().lower()
    META, CO_SUM = {}, {}
    sum_f = os.path.join(OUT, "Summary_CladeOrigin.csv")
    if os.path.exists(sum_f):
        so = pd.read_csv(sum_f, dtype=str, keep_default_na=False, encoding="utf-8-sig")
        for _, q in so.iterrows():
            ref = norm(q.Source.split(",")[0])
            META[(q.Clade, q.Parameter, ref)] = (q.get("Data_type", ""), q.get("Inference_Method", ""))
            CO_SUM[(q.Clade, q.Parameter, ref)] = q          # copied clades: Source and Note as in the clade-origin Summary
    # sources new to the mass-extinction analysis: Data_type and Inference_Method (vocabulary as in
    # Summary_CladeOrigin.csv). Crinoidea_Diplobathrida disparity: the series used is the SOV curve (Cole & Hopkins
    # 2021 Fig. 2B), so SOV.
    NEW_META = {("Ammonoidea", "phenotypic_disparity", "whalen_2020"): ("LMM (shell)", "MPD"),
                ("Ammonoidea", "taxonomic_diversity", "whalen_2020"): ("Occur", "RT"),
                ("Crinoidea_Diplobathrida", "phenotypic_disparity", "cole_hopkins_2021"): ("Discrete (WB)", "SOV"),
                ("Crinoidea_Diplobathrida", "taxonomic_diversity", "cole_hopkins_2021"): ("Occur", "RT"),
                ("Eurypterida", "phenotypic_disparity", "lamsdell_selden_2016"): ("Discrete (WB)", "SOV"),
                ("Eurypterida", "taxonomic_diversity", "lamsdell_selden_2016"): ("Occur", "RT"),
                ("Placentalia", "phenotypic_disparity", "halliday_goswami_2016"): ("Discrete & Cont", "MPD")}
    for cl in ("Brachiopoda_Athyridida", "Brachiopoda_Rhynchonellida", "Brachiopoda_Spiniferida", "Brachiopoda_Terebratuliida"):
        NEW_META[(cl, "phenotypic_disparity", "guo_2024")] = ("Discrete", "SOV")
        NEW_META[(cl, "taxonomic_diversity", "guo_2024")] = ("Occur", "RT")
    for k, v in NEW_META.items(): META.setdefault(k, v)
    # Crinoidea_Eucladida diversity (Wright 2017 Fig. S10) classed as RT
    META[("Crinoidea_Eucladida", "taxonomic_diversity", "wright_2017_si")] = ("Phylo & Occur", "RT")
    sm = []
    for _, r in T.iterrows():
        ref = str(r.Reference) if pd.notna(r.Reference) else "na"
        src = ", ".join(x for x in (ref, str(r.Figure)) if x and x != "nan")
        om, fm = str(r.Onset_method), str(r.Offset_method)
        is_na = om == "na" or ref.lower() == "na"           # no published source
        note = []
        if om == "too few points": note.append("Too few data points after the extinction (not analysed)")
        if om == "no peak": note.append("No clear peak (becomes na)")
        if str(r.Decision_onset).startswith("clade origin") and r.Clade not in COPY:   # copied clades: no note
            note.append("Clade-origin results (same trajectory as in the clade-origin analysis)")
        if fm in ("no offset", "none"): note.append("No offset (no decline): offset set to the end of the time series")
        if "Manual range" in (om, fm): note.append("Adjusted range (see exception)")   # copied clades too
        elif "(Exception)" in str(r.Onset_PE_basis) + str(r.Offset_PE_basis): note.append("Exception on method choice")
        if "SEG used for both" in str(r.Combination_notes): note.append("PELT offset older than SEG onset: SEG for both")
        if "capped at the extinction age" in str(r.Combination_notes): note.append("Older bound capped at the extinction age")
        dt, im = META.get((r.Clade, r.Parameter, norm(ref)), ("", ""))
        if r.Clade in COPY and not is_na:                    # copied clades: everything as in the clade-origin Summary,
            q = CO_SUM[(r.Clade, r.Parameter, norm(ref))]    # incl. its Source and Note
            src, note = q.Source, [q.Note] if q.Note else []
        if is_na: dt, im = "na", "na"
        sm.append(dict(Group=r.Group, Clade=r.Clade, Extinction_event=r.Extinction_event, Extinction_age_Ma=f"{r.Extinction_age_Ma:g}",
                       Parameter=r.Parameter, Data_type=dt, Inference_Method=im, Source="na" if is_na else src,
                       Onset_PE=f(r.Onset_PE), Onset_Max=f(r.Onset_Max), Onset_Min=f(r.Onset_Min),
                       Onset_method="na" if (om == "na" or ref.lower() == "na") else om,
                       Offset_PE=f(r.Offset_PE), Offset_Max=f(r.Offset_Max), Offset_Min=f(r.Offset_Min),
                       Offset_method="na" if (fm == "na" or ref.lower() == "na") else ("none" if fm == "none" else fm),
                       Onset_after_extinction_Myr=f(r.Onset_after_extinction_Myr), Offset_after_extinction_Myr=f(r.Offset_after_extinction_Myr),
                       Note="; ".join(note)))
    pd.DataFrame(sm).to_csv(os.path.join(OUT, "Summary_MassExtinction.csv"), index=False, encoding="utf-8")
    print("written", os.path.join(OUT, "Summary_MassExtinction.csv"))
