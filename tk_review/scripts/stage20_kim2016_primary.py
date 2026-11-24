#!/usr/bin/env python3
"""Kim 2016 Table 4 read at source: three rat renal clearances (Stage 20).

    python3 -I scripts/stage20_kim2016_primary.py

Writes db/primary_2026/kim2016_table4_rat.csv, db/renal_only_v2.csv and
db/gfr_sensitivity.csv. Idempotent.

WHY THIS STAGE EXISTS. Stage 18 raised an open provenance flag: the rat
PFHxS renal clearance of 0.76 mL/kg/day, one of only three rat inputs to
equation (5), had been recorded in db/vd_clearance.csv from ATSDR and EPA
secondary sources with no pmid or doi, and the Kim 2016 abstract reports
plasma and tissue distribution without mentioning urine. The PDF has now been
read.

FLAG RESOLVED, FAVOURABLY. Kim 2016 collected urine and feces, defines
CL_R = Xu(t)/AUC(0-inf) in its Methods, and reports CL_R in Table 4 for all
three compounds in both sexes after IV dosing (it is NA after oral dosing).
The 0.76 is genuine: PFHxS, male, IV. Equation (5)'s PFHxS row stands.

AND STAGE 18 WAS WRONG ABOUT PFOS. Stage 18 concluded that no rat PFOS renal
clearance existed anywhere and admitted PFOS to equation (5) only through the
loose bound CL_R <= CL_tot. In fact Kim 2016 Table 4 reports a MEASURED rat
PFOS renal clearance, 1.92 mL/kg/day in the male. It was in the very paper
whose provenance was being questioned, which had been judged from its
abstract rather than read. The bound is superseded by the measurement, and
the measurement is 4.8-fold below the tightest bound Stage 18 could offer.

A FOURTH RAT PFOA RENAL CLEARANCE, AND IT DISAGREES. Kim's male rat PFOA
CL_R is 4.26 mL/kg/day, against 18.2 (Kemper via OEHHA), 21.2 (Ohmori Fig. 2
digitised by Han) and 46.08 (Kudo 2002). The four span 10.8-fold. That is a
material widening of a disagreement the paper already carried at 2.5-fold.

THE TABLE IS INTERNALLY VALIDATED. Kim's Results text states urinary
percentages of the IV dose independently of Table 4, and CL_R/CL reproduces
them: PFOA male 8.99% against a stated 9.0%, PFHxS female 28.02% against a
stated 28.02%, PFHxS male 8.44% against 8.26%, PFOS female 17.5% against
17.7%. Four independent checks, worst discrepancy 0.2 percentage points. The
values below are therefore read correctly despite Table 4 being rotated and
scrambled by text extraction; each was taken from the page image and
confirmed against the text layer.

GFR SENSITIVITY IS INCLUDED HERE because it cannot be separated from these
numbers. The project holds TWO rat GFR values differing 1.82-fold -- 14400
mL/kg/day (OEHHA/Han) and 7920 (Ohmori/Shim) -- and Harada 2005 implies a
third near 12200. The paper uses 14400 throughout without comment. GFR
multiplies the glomerular term and divides the tubular one, so the choice
moves 0.26 log units from one to the other and changes nothing else.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

GFR_HUMAN = 2570.0
GFR_RAT = {"OEHHA/Han": 14400.0, "Ohmori/Shim": 7920.0, "Harada 2005": 12200.0}

# Kim 2016 Table 4, read from the page image and confirmed against the text
# layer. CL_R is reported for IV only; oral rows carry NA.
KIM_TABLE4 = [
    # chem,    sex,      route, cl_renal, sd,   cl_tot, cl_tot_sd
    ("PFOA",   "male",   "IV",   4.26,   0.03,   47.39, 3.40),
    ("PFOA",   "female", "IV", 323.60,   7.94,  612.84, 12.54),
    ("PFOS",   "male",   "IV",   1.92,   0.01,    9.24, 0.37),
    ("PFOS",   "female", "IV",   1.72,   0.03,    9.82, 0.21),
    ("PFHxS",  "male",   "IV",   0.76,   0.01,    9.01, 0.05),
    ("PFHxS",  "female", "IV",  63.87,   1.96,  227.93, 6.73),
]
# Urinary % of IV dose stated in Kim's Results text, independent of Table 4.
KIM_NARRATIVE_PCT = {
    ("PFOA", "male"): 9.0,
    ("PFHxS", "male"): 8.26,
    ("PFHxS", "female"): 28.02,
    ("PFOS", "female"): 17.7,
}

# Human male renal clearances already established by the project.
HUMAN_CLR = {
    "PFHpA": [(0.674, "Fujii 2015 Table 2")],
    "PFOA":  [(0.030, "Harada 2005 Sec 3.3"), (0.044, "Fujii 2015 Table 2")],
    "PFHxS": [(0.027, "Zhang 2013 Table 2")],
    "PFOS":  [(0.012, "Harada 2005 Table 5"), (0.015, "Harada 2005 Sec 3.3"),
              (0.024, "Zhang 2013 median"), (0.037, "Zhang 2013 mean")],
}
# Male rat renal clearances, now including Kim read at source.
RAT_CLR = {
    "PFHpA": [(983.1, "Ohmori 2003 Fig 2, digitised")],
    "PFOA":  [(4.26, "Kim 2016 Table 4, IV, READ AT SOURCE"),
              (18.2, "Kemper 2003 via OEHHA PHG A6.4"),
              (21.2, "Ohmori 2003 Fig 2, digitised by Han 2012"),
              (46.08, "Kudo 2002 Table 3")],
    "PFHxS": [(0.76, "Kim 2016 Table 4, IV, READ AT SOURCE")],
    "PFOS":  [(1.92, "Kim 2016 Table 4, IV, READ AT SOURCE")],
}
NC = {"PFHpA": 7, "PFOA": 8, "PFHxS": 6, "PFOS": 8}


def fu_table():
    out = {}
    for r in csv.DictReader((HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r["experiment"] == "Ryu 2024" and r["fu"]:
            sp = r["species"].split(" (")[0]
            if sp in ("human", "rat"):
                out[(sp, r["chemical"])] = float(r["fu"])
    return out


def main():
    fu = fu_table()
    bar = "=" * 78

    print(bar)
    print("KIM 2016 TABLE 4, READ AT SOURCE  (IV dosing; oral rows are NA)")
    print(bar)
    print(f"  {'cmpd':7}{'sex':8}{'CL_R':>9}{'CL_tot':>9}{'renal %':>9}"
          f"{'stated':>8}  check")
    prim = []
    for chem, sex, route, clr, clr_sd, clt, clt_sd in KIM_TABLE4:
        pct = clr / clt * 100
        stated = KIM_NARRATIVE_PCT.get((chem, sex))
        if stated is None:
            mark, shown = "not stated separately in the text", "--"
        else:
            mark, shown = f"agrees to {abs(pct-stated):.2f} pt", f"{stated:.2f}"
        print(f"  {chem:7}{sex:8}{clr:>9.2f}{clt:>9.2f}{pct:>9.2f}"
              f"{shown:>8}  {mark}")
        prim.append(dict(chemical=chem, sex=sex, route=route,
                         cl_renal_mL_kg_day=clr, cl_renal_sem=clr_sd,
                         cl_total_mL_kg_day=clt, cl_total_sem=clt_sd,
                         renal_pct_of_total=round(pct, 2),
                         narrative_pct=(stated if stated else ""),
                         source="Kim 2016 Table 4, Food Chem Toxicol "
                                "97:243-255, doi:10.1016/j.fct.2016.09.017"))
    print()
    print("  Four independent checks against Kim's own Results text; worst")
    print("  discrepancy 0.2 percentage points. The reading is sound.")
    print()

    print(bar)
    print("THE SEX DIFFERENCE IS IN THE TUBULE, AND IT CHANGES T's SIGN")
    print(bar)
    print(f"  {'cmpd':7}{'T male':>10}{'T female':>11}{'F/M':>9}  verdict")
    sexrows = []
    for chem in ("PFOA", "PFOS", "PFHxS"):
        m = [r for r in KIM_TABLE4 if r[0] == chem and r[1] == "male"][0]
        f_ = [r for r in KIM_TABLE4 if r[0] == chem and r[1] == "female"][0]
        filt = fu[("rat", chem)] * GFR_RAT["OEHHA/Han"]
        tm, tf = m[3] / filt, f_[3] / filt
        v = ("T crosses 1: male REABSORBS, female SECRETES"
             if tm < 1 < tf else
             "both below 1: net reabsorption in both sexes"
             if max(tm, tf) < 1 else "both above 1")
        print(f"  {chem:7}{tm:>10.4f}{tf:>11.4f}{tf/tm:>9.1f}  {v}")
        sexrows.append((chem, tm, tf))
    print()
    print("  Kim measured both sexes in ONE study, ONE strain, ONE assay, and")
    print("  reports no significant sex difference in V_d for any of the")
    print("  three. So the 76- to 84-fold sex difference in renal clearance")
    print("  for PFOA and PFHxS is not a distribution effect; it sits in the")
    print("  tubular term. PFOS, which has no sex difference in renal")
    print("  clearance, is the control that shows the contrast is real.")
    print()
    print("  THIS BEARS ON REABSORPTION vs SECRETION. For PFOA and PFHxS the")
    print("  male tubular factor is far below 1 and the female far above it,")
    print("  so within one study the same compound is net-reabsorbed in one")
    print("  sex and net-secreted in the other. A secretion-only account")
    print("  would need females to secrete 76-fold harder; a reabsorption")
    print("  account needs males to reabsorb what females do not, which is")
    print("  what androgen-regulated apical transport predicts, and which is")
    print("  the mechanism Kim's own discussion invokes.")
    print()

    print(bar)
    print("EQUATION (5), CORRECTED AND EXTENDED  (male, rat/human)")
    print(bar)
    print("  PFOS is now MEASURED, not bounded. Stage 18 gave T_rat <= 0.59")
    print(f"  at best; the measurement is "
          f"{1.92/(fu[('rat','PFOS')]*GFR_RAT['OEHHA/Han']):.4f}, "
          f"{0.5887/(1.92/(fu[('rat','PFOS')]*GFR_RAT['OEHHA/Han'])):.1f}-fold"
          f" below that bound.")
    print()
    rows = []
    for gname, grat in GFR_RAT.items():
        for chem in ("PFHpA", "PFOA", "PFHxS", "PFOS"):
            fh, fr = fu[("human", chem)], fu[("rat", chem)]
            fb, fg = fr / fh, grat / GFR_HUMAN
            for vh, sh in HUMAN_CLR[chem]:
                th = vh / (fh * GFR_HUMAN)
                for vr, sr in RAT_CLR[chem]:
                    tr = vr / (fr * grat)
                    ratio, ft = vr / vh, tr / th
                    lg = math.log10(ratio)
                    pb = math.log10(fb) / lg * 100
                    pg = math.log10(fg) / lg * 100
                    pt = math.log10(ft) / lg * 100
                    rows.append(dict(
                        gfr_rat_source=gname, gfr_rat=grat, chemical=chem,
                        n_carbons=NC[chem], cl_renal_human=vh,
                        human_source=sh, cl_renal_rat=vr, rat_source=sr,
                        cl_renal_ratio=round(ratio, 2), T_human=round(th, 6),
                        T_rat=round(tr, 5), factor_binding=round(fb, 4),
                        factor_glomerular=round(fg, 4),
                        factor_tubule=round(ft, 3),
                        pct_binding=round(pb, 1),
                        pct_glomerular=round(pg, 1),
                        pct_tubule=round(pt, 1)))
    base = [r for r in rows if r["gfr_rat_source"] == "OEHHA/Han"]
    print(f"  {'cmpd':7}{'CL_R rat':>10}{'ratio':>9}{'T rat':>9}"
          f"{'tubule':>10}{'tub %':>8}  rat source")
    for chem in ("PFHpA", "PFOA", "PFHxS", "PFOS"):
        sub = [r for r in base if r["chemical"] == chem]
        seen = set()
        for r in sub:
            if r["rat_source"] in seen:
                continue
            seen.add(r["rat_source"])
            print(f"  {chem:7}{r['cl_renal_rat']:>10.2f}"
                  f"{r['cl_renal_ratio']:>9.0f}{r['T_rat']:>9.4f}"
                  f"{r['factor_tubule']:>10.1f}{r['pct_tubule']:>8.0f}"
                  f"  {r['rat_source'][:34]}")
    print()
    pfoa = [r for r in base if r["chemical"] == "PFOA"]
    tfs = [r["factor_tubule"] for r in pfoa]
    pts = [r["pct_tubule"] for r in pfoa]
    print(f"  PFOA's four rat sources span {max(tfs)/min(tfs):.1f}-fold in the")
    print(f"  tubular FACTOR ({min(tfs):.1f} to {max(tfs):.1f}) but only "
          f"{min(pts):.0f}-{max(pts):.0f}% in its SHARE.")
    print("  The magnitude is badly determined; the conclusion is not.")
    print()

    print(bar)
    print("RAT GFR SENSITIVITY: 1.82-FOLD, AND WHAT IT MOVES")
    print(bar)
    print(f"  {'GFR source':14}{'GFR':>8}   mean share over compounds/sources")
    print(f"  {'':14}{'':8}   {'binding':>9}{'glomerular':>12}{'tubule':>9}")
    gfr_rows = []
    for gname, grat in GFR_RAT.items():
        sub = [r for r in rows if r["gfr_rat_source"] == gname]
        mb = sum(r["pct_binding"] for r in sub) / len(sub)
        mg = sum(r["pct_glomerular"] for r in sub) / len(sub)
        mt = sum(r["pct_tubule"] for r in sub) / len(sub)
        print(f"  {gname:14}{grat:>8.0f}   {mb:>9.1f}{mg:>12.1f}{mt:>9.1f}")
        gfr_rows.append(dict(gfr_rat_source=gname, gfr_rat=grat,
                             n_combinations=len(sub),
                             mean_pct_binding=round(mb, 1),
                             mean_pct_glomerular=round(mg, 1),
                             mean_pct_tubule=round(mt, 1)))
    print()
    print("  GFR trades directly against the tubule and touches nothing else:")
    print("  it multiplies the glomerular term and divides the tubular one by")
    print("  the same factor, so log10(1.82) = 0.26 units move from one to")
    print("  the other. Binding is untouched, which is why the paper's")
    print("  central claim is insensitive to the choice -- but the glomerular")
    print("  and tubular shares must be reported as a range, not a number.")
    print()

    out = []
    p = HERE / "db" / "primary_2026" / "kim2016_table4_rat.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prim[0].keys()))
        w.writeheader()
        for r in sorted(prim, key=lambda x: (x["chemical"], x["sex"])):
            w.writerow(r)
    out.append((p, len(prim)))

    p = HERE / "db" / "renal_only_v2.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["gfr_rat_source"],
                                             x["chemical"],
                                             x["cl_renal_human"],
                                             x["cl_renal_rat"])):
            w.writerow(r)
    out.append((p, len(rows)))

    p = HERE / "db" / "gfr_sensitivity.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gfr_rows[0].keys()))
        w.writeheader()
        for r in sorted(gfr_rows, key=lambda x: x["gfr_rat"]):
            w.writerow(r)
    out.append((p, len(gfr_rows)))

    for p, n in out:
        print(f"wrote {p.relative_to(HERE)}  ({n} rows)")


if __name__ == "__main__":
    main()
