#!/usr/bin/env python3
"""Renal-only decomposition: no residual term (Stage 17).

    python3 -I scripts/stage17_renal_only_decomposition.py

Writes db/renal_only_decomposition.csv. Idempotent.

WHY THIS IS DIFFERENT FROM STAGE 15. That decomposition used TOTAL clearance,
so its fourth factor was a residual absorbing both tubular handling and any
species difference in non-renal elimination. With direct renal clearances now
in hand for both species, the renal clearance ratio decomposes EXACTLY into
three measured terms and no residual at all:

    CL_R,rat     fu_rat      GFR_rat      T_rat
    ---------  = -------  x  ---------  x  -----
    CL_R,human   fu_human    GFR_human    T_human
                 BINDING     GLOMERULAR   TUBULE

Here T is genuinely the tubular factor, not a catch-all. Nothing is inferred
by subtraction.

THE PRICE IS COVERAGE. Only three compounds have a DIRECT renal clearance in
male human and male rat: PFHpA, PFOA and PFHxS. PFNA has one in the human and
only in the FEMALE rat, and pairing across sexes is not admissible when sex
changes the sign of the tubular term in this species.

PFDA CANNOT BE ADDED, AND THE ATTEMPT IS RECORDED BECAUSE IT FAILED TWICE.
Ohmori 2003 Figure 2 was re-examined at 300 dpi: the PFDA bars for both sexes
sit flat on the baseline, carrying only their significance letters, exactly
as Stage 10 found. Ohmori's own Figure 3 caption states that PFNA-male,
PFDA-male and PFDA-female "are not distinguished on the graph". Its published
regression, CL_R = 0.48*CL_tot + 15.63 (r2 = 0.981), cannot rescue them
either: at PFDA's total clearance of 5.2 mL/kg/day it predicts a renal
clearance of 18.1, which exceeds the total by 3.5-fold and is impossible. The
intercept is fitted across totals up to 3000 and swamps small values. PFDA's
renal clearance is not recoverable from this paper by any route.
"""
import csv, math, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
GFR = {"human": 2570.0, "rat": 14400.0}

# Direct renal clearances, mL/kg/day, MALE in both species.
# Several rat PFOA sources disagree, so all three are carried.
HUMAN_CLR = {
    "PFHpA": [(0.674, "Fujii 2015 Table 2, urinary, n=10")],
    "PFOA":  [(0.030, "Harada 2005 Section 3.3, n=20"),
              (0.044, "Fujii 2015 Table 2, urinary, n=10")],
    "PFHxS": [(0.027, "Zhang 2013 Table 2, male+older female, n=64")],
}
RAT_CLR = {
    "PFHpA": [(983.1, "Ohmori 2003 Fig 2, digitised")],
    "PFOA":  [(18.2, "Kemper 2003 via OEHHA PHG A6.4"),
              (21.2, "Ohmori 2003 Fig 2, digitised by Han 2012"),
              (46.08, "Kudo 2002 Table 3")],
    "PFHxS": [(0.76, "Kim 2016")],
}
# Renal share of total clearance, from SAME-STUDY pairs where possible.
F_REN = {
    ("human", "PFHpA"): (100.0, "Fujii 2015: urinary 0.674 of total 0.674"),
    ("human", "PFOA"):  (45.8, "Fujii 2015: 0.044 of 0.096"),
    ("human", "PFOA2"): (25.1, "Harada 2005 Table 5, male 61 kg"),
    ("human", "PFHxS"): (50.9, "CROSS-STUDY: Zhang 0.027 of Abraham 0.053"),
    ("rat", "PFHpA"):   (61.3, "Ohmori 2003: CL_R 983.1 of CL_tot 1603.9"),
    ("rat", "PFOA"):    (91.2, "Kudo 2002 / Ohmori 2003: 46.08 of 50.5"),
    ("rat", "PFOA2"):   (42.0, "Ohmori Fig 2 digitised: 21.2 of 50.5"),
    ("rat", "PFHxS"):   (8.4, "Kim 2016: 0.76 of 9.01"),
}
NC = {"PFHpA": 7, "PFOA": 8, "PFHxS": 6}


def fu_table():
    out = {}
    for r in csv.DictReader(
            (HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r["experiment"] == "Ryu 2024" and r["fu"]:
            sp = r["species"].split(" (")[0]
            if sp in ("human", "rat"):
                out[(sp, r["chemical"])] = float(r["fu"])
    return out


def main():
    fu = fu_table()
    bar = "=" * 78
    print(bar)
    print("TUBULAR FACTORS FROM DIRECT RENAL CLEARANCE, MALE, BOTH SPECIES")
    print(bar)
    print(f"  {'cmpd':7}{'species':8}{'CL_R':>9}{'fu':>10}{'filtered':>10}"
          f"{'T':>9}  source")
    T = {}
    for chem in ("PFHpA", "PFOA", "PFHxS"):
        for sp, table in (("human", HUMAN_CLR), ("rat", RAT_CLR)):
            f = fu[(sp, chem)]
            filt = f * GFR[sp]
            for v, src in table[chem]:
                t = v / filt
                T.setdefault((sp, chem), []).append((t, v, src))
                print(f"  {chem:7}{sp:8}{v:>9.3f}{f:>10.5f}{filt:>10.2f}"
                      f"{t:>9.4f}  {src}")
        print()

    print(bar)
    print("THE RENAL CLEARANCE RATIO, DECOMPOSED EXACTLY  (rat / human)")
    print(bar)
    print("  Three measured terms, no residual. T is genuinely tubular.")
    print()
    print(f"  {'cmpd':7}{'CL_R ratio':>12}{'= binding':>11}{'x glom':>9}"
          f"{'x tubule':>11}   share of log ratio (%)")
    print(f"  {'':7}{'':12}{'':11}{'':9}{'':11}   "
          f"{'bind':>6}{'glom':>7}{'TUBULE':>9}")
    rows = []
    for chem in ("PFHpA", "PFOA", "PFHxS"):
        fh, fr = fu[("human", chem)], fu[("rat", chem)]
        fb, fg = fr / fh, GFR["rat"] / GFR["human"]
        for th, vh, sh in T[("human", chem)]:
            for tr, vr, sr in T[("rat", chem)]:
                ratio = vr / vh
                ft = tr / th
                tot = math.log10(ratio)
                pb, pg, pt = (math.log10(fb) / tot * 100,
                              math.log10(fg) / tot * 100,
                              math.log10(ft) / tot * 100)
                rows.append(dict(chem=chem, ratio=ratio, fb=fb, fg=fg, ft=ft,
                                 pb=pb, pg=pg, pt=pt, th=th, tr=tr,
                                 vh=vh, vr=vr, sh=sh, sr=sr))
        # print one line per human source, using the rat median source
        for r in [x for x in rows if x["chem"] == chem][-len(T[("rat", chem)]):]:
            print(f"  {chem:7}{r['ratio']:>12.0f}{r['fb']:>11.2f}"
                  f"{r['fg']:>9.2f}{r['ft']:>11.1f}   "
                  f"{r['pb']:>6.0f}{r['pg']:>7.0f}{r['pt']:>9.0f}")
    print()
    n = len(rows)
    print(f"  mean over all {n} source combinations:"
          f"   binding {sum(r['pb'] for r in rows)/n:>5.0f}%"
          f"   glomerular {sum(r['pg'] for r in rows)/n:>5.0f}%"
          f"   TUBULE {sum(r['pt'] for r in rows)/n:>5.0f}%")
    print()
    print("  Compare Stage 15, where the same partition computed from TOTAL")
    print("  clearance gave binding 8% and a residual of 65%. Removing the")
    print("  non-renal routes does not rescue binding.")
    print()

    print(bar)
    print("HOW MUCH OF THE HALF-LIFE THIS ACCOUNTS FOR")
    print(bar)
    print(f"  {'species':8}{'cmpd':8}{'renal % of total':>18}  basis")
    for (sp, key), (pct, src) in sorted(F_REN.items()):
        print(f"  {sp:8}{key:8}{pct:>18.1f}  {src}")
    print()
    print("  The renal route is a minority of human elimination for PFOA and")
    print("  PFHxS and the whole of it for PFHpA, while in the rat it is the")
    print("  majority for PFOA and PFHpA and a small minority for PFHxS.")
    print("  So the renal decomposition above explains all of the PFHpA")
    print("  difference, roughly half of the PFOA difference, and only part")
    print("  of the PFHxS difference; the remainder is a difference in route")
    print("  partition rather than in renal handling.")
    print()

    p = HERE / "db" / "renal_only_decomposition.csv"
    with p.open("w", newline="") as fh_:
        w = csv.writer(fh_)
        w.writerow(["chemical", "n_carbons", "cl_renal_human",
                    "human_source", "cl_renal_rat", "rat_source",
                    "cl_renal_ratio", "T_human", "T_rat", "factor_binding",
                    "factor_glomerular", "factor_tubule", "pct_binding",
                    "pct_glomerular", "pct_tubule"])
        for r in rows:
            w.writerow([r["chem"], NC[r["chem"]], r["vh"], r["sh"], r["vr"],
                        r["sr"], round(r["ratio"], 1), round(r["th"], 5),
                        round(r["tr"], 4), round(r["fb"], 4),
                        round(r["fg"], 4), round(r["ft"], 3),
                        round(r["pb"], 1), round(r["pg"], 1),
                        round(r["pt"], 1)])
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} source combinations)")


if __name__ == "__main__":
    main()
