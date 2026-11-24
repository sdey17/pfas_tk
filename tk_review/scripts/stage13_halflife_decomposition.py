#!/usr/bin/env python3
"""Decompose the human/rat half-life ratio into its four physical factors.

    python3 -I scripts/stage13_halflife_decomposition.py

Writes db/halflife_decomposition.csv. Idempotent.

THE QUESTION THIS ANSWERS. The project has established that binding differs
little between species while the tubule differs greatly. That is an argument
about mechanism. This script turns it into an argument about OUTCOME: of the
observed difference in elimination half-life between human and rat, how much
is attributable to distribution, how much to binding, how much to glomerular
filtration, and how much to the tubule?

THE ALGEBRA. For a one-compartment system,

        t_half  =  ln2 * Vd / CL

so for one compound in two species the half-life ratio is

        t_h / t_r  =  (Vd_h / Vd_r)  *  (CL_r / CL_h)

and substituting CL = fu * GFR * T,

        t_h / t_r  =  (Vd_h/Vd_r) * (fu_r/fu_h) * (GFR_r/GFR_h) * (T_r/T_h)
                      \__________/  \__________/  \____________/  \_______/
                       DISTRIBUTION    BINDING      GLOMERULAR     TUBULE

Each factor is a multiplicative contribution to the half-life difference, and
on a log scale they ADD, so the percentage of the log ratio carried by each
is well defined.

WHAT T MEANS HERE, AND THE HONEST CAVEAT. Both clearances are TOTAL, not
renal, so T is the apparent tubular factor CL_tot/(fu*GFR). The fourth factor
is therefore a RESIDUAL: it contains the tubule plus any species difference
in non-renal elimination. In the human, 58-61% of PFOA elimination is not
renal, so the residual is not purely tubular and must not be described as if
it were. The identity is exact BY CONSTRUCTION -- the residual absorbs
whatever the first three factors do not explain -- so the result is not a
test of the identity but a partition of explanatory burden among terms that
are each independently measured except the last.

SOURCES. Human: Abraham 2024, one male volunteer, labelled single oral doses,
Vd and CL both MEASURED, and the identity ln2*Vd/CL closes at a consistent
implied body weight of 82.9 kg across all 15 compounds -- an internal
consistency check the rat data cannot match. Rat: male, measured, one study
per compound, named per row; the rat side is unavoidably multi-source and
that is a limitation, not a choice.
"""
import csv, math, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
GFR = {"human": 2570.0, "rat": 14400.0}

# Human: Abraham 2024, male volunteer. Vd mL/kg, CL mL/kg/day, both measured.
HUMAN = {
    "PFBA":  (132.0, 21.89), "PFPeA": (110.0, 146.6), "PFHxA": (119.0, 56.89),
    "PFHpA": (128.0, 0.584), "PFOA":  (121.0, 0.0417), "PFNA": (124.0, 0.0659),
    "PFDA":  (184.0, 0.148), "PFBS":  (137.0, 1.877), "PFHxS": (125.0, 0.0530),
    "PFOS":  (152.0, 0.0870),
}
# Rat: male, MEASURED Vd and CL. One study per compound, named.
# Ohmori's PFOA row in vd_clearance.csv carries the ATSDR transcription error
# (Vd 196.2 is the paper's PFHpA value, CL 135 is 2.7x the paper's 50.5), so
# the rat PFOA entry uses Kim 2016 instead -- see db/primary_2026/.
RAT = {
    "PFHpA": (196.2, 1604.0, "Ohmori 2003, Wistar, IV, 2-compartment Vss"),
    "PFHxA": (930.0, 2784.0, "Chengelis 2009, SD, IV, Vd-beta"),
    "PFOA":  (112.0, 47.39,  "Kim 2016, SD, IV, non-compartmental"),
    "PFOS":  (383.0, 9.24,   "Kim 2016, SD, IV, non-compartmental"),
    "PFHxS": (224.0, 4.73,   "Huang 2019 NTP, SD, IV, 2-compartment"),
    "PFBS":  (188.0, 828.0,  "Huang 2019 NTP, SD, IV, 2-compartment"),
}
NC = {"PFBA": 4, "PFPeA": 5, "PFHxA": 6, "PFHpA": 7, "PFOA": 8, "PFNA": 9,
      "PFDA": 10, "PFBS": 4, "PFHxS": 6, "PFOS": 8}


def fu_table():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] == "Ryu 2024" and r["fu"]:
            sp = r["species"].split(" (")[0]
            if sp in ("human", "rat"):
                out[(sp, r["chemical"])] = float(r["fu"])
    return out


def main():
    fu = fu_table()
    rows = []
    for chem in sorted(set(HUMAN) & set(RAT), key=lambda c: NC[c]):
        vh, clh = HUMAN[chem]
        vr, clr, src = RAT[chem]
        fh, fr = fu.get(("human", chem)), fu.get(("rat", chem))
        if not (fh and fr):
            continue
        th = math.log(2) * vh / clh          # days
        tr = math.log(2) * vr / clr
        ratio = th / tr
        f_dist = vh / vr
        f_bind = fr / fh
        f_glom = GFR["rat"] / GFR["human"]
        f_tub = (clr / clh) / (f_bind * f_glom)   # residual, exact by construction
        check = f_dist * f_bind * f_glom * f_tub
        tot = math.log10(ratio)
        rows.append(dict(
            chemical=chem, nc=NC[chem], t_h=th, t_r=tr, ratio=ratio,
            f_dist=f_dist, f_bind=f_bind, f_glom=f_glom, f_tub=f_tub,
            check=check, vh=vh, vr=vr, clh=clh, clr=clr, fh=fh, fr=fr,
            src=src,
            p_dist=math.log10(f_dist) / tot * 100 if tot else float("nan"),
            p_bind=math.log10(f_bind) / tot * 100 if tot else float("nan"),
            p_glom=math.log10(f_glom) / tot * 100 if tot else float("nan"),
            p_tub=math.log10(f_tub) / tot * 100 if tot else float("nan")))

    bar = "=" * 78
    print(bar)
    print("HUMAN / RAT HALF-LIFE RATIO, DECOMPOSED  (male, measured Vd and CL)")
    print(bar)
    print(f"  {'chem':7}{'nC':>3}{'t_h (d)':>10}{'t_r (d)':>9}{'ratio':>10}"
          f"{'  =  dist':>10}{'bind':>8}{'glom':>7}{'TUBULE':>10}")
    for r in rows:
        print(f"  {r['chemical']:7}{r['nc']:>3}{r['t_h']:>10.1f}{r['t_r']:>9.2f}"
              f"{r['ratio']:>10.0f}{r['f_dist']:>10.2f}{r['f_bind']:>8.2f}"
              f"{r['f_glom']:>7.2f}{r['f_tub']:>10.1f}")
    print()
    print("  identity check (product should equal the ratio):")
    for r in rows:
        ok = abs(r["check"] - r["ratio"]) / r["ratio"] < 1e-9
        print(f"    {r['chemical']:7} {r['check']:>12.1f} vs {r['ratio']:>12.1f}"
              f"  {'exact' if ok else 'MISMATCH'}")
    print()
    print(bar)
    print("SHARE OF THE log10 HALF-LIFE RATIO CARRIED BY EACH FACTOR (%)")
    print(bar)
    print(f"  {'chem':7}{'distribution':>14}{'binding':>10}{'glomerular':>12}"
          f"{'TUBULE':>10}")
    for r in rows:
        print(f"  {r['chemical']:7}{r['p_dist']:>14.0f}{r['p_bind']:>10.0f}"
              f"{r['p_glom']:>12.0f}{r['p_tub']:>10.0f}")
    n = len(rows)
    for k, lbl in (("p_dist", "distribution"), ("p_bind", "binding"),
                   ("p_glom", "glomerular"), ("p_tub", "TUBULE")):
        print(f"  {'mean':7}{sum(r[k] for r in rows)/n:>14.0f}" if k == "p_dist"
              else f"{sum(r[k] for r in rows)/n:>10.0f}", end="")
    print()
    print()
    print("  A NEGATIVE share means the factor works AGAINST the observed")
    print("  difference -- the human would be eliminating FASTER than the rat")
    print("  on that factor alone, and the others must overcome it.")
    print()
    print(bar)
    print("IS THE BINDING TERM EVEN CAPABLE OF EXPLAINING THIS?")
    print(bar)
    rs = sorted(((c, fu[("rat", c)] / fu[("human", c)])
                 for c in {x for sp, x in fu if sp == "human"}
                 if ("rat", c) in fu), key=lambda x: -x[1])
    lo, hi = rs[-1][1], rs[0][1]
    print(f"  rat/human fu ratio across all {len(rs)} compounds measured in")
    print(f"  one laboratory (Ryu 2024a): {lo:.2f} to {hi:.2f}")
    print(f"      widest:  {rs[0][0]} {hi:.2f}      narrowest: {rs[-1][0]} {lo:.2f}")
    print(f"  half-life ratios to be explained: "
          f"{min(r['ratio'] for r in rows):.0f} to "
          f"{max(r['ratio'] for r in rows):.0f}")
    print()
    print("  THE BINDING TERM IS BOUNDED AT ABOUT TWO-FOLD BY DIRECT")
    print("  MEASUREMENT. It is arithmetically incapable of accounting for a")
    print("  thousand-fold difference in half-life, whichever direction it")
    print("  points. No refinement of binding data can change this; the")
    print("  bound is on the quantity itself, not on our knowledge of it.")
    print()

    print(bar)
    print("SENSITIVITY: DOES THIS DEPEND ON ABRAHAM'S SINGLE SUBJECT?")
    print(bar)
    ALT = {"PFOA": (430.0, 0.095, "Chiu 2022 Bayesian hierarchical"),
           "PFOS": (320.0, 0.066, "Chiu 2022 Bayesian hierarchical"),
           "PFHxS": (290.0, 0.025, "Chiu 2022 Bayesian hierarchical")}
    print(f"  {'chem':7}{'human source':32}{'t_h (y)':>9}{'ratio':>8}"
          f"{'bind%':>8}{'TUBULE%':>9}")
    for chem in ("PFOA", "PFOS", "PFHxS"):
        vr, clr, _ = RAT[chem]
        fh, fr = fu[("human", chem)], fu[("rat", chem)]
        for vh, clh, lbl in ((HUMAN[chem][0], HUMAN[chem][1],
                              "Abraham 2024 (n=1, measured)"),
                             ALT[chem]):
            th = math.log(2) * vh / clh
            ratio = th / (math.log(2) * vr / clr)
            tot = math.log10(ratio)
            fb = fr / fh
            ft = (clr / clh) / (fb * GFR["rat"] / GFR["human"])
            print(f"  {chem:7}{lbl:32}{th/365:>9.1f}{ratio:>8.0f}"
                  f"{math.log10(fb)/tot*100:>8.0f}"
                  f"{math.log10(ft)/tot*100:>9.0f}")
        print()
    print("  The partition shifts but the conclusion does not: binding under")
    print("  10% on every combination, the tubular residual 57-80%.")
    print()

    p = HERE / "db" / "halflife_decomposition.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "t_half_human_d", "t_half_rat_d",
                    "halflife_ratio_h_over_r", "factor_distribution",
                    "factor_binding", "factor_glomerular",
                    "factor_tubule_residual", "identity_check",
                    "pct_log_distribution", "pct_log_binding",
                    "pct_log_glomerular", "pct_log_tubule",
                    "vd_human_mL_kg", "vd_rat_mL_kg", "cl_human_mL_kg_d",
                    "cl_rat_mL_kg_d", "fu_human", "fu_rat", "rat_source"])
        for r in rows:
            w.writerow([r["chemical"], r["nc"], round(r["t_h"], 2),
                        round(r["t_r"], 3), round(r["ratio"], 1),
                        round(r["f_dist"], 4), round(r["f_bind"], 4),
                        round(r["f_glom"], 4), round(r["f_tub"], 3),
                        round(r["check"], 1), round(r["p_dist"], 1),
                        round(r["p_bind"], 1), round(r["p_glom"], 1),
                        round(r["p_tub"], 1), r["vh"], r["vr"], r["clh"],
                        r["clr"], r["fh"], r["fr"], r["src"]])
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} compounds)")


if __name__ == "__main__":
    main()
