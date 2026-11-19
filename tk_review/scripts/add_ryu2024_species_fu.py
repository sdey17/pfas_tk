#!/usr/bin/env python3
"""Ryu 2024: measured plasma fu in human, rat and mouse. Closes the rat gap.

    python3 -I scripts/add_ryu2024_species_fu.py

Appends to db/plasma_unbound_fraction.csv and prints the recomputation of
OEHHA PHG Table A6.4 with measured values for all three species.

    Ryu S, Burchett W, Zhang S, Modaresi SMS, Agudelo Areiza J, Kaye E,
    Fischer FC, Slitt AL. "Species-Specific Unbound Fraction Differences in
    Highly Bound PFAS: A Comparative Study across Human, Rat, and Mouse
    Plasma and Albumin." Toxics 2024;12(4):253.
    doi:10.3390/toxics12040253  PMID 38668476, PMC11054487.

Presaturation equilibrium dialysis, 14 PFAS at 5 uM, in plasma of human, rat,
mouse CD-1 and mouse C57BL/6, and in the three serum albumins.

THIS IS THE MEASUREMENT THE ROADMAP HAS BEEN ASKING FOR. The rat value for
PFOA exists and is 0.00298.

AND IT CHANGES THE ANSWER IN BOTH DIRECTIONS AT ONCE.

OEHHA's assumed fu of 0.02 is wrong in VALUE -- eight-fold too high for the
rat and PFOA, six-fold too high for the human -- so the reabsorption
percentages the PHG prints are badly wrong. Recomputed with measured values:

                        OEHHA assumed      measured (Ryu 2024)
    human                   99.88%              99.05%
    male rat                93.68%              57.6%
    male mouse              97.01%              95.8%

The male rat does not reabsorb 94% of filtered PFOA. It reabsorbs about 58%.

But the assumption is approximately right in STRUCTURE for this pair, which
the earlier sensitivity analysis could not know. Rat and human plasma fu for
PFOA differ only 1.22-fold (0.00298 against 0.00245), and the paper's own
equivalence testing puts human and rat within two-fold with a mean absolute
fold difference of 1.32. So holding fu equal between rat and human, which is
what OEHHA does, happens to be close to true FOR THIS PAIR AND THIS COMPOUND.

The decomposition of the observed 303-fold difference therefore survives:

    assumed fu     303x  =  5.6x glomerular  x  54x tubular
    measured fu    303x  =  6.8x glomerular  x  44x tubular

THE MOUSE IS WHERE THE ASSUMPTION FAILS STRUCTURALLY, not the rat. Mouse
plasma fu for PFOA is 0.0142-0.0152, which is 5.8 to 6.2-fold the human
value, and Argoul 2025 independently found 6.7-fold by a different method in
a different laboratory. For mouse against human, holding fu equal moves the
difference from glomerular to tubular, as the previous version of this
analysis showed.

SO THE CORRECTED POSITION IS NARROWER AND FIRMER THAN EITHER EARLIER ONE.
The species contrast in renal reabsorption between human and rat IS
established, and it is larger in character than the PHG states: the human
reabsorbs essentially all filtered PFOA while the male rat reabsorbs a little
over half. The percentages OEHHA prints are not reliable. The attribution of
the difference to the tubule is.

TWO LABORATORIES AGREE, with a consistent offset. For PFOA, Ryu gives human
0.00245 and mouse 0.0142-0.0152; Argoul 2025 gives human 0.0013-0.0017 and
mouse 0.0087-0.012. Ryu runs 1.3 to 1.9-fold higher throughout, which is
within the spread of equilibrium-dialysis methods, and the species ORDERING
and RATIOS agree closely.

Idempotent.
"""
import csv, io, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "plasma_unbound_fraction.csv")
SRC = ("Ryu 2024, Toxics 12:253, doi:10.3390/toxics12040253, Table 1")
METHOD = ("presaturation equilibrium dialysis, 5 uM PFAS in species-specific "
          "plasma with 1% DMSO")
FIELDS = ["chemical", "mw_g_mol", "species", "experiment", "fu_percent",
          "fu", "range_low_percent", "range_high_percent", "n_replicates",
          "method", "source"]

#  chemical, (human, cv), (rat, cv), (mouse CD-1, cv), (mouse C57BL/6, cv)
T1 = [
    ("PFBA", 0.229, 0.37, 0.757, 0.767),
    ("PFPeA", 0.0888, 0.18, 0.609, 0.577),
    ("PFHxA", 0.038, 0.0358, 0.169, 0.182),
    ("PFHpA", 0.00445, 0.00454, 0.0331, 0.0323),
    ("PFOA", 0.00245, 0.00298, 0.0152, 0.0142),
    ("PFNA", 0.00187, 0.0021, 0.00455, 0.00438),
    ("PFDA", 0.000618, 0.001, 0.0015, 0.00127),
    ("PFUnDA", 0.000419, 0.000349, 0.000524, 0.00127),
    ("PFDoDA", 0.000105, 0.0000796, 0.000222, 0.000205),
    ("PFBS", 0.0365, 0.0542, 0.261, 0.28),
    ("PFHxS", 0.000695, 0.00076, 0.0128, 0.013),
    ("PFOS", 0.000753, 0.00109, 0.0128, 0.00256),
    ("PFOSA", 0.0019, 0.00197, 0.00326, 0.0023),
    ("6:2 FTS", 0.0197, 0.0236, 0.0255, 0.0253),
]
SPECIES = ("human", "rat", "mouse (CD-1)", "mouse (C57BL/6)")

GFR = {"human": 2570, "rat": 14400, "mouse": 16700}
CLR = {"human": 0.06, "rat male": 18.2, "mouse male": 10}


def main():
    existing = ""
    if os.path.exists(OUT):
        existing = open(OUT, encoding="utf-8").read()
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    rows = []
    for rec in T1:
        chem = rec[0]
        for sp, fu in zip(SPECIES, rec[1:]):
            rows.append({
                "chemical": chem, "mw_g_mol": "", "species": sp,
                "experiment": "Ryu 2024", "fu_percent": round(fu * 100, 5),
                "fu": fu, "range_low_percent": "", "range_high_percent": "",
                "n_replicates": "", "method": METHOD, "source": SRC})
    for r in rows:
        w.writerow(r)
    block = s.getvalue()
    if block in existing:
        print("no change needed (Ryu 2024 already appended)")
    else:
        with open(OUT, "a", encoding="utf-8") as fh:
            fh.write(block)
        print(f"appended {len(rows)} Ryu 2024 rows to "
              f"db/plasma_unbound_fraction.csv")

    fu = {c[0]: dict(zip(SPECIES, c[1:])) for c in T1}["PFOA"]
    print()
    print("=" * 72)
    print("OEHHA PHG TABLE A6.4 RECOMPUTED FOR PFOA, measured fu")
    print("=" * 72)
    print(f"  {'species':14s} {'fu':>9s} {'filtered':>9s} {'CL_R':>7s} "
          f"{'reabsorbed':>11s}   OEHHA printed")
    cases = (("human", fu["human"], GFR["human"], CLR["human"], 99.88),
             ("rat male", fu["rat"], GFR["rat"], CLR["rat male"], 93.68),
             ("mouse male", fu["mouse (C57BL/6)"], GFR["mouse"],
              CLR["mouse male"], 97.01))
    out = {}
    for name, f, g, c, printed in cases:
        filt = f * g
        pct = 100 * (1 - c / filt)
        out[name] = (filt, pct)
        print(f"  {name:14s} {f:9.5f} {filt:9.2f} {c:7.2f} {pct:10.1f}%   "
              f"{printed:.2f}%")
    print()
    print("  The male rat does not reabsorb 94% of filtered PFOA. It")
    print(f"  reabsorbs about {out['rat male'][1]:.0f}%.")
    print()
    print("=" * 72)
    print("BUT THE RAT/HUMAN STRUCTURE SURVIVES, AND THE MOUSE DOES NOT")
    print("=" * 72)
    print(f"  rat/human plasma fu for PFOA   {fu['rat'] / fu['human']:.2f}x  "
          f"-- close to the 1.0 OEHHA assumes")
    print(f"  mouse/human                    "
          f"{fu['mouse (C57BL/6)'] / fu['human']:.2f}x  -- not close at all")
    print()
    clr_ratio = CLR["rat male"] / CLR["human"]
    g_assumed = GFR["rat"] / GFR["human"]
    g_meas = out["rat male"][0] / out["human"][0]
    print(f"  the observed rat/human renal clearance ratio is {clr_ratio:.0f}x")
    print(f"     with fu assumed equal   {g_assumed:.1f}x glomerular x "
          f"{clr_ratio / g_assumed:.0f}x tubular")
    print(f"     with fu measured        {g_meas:.1f}x glomerular x "
          f"{clr_ratio / g_meas:.0f}x tubular")
    print()
    print("  So OEHHA's assumption is wrong in its VALUE and approximately")
    print("  right in its STRUCTURE for rat against human. The percentages it")
    print("  prints are not reliable; the attribution of the difference to")
    print("  the tubule is.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
