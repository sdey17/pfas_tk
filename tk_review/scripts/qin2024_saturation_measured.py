#!/usr/bin/env python3
"""The non-linearity of PFAS protein binding, measured across 11 compounds.

Section 7 of MBM_REVIEW argues that the assays run far above the binding
regime a linear model assumes.  Until now that argument ran through Qin et
al. 2023's saturable model fitted to four compounds.  Qin et al. 2024 Table
S6 reports the concentration dependence DIRECTLY, as Freundlich-type fits

    log D = slope x log Cw + intercept          Cw in umol/L of water

for eleven PFAS in two matrices.  The slope is the whole argument: a linear
sorption model asserts slope = 0.

Run:  python3 tk_review/scripts/qin2024_saturation_measured.py
"""
import csv
import math
import os
import re
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
ORDER = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "HFPO-DA",
         "PFHxS", "PFOS", "6:2 FTSA", "PFOSA"]
# Qin et al. 2023 Table 1, the four compounds the review has been using
QIN23 = {"PFBA": -0.391, "PFOA": -0.381, "PFHxS": -0.480, "PFOS": -0.381}


# The database now holds five Qin 2024 matrices -- BSA, the PPARgamma
# medium, human plasma, trout plasma, and the cell and structural-protein
# ratios -- across two papers. This script is about the first two only, so
# the matrix is selected explicitly rather than by "albumin or else".
# Matching on the protein string alone would sweep the plasma rows into
# "medium" and silently overwrite it.
MATRIX_OF = {
    "bovine serum albumin (BSA)": "BSA",
    "PPARgamma assay medium, 98% OptiMEM + 2% csFBS (protein + lipid)":
        "medium",
}
# Accepts both the equation wording this file used to write ("log D = ...")
# and the matrix-tagged wording ("log D_BSA/w = ...").
EQ = re.compile(r"log D(?:_[^\s]+)? = ([+-]?[\d.]+) x log Cw \+ ([\d.]+)")


def load():
    out = {}
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if r["study"] != "Qin et al." or r["year"] != "2024":
            continue
        mat = MATRIX_OF.get(r["protein"])
        if mat is None or "Freundlich-type fit" not in r["ka_units"]:
            continue
        m = EQ.search(r["source_table"])
        slope = float(m.group(1)) if m else 0.0
        out.setdefault(r["chemical"], {})[mat] = (slope, float(r["ka"]))
    missing = [nm for nm in ORDER
               if nm not in out or set(out[nm]) != {"BSA", "medium"}]
    if missing:
        raise SystemExit(f"missing BSA or medium fits for: {missing}")
    return out


def main():
    print(__doc__.split("Run:")[0].rstrip())
    D = load()

    print("\n" + "=" * 76)
    print("1.  THE MEASURED SLOPES.  A LINEAR MODEL ASSERTS ZERO.")
    print("=" * 76)
    print(f"\n  {'PFAS':10s}{'BSA slope':>11s}{'medium slope':>14s}"
          f"{'D_BSA/w falls per':>19s}")
    print(f"  {'':10s}{'':11s}{'':14s}{'decade of C_w':>19s}")
    sb = []
    for nm in ORDER:
        b, m = D[nm]["BSA"][0], D[nm]["medium"][0]
        sb.append(b)
        fall = f"{10 ** -b:.1f}x" if b else "flat"
        print(f"  {nm:10s}{b:>+11.3f}{m:>+14.3f}{fall:>19s}")
    nz = [x for x in sb if x]
    print(f"""
  Every compound with a fitted slope has a NEGATIVE one: {len(nz)} of
  {len(sb)}, from {min(nz):+.3f} to {max(nz):+.3f}, mean {st.mean(nz):+.3f}.
  The apparent albumin affinity therefore falls by {10 ** -st.mean(nz):.1f}x for every
  decade the free concentration rises.  A linear spLFER treats it as a
  constant.""")

    print("\n" + "=" * 76)
    print("2.  AGAINST THE FOUR COMPOUNDS THE REVIEW HAS BEEN USING")
    print("=" * 76)
    print(f"\n  {'PFAS':10s}{'Qin 2023':>11s}{'Qin 2024':>11s}{'diff':>9s}")
    for nm, s23 in QIN23.items():
        s24 = D[nm]["BSA"][0]
        print(f"  {nm:10s}{s23:>+11.3f}{s24:>+11.3f}{s24 - s23:>+9.3f}")
    print("""
  The two sulfonates reproduce almost exactly; the two carboxylates come out
  shallower in 2024.  Either way the sign and the order of magnitude hold, so
  the four-compound argument section 7 already makes is not overturned -- it
  is extended to eleven and to three more structural classes.""")

    print("\n" + "=" * 76)
    print("3.  WHAT IT COSTS ACROSS A REAL DOSE SERIES")
    print("=" * 76)
    print("""
  The three uptake datasets in this review dose over roughly 1.2 decades
  (Fragki 6.25-100 uM, Rosenmai 10-100 uM).  Free concentration spans less
  than nominal, but the measured slopes give the fall in apparent affinity
  directly.
""")
    print(f"  {'PFAS':10s}{'over 1.2 decades':>18s}{'over 3 decades':>17s}")
    for nm in ORDER:
        b = D[nm]["BSA"][0]
        if not b:
            print(f"  {nm:10s}{'flat':>18s}{'flat':>17s}")
            continue
        print(f"  {nm:10s}{10 ** (-b * 1.2):>17.1f}x{10 ** (-b * 3):>16.1f}x")
    print(f"""
  Across one dose series the apparent albumin affinity of a typical PFAS
  falls by about {10 ** (-st.mean(nz) * 1.2):.1f}x.  Across the three decades separating an
  assay from human serum it falls by {10 ** (-st.mean(nz) * 3):.1f}x on the mean slope, and by
  {10 ** (-min(nz) * 3):.0f}x for the steepest compound (HFPO-DA).

  That is the quantity a linear model has to hold constant, and it is now a
  measurement rather than an inference from a fitted saturable model.""")

    print("\n" + "=" * 76)
    print("4.  HEAD GROUP")
    print("=" * 76)
    grp = {"carboxylate": ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"],
           "sulfonate": ["PFHxS", "PFOS"],
           "other": ["HFPO-DA", "6:2 FTSA", "PFOSA"]}
    for g, ms in grp.items():
        v = [D[m]["BSA"][0] for m in ms if D[m]["BSA"][0]]
        print(f"  {g:14s} n={len(v)}  mean slope {st.mean(v):+.3f}"
              f"   range {min(v):+.3f} to {max(v):+.3f}")
    print("""
  PFUnA is the one compound whose binding is concentration-INDEPENDENT over
  the range tested (R2 0.96 on a constant).  It is also the longest chain and
  the one dosed at the lowest concentration, so it is the compound furthest
  into the linear low-occupancy limb -- which is the regime the whole of
  section 7 says the assays do not otherwise occupy.""")


if __name__ == "__main__":
    sys.exit(main())
