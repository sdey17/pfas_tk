#!/usr/bin/env python3
"""Chen 2025: two wrong notes, and the loading claim reversed.

Checked against the full text (PMC11740993, doi 10.1021/acs.est.4c06734).

DEFECT 1 -- A 1000x UNIT ERROR IN THE ASSAY NOTE.

add_chen2025_tableS5.py set protein_conc to 1.0 uM, which is right,
and then wrote an assay note saying "100 uL of 1.0 mM HSA".  The
paper's Methods say "100 uL of 1.0 uM HSA was added to the dialysis
cup".  Micromolar.  The structured field was never wrong; the prose
beside it was, by three orders of magnitude, and prose is what a
reader checks the number against.

DEFECT 2 -- THE LOADING CLAIM IS BACKWARDS.

The note called this "the closest to physiological of any
wide-coverage dataset here".  It is nearly the furthest.  The paper
states the design plainly: "PFAS to HSA molar ratios of 1:16, 1:8,
1:4, 1:2, 1:1, and 2:1".  The top point is TWO ligands per albumin.

Physiological occupancy is nu < 0.01.  Solving 1:1 binding with
depletion at the top point ([HSA] 1 uM, PFAS 2 uM) gives nu from
0.023 for the weakest binder to 0.84 for PFMBA -- all 60 compounds
above the physiological band, most of them by one to two orders.
Ranked against the other sources in this file by median
ligand:protein ratio, Chen 2025 sits at 2, alongside Gao 2019, while
Bischel 2011 (0.002), Fischer 2024 (0.004) and Allendorf 2019
(0.009) are 200-1000x closer.

The note was comparing against the DSF and fluorescence datasets,
where loading is higher still, and that narrower comparison is true.
Stated without it, the claim inverts the ranking.

WHAT IS NOT A DEFECT, and is recorded rather than corrected: 49 of
the 60 K_D values lie ABOVE the highest concentration tested (2 uM),
PFBA's 84.9 uM being 42x above it.  That is not a transcription
problem and not a criticism of the paper, which says so itself --
"we extrapolated the binding curves, assuming single-site binding".
Below the K_D a Langmuir fit constrains the initial slope, which is
B_max/K_D; it does not separate the two.  For a partition
coefficient, which is what this paper is deriving, the initial slope
is the right quantity and the regime is the correct one to work in.
The K_D is a re-parameterisation of it that inherits whatever n was
assumed, so it should not be read as a site affinity, and it is that
reading the new note blocks.

Run:  python3 tk_review/scripts/fix_chen2025_notes.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

NOTE = ("1000 uL PFAS solution in the vial, 100 uL of 1.0 uM HSA in "
        "the Slide-A-Lyzer vessel (Methods: '100 uL of 1.0 uM HSA'). "
        "PFAS:HSA molar ratios 1:16, 1:8, 1:4, 1:2, 1:1, 2:1, so the "
        "top point is TWO ligands per albumin")

RELEVANCE = (
    "K_D measured across PFAS:HSA molar ratios of 1:16 to 2:1 (HSA "
    "fixed at 1.0 uM, PFAS 0.0625-2 uM). NOT a physiological "
    "loading: solving 1:1 binding with depletion at the top point "
    "gives occupancy nu 0.023-0.84 against a physiological nu < "
    "0.01, so all 60 compounds sit above the physiological band. By "
    "median ligand:protein ratio this source sits at 2, with Gao "
    "2019; Bischel 2011, Fischer 2024 and Allendorf 2019 are "
    "200-1000x closer to physiological. K_D IS AN EXTRAPOLATION for "
    "49 of 60 compounds, whose K_D exceeds the 2 uM top spike (PFBA "
    "84.9 uM is 42x above it); the paper says so -- 'we extrapolated "
    "the binding curves, assuming single-site binding'. Below K_D a "
    "Langmuir fit fixes the initial slope B_max/K_D, not K_D and "
    "B_max separately, so these are partition coefficients "
    "re-expressed as affinities and should not be read as site "
    "affinities. Paired log K_MW lets the albumin and phospholipid "
    "sinks be compared directly: protein-water distribution exceeds "
    "membrane-water partitioning at low PFAS concentration and the "
    "balance shifts to the membrane as concentration rises")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"].startswith("Chen")
           and r["year"] == "2025"
           and "albumin" in r["protein"].lower()]
    n1 = n2 = 0
    for r in hit:
        if r["assay_volume_note"] != NOTE:
            if "1.0 mM HSA" in r["assay_volume_note"]:
                n1 += 1
            r["assay_volume_note"] = NOTE
        if r["tk_relevance_noted"] != RELEVANCE:
            n2 += 1
            r["tk_relevance_noted"] = RELEVANCE

    print(f"  {len(hit)} Chen 2025 HSA rows")
    print(f"  defect 1  assay note said 1.0 mM, paper says 1.0 uM: "
          f"{n1} rows")
    print(f"  defect 2  loading claim rewritten:                   "
          f"{n2} rows")
    bad = [r for r in hit if r["protein_conc"] != "1.0"
           or r["protein_conc_units"] != "uM"]
    print(f"  structured protein_conc wrong on {len(bad)} rows "
          f"(expected 0 -- the field was always right)")

    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
