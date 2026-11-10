#!/usr/bin/env python3
"""Qin 2024: the real per-compound BSA level, and what the ratio means.

Two defects, found by asking where the plotted Qin values came from.

DEFECT 1 -- THE PROTEIN CONCENTRATION WAS A PARSE, NOT A VALUE.
protein_conc held the prose "0.1-10 mg/mL BSA", and the condition
parser took the TOP of that range for every compound: 10 mg/mL =
10,000 mg/L = 150.46 uM, identical on all eleven albumin rows.

SI Table S5 gives the level actually used, and it is per compound,
tuned to keep the bound fraction measurable:

    PFHxS, PFOS        100 mg/L  =  1.50 uM
    PFHpA, PFOA, PFNA, PFUnDA,
    6:2 FTSA, PFOSA    500 mg/L  =  7.52 uM
    PFHxA, HFPO-DA    1500 mg/L  = 22.57 uM
    PFBA              5000 mg/L  = 75.23 uM

So the parsed 150.46 uM was wrong by 2x to 100x, and 10,000 mg/L is
not even a level these compounds were run at -- it is the FTOH level
from the same table.

DEFECT 2 is NOT fixed here, because it cannot be: the ratio and the
log K describe different conditions, and the condition the log K
belongs to is not in this file.

A Qin log K is not a fit over the whole isotherm.  It is SI Table
S6's Freundlich equation EVALUATED at one free concentration -- the
Cw reached from 10 uM nominal in 10% FBS -- and that Cw is applied
downstream, by albumin_conc_vs_logk.QIN_CW, not stored here.
ligand_over_protein_max meanwhile keeps its project-wide meaning:
the highest ratio the ASSAY reached, which for Qin is Table S5's
Ctot/BSA of 10 to 37.  Both are true of different things, which is
why dividing the concentration columns never reproduced the ratio
and why these rows sit on the working-table audit's ratio list.

The loading at the operating point follows from the coefficient
itself, nu = D*Cw*M/1e6, and comes out near 0.15 for all eleven --
about 100x below the isotherm maximum.  That calculation belongs
where the Cw is, so make_route_and_loading_figures.py does it and
plots Qin there instead of at the isotherm maximum.

Only the BSA rows are touched.  The plasma arms of the same study
were run at their own volume concentration -- 100, 30, 10 or 2 mL/L
of plasma -- and a BSA molarity does not describe them.

Run:  python3 tk_review/scripts/fix_qin2024_protein_conc.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
DOI = "10.1021/acs.est.3c10824"
MW_BSA = 66463.0

BSA_MGL = {"PFBA": 5000, "PFHxA": 1500, "HFPO-DA": 1500, "HFPO-DA (GenX)": 1500,
           "PFHpA": 500, "PFOA": 500, "PFNA": 500, "PFUnA": 500,
           "PFUnDA": 500, "6:2 FTSA": 500, "6:2 FTS": 500, "PFOSA": 500,
           "FOSA": 500, "hexaflumuron": 500, "flubendiamide": 500,
           "PFHxS": 100, "PFOS": 100,
           "6:2 FTOH": 10000, "8:2 FTOH": 10000, "10:2 FTOH": 10000}

NOTE = ("PROTEIN CONCENTRATION is SI Table S5's per-compound BSA level, "
        "{mg:g} mg/L = {um:.4g} uM. It replaces a parse of the prose "
        "'0.1-10 mg/mL BSA', which took the top of that range and put "
        "150.46 uM (10,000 mg/L) on every albumin row -- wrong by 2x to "
        "100x here, and a level only the FTOHs were run at")
MARKS = ("PROTEIN CONCENTRATION is SI Table S5",)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)
    hit = [r for r in rows if r["doi"] == DOI]
    alb = [r for r in hit if "albumin" in r["protein"].lower()]
    print(f"  {len(hit)} plasma-paper rows, {len(alb)} on albumin\n")
    print(f"  {'PFAS':14}{'BSA mg/L':>10}{'was (parsed)':>12}{'now uM':>11}")
    n = 0
    for r in rows:
        if r["doi"] != DOI:
            continue
        # BSA rows ONLY. Table S5's BSA column describes the protein
        # binding assay; the plasma arms were run at their own volume
        # concentration (100/30/10/2 mL/L) and writing a BSA molarity
        # onto them would be wrong.
        if "albumin" not in r["protein"].lower():
            continue
        mg = BSA_MGL.get(r["chemical"]) or BSA_MGL.get(
            r["chemical"].split(" (")[0])
        if not mg:
            continue
        um = mg / MW_BSA * 1000.0
        old = r["protein_conc"]
        r["protein_conc"] = f"{um:.4g}"
        r["protein_conc_units"] = "uM"
        note = NOTE.format(mg=mg, um=um)
        if old != r["protein_conc"]:
            print(f"  {r['chemical'][:13]:14}{mg:>10g}"
                  f"{('150.46' if '0.1-10' in old else old)[:9]:>12}"
                  f"{um:>11.4g}")
        parts = [x for x in r["assay_volume_note"].split(" | ")
                 if not any(m in x for m in MARKS)]
        r["assay_volume_note"] = " | ".join(parts + [note])
        n += 1
    print(f"\n  {n} rows updated")

    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader(); w.writerows(rows)
    print(f"  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
