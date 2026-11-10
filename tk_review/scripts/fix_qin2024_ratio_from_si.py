#!/usr/bin/env python3
"""Qin 2024: ligand:protein ratio and bound fraction, from SI Table S5.

Qin W-P, et al.  "Species Difference? Bovine, Trout, and Human Plasma
Protein Binding of Per- and Polyfluoroalkyl Substances."  Environ Sci
Technol 2024;58:9954.  doi 10.1021/acs.est.3c10824, PMID 38804966.
The SI was supplied directly.

WHAT WAS WRONG.  136 rows, 112 with a protein concentration and 111
with a ligand range, and ZERO ratios -- because protein_conc held
prose ("0.1-10 mg/mL BSA, chosen per compound", "plasma diluted to
the matching protein level") rather than a number, so nothing could
divide.  unbound_fraction was empty on every row.

SI TABLE S5 SUPPLIES BOTH, per compound, for all 16: the total PFAS
concentration in the assay, the BSA mass concentration actually used
for that compound, and the bound fraction in BSA, human plasma and
fish plasma.  The BSA level is NOT constant -- it was tuned per
compound to keep the bound fraction measurable, spanning 100 mg/L to
10,000 mg/L, a hundredfold:

    PFBA      5000 mg/L     PFHxS      100 mg/L
    PFHxA     1500 mg/L     PFOS       100 mg/L
    HFPO-DA   1500 mg/L     6:2 FTOH   10000 mg/L
    most C7+   500 mg/L     8:2 FTOH   10000 mg/L

That tuning is the reason a single ratio could never have been
inferred from the paper's text, and the reason the ratios that come
out span four orders of magnitude (0.0003 to 36) within ONE study.

THE RATIO IS TABLE S5's OWN PAIRING.  Ctot and the BSA level sit in
the same row of the same table and describe the same assay, so
Ctot / [BSA in uM] is internally consistent whatever relation Ctot
bears to the dilution series recorded in ligand_conc_range.  BSA is
converted at 66,463 g/mol.

BOUND FRACTIONS are written into unbound_fraction as the BSA value
for albumin rows and the human-plasma value for plasma rows, each
labelled; the fish-plasma column is kept in the note because this
database has no trout rows to put it on.  Table S5 states these are
"average values of all concentration points in the binding
isotherms", not single-point measurements, and that is recorded.

Run:  python3 tk_review/scripts/fix_qin2024_ratio_from_si.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MW_BSA = 66463.0

# Table S5: name -> (Ctot uM, BSA mg/L, BSA fbound %, plasma mL/L,
#                    human plasma fbound %, fish plasma fbound %)
S5 = {
 "PFBA": (2712, 5000, 46, 100, 35, 2),
 "PFHxA": (731, 1500, 59, 30, 58, 17),
 "PFHpA": (78, 500, 65, 10, 65, 34),
 "PFOA": (100, 500, 87, 10, 77, 72),
 "PFNA": (100, 500, 86, 10, 79, 86),
 "PFUnA": (76, 500, 87, 10, 84, 89),
 "PFUnDA": (76, 500, 87, 10, 84, 89),
 "HFPO-DA": (203, 1500, 58, 30, 42, 52),
 "HFPO-DA (GenX)": (203, 1500, 58, 30, 42, 52),
 "PFHxS": (54, 100, 57, 2, 53, 31),
 "PFOS": (55, 100, 72, 2, 64, 71),
 "6:2 FTSA": (106, 500, 65, 10, 62, 67),
 "6:2 FTS": (106, 500, 65, 10, 62, 67),
 "PFOSA": (96, 500, 78, 10, 63, 64),
 "FOSA": (96, 500, 78, 10, 63, 64),
 "hexaflumuron": (29, 500, 69, 10, 77, 86),
 "flubendiamide": (37, 500, 44, 10, 51, 65),
 "6:2 FTOH": (4.5, 10000, 42, 100, 35, 22),
 "8:2 FTOH": (0.45, 10000, 62, 100, 49, 58),
 "10:2 FTOH": (0.045, 10000, 79, 100, 51, 70),
}

NOTE = ("SI Table S5: Ctot {ct:g} uM against BSA {mg:g} mg/L = {um:.3g} uM "
        "(66,463 g/mol), so ligand:protein = {r:.4g}. The BSA level is NOT "
        "constant across this study -- it was tuned per compound to keep "
        "the bound fraction measurable, from 100 to 10,000 mg/L -- which is "
        "why no single ratio could be inferred from the text and why the "
        "ratios span four orders of magnitude within one study. Bound "
        "fractions from the same table: BSA {fb}%, human plasma {hp}%, fish "
        "plasma {fp}%; Table S5 states these are 'average values of all "
        "concentration points in the binding isotherms', not single points. "
        "Plasma was used at {pl:g} mL/L")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"].startswith("Qin") and r["year"] == "2024"]
    done = miss = 0
    seen = {}
    for r in hit:
        d = S5.get(r["chemical"]) or S5.get(r["chemical"].split(" (")[0])
        if not d:
            miss += 1
            continue
        ct, mg, fb, pl, hp, fp = d
        um = mg / MW_BSA * 1000.0
        ratio = ct / um
        r["ligand_over_protein_max"] = f"{ratio:.4g}"
        is_plasma = "plasma" in r["protein"].lower()
        r["unbound_fraction"] = (
            f"{100-hp}% unbound ({hp}% bound, human plasma, SI Table S5)"
            if is_plasma else
            f"{100-fb}% unbound ({fb}% bound, BSA, SI Table S5)")
        note = NOTE.format(ct=ct, mg=mg, um=um, r=ratio, fb=fb, hp=hp,
                           fp=fp, pl=pl)
        if "SI Table S5: Ctot" not in r["assay_volume_note"]:
            r["assay_volume_note"] = (r["assay_volume_note"] + " | " + note
                                      if r["assay_volume_note"].strip()
                                      else note)
        done += 1
        seen.setdefault(r["chemical"], (ct, mg, um, ratio, fb, hp, fp))

    print(f"  {len(hit)} Qin 2024 rows: {done} filled, {miss} with no "
          f"Table S5 entry\n")
    print(f"  {'PFAS':16}{'Ctot uM':>10}{'BSA mg/L':>10}{'BSA uM':>9}"
          f"{'ratio':>9}{'fb BSA':>8}{'fb HP':>7}{'fb FP':>7}")
    for c, (ct, mg, um, rt, fb, hp, fp) in sorted(seen.items(),
                                                  key=lambda kv: -kv[1][3]):
        print(f"  {c[:15]:16}{ct:>10g}{mg:>10g}{um:>9.3g}{rt:>9.4g}"
              f"{fb:>7}%{hp:>6}%{fp:>6}%")
    vals = sorted(float(r["ligand_over_protein_max"]) for r in hit
                  if r["ligand_over_protein_max"])
    if vals:
        print(f"\n  ratio spans {vals[0]:.4g} to {vals[-1]:.4g} "
              f"({vals[-1]/vals[0]:,.0f}-fold) within one study")
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
