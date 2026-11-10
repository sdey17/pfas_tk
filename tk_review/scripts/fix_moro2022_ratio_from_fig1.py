#!/usr/bin/env python3
"""Moro 2022: the ligand:protein ratio read off Figure 1, per compound.

An ITC isotherm is PLOTTED against molar ratio, so the ratio this
assay reached is printed on the x-axis of Figure 1b and 1c.  That is
better than any derivation: the earlier value rested on the PEAQ-ITC's
200 uL cell volume, which the paper never states, and it assumed one
protocol for all four compounds.

MEASURED, from tx2c00211_0002.jpg in the open-access PMC copy
(PMC9682524), by locating each panel's plot frame and the rightmost
data marker and mapping it onto the printed axis:

    compound    298 K   310 K      derived from the stated quantities
    HFPO-DA       6.9     6.9      6.9 simple dilution
    C6O4          7.6     6.7      7.6 with overflow correction
    PFHxA         3.8     3.8      <-- HALF. Not the same titration.
    PFOA          7.5     7.2

Three of the four land between 6.7 and 7.6, which confirms the
derivation and kills the 36.4 that was reaching the workbook: 4000 uM
is the syringe and 110 uM the cell, and dividing them ignores that
only 37.9 uL goes into a 200 uL cell.

PFHxA IS THE FINDING.  It reaches 3.8 at BOTH temperatures, about half
the others, and its axis is drawn 0-4 where the others are 0-7 or 0-8,
so this is the data ending and not a cropped plot.  Section 2.3 says
every ligand was 4 mM and dhSA 110 uM with the same injection
schedule, which cannot produce half the molar ratio for one compound.
Something about the PFHxA titration differed -- a 2 mM syringe, or
double the protein -- and the paper does not say which.  Recorded as
measured, with the discrepancy noted, rather than normalised to match
the Methods.

This does not change any constant.  It changes what each constant can
be compared with, and it is the per-compound number a loading plot
should show instead of one value stamped on all 12 rows.

Run:  python3 tk_review/scripts/fix_moro2022_ratio_from_fig1.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

# (compound, T_K) -> molar ratio at the last plotted point, Figure 1b/1c
MEASURED = {
    ("HFPO-DA (GenX)", 298): "6.9", ("HFPO-DA (GenX)", 310): "6.9",
    ("C6O4", 298): "7.6",           ("C6O4", 310): "6.7",
    ("PFHxA", 298): "3.8",          ("PFHxA", 310): "3.8",
    ("PFOA", 298): "7.5",           ("PFOA", 310): "7.2",
}

BASIS = ("measured off the x-axis of Figure 1{panel} (an ITC isotherm is "
         "plotted against molar ratio, so the paper states this number "
         "itself): panel frame located and the rightmost data marker "
         "mapped onto the printed 0-{amax} axis, giving {v} ligands per "
         "albumin at the end of the titration. This REPLACES a derived "
         "7.6 that assumed the PEAQ-ITC's 200 uL cell, a volume the "
         "paper does not state, and assumed one protocol for all four "
         "compounds")
PFHXA = ("PFHxA reaches only {v}, about half the other three, at BOTH "
         "temperatures, and its axis is drawn 0-4 where the others are "
         "0-7 or 0-8, so this is where the data ends rather than a "
         "cropped plot. Section 2.3 states 4 mM ligand and 110 uM dhSA "
         "with one injection schedule for every compound, which cannot "
         "give half the molar ratio for one of them; the paper does not "
         "say what differed")
AXIS = {("HFPO-DA (GenX)", 298): 7, ("HFPO-DA (GenX)", 310): 7,
        ("C6O4", 298): 8, ("C6O4", 310): 7,
        ("PFHxA", 298): 4, ("PFHxA", 310): 4,
        ("PFOA", 298): 8, ("PFOA", 310): 8}
OLD = ("LOADING: ligand_over_protein_max is the molar ratio reached in "
       "the ITC CELL")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"] == "Moro et al." and r["year"] == "2022"]
    print(f"  {len(hit)} Moro 2022 rows\n")
    print(f"  {'PFAS':16}{'T':>5}{'was':>7}{'now':>7}   basis")
    n = 0
    for r in sorted(hit, key=lambda r: (r["chemical"],
                                        r["temperature_c"])):
        T = int(r["temperature_c"].split()[0]) + 273
        key = (r["chemical"], T)
        if key not in MEASURED:
            print(f"  {r['chemical']:16}{T:>5}   no measurement, SKIPPED")
            continue
        v = MEASURED[key]
        was = r["ligand_over_protein_max"]
        if was != v:
            n += 1
        r["ligand_over_protein_max"] = v

        panel = "b" if T == 298 else "c"
        note = BASIS.format(panel=panel, amax=AXIS[key], v=v)
        if r["chemical"] == "PFHxA":
            note += ". " + PFHXA.format(v=v)
        # Replace the superseded derivation rather than stacking on it.
        parts = [p for p in r["assay_volume_note"].split(" | ")
                 if not p.startswith(OLD) and not p.startswith("LOADING:")]
        r["assay_volume_note"] = " | ".join(parts + ["LOADING: " + note])
        print(f"  {r['chemical']:16}{T:>5}{was:>7}{v:>7}   Fig 1{panel}, "
              f"axis 0-{AXIS[key]}")

    vals = sorted(float(r["ligand_over_protein_max"]) for r in hit)
    print(f"\n  {n} rows changed; ratio now spans {vals[0]:g} to {vals[-1]:g} "
          f"instead of a single 7.6 on every row")

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
