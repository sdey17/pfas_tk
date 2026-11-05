#!/usr/bin/env python3
"""Peng 2024's two rows per compound are two TEMPERATURES, not two studies.

Peng Y, Xu Y, Gao Y, et al.  Toxics 2024;12(1):43.
doi 10.3390/toxics12010043, PMID 38250999.  Open access.

WHAT WAS WRONG HERE.  This database held two albumin rows per compound
for Peng and treated them as two independent candidates -- they were
inflating the admissible spread in fu_physiological.py and, for PFNA,
producing a 509-fold range on their own.  They are the same binding
constant measured at two temperatures.

Peng's Table 2 reports Kb and n at 298, 304 and 310 K.  Table 3's
"no probe" row is the SAME 298 K number re-listed as the baseline for
the site-competition experiment -- identical to four significant
figures for PFNA, PFOA, PFO3DA, PFHpA and DFSA.  It is not a second
measurement.

Three faults, all in the ingestion, none in the paper:

  1  the temperature column carried all three values as one lump,
     "25, 31, 37 (298, 304, 310 K)", on every row, so neither row
     said which temperature it was.
  2  the 298 K row lost its n.  Table 2 gives n at every temperature;
     without it the conversion to log K assumed n = 1, which is wrong
     by log(1.51) = 0.18 for PFNA and by log(0.96) = -0.02 for PFHpA.
  3  nothing recorded that the two rows are one experiment, so a
     selection rule could pick between them as though they were
     independent evidence.

WHICH ONE A PHYSIOLOGICAL CALCULATION WANTS: 310 K.  That is 37 C.
298 K is 25 C and is not a body temperature.

A CAUTION THAT SURVIVES THE FIX.  The temperature dependence Peng
reports is extreme and not monotonic in one direction across
compounds: PFNA falls 79-fold from 298 to 310 K while HFPO-TA rises
624-fold and PFHpA rises 165-fold.  The van't Hoff enthalpies the
paper derives from this are +412 kJ/mol for HFPO-TA and +325 for
PFHpA.  Non-covalent protein-ligand binding is normally within about
+/- 50 kJ/mol; 412 kJ/mol is covalent-bond territory.  The numbers are
transcribed as reported and the temperatures are now right, but the
series should not be read as a reliable thermodynamic result.

Run:  python3 tk_review/scripts/fix_peng2024_temperatures.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR, DOI = "Peng et al.", "2024", "10.3390/toxics12010043"

# Table 2, transcribed as printed: Ka (L/mol) -> (temperature K, n)
# at 298 and 310 K. The 304 K series is in the paper and is not held
# in this database.
BY_KA = {
    "7.81e6": ("298", "1.51"), "9.9e4": ("310", "1.10"),
    "3.7e6":  ("298", "1.49"), "2.31e9": ("310", "2.09"),
    "2.27e5": ("298", "1.26"), "1.98e6": ("310", "1.47"),
    "1.59e5": ("298", "1.25"), "2.99e4": ("310", "1.10"),
    "4.53e3": ("298", "0.96"), "7.47e5": ("310", "1.50"),
    "1.52e3": ("298", "0.88"), "2.07e4": ("310", "1.13"),
}
TEMP_C = {"298": "25", "310": "37"}
NOTE = ("Table 2 (Kb and n at 298, 304 and 310 K). The 298 K value is "
        "also Table 3's no-probe baseline for the site-competition "
        "experiment, re-listed rather than measured again, so the two "
        "rows this database holds for each compound are ONE experiment "
        "at TWO temperatures and are not independent. A physiological "
        "calculation wants the 310 K row (37 C). Peng's van't Hoff "
        "enthalpies from this series reach +412 kJ/mol, far outside "
        "the normal range for non-covalent binding; treat the "
        "temperature dependence with caution.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, list(reader)

    fixed = 0
    for r in rows:
        if not (r["study"] == STUDY and r["year"] == YEAR
                and r["doi"] == DOI
                and "albumin" in r["protein"].lower()):
            continue
        hit = BY_KA.get(r["ka"].strip())
        if hit is None:
            continue
        K, n = hit
        r["temperature_c"] = f"{TEMP_C[K]} ({K} K)"
        if not r["n_binding_sites"].strip():
            r["n_binding_sites"] = n
        base = r["source_table"].split(" | ")[0]
        if "ONE experiment" not in r["source_table"]:
            r["source_table"] = f"{base} | {K} K | {NOTE}"
        fixed += 1
        print(f"  {r['chemical'][:22]:24}Ka {r['ka']:10}"
              f"-> {TEMP_C[K]} C ({K} K), n = {r['n_binding_sites']}")

    print(f"\n{fixed} rows corrected")
    if args.dry_run:
        print("--dry-run: nothing written")
        return
    with open(CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {os.path.relpath(CSV, os.path.join(HERE, '..'))}")


if __name__ == "__main__":
    main()
