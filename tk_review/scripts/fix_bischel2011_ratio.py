#!/usr/bin/env python3
"""Bischel 2011: give each compound its own ligand:protein ratio.

All 11 rows carry ligand_over_protein_max = 0.0034.  That is the
maximum across the whole study, applied uniformly, and it overstates
10 of the 11 compounds.

Nothing new has to be read to fix it.  Bischel holds BSA at 200 uM
for every compound and reports a PER-COMPOUND reservoir range, which
this database already has: PFOS ran 0.013-0.26 uM, PFDoDA 0.031-0.68.
The ratio is the top of that range over 200 uM, and it spans 0.00130
to 0.00340 -- a 2.6-fold difference between compounds.

This changes no conclusion.  Every value stays three orders of
magnitude below unit loading, which is the point of the study and
why Bischel is the most physiologically loaded source in the file.
It is corrected because a column that reports the same number for
eleven different experiments is not reporting a measurement, and
because the per-compound value is the one a loading plot should
show.

Run:  python3 tk_review/scripts/fix_bischel2011_ratio.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")


def top_uM(text):
    """Highest reservoir concentration from 'a to b'."""
    m = re.match(r"\s*([\d.eE+-]+)\s*(?:to|-|–)\s*([\d.eE+-]+)", text)
    return max(float(m.group(1)), float(m.group(2))) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"].startswith("Bischel")
           and r["year"] == "2011"]
    print(f"  {len(hit)} Bischel 2011 rows\n")
    print(f"  {'PFAS':10}{'ligand uM':>16}{'[BSA] uM':>10}"
          f"{'was':>9}{'now':>10}")
    n = 0
    for r in sorted(hit, key=lambda r: r["chemical"]):
        if r["ligand_conc_units"].strip() != "uM" \
                or r["protein_conc_units"].strip() != "uM":
            print(f"  {r['chemical']:10} SKIPPED, units not uM")
            continue
        hi = top_uM(r["ligand_conc_range"])
        P = float(r["protein_conc"])
        if hi is None or not P:
            print(f"  {r['chemical']:10} SKIPPED, range did not parse")
            continue
        v = f"{hi / P:.5f}"
        was = r["ligand_over_protein_max"]
        if was != v:
            n += 1
        r["ligand_over_protein_max"] = v
        print(f"  {r['chemical']:10}{r['ligand_conc_range']:>16}{P:>10.0f}"
              f"{was:>9}{v:>10}")

    vals = sorted(float(r["ligand_over_protein_max"]) for r in hit)
    print(f"\n  {n} rows changed; ratio now spans {vals[0]:.5f} to "
          f"{vals[-1]:.5f} ({vals[-1]/vals[0]:.1f}-fold)")
    print(f"  all still ~{1/vals[-1]:.0f}x below unit loading")

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
