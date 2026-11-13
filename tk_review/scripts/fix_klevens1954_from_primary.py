#!/usr/bin/env python3
"""Correct the Klevens & Ellenbogen 1954 row against the primary.

The row was ingested from a secondary (Starnes & Belcher 2026 Table 1)
and its source_table still said "primary not retrieved".  The primary
is Klevens HB, Ellenbogen E, "Protein-fluoroacid interaction: bovine
serum albumin-perfluoro-octanoic acid", Discuss Faraday Soc
1954;18:277-288, doi 10.1039/df9541800277.  Reading it settles four
things and leaves one open.

CONFIRMED.  Fig. 3's caption: "Value of n, number of binding sites
(extrapolated) = 85 and dissociation constant K = 3.1 x 10-3."  So
Kd = 3.1 mM, which is what the row carries (3106 uM, from 1/3.22e2).

WRONG, AND CORRECTED HERE.

  pH.  The row said "pH 3.54-4.3 -- THREE pH UNITS BELOW BLOOD and
  below BSA's isoelectric point".  The paper says "only one buffer,
  acetate, has been used in the region between pH 4.3-5.8", and
  Table 1 runs 4.30, 4.60, 5.03, 5.35, 5.55, 5.80.  There is no 3.54
  anywhere in it.  The number comes from a text-extraction artefact:
  "pH 4.3-5.8" loses its middle dots and reads as "pH 4-3-54".
  The correction matters in both directions -- the range is 1.6 to
  3.1 pH units below blood rather than three, and most of the data,
  Figs. 3 to 5 included, sit ABOVE BSA's isoelectric point near 4.7
  rather than below it.  The paper itself splits its curves on that
  point: "All points fall on one curve for values of pH above the
  isoelectric point; on another curve for values of pH below."

  Temperature.  The row said "20 (sedimentation reference)".  20 C is
  the reference medium for the SEDIMENTATION calculation, not a
  binding temperature.  The binding section says change in
  temperature "has little effect on the binding capacity ... for the
  binding of the perfluorooctanoate ion at pH 5.8 between 4 C and
  37 C", with "a slight increase in binding at lower temperatures".
  So 37 C is covered, and the constant is not strongly temperature
  dependent -- worth more than the sedimentation reference was.

  Provenance.  source_table now records the primary.

LEFT OPEN, DELIBERATELY: the stoichiometry.  n is stored as
"85 (extrapolated); curve breaks at a mole ratio of 28-29", which
site_count() cannot parse, so the converter falls through to n = 1 and
the sheet prints "n = 1 assumed by the reader, THE PAPER REPORTS
NONE".  That last clause is false -- the paper reports n in a figure
caption -- but it is produced by a parse failure, not by a judgement,
and making the string parseable would move this row by 1.9 log.  The
choice is not obvious either, because n here is not a site count in
the modern sense:

    0.5% BSA   n about 150        4% BSA    n about 45
    1%   BSA   n about 85         2% BSA    n about 60

  n falls with protein concentration because the protein aggregates --
  "a marked loss of potential binding sites with increase in protein
  concentration ... association of BSA molecules occurs".  The 85 in
  Fig. 3 belongs to the 1% solutions, i.e. 147 uM BSA, while the row
  stores 588 uM (4%), whose own n is about 45.  For comparison Han
  2003 counts 6 to 9 sites and Beesoon 2015 counts 2.

  So this script does NOT touch n.  See the decision note it prints.

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/fix_klevens1954_from_primary.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "primary read 2026-11-13"

NEW_BUFFER = (
    "acetate, ionic strength 0.1, pH 4.3-5.8 (Table 1: 4.30, 4.60, "
    "5.03, 5.35, 5.55, 5.80) -- 1.6 to 3.1 pH units below blood. "
    "BSA's isoelectric point is near 4.7, and the paper splits its "
    "curves on it: 'All points fall on one curve for values of pH "
    "above the isoelectric point; on another curve for values of pH "
    "below the isoelectric point.' Figs. 3-5, which carry the "
    "constant, are at the upper end (Fig. 4 and 5 state pH 5.8), so "
    "this row is ABOVE the isoelectric point, not below it. An "
    "earlier version of this field said pH 3.54-4.3; that number is "
    "a text-extraction artefact of 'pH 4.3-5.8' and appears nowhere "
    "in the paper")

NEW_TEMP = (
    "4 to 37 (binding); the paper reports that temperature 'has "
    "little effect on the binding capacity' for perfluorooctanoate "
    "with BSA at pH 5.8 between 4 and 37 C, with a slight increase "
    "at lower temperatures. The 20 C recorded here previously is the "
    "reference medium for the SEDIMENTATION calculation, not a "
    "binding temperature")

NEW_SOURCE = (
    "Fig. 3 caption, 'Value of n, number of binding sites "
    "(extrapolated) = 85 and dissociation constant K = 3.1 x 10-3' | "
    "read from the primary, Discuss Faraday Soc 1954;18:277-288 "
    "(" + MARK + "). First reported PFAS-albumin constant and the "
    "weakest in the literature. The secondary that seeded this row "
    "(Starnes & Belcher 2026 Table 1) reproduces it exactly")


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())

    hit = 0
    for r in rows:
        if "Klevens" not in (r.get("study") or ""):
            continue
        if MARK in (r.get("source_table") or ""):
            print("  already corrected; nothing to do")
            return 0
        r["buffer"] = NEW_BUFFER
        r["temperature_c"] = NEW_TEMP
        r["source_table"] = NEW_SOURCE
        hit += 1

    if not hit:
        print("  no Klevens row found", file=sys.stderr)
        return 1

    with open(DB, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  corrected {hit} row: buffer/pH, temperature, provenance")
    print()
    print("  NOT CHANGED -- the stoichiometry, which needs a decision:")
    print("    the sheet says 'n = 1 assumed by the reader, the paper")
    print("    reports none'.  The paper DOES report n, in Fig. 3's")
    print("    caption, but the stored string is unparseable so the")
    print("    converter defaults to 1.  Options and what each gives:")
    import math
    for lbl, n in (("n = 1, as now", 1.0),
                   ("n = 85, Fig. 3 (1% BSA = 147 uM)", 85.0),
                   ("n = 45, the 4% BSA this row stores", 45.0)):
        print(f"      {lbl:38s} log K "
              f"{4.1774 + math.log10(n) - math.log10(3106.0):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
