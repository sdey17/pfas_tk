#!/usr/bin/env python3
"""Gao 2025: fill the ligand:protein ratio, blank on all 42 rows.

WHERE THE NUMBERS COME FROM.  Verified against the Supporting
Information (es5c08779_si_001, doi 10.1021/acs.est.5c08779).

  PROTEIN.  0.8 uM BSA, stated in the caption of Figure S9 --
  "Fluorescence emission profiles of 0.8 uM BSA with increasing
  concentration of (a) PFBA, (b) PFHxA, (c) PFHpA, (d) PFOA, (e)
  PFNA, (f) PFDA, (g) PFBS, (h) PFHxS, (i) PFOS, and (j) GenX in
  the absence of salt in 20 mM NaH2PO4-Na2HPO4 buffer (pH 7) at
  20 C."  The same 0.8 uM appears in the captions of Figures S5 and
  S6.  It is constant across every compound and every salt level.

  LIGAND.  The titration series is in the Figure S9 PANELS, one per
  compound.  There is NO TABLE of them anywhere in the paper or the
  SI: Table S1 is a literature comparison that does not include this
  study's own conditions, and Table 1 carries only the fitted
  constants.  So the series in ligand_conc_range was read off the
  figure.  This script uses only the TOP of each range, which is the
  endpoint of a plotted axis and the most reliable thing to read off
  a panel; the intermediate points are not needed for a maximum and
  are the part of that transcription to trust least.

WHAT THE RATIO TURNS OUT TO BE.  Between 12 and 18750 ligands per
albumin.  This is by a wide margin the most heavily loaded source in
the database -- Allendorf runs at 4e-4 to 3.8e-2, Chen 2025 at 2 --
and the top of it is roughly a million-fold above physiological
occupancy.  That is not a criticism of the measurement: a Trp
quenching titration has to drive the signal to a plateau, and
reaching one requires saturating the protein.  It is a statement
about what the resulting constant can be compared with, which is
why the ratio belongs in the database rather than being left blank.

A SEPARATE FINDING, recorded here because this is where it surfaced
and it is not ours to fix: Gao's own SI Table S1 lists Allendorf
2019's BSA concentrations as "3 ng/L, 1 ng/L, 0.1 ng/L, 0.05 ng/L,
0.025 ng/L".  Allendorf's SI Table SI-2.3 gives those same fourteen
numbers in g/L.  Gao's comparison table is off by nine orders of
magnitude on the protein concentration of the study it compares
itself against.  Our own Allendorf rows were verified directly
against Allendorf's SI and are unaffected.

Run:  python3 tk_review/scripts/fix_gao2025_ratio.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

BSA_UM = 0.8                      # Figure S9 caption
TO_UM = {"mM": 1000.0, "uM": 1.0, "nM": 1e-3}


def top_of_range(value, units):
    """Highest ligand concentration in uM, from '0-8 mM (Figure S9: ...)'."""
    # "0-8 mM (Figure S9: ...)" -> head "0-8 mM". The unit is repeated
    # inside the cell as well as in ligand_conc_units, so strip it from
    # the head rather than assume the head is bare digits.
    head = value.split("(")[0].strip()
    unit_in = re.search(r"\b(mM|uM|nM)\b", head)
    head = re.sub(r"\b(mM|uM|nM)\b", "", head).strip()

    m = re.match(r"\s*([\d.eE+-]+)\s*(?:-|to|–)\s*([\d.eE+-]+)\s*$", head)
    if m:
        hi = max(float(m.group(1)), float(m.group(2)))
    else:
        m2 = re.match(r"\s*([\d.eE+-]+)\s*$", head)
        if not m2:
            return None
        hi = float(m2.group(1))

    # If the cell names its own unit and the units column disagrees,
    # that is a transcription error and guessing which side is right
    # would be worse than refusing.
    declared = units.strip().split()[0] if units.strip() else ""
    if unit_in and declared and unit_in.group(1).lower() != declared.lower():
        return None
    f = TO_UM.get(unit_in.group(1) if unit_in else declared)
    return None if f is None else hi * f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"].startswith("Gao")
           and r["year"] == "2025"]
    seen, n, skipped = {}, 0, 0
    for r in hit:
        if float(r["protein_conc"]) != BSA_UM or r["protein_conc_units"] != "uM":
            skipped += 1
            continue
        top = top_of_range(r["ligand_conc_range"], r["ligand_conc_units"])
        if top is None:
            skipped += 1
            continue
        v = f"{top / BSA_UM:g}"
        if r["ligand_over_protein_max"] != v:
            n += 1
        r["ligand_over_protein_max"] = v
        seen.setdefault(r["chemical"], (top, v))

    print(f"  {len(hit)} Gao 2025 rows, BSA {BSA_UM} uM (Figure S9 caption)")
    print(f"  {n} ratios written, {skipped} skipped\n")
    print(f"  {'PFAS':16}{'top ligand uM':>15}{'ratio':>12}")
    for c, (top, v) in sorted(seen.items(), key=lambda kv: -float(kv[1][1])):
        print(f"  {c[:15]:16}{top:>15,.4g}{float(v):>12,.0f}")

    if args.dry_run:
        print("\n  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
