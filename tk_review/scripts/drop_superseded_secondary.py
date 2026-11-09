#!/usr/bin/env python3
"""Drop a secondary transcription once the primary source is in hand.

A "secondary transcription" row is one read out of a LATER paper's
comparison table, because the original had not been retrieved.  It is
a stopgap, and it carries a source_table note saying so.  When the
original is later ingested, the stopgap should go -- otherwise the
same measurement sits in the database twice and is counted twice by
anything that pools.

That has happened once.  Gao 2019's PFDoA/HSA Kd of 103 uM was read
out of Jia et al. 2022 Table 1, whose superscript 19 attributes it to
Gao; the note says "Gao et al. 2019 not retrieved".  Gao's own SI
Table S9 has since been ingested, 42 rows of it, and gives the same
compound the same value -- 103 +/- 12 uM.

It was invisible because of a NAME.  The primary came in as PFDoDA
and the stopgap as PFDoA, and those are the same chemical,
perfluorododecanoic acid.  A duplicate check keyed on the raw string
sees two different compounds; this one canonicalises first, through
scripts/pfas_names.py, which already knows the pair.

THE RULE APPLIED HERE, deliberately narrow:

  drop a secondary row only when a PRIMARY row exists for the same
  study, year, canonical compound and protein, AND reports the same
  value.

A disagreement is NOT a duplicate -- it means the later paper
misquoted the earlier one, or we misread one of them, and that is
worth looking at rather than silently resolving.  Such a row is
reported and kept.  Of the 12 secondary rows in the database only
this one is superseded at all; the other 11 are membrane partition
coefficients whose originals are still not in hand.

Run:  python3 tk_review/scripts/drop_superseded_secondary.py [--dry-run]
Idempotent.
"""
import argparse
import collections
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

_spec = importlib.util.spec_from_file_location(
    "pfas_names", os.path.join(HERE, "pfas_names.py"))
PN = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PN)

MARK = "secondary transcription"


def is_secondary(r):
    return MARK in r["source_table"].lower()


def key(r):
    return (r["study"], r["year"], PN.canon(r["chemical"]), r["protein"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    idx = collections.defaultdict(list)
    for r in rows:
        idx[key(r)].append(r)

    drop, conflict, orphan = [], [], 0
    for r in rows:
        if not is_secondary(r):
            continue
        prim = [x for x in idx[key(r)] if not is_secondary(x)]
        if not prim:
            orphan += 1
            continue
        if all(x["kd"] == r["kd"] and x["ka"] == r["ka"] for x in prim):
            drop.append(r)
        else:
            conflict.append((r, prim))

    print(f"  {sum(1 for r in rows if is_secondary(r))} secondary rows; "
          f"{orphan} with no primary yet\n")
    for r in drop:
        print(f"  DROP      {r['study'].replace(' et al.','')} {r['year']} "
              f"{r['chemical']} -> {PN.canon(r['chemical'])}, "
              f"{r['protein'][:30]}, Kd {r['kd']}")
        print(f"            superseded by the primary row, same value")
    for r, prim in conflict:
        print(f"  CONFLICT  {r['study'].replace(' et al.','')} {r['year']} "
              f"{r['chemical']}: secondary Kd {r['kd']} vs primary "
              f"{[x['kd'] for x in prim]}")
        print(f"            KEPT -- a disagreement is not a duplicate")
    if not drop and not conflict:
        print("  nothing superseded")

    if args.dry_run:
        print("\n  (dry run, nothing written)")
        return 0
    if drop:
        keep = [r for r in rows if r not in drop]
        with open(DB, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=flds)
            w.writeheader()
            w.writerows(keep)
        print(f"\n  {len(rows)} -> {len(keep)} rows")
    return 1 if conflict else 0


if __name__ == "__main__":
    sys.exit(main())
