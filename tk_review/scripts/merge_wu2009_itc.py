#!/usr/bin/env python3
"""Collapse Wu 2009's two ITC rows into the one binding event they are.

THE PAPER.  Wu L-L, Gao H-W, Gao N-Y, Chen F-F, Chen L, "Interaction of
perfluorooctanoic acid with human serum albumin", BMC Struct Biol
2009;9:31, doi 10.1186/1472-6807-9-31.  Table 1 fits ONE two-step
sequential model to ONE titration:

    step 1   N1 = 1    K_b = 24.7 +/- 6.5  x 10^3 M^-1   Kd   40.49 uM
    step 2   N2 = 12   K_b = 0.55 +/- 0.14 x 10^3 M^-1   Kd 1818.18 uM

These are two CLASSES of site on one protein, not two measurements.
The number that describes the protein is their sum,

    K0 = (1e9/M) x SUM_i n_i/Kd_i
       = (1e9/66463) x (1/40.49 + 12/1818.18)
       -> log K 2.673

Ingested as two rows they entered the PFOA pool twice, at log K 2.570
and 1.997, neither of which is a partition coefficient on its own:
2.570 describes the single buried site as though the other twelve did
not exist, and 1.997 the twelve surface sites as though the buried one
did not.

This script merges them into one row carrying both constants, which is
the form the converter's multi-class path already expects -- the same
shape as Woodcroft 2010's "3.1e-6, 26.2e-6, 52.6e-6".  The one thing
that had to change in the converter is that the old path assumed one
site per class; Wu's classes are 1 and 12, so site_list() now reads a
per-class list and the sum uses it.  Woodcroft's rows are untouched.

The Temkin row is NOT removed here.  It stays in the database as the
record of what the paper reports, and the converter now declines it on
its own terms (see the Temkin check in build_binding_affinity_table.py
and the working in scripts/audit_wu2009_temkin.py), which is the right
place for that refusal: the defect is in the quantity, not in scope.

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/merge_wu2009_itc.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "two classes merged 2026-11-13"


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())

    itc = [r for r in rows
           if "Wu" in (r.get("study") or "") and r.get("year") == "2009"
           and "isothermal titration" in (r.get("method") or "").lower()]
    if any(MARK in (r.get("source_table") or "") for r in itc):
        # Already merged. Repair the two fields that must stay cleared
        # for the row to resolve as Kd, then stop.
        r = itc[0]
        want = {"ka": "", "ka_units": "",
                "ligand_conc_range": ("2500 (syringe); fitted over "
                                      "PFOA:HSA 0 to 20 across both "
                                      "steps"),
                "ligand_over_protein_max": "20"}
        if any((r.get(k) or "") != v for k, v in want.items()):
            r.update(want)
            with open(DB, "w", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields)
                w.writeheader()
                w.writerows(rows)
            print("  already merged; repaired fields that must hold "
                  "for the row to resolve as Kd and to describe the "
                  "whole titration")
        else:
            print("  already merged; nothing to do")
        return 0
    if len(itc) != 2:
        print(f"  expected 2 ITC rows, found {len(itc)}", file=sys.stderr)
        return 1

    step1 = min(itc, key=lambda r: float(r["kd"]))
    step2 = max(itc, key=lambda r: float(r["kd"]))
    keep, drop = step1, step2

    keep["kd"] = f"{float(step1['kd']):g}, {float(step2['kd']):g}"
    keep["kd_units"] = "uM"
    # ka is deliberately CLEARED.  The converter resolves a row's
    # measure from whichever of ka/kd is populated, and prefers Ka; with
    # a list in ka it read measure = Ka, found no single value to invert
    # and gave up ("convertible = no").  The multi-class path reads the
    # kd cell, so the row has to present as Kd -- the same shape as
    # Woodcroft 2010's three-constant rows.  The K_b values are kept in
    # the note, where they are a record rather than an input.
    keep["ka"] = ""
    keep["ka_units"] = ""
    # The merged row spans the WHOLE titration, so it takes step 2's
    # loading. Step 1 was fitted below PFOA:HSA 1.5 and step 2 from 1.5
    # to 20 (Table 1's c_L0/c_M0 column); keeping step 1's 1.5 would
    # describe a third of the experiment the row now stands for.
    keep["ligand_conc_range"] = ("2500 (syringe); fitted over PFOA:HSA 0 "
                                 "to 20 across both steps")
    keep["ligand_over_protein_max"] = "20"
    keep["n_binding_sites"] = (
        "1, 12 (by ITC, two classes: Table 1 gives N1 = 1 buried in the "
        "hydrophobic intracavity and N2 = 12 on the exposed outer "
        "surface; Fig. 3B's fit box reads N1 0.892, N2 12.0 +/- 2.01)")
    keep["source_table"] = (
        "Table 1, both steps | read from the full text, primary "
        f"({MARK}). The two rows this replaces carried the steps "
        "separately at log K 2.570 and 1.997; neither is a partition "
        "coefficient alone, and together they double-counted the study "
        "in the PFOA pool. Summed per the project's multi-class rule, "
        "K0 = (1e9/M) x SUM n_i/Kd_i, they give log K 2.673")
    note = keep.get("assay_volume_note") or ""
    if MARK not in note:
        keep["assay_volume_note"] = (
            note + " | TWO CLASSES, ONE EVENT: Table 1's steps 1 and 2 "
            "are one two-step sequential fit to one titration (ITC cell "
            "1.4685 mL holding 0.010 mM HSA, 6.0 uL pulses of 2.50 mM "
            "PFOA, 37 C). Step 1 is fitted below PFOA:HSA 1.5 and step 2 "
            "over 1.5 to 20, which are the c_L0/c_M0 column of Table 1, "
            "not separate experiments.").strip(" |")

    rows = [r for r in rows if r is not drop]
    with open(DB, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  merged 2 ITC rows into 1: kd = {keep['kd']} uM, "
          f"n = 1, 12")
    print("  Temkin row left in place; the converter declines it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
