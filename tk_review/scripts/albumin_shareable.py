#!/usr/bin/env python3
"""The albumin coefficients, cut down to something you can hand someone.

WHY A SEPARATE FILE. db/albumin_kd_to_k.csv is 302 rows and 28 columns and
it is built to be argued with: every refusal, every diagnostic, every
provenance string. That is the right shape for the person who owns the
analysis and the wrong shape for the person who was sent it.

WHAT WAS CUT, AND THE ONE THING THAT WAS NOT. Rows: only the 64 that
reach a partition coefficient without the reader supplying an
assumption. Columns: 28 down to 8. What is NOT cut is the statement of
how solid each number is -- the five isotherm diagnostics collapse into
ONE plain-language column, `how_solid`. Stripping that too would hand
someone 60 numbers that all move one for one with an assumption they
cannot see, which is worse than sending nothing.

THREE TIERS, and the shape of them is the finding:

    both measured      4   Moro 2022, stoichiometry counted by ITC
    n assumed, Kd ok  11   Chen 2025, Kd inside the tested range
    n assumed, Kd ex  49   Chen 2025, Kd fitted beyond the range tested

Only four rows in the whole albumin literature carry a measured
stoichiometry. Anyone using these numbers should know that before they
use them, which is why the column is first-class and not a footnote.

Writes db/albumin_logK_shareable.xlsx (Notes + Coefficients).

Run:  python3 -I tk_review/scripts/albumin_shareable.py
"""
import collections
import csv
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_kd_to_k.csv")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")

FIELDS = [("pfas_canonical", "PFAS"), ("protein", "Protein"),
          ("species", "Species"), ("kd_uM", "Kd (uM)"),
          ("n_sites_used", "Sites per protein (n)"),
          ("logK", "log K (L water / kg protein)"),
          ("how_solid", "How solid"), ("source", "Source")]

BOTH = "both measured - n counted by calorimetry"
OK = "n assumed = 1 by the paper's single-site model; Kd inside the range tested"
EXTRAP = ("n assumed = 1 by the paper's single-site model; Kd fitted "
          "beyond the highest concentration tested")

NOTES = [
    ("What this is",
     "Protein-water partition coefficients for PFAS on serum albumin, "
     "converted from each paper's own dissociation constant. log K is "
     "log10 of litres of water per kg of protein, at zero loading."),
    ("The conversion",
     "log K = log10(1e9 / M) + log n - log Kd(uM), with M = 66,463 g/mol "
     "for albumin. It is exact algebra from the Langmuir isotherm, not a "
     "correlation or a fit."),
    ("Why n matters so much",
     "log K moves one for one with log n. An n of 2 rather than 1 shifts "
     "every value by +0.30 log. Papers that count sites on albumin report "
     "anywhere from 1 to 9, which is 0.95 log of spread - so the 'How "
     "solid' column is not a formality."),
    ("How solid",
     "Three values. 'both measured' means the paper supplies the affinity "
     "AND counts the sites. The other two mean the paper fitted a "
     "single-site model, which fixes n at 1 by construction rather than "
     "by measurement, and they differ in whether the fitted Kd sits "
     "inside the concentrations the experiment actually tested."),
    ("What is NOT here",
     "238 further albumin rows. 207 report an affinity with no usable "
     "site count, so a coefficient can only be had by assuming n = 1 "
     "yourself; 31 report no single affinity at all, 23 of those because "
     "the paper looked and found no measurable binding. Ask if you want "
     "them - they are kept, with reasons, in the full table."),
    ("How reproducible is any of this",
     "Two laboratories measuring the same albumin Kd by the same method "
     "differ by 1.45 log, 45x on PFOA. Gao et al. 2019 reports the "
     "published range for that constant as two orders of magnitude. Treat "
     "a single value as reproducible to about an order, not to its "
     "decimal places."),
    ("Where it comes from",
     "Built from the project's protein-binding database. The full table, "
     "with every diagnostic and every refusal, is db/albumin_kd_to_k.csv."),
]


def tier(r):
    if r["n_basis"] not in ("assumed", ""):
        return BOTH
    try:
        inside = float(r["kd_over_top_conc"]) <= 1.0
    except (TypeError, ValueError):
        inside = False
    return OK if inside else EXTRAP


def main():
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    rows = [r for r in csv.DictReader(open(SRC, encoding="utf-8"))
            if r["convertible"] == "yes"]
    for r in rows:
        r["how_solid"] = tier(r)
    order = {BOTH: 0, OK: 1, EXTRAP: 2}
    rows.sort(key=lambda r: (order[r["how_solid"]], r["pfas_canonical"]))

    wb = Workbook()
    wb.remove(wb.active)

    ws = wb.create_sheet("Notes")
    ws["A1"] = "PFAS - serum albumin partition coefficients"
    ws["A1"].font = Font(bold=True, size=14)
    for i, (k, v) in enumerate(NOTES, start=3):
        ws.cell(row=i, column=1, value=k).font = Font(bold=True)
        ws.cell(row=i, column=2, value=v).alignment = Alignment(
            wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 96

    ws = wb.create_sheet("Coefficients")
    ws.append([label for _, label in FIELDS])
    for r in rows:
        line = []
        for key, _ in FIELDS:
            v = r[key]
            if key in ("kd_uM", "logK", "n_sites_used") and str(v).strip():
                try:
                    v = float(v)
                except ValueError:
                    pass
            line.append(v)
        ws.append(line)
    for c, (_, label) in enumerate(FIELDS, start=1):
        cell = ws.cell(row=1, column=c)
        cell.fill = PatternFill("solid", fgColor="2A78D6")
        cell.font = Font(bold=True, color="FFFFFF")
        cell.alignment = Alignment(vertical="center", wrap_text=True)
        width = max([len(label)] + [len(str(r[FIELDS[c - 1][0]]))
                                    for r in rows])
        ws.column_dimensions[get_column_letter(c)].width = \
            min(max(width + 2, 11), 64)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    wb.save(OUT)

    print(f"{len(rows)} rows, {len(FIELDS)} columns -> "
          f"{os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    for t, n in collections.Counter(r["how_solid"] for r in rows).most_common():
        srcs = sorted({r["source"] for r in rows if r["how_solid"] == t})
        print(f"  {n:3d}  {t[:62]:64s} {', '.join(srcs)}")
    print(f"  {len({r['pfas_canonical'] for r in rows})} distinct compounds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
