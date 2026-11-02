#!/usr/bin/env python3
"""The albumin coefficients, cut down to something you can hand someone.

MEASURED ROWS ARE INCLUDED, and the first version of this file wrongly
left them out. It was built by filtering db/albumin_kd_to_k.csv, which
is a CONVERSION table and so holds only rows reporting a Kd or Ka. Four
sources never report an affinity because they measured the partition
coefficient directly -- Fischer 2024, Allendorf 2019, Bischel 2011 and
Qin 2024 -- so they have nothing to convert and are absent from that
table by design. That filter suits the conversion analysis and does not
suit a reader who simply wants the best available coefficient: it drops
the strongest data in the set, including the only values measured at a
genuinely negligible loading (Fischer holds PFAS:protein below 0.004).

WHY A SEPARATE FILE. db/albumin_kd_to_k.csv is 302 rows and 28 columns and
it is built to be argued with: every refusal, every diagnostic, every
provenance string. That is the right shape for the person who owns the
analysis and the wrong shape for the person who was sent it.

WHAT WAS CUT, AND THE ONE THING THAT WAS NOT. Rows: every albumin
coefficient that needs no assumption from the reader -- measured where a
paper measured one, converted where a paper supplied both an affinity
and a stoichiometry. Columns: 28 down to 8. What is NOT cut is the statement of
how solid each number is -- the five isotherm diagnostics collapse into
ONE plain-language column, `how_solid`. Stripping that too would hand
someone 60 numbers that all move one for one with an assumption they
cannot see, which is worse than sending nothing.

FOUR TIERS, and the shape of them is the finding:

    measured directly  90   Fischer, Allendorf, Bischel - no conversion
    both measured       4   Moro 2022, stoichiometry counted by ITC
    n assumed, Kd ok   11   Chen 2025, Kd inside the tested range
    n assumed, Kd ex   49   Chen 2025, Kd fitted beyond the range tested

Only four CONVERTED rows in the whole albumin literature carry a
measured stoichiometry. Anyone using those numbers should know it before
they use them, which is why the column is first-class and not a
footnote.

QIN 2024 IS THE ONE MEASURED SOURCE HELD BACK. It reports albumin in two
components, a specific and a non-specific term, and quotes the specific
one at half-saturation (its Eq. 15, K = C_max/2 K_d). Using it needs the
two summed in linear space after the specific term is lifted by
log10(2); the stored row is neither, and would travel about 0.3 log low.
The corrected value belongs to the comparison analysis, not to a table
someone reads off.

Writes db/albumin_logK_shareable.xlsx (Notes + Coefficients).

Run:  python3 -I tk_review/scripts/albumin_shareable.py
"""
import collections
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_kd_to_k.csv")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")

# Human only by default. Albumin's own species difference is tiny -- HSA
# and BSA are 66,437 and 66,398 g/mol by sequence, 0.00025 log apart --
# but the SOURCES are not interchangeable, and mixing them puts bovine
# equilibrium dialysis beside human single-site fits in one column with
# nothing but a species cell to separate them. Pass a species to change
# it, or "all" to keep every row.
SPECIES = "human"

# No Species column: the species filter leaves one, and `protein` names
# it on every row ("human serum albumin (HSA)"), so a Species column
# would be 92 identical cells.
FIELDS = [("pfas_canonical", "PFAS"), ("protein", "Protein"),
          ("kd_uM", "Kd (uM)"),
          ("n_sites_used", "Sites per protein (n)"),
          ("logK", "log K (L water / kg protein)"),
          ("how_solid", "How solid"), ("source", "Source")]

MEASURED = "measured directly - no conversion, no stoichiometry needed"
BOTH = "converted; both measured - n counted by calorimetry"
OK = "converted; n assumed = 1 by the paper's single-site model; Kd inside the range tested"
EXTRAP = ("converted; n assumed = 1 by the paper's single-site model; Kd fitted "
          "beyond the highest concentration tested")

NOTES = [
    ("What this is",
     "Protein-water partition coefficients for PFAS on serum albumin. "
     "log K is log10 of litres of water per kg of protein, at zero "
     "loading. Some were measured directly; the rest were converted from "
     "each paper's own dissociation constant. The 'How solid' column "
     "says which, for every row."),
    ("Measured or converted",
     "Prefer a measured row where one exists for your compound. Fischer "
     "2024 measured the coefficient itself and needed no stoichiometry, "
     "so none of the assumptions below apply to its rows. It held "
     "PFAS:protein under 0.004, which puts those values within 0.002 log "
     "of zero loading."),
    ("Species",
     "Human serum albumin only, so there is no species column. 62 bovine "
     "rows exist (Allendorf 2019, Bischel 2011, and Fischer again) and "
     "are left out: every compound they cover is also measured in human, "
     "so they add no coverage at all. They are a cross-check on the "
     "human values, not extra data - ask if you want them."),
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
     "Four values, best first. 'both measured' means the paper supplies the affinity "
     "AND counts the sites. The other two mean the paper fitted a "
     "single-site model, which fixes n at 1 by construction rather than "
     "by measurement, and they differ in whether the fitted Kd sits "
     "inside the concentrations the experiment actually tested."),
    ("What is NOT here",
     "238 further albumin rows that report an affinity. 207 of them give "
     "no usable site count, so a coefficient can only be had by assuming "
     "n = 1 yourself; 31 report no single affinity at all, 23 of those "
     "because the paper looked and found no measurable binding. Ask if "
     "you want them - they are kept, with reasons, in the full table."),
    ("One measured source held back",
     "Qin 2024 reports albumin in two components and quotes the specific "
     "one at half-saturation, so using it means summing the two and "
     "lifting the specific term by log10(2). The raw value travels about "
     "0.3 log low, so it is left out rather than shipped uncorrected."),
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


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def is_serum_albumin(protein):
    p = protein.lower()
    if "serum albumin" not in p:
        return False
    return not any(w in p for w in ("globulin", "plasma", "whole serum",
                                    "serum)", "immobilised", "immobilized"))


def measured():
    """The albumin coefficients that were measured, not converted.

    These never appear in the Kd-to-K table: they report no affinity,
    because the experiment produced the coefficient itself. Qin 2024 is
    excluded for the reason in the module docstring.
    """
    B = _load("baf", "build_binding_affinity_table.py")
    N = _load("pfas_names", "pfas_names.py")
    out = []
    for r in B.build(filter_compounds=False)[1]:
        if r["measure"] != "log D" or not r["logK"]:
            continue
        if not is_serum_albumin(r["protein"]):
            continue
        if r["source"].startswith("Qin 2024"):
            continue
        out.append({
            "pfas_canonical": N.canon(r["pfas"]), "protein": r["protein"],
            "species": r["species"], "kd_uM": "", "n_sites_used": "",
            "logK": r["logK"], "how_solid": MEASURED,
            "source": r["source"]})
    return out


def main():
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    rows = [r for r in csv.DictReader(open(SRC, encoding="utf-8"))
            if r["convertible"] == "yes"]
    for r in rows:
        r["how_solid"] = tier(r)
    rows = measured() + rows
    want = (sys.argv[1] if len(sys.argv) > 1 else SPECIES).lower()
    dropped = collections.Counter()
    if want != "all":
        for r in rows:
            if r["species"].lower() != want:
                dropped[(r["species"], r["source"])] += 1
        rows = [r for r in rows if r["species"].lower() == want]
    order = {MEASURED: 0, BOTH: 1, OK: 2, EXTRAP: 3}
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

    print(f"{len(rows)} rows, {len(FIELDS)} columns, species = {want} -> "
          f"{os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    for t, n in collections.Counter(r["how_solid"] for r in rows).most_common():
        srcs = sorted({r["source"] for r in rows if r["how_solid"] == t})
        print(f"  {n:3d}  {t[:62]:64s} {', '.join(srcs)}")
    print(f"  {len({r['pfas_canonical'] for r in rows})} distinct compounds")
    if dropped:
        print(f"  left out by the species filter:")
        for (sp, src), n in sorted(dropped.items()):
            print(f"     {n:3d}  {sp:8s} {src}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
