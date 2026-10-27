#!/usr/bin/env python3
"""The five best sources on PFAS binding to an INDIVIDUAL protein.

db/binding_affinity.xlsx carries every measurement and ranks its "Major
studies" tab by row count alone.  That ranking puts Ryu et al. 2024 first,
and Ryu reports unbound fractions in whole tissue homogenate -- total protein
plus phospholipid, no protein resolved.  It is a large and useful dataset and
it answers a different question.

This file ranks on the question actually being asked: how many rows report a
quantitative binding or partition measurement against a NAMED SINGLE PROTEIN.
Whole plasma, serum, tissue homogenate, culture medium, cells and liposomes
do not count towards the ranking, because none of them isolates a protein.

    source              named-protein rows   compounds   quantity
    Fischer 2024                        83          29   log D -- HSA, BSA, globulin
    Chen 2025                           57          57   Kd -- HSA
    Muensterman 2026                    56          56   log D -- HSA, biomimetic
    Gao 2025                            42          10   Ka -- BSA, NaCl series
    Qin 2024 (plasma)                   32          16   log D -- BSA, by concentration
    ---- cut ----
    Starnes 2024                        32           8   Kd -- albumin, 4 species
    Starnes 2025                        28          28   Ka -- HSA, DSF
    Allendorf 2019                      24          14   log D -- BSA, dialysis

Gao et al. 2025 replaces Ryu as the fifth.  It is the next-ranked source on
this basis and earns the place on its own: it is the only source in the whole
database that varies IONIC STRENGTH (0-500 mM added NaCl, including a
physiological 150 mM point), and it supplies a second independent affinity
constant where the set is otherwise short of them -- without it, Chen 2025
would be the only dissociation constant among the five.

THE FIFTH PLACE IS A TIE, and it is broken by a stated rule rather than by
accident.  Starnes 2024 and Qin 2024 (plasma) both supply 32 named-protein
rows.  The tie-break is distinct compounds, which Qin wins 16 to 8; Qin's
rows are also the only concentration-resolved albumin measurements in the
database, which is the quantity this project's section 7 turns on.  Ranking
by row count alone would have decided this by dictionary insertion order.

Two sources that just miss, worth knowing about: Starnes et al. 2025 covers
more compounds than anything here except Chen (28), but measures by
differential scanning fluorimetry at the protein's melting transition, which
sits about three orders of magnitude off the dialysis values -- a known and
documented offset, not a disagreement to average over.  Allendorf et al.
2019 is arguably the cleanest measurement of the lot (equilibrium dialysis at
37 C, ligand:protein below 0.1 mol/mol, two albumin batches) but reports a
partition coefficient, a quantity three of the five already supply.

ROWS ARE NOT FILTERED BY MATRIX.  Selection is by source; once a source is in,
all of its rows come with it, mixtures included, because dropping a paper's
plasma or serum rows would silently discard measurements the user did not ask
to lose.  The `protein` column says what each row was measured in.

COLUMNS.  The nine requested, plus the PFAS and protein concentrations each
measurement was made at.  Each concentration appears twice: a micromolar
number where unit arithmetic alone reaches one, and the paper's own wording
beside it.  A blank number next to a populated "as stated" cell means the
source gave a mass concentration, a dilution factor, a stationary phase or a
molar ratio -- something that cannot be converted without assuming a molecular
weight or a volume.  Nothing here is inferred.

Run:  python3 tk_review/scripts/build_binding_major5.py
"""
import collections
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "binding_major5.xlsx")

_spec = importlib.util.spec_from_file_location(
    "baf", os.path.join(HERE, "build_binding_affinity_table.py"))
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)

# The nine columns asked for, then the concentrations.
FIELDS = ["pfas", "species", "protein", "measure", "value", "units",
          "kd_uM", "logK", "source",
          "pfas_conc_uM", "pfas_conc_as_stated",
          "protein_conc_uM", "protein_conc_as_stated"]

# A named single protein. Albumin, globulin and the named binding proteins
# count; a matrix that merely CONTAINS protein does not.
PROTEIN_WORDS = ("albumin", "globulin", "transthyretin", "fabp",
                 "fatty acid-binding", "glycoprotein")


def is_named_protein(protein):
    p = protein.lower()
    if any(w in p for w in ("plasma", "serum albumin")) and \
            "serum albumin" not in p:
        return False
    return any(w in p for w in PROTEIN_WORDS)


def rank(out):
    """Named-protein rows per source, and the ordering over them.

    Rows first, distinct compounds as the tie-break. The tie is real and not
    hypothetical: Starnes 2024 and Qin 2024 (plasma) both supply 32, and
    they differ 8 compounds against 16. Leaving the order to Counter's
    insertion order would have decided a selection by accident.
    """
    n = collections.Counter(r["source"] for r in out
                            if is_named_protein(r["protein"]))
    cmpds = {s: len({r["pfas"] for r in out
                     if r["source"] == s and is_named_protein(r["protein"])})
             for s in n}
    order = sorted(n, key=lambda s: (-n[s], -cmpds[s], s))
    return n, cmpds, order


def main():
    rows, out, skipped, qualitative, dropped = B.build()
    n, cmpds, order = rank(out)
    top5 = order[:5]
    tot = collections.Counter(r["source"] for r in out)

    print("Ranked on rows measured against a named single protein,")
    print("ties broken on distinct compounds:")
    print(f"  {'source':22s}{'protein rows':>13s}{'all rows':>10s}"
          f"{'compounds':>11s}  quantities")
    for s in order[:9]:
        prot = [r for r in out
                if r["source"] == s and is_named_protein(r["protein"])]
        meas = collections.Counter(r["measure"] for r in prot)
        mark = "  <- kept" if s in top5 else ""
        print(f"  {s:22s}{n[s]:>13d}{tot[s]:>10d}{cmpds[s]:>11d}  "
              + ", ".join(f"{k} {v}" for k, v in meas.most_common(3)) + mark)
    print()

    cond = {c["id"]: c for c in B.conditions(out)}
    sel = []
    for r in out:
        if r["source"] not in top5:
            continue
        c = cond[r["id"]]
        sel.append({f: (r.get(f) if f in r else c.get(f, "")) for f in FIELDS})
    sel.sort(key=lambda r: (top5.index(r["source"]), r["protein"], r["pfas"]))

    by = collections.Counter(r["source"] for r in sel)
    print(f"{len(sel)} rows from {len(top5)} sources:")
    for s in top5:
        print(f"  {by[s]:5d}  {s}")
    print()
    for f in ("pfas_conc_uM", "protein_conc_uM", "pfas_conc_as_stated",
              "protein_conc_as_stated", "kd_uM", "logK"):
        got = sum(1 for r in sel if str(r[f]).strip())
        print(f"  {f:24s} populated {got:4d} / {len(sel)}")
    write(sel, top5, n, tot)
    print(f"\nwrote {os.path.abspath(OUT)}")
    return 0


NOTES = [
    ("What this file is", "The five sources that report PFAS binding to a "
     "NAMED SINGLE PROTEIN, ranked by how many such measurements each "
     "supplies. Whole plasma, serum, tissue homogenate, culture medium, "
     "cells and liposomes do not count towards the ranking."),
    ("Why not Ryu 2024", "It is the largest contributor to the full "
     "database (320 rows) and reports unbound fractions in whole tissue "
     "homogenate -- total protein plus phospholipid, no protein resolved. "
     "A row count is not a ranking of usefulness for this question."),
    ("Why Gao 2025", "Next-ranked on this basis, and the better of the two "
     "candidates at that rank: it is the only source in the database that "
     "varies ionic strength (0-500 mM added NaCl), and without it Chen 2025 "
     "would be the only dissociation constant among the five. Runner-up: "
     "Allendorf 2019, cleaner but reporting a partition coefficient that "
     "three of the five already supply."),
    ("Row selection", "By source, not by matrix. Once a source is in, all of "
     "its rows come with it, mixtures included; the `protein` column says "
     "what each row was measured in."),
    ("", ""),
    ("pfas", "Chemical as the paper names it."),
    ("species", "Species the protein came from."),
    ("protein", "Protein or matrix, as the paper describes it."),
    ("measure", "Kd, Ka, log D (protein-water), log K_membrane, f_unbound. "
     "A '(computed)' suffix means the source reports the value as predicted "
     "rather than measured."),
    ("value", "The number as printed, including ranges and bounds."),
    ("units", "Canonicalised unit for that value."),
    ("kd_uM", "Everything that converts EXACTLY to a dissociation constant "
     "in micromolar. Blank is a deliberate state, never a gap: log D is not "
     "converted to Kd, and a Kd is not derived from an affinity without a "
     "stoichiometry the paper reports."),
    ("logK", "log10 of the protein-water partition coefficient, L_water per "
     "kg protein. Two exact routes only: the paper measured it, or reported "
     "BOTH an affinity and a stoichiometry."),
    ("source", "First author and year. 'Qin 2024 (plasma)' is doi "
     "10.1021/acs.est.3c10824; the other Qin 2024 paper, the baseline-"
     "toxicity one, is a different study and is labelled separately."),
    ("", ""),
    ("pfas_conc_uM", "Dosed PFAS concentration in micromolar where unit "
     "arithmetic alone reaches one. A range is carried as its midpoint only "
     "where the source states a range in molar units."),
    ("pfas_conc_as_stated", "The paper's own wording, verbatim."),
    ("protein_conc_uM", "Protein concentration in micromolar, same rule."),
    ("protein_conc_as_stated", "The paper's own wording, verbatim. Blank "
     "numbers beside a populated cell here are the cases that cannot be "
     "converted without assuming a molecular weight, a volume or a dilution "
     "-- a mass concentration (mg/mL), a stationary phase, or a molar "
     "ratio. None is inferred."),
    ("", ""),
    ("Why concentration matters", "PFAS binding to albumin is not linear. "
     "The apparent affinity falls about 2x per decade of free concentration "
     "(Qin 2024 Table 1), so two constants measured at different "
     "ligand:protein ratios are not comparable even for the same compound "
     "and protein. These columns are what makes that check possible."),
]


def write(sel, top5, n, tot):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    wb.remove(wb.active)
    ws = wb.create_sheet("Binding")
    ws.append(FIELDS)
    for r in sel:
        row = []
        for f in FIELDS:
            v = r[f]
            if f in B.NUMERIC and str(v).strip():
                try:
                    v = float(v)
                except ValueError:
                    pass
            row.append(v)
        ws.append(row)
    for c, f in enumerate(FIELDS, start=1):
        cell = ws.cell(row=1, column=c)
        cell.fill = PatternFill("solid", fgColor="2A78D6")
        cell.font = Font(bold=True, color="FFFFFF")
        cell.alignment = Alignment(vertical="center")
        width = max([len(f)] + [len(str(r[f])) for r in sel[:500]])
        ws.column_dimensions[get_column_letter(c)].width = \
            min(max(width + 2, 10), 44)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions

    ws = wb.create_sheet("Sources")
    ws.append(["source", "named-protein rows", "rows in this file",
               "compounds", "quantities"])
    by = collections.Counter(r["source"] for r in sel)
    for s in top5:
        sub = [r for r in sel if r["source"] == s]
        meas = collections.Counter(r["measure"] for r in sub)
        ws.append([s, n[s], by[s], len({r["pfas"] for r in sub}),
                   ", ".join(f"{k} {v}" for k, v in meas.most_common())])
    for c in range(1, 6):
        ws.cell(row=1, column=c).fill = PatternFill("solid", fgColor="2A78D6")
        ws.cell(row=1, column=c).font = Font(bold=True, color="FFFFFF")
    for c, w in zip("ABCDE", (24, 20, 18, 12, 52)):
        ws.column_dimensions[c].width = w

    ws = wb.create_sheet("Notes")
    ws["A1"] = "binding_major5.xlsx -- what is in it and why"
    ws["A1"].font = Font(bold=True, size=13)
    for i, (k, v) in enumerate(NOTES, start=3):
        ws.cell(row=i, column=1, value=k).font = Font(bold=True)
        ws.cell(row=i, column=2, value=v).alignment = \
            Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 96
    wb.save(OUT)


if __name__ == "__main__":
    sys.exit(main())
