#!/usr/bin/env python3
"""The leading sources on PFAS binding to an INDIVIDUAL protein.

db/binding_affinity.xlsx carries every measurement and ranks its "Major
studies" tab by row count alone.  That ranking puts Ryu et al. 2024 first,
and Ryu reports unbound fractions in whole tissue homogenate -- total protein
plus phospholipid, no protein resolved.  It is a large and useful dataset and
it answers a different question.

This file ranks on the question actually being asked: how many rows report a
MEASURED binding or partition value against a NAMED SINGLE PROTEIN.  Whole
plasma, serum, tissue homogenate, culture medium, cells and liposomes do not
count towards the ranking, because none of them isolates a protein; nor do
docking scores and predicted values, because neither is a measurement.

    rank  source              measured  compounds  quantity
       1  Fischer 2024              83         29  log D -- HSA, BSA, globulin
       2  Chen 2025                 57         57  Kd -- HSA
       3  Muensterman 2026          56         56  log D -- HSA, biomimetic
       4  Gao 2025                  42         10  Ka -- BSA, NaCl series
       5  Qin 2024 (plasma)         32         16  log D -- BSA, by concentration
       6  Starnes 2024              32          8  Kd -- albumin, 4 species
       7  Starnes 2025              28         28  Ka -- HSA, DSF
       8  Allendorf 2019            24         14  log D -- BSA, dialysis
       9  Bischel 2011              22         11  log D + f_unbound -- BSA
      10  Wu 2024                   10         10  Ka -- HSA
    ---- cut at 10 measured rows ----
      11  Moro 2022                  8          4  Kd -- HSA
      12  Peng 2024                  8          4  Ka -- HSA

TWO RULES DECIDE THE SELECTION, both stated so neither is left to accident.

Ranking is by measured named-protein rows, with distinct compounds ordering
any source that ties.  The tie is not hypothetical: Starnes 2024 and Qin 2024
(plasma) both supply 32, and an earlier version of this script resolved that
by dictionary insertion order, which silently dropped one of them.

A tie AT THE CUT keeps everyone in it.  A rank boundary that splits two
sources with identical counts is an arbitrary decision dressed as a ranking,
so the cut widens to the whole tied group.  At TOP_N = 10 nothing ties at the
boundary -- rank 10 has 10 rows and rank 11 has 8 -- so the rule does not
bind here; it is in the code for the next time the database grows.

Counting only measured rows changes one line and no selection.  Wu et al.
2024 reports 10 affinities and 10 docking scores; on either count it ranks
tenth.  Its affinities are worth treating with care rather than pooling: its
PFOA value of 103 M-1 is three to four orders below every other measurement
of PFOA on albumin in this database.

ROWS ARE NOT FILTERED BY MATRIX.  Selection is by source; once a source is
in, all of its rows come with it, mixtures and computed values included,
because dropping a paper's plasma or serum rows would silently discard
measurements the user did not ask to lose.  The `protein` column says what
each row was measured in and `measure` marks anything computed.

COLUMNS.  The nine requested, plus the PFAS and protein concentrations each
measurement was made at.  Each concentration appears twice: a micromolar
number where unit arithmetic alone reaches one, and the paper's own wording
beside it.  A blank number next to a populated "as stated" cell means the
source gave a dilution factor, a stationary phase or a molar ratio --
something that cannot be converted without assuming a molecular weight or a
volume.  Albumin mass concentrations ARE converted, at 66,463 g/mol, which
is arithmetic with one published constant rather than an assumption about
the sample.  Nothing here is inferred.

Run:  python3 tk_review/scripts/build_binding_major_sources.py
"""
import collections
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP_N = 10          # how many sources to keep; ties at the cut are all kept
OUT = os.path.join(HERE, "..", "db", "binding_major_sources.xlsx")

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


def is_measured(rec):
    """A docking score or a predicted value is not a measurement of binding.

    Such rows stay in the file where their source qualifies -- the `measure`
    column marks them -- but they do not earn a source its place.
    """
    return "comput" not in rec["measure"]


def rank(out):
    """Named-protein rows per source, and the ordering over them.

    Measured rows decide the rank; distinct compounds order sources that
    tie. Both tie-breaks are stated rather than left to dictionary order,
    which on the first run of this script silently decided a selection:
    Starnes 2024 and Qin 2024 (plasma) both supply 32 rows.
    """
    n = collections.Counter(r["source"] for r in out
                            if is_named_protein(r["protein"])
                            and is_measured(r))
    with_comp = collections.Counter(r["source"] for r in out
                                    if is_named_protein(r["protein"]))
    cmpds = {s: len({r["pfas"] for r in out if r["source"] == s
                     and is_named_protein(r["protein"]) and is_measured(r)})
             for s in n}
    order = sorted(n, key=lambda s: (-n[s], -cmpds[s], s))
    return n, with_comp, cmpds, order


def select(n, order, top_n=TOP_N):
    """The top N sources, plus everything tied with the last of them.

    A rank cut that splits two sources with identical counts is an arbitrary
    decision dressed as a ranking, so the cut is widened to include the
    whole tied group.
    """
    if len(order) <= top_n:
        return list(order), []
    cut = n[order[top_n - 1]]
    kept = [s for s in order if n[s] >= cut]
    tied = [s for s in kept[top_n:]] if len(kept) > top_n else []
    return kept, tied


def main():
    rows, out, skipped, qualitative, dropped = B.build()
    n, with_comp, cmpds, order = rank(out)
    kept, tied_in = select(n, order)
    tot = collections.Counter(r["source"] for r in out)

    print(f"Ranked on MEASURED rows against a named single protein. Top "
          f"{TOP_N}, plus any source tied with the {TOP_N}th.")
    print(f"  {'':4s}{'source':22s}{'measured':>9s}{'+computed':>10s}"
          f"{'all rows':>9s}{'cmpds':>7s}  quantities")
    for i, s in enumerate(order[:TOP_N + 4], 1):
        prot = [r for r in out if r["source"] == s
                and is_named_protein(r["protein"]) and is_measured(r)]
        meas = collections.Counter(r["measure"] for r in prot)
        mark = ("  <- kept (tied)" if s in tied_in else
                "  <- kept" if s in kept else "")
        print(f"  {i:>2d}  {s:22s}{n[s]:>9d}{with_comp[s]:>10d}{tot[s]:>9d}"
              f"{cmpds[s]:>7d}  "
              + ", ".join(f"{k} {v}" for k, v in meas.most_common(3)) + mark)
    cut = n[order[TOP_N - 1]]
    print()
    print(f"  cut at {cut} measured rows; "
          + (f"{len(tied_in)} further source(s) tied there and kept"
             if tied_in else "nothing tied at the cut this time"))
    print()

    cond = {c["id"]: c for c in B.conditions(out)}
    sel = []
    for r in out:
        if r["source"] not in kept:
            continue
        c = cond[r["id"]]
        sel.append({f: (r.get(f) if f in r else c.get(f, "")) for f in FIELDS})
    sel.sort(key=lambda r: (kept.index(r["source"]), r["protein"], r["pfas"]))

    by = collections.Counter(r["source"] for r in sel)
    print(f"{len(sel)} rows from {len(kept)} sources:")
    for s in kept:
        print(f"  {by[s]:5d}  {s}")
    print()
    for f in ("pfas_conc_uM", "protein_conc_uM", "pfas_conc_as_stated",
              "protein_conc_as_stated", "kd_uM", "logK"):
        got = sum(1 for r in sel if str(r[f]).strip())
        print(f"  {f:24s} populated {got:4d} / {len(sel)}")
    write(sel, kept, n, tot)
    print(f"\nwrote {os.path.abspath(OUT)}")
    return 0


NOTES = [
    ("What this file is", "The leading sources that report PFAS binding to "
     "a NAMED SINGLE PROTEIN, ranked by how many MEASURED such values each "
     "supplies. Whole plasma, serum, tissue homogenate, culture medium, "
     "cells and liposomes do not count towards the ranking, nor do docking "
     "scores or predicted values."),
    ("How many", "The top TOP_N sources (currently 10), plus any source "
     "tied with the last of them. A rank boundary that splits two sources "
     "with identical counts is an arbitrary decision dressed as a ranking. "
     "At 10 nothing ties at the boundary, so the rule does not bind here."),
    ("Ties inside the list", "Broken on distinct compounds. Starnes 2024 "
     "and Qin 2024 (plasma) both supply 32 rows; Qin covers 16 compounds "
     "against 8 and so is listed first. Both are kept either way."),
    ("Why not Ryu 2024", "It is the largest contributor to the full "
     "database (320 rows) and reports unbound fractions in whole tissue "
     "homogenate -- total protein plus phospholipid, no protein resolved. "
     "A row count is not a ranking of usefulness for this question."),
    ("Row selection", "By source, not by matrix. Once a source is in, all "
     "of its rows come with it, mixtures and computed values included; the "
     "`protein` column says what each row was measured in and `measure` "
     "marks anything computed."),
    ("Treat with care", "Wu 2024's affinities sit three to four orders "
     "below every other measurement of the same compound on albumin (its "
     "PFOA value is 103 M-1), and half its rows are docking scores. "
     "Starnes 2025 measures by differential scanning fluorimetry at the "
     "protein's melting transition, about three orders off the dialysis "
     "values -- a documented offset, not a disagreement to average over. "
     "Muensterman 2026 is biomimetic chromatography calibrated against "
     "pharmaceuticals, a surrogate rather than an equilibrium measurement, "
     "and the paper does not state its basis."),
    ("", ""),
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


def write(sel, kept, n, tot):
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
    for s in kept:
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
    ws["A1"] = "binding_major_sources.xlsx -- what is in it and why"
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
