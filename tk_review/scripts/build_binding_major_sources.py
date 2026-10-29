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
       9  Jia 2022                  13          9  Kd -- hL-FABP + HSA, MST
      10  Bischel 2011              11         11  log D -- BSA
    ---- cut at 11 measured rows ----
      11  Wu 2024                   10         10  Ka -- HSA
      12  Khazaee 2021              10          6  Kd -- L-FABP, I-FABP, PPARs
      13  Moro 2022                  8          4  Kd -- HSA

Two things about the bottom of that list are worth stating plainly rather
than leaving in the numbers.  Wu 2024 and Khazaee 2021 both supply 10 rows
and both fall ONE row short of the cut.  They tie with each other, not with
the tenth source, so the tie rule below does not fire -- but a one-row
margin is not a judgement about quality, and Khazaee is the only source in
the database measuring PFAS against the PPAR receptors and against
intestinal FABP.  Khazaee's rank also understates its content: 7 of its 17
reported pairs are non-detections ("ND", no binding up to the highest
concentration tested).  Those carry no number, so they are filed as
qualitative rows and cannot count here, yet a clean non-detection against a
named receptor is a real result.  Read it from the full database, not from
this file.

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

ROWS ARE NOT FILTERED BY MATRIX, ONLY BY QUANTITY.  Selection is by source; once a source is
in, all of its rows come with it, mixtures and computed values included,
because dropping a paper's plasma or serum rows would silently discard
measurements the user did not ask to lose.  The `protein` column says what
each row was measured in and `measure` marks anything computed.

COLUMNS.  The nine requested, plus the ligand-to-protein ratio the assay
reached, the stoichiometry each converted partition coefficient used, and a
`condition_id` key.  The concentrations, temperature and buffer used to
sit on every row; they are properties of the EXPERIMENT rather than of the
individual number, so they moved to a Conditions sheet that holds one row per
assay and the measurement rows carry the key instead.  That is mostly, but not
always, one row per study: 37 of the 50 sources in the full database use a
single setting, while Qin 2024, Bischel 2011 and Gao 2025 deliberately vary
concentration, dilution and ionic strength within the paper, which is the
result rather than untidiness.  A key into a block table keeps both cases
honest where a per-study table would flatten the second.

On the Conditions sheet each concentration still appears twice: a micromolar
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

# The nine columns asked for, then the two derived numbers worth scanning
# down a column and the key into Conditions. This is the reading copy, so it
# stays narrow: `binding_site` came out (16 of 578 rows carried one, and it
# is in db/binding_affinity.xlsx and db/combined/binding.csv), and in its
# place sit:
#
#   ligand_protein_ratio_max -- whether a constant was measured in the
#   specific-site regime or the partitioning one, which is what decides
#   whether two numbers are comparable at all;
#
#   n_used and n_basis -- the stoichiometry each converted logK used, and
#   whether the paper FITTED it or ASSUMED it. The pair is inseparable:
#   an assumed 1 and a fitted 1.10 look like the same kind of number in a
#   column and are not. Sparse
#   here on purpose, and the sparsity is the point: of the ten sources in
#   this file only Chen 2025 reaches a logK through the Kd-and-n route, and
#   it does so at n = 1 for every one of its rows. That single assumption is
#   what the -0.86 log median offset in the converted-against-measured
#   comparison rests on, so a reader of THIS file should be able to see it
#   without opening another one. A blank means the row's logK was read off
#   the paper directly and no stoichiometry entered it -- not that one is
#   missing.
#
# The per-row conversion basis is still left out: it is prose, it is the
# widest text in the database, and binding_affinity.xlsx carries it through
# a basis_id lookup for anyone who needs the full reason.
FIELDS = ["pfas", "pfas_canonical", "species", "protein", "measure",
          "value", "units", "kd_uM", "logK", "n_used", "n_basis",
          "ligand_protein_ratio_max", "condition_id", "source"]

# One name per compound, decided for the whole project in pfas_names.py.
# `pfas` is the source's own spelling and stays untouched; group and
# count on `pfas_canonical`.
_ns = importlib.util.spec_from_file_location(
    "pfas_names", os.path.join(HERE, "pfas_names.py"))
N = importlib.util.module_from_spec(_ns)
_ns.loader.exec_module(N)

# A named single protein. Albumin, globulin and the named binding proteins
# count; a matrix that merely CONTAINS protein does not.
PROTEIN_WORDS = ("albumin", "globulin", "transthyretin", "fabp",
                 "fatty acid-binding", "glycoprotein",
                 # Nuclear receptors are named single proteins too. They
                 # were missing from this list, so Khazaee 2021's ten PPAR
                 # rows did not count towards its rank on first ingestion.
                 "ppar", "peroxisome proliferator")


# A matrix can CONTAIN a named protein, or be named after an assay that
# uses one, without being that protein. "PPARgamma assay medium" and
# "HEK293H (PPARgamma-GeneBLAzer) whole cell" both carry the string "ppar";
# neither is a purified receptor. Matrix wins over the protein word.
# NOT in this list: "stationary phase". Muensterman 2026's albumin is
# immobilised on a column, which changes the METHOD (biomimetic, a
# surrogate) but not the material -- it is still human serum albumin. Its
# IAM rows are excluded by "membrane" instead, which is what they are.
MATRIX_WORDS = ("plasma", "whole cell", "cells", "medium", "homogenate",
                "lysate", "tissue", "liposome", "bilayer", "membrane",
                "milk", "whole human serum")


def is_named_protein(protein):
    """One purified, named protein -- not a matrix that contains proteins."""
    p = protein.lower()
    if any(w in p for w in MATRIX_WORDS):
        return False
    return any(w in p for w in PROTEIN_WORDS)


# This file holds experimental binding constants and partition coefficients
# and nothing else. Two kinds of row are excluded outright, not merely left
# out of the ranking:
#   - anything computed: docking scores, QSAR predictions. Not measurements.
#   - unbound fractions. Experimental, but neither a Kd nor a K, and not
#     convertible to either without the binding-site molarity that almost no
#     paper reports. They remain in db/binding_affinity.xlsx.
MEASURES = {"Kd", "Ka", "log D", "log K_membrane"}


def is_measured(rec):
    """An experimental dissociation constant, affinity or partition ratio."""
    return rec["measure"] in MEASURES


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
    cmpds = {s: len({N.canon(r["pfas"]) for r in out if r["source"] == s
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
    # Numbered over the unfiltered table so a condition_id means the same
    # assay here as in db/binding_affinity.xlsx.
    unfiltered = B.build(filter_compounds=False)[1]
    blocks = B.condition_index(unfiltered)
    B.stamp_conditions(out, blocks)
    # n_used rides with the basis stamp, which also sets a basis_id this
    # file does not print; the index is built over the unfiltered set for
    # the same reason the condition blocks are.
    B.stamp_basis(out, B.basis_index(unfiltered))
    n, with_comp, cmpds, order = rank(out)
    excluded = collections.Counter(
        r["measure"] for r in out
        if is_named_protein(r["protein"]) and not is_measured(r))
    kept, tied_in = select(n, order)
    tot = collections.Counter(r["source"] for r in out)

    print(f"Ranked on MEASURED rows against a named single protein. Top "
          f"{TOP_N}, plus any source tied with the {TOP_N}th.")
    print(f"  {'':4s}{'source':22s}{'Kd/Ka/K':>9s}{'+other':>10s}"
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
    if excluded:
        print("  excluded from both the ranking and the file, as neither an "
              "experimental Kd nor a K:")
        for k, v in excluded.most_common():
            print(f"    {v:4d}  {k}")
        print()
    print(f"  cut at {cut} measured rows; "
          + (f"{len(tied_in)} further source(s) tied there and kept"
             if tied_in else "nothing tied at the cut this time"))
    print()

    sel = []
    for r in out:
        if r["source"] not in kept or not is_measured(r):
            continue
        row = {f: r.get(f, "") for f in FIELDS}
        row["pfas_canonical"] = N.canon(r["pfas"])
        sel.append(row)
    sel.sort(key=lambda r: (kept.index(r["source"]), r["protein"], r["pfas"]))

    by = collections.Counter(r["source"] for r in sel)
    print(f"{len(sel)} rows from {len(kept)} sources:")
    for s in kept:
        print(f"  {by[s]:5d}  {s}")
    print()
    for f in ("kd_uM", "logK", "n_used", "n_basis",
              "ligand_protein_ratio_max", "condition_id"):
        got = sum(1 for r in sel if str(r[f]).strip())
        print(f"  {f:24s} populated {got:4d} / {len(sel)}")
    used = {r["condition_id"] for r in sel}
    cond_rows = [c for c in B.condition_table(blocks)
                 if c["condition_id"] in used]
    print(f"  {'assays behind them':24s}           {len(cond_rows)}")
    for f in ("pfas_conc_uM", "protein_conc_uM", "pfas_conc_as_stated",
              "protein_conc_as_stated"):
        got = sum(1 for c in cond_rows if str(c[f]).strip())
        print(f"    Conditions.{f:22s} {got:3d} / {len(cond_rows)}")
    write(sel, cond_rows, kept, n, tot)
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
     "At 10 nothing ties AT the boundary, so the rule does not bind here "
     "-- but Wu 2024 and Khazaee 2021 tie with each other at 10 rows, one "
     "row below the cut of 11. A one-row margin is not a judgement about "
     "quality."),
    ("Missing by one row", "Khazaee 2021 is the only source in the "
     "database measuring PFAS against PPAR-alpha, -gamma and -delta and "
     "against intestinal FABP, and its rank understates it: 7 of its 17 "
     "reported pairs are non-detections, which carry no number and so are "
     "filed as qualitative rows. Read those from db/protein_binding.csv."),
    ("Ties inside the list", "Broken on distinct compounds. Starnes 2024 "
     "and Qin 2024 (plasma) both supply 32 rows; Qin covers 16 compounds "
     "against 8 and so is listed first. Both are kept either way."),
    ("Why not Ryu 2024", "It is the largest contributor to the full "
     "database (320 rows) and reports unbound fractions in whole tissue "
     "homogenate -- total protein plus phospholipid, no protein resolved. "
     "A row count is not a ranking of usefulness for this question."),
    ("Row selection", "By source, not by matrix: once a source is in, its "
     "rows come with it and the `protein` column says what each was "
     "measured in. But only EXPERIMENTAL Kd, Ka, log D and log K_membrane "
     "rows are here. Docking scores and predicted values are excluded as "
     "not measurements; unbound fractions are excluded as neither a Kd nor "
     "a K, and not convertible to either without a binding-site molarity "
     "almost no paper reports. Both remain in db/binding_affinity.xlsx."),
    ("Treat with care", "Wu 2024's affinities sit three to four orders "
     "below every other measurement of the same compound on albumin (its "
     "PFOA value is 103 M-1); its docking scores are excluded here. "
     "Starnes 2025 measures by differential scanning fluorimetry at the "
     "protein's melting transition, about three orders off the dialysis "
     "values -- a documented offset, not a disagreement to average over. "
     "Muensterman 2026 is biomimetic chromatography calibrated against "
     "pharmaceuticals, a surrogate rather than an equilibrium measurement, "
     "and the paper does not state its basis. Jia 2022 measures by "
     "microscale thermophoresis; on hL-FABP its PFOA Kd is 181x the "
     "equilibrium-dialysis value of Khazaee 2021 and the fluorescence "
     "value is 509x it, so a liver-FABP affinity is only comparable "
     "within one method."),
    ("", ""),
    ("", ""),
    ("pfas", "Chemical as the paper names it, untouched."),
    ("pfas_canonical", "One name per compound across the whole database. "
     "The same PFAS is spelled differently by different groups -- PFUnA "
     "and PFUnDA, 6:2 FtS and 6:2 FTSA, FOSA and PFOSA, PFECHS and "
     "PFEtCHxS -- and counting raw names inflates every compound total. "
     "GROUP AND COUNT ON THIS COLUMN, not on `pfas`. The full mapping, "
     "with the evidence for each merge and the groups deliberately left "
     "apart, is db/pfas_name_map.csv."),
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
    ("n_basis", "WHICH EXPERIMENT produced that n: `ITC` or "
     "`displacement` (a counted site), `assumed` (n = 1 because the paper "
     "fitted a single-site model, which fixes it by construction), or "
     "`quenching fit` (the slope of a double-logarithmic, Scatchard or "
     "modified Stern-Volmer plot, which is a fitted exponent and not a "
     "site count -- those rows are NOT converted to a logK anywhere in "
     "this database). Every n_used in this file reads `assumed`, because "
     "Chen 2025 is its only Kd-and-n source and it fitted a single-site "
     "model. The word still matters here: it says the 57 converted "
     "coefficients below rest on a modelling choice, not a measurement."),
    ("n_used", "The stoichiometry that this row's converted logK actually "
     "used, as a number. Blank means no stoichiometry entered the value -- "
     "the paper reported a partition coefficient directly, or no logK was "
     "reached at all -- NOT that one is missing. Of the ten sources here "
     "only Chen 2025 converts through the Kd-and-n route, and it assumes "
     "n = 1 throughout: log K fixes only the ratio n/Kd, so every one of "
     "those values moves one-for-one with that assumption, and the "
     "-0.86 log median gap between converted and measured coefficients is "
     "measured against it. Read the Chen rows knowing that; the other nine "
     "sources' coefficients carry no such assumption."),
    ("ligand_protein_ratio_max", "The highest ratio this assay reached: "
     "most ligand over least protein. Below about 0.01 a measurement sees "
     "the specific site; above about 1 it sees non-specific partitioning, "
     "and the two differ by orders of magnitude, so two constants measured "
     "at different ratios are not comparable even for the same compound "
     "and protein. It is a property of the assay and also sits on the "
     "Conditions sheet; it is repeated here because it is the column you "
     "filter on, and having to do a lookup to learn which regime a number "
     "came from is a lookup too many."),
    ("condition_id", "Key into the Conditions sheet. The join is MANY "
     "measurements to ONE assay, so use a lookup, not a row-for-row "
     "alignment. The ids are numbered over the whole database, so the same "
     "id means the same assay in db/binding_affinity.xlsx and "
     "db/combined/binding_conditions.csv."),
    ("source", "First author and year. 'Qin 2024 (plasma)' is doi "
     "10.1021/acs.est.3c10824; the other Qin 2024 paper, the baseline-"
     "toxicity one, is a different study and is labelled separately."),
    ("", ""),
    ("Conditions", "One row per ASSAY. These columns used to repeat on "
     "every measurement row; they describe the experiment, not the "
     "number, so they live here once. For 37 of the 50 sources in the full "
     "database the setting -- concentrations, temperature, buffer -- is a "
     "single study constant and the sheet repeats it only where the study "
     "measured more than one protein or method. For 13 it is not: Qin 2024 "
     "runs a concentration series, Bischel 2011 a dilution series, Gao "
     "2025 a NaCl series. That is why the key points at an assay and not "
     "at a study."),
    ("rows / compounds", "How many measurements that assay produced, and "
     "over how many distinct compounds."),
    ("assays_in_source", "How many distinct assays the source reports in "
     "all, counting a change of protein or method as well as a change of "
     "setting. 1 or 2 means the conditions are effectively a study "
     "constant."),
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
    ("", "Moving them to their own sheet does not weaken that check -- the "
     "ratio sits in both places: beside the concentrations it is computed "
     "from, on Conditions, and on every measurement row, where it is "
     "actually used."),
    ("What is NOT here", "The binding site, and the per-row basis for each "
     "converted value -- the latter being prose, and the widest text in "
     "the database. This file is the reading copy and stays narrow. Both "
     "are in db/binding_affinity.xlsx and in db/combined/binding.csv: the "
     "site in its own column, the basis through a `basis_id` key into a "
     "Conversion basis table that says how every kd_uM and logK was "
     "reached, or why either is blank. `n_used` is the one part of that "
     "basis lifted out into this file, because it is a number and because "
     "it is the assumption the converted values turn on."),
]


def write(sel, cond_rows, kept, n, tot):
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

    ws = wb.create_sheet("Conditions")
    ws.append(B.COND_FIELDS)
    for c in cond_rows:
        row = []
        for f in B.COND_FIELDS:
            v = c[f]
            if f in B.NUMERIC and str(v).strip():
                try:
                    v = float(v)
                except ValueError:
                    pass
            row.append(v)
        ws.append(row)
    for c, f in enumerate(B.COND_FIELDS, start=1):
        cell = ws.cell(row=1, column=c)
        cell.fill = PatternFill("solid", fgColor="2A78D6")
        cell.font = Font(bold=True, color="FFFFFF")
        cell.alignment = Alignment(vertical="center")
        width = max([len(f)] + [len(str(r[f])) for r in cond_rows])
        ws.column_dimensions[get_column_letter(c)].width = \
            min(max(width + 2, 10), 44)
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = ws.dimensions

    ws = wb.create_sheet("Sources")
    ws.append(["source", "named-protein rows", "rows in this file",
               "compounds", "quantities"])
    by = collections.Counter(r["source"] for r in sel)
    for s in kept:
        sub = [r for r in sel if r["source"] == s]
        meas = collections.Counter(r["measure"] for r in sub)
        ws.append([s, n[s], by[s], len({N.canon(r["pfas"]) for r in sub}),
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
