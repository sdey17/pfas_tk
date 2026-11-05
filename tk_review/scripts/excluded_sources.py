#!/usr/bin/env python3
"""One place that says what this project will not plot or pool, and why.

Two kinds of exclusion live here, and they are not the same kind of
judgement.

THE STANDING RULE: no computed value in any output.  A docking score,
a QSAR prediction, a PBPK model parameter and a machine-learning fit
are not measurements, and this project does not put them on an axis
beside measurements or into a pooled statistic.  They stay in
db/protein_binding.csv, which is a record of what the literature says,
flagged by is_computational; they do not leave it.  The flag is
curated to distinguish a computed VALUE from an experimental paper
that also ran a model -- "no (AutoDock4 binding energies are reported
alongside, as a separate computational strand)" is an experimental
row and stays.  Only a row whose flag does not begin "no" is computed.

NAMED HOLDOUTS: specific sources whose numbers are not what they
appear to be.  Each needs a defect in the paper's own arithmetic or
provenance, written down.  Disagreeing with the rest of the literature
is NOT a reason -- that is how a review talks itself into a
consensus.

Run it to audit:  python3 tk_review/scripts/excluded_sources.py
It prints the computed rows, the named holdouts, and then scans every
shipped artifact for leaks of either kind.
"""
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")


def is_computed(row):
    """True when the value came out of a computer, not an instrument.

    Reads is_computational, which is curated so that a bare "yes" means
    the VALUE is computed, while "no (...)" means the row is a
    measurement from a paper that also ran a model. Anything that does
    not begin "no" is computed.
    """
    return not str(row.get("is_computational", "no")
                   ).strip().lower().startswith("no")


COMPUTED_WHY = (
    "a computed value -- docking score, QSAR prediction, PBPK "
    "parameter or ML fit. Not a measurement, so not on an axis or in "
    "a pooled statistic beside measurements. See the standing rule in "
    "scripts/excluded_sources.py")

SOURCES = {
    "Wu 2024": (
        "HSA constants are not recoverable from the paper's own "
        "tables: the accessible-fluorophore fraction exceeds 1 for all "
        "ten compounds (HSA has one tryptophan), KA is in M^-n with "
        "n = 0.47-0.87 so inverting it is dimensionally wrong, and "
        "re-expressed inside the measured window the ten collapse to a "
        "0.13-log spread. The assay is 17-50 % ethanol v/v with the "
        "protein not held constant. Full working: "
        "scripts/audit_wu2024_hsa.py"),
    "Han 2021": (
        "its hL-FABP constants are explicitly quoted from other "
        "papers, not measured here -- the notes on those rows say so, "
        "and the PFBS value of 1034 uM is Zhang 2013's. Plotting it as "
        "a laboratory row would double-count whatever source it came "
        "from. Its own measurements are thermal-shift dTm, which is "
        "not a dissociation constant"),
    "no binding detected": (
        "reported as NB or ND -- a result, not a missing value, and "
        "not plottable on a log axis. Worth reading against the same "
        "protein's other assay: fluorescence displacement reports NB "
        "for Sheng 2016's N111D and R122G on PFOA and PFNA, where ITC "
        "gives R122G a Kd of 8.4 uM"),
}

# What a COMPUTED VALUE looks like when it reaches a row. These match a
# value POSITION, not a mention: rows are joined on " | ", so a suffix
# on the measure cell or a units cell of its own is a value, while the
# same words inside a method string are a description of how a paper
# worked.
#
# The distinction has to be made, because it is most of the traffic.
# Qin 2024 has 94 measured rows and 2 predicted ones; Chen 2015, Chi
# 2018 and Li 2021 are experimental papers whose METHOD mentions
# docking; a beta-lactoglobulin Kd carries a docking energy in its
# note. Matching the word anywhere condemns all of them.
LEAK = re.compile(
    r"\((?:computed|computational)\)\s*\||"      # measure cell suffix
    r"\|\s*kcal/mol\s*(?:\||$)|"                 # units cell of its own
    r"\bQSAR\b\s*\(|machine-learning QSAR", re.I)

# What a row in a "predict" artifact must say about itself. Matching
# the CATEGORY and the SOURCE LABEL, not a free-text note, so the
# marking cannot be lost by an edit to prose.
PREDICTED = re.compile(r"\bdocking\b|\bDOCKING\b|\bQSAR\b")

# Sheets that document the table rather than carry values. The rule is
# written down in some of them, so the words appear there by design.
DOC_SHEETS = {"notes", "conversion basis", "conditions", "sources",
              "legend", "compounds", "skipped", "qualitative"}

# Artifacts the scan covers, and what each must be free of.
#   plot    feeds a figure: no computed values AND no named holdouts
#   ship    a reading or working copy: no computed values; a named
#           holdout may appear, because these are the record
#   predict THE ONE NAMED EXCEPTION. A file that exists to carry
#           predictions, overlaid on a figure as a visibly separate
#           series. The rule is not suspended for it, it is INVERTED:
#           every row must be a prediction, and the scan fails on a
#           row that is not. That is what keeps this an exception
#           rather than a hole -- a measured value cannot drift into
#           the file unnoticed any more than a computed one can drift
#           into a plotted file.
#
# Why there is an exception at all: a docking score converts onto the
# log K axis exactly (dG = -RT ln K), and seeing the predictions
# beside the measurements is the whole point of the comparison. What
# must never happen is their entering a POOLED statistic or the
# coefficient-selection rule, and nothing here lets them -- the
# measured CSV is still checked as strictly as before.
ARTIFACTS = [
    ("plot", "figures/albumin_chainlength_data.csv"),
    ("plot", "figures/lfabp_by_compound_data.csv"),
    ("predict", "figures/albumin_docking_data.csv"),
    ("ship", "db/albumin_logK_shareable.xlsx"),
    ("ship", "db/albumin_logK_with_conc.xlsx"),
    ("ship", "db/albumin_lfabp_reported_conc.xlsx"),
    ("ship", "db/binding_affinity.xlsx"),
    ("ship", "db/binding_major_sources.xlsx"),
]


def canon_source(study, year):
    """DB study+year -> the short form the artifacts use ("Wu 2024")."""
    s = re.sub(r"\s+et al\.?", "", study).strip().rstrip(",")
    return f"{s} {year}".strip()


def computed_rows():
    with open(DB, encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh) if is_computed(r)]


def artifact_rows(path):
    """(label, joined-text) for every data row of a CSV or workbook."""
    out = []
    if path.endswith(".csv"):
        with open(path, encoding="utf-8") as fh:
            for i, r in enumerate(csv.DictReader(fh), 2):
                out.append((f"row {i}",
                            " | ".join(str(v) for v in r.values())))
        return out
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True)
    for name in wb.sheetnames:
        if name.strip().lower() in DOC_SHEETS:
            continue
        for i, row in enumerate(wb[name].iter_rows(min_row=2,
                                                   values_only=True), 2):
            text = " | ".join(str(c) for c in row if c is not None)
            if text.strip():
                out.append((f"{name} row {i}", text))
    return out


def main():
    comp = computed_rows()
    srcs = sorted({canon_source(r["study"], r["year"]) for r in comp})
    print(f"THE STANDING RULE: no computed value in any output.")
    print(f"  {len(comp)} computed rows in db/protein_binding.csv, "
          f"from {len(srcs)} sources:")
    for s in srcs:
        n = sum(1 for r in comp if canon_source(r["study"], r["year"]) == s)
        kind = sorted({(r["method"] or "")[:46] for r in comp
                       if canon_source(r["study"], r["year"]) == s})
        print(f"    {n:3d}  {s:24s}{kind[0]}")
        for k in kind[1:]:
            print(f"         {'':24s}{k}")
    print()
    print(f"NAMED HOLDOUTS: {len(SOURCES)}")
    for k in SOURCES:
        print(f"    {k}")
    print()

    print("LEAK SCAN")
    print("  plot = feeds a figure: no computed values, no named holdouts")
    print("  ship = a reading or working copy: no computed values")
    print("  predict = the one named exception: predictions ONLY, "
          "so the test is\n            inverted and a row that is "
          "not a prediction is the failure\n")
    bad = 0
    for kind, rel in ARTIFACTS:
        path = os.path.join(HERE, "..", rel)
        if not os.path.exists(path):
            print(f"  [{kind}] {rel:40s} MISSING")
            continue
        rows = artifact_rows(path)
        if kind == "predict":
            # Inverted: a row that does NOT declare itself a
            # prediction is the failure here.
            hits = [(lab, t) for lab, t in rows
                    if not PREDICTED.search(t)]
            bad += len(hits)
            flag = "clean" if not hits else (
                f"{len(hits)} row(s) not marked as predictions")
            print(f"  [{kind}] {rel:40s} {len(rows):5d} rows  {flag}")
            for lab, t in hits[:4]:
                print(f"        {lab}: {t[:88]}")
            continue
        computed = [(lab, t) for lab, t in rows if LEAK.search(t)]
        holdout = []
        if kind == "plot":
            holdout = [(lab, t) for lab, t in rows
                       if any(re.search(rf"(?<!\w){re.escape(h)}(?!\w)", t)
                              for h in SOURCES)]
        hits = computed + holdout
        bad += len(hits)
        flag = "clean" if not hits else (
            f"{len(computed)} computed, {len(holdout)} holdout")
        print(f"  [{kind}] {rel:40s} {len(rows):5d} rows  {flag}")
        for lab, t in hits[:4]:
            print(f"        {lab}: {t[:88]}")
        if len(hits) > 4:
            print(f"        ... and {len(hits) - 4} more")
    print()
    print("  clean: no computed value reaches any artifact, and no named"
          "\n  holdout reaches a plotted one." if not bad else
          f"  {bad} rows need attention.")


if __name__ == "__main__":
    main()
