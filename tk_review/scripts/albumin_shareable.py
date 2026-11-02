#!/usr/bin/env python3
"""Albumin partition coefficients assembled for comparison against Qin 2024.

WHAT CHANGED AND WHY. An earlier version of this file was a handout: 92
human rows on one comparable scale, everything questionable removed. The
purpose is now different. Qin 2024 is being used in the study and the
question is what the rest of the albumin literature says beside it, so
breadth beats tidiness and a wide spread is the point rather than a
defect. Bovine is in, because Qin is bovine. The fluorescence work is in,
converted with each paper's own fitted exponent rather than refused. The
separative-method block is in at an assumed n = 1. Every row says which
of those it is, so a later cut can be made on any of them.

WHAT THIS MEANS FOR THE NUMBERS. PFOA on human albumin now appears at
1.64, 2.12, 2.67, 3.77, 4.48 and 4.70. That thousand-fold spread is
real and documented -- two laboratories measuring the same albumin Kd by
the same method differ by 1.45 log, and Gao 2019 puts the published
range for that constant at two orders. The spread is the finding, not an
error to be averaged away. Do not take a mean down a column.

SEVEN ROUTES, and `route` names the one each row took:

  measured directly              the paper measured the coefficient
  measured, Qin correction       Qin's two terms summed, the specific one
                                 lifted log10(2) off half-saturation
  converted, n by ITC            stoichiometry counted by calorimetry
  converted, n = 1 (paper)       a single-site model fixes n by
                                 construction; split on whether the
                                 fitted Kd sits inside the tested range
  converted, n = 1 (reader)      the paper reports no n at all
  converted, quenching exponent  the paper's own double-logarithmic,
                                 Scatchard or Stern-Volmer slope, used as
                                 n. It is a fitted exponent, not a
                                 counted site, and the main database
                                 refuses it; it is used here because the
                                 comparison was asked for.

STILL OUT, and cheap to add on request: differential scanning
fluorimetry (82 rows over human, bovine, rat and porcine -- its
constants sit about three orders off dialysis), microscale
thermophoresis (6), and Muensterman 2026's immobilised-HSA
chromatography (82 rows, 83 compounds, on a basis the paper never
states). Rat and porcine rows are dropped with the species filter.

Writes db/albumin_logK_shareable.xlsx.

Run:  python3 -I tk_review/scripts/albumin_shareable.py [species ...]
"""
import collections
import csv
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_kd_to_k.csv")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")

SPECIES = ("human", "bovine", "human; bovine")
ANCHOR = "Qin 2024 (plasma)"

FIELDS = [("pfas_canonical", "PFAS"), ("species", "Species"),
          ("kd_uM", "Kd (uM)"), ("n_used", "n used"),
          ("logK", "log K (L water / kg protein)"),
          ("route", "Route"), ("source", "Source")]

MEASURED = "measured directly"
QIN = "measured; Qin Table S6 equation at"
ITC = "converted; n counted by ITC"
PAPER_OK = "converted; n = 1 fixed by the paper's single-site model, Kd inside the range tested"
PAPER_EX = "converted; n = 1 fixed by the paper's single-site model, Kd beyond the range tested"
READER = "converted; n = 1 assumed by the reader, the paper reports none"
QUENCH = "converted; n is the paper's fitted quenching-plot exponent, not a counted site"
ORDER = [MEASURED, QIN, ITC, PAPER_OK, PAPER_EX, READER, QUENCH]


def rank(route):
    """Sort key. A route may carry its own suffix, e.g. the Cw it was
    evaluated at, so match on the stem rather than the whole string."""
    for i, stem in enumerate(ORDER):
        if route.startswith(stem):
            return i
    return len(ORDER)

# Qin W, Henneberger L, Gluge J, Konig M, Escher BI. "Baseline Toxicity
# Model to Identify the Specific and Nonspecific Effects of PFAS in
# Cell-Based Bioassays." Environ Sci Technol 2024;58:5727-5738,
# doi 10.1021/acs.est.3c09950. Main-text Table 1, log D_BSA/w column.
# Its footnote b says each value is the Freundlich equation of SI Table
# S6 evaluated at IC10,free; that reproduces to within 0.008 log for all
# eleven, and those equations are identical to the plasma paper's Table
# 1 - slope, intercept and R2 to the last digit.
# (compound, IC10,free in mol/L, log D_BSA/w per L of protein)
OPERATING = (
    "10 uM nominal, 10% FBS medium"
)
# log D_BSA/w per kg of protein, and the free concentration it was
# evaluated at, from Qin's SI Table S6 equations solved against his eq 6
# mass balance for the generic 10% FBS medium (Table S7: 3.0 mL/L protein,
# 0.07 mL/L lipid) plus the HEK293H cell term (Table S8). The same
# calculation run on Qin's own 2% csFBS composition reproduces his
# measured IC10,free with a median of +0.03 log and |median| 0.22 log,
# worst -0.67 (PFOS) -- so read these as good to a few tenths, not better.
# (compound, free Cw in uM at 10 uM nominal, log D_BSA/w per kg)
QIN_AT_OP = [
    ("PFBA", 3.79, 2.60), ("PFHxA", 1.134, 3.28),
    ("PFHpA", 0.2344, 4.01), ("PFOA", 0.05577, 4.64),
    ("PFNA", 0.06658, 4.56), ("PFUnA", 0.05805, 4.62),
    ("HFPO-DA", 1.141, 3.28), ("PFHxS", 0.03635, 4.82),
    ("PFOS", 0.01085, 5.35), ("6:2 FTSA", 0.3965, 3.76),
    ("PFOSA", 0.1278, 4.24),
]
_UNUSED_BASELINE_T1 = [
    ("PFBA", 1.78e-3, 1.94), ("PFHxA", 4.72e-4, 2.70),
    ("PFHpA", 6.76e-5, 3.40), ("PFOA", 1.84e-5, 3.98),
    ("PFNA", 7.72e-6, 4.39), ("PFUnA", 3.83e-7, 4.75),
    ("HFPO-DA", 2.20e-4, 2.28), ("PFHxS", 4.22e-5, 3.52),
    ("PFOS", 3.89e-6, 4.52), ("6:2 FTSA", 4.25e-5, 3.71),
    ("PFOSA", 3.28e-7, 4.33),
]

FAMILIES = [
    (("scanning fluorimetry", "(dsf)"), "thermal shift (DSF)"),
    (("thermophoresis",), "thermophoresis (MST)"),
    (("quench", "fluorescence", "fluorimetr"), "fluorescence"),
    (("dialysis", "ultrafiltration", "size exclusion", "filtration",
      "desalting", "chromatograph", "hplc"), "separative"),
    (("isothermal titration",), "calorimetry (ITC)"),
]


def family(method):
    m = (method or "").lower()
    for words, label in FAMILIES:
        if any(w in m for w in words):
            return label
    return "other"


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
                                    "serum)", "immobilis", "immobiliz"))


def row(r, route, logK, n_used="", at_cw=""):
    return {"at_cw_uM": at_cw,
            "pfas_canonical": r.get("pfas_canonical") or r["pfas"],
            "species": r["species"], "protein": r["protein"],
            "method_short": family(r.get("method") or r.get("_method")),
            "kd_uM": r.get("kd_uM", ""), "n_used": n_used,
            "logK": logK, "route": route, "source": r["source"]}


def collect():
    B = _load("baf", "build_binding_affinity_table.py")
    L = _load("lkc", "logk_converted_vs_measured.py")
    N = _load("pfas_names", "pfas_names.py")
    recs = B.build(filter_compounds=False)[1]
    out = []

    # 1. measured, and 2. the Qin anchor with its correction applied
    for r in recs:
        if r["measure"] != "log D" or not r["logK"]:
            continue
        if not is_serum_albumin(r["protein"]):
            continue
        if r["source"] == ANCHOR:
            continue             # replaced by the operating-point row
        rr = dict(r, pfas_canonical=N.canon(r["pfas"]), method=r["_method"])
        out.append(row(rr, MEASURED, float(r["logK"])))
    for name, cw, logk in QIN_AT_OP:
        rr = {"pfas": name, "pfas_canonical": N.canon(name),
              "species": "bovine",
              "protein": "bovine serum albumin (BSA)",
              "method": "BioSPME 96-pin, Freundlich isotherm",
              "kd_uM": "", "source": "Qin 2024"}
        out.append(row(rr, f"{QIN} Cw = {cw:.3g} uM ({OPERATING})",
                       logk))

    # 6. the quenching conversions, using each paper's own exponent
    quenched = set()
    for r in L.revived(recs):
        if not is_serum_albumin(r["protein"]):
            continue
        n = r["logK_basis"].split("n = ")[1].split(",")[0].strip()
        quenched.add((r["source"], r["pfas"], r["kd_uM"]))
        rr = dict(r, pfas_canonical=N.canon(r["pfas"]), method=r["_method"])
        out.append(row(rr, QUENCH, float(r["logK"]), n))

    # 3-5. everything the flat Kd-to-K table already classified
    for r in csv.DictReader(open(SRC, encoding="utf-8")):
        if r["convertible"] == "yes":
            if r["n_basis"] not in ("assumed", ""):
                route = ITC
            else:
                try:
                    inside = float(r["kd_over_top_conc"]) <= 1.0
                except (TypeError, ValueError):
                    inside = False
                route = PAPER_OK if inside else PAPER_EX
            out.append(row(r, route, float(r["logK"]), r["n_sites_used"]))
        elif r["convertible"] == "assume n = 1":
            if family(r["method"]) not in ("fluorescence", "separative"):
                continue               # DSF, MST and the rest stay out
            if (r["source"], r["pfas"], r["kd_uM"]) in quenched:
                continue               # already in, with its own exponent
            out.append(row(r, READER, float(r["logK_if_n_assumed_1"]), "1"))

    out = [r for r in out if r["species"] in SPECIES]
    out.sort(key=lambda r: (rank(r["route"]), r["pfas_canonical"],
                            r["source"]))
    return out


NOTES = [
    ("What this is",
     "Protein-water partition coefficients for PFAS on serum albumin, "
     "human and bovine. log K is log10 of litres of water per kg of "
     "protein."),
    ("Qin 2024 is at ONE operating point",
     "Its rows are the SI Table S6 Freundlich equations evaluated at 10 "
     "uM nominal in a 10% FBS medium. D_BSA/w is not a constant - it "
     "falls roughly 0.3 log per decade of free concentration - so each "
     "Qin row names the free Cw it was evaluated at, in its Route cell. "
     "At another nominal dose or another medium these values change; "
     "re-evaluate the equations rather than reusing the numbers."),
    ("How good is that Cw",
     "Free concentrations come from Qin's own mass balance (his eq 6) "
     "with the generic 10% FBS composition of his Table S7 and the "
     "HEK293H cell term of Table S8. Run on his own 2% csFBS medium the "
     "same calculation reproduces his measured IC10,free with a median "
     "of +0.03 log and a typical error of 0.22 log, worst 0.67. Treat "
     "the Qin rows as good to a few tenths of a log, not better."),
    ("Read the Route column",
     "Rows reach their number by seven routes, from a direct "
     "measurement to a conversion resting on an assumed stoichiometry. "
     "The sheet is sorted best-route first. Two values in the same "
     "column are not necessarily the same kind of number."),
    ("The spread is real",
     "Two laboratories measuring the same albumin Kd by the same method "
     "differ by 1.45 log, and Gao et al. 2019 puts the published range "
     "for that constant at two orders of magnitude. Do not take a mean "
     "down a column."),
    ("The conversion, where one was needed",
     "log K = log10(1e9 / M) + log n - log Kd(uM), M = 66,463 g/mol. "
     "Exact algebra from the Langmuir isotherm. log K moves one for one "
     "with log n, so every converted row turns on the stoichiometry its "
     "route names. Quenching rows use the paper's own fitted exponent, "
     "which is not a counted site."),
]


def write(rows):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    wb.remove(wb.active)
    blue = PatternFill("solid", fgColor="2A78D6")

    ws = wb.create_sheet("Notes")
    ws["A1"] = "PFAS - serum albumin partition coefficients, for comparison"
    ws["A1"].font = Font(bold=True, size=14)
    for i, (k, v) in enumerate(NOTES, start=3):
        ws.cell(row=i, column=1, value=k).font = Font(bold=True)
        ws.cell(row=i, column=2, value=v).alignment = Alignment(
            wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 28
    ws.column_dimensions["B"].width = 96

    def table(title, data, fields, widths=None):
        sh = wb.create_sheet(title)
        sh.append([lab for _, lab in fields])
        for r in data:
            line = []
            for key, _ in fields:
                v = r.get(key, "")
                if key in ("kd_uM", "logK", "n_used") and str(v).strip():
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                line.append(v)
            sh.append(line)
        for c, (key, lab) in enumerate(fields, start=1):
            cell = sh.cell(row=1, column=c)
            cell.fill = blue
            cell.font = Font(bold=True, color="FFFFFF")
            cell.alignment = Alignment(vertical="center", wrap_text=True)
            w = max([len(lab)] + [len(str(r.get(key, ""))) for r in data])
            sh.column_dimensions[get_column_letter(c)].width = \
                min(max(w + 2, 11), 58)
        sh.freeze_panes = "A2"
        sh.auto_filter.ref = sh.dimensions
        return sh

    table("Coefficients", rows, FIELDS)

    anchor = {r["pfas_canonical"] for r in rows
              if r["route"].startswith(QIN)}
    beside = [r for r in rows if r["pfas_canonical"] in anchor]
    beside.sort(key=lambda r: (r["pfas_canonical"],
                               not r["route"].startswith(QIN),
                               rank(r["route"]), r["source"]))
    table("vs Qin 2024", beside, FIELDS)
    wb.save(OUT)


def main():
    rows = collect()
    if len(sys.argv) > 1:
        want = {a.lower() for a in sys.argv[1:]}
        rows = [r for r in rows if r["species"].lower() in want]
    write(rows)

    print(f"{len(rows)} rows, {len(FIELDS)} columns -> "
          f"{os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print(f"  {len({r['pfas_canonical'] for r in rows})} compounds, "
          f"{len({r['source'] for r in rows})} sources, species "
          f"{dict(collections.Counter(r['species'] for r in rows))}")
    print()
    for route in ORDER:
        sub = [r for r in rows if r["route"].startswith(route)]
        if not sub:
            continue
        print(f"  {len(sub):4d}  {route[:58]:60s} "
              f"{', '.join(sorted({r['source'] for r in sub}))[:58]}")
    anchor = sorted({r["pfas_canonical"] for r in rows
                     if r["route"].startswith(QIN)})
    print()
    print(f"  'vs Qin 2024' sheet: {len(anchor)} compounds Qin covers "
          f"({', '.join(anchor)})")
    for c in anchor:
        v = sorted(float(r["logK"]) for r in rows
                   if r["pfas_canonical"] == c)
        q = [float(r["logK"]) for r in rows
             if r["pfas_canonical"] == c
             and r["route"].startswith(QIN)][0]
        print(f"     {c:14s} Qin {q:5.2f} | {len(v)} values, "
              f"{min(v):.2f} to {max(v):.2f}  (spread {max(v)-min(v):.2f} log)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
