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
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_kd_to_k.csv")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")

SPECIES = ("human", "bovine", "human; bovine")
ANCHOR = "Qin 2024 (plasma)"

FIELDS = [("pfas_canonical", "PFAS"), ("species", "Species"),
          ("protein", "Protein"), ("method_short", "Method"),
          ("kd_uM", "Kd (uM)"), ("n_used", "n used"),
          ("logK", "log K (L water / kg protein)"),
          ("route", "Route"), ("source", "Source")]

MEASURED = "measured directly"
QIN = "measured; Qin's two terms summed, specific term lifted off half-saturation"
QIN_RAW = "measured; Qin's specific column AS PUBLISHED, at half-saturation, uncorrected"
ITC = "converted; n counted by ITC"
PAPER_OK = "converted; n = 1 fixed by the paper's single-site model, Kd inside the range tested"
PAPER_EX = "converted; n = 1 fixed by the paper's single-site model, Kd beyond the range tested"
READER = "converted; n = 1 assumed by the reader, the paper reports none"
QUENCH = "converted; n is the paper's fitted quenching-plot exponent, not a counted site"
ORDER = [MEASURED, QIN, QIN_RAW, ITC, PAPER_OK, PAPER_EX, READER, QUENCH]

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


def row(r, route, logK, n_used=""):
    return {"pfas_canonical": r.get("pfas_canonical") or r["pfas"],
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
            continue                       # folded below instead
        rr = dict(r, pfas_canonical=N.canon(r["pfas"]), method=r["_method"])
        out.append(row(rr, MEASURED, float(r["logK"])))
    for r in L.totals(recs):
        if r["source"] == ANCHOR and is_serum_albumin(r["protein"]):
            rr = dict(r, pfas_canonical=N.canon(r["pfas"]),
                      method=r.get("_method", ""), kd_uM="")
            out.append(row(rr, QIN, float(r["logK"])))
    # The same values as Qin printed them, so the correction can be seen
    # rather than taken on trust. It moves 3 of the 7 compounds from
    # inside the measured range to above all of it.
    spec, _nons = L.qin_split(recs)
    for r in recs:
        if r["source"] != ANCHOR or not is_serum_albumin(r["protein"]):
            continue
        if "specific column" not in r["location"] or \
                "nonspecific" in r["location"]:
            continue
        rr = dict(r, pfas_canonical=N.canon(r["pfas"]), method=r["_method"],
                  kd_uM="")
        out.append(row(rr, QIN_RAW, round(float(r["logK"]), 3)))

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
    out.sort(key=lambda r: (ORDER.index(r["route"]), r["pfas_canonical"],
                            r["source"]))
    return out


NOTES = [
    ("What this is",
     "Protein-water partition coefficients for PFAS on serum albumin, "
     "human and bovine, assembled so they can be compared against Qin "
     "2024. log K is log10 of litres of water per kg of protein, at zero "
     "loading."),
    ("Read the Route column first",
     "Rows reach their number by seven different routes, from a direct "
     "measurement to a conversion resting on an assumed stoichiometry. "
     "The column says which, on every row, and the sheet is sorted "
     "best-route first. Two values in the same column are not "
     "necessarily the same kind of number."),
    ("The spread is real",
     "PFOA on human albumin appears here at 1.64, 2.12, 2.67, 3.77, 4.48 "
     "and 4.70. Two laboratories measuring the same albumin Kd by the "
     "same method differ by 1.45 log, and Gao et al. 2019 puts the "
     "published range for that constant at two orders of magnitude. The "
     "spread is the finding. Do not take a mean down a column."),
    ("The conversion",
     "log K = log10(1e9 / M) + log n - log Kd(uM), with M = 66,463 g/mol "
     "for albumin. Exact algebra from the Langmuir isotherm. Because log "
     "K moves one for one with log n, every converted row turns on the "
     "stoichiometry its route names."),
    ("Qin 2024, twice",
     "Bovine, and the reason this file exists. It reports albumin in two "
     "components and quotes the specific one at half-saturation (its Eq. "
     "15, K = C_max / 2 K_d). It appears on TWO rows per compound: as "
     "published, and corrected - the two terms summed in linear space "
     "after the specific one is lifted by log10(2). Decide which you "
     "want rather than inheriting the choice."),
    ("What the correction does",
     "It is worth 0.30 to 0.34 log and it changes where Qin sits. As "
     "published, Qin falls inside the range other laboratories measured "
     "for 6 of its 7 compounds. Corrected, it sits ABOVE every measured "
     "value for PFOA, PFOS and PFHxS - the three that matter most - and "
     "inside for GenX, PFBA, PFHpA and PFHxA. PFOS is the one compound "
     "where Qin disagrees with the others either way."),
    ("Quenching rows",
     "The n used is the slope of a double-logarithmic, Scatchard or "
     "modified Stern-Volmer plot. That is a fitted exponent, not a "
     "counted site, and the project's main database refuses to convert "
     "on it. These rows are included because the comparison was asked "
     "for; weigh them accordingly."),
    ("Species",
     "Human and bovine. Albumin itself barely differs - HSA and BSA are "
     "66,437 and 66,398 g/mol by sequence, 0.00025 log apart - so a "
     "human-bovine gap in this file is a laboratory and method "
     "difference, not a species one."),
    ("Still out, available on request",
     "Differential scanning fluorimetry (82 rows; its constants sit "
     "about three orders off dialysis values), microscale thermophoresis "
     "(6 rows), and Muensterman 2026's immobilised-albumin "
     "chromatography (82 rows over 83 compounds, on a basis the paper "
     "never states). Rat and porcine rows are also out."),
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

    anchor = {r["pfas_canonical"] for r in rows if r["route"] == QIN}
    beside = [r for r in rows if r["pfas_canonical"] in anchor]
    beside.sort(key=lambda r: (r["pfas_canonical"], r["route"] != QIN,
                               ORDER.index(r["route"]), r["source"]))
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
        sub = [r for r in rows if r["route"] == route]
        if not sub:
            continue
        print(f"  {len(sub):4d}  {route[:58]:60s} "
              f"{', '.join(sorted({r['source'] for r in sub}))[:58]}")
    anchor = sorted({r["pfas_canonical"] for r in rows if r["route"] == QIN})
    print()
    print(f"  'vs Qin 2024' sheet: {len(anchor)} compounds Qin covers "
          f"({', '.join(anchor)})")
    for c in anchor:
        v = sorted(float(r["logK"]) for r in rows
                   if r["pfas_canonical"] == c)
        q = [float(r["logK"]) for r in rows
             if r["pfas_canonical"] == c and r["route"] == QIN][0]
        print(f"     {c:14s} Qin {q:5.2f} | {len(v)} values, "
              f"{min(v):.2f} to {max(v):.2f}  (spread {max(v)-min(v):.2f} log)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
