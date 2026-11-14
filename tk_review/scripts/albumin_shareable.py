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

# The default species set. Rat is in it because the per-compound and
# chain-length figures all facet on human / bovine / rat, and Fan 2025
# and Han 2003 supply real rat serum-albumin measurements. Leaving it
# out meant a bare re-run of this script silently emptied the rat rows
# out of every downstream figure.
SPECIES = ("human", "bovine", "rat", "human; bovine")

# Sentinel for "no species filter at all". It cannot be None, because
# None is the default and means SPECIES.
ALL = "\0all"
# Rows the sheet must NOT also take from the "measured" pass, because the
# anchor module supplies them: both Qin 2024 labels, whichever is active.
# Dropping only the active one would leave the other paper's albumin rows
# in the sheet beside it, which is the cross-contamination this guards.
ANCHOR = ("Qin 2024", "Qin 2024 (plasma)")

FIELDS = [("pfas_canonical", "PFAS"), ("species", "Species"),
          ("kd_uM", "Kd (uM)"), ("n_used", "n used"),
          ("logK", "log K (L water / kg protein)"),
          ("route", "Route"), ("source", "Source")]

MEASURED = "measured directly"
QIN = "measured; Qin Table S6 equation at"
ITC = "converted; n counted by ITC"
FITTED_OK = "converted; n fitted by the paper from its own binding isotherm, Kd inside the range tested"
FITTED_EX = "converted; n fitted by the paper from its own binding isotherm, Kd beyond the range tested"
PAPER_OK = "converted; n = 1 fixed by the paper's single-site model, Kd inside the range tested"
PAPER_EX = "converted; n = 1 fixed by the paper's single-site model, Kd beyond the range tested"
# Sources whose method family sits outside the default two but which
# have been read from the primary and admitted DELIBERATELY, one at a
# time.  The default gate gives the reader's n = 1 route only to
# "fluorescence" and "separative", which silently excluded every DSF,
# MST and NMR row -- 98 of them, all already carrying a Kd.  That was
# a policy acting as a verdict, so admission is now explicit and
# per-source rather than implied by a family list.
#
# Starnes 2024 (DSF): admitted 2026-11-14 after reading Toxicol Sci
# 2024;199(1):132-148.  Its Table 2 is the only four-species albumin
# comparison in this compilation.  Its Kd are MILLIMOLAR and sit a
# median 2.05 log BELOW the per-compound, per-species pools -- every
# one of the twenty comparisons that has a pool is negative.  That is
# a method offset, not scatter: DSF reports the concentration at which
# a ligand measurably stabilises the protein against thermal
# unfolding, which needs far more occupancy than half-saturation, and
# the Kd comes off a 50-90 C ramp rather than one temperature.  The
# rows carry that in their note.  Full working:
# scripts/audit_starnes2024_dsf.py
#
# NOT admitted yet, and each needs its own primary read first:
# Starnes 2025 (49 DSF rows), Alesio 2024 (11 DSF), Jia 2022 (6 MST),
# Jackson 2021 (1 DSF), Fedorenko 2021 and MacManus-Spencer 2010 (NMR).
METHOD_ADMIT = {"Starnes 2024"}


READER = "converted; n = 1 assumed by the reader, the paper reports none"

# Sources where the paper DOES report a stoichiometry but the reader
# still uses n = 1, by decision.  READER's "the paper reports none"
# would be a false statement about these papers, so they carry their
# own wording.  Every replacement must keep the phrase "assumed by the
# reader" verbatim: what_drives_albumin_spread.py,
# make_route_and_loading_figures.py, build_notebook.py and
# fu_physiological.py all bucket this route by that substring, and a
# row that loses it silently leaves the "reader n=1" class.
READER_DECLINED = {
    "Klevens 1954":
        "converted; n = 1 assumed by the reader BY DECISION -- the "
        "paper does report n = 85 (Fig. 3, extrapolated), not used: "
        "n falls from ~150 to ~45 as BSA goes 0.5% to 4% because the "
        "protein aggregates, so it is a surfactant-adsorption count, "
        "not a site count. See fix_klevens1954_from_primary.py",
}


def reader_route(source):
    return READER_DECLINED.get(source, READER)
QUENCH = "converted; n is the paper's fitted quenching-plot exponent, not a counted site"
ORDER = [MEASURED, QIN, ITC, FITTED_OK, FITTED_EX, PAPER_OK,
         PAPER_EX, READER, QUENCH]


def rank(route):
    """Sort key. A route may carry its own suffix, e.g. the Cw it was
    evaluated at, so match on the stem rather than the whole string."""
    for i, stem in enumerate(ORDER):
        if route.startswith(stem):
            return i
    return len(ORDER)

# The Qin anchor lives in qin2024_anchor.py, which also holds the other
# paper's version of it and the audit that tells them apart.  ACTIVE there
# decides which one this sheet, the concentration columns and both figures
# all use.  Nothing about Qin is hard-coded here any more.
ANCHOR_MOD = "qin2024_anchor.py"

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


def collect(keep=None):
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
        if r["source"] in ANCHOR:
            continue             # replaced by the operating-point row
        rr = dict(r, pfas_canonical=N.canon(r["pfas"]), method=r["_method"])
        out.append(row(rr, MEASURED, float(r["logK"])))
    for a in _load("qin_anchor", ANCHOR_MOD).anchor():
        rr = {"pfas": a["pfas"], "pfas_canonical": N.canon(a["pfas"]),
              "species": "bovine",
              "protein": "bovine serum albumin (BSA)",
              "method": "BioSPME 96-pin, Freundlich isotherm",
              "kd_uM": "", "source": a["label"]}
        out.append(row(rr, f"{QIN} {a['where']}", a["logK"]))

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
            # n_basis says HOW the stoichiometry was obtained, and the
            # three cases are not interchangeable. "ITC" is counted by
            # calorimetry. "fitted" is read off the paper's own binding
            # isotherm -- equilibrium dialysis and ultrafiltration work,
            # never calorimetry; it used to fall into the ITC bucket here,
            # which mislabelled 25 rows over five studies. "assumed" (or
            # blank) is a single-site model that fixes n = 1 by
            # construction. The last two also split on whether the fitted
            # Kd lands inside the concentration actually tested.
            try:
                inside = float(r["kd_over_top_conc"]) <= 1.0
            except (TypeError, ValueError):
                inside = False
            if r["n_basis"] == "ITC":
                route = ITC
            elif r["n_basis"] == "fitted":
                route = FITTED_OK if inside else FITTED_EX
            else:
                route = PAPER_OK if inside else PAPER_EX
            out.append(row(r, route, float(r["logK"]), r["n_sites_used"]))
        elif r["convertible"] == "assume n = 1":
            if (family(r["method"]) not in ("fluorescence", "separative")
                    and r["source"] not in METHOD_ADMIT):
                continue               # see METHOD_ADMIT above
            if (r["source"], r["pfas"], r["kd_uM"]) in quenched:
                continue               # already in, with its own exponent
            out.append(row(r, reader_route(r["source"]),
                           float(r["logK_if_n_assumed_1"]), "1"))

    if keep is not ALL:
        out = [r for r in out if r["species"] in (keep or SPECIES)]
    out.sort(key=lambda r: (rank(r["route"]), r["pfas_canonical"],
                            r["source"]))
    return out


NOTES = [
    ("What this is",
     "Protein-water partition coefficients for PFAS on serum albumin, "
     "human and bovine. log K is log10 of litres of water per kg of "
     "protein."),
    ("Which Qin 2024 this is",
     "Two Qin 2024 papers in Environ Sci Technol both print a log D_BSA/w "
     "column holding different numbers. These rows are the BASELINE "
     "TOXICITY paper, 58(13):5727-5738, doi 10.1021/acs.est.3c09950, "
     "Table 1 as printed. The plasma-binding paper 58(23):9954-9966 is in "
     "the database but excluded from this sheet."),
    ("Qin 2024 is at ONE concentration, and it is a high one",
     "Each value is the Freundlich fit read at IC10,free, the free "
     "concentration causing 10% cytotoxicity. D_BSA/w falls about 0.3 log "
     "per decade of free concentration, so the Route cell names the Cw. "
     "At IC10,free the loading is 0.3 to 12 PFAS per albumin -- at or past "
     "the site count, and three orders above the physiological nu < 0.01. "
     "Read these as saturated-carrier values, not trace-limit constants."),
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
    # An argument REPLACES the default species set rather than narrowing
    # it, so a species outside SPECIES (rat, mouse) can be asked for.
    # "all" keeps every species the database holds.
    keep = None
    if len(sys.argv) > 1:
        want = {a.lower() for a in sys.argv[1:]}
        keep = ALL if "all" in want else want
    rows = collect(keep=keep)
    if keep and keep is not ALL:
        rows = [r for r in rows if r["species"].lower() in keep]
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
