#!/usr/bin/env python3
"""Append the L-FABP binding data to the reported-concentration extract.

Takes the output of albumin_conc_reported.py and adds the liver
fatty-acid-binding protein rows under the same rule: every concentration
cell is the paper's own string, nothing converted, nothing computed.

Two columns are added to make the combined file readable:

  Protein  -- the albumin extract assumed one protein throughout. Mixing
              L-FABP in without naming the protein would be wrong.
  Note     -- carries reported outcomes that are not numbers ("no
              binding detected") and provenance the row cannot be read
              without (Han's Kd values are quoted from Zhang 2013).

The log K column stays EMPTY for every L-FABP row. The bridge this
project uses, log K = log10(1e9/M) + log n - log Kd, takes M = 66,463
g/mol, which is albumin. L-FABP is roughly a fifth of that mass, so the
same arithmetic would be wrong by about 0.7 log, and no L-FABP paper
here reports a partition coefficient to use instead.

Run:  python3 scripts/add_lfabp_reported_conc.py [workbook.xlsx]
"""
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
DEFAULT = os.path.join(HERE, "..", "db", "albumin_logK_reported_conc.xlsx")
OUT = os.path.join(HERE, "..", "db", "albumin_lfabp_reported_conc.xlsx")

# Method -> the Route label, matched on a substring of the method string.
ROUTES = [
    ("microscale thermophoresis", "measured Kd (MST)"),
    ("equilibrium dialysis", "measured Kd (equilibrium dialysis)"),
    ("thermal shift", "thermal shift (QITSA) - dTm, not a Kd"),
    ("qitsa", "thermal shift (QITSA) - dTm, not a Kd"),
    ("isothermal titration", "measured Kd (ITC + fluorescence)"),
    ("fluorescence", "fluorescence displacement"),
    ("displacement", "displacement assay"),
    ("binding affinity assay", "binding affinity assay (rank order)"),
]
_PLAIN = re.compile(r"^(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(\s+\S.*)$")
# Word-boundary matching: the substring "ratio" also occurs inside
# "titration", which would falsely tag every MST row as reporting a
# ligand:protein ratio.
RATIO_RE = re.compile(r"mol/mol|\bmolar ratio\b|\bmole ratio\b"
                      r"|\bratios?\b|\d+\s*:\s*1\b", re.I)


# A cell can hold the paper's words and still report no number: "see Text
# S7", "not quantified in molar terms". A cell counts as reporting a
# CONCENTRATION only when it carries a number against a concentration
# unit -- mol/mol is a ratio, not a concentration.
_UNIT = r"(?:[munp]?M|mol/L|g/L|mg/mL|[mun]g/L|[mun]g/mL|mg/L)"
_HAS_NUM = re.compile(r"\d")
_HAS_UNIT = re.compile(rf"(?<![a-zA-Z/]){_UNIT}(?![a-zA-Z/])")


def has_conc(text):
    t = str(text or "")
    return bool(_HAS_NUM.search(t) and _HAS_UNIT.search(t))


# ---------------------------------------------------------------------
# Molar masses, from UniProt, for the Kd -> partition conversion.
# Each is the mature chain of the sequence the papers name, read from
# the UniProt record rather than quoted from memory:
#   P07148  FABPL_HUMAN   127 aa   14,208 Da
#   P02692  FABPL_RAT     127 aa   14,273 Da
#   P12710  FABPL_MOUSE   127 aa   14,246 Da
#   P12104  FABP2_HUMAN   132 aa   15,237 Da  (intestinal)
# Zebrafish L-FABP is deliberately absent: no zebrafish row carries a
# Kd, so no mass is needed and none is guessed.
FABP_MASS = {
    ("human", "liver"): (14208.0, "UniProt P07148"),
    ("rat", "liver"): (14273.0, "UniProt P02692"),
    ("mouse", "liver"): (14246.0, "UniProt P12710"),
    ("human", "intestinal"): (15237.0, "UniProt P12104"),
}


def fabp_mass(species, protein):
    kind = "intestinal" if "intestinal" in protein.lower() else "liver"
    return FABP_MASS.get((species.strip().lower(), kind), (None, ""))


# The same bridge the albumin table uses, written so it also covers a
# multi-site fit:  K = (1e9 / M) * sum(1 / Kd_i),  Kd_i in uM.
# One site with n copies collapses to log10(1e9/M) + log n - log Kd, and
# a list of per-site Kd sums the reciprocals instead of assuming they
# are equal.
def log_k(kds_uM, mass):
    if not kds_uM or not mass:
        return None
    return math.log10(1e9 / mass) + math.log10(sum(1.0 / k for k in kds_uM))


_KD_LIST = re.compile(r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")


def kd_values_uM(raw, units):
    """Every Kd the cell reports, in uM, or [] if it reports none.

    A unit conversion inside a derived column is not the same as putting
    a converted number in a column labelled 'as reported', so M and nM
    are accepted here and scaled.
    """
    u = (units or "").strip().lower()
    if raw.strip().startswith(">") or raw.strip().startswith("<"):
        return []                      # a bound, not a value
    nums = [float(x) for x in _KD_LIST.findall(raw or "")]
    if not nums:
        return []
    if u.startswith("um"):
        return nums
    if u.startswith("nm"):
        return [n / 1e3 for n in nums]
    if u.startswith("m") and not u.startswith("mol/mol"):
        return [n * 1e6 for n in nums]
    return []


def tidy(text):
    if not text:
        return ""
    if text.lower().startswith("not stated"):
        return "not stated"
    m = _PLAIN.match(text)
    if not m:
        return text
    v = float(m.group(1))
    out = f"{v:.10f}".rstrip("0").rstrip(".") if abs(v) < 1e7 else str(v)
    return out + m.group(2)


def as_stated(value, units):
    v, u = (value or "").strip(), (units or "").strip()
    return tidy(f"{v} {u}".strip() if v else "")


def stated_ratio(*texts):
    for t in texts:
        if RATIO_RE.search(t):
            return t
    return ""


def range_index():
    """(Source, PFAS) -> the range check's verdict for that row.

    A Kd is a property of the curve at half occupancy. If the
    titration never went near that concentration, the number is an
    extrapolation, and that belongs in the same row as the number
    rather than in a separate file nobody opens. Produced by
    scripts/kd_within_tested_range.py, which is run first.
    """
    path = os.path.join(HERE, "..", "db", "kd_range_check.csv")
    if not os.path.exists(path):
        return {}
    idx = {}
    for row in csv.DictReader(open(path, encoding="utf-8")):
        src = f"{row['study'].replace(' et al.', '')} {row['year']}"
        bits = [row["range_flag"]]
        if row["fold_outside"]:
            bits[0] += f" ({row['fold_outside']}x beyond)"
        if row["identifiability"]:
            bits.append(row["identifiability"])
        idx.setdefault((src, row["chemical"]), " - ".join(bits))
    return idx


def route_for(method):
    low = (method or "").lower()
    for key, label in ROUTES:
        if key in low:
            return label
    return "reported, method not classified"


def lfabp_rows():
    out = []
    for r in csv.DictReader(open(DB)):
        if "fabp" not in r["protein"].lower():
            continue
        if r["is_computational"].strip().lower() in ("yes", "true", "1"):
            continue
        kd_raw, kd_units = r["kd"].strip(), r["kd_units"].strip()
        # A number in the Kd column ONLY when the paper reports one in uM.
        kd = None
        note = []
        try:
            val = float(kd_raw)
            kd = val if kd_units.lower().startswith("um") else None
        except ValueError:
            val = None
        if kd is None and kd_raw:
            note.append(f"Kd as reported: {kd_raw}"
                        + (f" {kd_units}" if kd_units else ""))
        elif kd is None and kd_units:
            note.append(kd_units)
        if "NOT measured in this paper" in kd_units:
            note.append("value quoted from another study, not measured here")
        lig = as_stated(r["ligand_conc_range"], r["ligand_conc_units"])
        prot = as_stated(r["protein_conc"], r["protein_conc_units"])

        # --- the conversion -------------------------------------------
        mass, mass_src = fabp_mass(r["protein_species"], r["protein"])
        kds = kd_values_uM(kd_raw, kd_units)
        lk = log_k(kds, mass)
        n_used = None
        if lk is not None:
            n_used = len(kds)
            note.append(f"log K from {len(kds)} Kd value"
                        f"{'s' if len(kds) > 1 else ''} "
                        f"({', '.join(f'{k:g}' for k in kds)} uM) at "
                        f"M = {mass:,.0f} g/mol, {mass_src}")
            if len(kds) == 1 and "single-site" not in r["n_binding_sites"]:
                note.append("stoichiometry not reported; n = 1 assumed")
            elif len(kds) > 1:
                note.append("multi-site fit; reciprocal Kd summed over "
                            "sites, not assumed equal")
        elif kd_raw and mass is None:
            note.append("no molar mass resolved for this species")

        out.append({
            "PFAS": r["chemical"],
            "Protein": r["protein"],
            "Species": r["protein_species"],
            "Kd (uM)": kd,
            "n used": n_used,
            "log K (L water / kg protein)":
                round(lk, 2) if lk is not None else None,
            "Route": route_for(r["method"]),
            "Source": f"{r['study'].replace(' et al.', '')} {r['year']}",
            "PFAS conc (as reported)": lig or None,
            "Protein conc (as reported)": prot or None,
            "ligand:protein ratio (as reported)":
                stated_ratio(lig, prot) or None,
            "Note": "; ".join(note) or None,
        })
    return out


NOTES = [
    ("L-FABP log K uses L-FABP's own molar mass",
     "The albumin rows convert at 66,463 g/mol. The L-FABP rows convert "
     "at the mature chain mass of the protein each paper names, read "
     "from UniProt: human P07148 14,208 Da, rat P02692 14,273 Da, mouse "
     "P12710 14,246 Da, human intestinal P12104 15,237 Da. Because "
     "L-FABP is about a fifth of albumin's mass, its bridge constant "
     "log10(1e9/M) is 4.85 against albumin's 4.18 -- an L-FABP Kd and "
     "an albumin Kd of the same value give partition coefficients 0.67 "
     "log apart. No zebrafish mass is included and no zebrafish row "
     "carries a Kd, so none is guessed."),
    ("Only one L-FABP study reports both concentrations",
     "Khazaee 2021 states 1 uM protein, 0.33-153.5 ng/mL ligand and a "
     "ligand:protein mole ratio of 0.05-5. Han 2021 states 100 uM ligand "
     "but says the protein is endogenous in a lysate and 'not quantified "
     "in molar terms'. Jia 2022 points to Text S7 and 'MST titration "
     "series' without numbers. The rest report neither."),
    ("Han 2021's Kd values are not Han's",
     "The database records them as taken from Zhang et al. 2013, not "
     "measured in that paper. Those rows are flagged in the Note column "
     "and should not be counted as an independent measurement."),
]


# The five compounds under active comparison. A dedicated sheet puts
# albumin and L-FABP side by side for these, because the combined sheet
# is 301 rows and the question is a per-compound one.
FOCUS = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]


def focus_sheet(wb, hdr, rows_all):
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    if "Five compounds" in wb.sheetnames:
        del wb["Five compounds"]
    sh = wb.create_sheet("Five compounds")
    cols = ["PFAS", "Protein", "Species", "Kd (uM)", "n used",
            "log K (L water / kg protein)", "Route", "Source",
            "PFAS conc (as reported)", "Protein conc (as reported)",
            "ligand:protein ratio (as reported)",
            "Kd within tested range?", "Note"]
    sh.append(cols)
    idx = {c: hdr.index(c) for c in cols if c in hdr}
    picked = []
    for name in FOCUS:
        for r in rows_all:
            if (r[hdr.index("PFAS")] or "") != name:
                continue
            picked.append([r[idx[c]] if c in idx else None for c in cols])
    for row in sorted(picked, key=lambda x: (FOCUS.index(x[0]),
                                             "FABP" in str(x[1]),
                                             str(x[7]))):
        sh.append(row)

    bold, wrap = Font(bold=True), Alignment(vertical="top", wrap_text=True)
    for j, name in enumerate(cols):
        sh.cell(row=1, column=j + 1).font = bold
        wide = ("as reported" in name
                or name in ("Protein", "Note", "Route", "Source"))
        sh.column_dimensions[get_column_letter(j + 1)].width = (
            46 if "as reported" in name else 34 if wide else 16)
        if wide:
            for row in sh.iter_rows(min_row=2):
                row[j].alignment = wrap
    sh.freeze_panes = "A2"
    return sh


def main():
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    src = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    wb = load_workbook(src)
    ws = wb["Coefficients"]
    hdr = [c.value for c in ws[1]]

    # Protein goes after Species; Note goes last.
    ws.insert_cols(hdr.index("Species") + 2)
    ws.cell(row=1, column=hdr.index("Species") + 2, value="Protein")
    for row in ws.iter_rows(min_row=2):
        row[hdr.index("Species") + 1].value = "serum albumin"
    hdr = [c.value for c in ws[1]]
    ws.cell(row=1, column=len(hdr) + 1, value="Note")
    hdr = [c.value for c in ws[1]]
    ws.cell(row=1, column=len(hdr) + 1, value="concentrations reported")
    hdr = [c.value for c in ws[1]]
    ws.cell(row=1, column=len(hdr) + 1,
            value="Kd within tested range?")
    hdr = [c.value for c in ws[1]]

    n_alb = ws.max_row - 1
    rows = lfabp_rows()
    for rec in rows:
        ws.append([rec.get(c) for c in hdr])

    bold = Font(bold=True)
    wrap = Alignment(vertical="top", wrap_text=True)
    for j, name in enumerate(hdr):
        ws.cell(row=1, column=j + 1).font = bold
        wide = ("as reported" in name or name in ("Protein", "Note",
                                                  "Route", "Source"))
        ws.column_dimensions[get_column_letter(j + 1)].width = (
            46 if "as reported" in name else 34 if wide else 16)
        if wide:
            for row in ws.iter_rows(min_row=2):
                row[j].alignment = wrap
    ws.freeze_panes = "A2"

    # A concentration cell can carry the paper's own words and still
    # report no number ("see Text S7", "not quantified in molar
    # terms", "molar ratio held below 0.1"). This column says, per
    # row, which of the two are actually NUMBERS, so the file can be
    # filtered to rows whose concentrations the paper really gives.
    ci = hdr.index("concentrations reported")
    pi = hdr.index("PFAS conc (as reported)")
    qi = hdr.index("Protein conc (as reported)")
    tally = {}
    for row in ws.iter_rows(min_row=2):
        a, b = has_conc(row[pi].value), has_conc(row[qi].value)
        verdict = ("both" if a and b else "PFAS only" if a
                   else "protein only" if b else "neither")
        row[ci].value = verdict
        tally[verdict] = tally.get(verdict, 0) + 1

    # Where each Kd sits relative to the concentrations actually
    # tested. Blank means the paper does not state enough to check,
    # which is not the same as passing.
    ridx = range_index()
    ri = hdr.index("Kd within tested range?")
    si, fi = hdr.index("Source"), hdr.index("PFAS")
    rtally = {}
    for row in ws.iter_rows(min_row=2):
        src = str(row[si].value or "").split(" Cw")[0].strip()
        v = ridx.get((src, row[fi].value))
        row[ri].value = v
        key = (v or "not checkable").split(" - ")[0].split(" (")[0]
        rtally[key] = rtally.get(key, 0) + 1

    all_rows = list(ws.iter_rows(min_row=2, values_only=True))
    sh = focus_sheet(wb, hdr, all_rows)

    if "Notes" in wb.sheetnames:
        nt = wb["Notes"]
        for head, body in NOTES:
            nt.cell(row=nt.max_row + 1, column=1, value=head).font = bold
            nt.cell(row=nt.max_row, column=2, value=body).alignment = wrap

    wb.save(OUT)
    both = sum(1 for r in rows if has_conc(r["PFAS conc (as reported)"])
               and has_conc(r["Protein conc (as reported)"]))
    kd = sum(1 for r in rows if r["Kd (uM)"] is not None)
    print(f"wrote {OUT}")
    print("  Kd within the tested concentration range:")
    for k in sorted(rtally, key=lambda k: -rtally[k]):
        print(f"    {rtally[k]:5d}  {k}")
    print("  concentrations reported as NUMBERS, over every row:")
    for k in ("both", "PFAS only", "protein only", "neither"):
        if k in tally:
            print(f"     {k:<14}{tally[k]:>4}")
    print(f'  "Five compounds" sheet: {sh.max_row - 1} rows')
    print(f"  {n_alb} albumin rows + {len(rows)} L-FABP rows "
          f"= {ws.max_row - 1}\n")
    print(f"  L-FABP rows with a numeric Kd in uM      {kd:>3} / {len(rows)}")
    print(f"  L-FABP rows reporting BOTH as numbers  {both:>4} / {len(rows)}")
    src_counts = {}
    for r in rows:
        e = src_counts.setdefault(r["Source"], [0, 0, 0])
        e[0] += 1
        e[1] += has_conc(r["PFAS conc (as reported)"])
        e[2] += has_conc(r["Protein conc (as reported)"])
    print(f"\n  {'source':<18}{'rows':>5}{'PFAS conc':>11}{'prot conc':>11}")
    for s, (n, a, b) in sorted(src_counts.items()):
        print(f"  {s:<18}{n:>5}{a:>11}{b:>11}")


if __name__ == "__main__":
    main()
