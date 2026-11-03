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
        out.append({
            "PFAS": r["chemical"],
            "Protein": r["protein"],
            "Species": r["protein_species"],
            "Kd (uM)": kd,
            "n used": None,
            "log K (L water / kg protein)": None,
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
    ("L-FABP rows carry no log K",
     "The Kd-to-partition bridge in this project uses albumin's molar "
     "mass, 66,463 g/mol. L-FABP is about a fifth of that, so the same "
     "arithmetic applied to an L-FABP Kd would be wrong by roughly 0.7 "
     "log. No L-FABP source here reports a partition coefficient, so the "
     "column is left empty rather than filled with the wrong constant."),
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
