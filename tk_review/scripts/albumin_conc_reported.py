#!/usr/bin/env python3
"""Assay concentrations for the albumin extract, AS REPORTED only.

The companion script albumin_conc_vs_logk.py adds three DERIVED columns:
a PFAS concentration in uM (the top of a stated range, unit-converted), a
protein concentration in uM (a stated mass concentration divided by
66,463 g/mol) and a ligand:protein ratio (mostly computed as one over the
other).  Every one of those is arithmetic this project did, not a number
a paper printed.

This script carries the papers' own strings instead.  Nothing is
converted, nothing is reduced to a single number, and a ratio appears
only where the paper states one rather than where it can be computed.
Where a paper reports nothing, the cell is empty -- which is itself the
finding for those rows.

Run:  python3 scripts/albumin_conc_reported.py [workbook.xlsx]
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_reported_conc.xlsx")

NEW = ["PFAS conc (as reported)", "Protein conc (as reported)",
       "ligand:protein ratio (as reported)"]

# A ratio counts as REPORTED only when the paper states one.  Allendorf
# writes "molar ratio held below 0.1" and Fischer "PFAS:protein molar
# ratio <=0.004"; both carry mol/mol in the units field.  Bischel 2011's
# 0.0034, Peng's 18.00 and Wu's 1.00 are this project dividing one stated
# concentration by the other, so they stay blank here.
RATIO_MARKS = ("mol/mol", "molar ratio", ":1 ", "ratio")


NOTES = [
    ("The three concentration columns are verbatim",
     "Each cell is what the paper itself prints, with its own units and "
     "its own wording. Nothing is converted to a common unit, no range is "
     "reduced to a single number, and no concentration is inferred from "
     "another. Where a value is written in scientific notation in the "
     "source database it is rendered plainly (3e0 g/L -> 3 g/L); that "
     "changes the formatting, not the number."),
    ("An empty cell means the paper reports nothing",
     "It does not mean zero and it does not mean the assay had no "
     "concentration. Twenty-two of 261 rows have no PFAS concentration "
     "and twenty-four have no protein concentration, because those "
     "papers never state one."),
    ("The ratio column is blank wherever a ratio would have to be computed",
     "Allendorf and Fischer state a ligand:protein ratio as a design "
     "constraint, and Bischel 2010 states one in words; those appear. "
     "Bischel 2011, Chen, Peng, Starnes and Wu state two concentrations "
     "and no ratio, so dividing one by the other would be this project's "
     "arithmetic, not theirs, and the cell stays empty."),
    ("What this file does NOT change",
     "The log K column still contains converted values wherever the "
     "Route column says so, and the Qin 2024 rows are still quoted at a "
     "free concentration this project computed for 10 uM nominal in 10% "
     "FBS, which Qin does not report. Only the concentration columns "
     "were restricted to reported numbers."),
]


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


# The database stores some transcribed values in scientific notation
# ("3e0 g/L"), which is this project's formatting, not the paper's -- the
# paper prints "3 g/L".  Rendering the same number plainly changes no
# value.  Only a field that is one number followed by a unit is touched;
# a range, a list or a sentence is left exactly as it stands.
_PLAIN = re.compile(r"^(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(\s+\S.*)$")


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


def stated_ratio(lig, prot):
    """The paper's own ratio statement, or '' if it never states one."""
    for text in (lig, prot):
        low = text.lower()
        if any(m in low for m in RATIO_MARKS):
            return text
    return ""


def conditions():
    """(compound, species, source) -> the three as-reported strings."""
    B = _load("baf", "build_binding_affinity_table.py")
    N = _load("pfas_names", "pfas_names.py")
    recs = B.build(filter_compounds=False)[1]
    blocks = B.condition_index(recs)
    B.stamp_conditions(recs, blocks)
    out = {}
    for r in recs:
        if not is_serum_albumin(r["protein"]):
            continue
        if r["measure"] not in ("Kd", "Ka", "log D"):
            continue
        c = blocks[B.block_key(B.condition_row(r))]
        k = (N.canon(r["pfas"]), r["species"], r["source"])
        lig = tidy((c["pfas_conc_as_stated"] or "").strip())
        prot = tidy((c["protein_conc_as_stated"] or "").strip())
        cand = (lig, prot, stated_ratio(lig, prot))
        # Prefer the setting that actually says something.
        if k not in out or (not any(out[k]) and any(cand)):
            out[k] = cand
    return out


def main():
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    src = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    wb = load_workbook(src)
    ws = wb[[s for s in wb.sheetnames if s.lower() != "notes"][0]]
    hdr = [c.value for c in ws[1]]
    i_pfas, i_sp, i_src = (hdr.index("PFAS"), hdr.index("Species"),
                           hdr.index("Source"))

    cond = conditions()
    base = len(hdr)
    for j, name in enumerate(NEW):
        ws.cell(row=1, column=base + j + 1, value=name)

    filled = [0] * 3
    miss = set()
    for row in ws.iter_rows(min_row=2):
        pfas = row[i_pfas].value
        sp = (row[i_sp].value or "").strip()
        source = (row[i_src].value or "").split(" Cw =")[0].strip()
        key = (pfas, sp, source)
        got = cond.get(key)
        if got is None:
            # The database qualifies some sources that the extract names
            # plainly ("Qin 2024 (plasma)" against "Qin 2024"). Match on
            # the stem, but only when it resolves to exactly one entry.
            hits = {v for (p, s, src), v in cond.items()
                    if p == pfas and s == sp and src.startswith(source)}
            got = hits.pop() if len(hits) == 1 else None
        if got is None:
            miss.add(key)
            got = ("", "", "")
        for j, v in enumerate(got):
            if v:
                filled[j] += 1
            ws.cell(row=row[0].row, column=base + j + 1, value=v or None)

    bold = Font(bold=True)
    wrap = Alignment(vertical="top", wrap_text=True)
    for j in range(len(hdr) + len(NEW)):
        col = get_column_letter(j + 1)
        ws.cell(row=1, column=j + 1).font = bold
        width = 46 if j >= base else 18
        ws.column_dimensions[col].width = width
        if j >= base:
            for row in ws.iter_rows(min_row=2):
                row[j].alignment = wrap
    ws.freeze_panes = "A2"

    if "Notes" in wb.sheetnames:
        nt = wb["Notes"]
        for head, body in NOTES:
            nt.cell(row=nt.max_row + 1, column=1, value=head).font = bold
            nt.cell(row=nt.max_row, column=2, value=body).alignment = wrap

    wb.save(OUT)
    n = ws.max_row - 1
    print(f"wrote {OUT}")
    print(f"{n} rows\n")
    for name, k in zip(NEW, filled):
        print(f"  {name:<38} {k:>4} / {n} rows reported "
              f"({100 * k / n:.0f}%)")
    if miss:
        print(f"\n{len(miss)} (compound, species, source) keys unmatched:")
        for k in sorted(miss)[:12]:
            print("   ", k)


if __name__ == "__main__":
    main()
