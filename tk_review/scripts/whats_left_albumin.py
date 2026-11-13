#!/usr/bin/env python3
"""What albumin data exists in the database but never reaches the sheet.

Run this to answer "are we done?".  It does not judge; it reports what
is held back and by which rule, so the remaining decisions are visible
rather than implicit.

Three blocks dominate, and none of them is a missing paper -- every
one is already ingested and sitting in db/protein_binding.csv.

  THE METHOD-FAMILY GATE.  albumin_shareable.py admits the reader's
  n = 1 route only for method families "fluorescence" and
  "separative".  Differential scanning fluorimetry, microscale
  thermophoresis and 19F NMR are none of those, so ~98 rows that
  carry a Kd and a computed logK_if_n_assumed_1 stop there.  Among
  them are the only PORCINE albumin measurements in the whole
  database, which would add a fourth species to a sheet that has
  human, bovine and rat.

  MUENSTERMAN 2026.  165 rows over 83 compounds on an immobilised-HSA
  column (biomimetic chromatography).  Never converted, because its
  log K_HSA has no stated basis -- the ka_units field says so.  It is
  also a different quantity from a solution binding constant.  For
  scale: the sheet carries 68 distinct compounds in total.

  LI 2021.  26 equilibrium-dialysis rows refused as "no single
  dissociation constant reported".

Run:  python3 -I tk_review/scripts/whats_left_albumin.py
"""
import collections
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
K2K = os.path.join(HERE, "..", "db", "albumin_kd_to_k.csv")
SHEET = os.path.join(HERE, "..", "albumin_logK_Kd_revised.xlsx")


def fam(method):
    m = (method or "").lower()
    if "scanning fluorimetry" in m or "(dsf)" in m:
        return "DSF"
    if "thermophoresis" in m:
        return "MST"
    if "nmr" in m:
        return "19F NMR"
    return None


def main():
    import openpyxl
    ws = openpyxl.load_workbook(SHEET)["Coefficients"]
    h = [c.value for c in ws[1]]
    rows = [dict(zip(h, [c.value for c in r]))
            for r in ws.iter_rows(min_row=2)]
    insheet = {str(r["Source"]) for r in rows}
    print(f"SHEET: {len(rows)} rows, {len(insheet)} sources, "
          f"{len({str(r['PFAS']) for r in rows})} compounds, "
          f"species {sorted({str(r['Species']) for r in rows})}\n")

    print("A.  Rows with a Kd, held back only by the method-family gate")
    blocked = collections.Counter()
    for r in csv.DictReader(open(K2K, encoding="utf-8")):
        f = fam(r.get("method"))
        if f and r.get("kd_uM"):
            blocked[(f, str(r.get("source")))] += 1
    for (f, s), v in sorted(blocked.items(), key=lambda x: -x[1]):
        print(f"      {v:4d}  {f:8s} {s}")
    print(f"      {sum(blocked.values()):4d}  TOTAL\n")

    print("B.  Species the database has and the sheet does not")
    sp = collections.Counter()
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if "albumin" in (r.get("protein") or "").lower():
            sp[r.get("protein_species") or "?"] += 1
    sheet_sp = {str(r["Species"]) for r in rows}
    for k, v in sorted(sp.items(), key=lambda x: -x[1]):
        if k not in sheet_sp and v > 1:
            print(f"      {v:4d}  {k}   <-- absent from the sheet")
    print()

    print("C.  Albumin sources with a constant that reach no sheet row")
    db = collections.Counter()
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if ("albumin" in (r.get("protein") or "").lower()
                and (r.get("ka") or r.get("kd"))):
            base = re.sub(r"\s+et al\.?", "", r["study"] or "").strip()
            db[base.rstrip(",") + " " + (r["year"] or "")] += 1
    for k, v in sorted(db.items(), key=lambda x: -x[1]):
        if k not in insheet and v >= 2:
            print(f"      {v:4d}  {k}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
