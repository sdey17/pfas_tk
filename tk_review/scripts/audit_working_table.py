#!/usr/bin/env python3
"""Audit a working spreadsheet against the database, and list what it lacks.

The reviewer keeps a working copy of the albumin table in Excel and
edits it by hand.  It drifts from db/albumin_logK_with_conc.xlsx in
both directions: the workbook gains corrected rows the copy has not
picked up, and the copy holds choices the workbook does not make
(one salt level per compound from Gao 2025, one representative row
per compound from Moro 2022).  Re-deriving the differences by hand
each time is how a duplicate survives several passes.

WHAT IT CHECKS, and why each one earns its place:

  1. log K = BRIDGE + log10(n) - log10(Kd) on every row carrying
     both.  This is the one column nothing else can validate, and a
     single wrong digit here moves a compound an order of magnitude.

  2. ligand:protein ratio = PFAS conc / protein conc.  Reported, NOT
     failed, because it is legitimately false for some rows: in
     Allendorf 2019 the ratio column is the paper's own reported v,
     bound per protein, which is not the quotient of the two
     concentration columns -- PFECHS differs 8x, and that gap is in
     Allendorf's own SI.  A row that fails here needs a reason, not
     a correction.

  3. Duplicate (PFAS, species, source, log K).  Chases the failure
     mode this project has already hit twice: the same measurement
     entering once from a primary and once from a later paper's
     comparison table, invisible because the two spell the compound
     differently (PFDoA vs PFDoDA).

  4. Sources and row counts in each direction.  Separates "missing
     because nobody added it" from "absent on purpose".

It writes the rows the working copy lacks, in the working copy's own
column order, ready to paste.  It changes nothing, reads both files
read-only, and makes no judgement about which side is right.

Run:  python3 tk_review/scripts/audit_working_table.py <working.xlsx>
      [--sheet Coefficients] [-o missing.tsv]
"""
import argparse
import collections
import math
import os
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
OURS = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
BRIDGE = 4.1774                      # log10(1e9 / 66463), albumin
TOL_LOGK, TOL_RATIO = 0.02, 0.02     # absolute log units; relative

# The working copy's headers, and the workbook column each maps from
# where the two differ.
ALIAS = {"ligand:protein ratio": "ligand:protein ratio (max)",
         "PFAS conc (uM)": "PFAS conc (uM)",
         "PFAS conc (μM)": "PFAS conc (uM)",
         "Protein conc (uM)": "Protein conc (uM)",
         "Protein conc (μM)": "Protein conc (uM)"}


def col(row, name):
    """A column by the working copy's name, tolerating the micro-sign
    and the workbook's own spelling."""
    for k in (name, ALIAS.get(name, name)):
        if k in row:
            return row[k]
    want = name.replace("μ", "u")
    for k, v in row.items():
        if str(k).replace("μ", "u") == want:
            return v
    return None


def num(x):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return None


def load(path, sheet=None):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet] if sheet and sheet in wb.sheetnames else (
        wb["Coefficients"] if "Coefficients" in wb.sheetnames
        else wb.worksheets[0])
    hdr = [c.value for c in ws[1]]
    rows = []
    for n, r in enumerate(ws.iter_rows(min_row=2, values_only=True), 2):
        if any(v is not None and str(v).strip() for v in r):
            d = dict(zip(hdr, r))
            d["_row"] = n
            rows.append(d)
    return ws.title, hdr, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("working", help="path to the working .xlsx")
    ap.add_argument("--sheet", default=None)
    ap.add_argument("-o", "--out", default="", help="'-' for stdout")
    args = ap.parse_args()

    title, uhdr, U = load(args.working, args.sheet)
    _, _, O = load(OURS)
    print(f"  working copy : {os.path.basename(args.working)} "
          f"[{title}]  {len(U)} rows")
    print(f"  database     : {os.path.relpath(OURS, os.path.join(HERE,'..'))}"
          f"  {len(O)} rows\n")

    print("  1. log K vs BRIDGE + log10(n) - log10(Kd)")
    bad = []
    for r in U:
        kd, n, lk = (num(col(r, "Kd (uM)")), num(col(r, "n used")),
                     num(col(r, "log K (L water / kg protein)")))
        if kd and n and lk and kd > 0 and n > 0:
            e = BRIDGE + math.log10(n) - math.log10(kd)
            if abs(e - lk) > TOL_LOGK:
                bad.append((r["_row"], col(r, "PFAS"), col(r, "Source"), lk, e))
    print(f"     {len(bad)} mismatch(es)" + ("" if bad else "  -- clean"))
    for b in bad[:20]:
        print(f"       row {b[0]:>4} {str(b[1])[:20]:22}{str(b[2])[:18]:20}"
              f"file {b[3]}  expected {b[4]:.3f}")

    print("\n  2. ratio vs PFAS conc / protein conc  (reported, not failed)")
    off = collections.defaultdict(list)
    for r in U:
        p, pr = num(col(r, "PFAS conc (μM)")), num(col(r, "Protein conc (μM)"))
        rt = num(col(r, "ligand:protein ratio"))
        if p and pr and rt and pr > 0:
            e = p / pr
            if abs(e - rt) > max(TOL_RATIO * max(e, rt), 1e-9):
                off[str(col(r, "Source"))].append(r["_row"])
    if not off:
        print("     every row divides")
    for s, rr in sorted(off.items()):
        print(f"     {s:22}{len(rr):>4} row(s), e.g. {rr[:4]}")
        print(f"       -> check whether this column is the paper's own "
              f"reported v rather than the quotient")

    print("\n  3. duplicate (PFAS, species, source, log K)")
    seen = collections.defaultdict(list)
    for r in U:
        seen[(str(col(r, "PFAS")), str(col(r, "Species")),
              str(col(r, "Source")),
              str(col(r, "log K (L water / kg protein)")))].append(r["_row"])
    dup = {k: v for k, v in seen.items() if len(v) > 1}
    print(f"     {len(dup)} duplicated key(s)" + ("" if dup else "  -- clean"))
    for k, v in sorted(dup.items()):
        print(f"       rows {v}  {k[0][:22]:24}{k[2][:18]:20}log K {k[3]}")

    print("\n  4. sources")
    uc = collections.Counter(str(col(r, "Source")) for r in U)
    oc = collections.Counter(str(col(r, "Source")) for r in O)
    miss = [s for s in oc if s not in uc]
    print(f"     in the database, absent from the working copy: {len(miss)}")
    for s in sorted(miss, key=lambda s: -oc[s]):
        sp = sorted({str(col(r, "Species")) for r in O
                     if str(col(r, "Source")) == s})
        print(f"       {s:20}{oc[s]:>4} rows  {sp}")
    diff = [s for s in set(uc) & set(oc) if uc[s] != oc[s]]
    print(f"     present in both with different counts: {len(diff)}")
    for s in sorted(diff):
        print(f"       {s:20} working {uc[s]:>4}   database {oc[s]:>4}"
              f"   ({oc[s]-uc[s]:+d})")
    only = [s for s in uc if s not in oc]
    if only:
        print(f"     in the working copy only: {sorted(only)}")

    out = []
    for r in O:
        if str(col(r, "Source")) in miss:
            out.append([("" if col(r, h) is None else str(col(r, h)))
                        for h in uhdr])
    out.sort(key=lambda x: (x[uhdr.index("Source")], x[0]))
    body = "\t".join(str(h) for h in uhdr) + "\n" + \
           "\n".join("\t".join(c) for c in out) + "\n"
    if args.out == "-":
        sys.stdout.write(body)
    else:
        path = args.out or os.path.join(HERE, "..",
                                        "missing_rows_for_user_table.tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(body)
        print(f"\n  wrote {len(out)} missing rows to "
              f"{os.path.relpath(path, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
