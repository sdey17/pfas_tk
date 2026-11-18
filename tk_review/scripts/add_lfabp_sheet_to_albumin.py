#!/usr/bin/env python3
"""Add the L-FABP compilation to the reviewer's albumin workbook as extra sheets.

    python3 -I scripts/add_lfabp_sheet_to_albumin.py

albumin_logK_Kd_FINAL.xlsx is the reviewer's own hand-edited file and the
record for albumin, so this script is deliberately conservative about it:

  * their sheets are NOT renamed and their cells are NOT touched. The albumin
    data stays on 'Coefficients' and their (empty) 'Notes' stays empty.
  * two sheets are appended, 'L-FABP Coefficients' and 'L-FABP Notes', copied
    from lfabp_logK_Kd.xlsx.
  * the file is verified after writing: every albumin cell is compared against
    the pre-write copy, and the write is rolled back if any differ.

WHAT COULD GO WRONG, and what is checked. Their workbook carries an Excel
Table ('Table1') that Excel has extended across all 16,384 columns, which is
also why a 288-row file is 650 KB. openpyxl round-trips tables, but a silent
loss there would break sorting and filtering in their copy, so the check below
asserts the table is still present with the same ref after the save.

Re-running is safe: the two L-FABP sheets are removed and rewritten.
"""
import os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ALB = os.path.join(HERE, "..", "albumin_logK_Kd_FINAL.xlsx")
LF = os.path.join(HERE, "..", "lfabp_logK_Kd.xlsx")
NEW = {"Coefficients": "L-FABP Coefficients", "Notes": "L-FABP Notes"}


def snapshot(path):
    """Every non-empty cell of every pre-existing sheet, as a comparable dict."""
    import openpyxl
    wb = openpyxl.load_workbook(path)
    out = {}
    for ws in wb.worksheets:
        if ws.title in NEW.values():
            continue
        out[ws.title] = {(c.row, c.column): c.value
                         for row in ws.iter_rows() for c in row
                         if c.value is not None}
    return out


def main():
    import openpyxl
    from openpyxl.styles import Font
    from openpyxl.utils import get_column_letter

    before = snapshot(ALB)
    tables_before = {}
    wb = openpyxl.load_workbook(ALB)
    for ws in wb.worksheets:
        for name, ref in getattr(ws, "tables", {}).items():
            tables_before[(ws.title, name)] = getattr(ref, "ref", ref)

    src = openpyxl.load_workbook(LF)
    for old in NEW.values():                      # idempotent
        if old in wb.sheetnames:
            del wb[old]

    for s_name, d_name in NEW.items():
        s = src[s_name]
        d = wb.create_sheet(d_name)
        for row in s.iter_rows():
            for c in row:
                if c.value is None:
                    continue
                t = d.cell(row=c.row, column=c.column, value=c.value)
                if c.font and c.font.bold:
                    t.font = Font(bold=True)
        if s.freeze_panes:
            d.freeze_panes = s.freeze_panes
        for i in range(1, (s.max_column or 1) + 1):
            L = get_column_letter(i)
            if L in s.column_dimensions and s.column_dimensions[L].width:
                d.column_dimensions[L].width = s.column_dimensions[L].width
        print(f"  added {d_name!r}: {sum(1 for r in s.iter_rows() if r[0].value is not None)} rows")

    bak = ALB + ".bak"
    shutil.copy2(ALB, bak)
    wb.save(ALB)

    # --- verify, and roll back on any difference -------------------------
    after = snapshot(ALB)
    problems = []
    if set(before) != set(after):
        problems.append(f"sheet set changed: {sorted(before)} -> {sorted(after)}")
    for name in before:
        if before[name] != after.get(name):
            diff = {k for k in set(before[name]) | set(after.get(name, {}))
                    if before[name].get(k) != after.get(name, {}).get(k)}
            problems.append(f"{name}: {len(diff)} cell(s) changed, e.g. "
                            f"{sorted(diff)[:3]}")
    wb2 = openpyxl.load_workbook(ALB)
    tables_after = {}
    for ws in wb2.worksheets:
        for n, ref in getattr(ws, "tables", {}).items():
            tables_after[(ws.title, n)] = getattr(ref, "ref", ref)
    if tables_before != tables_after:
        problems.append(f"Excel tables changed: {tables_before} -> {tables_after}")
    for d_name in NEW.values():
        if d_name not in wb2.sheetnames:
            problems.append(f"{d_name} missing after save")

    if problems:
        shutil.copy2(bak, ALB)
        print("\nROLLED BACK, workbook restored. Problems:", file=sys.stderr)
        for p in problems:
            print("   " + p, file=sys.stderr)
        return 1

    os.remove(bak)
    alb_rows = len(set(r for r, _ in before["Coefficients"]) - {1})
    print(f"\nverified: {alb_rows} albumin data rows unchanged cell for cell, "
          f"Excel table intact")
    print(f"sheets now: {wb2.sheetnames}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
