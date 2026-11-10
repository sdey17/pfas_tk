#!/usr/bin/env python3
"""Export the albumin working table as CSV for the Jupyter notebook.

The project's artifact is an .xlsx with two sheets.  A notebook wants a
plain CSV it can read with the standard library, so this writes one file
per sheet into tk_review/notebook/ and leaves the workbook untouched.

Run:  python3 -I tk_review/scripts/export_notebook_bundle.py
"""
import csv
import os

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
SRC = os.path.join(ROOT, "albumin_logK_Kd_revised.xlsx")
OUT = os.path.join(ROOT, "notebook")

SHEETS = {"Coefficients": "albumin_logK_Kd_revised.csv",
          "Notes": "albumin_logK_Kd_revised_notes.csv"}


def main():
    os.makedirs(OUT, exist_ok=True)
    wb = openpyxl.load_workbook(SRC, data_only=True)
    for sheet, fn in SHEETS.items():
        if sheet not in wb.sheetnames:
            print(f"  (no sheet {sheet!r}, skipped)")
            continue
        ws = wb[sheet]
        path = os.path.join(OUT, fn)
        n = 0
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            for row in ws.iter_rows(values_only=True):
                if all(c is None for c in row):
                    continue
                w.writerow(["" if c is None else c for c in row])
                n += 1
        print(f"  wrote notebook/{fn}  ({n - 1} rows)")


if __name__ == "__main__":
    main()
