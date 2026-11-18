#!/usr/bin/env python3
"""Dump lfabp_logK_Kd.xlsx's Coefficients sheet as TSV for pasting into Excel.

Read-only; the workbook is the source. Run after
build_lfabp_working_table.py.
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
XL = os.path.join(HERE, "..", "lfabp_logK_Kd.xlsx")
OUT = os.path.join(HERE, "..", "lfabp_rows_for_excel.tsv")


def main():
    import openpyxl
    ws = openpyxl.load_workbook(XL)["Coefficients"]
    n = 0
    with open(OUT, "w", encoding="utf-8") as fh:
        for row in ws.iter_rows():
            if row[0].value is None:
                continue
            fh.write("\t".join("" if c.value is None else str(c.value)
                               for c in row) + "\n")
            n += 1
    print(f"wrote {OUT}  ({n - 1} data rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
