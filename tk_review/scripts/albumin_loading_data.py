#!/usr/bin/env python3
"""Build the CSV the loading plot reads. Workbook in, tidy table out.

Separated from the plotting so the plot script can live in a notebook
with no paths in it: this is the part that needs to know where the
project keeps things, and it runs in the repo.

Writes figures/albumin_loading_vs_logk_data.csv with one row per
value: pfas, n_cf, log_ratio, log_k, category, source. `category` is
how the number was arrived at, which is what the plot colours by.

Allendorf's log_ratio is the analyte:BSA ratio its SI Table SI-2.3
measured, not the "below 0.1" ceiling its main text quotes. SI-2.3
covers batch A3803 only and the A7906 row inherits it, which is an
assumption rather than a reading.

Run:  python3 tk_review/scripts/albumin_loading_data.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
OUT_CSV = os.path.join(HERE, "..", "figures",
                       "albumin_loading_vs_logk_data.csv")

WANT = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]        # chain-length order
NCF = {"PFBA": 3, "PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}
FIELDS = ["pfas", "n_cf", "log_ratio", "log_k", "category", "source"]


def from_workbook():
    from openpyxl import load_workbook
    ws = load_workbook(SRC)["Coefficients"]
    head = [c.value for c in ws[1]]
    ix = {h: i for i, h in enumerate(head)}
    rows = []
    for r in ws.iter_rows(min_row=2, values_only=True):
        name = str(r[ix["PFAS"]])
        if name not in WANT:
            continue
        route = str(r[ix["Route"]])
        src = str(r[ix["Source"]])
        cat = ("qin" if src.startswith("Qin 2024") else
               "measured" if route == "measured directly" else
               "quenching" if "quenching" in route else "converted")
        try:
            rat = float(str(r[ix["ligand:protein ratio (max)"]]).strip())
        except (TypeError, ValueError):
            rat = None
        rows.append({
            "pfas": name, "n_cf": NCF[name],
            "log_ratio": "" if not rat or rat <= 0 else round(
                math.log10(rat), 4),
            "log_k": float(r[ix["log K (L water / kg protein)"]]),
            "category": cat, "source": src})
    return rows


def main():
    rows = from_workbook()
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["n_cf"], r["category"],
                                                r["log_k"])))
    print(f"wrote {OUT_CSV}")
    print(f"  {len(rows)} rows over "
          f"{len({r['pfas'] for r in rows})} compounds")
    from collections import Counter
    for k, v in sorted(Counter(r["category"] for r in rows).items()):
        n = sum(1 for r in rows if r["category"] == k
                and r["log_ratio"] != "")
        print(f"  {k:<12}{v:>4} rows, {n} with a ratio")


if __name__ == "__main__":
    main()
