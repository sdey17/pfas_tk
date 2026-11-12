#!/usr/bin/env python3
"""Stamp the K_b/K_SV finding onto Peng 2024's rows.

audit_peng2024_kb_vs_ksv.py shows that Peng's K_b diverges from the
Stern-Volmer constant of the same titration by a factor that is a clean
function of n (r2 ~0.95), because K_b is the intercept of a
double-logarithmic plot whose slope is n -- so K_b carries units of
(L/mol)^n and is not a binding constant where n is not 1.

The rows stay: the user wants the compilation exhaustive, and the Route
column already says these rest on a quenching exponent.  But the finding
belongs on the row, not only in a script nobody runs, so it goes in
assay_volume_note.

Idempotent: re-running changes nothing.

Run:  python3 -I tk_review/scripts/flag_peng2024_kb_units.py
"""
import csv
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
DOI = "10.3390/toxics12010043"
MARK = "K_b is the INTERCEPT of a double-logarithmic plot"
NOTE = (" | " + MARK + " whose slope is n, so it carries units of "
        "(L/mol)^n and is not a binding constant where n is not 1. "
        "Table 2's own Stern-Volmer constant, from the same titration, "
        "is 3e3-2.5e4 L/mol on every row; K_b/K_SV runs 0.3 to 116,667 "
        "and log10 of that ratio is linear in n (r2 ~0.95), so the "
        "divergence IS n. Worst case HFPO-TA at 310 K: n = 2.09, "
        "K_b = 2.31e9 (L/mol)^2.09, converting to log K 7.86 -- the "
        "highest in this compilation, against ~5.5 for anything measured "
        "directly. See audit_peng2024_kb_vs_ksv.py")


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        fields, rows = rd.fieldnames, list(rd)

    hit = [r for r in rows if r["doi"] == DOI]
    todo = [r for r in hit if MARK not in r["assay_volume_note"]]
    print(f"  {len(hit)} rows with doi {DOI}; {len(todo)} need the note")
    if not todo:
        print("  already stamped; nothing to do")
        return
    shutil.copy(DB, DB + ".bak")
    for r in todo:
        r["assay_volume_note"] = (r["assay_volume_note"] or "").rstrip() + NOTE
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  stamped {len(todo)} rows")


if __name__ == "__main__":
    main()
