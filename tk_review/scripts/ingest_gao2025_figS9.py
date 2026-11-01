#!/usr/bin/env python3
"""Gao et al. 2025 Figure S9 -- the per-compound titration series.

WHAT THIS CLOSES. 42 rows from Gao 2025 carried "titration to
fluorescence saturation; per-compound series in Figures S5-S12" and
reached no molar concentration, so the isotherm-position columns could
not be computed for the second-largest albumin source in the database.
The supporting information was obtained and Figure S9 prints the series
outright, in its legend, for all ten compounds.

THE SERIES IS NOT ONE SERIES. It is chosen per compound and it spans
three orders of magnitude between them -- PFOS tops out at 0.01 mM and
PFBS at 15 mM -- which is exactly what a titration to saturation
requires, since the affinities differ by as much. A single "0-15 mM"
for the source would have been useless; the per-compound number is the
point.

SCOPE, STATED PRECISELY. Figure S9 is the NO-SALT condition: 0.8 uM BSA
in 20 mM NaH2PO4-Na2HPO4 (pH 7) at 20 C. The 42 rows are a NaCl series
at 0, 10, 30, 80, 150, 300 and 500 mM. The series is carried to the
salt rows as well, on the evidence of Figure S5, which shows PFHpA
titrated across eight different salt conditions with an identical
legend -- the same 14 concentrations, 0 to 8 mM, as Figure S9(c) gives
it with no salt. That is direct evidence for one compound and an
inference for the other nine, within one paper and one experimental
design, and each row says so rather than presenting it as read.

Writes db/gao2025_figS9_titration_series.csv.

Run:  python3 -I tk_review/scripts/ingest_gao2025_figS9.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "gao2025_figS9_titration_series.csv")

SOURCE = ("Gao et al. 2025, Environ Sci Technol, doi "
          "10.1021/acs.est.5c08779, supporting information Figure S9 "
          "(legend); no-salt condition, 0.8 uM BSA in 20 mM phosphate "
          "pH 7 at 20 C")
CARRIED = ("Figure S5 shows PFHpA titrated across eight salt conditions "
           "with the same 14 concentrations Figure S9(c) gives with no "
           "salt, so the series is taken to be invariant to salt")

# (compound, the legend's concentrations in mM, verbatim and in order)
SERIES = [
    ("PFBA",  [0, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10]),
    ("PFHxA", [0, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 8, 10]),
    ("PFHpA", [0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
               1, 2, 5, 8]),
    ("PFOA",  [0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02,
               0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5]),
    ("PFNA",  [0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]),
    ("PFDA",  [0, 0.0001, 0.0002, 0.0005, 0.001, 0.0015, 0.002, 0.005,
               0.01, 0.02, 0.1]),
    ("PFBS",  [0, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15]),
    ("PFHxS", [0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02,
               0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6]),
    ("PFOS",  [0, 0.00005, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001,
               0.002, 0.004, 0.01]),
    ("HFPO-DA (GenX)", [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4,
                        5]),
]

FIELDS = ["pfas", "n_levels", "lowest_nonzero_mM", "highest_mM",
          "highest_uM", "series_mM", "protein_conc_uM", "source",
          "carried_to_salt_rows_because"]


def _fmt(x):
    return ("%g" % x)


def rows():
    out = []
    for name, series in SERIES:
        nz = [c for c in series if c > 0]
        out.append({
            "pfas": name, "n_levels": len(series),
            "lowest_nonzero_mM": _fmt(min(nz)),
            "highest_mM": _fmt(max(series)),
            "highest_uM": _fmt(max(series) * 1000.0),
            "series_mM": "; ".join(_fmt(c) for c in series),
            "protein_conc_uM": "0.8", "source": SOURCE,
            "carried_to_salt_rows_because": CARRIED})
    return out


def main():
    recs = rows()
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(recs)
    print(f"{len(recs)} compounds from Figure S9, 0.8 uM BSA throughout")
    print(f"  {'compound':16s}{'levels':>7s}{'lowest':>11s}{'highest':>11s}"
          f"   ligand:protein at the top")
    for r in recs:
        top = float(r["highest_uM"])
        print(f"  {r['pfas']:16s}{r['n_levels']:>7d}"
              f"{r['lowest_nonzero_mM'] + ' mM':>11s}"
              f"{r['highest_mM'] + ' mM':>11s}   {top / 0.8:>10,.0f} : 1")
    tops = [float(r["highest_uM"]) for r in recs]
    print(f"\n  the series spans {min(tops):g} to {max(tops):g} uM across "
          f"compounds -- {max(tops) / min(tops):,.0f}x, which is why one "
          f"range for the source would say nothing")
    print(f"\nwrote {os.path.abspath(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
