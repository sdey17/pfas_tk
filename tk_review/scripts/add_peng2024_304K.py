#!/usr/bin/env python3
"""Peng 2024's missing middle temperature: the 304 K series.

Peng Y, et al.  Toxics 2024;12(1):43.  doi 10.3390/toxics12010043,
PMID 38250999.  Open access.  Table 2, transcribed as printed.

This database held two of Peng's three temperatures.  Table 2 reports
Kb, n, KSV and Kq at 298, 304 and 310 K; 304 K was never ingested, so
the temperature dependence could be described but not plotted, and
the paper's van't Hoff enthalpies -- which reach +412 kJ/mol -- could
not be checked against the series they came from.

With three points the check becomes possible: a van't Hoff plot is a
straight line in (1/T, ln K) only if a single enthalpy describes the
range.  See plot_peng2024_vanthoff.py.

Run:  python3 tk_review/scripts/add_peng2024_304K.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Peng et al.", "2024"
DOI, PMID = "10.3390/toxics12010043", "38250999"
PROTEIN = "human serum albumin (HSA)"
METHOD = ("fluorescence quenching (static; Kb from modified "
          "Stern-Volmer / double-log plot) with warfarin, ibuprofen "
          "and lidocaine site-marker competition")

# Table 2, 304 K row: (chemical, nCF or None, head group, Kb, n, KSV)
T304 = [
    ("PFNA", 8, "carboxylate", "2.02e5", "1.18", "2.35"),
    ("HFPO-TA", None, "ether carboxylate", "8.54e7", "1.78", "1.95"),
    ("PFOA", 7, "carboxylate", "2.23e6", "1.31", "1.46"),
    ("PFO3DA (perfluoro-3,6,9-trioxadecanoic acid)", None,
     "ether carboxylate", "4.55e4", "1.33", "1.03"),
    ("PFHpA", 6, "carboxylate", "5.35e5", "1.42", "0.53"),
    ("DFSA (dodecafluorosuberic acid)", None, "dicarboxylate",
     "4.93e3", "0.98", "0.51"),
]

NOTE = ("Table 2, 304 K row | Kb and n at 304 K, the middle of Peng's "
        "three temperatures (298, 304, 310 K) | this row completes the "
        "series so the paper's van't Hoff enthalpies can be checked "
        "against it; see plot_peng2024_vanthoff.py | the three "
        "temperatures are ONE experiment, not three studies, and a "
        "physiological calculation wants the 310 K row")
TK = ("partly - discussed as driver of bioaccumulation/transport; no "
      "explicit half-life or Vd linkage")


def blank(fields):
    return dict.fromkeys(fields, "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, list(reader)

    have = {(r["chemical"], r["ka"]) for r in rows
            if r["study"] == STUDY and r["year"] == YEAR}
    new = []
    for chem, ncf, head, kb, n, ksv in T304:
        if (chem, kb) in have:
            continue                      # already ingested
        r = blank(fields)
        r.update(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            n_perfluorinated_carbons="" if ncf is None else str(ncf),
            head_group=head, protein=PROTEIN, protein_species="human",
            method=METHOD, ka=kb, ka_units="L/mol (= M^-1)",
            n_binding_sites=n, protein_conc="1.0",
            protein_conc_units="uM", ligand_conc_range="0-18",
            ligand_conc_units="uM", temperature_c="31 (304 K)",
            is_computational="no", tk_relevance_noted=TK,
            source_table=f"{NOTE} | KSV {ksv} x 10^4 L/mol")
        new.append(r)

    print(f"{len(new)} rows to add"
          + (" (already present)" if not new else ""))
    for r in new:
        print(f"  {r['chemical'][:40]:42}Kb {r['ka']:9} n {r['n_binding_sites']}")
    if args.dry_run or not new:
        if args.dry_run:
            print("--dry-run: nothing written")
        return
    with open(CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows + new)
    print(f"wrote {os.path.relpath(CSV, os.path.join(HERE, '..'))}")


if __name__ == "__main__":
    main()
