#!/usr/bin/env python3
"""Add Abraham et al. 2024 to db/citations.csv. Idempotent.

    python3 -I scripts/add_abraham2024_citation.py

This study supplies EVERY human volume of distribution and total clearance
used in the half-life decomposition, and it was absent from the bibliography
while being cited by name in four scripts and three reports. Row-level
provenance was correct ("Abraham 2024 Tables 2-4"); only the citation table
was short.
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "citations.csv")
ROW = {
    "study": "Abraham et al.", "year": "2024", "pmid": "39476597",
    "doi": "10.1016/j.envint.2024.109047",
    "citation": ("Abraham K; Mertens H; Richter L; Mielke H; Schwerdtle T; "
                 "Monien BH. Kinetics of 15 per- and polyfluoroalkyl "
                 "substances (PFAS) after single oral application as a "
                 "mixture - A pilot investigation in a male volunteer. "
                 "Environ Int 2024;193:109047."),
    "authors": ("Abraham K; Mertens H; Richter L; Mielke H; Schwerdtle T; "
                "Monien BH"),
    "title": ("Kinetics of 15 per- and polyfluoroalkyl substances (PFAS) "
              "after single oral application as a mixture - A pilot "
              "investigation in a male volunteer"),
    "journal": "Environ Int", "pub_year": "2024", "volume": "193",
    "issue": "", "pages": "109047", "authority": "PubMed",
}


def main():
    with open(CSV) as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    if any(r.get("doi") == ROW["doi"] for r in rows):
        print("already present; nothing written")
        return
    rows.append({k: ROW.get(k, "") for k in fields})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"appended Abraham 2024 -> {len(rows)} rows")


if __name__ == "__main__":
    main()
