#!/usr/bin/env python3
"""Add the Ryu 2024 PLASMA paper to db/citations.csv. Idempotent.

    python3 -I scripts/add_ryu2024_toxics_citation.py

WHY THIS WAS MISSING, AND WHY IT MATTERS. There are TWO Ryu et al. 2024
papers from the same group, and the project uses both:

  Ryu 2024a  Toxics 12(4):253, doi 10.3390/toxics12040253, PMID 38668476
             "Species-Specific Unbound Fraction Differences in Highly Bound
             PFAS: A Comparative Study across Human, Rat, and Mouse PLASMA
             and Albumin."  14 PFAS, C4-C12, presaturation equilibrium
             dialysis, plasma + albumin of human, rat, CD-1 and C57BL/6.
             -> db/plasma_unbound_fraction.csv

  Ryu 2024b  Environ Sci Technol 58(33):14641-14650,
             doi 10.1021/acs.est.4c04050, PMID 39161261
             "Unbound Fractions of PFAS in Human and Rodent TISSUES: Rat
             Liver a Suitable Proxy for Evaluating Emerging PFAS?"  16 PFAS,
             eta_pfc 3-13, HTD96 dialysis, liver/brain/heart/kidney/lung.
             -> db/protein_binding.csv

Only 2024b was in citations.csv. 2024a is the one that supplies EVERY plasma
unbound fraction in the tubular-factor pipeline -- the rat 0.00298 and human
0.00245 that the whole species comparison rests on -- so the project's
bibliography was missing its most load-bearing source while carrying the
companion paper that the T calculation never uses.

The row-level provenance was correct throughout: the `source` field of
plasma_unbound_fraction.csv already reads "Ryu 2024, Toxics 12:253,
doi:10.3390/toxics12040253, Table 1". Only the bibliography was short.
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "citations.csv")

ROW = {
    "study": "Ryu et al. (plasma)",
    "year": "2024",
    "pmid": "38668476",
    "doi": "10.3390/toxics12040253",
    "citation": ("Ryu S; Burchett W; Zhang S; Modaresi SMS; Agudelo Areiza J; "
                 "Kaye E; Fischer FC; Slitt AL. Species-Specific Unbound "
                 "Fraction Differences in Highly Bound PFAS: A Comparative "
                 "Study across Human, Rat, and Mouse Plasma and Albumin. "
                 "Toxics 2024;12(4):253."),
    "authors": ("Ryu S; Burchett W; Zhang S; Modaresi SMS; Agudelo Areiza J; "
                "Kaye E; Fischer FC; Slitt AL"),
    "title": ("Species-Specific Unbound Fraction Differences in Highly Bound "
              "PFAS: A Comparative Study across Human, Rat, and Mouse Plasma "
              "and Albumin"),
    "journal": "Toxics",
    "pub_year": "2024",
    "volume": "12",
    "issue": "4",
    "pages": "253",
    "authority": "PubMed",
}


def main():
    with open(CSV) as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    if any(r.get("doi") == ROW["doi"] for r in rows):
        print("already present; nothing written")
    else:
        rows.append({k: ROW.get(k, "") for k in fields})
        with open(CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"appended Ryu 2024a (Toxics) -> {len(rows)} rows")
    for r in rows:
        if "ryu" in (r.get("study") or "").lower():
            print(f"  {r['study']:22} {r['journal']:26} {r['doi']}")


if __name__ == "__main__":
    main()
