#!/usr/bin/env python3
"""Add Harada 2005, Zhang 2013 and Fujii 2015 to db/citations.csv. Idempotent.

    python3 -I scripts/add_primary_renal_citations.py

These three supply every DIRECT human renal clearance used in the analysis,
and reading them resolved the apparent self-contradiction in EFSA's
tabulation: the per-compound table reproduces Zhang correctly and transcribes
Fujii's TOTAL clearance column as renal clearance.
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "citations.csv")
ROWS = [
    {"study": "Harada et al.", "year": "2005", "pmid": "16194675",
     "doi": "10.1016/j.envres.2004.12.003",
     "authors": "Harada K; Inoue K; Morikawa A; Yoshinaga T; Saito N; Koizumi A",
     "title": ("Renal clearance of perfluorooctane sulfonate and "
               "perfluorooctanoate in humans and their species-specific "
               "excretion"),
     "journal": "Environ Res", "pub_year": "2005", "volume": "99",
     "issue": "2", "pages": "253-261", "authority": "PubMed"},
    {"study": "Zhang et al.", "year": "2013", "pmid": "23980546",
     "doi": "10.1021/es401905e",
     "authors": "Zhang Y; Beesoon S; Zhu L; Martin JW",
     "title": ("Biomonitoring of perfluoroalkyl acids in human urine and "
               "estimates of biological half-life"),
     "journal": "Environ Sci Technol", "pub_year": "2013", "volume": "47",
     "issue": "18", "pages": "10619-10627", "authority": "PubMed"},
    {"study": "Fujii et al.", "year": "2015", "pmid": "25422127",
     "doi": "10.1539/joh.14-0136-OA",
     "authors": ("Fujii Y; Niisoe T; Harada KH; Uemoto S; Ogura Y; "
                 "Takenaka K; Koizumi A"),
     "title": ("Toxicokinetics of perfluoroalkyl carboxylic acids with "
               "different carbon chain lengths in mice and humans"),
     "journal": "J Occup Health", "pub_year": "2015", "volume": "57",
     "issue": "1", "pages": "1-12", "authority": "PubMed"},
]
for r in ROWS:
    r["citation"] = (f"{r['authors']}. {r['title']}. {r['journal']} "
                     f"{r['pub_year']};{r['volume']}"
                     f"({r['issue']}):{r['pages']}.")


def main():
    with open(CSV) as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    have = {r.get("doi") for r in rows}
    added = [r for r in ROWS if r["doi"] not in have]
    if not added:
        print("all three already present; nothing written")
        return
    for r in added:
        rows.append({k: r.get(k, "") for k in fields})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"appended {len(added)} -> {len(rows)} rows: "
          + ", ".join(r["study"] + " " + r["year"] for r in added))


if __name__ == "__main__":
    main()
