#!/usr/bin/env python3
"""Add Kudo 2001, Kudo 2002 and Ohmori 2003 to db/citations.csv. Idempotent.

    python3 -I scripts/add_kudo_citations.py

All three are cited by name across the project's scripts and reports -- Kudo
2001 for urinary and faecal fractions of dose, Kudo 2002 Tables 2 and 3 for
rat PFOA renal clearance by sex and hormonal state, Ohmori 2003 Table 1 and
Figure 2 for chain-length toxicokinetics -- and none was in the bibliography.

A NOTE ON ONE CITATION CHAIN. Harada 2005 Table 5 footnote b attributes its
rat half-lives and clearances to "Kudo et al. (2002)", and the values it
prints (male t1/2 5.63 d, female 0.08 d) are those reported in the ABSTRACT
of Ohmori 2003. Ohmori 2003 appeared in Toxicology volume 184 dated March
2003 but was available online in 2002, so this is almost certainly an
online-first citation rather than an error; Kudo is second author of it. The
numbers themselves agree with what this project already holds as "Kudo 2002
Table 3" (male PFOA CL_R 46.08, female 1054.1) against Harada's 46.1 and
1054.1, so the two routes to those values are consistent.
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "citations.csv")
ROWS = [
    {"study": "Kudo et al.", "year": "2001", "pmid": "11311214",
     "doi": "10.1016/s0009-2797(01)00155-7",
     "authors": "Kudo N; Suzuki E; Katakura M; Ohmori K; Noshiro R; Kawashima Y",
     "title": ("Comparison of the elimination between perfluorinated fatty "
               "acids with different carbon chain length in rats"),
     "journal": "Chem Biol Interact", "pub_year": "2001", "volume": "134",
     "issue": "2", "pages": "203-216", "authority": "PubMed"},
    {"study": "Kudo et al.", "year": "2002", "pmid": "11879818",
     "doi": "10.1016/s0009-2797(02)00006-6",
     "authors": "Kudo N; Katakura M; Sato Y; Kawashima Y",
     "title": "Sex hormone-regulated renal transport of perfluorooctanoic acid",
     "journal": "Chem Biol Interact", "pub_year": "2002", "volume": "139",
     "issue": "3", "pages": "301-316", "authority": "PubMed"},
    {"study": "Ohmori et al.", "year": "2003", "pmid": "12499116",
     "doi": "10.1016/s0300-483x(02)00573-5",
     "authors": "Ohmori K; Kudo N; Katayama K; Kawashima Y",
     "title": ("Comparison of the toxicokinetics between perfluorocarboxylic "
               "acids with different carbon chain length"),
     "journal": "Toxicology", "pub_year": "2003", "volume": "184",
     "issue": "2-3", "pages": "135-140", "authority": "PubMed"},
]
for r in ROWS:
    r["citation"] = (f"{r['authors']}. {r['title']}. {r['journal']} "
                     f"{r['pub_year']};{r['volume']}({r['issue']}):"
                     f"{r['pages']}.")


def main():
    with open(CSV) as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    have = {r.get("doi") for r in rows}
    add = [r for r in ROWS if r["doi"] not in have]
    if not add:
        print("all present; nothing written")
        return
    for r in add:
        rows.append({k: r.get(k, "") for k in fields})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"appended {len(add)} -> {len(rows)} rows: "
          + "; ".join(f"{r['study']} {r['year']}" for r in add))


if __name__ == "__main__":
    main()
