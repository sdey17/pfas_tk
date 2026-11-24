#!/usr/bin/env python3
"""Add the Thakur 2024 kidney proteomics citation.

    python3 -I scripts/add_thakur_citation.py

Appends to db/citations.csv, keyed on PMID. Idempotent.

WHY IT WAS MISSING, AND WHAT THE DRAFT GOT WRONG. Thakur 2024 is the single
source for the transporter-abundance section of the paper -- the claim that
human OAT4 is undetectable in rat and mouse and rodent Oatp1a1 undetectable
in human rests on it -- yet it had no row in db/citations.csv, and the
reference entry in report/PAPER_DRAFT.md carried a DESCRIPTIVE TITLE THAT IS
NOT THE PAPER'S TITLE ("Quantitative proteomic analysis of renal transporters
in human, rat and mouse kidney") together with a single page number, 235,
which had been taken for a meeting-abstract page. PubMed metadata shows the
real title is "Sex and the Kidney Drug-Metabolizing Enzymes and Transporters:
Are Preclinical Drug Disposition Data Translatable to Humans?" and the real
pagination 235-246, a full research paper. Both are now corrected from
metadata rather than memory.

The paper is open access as PMC11218045, so the abundance values the project
relies on can be verified against the source directly.
"""
import csv
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

NEW = dict(
    study="Thakur et al.", year="2024", pmid="38711199",
    doi="10.1002/cpt.3277",
    authors="Thakur A; Yue G; Ahire D; Mettu VS; Al Maghribi A; Ford K; "
            "Peixoto L; Leeder JS; Prasad B",
    title="Sex and the kidney drug-metabolizing enzymes and transporters: "
          "are preclinical drug disposition data translatable to humans?",
    journal="Clinical Pharmacology and Therapeutics", pub_year="2024",
    volume="116", issue="1", pages="235-246", authority="PubMed",
)


def main():
    p = HERE / "db" / "citations.csv"
    rows = list(csv.DictReader(p.open()))
    cols = list(rows[0].keys())

    if NEW["pmid"] in {r["pmid"] for r in rows if r["pmid"]}:
        print(f"db/citations.csv: {len(rows)} rows "
              f"(Thakur 2024 already present -- nothing to do)")
        return

    d = dict(NEW)
    d["citation"] = (f"{d['authors']}. {d['title']}. {d['journal']} "
                     f"{d['pub_year']};{d['volume']}({d['issue']}):"
                     f"{d['pages']}.")
    rows.append({c: d.get(c, "") for c in cols})
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"db/citations.csv: {len(rows)} rows (+ Thakur et al. 2024)")
    print(f"  9 authors, PMID {d['pmid']}, PMC11218045 (open access)")


if __name__ == "__main__":
    main()
