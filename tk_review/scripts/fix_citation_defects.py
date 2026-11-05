#!/usr/bin/env python3
"""Two defects in db/citations.csv, both found by the identity check.

Run scripts/check_citation_identity.py before and after this.

DEFECT 1 -- A WRONG PMID, so a citation describing another paper.

Chen & Guo 2009 was filed under PMID 19002450.  That PMID is real and
resolves cleanly; it is just somebody else's paper.  Its DOI is
10.1007/s00247-008-1049-9 -- `s00247` is Pediatric Radiology, while
`s00204` is Archives of Toxicology, one digit apart -- and the
citation built from it was a review of splenic function in sickle
cell disease.

The DOI already stored against the 20 database rows,
10.1007/s00204-008-0359-x, was correct throughout, which is why the
DOI-series check sees the same error independently.  The correct
record is PMID 18854981: Chen YM, Guo LH, "Fluorescence study on
site-specific binding of perfluoroalkyl acids to human serum
albumin", Arch Toxicol 2009;83(3):255-61.  Confirmed against PubMed
by DOI and by author search; the abstract's constants (2.7e5 M^-1
PFOA, 2.2e4 M^-1 PFOS at Trp-214; 1.0e6 / 2.2e6 at site I; 6.5e6 /
1.2e6 at site II) match what this database already holds from
Table 1, so the ROWS were transcribed from the right paper all along.
Only the identifier was wrong.

This fixes the PMID in BOTH files it appears in -- citations.csv and
the 20 rows of protein_binding.csv -- so nothing downstream can pick
the stale one up again.

DEFECT 2 -- THE YEAR IN THE CITATION STRING IS THE ONLINE-FIRST YEAR.

PubMed's publication_date is, for Springer and some other
publishers, the date the paper first appeared online, not the date of
the issue it ended up in.  The citations were built from that field,
so nine rows pair one year with another year's volume: Sheng 2016
reads "Arch Toxicol 2014;90(1):217-27", and Arch Toxicol volume 90
issue 1 is January 2016.

This project's own `year` column is the issue year, and it is correct
on all nine -- each volume and issue belongs to the year the column
already gives.  So the fix is to make the citation's year agree with
the volume it is printed beside, which is what a reference list
needs.  Volume, issue and pages are untouched.

SAFETY.  The original citation strings were assembled by a script no
longer in the repository, so before rewriting any string this one
re-derives the string for every row it is NOT changing and checks it
reproduces the stored text byte for byte.  If the formatter here
disagrees with the one that built the file, the run aborts rather
than reformat 71 citations into a slightly different house style.

Run:  python3 tk_review/scripts/fix_citation_defects.py
Idempotent; rewrites db/citations.csv and db/protein_binding.csv.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CITES = os.path.join(HERE, "..", "db", "citations.csv")
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

BAD_PMID, GOOD_PMID = "19002450", "18854981"

# The correct record, from PubMed, for the study the label names.
CHEN_GUO = dict(
    pmid=GOOD_PMID,
    doi="10.1007/s00204-008-0359-x",
    authors="Chen YM; Guo LH",
    title="Fluorescence study on site-specific binding of "
          "perfluoroalkyl acids to human serum albumin",
    journal="Arch Toxicol",
    pub_year="2009",          # issue year; PubMed's e-pub is 2008-10-15
    volume="83", issue="3", pages="255-61",
    authority="PubMed 18854981, resolved by DOI and by author search "
              "after the stored PMID was found to describe a "
              "different paper")


def fmt(r):
    """Rebuild a citation string in the file's existing house style."""
    vol = r["volume"] + (f"({r['issue']})" if r["issue"] else "")
    tail = f"{r['journal']} {r['pub_year']}" + (f";{vol}" if vol else "")
    tail += f":{r['pages']}" if r["pages"] else ""
    return f"{r['authors']}. {r['title']}. {tail}."


def main():
    with open(CITES, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fields, rows = list(reader.fieldnames), list(reader)

    # ---- which rows this run touches -------------------------------
    chen = [r for r in rows if r["study"] == "Chen & Guo"]
    offyear = []
    for r in rows:
        if not r["citation"] or r in chen:
            continue
        try:
            if int(r["year"]) != int(r["pub_year"]):
                offyear.append(r)
        except (ValueError, TypeError):
            pass
    touched = {id(r) for r in chen + offyear}

    # ---- SAFETY: the formatter must match the one that built the
    # ---- file, proved on every row this run leaves alone.
    checked = mismatch = 0
    for r in rows:
        if id(r) in touched or not r["citation"]:
            continue
        checked += 1
        if fmt(r) != r["citation"]:
            if mismatch < 3:
                print("  formatter disagrees on an untouched row:")
                print(f"    stored  {r['citation']!r}")
                print(f"    rebuilt {fmt(r)!r}")
            mismatch += 1
    if mismatch:
        print(f"\nABORT: {mismatch} of {checked} untouched citations "
              f"do not round-trip.")
        print("The house style here is not the one that built the "
              "file; nothing was written.")
        return 1
    print(f"formatter round-trips on all {checked} untouched "
          f"citations -- safe to rebuild strings\n")

    # ---- DEFECT 1 --------------------------------------------------
    for r in chen:
        was = r["pmid"]
        r.update(CHEN_GUO)
        r["citation"] = fmt(r)
        print(f"DEFECT 1  Chen & Guo {r['year']}: pmid {was} -> "
              f"{r['pmid']}")
        print(f"          {r['citation']}")
    if not chen:
        print("DEFECT 1  no Chen & Guo row found")

    # ---- DEFECT 2 --------------------------------------------------
    print()
    for r in offyear:
        was = r["pub_year"]
        r["pub_year"] = r["year"]
        r["citation"] = fmt(r)
        print(f"DEFECT 2  {r['study']} {r['year']}: citation year "
              f"{was} -> {r['pub_year']} "
              f"(volume {r['journal']} {r['volume']} is "
              f"{r['year']})")
    if not offyear:
        print("DEFECT 2  no row carries an online-first year")

    with open(CITES, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # ---- the same PMID, in the database rows -----------------------
    with open(DB, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        dbf, dbr = list(reader.fieldnames), list(reader)
    n = 0
    for r in dbr:
        if r["pmid"] == BAD_PMID:
            r["pmid"] = GOOD_PMID
            n += 1
    if n:
        with open(DB, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=dbf)
            w.writeheader()
            w.writerows(dbr)
    print(f"\nprotein_binding.csv: {n} row(s) moved off the stale "
          f"PMID")
    print(f"citations.csv: {len(chen)} identity fix, "
          f"{len(offyear)} year fix")
    return 0


if __name__ == "__main__":
    sys.exit(main())
