#!/usr/bin/env python3
"""Fill the three resolvable missing identifiers in db/protein_binding.csv.

Eight rows reached the file through the OEHHA PHG Table A6.2 compilation with
neither a PMID nor a DOI. Three of the seven studies resolve to a single
PubMed record on an author + topic search, and each was confirmed by reading
the returned title and author list against the row it belongs to:

  Jones et al. 2003        PMID 14587903  10.1897/02-553
    "Binding of perfluorinated fatty acids to serum proteins", Environ Toxicol
    Chem 22(11):2639-49 -- the steroid-binding-protein displacement work the
    SHBP and CBG rows come from.
  Salvalaglio et al. 2010  PMID 21028884  10.1021/jp106584b
    "Determination of energies and sites of binding of PFOA and PFOS to human
    serum albumin", J Phys Chem B 114(46):14860-74 -- the in-silico rows.
  Woodcroft et al. 2010    PMID 20821618  10.1002/etc.199
    "Experimental characterization of the mechanism of perfluorocarboxylic
    acids' liver protein bioaccumulation", Environ Toxicol Chem 29(8):1669-77
    -- the L-FABP rows.

The remaining four are left blank rather than guessed: Messina 2005 returns no
PubMed record on an author+topic search, Gao 2019 and Li 2009 are ambiguous on
surname and year alone, and Kerstner-Wood 2003 is a Southern Research
Institute laboratory report, which PubMed does not index.

Identifiers retrieved from PubMed (https://pubmed.ncbi.nlm.nih.gov/).

Run:  python3 scripts/fill_binding_dois.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

RESOLVED = {
    ("Jones et al.", "2003"): ("14587903", "10.1897/02-553"),
    ("Salvalaglio et al.", "2010"): ("21028884", "10.1021/jp106584b"),
    ("Woodcroft et al.", "2010"): ("20821618", "10.1002/etc.199"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, list(reader)

    n = 0
    for r in rows:
        key = (r["study"], r["year"])
        if key not in RESOLVED:
            continue
        if r["pmid"].strip() or r["doi"].strip():
            continue
        pmid, doi = RESOLVED[key]
        r["pmid"], r["doi"] = pmid, doi
        n += 1
        print(f"  {key[0]} {key[1]}  {r['chemical'][:20]:20s} "
              f"{r['protein'][:34]:34s} -> {pmid} / {doi}")

    if args.dry_run:
        print(f"\n(dry run) {n} rows would be filled")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nfilled {n} rows in {os.path.abspath(CSV)}")

    blank = sum(1 for r in rows if not (r["pmid"].strip() or r["doi"].strip()))
    print(f"still without an identifier: {blank} rows "
          f"(Messina 2005, Gao 2019, Li 2009, Kerstner-Wood 2003, EFSA 2018)")


if __name__ == "__main__":
    main()
