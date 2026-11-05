#!/usr/bin/env python3
"""The eight sources PubMed has no record of, resolved by hand.

PubMed indexes 63 of this database's 71 studies.  The other eight are
reports, a 1954 Faraday Discussion, and journals PubMed does not
index, so each was resolved from a named authority rather than from
memory, and the authority is recorded beside it:

  OEHHA's own reference list, from the PHG text held in papers/ --
  this is how the four sources that reach this database only through
  OEHHA Table A6.2 are cited, so the citation matches what OEHHA
  used.  Messina, Klevens & Ellenbogen, Li, and EFSA.
  Crossref, by DOI, where the database already held one.
  EPA HERO record 4771364 for the Kerstner-Wood report, which gives
  the study ID, the authoring organisation and the TSCA submission
  numbers that OEHHA's table does not.

Run:  python3 tk_review/scripts/fill_missing_citations.py
Rewrites db/citations.csv in place; idempotent.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CITES = os.path.join(HERE, "..", "db", "citations.csv")

FILLED = {
 ("EFSA", "2018"): dict(
    authors="EFSA Panel on Contaminants in the Food Chain (CONTAM)",
    title="Risk to human health related to the presence of "
          "perfluorooctane sulfonic acid and perfluorooctanoic acid "
          "in food",
    journal="EFSA J", pub_year="2018", volume="16", issue="12",
    pages="e05194", doi="10.2903/j.efsa.2018.5194",
    authority="OEHHA 2024 PHG reference list"),
 ("Fischer et al.", "2025"): dict(
    authors="Fischer FC; Thackray C; Ferguson N; Chicoine C; "
            "Skende O; Hu Z; Zhu Y; Slitt A; Sunderland EM",
    title="Understanding Mechanisms of PFAS Absorption, "
          "Distribution, and Elimination Using a Physiologically "
          "Based Toxicokinetic Model",
    journal="Environ Sci Technol", pub_year="2025", volume="59",
    issue="26", pages="13240-13250",
    doi="10.1021/acs.est.5c05473", authority="Crossref, by DOI"),
 ("Kerstner-Wood et al.", "2003"): dict(
    authors="Kerstner-Wood C; Coward L; Gorman G",
    title="Protein binding of perfluorobutane sulfonate, "
          "perfluorohexane sulfonate, perfluorooctane sulfonate and "
          "perfluorooctanoate to plasma (human, rat, and monkey), "
          "and various human-derived plasma protein fractions. "
          "Study ID 9921.7",
    journal="Southern Research Institute for 3M Company, St. Paul "
            "MN; UNPUBLISHED report, US EPA TSCA submission "
            "8EHQ-04-15845A / 88040000364",
    pub_year="2003", volume="", issue="", pages="", doi="",
    authority="EPA HERO record 4771364"),
 ("Klevens & Ellenbogen", "1954"): dict(
    authors="Klevens HB; Ellenbogen E",
    title="Protein-fluoroacid interaction. Bovine serum "
          "albumin-perfluoro-octanoic acid",
    journal="Discuss Faraday Soc", pub_year="1954", volume="18",
    issue="", pages="277-288", doi="10.1039/df9541800277",
    authority="OEHHA 2024 PHG reference list; Crossref, by DOI"),
 ("Li et al.", "2009"): dict(
    authors="Li L; Xu ZS; Song GW",
    title="Study on the Langmuir aggregation of fluorinated "
          "surfactants on protein",
    journal="J Fluor Chem", pub_year="2009", volume="130", issue="2",
    pages="225-230", doi="10.1016/j.jfluchem.2008.10.007",
    authority="OEHHA 2024 PHG reference list; Crossref, by title"),
 ("Messina et al.", "2005"): dict(
    authors="Messina P; Prieto G; Dodero V; Ruso JM; Schulz P; "
            "Sarmiento F",
    title="Ultraviolet-circular dichroism spectroscopy and "
          "potentiometric study of the interaction between human "
          "serum albumin and sodium perfluorooctanoate",
    journal="Biopolymers", pub_year="2005", volume="79", issue="6",
    pages="300-309", doi="10.1002/bip.20353", pmid="16113999",
    authority="OEHHA 2024 PHG reference list, confirmed against "
              "PubMed 16113999"),
 ("Peng et al.", "2025"): dict(
    authors="Peng S; Yang Y; Tian R; Lu N",
    title="Critical new insights into the interactions of "
          "hexafluoropropylene oxide-dimer acid (GenX or HFPO-DA) "
          "with albumin at molecular and cellular levels",
    journal="J Environ Sci", pub_year="2025", volume="149", issue="",
    pages="88-98", doi="10.1016/j.jes.2024.02.011",
    authority="Crossref, by DOI"),
 ("Wu et al.", "2009"): dict(
    authors="Wu L; Gao H; Gao N; Chen F; Chen L",
    title="Interaction of perfluorooctanoic acid with human serum "
          "albumin",
    journal="BMC Struct Biol", pub_year="2009", volume="9", issue="",
    pages="31", doi="10.1186/1472-6807-9-31",
    authority="Crossref, by DOI"),
}


def main():
    with open(CITES, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fields, rows = list(reader.fieldnames), list(reader)
    if "authority" not in fields:
        fields.append("authority")

    n = 0
    for r in rows:
        f = FILLED.get((r["study"], r["year"]))
        if not f or r["citation"]:
            r.setdefault("authority", "PubMed" if r["citation"] else "")
            continue
        for k, v in f.items():
            if k != "authority":
                r[k] = v or r.get(k, "")
        v = f["volume"] + (f"({f['issue']})" if f["issue"] else "")
        tail = f"{f['journal']} {f['pub_year']}" + (f";{v}" if v else "")
        tail += f":{f['pages']}" if f["pages"] else ""
        r["citation"] = f"{f['authors']}. {f['title']}. {tail}."
        r["authority"] = f["authority"]
        n += 1
        print(f"  filled  {r['study']} {r['year']}")
        print(f"          via {f['authority']}")

    with open(CITES, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    done = sum(1 for r in rows if r["citation"])
    print(f"\n{n} filled this run; {done}/{len(rows)} studies now have "
          f"a full citation")


if __name__ == "__main__":
    main()
