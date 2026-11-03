#!/usr/bin/env python3
"""Zhao et al. 2023: L-FABP docking + machine-learning QSAR for 1599 PFAS.

Zhao J, Shi X, Wang Z, Xiong S, Lin Y, Wei X, Li Y, Tang X.
"Hepatotoxicity assessment investigations on PFASs targeting L-FABP
using binding affinity data and machine learning-based QSAR model."
Ecotoxicol Environ Saf 2023;262:115310.
doi 10.1016/j.ecoenv.2023.115310.  Open access, CC BY-NC-ND.
Read from the PDF.

THIS PAPER CONTAINS NO MEASURED BINDING CONSTANTS, and that is the
reason to be careful with it rather than a reason to leave it out.  Its
"binding affinity data" are molecular docking scores in kcal/mol,
computed with Glide (Schrodinger) against PDB 2LKK, for a library of
1599 PFAS taken from Liu et al. 2019.  Every row here is therefore
is_computational = yes and carries no kd or ka, so the albumin and
L-FABP extracts, which drop computational rows, will not pick it up.

It is a DIFFERENT PAPER from the Zhao et al. 2023 already in this
database, which is an Environment International review
(doi 10.1016/j.envint.2023.107951) with no primary data.  Both are
kept; the doi distinguishes them.

WHAT IT IS GOOD FOR:

  * The only QSAR-on-binding source in the collection.  Five regression
    algorithms (Ridge, Lasso, decision tree, SVM, random forest) over
    13 molecular fingerprints and 192 2D descriptors, 4:1 train/verify
    split, ten-fold cross-validation, plus Bayesian Kernel Machine
    Regression for joint descriptor effects.
  * Its only tie to experiment: a Pearson r = 0.79 between the docking
    energies and the Kd values of Yang et al. 2020, measured by in
    vitro fluorescence displacement.  That is the paper's whole
    validation, and r = 0.79 on a log-scale quantity is a weak anchor.
    Yang et al. 2020 is NOT yet in this database and is the measurement
    this paper rests on.
  * It resolves L-FABP's two sites as PORTAL and NON-PORTAL, which is
    the same two-site architecture Sheng et al. 2016 calls outer and
    inner.  The nomenclature differs between the two papers.

A CONTRADICTION IN THE PAPER, recorded rather than resolved.  Section
3.1 states the portal-site affinities "ranged from -8.92 to -2.44
kcal/mol with a mean of -6.33 kcal/mol, and from -9.42 to -4.56
kcal/mol with a mean of -5.38 kcal/mol to the non-portal site", then
concludes "PFASs showed a higher affinity for the non-portal site than
the portal site.  The difference of about 1 kcal/mol".  As printed the
portal mean is the more negative of the two, so the means and the
conclusion disagree about which site is preferred.  Both are stored
verbatim; neither is silently corrected.

Run:  python3 scripts/add_zhao2023_lfabp_qsar.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Zhao et al.", "2023"
DOI = "10.1016/j.ecoenv.2023.115310"
PROTEIN = ("liver fatty acid-binding protein (L-FABP), PDB 2LKK "
           "receptor model")

DOCK = ("molecular docking, Glide (Schrodinger Maestro), OPLS4 force "
        "field, receptor relaxed at 310 K, grid 25 x 25 x 25 A centred "
        "on oleic acid, flexible induced-fit redocking; ionizable "
        "states set at pH 7.0 +/- 0.2")
QSAR = ("machine-learning QSAR on the docking scores: Ridge, Lasso, "
        "decision tree, SVM and random forest over 13 molecular "
        "fingerprints and 192 2D descriptors, 4:1 train/verify split, "
        "ten-fold cross-validation and external validation; plus "
        "Bayesian Kernel Machine Regression (R 4.2.1, 25,000 MCMC "
        "iterations) for joint descriptor effects")

ROWS = [
    ("1599 PFAS library (Liu et al. 2019), portal site", DOCK,
     "docking affinity, portal region: 766 of 1599 PFAS docked "
     "successfully | range -8.92 to -2.44 kcal/mol | mean -6.33 "
     "kcal/mol | natural ligand oleic acid -4.34 kcal/mol"),
    ("1599 PFAS library (Liu et al. 2019), non-portal site", DOCK,
     "docking affinity, non-portal region: 653 of 1599 PFAS docked "
     "successfully, 621 in both regions | range -9.42 to -4.56 "
     "kcal/mol | mean -5.38 kcal/mol | natural ligand oleic acid "
     "-3.46 kcal/mol | the text concludes PFAS prefer the non-portal "
     "site by about 1 kcal/mol, about six-fold in IC50 at the median, "
     "which CONTRADICTS the means as printed -- see the script "
     "docstring"),
    ("1599 PFAS library (Liu et al. 2019), docking validation", DOCK,
     "the paper's only tie to measurement: Pearson r = 0.79 between "
     "docking energy and the Kd values of Yang et al. 2020 (in vitro "
     "fluorescence displacement). Yang et al. 2020 is not yet in this "
     "database and is the experimental anchor this paper rests on"),
    ("1599 PFAS library (Liu et al. 2019), QSAR model", QSAR,
     "predicted binding affinity to L-FABP from structure; molecular "
     "flexibility identified by BKMR as the determinant of "
     "PFAS-induced hepatotoxicity, acting both on its own and through "
     "its interaction with other descriptors. Class comparison (Fig. "
     "3): polyfluoroalkyl sulfates bind the portal region more "
     "strongly than perfluoroalkyl carboxylates (p < 0.05), and "
     "perfluoroalkyl sulfonamides the non-portal region (p < 0.01)"),
]

TK = ("indirect - no measured constant. The paper's premise is that "
      "L-FABP binding drives hepatic accumulation, citing a positive "
      "correlation between -log Kd for L-FABP and PFAS accumulation in "
      "liver (Cao 2019, Jia 2022). Useful to this review as the only "
      "QSAR-on-binding source, not as a data point")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []
    for chem, method, note in ROWS:
        r = blank(fields)
        r.update(
            study=STUDY, year=YEAR, pmid="", doi=DOI, chemical=chem,
            protein=PROTEIN, protein_species="human", method=method,
            ka="", ka_units="", kd="",
            kd_units="NOT a dissociation constant: docking score in "
                     "kcal/mol" if "kcal" in note else "",
            n_binding_sites="2 (portal and non-portal; the same "
                            "two-site architecture Sheng et al. 2016 "
                            "calls outer and inner)",
            binding_site_identified="portal and non-portal regions; "
                                    "Lys31/Ser56 in the portal site, "
                                    "Arg122/Ser39/Ser124 in the "
                                    "non-portal site",
            is_computational="yes", tk_relevance_noted=TK,
            source_table=note)
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, existing = reader.fieldnames, list(reader)

    before = len(existing)
    # Keyed on the DOI, NOT on study+year: the other Zhao et al. 2023
    # in this database is a different paper and must survive.
    existing = [r for r in existing if r["doi"] != DOI]
    print(f"removed {before - len(existing)} existing rows for {DOI}")
    other = [r for r in existing
             if r["study"] == STUDY and r["year"] == YEAR]
    print(f"kept    {len(other)} other '{STUDY} {YEAR}' rows "
          f"(doi {other[0]['doi'] if other else '-'})")

    new = rows(fields)
    print(f"adding  {len(new)} rows, all is_computational=yes "
          f"-> {len(existing) + len(new)} total")
    for r in new:
        print(f"   {r['chemical'].split(', ')[-1]:<22}"
              f"{r['source_table'][:62]}")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(existing + new)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    main()
