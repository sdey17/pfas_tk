#!/usr/bin/env python3
"""How much of each binding protein a tissue actually contains.

db/protein_binding.csv holds affinities. An affinity ranks SITES. It
cannot rank COMPARTMENTS, because the mass a tissue holds is

    bound per litre of phase = C_protein x n x theta,   theta -> Cw/Kd

so the compartment comparison needs C_protein, which no binding paper
reports. This file supplies it, with provenance, and nothing else.

ONE ROW PER MEASURED ABUNDANCE. Every row says what was measured, in what
units the source reported it, what conversion was applied, and what the
source assumed to make that conversion. Nothing is averaged across
sources.

Writes db/protein_abundance.csv.

Run:  python3 tk_review/scripts/add_protein_abundance.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "db", "protein_abundance.csv")

FIELDS = ["protein", "species", "tissue", "basis", "value_uM", "low_uM",
          "high_uM", "cv_percent", "n_subjects", "method", "as_reported",
          "conversion", "study", "year", "doi", "pmid", "peer_reviewed",
          "location", "note"]

ROWS = [
    # ---------------------------------------------------------- FABP1
    dict(protein="liver fatty acid-binding protein (L-FABP / FABP1)",
         species="human", tissue="liver",
         basis="per litre of intracellular (cytosolic) water",
         value_uM="399", low_uM="211", high_uM="760", cv_percent="24.8",
         n_subjects="61",
         method="surrogate-peptide LC-MS/MS (signature peptide FTITAGSK) on "
                "human liver S9 from the University of Washington human "
                "liver bank; each liver measured in singlet on three "
                "separate days",
         as_reported="nmol FABP1 per mg S9 protein, extrapolated to liver "
                     "expression using the HLS9 yield",
         conversion="the SOURCE's own: tissue expression divided by the "
                    "liver fraction of intracellular water, fICW = 0.573, "
                    "with tissue density taken as 1 mL/g. Not recomputed "
                    "here",
         study="Wen & Isoherranen", year="2026",
         doi="10.1016/j.jpet.2026.104309", pmid="", peer_reviewed="yes",
         location="Results, prediction of individual human liver Kp "
                  "values; Figure 5C and Table S3-1",
         note="The largest direct quantification of human hepatic FABP1 "
              "found. Independent corroboration of the order: McCarty & "
              "Guengerich 2025 (J Biol Chem 301:108168) state FABP1 is "
              "7-11% of liver cytosolic protein, reaching up to 1 mM, "
              "which brackets this mean and maximum. A preprint of the "
              "same work is at doi 10.64898/2026.01.17.700130; the "
              "peer-reviewed JPET version is cited here"),
    # ---------------------------------------------------------- albumin
    dict(protein="serum albumin",
         species="human", tissue="plasma",
         basis="per litre of plasma",
         value_uM="", low_uM="527", high_uM="752", cv_percent="",
         n_subjects="",
         method="clinical reference interval, not a single study",
         as_reported="35-50 g/L",
         conversion="divided by 66,463 g/mol, the albumin molar mass used "
                    "throughout this project",
         study="clinical reference range", year="",
         doi="", pmid="", peer_reviewed="n/a",
         location="standard clinical chemistry reference interval",
         note="DELIBERATELY A RANGE, NOT A POINT. No single study is "
              "cited because none is needed: the comparison this feeds is "
              "reported across the whole interval so the conclusion can "
              "be read against it. A point estimate here would be false "
              "precision in the one input that varies most between people"),
]


def main():
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in ROWS:
            w.writerow({f: r.get(f, "") for f in FIELDS})
    print(f"wrote {os.path.abspath(OUT)}  ({len(ROWS)} rows)")
    for r in ROWS:
        rng = (f"{r['low_uM']}-{r['high_uM']} uM" if r["low_uM"] else "")
        print(f"  {r['protein'][:44]:46s}{r['species']:7s}{r['tissue']:8s}"
              f"{r['value_uM'] or '-':>6s} uM   {rng}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
