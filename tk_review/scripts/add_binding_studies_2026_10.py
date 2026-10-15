"""
Append three binding studies found 2026-10-15 to db/protein_binding.csv.

Only one of the three adds experimental binding measurements, and none of the
three adds a NUMBER, for reasons recorded per row rather than smoothed over:

  Gao 2025        Experimental, BSA, 10 PFAS resolved by chain length and head
                  group. Not in PMC and the full text was not obtained, so the
                  per-compound affinities are not here. The rows carry the
                  design and the qualitative direction, with ka left empty.

  Gallagher 2023  QSAR + read-across + docking. Its underlying affinities are
                  "collected from the literature" -- 24 PFAS, EC50 in mM, 4:2
                  and 6:2 FTOH non-binding, H/L cut at 1.45 mM. That is the
                  Starnes dataset, which db/protein_binding.csv already holds
                  (33 rows, Kd in mM, DSF, 6:2 FTOH nonbinding). So it is
                  logged as COMPUTATIONAL over existing data, not as new
                  measurements -- adding its values would double-count Starnes.

  Zhao 2023       Review. No primary data. One row, flagged as such, kept for
                  the citation trail.

Run:  python3 scripts/add_binding_studies_2026_10.py
"""

import csv
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "db", "protein_binding.csv")

GAO_NOTE = ("abstract only, full text not obtained (not in PMC): per-compound "
            "binding affinities NOT extracted. Qualitative findings: affinity "
            "depends on chain length and head group; weakens with increasing "
            "ionic strength and with decreasing cation valency; pH-dependent, "
            "salt effects stronger at lower pH; thermodynamics indicate a "
            "hydrophobic-dominant interaction. PFAS structure did NOT change "
            "how ionic strength affected affinity.")

# (chemical, n perfluorinated carbons, head group)
GAO = [("PFBA", "3", "carboxylate"), ("PFHxA", "5", "carboxylate"),
       ("PFHpA", "6", "carboxylate"), ("PFOA", "7", "carboxylate"),
       ("PFNA", "8", "carboxylate"), ("PFDA", "9", "carboxylate"),
       ("PFBS", "4", "sulfonate"), ("PFHxS", "6", "sulfonate"),
       ("PFOS", "8", "sulfonate"), ("GenX (HFPO-DA)", "4", "ether carboxylate")]


def main():
    with open(OUT, newline="") as fh:
        fields = next(csv.reader(fh))
    rows = []

    for chem, n, hg in GAO:
        rows.append(dict(
            study="Gao et al.", year="2025", pmid="41344664",
            doi="10.1021/acs.est.5c08779", chemical=chem,
            n_perfluorinated_carbons=n, head_group=hg,
            protein="bovine serum albumin (BSA)", protein_species="bovine",
            method="not stated in abstract; binding affinity vs salt "
                   "composition, ionic strength and pH",
            is_computational="no",
            tk_relevance_noted="no - framed as environmental fate and "
                               "wastewater partitioning, not toxicokinetics",
            source_table=GAO_NOTE))

    rows.append(dict(
        study="Gallagher et al.", year="2023", pmid="37513249",
        doi="10.3390/molecules28145375",
        chemical="24 PFAS (PFSA C4-C8, PFCA C4-C12, mono- and polyether "
                 "perfluoroalkyl ether acids, fluorotelomers)",
        n_perfluorinated_carbons="3-12", head_group="mixed",
        protein="human serum albumin (HSA)", protein_species="human",
        method="QSAR (classification + regression), GenRA read-across, and "
               "Glide XP docking onto PDB 7AAI",
        binding_site_identified="high-affinity PFAS form an H-bond with Arg348 "
                                "AND salt bridges to Arg348/Arg485; "
                                "low-affinity PFAS form neither",
        is_computational="yes",
        tk_relevance_noted="indirect - binding affinity as a bioaccumulation "
                           "proxy, no kinetic parameter",
        source_table="NO NEW MEASUREMENTS. Underlying EC50 values (mM, "
                     "4-parameter fit, 24 PFAS, 4:2 and 6:2 FTOH non-binding, "
                     "H/L threshold 1.45 mM) are collected from the literature "
                     "and match the Starnes dataset already in this file; "
                     "logged to avoid double-counting it. Useful finding: PFAA "
                     "with C7-C9 bind HSA most tightly and >C10 bind less, "
                     "which corroborates the non-monotonic chain-length "
                     "relationship seen independently in Fischer 2024."))

    rows.append(dict(
        study="Zhao et al.", year="2023", pmid="37126916",
        doi="10.1016/j.envint.2023.107951",
        chemical="legacy PFAS and alternatives (review)",
        head_group="mixed",
        protein="review across serum albumin, liver FABP, organic anion "
                "transporters, nuclear receptors, cell and mitochondrial "
                "membranes",
        protein_species="multiple species compared (interspecies variation)",
        method="literature review; no primary measurements",
        is_computational="yes",
        tk_relevance_noted="yes - reviews binding as a determinant of "
                           "distribution and accumulation",
        source_table="REVIEW, no primary data. Mechanistic summary: the "
                     "hydrophobic fluorinated chain occupies the protein "
                     "binding cavity while the acid head group hydrogen-bonds "
                     "to residues; chlorine substitution, ether oxygen "
                     "insertion and branching change chain length and "
                     "hydrophobicity and so change affinity. Recommends "
                     "integrating in vitro work with in silico QSAR."))

    with open(OUT, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        for r in rows:
            w.writerow(r)
    print(f"appended {len(rows)} rows to {os.path.relpath(OUT, HERE)}")
    print(f"  Gao 2025        {len(GAO)} experimental BSA rows, ka empty "
          f"(full text not obtained)")
    print("  Gallagher 2023   1 row, computational, flagged as Starnes-derived")
    print("  Zhao 2023        1 row, review, no primary data")


if __name__ == "__main__":
    main()
