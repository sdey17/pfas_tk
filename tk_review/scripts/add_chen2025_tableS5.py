#!/usr/bin/env python3
"""Chen et al. 2025 SI Table S5: K_D and log K_MW for 60 PFAS.

Equilibrium dialysis on Slide-A-Lyzer devices, HSA held at 1.0 uM, PFAS at
0.0625-2 uM, so the ligand:protein molar ratio runs 0.06-2.  That is far
closer to the physiological regime than most of this literature, and it makes
these the most directly usable dissociation constants in the file.

Membrane-water partition coefficients (log K_MW, solid-supported lipid
membrane, pH 7.0) are carried on the same rows as context, because the
paper's point is the comparison between them.

Replaces the single study-level row written by add_binding_qin2024_chen2025.py.

Run:  python3 tk_review/scripts/add_chen2025_tableS5.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

# (PFAS, log K_MW +/- SE, K_D uM +/- SE, nCF, head group)
T5 = [
    ("PFBA", "1.63 +/- 0.193", 84.9, 15.95, 3, "carboxylate"),
    ("PFPeA", "2.02 +/- 0.149", 16.27, 3.58, 4, "carboxylate"),
    ("PFHxA", "2.42 +/- 0.156", 13.98, 2.79, 5, "carboxylate"),
    ("PFHpA", "2.85 +/- 0.115", 6.2, 1.13, 6, "carboxylate"),
    ("PFOA", "3.28 +/- 0.135", 2.57, 0.45, 7, "carboxylate"),
    ("PFNA", "3.75 +/- 0.163", 3.75, 0.77, 8, "carboxylate"),
    ("PFDA", "4.18 +/- 0.213", 3.64, 0.85, 9, "carboxylate"),
    ("PFUnDA", "4.5 +/- 0.215", 4.73, 1.08, 10, "carboxylate"),
    ("PFDoDA", "4.82 +/- 0.1", 15.76, 3.34, 11, "carboxylate"),
    ("PFTrDA", "5.11 +/- 0.187", 31.44, 5.09, 12, "carboxylate"),
    ("PFTeDA", "5.36 +/- 0.241", 13.09, 3.01, 13, "carboxylate"),
    ("PFHxDA", "5.73 +/- 0.25", 10.45, 1.71, 15, "carboxylate"),
    ("PFBS", "2.72 +/- 0.308", 15.6, 2.68, 4, "sulfonate"),
    ("PFPeS", "3.16 +/- 0.127", 17.15, 2.58, 5, "sulfonate"),
    ("PFHxS", "3.65 +/- 0.115", 1.31, 0.25, 6, "sulfonate"),
    ("PFHpS", "4.09 +/- 0.174", 1.68, 0.36, 7, "sulfonate"),
    ("PFOS", "4.5 +/- 0.145", 0.94, 0.23, 8, "sulfonate"),
    ("PFNS", "4.82 +/- 0.179", 6.56, 1.36, 9, "sulfonate"),
    ("PFDS", "5.15 +/- 0.224", 1.88, 0.4, 10, "sulfonate"),
    ("PFDoS", "5.66 +/- 0.121", 1.29, 0.21, 12, "sulfonate"),
    ("4:2 FTS", "2.21 +/- 0.132", 18.96, 3.59, 4, "sulfonate (fluorotelomer)"),
    ("6:2 FTS", "3.12 +/- 0.135", 2.42, 0.6, 6, "sulfonate (fluorotelomer)"),
    ("8:2 FTS", "4.07 +/- 0.264", 4.49, 0.84, 8, "sulfonate (fluorotelomer)"),
    ("10:2 FTS", "4.57 +/- 0.071", 2.31, 0.51, 10, "sulfonate (fluorotelomer)"),
    ("PFOSA", "3.91 +/- 0.162", 11.18, 1.92, 8, "sulfonamide"),
    ("NMeFOSA", "4.17 +/- 0.172", 12.22, 1.94, 8, "sulfonamide (N-alkyl)"),
    ("NEtFOSA", "4.41 +/- 0.23", 12.98, 2.06, 8, "sulfonamide (N-alkyl)"),
    ("N-MeFOSAA", "4.21 +/- 0.179", 6.63, 1.49, 8, "sulfonamidoacetic acid"),
    ("N-EtFOSAA", "4.4 +/- 0.153", 5.44, 1.28, 8, "sulfonamidoacetic acid"),
    ("FOSAA", "4.28 +/- 0.257", 7.4, 1.69, 8, "sulfonamidoacetic acid"),
    ("NMeFOSE", "3.96 +/- 0.12", 5.88, 0.97, 8, "sulfonamidoethanol"),
    ("NEtFOSE", "4.21 +/- 0.089", 6.09, 1.39, 8, "sulfonamidoethanol"),
    ("HFPO-DA (GenX)", "2.07 +/- 0.129", 6.57, 1.49, 5, "ether carboxylate"),
    ("ADONA", "2.07", 5.68, 0.89, 6, "ether carboxylate"),
    ("PFMPA", "2 +/- 0.223", 47.55, 7.91, 3, "ether carboxylate"),
    ("PFMBA", "2.04 +/- 0.144", 0.22, 0.04, 4, "ether carboxylate"),
    ("NFDHA", "2", 17.34, 3.01, 4, "ether carboxylate"),
    ("9Cl-PF3ONS", "4.38 +/- 0.191", 1.58, 0.39, 8, "ether sulfonate (chlorinated)"),
    ("11Cl-PF3OUdS", "4.75 +/- 0.263", 3.52, 0.77, 10, "ether sulfonate (chlorinated)"),
    ("PFEESA", "2.49 +/- 0.246", 5.15, 1.22, 4, "ether sulfonate"),
    ("3:3 FTCA", "2.14 +/- 0.135", 70.57, 14.94, 3, "carboxylate (fluorotelomer)"),
    ("5:3 FTCA", "2.47 +/- 0.1", 10.75, 2.0, 5, "carboxylate (fluorotelomer)"),
    ("7:3 FTCA", "3.4 +/- 0.208", 1.82, 0.3, 7, "carboxylate (fluorotelomer)"),
    ("6:2 FTCA", "2.34 +/- 0.103", 3.64, 0.74, 6, "carboxylate (fluorotelomer)"),
    ("8:2 FTCA", "3.23 +/- 0.203", 1.99, 0.37, 8, "carboxylate (fluorotelomer)"),
    ("10:2 FTCA", "4.1 +/- 0.083", 9.06, 1.61, 10, "carboxylate (fluorotelomer)"),
    ("8Cl-PFOS", "4.73 +/- 0.198", 0.87, 0.17, 8, "sulfonate (chlorinated)"),
    ("PFEtCHxS", "3.81 +/- 0.156", 1.62, 0.34, 8, "sulfonate (cyclic)"),
    ("FBSA", "2.07 +/- 0.243", 5.58, 0.92, 4, "sulfonamide"),
    ("FHxSA", "3 +/- 0.203", 3.07, 0.57, 6, "sulfonamide"),
    ("6:2 FTUCA", "2.16 +/- 0.023", 4.23, 0.89, 6, "carboxylate (unsaturated fluorotelomer)"),
    ("8:2 FTUCA", "3.13 +/- 0.196", 2.26, 0.37, 8, "carboxylate (unsaturated fluorotelomer)"),
    ("6:2 diPAP", "4.62 +/- 0.18", 2.71, 0.41, 6, "phosphate diester"),
    ("8:2 diPAP", "5.17 +/- 0.15", 14.4, 3.11, 8, "phosphate diester"),
    ("diSAmPAP", "5.06 +/- 0.072", 5.94, 1.1, 8, "phosphate diester (sulfonamido)"),
    ("PFHxSaAm", "3.52", 6.19, 1.01, 6, "polyfluoroalkyl amine"),
    ("6:2 FtSaB", "2.98 +/- 0.2", 16.97, 3.52, 6, "betaine"),
    ("N-TAmP-FHxSA", "3.14", 33.68, 6.0, 6, "polyfluoroalkyl amine"),
    ("5:3 FTB", "1.63 +/- 0.1", 30.44, 5.78, 5, "betaine"),
    ("5:1:2 FTB", "2.01 +/- 0.1", 49.32, 10.1, 5, "betaine"),
]

TK = (
    "yes - K_D measured at a ligand:protein molar ratio of 0.06-2 (HSA fixed "
    "at 1.0 uM, PFAS 0.0625-2 uM), the closest to physiological of any "
    "wide-coverage dataset here. Paired log K_MW lets the albumin and "
    "phospholipid sinks be compared directly: protein-water distribution "
    "exceeds membrane-water partitioning at low PFAS concentration and the "
    "balance shifts to the membrane as concentration rises"
)


def rows():
    out = []
    for (chem, kmw, kd, kdse, ncf, head) in T5:
        ka = 1.0 / (kd * 1e-6)
        note = ""
        if chem == "PFMBA":
            note = (" | NOTE: 0.22 uM is the tightest K_D in the table by 4x "
                    "and sits far off the chain-length trend for a 4-carbon "
                    "ether acid; recorded as printed")
        out.append(dict(
            study="Chen et al.", year="2025", pmid="", doi="10.1021/acs.est.4c06734",
            chemical=chem, n_perfluorinated_carbons=str(ncf), head_group=head,
            protein="human serum albumin (HSA)", protein_species="human",
            method=("equilibrium dialysis, Slide-A-Lyzer 10 kDa MWCO, 48 h at "
                    "50 rpm; single-site Langmuir fit to the specific binding "
                    "curve; LC-MS/MS"),
            ka=f"{ka:.4g}", ka_units="M^-1 (computed as 1/K_D; the paper reports K_D)",
            kd=f"{kd:g}", kd_units="uM",
            n_binding_sites="B_max fitted per compound (one site assumed)",
            binding_site_identified="",
            unbound_fraction="",
            protein_conc="1.0", protein_conc_units="uM",
            ligand_conc_range="0.0625-2",
            ligand_conc_units="uM (0.0625, 0.125, 0.25, 0.5, 1, 2)",
            temperature_c="room temperature (not stated more precisely)",
            buffer="10 mM ammonium acetate, pH 7.4",
            is_computational="no",
            tk_relevance_noted=TK,
            source_table=(f"SI Table S5 | K_D {kd} +/- {kdse} uM | "
                          f"log K_MW (phospholipid membrane-water, pH 7.0) {kmw}"
                          + note),
            assay_volume_note=(
                "1000 uL PFAS solution in the vial, 100 uL of 1.0 mM HSA in the "
                "Slide-A-Lyzer vessel; ligand:protein molar ratio 0.06-2, which "
                "is 1-2 orders of magnitude closer to physiological (nu < 0.01) "
                "than the DSF and fluorescence datasets in this file"),
            ligand_over_protein_max="2",
        ))
    return out


def main():
    with open(CSV) as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames
        existing = list(rdr)
    before = len(existing)
    existing = [r for r in existing
                if (r["study"], r["year"]) != ("Chen et al.", "2025")]
    new = rows()
    for r in new:
        assert set(r) <= set(fields), set(r) - set(fields)
        existing.append({f: r.get(f, "") for f in fields})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(existing)
    kds = [t[2] for t in T5]
    print(f"Chen 2025 Table S5: {len(new)} compounds")
    print(f"  K_D range {min(kds)} - {max(kds)} uM "
          f"(tightest {min(T5, key=lambda t: t[2])[0]}, "
          f"weakest {max(T5, key=lambda t: t[2])[0]})")
    print(f"{before} -> {len(existing)} rows")


if __name__ == "__main__":
    main()
