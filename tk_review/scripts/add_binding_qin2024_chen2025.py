#!/usr/bin/env python3
"""Add Qin et al. 2024 (ES&T 58:9954) and Chen et al. 2025 (ES&T 59:82).

Qin 2024 is the same group as the Qin 2023 mass-balance model and measures,
by solid-phase microextraction over four orders of magnitude in concentration,
exactly the quantity that model's cell term is built from: the SPECIFIC and
NONSPECIFIC protein-water distribution ratios of 16 PFAS, in BSA, human plasma
and trout plasma.  Separating the two is the point -- it is the direct
measurement of the ligand:protein regime effect that sections 6.1a and 6.1b
could only infer.

Chen 2025 derives membrane-water and protein-water partition coefficients for
60 PFAS and, in its Eq. 4, writes the same Kd -> partition coefficient bridge
derived independently in section 6.1b.

Run:  python3 tk_review/scripts/add_binding_qin2024_chen2025.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

# Qin 2024 Table 2: log D (pH 7.4), volume basis (L_water per L_protein).
# (chem, nCF, head, logD_lipw, BSA_spec, BSA_nonspec, HP_spec, HP_nonspec,
#  TP_spec, TP_nonspec)   "" = not resolvable;  * values are concentration
#  averages with n_Freundlich fixed to 1.
T2 = [
    ("PFBA", 3, "carboxylate", "1.00",
     "2.44 (2.10-2.63)", "1.97 (1.90-2.04)", "2.20 (0.87-2.49)",
     "1.35 (1.23-1.44)", "", "1.43 (1.28-1.58)*"),
    ("PFHxA", 5, "carboxylate", "2.32",
     "3.26 (2.98-3.42)", "2.69 (2.66-2.72)", "3.19 (3.01-3.31)",
     "2.61 (2.56-2.66)", "", "2.49 (2.41-2.57)*"),
    ("PFHpA", 6, "carboxylate", "2.91",
     "4.12 (3.87-4.27)", "3.47 (3.42-3.52)", "4.10 (3.98-4.20)",
     "3.20 (3.15-3.24)", "", "3.32 (3.25-3.39)*"),
    ("PFOA", 7, "carboxylate", "3.52",
     "4.58 (4.47-4.67)", "3.88 (3.86-3.90)", "4.45 (4.34-4.54)",
     "3.65 (3.59-3.70)", "", "4.18 (4.14-4.21)*"),
    ("PFNA", 8, "carboxylate", "4.25",
     "", "4.61 (4.56-4.65)*", "", "4.19 (4.14-4.24)*", "", "4.71 (4.65-4.76)*"),
    ("PFUnDA", 10, "carboxylate", "4.54",
     "", "4.75 (4.71-4.80)*", "", "4.54 (4.49-4.59)*", "", "4.99 (4.96-5.03)*"),
    ("HFPO-DA (GenX)", 5, "ether carboxylate", "2.41",
     "3.31 (2.95-3.50)", "2.17 (2.05-2.26)", "3.35 (2.77-3.58)",
     "1.99 (1.93-2.04)", "3.56 (3.08-3.78)", "1.89 (1.67-2.04)"),
    ("PFHxS", 6, "sulfonate", "4.13",
     "5.02 (4.90-5.11)", "3.58 (3.50-3.64)", "4.98 (4.80-5.11)",
     "2.92 (2.78-3.03)", "", "4.00 (3.91-4.10)*"),
    ("PFOS", 8, "sulfonate", "4.89",
     "5.27 (4.84-5.48)", "4.17 (4.14-4.20)", "4.82 (4.41-5.03)",
     "4.07 (4.04-4.11)", "5.49 (4.88-5.74)", "4.57 (4.54-4.60)"),
    ("6:2 FTSA", 6, "sulfonate (fluorotelomer)", "3.87",
     "", "3.87 (3.84-3.90)*", "", "3.67 (3.63-3.70)*", "", "4.11 (4.04-4.19)*"),
    ("PFOSA", 8, "sulfonamide", "4.94",
     "", "4.32 (4.28-4.37)*", "", "3.90 (3.86-3.93)*", "", "4.33 (4.27-4.40)*"),
    ("hexaflumuron", 2, "benzoylurea pesticide (98% neutral at pH 7.4)", "4.58",
     "", "3.96 (3.93-4.00)*", "", "4.29 (4.24-4.35)*", "", "4.61 (4.53-4.69)*"),
    ("flubendiamide", 3, "phthalamide pesticide (98% neutral at pH 7.4)", "3.28",
     "", "3.98 (3.93-4.02)*", "", "3.88 (3.85-3.91)*", "", "4.46 (4.42-4.50)*"),
    ("6:2 FTOH", 6, "alcohol (fluorotelomer, neutral)", "4.38",
     "", "2.67 (2.54-2.77)*", "", "2.73 (2.49-2.89)*", "", "2.77 (2.49-2.94)*"),
    ("8:2 FTOH", 8, "alcohol (fluorotelomer, neutral)", "5.85",
     "", "4.61 (4.48-4.72)*", "", "4.55 (4.49-4.61)*", "", "5.01 (4.90-5.10)*"),
    ("10:2 FTOH", 10, "alcohol (fluorotelomer, neutral)", "7.73",
     "", "6.72 (6.44-6.89)*", "", "6.21 (6.11-6.28)*", "", "6.86 (6.68-6.99)*"),
]

MATRIX = [
    (4, "specific", "bovine serum albumin (BSA)", "bovine",
     "BSA 0.1-10 mg/mL in PBS pH 7.4, 37 C"),
    (5, "nonspecific", "bovine serum albumin (BSA)", "bovine",
     "BSA 0.1-10 mg/mL in PBS pH 7.4, 37 C"),
    (6, "specific", "human plasma (protein + lipid)", "human",
     "human plasma diluted in PBS to match the BSA protein level, pH 7.4, 37 C"),
    (7, "nonspecific", "human plasma (protein + lipid)", "human",
     "human plasma diluted in PBS to match the BSA protein level, pH 7.4, 37 C"),
    (8, "specific", "trout plasma (protein + lipid)", "rainbow trout",
     "trout plasma diluted in PBS to match the BSA protein level, pH 7.4, 25 C"),
    (9, "nonspecific", "trout plasma (protein + lipid)", "rainbow trout",
     "trout plasma diluted in PBS to match the BSA protein level, pH 7.4, 25 C"),
]

TK_Q = (
    "yes - D_protein/water offered directly as a PBTK input for free PFAS in "
    "plasma. The paper's central result for this review: at physiological "
    "PFAS:protein molar ratios (nu < 0.01) SPECIFIC binding dominates, while "
    "in vitro assays dose in excess (nu > 1) where NONSPECIFIC binding takes "
    "over - the two differ by up to 10x for the same compound and protein. "
    "Also: BSA is a good surrogate for human plasma, and the human/trout "
    "difference is explained by lipid:protein content, not by protein identity"
)


def qin2024_rows():
    out = []
    for row in T2:
        chem, nCF, head, lipw = row[0], row[1], row[2], row[3]
        for idx, kind, prot, sp, buf in MATRIX:
            val = row[idx]
            if not val:
                continue
            avg = val.endswith("*")
            point = val.split(" ")[0]
            ci = val[val.find("(") + 1:val.find(")")] if "(" in val else ""
            out.append(dict(
                study="Qin et al.", year="2024", pmid="38804966",
                doi="10.1021/acs.est.3c10824",
                chemical=chem, n_perfluorinated_carbons=str(nCF), head_group=head,
                protein=prot, protein_species=sp,
                method=(
                    "solid-phase microextraction (BioSPME 96-pin, C18) + LC-MS/MS, "
                    "binding isotherm over ~4 orders of magnitude; specific and "
                    "nonspecific terms separated by a combined "
                    "binding/partitioning model"
                    + ("; HS-SPME + GC-MS for the volatile FTOHs"
                       if "FTOH" in chem else "")),
                ka=point,
                ka_units=(
                    f"log D_protein/water at pH 7.4, {kind} binding "
                    f"(L_water per L_protein); 95% CI {ci}"
                    + ("; concentration-average with n_Freundlich fixed to 1, "
                       "i.e. the two terms could not be separated" if avg else "")),
                kd="", kd_units="", n_binding_sites="",
                binding_site_identified=(
                    "one saturable site assumed for the specific term"
                    if kind == "specific" else ""),
                unbound_fraction="",
                protein_conc=("0.1-10 mg/mL BSA" if "albumin" in prot
                              else "plasma diluted to the matching protein level"),
                protein_conc_units="",
                ligand_conc_range="9-step 2-fold dilution series spanning ~4 orders of magnitude",
                ligand_conc_units="uM",
                temperature_c="25" if "trout" in prot else "37",
                buffer=buf,
                is_computational="no",
                tk_relevance_noted=TK_Q,
                source_table=(
                    f"Table 2, {kind} column, {prot.split(' (')[0]}"
                    f" | log D_lip/w {lipw} (membrane lipid, from Droge/Ebert, "
                    "not measured here)"),
                assay_volume_note=(
                    "Plasma composition measured in the same paper (Table 3): "
                    "human 42.25 mL/L protein and 4.46 mL/L lipid, trout 15.46 "
                    "and 7.08. Protein density taken as 1.36 kg/L, lipid 1.00"),
                ligand_over_protein_max="",
            ))
    return out


CHEN = dict(
    study="Chen et al.", year="2025", pmid="", doi="10.1021/acs.est.4c06734",
    chemical="60 structurally diverse PFAS across 11 structural groups",
    n_perfluorinated_carbons="3-16",
    head_group="carboxylate, sulfonate, sulfonamide, sulfonamidoethanol, ether acid, fluorotelomer, phosphate diester, betaine",
    protein="human serum albumin (HSA)", protein_species="human",
    method=("equilibrium dialysis (Slide-A-Lyzer, 10 kDa MWCO, 48 h) at PFAS:HSA "
            "molar ratios 1:16, 1:8, 1:4, 1:2, 1:1 and 2:1, with AutoDock Vina "
            "docking alongside; membrane-water partitioning by solid-supported "
            "lipid membrane (TRANSIL)"),
    ka="", ka_units="log K_A (M^-1); per-congener values are in the figures and SI, not extracted",
    kd="", kd_units="",
    n_binding_sites="B_max fitted per compound (single-site Langmuir)",
    binding_site_identified=("docking places the highest-affinity PFAS in HSA fatty-acid "
                             "sites 1, 2 and 6; PFOA and PFOS in FA site 6"),
    unbound_fraction="",
    protein_conc="10", protein_conc_units="mg/mL (MW 66.4 kDa)",
    ligand_conc_range="PFAS:HSA molar ratio 1:16 to 2:1",
    ligand_conc_units="mol/mol",
    temperature_c="not stated for the dialysis step",
    buffer="ammonium acetate / ammonium hydroxide, pH 6 and pH 7 for K_MW; 10 mM ammonium acetate for SSLM",
    is_computational="no (experimental; the AutoDock Vina arm is computational)",
    tk_relevance_noted=(
        "yes - the paper's Eq. 4 writes the Kd-to-partition-coefficient bridge "
        "independently of this review's section 6.1b: D_HSA/W = B_max * rho_HSA "
        "/ [(K_D + C_aq/MW_PFAS) * MW_HSA], which is the same Langmuir "
        "expression with the albumin concentration cancelled. Headline numbers: "
        "K_MW rises 0.36 +/- 0.01 log units per CF2 for PFCAs (FC3-FC16) and "
        "0.37 +/- 0.02 for PFSAs; HSA affinity peaks between FC6 and FC10; "
        "subclass order PFSAs ~ ESAs > PFCAs > FTUCAs ~ diPAPs ~ PFECAs. "
        "Protein-water distribution EXCEEDS membrane-water partitioning at low "
        "PFAS concentrations and the balance shifts toward the membrane as "
        "concentration rises"),
    source_table="Figures 1-4 and Table 1; per-congener K_MW and K_A are in the SI, which was not retrieved",
    assay_volume_note=("1.0 mL PFAS solution in 1.5 mL cups, 100 uL into the dialysis "
                       "cup; SSLM lipid volumes 0.067-2.166 uL; methanol rinse step "
                       "added for the most hydrophobic compounds to recover "
                       "container-surface sorption"),
    ligand_over_protein_max="2",
)


def main():
    with open(CSV) as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames
        rows = list(rdr)
    before = len(rows)
    drop = {("Qin et al.", "2024"), ("Chen et al.", "2025")}
    rows = [r for r in rows if (r["study"], r["year"]) not in drop]

    new = qin2024_rows() + [CHEN]
    for r in new:
        unknown = set(r) - set(fields)
        if unknown:
            raise SystemExit(f"unknown fields: {unknown}")
        rows.append({f: r.get(f, "") for f in fields})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Qin 2024 Table 2 : {len(qin2024_rows())} rows "
          f"({len(T2)} PFAS x BSA / human plasma / trout plasma, "
          f"specific and nonspecific)")
    print(f"Chen 2025        : 1 study-level row (per-congener values need the SI)")
    print(f"{before} -> {len(rows)} rows")


if __name__ == "__main__":
    main()
