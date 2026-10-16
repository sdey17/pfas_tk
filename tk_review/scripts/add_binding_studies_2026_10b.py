#!/usr/bin/env python3
"""Ingest five primary binding papers supplied 2026-10-16 into db/protein_binding.csv.

  Chen & Guo 2009  Arch Toxicol 83:255-261   site-resolved Ka for 5 PFAA on HSA
  Gao et al. 2025  ES&T 59:26716             Ka_Hill for 10 PFAS on BSA vs added NaCl
  Han et al. 2021  ES&T 55:9012              QITSA dTm, hL-FABP panel + 4 species
  Ryu et al. 2024  ES&T 58:14641             f_unbound in 5 whole tissues, 4 organisms
  Starnes 2025     Toxicol Sci 203:67        DSF Ka/dTm for 65 PFAS + 3 fatty acids

Every row replaced here previously held an abstract-only or secondary-source
placeholder.  Rows are keyed on (study, year, protein, chemical, source_table)
so a re-run is idempotent.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")
# Starnes 2025 SI Table S4, extracted from the 135 MB supplementary workbook
# (openpyxl with explicit unmerge; LibreOffice will not open it).
TABLE_S4 = (sys.argv[1] if len(sys.argv) > 1
            else os.path.join(HERE, "..", "db", "starnes2025_tableS4.csv"))

# --------------------------------------------------------------------------
# 1. Chen & Guo 2009 -- three site-specific fluorescence probes, 15 cells
# --------------------------------------------------------------------------
CG_COMMON = dict(
    study="Chen & Guo",
    year="2009",
    pmid="19002450",
    doi="10.1007/s00204-008-0359-x",
    protein="human serum albumin (HSA)",
    protein_species="human",
    ka_units="M^-1",
    protein_conc="0.1",
    protein_conc_units="mM",
    ligand_conc_units="mM",
    temperature_c="not stated (ambient; Perkin-Elmer LS55, no thermostat reported)",
    buffer="100 mM NaH2PO4/Na2HPO4, pH 7.4; HSA fraction V, essentially fatty-acid-free",
    is_computational="no",
    tk_relevance_noted=(
        "yes - authors argue that at ng/mL human serum levels PFAAs are too dilute to "
        "displace fatty acids, but could displace weaker or scarcer endogenous ligands; "
        "binding framed as lowering free PFAA and so reducing bioavailability"
    ),
    source_table="Table 1",
    assay_volume_note=(
        "cuvette fluorescence titration, 1 cm path; Trp214 excited at 295 nm, probes at "
        "350 nm (dansylamide) / 375 nm (dansyl-L-proline); 15 nm slits both sides"
    ),
)

# (chemical, n_pfc, head_group, site, ka, n_sites, ligand_range, method_tail, note)
CG_ROWS = [
    # ---- tryptophan (Trp214) site, intrinsic fluorescence quenching -------
    ("PFBA", 3, "carboxylate", "Trp214", "", "", "0-6",
     "fluorescence quenching of Trp-214 (intrinsic, lambda_ex 295 nm)",
     "no binding - PFBA induced no change in tryptophan fluorescence"),
    ("PFBS", 4, "sulfonate", "Trp214", "", "", "0-6",
     "fluorescence quenching of Trp-214 (intrinsic, lambda_ex 295 nm)",
     "no binding - PFBS induced no change in tryptophan fluorescence"),
    ("PFOA", 7, "carboxylate", "Trp214", "2.7e5", "5.8", "0-0.8",
     "fluorescence quenching of Trp-214 (intrinsic, lambda_ex 295 nm) + Scatchard plot",
     "Ka 2.7e5 +/- 0.24e5 M^-1, binding ratio n = 5.8; emission blue-shifted 341->330 nm, "
     "plateau above 0.8 mM. This is the constant OEHHA 2024 Table A6.2 lists for PFOA "
     "(1/2.7e5 = 3.7e-6 M) and it is the WEAKEST of the three sites in this paper"),
    ("PFOS", 8, "sulfonate", "Trp214", "2.2e4", "14", "0-3",
     "fluorescence quenching of Trp-214 (intrinsic, lambda_ex 295 nm) + Scatchard plot",
     "Ka 2.2e4 +/- 0.43e4 M^-1, binding ratio n = 14; emission blue-shifted 341->331 nm, "
     "plateau above 3 mM. This is the constant OEHHA 2024 Table A6.2 lists for PFOS "
     "(1/2.2e4 = 4.5e-5 M) and it is the WEAKEST of the three sites in this paper"),
    ("PFDoA", 11, "carboxylate", "Trp214", "", "", "0-2",
     "fluorescence quenching of Trp-214 (intrinsic, lambda_ex 295 nm)",
     "binds but constant not obtained - solubility-limited at 2 mM so no plateau reached; "
     "from the titration slope, weaker than PFOS at this site"),
    # ---- Sudlow site I, dansylamide displacement -------------------------
    ("PFBA", 3, "carboxylate", "Sudlow site I", "1.1e6", "", "0-6",
     "competitive displacement of dansylamide (site I probe), Ka from IC50 via Cheng-Prusoff",
     "Ka 1.1e6 +/- 0.1e6 M^-1. 25x tighter than the Trp-site constant OEHHA used for PFOA; "
     "PFBA binds HSA at site I only"),
    ("PFBS", 4, "sulfonate", "Sudlow site I", "2.2e6", "", "0-6",
     "competitive displacement of dansylamide (site I probe), Ka from IC50 via Cheng-Prusoff",
     "Ka 2.2e6 +/- 0.1e6 M^-1"),
    ("PFOA", 7, "carboxylate", "Sudlow site I", "", "", "0-6",
     "competitive displacement of dansylamide (site I probe)",
     "binds but constant not obtained - multiphasic fluorescence transition"),
    ("PFOS", 8, "sulfonate", "Sudlow site I", "", "", "0-4",
     "competitive displacement of dansylamide (site I probe)",
     "binds but constant not obtained; signal rose rather than fell, interpreted as a "
     "conformational change raising probe fluorescence"),
    ("PFDoA", 11, "carboxylate", "Sudlow site I", "", "", "0-2",
     "competitive displacement of dansylamide (site I probe)",
     "no binding at site I - signal rose then levelled without displacing the probe"),
    # ---- Sudlow site II, dansyl-L-proline displacement -------------------
    ("PFBA", 3, "carboxylate", "Sudlow site II", "", "", "0-6",
     "competitive displacement of dansyl-L-proline (site II probe)",
     "no binding - DP fluorescence essentially unchanged with PFBA"),
    ("PFBS", 4, "sulfonate", "Sudlow site II", "6.5e6", "", "0-6",
     "competitive displacement of dansyl-L-proline (site II probe), Ka from IC50 via Cheng-Prusoff",
     "Ka 6.5e6 +/- 0.7e6 M^-1 - the tightest constant in the paper after PFOS site II"),
    ("PFOA", 7, "carboxylate", "Sudlow site II", "", "", "0-6",
     "competitive displacement of dansyl-L-proline (site II probe)",
     "binds but constant not obtained - multiphasic, like the site I curve"),
    ("PFOS", 8, "sulfonate", "Sudlow site II", "7.6e6", "", "0-6",
     "competitive displacement of dansyl-L-proline (site II probe), Ka from IC50 via Cheng-Prusoff",
     "Ka 7.6e6 +/- 0.3e6 M^-1 for the first (high-affinity) of two phases, i.e. Kd 0.13 uM. "
     "345x tighter than the 2.2e4 M^-1 Trp-site value OEHHA 2024 Table A6.2 lists for PFOS; "
     "the site choice, not the chemistry, is what makes PFOS look like a weak binder"),
    ("PFDoA", 11, "carboxylate", "Sudlow site II", "1.2e6", "", "0-2",
     "competitive displacement of dansyl-L-proline (site II probe), Ka from IC50 via Cheng-Prusoff",
     "Ka 1.2e6 +/- 0.04e6 M^-1"),
]

# Fatty-acid reference constants the same authors report in the discussion for
# the same three probes.  Flagged as unpublished-at-the-time.
CG_FA_ROWS = [
    ("butyric acid (C4)", "Sudlow site I", "1.8e5"),
    ("butyric acid (C4)", "Trp214", "9.2e3"),
    ("octanoic acid (C8)", "Sudlow site I", "3.1e6"),
    ("octanoic acid (C8)", "Sudlow site II", "1.1e7"),
    ("octanoic acid (C8)", "Trp214", "7.0e4"),
]


def chen_guo_rows():
    out = []
    for (chem, npfc, head, site, ka, n, lrange, method, note) in CG_ROWS:
        r = dict(CG_COMMON)
        r.update(
            chemical=chem,
            n_perfluorinated_carbons=str(npfc),
            head_group=head,
            method=method,
            ka=ka,
            kd="" if not ka else f"{1.0 / float(ka) * 1e6:.4g}",
            kd_units="" if not ka else "uM",
            n_binding_sites=n,
            binding_site_identified=site,
            unbound_fraction="",
            ligand_conc_range=lrange,
            # 0.1 mM HSA, ligand up to lrange high end
            ligand_over_protein_max=f"{float(lrange.split('-')[1]) / 0.1:.0f}",
        )
        r["source_table"] = "Table 1" + ("" if ka else " (qualitative cell)")
        r["assay_volume_note"] = CG_COMMON["assay_volume_note"] + " | " + note
        out.append(r)
    for (chem, site, ka) in CG_FA_ROWS:
        r = dict(CG_COMMON)
        r.update(
            chemical=chem,
            n_perfluorinated_carbons="0",
            head_group="carboxylate (non-fluorinated fatty acid; reference compound)",
            method="same three site-specific fluorescence probes as the PFAA series",
            ka=ka,
            kd=f"{1.0 / float(ka) * 1e6:.4g}",
            kd_units="uM",
            n_binding_sites="",
            binding_site_identified=site,
            unbound_fraction="",
            ligand_conc_range="not stated",
            ligand_over_protein_max="",
            source_table="Discussion text, cited as 'Chen and Guo, unpublished work'",
        )
        r["assay_volume_note"] = (
            "REFERENCE COMPOUND, not a PFAS. Included because the displacement argument "
            "turns on whether PFAS bind as tightly as the fatty acids they would have to "
            "compete with: at site II octanoic acid (1.1e7 M^-1) is tighter than every "
            "PFAS measured here. Value is from the discussion, attributed to unpublished "
            "work by the same authors, not from Table 1"
        )
        out.append(r)
    return out


# --------------------------------------------------------------------------
# 2. Gao et al. 2025 -- Ka_Hill vs added NaCl, Table 1
# --------------------------------------------------------------------------
# chemical -> n_pfc, head_group
GAO_CHEM = {
    "PFBA": (3, "carboxylate"),
    "PFHxA": (5, "carboxylate"),
    "PFHpA": (6, "carboxylate"),
    "PFOA": (7, "carboxylate"),
    "PFNA": (8, "carboxylate"),
    "PFDA": (9, "carboxylate"),
    "PFBS": (4, "sulfonate"),
    "PFHxS": (6, "sulfonate"),
    "PFOS": (8, "sulfonate"),
    "HFPO-DA (GenX)": (5, "ether carboxylate"),
}
# (chemical, added NaCl mM, Ka_Hill x1e-4 L/mol, SEM, Kas1 x1e-4, Kas_app x1e-4, R2)
GAO_T1 = [
    ("PFBA", 0, 0.03, 0.01, "0.04 +/- 0.003", 0.03, 0.86),
    ("PFBA", 150, 0.02, 0.001, "", 0.02, ""),
    ("PFHxA", 0, 0.22, 0.03, "0.26 +/- 0.01", 0.20, 0.83),
    ("PFHxA", 10, 0.18, 0.02, "", 0.19, ""),
    ("PFHxA", 30, 0.17, 0.01, "", 0.18, ""),
    ("PFHxA", 80, 0.13, 0.01, "", 0.14, ""),
    ("PFHpA", 0, 2.43, 0.19, "2.73 +/- 0.16", 2.13, 0.94),
    ("PFHpA", 10, 1.96, 0.14, "", 2.03, ""),
    ("PFHpA", 30, 1.92, 0.18, "", 1.85, ""),
    ("PFHpA", 80, 1.45, 0.10, "", 1.52, ""),
    ("PFHpA", 150, 1.22, 0.14, "", 1.22, ""),
    ("PFHpA", 300, 0.88, 0.07, "", 0.85, ""),
    ("PFHpA", 500, 0.71, 0.04, "", 0.61, ""),
    ("PFOA", 0, 30.94, 3.22, "32.68 +/- 2.97", 25.47, 0.83),
    ("PFOA", 10, 23.22, 2.21, "", 24.26, ""),
    ("PFOA", 30, 20.37, 1.24, "", 22.15, ""),
    ("PFOA", 80, 14.11, 0.66, "", 18.21, ""),
    ("PFNA", 0, 99.42, 2.09, "114.62 +/- 5.61", 89.33, 0.81),
    ("PFNA", 10, 84.92, 3.31, "", 85.08, ""),
    ("PFNA", 30, 71.73, 2.55, "", 77.70, ""),
    ("PFNA", 80, 57.20, 2.26, "", 63.85, ""),
    ("PFDA", 0, 105.21, 6.43, "131.10 +/- 2.19", 102.17, 0.95),
    ("PFDA", 10, 93.10, 7.10, "", 97.32, ""),
    ("PFDA", 30, 89.55, 4.69, "", 88.87, ""),
    ("PFDA", 80, 73.57, 8.25, "", 73.03, ""),
    ("PFBS", 0, 0.14, 0.03, "0.25 +/- 0.001", 0.15, 0.99),
    ("PFBS", 150, 0.09, 0.01, "", 0.09, ""),
    ("PFHxS", 0, 19.10, 1.31, "19.80 +/- 1.867", 15.43, 0.67),
    ("PFHxS", 10, 13.19, 0.73, "", 14.70, ""),
    ("PFHxS", 30, 12.29, 0.41, "", 13.42, ""),
    ("PFHxS", 80, 9.29, 0.10, "", 11.03, ""),
    ("PFOS", 0, 143.44, 10.92, "237.19 +/- 7.41", 143.30, 0.90),
    ("PFOS", 10, 144.22, 9.46, "", 136.49, ""),
    ("PFOS", 30, 110.49, 8.81, "", 124.65, ""),
    ("PFOS", 80, 93.47, 5.40, "", 102.42, ""),
    ("PFOS", 150, 92.25, 11.22, "", 81.97, ""),
    ("PFOS", 300, 58.27, 3.97, "", 57.40, ""),
    ("PFOS", 500, 58.40, 2.77, "", 41.01, ""),
    ("HFPO-DA (GenX)", 0, 0.24, 0.01, "0.30 +/- 0.01", 0.23, 0.84),
    ("HFPO-DA (GenX)", 10, 0.22, 0.02, "", 0.22, ""),
    ("HFPO-DA (GenX)", 30, 0.21, 0.02, "", 0.20, ""),
    ("HFPO-DA (GenX)", 80, 0.14, 0.02, "", 0.16, ""),
]


def gao_rows():
    out = []
    for (chem, nacl, ka4, sem4, kas1, kasapp4, r2) in GAO_T1:
        npfc, head = GAO_CHEM[chem]
        ka = ka4 * 1e4
        kd_uM = 1.0 / ka * 1e6
        phys = " (physiological ionic strength)" if nacl == 150 else ""
        out.append(dict(
            study="Gao et al.",
            year="2025",
            pmid="41344664",
            doi="10.1021/acs.est.5c08779",
            chemical=chem,
            n_perfluorinated_carbons=str(npfc),
            head_group=head,
            protein="bovine serum albumin (BSA)",
            protein_species="bovine",
            method=(
                "intrinsic Trp fluorescence quenching (lambda_ex 280 nm, em 300-500 nm, "
                "5 nm slits), Ka_Hill from nonlinear least-squares fit of the Hill equation"
            ),
            ka=f"{ka:.4g}",
            ka_units="M^-1 (Ka_Hill; equals Ka only if n_Hill = 1)",
            kd=f"{kd_uM:.4g}",
            kd_units="uM (1/Ka_Hill)",
            n_binding_sites=(
                "n_Hill 0.69-2.5 overall; ~1 site for PFHxA, PFHpA, PFOA, PFHxS and GenX, "
                "~2 for PFNA and PFDA; 1.5-2 for PFBA and the PFSAs"
            ),
            binding_site_identified="",
            unbound_fraction="",
            protein_conc="0.8",
            protein_conc_units="uM",
            ligand_conc_range="titration to fluorescence saturation; per-compound series in Figures S5-S12",
            ligand_conc_units="uM",
            temperature_c="20",
            buffer=(
                f"20 mM sodium phosphate, pH 7, plus {nacl} mM added NaCl{phys}"
            ),
            is_computational="no (experimental; a competitive-Langmuir salt model is fitted to the data)",
            tk_relevance_noted=(
                "no - framed as environmental fate and wastewater partitioning. Relevant here "
                "as the only chain-length series measured across ionic strength: affinity falls "
                "2-3.4x from 0 to 500 mM NaCl, so constants measured in low-salt buffer "
                "overstate binding at the 150 mM of plasma"
            ),
            source_table=(
                f"Table 1 (Kas1 = {kas1 if kas1 else 'n/a'} x1e4 L/mol; "
                f"model Kas_app = {kasapp4} x1e4 L/mol; R2 = {r2 if r2 != '' else 'see 0 mM row'})"
            ),
            assay_volume_note=(
                "1 cm quartz cuvette; BSA stock made fresh daily. Note pH 7 and 20 C, not "
                "pH 7.4 and 37 C: binding affinity rose with temperature in this study, so "
                "20 C values are a lower bound for body temperature"
            ),
            ligand_over_protein_max="",
        ))
    return out


# --------------------------------------------------------------------------
# 3. Han et al. 2021 -- QITSA thermal shift, L-FABP
# --------------------------------------------------------------------------
HAN_COMMON = dict(
    study="Han et al.",
    year="2021",
    pmid="34133149",
    doi="10.1021/acs.est.1c00509",
    ka="",
    ka_units="",
    protein_conc="endogenous L-FABP in cell lysate or liver tissue extract; not quantified in molar terms",
    protein_conc_units="",
    ligand_conc_range="100",
    ligand_conc_units="uM",
    temperature_c="thermal ramp; L-FABP Tm 73-76 C in the absence of PFOS",
    buffer="cell lysis / tissue extraction buffer, PFAS added from DMSO; vehicle control = DMSO",
    is_computational="no (AutoDock4 binding energies are reported alongside, as a separate computational strand)",
    unbound_fraction="",
    n_binding_sites="",
    binding_site_identified="",
)

# (chemical, n_pfc, head, kd_uM, dTm, sd, note)
HAN_PANEL = [
    ("PFBA", 3, "carboxylate", ">1000", "-0.67", "0.68", ""),
    ("PFBS", 4, "sulfonate", "1034", "0.35", "0.66", ""),
    ("PFHxS", 6, "sulfonate", "85.7", "4.43", "0.74",
     "the paper's running text labels this third dTm 'PFHxA', but Figure 3A lists the panel "
     "as PFBA, PFOA, PFNA, PFBS, PFHxS and PFOS, and 4.43 C is the value that sits "
     "monotonically between PFBS (Kd 1034 uM) and PFOA (Kd 50.4 uM) on the Kd-dTm line. "
     "Recorded as PFHxS"),
    ("PFOA", 7, "carboxylate", "50.4", "5.42", "0.43", ""),
    ("PFOS", 8, "sulfonate", "18.5", "6.37", "0.67", ""),
    ("PFNA", 8, "carboxylate", "16.2", "8.04", "1.81", ""),
]

# (species label, protein label, dTm pre-dialysis, sd, dTm post-dialysis, sd, matrix)
HAN_SPECIES = [
    ("human", "liver fatty acid-binding protein (hL-FABP)", "7.32", "1.81", "14.71", "2.09",
     "endogenous protein in HepG2 cell lysate"),
    ("mouse", "liver fatty acid-binding protein (mL-FABP)", "2.01", "0.61", "2.32", "1.21",
     "endogenous protein in mouse liver tissue extract"),
    ("rat", "liver fatty acid-binding protein (rL-FABP)", "1.90", "0.73", "1.92", "0.89",
     "endogenous protein in rat liver tissue extract"),
    ("zebrafish", "liver fatty acid-binding protein (zL-FABP)", "-0.92", "0.62", "-1.00", "0.90",
     "endogenous protein in zebrafish liver tissue extract"),
]

HAN_RECOMB = [
    ("rat", "recombinant liver fatty acid-binding protein (rL-FABP)", "3.91", "1.32"),
    ("mouse", "recombinant liver fatty acid-binding protein (mL-FABP)", "5.61", "0.56"),
]


def han_rows():
    out = []
    for (chem, npfc, head, kd, dtm, sd, extra) in HAN_PANEL:
        r = dict(HAN_COMMON)
        note = (
            f"QITSA dTm = {dtm} +/- {sd} C at 100 uM. Kd is NOT measured here: it is quoted "
            "from Zhang et al. (ref 34), a fluorescence-displacement assay on recombinant "
            "hL-FABP, and used as the independent x-axis against which dTm is calibrated "
            "(log Kd vs dTm, R2 = 0.98)"
        )
        if extra:
            note += ". " + extra
        r.update(
            chemical=chem,
            n_perfluorinated_carbons=str(npfc),
            head_group=head,
            protein="liver fatty acid-binding protein (hL-FABP), endogenous in HepG2",
            protein_species="human",
            method=(
                "quantitative interspecies thermal shift assay (QITSA): PFAS-induced dTm of "
                "endogenous protein measured by targeted LC-MS/MS proteomics after thermal "
                "denaturation and centrifugation (n = 2)"
            ),
            kd=kd,
            kd_units="uM (from Zhang et al. ref 34, NOT measured in this paper)",
            tk_relevance_noted=(
                "yes - the paper's thesis is that L-FABP binding differences track interspecies "
                "differences in PFOS elimination half-life. Note Ryu et al. 2024 reaches the "
                "opposite conclusion from whole-tissue f_unbound (liver binding no stronger than "
                "other tissues despite L-FABP enrichment)"
            ),
            source_table="Results text, Figures 3B and 4B",
            assay_volume_note=note,
        )
        out.append(r)
    for (sp, prot, pre, presd, post, postsd, matrix) in HAN_SPECIES:
        r = dict(HAN_COMMON)
        r.update(
            chemical="PFOS",
            n_perfluorinated_carbons="8",
            head_group="sulfonate",
            protein=prot,
            protein_species=sp,
            method=(
                "QITSA on liver tissue extract / cell lysate, after dialysis to strip "
                "endogenous ligands; shotgun then PRM-targeted proteomics for L-FABP"
            ),
            kd="",
            kd_units="",
            tk_relevance_noted=(
                "yes - ranked thermal enhancement human > mouse > rat > zebrafish is offered as "
                "a partial explanation of interspecies PFOS half-life differences, alongside "
                "organic anion transporters"
            ),
            source_table="Results text, Figure 5A",
            assay_volume_note=(
                f"{matrix}. dTm = {post} +/- {postsd} C after dialysis; {pre} +/- {presd} C "
                "before dialysis. Dialysis matters because endogenous fatty acids occupy the "
                "site: for human L-FABP it doubles the apparent shift. Order is the same either "
                "way. Interspecies spread is 7.7x (human vs mouse) on dTm, which is a thermal "
                "shift and not a Kd"
            ),
        )
        out.append(r)
    for (sp, prot, dtm, sd) in HAN_RECOMB:
        r = dict(HAN_COMMON)
        r.update(
            chemical="PFOS",
            n_perfluorinated_carbons="8",
            head_group="sulfonate",
            protein=prot,
            protein_species=sp,
            method="QITSA on recombinant protein expressed in E. coli (validation arm)",
            kd="",
            kd_units="",
            protein_conc="recombinant protein, 1 ug aliquots",
            tk_relevance_noted=(
                "yes - recombinant arm shows the weak rodent response is intrinsic to the "
                "protein, not an artefact of endogenous ligand occupancy"
            ),
            source_table="Results text, Figure 5C",
            assay_volume_note=(
                f"dTm = {dtm} +/- {sd} C, above the tissue-lysate value for the same species but "
                "still far below human L-FABP (14.71 C). zL-FABP was not available as recombinant"
            ),
        )
        out.append(r)
    r = dict(HAN_COMMON)
    r.update(
        chemical="palmitic acid (C16:0)",
        n_perfluorinated_carbons="0",
        head_group="carboxylate (endogenous fatty acid; reference compound)",
        protein="liver fatty acid-binding protein (hL-FABP), endogenous in HepG2",
        protein_species="human",
        method="QITSA, same protocol as the PFAS panel",
        kd="0.63",
        kd_units="uM (from Zhang et al. ref 34, NOT measured in this paper)",
        tk_relevance_noted=(
            "yes - the endogenous ligand binds hL-FABP 28x tighter than PFOS, so PFAS do not "
            "outcompete fatty acids on this protein"
        ),
        source_table="Results text, Figure 4C",
        assay_volume_note=(
            "dTm = 2.86 +/- 0.44 C, 2.2x SMALLER than PFOS despite a 28x tighter Kd. The "
            "authors attribute this to a different unfolding path (two sequential palmitate "
            "bindings), which is the paper's own warning that dTm is not a monotone proxy for "
            "affinity across binding modes"
        ),
    )
    out.append(r)
    return out


# --------------------------------------------------------------------------
# 4. Ryu et al. 2024 -- rewrite the single row
# --------------------------------------------------------------------------
def ryu_row():
    return dict(
        study="Ryu et al.",
        year="2024",
        pmid="39161261",
        doi="10.1021/acs.est.4c04050",
        chemical=(
            "16 PFAS: PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFNA, PFDA, PFUnDA, PFDoDA, PFTriDA, "
            "PFTeDA, PFBS, PFHxS, PFOS, PFOSA, 6:2 FtS"
        ),
        n_perfluorinated_carbons="3-13",
        head_group="11 carboxylates, 3 sulfonates, 1 sulfonamide, 1 fluorotelomer sulfonate",
        protein=(
            "whole tissue homogenate (brain, heart, kidney, liver, lung) - total protein plus "
            "phospholipid binding, not an isolated protein. NO plasma or serum in this study"
        ),
        protein_species=(
            "human (pooled donors; brain = 2M + 2F from NIH NeuroBioBank), Wistar Han rat, "
            "CD-1 mouse, C57BL/6J mouse; rodent tissues pooled by sex"
        ),
        method=(
            "rapid equilibrium dialysis (HTD96, 12-14 kDa cutoff cellulose), 18 h at 37 C, "
            "200 rpm, 5% CO2, 75% RH, n = 4; LC-MS/MS quantification; diluted f_unbound "
            "corrected to undiluted with the standard dilution-factor equation"
        ),
        ka="",
        ka_units="",
        kd="",
        kd_units="",
        n_binding_sites="",
        binding_site_identified="",
        unbound_fraction=(
            "6.8e-6 to 0.3 across the 16 PFAS x 5 tissues x 4 organisms; 6 of 16 PFAS "
            ">99.999% bound. Falls with molecular weight and log D. K_tissue/w highest in "
            "brain for 14/16 PFAS; lung lowest, especially for C10-C13 carboxylates. "
            "Per-congener values are in the Supporting Information, NOT in the main text"
        ),
        protein_conc="tissue homogenate diluted 1:5 in DPBS, or 1:10 for brain",
        protein_conc_units="dilution factor",
        ligand_conc_range="5",
        ligand_conc_units="uM (single concentration, 1% DMSO; chosen to stay below albumin-type saturation)",
        temperature_c="37",
        buffer="Dulbecco's PBS, pH 7.4",
        is_computational="no (experimental; PLSR models on the data are a separate strand)",
        tk_relevance_noted=(
            "yes - f_unbound and K_tissue/w are offered directly as PBTK inputs. Three findings "
            "bear on species extrapolation: (1) log f_unbound correlates across all species and "
            "tissue pairs, R2 >= 0.97, Spearman 0.96-1, so tissue binding is NOT a species "
            "difference; (2) rat liver is the best single surrogate for human tissues, R2 >= 0.98, "
            "with the scaling equations in Table 1; (3) liver binding was within 2-fold of other "
            "tissues for 14/16 PFAS despite L-FABP enrichment, which argues L-FABP matters less "
            "for retention than Han et al. 2021 proposes. The paper's own half-life claim is "
            "weaker than it reads: PFOS f_unbound is 5.8-11.3x below PFOA, but the human "
            "half-lives it cites differ only 1.4x (3.8 vs 5.4 y)"
        ),
        source_table=(
            "Table 1 gives rat-liver-to-tissue scaling of the form "
            "f_unbound = 1/(1 + k*(1/f_unbound,rat,liver - 1)), with k = 1.78 (human brain), "
            "1.68 (rat brain), 1.64 (CD-1 mouse brain), 1.75 (C57BL/6 mouse brain); "
            "0.462 / 0.589 / 0.761 for kidney; 0.468 / 0.359 / 0.416 for lung; heart and liver "
            "take f_unbound,rat,liver unchanged. Per-congener f_unbound: SI (not extracted)"
        ),
        assay_volume_note=(
            "150 uL spiked homogenate donor / 150 uL DPBS receiver; membranes presaturated "
            "twice; matrix-matched samples built for LC-MS/MS; tolbutamide and carbamazepine "
            "as internal standards"
        ),
        ligand_over_protein_max="",
    )


# --------------------------------------------------------------------------
# 5. Starnes et al. 2025 -- Table S4, 68 chemicals
# --------------------------------------------------------------------------
# class -> (n_pfc placeholder rule, head group text)
S25_HEAD = {
    "Fatty acids": "carboxylate (non-fluorinated fatty acid; reference compound)",
    "PFCA": "carboxylate",
    "PFSA": "sulfonate",
    "PFECA": "ether carboxylate",
    "PFESA": "ether sulfonate",
    "Polyfluorinated carboxylic acids": "carboxylate (H-substituted, not fully fluorinated)",
    "FTCA": "carboxylate (fluorotelomer)",
    "FTSA": "sulfonate (fluorotelomer)",
    "FTOH": "alcohol (fluorotelomer)",
    "PFAE": "ether (neutral)",
    "Polyfluoroalkyl amines": "amine",
    "Perfluoroalkyl amides": "amide",
    "N-alkyl FASAs and FASEs": "sulfonamide / sulfonamidoethanol",
    "Per and polyfluoroalkyl trifluoromethane-sulfonates": "trifluoromethanesulfonate ester",
    "Per and polyfluoroalkyl methyl esters": "methyl ester",
    "Per and polyfluoroalkyl diols": "diol",
    "Trifluoro-crotonic acids and crotonates": "unsaturated carboxylate / ester",
    "FTACs and MFTACs": "acrylate / methacrylate (fluorotelomer)",
    "Polyfluoralkyl tetrols": "tetrol",
    "Polyfluoroalkyl aldehydes": "aldehyde",
    "Hydrofluoroolefins": "olefin (neutral)",
}
# Abbreviations in Table S4 are normalised to the convention already used in
# protein_binding.csv, so the new rows join to the existing ones.
S25_RENAME = {
    "PFunDA": "PFUnDA",
    "HFPO-DA, GenX": "HFPO-DA (GenX)",
    "EtFOSA, Sulfluramid, Finatron": "EtFOSA",
    "Nafion Byproduct 2": "Nafion byproduct 2",
}


def starnes2025_rows(path, npfc_lookup):
    raw = list(csv.reader(open(path)))[3:]
    out = []
    for r in raw:
        if len(r) < 9:
            continue
        cls, cas, name, abbr, bound, ka, r2, dtm = (x.strip() for x in r[1:9])
        if not name or cls.startswith("*"):
            continue
        head = S25_HEAD.get(cls.rstrip("*"), cls)
        flags = []
        # The "**" outlier marker sits on the chemical name, but for one row the
        # spreadsheet merge put it in the class cell (PFECA** / PFMOAA).
        marked = name + ("**" if cls.endswith("**") else "")
        if "*" in marked and "**" not in marked:
            flags.append("flagged as an outlier by one of two detection methods (KNN / isolation forest)")
        if "**" in marked:
            flags.append("flagged as an outlier by BOTH methods and removed from the final ML models")
        if "☨" in name or "†" in name:
            flags.append("binding destabilised HSA")
        if "✧" in name:
            flags.append("machine-learning test set (rest were training set)")
        if "➰" in name:
            flags.append("excluded from the ML dataset for want of molecular-descriptor data")
        clean = name.replace("*", "").replace("☨", "").replace("✧", "").replace("➰", "").strip()
        label = S25_RENAME.get(abbr or clean, abbr or clean)
        npfc = "0" if cls == "Fatty acids" else npfc_lookup.get(label, "")
        binds = bound.lower().startswith("y")
        if binds:
            ka_mM = float(ka.split("+/-")[0].split("±")[0].strip())
            ka_M = ka_mM * 1000.0
            kd_uM = 1.0 / ka_mM * 1000.0
            ka_txt = f"{ka_M:.4g}"
            ka_units = "M^-1 (paper reports mM^-1: " + ka.replace("±", "+/-") + ")"
            kd_txt = f"{kd_uM:.4g}"
            kd_units = "uM (1/Ka)"
        else:
            ka_txt, ka_units, kd_txt, kd_units = "", "no binding detected (Ka = 0)", "", ""
        out.append(dict(
            study="Starnes et al.",
            year="2025",
            pmid="39298512",
            doi="10.1093/toxsci/kfae124",
            chemical=label,
            n_perfluorinated_carbons=npfc,
            head_group=head,
            protein="human serum albumin (HSA)",
            protein_species="human",
            method=(
                "differential scanning fluorimetry (DSF), GloMelt dye with ROX reference; "
                "Ka from the concentration-response of dTm by nonlinear regression"
            ),
            ka=ka_txt,
            ka_units=ka_units,
            kd=kd_txt,
            kd_units=kd_units,
            n_binding_sites="",
            binding_site_identified="",
            unbound_fraction="",
            protein_conc="0.13",
            protein_conc_units="mM",
            ligand_conc_range=(
                "concentration-response series per chemical (Table S3, not extracted); "
                "positive controls octanoic acid 1-3 mM and PFHxA 0.3-3 mM"
            ),
            ligand_conc_units="mM",
            temperature_c="thermal ramp - Ka is an apparent affinity at the melting transition, not a 37 C equilibrium constant",
            buffer="HEPES-buffered saline: 140 mM NaCl, 50 mM HEPES, 0.38 mM Na2HPO4, pH 7.4",
            is_computational="no (experimental DSF; the ML models built on these data are computational)",
            tk_relevance_noted=(
                "yes - protein binding framed as a determinant of PFAS bioaccumulation; "
                "the dataset is the training set for binary-binding and Ka regression models"
            ),
            source_table=(
                "Table S4" + (f" | CAS {cas}" if cas else "") +
                (f" | R2 {r2}" if r2 and r2 != "N/A" else "") +
                (f" | max dTm {dtm} C" if dtm and dtm != "N/A" else "") +
                (" | " + "; ".join(flags) if flags else "")
            ),
            assay_volume_note=(
                "20 uL in sealed optical 96-well plates, >=2 independent plates, 3 replicates "
                "per concentration per plate; HSA fraction V, purity >=95%. Ka here is a "
                "thermal-shift apparent affinity: it ranks congeners but is not interchangeable "
                "with an equilibrium-dialysis Kd. Per-chemical top concentrations are in "
                "Table S3, which was not extracted, so ligand:protein could not be computed; "
                "the positive controls alone reach 23:1 over the 0.13 mM HSA"
            ),
            ligand_over_protein_max="",
        ))
    return out


# --------------------------------------------------------------------------
def key(r):
    return (r["study"], r["year"], r["protein"], r["chemical"],
            r["binding_site_identified"], r["buffer"], r["source_table"][:40])


def main():
    with open(CSV) as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames
        rows = list(rdr)
    before = len(rows)

    def drop(pred, label):
        nonlocal rows
        n = sum(1 for r in rows if pred(r))
        rows = [r for r in rows if not pred(r)]
        print(f"  dropped {n:3d} placeholder rows: {label}")

    print("Replacing placeholders")
    drop(lambda r: r["study"] == "Chen & Guo", "Chen & Guo 2009 (abstract-only)")
    drop(lambda r: r["study"] == "Gao et al." and r["year"] == "2025",
         "Gao 2025 (abstract-only, no constants)")
    drop(lambda r: r["study"] == "Han et al." and r["year"] == "2021",
         "Han 2021 (abstract-only)")
    drop(lambda r: r["study"] == "Ryu et al.", "Ryu 2024 (superseded)")
    drop(lambda r: r["study"] == "Starnes et al." and r["year"] == "2025",
         "Starnes 2025 (study-level summary)")

    new = []
    new += chen_guo_rows()
    print(f"  + {len(new):3d} Chen & Guo 2009")
    n0 = len(new)
    new += gao_rows()
    print(f"  + {len(new)-n0:3d} Gao 2025")
    n0 = len(new)
    new += han_rows()
    print(f"  + {len(new)-n0:3d} Han 2021")
    n0 = len(new)
    new.append(ryu_row())
    print(f"  + {len(new)-n0:3d} Ryu 2024")
    n0 = len(new)
    if TABLE_S4 and os.path.exists(TABLE_S4):
        npfc_lookup = {}
        for r in rows:
            c, n = r["chemical"], r["n_perfluorinated_carbons"]
            if n and n.isdigit() and c not in npfc_lookup:
                npfc_lookup[c] = n
        new += starnes2025_rows(TABLE_S4, npfc_lookup)
        print(f"  + {len(new)-n0:3d} Starnes 2025 Table S4")
    else:
        print("  !! Table S4 path not given or missing; Starnes 2025 rows NOT added")

    seen = {key(r) for r in rows}
    added = 0
    for r in new:
        full = {f: r.get(f, "") for f in fields}
        unknown = set(r) - set(fields)
        if unknown:
            raise SystemExit(f"unknown fields: {unknown}")
        if key(full) in seen:
            continue
        seen.add(key(full))
        rows.append(full)
        added += 1

    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\n{before} -> {len(rows)} rows ({added} added)")


if __name__ == "__main__":
    main()
