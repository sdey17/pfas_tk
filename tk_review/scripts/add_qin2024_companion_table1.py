#!/usr/bin/env python3
"""Qin et al. 2024 (the plasma paper) Table 1: BSA, human plasma, trout plasma.

Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L.  "Species Difference?
Bovine, Trout, and Human Plasma Protein Binding of Per- and Polyfluoroalkyl
Substances."  Environ Sci Technol 2024;58(23):9954-9966.
doi:10.1021/acs.est.3c10824.  PMID 38804966, PMC11171458 (open access).
Retrieved via PubMed Central.

THIS SCRIPT DOES TWO THINGS.

1. CORRECTS 11 ROWS THAT WERE FILED UNDER THE WRONG PAPER.

add_qin2024_tableS6.py ingested the log D_BSA/w column of SI Table S6 of the
BASELINE-TOXICITY paper (doi 10.1021/acs.est.3c09950) as if that paper had
measured it.  It did not.  That column is Table 1 of THIS paper, reprinted
verbatim -- all 11 slopes, intercepts and R2 agree to the last digit, which
qin2024_tableS6_provenance.py verifies cell by cell.  The baseline paper does
not cite it; the companion appeared in issue 23 of volume 58, the baseline
paper in issue 13.

Attribution was not the only thing wrong.  Every assay field on those rows
described the CYTOTOXICITY experiment instead of the BSA isotherm:

    field              was (the bioassay)            is (the BSA isotherm)
    -----------------  ----------------------------  ----------------------
    method             C18-SPME                      BioSPME 96-pin, pooled
                                                     with C18-SPME single
                                                     fibre for 3 compounds
    protein_conc       2% csFBS in assay medium       0.1-10 mg/mL BSA in PBS
    ligand_conc_range  from SI Table S2 C_tot         9-step 2-fold series
    buffer             cell culture medium            PBS pH 7.4
    assay_volume_note  220 uL, 96-well, 25000 cells   600 uL, 96-deep-well,
                                                     no cells

The BSA isotherm has no cells in it at all, so a cell count on the row was
plainly wrong.  All five fields are rewritten here.

2. ADDS THE OTHER TWO COLUMNS OF THE SAME TABLE, which nothing had ingested:
concentration-dependent Freundlich fits for HUMAN PLASMA and TROUT PLASMA,
13 chemicals each.  These are the species-resolved plasma counterpart of the
albumin column and they carry the paper's central result directly: the same
compound saturates differently in the two plasmas because the lipid:protein
ratio differs, not because the proteins differ.

WHAT THE FIT IS, from the paper's Data Evaluation section.  An empirical
Freundlich isotherm is fitted to bound against free concentration; after log
transformation that is a straight line of log C_bound against log Cw, giving
log K_Fr and the exponent n_Fr, from which log D at any Cw follows.  Whether
a chemical gets an equation or a flat constant is decided by n_Fr: below 0.9
the isotherm counts as concentration-dependent and the regression is
reported; at 0.90-1 (weakly dependent) or at or above 1 the authors FIX
n_Fr = 1 and report the average log D.  That rule, not chemistry, is why
PFUnA is a constant 4.75 between two neighbours that carry slopes.

The stored `ka` is the distribution ratio at Cw = 1 umol/L, i.e. 10^intercept,
because a Freundlich fit has no single constant.  The slope, R2 and the
equation verbatim are in `source_table`.

Idempotent: re-running changes nothing.

Run:  python3 tk_review/scripts/add_qin2024_companion_table1.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

DOI = "10.1021/acs.est.3c10824"
PMID = "38804966"
BASELINE_DOI = "10.1021/acs.est.3c09950"

# chemical: (nCF, head group,
#            (BSA slope, BSA intercept, BSA R2),
#            (human plasma slope, intercept, R2),
#            (trout plasma slope, intercept, R2))
# A slope of None is a value the paper reports as concentration-independent,
# its n_Fr having been fixed to 1.
T1 = {
    "PFBA":     (3,  "carboxylate",
                 (-0.298, 2.91, 0.60), (-0.402, 2.80, 0.68),
                 (None,   1.43, 0.74)),
    "PFHxA":    (5,  "carboxylate",
                 (-0.271, 3.43, 0.88), (-0.287, 3.29, 0.84),
                 (None,   2.49, 0.91)),
    "PFHpA":    (6,  "carboxylate",
                 (-0.305, 3.95, 0.81), (-0.389, 3.79, 0.89),
                 (None,   3.32, 0.93)),
    "PFOA":     (7,  "carboxylate",
                 (-0.314, 4.38, 0.93), (-0.410, 4.02, 0.86),
                 (None,   4.18, 0.98)),
    "PFNA":     (8,  "carboxylate",
                 (-0.147, 4.52, 0.61), (-0.183, 4.13, 0.91),
                 (None,   4.71, 0.98)),
    "PFUnA":    (10, "carboxylate",
                 (None,   4.75, 0.96), (None,   4.54, 0.95),
                 (None,   4.99, 0.98)),
    "HFPO-DA":  (5,  "ether carboxylate",
                 (-0.493, 3.44, 0.86), (-0.633, 3.41, 0.95),
                 (-0.646, 3.20, 0.97)),
    "PFHxS":    (6,  "sulfonate",
                 (-0.472, 4.28, 0.92), (-0.677, 3.95, 0.96),
                 (-0.236, 3.96, 0.68)),
    "PFOS":     (8,  "sulfonate",
                 (-0.379, 4.74, 0.92), (-0.336, 4.46, 0.85),
                 (-0.303, 4.91, 0.92)),
    "6:2 FTSA": (6,  "fluorotelomer sulfonate",
                 (-0.092, 3.86, 0.54), (-0.144, 3.66, 0.71),
                 (None,   4.11, 0.91)),
    "PFOSA":    (8,  "sulfonamide",
                 (-0.105, 4.28, 0.50), (-0.103, 3.87, 0.52),
                 (-0.164, 4.26, 0.60)),
    "hexaflumuron": (2, "fluorinated pesticide (98% neutral at pH 7.4)",
                     (None, 3.96, 0.97), (None, 4.29, 0.94),
                     (None, 4.61, 0.91)),
    "flubendiamide": (3, "fluorinated pesticide (98% neutral at pH 7.4)",
                      (None, 3.98, 0.96), (None, 3.88, 0.99),
                      (None, 4.46, 0.97)),
}

METHOD = ("BioSPME 96-Pin Device (C18, ~80 nL coating) on a Hamilton Star "
          "robot + LC-MS/MS; empirical Freundlich isotherm fitted to bound "
          "against free concentration, log-transformed to a linear fit of "
          "log C_bound against log Cw. For PFBA, PFOA, PFHxS and PFOS the "
          "fit POOLS these data with the single-fibre C18-SPME data of Qin "
          "et al. 2023 (doi 10.1016/j.envint.2023.107857), which overlapped; "
          "for PFOS several low-concentration points were excluded")

NFR_RULE = ("whether a chemical gets an equation or a flat constant is set by "
            "n_Fr: below 0.9 the regression is reported, at 0.90-1 or above "
            "n_Fr is FIXED to 1 and the average log D reported instead")

BSA_NOTE = ("yes - the Freundlich slope IS the saturation: a linear sorption "
            "model assumes slope 0 and these are negative, so the apparent "
            "distribution ratio falls as the assay dose rises. Reprinted as "
            "SI Table S6 of the baseline-toxicity paper (doi "
            "10.1021/acs.est.3c09950), where it supplies the protein term of "
            "that paper's mass-balance model")

PLASMA_NOTE = ("yes - D_plasma/water offered directly as a PBTK input, and "
               "the paper's species result rests on it: PFBA, PFHxA, PFHpA, "
               "PFOA and PFHxS differ between trout and human plasma at LOW "
               "concentration and converge at high, because human plasma is "
               "protein-rich (42.25 vs 4.46 mL/L) so specific protein "
               "binding dominates there, while trout plasma carries "
               "proportionally more lipid (15.46 vs 7.08 mL/L) and lipid "
               "binding masks the specific term. The species difference is "
               "composition, not protein identity")

MATRIX = {
    "BSA": dict(
        protein="bovine serum albumin (BSA)", species="bovine",
        units=("log D_BSA/w at Cw = 1 umol/L (L_water per L_protein); "
               "Freundlich-type fit, NOT an association constant"),
        conc="0.1-10 mg/mL BSA in PBS, chosen per compound to hold the "
             "bound fraction between 30 and 90%",
        conc_units="mg/mL", temp="37",
        buffer="BSA in PBS, pH 7.4, 37 C, 250 rpm overnight",
        col="BSA", note=BSA_NOTE),
    "HP": dict(
        protein="human plasma (protein + lipid)", species="human",
        units=("log D_plasma/w at Cw = 1 umol/L (L_water per L_prot+lip); "
               "Freundlich-type fit, NOT an association constant"),
        conc="plasma diluted with PBS to match the BSA protein level",
        conc_units="v/v plasma", temp="37",
        buffer="human plasma diluted in PBS, pH 7.4, 37 C, 250 rpm overnight",
        col="human plasma", note=PLASMA_NOTE),
    "TP": dict(
        protein="trout plasma (protein + lipid)", species="rainbow trout",
        units=("log D_plasma/w at Cw = 1 umol/L (L_water per L_prot+lip); "
               "Freundlich-type fit, NOT an association constant"),
        conc="plasma diluted with PBS to match the BSA protein level",
        conc_units="v/v plasma", temp="25",
        buffer="trout plasma diluted in PBS, pH 7.4, 25 C (room "
               "temperature), 250 rpm overnight",
        col="trout plasma", note=PLASMA_NOTE),
}

VOLNOTE = ("600 uL aliquots in 96-deep-well plates (glass-coated for "
           "PFUnA, PFOSA, hexaflumuron and flubendiamide to limit sorption "
           "to plastic, plastic for the other 9); NO CELLS in this assay. "
           "Plasma composition measured in the same paper (Table 3): human "
           "42.25 mL/L protein and 4.46 mL/L lipid, trout 15.46 and 7.08; "
           "protein density 1.36 kg/L, lipid 1.00")


def eq(matrix, slope, inter, r2):
    if slope is None:
        return (f"log D_{matrix}/w = {inter:.2f} (concentration-independent, "
                f"n_Fr fixed to 1) | R2 {r2}")
    return (f"log D_{matrix}/w = {slope:+.3f} x log Cw + {inter:.2f} | "
            f"R2 {r2}")


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    fields = list(rows[0])
    before = len(rows)

    # --- 1. correct the mis-filed BSA rows in place ---------------------
    fixed = 0
    for r in rows:
        if (r["study"], r["year"]) != ("Qin et al.", "2024"):
            continue
        if r["doi"] != BASELINE_DOI or not r["ka_units"].startswith(
                "log D_BSA/w"):
            continue
        nm = r["chemical"]
        if nm not in T1:
            continue
        slope, inter, r2 = T1[nm][2]
        m = MATRIX["BSA"]
        r.update(pmid=PMID, doi=DOI, method=METHOD,
                 protein_conc=m["conc"], protein_conc_units=m["conc_units"],
                 ligand_conc_range="9-step 2-fold dilution series",
                 ligand_conc_units="uM (per-compound range, SI)",
                 buffer=m["buffer"], tk_relevance_noted=m["note"],
                 source_table=(f"Table 1, BSA column | "
                               f"{eq('BSA', slope, inter, r2)} | Cw in "
                               f"umol/L | {NFR_RULE} | reprinted verbatim as "
                               f"SI Table S6 of doi {BASELINE_DOI}, which is "
                               f"where this row was previously filed"),
                 assay_volume_note=VOLNOTE)
        fixed += 1

    # --- 2. add the two plasma columns ---------------------------------
    key = {(r["study"], r["year"], r["chemical"], r["protein"],
            r["ka_units"][:48]) for r in rows}
    added = 0
    for nm, (ncf, hg, bsa, hp, tp) in T1.items():
        for tag, (slope, inter, r2) in (("HP", hp), ("TP", tp)):
            m = MATRIX[tag]
            k = ("Qin et al.", "2024", nm, m["protein"], m["units"][:48])
            if k in key:
                continue
            key.add(k)
            rows.append({f: "" for f in fields} | dict(
                study="Qin et al.", year="2024", pmid=PMID, doi=DOI,
                chemical=nm, n_perfluorinated_carbons=str(ncf),
                head_group=hg, protein=m["protein"],
                protein_species=m["species"], method=METHOD,
                ka=f"{10 ** inter:.3g}", ka_units=m["units"],
                protein_conc=m["conc"], protein_conc_units=m["conc_units"],
                ligand_conc_range="9-step 2-fold dilution series",
                ligand_conc_units="uM (per-compound range, SI)",
                temperature_c=m["temp"], buffer=m["buffer"],
                is_computational="no", tk_relevance_noted=m["note"],
                source_table=(f"Table 1, {m['col']} column | "
                              f"{eq(m['col'].split()[0], slope, inter, r2)} | "
                              f"Cw in umol/L | {NFR_RULE}"),
                assay_volume_note=VOLNOTE))
            added += 1

    # --- 3. say what basis the "specific" column is on -----------------
    # Qin's specific binding constant is Eq. 15, K_spec = Cmax/(2 Kd): the
    # distribution ratio AT HALF-SATURATION, which is exactly log10(2) =
    # 0.301 BELOW the Cw -> 0 plateau of the specific term. Verified against
    # Qin 2023 Table 1 for all four compounds there to within 0.003 log
    # (mbm_qin_vs_armitage.py, "Qin Eq15" check). A modeller reading the
    # column as a low-concentration plateau is 2x low, so the row says so.
    BASIS = (" BASIS: this is Qin Eq. 15, K_spec = Cmax/(2 Kd), the ratio AT "
             "HALF-SATURATION -- 0.301 log BELOW the Cw->0 plateau of the "
             "specific term, which is Cmax/Kd. Not a plateau")
    annotated = 0
    for r in rows:
        if r["doi"] != DOI or "specific binding" not in r["ka_units"]:
            continue
        if "nonspecific" in r["ka_units"] or BASIS in r["ka_units"]:
            continue
        r["ka_units"] += ";" + BASIS
        annotated += 1

    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"{before} -> {len(rows)} rows (+{added}); "
          f"{fixed} BSA rows re-attributed to {DOI}; "
          f"{annotated} specific-binding rows given their basis")
    return 0


if __name__ == "__main__":
    sys.exit(main())
