#!/usr/bin/env python3
"""Add Qin et al. 2024 Table S6 to db/protein_binding.csv.

Qin W, Henneberger L, Gluge J, Konig M, Escher BI.  Environ Sci Technol
2024;58(13):5727-5738.  doi:10.1021/acs.est.3c09950.  SI Table S6.

WHAT THIS ADDS.  Concentration-dependent (Freundlich-type) distribution
ratios for ELEVEN PFAS, against the four in Qin et al. 2023:

    log D = slope x log Cw + intercept          Cw in umol/L of water

in two matrices -- bovine serum albumin (L_w per L_protein) and the whole
PPARgamma assay medium, 98% OptiMEM + 2% csFBS (L_w per L_prot+lip).

WHY IT MATTERS HERE.  The slope IS the saturation.  A linear sorption model
assumes slope = 0; every value below is negative, so the apparent affinity
falls as the assay dose rises.  This is the measured version of the argument
section 7 of MBM_REVIEW makes from a model, now across eleven compounds and
both head groups rather than four compounds.

PROVENANCE CAUTION, recorded on every row.  The SI describes assays for the
medium (BioSPME), the cells and structural protein, but NO dedicated BSA
binding assay, and its Figure S8b shows D_BSA/w for only four compounds
(PFHxS, PFOS, 6:2 FTSA, PFOSA).  For the four compounds overlapping Qin 2023
the parameters are close but not identical (PFHxS -0.480 -> -0.472, PFOS
-0.381 -> -0.379, but PFBA -0.391 -> -0.298 and PFOA -0.381 -> -0.314), so
these are a refit or a re-measurement rather than a straight reprint.  The
rows carry the equation verbatim and say so.

The stored `ka` is the distribution ratio AT Cw = 1 umol/L, i.e. 10^intercept,
because a Freundlich fit has no single constant.  The slope, R2 and the full
equation are in `source_table`.

Run:  python3 tk_review/scripts/add_qin2024_tableS6.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

# name: (nCF, head group, slope_BSA, int_BSA, R2_BSA,
#                         slope_med, int_med, R2_med, Ctot_top_M)
S6 = {
    "PFBA":     (3,  "carboxylate",          -0.298, 2.91, 0.60,
                 None,   2.60, 0.94, 8.65e-3),
    "PFHxA":    (5,  "carboxylate",          -0.271, 3.43, 0.88,
                 None,   2.78, 0.90, 3.50e-3),
    "PFHpA":    (6,  "carboxylate",          -0.305, 3.95, 0.81,
                 None,   3.21, 0.91, 1.53e-3),
    "PFOA":     (7,  "carboxylate",          -0.314, 4.38, 0.93,
                 -0.172, 3.88, 0.60, 1.52e-3),
    "PFNA":     (8,  "carboxylate",          -0.147, 4.52, 0.61,
                 -0.197, 4.31, 0.67, 1.16e-3),
    "PFUnA":    (10, "carboxylate",           None,  4.75, 0.96,
                 None,   4.74, 0.97, 2.89e-4),
    "HFPO-DA":  (5,  "ether carboxylate",    -0.493, 3.44, 0.86,
                 -0.167, 3.57, 0.79, 2.01e-3),
    "PFHxS":    (6,  "sulfonate",            -0.472, 4.28, 0.92,
                 -0.389, 3.77, 0.89, 6.26e-4),
    "PFOS":     (8,  "sulfonate",            -0.379, 4.74, 0.92,
                 -0.343, 4.25, 0.86, 2.53e-4),
    "6:2 FTSA": (6,  "fluorotelomer sulfonate", -0.092, 3.86, 0.54,
                 -0.207, 3.79, 0.66, 7.17e-4),
    "PFOSA":    (8,  "sulfonamide",          -0.105, 4.28, 0.50,
                 -0.162, 4.11, 0.70, None),
}
DOI = "10.1021/acs.est.3c09950"
PMID = "38394616"
NOTE = ("yes - the Freundlich slope IS the saturation: a linear sorption model "
        "assumes slope 0, and every measured slope here is negative, so the "
        "apparent distribution ratio falls as the assay dose rises. Eleven "
        "compounds and both head groups, against the four in Qin 2023")
CAUTION = ("The SI describes no dedicated BSA binding assay and its Figure S8b "
           "shows D_BSA/w for only 4 compounds; for the 4 overlapping Qin 2023 "
           "the parameters are close but not identical, so this is a refit or "
           "re-measurement, not a reprint")


def eq(slope, intercept):
    if slope is None:
        return f"log D = {intercept:.2f} (concentration-independent)"
    return f"log D = {slope:+.3f} x log Cw + {intercept:.2f}"


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    fields = list(rows[0])
    before = len(rows)
    key = {(r["study"], r["year"], r["chemical"], r["protein"],
            r["ka_units"][:40]) for r in rows}
    added = 0
    for nm, (ncf, hg, s_b, i_b, r_b, s_m, i_m, r_m, ctot) in S6.items():
        for matrix, slope, inter, r2, unit, meth in (
            ("bovine serum albumin (BSA)", s_b, i_b, r_b,
             "log D_BSA/w at Cw = 1 umol/L (L_water per L_protein); "
             "Freundlich-type fit, NOT an association constant",
             "C18-SPME, Freundlich-type fit of log D against log Cw"),
            ("PPARgamma assay medium, 98% OptiMEM + 2% csFBS (protein + lipid)",
             s_m, i_m, r_m,
             "log D_medium/w at Cw = 1 umol/L (L_water per L_prot+lip); "
             "Freundlich-type fit, NOT an association constant",
             "BioSPME 96-pin, Freundlich-type fit of log D against log Cw"),
        ):
            k = ("Qin et al.", "2024", nm, matrix, unit[:40])
            if k in key:
                continue
            rng = ("" if ctot is None else
                   f"{ctot * 1e6 / 512:.3g}-{ctot * 1e6:.3g}")
            rows.append({f: "" for f in fields} | dict(
                study="Qin et al.", year="2024", pmid=PMID, doi=DOI,
                chemical=nm, n_perfluorinated_carbons=str(ncf), head_group=hg,
                protein=matrix,
                protein_species="bovine" if "albumin" in matrix else "mixed",
                method=meth, ka=f"{10 ** inter:.3g}", ka_units=unit,
                is_computational="no", tk_relevance_noted=NOTE,
                protein_conc="2% csFBS in the assay medium; see SI Tables S7",
                protein_conc_units="v/v serum",
                ligand_conc_range=rng, ligand_conc_units="uM (10-point, 2x dilution series)",
                temperature_c="37", buffer="cell culture medium, pH 7.4",
                source_table=(f"SI Table S6 | {eq(slope, inter)} | R2 {r2} | "
                              f"Cw in umol/L | {CAUTION}"),
                assay_volume_note=("220 uL per well, 96-well plate, 25000 "
                                   "HEK293H cells per well"),
            ))
            added += 1
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"{before} -> {len(rows)} rows (+{added})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
