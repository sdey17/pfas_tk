#!/usr/bin/env python3
"""Fill the per-compound rows behind Qin et al. 2024 Table 1.

Qin W, Henneberger L, Gluge J, Konig M, Escher BI.  "Baseline Toxicity Model
to Identify the Specific and Nonspecific Effects of Per- and Polyfluoroalkyl
Substances in Cell-Based Bioassays."  Environ Sci Technol 2024;58(13):
5727-5738.  doi:10.1021/acs.est.3c09950.  Table 1.

WHY THIS SCRIPT EXISTS.  db/protein_binding.csv carried Ebert et al. 2020 as
a single narrative row -- "6 PFCAs, 3 PFSAs, 4 alternatives", source_table
"abstract (per-congener values not retrieved)".  The congener values were
never retrieved from Ebert's own paper, but Qin 2024 Table 1 reprints eight
of them, and that is a citable secondary source.  This fills them in.

ATTRIBUTION IS NOT UNIFORM ACROSS THAT COLUMN, and the footnotes say so:

    footnote d   PFBA 1.00                    -> Droge 2019 (ref 31)
    footnote e   8 compounds, PFHxA..PFOS     -> Ebert et al. 2020 (ref 25)
    footnote f   6:2 FTSA 3.87, PFOSA 4.94    -> PREDICTED by Qin from eq S10

So the eleven values become rows under three different studies, and the two
predicted ones are flagged is_computational=yes.  Treating the column as
"Ebert's measurements" would misattribute three of eleven.

Every D_lip/w row records that the number was read from Qin's table and not
from the measuring paper's own table -- a secondary transcription.  Ebert
reports membrane/water partitioning per kg of membrane; Qin states the
column as L_water per L_lipid.  With lipid density near 1 kg/L the two are
numerically close, but they are not the same basis, and the row says which
one it carries (Qin's).

ALSO ADDED, from the same table: the distribution ratios Qin MEASURED in
this paper -- D_cell/w (HEK293H PPARgamma cells) and D_SP/w (structural
protein), eleven compounds each.  Neither was in the database.  D_cell/w is
the quantity section 6.7 of MBM_REVIEW tests the published mass-balance
models against, so having the measured values beside the modelled ones is
the point of the exercise.  Table 1's other two columns, D_medium/w and
D_BSA/w, are the Table S6 equations evaluated at IC_10,free; those equations
are already in the database (add_qin2024_tableS6.py), so the IC_10 point
estimates are recorded in the note rather than as separate rows.

Idempotent: re-running changes nothing.

Run:  python3 tk_review/scripts/add_qin2024_table1.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

QIN_DOI = "10.1021/acs.est.3c09950"
QIN_PMID = "38394616"
EBERT = ("Ebert et al.", "2020", "32212724", "10.1021/acs.est.0c00175")
DROGE = ("Droge", "2019", "30572703", "10.1021/acs.est.8b05052")
QIN = ("Qin et al.", "2024", QIN_PMID, QIN_DOI)

# chemical: (nCF, head group, IC10_nom_M, IC10_free_M, logD_medium_w,
#            logD_cell_w, logD_BSA_w, logD_SP_w, logD_lip_w, footnote)
T1 = {
    "PFBA":     (3,  "carboxylate",             3.19e-3, 1.78e-3,
                 2.60, 2.40, 1.94, 2.11, 1.00, "d"),
    "PFHxA":    (5,  "carboxylate",             9.95e-4, 4.72e-4,
                 2.78, 3.45, 2.70, 2.45, 2.32, "e"),
    "PFHpA":    (6,  "carboxylate",             3.75e-4, 6.76e-5,
                 3.21, 3.69, 3.40, 2.65, 2.91, "e"),
    "PFOA":     (7,  "carboxylate",             1.30e-4, 1.84e-5,
                 3.67, 3.94, 3.98, 3.14, 3.52, "e"),
    "PFNA":     (8,  "carboxylate",             9.07e-5, 7.72e-6,
                 4.14, 3.74, 4.39, 3.39, 4.25, "e"),
    "PFUnA":    (10, "carboxylate",             2.16e-5, 3.83e-7,
                 4.74, 3.78, 4.75, 3.97, 4.54, "e"),
    "HFPO-DA":  (5,  "ether carboxylate",       6.54e-4, 2.20e-4,
                 3.18, 3.09, 2.28, 2.77, 2.41, "e"),
    "PFHxS":    (6,  "sulfonate",               1.46e-4, 4.22e-5,
                 3.14, 3.56, 3.52, 2.95, 4.13, "e"),
    "PFOS":     (8,  "sulfonate",               4.78e-5, 3.89e-6,
                 4.04, 4.32, 4.52, 3.94, 4.89, "e"),
    "6:2 FTSA": (6,  "fluorotelomer sulfonate", 3.07e-4, 4.25e-5,
                 3.45, 3.90, 3.71, 3.39, 3.87, "f"),
    "PFOSA":    (8,  "sulfonamide",             7.85e-6, 3.28e-7,
                 4.18, 3.81, 4.33, 3.52, 4.94, "f"),
}

FOOTNOTE = {
    "d": (DROGE, "no",
          "Table 1 footnote d: log D_lip/w from Droge 2019"),
    "e": (EBERT, "no",
          "Table 1 footnote e: log D_lip/w from Ebert et al. 2020"),
    "f": (QIN, "yes",
          "Table 1 footnote f: log D_lip/w PREDICTED by Qin from the linear "
          "relationship of experimental log D_lip/w against the number of "
          "fluorinated carbons (Figure S8, eq S10), not measured"),
}

SECONDARY = ("secondary transcription: read from Qin et al. 2024 "
             "(doi 10.1021/acs.est.3c09950) Table 1, not from the measuring "
             "paper's own table. Qin states the column as L_water per "
             "L_lipid; Ebert 2020 reports membrane/water partitioning per kg "
             "of membrane, a near-equal but not identical basis")

LIP_NOTE = ("yes - the liposome-water partition coefficient is the lipid "
            "sorption term in every published in vitro mass-balance model "
            "(Armitage 2014, Fischer 2017, Qin 2023/2024) and the input "
            "section 6.7 of MBM_REVIEW varies. Measured values exist for 9 "
            "of these 11 compounds; the other 2 are predicted from chain "
            "length, which is exactly the regression that section tests")

CELL_NOTE = ("yes - D_cell/w measured in the bioassay this paper runs. "
             "Section 6.7 of MBM_REVIEW compares the published mass-balance "
             "models against measured cell-water ratios; these are the "
             "measured side of that comparison, for 11 compounds spanning "
             "both head groups")

SP_NOTE = ("yes - structural protein is the cellular protein phase the "
           "mass-balance models sorb into. Qin measures it separately from "
           "albumin, so the two protein surrogates Armitage 2014 and Fischer "
           "2017 treat as interchangeable can be told apart here")

IC_NOTE = ("Table 1 also gives log D_medium/w and log D_BSA/w at IC_10,free; "
           "those are the SI Table S6 equations evaluated at that one "
           "concentration and the equations are already in this database")


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    fields = list(rows[0])
    before = len(rows)
    key = {(r["study"], r["year"], r["chemical"], r["protein"],
            r["ka_units"][:48]) for r in rows}
    added = 0

    def add(study, chem, ncf, hg, protein, species, method, val, units,
            comp, note, src, assay=""):
        nonlocal added
        k = (study[0], study[1], chem, protein, units[:48])
        if k in key:
            return
        key.add(k)
        rows.append({f: "" for f in fields} | dict(
            study=study[0], year=study[1], pmid=study[2], doi=study[3],
            chemical=chem, n_perfluorinated_carbons=str(ncf), head_group=hg,
            protein=protein, protein_species=species, method=method,
            ka=f"{val:.2f}", ka_units=units, is_computational=comp,
            tk_relevance_noted=note, source_table=src,
            temperature_c="37", assay_volume_note=assay))
        added += 1

    for chem, (ncf, hg, ic_nom, ic_free, d_med, d_cell, d_bsa, d_sp,
               d_lip, fn) in T1.items():
        study, comp, fn_text = FOOTNOTE[fn]

        # 1. liposome-water, attributed by footnote
        meth = ("solid-supported lipid membrane / liposome-water "
                "partitioning (measuring paper)") if comp == "no" else \
               ("predicted from the number of fluorinated carbons, "
                "Qin 2024 eq S10 (log D_lip/w = 0.535 N + 0.658, R2 0.977)")
        add(study, chem, ncf, hg,
            "phospholipid bilayer membrane (not a protein)",
            "model lipid bilayer (liposome)", meth, d_lip,
            "log D_lip/w at pH 7.4 (L_water per L_lipid)",
            comp, LIP_NOTE,
            f"Qin et al. 2024 Table 1, log D_lip/w column | {fn_text} | "
            f"{SECONDARY}")

        # 2. cell-water, measured here
        add(QIN, chem, ncf, hg,
            "HEK293H (PPARgamma-GeneBLAzer) whole cell, protein + lipid",
            "human", "measured at PFAS concentrations near IC_10,nom "
            "(SI Table S2); solid-phase microextraction + LC-MS/MS",
            d_cell, "log D_cell/w at pH 7.4 (L_water per L_prot+lip)",
            "no", CELL_NOTE,
            f"Qin et al. 2024 Table 1, log D_cell/w column | measured in "
            f"this paper near IC_10,nom = {ic_nom:.3g} mol/L | {IC_NOTE}",
            "220 uL per well, 96-well plate, 25000 HEK293H cells per well")

        # 3. structural protein-water, measured here
        add(QIN, chem, ncf, hg,
            "structural protein of HEK293H cells",
            "human", "measured at PFAS concentrations near IC_10,nom "
            "(SI Table S2); solid-phase microextraction + LC-MS/MS",
            d_sp, "log D_SP/w at pH 7.4 (L_water per L_protein)",
            "no", SP_NOTE,
            f"Qin et al. 2024 Table 1, log D_SP/w column | measured in this "
            f"paper near IC_10,nom = {ic_nom:.3g} mol/L | IC_10,free "
            f"{ic_free:.3g} mol/L; log D_medium/w {d_med:.2f} and "
            f"log D_BSA/w {d_bsa:.2f} at IC_10,free | {IC_NOTE}",
            "220 uL per well, 96-well plate, 25000 HEK293H cells per well")

    # Point the surviving Ebert narrative row at the per-congener rows it
    # used to stand in for, rather than deleting a row that carries the
    # COSMO-RS and permeability findings nothing else does.
    for r in rows:
        if r["study"] == "Ebert et al." and \
                r["chemical"].startswith("6 PFCAs"):
            r["source_table"] = (
                "abstract (per-congener values not retrieved from this paper "
                "itself) | 8 of the congener log D_lip/w values are now in "
                "this database as separate Ebert rows, transcribed from Qin "
                "et al. 2024 Table 1 footnote e; this row is kept for the "
                "COSMO-RS and P_ion findings, which are the paper's and are "
                "not in Qin's table")

    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"{before} -> {len(rows)} rows (+{added})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
