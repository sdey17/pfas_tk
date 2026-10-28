#!/usr/bin/env python3
"""Two studies of PFAS binding to LIVER proteins, not albumin.

(1) Khazaee M, Christie E, Cheng W, Michalsen M, Field J, Ng C.
    "Perfluoroalkyl Acid Binding with Peroxisome Proliferator-Activated
    Receptors alpha, gamma, and delta, and Fatty Acid Binding Proteins by
    Equilibrium Dialysis with a Comparison of Methods."
    Toxics 2021;9(3):45. doi:10.3390/toxics9030045. PMID 33652875,
    PMC7996760, open access (CC BY). Table 3, read from the publisher's
    HTML because the PMC full text strips tables.

    17 protein-PFAS pairs by EQUILIBRIUM DIALYSIS: L-FABP, I-FABP and the
    three PPAR isoforms. First experimental data for I-FABP and PPAR-delta
    with any PFAS, and the first sub-micromolar Kd reported for short-chain
    PFAS with a PPAR.

(2) Jia Y, Zhu Y, Xu D, Feng X, Yu X, Shan G, Zhu L.
    "Insights into the Competitive Mechanisms of Per- and Polyfluoroalkyl
    Substances Partition in Liver and Blood."
    Environ Sci Technol 2022;56(10):6192-6200.
    doi:10.1021/acs.est.1c08493. PMID 35436088. Table 1, from the PDF.

    11 PFAS against hL-FABP by MICROSCALE THERMOPHORESIS at 25 C, with
    HSA alongside -- the first study in this database to put an L-FABP and
    an albumin constant for the same compound in the same table.

    Six of the eleven HSA values are NOT this paper's. They carry a
    superscript 19 in Table 1, pointing at Gao et al. 2019 (Environ Sci
    Technol 53:6529, doi 10.1021/acs.est.9b00715), which measured them by
    equilibrium dialysis. Those six are filed HERE UNDER GAO, with the
    transcription route recorded, exactly as the Ebert liposome values were
    filed under Ebert rather than under the Qin paper that reprinted them.
    Jia keeps the five HSA values it measured and all eleven hL-FABP ones.

WHY THESE TWO TOGETHER. They measure the same quantity -- Kd for PFAS on
human L-FABP -- by two different methods, and a third already in the
database (fluorescence displacement, via Han 2021 citing Zhang 2013) covers
the same compounds. The three disagree by about two orders of magnitude and
in a consistent direction:

    PFOA on hL-FABP   EqD 0.099 uM   MST 17.95 uM   fluorescence 50.4 uM
    PFOS on hL-FABP   EqD 0.18  uM   MST 16.45 uM   fluorescence 18.5 uM

Khazaee's own conclusion, from a 37-study comparison, is that this is
systematic rather than particular: "EqD consistently generates lower KD
values (stronger binding affinities) than other approaches", and "binding
affinities of PFAS to proteins considered here ... have been historically
underestimated". Every row below carries its method in `method` so the
comparison cannot be made by accident.

Idempotent: re-running changes nothing.

Run:  python3 tk_review/scripts/add_khazaee2021_jia2022.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

# ---------------------------------------------------------------- Khazaee
K_DOI, K_PMID = "10.3390/toxics9030045", "33652875"
K_METHOD = ("equilibrium dialysis (Slide-A-Lyzer mini, 10 kDa MWCO), 36 h "
            "on a rocker at room temperature, ligand:protein mole ratios "
            "0.05/0.1/0.5/1/5; LC-MS/MS of the dialysate; single-site "
            "nonlinear fit. Purified human proteins (Novus Biologicals)")
K_BUFFER = "18.1 mS/cm phosphate-buffered saline, pH 7.4"
K_NOTE = ("yes - equilibrium dialysis is the comparison arm of this paper's "
          "central methodological result: across 37 studies it consistently "
          "returns LOWER Kd (stronger binding) than fluorescence "
          "displacement, by up to 3-4 orders of magnitude, for FABPs, PPARs "
          "and both serum albumins. The authors leave open whether EqD "
          "measures binding the other methods miss, in which case published "
          "affinities are systematically too weak")
# (protein label, species, protein conc uM, [(chemical, nCF, head, Kd uM, SE)])
# Kd None = "ND" in Table 3: no constant could be determined, i.e. low to no
# binding. That is a result and is carried as one, not as a blank.
KHAZAEE = [
    ("liver fatty acid-binding protein (L-FABP)", 1.0, [
        ("PFHxA", 5, "carboxylate", None, None),
        ("PFOA", 7, "carboxylate", 0.099, 0.015),
        ("PFBS", 4, "sulfonate", None, None),
        ("PFHxS", 6, "sulfonate", 1.7, 0.031),
        ("PFOS", 8, "sulfonate", 0.18, 0.032)]),
    ("intestinal fatty acid-binding protein (I-FABP)", 1.0, [
        ("PFHpA", 6, "carboxylate", None, None),
        ("PFNA", 8, "carboxylate", None, None)]),
    ("peroxisome proliferator-activated receptor alpha (PPAR-alpha)", 0.48, [
        ("PFBA", 3, "carboxylate", None, None),
        ("PFHxA", 5, "carboxylate", 0.097, 0.070),
        ("PFHpA", 6, "carboxylate", None, None),
        ("PFNA", 8, "carboxylate", 0.083, 0.028)]),
    ("peroxisome proliferator-activated receptor gamma (PPAR-gamma)", 1.0, [
        ("PFOA", 7, "carboxylate", 0.057, 0.027),
        ("PFOS", 8, "sulfonate", 8.5, 0.46)]),
    ("peroxisome proliferator-activated receptor delta (PPAR-delta)", 0.48, [
        ("PFBA", 3, "carboxylate", 0.044, 0.013),
        ("PFBS", 4, "sulfonate", None, None),
        ("PFHxS", 6, "sulfonate", 0.035, 0.0020),
        ("PFOS", 8, "sulfonate", 0.69, 0.33)]),
]

# ---------------------------------------------------------------- Jia
J_DOI, J_PMID = "10.1021/acs.est.1c08493", "35436088"
J_METHOD_FABP = ("microscale thermophoresis (MST), Monolith NT.115; Kd over "
                 "21-40 C fitted by the van't Hoff equation, the value here "
                 "being 25 C")
J_METHOD_HSA = ("microscale thermophoresis (MST), Monolith NT.115, 25 C")
J_NOTE = ("yes - the first source in this database to report an L-FABP and "
          "an albumin dissociation constant for the same compounds in the "
          "same table, which is what the Kd-to-partition bridge needs and "
          "has never had. The paper's result is that liver-blood "
          "partitioning R_L/S is set by COMPETITION between the two: R_L/S "
          "correlates with log Kd(hL-FABP)/log Kd(HSA) at r_s 0.89-0.97, "
          "better than with either constant alone")
# chemical: (nCF, head group, hL-FABP Kd uM, HSA Kd uM, HSA from ref 19?)
# Table 1 marks six HSA values with a superscript 19 = Gao et al. 2019,
# doi 10.1021/acs.est.9b00715. Those are NOT this paper's measurements.
JIA = {
    "PFHxA":        (5,  "carboxylate",              538.07, 195.0,  True),
    "PFOA":         (7,  "carboxylate",               17.95, 115.0,  True),
    "PFNA":         (8,  "carboxylate",                6.62,  73.0,  True),
    "PFDoA":        (11, "carboxylate",                3.45, 103.0,  True),
    "HFPO-DA":      (5,  "ether carboxylate",        267.78, 1156.0, False),
    "HFPO-TA":      (8,  "ether carboxylate",         19.94,   84.0, False),
    "HFPO-TeA":     (11, "ether carboxylate",          1.22,  179.0, False),
    "PFOS":         (8,  "sulfonate",                 16.45,   38.0, True),
    "PFEESA":       (4,  "ether sulfonate",         3138.20,  874.0, False),
    "6:2 Cl-PFESA": (8,  "chlorinated ether sulfonate", 12.47, 67.0, True),
    "8:2 Cl-PFESA": (10, "chlorinated ether sulfonate",  5.39, 31.0, False),
}
GAO19 = "10.1021/acs.est.9b00715"


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    fields = list(rows[0])
    before = len(rows)
    key = {(r["study"], r["year"], r["chemical"], r["protein"],
            (r["kd_units"] or r["ka_units"])[:44]) for r in rows}
    added = 0

    def add(**kw):
        nonlocal added
        k = (kw["study"], kw["year"], kw["chemical"], kw["protein"],
             (kw.get("kd_units") or kw.get("ka_units", ""))[:44])
        if k in key:
            return
        key.add(k)
        rows.append({f: "" for f in fields} | kw)
        added += 1

    # ---- Khazaee 2021 -------------------------------------------------
    for protein, pconc, panel in KHAZAEE:
        for chem, ncf, hg, kd, se in panel:
            detected = kd is not None
            add(study="Khazaee et al.", year="2021", pmid=K_PMID, doi=K_DOI,
                chemical=chem, n_perfluorinated_carbons=str(ncf),
                head_group=hg, protein=protein, protein_species="human",
                method=K_METHOD,
                kd=(f"{kd:g}" if detected else ""),
                kd_units=(f"uM (+/- {se:g} SE)" if detected
                          else "no binding detected; reported as ND, no "
                               "dissociation constant determinable"),
                n_binding_sites="single-site model assumed in the fit",
                protein_conc=f"{pconc:g}", protein_conc_units="uM",
                ligand_conc_range="0.33-153.5 ng/mL (ligand:protein mole "
                                  "ratio 0.05-5)",
                ligand_conc_units="ng/mL",
                temperature_c="", buffer=K_BUFFER + ", room temperature",
                is_computational="no", tk_relevance_noted=K_NOTE,
                source_table="Table 3 (equilibrium dialysis Kd +/- SE), read "
                             "from the publisher HTML; the PMC full text "
                             "strips tables",
                ligand_over_protein_max="5")

    # ---- Jia 2022 -----------------------------------------------------
    for chem, (ncf, hg, kd_fabp, kd_hsa, from_gao) in JIA.items():
        add(study="Jia et al.", year="2022", pmid=J_PMID, doi=J_DOI,
            chemical=chem, n_perfluorinated_carbons=str(ncf), head_group=hg,
            protein="liver fatty acid-binding protein (hL-FABP)",
            protein_species="human", method=J_METHOD_FABP,
            kd=f"{kd_fabp:g}", kd_units="uM (25 C)",
            protein_conc="see Text S7", protein_conc_units="",
            ligand_conc_range="MST titration series", ligand_conc_units="",
            temperature_c="25", buffer="Tris-HCl",
            is_computational="no", tk_relevance_noted=J_NOTE,
            source_table="Table 1, hL-FABP Kd column (25 C); the same table "
                         "gives dS, dH and dG from the van't Hoff fit over "
                         "21-40 C")
        if from_gao:
            # Gao's measurement, read from Jia's table. Filed under Gao.
            add(study="Gao et al.", year="2019", pmid="31099564", doi=GAO19,
                chemical=chem, n_perfluorinated_carbons=str(ncf),
                head_group=hg, protein="human serum albumin (HSA)",
                protein_species="human",
                method="equilibrium dialysis (Gao et al. 2019); value not "
                       "measured by the paper it was transcribed from",
                kd=f"{kd_hsa:g}", kd_units="uM",
                temperature_c="", buffer="",
                is_computational="no",
                tk_relevance_noted="yes - albumin affinity paired against "
                                   "L-FABP affinity for the same compound, "
                                   "which is what a liver-blood competition "
                                   "model needs",
                source_table=f"secondary transcription: read from Jia et al. "
                             f"2022 (doi {J_DOI}) Table 1, HSA column, where "
                             f"it carries superscript 19 attributing it to "
                             f"this paper. Gao et al. 2019 not retrieved at "
                             f"source; its own table may carry more "
                             f"compounds and a stated uncertainty")
        else:
            add(study="Jia et al.", year="2022", pmid=J_PMID, doi=J_DOI,
                chemical=chem, n_perfluorinated_carbons=str(ncf),
                head_group=hg, protein="human serum albumin (HSA)",
                protein_species="human", method=J_METHOD_HSA,
                kd=f"{kd_hsa:g}", kd_units="uM (25 C)",
                protein_conc="see Text S7", protein_conc_units="",
                ligand_conc_range="MST titration series",
                ligand_conc_units="", temperature_c="25", buffer="Tris-HCl",
                is_computational="no", tk_relevance_noted=J_NOTE,
                source_table="Table 1, HSA Kd column, measured in this paper")

    # ---- point the old aggregate Jia row at the per-compound rows -----
    fixed = 0
    for r in rows:
        if r["doi"] == J_DOI and r["chemical"].startswith("HFPO-DA (GenX),"):
            r["source_table"] = (
                "secondary (Starnes & Belcher 2026 Table 1): an aggregated "
                "Ka RANGE over the five HSA values Jia measured himself. "
                "Superseded by the per-compound Table 1 rows now in this "
                "database, which carry the primary Kd for all 11 compounds "
                "against hL-FABP and HSA. Kept because it records how the "
                "value reached this project")
            fixed += 1

    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"{before} -> {len(rows)} rows (+{added}); "
          f"{fixed} aggregate Jia row(s) annotated as superseded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
