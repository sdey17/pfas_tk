#!/usr/bin/env python3
"""Replace the Sheng 2020 secondary stub with its two measured rows.

THE PAPER.  Sheng N, Wang J, Guo Y, Wang J, Dai J.  "Interactions of
Perfluorooctanesulfonate and 6:2 Chlorinated Polyfluorinated Ether
Sulfonate with Human Serum Albumin: A Comparative Study."  Chem Res
Toxicol 2020;33(6):1478-1486.  doi 10.1021/acs.chemrestox.0c00075,
PMID 32423201.

WHAT THE STUB SAID, AND WHY IT WAS WRONG.  One row, from a secondary
(Starnes & Belcher 2026 Table 1) whose source_table read "primary full
text not obtained", carrying chemical "PFOS and 6:2 Cl-PFESA (F-53B)"
and ka "3.26e4-5.99e4 M^-1".  That apparent RANGE is not a range: it
is the two compounds' constants, 1/30.7 uM = 3.26e4 and 1/16.7 uM =
5.99e4.  The stub took their geometric midpoint, 22.63 uM, and with
"n = 1 assumed by the reader, the paper reports none" produced
log K 2.823 -- a number belonging to neither compound, on a
stoichiometry the paper does print.

WHAT THE PAPER ACTUALLY REPORTS.  Figure 4, Scatchard plots from
ultrafiltration centrifugation, one point per experiment:

                 printed fit            n          Kd        log K
  PFOS           y=0.264x+8.11  R2 .994  3.8+/-0.03  30.7 uM   3.270
  6:2 Cl-PFESA   y=0.301x+5.02  R2 .988  3.3+/-0.05  16.7 uM   3.473

Both boxes follow exactly from their own fits under the Scatchard form
L_F/r = L_F/n + Kd/n, so n = 1/slope and Kd = intercept/slope:
1/0.264 = 3.788 and 8.11/0.264 = 30.72; 1/0.301 = 3.322 and
5.02/0.301 = 16.68.  The abstract states the same pair.

UNITS ARE UNAMBIGUOUS HERE.  Both axes are labelled (uM) and Kd comes
out in uM.  That is worth recording because this paper follows
Beesoon & Martin's method explicitly (its ref 29), and Beesoon's own
figures label their axes (mM) while reporting Kd in M -- the 1000x
problem in scripts/audit_beesoon2015_scatchard.py.  Sheng does not
inherit it.

MASS BALANCE CLOSES.  HSA 5 uM with PFAS dosed 1-50 uM; at the top
plotted points the totals are 33.5 uM (PFOS) and 40.7 uM
(6:2 Cl-PFESA), inside the dosed range.

CONDITIONS, all from Materials and Methods, "Calculation of
Dissociation Constants ... by Ultrafiltration Centrifugation":
fatty-acid-free HSA at 5 uM, PFOS or 6:2 Cl-PFESA at 1-50 uM, 1 h at
37 C, then 700 uL through a Centrifree device, 10 kDa cutoff, 1500 g
for 30 min at 37 C.  Maximum loading 10.

pH IS NOT STATED for this assay.  The HSA was "diluted with
HPLC-grade water"; the only pH values in the paper are the PBS the
cells were harvested in and the Tris-HCl pH 8.0 of the limited
proteolysis, neither of which is the binding experiment.  Recorded as
not stated rather than assumed to be 7.4.

NOT INGESTED: Table 1's CDOCKER docking energies over the seven
fatty-acid sites, which are computed values and barred by the standing
rule in scripts/excluded_sources.py; and the 1,8-ANS fluorescence
displacement Kd on the two purified core peptides (60.7 and 48.9 uM
for PFOS), which are peptide fragments, not albumin.

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/add_sheng2020_primary.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "Fig. 4 | read from the full text, primary 2026-11-13"

METHOD = ("ultrafiltration centrifugation (Centrifree, 10 kDa cutoff), "
          "Scatchard plot L_F/r vs L_F: slope 1/n, intercept Kd/n")

COMMON = {
    "study": "Sheng et al.", "year": "2020", "pmid": "32423201",
    "doi": "10.1021/acs.chemrestox.0c00075",
    "protein": "human serum albumin (HSA), fatty acid free",
    "protein_species": "human", "method": METHOD,
    "ka_units": "M^-1, computed as 1/K_D",
    "kd_units": "uM", "protein_conc": "5", "protein_conc_units": "uM",
    "ligand_conc_range": "1 to 50", "ligand_conc_units": "uM",
    "ligand_over_protein_max": "10", "temperature_c": "37",
    "buffer": ("HPLC-grade water, pH NOT STATED for this assay. The HSA "
               "was 'diluted with HPLC-grade water and then spiked with "
               "different concentrations'; the only pH values anywhere "
               "in the paper are the PBS the JAR cells were harvested "
               "in and the Tris-HCl pH 8.0 of the limited proteolysis, "
               "neither of which is the binding experiment. The method "
               "follows Beesoon & Martin (ref 29), whose equivalent "
               "assay was HPLC water at pH 6.4"),
    "is_computational": ("no (CDOCKER docking energies over seven "
                         "fatty-acid sites are reported alongside in "
                         "Table 1, as a separate computational strand)"),
    "tk_relevance_noted": ("yes - the paper's framing is that the "
                           "stronger HSA affinity of 6:2 Cl-PFESA 'may "
                           "contribute to its higher bioaccumulation "
                           "potential than PFOS'"),
}

NOTE = (
    "Materials and Methods, 'Calculation of Dissociation Constants of "
    "PFOS and 6:2 Cl-PFESA to HSA by Ultrafiltration Centrifugation': "
    "fatty-acid-free HSA diluted in HPLC-grade water and spiked to 1-50 "
    "uM PFAS, albumin held at 5 uM, 1 h at 37 C, then 700 uL through a "
    "Centrifree device (10 kDa) at 1500 g and 37 C for 30 min. The "
    "authors chose ultrafiltration over equilibrium dialysis "
    "deliberately: 'it's hard for PFOS and 6:2 Cl-PFESA to reach the "
    "required concentrations (hundreds to more than 1000 uM) in "
    "equilibrium dialysis assay due to their low solubilities in "
    "water'. | BOTH Fig. 4 BOXES REPRODUCE from their own printed fits "
    "under L_F/r = L_F/n + Kd/n: PFOS y=0.264x+8.11 gives n 3.788 and "
    "Kd 30.72; 6:2 Cl-PFESA y=0.301x+5.02 gives n 3.322 and Kd 16.68. "
    "Axes are labelled (uM) and Kd comes out in uM, so this paper does "
    "NOT carry the mM/M ambiguity of Beesoon 2015 whose method it "
    "follows. | REPLACES a secondary stub that read chemical 'PFOS and "
    "6:2 Cl-PFESA (F-53B)' with ka '3.26e4-5.99e4 M^-1'. That was never "
    "a range: it is these two constants, 1/30.7 uM and 1/16.7 uM, "
    "collapsed into one cell, from which the stub took a geometric "
    "midpoint of 22.63 uM and assumed n = 1 'the paper reports none' -- "
    "while the paper prints n = 3.8 and 3.3 in Fig. 4.")

ROWS = [
    dict(COMMON, chemical="PFOS", n_perfluorinated_carbons="8",
         head_group="sulfonate", ka="3.26e4", kd="30.7",
         n_binding_sites="3.8 (Scatchard fit to the ultrafiltration "
                         "isotherm; Fig. 4 gives 3.8 +/- 0.03)",
         binding_site_identified=(
             "Sudlow site II / the Trp214 region is the optimal site for "
             "PFOS; limited proteolysis plus docking (Table 1) put its "
             "best CDOCKER interaction energy at fatty-acid site 6"),
         source_table="Figure 4, PFOS panel | " + MARK,
         assay_volume_note=NOTE),
    dict(COMMON, chemical="6:2 Cl-PFESA (F-53B)",
         n_perfluorinated_carbons="8",
         head_group="chlorinated ether sulfonate", ka="5.99e4", kd="16.7",
         n_binding_sites="3.3 (Scatchard fit to the ultrafiltration "
                         "isotherm; Fig. 4 gives 3.3 +/- 0.05)",
         binding_site_identified=(
             "prefers Sudlow site I over Trp214 / Sudlow site II -- the "
             "opposite preference to PFOS; the paper attributes it to "
             "the chlorine and ether oxygen. Best CDOCKER interaction "
             "energy at fatty-acid site 2"),
         source_table="Figure 4, 6:2 Cl-PFESA panel | " + MARK,
         assay_volume_note=NOTE),
]


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())

    if any(MARK in (r.get("source_table") or "") for r in rows):
        print("  already ingested; nothing to do")
        return 0

    before = len(rows)
    rows = [r for r in rows
            if not ("Sheng" in (r.get("study") or "")
                    and r.get("year") == "2020")]
    removed = before - len(rows)
    for spec in ROWS:
        rows.append({f: spec.get(f, "") for f in fields})

    with open(DB, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  removed {removed} secondary stub row, added {len(ROWS)}:")
    for spec in ROWS:
        print(f"    {spec['chemical']:22s} Kd {spec['kd']:>5s} uM  "
              f"n {spec['n_binding_sites'].split()[0]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
