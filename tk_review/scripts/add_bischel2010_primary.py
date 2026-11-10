#!/usr/bin/env python3
"""Bischel 2010 read at source, replacing six abstract-level rows.

Bischel HN, MacManus-Spencer LA, Luthy RG.  "Noncovalent interactions
of long-chain perfluoroalkyl acids with serum albumin."  Environ Sci
Technol 2010;44(13):5263-5269.  doi 10.1021/es101334s, PMID 20540534.
Closed access; the full text and SI were supplied directly.

WHAT WAS THERE.  Six rows, all sourced "abstract", carrying order-of-
magnitude constants ("~1e6", "~1e5", "~1e4") with no stoichiometry,
no concentrations and no ratio.  Downstream that became Kd = 1 uM for
PFOA and PFNA -- "Ka ~1e6" rounded.

THE PAPER MEASURES log K_PW DIRECTLY, which almost nothing else in
this file does.  Table 1 reports the protein-water partition
coefficient for PFNA at 500 uM albumin, computed from the measured
bound fraction with a partial specific volume of 0.733 mL/g -- the
same quantity the bridge equation elsewhere has to INFER from Kd and
n.  These rows need no conversion and no assumed stoichiometry.

TABLE 1, and it carries two comparisons this project has been short of:

   albumin  conditions                          % bound   log K_PW   [L]:[P]
   HSA      50 mM Na phosphate pH 7.4          99.95     4.93+-0.05  5e-4
   BSA      50 mM Na phosphate pH 7.4          99.92     4.74+-0.05  6e-4
   BSA      same + 9 g/L NaCl                  99.89     4.56+-0.05  10e-4

  SPECIES, within one study and one method: HSA 4.93 against BSA
  4.74, a difference of 0.19 log.
  SALT, within one study: adding 9 g/L NaCl -- physiological saline
  -- drops log K_PW by 0.18, from 4.74 to 4.56.  That is a direct
  measurement of the salt effect on albumin binding at a
  PHYSIOLOGICAL ligand:protein ratio, which is what Gao 2025's salt
  series cannot offer because it runs at thousands of ligands per
  albumin.
  FATTY ACIDS: essentially fatty-acid-free BSA gave 4.80 +- 0.08
  against 4.74 for the gamma-globulin-free material, so "trace fatty
  acids had little influence on the binding affinities" HERE -- a
  useful counterweight to the defatting caveat that Moro 2022 and
  Han 2003 raise.

THE RATIO IS THE POINT OF THE PAPER.  These are the only rows in the
file measured at a genuinely PHYSIOLOGICAL loading: 5-10 x 10^-4
ligands per albumin, where the paper finds >99.9% bound.  Its thesis
is that earlier studies reporting weak binding did so because they
worked at high ratios: "Many studies have used relatively high
PFAA:albumin mole ratios ... and may explain weaker observed
binding."  Table 2 demonstrates it inside this one paper -- the same
PFOA data fitted over 0.04-5 gives Ka 0.2e6, and over 0.04-70 gives
1.4e6.

TABLE 2, equilibrium dialysis against 1 uM BSA:

   PFOA  Ka 0.20(+-0.14)e6  n 4.3(+-2.0)  R2 0.772  [L]:[P] 0.04-5    one-class
   PFOA  Ka 1.4 (+-1.9) e6  n 1.4(+-0.5)  R2 0.913  [L]:[P] 0.04-70   two-class
   PFNA  Ka 1.1 (+-0.2) e6  n 4.6(+-0.3)  R2 0.943  [L]:[P] 0.02-5    one-class
   PFNA  Ka 3.3 (+-3.2) e6  n 2.9(+-0.7)  R2 0.953  [L]:[P] 0.02-120  two-class

Both fits of each compound are kept.  They are not competing
estimates of one number: they are the same data read over different
loading windows, and the pair IS the paper's argument.  The
two-class PFOA constant has a 136% relative error and the one-class
fit an R2 of 0.772, so neither is precise, and that is recorded
rather than hidden by picking one.

nanoESI-MS is NOT ingested as a constant.  It is a gas-phase method
at 50 uM BSA in 9 mM ammonium acetate, reporting up to eight occupied
sites at a 4:1 mole ratio and Ka on the order of 1e5 -- an order of
magnitude below the solution values from the same paper.  It is
recorded as a qualitative row, because a gas-phase association
constant is not a solution binding constant.

Run:  python3 tk_review/scripts/add_bischel2010_primary.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
STUDY, YEAR = "Bischel et al.", "2010"
PMID, DOI = "20540534", "10.1021/es101334s"
BUF = "50 mM sodium phosphate, pH 7.4"

DIAL = ("equilibrium dialysis + LC-MS/MS, Spectra/Por membrane, 48 h to "
        "equilibrium; log K_PW from the measured bound fraction with a "
        "partial specific volume of 0.733 mL/g")
TITR = ("equilibrium dialysis + LC-MS/MS, binding isotherm fitted to a "
        "stepwise-equilibrium model")
MS = "nanoelectrospray ionization MS (nanoESI-MS), gas phase"

T1_NOTE = ("Table 1, measured at 500 uM albumin and a PHYSIOLOGICAL "
           "ligand:protein ratio. log K_PW is MEASURED, not converted: it "
           "comes from the bound fraction and a partial specific volume of "
           "0.733 mL/g, so it needs no Kd and no assumed stoichiometry. "
           "Initial dialysis-bag PFAA ran 1.6 to 2700 uM before "
           "equilibration; mass-balance controls recovered 94-107%, so "
           "PFAAs did not bind the membrane or vessels. Within-study "
           "comparisons: HSA 4.93 vs BSA 4.74 (species, 0.19 log) and BSA "
           "4.74 vs 4.56 with 9 g/L NaCl (salt, 0.18 log). Fatty-acid-free "
           "BSA gave 4.80 +/- 0.08, so trace fatty acids had little "
           "influence here")
T2_NOTE = ("Table 2, equilibrium dialysis against 1 uM BSA over a wide "
           "loading range; initial bag PFAA 1.6-2700 uM. BOTH fits of each "
           "compound are kept: they are the same data read over different "
           "[L]:[P] windows, and the pair is the paper's central argument "
           "-- that studies reporting weak PFAA-albumin binding did so "
           "because they worked at high ratios. PFOA over 0.04-5 gives "
           "Ka 0.2e6, the same PFOA over 0.04-70 gives 1.4e6")

ROWS = [
 # Table 1: measured log K_PW. Stored as a log D so no conversion happens.
 dict(chem="PFNA", prot="human serum albumin (HSA)", sp="human", meth=DIAL,
      ka="4.93", kau="log K_PW (L_water/kg_protein), MEASURED",
      kd="", n="", prot_c="500", lig="", ratio="0.0005",
      fb="0.05% unbound (99.95 +/- 0.01% bound)", buf=BUF,
      src="Table 1, HSA row", note=T1_NOTE),
 dict(chem="PFNA", prot="bovine serum albumin (BSA)", sp="bovine", meth=DIAL,
      ka="4.74", kau="log K_PW (L_water/kg_protein), MEASURED",
      kd="", n="", prot_c="500", lig="", ratio="0.0006",
      fb="0.08% unbound (99.92 +/- 0.01% bound)", buf=BUF,
      src="Table 1, BSA row", note=T1_NOTE),
 dict(chem="PFNA", prot="bovine serum albumin (BSA)", sp="bovine", meth=DIAL,
      ka="4.56", kau="log K_PW (L_water/kg_protein), MEASURED",
      kd="", n="", prot_c="500", lig="", ratio="0.001",
      fb="0.11% unbound (99.89 +/- 0.01% bound)",
      buf=BUF + "; plus 9 g/L NaCl (physiological saline)",
      src="Table 1, BSA + NaCl row", note=T1_NOTE +
      " | THE SALT ROW: 9 g/L NaCl lowers log K_PW from 4.74 to 4.56, a "
      "0.18 log drop measured at physiological loading"),
 # Table 2: association constants over wide loading.
 dict(chem="PFOA", prot="bovine serum albumin (BSA)", sp="bovine", meth=TITR,
      ka="2.0e5", kau="M^-1 (Ka,1)", kd="5", n="4.3",
      prot_c="1", lig="0.04 to 5 (as [L]:[P])", ratio="5",
      fb="", buf=BUF, src="Table 2, PFOA one-class fit, [L]:[P] 0.04-5",
      note=T2_NOTE + " | one-class fit, R2 0.772, n 4.3 +/- 2.0, "
                     "Ka 0.20 +/- 0.14 x1e6"),
 dict(chem="PFOA", prot="bovine serum albumin (BSA)", sp="bovine", meth=TITR,
      ka="1.4e6", kau="M^-1 (Ka,1)", kd="0.7143", n="1.4",
      prot_c="1", lig="0.04 to 70 (as [L]:[P])", ratio="70",
      fb="", buf=BUF, src="Table 2, PFOA two-class fit, [L]:[P] 0.04-70",
      note=T2_NOTE + " | two-class fit, R2 0.913, n 1.4 +/- 0.5, "
                     "Ka 1.4 +/- 1.9 x1e6 -- a 136% relative error"),
 dict(chem="PFNA", prot="bovine serum albumin (BSA)", sp="bovine", meth=TITR,
      ka="1.1e6", kau="M^-1 (Ka,1)", kd="0.9091", n="4.6",
      prot_c="1", lig="0.02 to 5 (as [L]:[P])", ratio="5",
      fb="", buf=BUF, src="Table 2, PFNA one-class fit, [L]:[P] 0.02-5",
      note=T2_NOTE + " | one-class fit, R2 0.943, n 4.6 +/- 0.3, "
                     "Ka 1.1 +/- 0.2 x1e6 -- the best-constrained fit here"),
 dict(chem="PFNA", prot="bovine serum albumin (BSA)", sp="bovine", meth=TITR,
      ka="3.3e6", kau="M^-1 (Ka,1)", kd="0.303", n="2.9",
      prot_c="1", lig="0.02 to 120 (as [L]:[P])", ratio="120",
      fb="", buf=BUF, src="Table 2, PFNA two-class fit, [L]:[P] 0.02-120",
      note=T2_NOTE + " | two-class fit, R2 0.953, n 2.9 +/- 0.7, "
                     "Ka 3.3 +/- 3.2 x1e6 -- a 97% relative error"),
 # nanoESI-MS, qualitative only.
 dict(chem="PFOA", prot="bovine serum albumin (BSA)", sp="bovine", meth=MS,
      ka="", kau="", kd="", n="up to 8 occupied sites at a 4:1 mole ratio",
      prot_c="50", lig="200 (4:1 mole ratio)", ratio="4", fb="",
      buf="9 mM ammonium acetate, pH 7",
      src="Figure 2 and Results, nanoESI-MS",
      note=("GAS PHASE, no solution constant taken. nanoESI-MS at 50 uM BSA "
            "resolves up to eight occupied sites at a 4:1 PFAA:albumin mole "
            "ratio and gives Ka on the order of 1e5 M^-1 -- an order of "
            "magnitude BELOW the solution values in Table 2 from the same "
            "paper. A gas-phase association constant is not a solution "
            "binding constant, so the value is recorded here and not in "
            "the ka column")),
 dict(chem="PFNA", prot="bovine serum albumin (BSA)", sp="bovine", meth=MS,
      ka="", kau="", kd="", n="up to 8 occupied sites at a 4:1 mole ratio",
      prot_c="50", lig="200 (4:1 mole ratio)", ratio="4", fb="",
      buf="9 mM ammonium acetate, pH 7",
      src="Figure 2 and Results, nanoESI-MS",
      note=("GAS PHASE, no solution constant taken -- see the PFOA "
            "nanoESI-MS row. Ka on the order of 1e5 M^-1, an order of "
            "magnitude below this paper's own solution values")),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)
    old = [r for r in rows if r["study"] == STUDY and r["year"] == YEAR]
    keep = [r for r in rows if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"  removing {len(old)} abstract-level Bischel 2010 row(s)")
    new = []
    for d in ROWS:
        new.append(dict(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=d["chem"],
            n_perfluorinated_carbons={"PFOA": "7", "PFNA": "8"}[d["chem"]],
            head_group="carboxylate", protein=d["prot"],
            protein_species=d["sp"], method=d["meth"],
            ka=d["ka"], ka_units=d["kau"], kd=d["kd"],
            kd_units=("uM" if d["kd"] else ""), n_binding_sites=d["n"],
            binding_site_identified="", unbound_fraction=d["fb"],
            protein_conc=d["prot_c"], protein_conc_units="uM",
            ligand_conc_range=d["lig"],
            ligand_conc_units=("uM" if d["lig"] and "[L]:[P]" not in d["lig"]
                               else "mole ratio"),
            temperature_c="not stated", buffer=d["buf"], is_computational="no",
            tk_relevance_noted=(
                "yes - the paper's thesis is that PFAA-albumin binding must "
                "be measured at PHYSIOLOGICAL ligand:protein ratios (<1e-3), "
                "where >99.9% is bound, and that studies reporting weak "
                "binding worked at high ratios"),
            source_table=f"{d['src']} | read from the full text, primary",
            assay_volume_note=d["note"],
            ligand_over_protein_max=d["ratio"]))
    print(f"  adding {len(new)} rows:")
    print(f"    {'PFAS':6}{'protein':10}{'method':12}{'value':>22}"
          f"{'n':>8}{'ratio':>9}")
    for r in new:
        m = ("dialysis K_PW" if "K_PW" in r["ka_units"] else
             "nanoESI-MS" if "nanoESI" in r["method"] else "dialysis Ka")
        v = (f"log K_PW {r['ka']}" if "K_PW" in r["ka_units"]
             else (f"Ka {r['ka']}" if r["ka"] else "qualitative"))
        print(f"    {r['chemical']:6}{r['protein_species']:10}{m:12}{v:>22}"
              f"{(r['n_binding_sites'] or '-')[:7]:>8}"
              f"{r['ligand_over_protein_max']:>9}")
    if args.dry_run:
        print("\n  (dry run, nothing written)")
        return 0
    out = keep + new
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader(); w.writerows(out)
    print(f"\n  {len(rows)} -> {len(out)} rows; wrote "
          f"{os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
