#!/usr/bin/env python3
"""Concentrations and ratios for three studies whose full texts arrived.

Beesoon & Martin 2015, Klevens & Ellenbogen 1954 and Li 2021 all sat
in the database with the ligand:protein ratio blank.  Two of them
were sourced from the ABSTRACT alone and one from an SI spreadsheet
that carried constants but not conditions.  The full texts were
supplied directly, so the conditions are now read at source.

No binding constant changes in any of the three.  What changes is the
loading each constant was measured at, which is what decides whether
it describes a site or a coated protein.

----------------------------------------------------------------------
BEESOON & MARTIN 2015  (Environ Sci Technol 49:5722, doi
10.1021/es505399w, PMID 25826685).  "Incubation and Ultrafiltration
of Pure Isomers with HSA": HSA 0.006 mM, isomers 0.002-0.012 mM, so
the ratio runs 0.33 to 2.0.  That span is DELIBERATE, and the paper
says why: "This particular concentration of HSA was chosen such that
it falls within the range of PFOS and PFOA concentrations used in
this experiment, thus a 1:1 molar ratio was achieved at a given point
within the range of selected concentrations."  Almost no other study
in this file designed its loading on purpose.

Two conditions worth carrying: the incubation is in HPLC-grade WATER
at pH 6.4, not a pH 7.4 buffer, and residual methanol runs 2-10%.
37 C for 1 h, then 700 uL through a 30 kDa Centrifree device at
1500g.  Scatchard also gives a stoichiometry the database lacked --
close to 2 sites for linear PFOS and linear PFOA, about 1 for every
branched isomer -- so n is filled too.

----------------------------------------------------------------------
KLEVENS & ELLENBOGEN 1954  (Discuss Faraday Soc 18:277, doi
10.1039/df9541800277).  The database had one secondary row, Ka 3.22e2
via Starnes & Belcher 2026.  The primary confirms it exactly: Klevens
reports a dissociation constant of 3.1e-3 M, and 1/3.1e-3 = 3.2e2.

It also supplies what the secondary could not.  n = 85 sites
(extrapolated).  BSA at 0.5, 1 and 4% w/v, which at the paper's own
assumed MW of 68,000 is 74, 147 and 588 uM.  Mole ratios of PFOA to
BSA run past the break at 28-29 that the authors match to "27-28
active anion binding sites".

THE CONDITION THAT DISQUALIFIES IT FROM PHYSIOLOGICAL COMPARISON:
acetate buffer at pH 3.54-4.3, ionic strength 0.1.  That is three pH
units below blood and below BSA's isoelectric point, and the paper is
about flocculation and protein-protein aggregation -- "a marked loss
of potential binding sites with increase in protein concentration".
The constant is real; it is not a pH 7.4 albumin affinity.

----------------------------------------------------------------------
LI 2021  (Toxics 9:63).  26 rows already ingested from the SI
spreadsheet, with protein 600 uM but the ligand "not stated for the
26 standards" and no ratio.  Section 2.2 states the rest: 7.97 mg HSA
in PBS pH 7.4 to 600 uM in one chamber of a 96-well cell split by a
10 kDa membrane, 200 uL of either an AFFF dilution (4000-16,000x) or
the 26-PFAS mix in the other, six concentration levels, four
replicates, 37 C, 30 RPM, 108 h.

7.97 mg of HSA is 120 nmol, and 120 nmol at 600 uM is 200 uL, so both
chambers are 200 uL and the cell holds 120 nmol of albumin.  Section
3 gives "the total spiked mass (80 ng) of each of these compounds"
for PFBS, PFOA and PFOS; 80 ng of PFOA is 0.19 nmol, so the ratio is
0.19/120 = 0.0016.

That puts Li 2021 among the most lightly loaded sources in the file,
beside Bischel 2011 at 0.0013-0.0034, and far below the titration
studies.  It rests on the 80 ng applying per cell, which the paper
states for three compounds and not for the other 23, so the ratio is
written ONLY on those three and the rest keep the conditions without
a ratio.

Run:  python3 tk_review/scripts/fill_conditions_from_fulltext.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

MW = {"PFOA": 414.07, "PFOS": 500.13, "PFBS": 300.10}

BEESOON = dict(
    protein_conc="6", protein_conc_units="uM",
    ligand_conc_range="2 to 12", ligand_conc_units="uM",
    ligand_over_protein_max="2",
    temperature_c="37",
    buffer=("HPLC-grade WATER at pH 6.4, NOT a pH 7.4 buffer; residual "
            "methanol 2-10% from the isomer standards"),
    note=("ultrafiltration: 900-1000 uL of 0.006 mM HSA spiked with pure "
          "isomer standards to 1000 uL, incubated 37 C for 1 h, then "
          "700 uL through a 30 kDa Centrifree device at 1500g for 30 min "
          "at 37 C; free isomer by LC-MS/MS against the uncentrifuged "
          "total. THE LOADING IS DELIBERATE: 'This particular "
          "concentration of HSA was chosen such that it falls within the "
          "range of PFOS and PFOA concentrations used in this "
          "experiment, thus a 1:1 molar ratio was achieved at a given "
          "point within the range of selected concentrations' -- the "
          "ratio spans 0.33 to 2.0 by design, which few studies in this "
          "file can say"))
BEESOON_N = {"linear PFOS (n-PFOS)": "2 (Scatchard)",
             "linear PFOA (n-PFOA)": "2 (Scatchard)",
             "branched PFOS isomers (range)": "~1 (Scatchard)",
             "branched PFOA isomers (range)": "~1 (Scatchard)"}

KLEVENS = dict(
    protein_conc="74 to 588", protein_conc_units="uM",
    ligand_conc_range="mole ratios past the break at 28-29",
    ligand_conc_units="mole ratio PFOA:BSA",
    ligand_over_protein_max="29",
    temperature_c="20 (sedimentation reference)",
    buffer=("acetate, ionic strength 0.1, pH 3.54-4.3 -- THREE pH UNITS "
            "BELOW BLOOD and below BSA's isoelectric point"),
    note=("equilibrium dialysis in Cellophane casings against 10-20 "
          "volumes of buffered perfluoroacid, rocked 24 h; bound PFOA by "
          "surface tension, conductivity and a spectral-dye method. BSA "
          "at 0.5, 1 and 4% w/v = 74, 147 and 588 uM at the paper's own "
          "assumed MW of 68,000 with a 5.9% moisture correction. n = 85 "
          "sites extrapolated; the binding curve breaks at a mole ratio "
          "of 28-29, which the authors match to '27-28 active anion "
          "binding sites'. NOT A PHYSIOLOGICAL AFFINITY: the study is "
          "about flocculation and protein-protein aggregation ('a marked "
          "loss of potential binding sites with increase in protein "
          "concentration'), at pH 3.5-4.3. The primary confirms the "
          "secondary value exactly: Kd 3.1e-3 M, and 1/3.1e-3 = 3.2e2"))
KLEVENS_N = "85 (extrapolated); curve breaks at a mole ratio of 28-29"

LI_NOTE = ("equilibrium dialysis, 96-well, two chambers split by a 10 kDa "
           "regenerated-cellulose membrane: 7.97 mg HSA in PBS pH 7.4 to "
           "600 uM in one, 200 uL of an AFFF dilution (4000-16,000x) or "
           "the 26-PFAS mix in the other; six concentration levels, four "
           "replicates, 37 C, 30 RPM, 108 h. 7.97 mg HSA is 120 nmol and "
           "120 nmol at 600 uM is 200 uL, so both chambers are 200 uL and "
           "the cell holds 120 nmol of albumin")
LI_RATIO_NOTE = ("ratio from the paper's stated 'total spiked mass (80 ng) "
                 "of each of these compounds', given for PFBS, PFOA and "
                 "PFOS only: {ng} ng / {mw} g/mol = {nmol:.3g} nmol against "
                 "120 nmol HSA = {ratio:.4f}. Among the most lightly loaded "
                 "sources in this file, beside Bischel 2011. The other 23 "
                 "standards get the conditions but NO ratio, because the "
                 "80 ng is not stated for them")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    def touch(r, d, nmap=None, nval=None):
        for k, v in d.items():
            if k == "note":
                if v not in r["assay_volume_note"]:
                    r["assay_volume_note"] = (
                        (r["assay_volume_note"] + " | " if
                         r["assay_volume_note"].strip() else "") + v)
            else:
                r[k] = v
        if nmap is not None and r["chemical"] in nmap:
            r["n_binding_sites"] = nmap[r["chemical"]]
        if nval is not None:
            r["n_binding_sites"] = nval

    n = 0
    print("  BEESOON & MARTIN 2015")
    for r in rows:
        if r["study"] == "Beesoon & Martin" and r["year"] == "2015":
            if "technical mixtures" in r["chemical"]:
                print(f"    {r['chemical'][:46]:48} SKIPPED, whole-serum arm")
                continue
            touch(r, BEESOON, nmap=BEESOON_N); n += 1
            print(f"    {r['chemical'][:46]:48} HSA 6 uM, lig 2-12 uM, "
                  f"ratio 0.33-2.0, n={r['n_binding_sites']}")

    print("\n  KLEVENS & ELLENBOGEN 1954")
    for r in rows:
        if r["study"] == "Klevens & Ellenbogen" and r["year"] == "1954":
            touch(r, KLEVENS, nval=KLEVENS_N); n += 1
            print(f"    {r['chemical']:48} BSA 74-588 uM, ratio to 29, n=85")

    print("\n  LI 2021")
    for r in rows:
        if r["study"].startswith("Li") and r["year"] == "2021":
            if LI_NOTE not in r["assay_volume_note"]:
                r["assay_volume_note"] = (r["assay_volume_note"] + " | "
                                          + LI_NOTE)
            r["temperature_c"] = r["temperature_c"] or "37"
            if r["chemical"].startswith("PFOS ("):
                tag = ("no ratio: the paper's 80 ng is stated for PFOS as "
                       "one compound, but this database splits it into "
                       "linear and branched rows and the paper does not "
                       "say how that mass divides, so writing it on both "
                       "would double-count")
                if tag not in r["assay_volume_note"]:
                    r["assay_volume_note"] += " | " + tag
                print(f"    {r['chemical']:18} 80 ng stated but SPLIT "
                      f"linear/branched -> no ratio")
                n += 1
                continue
            mw = MW.get(r["chemical"])
            if mw:
                nmol = 80e-9 / mw * 1e9
                ratio = nmol / 120.0
                r["ligand_conc_range"] = "80 ng spiked per cell (stated)"
                r["ligand_over_protein_max"] = f"{ratio:.4f}"
                msg = LI_RATIO_NOTE.format(ng=80, mw=mw, nmol=nmol,
                                           ratio=ratio)
                if msg not in r["assay_volume_note"]:
                    r["assay_volume_note"] += " | " + msg
                print(f"    {r['chemical']:10} 80 ng = {nmol:.3g} nmol / "
                      f"120 nmol HSA = ratio {ratio:.4f}")
            n += 1
    print(f"    {sum(1 for r in rows if r['study'].startswith('Li') and r['year']=='2021' and not r['ligand_over_protein_max'])}"
          f" of the 26 keep conditions but no ratio (80 ng not stated)")

    print(f"\n  {n} rows updated")
    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader(); w.writerows(rows)
    print(f"  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
