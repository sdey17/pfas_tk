#!/usr/bin/env python3
"""Khazaee 2021 again, with the Supporting Information in hand. I was wrong.

Khazaee M, Christie E, Cheng W, Michalsen M, Field J, Ng C. "Perfluoroalkyl
Acid Binding with Peroxisome Proliferator-Activated Receptors alpha, gamma,
and delta, and Fatty Acid Binding Proteins by Equilibrium Dialysis with a
Comparison of Methods." Toxics 2021;9(3):45. doi:10.3390/toxics9030045
Open access. SI retrieved from https://www.mdpi.com/2305-6304/9/3/45/s1

THIS SUPERSEDES fix_khazaee2021_protein_conc.py, WHICH WAS WRONG.

That script changed protein_conc from 1 uM to 0.03-0.075 uM. Its reasoning was
that the paper states the ligand series as "0.33-153.5 ng/mL across
ligand:protein mole ratios of 0.05-5", which inverts to tens of nanomolar, and
that 1 uM protein would put PFOA at the top ratio at 2070 ng/mL against a
stated maximum of 153.5. The arithmetic was right. The conclusion was not,
because it resolved an internal inconsistency in the paper by discarding the
WRONG side of it.

The methods section states the protein concentration outright:

    "A known volume of protein in buffer (20 to 50 uL) was added to reach a
     1 uM concentration for L-FABP, I-FABP, and PPAR-gamma, and 0.48 uM for
     PPAR-delta and PPAR-alpha."

and SI Figure S5 settles which statement to keep. Its abscissa is free ligand
concentration in uM, and the L-FABP curves run out to 3.75 uM (PFOA, S5A) and
3.1 uM (PFHxS, S5D). At 1 uM protein the top mole ratio of 5 is 5 uM total
ligand, and a free concentration of 3.75 uM sits comfortably inside that. At
the 0.075 uM protein I derived, the top ratio would be 0.375 uM TOTAL, which
the figure's 3.75 uM FREE exceeds tenfold. The figure refutes the derivation.

So protein_conc goes back to 1 uM, which is what it was before I touched it,
and the ng/mL sentence is recorded as the inconsistent statement it is.

CONSEQUENCE 1: PFHxS WAS WRONGLY EXCLUDED. The working table dropped it under
"Kd outside the concentration series tested", on the derived ligand ceiling of
0.384 uM. The real series reaches 3.1 uM free, which is 1.8x ABOVE the Kd of
1.7 uM, so the fit is properly bracketed. The gate is removed in
build_lfabp_working_table.py and the row returns.

CONSEQUENCE 2: THE CAPACITY IS IN THE FIGURE, AND IT EXPLAINS THE OUTLIER.
Figure S5's ordinate is "Specific binding (umol PFAS/umol protein)", so each
curve carries Bmax as well as Kd. The paper fits a single-site model in
GraphPad Prism and reports only Kd. Two plotted points plus a one-site model
determine both parameters exactly, and the recovered Kd checks against the
published value:

    PFOA  (S5A)  points (0.33, 0.027) and (3.75, 0.034)
                 -> Kd 0.096 uM  vs published 0.099    Bmax 0.035 mol/mol
    PFHxS (S5D)  points (0.25, 0.12 ) and (3.1 , 0.65 )
                 -> Kd 1.96  uM  vs published 1.7      Bmax 1.061 mol/mol

PFHxS saturates at one site per protein, as it should. PFOA saturates at
0.035 -- about one PFOA per thirty L-FABP molecules. Khazaee's PFOA Kd of
0.099 uM is the tightest value in the whole L-FABP compilation, 97x below the
median of every other source, and it is the number Modaresi 2025's knockout
dose design rests on. A capacity of 3.5% of a site says that constant
describes a small competent sub-population, not the bulk protein. Carried
through the bridge it would move log K from 5.85 to 4.39, which is no longer
an outlier at all.

NOT APPLIED TO THE SHEET. n stays 1, because the albumin sheet's rule is to
use the n a paper reports and 1 where its model fixes 1, and Khazaee reports
no Bmax. A value read off two points of an SI figure is not the same kind of
number as a printed one. It is recorded here and in audit_lfabp_anchors.py so
that the outlier is explained rather than silently carried.

CONDITIONS, now stated rather than derived: 1 uM protein, mole ratios 0.05,
0.1, 0.5, 1 and 5, so the top total ligand is 5 uM, at pH 7.4 and ionic
strength 18.1 mS/cm. Single-site nonlinear regression in GraphPad Prism V8.1.2.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "protein_binding.csv")
DOI = "10.3390/toxics9030045"

PROT = ("1 uM, as the methods state: \"A known volume of protein in buffer "
        "(20 to 50 uL) was added to reach a 1 uM concentration for L-FABP, "
        "I-FABP, and PPAR-gamma\". The paper also says its mole ratios of "
        "0.05-5 correspond to 0.33-153.5 ng/mL, which is not consistent with "
        "1 uM protein (PFOA at a ratio of 5 would be 2070 ng/mL); SI Figure "
        "S5 settles it in favour of 1 uM, its L-FABP curves running to 3.75 "
        "uM FREE ligand, tenfold beyond what 0.075 uM protein could reach")
LIG = ("0.05-5 uM total, from mole ratios of 0.05, 0.1, 0.5, 1 and 5 against "
       "1 uM protein; SI Figure S5 plots free ligand to 3.75 uM for PFOA and "
       "3.1 uM for PFHxS")
NOTE = ("SI Figure S5 ordinate is specific binding in umol PFAS/umol protein, "
        "so it carries Bmax. Recovered from two plotted points with the "
        "paper's one-site model: PFOA Bmax 0.035 mol/mol (recovered Kd 0.096 "
        "vs published 0.099), PFHxS Bmax 1.06 mol/mol (recovered Kd 1.96 vs "
        "published 1.7). The paper reports no Bmax. PFOA's 0.035 means its "
        "0.099 uM describes about one competent protein in thirty")


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh); fields, rows = rdr.fieldnames, list(rdr)
    changed = 0
    for r in rows:
        if (r.get("doi") or "") != DOI:
            continue
        for k, v in (("protein_conc", PROT), ("protein_conc_units", "uM"),
                     ("ligand_conc_range", LIG), ("ligand_conc_units", "uM"),
                     ("ligand_over_protein_max", "5"),
                     ("temperature_c", ""), ("assay_volume_note", NOTE)):
            if k == "temperature_c":
                continue
            if (r.get(k) or "") != v:
                r[k] = v; changed += 1
    if not changed:
        print("no change needed (Khazaee 2021 already corrected from the SI)")
        return 0
    shutil.copy2(DB, DB + ".bak")
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"updated {changed} cell(s) across the Khazaee 2021 rows")
    print("  protein_conc 0.03-0.075 uM -> 1 uM (the paper states it; SI Fig S5 confirms)")
    print("  ligand range recorded as 0.05-5 uM total, from the paper's own mole ratios")
    print("  Bmax recovered from SI Fig S5 recorded in assay_volume_note")
    print("  -> PFHxS is now inside the tested series and must be re-admitted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
