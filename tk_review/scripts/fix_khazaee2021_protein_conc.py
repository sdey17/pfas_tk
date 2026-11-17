#!/usr/bin/env python3
"""Khazaee 2021 checked against the primary: values correct, protein conc wrong.

Khazaee M, Christie E, Cheng W, Michalsen M, Field J, Ng C. "Perfluoroalkyl
Acid Binding with Peroxisome Proliferator-Activated Receptors alpha, gamma, and
delta, and Fatty Acid Binding Proteins by Equilibrium Dialysis with a
Comparison of Methods." Toxics 2021;9(3):45. doi:10.3390/toxics9030045
Open access; Table 3 read from the publisher HTML.

WHY IT WAS CHECKED. Khazaee is the tightest-binding source in the L-FABP
compilation by about two orders of magnitude: 3 rows, median log K 5.59 against
2.54-4.03 for the other six sources. Its PFOS value of 0.18 uM is also the
figure Modaresi 2025's knockout dose design rests on, so whether it is sound
decides whether that design tested saturating conditions.

TRANSCRIPTION IS CORRECT. Table 3 verbatim for L-FABP: PFHxA ND, PFOA
0.099 +/- 0.015, PFBS ND, PFHxS 1.7 +/- 0.031, PFOS 0.18 +/- 0.032 uM. All five
of our rows match.

THE ERROR IS OURS. protein_conc was stored as 1 uM. The paper does not print a
protein concentration, but it states the ligand series as 0.33-153.5 ng/mL
across ligand:protein mole ratios of 0.05-5, and that inverts:

    153.5 ng/mL / 414.07 g/mol (PFOA) / 5  =  74 nM
    0.33 ng/mL / 214.04 g/mol (PFBA) / 0.05 =  31 nM

so the working protein concentration is TENS OF NANOMOLAR, around 30-75 nM
depending on the compound's mass. The stored 1 uM is wrong by more than an
order of magnitude, and it is refutable directly: at 1 uM protein, PFOA at a
ratio of 5 would be 5 uM, i.e. 2070 ng/mL, against a stated maximum of 153.5.

WHAT THE CORRECTED CONCENTRATION DOES TO THE VALUES. It makes two of the three
better than they looked and one worse.

  PFOA  Kd 0.099 uM is 1.3x the protein concentration, and sits inside the
        tested ligand series (3.7-371 nM). Properly bracketed.
  PFOS  Kd 0.18 uM is 2.9x the protein concentration, also inside the series.
        Properly bracketed.
  PFHxS Kd 1.7 uM is 22x the protein concentration and 4.4x ABOVE the top of
        the ligand series (384 nM). That value is extrapolated beyond the data
        and should not be used as a measurement.

So the suspicion that these are an artefact of fitting below the protein
concentration -- the usual reason a Kd comes out implausibly tight -- does NOT
hold for PFOA and PFOS. They are legitimate fits at genuinely trace loading.

WHAT REMAINS WRONG WITH THE SOURCE, recorded but not acted on:

  Recovery. Table 2's sorption test spiked 2000 ng/L and recovered 101-285% in
  series 1 and 75-130% in series 2; the same spike gave 5700 and 2500 ng/L in
  two replicates. The text reports this range as "75-235%" and reads it as
  showing no sorption. A 2.3-fold spread between replicates of one spike is
  not compatible with the +/-15% standard error reported on a 0.099 uM Kd.

  Membrane. A 10 kDa MWCO device against a 14.2 kDa protein is a 1.4-fold
  margin where 3-fold is conventional. Protein crossing the membrane inflates
  the apparent free concentration and would bias Kd WEAKER, so it does not
  explain the tight values, but it is uncontrolled.

  The paper's own claim is narrower than its use elsewhere. It says its results
  "generally agreed with previous observations in terms of RELATIVE
  affinities", and validates the absolute values against its own molecular
  dynamics predictions rather than against prior measurements. It also reports
  ND for PFHxA and PFBS, which Gao 2019 and Yang 2020 both measure.

NOT EXPLAINED. The 100-fold offset from the other six sources is not accounted
for by loading. Only three L-FABP sources report a ligand:protein ratio, and
across them the relationship between loading and apparent affinity runs the
wrong way for a saturation explanation (slope +1.79 log per decade, r = +0.48
on n = 3). Protein concentration does not explain it either: Gao 2019 used
6 nM protein, lower than Khazaee's, and reports the WEAKEST binding in the
compilation (median Kd 203 uM). The offset stands unexplained.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.3390/toxics9030045'

CONC = ('0.03-0.075 (DERIVED, not printed: the paper gives the ligand series as '
        '0.33-153.5 ng/mL across ligand:protein mole ratios 0.05-5, which '
        'inverts to 31 nM for PFBA at the bottom and 74 nM for PFOA at the '
        'top. Previously stored as 1 uM, which the paper refutes: at 1 uM '
        'protein, PFOA at a ratio of 5 would be 2070 ng/mL against a stated '
        'maximum of 153.5)')

def main():
    with open(DB, newline='', encoding='utf-8') as fh:
        rdr = csv.DictReader(fh); fields, rows = rdr.fieldnames, list(rdr)
    changed = 0
    for r in rows:
        if (r.get('doi') or '') != DOI:
            continue
        if (r.get('protein_conc') or '') != CONC:
            r['protein_conc'] = CONC; changed += 1
        if (r.get('protein_conc_units') or '') != 'uM':
            r['protein_conc_units'] = 'uM'; changed += 1
    if not changed:
        print("no change needed (Khazaee 2021 protein concentration already corrected)")
        return 0
    shutil.copy2(DB, DB + '.bak')
    with open(DB, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"updated {changed} cell(s) across the Khazaee 2021 rows")
    print("  protein_conc 1 uM -> 0.03-0.075 uM, recorded as derived")
    return 0

if __name__ == '__main__':
    sys.exit(main())
