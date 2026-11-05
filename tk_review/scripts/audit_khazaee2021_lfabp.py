#!/usr/bin/env python3
"""Khazaee 2021's L-FABP constants, and why the flag belongs elsewhere.

Khazaee M, Christie E, Cheng W, Michalsen M, Field J, Ng C.
Perfluoroalkyl Acid Binding with Peroxisome Proliferator-Activated
Receptors alpha, gamma, and delta, and Fatty Acid Binding Proteins by
Equilibrium Dialysis with a Comparison of Methods.  Toxics
2021;9(3):45.  doi 10.3390/toxics9030045, PMID 33652875, PMC7996760.
Open access.  Full text read via PubMed Central.

WHAT PROMPTED THIS.  On the L-FABP figure Khazaee's PFOA sits at
log K 5.85, more than two logs above the next human value, from a
single Kd of 0.099 uM.  It looked like the L-FABP set's Wu 2024.
It is not.  The arithmetic says the opposite, and the flag belongs on
the study it disagrees with.

THE PAPER'S OWN ANSWER IS REAL BUT DOES NOT COVER THIS CASE.  Section
3.3 argues that equilibrium dialysis systematically returns lower Kd
than fluorescence displacement, and backs it with a literature search
over 37 studies: "the low Kd values measured here are not an artifact
of this study but rather a consistent outcome of the EqD approach."
That is a methodological claim about EqD versus FD.  It cannot explain
the disagreement with Gao 2019, because GAO 2019 IS ALSO EQUILIBRIUM
DIALYSIS, on the same protein, for the same compound, 1990-fold apart.

THE REAL DIFFERENCE IS WHERE EACH TITRATION SAT.  A dissociation
constant is only identifiable if the ligand concentrations bracket it.

    Khazaee  ligand 0.05-5 uM against its reported Kd of 0.099 uM
             = 0.5x to 50x the Kd.  Brackets it.
    Gao      ligand 2-10 uM against its reported Kd of 197 uM
             = 0.010x to 0.051x the Kd.  Never gets within 20-fold.

Run the identifiability test below and the consequence is stark.  At
Khazaee's concentrations a Kd ten times tighter would change the
measured free concentration by up to 89 %, against a stated
whole-method precision of 4-18 % -- so the data can see the Kd.  At
Gao's concentrations the same ten-fold change moves the free
concentration by 15-20 %, barely outside that precision, and only
about 3 % of the ligand is bound at all.  Gao's Kd is an
extrapolation roughly twenty-fold beyond its own data.

CONCLUSION.  Khazaee's L-FABP constants stay, unflagged.  The
better-designed titration of the two is the one that disagrees with
the pack.  What this audit actually found is that GAO 2019's L-FABP
Kd values are weakly identified, and progressively worse the shorter
the chain -- because the titration is fixed at 2-10 uM while the
reported Kd climbs from 77 uM at PFDA to 879 uM at PFBA.  Gao supplies
the only published PFBA value for L-FABP in this database, and it is
the least identifiable number in the set.

WHAT IS STILL WORTH KNOWING ABOUT KHAZAEE, none of it fatal:
  - duplicates only.  "All dialysis tests were performed in duplicate",
    five ligand:protein ratios, no replication beyond n = 2.
  - negative binding coefficients.  Some dialysate concentrations came
    out ABOVE the non-binding control, which is physically impossible;
    the paper attributes this to 15- to 3000-fold dilutions and says a
    decision tree was used to handle them.  The tree is not in the
    main text.
  - the sorption check does not establish what it is cited for.
    Recovery of 75-235 % is quoted as showing "no level of detectable
    sorption".  A window that wide cannot exclude a two-fold loss.
  - the molecular-dynamics screen that chose the pairs did not predict
    the outcomes: PFBA was predicted to bind PPAR-alpha and did not,
    PFBS was predicted for PPAR-delta and did not, PFBA bound
    PPAR-delta without being predicted, and for I-FABP both predicted
    binders showed no binding at all.  Relevant to the QSAR aim and
    recorded in db/insilico_models.csv.

Run:  python3 tk_review/scripts/audit_khazaee2021_lfabp.py
"""
import math

M_LFABP = 14208.0                      # g/mol, UniProt P07148
BRIDGE = math.log10(1e9 / M_LFABP)     # 4.8475

# (label, [protein] uM, ligand concentrations uM, reported Kd uM)
RUNS = [
    ("Khazaee 2021  human L-FABP, PFOA", 1.0,
     [0.05, 0.1, 0.5, 1.0, 5.0], 0.099),
    ("Khazaee 2021  human L-FABP, PFOS", 1.0,
     [0.05, 0.1, 0.5, 1.0, 5.0], 0.18),
    ("Khazaee 2021  human L-FABP, PFHxS", 1.0,
     [0.05, 0.1, 0.5, 1.0, 5.0], 1.7),
    ("Khazaee 2021  PPAR-delta, PFHxS (tightest reported)", 0.48,
     [0.024, 0.048, 0.24, 0.48, 2.4], 0.035),
    ("Gao 2019      human L-FABP, PFOA", 6.0,
     [2.0, 4.0, 6.0, 8.0, 10.0], 197.0),
    ("Gao 2019      human L-FABP, PFBA (the only one in the literature)",
     6.0, [2.0, 4.0, 6.0, 8.0, 10.0], 879.0),
]
PRECISION = (4.0, 18.0)   # Khazaee whole-method RSD, per cent


def free(L, P, Kd):
    """Solve L = f + P f/(Kd+f) for the free concentration."""
    b = Kd + P - L
    return (-b + math.sqrt(b * b + 4.0 * Kd * L)) / 2.0


def run(label, P, Ls, Kd, factor=10.0):
    print(f"\n{label}")
    print(f"  [protein] {P} uM,  reported Kd {Kd} uM,  ligand spans "
          f"{min(Ls)/Kd:.3g}x to {max(Ls)/Kd:.3g}x the Kd")
    print(f"  {'L_tot':>9}{'free':>11}{'bound %':>10}"
          f"{'free if Kd/10':>15}{'shift':>9}")
    worst = 0.0
    for L in Ls:
        f0 = free(L, P, Kd)
        ft = free(L, P, Kd / factor)
        shift = abs(f0 - ft) / f0 * 100.0
        worst = max(worst, shift)
        print(f"  {L:9.3f}{f0:11.4g}{100*(1-f0/L):10.1f}"
              f"{ft:15.4g}{shift:8.1f}%")
    # Three bands, not two. Clearing the precision ceiling by a hair
    # is not the same as clearing it by a factor of five, and calling
    # both "identifiable" would overstate the weaker case.
    verdict = ("IDENTIFIABLE" if worst > 2 * PRECISION[1] else
               "MARGINAL -- clears the precision ceiling only just"
               if worst > PRECISION[1] else "NOT IDENTIFIABLE")
    print(f"  A Kd ten times tighter moves the measured free "
          f"concentration by at most {worst:.1f} %.")
    print(f"  Method precision is {PRECISION[0]}-{PRECISION[1]} % RSD."
          f"   ->  {verdict}")
    return worst


def main():
    print("=" * 74)
    print("IDENTIFIABILITY: can each titration see the Kd it reports?")
    print("=" * 74)
    print("""
Equilibrium dialysis measures the FREE ligand concentration. The Kd is
identifiable only if changing it changes that measurement by more than
the assay's precision. The test below asks what a ten-fold tighter Kd
would do to the number the instrument actually sees.""")
    for args in RUNS:
        run(*args)

    print()
    print("=" * 74)
    print("THE SAME PAIR, TWO LABORATORIES, THE SAME METHOD")
    print("=" * 74)
    k, g = 0.099, 197.0
    print(f"  Khazaee 2021  PFOA / human L-FABP, equilibrium dialysis, "
          f"Kd {k} uM  ->  log K {BRIDGE - math.log10(k):.2f}")
    print(f"  Gao 2019      PFOA / human L-FABP, equilibrium dialysis, "
          f"Kd {g} uM  ->  log K {BRIDGE - math.log10(g):.2f}")
    print(f"  {g/k:,.0f}-fold apart, {math.log10(g/k):.2f} log.")
    print()
    print("  Khazaee's defence -- that equilibrium dialysis returns")
    print("  lower Kd than fluorescence displacement, shown across 37")
    print("  studies -- does not reach this comparison, because both")
    print("  of these are equilibrium dialysis.")


if __name__ == "__main__":
    main()
