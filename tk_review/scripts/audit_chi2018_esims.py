#!/usr/bin/env python3
"""Chi 2018: 24 constants behind a two-number "range", and they are good.

THE PAPER.  Chi Q, Li Z, Huang J, Ma J, Wang X.  "Interactions of
perfluorooctanoic acid and perfluorooctanesulfonic acid with serum
albumins by native mass spectrometry, fluorescence and molecular
docking."  Chemosphere 2018;198:442-449.
doi 10.1016/j.chemosphere.2018.01.152, PMID 29425944.

THE STUB said chemical "PFOA and PFOS", protein_species "human;
bovine", ka "0.52e5-9.90e5 M^-1" from a secondary, geometric midpoint
Kd 4.407 uM, n = 1 assumed, log K 3.533 -- one row for everything.

Table 1 holds TWENTY-FOUR constants: 2 compounds x 2 albumins x up to
3 stepwise steps x 3 mole ratios.  The stub's two numbers are exactly
their minimum (HSA-PFOA K_a,3 at 1:2, 5.2e4) and maximum (BSA-PFOS
K_a,1 at 1:0.5, 9.9e5).  This is the most extreme collapse of the
four stubs read so far.

TABLE 1, association constants from deconvolved native ESI-MS spectra
(L mol^-1):

  P:L     step    BSA-PFOA   BSA-PFOS   HSA-PFOA   HSA-PFOS
  1:0.5   K_a,1    8.3e5      9.9e5      1.3e5      2.0e5
  1:1     K_a,1    2.6e5      9.1e5      1.2e5      1.4e5
          K_a,2    1.5e5      2.7e5      6.3e4      6.9e4
  1:2     K_a,1    2.1e5      7.6e5      1.1e5      1.1e5
          K_a,2    1.1e5      2.4e5      6.1e4      7.6e4
          K_a,3    5.5e4      1.2e5      5.2e4      5.4e4

WHY THIS SOURCE IS DIFFERENT FROM THE LAST THREE STUBS.  Peng 2025,
Qin 2010 and Chen 2020 all turned on a quenching constant that is not
a binding constant.  Chi's are mass-action stepwise association
constants read from counted complexes: native ESI-MS resolves PL_1,
PL_2 ... PL_8 as separate peaks, and eq 2-4 convert their relative
abundances into K_a,j.  The stoichiometry is observed, not fitted to
an optical model.

THE ONE REAL PROBLEM: K_a,1 DRIFTS WITH THE LOADING IT WAS MEASURED
AT.

    BSA-PFOA   8.3e5 -> 2.6e5 -> 2.1e5   4.0x across 1:0.5 to 1:2
    HSA-PFOS   2.0e5 -> 1.4e5 -> 1.1e5   1.8x
    BSA-PFOS   9.9e5 -> 9.1e5 -> 7.6e5   1.3x
    HSA-PFOA   1.3e5 -> 1.2e5 -> 1.1e5   1.2x

A mass-action constant should not depend on the dose used to measure
it.  The paper reports the stepwise DECREASE within a ratio (K_a,1 >
K_a,2 > K_a,3) and reads it as negative cooperativity, but does not
remark on K_a,1 itself moving between ratios.  The 1:0.5 column is the
lowest loading and so the closest to the trace limit a partition
coefficient wants, which is the argument for taking it.

WHERE THE 1:0.5 VALUES LAND, n = 1, against each species' own pool:

    BSA-PFOA   Kd 1.205 uM   log K 4.096   +0.04 vs median 4.059
    BSA-PFOS   Kd 1.010 uM   log K 4.173   -0.26 vs median 4.433
    HSA-PFOA   Kd 7.692 uM   log K 3.291   +0.12 vs median 3.167
    HSA-PFOS   Kd 5.000 uM   log K 3.478   +0.01 vs median 3.468

Three of four inside 0.12 log of consensus.  No other source read in
this review sits that close.

AND IT CARRIES A WITHIN-STUDY SPECIES COMPARISON, same method, same
experiment, same day:

    PFOA   BSA 8.3e5 vs HSA 1.3e5   BSA 6.4x tighter, +0.81 log
    PFOS   BSA 9.9e5 vs HSA 2.0e5   BSA 5.0x tighter, +0.69 log

The sheet has only three such sources now (Bischel 2010, Fischer 2024,
Han 2003), and this is the only one that measures BOTH compounds on
BOTH albumins in one run.

CONDITIONS: albumin 5 uM (BSA 66,430 Da, HSA 66,550 Da as the paper
states), PFOA/PFOS 0, 2.5, 5, 10, 20, 30, 40 uM so loadings 0 to 8,
45 min at 37 C.  BUFFER IS 10 mM AMMONIUM ACETATE at pH 7.4, not
phosphate: ESI-MS needs a volatile buffer, so the ionic strength is
about a fifteenth of blood's and the electrostatics are not
physiological.  Fatty-acid-free and IgG-free albumin.

NOT FOR INGESTION: the molecular docking (computational), and the
fluorescence strand, whose only numeric output is a Stern-Volmer slope
for BSA-PFOA, F0/F = 1 + 7.781e4 [Q].  The paper itself treats that as
corroboration rather than a constant -- "the slope of the curve of
BSA-PFOA can be considered as association constant which is slightly
lower than K_a,1" -- and PFOS's plot is explicitly nonlinear, so no
constant is drawn from it.

Nothing is changed here.

Run:  python3 -I tk_review/scripts/audit_chi2018_esims.py
"""
import math

BR = 4.1774
RATIOS = ["1:0.5", "1:1", "1:2"]
SYS = ["BSA-PFOA", "BSA-PFOS", "HSA-PFOA", "HSA-PFOS"]
T1 = {
    "1:0.5": {"Ka1": dict(zip(SYS, [8.3e5, 9.9e5, 1.3e5, 2.0e5]))},
    "1:1": {"Ka1": dict(zip(SYS, [2.6e5, 9.1e5, 1.2e5, 1.4e5])),
            "Ka2": dict(zip(SYS, [1.5e5, 2.7e5, 6.3e4, 6.9e4]))},
    "1:2": {"Ka1": dict(zip(SYS, [2.1e5, 7.6e5, 1.1e5, 1.1e5])),
            "Ka2": dict(zip(SYS, [1.1e5, 2.4e5, 6.1e4, 7.6e4])),
            "Ka3": dict(zip(SYS, [5.5e4, 1.2e5, 5.2e4, 5.4e4]))}}
MEDIAN = {"BSA-PFOA": 4.059, "BSA-PFOS": 4.433,
          "HSA-PFOA": 3.167, "HSA-PFOS": 3.468}


def main():
    vals = [v for r in T1.values() for s in r.values()
            for v in s.values()]
    print(f"A.  Table 1 holds {len(vals)} constants, "
          f"{min(vals):.3g} to {max(vals):.3g} M^-1")
    print("    The stub's range, 0.52e5-9.90e5, is exactly that minimum")
    print("    and maximum: 2 compounds x 2 albumins x up to 3 stepwise")
    print("    steps x 3 mole ratios, in one cell.\n")

    print("B.  K_a,1 drifts with the loading it was measured at")
    print(f"  {'system':10s} " + " ".join(f"{r:>9s}" for r in RATIOS)
          + "   fold")
    for s in SYS:
        row = [T1[r]["Ka1"][s] for r in RATIOS]
        print(f"  {s:10s} " + " ".join(f"{v:9.3g}" for v in row)
              + f"   {row[0]/row[-1]:4.1f}x")
    print("  A mass-action constant should not depend on the dose used")
    print("  to measure it.  The paper remarks on the stepwise decrease")
    print("  within a ratio, not on K_a,1 moving between ratios.\n")

    print("C.  The 1:0.5 column, lowest loading, n = 1")
    for s in SYS:
        ka = T1["1:0.5"]["Ka1"][s]
        kd = 1e6 / ka
        lk = BR - math.log10(kd)
        print(f"    {s:10s} Kd {kd:6.3f} uM -> log K {lk:.3f}"
              f"   ({lk - MEDIAN[s]:+.2f} vs its species pool)")
    print("    Three of four inside 0.12 log.  No other source read in")
    print("    this review sits that close to consensus.\n")

    print("D.  Within-study species comparison, one experiment")
    for c in ("PFOA", "PFOS"):
        b, h = (T1["1:0.5"]["Ka1"]["BSA-" + c],
                T1["1:0.5"]["Ka1"]["HSA-" + c])
        print(f"    {c}: BSA {b:.3g} vs HSA {h:.3g} -> BSA {b/h:.1f}x "
              f"tighter, {math.log10(b/h):+.2f} log")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
