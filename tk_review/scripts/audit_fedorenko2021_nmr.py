#!/usr/bin/env python3
"""Fedorenko 2021: two Kd per compound, and n is reported after all.

THE PAPER.  Fedorenko M, Alesio J, Fedorenko A, Slitt A, Bothun GD.
"Dominant entropic binding of perfluoroalkyl substances (PFASs) to
albumin protein revealed by 19F NMR."  Chemosphere 2021;263:128083.
doi 10.1016/j.chemosphere.2020.128083, PMID 33297081.  PFOA, PFNA,
PFHxS and PFOS against BSA.

IT REACHES NO SHEET ROW TODAY, and that is a standing policy rather
than a judgement about this paper.  albumin_shareable.py admits the
reader's n = 1 route only for method families "fluorescence" and
"separative"; 19F NMR is neither, so Han 2003's NMR row is excluded
the same way.  Worth deciding deliberately, because this is the same
method Han uses -- Fedorenko's eq 3 and Han's eq 1 are the same
chemical-shift-perturbation form,

    [L]_T = (n [P]_T / d_delta) x d_delta_Bapp - K_d

with K_d read off as the negative y-intercept.

THE STUB said chemical "PFOA, PFNA, PFHxS, PFOS", ka "0.17e5-3.57e5
M^-1" from a secondary, and "no stoichiometry reported".  Two of those
three claims are wrong.

WHERE THE RANGE COMES FROM.  Page 6 states, in words: "we report K_a
values of 17 mM-1 (a) and 31 mM-1 (w) for PFHxS; 18 mM-1 (a) and
155 mM-1 (w) for PFOA; and 18 mM-1 (a) and 357 mM-1 (w) for PFOS."
Those six numbers span 17 to 357 mM^-1, i.e. 1.7e4 to 3.57e5 M^-1 --
exactly the stub's range.  PFNA has no numeric K_a in the text; it
appears only in Fig. 3C.

THE REAL STRUCTURE: TWO Kd PER COMPOUND, FROM TWO FLUORINE PROBES.
Alpha is the -CF2- adjacent to the anionic headgroup; omega is the
terminal -CF3.  These are not two binding events but two reporters on
the same molecule, and they disagree:

  PFAS    Kd alpha   Kd omega   ratio     log K (n=1)  alpha / omega
  PFHxS   58.82 uM   32.26 uM    1.8x       2.408   2.669
  PFOA    55.56 uM    6.45 uM    8.6x       2.433   3.368
  PFOS    55.56 uM    2.80 uM   19.8x       2.433   3.730

The alpha values are nearly identical across all three compounds
(17-18 mM^-1) while omega spans 21-fold.  The paper's reading is that
the hydrophobic tail drives binding and the headgroup barely
participates: "The lowest K_d values (highest protein affinity) were
also observed for the -CF3 omega position, which is more hydrophobic
itself".  So omega is the number that tracks affinity and alpha is
close to a constant -- which is itself the paper's finding, not noise.

n IS REPORTED.  "For PFOA, PFNA, and PFHxS this ratio was
[PFAS]:[BSA] ~ 12-14, denoting the shift from binding site to surface
adsorption", and the highlights say "PFAS bound per protein ranged
from 12 to 14".  It is a saturation stoichiometry of the same kind as
Klevens' break at 28-29, not a count from the Kd fit -- eq 3 tangles n
with d_delta_Bapp in the slope, so the fit cannot separate them.  But
"no stoichiometry reported" is false.

WHAT THE CHOICES COST, against each compound's bovine pool median:

                 alpha n=1  omega n=1  alpha n=13  omega n=13
  PFHxS  (4.165)   -1.76      -1.50      -0.64       -0.38
  PFOA   (4.096)   -1.66      -0.73      -0.55       +0.39
  PFOS   (4.409)   -1.98      -0.68      -0.86       +0.43

omega with n = 13 is the only combination that straddles consensus.

CONDITIONS, and they are good: BSA held at 10 uM, PFAS 10 uM to 1 mM,
PBS at pH 7.4 made to a stated physiological recipe (2.68 mM KCl,
1.47 mM KH2PO4, 136.89 mM NaCl, 8.06 mM Na2HPO4.7H2O), three
temperatures 298, 304 and 310 K, triplicate.  The fit is restricted to
the linear range below [PFAS]:[BSA] = 25:1.  Unusually for this
compilation, Kd DECREASE with rising temperature -- binding is
entropically driven, which is the paper's headline.

Nothing is changed here.

Run:  python3 -I tk_review/scripts/audit_fedorenko2021_nmr.py
"""
import math

BR = 4.1774
KA = {"PFHxS": (17, 31), "PFOA": (18, 155), "PFOS": (18, 357)}
MED = {"PFHxS": 4.165, "PFOA": 4.096, "PFOS": 4.409}
NSAT = 13.0


def main():
    vals = [v for p in KA.values() for v in p]
    print("A.  The stub's range is the six K_a stated in words on p6")
    print(f"    {min(vals)} to {max(vals)} mM^-1 = {min(vals)*1e3:.3g} "
          f"to {max(vals)*1e3:.3g} M^-1; the stub said 0.17e5-3.57e5.\n")

    print("B.  Two Kd per compound: alpha = CF2 by the headgroup,")
    print("    omega = terminal CF3.  Two probes, one molecule.\n")
    print(f"  {'PFAS':6s} {'Kd a':>8s} {'Kd w':>8s} {'a/w':>6s} | "
          f"{'logK a':>7s} {'logK w':>7s}")
    print("  " + "-" * 48)
    for c, (a, w) in KA.items():
        ka, kw = 1000 / a, 1000 / w
        print(f"  {c:6s} {ka:8.2f} {kw:8.2f} {ka/kw:5.1f}x | "
              f"{BR-math.log10(ka):7.3f} {BR-math.log10(kw):7.3f}")
    print("\n    alpha is nearly constant (17-18 mM^-1) across all three;")
    print("    omega spans 21-fold.  That IS the paper's finding -- the")
    print("    tail drives binding, the headgroup barely participates.")

    print("\nC.  n is reported: [PFAS]:[BSA] ~ 12-14 at the shift from")
    print("    site binding to surface adsorption.  'No stoichiometry")
    print("    reported' in the stub is false.\n")
    print(f"  {'PFAS':6s} {'median':>7s} | {'a n=1':>7s} {'w n=1':>7s} "
          f"{'a n=13':>7s} {'w n=13':>7s}")
    print("  " + "-" * 50)
    for c, (a, w) in KA.items():
        ka, kw = 1000 / a, 1000 / w
        m = MED[c]
        row = [BR - math.log10(ka), BR - math.log10(kw),
               BR + math.log10(NSAT) - math.log10(ka),
               BR + math.log10(NSAT) - math.log10(kw)]
        print(f"  {c:6s} {m:7.3f} | "
              + " ".join(f"{v-m:+7.2f}" for v in row))
    print("\n    omega with n = 13 is the only combination that")
    print("    straddles consensus rather than sitting under it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
