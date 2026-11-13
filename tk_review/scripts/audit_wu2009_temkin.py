#!/usr/bin/env python3
"""Wu 2009's three constants, and why one of them is not a Kd.

THE PAPER.  Wu L-L, Gao H-W, Gao N-Y, Chen F-F, Chen L.  "Interaction
of perfluorooctanoic acid with human serum albumin."  BMC Struct Biol
2009;9:31.  doi 10.1186/1472-6807-9-31, PMID 19442292, PMC2688008,
open access.  PFOA on HSA, pH 7.40 Britton-Robinson buffer with
0.15 M NaCl, 37 C -- the paper's own "normal physiological condition".
Not to be confused with Wu 2024 (doi 10.3390/toxics12010046), a
different group and a named holdout here.

WHAT IT REPORTS.  Three constants from two independent experiments,
all verified below against Table 1 and the Fig. 3B fit box:

  ITC, two-site sequential   step 1  N=1   K_b 2.47e4 M^-1
                             step 2  N=12  K_b 5.5e2 M^-1
  equilibrium dialysis,
  Temkin isotherm                    N=13  K_a 3.12e4 M^-1

The ITC fit is a proper mass-action model and its two steps are sound.
The Temkin constant is a different kind of quantity, and the sheet has
been treating it as a Kd.

WHY THE TEMKIN ROW IS NOT CONVERTIBLE.  The model the paper fits is

    gamma/N = -(RT/dQ) ln(K_a c_L)

which is logarithmic in c_L, not hyperbolic.  Three consequences, each
checked numerically below:

  1. 1/K_a is a ZERO-coverage threshold, not a half-saturation.
     Setting gamma/N = 0 gives c_L = 1/K_a = 32.05 uM, the
     concentration at which the model says NOTHING is bound.  A
     Langmuir Kd is where HALF the sites are filled.  Reading 1/K_a
     as a Kd misreads the one number by its definition.

  2. There is no trace limit.  Below 1/K_a the fitted coverage is
     NEGATIVE (-29.8 bound per HSA at 1 uM).  A protein-water
     partition coefficient is a low-occupancy quantity and this model
     has no linear region in which to take that limit, so no
     conversion exists -- not a wrong conversion, no conversion.

  3. The isotherm has no asymptote, so it passes its own saturation
     number.  At 300 uM it claims 148% of N, at 1 mM 228%, and Fig. 2
     runs to c_L0/c_M0 = 30, about 1 mM.  N = 13 is from the same
     paper.

WHAT IT COSTS.  Converted anyway, as the sheet does now, the Temkin
row gives log K 3.785.  The ITC from the same paper, aggregated over
its two classes the way this project aggregates multi-class fits
(K0 = (1e9/M) x SUM_i n_i/Kd_i), gives 2.673.  That is a 1.11 log,
13-fold disagreement between two methods in one paper at one pH and
one temperature.

The paper does claim the methods agree: "By comparison of N obtained
by equilibrium dialysis (Fig. 2) and that (N = N1 + N2 = 13) obtained
by ITC (Fig. 3B), two methods achieves the same result."  That
agreement is about N, the site COUNT.  On affinity the two are never
compared, and they differ 13-fold.

A SECOND, SEPARATE ISSUE: the ITC is stored as TWO ROWS.  Steps 1 and
2 are two classes of one binding event, so the single number that
describes the protein is the aggregate, 2.673.  As two rows they enter
the PFOA pool twice, at 2.570 and 1.997, neither of which is a
partition coefficient on its own.  The converter already knows how to
add classes, but only when they share one row; this study was ingested
as two.

Nothing is repaired here.

Run:  python3 -I tk_review/scripts/audit_wu2009_temkin.py
"""
import math
BR = 4.1774
R, T = 8.314, 310.15          # J/mol/K, 37 C

print("A.  The three constants, against the paper")
print("    Table 1 'Thermodynamic property of the HSA-PFOA interaction")
print("    at pH 7.40 at 37 C', and Fig. 3B's fit box\n")
ITC = [("step 1", 1,  24.7e3), ("step 2", 12, 0.55e3)]
for lbl, n, Kb in ITC:
    kd = 1e6/Kb
    print(f"  ITC {lbl}: N={n:<3} K_b={Kb:.3g} M^-1  -> Kd {kd:8.2f} uM"
          f"  log K {BR+math.log10(n)-math.log10(kd):.3f}")
KA_T, N_T = 3.12e4, 13
kd_t = 1e6/KA_T
print(f"  Temkin   : N={N_T}  K_a={KA_T:.3g} M^-1  -> Kd {kd_t:8.2f} uM"
      f"  log K {BR+math.log10(N_T)-math.log10(kd_t):.3f}   <- the sheet")

print("\nB.  Does the Temkin fit reproduce its own K_a?")
print("    Fig. 2 curve 3, printed fit: y = 0.662x + 6.85, R = 0.993")
print("    model:  gamma/N = -(RT/dQ) ln(K_a c_L)")
slope, icept = 0.662, 6.85
dQ = -R*T/slope
lnKa = icept/slope
print(f"      slope -> dQ    = {dQ/1000:.2f} kJ/mol   (paper: -3.89)")
print(f"      icept -> K_a   = {math.exp(lnKa):.3g} M^-1  (paper: 3.12e4)")
print("    Reproduces.  The fit is sound; the question is what K_a MEANS.")

print("\nC.  1/K_a is NOT a half-saturation")
print("    Setting gamma/N = 0 gives K_a c_L = 1, i.e. c_L = 1/K_a.")
print(f"    So at c_L = {kd_t:.2f} uM the Temkin model says coverage is")
print("    ZERO, not one half.  A Langmuir Kd is the concentration at")
print("    HALF saturation.  They are different quantities.")

print("\nD.  Temkin has no trace limit, so it yields no partition coefficient")
for c in (1, 10, 32.05, 100, 1000):
    th = slope*math.log(KA_T*c*1e-6)
    itc = sum(n*(c/(1e6/Kb))/(1+c/(1e6/Kb)) for _, n, Kb in ITC)
    print(f"    c_L {c:7.2f} uM   Temkin bound/HSA {th*N_T:7.2f}"
          f"   ITC bound/HSA {itc:6.3f}")
print("    Below 1/K_a the Temkin coverage goes NEGATIVE.  A protein-water")
print("    partition coefficient is a trace-limit quantity, and this model")
print("    has no linear region to take that limit in.")

print("\nE.  The two methods on affinity, properly compared")
s = sum(n/(1e6/Kb) for _, n, Kb in ITC)
print(f"    Multi-class sum for ITC:  SUM n_i/Kd_i = "
      f"1/{1e6/ITC[0][2]:.2f} + 12/{1e6/ITC[1][2]:.1f} = {s:.6f} uM^-1")
print(f"    log K (ITC, aggregated)   = {BR+math.log10(s):.3f}")
print(f"    log K (Temkin row in sheet) = "
      f"{BR+math.log10(N_T)-math.log10(kd_t):.3f}")
gap = (BR+math.log10(N_T)-math.log10(kd_t)) - (BR+math.log10(s))
print(f"    gap = {gap:.2f} log  ({10**gap:.0f}x) between two methods in")
print("    ONE paper, same buffer, same pH, same temperature.")
print("\n    The paper says the methods agree -- 'By comparison of N")
print("    obtained by equilibrium dialysis and that (N = N1 + N2 = 13)")
print("    obtained by ITC, two methods achieves the same result.'")
print("    That agreement is about N, the site COUNT.  On affinity they")
print("    differ 13-fold, and the paper never compares them.")
