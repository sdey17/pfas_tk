#!/usr/bin/env python3
"""What does suppressing bile acid synthesis do to cholesterol?

A deliberately small hepatic cholesterol balance, built to answer one
question: if PFAS lowers CYP7A1 (bile acid synthesis), does serum
cholesterol go up?

THE NAIVE ARGUMENT
    Bile acid synthesis is the main route by which the body gets RID of
    cholesterol -- the carbon skeleton is not broken down, it is converted to
    bile acids and excreted.  Shut that route and cholesterol must accumulate.
    Hepatic cholesterol rises, SREBP-2 is switched off, LDL receptors fall,
    LDL stops being cleared from blood, serum LDL-C rises.

    Human genetics supports this: people homozygous for CYP7A1 loss of
    function have high LDL-C, gallstones, and poor statin response
    (Pullinger 2002).  So does pharmacology in reverse: bile acid
    sequestrants interrupt bile acid recycling, DRIVE UP CYP7A1, and lower
    LDL-C.  That is an approved drug class.

WHY IT IS NOT THAT SIMPLE HERE
    In the HepaRG experiments behind the Fragki QIVIVE, PFAS did not only
    lower CYP7A1 and CYP8B1.  They lowered HMGCR and LSS as well -- the
    rate-limiting step and a downstream step of cholesterol SYNTHESIS.  Both
    the inflow and the outflow of the hepatic cholesterol pool go down
    together.  And consistently, intracellular cholesterol "hardly changed"
    (Behr 2020; Louisse 2020).

    So the sign of the whole-body effect is not set by the fact that bile
    acid synthesis falls.  It is set by WHICH ARM FALLS FURTHER, in absolute
    flux.  This script works out where that crossover sits.

THE MODEL (dimensionless; baseline total flux = 1.0)
    liver free cholesterol C, plasma LDL-C pool L

    in : de novo synthesis (SREBP-2 controlled, HMGCR-scaled)
         LDL receptor uptake (SREBP-2 controlled, LDLR-scaled)
         dietary / remnant delivery (constant)
    out: bile acid synthesis (CYP7A1-scaled)
         biliary cholesterol secretion (rises with C)
         VLDL export (rises with C)

    SREBP-2 activity s(C) = 1 / (1 + (C/C0)^n) -- sterol-sensing, so a fall
    in hepatic cholesterol DERESPRESSES synthesis and LDLR, and a rise
    represses them.  That feedback is the whole mechanism linking hepatic
    cholesterol to serum LDL-C.

Run:  python3 tk_learning/cholesterol/bile_acid_cholesterol.py
"""

# ---- baseline flux shares of the hepatic cholesterol pool (sum to 1.0) ----
# Chosen so that EVERY compartment balances at C = L = 1:
#   liver  in  0.45 synthesis + 0.35 LDLR uptake + 0.20 diet      = 1.00
#   liver  out 0.35 bile acids + 0.40 biliary + 0.25 VLDL         = 1.00
#   plasma in  0.25 VLDL + 0.10 peripheral return                 = 0.35
#   plasma out 0.35 LDLR uptake                                   = 0.35
#   body   in  0.45 hepatic + 0.10 peripheral synthesis + 0.20 diet = 0.75
#   body   out 0.35 bile acids + 0.40 net biliary loss             = 0.75
# Roughly human in proportion: whole-body synthesis ~1 g/d, bile acid
# synthesis ~0.5 g/d, dietary absorption ~0.3 g/d.  Treat as a toy.
F_SYN, F_UP, F_DIET = 0.45, 0.35, 0.20        # liver in
F_BA, F_BILE, F_VLDL = 0.35, 0.40, 0.25       # liver out
F_PERIPH = 0.10                                # reverse transport into plasma
HILL = 4.0                                     # sterol-sensing sharpness


def srebp(C):
    """SREBP-2 activity, normalised to 1.0 at baseline (C = 1)."""
    return (1.0 / (1.0 + C ** HILL)) / (1.0 / (1.0 + 1.0 ** HILL))


def steady_state(a_cyp7a1=1.0, a_hmgcr=1.0, a_ldlr=1.0,
                 iters=400000, dt=0.002):
    """Integrate to steady state. a_* are activity multipliers, 1.0 = normal."""
    C, L = 1.0, 1.0
    for _ in range(iters):
        s = srebp(C)
        v_syn = F_SYN * a_hmgcr * s
        v_up = F_UP * a_ldlr * s * L
        v_ba = F_BA * a_cyp7a1
        v_bile = F_BILE * C
        v_vldl = F_VLDL * C
        dC = v_syn + v_up + F_DIET - v_ba - v_bile - v_vldl
        # plasma: VLDL export + peripheral return in, LDLR uptake out
        dL = v_vldl + F_PERIPH - v_up
        C += dt * dC
        L += dt * dL
        C, L = max(C, 1e-9), max(L, 1e-9)
    return C, L, srebp(C)


def pct(x):
    return f"{100 * (x - 1):+.1f}%"


def main():
    print(__doc__.split("Run:")[0].rstrip())
    base_C, base_L, _ = steady_state()
    print("\n" + "=" * 74)
    print("CASE A -- knock down CYP7A1 ALONE (the naive thought experiment)")
    print("=" * 74)
    print(f"\n  {'CYP7A1 left':>12s}{'liver chol':>13s}{'SREBP-2':>11s}"
          f"{'LDLR flux':>12s}{'plasma LDL':>13s}")
    for a in (1.0, 0.8, 0.6, 0.5, 0.3):
        C, L, s = steady_state(a_cyp7a1=a)
        print(f"  {a:>11.0%}{pct(C / base_C):>13s}{pct(s):>11s}"
              f"{pct(s):>12s}{pct(L / base_L):>13s}")
    print("""
  Exactly the textbook chain: less cholesterol leaves as bile acid, the liver
  pool rises, SREBP-2 is repressed, LDL receptors fall, plasma LDL-C rises.
  This is the CYP7A1-deficiency phenotype.""")

    print("\n" + "=" * 74)
    print("CASE B -- knock down CYP7A1 AND HMGCR together (what PFAS did)")
    print("=" * 74)
    print(f"\n  {'both left at':>14s}{'liver chol':>13s}{'SREBP-2':>11s}"
          f"{'plasma LDL':>13s}")
    for a in (1.0, 0.8, 0.6, 0.5, 0.3):
        C, L, s = steady_state(a_cyp7a1=a, a_hmgcr=a)
        print(f"  {a:>13.0%}{pct(C / base_C):>13s}{pct(s):>11s}"
              f"{pct(L / base_L):>13s}")
    print("""
  The sign REVERSES.  Suppressing both arms equally removes more inflow
  (synthesis, 0.45 of the budget) than outflow (bile acids, 0.35), so the
  liver pool falls slightly, SREBP-2 is DEREPRESSED, LDL receptors go UP, and
  plasma LDL-C goes DOWN.  Same CYP7A1 suppression, opposite answer.""")

    print("\n" + "=" * 74)
    print("WHERE THE SIGN FLIPS")
    print("=" * 74)
    print(f"\n  rows: CYP7A1 remaining   cols: HMGCR remaining"
          f"   cell: change in plasma LDL-C\n")
    cols = [1.0, 0.8, 0.6, 0.5, 0.4]
    print("        " + "".join(f"{c:>10.0%}" for c in cols))
    for ac in (1.0, 0.8, 0.6, 0.5, 0.4):
        row = []
        for ah in cols:
            _, L, _ = steady_state(a_cyp7a1=ac, a_hmgcr=ah)
            row.append(pct(L / base_L))
        print(f"  {ac:>5.0%} " + "".join(f"{v:>10s}" for v in row))
    print(f"""
  The diagonal is near zero.  Serum cholesterol rises only BELOW the diagonal
  -- where CYP7A1 is suppressed harder than HMGCR -- and falls above it.
  The crossover is at

      (1 - a_CYP7A1) x {F_BA}  =  (1 - a_HMGCR) x {F_SYN}

  i.e. the CYP7A1 effect has to be about {F_SYN / F_BA:.1f}x deeper than the HMGCR
  effect before the liver pool rises at all.""")

    print("\n" + "=" * 74)
    print("WHAT FRAGKI'S OWN POTENCIES SAY ABOUT WHICH SIDE WE ARE ON")
    print("=" * 74)
    bmd = {"PFOA": (0.7, 0.8), "PFNA": (1.2, 0.8),
           "PFHxS": (3.6, 3.5), "PFOS": (1.2, 1.4)}
    print(f"\n  {'PFAS':8s}{'HMGCR BMD50':>13s}{'CYP7A1 BMD50':>14s}"
          f"{'ratio':>9s}   more sensitive")
    for k, (h, c) in bmd.items():
        print(f"  {k:8s}{h:>13.1f}{c:>14.1f}{h / c:>9.2f}   "
              f"{'CYP7A1' if c < h else 'HMGCR'}")
    print(f"""
  (BMD50 in ng/kg bw/d, scenario 3, Fragki et al. 2022 Table 4.  A LOWER BMD
  means the gene responds at a LOWER dose, i.e. is more sensitive.)

  The two genes are within 1.5x of each other for every compound, and which
  one is more sensitive flips between compounds.  The crossover above needs a
  factor of about {F_SYN / F_BA:.1f}.  So these data sit ON the crossover: they do not
  determine the sign of the serum cholesterol change, in either direction.""")


if __name__ == "__main__":
    main()
