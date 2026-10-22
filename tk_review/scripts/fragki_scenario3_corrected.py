#!/usr/bin/env python3
"""What do Fragki's scenario-3 effect doses become under the mass-balance models?

THE CHAIN, AND WHY IT IS LINEAR
  Fragki et al. 2022 scenario 3 takes the MEASURED cell-associated
  concentration in HepaRG as the human liver concentration, then runs a PBK
  model backwards to the oral dose that would produce it.  So

      oral equivalent effect dose  proportional to  assumed C_liver at the POD

  at doses far below saturation of renal resorption, which ng/kg bw/d
  certainly is.  Nothing else in the chain moves: the in vitro
  concentration-response curve, the BMD analysis and the PBK model are
  untouched.  So a factor F on the cellular concentration is a factor F on
  the effect dose, and the comparison with EFSA exposure shifts by F.

  DIRECTION.  If the cells really hold MORE than was measured, the same in
  vitro POD corresponds to a HIGHER liver concentration, which takes a HIGHER
  oral dose to reach.  Effect doses go UP and the margin against dietary
  exposure WIDENS.  Correcting upward is reassuring, not alarming.

WHAT F IS
  Both mass-balance models predict a larger fraction of dose in the cells
  than Fragki measured.  F is that ratio, per compound, over Fragki's own
  dose series.  This is emphatically a 'if the model is right' sensitivity
  analysis: the review's own conclusion is that the models are wrong in ways
  it localises, and no result here says the measurement is at fault.  The
  point is to size the uncertainty, not to assign blame.

Run:  python3 tk_review/scripts/fragki_scenario3_corrected.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P      # noqa: E402
import mbm_qin_vs_armitage as Q          # noqa: E402

# Fragki Table 4, scenario 3: BMD50 oral equivalent effect doses, ng/kg bw/d
GENES = {
    "PFOA":  dict(HMGCR=.7, LSS=.1, CYP7A1=.8, CYP8B1=.4, CXCL10=.2, PKD4=.06,
                  ANGPTL4=.09, OAT5=.08, ATF4=.6, SLC7A11=.5, YARS1=.5, THRSP=.5),
    "PFNA":  dict(HMGCR=1.2, LSS=.8, CYP7A1=.8, CYP8B1=1.0, CXCL10=.7, PKD4=.4,
                  ANGPTL4=.7, OAT5=.2, ATF4=1.4, SLC7A11=.9, YARS1=1.2, THRSP=.7),
}
TG = {"PFOA": 0.8, "PFNA": 3.2}
EFSA = {"PFOA": dict(LB=0.18, UB=4.18), "PFNA": dict(LB=0.04, UB=3.68)}
NCF = {"PFOA": 7, "PFNA": 8}


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def factors():
    """Model-predicted / Fragki-measured fraction of dose in cells."""
    a = P.ASSAYS["Fragki"]
    V = P.volumes(a)
    out = {}
    for cp, nCF in NCF.items():
        lk = 5.68 if cp == "PFOA" else Q.KOW_OPERA(nCF)
        q, ar = [], []
        for d in a["doses"]:
            m = P.measured_x_cell("Fragki", cp, d, V, a)
            q.append(P.qin_run(cp, d, V)["x_cell"] / m)
            ar.append(P.armitage_run(cp, lk, d, V)["x_cell"] / m)
        out[cp] = dict(qin=gm(q), arm=gm(ar),
                       qin_rng=(min(q), max(q)), arm_rng=(min(ar), max(ar)))
    return out


def main():
    print(__doc__.split("Run:")[0].rstrip())
    F = factors()

    print("\n" + "=" * 76)
    print("THE CORRECTION FACTORS, ON FRAGKI'S OWN REPORTED UPTAKE")
    print("=" * 76)
    print(f"\n  {'PFAS':7s}{'Armitage':>11s}{'range':>16s}{'Qin':>9s}{'range':>16s}")
    for cp, f in F.items():
        ar_r = "%.1f-%.1fx" % f["arm_rng"]
        q_r = "%.1f-%.1fx" % f["qin_rng"]
        print(f"  {cp:7s}{f['arm']:>10.1f}x{ar_r:>16s}"
              f"{f['qin']:>8.1f}x{q_r:>16s}")
    print("""
  Both models say the cells hold several times more than was measured.  The
  spread across the dose series is the models' dose-dependence, not scatter
  in the measurement: Armitage is linear so its predicted FRACTION is
  dose-independent, while the measured fraction rises with dose, so its
  factor falls as dose rises.""")

    print("\n" + "=" * 76)
    print("SCENARIO-3 EFFECT DOSES, AS PUBLISHED AND AS CORRECTED (ng/kg bw/d)")
    print("=" * 76)
    for cp in ("PFOA", "PFNA"):
        f = F[cp]
        v = sorted(GENES[cp].values())
        lo, hi = v[0], v[-1]
        e = EFSA[cp]
        print(f"\n  {cp}   EFSA mean dietary exposure  LB {e['LB']}  UB {e['UB']}")
        print(f"    {'parameterisation':28s}{'gene effect doses':>22s}"
              f"{'TG':>8s}{'genes under UB':>17s}")
        rows = [("Fragki as published", 1.0),
                (f"x Armitage factor ({f['arm']:.1f}x)", f["arm"]),
                (f"x Qin factor ({f['qin']:.1f}x)", f["qin"])]
        for lab, k in rows:
            n_ub = sum(1 for x in v if x * k < e["UB"])
            print(f"    {lab:28s}{f'{lo*k:.2f} - {hi*k:.2f}':>22s}"
                  f"{TG[cp]*k:>8.1f}{f'{n_ub}/12':>17s}")

    print("\n" + "=" * 76)
    print("WHAT IT DOES TO THE HEADLINE COMPARISON")
    print("=" * 76)
    print(f"\n  {'PFAS':6s}{'bound':7s}{'as published':>26s}{'x Armitage':>22s}"
          f"{'x Qin':>22s}")
    for cp in ("PFOA", "PFNA"):
        f, v, e = F[cp], sorted(GENES[cp].values()), EFSA[cp]
        for b in ("LB", "UB"):
            cells = []
            for k in (1.0, f["arm"], f["qin"]):
                lo, hi = v[0] * k, v[-1] * k
                x = e[b]
                if lo <= x <= hi:
                    cells.append("INSIDE the span")
                elif x < lo:
                    cells.append(f"below by {lo/x:.0f}x")
                else:
                    cells.append(f"above by {x/hi:.0f}x")
            print(f"  {cp:6s}{b:7s}{cells[0]:>26s}{cells[1]:>22s}"
                  f"{cells[2]:>22s}")
    print("""
  Both corrections push the effect doses UP, away from dietary exposure.
  The LB bound, which sat inside the PFOA gene span as published, drops below
  it.  The UB bound, which sat above the span for both compounds, moves into
  it or below it -- so the paper's 'in the same range, in particular for UB'
  becomes MORE defensible as a statement about UB, while the exceedance that
  drives the conclusion gets weaker.""")


if __name__ == "__main__":
    main()
