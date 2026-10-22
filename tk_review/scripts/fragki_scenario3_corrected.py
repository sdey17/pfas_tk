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
    "PFHxS": dict(HMGCR=3.6, LSS=2.4, CYP7A1=3.5, CYP8B1=5.2, CXCL10=1.0,
                  PKD4=.5, ANGPTL4=1.1, OAT5=.7, ATF4=4.2, SLC7A11=3.4,
                  YARS1=3.8, THRSP=3.0),
    "PFOS":  dict(HMGCR=1.2, LSS=.3, CYP7A1=1.4, CYP8B1=1.7, CXCL10=1.1,
                  PKD4=2.0, ANGPTL4=2.2, OAT5=.2, ATF4=1.9, SLC7A11=1.3,
                  YARS1=1.9, THRSP=1.0),
}
TG = {"PFOA": 0.8, "PFNA": 3.2, "PFHxS": 3.0, "PFOS": 1.5}
EFSA = {"PFOA": dict(LB=0.18, UB=4.18), "PFNA": dict(LB=0.04, UB=3.68),
        "PFHxS": dict(LB=0.08, UB=3.45), "PFOS": dict(LB=0.58, UB=4.47)}
NCF = {"PFOA": 7, "PFNA": 8, "PFHxS": 6, "PFOS": 8}
ORDER = ("PFOA", "PFNA", "PFHxS", "PFOS")
# log K_OW for the Armitage limb.  Carboxylates from Q.COMPOUNDS / the OPERA
# fill; sulfonates from Qin Table 1, with COSMOtherm carried where it differs.
LKOW = {"PFOA": (5.68, None), "PFNA": (None, None),
        "PFHxS": (3.69, 5.29), "PFOS": (5.77, None)}


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def logkow(cp):
    """Primary log K_OW, and the alternative where the sources disagree."""
    pri, alt = LKOW[cp]
    if pri is None:
        pri = Q.KOW_OPERA(NCF[cp])          # carboxylate OPERA fill
    return pri, alt


def factors():
    """Model-predicted / Fragki-measured fraction of dose in cells."""
    a = P.ASSAYS["Fragki"]
    V = P.volumes(a)
    out = {}
    for cp in ORDER:
        lk, alt = logkow(cp)
        q, ar, ar2 = [], [], []
        for d in a["doses"]:
            m = P.measured_x_cell("Fragki", cp, d, V, a)
            q.append(P.qin_run(cp, d, V)["x_cell"] / m)
            ar.append(P.armitage_run(cp, lk, d, V)["x_cell"] / m)
            if alt is not None:
                ar2.append(P.armitage_run(cp, alt, d, V)["x_cell"] / m)
        out[cp] = dict(qin=gm(q), arm=gm(ar),
                       qin_rng=(min(q), max(q)), arm_rng=(min(ar), max(ar)),
                       arm_alt=gm(ar2) if ar2 else None,
                       qin_src=Q.qin_params(cp, NCF[cp])["source"])
    return out


def main():
    print(__doc__.split("Run:")[0].rstrip())
    F = factors()

    print("\n" + "=" * 76)
    print("THE CORRECTION FACTORS, ON FRAGKI'S OWN REPORTED UPTAKE")
    print("=" * 76)
    print(f"\n  {'PFAS':7s}{'Armitage':>11s}{'range':>16s}{'Qin':>9s}"
          f"{'range':>16s}   Qin params")
    for cp in ORDER:
        f = F[cp]
        ar_r = "%.1f-%.1fx" % f["arm_rng"]
        q_r = "%.1f-%.1fx" % f["qin_rng"]
        print(f"  {cp:7s}{f['arm']:>10.1f}x{ar_r:>16s}"
              f"{f['qin']:>8.1f}x{q_r:>16s}   {f['qin_src']}")
    alt = F["PFHxS"]["arm_alt"]
    print(f"""
  PFHxS carries a second, larger uncertainty that the others do not.  Its
  log K_OW is 3.69 by OPERA and 5.29 by COSMOtherm (Qin Table 1), a factor of
  40, and the Armitage limb rides on it directly: the factor above is
  {F['PFHxS']['arm']:.1f}x on OPERA and {alt:.1f}x on COSMOtherm.  Nothing in this
  analysis can choose between them.""")
    print("""
  The two models do NOT agree about what needs correcting, and the split is
  by head group, not by chain length:

    Qin over-predicts all four, by 4.4x to 14.5x.
    Armitage over-predicts the two CARBOXYLATES by 6-7x and is essentially
    unbiased on the two SULFONATES (0.9x, 1.1x).

  The reason is visible in the raw numbers.  Fragki measured roughly seven
  times more PFOS than PFOA in the cells (3.9-10.9% of dose against
  0.6-1.6%), while Armitage predicts almost the same fraction for both (7.4%
  against 6.6%).  Its log K_OW-driven prediction does not separate the head
  groups nearly as much as the measurement does, so it happens to land on the
  sulfonates and misses the carboxylates.  That is a coincidence of level,
  not evidence that the sulfonate case is better modelled.

  The spread across the dose series is the models' dose-dependence, not
  scatter in the measurement: Armitage is linear, so its predicted FRACTION
  is dose-independent while the measured fraction rises with dose, and its
  factor therefore falls as dose rises.""")

    print("\n" + "=" * 76)
    print("SCENARIO-3 EFFECT DOSES, AS PUBLISHED AND AS CORRECTED (ng/kg bw/d)")
    print("=" * 76)
    for cp in ORDER:
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
    for cp in ORDER:
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
    tot = {}
    for k in ("pub", "arm", "qin"):
        tot[k] = sum(
            1 for cp in ORDER for x in GENES[cp].values()
            if x * (1.0 if k == "pub" else F[cp][k]) < EFSA[cp]["UB"])
    n = sum(len(GENES[cp]) for cp in ORDER)
    print(f"""
  Across all four compounds, genes whose effect dose falls BELOW the UB mean
  dietary exposure -- the exceedance that carries the paper's conclusion:

      as published      {tot['pub']:>2d} / {n}
      x Armitage        {tot['arm']:>2d} / {n}
      x Qin             {tot['qin']:>2d} / {n}

  The two models do not agree on what the correction means.  Under Qin's
  factors the exceedance largely collapses.  Under Armitage's it survives,
  because Armitage needs no correction on the sulfonates -- and for PFHxS its
  factor is slightly BELOW one, so the correction moves that compound the
  other way, from 7/12 genes to 10/12.

  So the honest summary is not "the effect doses go up".  It is that the in
  vitro dosimetry step carries a 0.9x to 14.5x uncertainty whose size and
  even whose DIRECTION depend on which mass-balance model is believed, and
  the two models disagree most on exactly the two compounds (PFHxS, PFOS) for
  which Qin's own parameters are fully measured and need no interpolation.
  That is a strong argument for measuring the cellular concentration rather
  than modelling it -- which is what Fragki's scenario 3 already does, and is
  the reason it remains the most defensible of their three scenarios.""")


if __name__ == "__main__":
    main()
