#!/usr/bin/env python3
"""What changed in the IV-MBM cell term between v1.0 (2014) and v2.0 (2021)?

The 2014 paper (Armitage, Wania & Arnot, Environ Sci Technol 48:9770-9779) is
now in hand.  Its Eq. 1 is

  C_W = M_T / (K_AW*V_A + V_W + K_SaW*V_Sa + K_SlW*V_Sl + K_DW*V_D + K_CW*V_C)

and the cell term is ONE lumped phase:

  "In the absence of user-specified values, K_CW is estimated as K_OW
   multiplied by the volume fraction of lipid equivalents (i.e., pseudo-octanol
   content) (f_L)"                                         -- citing Mackay 2001

with f_L for cells/tissue "assumed to be 0.05 for the results displayed" and
varied over 0.01-0.05 in the sensitivity analysis.  Three things follow, and
they settle questions this project had left open:

  1. BULK CELL PROTEIN HAS NO SEPARATE TERM IN v1.0.  It is not that the 3.5%
     octanol-equivalence convention is documented there with a different
     value; there is no protein term at all, and no membrane-lipid term and no
     cell-water term either.  The whole cell is K_OW * f_L.  So v2.0's
     "expanded to include ... bulk proteins" is literally a new phase, and the
     coefficient attached to it cannot be recovered from the 2014 record.
     The 3.5% does not appear in the 2014 paper (its only "3.5%" is a
     supersaturation statistic).

  2. v1.0 EXPLICITLY EXCLUDES PFAS.  "The model is intended for neutral
     organic chemicals and should only be applied to ionogenic organic
     chemicals (IOCs) if the compound is predominantly neutral at the test
     pH."  The paper's own future-work list names the gap: "An additional task
     is to expand the approach for neutral organic chemicals to ionogenic
     organic chemicals that are predominantly charged in the test system."
     v2.0's IOC machinery -- the d scaling factors, anionic phospholipids,
     lysosomal trapping -- is that extension.  Anionic PFAS are the first
     application of it, not a case the formalism was calibrated on.

  3. THE TWO CELL TERMS DISAGREE BY A LARGE FACTOR FOR AN ANION, and the
     direction is not obvious.  This script computes it.

Run:  python3 tk_review/scripts/mbm_armitage_v1_vs_v2.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P      # noqa: E402
import mbm_qin_vs_armitage as Q          # noqa: E402

FL_V1 = (0.01, 0.05)        # the range the 2014 paper uses; 0.05 is its default


def main():
    print(__doc__.split("Run:")[0].rstrip())
    comp = P.CELL_COMP["AREc2"]
    print("\n" + "=" * 78)
    print("CELL-WATER DISTRIBUTION RATIO, v1.0 AGAINST v2.0")
    print("=" * 78)
    print(f"""
  v1.0   K_CW = K_OW * f_L                                   (f_L = 0.05)
  v2.0   D_CW = sl*D_OW + ml*D_MW + 0.035*prot*D_OW + water
         sl {comp['sl']}  ml {comp['ml']}  prot {comp['prot']}  water {comp['water']}
         (AREc2, Armitage 2021 SI Table S5)

  Both are evaluated on the ANION at pH 7.4, i.e. with v2.0's own scaling
  factors d_OW = 3.1 and d_MW = 1.0.  v1.0 has no scaling factors, so its
  K_OW is the neutral value -- which is the point of item 2 above: applying
  v1.0 to an anion at all means ignoring speciation.  Both readings are shown.
""")
    hdr = (f"  {'cmpd':8s}{'nCF':>4s}{'logKow':>8s} | {'v1 f=.05':>10s}"
           f"{'v1 anion':>10s} | {'v2 D_CW':>10s} | {'v2/v1n':>9s}{'v2/v1a':>9s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
        if nCF > 8:
            continue
        lk = ko if ko is not None else Q.KOW_OPERA(nCF)
        D = P.armitage_D(lk, comp)
        v1_neutral = 10 ** lk * FL_V1[1]
        v1_anion = 10 ** (lk - P.SI["d_OW"]) * FL_V1[1]
        print(f"  {nm:8s}{nCF:>4d}{lk:>8.2f} | {v1_neutral:>10.3g}"
              f"{v1_anion:>10.3g} | {D['Dcw']:>10.3g} |"
              f" {D['Dcw'] / v1_neutral:>9.3f}{D['Dcw'] / v1_anion:>9.1f}")

    print("""
  Read the last column.  For an anion, v2.0's cell sink is 1-2 orders of
  magnitude LARGER than v1.0's would be at the same f_L, and the whole of that
  difference is the membrane-lipid phase: octanol carries an anion penalty of
  3.1 log units in v2.0 and the membrane carries only 1.0, so the phase that
  v1.0 did not have is also the phase that is least penalised by the charge.
  Against the NEUTRAL reading of v1.0, v2.0's cell is smaller -- the familiar
  direction, and the one that would be inferred from f_L alone.""")

    print("\n" + "=" * 78)
    print("WHAT CARRIES THE v2.0 CELL TERM")
    print("=" * 78)
    print(f"\n  {'cmpd':8s}{'nCF':>4s} | {'stor.lip':>9s}{'mem.lip':>9s}"
          f"{'protein':>9s}{'water':>9s}")
    for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
        if nCF > 8:
            continue
        lk = ko if ko is not None else Q.KOW_OPERA(nCF)
        D = P.armitage_D(lk, comp)
        parts = [comp["sl"] * D["Dow"], comp["ml"] * D["Dmw"],
                 P.P_NLOM * comp["prot"] * D["Dow"], comp["water"]]
        tot = sum(parts)
        print(f"  {nm:8s}{nCF:>4d} | " +
              "".join(f"{100 * p / tot:>8.2f}%" for p in parts))
    print("""
  The undocumented coefficient (0.035, applied to 5.14% protein) contributes
  at most 0.25% of the v2.0 cell term AT ITS OWN VALUE, because the
  membrane-lipid phase dominates from nCF 5 upward.  So the published model's
  output is insensitive to that number over any plausible range, which bounds
  the consequence of its being undocumented.

  That is NOT the same as saying the protein phase is unimportant.  It stands
  in for a modelling choice -- octanol-equivalence -- and replacing that choice
  with a measured structural-protein affinity has large leverage:""")
    print(f"  {'cmpd':8s}{'nCF':>4s} | {'mem.lip':>10s}{'prot(.035)':>12s}"
          f"{'prot(meas)':>12s} | {'D_CW x':>8s}")
    for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
        if nCF > 8:
            continue
        lk = ko if ko is not None else Q.KOW_OPERA(nCF)
        D = P.armitage_D(lk, comp)
        mem = comp["ml"] * D["Dmw"]
        pr_oct = P.P_NLOM * comp["prot"] * D["Dow"]
        dns = 10 ** Q.qin_params(nm, nCF)["logDnonspec"]
        pr_meas = comp["prot"] * dns
        sub = (D["Dcw"] - pr_oct + pr_meas) / D["Dcw"]
        print(f"  {nm:8s}{nCF:>4d} | {mem:>10.3g}{pr_oct:>12.3g}"
              f"{pr_meas:>12.3g} | {sub:>7.2f}x")
    print("""
  Substituting Qin's measured non-specific protein-water ratio for the octanol
  surrogate makes the protein term comparable to the membrane term and enlarges
  the cell sink, which is why that substitution moves the aggregate bias so far
  in the review's parameter-substitution table.  The precise statement is that
  the model is insensitive to the VALUE of the undocumented coefficient and
  highly sensitive to the CHOICE it encodes.""")


if __name__ == "__main__":
    main()
