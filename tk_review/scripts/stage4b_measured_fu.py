#!/usr/bin/env python3
"""Stage 4b, with measured unbound fractions for mouse and human.

    python3 -I scripts/stage4b_measured_fu.py

Read-only. Supersedes the sensitivity analysis in
scripts/stage4_fu_sensitivity.py for the two species where fu is now
measured, and corrects an inference made there before the paper arrived.

    Argoul CML, Toutain PL, Berard S, Picard-Hagen N, Gayrard V, Lacroix MZ.
    "Structure-related differences in plasma protein binding of per- and
    polyfluoroalkyl substances in mice and humans." Environ Toxicol Pharmacol
    2025;120:104879. doi:10.1016/j.etap.2025.104879  PMID 41314320.

Equilibrium dialysis on the Dianorm system, chosen over a 96-well dialyser
and over RED after an equilibrium-time study; 15 PFAS, human plasma pooled
from 77 donors and C57BL/6N mouse plasma, two independent experiments. Table
2 gives both experiments separately, which is what is used here: no averaging
where the paper reports a range.

THE INFERENCE THAT WAS WRONG. The roadmap argued, from this paper's abstract,
that interspecies fu differences vanish for long-chain PFAS, so a rat fu near
the human value was the defensible default. Table 2 shows the convergence is
COMPOUND-SPECIFIC and PFOA IS NOT ONE OF THE CONVERGED ONES:

    PFOS    mouse 0.21-0.25%   human 0.16%        1.4x   converged
    PFDA    mouse 0.36-0.42%   human 0.25-0.38%   1.2x   converged
    PFNA    mouse 0.35-0.36%   human 0.13-0.17%   2.4x
    PFOA    mouse 0.87-1.2%    human 0.13-0.17%   6.7x   NOT converged
    PFHxS   mouse 1.1-1.3%     human 0.06-0.08%    17x
    PFBS    mouse 25-29%       human 1.7%          16x

The paper says it plainly: mouse values "were at least 2 times higher than in
humans for all PFAS excepted PFOS and PFDA". PFOA is 6.7-fold apart. So a rat
fu for PFOA cannot be assumed equal to the human one either, and the
break-even argument built on that assumption does not hold.

WHAT THE MOUSE NOW SHOWS, AND IT IS THE POINT OF THE WHOLE SECTION. For mouse
and human the assumption can be removed entirely, because both are measured.
Doing so moves most of the species difference from the tubule to the
glomerulus. That is the circularity demonstrated rather than argued.
"""
import sys

GFR = {"human": 2570, "mouse": 16700, "rat": 14400}   # mL/kg-day, OEHHA A6.4
CLR = {"human": 0.06, "mouse male": 10, "mouse female": 16, "rat male": 18.2}
# Argoul 2025 Table 2, PFOA, the two experiments reported separately
FU_PFOA = {"human": (0.0017, 0.0013), "mouse": (0.0087, 0.012)}
OEHHA = 0.02


def reab(fu, gfr, clr):
    f = fu * gfr
    return f, (1 - clr / f)


def main():
    print("=" * 76)
    print("1. OEHHA's MOUSE ROWS, REDONE WITH MEASURED fu (PFOA)")
    print("=" * 76)
    print(f"  {'basis':34s} {'fu':>8s} {'filtered':>9s} {'male':>8s} "
          f"{'female':>8s}")
    rows = [("OEHHA assumption", OEHHA)]
    rows += [(f"MEASURED, Argoul 2025 exp {i + 1}", v)
             for i, v in enumerate(FU_PFOA["mouse"])]
    for label, fu in rows:
        f, m = reab(fu, GFR["mouse"], CLR["mouse male"])
        _, w = reab(fu, GFR["mouse"], CLR["mouse female"])
        print(f"  {label:34s} {fu:8.4f} {f:9.1f} {100 * m:7.1f}% "
              f"{100 * w:7.1f}%")
    be = CLR["mouse male"] / GFR["mouse"]
    print(f"\n  Mouse break-even is fu = {be:.5f}, which is "
          f"{min(FU_PFOA['mouse']) / be:.0f}-{max(FU_PFOA['mouse']) / be:.0f}"
          f"-fold below")
    print("  the measured values. The mouse is a strong reabsorber of PFOA on")
    print("  any reading, but 93-95% rather than OEHHA's 97%.")

    print("\n" + "=" * 76)
    print("2. THE HUMAN ROW, REDONE THE SAME WAY")
    print("=" * 76)
    for label, fu in [("OEHHA assumption", OEHHA)] + [
            (f"MEASURED, exp {i + 1}", v)
            for i, v in enumerate(FU_PFOA["human"])]:
        f, m = reab(fu, GFR["human"], CLR["human"])
        print(f"  {label:34s} {fu:8.4f} {f:9.2f} {100 * m:7.2f}%")
    print("\n  98.4-98.6% measured against 99.88% assumed. The human")
    print("  conclusion was never at risk.")

    print("\n" + "=" * 76)
    print("3. MOUSE AGAINST HUMAN, WITH AND WITHOUT THE ASSUMPTION")
    print("=" * 76)
    clr_ratio = CLR["mouse male"] / CLR["human"]
    print(f"  Observed: the male mouse clears PFOA renally "
          f"{clr_ratio:.0f}x faster.")
    print("  That is a fact and does not depend on fu. The split does:\n")
    print(f"  {'basis':30s} {'filtered load':>14s} {'tubular escape':>15s}")
    fh_a, fm_a = OEHHA * GFR["human"], OEHHA * GFR["mouse"]
    print(f"  {'fu assumed equal (OEHHA)':30s} {fm_a / fh_a:13.1f}x "
          f"{clr_ratio / (fm_a / fh_a):14.1f}x")
    for i in range(2):
        fh = FU_PFOA["human"][i] * GFR["human"]
        fm = FU_PFOA["mouse"][i] * GFR["mouse"]
        print(f"  {'MEASURED, experiment ' + str(i + 1):30s} "
              f"{fm / fh:13.1f}x {clr_ratio / (fm / fh):14.1f}x")
    print()
    print("  Removing the assumption moves the mouse-human difference from")
    print("  mostly tubular to mostly GLOMERULAR. Under OEHHA's equal fu,")
    print(f"  {clr_ratio / (fm_a / fh_a):.0f}-fold of the {clr_ratio:.0f} is "
          f"attributed to the tubule; with the")
    print("  measured values it is about 4-fold, and the filtered load")
    print("  carries the rest. The conclusion that renal reabsorption")
    print("  explains the species difference does not survive for the mouse.")
    print("\n  This is no longer a hypothetical. It is the same arithmetic")
    print("  OEHHA applies, with two measured numbers in place of one")
    print("  assumed one, and it reverses which compartment the difference")
    print("  lives in.")

    print("\n" + "=" * 76)
    print("4. WHAT THIS IMPLIES FOR THE RAT, HONESTLY")
    print("=" * 76)
    print(f"  {'assumed rat fu':34s} {'filtered':>9s} {'reabsorbed':>12s}")
    for label, fu in (("OEHHA's 0.02", OEHHA),
                      ("= the mouse's measured PFOA fu", 0.010),
                      ("= the human's measured PFOA fu", 0.0015),
                      ("break-even", CLR["rat male"] / GFR["rat"])):
        f, m = reab(fu, GFR["rat"], CLR["rat male"])
        tag = f"{100 * m:11.1f}%" if m > 1e-9 else "        0.0%"
        print(f"  {label:34s} {f:9.1f} {tag}")
    print()
    print("  The rat is not measured, and the honest range is wide. If the")
    print("  rat resembles the MOUSE -- the closer rodent -- it reabsorbs")
    print("  about 87% of filtered PFOA, lower than OEHHA's 93.7% but still")
    print("  clearly a reabsorber, and the break-even argument in the")
    print("  previous version of this section was too strong. If it")
    print("  resembles the HUMAN it sits near break-even. Nothing in this")
    print("  paper decides between those, because PFOA is precisely the")
    print("  compound whose fu does NOT converge across the two species")
    print("  that were measured.")
    print()
    print("  So the request stands and is now better motivated: a rat plasma")
    print("  fu for PFOA, by this same Dianorm method, from this same group.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
