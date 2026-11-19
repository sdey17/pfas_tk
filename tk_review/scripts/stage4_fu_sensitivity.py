#!/usr/bin/env python3
"""Stage 4b: the species comparison WITHOUT assuming a constant unbound fraction.

    python3 -I scripts/stage4_fu_sensitivity.py

Read-only.

THE OBJECTION. OEHHA PHG Table A6.4, and Han et al. 2012 before it, computes
net renal tubular reabsorption as fu x GFR - CLR with fu ASSUMED TO BE 0.02
FOR EVERY SPECIES. The PHG states it in one clause -- "Assuming the unbound
fraction in blood (fu) as 0.02" -- with no citation and no species
justification, and then draws a headline conclusion from it: "According to
this calculation, human kidney reabsorbed 99.8% of PFOA. This mechanism would
explain the exceptionally long human T1/2 relative to animals."

The reviewer's objection is that a single fu cannot be right across species.
It is worse than that. A single fu is not right across COMPOUNDS within one
species. Argoul 2026 Table 1 measures human plasma fu for nine PFAS and it
spans fifty-two-fold, from 0.08% for PFHxS to 4.2% for GenX. The assumed 0.02
is close to the measured value for the SHORT-chain compounds (PFHxA 2.0%,
PFBS 1.7%) and 12 to 25-fold too high for the long-chain ones -- PFOA, PFOS,
PFNA, PFHxS -- which are the compounds the PHG is about.

THE CIRCULARITY. With fu held equal between species, the ratio of filtered
loads is forced to equal the ratio of GFRs, and every remaining part of the
species difference is attributed to the tubule by construction. The
assumption does not merely add uncertainty to the conclusion; in part it
produces it.

WHAT SURVIVES AND WHAT DOES NOT, below. The human conclusion is robust: human
renal clearance of PFOA is so small that essentially any plausible fu leaves
reabsorption above 98%. The RAT conclusion is not robust at all, and that is
the finding.
"""
import sys

GFR = {"human": 2570, "rat": 14400, "mouse": 16700}   # mL/kg-day, OEHHA A6.4
CLR = {"human (OEHHA)": ("human", 0.06),
       "human (Argoul)": ("human", 0.042),
       "rat male (Kemper)": ("rat", 18.2),
       "rat male (Ohmori)": ("rat", 21.2),
       "rat female (Kemper)": ("rat", 666),
       "mouse male": ("mouse", 10)}
# Argoul 2026 Table 1, measured human plasma unbound fraction
FU_HUMAN = {"PFHxS": 0.0008, "PFOS": 0.0016, "PFOA": 0.0017, "PFNA": 0.0017,
            "PFDA": 0.0038, "PFO2OA": 0.0047, "PFBS": 0.017, "PFHxA": 0.020,
            "HFPO-DA (GenX)": 0.042}


def reab(fu, gfr, clr):
    f = fu * gfr
    return 1 - clr / f, f


def main():
    print("=" * 76)
    print("1. THE ASSUMED fu IS A SHORT-CHAIN VALUE USED ON LONG-CHAIN PFAS")
    print("=" * 76)
    print(f"  {'PFAS':16s} {'measured human fu':>18s} {'OEHHA 0.02 is':>16s}")
    for c, fu in sorted(FU_HUMAN.items(), key=lambda x: x[1]):
        print(f"  {c:16s} {fu:18.5f} {0.02 / fu:14.1f}x high")
    print(f"\n  Range within one species: "
          f"{max(FU_HUMAN.values()) / min(FU_HUMAN.values()):.0f}-fold.")

    print("\n" + "=" * 76)
    print("2. THE HUMAN CONCLUSION SURVIVES. PFOA, varying fu over 100-fold:")
    print("=" * 76)
    print(f"  {'fu':>8s} {'filtered':>9s} {'reabsorbed':>11s}   source of this fu")
    for fu, src in ((0.02, "OEHHA assumption"),
                    (0.0047, "measured, PFO2OA"),
                    (0.0017, "MEASURED, PFOA (Argoul 2026)"),
                    (0.0008, "measured, PFHxS - the lowest seen"),
                    (0.0002, "4x below anything measured")):
        pct, f = reab(fu, GFR["human"], 0.06)
        print(f"  {fu:8.4f} {f:9.2f} {100 * pct:10.2f}%   {src}")
    print("\n  Human renal clearance of PFOA is 0.06 mL/kg-day. Reabsorption")
    print("  stays above 98% for every fu down to 0.0008, and the break-even")
    print(f"  point -- where the human stops reabsorbing -- is fu = "
          f"{0.06 / GFR['human']:.6f},")
    print("  thirty times below the lowest value ever measured. The 99%")
    print("  figure is not an artefact of the assumption.")

    print("\n" + "=" * 76)
    print("3. THE RAT CONCLUSION DOES NOT SURVIVE. Male rat, CL_R 18.2:")
    print("=" * 76)
    print(f"  {'fu':>8s} {'filtered':>9s} {'reabsorbed':>12s}")
    be = 18.2 / GFR["rat"]
    for fu in sorted((0.02, 0.01, 0.005, 0.0025, be, 0.0017, 0.0008),
                     reverse=True):
        pct, f = reab(fu, GFR["rat"], 18.2)
        tag = (f"{100 * pct:11.1f}%" if pct > 1e-9 else
               "        0.0%" if pct > -1e-9 else "  SECRETION")
        extra = ""
        if abs(fu - be) < 1e-9:
            extra = "  <- break-even: no net reabsorption"
        elif fu == 0.0017:
            extra = "  <- if the rat's fu matched the human's"
        elif fu == 0.02:
            extra = "  <- the OEHHA assumption"
        print(f"  {fu:8.5f} {f:9.1f} {tag}{extra}")
    print(f"\n  The male rat stops being a net reabsorber at fu = {be:.5f}.")
    print("  That is only 1.4-fold below the human's MEASURED PFOA fu of")
    print("  0.0017. So whether the male rat reabsorbs PFOA at all, or")
    print("  secretes it like the female, turns on a rat fu that nobody in")
    print("  this database has measured -- and the whole 93.7% rests on an")
    print("  assumed value 16-fold above that break-even.")

    print("\n" + "=" * 76)
    print("4. THE DECOMPOSITION IS FORCED BY THE ASSUMPTION")
    print("=" * 76)
    print("  The observed 303-fold difference in renal clearance is a fact and")
    print("  does not depend on fu. How it SPLITS does:")
    print()
    print(f"  {'fu_rat / fu_human':>18s} {'filtered load':>14s} "
          f"{'tubular escape':>15s}")
    for ratio in (0.5, 1.0, 2.0, 5.0, 11.8):
        filt = (GFR["rat"] / GFR["human"]) * ratio
        print(f"  {ratio:18.1f} {filt:13.1f}x {303 / filt:14.1f}x")
    print()
    print("  OEHHA sets the first column to 1.0 by assumption, which fixes")
    print("  the filtered-load term at 5.6x and sends the remaining 54x to")
    r = 11.8
    filt = (GFR["rat"] / GFR["human"]) * r
    print(f"  the tubule. If the rat's fu is instead {r}-fold the human's --")
    print("  which is what you get if the rat really is near 0.02 while the")
    print(f"  human is measured at 0.0017 -- then {filt:.0f}x of the 303x is")
    print(f"  GLOMERULAR and only {303 / filt:.1f}x is tubular. The conclusion")
    print("  that renal reabsorption explains the human half-life would not")
    print("  follow.")
    print()
    print("  This is the circularity: assuming equal fu does not merely add")
    print("  uncertainty, it assigns the difference to the tubule.")

    print("\n" + "=" * 76)
    print("5. WHAT TO DO ABOUT IT")
    print("=" * 76)
    print("  - Quote the human numbers with the MEASURED fu, per compound,")
    print("    from Argoul 2026 Table 1. They are already in")
    print("    db/primary_2026/argoul2026_table1_human_reabsorption.csv and")
    print("    they give human PFOA 99.0% reabsorbed rather than 99.8%.")
    print("  - Do NOT quote the rat, mouse, dog, macaque or rabbit")
    print("    percentages from Table A6.4 without the assumption attached.")
    print("    They are conditional on fu = 0.02 and the rat's is close to a")
    print("    sign change.")
    print("  - The single measurement that would fix this is a rat plasma fu")
    print("    for PFOA, by the method Argoul used for human and mouse.")
    print("  - Until then the honest statement of aim 1 is: the human")
    print("    reabsorbs PFOA almost completely, which is robust; the")
    print("    species CONTRAST in reabsorption is not established, because")
    print("    the rat side of it rests on an assumed binding value.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
