#!/usr/bin/env python3
"""Does liver out-compete plasma for PFAS? Affinity cannot answer it.

A tighter Kd on L-FABP than on albumin does NOT mean the liver holds more.
The mass a phase holds at the trace limit is

    bound per litre of phase = C_protein x n x Cw / Kd

so the comparison between two compartments is a ratio of C x n / Kd, and
the concentration term is the one no binding paper reports. Until this
file, the database had every Kd and no C, which is why the question could
be asked and not answered.

WHAT MAKES THE COMPARISON CLEAN. Two rules have to hold at once, and
until now nothing in the database satisfied both:

  ONE METHOD. Khazaee et al. 2021 measure L-FABP by equilibrium dialysis
  and Chen et al. 2025 measure albumin by equilibrium dialysis, but they
  are different laboratories, and part 8 of logk_converted_vs_measured.py
  shows two dialysis laboratories differ on the SAME albumin Kd by a
  median 1.45 log. A cross-laboratory ratio is therefore uninterpretable
  however well the methods match.

  ONE PAPER. Jia et al. 2022 measures BOTH proteins, by microscale
  thermophoresis, in one table. Five compounds carry an L-FABP and an HSA
  constant from that single run. Those five are the whole of the usable
  evidence, and they are all alternative PFAS -- the legacy compounds in
  Jia's table take their HSA values from Gao et al. 2019 and so cannot be
  used here.

THE CONCENTRATIONS, from db/protein_abundance.csv:
  FABP1   399 uM in human liver cytosol (211-760 over 61 livers),
          Wen & Isoherranen 2026, LC-MS/MS
  albumin 527-752 uM in plasma, the 35-50 g/L clinical interval

WHAT IS ASSUMED, stated once. The ratio carries n_albumin / n_LFABP, and
MST reports no stoichiometry for either protein, so the headline figures
take them equal and the sensitivity to that is computed explicitly. It is
the weakest step and it is not hidden: Maso et al. 2021 co-crystallised
four PFOA per albumin, and if albumin really carries four sites to
L-FABP's one, every ratio below falls four-fold.

Run:  python3 tk_review/scripts/lfabp_capacity_vs_albumin.py
"""
import csv
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
ABUND = os.path.join(os.path.dirname(HERE), "db", "protein_abundance.csv")

JIA_DOI = "10.1021/acs.est.1c08493"
F_ICW = 0.573          # liver fraction of intracellular water, Wen 2026


def rule(ch="=", n=78):
    print(ch * n)


def load_abundance():
    out = {}
    for r in csv.DictReader(open(ABUND, encoding="utf-8")):
        out[(r["protein"].split("(")[0].strip(), r["tissue"])] = r
    return out


def load_jia():
    """Jia's own MST constants, L-FABP and HSA, keyed by compound."""
    lf, hsa = {}, {}
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if r["doi"] != JIA_DOI or not r["kd"].strip():
            continue
        try:
            v = float(r["kd"].split()[0])
        except ValueError:
            continue
        p = r["protein"].lower()
        if "fatty acid" in p:
            lf[r["chemical"]] = v
        elif "serum albumin" in p:
            hsa[r["chemical"]] = v
    return lf, hsa


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    ab = load_abundance()
    lf, hsa = load_jia()
    shared = sorted(set(lf) & set(hsa), key=lambda c: lf[c])

    fab = ab[("liver fatty acid-binding protein", "liver")]
    alb = ab[("serum albumin", "plasma")]
    C_FAB = float(fab["value_uM"])
    C_FAB_LO, C_FAB_HI = float(fab["low_uM"]), float(fab["high_uM"])
    C_ALB_LO, C_ALB_HI = float(alb["low_uM"]), float(alb["high_uM"])
    C_ALB = math.sqrt(C_ALB_LO * C_ALB_HI)

    rule()
    print("1. THE TWO CONCENTRATIONS, WHICH ARE THE NEW INFORMATION")
    rule()
    print(f"  FABP1 in human liver cytosol   {C_FAB:.0f} uM   "
          f"({C_FAB_LO:.0f}-{C_FAB_HI:.0f}, CV {fab['cv_percent']}%, "
          f"n = {fab['n_subjects']} livers)")
    print(f"  albumin in human plasma        {C_ALB:.0f} uM   "
          f"({C_ALB_LO:.0f}-{C_ALB_HI:.0f}, the 35-50 g/L interval)")
    print()
    print("  THAT IS THE FIRST SURPRISE. These are the same order. The")
    print("  intuition that a nuclear-receptor-like protein is too scarce")
    print("  to matter does not apply to L-FABP at all: per litre of")
    print(f"  cytosol it is {C_FAB / C_ALB:.2f}x the molarity of albumin "
          f"per litre of plasma.")
    print()
    print(f"  Per litre of LIVER rather than of cytosol, multiply by the")
    print(f"  intracellular water fraction {F_ICW}: {C_FAB * F_ICW:.0f} uM.")
    print("  Both bases are carried below, because which one is right")
    print("  depends on whether the question is about a cell or an organ.")
    print()

    rule()
    print("2. THE AFFINITIES, FROM ONE PAPER AND ONE RUN")
    rule()
    print("  Jia et al. 2022, microscale thermophoresis, 25 C. Five")
    print("  compounds carry both constants from this paper's own hands.")
    print()
    print(f"  {'compound':16s}{'L-FABP Kd':>11s}{'HSA Kd':>10s}"
          f"{'HSA/L-FABP':>12s}   which binds tighter")
    for c in shared:
        r = hsa[c] / lf[c]
        print(f"  {c:16s}{lf[c]:>11.2f}{hsa[c]:>10.1f}{r:>12.1f}   "
              + ("L-FABP" if r > 1 else "albumin"))
    print()
    tighter = sum(1 for c in shared if hsa[c] > lf[c])
    print(f"  L-FABP binds tighter for {tighter} of {len(shared)}.")
    print("  Note these are all ALTERNATIVE PFAS. Jia's legacy-compound")
    print("  HSA values are transcribed from Gao et al. 2019 and are a")
    print("  different laboratory and method, so they are excluded.")
    print()

    rule()
    print("3. CAPACITY, WHICH IS WHAT THE QUESTION ACTUALLY ASKED")
    rule()
    print("  ratio = (C_FABP1 / Kd_LFABP) / (C_albumin / Kd_HSA),")
    print("  equal stoichiometry assumed on both sides (see the header).")
    print()
    print(f"  {'compound':16s}{'per L cytosol':>15s}{'per L liver':>13s}"
          f"{'  over the albumin interval (per L cytosol)':>44s}")
    ratios = []
    for c in shared:
        aff = hsa[c] / lf[c]
        cyt = (C_FAB / lf[c]) / (C_ALB / hsa[c])
        liv = cyt * F_ICW
        lo = (C_FAB / lf[c]) / (C_ALB_HI / hsa[c])
        hi = (C_FAB / lf[c]) / (C_ALB_LO / hsa[c])
        ratios.append(cyt)
        print(f"  {c:16s}{cyt:>14.1f}x{liv:>12.1f}x"
              f"{'':>22s}{lo:>8.1f} to {hi:<8.1f}")
    print()
    gm = math.exp(statistics.fmean(math.log(r) for r in ratios))
    print(f"  geometric mean {gm:.1f}x per litre of cytosol, "
          f"{gm * F_ICW:.1f}x per litre of liver.")
    print(f"  range {min(ratios):.2f}x to {max(ratios):.1f}x -- "
          f"{sum(1 for r in ratios if r > 1)} of {len(ratios)} favour liver.")
    print()
    print("  SO THE ANSWER IS COMPOUND-SPECIFIC, AND THE SPREAD IS THE")
    print(f"  FINDING: {max(ratios) / min(ratios):.0f}-fold between the "
          f"most and least liver-favouring")
    print("  of five compounds measured in one run. 'Liver out-competes")
    print("  plasma' is true for four of them and false for the fifth.")
    print()

    rule()
    print("4. HOW MUCH OF THAT IS THE CONCENTRATION AND HOW MUCH THE Kd")
    rule()
    print("  The ratio factorises exactly:")
    print()
    print("      capacity ratio  =  (C_FABP1 / C_albumin)  x  (Kd_HSA / Kd_LFABP)")
    print("                      =   concentration term    x   affinity term")
    print()
    conc = C_FAB / C_ALB
    print(f"  The concentration term is the same {conc:.2f} for every")
    print("  compound. Everything that varies is affinity.")
    print()
    print(f"  {'compound':16s}{'concentration':>14s}{'affinity':>11s}"
          f"{'product':>10s}")
    for c in shared:
        aff = hsa[c] / lf[c]
        print(f"  {c:16s}{conc:>14.2f}{aff:>11.1f}{conc * aff:>10.1f}")
    print()
    print("  Which is worth stating plainly: the concentration that was")
    print("  missing turns out NOT to be what decides the comparison. It")
    print(f"  is a flat {conc:.2f}x handicap on the liver side. The answer is")
    print("  carried by affinity after all -- but only because the two")
    print("  concentrations happen to be within a factor of two, and that")
    print("  could not be known until the concentration was in hand.")
    print()

    rule()
    print("5. WHAT WOULD OVERTURN THIS")
    rule()
    print("  1. STOICHIOMETRY. The ratio carries n_HSA / n_LFABP and MST")
    print("     reports neither. Maso et al. 2021 co-crystallised four")
    print("     PFOA per albumin against L-FABP's usual one or two. At")
    print("     n_HSA = 4, n_LFABP = 1 every ratio above falls 4x:")
    for c in shared:
        cyt = (C_FAB / lf[c]) / (C_ALB / hsa[c])
        print(f"       {c:16s}{cyt:>8.1f}x  ->{cyt / 4:>8.2f}x")
    print(f"     {sum(1 for c in shared if (C_FAB / lf[c]) / (C_ALB / hsa[c]) / 4 > 1)}"
          f" of {len(shared)} would still favour liver. The sign of the")
    print("     answer survives for the strongest binders and not for the")
    print("     rest, so this is a real sensitivity, not a quibble.")
    print()
    print("  2. ALBUMIN IN THE LIVER. Liver is not only cytosol. It")
    print("     carries interstitial albumin that this comparison ignores")
    print("     entirely, which can only push the ratio DOWN. What is")
    print("     computed above is cytosolic L-FABP against plasma")
    print("     albumin, not liver against plasma.")
    print()
    print("  3. METHOD. Every number in part 2 is microscale")
    print("     thermophoresis. Khazaee's equilibrium dialysis puts PFOA")
    print("     on L-FABP at 0.099 uM where Jia's MST puts it at 17.95 --")
    print("     181x. If MST is biased the same way on both proteins the")
    print("     ratio survives, and if it is not, it does not. One run")
    print("     cannot tell, and no second paper measures both.")
    print()
    print("  4. FREE CONCENTRATION. All of this is the trace limit. The")
    print("     two proteins saturate at different places, so at a real")
    print("     hepatic exposure the ratio moves toward whichever is")
    print("     further from saturation.")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
