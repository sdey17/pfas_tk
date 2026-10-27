#!/usr/bin/env python3
"""Converting between a PFAS-albumin Kd and an albumin-water partition
coefficient: the formula, what it costs you, and when you may extrapolate.

kd_to_partition.py derives the bridge and validates it inside Qin 2023.
This file answers the practical question -- can I take a Kd from one paper
and a partition coefficient from another and move between them? -- and the
short answer is: the algebra is exact, the application is not, and the
reason is diagnosable rather than mysterious.

Run:  python3 tk_review/scripts/kd_k_interconversion.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

M_ALB = 66463.0           # g/mol, serum albumin
RHO_ALB = 1.0 / 0.733     # kg/L, from the 0.733 mL/g partial specific volume
# The albumin molarity inside the pure protein phase, in micromoles per kg
# (or per L) of protein. 1 kg of albumin is 1e9/M_alb umol; one litre is
# 1.364x that. Kd is carried in uM throughout, so the 1e9 and not 1e6 is
# what makes the two sides of the conversion commensurate -- getting this
# wrong is a silent factor of 1000.
C_MASS = math.log10(1e9 / M_ALB)              # 4.1774, per kg of protein
C_VOL = math.log10(1e9 * RHO_ALB / M_ALB)     # 4.3123, per L of protein

PFCA = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
        "PFUnDA", "PFDoDA", "PFTrDA", "PFTeDA"]
PFSA = ["PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS"]


def rule(ch="=", n=78):
    print(ch * n)


def load():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    ncf = {}
    for r in rows:
        v = r["n_perfluorinated_carbons"].strip()
        if v.isdigit():
            ncf.setdefault(r["chemical"], int(v))
    chen = {r["chemical"]: float(r["kd"]) for r in rows
            if r["doi"] == "10.1021/acs.est.4c06734"
            and "albumin" in r["protein"].lower() and r["kd"].strip()}
    fis = {r["chemical"]: float(r["ka"]) for r in rows
           if r["study"] == "Fischer et al." and r["year"] == "2024"
           and r["protein"] == "human serum albumin (HSA)" and r["ka"].strip()}
    return ncf, chen, fis


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((a - mx) ** 2 for a in xs)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / sxx, sxy / math.sqrt(sxx * syy)


# =====================================================================
def part1():
    rule()
    print("1. THE TWO DIRECTIONS")
    rule()
    print("  One class of n equivalent sites, dissociation constant Kd, and")
    print("  Langmuir occupancy. The albumin concentration CANCELS, which is")
    print("  why a partition coefficient can be written from an affinity at")
    print("  all. In the trace limit Cw << Kd (derivation in")
    print("  kd_to_partition.py):")
    print()
    print("      FORWARD   Kd  ->  K")
    print(f"        log K [L_w/kg protein] = {C_MASS:.4f} + log n "
          f"- log Kd(uM)")
    print(f"        log D [L_w/L protein]  = {C_VOL:.4f} + log n "
          f"- log Kd(uM)")
    print()
    print("      REVERSE   K  ->  Kd")
    print(f"        log Kd(uM) = {C_MASS:.4f} + log n - log K [L_w/kg]")
    print(f"        log Kd(uM) = {C_VOL:.4f} + log n - log D [L_w/L]")
    print()
    print("  The two constants differ by log10(1.364) = 0.135, the density of")
    print("  the pure albumin phase. Mixing the bases is a 1.4x error, which")
    print("  is small enough to hide and large enough to matter.")
    print()
    print("  Worked both ways, PFOA, Chen et al. 2025 Kd = 2.57 uM, n = 1:")
    kd = 2.57
    lk = C_MASS - math.log10(kd)
    print(f"    forward:  log K = {C_MASS:.4f} - log10({kd}) = {lk:.2f} "
          f"[L_w/kg]")
    print(f"    reverse:  log Kd = {C_MASS:.4f} - {lk:.2f} = "
          f"{C_MASS - lk:.2f}  ->  Kd = {10 ** (C_MASS - lk):.2f} uM")
    print()
    print("  Note what the forward conversion is NOT: it is not a fit, and")
    print("  there is no empirical coefficient in it. If both numbers are")
    print("  right and describe the same sites, it holds exactly.")
    print()


def part2():
    rule()
    print("2. THE THREE THINGS YOU MUST SUPPLY, AND WHAT EACH COSTS")
    rule()
    print("  (a) n, the stoichiometry. It enters as log n, so it is a pure")
    print("      multiplier on the answer and almost no paper reports it.")
    print(f"        n = 1  ->  +0.00 log        n = 6  -> "
          f"{math.log10(6):+.2f} log ({6:.0f}x)")
    print(f"        n = 2  -> {math.log10(2):+.2f} log        n = 20 -> "
          f"{math.log10(20):+.2f} log ({20:.0f}x)")
    print("      Qin 2023 is the rare study that fits it, and gets 0.39 to")
    print("      1.24 sites per albumin -- SUB-stoichiometric, so n is an")
    print("      effective capacity, not a count of pockets.")
    print()
    print("  (b) The concentration regime of BOTH measurements. The trace")
    print("      limit is an assumption, and the exact correction is")
    print()
    print("          D(Cw) = D0 / (1 + Cw/Kd)")
    print()
    print("      so a coefficient measured at Cw reads low by that factor:")
    print(f"        {'Cw/Kd':>8s}{'factor':>10s}   {'reads':>12s}")
    for ratio in (0.1, 1, 3, 10, 100):
        print(f"        {ratio:>8g}{1 + ratio:>10.1f}x   "
              f"{1 / (1 + ratio):>11.3f}x of trace")
    print("      A bioassay at 10 uM against a Kd of 2.6 uM is already 4.9x")
    print("      below the trace value. This is the single commonest reason")
    print("      two literature constants for one compound disagree.")
    print()
    print("  (c) The non-specific term. Real albumin sorbs PFAS by a route")
    print("      that does NOT saturate, on top of the saturable site:")
    print()
    print("          D_total(Cw) = n*M_alb/(Kd + Cw)  +  K_nonspecific")
    print()
    print("      A Kd describes only the first term. Part 4 shows that the")
    print("      second is what breaks the cross-study conversion.")
    print()


def part3(ncf, chen, fis):
    rule()
    print("3. THE TEST: does it work between two independent studies?")
    rule()
    print("  Both measured HUMAN serum albumin, so no species is mixed.")
    print("  Chen et al. 2025 gives Kd by equilibrium dialysis; Fischer et")
    print("  al. 2024 gives log D by SPME at a PFAS:protein ratio below")
    print("  0.004 mol/mol, which IS the trace regime the formula wants.")
    print("  So assumption (b) is satisfied by construction. Solve for the n")
    print("  that would reconcile them:")
    print()
    print(f"  {'cmpd':10s}{'nCF':>4s}{'Kd uM':>9s}{'logD meas':>10s}"
          f"{'logD n=1':>9s}{'n implied':>10s}")
    ns = []
    for c in PFCA + PFSA:
        if c not in chen or c not in fis:
            continue
        pred = C_MASS - math.log10(chen[c])
        n = 10 ** (fis[c] - pred)
        ns.append((ncf[c], n, c))
        print(f"  {c:10s}{ncf[c]:>4d}{chen[c]:>9.3g}{fis[c]:>10.2f}"
              f"{pred:>9.2f}{n:>10.2f}")
    lo = min(ns, key=lambda t: t[1])
    hi = max(ns, key=lambda t: t[1])
    print()
    print(f"  n runs from {lo[1]:.2f} ({lo[2]}) to {hi[1]:.1f} ({hi[2]}) -- a "
          f"{hi[1] / lo[1]:.0f}-fold span.")
    print("  No protein has between 0.2 and 120 sites. The conversion has")
    print("  failed, and n is absorbing something that is not stoichiometry.")
    sl, r = fit([t[0] for t in ns], [math.log10(t[1]) for t in ns])
    print(f"  And the failure is ORDERED: implied n rises {sl:+.3f} log per")
    print(f"  CF2 (r {r:+.2f}). Ordered residuals are a missing mechanism,")
    print("  not measurement scatter.")
    print()


def part4(ncf, chen, fis):
    rule()
    print("4. THE DIAGNOSIS: the two quantities have different chain-length")
    print("   slopes, so no single n can reconcile them")
    rule()
    print("  The formula says log K and -log Kd differ by the constant")
    print("  log n. Differentiate along the chain-length series: if n does")
    print("  not depend on chain length, the two slopes must be EQUAL.")
    print()
    for lab, ser in (("carboxylates", PFCA), ("sulfonates", PFSA)):
        xk = [(ncf[c], -math.log10(chen[c])) for c in ser if c in chen]
        xd = [(ncf[c], fis[c]) for c in ser if c in fis]
        sk, rk = fit([a for a, _ in xk], [b for _, b in xk])
        sd, rd = fit([a for a, _ in xd], [b for _, b in xd])
        print(f"  {lab}")
        print(f"    -log Kd vs nCF  {sk:+.3f} per CF2  (r {rk:+.2f}, "
              f"n={len(xk)})   affinity")
        print(f"     log D  vs nCF  {sd:+.3f} per CF2  (r {rd:+.2f}, "
              f"n={len(xd)})   capacity")
        print(f"    gap             {sd - sk:+.3f} per CF2  -> "
              f"{10 ** (sd - sk):.2f}x per CF2 of unexplained capacity")
    print()
    print("  Read the affinity row first. Chen's Kd is nearly FLAT in chain")
    print("  length for the carboxylates, and it is not monotone: it falls")
    print("  to an optimum near PFOA-PFDA and rises again for the longest")
    print("  chains (PFDoDA 15.8 uM, PFTrDA 31.4 uM, against PFOA 2.6 uM).")
    print("  Affinity peaks and then weakens. Partitioning does not -- it")
    print("  climbs steadily across the whole series.")
    print()
    print("  Those two facts are only compatible if what the partition")
    print("  coefficient measures is not only the saturable site. It is")
    print("  the sum in part 2(c), and the non-specific term grows with")
    print("  hydrophobicity while the specific site does not. That is the")
    print("  gap above, and it is exactly the split Qin 2023 and Qin 2024")
    print("  fit explicitly as D_specific and D_nonspecific.")
    print()
    print("  So the implied n of part 3 is not a count of sites. It is")
    print()
    print("      n_implied = (total capacity) / (one-site capacity)")
    print()
    print("  which is near 1 only where the specific site dominates -- the")
    print("  short carboxylates -- and grows without bound as the")
    print("  non-specific term takes over.")
    print()


def part5():
    rule()
    print("5. SO WHAT CAN YOU ACTUALLY DO")
    rule()
    print("  WORKS -- converting within one study that reports both.")
    print("    Qin 2023 reports Kd AND Cmax, so n is known and the")
    print("    conversion is exact; kd_to_partition.py part 2 recovers the")
    print("    albumin molarity to within 1% for all four compounds. If a")
    print("    paper gives you an affinity and a capacity, use them.")
    print()
    print("  WORKS -- converting a Kd into the SPECIFIC term only.")
    print("    log D_specific = C + log n - log Kd is sound. Just do not")
    print("    call the result the albumin-water partition coefficient: it")
    print("    is one of two terms, and for anything above about C8 it is")
    print("    the smaller one.")
    print()
    print("  WORKS -- correcting a coefficient for its measurement regime.")
    print("    Multiply by (1 + Cw/Kd) to reach the trace limit. You need")
    print("    the Cw the constant was measured at, which is why the")
    print("    concentration columns were added to db/binding_major_"
          "sources.xlsx.")
    print()
    print("  DOES NOT WORK -- a Kd from study A into a total partition")
    print("    coefficient in study B. Part 3 is that experiment, done")
    print("    under the best conditions available in this database -- one")
    print("    species, both in the trace regime -- and it fails by up to")
    print("    120x with a chain-length-ordered residual.")
    print()
    print("  DOES NOT WORK -- extrapolating a Kd along chain length to get")
    print("    a partition coefficient. The two have different slopes")
    print("    (part 4), so the error grows about 1.4x per CF2. Over PFBA")
    print("    to PFDoDA that is two orders of magnitude.")
    print()
    print("  WHAT TO DO INSTEAD, if you need a partition coefficient for a")
    print("  compound nobody has measured one for:")
    print("    1. Regress log D on chain length DIRECTLY, within a head-")
    print("       group class and within one study's basis. The slopes are")
    print("       +0.21 (carboxylates) and +0.26 (sulfonates) per CF2 in")
    print("       Fischer's human albumin series, and the same measurement")
    print("       exists in Allendorf 2019 (bovine) and Qin 2024 (bovine).")
    print("    2. Interpolate, do not extrapolate past the measured range.")
    print("       Affinity turns over near C10 and capacity may too; the")
    print("       Fischer series itself drops from 5.18 at PFUnDA to 3.87")
    print("       at PFTeDA.")
    print("    3. Use the Kd only for what it is good for -- telling you")
    print("       WHERE the saturable site runs out, which is what decides")
    print("       whether any of this is in the linear regime at all.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    ncf, chen, fis = load()
    part1()
    part2()
    part3(ncf, chen, fis)
    part4(ncf, chen, fis)
    part5()
    return 0


if __name__ == "__main__":
    sys.exit(main())
