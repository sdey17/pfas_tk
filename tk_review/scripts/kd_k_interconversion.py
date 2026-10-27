#!/usr/bin/env python3
"""Converting between a PFAS-albumin Kd and an albumin-water partition
coefficient: the formula, the procedure, and how many rows can take it.

kd_to_partition.py derives the bridge and validates it inside Qin 2023.
This file is the operating manual: the two directions, the three inputs
each needs, a worked seven-step conversion on a real row, and an inventory
of how much of the database can actually be converted.

THE RULE THIS FILE ENFORCES: convert only within one paper, which must
supply the affinity AND the stoichiometry. Part 5 shows what happens when
you take a Kd from one study and a partition coefficient from another --
the reconciling n spans 752-fold and is ordered by chain length -- so that
is a diagnosis of why the rule exists, not a method to use.

Part 6 keeps the alternative that avoids the conversion entirely:
regressing log K on chain length within a single study.

THE COUNT, for the impatient: 82 of 532 albumin rows can be converted
forward, all of them because one paper reported n. Nought can be converted
in reverse, and that zero is structural rather than a gap in extraction --
a trace-regime partition coefficient constrains only the ratio n/Kd and
cannot separate them in principle.

Run:  python3 tk_review/scripts/kd_k_interconversion.py
"""
import collections
import csv
import math
import re
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
    return ncf, chen, fis, rows


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


# ---------------------------------------------------------------------
# Which rows carry a usable stoichiometry, and how it is stated.
# ---------------------------------------------------------------------
NUMERIC_N = re.compile(r"^\s*\d+(?:\.\d+)?\s*$")


def stated_n(row):
    """(n, how it was stated) for one database row, or (None, None).

    Three ways a paper can give you n, in descending order of strength:
      - a fitted number in the row;
      - a model statement that fixes it ("one site assumed"), which is the
        paper's own assumption and may be carried as such;
      - a per-compound figure in prose. Gao 2025 states Hill coefficients
        and reads sites off them, which is weaker: a Hill coefficient is a
        cooperativity, and equating it with a site count is the paper's
        interpretation, not a measurement. Kept separate for that reason.
    """
    v = row["n_binding_sites"].strip()
    if NUMERIC_N.match(v):
        return float(v), "fitted"
    if "one site assumed" in v:
        return 1.0, "model"
    if "site for" in v:
        return None, "prose"
    return None, None


def is_affinity(row):
    return (bool(row["kd"].strip() or row["ka"].strip())
            and not row["ka_units"].startswith("log"))


def is_partition(row):
    return row["ka_units"].startswith("log") and bool(row["ka"].strip())


def part3():
    rule()
    print("3. THE PROCEDURE, STEP BY STEP, ON ONE ROW")
    rule()
    print("  Worked on Chen et al. 2025, PFOA on human serum albumin,")
    print("  Kd = 2.57 uM, equilibrium dialysis, 0.0625-2 uM ligand against")
    print("  1.0 uM HSA. Every step is a check you can fail; failing one")
    print("  means stop, not adjust.")
    print()
    kd, n, cw_max, prot = 2.57, 1.0, 2.0, 1.0

    print("  STEP 1  Fix the material and the species.")
    print("          Human serum albumin. The output is a coefficient for")
    print("          HUMAN albumin and for nothing else. Do not carry it to")
    print("          bovine, to plasma, or to 'protein'.")
    print()
    print("  STEP 2  Put Kd in uM.")
    print(f"          Kd = {kd} uM. If the paper gives Ka in M^-1, "
          f"Kd(uM) = 1e6/Ka.")
    print("          If it gives Kd in mM, multiply by 1000. This is where")
    print("          most of the published disagreement actually comes from.")
    print()
    print("  STEP 3  Get n FROM THE SAME PAPER. Do not supply it yourself.")
    print(f"          Chen fits B_max per compound with one site assumed, so")
    print(f"          n = {n:.0f} is the paper's own model. log n = "
          f"{math.log10(n):+.2f}.")
    print("          No n in the paper -> the conversion cannot be done. It")
    print("          is not optional and it is not 1 by default.")
    print()
    print("  STEP 4  Check the regime. The formula is a trace-limit result.")
    ratio = cw_max / kd
    print(f"          Highest Cw dosed {cw_max} uM, Kd {kd} uM, so "
          f"Cw/Kd = {ratio:.2f}.")
    print(f"          Correction factor 1 + Cw/Kd = {1 + ratio:.2f}.")
    print(f"          At the top of Chen's series the apparent coefficient")
    print(f"          is {1 / (1 + ratio):.2f} of its trace value. Below "
          f"Cw/Kd ~ 0.1 you")
    print("          may ignore this; here you may not.")
    print()
    print("  STEP 5  Apply the conversion.")
    lk = C_MASS + math.log10(n) - math.log10(kd)
    print(f"          log K = {C_MASS:.4f} + {math.log10(n):.2f} - "
          f"log10({kd}) = {lk:.2f}  [L_w/kg protein]")
    print(f"          log D = {C_VOL:.4f} + {math.log10(n):.2f} - "
          f"log10({kd}) = {lk + 0.1349:.2f}  [L_w/L protein]")
    print()
    print("  STEP 6  State the basis on the output, always.")
    print("          Per KILOGRAM and per LITRE of protein differ by 0.135")
    print("          log. Fischer and Allendorf report per kg; Qin reports")
    print("          per L. A number without its basis is unusable.")
    print()
    print("  STEP 7  Sanity-check against a direct measurement of the SAME")
    print("          compound in the SAME species -- and if they disagree,")
    print("          report the disagreement. Do not tune n until they")
    print("          agree; n is a measurement, not a free parameter.")
    print(f"          converted   {lk:.2f}")
    print(f"          Fischer 2024 measured, human HSA, per kg   4.48")
    print(f"          -> {10 ** (4.48 - lk):.1f}x apart. That gap is real "
          f"and is diagnosed")
    print("             in part 5; it is not an error in the arithmetic.")
    print()
    print("  THE REVERSE DIRECTION, K -> Kd, is the same formula rearranged")
    print("  and needs the same n from the same paper. Part 4 shows that no")
    print("  row in this database can take it.")
    print()


def part4(rows):
    rule()
    print("4. HOW MANY ROWS CAN ACTUALLY BE CONVERTED")
    rule()
    alb = [r for r in rows if "albumin" in r["protein"].lower()]
    fwd = collections.defaultdict(list)
    rev = collections.defaultdict(list)
    no_n_aff = no_n_part = 0
    prose = 0
    for r in alb:
        n, how = stated_n(r)
        src = f"{r['study']} {r['year']}"
        if is_affinity(r):
            if n is not None:
                fwd[(src, how)].append(r)
            elif how == "prose":
                prose += 1
            else:
                no_n_aff += 1
        elif is_partition(r):
            if n is not None:
                rev[(src, how)].append(r)
            else:
                no_n_part += 1
    nfwd = sum(len(v) for v in fwd.values())
    nrev = sum(len(v) for v in rev.values())
    print(f"  {len(alb)} albumin rows in db/protein_binding.csv.")
    print()
    print(f"  FORWARD, Kd or Ka -> partition coefficient: {nfwd} rows")
    for (src, how), v in sorted(fwd.items(), key=lambda kv: -len(kv[1])):
        print(f"    {len(v):4d}  {src:24s} n {how}")
    print(f"    ---- not counted in the {nfwd}, decide for yourself ----")
    print(f"    {prose:4d}  Gao et al. 2025          n from prose Hill "
          f"coefficients -- weaker,")
    print(f"          counted separately: a Hill coefficient is a "
          f"cooperativity,")
    print(f"          and reading a site count off it is the paper's "
          f"interpretation")
    print(f"    {no_n_aff:4d}  affinity rows with NO n anywhere -- cannot be "
          f"converted")
    print()
    print(f"  REVERSE, partition coefficient -> Kd: {nrev} rows")
    print(f"    {no_n_part:4d}  partition rows, and not one of them carries "
          f"an n.")
    print()
    print("  The reverse zero is STRUCTURAL, not a gap in our extraction.")
    print("  In the trace limit the measurement is")
    print()
    print("      D0 = n x M_alb / Kd")
    print()
    print("  so a partition coefficient determines only the RATIO n/Kd. It")
    print("  cannot separate them, in principle, however well it is made. A")
    print("  paper that works in the trace regime -- which is what makes its")
    print("  coefficient a constant at all -- has thrown away the")
    print("  information needed to recover a Kd. You need an experiment that")
    print("  reaches saturation, and then you are measuring a Kd directly")
    print("  and do not need the conversion.")
    print()
    print("  Restricted to the ten sources in db/binding_major_sources.xlsx:")
    keep = {"Chen et al. 2025", "Wu et al. 2024", "Gao et al. 2025"}
    infile = sum(len(v) for (src, how), v in fwd.items() if src in keep)
    print(f"    {infile} rows convertible forward, {prose} more on the "
          f"weaker prose route,")
    print(f"    0 convertible in reverse.")
    print()


def part5(ncf, chen, fis):
    rule()
    print("5. WHY THE CONVERSION MUST STAY INSIDE ONE PAPER")
    rule()
    print("  Step 7 left a gap: Chen's Kd converts to 3.77 where Fischer")
    print("  measured 4.48 on the same protein and species, 5.2x apart. It")
    print("  is worth knowing why, because the reason says when the")
    print("  conversion is safe at all.")
    print()
    print("  Solve for the n that would reconcile the two, compound by")
    print("  compound. If the conversion were transferable it would come")
    print("  out constant.")
    print()
    print(f"  {'cmpd':10s}{'nCF':>4s}{'Kd uM':>9s}{'logK meas':>10s}"
          f"{'logK n=1':>9s}{'n implied':>10s}")
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
    sl, r = fit([t[0] for t in ns], [math.log10(t[1]) for t in ns])
    print()
    print(f"  It is not constant: {lo[1]:.2f} ({lo[2]}) to {hi[1]:.1f} "
          f"({hi[2]}), a {hi[1] / lo[1]:.0f}-fold span,")
    print(f"  and it rises {sl:+.3f} log per CF2 (r {r:+.2f}). An ordered")
    print("  residual is a missing mechanism, not scatter. The mechanism:")
    print()
    for lab, ser in (("carboxylates", PFCA), ("sulfonates", PFSA)):
        xk = [(ncf[c], -math.log10(chen[c])) for c in ser if c in chen]
        xd = [(ncf[c], fis[c]) for c in ser if c in fis]
        sk, rk = fit([a for a, _ in xk], [b for _, b in xk])
        sd, rd = fit([a for a, _ in xd], [b for _, b in xd])
        print(f"    {lab:14s} -log Kd {sk:+.3f}/CF2 (r {rk:+.2f})   "
              f"log K {sd:+.3f}/CF2 (r {rd:+.2f})   gap {sd - sk:+.3f}")
    print()
    print("  The formula says the two must have the SAME chain-length slope")
    print("  if n does not depend on chain length. They do not. Affinity is")
    print("  nearly flat and not even monotone -- Chen's Kd reaches an")
    print("  optimum near PFOA-PFDA and weakens again (PFDoDA 15.8 uM")
    print("  against PFOA 2.6 uM) -- while partitioning climbs steadily.")
    print()
    print("  Those are compatible only if the partition coefficient includes")
    print("  something the Kd does not: a non-saturating term that grows")
    print("  with hydrophobicity. So the implied n is not a site count, it")
    print("  is (total capacity)/(one-site capacity), and it is near 1 only")
    print("  for the short carboxylates where the specific site dominates.")
    print()
    print("  CONCLUSION. Convert only where one paper gives you the affinity")
    print("  and the stoichiometry together, and label the output what it")
    print("  is: the SPECIFIC term, not the albumin-water partition")
    print("  coefficient. Above about C8 the specific term is the smaller")
    print("  of the two contributions.")
    print()


def part6(ncf, fis):
    rule()
    print("6. THE OTHER OPTION: REGRESS log K ON CHAIN LENGTH DIRECTLY")
    rule()
    print("  This never touches a Kd, so none of the above applies. Fit")
    print("  within ONE study, ONE head-group class and ONE basis, and")
    print("  interpolate inside the measured range.")
    print()
    print("  Fischer et al. 2024, human serum albumin, L_w per kg protein:")
    for lab, ser in (("carboxylates", PFCA), ("sulfonates", PFSA)):
        pts = [(ncf[c], fis[c], c) for c in ser if c in fis]
        sl, r = fit([p[0] for p in pts], [p[1] for p in pts])
        xs = [p[0] for p in pts]
        inter = (sum(p[1] for p in pts) / len(pts)
                 - sl * sum(xs) / len(xs))
        print(f"    {lab:14s} log K = {sl:+.3f} x nCF {inter:+.2f}   "
              f"r {r:+.2f}, n = {len(pts)}, nCF {min(xs)}-{max(xs)}")
    print()
    print("  Two cautions that the r values hide.")
    print("    - The series TURNS OVER. Fischer's own carboxylates run 5.18")
    print("      at PFUnDA down to 3.87 at PFTeDA, so a line fitted across")
    print("      the whole range understates the middle and overstates both")
    print("      ends. Interpolate; do not extrapolate past the last")
    print("      measured chain length.")
    print("    - The slope is a property of the study's basis. Allendorf")
    print("      2019 (bovine) and Qin 2024 (bovine) measure the same")
    print("      quantity and would give their own intercepts; do not mix")
    print("      an intercept from one with a slope from another.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    ncf, chen, fis, rows = load()
    part1()
    part2()
    part3()
    part4(rows)
    part5(ncf, chen, fis)
    part6(ncf, fis)
    return 0


if __name__ == "__main__":
    sys.exit(main())
