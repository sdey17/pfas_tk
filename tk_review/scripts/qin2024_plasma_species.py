#!/usr/bin/env python3
"""The human/trout plasma species difference, from measured isotherms.

Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L.  Environ Sci Technol
2024;58(23):9954-9966.  doi:10.1021/acs.est.3c10824.  PMC11171458.

Table 1 gives Freundlich-type fits of log D_plasma/w against log Cw for 13
chemicals in BOTH human and rainbow trout plasma, and Table 3 gives each
plasma's measured protein and lipid content.  That combination lets the
species difference be examined as a function of concentration rather than at
one point, which matters because the two plasmas do not saturate alike.

THE PAPER'S CLAIM, in its own words: "Differences in the binding behavior of
PFAS between human and trout plasma were explained by differences in lipid
and protein contents."  Composition, not protein identity.  Its own test of
that is a prediction "within a factor of 10" at single concentrations
(Figure 4).  With the fits in hand the same claim can be tested about three
decades wide and to much better than a factor of 10.

THE FINDING.  Composition is necessary and it is not sufficient, and the
test that shows it needs no partition coefficient at all.  With one shared
protein coefficient and one shared lipid coefficient, the human/trout ratio
is algebraically confined to 0.30x-1.32x by the two measured compositions
alone -- it is monotone in D_protein/D_lipid, so it lies between the two
weight ratios.  Measured, the ratio runs from 23x (PFBA) to 0.26x (PFNA) at
1 umol/L, and 5 of 11 PFAS fall outside a bound no choice of coefficients
can widen; extrapolated to 1 nmol/L it is 7 of 11.  The residual is ordered
by chain length, not scattered: human plasma binds the short compounds
harder and trout the long ones, with the handover near nCF 8 and a slope of
-0.30 to -0.46 log units per CF2 (r -0.94 to -0.97).  The shape disagrees
too -- trout plasma shows no detectable saturation for 9 of 13 chemicals
where human plasma saturates for 10 of 13, and a weighted mean of a curved
term and a constant is still curved.

The alternatives this dataset cannot separate are set out in part 4; the
temperature difference between the two arms (25 C for trout, 37 C for human)
is a real confound in the bound, not a quibble.

WHY THIS PROJECT CARES.  A species difference that depends on dose cannot be
carried across species as a single number.  If the ratio between human and
trout plasma binding is 3x at assay concentrations and 30x at physiological
ones, then a free fraction measured in one species at bioassay doses is not
transferable to the other at human exposure levels -- which is exactly the
operation a cross-species toxicokinetic extrapolation performs.

BASIS.  Both plasma columns are [L_water / L_protein+lipid], so they are
directly comparable to each other.  The BSA column is [L_water / L_protein].
Table 3: human plasma 42.25 mL/L protein and 4.46 mL/L lipid; trout 15.46
and 7.08.  So of the SORPTIVE volume, protein is 90.5% in human plasma and
68.6% in trout -- the whole of the paper's mechanism is in those two numbers.

Run:  python3 tk_review/scripts/qin2024_plasma_species.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from add_qin2024_companion_table1 import T1     # Table 1, all three columns

# Table 3, measured. mL per L of plasma.
HP_PROT, HP_LIP = 42.25, 4.46
TP_PROT, TP_LIP = 15.46, 7.08
HP_FPROT = HP_PROT / (HP_PROT + HP_LIP)
TP_FPROT = TP_PROT / (TP_PROT + TP_LIP)
# Whole-plasma volume fractions, for the free-fraction calculation.
HP_VPL = (HP_PROT + HP_LIP) / 1000.0
TP_VPL = (TP_PROT + TP_LIP) / 1000.0

# Table 2, log D_lip/w (liposome-water) -- the lipid coefficient of the MBM.
# Measured by Droge 2019 / Ebert 2020 except 6:2 FTSA and PFOSA, predicted by
# Qin from chain length, and the two pesticides, predicted with COSMOtherm.
D_LIP = {"PFBA": 1.00, "PFHxA": 2.32, "PFHpA": 2.91, "PFOA": 3.52,
         "PFNA": 4.25, "PFUnA": 4.54, "HFPO-DA": 2.41, "PFHxS": 4.13,
         "PFOS": 4.89, "6:2 FTSA": 3.87, "PFOSA": 4.94,
         "hexaflumuron": 4.58, "flubendiamide": 3.28}
PREDICTED_LIP = {"6:2 FTSA", "PFOSA", "hexaflumuron", "flubendiamide"}
PFAS = [n for n in T1 if n not in ("hexaflumuron", "flubendiamide")]

# Concentrations to evaluate at. 1 uM sits inside the measured series for
# every compound; the others are stated as extrapolation where they are.
CW = [("1 nM", 1e-3), ("100 nM", 0.1), ("1 uM", 1.0), ("100 uM", 100.0)]


def logD(fit, cw):
    slope, inter, _ = fit
    return inter if slope is None else slope * math.log10(cw) + inter


def bsa(nm):
    return T1[nm][2]


def hp(nm):
    return T1[nm][3]


def tp(nm):
    return T1[nm][4]


def rule(ch="=", n=78):
    print(ch * n)


def gm(vals):
    return 10 ** (sum(vals) / len(vals))


# =====================================================================
def part1():
    rule()
    print("1. SATURATION IS A PROTEIN PHENOMENON -- count the flat fits")
    rule()
    print("  Qin's rule: a fitted Freundlich exponent at or above 0.9 is")
    print("  forced to 1 and the result reported as a concentration-")
    print("  independent constant. So a flat entry means the isotherm was")
    print("  not detectably curved over the range tested.")
    print()
    print(f"  {'chemical':15s} {'BSA':>8s} {'human':>8s} {'trout':>8s}"
          f"   human-trout")
    nflat = {"BSA": 0, "human": 0, "trout": 0}
    diffs = []
    for nm in T1:
        f = lambda x: "flat" if x[0] is None else f"{x[0]:+.3f}"
        for lab, fit in (("BSA", bsa(nm)), ("human", hp(nm)),
                         ("trout", tp(nm))):
            if fit[0] is None:
                nflat[lab] += 1
        sh = hp(nm)[0] or 0.0
        st = tp(nm)[0] or 0.0
        diffs.append(sh - st)
        print(f"  {nm:15s} {f(bsa(nm)):>8s} {f(hp(nm)):>8s} "
              f"{f(tp(nm)):>8s}   {sh - st:+.3f}")
    n = len(T1)
    print()
    for lab in ("BSA", "human", "trout"):
        print(f"  {lab:>6s}: {n - nflat[lab]:2d} of {n} concentration-"
              f"dependent, {nflat[lab]:2d} flat")
    print()
    print("  That is the species result in its simplest form. Human plasma")
    print("  saturates for 10 of 13 chemicals, trout plasma for 4. Pure")
    print("  albumin saturates for 10 of 13, so trout is the outlier, not")
    print("  human. Saturation requires a saturable site, the site is on")
    print("  protein, and trout plasma is the one matrix here where protein")
    print(f"  carries only {TP_FPROT * 100:.0f}% of the sorptive volume "
          f"against {HP_FPROT * 100:.0f}% in human.")
    print()
    mean_h = sum(abs(hp(nm)[0] or 0.0) for nm in T1) / n
    mean_t = sum(abs(tp(nm)[0] or 0.0) for nm in T1) / n
    print(f"  Mean |exponent| over all 13: human {mean_h:.3f}, trout "
          f"{mean_t:.3f} -- a {mean_h / mean_t:.1f}-fold")
    print("  difference in how strongly the apparent affinity falls with")
    print("  dose. The dilution a bioassay imposes on measured binding is")
    print("  therefore species-specific, before any chemistry is invoked.")
    print()


def part2():
    rule()
    print("2. THE SPECIES RATIO DEPENDS ON DOSE")
    rule()
    print("  D_human/D_trout, from the two fits, at four concentrations.")
    print("  A single cross-species factor would be constant across the row.")
    print()
    hdr = "".join(f"{lab:>10s}" for lab, _ in CW)
    print(f"  {'chemical':15s}{hdr}   {'span':>7s}  crossover")
    spans = []
    for nm in T1:
        vals = [logD(hp(nm), c) - logD(tp(nm), c) for _, c in CW]
        span = 10 ** (max(vals) - min(vals))
        spans.append(span)
        # Cw where the two curves cross, if the slopes differ.
        dsl = (hp(nm)[0] or 0.0) - (tp(nm)[0] or 0.0)
        dint = hp(nm)[1] - tp(nm)[1]
        if abs(dsl) > 1e-9:
            cross = 10 ** (-dint / dsl)
            # Only report a crossover inside a range any experiment could
            # occupy. Where the two exponents are nearly equal the algebra
            # returns a number like 1e-17 uM, which is not a crossover, it
            # is two near-parallel lines.
            cs = (f"{cross:.3g} uM" if 1e-5 <= cross <= 1e5
                  else "none in range")
        else:
            cs = "parallel"
        cells = "".join(f"{10 ** v:>9.2f}x" for v in vals)
        print(f"  {nm:15s}{cells}   {span:>6.1f}x  {cs}")
    print()
    print(f"  The ratio moves by a median of "
          f"{sorted(spans)[len(spans) // 2]:.1f}x across those four")
    print("  concentrations. It is not a species constant.")
    print()
    print("  The sign is not a property of the species pair. It is ordered")
    print("  by chain length, and it reverses. Over the carboxylate series,")
    print("  where chain length is the only thing varying:")
    print()
    print(f"  {'chemical':10s}{'nCF':>4s}{'ratio @1 nM':>13s}"
          f"{'ratio @1 uM':>13s}")
    carb = [n for n in PFAS if T1[n][1] == "carboxylate"]
    pts = []
    for nm in carb:
        r_lo = logD(hp(nm), 1e-3) - logD(tp(nm), 1e-3)
        r_hi = logD(hp(nm), 1.0) - logD(tp(nm), 1.0)
        pts.append((T1[nm][0], r_lo, r_hi))
        print(f"  {nm:10s}{T1[nm][0]:>4d}{10 ** r_lo:>12.2f}x"
              f"{10 ** r_hi:>12.2f}x")
    for lab, idx in (("1 nM", 1), ("1 uM", 2)):
        xs = [p[0] for p in pts]
        ys = [p[idx] for p in pts]
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((a - mx) ** 2 for a in xs)
        sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
        sl = sxy / sxx
        syy = sum((b - my) ** 2 for b in ys)
        r = sxy / math.sqrt(sxx * syy)
        cross = mx - my / sl
        print(f"  at {lab:5s}: {sl:+.3f} log units per CF2, r {r:+.3f}, "
              f"sign reverses at nCF {cross:.1f}")
    print()
    print("  Human plasma binds the short compounds harder, trout plasma the")
    print("  long ones, and the handover is near nCF 8. That is the paper's")
    print("  mechanism made quantitative: the short compounds depend on the")
    print("  saturable protein site, which human plasma has more of, and the")
    print("  long ones partition into lipid, which trout plasma has more of.")
    print("  A single species factor would be a horizontal line here; the")
    print("  measured one falls through two orders of magnitude and crosses")
    print("  unity inside the series.")
    print()


def part3():
    rule()
    print("3. IS IT COMPOSITION, OR PROTEIN IDENTITY?")
    rule()
    print("  The test. Take ONE protein coefficient -- bovine albumin, the")
    print("  same for both species, concentration-dependent as measured --")
    print("  and ONE lipid coefficient, and weight them by each plasma's own")
    print("  measured composition:")
    print()
    print("      D_plasma/w = D_BSA/w(Cw) x f_prot + D_lip/w x f_lip")
    print()
    print(f"      human  f_prot {HP_FPROT:.3f}  f_lip {1 - HP_FPROT:.3f}")
    print(f"      trout  f_prot {TP_FPROT:.3f}  f_lip {1 - TP_FPROT:.3f}")
    print()
    print("  If the species difference is composition, this reproduces both")
    print("  measured curves with no species-specific protein term anywhere.")
    print("  Residuals are predicted minus measured, log units.")
    print()
    for lab, fprot, getter in (("human", HP_FPROT, hp),
                               ("trout", TP_FPROT, tp)):
        print(f"  -- {lab} plasma")
        print(f"  {'chemical':15s}" + "".join(f"{c:>10s}" for c, _ in CW)
              + f"   {'bias@1uM':>9s}")
        res1 = []
        for nm in PFAS:
            row = []
            for _, c in CW:
                pred = math.log10(10 ** logD(bsa(nm), c) * fprot
                                  + 10 ** D_LIP[nm] * (1 - fprot))
                row.append(pred - logD(getter(nm), c))
            res1.append(row[2])
            star = "*" if nm in PREDICTED_LIP else " "
            print(f"  {nm:15s}" + "".join(f"{v:>+10.2f}" for v in row)
                  + f"   {10 ** row[2]:>8.2f}x{star}")
        amax = max(res1, key=abs)
        print(f"  {'geometric mean':15s}" + " " * 40
              + f"   {gm(res1):>8.2f}x")
        print(f"  {'mean |residual|':15s}" + " " * 40
              + f"   {10 ** (sum(abs(v) for v in res1) / len(res1)):>8.2f}x"
              f"   worst {10 ** amax:.2f}x")
        print()
    print("  * lipid coefficient predicted from chain length, not measured")
    print()
    print("  Both species come out near 1.7x on the geometric mean, inside")
    print("  the factor of 10 the paper claims. But the two means mean")
    print("  different things: every human residual but one is POSITIVE, so")
    print("  1.69x there is a systematic over-prediction, while the trout")
    print("  residuals change sign across the series and their mean is")
    print("  cancellation rather than accuracy -- hence the second row.")
    print("  Note also that the trout residuals run from +1.32 at PFBA to")
    print("  -0.30 at PFUnA: ordered by chain length, not scattered.")
    print()


def part4():
    rule()
    print("4. A PARAMETER-FREE BOUND ON WHAT COMPOSITION CAN EXPLAIN")
    rule()
    print("  The claim is that the species difference is composition, with")
    print("  albumin standing in for the protein of both plasmas. Write it")
    print("  out, with D_p the shared protein coefficient and D_l the shared")
    print("  lipid one:")
    print()
    print("      D_human     a_h D_p + b_h D_l")
    print("      -------  =  -----------------")
    print("      D_trout     a_t D_p + b_t D_l")
    print()
    print(f"      a_h {HP_FPROT:.4f}  b_h {1 - HP_FPROT:.4f}   (human)")
    print(f"      a_t {TP_FPROT:.4f}  b_t {1 - TP_FPROT:.4f}   (trout)")
    print()
    hi = HP_FPROT / TP_FPROT
    lo = (1 - HP_FPROT) / (1 - TP_FPROT)
    print("  That expression is monotone in D_p/D_l -- the derivative has")
    print(f"  sign a_h b_t - a_t b_h = "
          f"{HP_FPROT * (1 - TP_FPROT) - TP_FPROT * (1 - HP_FPROT):+.4f} > 0 "
          f"-- so it is")
    print("  confined between its two limits, whatever the coefficients are:")
    print()
    print(f"      D_p >> D_l  ->  a_h/a_t = {hi:.3f}x  ({math.log10(hi):+.3f} "
          f"log)")
    print(f"      D_l >> D_p  ->  b_h/b_t = {lo:.3f}x  ({math.log10(lo):+.3f} "
          f"log)")
    print()
    print(f"  So composition alone cannot produce a human/trout ratio outside")
    print(f"  {lo:.2f}x to {hi:.2f}x, for ANY compound at ANY concentration.")
    print("  No partition coefficient enters the bound; it is Table 3 and")
    print("  nothing else. Against the measured ratios:")
    print()
    print(f"  {'chemical':15s}{'nCF':>4s}" +
          "".join(f"{c:>11s}" for c, _ in CW))
    tally = {lab: [0, 0] for lab, _ in CW}
    for nm in PFAS:
        cells = []
        for lab, c in CW:
            rt = 10 ** (logD(hp(nm), c) - logD(tp(nm), c))
            out = rt > hi or rt < lo
            tally[lab][0] += 1
            tally[lab][1] += 1 if out else 0
            cells.append(f"{rt:>9.2f}x" + ("!" if out else " "))
        print(f"  {nm:15s}{T1[nm][0]:>4d}" + "".join(cells))
    print()
    print("  ! = outside what composition can explain")
    for lab, _ in CW:
        n, o = tally[lab]
        print(f"    at {lab:7s} {o} of {n} compounds outside the bound")
    print()
    print("  THE VERDICT. Composition is necessary and it is not sufficient.")
    print("  At 1 uM, within the measured range, 5 of 11 compounds sit")
    print("  outside a bound that no choice of partition coefficients can")
    print("  widen; PFBA's measured ratio of 23x exceeds the 1.32x ceiling")
    print(f"  by {23.4 / hi:.0f}-fold. Extrapolated down to 1 nM it is 7 of 11.")
    print("  The residuals are ordered by chain length rather than scattered")
    print("  (part 3), which is the signature of a missing mechanism, not of")
    print("  measurement error.")
    print()
    print("  What the paper's own test could not see. Its Figure 4 compares")
    print("  predicted against measured plasma partitioning at single")
    print("  concentrations and reports agreement 'within a factor of 10'.")
    print("  Our part 3 reproduces that: 1.69x for human, 1.67x for trout,")
    print("  comfortably inside. A factor of 10 on each ABSOLUTE level is")
    print("  simply too loose to test a RATIO whose whole admissible range")
    print(f"  is {lo:.2f}x to {hi:.2f}x -- a span of {hi / lo:.1f}x. The claim and")
    print("  the test were not on the same quantity.")
    print()
    print("  WHAT ELSE IT COULD BE, and this dataset cannot separate them:")
    print("    - temperature. Trout plasma was equilibrated at 25 C and")
    print("      human at 37 C. That breaks the shared-coefficient premise")
    print("      the bound rests on, so it is a live alternative and not a")
    print("      quibble: the bound tests 'composition at equal chemistry',")
    print("      and the two arms are not at equal temperature.")
    print("    - protein identity, which is what the paper argues against.")
    print("      Trout plasma binding is FLAT for 9 of 13 chemicals where")
    print("      human is not (part 1), and a weighted mean of a curved")
    print("      term and a constant is still curved -- so the shape, not")
    print("      just the level, needs something the composition model does")
    print("      not have. Fewer or weaker saturable sites on trout albumin")
    print("      would do it.")
    print("    - a sorptive phase in trout plasma that is neither of these")
    print("      two, which the same argument as section 6.7.3 would reach:")
    print("      a mixture outside the range of its components is not a")
    print("      mixing problem.")
    print()
    print("  Note which way this cuts. It does NOT weaken the paper's")
    print("  negative result -- that you cannot read a species difference")
    print("  off protein identity. It strengthens the positive one the")
    print("  project needs: a cross-species plasma-binding factor is not a")
    print("  constant, and composition alone will not supply it.")
    print()


def part5():
    rule()
    print("5. THE CONSEQUENCE FOR CROSS-SPECIES EXTRAPOLATION")
    rule()
    print("  Free fraction in UNDILUTED plasma, from the measured fits:")
    print()
    print("      f_u = V_w / (V_w + D_plasma/w x V_prot+lip)")
    print()
    print(f"      human  V_prot+lip {HP_VPL:.5f} L/L, V_w {1 - HP_VPL:.4f}")
    print(f"      trout  V_prot+lip {TP_VPL:.5f} L/L, V_w {1 - TP_VPL:.4f}")
    print()
    print("  Evaluated at a bioassay concentration and at a physiological")
    print("  one. The 1 nM column is an EXTRAPOLATION three decades below")
    print("  the fitted series and is shown to make the direction visible,")
    print("  not as a measurement.")
    print()
    print(f"  {'chemical':15s}{'human 1uM':>10s}{'trout 1uM':>10s}"
          f"{'ratio':>7s}  |{'human 1nM':>10s}{'trout 1nM':>10s}{'ratio':>7s}")
    for nm in PFAS:
        row = []
        for c in (1.0, 1e-3):
            fu_h = (1 - HP_VPL) / (1 - HP_VPL + 10 ** logD(hp(nm), c) * HP_VPL)
            fu_t = (1 - TP_VPL) / (1 - TP_VPL + 10 ** logD(tp(nm), c) * TP_VPL)
            row.append((fu_h, fu_t))
        print(f"  {nm:15s}{row[0][0]:>10.4f}{row[0][1]:>10.4f}"
              f"{row[0][1] / row[0][0]:>7.2f}  |"
              f"{row[1][0]:>10.5f}{row[1][1]:>10.5f}"
              f"{row[1][1] / row[1][0]:>7.2f}")
    print()
    print("  Read the two ratio columns against each other. The trout/human")
    print("  free-fraction ratio is not the same number at the two")
    print("  concentrations, and for several compounds it changes by more")
    print("  than the ratio itself is worth. A free fraction measured in one")
    print("  species at bioassay doses does not transfer to the other at")
    print("  exposure doses, and the error is not a constant that can be")
    print("  folded into an uncertainty factor -- it is a function of dose.")
    print()
    print("  This is the plasma-side counterpart of section 7's regime gap.")
    print("  Section 7 shows a single assay's own dose series moves the")
    print("  apparent albumin affinity by 2-4x. Here the same mechanism")
    print("  makes a SPECIES COMPARISON move with dose, which is worse,")
    print("  because an allometric or body-weight scaling step cannot see it.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    part1()
    part2()
    part3()
    part4()
    part5()
    return 0


if __name__ == "__main__":
    sys.exit(main())
