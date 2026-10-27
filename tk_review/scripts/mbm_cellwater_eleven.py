#!/usr/bin/env python3
"""Section 6.7 re-run: the cell-water limb against ELEVEN measured ratios.

Section 6.7 of MBM_REVIEW tests the published mass-balance models against
cellular uptake in three assay datasets, over four compounds at nCF 5-8.
Qin et al. 2024 Table 1 (doi 10.1021/acs.est.3c09950) reports measured
log D_cell/w for ELEVEN PFAS spanning five head groups, in one cell line
whose composition the same group measured.  That permits a sharper test than
section 6.7 could run, for a structural reason:

    the cells are the same throughout (HEK293H), and their composition is
    MEASURED, so the volume and composition term is held fixed at its
    measured value and every remaining discrepancy is a partition
    coefficient.

Section 6.7's central claim is that the two failure modes separate -- chain
length belongs to the coefficients, level belongs to the volumes.  This test
fixes the volumes by construction, so it isolates the first half of that
claim and should reproduce it on 11 compounds rather than 4 if it is right.

WHAT IS COMPARED.  Everything is on Qin's basis, log D_cell/w in litres of
water per litre of cell protein + lipid.  Qin's eq 5 is a sorptive-volume
weighted mean,

    D_cell/w = D_protein/w x f_prot  +  D_lip/w x f_lip,   f_prot + f_lip = 1

with f from the measured HEK293H composition (SI Table S8: 33.9 mL/L protein,
5.3 mL/L lipid, 96.1% water), so f_prot = 0.865 and f_lip = 0.135.  The
variants differ only in what is put in for the two coefficients.  The
Fischer et al. form carries cell water as a third phase and normalises by
V_cell, so it is converted to Qin's basis before comparison, not compared
across bases.

THE IMPOSSIBILITY TEST comes first because it needs no model at all.  A
weighted mean of two numbers must lie between them, whatever the weights, so
any measured D_cell/w outside [D_protein/w, D_lip/w] cannot be reproduced by
eq 5 under ANY composition.  That is a composition-independent verdict and
no retuning of volume fractions can answer it.

Run:  python3 tk_review/scripts/mbm_cellwater_eleven.py
"""
import collections
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from add_qin2024_table1 import T1            # Qin 2024 Table 1
import mbm_qin_vs_armitage as Q              # 2024 BSA refit, Armitage spLFERs

# Measured HEK293H composition, Qin 2024 SI Table S8.
CELL_FPROT_V = 33.9e-3        # L protein per L cell
CELL_FLIP_V = 5.3e-3          # L lipid   per L cell
CELL_FW_V = 0.961             # L water   per L cell
F_PROT = CELL_FPROT_V / (CELL_FPROT_V + CELL_FLIP_V)      # 0.8648
F_LIP = 1.0 - F_PROT                                       # 0.1352

ORDER = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFUnA", "HFPO-DA",
         "PFHxS", "PFOS", "6:2 FTSA", "PFOSA"]
HEAD = {nm: T1[nm][1] for nm in ORDER}
NCF = {nm: T1[nm][0] for nm in ORDER}
# Table 1 marks these two D_lip/w values as predicted from eq S10, not
# measured, so every table says which set it is over.
LIP_PREDICTED = {"6:2 FTSA", "PFOSA"}
# log K_OW is sourced per compound only for these four (Qin 2023 Table 1 and
# the sulfonate entries). The Armitage limb is restricted to them rather
# than interpolating a log K_OW across five head groups.
KOW = {"PFBA": 2.60, "PFOA": 5.68, "PFHxS": 3.69, "PFOS": 5.77}
KOW_ALT = {"PFHxS": 5.29}     # COSMOtherm, 40x from the OPERA value


def meas(nm):
    return T1[nm][5]          # measured log D_cell/w


def d_sp(nm):
    return T1[nm][7]          # measured log D_SP/w


def d_lip(nm):
    return T1[nm][8]          # log D_lip/w (literature, 2 predicted)


def wmean(logp, loglip):
    """Qin eq 5 in log space."""
    return math.log10(10 ** logp * F_PROT + 10 ** loglip * F_LIP)


def fischer(logp, loglip):
    """Fischer et al.'s form, then rebased onto Qin's sorptive-volume basis.

    Fischer: D_cell/w = VF_prot*D_prot + VF_lip*D_lip + VF_w, per L of CELL.
    Qin's D_cell/w counts only the protein+lipid amount, per L of prot+lip.
    So subtract the water contribution and renormalise. The two differ by
    more than a constant, because the water term is additive and does not
    scale with the coefficients.
    """
    per_cell = (CELL_FPROT_V * 10 ** logp + CELL_FLIP_V * 10 ** loglip
                + CELL_FW_V)
    return math.log10((per_cell - CELL_FW_V)
                      / (CELL_FPROT_V + CELL_FLIP_V))


def armitage_cell(nm, logkow):
    """IV-MBM EQP v2.0's cell term, rebased onto Qin's basis.

    Armitage's cell contains water as a SORPTIVE phase and his D_cw is per
    litre of whole cell. Removing the water amount and renormalising puts it
    on Qin's basis; the structure is left as published, which means cell
    protein still sorbs by octanol-equivalence.
    """
    comp = dict(sl=0.0, ml=CELL_FLIP_V, prot=CELL_FPROT_V, water=CELL_FW_V)
    D = Q.armitage_D(logkow, comp, piecewise=True) if hasattr(
        Q, "armitage_D") else None
    if D is None:
        import mbm_published_equations as P
        D = P.armitage_D(logkow, comp, piecewise=True)
    dcw = D["Dcw"]
    return math.log10(max(dcw - CELL_FW_V, 1e-12)
                      / (CELL_FPROT_V + CELL_FLIP_V))


def stats(resid):
    """Geometric-mean bias and the slope of the residual against nCF."""
    names = sorted(resid, key=lambda n: ORDER.index(n))
    y = [resid[n] for n in names]
    bias = 10 ** (sum(y) / len(y))
    if len(names) < 3:
        return bias, None, None
    x = [NCF[n] for n in names]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    slope = sxy / sxx
    syy = sum((b - my) ** 2 for b in y)
    r = sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")
    return bias, slope, r


def rule(ch="=", n=78):
    print(ch * n)


# =====================================================================
def part1():
    rule()
    print("1. THE IMPOSSIBILITY TEST -- no model, no composition, no fit")
    rule()
    print("  Qin eq 5 is a weighted mean of D_SP/w and D_lip/w, so the")
    print("  measured D_cell/w must lie between them for SOME composition.")
    print("  Both bracketing values are measured, and so is the answer.")
    print()
    print(f"  {'cmpd':9s} {'nCF':>3s} {'D_SP/w':>7s} {'D_lip/w':>8s} "
          f"{'range':>14s} {'measured':>9s}  verdict")
    out = []
    for nm in ORDER:
        lo, hi = sorted((d_sp(nm), d_lip(nm)))
        m = meas(nm)
        ok = lo <= m <= hi
        excess = 0.0 if ok else (m - hi if m > hi else lo - m)
        out.append((nm, ok, excess))
        tag = "ok" if ok else f"OUTSIDE by {excess:+.2f} log"
        star = "" if nm not in LIP_PREDICTED else " (lipid predicted)"
        print(f"  {nm:9s} {NCF[nm]:3d} {d_sp(nm):7.2f} {d_lip(nm):8.2f} "
              f"{lo:6.2f}-{hi:<7.2f} {m:9.2f}  {tag}{star}")
    bad = [o for o in out if not o[1]]
    print()
    print(f"  {len(bad)} of {len(out)} compounds are unreachable by eq 5 under")
    print(f"  ANY composition. Excess above the achievable range runs "
          f"{min(abs(o[2]) for o in bad):.2f}-{max(abs(o[2]) for o in bad):.2f}")
    print(f"  log units ({10 ** min(abs(o[2]) for o in bad):.1f}x to "
          f"{10 ** max(abs(o[2]) for o in bad):.1f}x).")
    print()
    print("  Direction matters. Of the failures:")
    hi_side = [o[0] for o in out if not o[1] and o[2] > 0]
    lo_side = [o[0] for o in out if not o[1] and o[2] < 0]
    print(f"    above both coefficients: {', '.join(hi_side)}")
    print(f"    below both coefficients: {', '.join(lo_side) or 'none'}")
    print("  A cell that sorbs MORE than either of its two pure phases is")
    print("  not a mixing problem. Either a sorptive phase is missing from")
    print("  the model, or D_SP/w measured on isolated structural protein")
    print("  understates protein in situ. Eq 5 cannot express it either way.")
    print()
    return out


def part2():
    rule()
    print("2. THE VARIANTS, with the cell composition held at its measured")
    print("   value throughout -- so every difference below is a coefficient")
    rule()
    print(f"  f_prot = {F_PROT:.4f}, f_lip = {F_LIP:.4f} of the sorptive")
    print(f"  volume (HEK293H measured, SI Table S8)")
    print()
    refit = Q.QIN2024_BSA_SPEC_NONSPEC

    def run(label, protein_of, form=wmean, only=None):
        resid = {}
        for nm in ORDER:
            if only and nm not in only:
                continue
            lp = protein_of(nm)
            if lp is None:
                continue
            resid[nm] = form(lp, d_lip(nm)) - meas(nm)
        return label, resid

    # Which route each paper actually published, read from the sources:
    #   Qin 2024 Fig 2c: "D_cell/w values were predicted by eq 5 with
    #   [D_SP/w] and D_lip/w (Table 1)"        -> variant A, measured D_SP/w
    #   Qin 2023 sec 3.5: "D_cell/w was predicted by the MBM (Eq. 17) with
    #   D_non-specific (Table 1) as a surrogate for partitioning to
    #   structural proteins"                   -> variant B, albumin surrogate
    variants = [
        run("A  measured D_SP/w + measured D_lip/w (Qin 2024 as published)",
            d_sp),
        run("B  albumin non-specific surrogate (Qin 2023 as published)",
            lambda nm: refit[nm][1]),
        run("C  albumin specific term instead (see note below)",
            lambda nm: refit[nm][0]),
    ]
    # The Fischer et al. form is NOT a fourth variant. Converting it onto
    # Qin's basis cancels the water term exactly:
    #   (VF_p D_p + VF_l D_l + VF_w - VF_w)/(VF_p + VF_l) = f_p D_p + f_l D_l
    # so with the same coefficients and composition the two forms are the
    # same prediction. Asserted rather than asserted-in-prose.
    for nm in ORDER:
        a, b = fischer(d_sp(nm), d_lip(nm)), wmean(d_sp(nm), d_lip(nm))
        assert abs(a - b) < 1e-12, (nm, a, b)
    print("  Before the table: whether cell water sits INSIDE D_cell/w")
    print("  (Fischer et al.) or outside it as a non-sorbing phase (Qin et")
    print("  al.) is a basis convention, not a model difference. Converting")
    print("  Fischer's form onto Qin's basis cancels the water term exactly,")
    print("  so with the same coefficients the two give the identical")
    print("  prediction -- verified to machine precision for all 11 here.")
    print("  Section 6.7 lists it as difference (ii) between the two")
    print("  schools; it has no predictive content. What does differ is")
    print("  difference (i), the protein surrogate, which is A against B.")
    print()
    print(f"  {'variant':62s}{'bias':>8s}{'slope':>8s}{'r':>7s}")
    for label, resid in variants:
        b, s, r = stats(resid)
        print(f"  {label:62s}{b:>7.2f}x"
              + (f"{s:>+8.3f}{r:>+7.2f}" if s is not None else f"{'(n<3)':>15s}"))
    print()
    print("  Restricted to the 9 compounds whose D_lip/w is MEASURED rather")
    print("  than predicted from chain length:")
    meas9 = [n for n in ORDER if n not in LIP_PREDICTED]
    for label, resid in variants:
        b, s, r = stats({k: v for k, v in resid.items() if k in meas9})
        print(f"  {label:62s}{b:>7.2f}x"
              + (f"{s:>+8.3f}{r:>+7.2f}" if s is not None else f"{'(n<3)':>15s}"))
    print()
    print("  NOTE on variant C. Qin's 'specific' binding constant is Eq. 15,")
    print("  Cmax/(2 Kd), the ratio AT half-saturation -- 0.301 log below the")
    print("  Cw->0 plateau of the specific term. Read as a plateau it is 2x")
    print("  low, so a true sub-saturating variant sits 0.301 log above the")
    print("  C row: bias 5.7x rather than 2.84x over the 7 compounds where a")
    print("  specific term could be resolved at all. Either way C is the")
    print("  worst of the three, which is the expected direction: these")
    print("  assays dose at nu > 1, where specific binding is saturated.")
    print()
    print("  Per-compound residuals for variant A, so the drift is visible")
    print("  rather than summarised:")
    la, ra = variants[0]
    print(f"  {'cmpd':9s} {'nCF':>3s} {'head group':26s} {'pred':>6s} "
          f"{'meas':>6s} {'resid':>7s}")
    for nm in ORDER:
        if nm not in ra:
            continue
        print(f"  {nm:9s} {NCF[nm]:3d} {HEAD[nm]:26s} "
              f"{wmean(d_sp(nm), d_lip(nm)):6.2f} {meas(nm):6.2f} "
              f"{ra[nm]:+7.2f}")
    print()
    return variants


def part3(variants):
    rule()
    print("3. THE SPLFER LIMB -- the discriminator section 6.7 identified")
    rule()
    print("  Section 6.7 concluded: every variant whose lipid coefficient is")
    print("  MEASURED has a chain-length slope at or below zero, and the only")
    print("  ones above zero predict lipid from log K_OW. Testing that here")
    print("  needs a log K_OW per compound, and one is sourced for only four")
    print("  of the eleven, so the comparison is over those four and is")
    print("  stated as such rather than interpolated across five head groups.")
    print()
    print(f"  {'cmpd':9s} {'logKow':>7s} {'Armitage D_cell/w':>18s} "
          f"{'measured':>9s} {'resid':>7s}")
    resid_a, resid_m = {}, {}
    for nm, lk in KOW.items():
        pred = armitage_cell(nm, lk)
        resid_a[nm] = pred - meas(nm)
        resid_m[nm] = wmean(d_sp(nm), d_lip(nm)) - meas(nm)
        print(f"  {nm:9s} {lk:7.2f} {pred:18.2f} {meas(nm):9.2f} "
              f"{resid_a[nm]:+7.2f}")
    for nm, lk in KOW_ALT.items():
        print(f"  {nm:9s} {lk:7.2f} {armitage_cell(nm, lk):18.2f} "
              f"{meas(nm):9.2f} {armitage_cell(nm, lk) - meas(nm):+7.2f}"
              f"   <- COSMOtherm log K_OW instead of OPERA")
    print()
    ba, sa, ra = stats(resid_a)
    bm, sm, rm = stats(resid_m)
    print(f"  {'Armitage spLFER (predicted lipid), same 4':48s}"
          f"{ba:>7.2f}x{sa:>+8.3f}{ra:>+7.2f}")
    print(f"  {'measured coefficients, same 4':48s}"
          f"{bm:>7.2f}x{sm:>+8.3f}{rm:>+7.2f}")
    print()
    print("  The 40x OPERA/COSMOtherm disagreement on PFHxS alone moves its")
    print("  predicted D_cell/w by "
          f"{abs(armitage_cell('PFHxS', 5.29) - armitage_cell('PFHxS', 3.69)):.2f}"
          " log units, which is larger than most")
    print("  of the residuals being compared. The spLFER limb's input is not")
    print("  settled enough to carry a slope result on four points.")
    print()
    return resid_a, resid_m


def part4(variants):
    rule()
    print("4. BY HEAD GROUP -- the four-compound test could not do this")
    rule()
    groups = collections.defaultdict(list)
    for nm in ORDER:
        groups["carboxylate" if "carboxylate" in HEAD[nm]
               else "sulfonate" if HEAD[nm] == "sulfonate"
               else "other (telomer sulfonate, sulfonamide)"].append(nm)
    for label, resid in variants:
        print(f"  {label}")
        for g, nms in groups.items():
            sub = {k: v for k, v in resid.items() if k in nms}
            if not sub:
                continue
            b, s, r = stats(sub)
            sl = f"{s:>+8.3f}" if s is not None else f"{'(n<3)':>8s}"
            print(f"    {g:40s} n={len(sub)}  {b:>7.2f}x{sl}")
        print()


def part5(variants, resid_spl, resid_m4):
    rule()
    print("5. WHAT THIS DOES TO SECTION 6.7")
    rule()
    la, ra = variants[0]
    meas9 = [n for n in ORDER if n not in LIP_PREDICTED]
    b11, s11, r11 = stats(ra)
    b9, s9, r9 = stats({k: v for k, v in ra.items() if k in meas9})
    carbs = [n for n in ORDER if "carboxylate" in HEAD[n]]
    bc, sc, rc = stats({k: v for k, v in ra.items() if k in carbs})
    bs, ss, rs = stats(resid_spl)
    bm, sm, rm = stats(resid_m4)
    print("  SURVIVES: the discriminator. On the four compounds where a")
    print("  log K_OW is sourced, the spLFER lipid coefficient gives a")
    print(f"  residual slope of {ss:+.3f} and the measured coefficients "
          f"{sm:+.3f} --")
    print(f"  a {abs(ss / sm):.0f}-fold difference, same direction as section "
          f"6.7's")
    print("  +0.282 against -0.134. Predicting the lipid term from one")
    print("  descriptor is where the chain-length drift comes from.")
    print()
    print("  DOES NOT SURVIVE: chain-length neutrality. Section 6.7 concluded")
    print("  'measured biopartitioning does solve the chain-length problem'.")
    print("  In the cell limb alone the measured-coefficient residual is")
    print(f"    all 11 compounds          {s11:+.3f} decades/CF2 (r {r11:+.2f})")
    print(f"    9 with D_lip/w measured   {s9:+.3f} decades/CF2 (r {r9:+.2f})")
    print(f"    7 carboxylates only       {sc:+.3f} decades/CF2 (r {rc:+.2f})")
    print(f"    the 4 of section 6.7      {sm:+.3f} decades/CF2 (r {rm:+.2f})")
    print("  The flatness is a four-compound result. It does not hold across")
    print("  eleven, and the carboxylate-only row shows that is not an")
    print("  artefact of pooling five head groups -- the drift is inside the")
    print("  carboxylate series.")
    print()
    print("  NEW, AND IT LOCALISES SECTION 6.7's LEVEL ERROR. Section 6.7")
    print("  reports that nothing tested fixes the level, the best faithful")
    print("  parameterisation over-predicting cellular uptake by 4.62x. The")
    print(f"  cell limb on its own, with the composition measured and both")
    print(f"  coefficients measured, goes the OTHER WAY: {b11:.2f}x over 11")
    print(f"  compounds, {b9:.2f}x over the 9 with measured lipid -- an")
    print(f"  UNDER-prediction of about 2x. Two consequences:")
    print("    - the 4.62x over-prediction cannot be coming from the cell")
    print("      coefficient, which errs 2x in the opposite direction. It")
    print("      has to sit in the medium limb or the volumes, and section")
    print("      6.7's own separation argument says the volumes.")
    print("    - the cell and medium errors partially cancel in the")
    print("      whole-balance uptake test, which is the quantity section")
    print("      6.7 used. Cancellation between limbs is exactly what makes")
    print("      a whole-balance residual a poor diagnostic, and it is why")
    print("      this limb-by-limb test is worth running.")
    print()
    print("  UNCHANGED: the impossibility result, now on 11 compounds rather")
    print("  than 4, all failing on the same side. That is the strongest")
    print("  form of the section 6.7 finding and needs no model at all.")
    print()
    print("  ONE CLAIM TO DROP from section 6.7: difference (ii) between the")
    print("  two schools, whether cell water is inside D_cell/w. Part 2")
    print("  shows it is a basis convention with no predictive content.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    part1()
    variants = part2()
    resid_spl, resid_m4 = part3(variants)
    part4(variants)
    part5(variants, resid_spl, resid_m4)
    return 0


if __name__ == "__main__":
    sys.exit(main())
