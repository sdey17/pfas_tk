#!/usr/bin/env python3
"""How different is one Argoul female mouse from the next?

The side quest was to plot individual serum concentration-time curves for
Argoul 2026, the way we did for the CPHEA raw curves.  The answer is that
the design forbids it, and this script says exactly why, then recovers the
inter-animal variability that *is* computable from the published numbers.

Everything below comes from Argoul et al., Environmental Research 303 (2026)
124802, Tables 1 and 2 plus the Methods.  No digitising, no simulation of
unmeasured noise: section 5's spread is the paper's own covariate model
evaluated at the two ends of its own observed diuresis range.

Run:  python3 scripts/argoul_individual_variability.py
      python3 scripts/argoul_individual_variability.py --figure
"""
import argparse
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")

# ----------------------------------------------------------------------------
# Published inputs
# ----------------------------------------------------------------------------

# Methods, 2.2 / 2.3: 108 female C57BL/6, 12 weeks, 20.0 +- 0.8 g.
N_TOTAL, N_IV, N_ORAL, N_CONTROL = 108, 51, 49, 8

# Methods, 2.3, the designed plasma sampling times: "15 and 30 min; 1, 2, 4,
# 8, 24 and 48 h; 3, 8, 15, 30, 60 and 119 days" -- 2 + 6 + 6 = 14 times, with
# 3 to 5 mice at each. 24 h and 48 h are distinct from the day-3 entry, so
# nothing is double-counted.
TIMES_D = [15 / 1440, 30 / 1440, 1 / 24, 2 / 24, 4 / 24, 8 / 24, 1.0, 2.0,
           3.0, 8.0, 15.0, 30.0, 60.0, 119.0]
N_TIMES = len(TIMES_D)

# Methods, 2.5: "In plasma, LLOQ ranged from 0.5 to 10 ng/mL". The Fig. 1
# caption distinguishes an LLOQ1 for PFOA and PFOS from an LLOQ2 for the
# other long-chain PFAS but prints neither number, so the bracket is used.
LLOQ_LO_NG_ML, LLOQ_HI_NG_ML = 0.5, 10.0

# Methods: "A urine volume covariate (ranging from 0.19 to 1.64 mL/24 h)"
UV_MIN, UV_MAX, UV_REF = 0.19, 1.64, 0.50

# Table 1: MW, IV dose (mg/kg), CL (mL/kg/day), Vss (L/kg), MRT (d)
# Table 2: urinary share of CL (%), fecal share (%), theta_UV with 95% CI,
#          and the paper's own illustrative CLu at UV = 0.2 and 2.5 mL/day.
# "cpt" is the structural model the paper selected for that compound.
PFAS = [
    # name      MW      dose    CL     Vss    MRT   CLu%  CLf%  theta   lo      hi     CLu@0.2 CLu@2.5 cpt
    ("PFBA",   214.04, 1.41,   783.0, 0.14,  0.18, 100.0,  0.0, -0.134, -0.586, 0.799, None,  None,   2),
    ("PFHxA",  314.05, 1.28,  6830.0, 4.00,  0.63,  99.3,  0.7, -1.040, -1.598, -0.278, 19400, 1970,  3),
    ("PFOA",   414.07, 0.080,    4.54, 0.089, 19.5, 79.1, 20.9,  0.240,  0.038, 0.550,  2.79,  5.3,   1),
    ("PFNA",   464.08, 0.11,     2.10, 0.11,  52.8, 28.6, 71.4,  0.546,  0.214, 0.989,  0.36,  1.6,   1),
    ("PFDA",   514.08, 0.077,    1.78, 0.12,  68.1,  7.3, 92.7,  1.215,  1.178, 1.251,  0.043, 0.90,  1),
    ("PFHxS",  400.12, 0.094,    1.34, 0.062, 46.5, 44.0, 56.0,  0.699,  0.359, 0.928,  0.30,  1.37,  1),
    ("PFOS",   500.13, 0.038,    2.68, 0.11,  43.1, 46.3, 53.7,  0.543,  0.261, 0.857,  0.71,  2.98,  1),
    ("PFDS",   600.15, 0.019,    5.15, 0.28,  53.3, 70.0, 30.0,  0.260, -0.013, 0.855,  None,  None,  1),
    ("GenX",   330.05, 0.052,  274.0, 0.095, 0.35,  98.2,  1.8,  0.0839, 0.082, 0.085,  249,   307,   2),
    ("PFO2OA", 346.05, 0.17,  1570.0, 0.48,  0.32,  97.6,  2.4, -0.828, -1.523, 0.294,  None,  None,  3),
]
# PFBS excluded from the covariate work: Table 2 records NA for its slope
# (nonlinear kinetics, CL 2589 -> 22 mL/kg/day across the concentration range).

FIELDS = ("name mw dose cl vss mrt clu_pct clf_pct theta lo hi "
          "clu_lo_paper clu_hi_paper cpt").split()


def rows():
    for rec in PFAS:
        yield dict(zip(FIELDS, rec))


def significant(r):
    """theta_UV is significant when its bootstrap 95% CI excludes zero."""
    return (r["lo"] > 0) == (r["hi"] > 0)


def clu_at(r, uv):
    """Argoul Eq. 1:  CLu_i = theta_CLu_pop * (UV / 0.5) ** theta_UV."""
    clu_pop = r["cl"] * r["clu_pct"] / 100.0
    return clu_pop * (uv / UV_REF) ** r["theta"]


def banner(txt):
    print()
    print(txt)
    print("-" * len(txt))


# ----------------------------------------------------------------------------
# 1. Why there are no individual curves
# ----------------------------------------------------------------------------

def section1():
    banner("1. The design forbids individual concentration-time curves")
    print("Methods 2.3, verbatim:")
    print('  "blood samples were taken immediately antemortem by puncture of')
    print('   the submandibular vein using a lancet in anesthetized mouse at')
    print('   designated time points after administration: 15 and 30 min; 1, 2,')
    print('   4, 8, 24, and 48 h; 3, 8, 15, 30, 60, and 119 days, from 3 to 5')
    print('   mice at each sampling time."')
    print()
    print("Antemortem means the blood draw is the terminal act for that animal.")
    print("So the arithmetic has to close with one plasma observation per mouse:")
    print()
    print(f"  mice in the study            {N_TOTAL}")
    print(f"    intravenous arm            {N_IV}")
    print(f"    oral arm                   {N_ORAL}")
    print(f"    untreated controls         {N_CONTROL}")
    print(f"    accounted for              {N_IV + N_ORAL + N_CONTROL}")
    print(f"  designed sampling times      {N_TIMES}")
    print(f"  IV   mice per time point     {N_IV / N_TIMES:.2f}   (paper says 3-5)")
    print(f"  oral mice per time point     {N_ORAL / N_TIMES:.2f}   (paper says 3-5)")
    print()
    print("If any mouse had been bled twice, 51 animals could not cover 13")
    print("times at 3-5 per time; the arms are sized for exactly one draw each.")
    print()
    print("Two things follow, and they pull in opposite directions:")
    print()
    print("  (a) No mouse has a concentration-time curve. The paper's own Fig. 1")
    print("      is captioned 'individual plasma concentrations versus time' --")
    print("      those are individual OBSERVATIONS, 3-5 of them stacked at each")
    print("      time, not individual trajectories. Nothing can be joined up.")
    print()
    print("  (b) Each mouse got the 11-PFAS cocktail, so one animal yields 11")
    print("      simultaneous measurements. That is one point in 11 dimensions,")
    print("      not 11 points in time, so it still buys no curve -- but it is")
    print("      why the paper could fit one NLME per compound on shared")
    print("      per-animal covariates.")
    print()
    print("This is the textbook sparse / destructive-sampling case, and it is")
    print("exactly why the authors went to nonlinear mixed effects instead of")
    print("fitting each animal: individual parameters are not identifiable from")
    print("one observation, so they are shrunk toward a population mean.")


# ----------------------------------------------------------------------------
# 2. Check that Eq. 1 as transcribed reproduces the paper's own numbers
# ----------------------------------------------------------------------------

def section2():
    banner("2. Does Eq. 1 as transcribed reproduce Table 2?")
    print("Table 2 prints CLu at UV = 0.2 and 2.5 mL/day for seven compounds.")
    print("Recomputing them from Table 1's CL, Table 2's CLu% and theta_UV is")
    print("the check that the covariate model has been read correctly.")
    print()
    print(f"{'PFAS':9s} {'CLu@0.2 calc':>13s} {'paper':>9s} {'ratio':>7s}  "
          f"{'CLu@2.5 calc':>13s} {'paper':>9s} {'ratio':>7s}")
    agree, total = 0, 0
    for r in rows():
        if r["clu_lo_paper"] is None:
            continue
        c_lo, c_hi = clu_at(r, 0.2), clu_at(r, 2.5)
        p_lo, p_hi = r["clu_lo_paper"], r["clu_hi_paper"]
        rl, rh = c_lo / p_lo, c_hi / p_hi
        print(f"{r['name']:9s} {c_lo:13.4g} {p_lo:9.4g} {rl:7.2f}  "
              f"{c_hi:13.4g} {p_hi:9.4g} {rh:7.2f}")
        total += 2
        agree += sum(1 for x in (rl, rh) if 0.9 <= x <= 1.1)
    print()
    print(f"{agree} of {total} recomputed cells land within 10% of the printed value.")
    print()
    print("Both misses are the paper's own internal inconsistency, not a")
    print("transcription error, and both are worth naming:")
    print()
    imp = -math.log(19400 / 1970) / math.log(2.5 / 0.2)
    print(f"  PFHxA  the printed pair falls 19400 -> 1970, a factor "
          f"{19400 / 1970:.2f}, which")
    print(f"         inverts to theta_UV = {imp:.3f} over the 12.5x change in UV.")
    print("         The same row prints -1.040 as the median bootstrap slope,")
    print(f"         and Eq. 1 at -1.040 gives a "
          f"{(2.5 / 0.2) ** 1.040:.1f}x fall, not {19400 / 1970:.2f}x.")
    hxs = dict(zip(FIELDS, PFAS[5]))
    print(f"  PFHxS  the printed 1.37 is what Eq. 1 returns at UV = 1.64 "
          f"({clu_at(hxs, 1.64):.3f}),")
    print("         the largest urine volume actually observed in the study,")
    print(f"         not at the 2.5 the column header gives "
          f"({clu_at(hxs, 2.5):.3f}).")
    print()
    print("So Table 2's two illustration columns are not all evaluated at the")
    print("same UV, or not all with the median slope. Five of the seven")
    print("compounds reproduce to within a few percent, which is enough to")
    print("trust Eq. 1 as transcribed. Sections 3-5 use the median slope")
    print("throughout, that being the estimate the paper reports with a CI,")
    print("and stay inside the observed UV range rather than extrapolating")
    print("to 2.5 mL/24 h, which no mouse in the study reached.")


# ----------------------------------------------------------------------------
# 3. Individual renal clearance across the observed diuresis range
# ----------------------------------------------------------------------------

def section3():
    banner("3. How far apart are two real mice in urinary clearance?")
    print(f"Observed urine output ran {UV_MIN} to {UV_MAX} mL/24 h, a "
          f"{UV_MAX / UV_MIN:.2f}x spread")
    print("between the driest and the wettest animal in the study. Eq. 1 turns")
    print("that into a clearance spread of (UV_max/UV_min) ** theta_UV.")
    print()
    print("This is measured inter-animal variability, not a model assumption:")
    print("urine volume was weighed per mouse per 24-h collection.")
    print()
    print(f"{'PFAS':9s} {'theta_UV':>9s} {'95% CI':>18s} {'sig':>4s} "
          f"{'CLu@0.19':>10s} {'CLu@1.64':>10s} {'fold':>8s}")
    span = UV_MAX / UV_MIN
    for r in rows():
        sig = "yes" if significant(r) else "no"
        lo, hi = clu_at(r, UV_MIN), clu_at(r, UV_MAX)
        fold = hi / lo
        tag = f"{fold:7.2f}x" if fold >= 1 else f"1/{1 / fold:5.2f}x"
        ci = "[%.3f, %.3f]" % (r["lo"], r["hi"])
        print(f"{r['name']:9s} {r['theta']:9.3f} {ci:>18s} {sig:>4s} "
              f"{lo:10.4g} {hi:10.4g} {tag:>8s}")
    print()
    print("Read the significant rows:")
    for r in rows():
        if not significant(r):
            continue
        fold = clu_at(r, UV_MAX) / clu_at(r, UV_MIN)
        direction = "rises" if fold > 1 else "falls"
        print(f"  {r['name']:7s} CLu {direction} {max(fold, 1 / fold):5.2f}x "
              f"from the driest to the wettest mouse")
    print()
    print("PFDA is the extreme: theta_UV = 1.215, so renal clearance is")
    print("slightly supralinear in urine flow and a 8.6x diuresis spread")
    print(f"becomes a {span ** 1.215:.1f}x spread in CLu. Its CI, [1.178, 1.251], is")
    print("also by far the tightest in the table, so this is not a wobbly")
    print("estimate. PFHxA runs the other way: theta_UV is negative, the only")
    print("compound for which more urine means less renal clearance.")
    print()
    print("The three non-significant slopes (PFBA, PFDS, PFO2OA) are the")
    print("compounds for which diuresis explains none of the between-mouse")
    print("difference. For PFBA and PFO2OA that is unsurprising -- essentially")
    print("the whole dose is gone in urine within 24 h either way, so there is")
    print("no room for flow to matter.")


# ----------------------------------------------------------------------------
# 4. Damping: total clearance moves much less than renal clearance
# ----------------------------------------------------------------------------

def total_cl(r, uv):
    """CL_i = CLu_i(UV) + CLf, holding the fecal limb at its population value."""
    clf = r["cl"] * r["clf_pct"] / 100.0
    return clu_at(r, uv) + clf


def section4():
    banner("4. Fecal elimination damps the diuresis signal")
    print("Only the renal limb responds to urine flow. The fecal limb is a")
    print("constant, so the fold-range in TOTAL clearance is pulled toward 1 in")
    print("proportion to how much of the compound leaves in feces.")
    print()
    print(f"{'PFAS':9s} {'CLu share':>10s} {'CLu fold':>9s} {'CL@0.19':>9s} "
          f"{'CL@1.64':>9s} {'CL fold':>8s} {'damping':>8s}")
    for r in rows():
        if not significant(r):
            continue
        cf = clu_at(r, UV_MAX) / clu_at(r, UV_MIN)
        t_lo, t_hi = total_cl(r, UV_MIN), total_cl(r, UV_MAX)
        tf = t_hi / t_lo
        damp = (max(cf, 1 / cf) - 1) / (max(tf, 1 / tf) - 1) if abs(tf - 1) > 1e-9 else float("inf")
        print(f"{r['name']:9s} {r['clu_pct']:9.1f}% {cf:8.2f}x {t_lo:9.4g} "
              f"{t_hi:9.4g} {tf:7.2f}x {damp:7.1f}x")
    print()
    print("PFDA makes the point hardest: a 13.7x swing in renal clearance")
    print("shows up as a 1.30x swing in total clearance, because 92.7% of PFDA")
    print("leaves in feces. PFHxS, split 44/56, keeps a 2.0x swing. PFOA, 79%")
    print("renal, keeps 1.5x from a 1.7x renal swing -- almost undamped.")
    print()
    print("So the compound where diuresis matters most for CLu is NOT the")
    print("compound where it matters most for the body burden. Ranking by")
    print("effect on total clearance:")
    out = []
    for r in rows():
        if not significant(r):
            continue
        tf = total_cl(r, UV_MAX) / total_cl(r, UV_MIN)
        out.append((max(tf, 1 / tf), r["name"]))
    for fold, name in sorted(out, reverse=True):
        print(f"  {name:8s} {fold:5.2f}x")


# ----------------------------------------------------------------------------
# 5. The fan: individual plasma concentration-time curves, as far as the
#    published numbers license them
# ----------------------------------------------------------------------------

def conc_umol(r, t, uv):
    """One-compartment IV: C(t) = D/V1 * exp(-CL t / V1), in umol/L."""
    v1 = r["vss"]                       # L/kg; for a 1-cpt model V1 == Vss
    cl = total_cl(r, uv) / 1000.0       # mL/kg/day -> L/kg/day
    c_mg_l = r["dose"] / v1 * math.exp(-cl * t / v1)
    return c_mg_l / r["mw"] * 1000.0


def onecpt_significant():
    return [r for r in rows() if r["cpt"] == 1 and significant(r)]


def section5():
    banner("5. The individual curves the published numbers do license")
    print("For the six compounds the paper fitted with a one-compartment model")
    print("(PFHxS, PFOA, PFOS, PFNA, PFDA, PFDS) the central volume IS Vss, so")
    print("Table 1 plus Table 2 fully specify a per-mouse plasma curve once you")
    print("know that mouse's urine output:")
    print()
    print("    C_i(t) = D/Vss * exp(-[CLu_pop*(UV_i/0.5)^theta + CLf] t / Vss)")
    print()
    print("Five of the six have a significant diuresis slope. Evaluating at the")
    print("driest (0.19 mL/24 h) and wettest (1.64) animal gives the envelope")
    print("that measured diuresis alone opens up between two mice.")
    print()
    print("A clearance difference is a difference in SLOPE on a semilog plot, so")
    print("the gap is not a fixed offset -- it fans out exponentially in time:")
    print()
    probe = [1.0, 8.0, 15.0, 30.0, 60.0, 119.0]
    hdr = "".join(f"{f'd{int(t)}':>9s}" for t in probe)
    print(f"{'PFAS':9s} {'MRT_d':>6s}{hdr}")
    for r in onecpt_significant():
        cells = ""
        for t in probe:
            fold = conc_umol(r, t, UV_MIN) / conc_umol(r, t, UV_MAX)
            cells += f"{fold:8.2f}x"
        print(f"{r['name']:9s} {r['mrt']:6.1f}{cells}")
    print()
    print("(ratio = dry mouse / wet mouse. Day 119 IS a designed plasma")
    print(" sampling time -- it is the last of the 14 -- so the right-hand")
    print(" column is not extrapolation in time. Whether it is observable is a")
    print(" separate question, settled next.)")
    print()
    print("Is the day-119 gap above the assay floor, or is the dry mouse the")
    print("only one still quantifiable? Plasma LLOQ ran 0.5-10 ng/mL:")
    print()
    print(f"{'PFAS':9s} {'C(119) dry':>12s} {'C(119) wet':>12s}   both > "
          f"{LLOQ_LO_NG_ML} / {LLOQ_HI_NG_ML:.0f} ng/mL?")
    for r in onecpt_significant():
        dry = conc_umol(r, 119.0, UV_MIN) * r["mw"] / 1000.0 * 1000.0
        wet = conc_umol(r, 119.0, UV_MAX) * r["mw"] / 1000.0 * 1000.0
        v_lo = "yes" if wet > LLOQ_LO_NG_ML else "no"
        v_hi = "yes" if wet > LLOQ_HI_NG_ML else "no"
        print(f"{r['name']:9s} {dry:10.3g}   {wet:10.3g}        "
              f"{v_lo:>3s} / {v_hi:>3s}")
    print()
    print("(ng/mL, from the same one-compartment solution)")
    print()
    print("PFHxS, PFNA and PFDA clear even the pessimistic 10 ng/mL floor in")
    print("both mice, so for those three the day-119 spread is an observable")
    print("spread. PFOS clears the optimistic 0.5 floor in both but not the")
    print("pessimistic one, so its day-119 figure depends which LLOQ applied.")
    print("PFOA clears neither in the wet mouse. That ordering is consistent")
    print("with the paper reporting PFOA recovered in urine 'within 60 days'")
    print("while PFHxS and PFDS run to 119, and with the M3 method for")
    print("below-quantification values having had real work to do. PFOA's 13x")
    print("at day 119 is therefore a model statement, not a measured one; its")
    print("day-60 figure of 3.7x is the defensible one.")
    print()
    print("Answering the question as asked -- how different is one female mouse")
    print("from another -- at the last time each compound is still above the")
    print("pessimistic floor in both animals:")
    for r in onecpt_significant():
        t = 119.0
        while t > 1.0:
            wet = conc_umol(r, t, UV_MAX) * r["mw"]
            if wet > LLOQ_HI_NG_ML:
                break
            t -= 1.0
        fold = conc_umol(r, t, UV_MIN) / conc_umol(r, t, UV_MAX)
        print(f"  {r['name']:7s} at day {t:5.0f}: {fold:5.2f}x between the driest "
              f"and wettest mouse")
    print()
    print("PFHxS is the headline: two mice that differ in nothing but their")
    print("measured urine output are predicted to differ 7.5x in plasma PFHxS")
    print("four months after the same dose, and both values are well above")
    print("even the pessimistic assay floor. PFOS is next at 4.5x by day 102")
    print("on the pessimistic floor, 5.8x by day 119 on the optimistic one.")
    print()
    print("Note which compound fans FASTEST, because it is not the one with the")
    print("biggest clearance spread. The fan rate is d(ln ratio)/dt = dCL/V1:")
    print()
    print(f"{'PFAS':9s} {'dCL':>8s} {'V1':>8s} {'fan rate':>11s} "
          f"{'doubling':>10s}")
    for r in onecpt_significant():
        dcl = total_cl(r, UV_MAX) - total_cl(r, UV_MIN)
        v1_ml = r["vss"] * 1000.0
        rate = dcl / v1_ml
        print(f"{r['name']:9s} {dcl:8.3f} {v1_ml:8.1f} {rate:9.4f}/d "
              f"{math.log(2) / rate:9.1f} d")
    print()
    print("PFOA fans fastest per day -- 0.022/day, a doubling of the gap every")
    print("32 days -- because its small central volume turns a modest clearance")
    print("difference into a steep slope difference. PFDA fans slowest despite")
    print("having the largest CLu spread of the five, because 93% fecal")
    print("elimination damps the total and its volume is larger. Chain length")
    print("is not what orders this column; the renal share and V1 are.")
    print()
    print("Two caveats, both of which make this an UNDER-statement of the real")
    print("spread:")
    print("  (i)  only the explained part of the variability is here. The paper")
    print("       put random effects on CLu, CLf and V1; this fan is the CLu")
    print("       covariate alone, with CLf and V1 held at population values.")
    print("  (ii) the fan is dry-vs-wet, i.e. the two extremes of 51 animals,")
    print("       so it is a range, not an SD.")


# ----------------------------------------------------------------------------
# 6. What is genuinely out of reach
# ----------------------------------------------------------------------------

def section6():
    banner("6. What the paper has that we do not")
    print("The unexplained between-subject variability is the quantity that")
    print("would answer 'how different are these mice' without going through")
    print("diuresis. The paper reports it:")
    print()
    print('  "Three random effects were included in the model: urinary')
    print('   clearance, fecal clearance, and central volume of distribution')
    print('   (V1) ... Random effects were modeled using exponential models.')
    print('   The between-subject variability (BSV%) was calculated using the')
    print('   coefficient of variation formula."')
    print()
    print("With P_i = P_pop * exp(eta_i), eta ~ N(0, omega^2), the CV formula is")
    print("  BSV% = 100 * sqrt(exp(omega^2) - 1)")
    print("so a reported BSV maps straight back to omega:")
    print()
    print(f"{'BSV%':>6s} {'omega':>7s} {'5th-95th pct spread':>20s}")
    for bsv in (20, 30, 40, 50, 75, 100):
        om = math.sqrt(math.log1p((bsv / 100.0) ** 2))
        spread = math.exp(2 * 1.645 * om)
        print(f"{bsv:6d} {om:7.3f} {spread:19.1f}x")
    print()
    print("Those omegas are in Supplementary Tables S5-S26, one per compound")
    print("and matrix, which are not in hand -- the main text only states that")
    print("the relative standard errors are mostly below 20%, which is the")
    print("precision of the population estimates, a different quantity. Table 1")
    print("and Table 2 carry no omega.")
    print()
    print("Worth noting what the paper says it did with them: random effects")
    print('were "retained or excluded based on biological relevance and')
    print('acceptable shrinkage values", and incorporating urinary volume')
    print('"improved the model and reduced inter-individual variability for all')
    print('PFAS, with the exception of PFBS". That last sentence is the direct')
    print("statement that the fan in section 5 is a real slice of the total:")
    print("adding the UV covariate moved variance out of omega and into an")
    print("explained covariate effect for 10 of the 11 compounds.")
    print()
    print("If the SI turns up, the one line to pull is omega for CLu, CLf and")
    print("V1 per compound. Combined with section 5 that would split the")
    print("between-mouse spread into a diuresis-explained part and a residual,")
    print("which is the quantity a PBPK model would need for a Monte Carlo.")


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------

def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    BLUE, ORANGE = "#2a78d6", "#eb6834"
    SURFACE = "#fcfcfb"
    INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
    RULE = "#d8d7d3"
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.edgecolor": RULE, "axes.linewidth": 0.8,
        "axes.labelcolor": INK2, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.fontsize": 8.5,
    })

    comps = onecpt_significant()
    fig, axes = plt.subplots(1, len(comps), figsize=(2.85 * len(comps), 3.6))
    grid = [i / 8.0 for i in range(1, 8 * 120)]

    for ax, r in zip(axes, comps):
        dry = [conc_umol(r, t, UV_MIN) for t in grid]
        wet = [conc_umol(r, t, UV_MAX) for t in grid]
        mid = [conc_umol(r, t, UV_REF) for t in grid]

        # assay floor, converted from ng/mL to umol/L for this compound
        f_lo = LLOQ_LO_NG_ML / r["mw"]
        f_hi = LLOQ_HI_NG_ML / r["mw"]
        ax.axhspan(f_lo, f_hi, color=MUTED, alpha=0.18, lw=0, zorder=0)

        ax.fill_between(grid, wet, dry, color=BLUE, alpha=0.16, lw=0, zorder=2)
        ax.plot(grid, dry, color=BLUE, lw=1.7, zorder=3)
        ax.plot(grid, wet, color=ORANGE, lw=1.7, zorder=3)
        ax.plot(grid, mid, color=MUTED, lw=1.0, ls=(0, (4, 2)), zorder=3)

        ax.set_yscale("log")
        ax.set_xlim(-4, 124)
        ax.set_ylim(min(f_lo * 0.35, min(wet)), max(dry) * 2.2)

        # the 14 designed plasma sampling times, as a rug at the axis foot
        y0 = ax.get_ylim()[0]
        for t in TIMES_D:
            ax.plot([t, t], [y0, y0 * 1.45], color=MUTED, lw=0.9, zorder=1)

        fold = conc_umol(r, 119.0, UV_MIN) / conc_umol(r, 119.0, UV_MAX)
        ax.set_title(f"{r['name']}\n"
                     f"CLu {r['clu_pct']:.0f}% of CL, "
                     rf"$\theta_{{UV}}$ = {r['theta']:.2f}",
                     fontsize=9.5, color=INK, pad=6)
        ax.annotate(f"{fold:.1f}x apart at day 119", xy=(0.96, 0.955),
                    xycoords="axes fraction", ha="right", va="top",
                    fontsize=8.5, color=INK2)
        ax.set_xlabel("days after IV dose")
        ax.set_xticks([0, 30, 60, 90, 119])

    axes[0].set_ylabel(r"plasma concentration ($\mu$mol/L)")

    handles = [
        plt.Line2D([], [], color=BLUE, lw=1.7,
                   label="driest mouse, 0.19 mL urine/24 h"),
        plt.Line2D([], [], color=ORANGE, lw=1.7,
                   label="wettest mouse, 1.64 mL/24 h"),
        plt.Line2D([], [], color=MUTED, lw=1.0, ls=(0, (4, 2)),
                   label="median mouse, 0.50 mL/24 h"),
        plt.Rectangle((0, 0), 1, 1, color=MUTED, alpha=0.18,
                      label="plasma LLOQ band, 0.5-10 ng/mL"),
        plt.Line2D([], [], color=MUTED, lw=0.9,
                   label="the 14 designed sampling times"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, -0.07))

    fig.suptitle(
        "Argoul 2026: the between-mouse spread that measured diuresis alone "
        "accounts for\n"
        "one-compartment PFAS, from Table 1 + Table 2 Eq. 1. No mouse in the "
        "study has more than one plasma point, so these are the individual "
        "curves the design permits, not observed trajectories.",
        fontsize=10, color=INK, y=1.04)
    fig.tight_layout()
    path = os.path.join(OUT, "fig16_argoul_individual_spread.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"wrote {os.path.abspath(path)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", action="store_true")
    args = ap.parse_args()
    if args.figure:
        figure()
        return
    print("=" * 78)
    print("ARGOUL 2026: INDIVIDUAL-MOUSE VARIABILITY, AND WHY THERE ARE NO")
    print("INDIVIDUAL CURVES")
    print("=" * 78)
    section1()
    section2()
    section3()
    section4()
    section5()
    section6()
    print()


if __name__ == "__main__":
    main()
