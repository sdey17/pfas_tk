#!/usr/bin/env python3
"""Rat renal clearance for PFBA, PFBS, PFDA and PFOS: search and bounds (Stage 18).

    python3 -I scripts/stage18_rat_renal_clearance_search.py

Writes db/rat_renal_clearance_search.csv, db/kudo_ohmori_crossvalidation.csv
and db/renal_only_bounds.csv. Idempotent.

WHAT WAS ASKED. Equation (5), the renal-only decomposition of Stage 17, is
exact but covers only three compounds (PFHpA, PFOA, PFHxS) because it needs a
DIRECT renal clearance in male human AND male rat. This stage asks whether
PFBA, PFBS, PFDA or PFOS can be added.

WHAT THE SEARCH FOUND. No directly measured rat renal clearance exists in the
accessible literature for any of the four. The reasons differ by compound and
are recorded individually in RAT_SEARCH below, because three of them are
retrievable in principle and one is not:

  PFDA  Ohmori 2003 Fig. 2 DID measure it, in both sexes, but both bars sit
        on the baseline. The figure yields only CL_R < 15 mL/kg/day, hence
        T < 1.04, which straddles unity and settles nothing.
  PFOS  No rat urine study found. The NTP study (Huang 2019) collected
        plasma, liver, kidney and brain only -- no urine, so no renal
        clearance can come from it at any dose. Harada 2005's abstract
        mentions "sex differences in renal clearance ... reported in rats",
        but its Table 5 carries rat rows for PFOA ONLY; PFOS appears there
        for cynomolgus monkey and human. That lead resolves to PFOA, which
        the project already holds.
  PFBA  Chang 2008 measured rat urine (its Methods list serum, liver, urine
        and feces) but reports no renal clearance in the abstract. PAYWALLED.
  PFBS  Olsen 2009 measured rat urine and states "urine appeared to be a
        major route of elimination", but reports no renal clearance in the
        abstract. PAYWALLED.

THE HUMAN SIDE IS THE TIGHTER CONSTRAINT, AND THIS WAS NOT EXPECTED. Even a
perfect rat number would not admit PFBA or PFBS to equation (5), because no
human renal clearance exists for either. PFDA and PFOS do have human renal
clearances (Fujii 2015, Zhang 2013, Harada 2005). So the two compounds whose
rat PDFs would help are the two that cannot be used, and the two that can be
used are limited by the rat.

WHAT CAN STILL BE DONE. For PFDA and PFOS the rat renal clearance can be
bounded ABOVE, and a one-sided bound on CL_R gives a one-sided bound on every
term of equation (5). Two different bounds are used:

  PFDA  Kudo 2001 reports 0.2% of an IP dose in urine over 120 h in male
        Wistar rats. Over a window tau, A_urine = CL_R * AUC(0,tau) and
        Dose = CL_tot * AUC(0,inf), so
            CL_R = (A_urine/Dose) * CL_tot / [AUC(0,tau)/AUC(0,inf)].
        With a one-compartment curve the bracket is 1-exp(-ln2*tau/t_half).
        For any multi-exponential curve with the SAME terminal half-life the
        bracket is LARGER -- it is a weighted average of (1-exp(-k*tau)) over
        rate constants all of which exceed the terminal one -- so dividing by
        the one-compartment value OVERSTATES CL_R. The result is an upper
        bound, which is the direction needed.
  PFOS  CL_R <= CL_tot identically. This is loose, because the rat excretes
        PFOS substantially in faeces, so the true renal clearance is well
        below the bound.

THE PFDA ROUTE IS CROSS-VALIDATED, NOT ASSUMED. Four compounds have BOTH a
Kudo 2001 urinary fraction and an Ohmori 2003 renal clearance, so the same
arithmetic can be run where the answer is already known. It reads high by
1.05x to 1.80x, and high in all four cases. A systematic overstatement is
exactly what an upper bound tolerates, so the validation supports the bound's
direction while denying it any claim to be a point estimate.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

GFR = {"human": 2570.0, "rat": 14400.0}

# ---------------------------------------------------------------- the search
# verdict is one of: measured_but_unresolvable, no_urine_study, paywalled
RAT_SEARCH = [
    dict(chemical="PFBA", verdict="paywalled",
         rat_clr="none published",
         measured_urine="yes (Chang 2008: serum, liver, urine, feces)",
         blocking_pdf="Chang SC, Das K, Ehresman DJ, et al. Comparative "
                      "pharmacokinetics of perfluorobutyrate in rats, mice, "
                      "monkeys, and humans and relevance to human exposure "
                      "via drinking water. Toxicol Sci 2008;104(1):40-53. "
                      "doi:10.1093/toxsci/kfn057",
         human_clr="NONE -- no human renal clearance for PFBA in the project",
         admissible_to_eq5="no, blocked on BOTH sides",
         note="Rat half-life is 6-9 h, so a 24 h urine collection spans many "
              "half-lives and the AUC correction would be ~1. If the PDF "
              "gives cumulative urinary recovery and total clearance, CL_R "
              "follows almost exactly. But the human side is empty, so this "
              "would extend the rat table, not equation (5)."),
    dict(chemical="PFBS", verdict="paywalled",
         rat_clr="none published",
         measured_urine="yes, in TWO papers. Olsen 2009: serum, liver, urine, "
                        "feces, 'urine appeared to be a major route of "
                        "elimination'. Chengelis 2009: urinalysis, and states "
                        "'there was no appreciable difference in the extent "
                        "or rate of urinary elimination between compounds or "
                        "genders'",
         blocking_pdf="Olsen GW, Chang SC, Noker PE, Gorman GS, Ehresman DJ, "
                      "Lieder PH, Butenhoff JL. A comparison of the "
                      "pharmacokinetics of perfluorobutanesulfonate (PFBS) "
                      "in rats, monkeys, and humans. Toxicology "
                      "2009;256(1-2):65-74. doi:10.1016/j.tox.2008.11.008"
                      "  ||  Chengelis CP, Kirkpatrick JB, Myers NR, "
                      "Shinohara M, Stetson PL, Sved DW. Comparison of the "
                      "toxicokinetic behavior of perfluorohexanoic acid "
                      "(PFHxA) and nonafluorobutane-1-sulfonic acid (PFBS) "
                      "in cynomolgus monkeys and rats. Reprod Toxicol "
                      "2009;27(3-4):400-406. "
                      "doi:10.1016/j.reprotox.2009.01.013",
         human_clr="NONE -- no human renal clearance for PFBS in the project",
         admissible_to_eq5="no, blocked on BOTH sides",
         note="Same situation as PFBA on equation (5), but the Chengelis "
              "abstract carries a mechanistic claim worth testing against "
              "its own table: serum clearance of PFBS was faster in females "
              "(the project holds the sex pair) while urinary elimination "
              "showed NO sex difference. If the tables bear that out, the "
              "PFBS sex difference is NOT renal, which cuts against the "
              "OAT-reabsorption account that explains PFOA. Three rat TOTAL "
              "clearances are already held (Huang 2019 828, Chengelis 2009 "
              "946, Olsen 2009 2856 mL/kg/day), disagreeing 3.4-fold."),
    dict(chemical="PFDA", verdict="measured_but_unresolvable",
         rat_clr="measured by Ohmori 2003 Fig 2 but < 15 mL/kg/day, below "
                 "the figure's resolution in both sexes",
         measured_urine="yes (Ohmori 2003 renal clearance; Kudo 2001 "
                        "urinary fraction 0.2% of dose over 120 h)",
         blocking_pdf="none -- both papers already in hand",
         human_clr="YES: Fujii 2015 0.015, Zhang 2013 mean 0.096 / "
                   "median 0.035 mL/kg/day",
         admissible_to_eq5="yes, as a ONE-SIDED BOUND via the Kudo route",
         note="The only one of the four that can be advanced with data "
              "already held."),
    dict(chemical="PFOS", verdict="no_urine_study",
         rat_clr="none found",
         measured_urine="NO -- Huang 2019 (NTP) collected plasma, liver, "
                        "kidney and brain only; no urine at any dose",
         blocking_pdf="none identified -- no rat PFOS urine study located",
         human_clr="YES: Harada 2005 0.015 (n=20) / 0.012 (male), "
                   "Zhang 2013 mean 0.037 / median 0.024 mL/kg/day",
         admissible_to_eq5="yes, as a LOOSE one-sided bound via CL_R <= CL_tot",
         note="Harada 2005's abstract claim of rat renal-clearance sex "
              "differences was traced to its Table 5, which has rat rows "
              "for PFOA only. PFOS there is monkey and human."),
]

# ------------------------------------------------- cross-validation inputs
# Kudo 2001: fraction of IP dose recovered in urine over 120 h, Wistar.
KUDO_URINARY_FRACTION = {
    ("male", "PFHpA"): 0.92, ("male", "PFOA"): 0.55,
    ("male", "PFNA"): 0.020, ("male", "PFDA"): 0.002,
    ("female", "PFOA"): 0.80, ("female", "PFNA"): 0.51,
}
KUDO_WINDOW_DAYS = 5.0  # 120 h

# ------------------------------------------------------- rat total clearance
# Male rat PFOS total clearances already held by the project, mL/kg/day.
PFOS_RAT_CLTOT = [
    (9.24, "Kim 2016"),
    (11.3, "Chang 2012"),
    (13.1, "Huang 2019 (NTP), IV 2 mg/kg"),
]

# ----------------------------------------------------- human renal clearance
HUMAN_CLR = {
    "PFDA": [(0.015, "Fujii 2015 Table 2, urinary, n=10"),
             (0.035, "Zhang 2013 Table 2, median, n=60"),
             (0.096, "Zhang 2013 Table 2, mean, n=60")],
    "PFOS": [(0.012, "Harada 2005 Table 5, male 61 kg"),
             (0.015, "Harada 2005 Section 3.3, n=20"),
             (0.024, "Zhang 2013 Table 2, median, n=66"),
             (0.037, "Zhang 2013 Table 2, mean, n=66")],
}
NC = {"PFDA": 10, "PFOS": 8}

# ------------------------------------------------------- provenance flags
# Turned up while checking which of the rat clearance sources had citations.
# Recorded rather than acted on, because it concerns an input to equation (5)
# and should be settled from the primary paper, not by inference.
PROVENANCE_FLAGS = [
    dict(flag="rat PFHxS renal clearance 0.76 attributed to Kim 2016 is a "
              "SECONDARY citation",
         affects="scripts/stage17_renal_only_decomposition.py RAT_CLR, one of "
                 "only three rat inputs to equation (5); also "
                 "db/rat_tubular_factor.csv, db/tubular_factor_all_species."
                 "csv, db/tubular_factor_evidence.csv, "
                 "scripts/stage6_consolidate_tubular.py",
         evidence="db/vd_clearance.csv row 70 records the value with source "
                  "'ATSDR 2021 Table 3-6; EPA 2023 PFHxA IRIS sec 5.2.1' and "
                  "leaves pmid and doi blank -- it was never read from the "
                  "primary paper. Kim 2016 is now identified as Kim SJ et "
                  "al., Food Chem Toxicol 97:243-255, PMID 27637925, whose "
                  "abstract reports plasma and TISSUE distribution and does "
                  "not mention urine or renal clearance at all.",
         why_it_matters="If the 0.76 is not in Kim 2016, the PFHxS row of "
                        "equation (5) rests on a misattribution, and PFHxS "
                        "is one of three compounds there. The project "
                        "already flags a separate 'MAJOR INTERNAL "
                        "CONTRADICTION' in EPA Appendix B over a 'Kim et al. "
                        "(2016b)' whose half-lives differ 13-157 fold "
                        "between EPA's own narrative and tables, with a "
                        "suspected chemical misattribution. Two independent "
                        "provenance problems around the same author-year is "
                        "reason to check rather than assume.",
         resolution="Needs the Kim 2016 PDF (doi:10.1016/j.fct.2016.09.017) "
                    "to confirm or retract. NOT changed in this stage.",
         status="open"),
]


def auc_fraction_one_compartment(window_days, halflife_days):
    """Fraction of total AUC accumulated by `window_days`, one compartment.

    For any multi-exponential curve with the same TERMINAL half-life this is
    an UNDERESTIMATE, so dividing by it overstates CL_R.
    """
    return 1.0 - math.exp(-math.log(2.0) * window_days / halflife_days)


def read_ohmori():
    """Ohmori 2003 Table 1 totals and Fig. 2 renal clearances, by sex."""
    tot, ren = {}, {}
    p = HERE / "db" / "primary_2026" / "ohmori2003_table1_primary.csv"
    for r in csv.DictReader(p.open()):
        tot[(r["sex"], r["chemical"])] = (float(r["cltot_mL_kg_day"]),
                                          float(r["halflife_d"]))
    p = HERE / "db" / "primary_2026" / "ohmori2003_fig2_renal_clearance.csv"
    for r in csv.DictReader(p.open()):
        if r["gfr_source"] != "OEHHA/Han":
            continue
        key = (r["sex"], r["chemical"])
        if key in ren:
            continue
        ren[key] = (r["cl_renal_mL_kg_day"], r["value_kind"])
    return tot, ren


def fu_table():
    out = {}
    p = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(p.open()):
        if r["experiment"] == "Ryu 2024" and r["fu"]:
            sp = r["species"].split(" (")[0]
            if sp in ("human", "rat"):
                out[(sp, r["chemical"])] = float(r["fu"])
    return out


def cross_validate(tot, ren):
    """Run the Kudo arithmetic where Ohmori already gives the answer."""
    rows = []
    for (sex, chem), frac in sorted(KUDO_URINARY_FRACTION.items()):
        if (sex, chem) not in tot or (sex, chem) not in ren:
            continue
        cltot, t12 = tot[(sex, chem)]
        clr_text, kind = ren[(sex, chem)]
        af = auc_fraction_one_compartment(KUDO_WINDOW_DAYS, t12)
        predicted = frac * cltot / af
        row = dict(sex=sex, chemical=chem, kudo_urinary_fraction=frac,
                   cltot=cltot, halflife_d=t12, auc_fraction=round(af, 4),
                   predicted_clr=round(predicted, 3),
                   ohmori_clr=clr_text, ohmori_kind=kind)
        if kind.startswith("digitised value"):
            obs = float(clr_text)
            row["ratio_predicted_over_observed"] = round(predicted / obs, 3)
        else:
            row["ratio_predicted_over_observed"] = ""
        rows.append(row)
    return rows


def pfda_rat_bound(tot):
    cltot, t12 = tot[("male", "PFDA")]
    frac = KUDO_URINARY_FRACTION[("male", "PFDA")]
    af = auc_fraction_one_compartment(KUDO_WINDOW_DAYS, t12)
    return frac * cltot / af, cltot, t12, af, frac


def main():
    fu = fu_table()
    tot, ren = read_ohmori()
    bar = "=" * 78

    # ------------------------------------------------------------- 1. search
    print(bar)
    print("RAT RENAL CLEARANCE SEARCH: PFBA, PFBS, PFDA, PFOS")
    print(bar)
    print("  Directly measured rat renal clearance found for: NONE of the four.")
    print()
    for d in RAT_SEARCH:
        print(f"  {d['chemical']}  [{d['verdict']}]")
        print(f"    rat CL_R   : {d['rat_clr']}")
        print(f"    rat urine  : {d['measured_urine']}")
        print(f"    human CL_R : {d['human_clr']}")
        print(f"    into eq(5) : {d['admissible_to_eq5']}")
        print()
    print("  THE HUMAN SIDE BLOCKS PFBA AND PFBS REGARDLESS OF THE RAT.")
    print("  The Chang 2008, Olsen 2009 and Chengelis 2009 PDFs would")
    print("  complete the rat row for those two compounds -- and Chengelis")
    print("  would test its own sex-difference claim -- but none of them can")
    print("  admit PFBA or PFBS to equation (5) while the human side is empty.")
    print()

    # ---------------------------------------------------- 2. cross-validation
    print(bar)
    print("CROSS-VALIDATION OF THE KUDO ROUTE, WHERE OHMORI GIVES THE ANSWER")
    print(bar)
    print("  predicted CL_R = (urinary fraction) x CL_tot / (AUC fraction)")
    print()
    print(f"  {'sex':8}{'cmpd':7}{'frac':>7}{'CL_tot':>9}{'t_half':>8}"
          f"{'AUCfrac':>9}{'predicted':>11}{'Ohmori':>10}{'ratio':>8}")
    xv = cross_validate(tot, ren)
    ratios = []
    for r in xv:
        rr = r["ratio_predicted_over_observed"]
        if rr != "":
            ratios.append(rr)
        print(f"  {r['sex']:8}{r['chemical']:7}"
              f"{r['kudo_urinary_fraction']:>7.3f}{r['cltot']:>9.1f}"
              f"{r['halflife_d']:>8.2f}{r['auc_fraction']:>9.4f}"
              f"{r['predicted_clr']:>11.3f}{str(r['ohmori_clr']):>10}"
              f"{(f'{rr:.2f}' if rr != '' else '--'):>8}")
    print()
    print(f"  The route reads HIGH on all {len(ratios)} testable compounds, by "
          f"{min(ratios):.2f}x to {max(ratios):.2f}x.")
    print("  A one-directional overstatement is what an upper bound needs. It")
    print("  also denies the route any claim to be a point estimate.")
    print()
    print("  NOTE WHERE THE DISAGREEMENT CANNOT BE BLAMED ON THE MODEL. For")
    print("  PFHpA-male and PFOA-female the half-lives are 0.10 and 0.08 d,")
    print("  so the 120 h window captures the entire AUC and the correction")
    print("  is exactly 1. The 1.50x and 1.80x gaps there are straight")
    print("  inter-study disagreements between Kudo's urinary recovery and")
    print("  Ohmori's digitised bars -- different doses and routes (Kudo IP,")
    print("  Ohmori IV) in the same strain and lab. That is the floor on how")
    print("  well any derived renal clearance can be trusted here.")
    print()

    # -------------------------------------------------------- 3. rat bounds
    print(bar)
    print("UPPER BOUNDS ON MALE RAT RENAL CLEARANCE")
    print(bar)
    clr_pfda, cltot_pfda, t12_pfda, af_pfda, frac_pfda = pfda_rat_bound(tot)
    print(f"  PFDA   Kudo route: {frac_pfda} x {cltot_pfda} / {af_pfda:.4f}"
          f" = {clr_pfda:.4f} mL/kg/day")
    print(f"         (t_half {t12_pfda} d, window {KUDO_WINDOW_DAYS} d, so the"
          f" 120 h collection")
    print(f"          captured only {af_pfda*100:.1f}% of the total AUC)")
    f = fu[("rat", "PFDA")]
    print(f"         f_u {f} x GFR {GFR['rat']:.0f} = {f*GFR['rat']:.2f}"
          f"  ->  T_rat <= {clr_pfda/(f*GFR['rat']):.5f}")
    print(f"         Compare the digitisation's bound, T < 1.04: this is"
          f" {1.04/(clr_pfda/(f*GFR['rat'])):.0f}-fold tighter")
    print("         and, unlike it, lands decisively below unity.")
    print()
    f = fu[("rat", "PFOS")]
    print(f"  PFOS   CL_R <= CL_tot, f_u {f} x GFR {GFR['rat']:.0f}"
          f" = {f*GFR['rat']:.2f}")
    for v, src in PFOS_RAT_CLTOT:
        print(f"         CL_tot {v:>6.2f}  ->  T_rat <= "
              f"{v/(f*GFR['rat']):.4f}   {src}")
    print("         Loose: the rat excretes PFOS substantially in faeces, so")
    print("         the true renal clearance is well below these bounds.")
    print()

    # ------------------------------------------- 4. bounded decomposition
    print(bar)
    print("EQUATION (5) EXTENDED AS ONE-SIDED BOUNDS  (rat / human, male)")
    print(bar)
    print("  CL_R,rat / CL_R,human  =  binding x glomerular x tubule")
    print("  Every rat input is an upper bound, so the ratio and the tubule")
    print("  factor are upper bounds too. Binding and glomerular are exact.")
    print()
    rat_bounds = {
        "PFDA": [(clr_pfda, "Kudo 2001 urinary fraction / Ohmori 2003 CL_tot,"
                            " AUC-corrected, UPPER BOUND")],
        "PFOS": [(v, f"CL_R <= CL_tot, {src}, UPPER BOUND")
                 for v, src in PFOS_RAT_CLTOT],
    }
    print(f"  {'cmpd':6}{'CL_R ratio':>12}{'binding':>9}{'glom':>7}"
          f"{'tubule':>10}  {'bind%':>6}{'glom%':>7}{'tub%':>7}  stable?")
    rows = []
    for chem in ("PFDA", "PFOS"):
        fh, fr = fu[("human", chem)], fu[("rat", chem)]
        fb, fg = fr / fh, GFR["rat"] / GFR["human"]
        for vh, sh in HUMAN_CLR[chem]:
            th = vh / (fh * GFR["human"])
            for vr, sr in rat_bounds[chem]:
                tr = vr / (fr * GFR["rat"])
                ratio, ft = vr / vh, tr / th
                lg = math.log10(ratio)
                # Shares are meaningless when the three log terms nearly
                # cancel. The symptom is a share exceeding the whole, so
                # test that directly rather than guessing a ratio threshold.
                pb = math.log10(fb) / lg * 100
                pg = math.log10(fg) / lg * 100
                pt = math.log10(ft) / lg * 100
                stable = max(abs(pb), abs(pg), abs(pt)) <= 100.0
                if stable:
                    shares = f"  {pb:>6.0f}{pg:>7.0f}{pt:>7.0f}  yes"
                else:
                    pb = pg = pt = ""
                    shares = f"  {'--':>6}{'--':>7}{'--':>7}  NO"
                print(f"  {chem:6}{ratio:>12.1f}{fb:>9.2f}{fg:>7.2f}"
                      f"{ft:>10.3f}{shares}")
                rows.append(dict(
                    chemical=chem, n_carbons=NC[chem],
                    cl_renal_human=vh, human_source=sh,
                    cl_renal_rat_bound=round(vr, 4), rat_source=sr,
                    cl_renal_ratio_bound=round(ratio, 2),
                    T_human=round(th, 6), T_rat_bound=round(tr, 6),
                    factor_binding=round(fb, 4),
                    factor_glomerular=round(fg, 4),
                    factor_tubule_bound=round(ft, 4),
                    pct_binding=(round(pb, 1) if stable else ""),
                    pct_glomerular=(round(pg, 1) if stable else ""),
                    pct_tubule=(round(pt, 1) if stable else ""),
                    shares_stable=("yes" if stable else "no"),
                    bound_direction="upper"))
    print()
    print("  'stable?' is no when any single share exceeds 100% of the whole.")
    print("  That is the symptom of the three log terms nearly cancelling: the")
    print("  percentages explode while the underlying factors stay perfectly")
    print("  well determined -- the same catastrophic-cancellation trap as the")
    print("  OEHHA difference form. Read the FACTORS in those rows, never the")
    print("  shares.")
    print()
    print(bar)
    print("WHAT THE TWO NEW COMPOUNDS SAY")
    print(bar)
    pfda = [r for r in rows if r["chemical"] == "PFDA"]
    pfos = [r for r in rows if r["chemical"] == "PFOS"]
    tb_pfda = [r["factor_tubule_bound"] for r in pfda]
    tb_pfos = [r["factor_tubule_bound"] for r in pfos]
    print("  PFDA breaks the pattern of the three compounds in equation (5).")
    print(f"  Its tubule factor is bounded above by {max(tb_pfda):.2f}, and")
    print(f"  below 1 against ALL {len(pfda)} human sources ({min(tb_pfda):.2f}"
          f" to {max(tb_pfda):.2f}): the rat")
    print("  tubule is NOT more secretory than the human's for PFDA. The")
    print("  renal clearance difference is carried by glomerular filtration")
    print("  (5.6x, a species constant) and binding (1.6x) alone. For PFHpA,")
    print("  PFOA and PFHxS the tubule term dominated instead.")
    print()
    sh_pfos = [r["pct_tubule"] for r in pfos if r["pct_tubule"] != ""]
    print("  PFOS keeps the familiar pattern. Its tubule factor is bounded")
    print(f"  above by {max(tb_pfos):.0f} and exceeds {min(tb_pfos):.0f} on"
          f" every one of the {len(pfos)} source")
    print(f"  combinations, carrying {min(sh_pfos):.0f}-{max(sh_pfos):.0f}% of"
          f" the log ratio, so tubular")
    print("  handling dominates there even at the bound.")
    print()

    # ------------------------------------------------- 5. provenance flags
    print(bar)
    print("PROVENANCE FLAG RAISED WHILE SEARCHING  (not acted on)")
    print(bar)
    for d in PROVENANCE_FLAGS:
        print(f"  [{d['status'].upper()}] {d['flag']}")
        print(f"    affects   : {d['affects']}")
        print(f"    evidence  : {d['evidence']}")
        print(f"    matters   : {d['why_it_matters']}")
        print(f"    to resolve: {d['resolution']}")
        print()

    # -------------------------------------------------------------- 6. write
    out = []
    p1 = HERE / "db" / "rat_renal_clearance_search.csv"
    cols = ["chemical", "verdict", "rat_clr", "measured_urine",
            "human_clr", "admissible_to_eq5", "blocking_pdf", "note"]
    with p1.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in sorted(RAT_SEARCH, key=lambda x: x["chemical"]):
            w.writerow({k: d[k] for k in cols})
    out.append((p1, len(RAT_SEARCH)))

    p2 = HERE / "db" / "kudo_ohmori_crossvalidation.csv"
    cols = ["sex", "chemical", "kudo_urinary_fraction", "cltot",
            "halflife_d", "auc_fraction", "predicted_clr", "ohmori_clr",
            "ohmori_kind", "ratio_predicted_over_observed"]
    with p2.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(xv, key=lambda x: (x["chemical"], x["sex"])):
            w.writerow(r)
    out.append((p2, len(xv)))

    p3 = HERE / "db" / "renal_only_bounds.csv"
    with p3.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["chemical"],
                                             x["cl_renal_human"],
                                             x["cl_renal_rat_bound"])):
            w.writerow(r)
    out.append((p3, len(rows)))

    p4 = HERE / "db" / "renal_clearance_provenance_flags.csv"
    cols = ["flag", "affects", "evidence", "why_it_matters", "resolution",
            "status"]
    with p4.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in sorted(PROVENANCE_FLAGS, key=lambda x: x["flag"]):
            w.writerow({k: d[k] for k in cols})
    out.append((p4, len(PROVENANCE_FLAGS)))

    for p, n in out:
        print(f"wrote {p.relative_to(HERE)}  ({n} rows)")


if __name__ == "__main__":
    main()
