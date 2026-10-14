"""
Argoul 2026 (eleven PFAS as one cocktail) against the single-compound studies.

Argoul carries a lot of weight in this review: it is the only source that
measures clearance and volume for many PFAS in one experiment, which is what
sections 3 and 10 need. Every other mouse value in the database comes from a
study that dosed ONE compound. So the question is whether dosing eleven at
once, at low dose, by mixed-effects fitting, reproduces what the
one-at-a-time literature found.

The answer has a shape. Clearance and volume are both systematically LOW --
and their ratio is not. That matters, because a half-life is their ratio.

All comparators are FEMALE mouse, matching Argoul; db/animal_halflife_measured.csv
and db/combined/tk_parameters.csv. Units harmonised to days and mL/kg.

Run:  python3 scripts/argoul_vs_single_compound.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LN2 = math.log(2)
OUT = os.path.join(HERE, "db", "argoul_vs_single_compound.csv")
PAIRS = []

# ---------------------------------------------------------------------------
# Provenance. Every comparator is named by a key here, so a value in the output
# CSV can be traced to a paper, a DOI and the table it was actually read from.
# "via" matters: several values reach this review through an agency compilation
# rather than the primary paper, and that is recorded rather than smoothed over.
# ---------------------------------------------------------------------------
REFS = {
 "lou2009": dict(
   ref="Lou I, Wambaugh JF, Lau C, Hanson RG, Lindstrom AB, Strynar MJ, Zehr RD, "
       "Setzer RW, Barton HA (2009). Modeling single and repeated dose "
       "pharmacokinetics of PFOA in mice. Toxicological Sciences 107(2):331-341.",
   doi="10.1093/toxsci/kfn234", pmid="19005225",
   via="primary full text; Lou 2009 Table 2 via db/primary_2026/lou2009_mouse_tk.csv"),
 "fujii2015": dict(
   ref="Fujii Y, Niisoe T, Harada KH, Uemoto S, Ogura Y, Takenaka K, Koizumi A "
       "(2015). Toxicokinetics of perfluoroalkyl carboxylic acids with different "
       "carbon chain lengths in mice and humans. Journal of Occupational Health "
       "57(1):1-12.",
   doi="10.1539/joh.14-0136-OA", pmid="25422127",
   via="OEHHA 2024 PHG Table 4.8.1; ATSDR 2021 Table 3-6 (compilations)"),
 "sundstrom2012": dict(
   ref="Sundstrom M, Chang S-C, Noker PE, Gorman GS, Hart JA, Ehresman DJ, "
       "Bergman A, Butenhoff JL (2012). Comparative pharmacokinetics of "
       "perfluorohexanesulfonate (PFHxS) in rats, mice, and monkeys. "
       "Reproductive Toxicology 33(4):441-451.",
   doi="10.1016/j.reprotox.2011.07.004", pmid="21856411",
   via="primary full text; Sundstrom 2012 Table 3 via "
       "db/primary_2026/sundstrom2012_pfhxs_three_species.csv"),
 "chang2008": dict(
   ref="Chang S-C, Das K, Ehresman DJ, Ellefson ME, Gorman GS, Hart JA, Noker PE, "
       "Tan Y-M, Lieder PH, Lau C, Olsen GW, Butenhoff JL (2008). Comparative "
       "pharmacokinetics of perfluorobutyrate in rats, mice, monkeys, and humans "
       "and relevance to human exposure via drinking water. Toxicological "
       "Sciences 104(1):40-53.",
   doi="10.1093/toxsci/kfn057", pmid="18353799",
   via="ATSDR 2021 Toxicological Profile for Perfluoroalkyls, Table 3-5 "
       "(compilation); Chang 2008 abstract"),
 "chang2012": dict(
   ref="Chang S-C, Noker PE, Gorman GS, Gibson SJ, Hart JA, Ehresman DJ, "
       "Butenhoff JL (2012). Comparative pharmacokinetics of "
       "perfluorooctanesulfonate (PFOS) in rats, mice, and monkeys. "
       "Reproductive Toxicology 33(4):428-440.",
   doi="10.1016/j.reprotox.2011.07.002", pmid="21889587",
   via="ATSDR 2021 Toxicological Profile for Perfluoroalkyls, Table 3-5 "
       "(compilation)"),
 "tatum2011": dict(
   ref="Tatum-Gibbs K, Wambaugh JF, Das KP, Zehr RD, Strynar MJ, Lindstrom AB, "
       "Delinsky A, Lau C (2011). Comparative pharmacokinetics of "
       "perfluorononanoic acid in rat and mouse. Toxicology 281(1-3):48-55.",
   doi="10.1016/j.tox.2011.01.003", pmid="21237237",
   via="Tatum-Gibbs 2011 abstract; ATSDR 2021 Table 3-5 (compilation)"),
 "epa2023pfhxa": dict(
   ref="US EPA (2023). IRIS Toxicological Review of Perfluorohexanoic Acid "
       "(PFHxA) and Related Salts. EPA/635/R-23/027Fa.",
   doi="", pmid="",
   via="EPA 2023 PFHxA IRIS sec 5.2.1 p.5-13; Vd,beta COMPUTED by EPA from "
       "Gannon et al. 2011, not measured"),
 "zurlinden2025": dict(
   ref="Zurlinden TJ, Dzierlenga MW, Kapraun DF, Ring C, Bernstein AS, "
       "Schlosser PM, Morozov V (2025). Estimation of species- and sex-specific "
       "PFAS pharmacokinetics in mice, rats, and non-human primates using a "
       "Bayesian hierarchical methodology. Toxicology and Applied Pharmacology "
       "499:117336.",
   doi="10.1016/j.taap.2025.117336", pmid="40210099",
   via="EPA CPHEA-Animal-PFAS-PK Bayesian hierarchical fit, Table 3; pools the "
       "single-compound studies above rather than adding new animals"),
 "argoul2026": dict(
   ref="Argoul CML, Toutain P-L, Picard-Hagen N, Mselli-Lakhal L, Dauwe Y, "
       "Roques BB, Lacroix MZ, Gayrard V (2026). Nonlinear mixed-effects "
       "modeling of the intravenous and oral kinetics of eleven perfluoroalkyl "
       "substances in female mice. Environmental Research 303:124802.",
   doi="10.1016/j.envres.2026.124802", pmid="",
   via="primary full text; Argoul 2026 Tables 1/2/4 via "
       "db/primary_2026/argoul2026_mouse_tk.csv"),
}

# Female-mouse comparators from single-compound studies, half-life in DAYS.
# Chang and Sundstrom report hours for serum; converted here.
# Each entry: (label, half-life days, reference key, dose string, lowest dose
# in mg/kg for the gap arithmetic, route). Dose is recorded because the dose
# mismatch against the cocktail is the first thing anyone asks about -- see
# scripts/argoul_dose_matching.py, which tests it and finds it insufficient.
SINGLE = {
    "PFBA":  [("Chang 2008, 10 mg/kg",  2.87 / 24, "chang2008", "10", 10.0, "oral"),
              ("Chang 2008, 30 mg/kg",  3.08 / 24, "chang2008", "30", 30.0, "oral"),
              ("Chang 2008, 100 mg/kg", 2.79 / 24, "chang2008", "100", 100.0, "oral")],
    "PFOA":  [("Lou 2009, 1-10 mg/kg", 15.6, "lou2009", "1, 10", 1.0, "oral")],
    "PFNA":  [("Tatum-Gibbs 2011, 1-10 mg/kg low", 25.8, "tatum2011",
               "1, 10", 1.0, "oral gavage"),
              ("Tatum-Gibbs 2011, 1-10 mg/kg high", 68.4, "tatum2011",
               "1, 10", 1.0, "oral gavage")],
    "PFHxS": [("Sundstrom 2012, 1 mg/kg", 24.8, "sundstrom2012", "1", 1.0, "oral"),
              ("Sundstrom 2012, 20 mg/kg", 26.8, "sundstrom2012", "20", 20.0, "oral")],
    "PFOS":  [("Chang 2012, 1 mg/kg", 907 / 24, "chang2012", "1", 1.0, "oral"),
              ("Chang 2012, 20 mg/kg", 731 / 24, "chang2012", "20", 20.0, "oral")],
    "PFBS":  [("Lau 2020, 30 mg/kg", 3.0 / 24, "", "30", 30.0, "oral"),
              ("Lau 2020, 300 mg/kg", 3.9 / 24, "", "300", 300.0, "oral")],
}
# The EPA Bayesian hierarchical fit pools the single-compound studies, so it is
# a second summary of the same literature rather than an independent study.
POOLED = {"PFBA": 0.186, "PFBS": 0.125, "PFHxA": 0.251, "PFHxS": 26.665,
          "PFNA": 59.439, "PFOA": 21.453, "PFOS": 32.810}

# Direct CL / Vd comparators, female mouse (mL/kg/d, mL/kg).
#
# The last field says whether that study's Vd is INDEPENDENT of its clearance.
# It usually is not: non-compartmental analysis gives Vz = CL/lambda_z and
# Vss = CL*MRT, and Lou's clearance is ke*Vd derived in this review, so in those
# studies "clearance low AND volume low" is one observation seen twice. Fujii is
# the exception -- Vd from Dose/C(0) by back-extrapolation, CL from Dose/AUC --
# and is therefore the only row where the two agree independently.
CLVD = {
    "PFOA":  [("Lou 2009", 5.99, 135, "lou2009", False, "1, 10", 1.0, "oral"),
              # 0.31 umol/kg x 414.07 g/mol = 0.1284 mg/kg -- near dose-matched
              # to the cocktail's 0.080, and the WORST-agreeing comparison.
              ("Fujii 2015", 11.8, 150, "fujii2015", True, "0.128", 0.1284, "IV")],
    "PFHxS": [("Sundstrom 2012, 1 mg/kg", 2.68, 96, "sundstrom2012", False,
               "1", 1.0, "oral"),
              ("Sundstrom 2012, 20 mg/kg", 3.79, 147, "sundstrom2012", False,
               "20", 20.0, "oral")],
    # EPA computed this Vd,beta from Gannon 2011, whose mouse arm ran 100 mg/kg.
    "PFHxA": [("US EPA PFHxA IRIS", None, 780, "epa2023pfhxa", False,
               "100 (Gannon 2011)", 100.0, "oral gavage")],
}


def check_refs_against_pdf():
    """Warn if a comparator's DOI here disagrees with the review's reference
    list, or is missing from it. A wrong citation in one file and not the other
    is exactly the error this comparison introduced once already."""
    import re
    try:
        pdf = open(os.path.join(HERE, "scripts", "make_summary_pdf.py")).read()
    except OSError:
        return
    listed = dict(re.findall(r'^"([a-z0-9]+)":\s*"(.*?)",\n', pdf, re.M | re.S))
    bad = []
    for k, r in REFS.items():
        if k not in listed:
            bad.append(f"{k}: not in the review's reference list")
        elif r["doi"] and r["doi"] not in listed[k]:
            bad.append(f"{k}: doi {r['doi']} not in the listed reference")
    if bad:
        print("\n  WARNING, references disagree with make_summary_pdf.py:")
        for b in bad:
            print(f"    {b}")
    else:
        print(f"  references: {len(REFS)} keys agree with the review's list")


def dose(a_row, dose_str="", dose_low=None, route=""):
    """Dose fields for one comparison. Argoul's IV and oral doses differ and the
    fit uses both, so both are carried; the gap is computed against the IV dose
    because every clearance here is an IV or IV-anchored estimate."""
    iv = float(a_row["iv_dose_mgkg"]) if a_row["iv_dose_mgkg"] else None
    orl = float(a_row["oral_dose_mgkg"]) if a_row["oral_dose_mgkg"] else None
    gap = (dose_low / iv) if (iv and dose_low) else ""
    return dict(argoul_dose_iv_mgkg=iv if iv else "",
                argoul_dose_oral_mgkg=orl if orl else "",
                single_dose_mgkg=dose_str, single_dose_route=route,
                dose_gap_single_over_argoul=round(gap, 1) if gap else "")


def cite(key):
    """Reference fields for one comparator, plus the Argoul side, so a row of
    the output CSV is traceable without the rest of the repository."""
    r = REFS.get(key, {})
    a = REFS["argoul2026"]
    return dict(comparator_ref=r.get("ref", ""), comparator_doi=r.get("doi", ""),
                comparator_pmid=r.get("pmid", ""), comparator_via=r.get("via", ""),
                argoul_ref=a["ref"], argoul_doi=a["doi"])


def argoul():
    p = os.path.join(HERE, "db", "primary_2026", "argoul2026_mouse_tk.csv")
    return {r["chemical"]: r for r in csv.DictReader(open(p))}


def main():
    a = argoul()
    print(__doc__.split("Run:")[0].rstrip())

    print("=" * 78)
    print("PART 1.  Clearance and volume, where a single-compound study reports both")
    print("=" * 78)
    print(f"\n  {'cmpd':<7s} {'source':<26s} {'CL':>8s} {'Vd':>7s} "
          f"{'CL ratio':>9s} {'Vd ratio':>9s}")
    cl_r, vd_r = [], []
    for c, comps in CLVD.items():
        r = a[c]
        acl = float(r["cl_mL_kg_day"]) if r["cl_mL_kg_day"] else None
        avd = float(r["vss_L_kg"]) * 1000
        print(f"  {c:<7s} {'Argoul 2026 (cocktail)':<26s} "
              f"{(f'{acl:g}' if acl else '-'):>8s} {avd:>7g} "
              f"{'-':>9s} {'-':>9s}")
        for name, ccl, cvd, key, indep, dstr, dlow, route in comps:
            dz = dose(r, dstr, dlow, route)
            rc = f"{acl/ccl:.2f}x" if (acl and ccl) else "-"
            rv = f"{avd/cvd:.2f}x" if cvd else "-"
            if acl and ccl: cl_r.append(acl / ccl)
            if cvd: vd_r.append(avd / cvd)
            PAIRS.append(dict(
                parameter="clearance" if (acl and ccl) else "", chemical=c,
                units="mL/kg/d", comparator=name, argoul=acl, single=ccl,
                ratio=(acl / ccl) if (acl and ccl) else "",
                vd_independent_of_cl="", **cite(key), **dz))
            PAIRS.append(dict(
                parameter="volume_of_distribution" if cvd else "", chemical=c,
                units="mL/kg", comparator=name, argoul=avd, single=cvd,
                ratio=(avd / cvd) if cvd else "",
                vd_independent_of_cl="yes" if indep else "no",
                **cite(key), **dz))
            print(f"  {'':<7s} {name:<26s} {(f'{ccl:g}' if ccl else '-'):>8s} "
                  f"{cvd:>7g} {rc:>9s} {rv:>9s}")
    gm = lambda v: math.exp(sum(math.log(x) for x in v) / len(v))
    vd_no_hxa = [v for v in vd_r if v < 3]
    print(f"""
  PFHxA's volume is the one outlier (5.1x HIGH) and is held out of the mean
  below; it is also the only compound Argoul puts in net tubular secretion,
  so its kinetics differ in kind.

  Geometric mean ratio, Argoul / single-compound:
      clearance                {gm(cl_r):.2f}x   ({1/gm(cl_r):.1f}x low)
      volume, excluding PFHxA  {gm(vd_no_hxa):.2f}x   ({1/gm(vd_no_hxa):.1f}x low)

  Both low, by similar factors. Hold on to that -- it is the whole result.

  One caveat on how much of a check that is. In most of these studies Vd is
  not measured apart from clearance: non-compartmental analysis gives
  Vz = CL/lambda_z and Vss = CL*MRT, and Lou's clearance is ke*Vd derived in
  this review. So for Lou, Sundstrom and the EPA PFHxA value, "clearance low
  AND volume low" is largely one observation seen twice. Fujii 2015 is the
  exception: Vd from Dose/C(0) by back-extrapolation, CL from Dose/AUC, which
  are independent. There CL is 0.38x and Vd 0.59x -- both low, independently.
  One row is thin evidence, so the defensible statement is the narrower one:
  the RATE CONSTANT agrees across ten comparisons and the SCALE does not.
""")

    print("=" * 78)
    print("PART 2.  Half-life: the ratio of the two, so the offsets should cancel")
    print("=" * 78)
    print("""
  Argoul reports MRT, not a terminal half-life. MRT = Vss/CL exactly in its
  table (checked: agrees within 6% for 9 of 10 compounds), so ln2 x MRT is the
  half-life a one-compartment system with that Vss and CL would show. Against a
  measured TERMINAL half-life that is an approximation, and for a compound with
  a deep second compartment it is a LOWER bound. Read the ratios as indicative.
""")
    print(f"  {'cmpd':<7s} {'Argoul ln2*MRT':>14s} {'single-compound':>16s} "
          f"{'ratio':>7s}   source")
    ratios, pooled_r = [], []
    for c in ["PFBA", "PFBS", "PFHxA", "PFOA", "PFHxS", "PFOS", "PFNA"]:
        r = a.get(c)
        if not r or not r["mrt_d"]:
            continue
        at = LN2 * float(r["mrt_d"])
        for name, t, key, dstr, dlow, route in SINGLE.get(c, []):
            ratios.append(at / t)
            print(f"  {c:<7s} {at:>14.3f} {t:>16.3f} {at/t:>6.2f}x   {name}")
            PAIRS.append(dict(parameter="half_life", chemical=c, units="d",
                              comparator=name, argoul=at, single=t,
                              ratio=at / t, vd_independent_of_cl="",
                              **cite(key), **dose(r, dstr, dlow, route)))
        if c in POOLED:
            # Kept for the note below, NOT written to the CSV: the pooled fit
            # re-analyses the very studies above, so a row of it next to them
            # would count the same animals twice.
            pooled_r.append(at / POOLED[c])
    both = ratios + pooled_r
    print(f"\n  {len(ratios)} comparisons, every one against a study that dosed that "
          f"compound\n  alone. Ratios span {min(ratios):.2f}x to {max(ratios):.2f}x, "
          f"geometric mean {gm(ratios):.2f}x.")
    print(f"""
  The EPA pooled fit [zurlinden2025] is NOT in the table above, and not in the
  output CSV. It is a Bayesian hierarchical re-analysis of the same single-
  compound studies, so putting it beside them counts the same animals twice,
  and it is a model output rather than a measurement -- for PFBA it returns
  {POOLED['PFBA']:.3f} d where every Chang arm sits at 0.116-0.128. For the record,
  its {len(pooled_r)} ratios would move the geometric mean to {gm(both):.2f}x
  (span {min(both):.2f}x to {max(both):.2f}x), which does not change the conclusion.

  One thing that costs: PFHxA has no single-compound female-mouse comparator,
  so dropping the pooled rows drops PFHxA from the half-life comparison
  entirely. It survives in the volume comparison only.

  PFBS and PFDS have no comparator either; Argoul reports no clearance for
  PFBS at all.""")
    print("""
  So the half-lives broadly agree while the clearances and volumes do not.
  That is the signature of a shared scaling on both terms rather than a
  disagreement about elimination: in t = ln2*Vd/CL, a factor that multiplies
  Vd and CL together cancels.""")

    print()
    print("=" * 78)
    print("PART 3.  What that leaves")
    print("=" * 78)
    print("""
  Three candidate causes, none separable with the data in hand:

   DOSE      The obvious candidate, and it does not survive testing. See
             scripts/argoul_dose_matching.py. The dose is tuned per compound
             rather than swept, and matched that way it fails to order the
             discrepancies: Fujii is nearly dose-matched (0.128 vs 0.080
             mg/kg) and disagrees MOST, Lou has a >=12x dose gap and agrees
             BEST. The +0.20 slope once quoted here was fitted through
             Argoul's own point, which carries 79% leverage on it; dropping
             it and extrapolating Sundstrom's two doses leaves ~1.6x
             unexplained. A contributor, not the explanation.

   COCKTAIL  Eleven PFAS dosed together could compete for the same
             reabsorptive transporters. That would RAISE clearance, not lower
             it, so it argues against the observed direction -- unless the
             competition is for plasma binding sites instead, which would
             raise the unbound fraction and lower the apparent volume.

   METHOD    Nonlinear mixed-effects fitting shrinks estimates toward the
             population mean; non-compartmental analysis does not. Shrinkage
             would compress the spread rather than shift the centre, so it is
             the weakest of the three.

  The useful conclusion for this review: Argoul's WITHIN-experiment contrasts
  (section 3's 5,254x clearance span against 7.7x in Vss) are unaffected,
  because a shared scaling cancels from a ratio. Its ABSOLUTE clearances are
  what section 10's endpoint R is built from, and those carry a roughly 2-3x
  low bias relative to the single-compound literature.""")

    fields = ["parameter", "chemical", "units", "argoul", "single", "ratio",
              "argoul_dose_iv_mgkg", "argoul_dose_oral_mgkg",
              "single_dose_mgkg", "single_dose_route",
              "dose_gap_single_over_argoul",
              "comparator", "comparator_ref", "comparator_doi",
              "comparator_pmid", "comparator_via", "vd_independent_of_cl",
              "argoul_ref", "argoul_doi"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in PAIRS:
            if row["parameter"]:
                w.writerow(row)
    print(f"\n  wrote {os.path.relpath(OUT, HERE)} "
          f"({sum(1 for r in PAIRS if r['parameter'])} comparison rows)")
    check_refs_against_pdf()


if __name__ == "__main__":
    main()
