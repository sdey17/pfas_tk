#!/usr/bin/env python3
"""The rat tubule, from rat studies only (Stage 4d).

    python3 -I scripts/stage4d_rat_within_species.py

Writes db/rat_tubular_factor.csv. Idempotent.

THE RULE THIS FILE OBEYS. Cross-study joins are allowed, within a species.
Every clearance here was measured in a rat; every fu here was measured in rat
plasma. No value is carried across species, in either direction.

TWO CONVENTIONS, STATED ONCE.

  1. THE TUBULAR FACTOR, PER SPECIES, NOT AS A RATIO.

         T  =  CL_R / (fu x GFR)

     T is the observed renal clearance as a fraction of the filtration
     clearance, and it has a hard reading:

         T  <  1    net tubular reabsorption, of (1 - T) of the filtered load
         T  =  1    filtration alone
         T  >  1    net tubular secretion

     Reporting T per species rather than only its between-species ratio is
     the change this file makes. A ratio of 45x hides whether both species
     reabsorb or one secretes, and those are different mechanisms.

  2. WHY THERE IS NO THIRD TERM. Writing FR for fractional reabsorption and
     S for secretory clearance,

         CL_R  =  fu x GFR x (1 - FR)  +  S        so   T  =  1 - FR + S/(fu x GFR)

     and only the COMBINATION is identifiable from in vivo data. Splitting it
     needs an independent measurement of FR or S in the same species, and
     scaled in vitro CL_int is not that measurement: Stage 3 showed the
     expression gate and the REF/RAF choice move it about five-fold, and the
     albumin effect is unquantified for PFAS. Fitting a three-term model on it
     would reintroduce exactly the failure of fu = 0.02 -- an unmeasured input
     treated as known. So T stays lumped, and the in vitro work constrains it
     rather than decomposing it.

TWO TIERS OF EVIDENCE, AND THE SECOND NEEDS NO ASSUMPTION AT ALL.

  DIRECT   a rat study reports renal clearance. T is computed outright.
  BOUNDED  only total clearance is reported. Since CL_R <= CL_total always,
           T <= CL_total / (fu x GFR). WHERE THAT UPPER BOUND IS BELOW ONE,
           NET REABSORPTION IS PROVEN, and the reabsorbed fraction is at
           least 1 - T_upper. No renal-fraction assumption is needed, which
           is why this tier is worth having. Where the bound exceeds one it
           says nothing, and is reported as uninformative rather than
           dropped.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
GFR_RAT = 14400.0          # mL/kg/day, OEHHA PHG Table A6.4 / Han 2012

# ---- rat renal clearance, mL/kg/day, from rat studies -------------------
# (chemical, sex, CL_R, source, note)
DIRECT = [
    ("PFOA", "male", 18.2,
     "Kemper 2003 (unpublished), as used by OEHHA PHG Table A6.4", ""),
    ("PFOA", "male", 21.2, "Ohmori 2003",
     "digitized from a graph by Han 2012"),
    ("PFOA", "male", 0.032 * 1440,
     "Kudo 2002 Table 3", "0.032 mL/min/kg, 5 h after IV; same table's "
     "castrated male 0.447 and probenecid 0.010 mL/min/kg"),
    ("PFOA", "female", 666.0,
     "Kemper 2003 (unpublished), as used by OEHHA PHG Table A6.4", ""),
    ("PFOA", "female", 1009.0, "Ohmori 2003", "digitized from a graph"),
    ("PFOA", "female", 0.732 * 1440, "Kudo 2002 Table 3",
     "0.732 mL/min/kg, 5 h after IV"),
    ("PFHxS", "male", 0.76, "Kim 2016",
     "reported alongside total CL 9.01, so renal is 8% of total"),
]

# ---- rat TOTAL clearance, mL/kg/day, from rat studies -------------------
TOTAL = [
    ("PFOA",  "male",   47.39, "Kim 2016", "IV"),
    ("PFOA",  "male",   50.4,  "Kudo 2002", "IV"),
    ("PFOA",  "male",  135.0,  "Ohmori 2003 via ATSDR 2021 Table 3-6",
     "inconsistent with its own Vd 196 and t1/2 5.63 d, which give 24"),
    ("PFHpA", "male", 1604.0,  "Ohmori 2003 via ATSDR 2021 Table 3-6", ""),
    ("PFNA",  "male",    6.9,  "Ohmori 2003 via ATSDR 2021 Table 3-6", ""),
    ("PFDA",  "male",  207.0,  "Ohmori 2003 via ATSDR 2021 Table 3-6",
     "ATSDR prints SD 0.054 for this, implausible; treat with caution"),
    ("PFOS",  "male",    9.24, "Kim 2016", "IV"),
    ("PFOS",  "male",   13.1,  "Huang 2019 (NTP)", "beta phase"),
    ("PFOS",  "male",    1.63, "Huang 2019", ""),
    ("PFOS",  "male",   11.3,  "Chang 2012", "2 mg/kg; falls to 4.9 at 15"),
    ("PFHxS", "male",    6.7,  "Sundstrom 2012", "IV, 10 mg/kg, 10 weeks"),
    ("PFHxS", "male",    9.01, "Kim 2016", "IV"),
    ("PFHxS", "male",    4.73, "Huang 2019 (NTP)", "beta phase"),
    ("PFBS",  "male",  828.0,  "Huang 2019 (NTP)", "beta phase"),
    ("PFBS",  "male",  946.0,  "Chengelis 2009", ""),
    ("PFBS",  "male", 2856.0,  "Olsen 2009", ""),
    ("PFBA",  "male",  851.0,  "Chang 2008", "IV; 494 oral"),
    ("PFHxA", "male", 2784.0,  "Chengelis 2009", ""),
]

# ---- rat plasma fu measured IN A RAT STUDY ------------------------------
# Kim 2016 measured it by ultrafiltration, by sex. Ryu 2024 by presaturation
# equilibrium dialysis, sex-pooled. Both are rat values.
FU_RAT_KIM = {("PFHxS", "male"): 0.00076, ("PFHxS", "female"): 0.00069}


def ryu_rat_fu():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] == "Ryu 2024" and r["species"] == "rat" and r["fu"]:
            out[r["chemical"]] = float(r["fu"])
    return out


def main():
    fu = ryu_rat_fu()
    rows = []

    for chem, sex, clr, src, note in DIRECT:
        f = fu.get(chem)
        if f is None:
            continue
        filt = f * GFR_RAT
        T = clr / filt
        rows.append(dict(chemical=chem, sex=sex, tier="DIRECT", fu=f,
                         filt=filt, cl=clr, T=T, src=src, note=note))

    for chem, sex, clt, src, note in TOTAL:
        f = fu.get(chem)
        if f is None:
            continue
        filt = f * GFR_RAT
        T = clt / filt
        rows.append(dict(chemical=chem, sex=sex, tier="BOUNDED", fu=f,
                         filt=filt, cl=clt, T=T, src=src, note=note))

    p = HERE / "db" / "rat_tubular_factor.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "sex", "tier", "fu_rat", "fu_source",
                    "gfr_mL_kg_day", "filtered_load_mL_kg_day",
                    "clearance_mL_kg_day", "clearance_kind", "T",
                    "verdict", "min_pct_reabsorbed", "clearance_source",
                    "note"])
        for r in rows:
            direct = r["tier"] == "DIRECT"
            if direct:
                verdict = ("net_reabsorption" if r["T"] < 1
                           else "NET_TUBULAR_SECRETION")
                mn = round(100 * (1 - r["T"]), 1) if r["T"] < 1 else ""
            elif r["T"] < 1:
                verdict = "net_reabsorption_PROVEN_by_upper_bound"
                mn = round(100 * (1 - r["T"]), 1)
            else:
                verdict = "UNINFORMATIVE_bound_exceeds_filtration"
                mn = ""
            w.writerow([r["chemical"], r["sex"], r["tier"], r["fu"],
                        "Ryu 2024 Toxics 12:253, rat plasma", GFR_RAT,
                        round(r["filt"], 2), r["cl"],
                        "renal" if direct else "total_upper_bound",
                        round(r["T"], 4), verdict, mn, r["src"], r["note"]])

    report(rows, fu)
    print(f"\nwrote {p.relative_to(HERE)}")


def report(rows, fu):
    bar = "=" * 78
    print(bar)
    print("A.  THE MALE RAT PFOA NUMBER IS NOT SETTLED, AND THE DISAGREEMENT")
    print("    STRADDLES THE BOUNDARY")
    print(bar)
    filt = fu["PFOA"] * GFR_RAT
    print(f"  measured rat fu (Ryu 2024)          {fu['PFOA']:.5f}")
    print(f"  filtration clearance fu x GFR       {filt:.2f} mL/kg/day")
    print("  -- so any renal clearance above that figure is net SECRETION --")
    print()
    print(f"  {'source':52}{'CL_R':>8}{'T':>8}  reading")
    for r in rows:
        if r["chemical"] != "PFOA" or r["sex"] != "male" or r["tier"] != "DIRECT":
            continue
        rd = (f"reabsorbs {100*(1-r['T']):.0f}%" if r["T"] < 1
              else "NET SECRETION")
        print(f"  {r['src'][:52]:52}{r['cl']:>8.1f}{r['T']:>8.3f}  {rd}")
    print()
    print("  Two primary sources put the male rat at 51-58% reabsorbed. The")
    print("  third, Kudo 2002, puts it marginally into net secretion. The")
    print("  three span 2.5-fold and the filtration clearance sits inside that")
    print("  span, so the SIGN of the male rat tubular term depends on which")
    print("  primary source is taken.")
    print()
    print("  Kudo 2002 is internally consistent: its renal clearance is 91% of")
    print("  its own total clearance (0.032 of 0.035 mL/min/kg), which is what")
    print("  a mostly-renal compound should give. Kemper 2003 is unpublished")
    print("  and reached this project only through OEHHA Table A6.4; Ohmori's")
    print("  value was digitized from a graph by Han 2012.")
    print()
    print("  THE ROBUST STATEMENT IS THE ONE THAT SURVIVES ALL THREE:")
    print("  T_rat(PFOA, male) lies between 0.42 and 1.07, so THE MALE RAT")
    print("  TUBULE RETURNS AT MOST HALF THE FILTERED LOAD AND POSSIBLY NONE.")
    print("  Section 16 reported 57.6% reabsorbed on OEHHA's source alone;")
    print("  that figure is the most reabsorptive of the three, not a")
    print("  consensus.")

    print()
    print(bar)
    print("B.  WHAT THAT DOES TO THE RAT / HUMAN TUBULAR RATIO")
    print(bar)
    # human side, both available sources
    FU_H, GFR_H = 0.00245, 2570.0
    filt_h = FU_H * GFR_H
    for hsrc, hclr in (("OEHHA PHG Table A6.4", 0.06),
                       ("Argoul 2026 Table 1", 0.042)):
        Th = hclr / filt_h
        lo = 18.2 / filt
        hi = 0.032 * 1440 / filt
        print(f"  human CL_R {hclr:5.3f} ({hsrc:22}) -> T_human {Th:.5f}")
        print(f"      T_rat/T_human spans {lo/Th:6.0f}x to {hi/Th:6.0f}x")
    print()
    print("  Four source combinations, and the ratio spans 45x to 160x. THE")
    print("  NUMBER IS NOT WELL DETERMINED; THE DOMINANCE IS. Nothing in that")
    print("  range comes near 1, and the glomerular factor is 6.8x, so the")
    print("  tubule carries the rat/human difference under every combination.")
    print("  Section 16's '45x' is the low end of the range, not the estimate.")

    print()
    print(bar)
    print("C.  THE BOUNDED TIER: FOUR MORE COMPOUNDS, NO ASSUMPTIONS")
    print(bar)
    print("  CL_R <= CL_total, so T <= CL_total/(fu x GFR). Where that is")
    print("  below one, net reabsorption is PROVEN whatever the renal share.")
    print()
    print("  Where several rat studies report total clearance they disagree")
    print("  -- PFOS spans 1.63 to 13.1 mL/kg/day across four. Taking the")
    print("  LOWEST would give the strongest claim and stake it on the lowest")
    print("  value being right, so the claim below is built on the HIGHEST")
    print("  total clearance for each compound: the loosest bound, the one")
    print("  that holds whichever source is correct.")
    print()
    print(f"  {'chem':7}{'fu':>9}{'filtered':>10}{'CL_tot range':>16}"
          f"{'T_up':>7}  conclusion on the loosest bound")
    grp = {}
    for r in rows:
        if r["tier"] != "BOUNDED" or r["sex"] != "male":
            continue
        grp.setdefault(r["chemical"], []).append(r)
    worst = {}
    for k, rs_ in grp.items():
        worst[k] = max(rs_, key=lambda r: r["cl"])
    for k in sorted(worst, key=lambda x: worst[x]["T"]):
        r = worst[k]
        cls = [x["cl"] for x in grp[k]]
        rng = (f"{min(cls):.1f}" if len(cls) == 1
               else f"{min(cls):.1f}-{max(cls):.1f}")
        con = (f"PROVEN: reabsorbs >= {100*(1-r['T']):.0f}%" if r["T"] < 1
               else "uninformative, either sign")
        print(f"  {k:7}{r['fu']:>9.5f}{r['filt']:>10.1f}{rng:>16}"
              f"{r['T']:>7.2f}  {con}")
    proven = sorted(k for k, r in worst.items() if r["T"] < 1)
    print()
    print(f"  NET REABSORPTION HOLDS IN THE MALE RAT FOR {len(proven)} COMPOUNDS")
    print(f"  ON THE LOOSEST BOUND AVAILABLE: {', '.join(proven)}.")
    print("  No renal-share assumption anywhere, and no reliance on which of")
    print("  several conflicting total clearances is right. On the tightest")
    print("  bounds the same four reach 90%, 84%, 77% and 57%.")
    print()
    print("  'Proven' is conditional on one thing only: Ryu's measured rat fu.")
    print("  Everything else is an inequality. If that fu were two-fold low")
    print("  the filtered load would double and all four bounds would only")
    print("  tighten, so the direction of that dependence is favourable --")
    print("  the conclusion fails only if the true fu is LOWER than measured.")
    print()
    print("  PFHxS also has a DIRECT renal clearance (Kim 2016, 0.76), giving")
    print("  T = 0.07 and 93% reabsorbed, consistent with its bound.")
    print()
    print("  Uninformative for PFHxA, PFHpA, PFBS, PFDA and -- on this tier --")
    print("  PFOA: total clearance exceeds the filtered load, so the bound")
    print("  permits either sign. Section 11's direction table predicts")
    print("  secretion for the short chains, which is consistent, but a bound")
    print("  that permits either sign is not evidence and they stay open.")
    print()
    print("  SO THE RAT LAYER GOES FROM 1 OF 14 COMPOUNDS TO 6 OF 14:")
    print("    PFOA   direct, three conflicting sources, T 0.42 to 1.07")
    print("    PFHxS  direct, T 0.07, and bounded consistently")
    print("    PFNA, PFBA, PFOS  bounded, net reabsorption holds")
    print("    PFDA, PFHpA, PFBS, PFHxA  bounded, uninformative")

    print()
    print(bar)
    print("D.  A RAT STUDY MEASURED RAT fu BY SEX, AND IT ANSWERS THE")
    print("    POOLING QUESTION FOR ONE COMPOUND")
    print(bar)
    print("  Section 16 flagged that Ryu's plasma is equal-sex-pooled while")
    print("  the rat clearances are male, and that no study measures whether")
    print("  rat plasma fu differs by sex. Kim 2016 does, for PFHxS, by")
    print("  ultrafiltration:")
    print()
    print(f"    male    {FU_RAT_KIM[('PFHxS','male')]*100:.3f}%")
    print(f"    female  {FU_RAT_KIM[('PFHxS','female')]*100:.3f}%"
          f"   ratio {FU_RAT_KIM[('PFHxS','male')]/FU_RAT_KIM[('PFHxS','female')]:.2f}x")
    print(f"    Ryu 2024, sex-pooled, same species: {fu['PFHxS']*100:.3f}%")
    print()
    print("  Rat plasma binding of PFHxS does NOT differ appreciably by sex --")
    print("  1.1-fold, against a 7-fold sex difference in that compound's")
    print("  total clearance (6.7 male against 53.4 female, Sundstrom 2012).")
    print("  So for PFHxS the sex difference is NOT a binding difference, and")
    print("  pooling the sexes for fu is defensible.")
    print()
    print("  The absolute agreement with Ryu is close enough to be worth")
    print("  noting and not close enough to lean on: ultrafiltration and")
    print("  equilibrium dialysis are different methods and ultrafiltration is")
    print("  the more suspect of the two for PFAS, through binding to the")
    print("  filter. The male/female RATIO is internal to Kim's method, so it")
    print("  holds even if the absolute value does not.")
    print()
    print("  ONE COMPOUND IS NOT THE CLASS. PFOA is the compound where the rat")
    print("  sex difference is largest (70-fold in half-life) and it is still")
    print("  unmeasured by sex.")


if __name__ == "__main__":
    main()
