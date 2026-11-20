#!/usr/bin/env python3
"""Ohmori 2003 from the paper itself, and a third evidence tier (Stage 4f).

    python3 -I scripts/stage4f_ohmori2003_primary.py

Writes db/primary_2026/ohmori2003_table1_primary.csv and
db/rat_tubular_factor_ohmori.csv. Idempotent.

    Ohmori K, Kudo N, Katayama K, Kawashima Y (2003). "Comparison of the
    toxicokinetics between perfluorocarboxylic acids with different carbon
    chain length." Toxicology 184(2-3):135-140. PMID 12499116,
    doi:10.1016/s0300-483x(02)00573-5

DESIGN. Male and female Wistar rats (SLC, Hamamatsu), single IV injection of
48.64 UMOL/kg -- the paper prints "mmol" in the extracted text because
pdftotext drops the micro sign, the same defect that corrupted Uchida 2016 in
this project. 48.64 mmol/kg of PFOA would be 20 g/kg; 48.64 umol/kg is 20
mg/kg, which is the dose OEHHA records. Two-compartment open model, CLtot =
dose / AUC(0-inf) by MULTI. CLR = amount excreted in urine 0 to t divided by
AUC(0-t), over 120 min for PFHpA, 240 (male) or 360 (female) for PFOA, and 360
for PFNA and PFDA, under saline infusion to maintain urine flow.

FOUR THINGS THE PAPER SETTLES.

1. TWO ATSDR VALUES ARE WRONG, ONE BY FORTY-FOLD.

      male PFOA CLtot      ATSDR 2021 Table 3-6: 135     Ohmori:   50.5
      male PFDA CLtot      ATSDR:                207     Ohmori:    5.2
      female PFDA CLtot    ATSDR:                140     Ohmori:    5.3

   The other five values ATSDR attributes to this paper are right. The PFDA
   error is the consequential one: 207 was the sole reason rat PFDA sat as
   "uninformative" in Stage 4e, and at the true 5.2 the compound closes.

2. AND ONE OF MY OWN CLAIMS WAS WRONG. Stage 4e argued ATSDR's 135 was an
   error because Ohmori's "own Vss of 196 and half-life of 5.63 d" give
   ln2 x Vd / t = 24.1. Both halves of that were wrong. 196.2 is Ohmori's Vss
   for PFHpA, not PFOA -- this project's vd_clearance.csv has it on the wrong
   row -- and PFOA's is 338.8. And ln2 x Vd / t assumes one compartment where
   the paper fits two, so the test did not apply. The conclusion survives
   (135 is wrong) but the reasoning did not, and the true value, 50.5, is
   twice my estimate. It agrees with Kudo 2002's 50.4 from the same
   laboratory, which is the cross-check that should have been used.

3. OHMORI'S RAT GFR IS 5.5 ml/(min kg), NOT OEHHA'S 10. Citing Shim 1984.
   That is 7920 mL/kg/day against 14400, a 1.82-fold fork in a quantity this
   project has been treating as known, and it moves every rat T upward by that
   factor. Both are reported below.

4. OHMORI CONCLUDES TUBULAR SECRETION, BY THE SAME TEST THIS PROJECT USES.
   "Some of the PFCAs are thought to be actively excreted into urine by renal
   tubular secretion because CLR of free PFCA was higher than the glomerular
   filtration rate ... if plasma protein binding in vivo is the same as that
   observed in vitro", supported by probenecid inhibition in Kudo 2002.
   CLR/fu > GFR is exactly T > 1. The paper reached the T framework in 2003.
   Its own binding is unusable quantitatively -- ">98% bound, data not shown",
   by ultrafiltration at 300 nmol/mL, sixty times Ryu's 5 uM and close to rat
   albumin saturation -- so the test is run here on Ryu's measured fu instead.

A THIRD TIER: A LOWER BOUND ON T FROM THE RENAL FRACTION.

Stage 4e had DIRECT and an upper bound from total clearance. The mirror also
works. With f_ren the fraction of elimination that is renal,

    CLR = f_ren x CLtot        so        T = f_ren x CLtot / (fu x GFR)

and T > 1 as soon as

    f_ren  >  (fu x GFR) / CLtot  ==  the THRESHOLD below

Where a study IN THE SAME SPECIES measures urinary recovery above that
threshold, NET SECRETION IS ESTABLISHED. Kudo 2001 measured it in Wistar rats,
from the same laboratory as Ohmori. Where the threshold exceeds 1 secretion is
impossible and reabsorption is certain, which recovers the Stage 4e upper
bound as a special case.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

# ---- Ohmori 2003 Table 1, read off the paper ---------------------------
# chemical: {sex: (CLtot mL/kg/day, SD, t_half d, SD, Vss mL/kg, SD)}
TABLE1 = {
    "PFHpA": {"male":   (1603.9, 557.6, 0.10, 0.05, 196.2, 18.6),
              "female": (3070.5, 781.3, 0.05, 0.01, 200.6, 29.0)},
    "PFOA":  {"male":   (  50.5,  17.0, 5.63, 1.20, 338.8, 68.1),
              "female": (2233.2, 805.0, 0.08, 0.03, 211.2, 28.2)},
    "PFNA":  {"male":   (   6.9,   0.6, 29.57, 2.28, 286.8, 12.5),
              "female": ( 105.7,  31.1, 2.44, 0.41, 243.1, 48.9)},
    "PFDA":  {"male":   (   5.2,   1.3, 39.92, 8.62, 347.7, 15.2),
              "female": (   5.3,   0.2, 58.57, 5.84, 441.1, 55.1)},
}
ATSDR = {("PFHpA", "male"): 1604.0, ("PFHpA", "female"): 3071.0,
         ("PFOA", "male"): 135.0,   ("PFOA", "female"): 2233.0,
         ("PFNA", "male"): 6.9,     ("PFNA", "female"): 106.0,
         ("PFDA", "male"): 207.0,   ("PFDA", "female"): 140.0}

# Kudo 2001, Wistar rats, same laboratory: % of dose in urine
# (the fraction of ELIMINATION that is renal, as measured)
F_REN = {("PFHpA", "male"): 0.92, ("PFOA", "male"): 0.55,
         ("PFNA", "male"): 0.020, ("PFDA", "male"): 0.002,
         ("PFOA", "female"): 0.80, ("PFNA", "female"): 0.51}
F_REN_SRC = "Kudo 2001, % of dose in urine, Wistar rats, same laboratory"

GFRS = [(14400.0, "OEHHA PHG Table A6.4 / Han 2012 (10 mL/min/kg)"),
        (7920.0, "Ohmori 2003, citing Shim 1984 (5.5 mL/min/kg)")]


def ryu_rat_fu():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] == "Ryu 2024" and r["species"] == "rat" and r["fu"]:
            out[r["chemical"]] = float(r["fu"])
    return out


def main():
    fu = ryu_rat_fu()

    # ---- the primary table ------------------------------------------
    p1 = HERE / "db" / "primary_2026" / "ohmori2003_table1_primary.csv"
    with p1.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "sex", "strain",
                    "cltot_mL_kg_day", "cltot_sd", "halflife_d",
                    "halflife_sd", "vss_mL_kg", "vss_sd",
                    "atsdr_2021_value", "atsdr_agrees", "source"])
        nc = {"PFHpA": 7, "PFOA": 8, "PFNA": 9, "PFDA": 10}
        for chem, bysex in TABLE1.items():
            for sex, (cl, clsd, t, tsd, v, vsd) in bysex.items():
                a = ATSDR[(chem, sex)]
                agree = "yes" if abs(a - cl) / cl < 0.05 else "NO"
                w.writerow([chem, nc[chem], sex, "Wistar", cl, clsd, t, tsd,
                            v, vsd, a, agree,
                            "Ohmori 2003 Table 1, Toxicology 184:135-140, "
                            "doi:10.1016/s0300-483x(02)00573-5"])

    # ---- the three-tier rat table -----------------------------------
    rows = []
    for gfr, gsrc in GFRS:
        for chem, bysex in TABLE1.items():
            f = fu.get(chem)
            if f is None:
                continue
            filt = f * gfr
            for sex, (cl, *_rest) in bysex.items():
                thresh = filt / cl
                fr = F_REN.get((chem, sex))
                T_up = cl / filt
                if thresh > 1:
                    verdict = "net_reabsorption_CERTAIN"
                    detail = f"reabsorbs >= {100*(1-T_up):.0f}%"
                elif fr is None:
                    verdict = "undetermined_no_urinary_recovery"
                    detail = f"secretion if renal fraction > {100*thresh:.1f}%"
                elif fr > thresh:
                    verdict = "NET_SECRETION_ESTABLISHED"
                    detail = (f"renal fraction {100*fr:.1f}% exceeds the "
                              f"{100*thresh:.1f}% threshold; T ~ "
                              f"{fr*cl/filt:.1f}")
                else:
                    verdict = "undetermined_recovery_below_threshold"
                    detail = (f"renal fraction {100*fr:.1f}% is below the "
                              f"{100*thresh:.1f}% threshold")
                rows.append(dict(chemical=chem, sex=sex, gfr=gfr, gsrc=gsrc,
                                 fu=f, filt=filt, cltot=cl, thresh=thresh,
                                 f_ren=fr, T_up=T_up, verdict=verdict,
                                 detail=detail))

    p2 = HERE / "db" / "rat_tubular_factor_ohmori.csv"
    with p2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "sex", "fu_rat", "fu_source",
                    "gfr_mL_kg_day", "gfr_source",
                    "filtered_load_mL_kg_day", "cltot_mL_kg_day",
                    "cltot_source", "secretion_threshold_renal_fraction",
                    "measured_renal_fraction", "renal_fraction_source",
                    "T_upper_bound", "verdict", "detail"])
        for r in rows:
            w.writerow([r["chemical"], r["sex"], r["fu"],
                        "Ryu 2024 Toxics 12:253, rat plasma", r["gfr"],
                        r["gsrc"], round(r["filt"], 2), r["cltot"],
                        "Ohmori 2003 Table 1", round(r["thresh"], 4),
                        r["f_ren"] if r["f_ren"] is not None else "",
                        F_REN_SRC if r["f_ren"] is not None else "",
                        round(r["T_up"], 4), r["verdict"], r["detail"]])

    report(rows, fu)
    print(f"\nwrote {p1.relative_to(HERE)}")
    print(f"wrote {p2.relative_to(HERE)}  ({len(rows)} rows)")


def report(rows, fu):
    bar = "=" * 78
    print(bar)
    print("A.  OHMORI 2003 TABLE 1 AGAINST ATSDR 2021 TABLE 3-6")
    print(bar)
    print(f"  {'':16}{'Ohmori':>10}{'ATSDR':>10}{'ratio':>8}")
    for chem, bysex in TABLE1.items():
        for sex, (cl, *_r) in bysex.items():
            a = ATSDR[(chem, sex)]
            flag = "" if abs(a - cl) / cl < 0.05 else "  <-- WRONG"
            print(f"  {chem+' '+sex:16}{cl:>10.1f}{a:>10.1f}"
                  f"{a/cl:>8.1f}{flag}")
    print()
    print("  Five of eight agree. Male PFOA is 2.7x too high in ATSDR and")
    print("  BOTH PFDA values are 26-40x too high. The PFDA error is the")
    print("  consequential one.")

    print()
    print(bar)
    print("B.  CORRECTING STAGE 4E'S OWN REASONING")
    print(bar)
    print("  Stage 4e said ATSDR's 135 must be wrong because Ohmori's Vss of")
    print("  196 and half-life of 5.63 d give ln2 x Vd / t = 24.1.")
    print()
    print("  Both halves of that were wrong:")
    print("    - 196.2 is Ohmori's Vss for PFHpA. PFOA's is 338.8. This")
    print("      project's vd_clearance.csv carries 196 on the PFOA row.")
    print("    - the paper fits a TWO-compartment model, so ln2 x Vd / t does")
    print("      not apply; CLtot = dose/AUC, and Vss/CLtot = MRT = 6.7 d,")
    print("      against t_half/ln2 = 8.1 d, which is self-consistent for two")
    print("      compartments.")
    print()
    print("  The conclusion held -- 135 is wrong -- but by luck. The true")
    print("  value, 50.5, is twice the estimate, and it matches Kudo 2002's")
    print("  50.4 from the same laboratory. THAT was the available")
    print("  cross-check and it should have been used instead of a")
    print("  one-compartment identity.")

    print()
    print(bar)
    print("C.  THE RAT GFR IS ALSO NOT AGREED")
    print(bar)
    for g, s in GFRS:
        print(f"  {g:>8.0f} mL/kg/day   {s}")
    print(f"  ratio {GFRS[0][0]/GFRS[1][0]:.2f}x, and every rat T scales "
          f"inversely with it")
    print()
    print("  Rat PFOA, male, under both:")
    for g, s in GFRS:
        filt = fu["PFOA"] * g
        print(f"    GFR {g:>6.0f}: filtered {filt:6.2f}  ->  "
              f"T = {18.2/filt:.2f} (Kemper), {21.2/filt:.2f} (Ohmori), "
              f"{46.08/filt:.2f} (Kudo)")
    print()
    print("  ON OHMORI'S OWN GFR THE MALE RAT BARELY REABSORBS PFOA AT ALL --")
    print("  at most about 23% on the lowest renal clearance, and clearly")
    print("  secreting on the highest. Section 17 said 'at most half,")
    print("  possibly none' on OEHHA's GFR; on Ohmori's it is 'at most a")
    print("  quarter, possibly net secretion'. The direction of section 16's")
    print("  conclusion is unaffected -- the human tubule reabsorbs over 99%")
    print("  on any of these -- but the rat figure is softer again.")

    print()
    print(bar)
    print("D.  THE THIRD TIER: SECRETION ESTABLISHED FROM THE RENAL FRACTION")
    print(bar)
    print("  T = f_ren x CLtot / (fu x GFR), so secretion needs")
    print("  f_ren > (fu x GFR)/CLtot. Kudo 2001 measured f_ren in Wistar")
    print("  rats in the same laboratory. On OEHHA's GFR:")
    print()
    print(f"  {'':14}{'CLtot':>9}{'filtered':>10}{'threshold':>11}"
          f"{'measured':>10}   verdict")
    for r in rows:
        if r["gfr"] != 14400.0:
            continue
        th = "impossible" if r["thresh"] > 1 else f"{100*r['thresh']:.1f}%"
        me = "-" if r["f_ren"] is None else f"{100*r['f_ren']:.1f}%"
        v = {"net_reabsorption_CERTAIN": "REABSORPTION certain",
             "NET_SECRETION_ESTABLISHED": "SECRETION established",
             "undetermined_no_urinary_recovery": "undetermined, no recovery",
             "undetermined_recovery_below_threshold":
                 "undetermined"}[r["verdict"]]
        print(f"  {r['chemical']+' '+r['sex']:14}{r['cltot']:>9.1f}"
              f"{r['filt']:>10.2f}{th:>11}{me:>10}   {v}")
    print()
    print("  SO THE PAPER CLOSES FOUR CELLS, IN BOTH DIRECTIONS:")
    print("    PFDA male and female   REABSORPTION certain, >= 64% and 63%,")
    print("                           because the true CLtot is 5.2 and 5.3")
    print("                           rather than ATSDR's 207 and 140")
    print("    PFHpA male             NET SECRETION, threshold 4.1% against a")
    print("                           measured 92% renal; T about 23")
    print("    PFOA female            NET SECRETION, threshold 1.9% against a")
    print("                           measured 80% renal; T about 42")
    print("    PFNA female            NET SECRETION, threshold 28.6% against")
    print("                           a measured 51% renal; T about 1.8")
    print()
    print("  AND THE MALE/FEMALE CONTRAST IS NOW IN THE SAME UNITS AS")
    print("  EVERYTHING ELSE. The female rat NET SECRETES PFOA and PFNA while")
    print("  the male rat does not: male PFOA sits at a threshold of 85%")
    print("  against a measured 55% renal fraction, so it fails the test, and")
    print("  male PFNA cannot secrete at all (threshold 438%). The 70-fold")
    print("  sex difference in PFOA half-life is a difference between a")
    print("  secreting kidney and a reabsorbing one, not a difference of")
    print("  degree within either.")
    print()
    print("  ROBUST TO THE fu POOLING. Ryu's rat fu is sex-pooled, which")
    print("  section 16 flagged. The female PFOA verdict needs only a 1.9%")
    print("  renal fraction against 80% measured, so a TEN-FOLD error in fu")
    print("  would raise the threshold to 19% and leave the verdict intact.")
    print("  Female PFNA is the tight one: threshold 28.6% against 51%, so a")
    print("  1.8-fold error in fu would overturn it.")


if __name__ == "__main__":
    main()
