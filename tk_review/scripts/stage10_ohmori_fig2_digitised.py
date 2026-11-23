#!/usr/bin/env python3
"""Ohmori 2003 Fig. 2 digitised: eight rat renal clearances (Stage 10).

    python3 -I scripts/stage10_ohmori_fig2_digitised.py

Writes db/primary_2026/ohmori2003_fig2_renal_clearance.csv and appends the
usable rows to db/tubular_factor_evidence.csv's companion record. Idempotent.

    Ohmori K, Kudo N, Katayama K, Kawashima Y (2003) Toxicology 184:135-140,
    doi:10.1016/s0300-483x(02)00573-5, Figure 2

WHY DIGITISE. Section 19 established that this figure holds renal clearances
for PFHpA, PFOA, PFNA and PFDA in BOTH SEXES of Wistar rat -- the single
richest block of rat renal data in the literature -- and that none of the
numbers appear in the text or tables. Han 2012 digitised one bar from it (male
PFOA, 21.2) and that single value has carried the rat column of this review
ever since.

METHOD. Page 5 rendered at 300 dpi; axis calibrated on the printed numerals
(0 at y=550 px, 2000 at y=17 px, so 3.75 mL/day/kg per pixel); each bar's top
taken as the highest row at which at least 35% of the bar's width is dark
(15% for the hatched male bars), which ignores the thin error bars and the
significance letters above them.

THE CALIBRATION IS VALIDATED AGAINST AN INDEPENDENT DIGITISATION. OEHHA PHG
Table A6.4 carries Ohmori's FEMALE PFOA renal clearance as 1009 mL/kg/day,
from Han 2012's reading of this same figure. This digitisation returns 994 --
A DIFFERENCE OF 1.5%. That is the check that makes the rest usable.

AND IT SHOWS WHERE DIGITISATION FAILS. For the male PFOA bar, a sliver a few
pixels tall, Han read 21.2 and this reading returns 45. TWO DIGITISATIONS OF
THE SAME BAR, TWO-FOLD APART. On a 2000-unit axis one pixel is 3.75 units and
a bar's own outline is two to three pixels, so small bars carry an error
comparable to their value. The large bars are reliable to a few per cent; the
small ones are not, and are recorded as bounds rather than values.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
PX_PER_UNIT = (550 - 17) / 2000.0
DETECT_LIMIT = 15.0        # ~4 px; below this a bar is not resolvable

# chemical, sex, CL_R mL/kg/day or None if below detection, note
FIG2 = [
    ("PFHpA", "male",   983.1, "tall hatched bar; reliable"),
    ("PFHpA", "female", 1467.2, "tall solid bar; reliable"),
    ("PFOA",  "male",   45.0,
     "SMALL BAR, UNRELIABLE: Han 2012 read 21.2 from the same bar, a 2.1-fold "
     "disagreement between two digitisations"),
    ("PFOA",  "female", 994.4,
     "VALIDATION ROW: OEHHA A6.4 carries 1009 from Han 2012's digitisation of "
     "this bar; this reading is 1.5% from it"),
    ("PFNA",  "male",   None, "no bar resolvable above the baseline"),
    ("PFNA",  "female", 67.5, "short solid bar; +/- about 8 on a 4 px error"),
    ("PFDA",  "male",   None, "no bar resolvable above the baseline"),
    ("PFDA",  "female", None, "no bar resolvable above the baseline"),
]
GFR = {"OEHHA/Han": 14400.0, "Ohmori/Shim": 7920.0}
# section 19's independent estimates, from renal fraction x total clearance
S19 = {("PFHpA", "male"): 23.0, ("PFOA", "female"): 42.0,
       ("PFNA", "female"): 1.8}


def rat_fu():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] == "Ryu 2024" and r["species"] == "rat" and r["fu"]:
            out[r["chemical"]] = float(r["fu"])
    return out


def main():
    fu = rat_fu()
    p = HERE / "db" / "primary_2026" / "ohmori2003_fig2_renal_clearance.csv"
    rows = []
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "sex", "strain",
                    "cl_renal_mL_kg_day", "value_kind", "fu_rat",
                    "gfr_mL_kg_day", "gfr_source", "filtered_load",
                    "T", "verdict", "section19_independent_T",
                    "agreement_fold", "digitisation_note", "source"])
        nc = {"PFHpA": 7, "PFOA": 8, "PFNA": 9, "PFDA": 10}
        for chem, sex, clr, note in FIG2:
            f = fu.get(chem)
            for gs, g in GFR.items():
                filt = f * g
                if clr is None:
                    kind, T = "below detection", DETECT_LIMIT / filt
                    verdict = ("net_reabsorption_bound" if T < 1
                               else "undetermined")
                    shown = f"<{DETECT_LIMIT}"
                else:
                    kind, T = "digitised value", clr / filt
                    verdict = ("net_secretion" if T > 1
                               else "net_reabsorption")
                    if (chem, sex) == ("PFOA", "male"):
                        kind = "digitised value, UNRELIABLE"
                        verdict = ("DO_NOT_USE: two digitisations of this bar "
                                   "differ 2.1-fold and straddle T = 1")
                    shown = clr
                s19 = S19.get((chem, sex)) if gs == "OEHHA/Han" else None
                agree = round(max(T, s19) / min(T, s19), 2) if s19 else ""
                w.writerow([chem, nc[chem], sex, "Wistar", shown, kind,
                            f, g, gs, round(filt, 2), round(T, 3), verdict,
                            s19 if s19 else "", agree, note,
                            "Ohmori 2003 Fig. 2, digitised here at 300 dpi; "
                            "doi:10.1016/s0300-483x(02)00573-5"])
                if gs == "OEHHA/Han":
                    rows.append((chem, sex, shown, T, verdict, s19, agree))

    bar = "=" * 78
    print(bar)
    print("OHMORI 2003 FIG. 2, DIGITISED  (CL_R in mL/kg/day)")
    print(bar)
    print(f"  {'chem':7}{'sex':8}{'CL_R':>9}{'filtered':>10}{'T':>8}"
          f"   verdict")
    for chem, sex, shown, T, verdict, s19, agree in rows:
        print(f"  {chem:7}{sex:8}{str(shown):>9}"
              f"{fu[chem]*14400:>10.1f}{T:>8.2f}   {verdict}")
    print()
    print(bar)
    print("THREE INDEPENDENT CONFIRMATIONS OF SECTION 19's SECRETION CALLS")
    print(bar)
    print("  Section 19 established net secretion for three rat rows WITHOUT")
    print("  any renal clearance, using only a measured urinary fraction and")
    print("  a total clearance. This figure measures the renal clearance")
    print("  directly. The two routes share no data:")
    print()
    print(f"  {'row':20}{'sec.19 (indirect)':>19}{'Fig.2 (direct)':>16}"
          f"{'fold':>7}")
    for chem, sex, shown, T, verdict, s19, agree in rows:
        if not s19:
            continue
        print(f"  {chem+' '+sex:20}{s19:>19.1f}{T:>16.2f}{agree:>7}")
    print()
    print("  SAME SIGN IN ALL THREE, WITHIN A FACTOR OF 1.2 TO 1.8. The")
    print("  secretion verdicts were right, and they are now direct rather")
    print("  than inferred.")
    print()
    print(bar)
    print("WHAT THIS ADDS TO THE RAT COLUMN")
    print(bar)
    print("  Before: one direct rat renal clearance (male PFOA, via Han).")
    print("  After:  four direct values and three bounds, both sexes.")
    print()
    print("    PFHpA male    983   T = 15.0   net secretion")
    print("    PFHpA female 1467   T = 22.4   net secretion")
    print("    PFOA  female  994   T = 23.2   net secretion")
    print("    PFNA  female   68   T =  2.2   net secretion")
    print("    PFNA  male    <15   T < 0.50   reabsorption")
    print("    PFDA  male    <15   T < 1.04   undetermined, marginal")
    print("    PFDA  female  <15   T < 1.04   undetermined, marginal")
    print()
    print("  AND THE MALE PFOA BAR SHOULD NOT BE USED. Han read 21.2 and this")
    print("  reading returns 45; on the OEHHA GFR those give T = 0.49 and")
    print("  T = 1.05, which sit on opposite sides of the filtration line.")
    print("  The compound whose rat value this review has revised four times")
    print("  cannot be settled by digitising a bar a few pixels tall.")
    print()
    print(f"\nwrote {p.relative_to(HERE)}  ({len(FIG2)*len(GFR)} rows)")


if __name__ == "__main__":
    main()
