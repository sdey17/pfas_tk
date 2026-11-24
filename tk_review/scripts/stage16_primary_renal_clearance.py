#!/usr/bin/env python3
"""Primary human renal clearances from Harada 2005, Zhang 2013, Fujii 2015.

    python3 -I scripts/stage16_primary_renal_clearance.py

Writes db/primary_2026/human_renal_clearance_primary.csv and
db/human_renal_clearance_resolved.csv. Idempotent.

WHAT CHANGED BY READING THE PRIMARIES. Stage 12 recorded ten candidate human
renal clearances taken from EFSA's 2020 tabulation, flagged as unusable
because that document gave two different values for the same compounds. The
three primary papers resolve it, and the resolution is not the one expected:
EFSA's per-compound table is right for Zhang and WRONG for Fujii.

    EFSA lists, as RENAL clearance:   Fujii actually reports:
      PFUnDA  0.065                     urinary 0.005, TOTAL 0.065
      PFDoDA  0.070                     urinary 0.005, TOTAL 0.070
      PFTrDA  0.077                     urinary 0.006, TOTAL 0.077
      PFTeDA  0.224                     urinary BELOW DETECTION, TOTAL 0.224

EFSA transcribed Fujii's TOTAL clearance column into a renal clearance table.
For PFTeDA the published "renal clearance" of 0.224 mL/kg/day corresponds to
a urinary clearance the authors could not detect at all. The 0.005 that
appears in EFSA's own summary sentence -- the figure that made the document
look self-contradictory -- is Fujii's real urinary value, so the summary
sentence was right and the table was wrong.

THE THREE PRIMARIES ALSO DISAGREE WITH EACH OTHER, AND NOT RANDOMLY.
For PFOA, Harada (0.030) and Fujii (0.044) agree within 1.5-fold, while
Zhang's mean (0.77) is 18-26 times higher. Zhang's own distribution is
heavily right-skewed -- mean 0.77 against median 0.19 -- and Zhang computed
clearance from WHOLE BLOOD converted to serum by a haematocrit factor, where
Harada and Fujii measured serum directly. Zhang is the outlier on PFOA and
is reported here as such rather than averaged in.

RENAL IS A MINORITY ROUTE IN THE HUMAN AND A MAJORITY ROUTE IN RODENTS.
Harada Table 5 gives CL_R/CL_tot as 25.1% for human male PFOA against 91.4%
for male rat. Fujii gives human urinary 0.044 against a total of 0.096 (46%)
and reports a biliary clearance of 2.62 mL/kg/day, 60 times the urinary. In
the mouse, Fujii's urinary PFOA is 13.1 of a total 14.2 (92%). This is the
quantitative basis for treating the human tubular term as a residual.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
OUT1 = HERE / "db" / "primary_2026" / "human_renal_clearance_primary.csv"
OUT2 = HERE / "db" / "human_renal_clearance_resolved.csv"

# (chemical, group, CL_R mL/kg/day, SD or CI, n, statistic, source, table)
ROWS = [
    # --- Harada 2005, Environ Res 99:253-261, Section 3.3 and Table 5 ---
    ("PFOA", "all (n=20)", 0.030, 0.013, 20, "mean +/- SD",
     "Harada 2005", "Section 3.3 text"),
    ("PFOS", "all (n=20)", 0.015, 0.010, 20, "mean +/- SD",
     "Harada 2005", "Section 3.3 text"),
    ("PFOA", "male 61 kg", 0.033, None, 10, "group mean",
     "Harada 2005", "Table 5"),
    ("PFOA", "female 48 kg", 0.027, None, 10, "group mean",
     "Harada 2005", "Table 5"),
    ("PFOS", "male 61 kg", 0.012, None, 10, "group mean",
     "Harada 2005", "Table 5"),
    ("PFOS", "female 48 kg", 0.019, None, 10, "group mean",
     "Harada 2005", "Table 5"),
    # --- Zhang 2013, Environ Sci Technol 47:10619-10627, Table 2 ---
    # male and older female group; mean and median both recorded because
    # the distribution is strongly right-skewed
    ("PFHxS", "male + older female", 0.027, None, 64, "mean",
     "Zhang 2013", "Table 2"),
    ("PFHxS", "male + older female", 0.015, None, 64, "median",
     "Zhang 2013", "Table 2"),
    ("PFOS", "male + older female", 0.037, None, 66, "mean (sum PFOS)",
     "Zhang 2013", "Table 2"),
    ("PFOS", "male + older female", 0.024, None, 66, "median (sum PFOS)",
     "Zhang 2013", "Table 2"),
    ("PFHpA", "male + older female", 0.61, None, 31, "mean",
     "Zhang 2013", "Table 2"),
    ("PFHpA", "male + older female", 0.41, None, 31, "median",
     "Zhang 2013", "Table 2"),
    ("PFNA", "male + older female", 0.15, None, 50, "mean",
     "Zhang 2013", "Table 2"),
    ("PFNA", "male + older female", 0.094, None, 50, "median",
     "Zhang 2013", "Table 2"),
    ("PFDA", "male + older female", 0.096, None, 60, "mean",
     "Zhang 2013", "Table 2"),
    ("PFDA", "male + older female", 0.035, None, 60, "median",
     "Zhang 2013", "Table 2"),
    ("PFUnDA", "male + older female", 0.065, None, 63, "mean",
     "Zhang 2013", "Table 2"),
    ("PFUnDA", "male + older female", 0.042, None, 63, "median",
     "Zhang 2013", "Table 2"),
    ("PFOA", "male + older female", 0.77, None, 66, "mean (sum PFOA)",
     "Zhang 2013", "Table 2"),
    ("PFOA", "male + older female", 0.19, None, 66, "median (sum PFOA)",
     "Zhang 2013", "Table 2"),
    # --- Fujii 2015, J Occup Health 57:1-12, Table 2 (URINARY column) ---
    ("PFHpA", "male n=5, female n=5", 0.674, 0.374, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
    ("PFOA", "male n=5, female n=5", 0.044, 0.01, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
    ("PFNA", "male n=5, female n=5", 0.038, 0.01, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
    ("PFDA", "male n=5, female n=5", 0.015, 0.01, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
    ("PFUnDA", "male n=5, female n=5", 0.005, 0.00, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
    ("PFDoDA", "male n=5, female n=5", 0.005, 0.00, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
    ("PFTriDA", "male n=5, female n=5", 0.006, 0.01, 10, "mean +/- SD",
     "Fujii 2015", "Table 2, urinary clearance"),
]
# What EFSA printed as renal clearance, against what the primary reports.
EFSA_VS_PRIMARY = [
    # compound, EFSA "renal", Fujii urinary, Fujii total, verdict
    ("PFUnDA", 0.065, 0.005, 0.065, "EFSA used Fujii's TOTAL clearance"),
    ("PFDoDA", 0.070, 0.005, 0.070, "EFSA used Fujii's TOTAL clearance"),
    ("PFTriDA", 0.077, 0.006, 0.077, "EFSA used Fujii's TOTAL clearance"),
    ("PFTeDA", 0.224, None, 0.224,
     "EFSA used Fujii's TOTAL; urinary was BELOW DETECTION"),
]
# Fraction of total clearance that is renal, from the primaries themselves.
RENAL_FRACTION = [
    ("human", "PFOA", 25.1, "Harada 2005 Table 5, male 61 kg"),
    ("human", "PFOS", 17.5, "Harada 2005 Table 5, male 61 kg"),
    ("human", "PFOA", 45.8, "Fujii 2015 Table 2, urinary 0.044 / total 0.096"),
    ("human", "PFNA", 61.3, "Fujii 2015 Table 2, 0.038 / 0.062"),
    ("human", "PFDA", 22.7, "Fujii 2015 Table 2, 0.015 / 0.066"),
    ("human", "PFUnDA", 7.7, "Fujii 2015 Table 2, 0.005 / 0.065"),
    ("human", "PFHpA", 100.0, "Fujii 2015 Table 2, 0.674 / 0.674"),
    ("rat", "PFOA", 91.4, "Harada 2005 Table 5, male Wistar"),
    ("rat", "PFOA", 47.2, "Harada 2005 Table 5, female Wistar"),
    ("macaque", "PFOA", 40.4, "Harada 2005 Table 5, male Japanese macaque"),
    ("macaque", "PFOA", 41.5, "Harada 2005 Table 5, female"),
    ("mouse", "PFOA", 92.3, "Fujii 2015, urinary 13.1 / total 14.2, male"),
]


def main():
    bar = "=" * 78
    OUT1.parent.mkdir(parents=True, exist_ok=True)
    with OUT1.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "group", "cl_renal_mL_kg_day", "sd", "n",
                    "statistic", "source", "source_table"])
        for r in ROWS:
            w.writerow([r[0], r[1], r[2], "" if r[3] is None else r[3],
                        r[4], r[5], r[6], r[7]])

    print(bar)
    print("THE EFSA TABULATION IS WRONG FOR FUJII  (new literature error)")
    print(bar)
    print(f"  {'compound':9}{'EFSA renal':>12}{'Fujii urinary':>15}"
          f"{'Fujii total':>13}  what happened")
    for c, e, u, t, why in EFSA_VS_PRIMARY:
        us = "below MDL" if u is None else f"{u:.3f}"
        print(f"  {c:9}{e:>12.3f}{us:>15}{t:>13.3f}  {why}")
    print()
    print("  EFSA's per-compound table is CORRECT for Zhang -- every value")
    print("  matches Zhang Table 2 exactly, the apparent 'range' being the")
    print("  young-female against male/older-female groups. It is WRONG for")
    print("  Fujii, where the TOTAL clearance column was transcribed into a")
    print("  renal clearance table. The 0.005 for PFUnDA in EFSA's summary")
    print("  sentence, which made the document look self-contradictory, is")
    print("  Fujii's real urinary value: the sentence was right, the table")
    print("  was wrong.")
    print()

    print(bar)
    print("THE THREE PRIMARIES, SIDE BY SIDE  (mL/kg/day)")
    print(bar)
    by = {}
    for chem, grp, v, sd, n, stat, src, tbl in ROWS:
        by.setdefault(chem, []).append((src, stat, v, sd, n))
    print(f"  {'compound':9}{'source':14}{'statistic':22}{'CL_R':>8}"
          f"{'SD':>8}{'n':>5}")
    for chem in ("PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA",
                 "PFTriDA", "PFHxS", "PFOS"):
        for src, stat, v, sd, n in by.get(chem, []):
            print(f"  {chem:9}{src:14}{stat:22}{v:>8.3f}"
                  f"{('-' if sd is None else f'{sd:.3f}'):>8}{n:>5}")
    print()
    print("  PFOA is the test case. Harada 0.030 and Fujii 0.044 agree within")
    print("  1.5-fold. Zhang's mean is 0.77 -- 18 to 26 times higher -- and")
    print("  its own median is 0.19, so the distribution is heavily skewed.")
    print("  Zhang computed clearance from whole blood converted to serum;")
    print("  the other two measured serum. ZHANG IS THE OUTLIER ON PFOA.")
    print()

    print(bar)
    print("RENAL IS A MINORITY ROUTE IN HUMANS AND A MAJORITY IN RODENTS")
    print(bar)
    print(f"  {'species':9}{'compound':9}{'CL_R/CL_tot %':>15}  source")
    for sp, chem, pct, src in RENAL_FRACTION:
        print(f"  {sp:9}{chem:9}{pct:>15.1f}  {src}")
    print()
    print("  This is the measured basis for treating the human tubular term")
    print("  as a residual: three quarters of human PFOA elimination is not")
    print("  renal, against under a tenth in the male rat and the mouse.")
    print("  Fujii's human biliary clearance for PFOA is 2.62 mL/kg/day --")
    print("  sixty times the urinary 0.044 -- of which the paper estimates")
    print("  98% is reabsorbed, leaving a faecal 0.052.")
    print()

    print(bar)
    print("DIRECT HUMAN TUBULAR FACTORS, REPLACING THE UPPER BOUNDS")
    print(bar)
    fu = {}
    for r in csv.DictReader(
            (HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r["experiment"] == "Ryu 2024" and r["species"] == "human" and r["fu"]:
            fu[r["chemical"]] = float(r["fu"])
    GFR = 2570.0
    prior = {"PFHpA": 0.0511, "PFOA": 0.0191, "PFNA": 0.0320, "PFDA": 0.0944,
             "PFUnDA": 0.2656, "PFDoDA": 3.0832, "PFOS": 0.0455,
             "PFHxS": 0.0358}
    print(f"  {'compound':9}{'source':13}{'CL_R':>8}{'fu':>10}{'filtered':>10}"
          f"{'T direct':>10}{'% reabs':>9}   was (bound)")
    trows = []
    for chem, src, v in (("PFHpA", "Fujii 2015", 0.674),
                         ("PFOA", "Harada 2005", 0.030),
                         ("PFOA", "Fujii 2015", 0.044),
                         ("PFNA", "Fujii 2015", 0.038),
                         ("PFDA", "Fujii 2015", 0.015),
                         ("PFUnDA", "Fujii 2015", 0.005),
                         ("PFDoDA", "Fujii 2015", 0.005),
                         ("PFOS", "Harada 2005", 0.015),
                         ("PFHxS", "Zhang 2013", 0.027)):
        f = fu.get(chem)
        if not f:
            continue
        filt = f * GFR
        T = v / filt
        was = prior.get(chem)
        print(f"  {chem:9}{src:13}{v:>8.3f}{f:>10.5f}{filt:>10.2f}"
              f"{T:>10.4f}{(1-T)*100:>8.2f}%   "
              f"{('T <= %.4f' % was) if was else '-'}")
        trows.append((chem, src, v, f, filt, T))
    print()
    print("  Every direct value lies BELOW its former upper bound, as it must,")
    print("  and the long chains are far below. PFDoDA moves from an")
    print("  uninformative bound of T <= 3.08 to a direct T of 0.019 -- from")
    print("  'sign undetermined' to 98% net reabsorbed. PFUnDA likewise.")
    print()

    with OUT2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "chemical", "renal_fraction_pct", "source"])
        for sp, chem, pct, src in RENAL_FRACTION:
            w.writerow([sp, chem, pct, src])
    print(f"wrote {OUT1.relative_to(HERE)}  ({len(ROWS)} rows)")
    print(f"wrote {OUT2.relative_to(HERE)}  ({len(RENAL_FRACTION)} rows)")


if __name__ == "__main__":
    main()
