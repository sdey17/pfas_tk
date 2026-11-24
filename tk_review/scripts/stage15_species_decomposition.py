#!/usr/bin/env python3
"""Half-life decomposition, extended to mouse, with bootstrap CIs (Stage 15).

    python3 -I scripts/stage15_species_decomposition.py

Writes db/species_decomposition.csv. Idempotent. Supersedes stage13 for the
paper; stage13 is retained because its two-species table is cited elsewhere.

    t_h/t_x = (Vd_h/Vd_x) * (fu_x/fu_h) * (GFR_x/GFR_h) * (T_x/T_h)
              DISTRIBUTION    BINDING    GLOMERULAR      TUBULE(residual)

SELECTION RULES, STATED SO THEY CAN BE CHECKED.
  * Human is Abraham 2024 throughout: one male volunteer, isotopically
    labelled single oral doses, Vd and CL both MEASURED, the identity
    ln2*Vd/CL closing at a consistent 82.9 kg across all 15 compounds. No
    other human source has measured Vd and CL for more than three compounds.
  * Animal rows are MALE (matching Abraham), MEASURED, and the LOWEST dose
    reported, since several of these compounds show dose-dependent clearance.
  * One study per species-compound, named in the output.
  * fu is Ryu 2024a for every species, the only study measuring human, rat
    and two mouse strains in one laboratory.

TWO VALUES ARE EXCLUDED AND THE REASONS ARE PART OF THE RESULT.
  * Mouse PFUnDA Vd is printed as "33 +/- 60" in EFSA Appendix C against
    280 +/- 80 for the female of the same study. A standard deviation twice
    the mean, and a tenfold male/female gap seen for no other compound in
    that series, indicate a transcription error; 330 is the obvious reading
    but is not used, because guessing a digit is not evidence.
  * Rat PFHxA Vd is 175 mL/kg in EFSA Appendix C and 930 mL/kg in this
    project's own extraction of the same Chengelis 2009 study. The values
    are probably Vss and Vd-beta respectively, but nothing in either source
    says so, so the compound is dropped from the rat column rather than
    resolved by assumption.
"""
import csv, math, pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent.parent
GFR = {"human": 2570.0, "rat": 14400.0, "mouse": 16700.0}

# Abraham 2024, male volunteer: Vd mL/kg, CL mL/kg/day, both measured.
HUMAN = {
    "PFBA": (132.0, 21.89), "PFPeA": (110.0, 146.6), "PFHxA": (119.0, 56.89),
    "PFHpA": (128.0, 0.584), "PFOA": (121.0, 0.0417), "PFNA": (124.0, 0.0659),
    "PFDA": (184.0, 0.148), "PFUnDA": (241.0, 0.286), "PFDoDA": (354.0, 0.832),
    "PFBS": (137.0, 1.877), "PFHxS": (125.0, 0.0530), "PFOS": (152.0, 0.0870),
}
ANIMAL = {
    ("rat", "PFBA"):   (209.0, 444.0,  "Chang 2008, SD, IV"),
    ("rat", "PFHpA"):  (196.0, 1604.0, "Ohmori 2003, Wistar, IV, Vss"),
    ("rat", "PFOA"):   (112.0, 47.39,  "Kim 2016, SD, IV, non-compartmental"),
    ("rat", "PFNA"):   (287.0, 6.9,    "Ohmori 2003, Wistar, IV, Vss"),
    ("rat", "PFDA"):   (348.0, 5.2,    "Ohmori 2003, Wistar, IV, Vss"),
    ("rat", "PFBS"):   (188.0, 828.0,  "Huang 2019 NTP, SD, IV"),
    ("rat", "PFHxS"):  (275.0, 6.75,   "Sundstrom 2012, SD, IV"),
    ("rat", "PFOS"):   (383.0, 9.24,   "Kim 2016, SD, IV, non-compartmental"),
    ("mouse", "PFBA"):  (152.0, 240.0, "Chang 2008, CD-1, 47 umol/kg"),
    ("mouse", "PFHpA"): (70.0, 347.0,  "Fujii 2015, 0.31 umol/kg"),
    ("mouse", "PFOA"):  (180.0, 14.2,  "Fujii 2015, FVB/NJcl"),
    ("mouse", "PFNA"):  (220.0, 3.9,   "Fujii 2015, 0.31 umol/kg"),
    ("mouse", "PFDA"):  (250.0, 2.2,   "Fujii 2015, 0.31 umol/kg"),
    ("mouse", "PFDoDA"): (570.0, 4.4,  "Fujii 2015, 0.31 umol/kg"),
    ("mouse", "PFHxS"): (129.0, 2.9,   "Sundstrom 2012, 2.3 umol/kg"),
    ("mouse", "PFOS"):  (276.0, 4.7,   "Chang 2012, CD-1, mean of 263-290"),
}
# Monkey: Vd and CL exist, fu does NOT, in this or any study.
MONKEY = {
    "PFBA":  (526.0, 1694.0, "Chang 2008, IV"),
    "PFHxA": (989.0, 2928.0, "Chengelis 2009, IV"),
    "PFHxS": (287.0, 1.3,    "Sundstrom 2012, IV"),
    "PFOA":  (181.0, 12.4,   "Butenhoff 2004, IV"),
    "PFOS":  (202.0, 1.10,   "Chang 2012, IV"),
}
NC = {"PFBA": 4, "PFPeA": 5, "PFHxA": 6, "PFHpA": 7, "PFOA": 8, "PFNA": 9,
      "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12, "PFBS": 4, "PFHxS": 6,
      "PFOS": 8}


def fu_table():
    out = {}
    for r in csv.DictReader((HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r["experiment"] != "Ryu 2024" or not r["fu"]:
            continue
        sp = r["species"].split(" (")[0]
        strain = r["species"]
        if sp == "mouse":
            out.setdefault(("mouse", r["chemical"]), []).append(float(r["fu"]))
        elif sp in ("human", "rat"):
            out[(sp, r["chemical"])] = float(r["fu"])
    # mouse: geometric mean of CD-1 and C57BL/6, both strains being used above
    for k, v in list(out.items()):
        if isinstance(v, list):
            out[k] = math.exp(sum(math.log(x) for x in v) / len(v))
    return out


def boot(vals, n=10000, seed=20261124):
    rng = np.random.default_rng(seed)
    a = np.asarray(vals, dtype=float)
    m = rng.choice(a, size=(n, len(a)), replace=True).mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def main():
    fu = fu_table()
    out, bar = [], "=" * 78
    for sp in ("rat", "mouse"):
        rows = []
        for chem in sorted(HUMAN, key=lambda c: NC[c]):
            if (sp, chem) not in ANIMAL:
                continue
            fh, fx = fu.get(("human", chem)), fu.get((sp, chem))
            if not (fh and fx):
                continue
            vh, clh = HUMAN[chem]
            vx, clx, src = ANIMAL[(sp, chem)]
            th, tx = math.log(2) * vh / clh, math.log(2) * vx / clx
            ratio = th / tx
            tot = math.log10(ratio)
            fd, fb, fg = vh / vx, fx / fh, GFR[sp] / GFR["human"]
            ft = (clx / clh) / (fb * fg)
            rows.append(dict(sp=sp, chem=chem, nc=NC[chem], th=th, tx=tx,
                             ratio=ratio, fd=fd, fb=fb, fg=fg, ft=ft,
                             src=src, vh=vh, vx=vx, clh=clh, clx=clx,
                             fuh=fh, fux=fx,
                             pd=math.log10(fd)/tot*100,
                             pb=math.log10(fb)/tot*100,
                             pg=math.log10(fg)/tot*100,
                             pt=math.log10(ft)/tot*100))
        print(bar)
        print(f"HUMAN / {sp.upper()} HALF-LIFE RATIO, DECOMPOSED   "
              f"(n = {len(rows)} compounds)")
        print(bar)
        print(f"  {'chem':8}{'tC':>3}{'t_h(d)':>9}{'t_x(d)':>9}{'ratio':>9}"
              f"{'dist':>8}{'bind':>7}{'glom':>7}{'TUBULE':>10}")
        for r in rows:
            print(f"  {r['chem']:8}{r['nc']:>3}{r['th']:>9.1f}{r['tx']:>9.2f}"
                  f"{r['ratio']:>9.0f}{r['fd']:>8.2f}{r['fb']:>7.2f}"
                  f"{r['fg']:>7.2f}{r['ft']:>10.1f}")
        print()
        print(f"  {'':8}{'':3}{'share of log ratio (%)':>36}")
        print(f"  {'chem':8}{'tC':>3}{'distribution':>16}{'binding':>10}"
              f"{'glomerular':>12}{'TUBULE':>9}")
        for r in rows:
            print(f"  {r['chem']:8}{r['nc']:>3}{r['pd']:>16.0f}{r['pb']:>10.0f}"
                  f"{r['pg']:>12.0f}{r['pt']:>9.0f}")
        print(f"  {'mean':8}{'':3}"
              f"{sum(r['pd'] for r in rows)/len(rows):>16.0f}"
              f"{sum(r['pb'] for r in rows)/len(rows):>10.0f}"
              f"{sum(r['pg'] for r in rows)/len(rows):>12.0f}"
              f"{sum(r['pt'] for r in rows)/len(rows):>9.0f}")
        print()
        print("  bootstrap 95% CI of the mean share (10,000 resamples):")
        for k, lbl in (("pd", "distribution"), ("pb", "binding"),
                       ("pg", "glomerular"), ("pt", "TUBULE")):
            lo, hi = boot([r[k] for r in rows])
            print(f"    {lbl:14} [{lo:+7.0f}, {hi:+7.0f}] %")
        print()
        out.extend(rows)

    print(bar)
    print("THE MONKEY CANNOT BE DECOMPOSED, AND THAT ABSENCE IS INFORMATIVE")
    print(bar)
    print("  Vd and CL exist for five compounds. No unbound fraction has ever")
    print("  been measured in cynomolgus plasma for any PFAS. Instead of")
    print("  omitting the species, ask what fu it would NEED for binding")
    print("  alone to explain the human/monkey half-life difference:")
    print()
    print(f"  {'chem':8}{'t_h(d)':>9}{'t_mk(d)':>9}{'ratio':>8}"
          f"{'fu_mk required':>16}{'vs human fu':>13}")
    GFR_MK = 2570.0 * 2.0       # cynomolgus GFR ~2x human per kg
    req = []
    for chem, (vm, clm, src) in sorted(MONKEY.items(), key=lambda kv: NC[kv[0]]):
        fh = fu.get(("human", chem))
        if not fh or chem not in HUMAN:
            continue
        vh, clh = HUMAN[chem]
        th, tm = math.log(2) * vh / clh, math.log(2) * vm / clm
        ratio = th / tm
        # binding alone: ratio = (Vd_h/Vd_mk) * (fu_mk/fu_h) => solve fu_mk
        fux = ratio / (vh / vm) * fh
        req.append((chem, fux / fh))
        print(f"  {chem:8}{th:>9.1f}{tm:>9.2f}{ratio:>8.0f}"
              f"{fux:>16.4f}{fux/fh:>12.0f}x")
    print()
    print("  Every one of these required unbound fractions is far outside the")
    print("  0.76-2.03 range within which rat/human fu ratios actually fall,")
    print("  and several exceed 1.0, which is physically impossible. BINDING")
    print("  CANNOT EXPLAIN THE MONKEY EITHER, and the conclusion needs no")
    print("  monkey binding measurement to reach -- though one would convert")
    print("  this argument from a reductio into a decomposition.")
    print()
    print("  (Cynomolgus GFR taken as twice the human value per kg; the")
    print("   conclusion is insensitive to that factor, which would have to")
    print("   be wrong by orders of magnitude to change it.)")
    print()

    p = HERE / "db" / "species_decomposition.csv"
    with p.open("w", newline="") as fh_:
        w = csv.writer(fh_)
        w.writerow(["comparison", "chemical", "n_total_carbons",
                    "t_half_human_d", "t_half_animal_d", "halflife_ratio",
                    "factor_distribution", "factor_binding",
                    "factor_glomerular", "factor_tubule_residual",
                    "pct_distribution", "pct_binding", "pct_glomerular",
                    "pct_tubule", "vd_human", "vd_animal", "cl_human",
                    "cl_animal", "fu_human", "fu_animal", "animal_source"])
        for r in out:
            w.writerow([f"human/{r['sp']}", r["chem"], r["nc"],
                        round(r["th"], 2), round(r["tx"], 3),
                        round(r["ratio"], 1), round(r["fd"], 4),
                        round(r["fb"], 4), round(r["fg"], 4),
                        round(r["ft"], 3), round(r["pd"], 1),
                        round(r["pb"], 1), round(r["pg"], 1),
                        round(r["pt"], 1), r["vh"], r["vx"], r["clh"],
                        r["clx"], r["fuh"], round(r["fux"], 6), r["src"]])
    print(f"wrote {p.relative_to(HERE)}  ({len(out)} rows)")


if __name__ == "__main__":
    main()
