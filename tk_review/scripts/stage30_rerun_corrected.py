#!/usr/bin/env python3
"""Human models re-run on the corrected renal T with Smeltz f_u (Stage 30).

    python3 -I scripts/stage30_rerun_corrected.py

Writes db/stage30_models.csv, db/stage30_summary.csv,
db/stage30_vs_stage28.csv, db/stage30_T_smeltz.csv. Idempotent.

TWO CHANGES FROM STAGE 28, BOTH ESTABLISHED IN STAGE 29.

1. T IS NOW RENAL-ONLY. Abraham 2024 Table 4 reports faecal clearance for
   the long-chain PFCAs at 17% to 61% of total, and urinary clearance as
   "n.d." throughout. T computed from total clearance credited the kidney
   with the faecal route, so it is replaced by

       T = CL_total * (1 - CL_fec/CL_tot) / (f_u * GFR)

   which is a renal-only UPPER BOUND, not a point estimate.

2. f_u IS SMELTZ 2023, NOT RYU 2024. Smeltz has the best worst-case error
   on the common compound set, the widest coverage (8 of 9), the largest
   series (67 PFAS), and is the f_u Abraham itself used in its Table 4 --
   so half-life, V_d and f_u become one published pairing.

   SMELTZ HAS NO PFDoDA. That compound takes Ryu's value, so the f_u
   column is NOT from a single source and every table below says so.

T IS DERIVED WITH THE SAME f_u THE MODEL DIVIDES BY. Mixing one study's
f_u into T and another's into the prediction was a bug in Stage 29's
first draft; it is not repeated here.

WHY THE RESULTS ARE STRATIFIED BY CHAIN LENGTH. Under Smeltz, T spans
84-fold across the eight compounds it covers but only 3.4-fold across the
five with eight or more carbons. The spread is almost entirely in the
short and mid chains -- Smeltz's PFHpA f_u of 0.0004 is the lowest value
any study reports for it, 11-fold below Ryu's, and it implies a tubular
factor of 0.57, i.e. almost no reabsorption. A single-constant T model
therefore cannot be judged on the pooled set: it is being asked to cover
two different regimes, and the stratification shows which one it fits.

THE QSAR IS UNCHANGED BY BOTH CORRECTIONS, because it never used f_u or
T. It is carried here as a fixed benchmark.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
LN2 = math.log(2.0)
GFR_HUMAN = 2570.0
ALBUMIN_M = 6.0e-4
SULFONATES = {"PFBS", "PFHxS", "PFOS"}

# Abraham 2024 Table 4: faecal / total clearance, both mL/min as published.
FECAL = {"PFOA": 0.0004 / 0.0024, "PFNA": 0.0022 / 0.0038,
         "PFDA": 0.0051 / 0.0084, "PFDoDA": 0.0286 / 0.0474}


def sig(x, n=6):
    if x is None:
        return None
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


def gmean(vs):
    return 10 ** (sum(math.log10(v) for v in vs) / len(vs))


def load_human():
    out = {}
    for r in csv.DictReader((HERE / "db" / "species_decomposition.csv").open()):
        c = r["chemical"]
        if c in out:
            continue
        out[c] = dict(t_half=float(r["t_half_human_d"]),
                      vd=float(r["vd_human"]), cl=float(r["cl_human"]),
                      fu_ryu=float(r["fu_human"]),
                      nc=int(r["n_total_carbons"]))
    return out


def load_smeltz():
    out = {}
    for r in csv.DictReader((HERE / "db" / "human_fu_comparison.csv").open()):
        if r.get("smeltz_UC"):
            out[r["chemical"]] = float(r["smeltz_UC"])
    return out


def load_albumin_ka():
    acc = {}
    for r in csv.DictReader((HERE / "db" / "protein_binding.csv").open()):
        if ("ALBUMIN" not in r["protein"].upper()
                or r["protein_species"] != "human"):
            continue
        u_ka, u_kd = (r["ka_units"] or "").lower(), (r["kd_units"] or "").lower()
        ka = None
        if r["ka"]:
            try:
                v = float(r["ka"])
                ka = 10 ** v if u_ka.startswith("log") else (
                    v if u_ka.startswith("m^-1") else None)
            except ValueError:
                ka = None
        elif r["kd"]:
            try:
                v = float(r["kd"])
                f = {"um": 1e-6, "mm": 1e-3, "m": 1.0}.get(u_kd)
                ka = 1.0 / (v * f) if f and v else None
            except ValueError:
                ka = None
        if ka and ka > 0:
            acc.setdefault(r["chemical"], []).append(ka)
    return {c: gmean(v) for c, v in acc.items()}


def ols(rows, y):
    p = len(rows[0])
    aug = [[sum(r[i] * r[j] for r in rows) for j in range(p)]
           + [sum(rows[k][i] * y[k] for k in range(len(rows)))]
           for i in range(p)]
    for i in range(p):
        piv = max(range(i, p), key=lambda k: abs(aug[k][i]))
        if abs(aug[piv][i]) < 1e-12:
            return None
        aug[i], aug[piv] = aug[piv], aug[i]
        for k in range(p):
            if k == i:
                continue
            f = aug[k][i] / aug[i][i]
            for j in range(i, p + 1):
                aug[k][j] -= f * aug[i][j]
    return [aug[i][p] / aug[i][i] for i in range(p)]


def main():
    hum, smeltz, ka = load_human(), load_smeltz(), load_albumin_ka()
    bar = "=" * 78
    chems = sorted(hum)

    # f_u actually used, with its provenance recorded per compound
    fu, prov = {}, {}
    for c in chems:
        if c in smeltz:
            fu[c], prov[c] = smeltz[c], "Smeltz 2023 (UC)"
        else:
            fu[c], prov[c] = hum[c]["fu_ryu"], "Ryu 2024 (Smeltz has none)"

    def T(c):
        """Renal-only upper-bound tubular factor, self-consistent with fu."""
        return hum[c]["cl"] * (1 - FECAL.get(c, 0.0)) / (fu[c] * GFR_HUMAN)

    print(bar)
    print("INPUTS AS CORRECTED")
    print(bar)
    print(f"  {'cmpd':9}{'nC':>3}{'f_u':>10}{'fecal':>8}{'T renal':>10}"
          f"   f_u source")
    trows = []
    for c in chems:
        f = FECAL.get(c)
        trows.append(dict(chemical=c, n_total_carbons=hum[c]["nc"],
                          fu_used=sig(fu[c]), fu_source=prov[c],
                          fecal_pct_of_total=sig(f * 100, 4) if f else "",
                          T_renal_upper_bound=sig(T(c), 4)))
        print(f"  {c:9}{hum[c]['nc']:>3}{fu[c]:>10.5f}"
              f"{(f'{f*100:.0f}%' if f else '-'):>8}{T(c):>10.5f}   {prov[c]}")
    print()
    mixed = [c for c in chems if "Ryu" in prov[c]]
    print(f"  f_u IS NOT FROM A SINGLE SOURCE: {', '.join(mixed)} takes Ryu's")
    print("  value because Smeltz does not report it. Stated here and in")
    print("  db/stage30_T_smeltz.csv rather than presented as uniform.")
    print()
    print("  ALSO ON THE RECORD, from db/human_fu_recommended.csv: Smeltz")
    print("  gives PFBS as 0.0046 in its Table 2 and 0.0128 in its Table 3")
    print("  for the SAME assay, a 2.8-fold internal discrepancy. The lower")
    print("  value is used, as the comparison table carries it.")
    print()

    # strata
    # PFDoDA is held out of one stratum on Stage 29's finding, not on fit:
    # its urinary clearance was not detectable, so its renal T is an upper
    # bound of undetermined value rather than a measurement.
    strata = [
        ("all 9 (f_u mixed)", chems),
        ("Smeltz coverage, 8", [c for c in chems if c in smeltz]),
        ("long chain, nC>=8", [c for c in chems if hum[c]["nc"] >= 8]),
        ("long chain excl PFDoDA", [c for c in chems
                                    if hum[c]["nc"] >= 8 and c != "PFDoDA"]),
        ("short/mid, nC<=7", [c for c in chems if hum[c]["nc"] <= 7]),
    ]
    print(bar)
    print("WHY STRATIFY: THE TUBULAR FACTOR IS NOT ONE POPULATION")
    print(bar)
    for name, cs in strata:
        ts = [T(c) for c in cs]
        print(f"  {name:24} n={len(cs)}  T {min(ts):.5f} to {max(ts):.5f}"
              f"  = {max(ts)/min(ts):>5.1f}-fold")
    print()
    spans = {name: max(T(c) for c in cs) / min(T(c) for c in cs)
             for name, cs in strata}
    tight = min(spans, key=lambda k: spans[k])
    print("  A model that fits ONE constant T can only work where T is")
    print(f"  actually near-constant. The TIGHTEST stratum is \"{tight}\"")
    print(f"  at {spans[tight]:.1f}-fold; the widest is "
          f"{max(spans, key=lambda k: spans[k])} at "
          f"{max(spans.values()):.0f}-fold.")
    print()
    pf = T("PFDoDA")
    rest = [T(c) for c in chems if c != "PFDoDA"]
    hi = max(rest)
    hic = max((c for c in chems if c != "PFDoDA"), key=lambda c: T(c))
    print(f"  ONE COMPOUND DRIVES MOST OF THE SPREAD. PFDoDA's T is {pf:.3f},")
    print(f"  {pf/hi:.1f}x above the next highest ({hic}) and {pf/min(rest):.0f}x "
          f"above the lowest. Stage 29")
    print("  established that it is an UPPER BOUND on an undetectable urinary")
    print("  clearance -- not a measurement. Dropping it takes the nC>=8")
    print(f"  stratum from {spans['long chain, nC>=8']:.0f}-fold to "
          f"{spans['long chain excl PFDoDA']:.1f}-fold. That is why it has its own")
    print("  stratum, and the reason is Stage 29's verdict, not the fit.")
    print()

    # ---------------------------------------------------------- models
    rows, summary = [], []

    def score(model, params, stratum, cs, preds, note=""):
        folds = []
        for c in cs:
            p = preds.get(c)
            if p is None or p <= 0:
                continue
            obs = hum[c]["t_half"]
            fold = max(p / obs, obs / p)
            folds.append(fold)
            rows.append(dict(
                model=model, stratum=stratum, chemical=c,
                n_total_carbons=hum[c]["nc"],
                t_half_measured_d=sig(obs), t_half_predicted_d=sig(p),
                fold_error=sig(fold, 4),
                direction="over" if p > obs else "under"))
        if not folds:
            return
        folds.sort()
        summary.append(dict(
            model=model, stratum=stratum, n=len(folds), n_params=params,
            median_fold=sig(folds[len(folds) // 2], 4),
            worst_fold=sig(folds[-1], 4), note=note))

    def qsar_loo(form, cs):
        def design(c):
            x = [1.0, float(hum[c]["nc"])]
            if "nC^2" in form:
                x.append(float(hum[c]["nc"]) ** 2)
            if "head" in form:
                x.append(1.0 if c in SULFONATES else 0.0)
            return x
        p = len(design(cs[0]))
        if len(cs) <= p + 1:
            return None, p
        out = {}
        for c in cs:
            tr = [x for x in cs if x != c]
            b = ols([design(x) for x in tr],
                    [math.log10(hum[x]["t_half"]) for x in tr])
            out[c] = 10 ** sum(a * d for a, d in zip(b, design(c))) if b else None
        return out, p

    for name, cs in strata:
        # M1a  T = 1
        score("M1a: T=1 (filtration alone)", 0, name, cs,
              {c: LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN) for c in cs},
              "parameter-free")
        # M1b  T = fitted constant, LOO
        pb = {}
        for c in cs:
            tr = [x for x in cs if x != c]
            if len(tr) < 2:
                continue
            pb[c] = LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN
                                          * gmean([T(x) for x in tr]))
        score("M1b: T = fitted constant", 1, name, cs, pb, "leave-one-out")
        # M1c  log T = a + b nC, LOO
        pc = {}
        for c in cs:
            tr = [x for x in cs if x != c]
            if len(tr) < 3:
                continue
            b = ols([[1.0, float(hum[x]["nc"])] for x in tr],
                    [math.log10(T(x)) for x in tr])
            if not b:
                continue
            That = 10 ** (b[0] + b[1] * hum[c]["nc"])
            pc[c] = LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN * That)
        score("M1c: log T = a + b*nC", 2, name, cs, pc, "leave-one-out")
        # M1d  f_u from albumin K_a, T fitted, LOO
        sh = [c for c in cs if c in ka]
        pd = {}
        for c in sh:
            tr = [x for x in sh if x != c]
            if len(tr) < 3:
                continue
            n = 2
            Ts = [hum[x]["cl"] * (1 - FECAL.get(x, 0.0))
                  * (1 + n * ALBUMIN_M * ka[x]) / GFR_HUMAN for x in tr]
            fub = 1.0 / (1 + n * ALBUMIN_M * ka[c])
            pd[c] = LN2 * hum[c]["vd"] / (fub * GFR_HUMAN * gmean(Ts))
        score("M1d: f_u from albumin K_a (n=2), T fitted", 1, name, sh, pd,
              "leave-one-out; albumin replaces measured f_u")
        # QSAR, unaffected by either correction
        for form in ("log t½ ~ nC", "log t½ ~ nC + nC^2",
                     "log t½ ~ nC + nC^2 + head"):
            pq, k = qsar_loo(form, cs)
            if pq:
                score("QSAR " + form, k, name, cs, pq,
                      "leave-one-out; uses no f_u and no T")

    for name, cs in strata:
        print(bar)
        print(f"{name.upper()}   (n={len(cs)}: {', '.join(cs)})")
        print(bar)
        print(f"  {'model':44}{'k':>3}{'median':>10}{'worst':>10}")
        for s in sorted((s for s in summary if s["stratum"] == name),
                        key=lambda s: s["median_fold"]):
            print(f"  {s['model']:44}{s['n_params']:>3}"
                  f"{s['median_fold']:>9.1f}x{s['worst_fold']:>9.1f}x")
        print()

    # ---------------------------------------------------------- vs stage 28
    print(bar)
    print("WHAT THE TWO CORRECTIONS BOUGHT, AGAINST STAGE 28")
    print(bar)
    old = {}
    for r in csv.DictReader((HERE / "db" / "human_tk_summary.csv").open()):
        old[r["model"].split(",")[0]] = (float(r["median_fold"]),
                                         float(r["worst_fold"]))
    for r in csv.DictReader((HERE / "db" / "human_qsar_summary.csv").open()):
        old["QSAR " + r["form"]] = (float(r["loo_median_fold"]),
                                    float(r["loo_worst_fold"]))
    print("  Stage 28: Ryu f_u, T from total clearance, all 9 compounds.")
    print("  Stage 30: Smeltz f_u, renal-only T, all 9 compounds.")
    print()
    print(f"  {'model':40}{'S28 med':>9}{'S30 med':>9}"
          f"{'S28 wst':>9}{'S30 wst':>9}")
    cmp_rows = []
    for s in sorted((s for s in summary if s["stratum"] == "all 9 (f_u mixed)"),
                    key=lambda s: s["model"]):
        key = s["model"].split(",")[0]
        o = old.get(key) or old.get(s["model"])
        if not o:
            continue
        cmp_rows.append(dict(
            model=s["model"], stage28_median=sig(o[0], 4),
            stage30_median=s["median_fold"], stage28_worst=sig(o[1], 4),
            stage30_worst=s["worst_fold"],
            median_change=sig(s["median_fold"] / o[0], 3),
            worst_change=sig(s["worst_fold"] / o[1], 3)))
        print(f"  {s['model'][:40]:40}{o[0]:>8.1f}x{s['median_fold']:>8.1f}x"
              f"{o[1]:>8.1f}x{s['worst_fold']:>8.1f}x")
    print()
    q = [r for r in cmp_rows if r["model"].startswith("QSAR")]
    if q and all(abs(r["median_change"] - 1) < 0.01 for r in q):
        print("  THE QSAR ROWS ARE IDENTICAL, AS THEY MUST BE. It uses neither")
        print("  f_u nor T, so a change to either cannot move it. That the")
        print("  numbers match exactly is a check on the rerun, not a result.")
    print()

    print(bar)
    print("WHAT TO TAKE FROM THIS")
    print(bar)
    LC, SC = "long chain excl PFDoDA", "short/mid, nC<=7"

    def best(stratum):
        return min((s for s in summary if s["stratum"] == stratum),
                   key=lambda s: s["median_fold"])

    def find(stratum, model):
        return next((s for s in summary if s["stratum"] == stratum
                     and s["model"] == model), None)

    bl, bs = best(LC), best(SC)
    print("  1. THE RANKING REVERSES BETWEEN THE TWO REGIMES, AND THAT IS")
    print("     THE RESULT WORTH REPORTING.")
    print(f"     {'stratum':26}{'best model':34}{'median':>9}{'worst':>9}")
    for st, b in ((LC, bl), (SC, bs)):
        print(f"     {st:26}{b['model'][:34]:34}"
              f"{b['median_fold']:>8.1f}x{b['worst_fold']:>8.1f}x")
    cross = [(find(SC, bl["model"]), bl["model"], SC),
             (find(LC, bs["model"]), bs["model"], LC)]
    print()
    print("     Each winner carried into the OTHER stratum:")
    for row, model, st in cross:
        if row:
            print(f"     {model[:34]:34} on {st:24}"
                  f"{row['median_fold']:>8.1f}x{row['worst_fold']:>8.1f}x")
    print()
    print("     So these are not competing models with one winner. The QSAR")
    print("     works where half-life is a smooth function of chain length")
    print("     (2011 -> 1305 -> 862 -> 295 d across PFOA to PFDoDA); the")
    print("     mechanistic model works where it is not, but T is near-")
    print("     constant. PFBA 4.2 d, PFBS 50.6 d, PFHxS 1635 d, PFHpA 152 d")
    print("     is no function of chain length at all, and the QSAR fails")
    print(f"     there ({bs['model'] != 'QSAR log t½ ~ nC' and find(SC, 'QSAR log t½ ~ nC') and find(SC, 'QSAR log t½ ~ nC')['median_fold'] or 0:.1f}x median) while a fitted T does not.")
    print()

    m1b_all = find("all 9 (f_u mixed)", "M1b: T = fitted constant")
    print("  2. THE CORRECTIONS MADE THE POOLED NUMBERS WORSE, NOT BETTER.")
    print(f"     M1b on all 9 goes 2.1x (Stage 28) -> "
          f"{m1b_all['median_fold']:.1f}x here; M1a 26.9x -> "
          f"{find('all 9 (f_u mixed)', 'M1a: T=1 (filtration alone)')['median_fold']:.1f}x.")
    print("     This is not a reason to revert. Stage 28's lower numbers came")
    print("     from a T that credited the kidney with faecal elimination and")
    print("     from Ryu f_u, whose T happens to be more uniform across")
    print("     compounds. A model is not improved by an input that is wrong")
    print("     in a convenient direction. The corrected inputs are right and")
    print("     the error they expose is real.")
    print()

    print("  3. SMELTZ'S ADVANTAGE IS CONFINED TO THE LONG CHAINS. Its PFHpA")
    print("     f_u of 0.0004 -- the lowest any study reports, 11-fold below")
    print(f"     Ryu's -- implies T = {T('PFHpA'):.2f}, i.e. almost no net")
    print("     reabsorption of a compound with a 152-day half-life. Either")
    print("     that f_u is too low or PFHpA is cleared by a route this model")
    print("     does not contain. It is the single largest contributor to the")
    print("     short-chain failure and should be flagged in the paper.")
    print()

    thin = [s for s in summary
            if s["stratum"] in (LC, SC) and s["n"] - s["n_params"] <= 2]
    print("  4. EVERY NUMBER IN THE TWO STRATA RESTS ON FOUR COMPOUNDS.")
    print(f"     {len(thin)} of the stratified fits have n - k <= 2, i.e. at most two")
    print("     residual degrees of freedom under leave-one-out. The reversal")
    print("     in point 1 is a real feature of these data, but its SIZE is")
    print("     not estimable from four points, and the paper should say so")
    print("     rather than quote 1.5x as an accuracy.")
    print()
    print("  5. AND THE WHOLE HUMAN ARM IS n=1. Every half-life and every V_d")
    print("     is Abraham's single volunteer. No fold error here is an")
    print("     estimate of population accuracy.")
    print()

    out = [
        ("stage30_models.csv", rows,
         ["model", "stratum", "chemical", "n_total_carbons",
          "t_half_measured_d", "t_half_predicted_d", "fold_error",
          "direction"]),
        ("stage30_summary.csv", summary,
         ["model", "stratum", "n", "n_params", "median_fold", "worst_fold",
          "note"]),
        ("stage30_vs_stage28.csv", cmp_rows,
         ["model", "stage28_median", "stage30_median", "median_change",
          "stage28_worst", "stage30_worst", "worst_change"]),
        ("stage30_T_smeltz.csv", trows,
         ["chemical", "n_total_carbons", "fu_used", "fu_source",
          "fecal_pct_of_total", "T_renal_upper_bound"]),
    ]
    for name, rr, cols in out:
        p = HERE / "db" / name
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in sorted(rr, key=lambda d: [str(d.get(x, "")) for x in cols]):
                w.writerow({x: r.get(x, "") for x in cols})
        print(f"wrote db/{name}  ({len(rr)} rows)")


if __name__ == "__main__":
    main()
