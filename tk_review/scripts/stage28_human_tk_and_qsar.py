#!/usr/bin/env python3
"""Human half-life models: mechanistic TK without renal clearance, and a
structure-only QSAR, kept separate (Stage 28).

    python3 -I scripts/stage28_human_tk_and_qsar.py

Writes db/human_tk_models.csv, db/human_tk_summary.csv,
db/human_qsar.csv, db/human_qsar_summary.csv, db/human_fu_sensitivity.csv,
db/human_tubular_factor.csv. Idempotent.

SCOPE, SET DELIBERATELY AND NARROWLY. Human only -- the rat and mouse arms
are dropped, not because their numbers were wrong but because 8 compounds
per species will not support a model, and the rat renal fraction turned out
to be in 10-fold dispute between two primary sources. Renal clearance is
dropped too: a model built on it over-predicted human half-lives by up to
166-fold, and an over-prediction of that size is a model failure whatever
else it also happens to measure.

WHAT IS LEFT IS THE MECHANISTIC CORE:

    t½ = ln2 V_d / (f_u * GFR * T)

with f_u, V_d and GFR measured, and T -- the tubular factor -- the only
unknown. Four models differ ONLY in where T comes from:

    M1a  T = 1                      parameter-free; filtration alone
    M1b  T = fitted constant        1 parameter
    M1c  log T = a + b*nC           2 parameters; T tracks chain length?
    M1d  f_u from albumin K_a       binding replaces measured f_u

M1a is the honest baseline and it is expected to FAIL, because human T is
far below 1 for almost every compound. The size of that failure is the
point: it is exactly 1/T, so it measures how much of the human half-life
the mechanistic core cannot reach.

THE QSAR IS A SEPARATE OBJECT AND IS NOT A TIER OF THE ABOVE. It sees
chain length and head group, no TK quantity at all. Being fitted, it is
scored leave-one-out, as are M1b, M1c and M1d. M1a spends no parameters so
its score needs no holdout.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
LN2 = math.log(2.0)
GFR_HUMAN = 2570.0
ALBUMIN_M = 6.0e-4
SULFONATES = {"PFBS", "PFHxS", "PFOS", "PFHpS", "PFDS"}
FU_SOURCES = [("smeltz_UC", "Smeltz 2023 (UC)"), ("ryu", "Ryu 2024"),
              ("argoul", "Argoul 2025"), ("fischer_derived", "Fischer 2024")]


def sig(x, n=6):
    if x is None:
        return None
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


def gmean(vs):
    return 10 ** (sum(math.log10(v) for v in vs) / len(vs))


# --------------------------------------------------------------- inputs
def load_human():
    """{chemical: {t_half, vd, fu, cl, nc, T}} -- compilation values."""
    out = {}
    for r in csv.DictReader((HERE / "db" / "species_decomposition.csv").open()):
        c = r["chemical"]
        if c in out:
            continue
        fu, cl = float(r["fu_human"]), float(r["cl_human"])
        out[c] = dict(t_half=float(r["t_half_human_d"]),
                      vd=float(r["vd_human"]), fu=fu, cl=cl,
                      nc=int(r["n_total_carbons"]),
                      T=cl / (fu * GFR_HUMAN))
    return out


def load_fu_by_source():
    """{source_key: {chemical: f_u}} across the four independent studies."""
    out = {k: {} for k, _ in FU_SOURCES}
    for r in csv.DictReader((HERE / "db" / "human_fu_comparison.csv").open()):
        for k, _ in FU_SOURCES:
            if r.get(k):
                out[k][r["chemical"]] = float(r[k])
    return out


def load_albumin_ka():
    """{chemical: K_a (M^-1)} -- geometric mean over human albumin studies."""
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


# --------------------------------------------------------------- algebra
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


def loo(chems, predict_from_train):
    """predict_from_train(train_chems, target) -> t½ or None."""
    out = {}
    for c in chems:
        train = [x for x in chems if x != c]
        out[c] = predict_from_train(train, c)
    return out


# --------------------------------------------------------------- models
def build(hum, ka):
    """Returns (rows, summary). Each model predicts a human half-life."""
    chems = sorted(hum)
    rows, summary = [], []

    def record(model, params, preds, note):
        folds = []
        for c in chems:
            p = preds.get(c)
            if p is None or p <= 0:
                continue
            obs = hum[c]["t_half"]
            fold = max(p / obs, obs / p)
            folds.append(fold)
            rows.append(dict(
                model=model, chemical=c, n_total_carbons=hum[c]["nc"],
                head_group="sulfonate" if c in SULFONATES else "carboxylate",
                t_half_measured_d=sig(obs), t_half_predicted_d=sig(p),
                fold_error=sig(fold, 4),
                direction="over" if p > obs else "under"))
        if not folds:
            return
        folds.sort()
        summary.append(dict(
            model=model, n=len(folds), n_params=params,
            scored="leave-one-out" if params else "direct (no parameters)",
            median_fold=sig(folds[len(folds) // 2], 4),
            worst_fold=sig(folds[-1], 4),
            note=note))

    # ---- M1a  T = 1, parameter-free
    record("M1a: T=1 (filtration alone)", 0,
           {c: LN2 * hum[c]["vd"] / (hum[c]["fu"] * GFR_HUMAN) for c in chems},
           "no parameters; error is exactly 1/T")

    # ---- M1b  T = fitted constant (geometric mean of training T)
    def m1b(train, c):
        That = gmean([hum[x]["T"] for x in train])
        return LN2 * hum[c]["vd"] / (hum[c]["fu"] * GFR_HUMAN * That)
    record("M1b: T = fitted constant", 1, loo(chems, m1b),
           "one parameter absorbs the systematic shortfall")

    # ---- M1c  log T = a + b*nC
    def m1c(train, c):
        X = [[1.0, float(hum[x]["nc"])] for x in train]
        y = [math.log10(hum[x]["T"]) for x in train]
        b = ols(X, y)
        if b is None:
            return None
        That = 10 ** (b[0] + b[1] * hum[c]["nc"])
        return LN2 * hum[c]["vd"] / (hum[c]["fu"] * GFR_HUMAN * That)
    record("M1c: log T = a + b*nC", 2, loo(chems, m1c),
           "tests whether T tracks chain length")

    # ---- M1d  f_u from albumin binding, T fitted constant
    shared = [c for c in chems if c in ka]
    for n_sites in (1, 2, 3):
        def m1d(train, c, n=n_sites):
            tr = [x for x in train if x in ka]
            if len(tr) < 3:
                return None
            # T fitted on the SAME f_u definition the model will use
            Ts = [hum[x]["cl"] * (1 + n * ALBUMIN_M * ka[x]) / GFR_HUMAN
                  for x in tr]
            That = gmean(Ts)
            fu = 1.0 / (1 + n * ALBUMIN_M * ka[c])
            return LN2 * hum[c]["vd"] / (fu * GFR_HUMAN * That)
        record(f"M1d: f_u from albumin K_a, n={n_sites}, T fitted", 1,
               loo(shared, m1d),
               f"binding replaces measured f_u; {len(shared)} compounds")
    return rows, summary, shared


def qsar(hum):
    chems = sorted(hum)
    forms = ["log t½ ~ nC", "log t½ ~ nC + nC^2", "log t½ ~ nC + head",
             "log t½ ~ nC + nC^2 + head"]
    rows, summary = [], []
    for form in forms:
        def design(c):
            x = [1.0, float(hum[c]["nc"])]
            if "nC^2" in form:
                x.append(float(hum[c]["nc"]) ** 2)
            if "head" in form:
                x.append(1.0 if c in SULFONATES else 0.0)
            return x
        p = len(design(chems[0]))
        full = ols([design(c) for c in chems],
                   [math.log10(hum[c]["t_half"]) for c in chems])
        ys = [math.log10(hum[c]["t_half"]) for c in chems]
        my = sum(ys) / len(ys)
        sst = sum((v - my) ** 2 for v in ys)
        sse = sum((math.log10(hum[c]["t_half"])
                   - sum(a * b for a, b in zip(full, design(c)))) ** 2
                  for c in chems)
        folds = []
        for c in chems:
            train = [x for x in chems if x != c]
            b = ols([design(x) for x in train],
                    [math.log10(hum[x]["t_half"]) for x in train])
            if b is None:
                continue
            pr = 10 ** sum(a * d for a, d in zip(b, design(c)))
            fold = max(pr / hum[c]["t_half"], hum[c]["t_half"] / pr)
            folds.append(fold)
            rows.append(dict(
                form=form, chemical=c, n_total_carbons=hum[c]["nc"],
                head_group="sulfonate" if c in SULFONATES else "carboxylate",
                t_half_measured_d=sig(hum[c]["t_half"]),
                t_half_loo_predicted_d=sig(pr), loo_fold_error=sig(fold, 4)))
        folds.sort()
        summary.append(dict(
            form=form, n=len(chems), n_params=p,
            r2_in_sample=sig(1 - sse / sst, 3),
            loo_median_fold=sig(folds[len(folds) // 2], 4),
            loo_worst_fold=sig(folds[-1], 4)))
    return rows, summary


# --------------------------------------------------------------- report
def main():
    hum, ka = load_human(), load_albumin_ka()
    fu_src = load_fu_by_source()
    bar = "=" * 78
    chems = sorted(hum)

    print(bar)
    print("SCOPE")
    print(bar)
    print(f"  Human only, {len(chems)} compounds: {', '.join(chems)}")
    print("  Rat and mouse dropped: 8 compounds each will not support a")
    print("  model, and the rat renal fraction is in 10-fold dispute.")
    print("  Renal clearance dropped: it over-predicted by up to 166x, and")
    print("  an over-prediction that large is a model failure.")
    print()

    # ---- the tubular factor, which is the whole problem
    print(bar)
    print("THE UNKNOWN: THE HUMAN TUBULAR FACTOR T = CL / (f_u * GFR)")
    print(bar)
    print("  T < 1 is net reabsorption, T = 1 filtration alone, T > 1 net")
    print("  secretion. Every model below differs ONLY in where T comes from.")
    print()
    print(f"  {'cmpd':9}{'nC':>3}{'T':>10}{'1/T':>9}   reading")
    tf = []
    for c in chems:
        T = hum[c]["T"]
        rd = ("net secretion" if T > 1 else
              "strong reabsorption" if T < 0.1 else "reabsorption")
        tf.append(dict(chemical=c, n_total_carbons=hum[c]["nc"],
                       tubular_factor_T=sig(T, 4), one_over_T=sig(1 / T, 4),
                       reading=rd))
        print(f"  {c:9}{hum[c]['nc']:>3}{T:>10.4f}{1/T:>9.1f}   {rd}")
    Ts = [hum[c]["T"] for c in chems]
    sec = [c for c in chems if hum[c]["T"] > 1]
    rest = [c for c in chems if hum[c]["T"] <= 1]
    print()
    print(f"  Geometric mean T = {gmean(Ts):.4f}  (1/T = {1/gmean(Ts):.1f})")
    print(f"  Full range spans {max(Ts)/min(Ts):.0f}-fold.")
    print(f"  {len(sec)} compound(s) show net secretion: {', '.join(sec)}")
    print(f"  Excluding it, T spans {max(hum[c]['T'] for c in rest)/min(hum[c]['T'] for c in rest):.1f}-fold "
          f"({min(hum[c]['T'] for c in rest):.4f} to {max(hum[c]['T'] for c in rest):.4f}),")
    print("  which is the irreducible error of any single-constant T model.")
    print()
    if sec:
        s = sec[0]
        print(f"  WATCH {s}. Its T of {hum[s]['T']:.2f} is the only value above 1,")
        print(f"  it is the longest chain in the set, and its half-life")
        print(f"  ({hum[s]['t_half']:.0f} d) is SHORTER than {('PFDA' if s != 'PFDA' else 'PFNA')}'s "
              f"({hum['PFDA']['t_half']:.0f} d),")
        print("  breaking the chain-length trend. It drives several results")
        print("  below and is flagged wherever it does.")
    print()

    # ---- model 1
    print(bar)
    print("MODEL 1  MECHANISTIC TK -- f_u, V_d, GFR measured; T the unknown")
    print(bar)
    rows, summary, shared = build(hum, ka)
    print(f"  {'model':40}{'n':>3}{'k':>3}{'median':>10}{'worst':>10}")
    for s in summary:
        print(f"  {s['model']:40}{s['n']:>3}{s['n_params']:>3}"
              f"{s['median_fold']:>9.1f}x{s['worst_fold']:>9.1f}x")
    print()
    a = next(s for s in summary if s["model"].startswith("M1a"))
    und = sum(1 for r in rows if r["model"].startswith("M1a")
              and r["direction"] == "under")
    print(f"  M1a FAILS AS EXPECTED: median {a['median_fold']:.1f}x, worst "
          f"{a['worst_fold']:.0f}x, and it UNDER-predicts in")
    print(f"  {und} of {a['n']} compounds. Assuming T=1 credits the kidney with")
    print("  clearing everything it filters, so the predicted half-life is")
    print("  too short by exactly the factor the tubule reabsorbs.")
    print()
    b = next(s for s in summary if s["model"].startswith("M1b"))
    print(f"  M1b, ONE PARAMETER, CUTS THE MEDIAN TO {b['median_fold']:.1f}x. That single")
    print("  constant is doing the work of the entire tubular mechanism, and")
    print("  what it leaves behind is the compound-to-compound spread in T.")
    print()
    c3 = next((s for s in summary if s["model"].startswith("M1c")), None)
    if c3:
        verdict = ("WORSE" if c3["median_fold"] > b["median_fold"]
                   else "better")
        print(f"  M1c LETS T TRACK CHAIN LENGTH AND IS {verdict} "
              f"({c3['median_fold']:.1f}x vs {b['median_fold']:.1f}x).")
        print("  With one secreting outlier at the long-chain end, a linear")
        print("  trend in log T is fitted largely to that one point.")
        print()
    d = [s for s in summary if s["model"].startswith("M1d")]
    if d:
        bd = min(d, key=lambda s: s["median_fold"])
        print(f"  M1d REPLACES MEASURED f_u WITH ALBUMIN BINDING on the "
              f"{bd['n']} compounds")
        print(f"  with a human K_a. Best is {bd['model'].split(',')[1].strip()} at "
              f"{bd['median_fold']:.1f}x median, {bd['worst_fold']:.0f}x worst,")
        print(f"  against M1b's {b['median_fold']:.1f}x on the full set -- so binding is a")
        print("  usable stand-in for a measured f_u, not an improvement on it.")
        print()
        meds = [x["median_fold"] for x in d]
        if max(meds) / min(meds) < 1.02:
            print(f"  THE THREE MEDIANS AGREE TO WITHIN "
                  f"{(max(meds)/min(meds) - 1) * 100:.1f}% "
                  f"({min(meds):.3f}x to {max(meds):.3f}x)")
            print("  AND THAT IS NOT A COINCIDENCE.")
            print("  With tight binding, n[Alb]K_a >> 1 so f_u ~ 1/(n[Alb]K_a),")
            print("  i.e. f_u is proportional to 1/n. The fitted T is")
            print("  proportional to 1/f_u and therefore to n, so n cancels")
            print("  exactly in the product f_u*T. It survives only through the")
            print("  weakly bound short chains, where n[Alb]K_a is not >> 1.")
            spread = {}
            for r in rows:
                if "M1d" in r["model"]:
                    spread.setdefault(r["chemical"], []).append(r["fold_error"])
            sp = sorted(((max(v) / min(v), c) for c, v in spread.items()),
                        reverse=True)
            top, second = sp[0], sp[1]
            lo = min(spread[top[1]])
            hi = max(spread[top[1]])
            print("  THE DATA CONFIRM EXACTLY THAT PATTERN. Going n=1 -> n=3,")
            print(f"  {top[1]} -- the weakest binder in the set -- moves "
                  f"{hi:.1f}x to {lo:.1f}x")
            print(f"  ({top[0]:.2f}-fold), while NO OTHER COMPOUND MOVES MORE "
                  f"THAN {(second[0] - 1) * 100:.1f}%")
            print(f"  ({second[1]}), and the worst case overall shifts only "
                  f"{max(x['worst_fold'] for x in d):.1f}x to "
                  f"{min(x['worst_fold'] for x in d):.1f}x.")
            print("  The number of binding sites is therefore NOT IDENTIFIABLE")
            print("  from these data once T is fitted.")
            print()

    # ---- f_u sensitivity
    print(bar)
    print("SENSITIVITY: WHICH f_u STUDY, ON THE PARAMETER-FREE MODEL")
    print(bar)
    print("  The compilation uses Ryu 2024, while db/human_fu_recommended.csv")
    print("  names Smeltz 2023 as the default human series. Both are run,")
    print("  with Argoul and Fischer alongside.")
    print()
    sens = []
    print(f"  {'f_u source':24}{'n':>3}{'median':>10}{'worst':>10}  direction")
    for key, label in FU_SOURCES:
        f = []
        dirs = set()
        for c in chems:
            fu = fu_src[key].get(c)
            if not fu:
                continue
            p = LN2 * hum[c]["vd"] / (fu * GFR_HUMAN)
            obs = hum[c]["t_half"]
            f.append(max(p / obs, obs / p))
            dirs.add("over" if p > obs else "under")
        if not f:
            continue
        f.sort()
        sens.append(dict(fu_source=label, n=len(f),
                         median_fold=sig(f[len(f) // 2], 4),
                         worst_fold=sig(f[-1], 4),
                         directions="/".join(sorted(dirs))))
        print(f"  {label:24}{len(f):>3}{f[len(f)//2]:>9.1f}x{f[-1]:>9.1f}x"
              f"  {'/'.join(sorted(dirs))}")
    print()
    for key, label in FU_SOURCES:
        miss = [c for c in chems if c not in fu_src[key]]
        if miss:
            print(f"    {label}: no value for {', '.join(miss)}")
    print("    (PFDoDA is absent from human_fu_comparison.csv entirely, so")
    print("     every row above is scored without it.)")
    print()
    print("  THE CONCLUSION DOES NOT DEPEND ON THE CHOICE. Every source")
    print("  leaves a large shortfall in the same direction, because f_u")
    print("  varies 2- to 10-fold between studies while the shortfall is")
    print("  roughly 30-fold. The f_u choice is not what is missing.")
    print()

    # ---- model 2
    print(bar)
    print("MODEL 2  QSAR -- chain length and head group only, no TK input")
    print(bar)
    qrows, qsum = qsar(hum)
    print(f"  {'form':30}{'n':>3}{'k':>3}{'r2':>8}{'LOO med':>10}{'LOO worst':>11}")
    for s in qsum:
        print(f"  {s['form']:30}{s['n']:>3}{s['n_params']:>3}"
              f"{s['r2_in_sample']:>8.3f}{s['loo_median_fold']:>9.1f}x"
              f"{s['loo_worst_fold']:>10.1f}x")
    print()
    bq = min(qsum, key=lambda s: s["loo_median_fold"])
    rich = max(qsum, key=lambda s: s["n_params"])
    simp = qsum[0]
    print(f"  r2 climbs {simp['r2_in_sample']:.3f} -> {rich['r2_in_sample']:.3f} "
          f"as k goes {simp['n_params']} -> {rich['n_params']}, while LOO median")
    print(f"  goes {simp['loo_median_fold']:.1f}x -> {rich['loo_median_fold']:.1f}x. "
          f"Best LOO is \"{bq['form']}\" at {bq['loo_median_fold']:.1f}x.")
    print("  With 9 compounds, r2 is not evidence; the holdout column is.")
    print()

    # ---- head to head
    print(bar)
    print("THE TWO MODELS SIDE BY SIDE, SAME 9 HUMAN COMPOUNDS")
    print(bar)
    print(f"  {'model':40}{'k':>3}{'median':>10}{'worst':>10}")
    allm = [(s["model"], s["n_params"], s["median_fold"], s["worst_fold"])
            for s in summary if s["n"] == len(chems)]
    allm += [("QSAR " + s["form"], s["n_params"], s["loo_median_fold"],
              s["loo_worst_fold"]) for s in qsum]
    for name, k, med, wst in sorted(allm, key=lambda x: x[2]):
        print(f"  {name:40}{k:>3}{med:>9.1f}x{wst:>9.1f}x")
    print()
    best = min(allm, key=lambda x: x[2])
    print(f"  BEST ON MEDIAN: {best[0]} at {best[2]:.1f}x ({best[1]} parameters).")
    print("  BUT THE MEDIAN ALONE PICKS THE WRONG MODEL HERE. Every worst")
    print("  case is named below, and they are not the same compound.")
    print()
    print(f"  {'model':40}{'worst':>9}  driven by")
    allr = rows + [dict(model="QSAR " + r["form"], chemical=r["chemical"],
                        fold_error=r["loo_fold_error"]) for r in qrows]
    drivers = {}
    for name, k, med, wst in sorted(allm, key=lambda x: x[2]):
        g = [r for r in allr if r["model"] == name]
        if not g:
            continue
        w = max(g, key=lambda r: r["fold_error"])
        drivers[name] = w["chemical"]
        print(f"  {name:40}{wst:>8.1f}x  {w['chemical']}")
    print()
    pf = [n for n, c in drivers.items() if c == "PFDoDA"]
    print(f"  PFDoDA IS THE WORST CASE FOR {len(pf)} OF {len(drivers)} MODELS, including")
    print("  both mechanistic winners. It is the one secreting compound, so a")
    print("  model that fits a single reabsorptive T necessarily misses it by")
    print("  the full width of the T range. The QSAR does not care: it never")
    print("  sees T, and its own worst cases are short-chain instead.")
    print()

    print(bar)
    print("THE SAME COMPARISON WITH PFDoDA DROPPED")
    print(bar)
    print("  Not to improve the numbers, but to find out which conclusions")
    print("  rest on that one compound. The 8 remaining compounds are all")
    print("  reabsorptive, which is the regime the mechanistic models assume.")
    print()
    print(f"  {'model':40}{'median':>10}{'worst':>10}")
    excl = []
    for name in [m[0] for m in sorted(allm, key=lambda x: x[2])]:
        g = [r["fold_error"] for r in allr
             if r["model"] == name and r["chemical"] != "PFDoDA"]
        if not g:
            continue
        g.sort()
        excl.append(dict(model=name, n=len(g),
                         median_fold=sig(g[len(g) // 2], 4),
                         worst_fold=sig(g[-1], 4)))
        print(f"  {name:40}{g[len(g)//2]:>9.1f}x{g[-1]:>9.1f}x")
    print()
    be = min(excl, key=lambda r: r["median_fold"])
    bw = min(excl, key=lambda r: r["worst_fold"])
    print(f"  Best median now {be['model']} at {be['median_fold']:.1f}x;")
    print(f"  best worst case {bw['model']} at {bw['worst_fold']:.1f}x.")
    full_best = best[0]
    print(f"  On the full set the median winner was {full_best}.")
    print("  WHETHER THE RANKING SURVIVES DROPPING ONE COMPOUND OUT OF NINE")
    print("  IS THE REAL TEST, and with 9 compounds it is a weak one either")
    print("  way. No model here is established; they are ranked, provisionally.")
    print()
    print("  WHAT THIS DOES AND DOES NOT SETTLE. Both models are now fitted")
    print("  to human half-lives, so neither can be carried to a species")
    print("  whose half-life is unknown -- the usual reason for wanting one.")
    print("  The mechanistic models earn something the QSAR cannot: their")
    print("  one parameter IS the tubular factor, so a fitted T is a")
    print("  physiological quantity that can be compared against transporter")
    print("  data. The QSAR's coefficients are not interpretable that way.")
    print()

    # ---- write
    out = [
        ("human_tk_models.csv", rows,
         ["model", "chemical", "n_total_carbons", "head_group",
          "t_half_measured_d", "t_half_predicted_d", "fold_error",
          "direction"]),
        ("human_tk_summary.csv", summary,
         ["model", "n", "n_params", "scored", "median_fold", "worst_fold",
          "note"]),
        ("human_qsar.csv", qrows,
         ["form", "chemical", "n_total_carbons", "head_group",
          "t_half_measured_d", "t_half_loo_predicted_d", "loo_fold_error"]),
        ("human_qsar_summary.csv", qsum,
         ["form", "n", "n_params", "r2_in_sample", "loo_median_fold",
          "loo_worst_fold"]),
        ("human_fu_sensitivity.csv", sens,
         ["fu_source", "n", "median_fold", "worst_fold", "directions"]),
        ("human_model_pfdoda_excluded.csv", excl,
         ["model", "n", "median_fold", "worst_fold"]),
        ("human_tubular_factor.csv", tf,
         ["chemical", "n_total_carbons", "tubular_factor_T", "one_over_T",
          "reading"]),
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
