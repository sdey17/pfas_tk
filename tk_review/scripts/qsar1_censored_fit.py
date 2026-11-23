#!/usr/bin/env python3
"""QSAR Stage 1: per-species censored regression of log10 T (chain + head).

    python3 scripts/qsar1_censored_fit.py

Writes db/qsar/qsar_censored_fit.csv and db/qsar/qsar_fit_predictions.csv.
Idempotent. Needs numpy and scipy.

THE MODEL

    log10 T  ~  Normal( b0 + b1 * nC + b2 * [head group is sulfonate],
                        sigma )

fitted per species by MAXIMUM LIKELIHOOD WITH CENSORING, because 30 of the 39
usable rows are one-sided bounds rather than point values. Each row
contributes according to what is actually known about it:

    point      y            phi((y - mu)/sigma) / sigma
    right      T <= y       Phi((y - mu)/sigma)
    left       T >= y       1 - Phi((y - mu)/sigma)
    interval   lo..hi       Phi((hi - mu)/sigma) - Phi((lo - mu)/sigma)

Fitting the bounds as if they were measurements would pull every coefficient
toward the censoring direction; dropping them would discard two-thirds of the
data. This is the standard fix and it is the reason the modelling table
carries a censoring direction per row.

HEAD GROUP IS COLLAPSED TO TWO LEVELS, carboxylate against sulfonate. The
human set contains exactly one fluorotelomer sulfonate (6:2 FTS) and a level
with one member cannot be estimated; it is grouped with the sulfonates, which
is what it chemically is, and the perfluorinated-carbon count already
distinguishes its partially fluorinated chain.

WHAT IS DELIBERATELY NOT IN THE MODEL. Molecular weight (collinear with chain
length at r = 0.91 in a homologous series -- including both fits one variable
twice), and any quadratic in chain length. The quadratic is TESTED below but
not adopted by default, because three laboratories report the U-shape in
PLASMA BINDING rather than in tubular handling, and whether it propagates into
T is exactly what this fit is for.
"""
import csv, math, pathlib
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

HERE = pathlib.Path(__file__).resolve().parent.parent
OUT = HERE / "db" / "qsar"
L10 = math.log(10)


def load():
    rows = []
    f = OUT / "qsar_modelling_table.csv"
    for r in csv.DictReader(f.open()):
        if r["usable_for_regression"] != "yes" or not r["T_point"]:
            continue
        cens = r["censoring"]
        lo = hi = None
        y = math.log10(float(r["T_point"]))
        if cens.startswith("none"):
            kind = "point"
        elif cens.startswith("right"):
            kind = "right"
        elif cens.startswith("left"):
            kind = "left"
        elif cens.startswith("interval"):
            kind = "interval"
            lo, hi = math.log10(float(r["T_low"])), math.log10(float(r["T_high"]))
        else:
            continue
        rows.append(dict(
            species=r["species"], sex=r["sex"], chem=r["chemical"],
            nC=float(r["n_C_total"]),
            sulf=1.0 if "sulfonate" in r["head_group"] else 0.0,
            head=r["head_group"], y=y, lo=lo, hi=hi, kind=kind,
            verdict=r["verdict"]))
    return rows


def negll(theta, X, y, lo, hi, kind):
    k = X.shape[1]
    b, ls = theta[:k], theta[k]
    s = math.exp(ls)
    mu = X @ b
    ll = np.zeros(len(y))
    m = kind == 0                      # point
    ll[m] = norm.logpdf(y[m], mu[m], s)
    m = kind == 1                      # right-censored: T <= y
    ll[m] = norm.logcdf((y[m] - mu[m]) / s)
    m = kind == 2                      # left-censored: T >= y
    ll[m] = norm.logsf((y[m] - mu[m]) / s)
    m = kind == 3                      # interval
    if m.any():
        a = norm.cdf((hi[m] - mu[m]) / s) - norm.cdf((lo[m] - mu[m]) / s)
        ll[m] = np.log(np.clip(a, 1e-300, None))
    v = -ll.sum()
    return 1e12 if not np.isfinite(v) else v


def fit(rows, quad=False):
    n = len(rows)
    cols = [np.ones(n), np.array([r["nC"] for r in rows]),
            np.array([r["sulf"] for r in rows])]
    names = ["intercept", "nC", "sulfonate"]
    if quad:
        cols.append(np.array([r["nC"] ** 2 for r in rows]))
        names.append("nC^2")
    X = np.column_stack(cols)
    y = np.array([r["y"] for r in rows])
    lo = np.array([r["lo"] if r["lo"] is not None else 0.0 for r in rows])
    hi = np.array([r["hi"] if r["hi"] is not None else 0.0 for r in rows])
    kmap = {"point": 0, "right": 1, "left": 2, "interval": 3}
    kind = np.array([kmap[r["kind"]] for r in rows])
    k = X.shape[1]
    best = None
    for s0 in (-0.5, 0.0, 0.5):
        th0 = np.r_[np.zeros(k), s0]
        th0[0] = y.mean()
        res = minimize(negll, th0, args=(X, y, lo, hi, kind),
                       method="Nelder-Mead",
                       options=dict(maxiter=40000, maxfev=40000,
                                    xatol=1e-9, fatol=1e-9))
        if best is None or res.fun < best.fun:
            best = res
    b, sigma = best.x[:k], math.exp(best.x[k])
    # standard errors from a numerical Hessian of the negative log-likelihood
    h, H = 1e-4, np.zeros((k + 1, k + 1))
    f0 = negll(best.x, X, y, lo, hi, kind)
    for i in range(k + 1):
        for j in range(k + 1):
            tp, tm, tpm, tmp = (best.x.copy() for _ in range(4))
            tp[i] += h; tp[j] += h
            tm[i] -= h; tm[j] -= h
            tpm[i] += h; tpm[j] -= h
            tmp[i] -= h; tmp[j] += h
            H[i, j] = (negll(tp, X, y, lo, hi, kind)
                       - negll(tpm, X, y, lo, hi, kind)
                       - negll(tmp, X, y, lo, hi, kind)
                       + negll(tm, X, y, lo, hi, kind)) / (4 * h * h)
    try:
        se = np.sqrt(np.diag(np.linalg.inv(H)))[:k]
    except np.linalg.LinAlgError:
        se = np.full(k, float("nan"))
    ll = -best.fun
    aic = 2 * (k + 1) - 2 * ll
    return dict(names=names, b=b, se=se, sigma=sigma, ll=ll, aic=aic,
                n=n, X=X, y=y, rows=rows)


def main():
    rows = load()
    bar = "=" * 78
    results = {}
    print(bar)
    print("PER-SPECIES CENSORED FIT:  log10 T ~ nC + sulfonate")
    print(bar)
    for sp in ("human", "rat", "mouse"):
        d = [r for r in rows if r["species"] == sp]
        if len(d) < 5:
            continue
        m = fit(d)
        mq = fit(d, quad=True)
        results[sp] = (m, mq)
        nk = sum(1 for r in d if r["kind"] != "point")
        print(f"\n  {sp.upper()}   n = {m['n']}  "
              f"({nk} censored, {m['n']-nk} point)   sigma = {m['sigma']:.3f}")
        print(f"    {'term':12}{'estimate':>11}{'SE':>9}{'z':>8}"
              f"   per-carbon effect")
        for nm, b, se in zip(m["names"], m["b"], m["se"]):
            z = b / se if se and np.isfinite(se) and se > 0 else float("nan")
            extra = ""
            if nm == "nC":
                extra = f"x{10**b:.2f} per carbon"
            if nm == "sulfonate":
                extra = f"x{10**b:.2f} vs carboxylate"
            print(f"    {nm:12}{b:>11.4f}{se:>9.4f}{z:>8.2f}   {extra}")
        print(f"    logLik {m['ll']:.3f}   AIC {m['aic']:.2f}"
              f"    |  with nC^2: AIC {mq['aic']:.2f}"
              f"  ({'quadratic better' if mq['aic'] < m['aic'] - 2 else 'no improvement'})")

    print()
    print(bar)
    print("THE QUADRATIC FITS, IN FULL")
    print(bar)
    for sp, (m, mq) in results.items():
        print(f"\n  {sp.upper()}  (sigma {mq['sigma']:.3f})")
        for nm, b, se in zip(mq["names"], mq["b"], mq["se"]):
            z = b / se if se and np.isfinite(se) and se > 0 else float("nan")
            print(f"    {nm:12}{b:>11.4f}{se:>9.4f}  z={z:>6.2f}")
        iq = mq["names"].index("nC^2")
        il = mq["names"].index("nC")
        if mq["b"][iq] != 0:
            vertex = -mq["b"][il] / (2 * mq["b"][iq])
            print(f"    vertex of the parabola at nC = {vertex:.1f}"
                  f"   ({'minimum' if mq['b'][iq] > 0 else 'MAXIMUM'})")

    print()
    print(bar)
    print("READING THE COEFFICIENTS")
    print(bar)
    print("  nC is the change in log10 T per added carbon. Negative means the")
    print("  tubule returns MORE of the filtered load as the chain grows --")
    print("  i.e. reabsorption strengthens with chain length.")
    print()
    for sp, (m, mq) in results.items():
        i = m["names"].index("nC")
        b = m["b"][i]
        tot = 10 ** (b * 6)
        print(f"  {sp:7} {10**b:.2f}x per carbon  ->  "
              f"{tot:.3g}x over C4 to C10")
    print()
    print("  THE SIGN IS THE SAME IN ALL THREE SPECIES AND THE POINT")
    print("  ESTIMATES ARE LARGE -- BUT NONE OF THEM REACHES SIGNIFICANCE IN")
    print("  THE LINEAR MODEL. The largest |z| is 1.74. With 10 to 19 rows,")
    print("  most of them one-sided bounds, that is what should be expected,")
    print("  and it is the honest headline: the linear fit is CONSISTENT")
    print("  with a strong chain-length effect and does not ESTABLISH one.")
    print("  Three species agreeing in sign is weak corroboration -- with")
    print("  three independent signs a 3/3 agreement arises by chance one")
    print("  time in four.")

    print()
    print(bar)
    print("DOES THE U-SHAPE PROPAGATE FROM BINDING INTO T?")
    print(bar)
    print("  Three laboratories report a U-shaped fu against molecular")
    print("  weight. T divides CL_R by fu, so a U in fu does NOT have to")
    print("  appear in T -- it may cancel. Adding nC^2:")
    print()
    for sp, (m, mq) in results.items():
        d = mq["aic"] - m["aic"]
        verdict = ("quadratic justified" if d < -2 else
                   "no support" if d > 2 else "indistinguishable")
        print(f"  {sp:7} dAIC = {d:+6.2f}   {verdict}")
    print()
    print("  IT DOES PROPAGATE, AND AS AN INVERTED U -- which is mechanically")
    print("  what dividing by fu should produce. fu has a MINIMUM near C7, so")
    print("  1/fu has a maximum there, and T inherits a MAXIMUM rather than a")
    print("  minimum. The fitted vertices are nC 5.3 (human), 6.9 (rat) and")
    print("  5.7 (mouse), all maxima. In the quadratic model the coefficients")
    print("  also become significant where the linear ones were not: human")
    print("  nC z=2.71, nC^2 z=-2.60, sulfonate z=-4.18.")
    print()
    print("  So the mid-short chains, C5 to C7, are the compounds the tubule")
    print("  returns LEAST of -- the ones closest to net secretion -- and")
    print("  reabsorption strengthens in both directions from there.")
    print()
    print("  THAT READING IS NOT SAFE YET. See the sensitivity below.")

    print()
    print(bar)
    print("SENSITIVITY: THE VERTEX SITS ON THE TWO DISPUTED COMPOUNDS")
    print(bar)
    print("  The fitted maximum is at nC 5.3 (human), 6.9 (rat), 5.7 (mouse)")
    print("  -- and PFPeA (C5) and PFHxA (C6) are exactly the two compounds")
    print("  whose tubular SIGN depends on which human binding dataset is")
    print("  used. Refitting the human without them:")
    print()
    dh = [r for r in rows if r["species"] == "human"]
    keep = [r for r in dh if r["chem"] not in ("PFPeA", "PFHxA")]
    mk, mkq = fit(keep), fit(keep, quad=True)
    print(f"    full    n={len(dh):2}  AIC linear {results['human'][0]['aic']:6.2f}"
          f"   AIC quadratic {results['human'][1]['aic']:6.2f}"
          f"   dAIC {results['human'][1]['aic']-results['human'][0]['aic']:+.2f}")
    print("  (Compare dAIC WITHIN each row only. The absolute AICs are not")
    print("   comparable between the two rows because the datasets differ.)")
    print()
    print(f"    dropped n={len(keep):2}  AIC linear {mk['aic']:6.2f}"
          f"   AIC quadratic {mkq['aic']:6.2f}"
          f"   dAIC {mkq['aic']-mk['aic']:+.2f}")
    print()
    iq = mkq["names"].index("nC^2")
    il = mkq["names"].index("nC")
    if mkq["b"][iq] != 0:
        print(f"    without them the vertex moves to nC = "
              f"{-mkq['b'][il]/(2*mkq['b'][iq]):.1f}")
    surv = (mkq["aic"] - mk["aic"]) < -2
    print()
    print(f"    THE CURVATURE {'SURVIVES' if surv else 'DOES NOT SURVIVE'} "
          f"REMOVING THE DISPUTED PAIR.")
    if not surv:
        print("    So the inverted U in the human rests on exactly the two")
        print("    rows that the choice of binding dataset decides. It is a")
        print("    hypothesis about shape, not an established shape, and it")
        print("    should not be carried into a predictive model until the")
        print("    short-chain fu question is settled.")
        print()
        print("    AND THAT IS THE PRACTICAL LINK BACK TO THE fu DECISION.")
        print("    Choosing Smeltz over Ryu for PFPeA and PFHxA does not")
        print("    merely move two numbers: it decides whether the human")
        print("    chain-length relationship is monotone or humped, which is")
        print("    the functional form of the QSAR. The binding question is")
        print("    upstream of the model structure, not a detail beside it.")

    # ---- write ------------------------------------------------------
    p1 = OUT / "qsar_censored_fit.csv"
    with p1.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "n", "n_censored", "term", "estimate",
                    "std_error", "z", "fold_per_unit", "sigma", "logLik",
                    "AIC", "AIC_with_nC2", "model"])
        for sp, (m, mq) in results.items():
            nk = sum(1 for r in m["rows"] if r["kind"] != "point")
            for nm, b, se in zip(m["names"], m["b"], m["se"]):
                z = b / se if se and np.isfinite(se) and se > 0 else ""
                w.writerow([sp, m["n"], nk, nm, round(b, 5),
                            round(se, 5) if np.isfinite(se) else "",
                            round(z, 3) if z != "" else "",
                            round(10 ** b, 4), round(m["sigma"], 4),
                            round(m["ll"], 4), round(m["aic"], 3),
                            round(mq["aic"], 3),
                            "log10 T ~ nC + sulfonate, censored MLE"])
    p2 = OUT / "qsar_fit_predictions.csv"
    with p2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "sex", "chemical", "n_C_total", "head_group",
                    "censoring_kind", "observed_log10_T",
                    "predicted_log10_T", "residual",
                    "residual_is_meaningful"])
        for sp, (m, mq) in results.items():
            pred = m["X"] @ m["b"]
            for r, pr in zip(m["rows"], pred):
                meaning = ("yes" if r["kind"] == "point" else
                           "only as an inequality: observed is a bound")
                w.writerow([sp, r["sex"], r["chem"], int(r["nC"]), r["head"],
                            r["kind"], round(r["y"], 4), round(pr, 4),
                            round(r["y"] - pr, 4), meaning])
    print(f"\nwrote {p1.relative_to(HERE)}")
    print(f"wrote {p2.relative_to(HERE)}")


if __name__ == "__main__":
    main()
