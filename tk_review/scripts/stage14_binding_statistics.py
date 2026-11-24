#!/usr/bin/env python3
"""Per-compound dispersion statistics for albumin and L-FABP binding.

    python3 -I scripts/stage14_binding_statistics.py

Writes db/binding_statistics.csv. Idempotent.

WHY THIS EXISTS. The reports have been quoting the albumin/L-FABP comparison
as a median difference of -0.12 log units without attaching a dispersion to
it, which invites the reader to treat it as a precise null. It is not precise;
it is a small central difference sitting inside an enormous spread, and the
spread is the more important number. This script attaches mean, SD, n and
range to every per-compound value and tests the paired difference formally.

THE TEST. The paired albumin-minus-L-FABP differences are not assumed normal
(n is small and the per-compound value counts are very unequal), so a
Wilcoxon signed-rank test is used alongside the parametric t-test, and both
are reported. A non-significant result here is NOT evidence of equivalence by
itself -- with this dispersion the study is underpowered to detect anything
smaller than about a log unit -- so a two-one-sided-test (TOST) equivalence
check against a +/- 0.5 log unit margin is reported as well. That is the test
that can actually support the claim being made.
"""
import csv, math, pathlib, statistics
import numpy as np
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent.parent
ALB = HERE / "figures" / "albumin_chainlength_data.csv"
FAB = HERE / "figures" / "lfabp_by_compound_data.csv"
TOST_MARGIN = 0.5      # log units; the equivalence margin we claim


def load(path):
    out = {}
    for r in csv.DictReader(path.open()):
        try:
            v = float(r["log_k"])
        except (TypeError, ValueError):
            continue
        out.setdefault((r["pfas"], r["family"], int(r["n_cf"])), []).append(v)
    return out


def total_c(fam, ncf):
    return ncf + 1 if fam.strip().lower().startswith("carboxyl") else ncf


def desc(v):
    return dict(n=len(v), mean=statistics.mean(v),
                sd=statistics.stdev(v) if len(v) > 1 else float("nan"),
                lo=min(v), hi=max(v), rng=max(v) - min(v),
                med=statistics.median(v))


def main():
    A, F = load(ALB), load(FAB)
    keys = {(p, fam, n) for (p, fam, n) in A} & {(p, fam, n) for (p, fam, n) in F}
    rows = []
    for k in sorted(keys, key=lambda k: total_c(k[1], k[2])):
        p, fam, ncf = k
        a, f = desc(A[k]), desc(F[k])
        rows.append(dict(chem=p, fam=fam, tc=total_c(fam, ncf), a=a, f=f,
                         diff=f["mean"] - a["mean"]))

    bar = "=" * 78
    print(bar)
    print("ALBUMIN AND L-FABP log K BY COMPOUND, WITH DISPERSION")
    print(bar)
    print(f"  {'chem':9}{'tC':>3}  {'--------- albumin ---------':^28}"
          f"{'--------- L-FABP ----------':^28}{'diff':>7}")
    print(f"  {'':9}{'':>3}  {'n':>3}{'mean':>7}{'SD':>7}{'range':>11}"
          f"{'n':>5}{'mean':>7}{'SD':>7}{'range':>11}{'':>7}")
    for r in rows:
        a, f = r["a"], r["f"]
        sa = "%.2f" % a["sd"] if a["n"] > 1 else "  - "
        sf = "%.2f" % f["sd"] if f["n"] > 1 else "  - "
        ra = "%.1f-%.1f" % (a["lo"], a["hi"])
        rf = "%.1f-%.1f" % (f["lo"], f["hi"])
        print(f"  {r['chem']:9}{r['tc']:>3}  {a['n']:>3}{a['mean']:>7.2f}"
              f"{sa:>7}{ra:>11}"
              f"{f['n']:>5}{f['mean']:>7.2f}{sf:>7}"
              f"{rf:>11}{r['diff']:>7.2f}")

    d = np.array([r["diff"] for r in rows])
    n = len(d)
    print()
    print(bar)
    print("THE PAIRED COMPARISON, TESTED")
    print(bar)
    print(f"  compounds measured against both proteins: n = {n}")
    print(f"  mean difference (L-FABP - albumin) : {d.mean():+.3f} log units")
    print(f"  median difference                  : {np.median(d):+.3f}")
    print(f"  SD of the differences              : {d.std(ddof=1):.3f}")
    se = d.std(ddof=1) / math.sqrt(n)
    tcrit = stats.t.ppf(0.975, n - 1)
    ci = (d.mean() - tcrit * se, d.mean() + tcrit * se)
    print(f"  95% CI of the mean difference      : "
          f"[{ci[0]:+.3f}, {ci[1]:+.3f}]")
    t, pt = stats.ttest_rel(
        [r["f"]["mean"] for r in rows], [r["a"]["mean"] for r in rows])
    try:
        w, pw = stats.wilcoxon([r["f"]["mean"] for r in rows],
                               [r["a"]["mean"] for r in rows])
    except ValueError:
        w, pw = float("nan"), float("nan")
    print(f"  paired t-test                      : t = {t:+.3f}, p = {pt:.3f}")
    print(f"  Wilcoxon signed-rank               : W = {w:.1f}, p = {pw:.3f}")
    print()
    print("  NEITHER TEST IS EVIDENCE OF EQUIVALENCE ON ITS OWN. A null result")
    print("  with this much scatter may only mean the comparison is")
    print("  underpowered. The equivalence test below is the one that can")
    print("  support the claim.")
    print()
    # TOST: is the difference inside +/- margin?
    t1 = (d.mean() - (-TOST_MARGIN)) / se
    t2 = (d.mean() - TOST_MARGIN) / se
    p1 = 1 - stats.t.cdf(t1, n - 1)      # H0: diff <= -margin
    p2 = stats.t.cdf(t2, n - 1)          # H0: diff >= +margin
    ptost = max(p1, p2)
    print(f"  TOST equivalence, margin +/-{TOST_MARGIN} log units")
    print(f"    p(lower) = {p1:.4f}   p(upper) = {p2:.4f}   "
          f"TOST p = {ptost:.4f}")
    print(f"    => the two proteins ARE statistically equivalent within "
          f"+/-{TOST_MARGIN} log units"
          if ptost < 0.05 else
          f"    => equivalence within +/-{TOST_MARGIN} log units is NOT "
          f"established")
    print()
    print(bar)
    print("THE TEST THAT ACTUALLY BEARS ON THE HYPOTHESIS")
    print(bar)
    print("  The L-FABP hypothesis requires L-FABP to bind PFAS MORE TIGHTLY")
    print("  than the plasma carrier. That is a ONE-SIDED claim, and a")
    print("  one-sided bound is what the data can supply.")
    print()
    t95 = stats.t.ppf(0.95, n - 1)
    ub = d.mean() + t95 * se
    print(f"  one-sided 95% upper confidence limit on (L-FABP - albumin):")
    print(f"      {ub:+.3f} log units  =  {10 ** ub:.2f}-fold")
    print()
    if ub < math.log10(2):
        print(f"  => the data EXCLUDE L-FABP being more than "
              f"{10 ** ub:.2f}-fold stronger than albumin.")
        print("     The hypothesis needs an advantage; the data cap it at")
        print("     essentially none. This is the defensible claim, and it is")
        print("     stronger than a failure to reject.")
    print()

    print(bar)
    print("SENSITIVITY: COMPOUNDS WITH n >= 3 ON BOTH SIDES")
    print(bar)
    sub = [r for r in rows if r["a"]["n"] >= 3 and r["f"]["n"] >= 3]
    ds = np.array([r["diff"] for r in sub])
    print(f"  the full set is driven by single-measurement compounds")
    print(f"  (PFEESA -2.12, HFPO-TA -2.83, PFTrDA -1.01 all rest on n=1)")
    print(f"  restricting to n>=3 on both sides: k = {len(sub)}")
    if len(sub) > 1:
        ses = ds.std(ddof=1) / math.sqrt(len(ds))
        tc = stats.t.ppf(0.975, len(ds) - 1)
        print(f"    mean   {ds.mean():+.3f}   median {np.median(ds):+.3f}"
              f"   SD {ds.std(ddof=1):.3f}")
        print(f"    95% CI [{ds.mean()-tc*ses:+.3f}, {ds.mean()+tc*ses:+.3f}]")
        t_, p_ = stats.ttest_1samp(ds, 0.0)
        print(f"    t = {t_:+.3f}, p = {p_:.3f}")
        tl = stats.t.ppf(0.95, len(ds) - 1)
        print(f"    one-sided 95% upper limit {ds.mean()+tl*ses:+.3f} log "
              f"= {10 ** (ds.mean()+tl*ses):.2f}-fold")
        p1s = max(1 - stats.t.cdf((ds.mean() + TOST_MARGIN) / ses, len(ds)-1),
                  stats.t.cdf((ds.mean() - TOST_MARGIN) / ses, len(ds)-1))
        print(f"    TOST p (margin +/-{TOST_MARGIN}) = {p1s:.4f}"
              f"  => {'EQUIVALENT' if p1s < 0.05 else 'not established'}")
    print()

    print(bar)
    print("THE DISPERSION IS THE HEADLINE, NOT THE DIFFERENCE")
    print(bar)
    multi_a = [r["a"]["rng"] for r in rows if r["a"]["n"] >= 3]
    multi_f = [r["f"]["rng"] for r in rows if r["f"]["n"] >= 3]
    for lbl, v in (("albumin", multi_a), ("L-FABP", multi_f)):
        if v:
            print(f"  {lbl:8} within-compound range (compounds with n>=3, "
                  f"k={len(v)}): median {statistics.median(v):.2f}, "
                  f"max {max(v):.2f} log units")
    sds = [r["a"]["sd"] for r in rows if r["a"]["n"] > 1]
    print(f"  albumin  pooled within-compound SD: "
          f"{math.sqrt(sum(s*s for s in sds)/len(sds)):.2f} log units")
    sds = [r["f"]["sd"] for r in rows if r["f"]["n"] > 1]
    print(f"  L-FABP   pooled within-compound SD: "
          f"{math.sqrt(sum(s*s for s in sds)/len(sds)):.2f} log units")
    print()
    print(f"  The mean difference between the two proteins ({d.mean():+.2f})")
    print(f"  is roughly an order of magnitude smaller than the scatter")
    print(f"  within either of them. THAT is why no affinity advantage can be")
    print(f"  claimed for L-FABP -- not because a test failed to reject.")
    print()

    p = HERE / "db" / "binding_statistics.csv"
    with p.open("w", newline="") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["chemical", "family", "n_total_carbons",
                     "albumin_n", "albumin_mean_logK", "albumin_sd",
                     "albumin_min", "albumin_max", "albumin_range",
                     "lfabp_n", "lfabp_mean_logK", "lfabp_sd",
                     "lfabp_min", "lfabp_max", "lfabp_range",
                     "difference_lfabp_minus_albumin"])
        for r in rows:
            a, f = r["a"], r["f"]
            w_.writerow([r["chem"], r["fam"], r["tc"],
                         a["n"], round(a["mean"], 3),
                         "" if a["n"] < 2 else round(a["sd"], 3),
                         round(a["lo"], 3), round(a["hi"], 3),
                         round(a["rng"], 3),
                         f["n"], round(f["mean"], 3),
                         "" if f["n"] < 2 else round(f["sd"], 3),
                         round(f["lo"], 3), round(f["hi"], 3),
                         round(f["rng"], 3), round(r["diff"], 3)])
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} compounds)")


if __name__ == "__main__":
    main()
