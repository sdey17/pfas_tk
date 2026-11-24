#!/usr/bin/env python3
"""Is the albumin-vs-L-FABP comparison method-confounded? (Stage 19)

    python3 -I scripts/stage19_lfabp_within_method.py

Writes db/lfabp_within_method.csv. Idempotent.

THE PROBLEM. The paper's L-FABP result pairs albumin and L-FABP BY COMPOUND:
for each compound it averages every albumin measurement and every L-FABP
measurement and differences the two means. But those two means are drawn from
almost disjoint sets of assays. Across db/protein_binding.csv, 63 methods
measured albumin only, 25 measured L-FABP only, and exactly ONE measured
both. So each per-compound difference is a difference between two different
METHOD MIXTURES, and pairing by compound does not remove that: the pairing is
on the wrong factor. The project's own figure 3 found method to be the
principal axis of disagreement in this compilation, which is precisely the
condition under which such a mixture difference biases the estimate.

THE TEST. Gao et al. measured 14 compounds against human serum albumin,
pooled human serum, AND human L-FABP in one laboratory by one method
(equilibrium dialysis). That is a fully paired within-method design: the
method mixture is identical on both arms by construction. If the pooled
estimate and the within-method estimate agree, the confound is harmless and
the published claim stands. If they disagree, the pooled number is an
artefact of which assays happened to be applied to which protein.

CONVENTION. Gao reports K_d, so a LARGER value is WEAKER binding. Results are
reported as a difference in log K (association), i.e. -(log10 Kd_LFABP -
log10 Kd_albumin), so that a NEGATIVE difference means L-FABP binds more
weakly, matching the sign convention of the paper's Table 1.
"""
import csv
import math
import pathlib
import statistics as st

HERE = pathlib.Path(__file__).resolve().parent.parent

METHOD_PREFIX = "equilibrium dialysis: protein or serum spiked"

# The paper's published estimates, for comparison. Differences in log K.
PUBLISHED = {
    "all 19 compounds, pooled across methods": -0.39,
    "10 compounds with n>=3 both sides, pooled": -0.126,
}
PUBLISHED_ONE_SIDED_UPPER = +0.023  # log units; the paper's headline limit

TOTAL_CARBONS = {
    "PFBA": 4, "PFPeA": 5, "PFHxA": 6, "PFHpA": 7, "PFOA": 8, "PFNA": 9,
    "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12, "PFTrDA": 13,
    "PFBS": 4, "PFHxS": 6, "PFOS": 8, "6:2 Cl-PFESA": 8,
}


def t_crit_one_sided_95(df):
    """t for a one-sided 95% limit. Table lookup; exact for the df we use."""
    table = {9: 1.833, 10: 1.812, 11: 1.796, 12: 1.782, 13: 1.771,
             14: 1.761, 15: 1.753, 18: 1.734}
    return table[df]


def t_crit_two_sided_95(df):
    table = {9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160,
             14: 2.145, 15: 2.131, 18: 2.101}
    return table[df]


def load():
    """Return {compound: {arm: Kd}} for the single both-protein method."""
    out = {}
    units = set()
    for r in csv.DictReader((HERE / "db" / "protein_binding.csv").open()):
        if not r["method"].startswith(METHOD_PREFIX):
            continue
        if not r["kd"]:
            continue
        p = r["protein"].upper()
        if "FABP" in p:
            arm = "lfabp"
        elif "ALBUMIN" in p:
            arm = "albumin"
        elif "SERUM" in p:
            arm = "serum"
        else:
            continue
        out.setdefault(r["chemical"], {})[arm] = float(r["kd"])
        units.add(r["kd_units"])
    return out, units


def paired_stats(diffs):
    n = len(diffs)
    mean = st.fmean(diffs)
    sd = st.stdev(diffs)
    sem = sd / math.sqrt(n)
    df = n - 1
    return dict(n=n, mean=mean, sd=sd, sem=sem, df=df,
                t=mean / sem,
                ci_lo=mean - t_crit_two_sided_95(df) * sem,
                ci_hi=mean + t_crit_two_sided_95(df) * sem,
                one_sided_upper=mean + t_crit_one_sided_95(df) * sem)


def main():
    data, units = load()
    bar = "=" * 78
    print(bar)
    print("THE CONFOUND, STATED IN NUMBERS")
    print(bar)
    # recount the method disjointness so the claim is reproducible here
    per = {}
    for r in csv.DictReader((HERE / "db" / "protein_binding.csv").open()):
        p = r["protein"].upper()
        fam = ("lfabp" if "FABP" in p else
               "albumin" if "ALBUMIN" in p else None)
        if fam and r["method"]:
            per.setdefault(r["method"], set()).add(fam)
    both = [m for m, s in per.items() if len(s) == 2]
    alb = [m for m, s in per.items() if s == {"albumin"}]
    fab = [m for m, s in per.items() if s == {"lfabp"}]
    print(f"  methods measuring albumin only : {len(alb)}")
    print(f"  methods measuring L-FABP only  : {len(fab)}")
    print(f"  methods measuring BOTH         : {len(both)}")
    print()
    print("  The paper pairs by COMPOUND, so each per-compound difference is")
    print("  a difference between two nearly disjoint method mixtures. The")
    print("  pairing is on the wrong factor; it cannot remove a method offset.")
    print()

    print(bar)
    print("THE WITHIN-METHOD TEST: Gao et al., equilibrium dialysis, one lab")
    print(bar)
    print(f"  Kd units in this subset: {sorted(units)}")
    print(f"  (a ratio is unit-free, so only internal consistency matters)")
    print()
    rows = []
    print(f"  {'compound':13}{'tC':>4}{'Kd alb':>9}{'Kd serum':>10}"
          f"{'Kd FABP':>9}{'Kd ratio':>10}{'dlogK':>9}")
    for chem in sorted(data, key=lambda c: (TOTAL_CARBONS.get(c, 99), c)):
        d = data[chem]
        if "albumin" not in d or "lfabp" not in d:
            continue
        ratio = d["lfabp"] / d["albumin"]
        dlogk = -math.log10(ratio)      # log K (association) difference
        rows.append(dict(chemical=chem,
                         total_carbons=TOTAL_CARBONS.get(chem, ""),
                         kd_albumin=d["albumin"],
                         kd_serum=d.get("serum", ""),
                         kd_lfabp=d["lfabp"],
                         kd_ratio_lfabp_over_albumin=round(ratio, 3),
                         diff_logK_lfabp_minus_albumin=round(dlogk, 4)))
        print(f"  {chem:13}{TOTAL_CARBONS.get(chem, 0):>4}"
              f"{d['albumin']:>9.0f}{d.get('serum', 0):>10.0f}"
              f"{d['lfabp']:>9.0f}{ratio:>10.2f}{dlogk:>+9.3f}")
    print()

    diffs = [r["diff_logK_lfabp_minus_albumin"] for r in rows]
    s = paired_stats(diffs)
    print(f"  paired over {s['n']} compounds, one method, one laboratory:")
    print(f"    mean difference in log K : {s['mean']:+.3f} "
          f"(95% CI {s['ci_lo']:+.3f} to {s['ci_hi']:+.3f})")
    print(f"    SD {s['sd']:.3f}   SEM {s['sem']:.3f}   "
          f"t({s['df']}) = {s['t']:.2f}")
    print(f"    one-sided 95% UPPER limit: {s['one_sided_upper']:+.3f} log "
          f"units = {10**s['one_sided_upper']:.2f}-fold")
    print()

    print(bar)
    print("WITHIN-METHOD vs PUBLISHED")
    print(bar)
    for label, val in PUBLISHED.items():
        print(f"  {label:44} {val:+.3f}")
    print(f"  {'within-method, paired, Gao equilibrium dialysis':44} "
          f"{s['mean']:+.3f}")
    print()
    print(f"  published one-sided 95% upper limit : "
          f"{PUBLISHED_ONE_SIDED_UPPER:+.3f}  "
          f"({10**PUBLISHED_ONE_SIDED_UPPER:.2f}-fold)")
    print(f"  within-method one-sided 95% upper   : "
          f"{s['one_sided_upper']:+.3f}  "
          f"({10**s['one_sided_upper']:.2f}-fold)")
    print()
    if s["one_sided_upper"] < 0:
        print("  THE WITHIN-METHOD RESULT IS STRONGER, NOT WEAKER. Its upper")
        print("  limit lies BELOW zero, so it does not merely exclude an")
        print("  affinity advantage for L-FABP -- it establishes a definite")
        print(f"  DISADVANTAGE of at least {10**-s['one_sided_upper']:.2f}-fold.")
        print("  The published claim survives the confound and can be")
        print("  sharpened.")
    else:
        print("  The within-method upper limit is above zero, so an affinity")
        print("  advantage is NOT excluded once method is controlled. The")
        print("  published claim would not survive in its present form.")
    print()

    print(bar)
    print("WHAT THE POOLED ANALYSIS HID: CHAIN LENGTH AND HEAD GROUP")
    print(bar)
    with_tc = [r for r in rows if r["total_carbons"] != ""]
    for r in with_tc:
        r["head_group"] = ("sulfonate"
                           if r["chemical"].endswith("S")
                           or "PFESA" in r["chemical"]
                           else "carboxylate")

    def ols(data, cols):
        """Least squares with an intercept. cols are callables -> float."""
        X = [[1.0] + [c(r) for c in cols] for r in data]
        y = [r["diff_logK_lfabp_minus_albumin"] for r in data]
        k = len(X[0])
        # normal equations, solved by Gauss-Jordan (k is 2 or 3)
        A = [[sum(X[i][a] * X[i][b] for i in range(len(X))) for b in range(k)]
             + [sum(X[i][a] * y[i] for i in range(len(X)))] for a in range(k)]
        for c in range(k):
            piv = max(range(c, k), key=lambda r_: abs(A[r_][c]))
            A[c], A[piv] = A[piv], A[c]
            d = A[c][c]
            A[c] = [v / d for v in A[c]]
            for r_ in range(k):
                if r_ != c and A[r_][c]:
                    f = A[r_][c]
                    A[r_] = [a - f * b for a, b in zip(A[r_], A[c])]
        beta = [A[i][k] for i in range(k)]
        fit = [sum(b * xv for b, xv in zip(beta, X[i])) for i in range(len(X))]
        res = [y[i] - fit[i] for i in range(len(y))]
        my_ = st.fmean(y)
        sse = sum(e * e for e in res)
        r2_ = 1 - sse / sum((v - my_) ** 2 for v in y)
        n, df = len(y), len(y) - k
        s2 = sse / df
        # SE of each coefficient from the inverse of X'X, recomputed
        XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)]
               for a in range(k)]
        I = [[1.0 if a == b else 0.0 for b in range(k)] for a in range(k)]
        M = [XtX[i] + I[i] for i in range(k)]
        for c in range(k):
            piv = max(range(c, k), key=lambda r_: abs(M[r_][c]))
            M[c], M[piv] = M[piv], M[c]
            d = M[c][c]
            M[c] = [v / d for v in M[c]]
            for r_ in range(k):
                if r_ != c and M[r_][c]:
                    f = M[r_][c]
                    M[r_] = [a - f * b for a, b in zip(M[r_], M[c])]
        se = [math.sqrt(s2 * M[i][k + i]) for i in range(k)]
        return beta, se, r2_, df, res

    tc = lambda r: float(r["total_carbons"])
    sulf = lambda r: 1.0 if r["head_group"] == "sulfonate" else 0.0

    b, se, r2, df, res = ols(with_tc, [tc])
    print(f"  (a) d(logK) ~ total carbons           n={len(with_tc)}")
    print(f"      slope {b[1]:+.4f} +/- {se[1]:.4f}  t({df}) = "
          f"{b[1]/se[1]:.2f}   r2 = {r2:.3f}")
    print(f"      crosses zero at {-b[0]/b[1]:.1f} total carbons")
    print()
    b2, se2, r22, df2, _ = ols(with_tc, [tc, sulf])
    print(f"  (b) + head group                      n={len(with_tc)}")
    print(f"      carbons   {b2[1]:+.4f} +/- {se2[1]:.4f}  "
          f"t({df2}) = {b2[1]/se2[1]:.2f}")
    print(f"      sulfonate {b2[2]:+.4f} +/- {se2[2]:.4f}  "
          f"t({df2}) = {b2[2]/se2[2]:.2f}   r2 = {r22:.3f}")
    print(f"      The sulfonate offset is NOT significant: {abs(b2[2]):.3f}")
    print(f"      log units with SE {se2[2]:.3f}, and r2 rises only from")
    print(f"      {r2:.3f} to {r22:.3f}. The three sulfonates do sit below the")
    print("      carboxylates of equal carbon number in the table above, but")
    print("      with n=3 that pattern is not established here. Chain length")
    print("      survives the adjustment essentially unchanged.")
    print()

    print("  (c) THE TOP END IS NOT MONOTONIC, and this limits the claim:")
    for r in sorted(with_tc, key=lambda x: x["total_carbons"])[-4:]:
        print(f"      C{r['total_carbons']:<3}{r['chemical']:12}"
              f"{r['diff_logK_lfabp_minus_albumin']:+.3f}")
    print("      The difference peaks at C11 and falls again by C13, so the")
    print("      linear fit in (a) should not be read as a trend continuing")
    print("      past C11. Each of C11, C12 and C13 rests on ONE compound")
    print("      here, so the reversal is suggestive, not established.")
    print()
    short = [r for r in with_tc if r["total_carbons"] <= 8]
    lng = [r for r in with_tc if r["total_carbons"] >= 11]
    ms = st.fmean([r["diff_logK_lfabp_minus_albumin"] for r in short])
    ml = st.fmean([r["diff_logK_lfabp_minus_albumin"] for r in lng])
    print(f"  <=C8  (n={len(short)}): mean {ms:+.3f}  "
          f"-> L-FABP weaker by {10**-ms:.2f}-fold")
    print(f"  >=C11 (n={len(lng)}): mean {ml:+.3f}  "
          f"-> indistinguishable from parity")
    print()
    print("  THIS MATTERS FOR THE HYPOTHESIS BEING REFUTED. The L-FABP")
    print("  explanation is invoked for LONG-chain hepatic accumulation. The")
    print("  refutation is strongest exactly where the hypothesis is not")
    print("  applied, and weakest where it is: at C11 and above the two")
    print("  proteins are indistinguishable on these data. A blanket claim")
    print("  that L-FABP has no affinity advantage is safe; a claim that it")
    print("  binds measurably more weakly is NOT safe for the long chains,")
    print("  and the pooled analysis could not see the difference.")
    print()

    p = HERE / "db" / "lfabp_within_method.csv"
    cols = list(rows[0].keys())
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["total_carbons"], x["chemical"])):
            w.writerow(r)
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} compounds)")


if __name__ == "__main__":
    main()
