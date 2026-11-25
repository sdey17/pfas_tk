#!/usr/bin/env python3
"""Half-life from albumin binding and nothing else (Stage 26).

    python3 -I scripts/stage26_albumin_only_model.py

Writes db/tk_albumin_only.csv. Idempotent.

THE CONSTRAINT, AND WHAT IT FORCES. The only input is an albumin binding
constant. No volume of distribution, no GFR, no unbound fraction, no
tubular factor. That rules out the mechanistic route of Stage 25, because a
binding constant has units of M^-1 and a half-life has units of time: with
nothing else to supply a volume and a clearance scale, the relationship can
only be CALIBRATED, not derived.

So this model is a one-variable regression:

    log10 t½ = a + b * log10 K_a

fitted within a species, on a single binding source, and scored by LEAVE-ONE-
OUT so that each compound is predicted from a fit that excluded it. That is
the fairest score available with nine points and two parameters.

WHAT THIS BUYS AND WHAT IT COSTS. It buys a model that cannot be wrong in the
way Stage 25's was: the calibration absorbs the tubular factor and the
V_d/GFR scale, so the systematic 128- to 1180-fold shortfall disappears by
construction. It costs the mechanism. The fitted constants a and b have no
physiological reading, they differ between species, and because they are
estimated FROM half-lives the model cannot be carried to a species whose
half-lives are unknown -- which is the main thing an interspecies model is
wanted for.

COVERAGE. Human albumin supports three sources covering all nine compounds.
Rat albumin supports one, Fan et al., which covers all eight rat compounds by
TWO methods that disagree by up to 63-fold within that single paper (PFBS
2.78 uM by ultrafiltration against 176.30 uM by calorimetry), so both are
run. MOUSE ALBUMIN DOES NOT EXIST IN THIS COMPILATION, so the mouse arm
cannot be built at all under this constraint.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

# (species, label, study, method substring or None)
ARMS = [
    ("human", "Li (fluorescence)", "Li et al.", None),
    ("human", "Chen (dialysis)", "Chen et al.", None),
    ("human", "Gao (dialysis)", "Gao et al.", None),
    ("rat", "Fan (ultrafiltration)", "Fan et al.", "ultrafiltration"),
    ("rat", "Fan (calorimetry)", "Fan et al.", "isothermal"),
]


def sig(x, n=6):
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


def ka_from_row(r):
    """K_a in M^-1, or None. Units differ by study; convert explicitly."""
    u_ka, u_kd = (r["ka_units"] or "").lower(), (r["kd_units"] or "").lower()
    if r["ka"]:
        try:
            v = float(r["ka"])
        except ValueError:
            return None
        if u_ka.startswith("log"):
            return 10 ** v
        if u_ka.startswith("m^-1"):
            return v
        return None
    if r["kd"]:
        try:
            v = float(r["kd"])
        except ValueError:
            return None
        return {"um": 1e-6, "mm": 1e-3, "m": 1.0}.get(u_kd) and \
            1.0 / (v * {"um": 1e-6, "mm": 1e-3, "m": 1.0}[u_kd])
    return None


def half_lives():
    """{species: {chemical: t_half_days}}"""
    out = {"human": {}, "rat": {}, "mouse": {}}
    for r in csv.DictReader(
            (HERE / "db" / "species_decomposition.csv").open()):
        out["human"][r["chemical"]] = float(r["t_half_human_d"])
        sp = r["comparison"].split("/")[1]
        out[sp][r["chemical"]] = float(r["t_half_animal_d"])
    return out


def binding(species, study, method):
    out = {}
    for r in csv.DictReader((HERE / "db" / "protein_binding.csv").open()):
        if ("ALBUMIN" not in r["protein"].upper()
                or r["protein_species"] != species or r["study"] != study):
            continue
        if method and method not in r["method"]:
            continue
        ka = ka_from_row(r)
        if ka:
            out.setdefault(r["chemical"], []).append(ka)
    return {c: 10 ** (sum(math.log10(v) for v in vs) / len(vs))
            for c, vs in out.items()}


def fit(xs, ys):
    """Least squares for y = a + b x. Returns (a, b) or None."""
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def main():
    hl = half_lives()
    bar = "=" * 78
    print(bar)
    print("THE MODEL: ALBUMIN BINDING AND NOTHING ELSE")
    print(bar)
    print("    log10 t½ = a + b * log10 K_a")
    print()
    print("  No V_d, no GFR, no f_u, no T. A binding constant is M^-1 and a")
    print("  half-life is time, so with nothing else to set the scale the")
    print("  relationship has to be fitted rather than derived. a and b are")
    print("  calibration constants with no physiological reading.")
    print()
    print(f"  MOUSE: not possible. The compilation holds no mouse albumin")
    print(f"  binding at all, so this constraint excludes the mouse entirely.")
    print()

    rows, summary = [], []
    for sp, label, study, method in ARMS:
        b = binding(sp, study, method)
        pairs = sorted((c, k) for c, k in b.items() if c in hl[sp])
        if len(pairs) < 4:
            print(f"  {label}: only {len(pairs)} compounds, skipped")
            continue
        chems = [c for c, _ in pairs]
        xs = [math.log10(k) for _, k in pairs]
        ys = [math.log10(hl[sp][c]) for c in chems]
        ab = fit(xs, ys)
        pred_in = [ab[0] + ab[1] * x for x in xs]
        my = sum(ys) / len(ys)
        r2 = 1 - (sum((ys[i] - pred_in[i]) ** 2 for i in range(len(ys)))
                  / sum((v - my) ** 2 for v in ys))
        folds = []
        for i in range(len(chems)):
            xi = [xs[j] for j in range(len(xs)) if j != i]
            yi = [ys[j] for j in range(len(ys)) if j != i]
            f2 = fit(xi, yi)
            if not f2:
                continue
            pr = f2[0] + f2[1] * xs[i]
            fold = 10 ** abs(pr - ys[i])
            folds.append(fold)
            rows.append(dict(
                species=sp, source=label, chemical=chems[i],
                ka_M_inv=sig(10 ** xs[i]),
                t_half_measured_d=hl[sp][chems[i]],
                t_half_loo_predicted_d=sig(10 ** pr),
                loo_fold_error=sig(fold, 4)))
        folds.sort()
        med = folds[len(folds) // 2]
        summary.append(dict(species=sp, source=label, n=len(chems),
                            slope=sig(ab[1], 4), intercept=sig(ab[0], 4),
                            r2=sig(r2, 3), loo_median=sig(med, 4),
                            loo_worst=sig(max(folds), 4)))
        print(f"  {sp.upper():6} {label:24} n={len(chems)}  "
              f"slope {ab[1]:+.2f}  r² {r2:.3f}   "
              f"LOO median {med:.1f}x   worst {max(folds):.1f}x")
    print()

    print(bar)
    print("READING THE SLOPES")
    print(bar)
    print("  A model in which tighter binding means a longer half-life needs")
    print("  a POSITIVE slope. The sign is therefore the first thing to look")
    print("  at, before any measure of fit.")
    print()
    print("  MECHANISM FIXES THE SLOPE AT +1, WHICH MAKES IT A DIAGNOSTIC.")
    print("  If V_d, GFR and T do not vary systematically with K_a, then")
    print("      t\u00bd \u221d 1/f_u,   f_u \u2248 1/(n[Alb]K_a) for tight binding")
    print("  so t\u00bd \u221d K_a and the slope of log t\u00bd on log K_a must be +1.")
    print("  A slope far from 1 means either that those three terms DO track")
    print("  binding, or that the binding numbers are not measuring what the")
    print("  model assumes.")
    print()
    pos = [s for s in summary if s["slope"] > 0]
    for s in summary:
        sl = s["slope"]
        if sl <= 0:
            d = "WRONG SIGN"
        elif 0.5 <= sl <= 2.0:
            d = "within 2x of the mechanistic +1"
        else:
            d = f"{max(sl, 1/sl):.1f}x from the mechanistic +1"
        print(f"    {s['species']:6} {s['source']:24} {sl:+7.2f}  {d}")
    print()
    print(f"  {len(pos)} of {len(summary)} arms give the expected sign, and the")
    print(f"  slopes span {min(s['slope'] for s in summary):+.2f} to "
          f"{max(s['slope'] for s in summary):+.2f} -- a range that includes a")
    print("  sign reversal BETWEEN TWO METHODS IN THE SAME PAPER (Fan, rat).")
    print("  The fitted relationship is therefore set more by which binding")
    print("  assay supplied K_a than by the binding itself.")
    print()

    print(bar)
    print("AGAINST THE ALTERNATIVES ALREADY ON THE TABLE")
    print(bar)
    print(f"  {'model':44}{'LOO median':>13}{'LOO worst':>12}")
    best = min(summary, key=lambda s: s["loo_median"])
    worst = max(summary, key=lambda s: s["loo_median"])
    print(f"  {'albumin binding only, best arm':44}"
          f"{best['loo_median']:>12.1f}x{best['loo_worst']:>11.1f}x")
    print(f"  {'albumin binding only, worst arm':44}"
          f"{worst['loo_median']:>12.1f}x{worst['loo_worst']:>11.1f}x")
    print(f"  {'chain length + head group, quadratic (Stage 25)':44}"
          f"{3.4:>12.1f}x{9.0:>11.1f}x")
    print()
    print("  The chain-length QSAR is the benchmark because it is the same")
    print("  kind of object: one or two descriptors, fitted, scored")
    print("  leave-one-out. It uses a number anyone can read off a formula;")
    print("  this model needs a measured binding constant per compound.")
    print()

    p = HERE / "db" / "tk_albumin_only.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["species"], x["source"],
                                             x["chemical"])):
            w.writerow(r)
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
