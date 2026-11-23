#!/usr/bin/env python3
"""Does excluding the weak bounds bias the fit? (Stage QSAR-2)

    python3 -I scripts/qsar2_weak_bound_sensitivity.py

Reads db/qsar/qsar_modelling_table.csv, writes db/qsar/qsar_weak_bound_
sensitivity.csv. Reuses the likelihood and fitter from qsar1. Idempotent.

THE QUESTION. qsar0 marks a row "usable_for_regression = no (bound permits
either sign)" when its upper bound on T sits above 1. For the VERDICT question
that rule is exactly right: a bound of T <= 1.29 establishes nothing about
direction. But qsar1 fits a CENSORED likelihood, and there a right-censored
observation at T <= 1.29 is perfectly good data -- it contributes
Phi((log10 1.29 - mu)/sigma), a real term. Dropping it is not required by the
statistics.

WORSE, THE DROPPING IS NOT RANDOM. A bound is weak when the filtered load is
small, and the filtered load is small when fu is small -- which happens at the
LONG-CHAIN end. So the exclusion removes long-chain rows preferentially, which
is selection on a variable correlated with the outcome, and that is the classic
way to manufacture a slope. This script measures the size of that effect
instead of assuming it is negligible.
"""
import csv, importlib.util, math, pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "q1", pathlib.Path(__file__).resolve().parent / "qsar1_censored_fit.py")
q1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q1)


def load_all():
    """Every row with a usable censoring type, INCLUDING the weak bounds."""
    rows = []
    f = HERE / "db" / "qsar" / "qsar_modelling_table.csv"
    for r in csv.DictReader(f.open()):
        if not r["T_point"]:
            continue
        cens = r["censoring"]
        kind = ("point" if cens.startswith("none") else
                "right" if cens.startswith("right") else
                "left" if cens.startswith("left") else
                "interval" if cens.startswith("interval") else None)
        if kind is None:
            continue
        lo = hi = None
        if kind == "interval":
            lo = math.log10(float(r["T_low"]))
            hi = math.log10(float(r["T_high"]))
        rows.append(dict(
            species=r["species"], sex=r["sex"], chem=r["chemical"],
            nC=float(r["n_C_total"]),
            sulf=1.0 if "sulfonate" in r["head_group"] else 0.0,
            head=r["head_group"], y=math.log10(float(r["T_point"])),
            lo=lo, hi=hi, kind=kind, verdict=r["verdict"],
            weak=r["usable_for_regression"].startswith("no")))
        rows[-1]["T_point"] = float(r["T_point"])
    return rows


def run(rows, label):
    """q1.fit returns a dict. No try/except here on purpose: a fit that
    fails should stop the script loudly, not vanish into a blank table."""
    f = q1.fit(rows)
    b, se = f["b"], f["se"]
    return dict(label=label, n=f["n"], sigma=f["sigma"], ll=f["ll"],
                aic=f["aic"], nC=b[1], nC_se=se[1],
                nC_z=b[1] / se[1] if se[1] else float("nan"),
                sulf=b[2], sulf_se=se[2])


def main():
    allrows = load_all()
    bar = "=" * 78
    print(bar)
    print("WHICH ROWS THE 'BOUND PERMITS EITHER SIGN' RULE DROPS")
    print(bar)
    weak = [r for r in allrows if r["weak"]]
    if not weak:
        print("  none -- every row has an informative bound.")
    print(f"  {'species':8}{'chem':9}{'nC':>3} {'sex':8}{'kind':10}{'T bound':>9}")
    for r in sorted(weak, key=lambda r: (r["species"], r["nC"], r["chem"])):
        print(f"  {r['species']:8}{r['chem']:9}{int(r['nC']):>3} {r['sex']:8}"
              f"{r['kind']:10}{r['T_point']:>9.3f}")
    print()
    bysp = {}
    for r in allrows:
        bysp.setdefault(r["species"], []).append(r)
    print(f"  {'species':8}{'all':>5}{'kept':>6}{'dropped':>9}"
          f"{'mean nC kept':>14}{'mean nC dropped':>17}")
    for sp, rs in sorted(bysp.items()):
        k = [r for r in rs if not r["weak"]]
        d = [r for r in rs if r["weak"]]
        mk = sum(r["nC"] for r in k) / len(k) if k else float("nan")
        md = sum(r["nC"] for r in d) / len(d) if d else float("nan")
        print(f"  {sp:8}{len(rs):>5}{len(k):>6}{len(d):>9}{mk:>14.2f}"
              f"{(f'{md:.2f}' if d else '-'):>17}")
    print()
    print("  IF THE DROPPED ROWS SIT AT LONGER CHAINS THAN THE KEPT ONES,")
    print("  the exclusion is removing exactly the observations that would")
    print("  flatten the chain-length slope, and the slope is overstated by")
    print("  however much the refit below moves.")
    print()

    res = []
    print(bar)
    print("REFIT WITH THE WEAK BOUNDS PUT BACK  (log10 T ~ nC + sulfonate)")
    print(bar)
    print(f"  {'species':8}{'fit':28}{'n':>4}{'sigma':>8}{'nC coef':>10}"
          f"{'SE':>8}{'z':>7}{'x/carbon':>10}")
    for sp, rs in sorted(bysp.items()):
        kept = [r for r in rs if not r["weak"]]
        a = run(kept, "excluding weak bounds (qsar1)")
        b = run(rs, "including weak bounds")
        for d in (a, b):
            if d is None:
                continue
            d["species"] = sp
            res.append(d)
            print(f"  {sp:8}{d['label']:28}{d['n']:>4}{d['sigma']:>8.3f}"
                  f"{d['nC']:>10.4f}{d['nC_se']:>8.4f}{d['nC_z']:>7.2f}"
                  f"{10**d['nC']:>10.3f}")
        if a and b:
            dz = b["nC_z"] - a["nC_z"]
            dc = b["nC"] - a["nC"]
            flag = ("UNCHANGED" if abs(dc) < 0.02 else
                    "MOVES" if abs(dc) < 0.08 else "MOVES MATERIALLY")
            print(f"  {'':8}{'-> delta':28}{'':4}{'':8}{dc:>10.4f}"
                  f"{'':8}{dz:>7.2f}   {flag}")
        print()

    print(bar)
    print("AND DOES THE CURVATURE CLAIM SURVIVE?  (nC^2 added)")
    print(bar)
    print("  The report rests its strongest shape claim on the RAT: an")
    print("  inverted U that survives dropping the digitised rows. If the")
    print("  weak bounds also bear on the curvature, that claim needs")
    print("  restating. Compare dAIC WITHIN a row only.")
    print()
    print(f"  {'species':8}{'fit':28}{'n':>4}{'AIC lin':>10}{'AIC quad':>10}"
          f"{'dAIC':>8}{'vertex':>9}")
    quad = []
    for sp, rs in sorted(bysp.items()):
        for label, sel in (("excluding weak bounds (qsar1)",
                            [r for r in rs if not r["weak"]]),
                           ("including weak bounds", rs)):
            fl, fq = q1.fit(sel), q1.fit(sel, quad=True)
            d = fq["aic"] - fl["aic"]
            bq = fq["b"]
            vx = -bq[1] / (2 * bq[3]) if bq[3] else float("nan")
            print(f"  {sp:8}{label:28}{len(sel):>4}{fl['aic']:>10.2f}"
                  f"{fq['aic']:>10.2f}{d:>8.2f}{vx:>9.1f}")
            quad.append(dict(species=sp, label=label, n=len(sel),
                             aic_lin=fl["aic"], aic_quad=fq["aic"],
                             daic=d, vertex=vx))
        print()
    print("  A dAIC that stays negative means the curvature is not an")
    print("  artefact of the exclusion rule.")
    print()

    pq = HERE / "db" / "qsar" / "qsar_weak_bound_curvature.csv"
    with pq.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "fit", "n", "aic_linear", "aic_quadratic",
                    "delta_aic", "vertex_nC"])
        for d in quad:
            w.writerow([d["species"], d["label"], d["n"],
                        round(d["aic_lin"], 3), round(d["aic_quad"], 3),
                        round(d["daic"], 3), round(d["vertex"], 2)])
    print(f"wrote {pq.relative_to(HERE)}  ({len(quad)} rows)")

    p = HERE / "db" / "qsar" / "qsar_weak_bound_sensitivity.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "fit", "n", "sigma", "loglik", "aic",
                    "nC_coef", "nC_se", "nC_z", "fold_per_carbon",
                    "sulfonate_coef", "sulfonate_se"])
        for d in res:
            w.writerow([d["species"], d["label"], d["n"], round(d["sigma"], 4),
                        round(d["ll"], 4), round(d["aic"], 3),
                        round(d["nC"], 5), round(d["nC_se"], 5),
                        round(d["nC_z"], 3), round(10 ** d["nC"], 5),
                        round(d["sulf"], 5), round(d["sulf_se"], 5)])
    print(f"wrote {p.relative_to(HERE)}  ({len(res)} rows)")


if __name__ == "__main__":
    main()
