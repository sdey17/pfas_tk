#!/usr/bin/env python3
"""Two models, built and scored separately (Stage 27).

    python3 -I scripts/stage27_tk_tiers_and_qsar.py

Writes db/tk_tiers.csv, db/tk_tier_summary.csv, db/qsar_models.csv,
db/qsar_summary.csv, db/renal_fraction_conflict.csv. Idempotent.

MODEL 1, THE TK MODEL, USES EVERY MEASURED INPUT THERE IS. It is built as a
LADDER: each rung adds one more measured quantity and nothing else, so the
drop in error from one rung to the next is the value of that measurement.

    tier 1   f_u, V_d, GFR, T assumed 1      t½ = ln2 V_d / (f_u GFR)
    tier 2   + measured renal clearance      t½ = ln2 V_d / CL_R
    tier 3   + measured renal fraction       t½ = ln2 V_d f_renal / CL_R
    tier 4   + compilation total clearance   t½ = ln2 V_d / CL_tot

NO TIER HAS A FITTED PARAMETER. Every rung is a forward calculation from
measurements, which is what separates it from Model 2.

TIER 4 IS AN IDENTITY, NOT A PREDICTION, and the script proves it rather
than asserting it: CL in db/species_decomposition.csv was computed as
ln2 V_d / t½, so inverting it returns the half-life it came from. It is
carried here only to mark the ceiling, and is excluded from every score.

Tiers 2 and 3 escape that circle because the renal clearances are measured
against urine, not fitted to a serum decay curve: Kim 2016 in the rat
(dose/AUC, with urinary recovery), Harada 2005, Zhang 2013 and Fujii 2015
in the human. They are the only genuinely independent clearance inputs in
the compilation.

ONE CAVEAT THAT APPLIES TO EVERY TIER. V_d is taken as reported. Where a
source study derived it as V_d,area = CL/lambda_z it is not independent of
that study's own half-life, which would make the ladder partly circular at
every rung. The compilation does not record the derivation per study, so
this cannot be resolved here; it is stated, not corrected.

MODEL 2, THE QSAR, IS KEPT ENTIRELY SEPARATE. It sees no f_u, no V_d, no
GFR, no clearance -- only what can be read off a structural formula: chain
length and head group. It is FITTED, so it is scored leave-one-out, and the
comparison against Model 1 is therefore between a parameter-free
calculation and a calibrated regression. Those are different kinds of
object and the write-up says so.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
LN2 = math.log(2.0)

GFR = {"human": 2570.0, "rat": 14400.0, "mouse": 16700.0}
GFR_SRC = {"human": "Chiu 2022 / OEHHA", "rat": "Han 2012 (OEHHA primary)",
           "mouse": "Chiu 2022 / OEHHA"}
SULFONATES = {"PFBS", "PFHxS", "PFOS"}


def sig(x, n=6):
    """Round to n significant figures; fixed decimals lose small values."""
    if x is None or x == 0:
        return 0.0 if x == 0 else None
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


def gmean(vs):
    return 10 ** (sum(math.log10(v) for v in vs) / len(vs))


# --------------------------------------------------------------- inputs
def load_tk():
    """{species: {chemical: dict}} with t½, V_d, f_u, CL from compilation."""
    out = {"human": {}, "rat": {}, "mouse": {}}
    for r in csv.DictReader((HERE / "db" / "species_decomposition.csv").open()):
        c, sp = r["chemical"], r["comparison"].split("/")[1]
        out["human"].setdefault(c, dict(
            t_half=float(r["t_half_human_d"]), vd=float(r["vd_human"]),
            fu=float(r["fu_human"]), cl=float(r["cl_human"]),
            nc=int(r["n_total_carbons"])))
        out[sp][c] = dict(
            t_half=float(r["t_half_animal_d"]), vd=float(r["vd_animal"]),
            fu=float(r["fu_animal"]), cl=float(r["cl_animal"]),
            nc=int(r["n_total_carbons"]), source=r["animal_source"])
    return out


def load_renal():
    """{(species, chemical): [(CL_R, source), ...]} -- per study, unpooled."""
    acc = {}
    for r in csv.DictReader(
            (HERE / "db" / "primary_2026"
             / "human_renal_clearance_primary.csv").open()):
        acc.setdefault(("human", r["chemical"]), []).append(
            (float(r["cl_renal_mL_kg_day"]), r["source"]))
    for r in csv.DictReader(
            (HERE / "db" / "primary_2026" / "kim2016_table4_rat.csv").open()):
        if r["sex"] != "male":            # male only: the paper's own arms
            continue                      # differ 75-fold by sex in PFOA
        acc.setdefault(("rat", r["chemical"]), []).append(
            (float(r["cl_renal_mL_kg_day"]), "Kim 2016 (male)"))
    return acc


def load_renal_fraction():
    """{(species, chemical): [(fraction, source), ...]}"""
    out = {}
    for r in csv.DictReader(
            (HERE / "db" / "human_renal_clearance_resolved.csv").open()):
        out.setdefault((r["species"], r["chemical"]), []).append(
            (float(r["renal_fraction_pct"]) / 100.0, r["source"]))
    # Kim 2016 reports CL_R and CL_tot in the same table, so its renal
    # fraction is internal to one study rather than assembled across two.
    for r in csv.DictReader(
            (HERE / "db" / "primary_2026" / "kim2016_table4_rat.csv").open()):
        if r["sex"] != "male":
            continue
        out.setdefault(("rat", r["chemical"]), []).append(
            (float(r["cl_renal_mL_kg_day"]) / float(r["cl_total_mL_kg_day"]),
             "Kim 2016 Table 4 (male), CL_R/CL_tot"))
    return out


def screen_renal(tk, renal):
    """Drop study values that exceed the CL_tot they would be paired with.

    CL_R is a COMPONENT of CL_tot, so CL_R > CL_tot is not a large error --
    it is an impossible pairing. Which side is wrong cannot be settled here:
    Zhang 2013 reports both a high renal clearance AND a short half-life for
    PFOA in young women, so its CL_R is self-consistent with ITS OWN
    half-life and inconsistent only with the compilation's. Either way the
    two cannot be multiplied together, so the pairing is void and recorded.
    """
    kept, audit = {}, []
    for (sp, c), vs in sorted(renal.items()):
        if c not in tk.get(sp, {}):
            continue
        cl_tot = tk[sp][c]["cl"]
        ok = []
        for v, src in sorted(vs):
            ratio = v / cl_tot
            good = ratio <= 1.0
            audit.append(dict(
                species=sp, chemical=c, cl_renal_mL_kg_day=v,
                cl_total_implied_mL_kg_day=sig(cl_tot),
                ratio_renal_to_total=sig(ratio, 4),
                verdict="usable" if good else "VOID: CL_R exceeds CL_tot",
                source=src))
            if good:
                ok.append((v, src))
        if ok:
            kept[(sp, c)] = (gmean([v for v, _ in ok]), len(ok),
                             min(v for v, _ in ok), max(v for v, _ in ok),
                             "; ".join(sorted({s for _, s in ok})))
    out = {"human": {}, "rat": {}}
    for (sp, c), v in kept.items():
        out[sp][c] = v
    return out, audit


# --------------------------------------------------------------- algebra
def ols(rows, y):
    """Least squares for y = X beta. rows include the intercept column."""
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


def design(form, nc, sulf):
    x = [1.0, float(nc)]
    if "nC^2" in form:
        x.append(float(nc) ** 2)
    if "head" in form:
        x.append(1.0 if sulf else 0.0)
    return x


QSAR_FORMS = ["log t½ ~ nC", "log t½ ~ nC + nC^2", "log t½ ~ nC + head",
              "log t½ ~ nC + nC^2 + head"]


# --------------------------------------------------------------- model 1
T1 = "tier 1: f_u + V_d + GFR, T=1"
T2 = "tier 2: + measured CL_R"
T3 = "tier 3: + measured renal fraction"
T4 = "tier 4: + compilation CL_tot"


def tk_ladder(tk, renal, frac, excluded):
    rows = []
    for sp in ("human", "rat", "mouse"):
        for c in sorted(tk[sp]):
            d = tk[sp][c]
            obs = d["t_half"]
            cand = [(T1, LN2 * d["vd"] / (d["fu"] * GFR[sp]), "prediction",
                     "", f"GFR {GFR_SRC[sp]}")]
            rr = renal.get(sp, {}).get(c)
            if rr:
                cand.append((T2, LN2 * d["vd"] / rr[0], "prediction", "",
                             f"CL_R {rr[4]} (n={rr[1]})"))
                for fr, src in frac.get((sp, c), []):
                    # The rodent arms are male throughout (see animal_source);
                    # Harada's own male and female rat PFOA clearances differ
                    # 75-fold, so a female fraction cannot be carried onto a
                    # male half-life. Recorded as excluded, not silently cut.
                    if sp != "human" and "female" in src.lower():
                        excluded.append(dict(
                            species=sp, chemical=c,
                            renal_fraction_pct=sig(fr * 100, 4), source=src,
                            reason="female fraction; t½ and CL_R are male"))
                        continue
                    cand.append((T3, LN2 * d["vd"] * fr / rr[0], "prediction",
                                 sig(fr * 100, 4), src))
            cand.append((T4, LN2 * d["vd"] / d["cl"], "IDENTITY", "",
                         "circular"))
            for name, pred, kind, fpct, prov in cand:
                rows.append(dict(
                    species=sp, chemical=c, tier=name, kind=kind,
                    renal_fraction_pct=fpct,
                    t_half_measured_d=sig(obs),
                    t_half_predicted_d=sig(pred),
                    fold_error=sig(max(pred / obs, obs / pred), 4),
                    direction="over" if pred > obs else "under",
                    provenance=prov))
    return rows


# --------------------------------------------------------------- model 2
def qsar(tk):
    rows, summary = [], []
    for sp in ("human", "rat", "mouse"):
        chems = sorted(tk[sp])
        for form in QSAR_FORMS:
            X = [design(form, tk[sp][c]["nc"], c in SULFONATES) for c in chems]
            y = [math.log10(tk[sp][c]["t_half"]) for c in chems]
            p = len(X[0])
            if len(chems) <= p:
                continue
            full = ols(X, y)
            if full is None:
                continue
            resid = [y[i] - sum(a * b for a, b in zip(full, X[i]))
                     for i in range(len(y))]
            my = sum(y) / len(y)
            sst = sum((v - my) ** 2 for v in y)
            r2 = 1 - sum(e ** 2 for e in resid) / sst if sst else float("nan")
            folds = []
            for i in range(len(chems)):
                Xi = [X[j] for j in range(len(X)) if j != i]
                yi = [y[j] for j in range(len(y)) if j != i]
                if len(Xi) <= p:
                    continue
                b = ols(Xi, yi)
                if b is None:
                    continue
                pr = sum(a * c2 for a, c2 in zip(b, X[i]))
                fold = 10 ** abs(pr - y[i])
                folds.append(fold)
                rows.append(dict(
                    species=sp, form=form, chemical=chems[i],
                    n_total_carbons=tk[sp][chems[i]]["nc"],
                    head_group="sulfonate" if chems[i] in SULFONATES
                    else "carboxylate",
                    t_half_measured_d=sig(tk[sp][chems[i]]["t_half"]),
                    t_half_loo_predicted_d=sig(10 ** pr),
                    loo_fold_error=sig(fold, 4)))
            if not folds:
                continue
            folds.sort()
            summary.append(dict(
                species=sp, form=form, n=len(chems), n_params=p,
                r2_in_sample=sig(r2, 3),
                loo_median_fold=sig(folds[len(folds) // 2], 4),
                loo_worst_fold=sig(folds[-1], 4)))
    return rows, summary


# --------------------------------------------------------------- report
def main():
    tk, frac = load_tk(), load_renal_fraction()
    renal, rscreen = screen_renal(tk, load_renal())
    excluded = []
    bar = "=" * 78

    bad = [r for r in rscreen if r["verdict"] != "usable"]
    print(bar)
    print("BEFORE ANY MODEL: SCREENING THE MEASURED RENAL CLEARANCES")
    print(bar)
    print("  CL_R is one component of CL_tot, so CL_R > CL_tot is impossible,")
    print("  not merely discrepant. Each measured value is therefore checked")
    print("  against the total clearance it would be paired with.")
    print()
    print(f"  {len(bad)} of {len(rscreen)} values FAIL:")
    for r in bad:
        print(f"    {r['species']} {r['chemical']:8} "
              f"CL_R {r['cl_renal_mL_kg_day']:>6} vs CL_tot "
              f"{r['cl_total_implied_mL_kg_day']:<8} "
              f"= {r['ratio_renal_to_total']:.2f}x   {r['source']}")
    print()
    print("  These are not dropped for being inconvenient. A renal clearance")
    print("  above total clearance means the two numbers come from")
    print("  populations with different half-lives -- Zhang 2013 pairs its")
    print("  high PFOA clearance with a ~1.2-year half-life in young women,")
    print("  against the ~2-4 years in the compilation -- so the product of")
    print("  one study's CL_R and another's V_d has no referent. WHICH SIDE")
    print("  IS WRONG IS NOT SETTLED HERE; only the pairing is void.")
    print()
    print("  ONE EXCEPTION, DELIBERATE. Tier 3 still uses Fujii's PFOA renal")
    print("  FRACTION (0.044/0.096 = 45.8%) even though its absolute CL_R was")
    print("  just voided. A fraction is dimensionless and internally")
    print("  consistent within one study, so it survives a disagreement about")
    print("  the absolute scale that the raw clearance does not.")
    print()

    print(bar)
    print("MODEL 1  THE TK LADDER -- every measured input, no fitted parameter")
    print(bar)
    print("  Each rung adds ONE measurement. The fall in error from one rung")
    print("  to the next is what that measurement is worth.")
    print()

    tiers = tk_ladder(tk, renal, frac, excluded)

    # the identity, demonstrated rather than asserted
    ident = [r for r in tiers if r["kind"] == "IDENTITY"]
    worst_id = max(r["fold_error"] for r in ident)
    print("  FIRST, TIER 4 IS DISQUALIFIED, AND HERE IS THE PROOF.")
    print(f"    {len(ident)} compound-species rows, worst fold error "
          f"{worst_id:.4f}x.")
    print("    A model that reproduces every half-life to the fourth decimal")
    print("    has not predicted them; CL in the compilation WAS ln2 V_d/t½,")
    print("    so tier 4 inverts the arithmetic that produced its own target.")
    print("    It is excluded from every score below.")
    print()

    print(f"  {'tier':40}{'n':>4}{'median':>10}{'worst':>10}  direction")
    rank = []
    for name in (T1, T2, T3):
        grp = [r for r in tiers if r["tier"] == name]
        if not grp:
            continue
        f = sorted(r["fold_error"] for r in grp)
        dirs = {r["direction"] for r in grp}
        rank.append((name, len(grp), f[len(f) // 2], f[-1]))
        print(f"  {name:40}{len(grp):>4}{f[len(f)//2]:>9.1f}x{f[-1]:>9.1f}x  "
              f"{'/'.join(sorted(dirs))}")
    print()
    if excluded:
        print("  EXCLUDED from tier 3, with the reason recorded:")
        for e in excluded:
            print(f"    {e['species']} {e['chemical']} "
                  f"{e['renal_fraction_pct']}%  -- {e['reason']}")
    print()

    # what each measurement buys, per species, on the compounds that have it
    print(bar)
    print("WHAT EACH MEASUREMENT BUYS (same compounds throughout)")
    print(bar)
    for sp in ("human", "rat"):
        have = sorted({r["chemical"] for r in tiers
                       if r["species"] == sp
                       and r["tier"] == T2})
        if not have:
            continue
        print(f"  {sp.upper()}  on the {len(have)} compounds with a measured "
              f"CL_R: {', '.join(have)}")
        for pref in (T1, T2, T3):
            grp = [r for r in tiers if r["species"] == sp
                   and r["chemical"] in have and r["tier"] == pref]
            if not grp:
                print(f"    {pref.split(':')[0]:8} unavailable")
                continue
            f = sorted(r["fold_error"] for r in grp)
            over = sum(1 for r in grp if r["direction"] == "over")
            print(f"    {pref.split(':')[0]:8} median {f[len(f)//2]:7.1f}x"
                  f"   worst {f[-1]:7.1f}x   ({over}/{len(grp)} over-predict)")
        print()

    print(bar)
    print("AUDIT: IS TIER 2's ERROR REALLY 1/f_renal?")
    print(bar)
    print("  The claim is arithmetic, so it is testable. If tier 2 omits every")
    print("  non-renal route then its fold error must equal CL_tot/CL_R. Both")
    print("  sides are computed independently below.")
    print()
    print(f"  {'':6}{'':9}{'tier 2 fold':>13}{'CL_tot/CL_R':>13}{'agree':>8}")
    bad = 0
    for r in sorted((x for x in tiers if x["tier"] == T2),
                    key=lambda x: (x["species"], x["chemical"])):
        sp, c = r["species"], r["chemical"]
        raw = tk[sp][c]["cl"] / renal[sp][c][0]
        implied = max(raw, 1.0 / raw)   # fold_error is symmetrised; match it
        ok = abs(implied - r["fold_error"]) / r["fold_error"] < 0.02
        bad += 0 if ok else 1
        print(f"  {sp:6}{c:9}{r['fold_error']:>12.2f}x{implied:>12.2f}x"
              f"{'yes' if ok else 'NO':>8}")
    print()
    if bad:
        print(f"  {bad} row(s) DISAGREE -- the tier-2 reading is wrong there.")
    else:
        print("  Every row agrees. Tier 2's error is not a modelling failure;")
        print("  it IS the reciprocal of the renal share, measured two ways.")
    print()
    o = max((x for x in tiers if x["tier"] == T2),
            key=lambda x: x["fold_error"])
    print(f"  The worst row, {o['species']} {o['chemical']} at "
          f"{o['fold_error']:.0f}x, is therefore not an outlier in the model.")
    print(f"  It says that compound's renal clearance is "
          f"1/{o['fold_error']:.0f} of its total.")
    print()
    print("  READING TIER 2. It uses a renal clearance and nothing about the")
    print("  other routes, so it must OVER-predict the half-life by exactly")
    print("  1/f_renal. The over-prediction is not an error in the model; it")
    print("  is a measurement of how much elimination is non-renal.")
    print()

    # ---- the rat renal fraction conflict, which tier 3 exposes
    print(bar)
    print("A CONFLICT TIER 3 FORCES INTO THE OPEN: THE RAT RENAL FRACTION")
    print(bar)
    conflict = []
    for (sp, c), vs in sorted(frac.items()):
        if len(vs) < 2:
            continue
        lo, hi = min(v for v, _ in vs), max(v for v, _ in vs)
        conflict.append(dict(
            species=sp, chemical=c, n_sources=len(vs),
            lowest_pct=sig(lo * 100, 4), highest_pct=sig(hi * 100, 4),
            spread_fold=sig(hi / lo, 4),
            sources=" | ".join(f"{v*100:.1f}% {s}" for v, s in sorted(vs))))
        print(f"  {sp} {c}: {lo*100:.1f}% to {hi*100:.1f}%  "
              f"({hi/lo:.1f}-fold apart)")
        for v, s in sorted(vs):
            print(f"      {v*100:6.1f}%   {s}")
    print()
    print("  The rat PFOA figure is the one that matters. Harada 2005 Table 5")
    print("  has renal clearance accounting for essentially all of it, which")
    print("  is the basis for the familiar claim that the rat clears PFOA")
    print("  renally while the human does not. Kim 2016 reports CL_R and")
    print("  CL_tot IN THE SAME TABLE, on the same animals, and gets 9%.")
    print("  Both cannot be right, and tier 3 for the rat swings by that")
    print("  same factor depending on which is used. This is a disagreement")
    print("  between two primary sources, not a modelling choice, and it is")
    print("  recorded rather than resolved.")
    print()

    # --------------------------------------------------------- model 2
    print(bar)
    print("MODEL 2  THE QSAR -- structure only, kept entirely separate")
    print(bar)
    print("  Inputs: chain length and head group. NOTHING ELSE -- no f_u, no")
    print("  V_d, no GFR, no clearance. It is fitted, so every number below")
    print("  is leave-one-out: each compound predicted from a fit without it.")
    print()
    qrows, qsum = qsar(tk)
    print(f"  {'species':8}{'form':30}{'n':>3}{'k':>3}{'r2':>7}"
          f"{'LOO med':>10}{'LOO worst':>11}")
    for s in qsum:
        print(f"  {s['species']:8}{s['form']:30}{s['n']:>3}{s['n_params']:>3}"
              f"{s['r2_in_sample']:>7.3f}{s['loo_median_fold']:>9.1f}x"
              f"{s['loo_worst_fold']:>10.1f}x")
    print()
    print("  THE OVERFITTING IS VISIBLE IN THE TABLE. Read r2 against LOO")
    print("  across the forms within a species: adding terms always raises")
    print("  r2, because it must, while LOO error can rise at the same time.")
    print("  With 8 or 9 compounds and up to 4 parameters that is the")
    print("  expected failure, and it is the reason for scoring LOO at all.")
    print()
    for sp in ("human", "rat", "mouse"):
        ss = [s for s in qsum if s["species"] == sp]
        if len(ss) < 2:
            continue
        simplest, best = ss[0], min(ss, key=lambda s: s["loo_median_fold"])
        richest = max(ss, key=lambda s: s["n_params"])
        print(f"    {sp:6} r2 climbs {simplest['r2_in_sample']:.3f} -> "
              f"{richest['r2_in_sample']:.3f} as k goes "
              f"{simplest['n_params']} -> {richest['n_params']}, "
              f"LOO median {simplest['loo_median_fold']:.1f}x -> "
              f"{richest['loo_median_fold']:.1f}x")
        print(f"    {'':6} best LOO is \"{best['form']}\" at "
              f"{best['loo_median_fold']:.1f}x")
    print()

    # --------------------------------------------------------- comparison
    print(bar)
    print("THE TWO MODELS SIDE BY SIDE")
    print(bar)
    print("  Scored on the SAME compounds -- the ones with a measured CL_R,")
    print("  so that tier 2 can appear at all.")
    print()
    for sp in ("human", "rat"):
        have = sorted({r["chemical"] for r in tiers if r["species"] == sp
                       and r["tier"] == T2})
        if not have:
            continue
        print(f"  {sp.upper()}  ({len(have)} compounds)")
        print(f"    {'model':46}{'median':>9}{'worst':>9}{'params':>8}")
        for pref in (T1, T2, T3):
            grp = [r for r in tiers if r["species"] == sp
                   and r["chemical"] in have and r["tier"] == pref]
            if not grp:
                continue
            f = sorted(r["fold_error"] for r in grp)
            lbl = "TK " + pref.split(":")[0] + " (parameter-free)"
            print(f"    {lbl:46}"
                  f"{f[len(f)//2]:>8.1f}x{f[-1]:>8.1f}x{0:>8}")
        for s in qsum:
            if s["species"] != sp:
                continue
            g = [r for r in qrows if r["species"] == sp
                 and r["form"] == s["form"] and r["chemical"] in have]
            if not g:
                continue
            f = sorted(r["loo_fold_error"] for r in g)
            print(f"    {'QSAR ' + s['form'] + ' (fitted, LOO)':46}"
                  f"{f[len(f)//2]:>8.1f}x{f[-1]:>8.1f}x"
                  f"{s['n_params']:>8}")
        print()
    print("  THE COMPARISON IS NOT LIKE FOR LIKE AND SHOULD NOT BE READ AS")
    print("  ONE. The TK tiers spend no degrees of freedom and need three to")
    print("  five measured inputs per compound; the QSAR spends two to four")
    print("  and needs a structural formula. A QSAR that wins on fold error")
    print("  has still been calibrated on the very half-lives it predicts,")
    print("  so it cannot be carried to a species whose half-lives are the")
    print("  unknown -- which is the usual reason for wanting the model.")
    print()

    # --------------------------------------------------------- write
    out = [
        ("tk_tiers.csv", tiers,
         ["species", "chemical", "tier", "kind", "renal_fraction_pct",
          "t_half_measured_d", "t_half_predicted_d", "fold_error",
          "direction", "provenance"]),
        ("tk_tier3_excluded.csv", excluded,
         ["species", "chemical", "renal_fraction_pct", "source", "reason"]),
        ("tk_tier_summary.csv",
         [dict(tier=n, n=k, median_fold=sig(m, 4), worst_fold=sig(w, 4))
          for n, k, m, w in rank],
         ["tier", "n", "median_fold", "worst_fold"]),
        ("qsar_models.csv", qrows,
         ["species", "form", "chemical", "n_total_carbons", "head_group",
          "t_half_measured_d", "t_half_loo_predicted_d", "loo_fold_error"]),
        ("qsar_summary.csv", qsum,
         ["species", "form", "n", "n_params", "r2_in_sample",
          "loo_median_fold", "loo_worst_fold"]),
        ("renal_clearance_consistency.csv", rscreen,
         ["species", "chemical", "cl_renal_mL_kg_day",
          "cl_total_implied_mL_kg_day", "ratio_renal_to_total", "verdict",
          "source"]),
        ("renal_fraction_conflict.csv", conflict,
         ["species", "chemical", "n_sources", "lowest_pct", "highest_pct",
          "spread_fold", "sources"]),
    ]
    for name, rows, cols in out:
        p = HERE / "db" / name
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in sorted(rows, key=lambda d: [str(d.get(c, "")) for c in cols]):
                w.writerow({c: r.get(c, "") for c in cols})
        print(f"wrote db/{name}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
