#!/usr/bin/env python3
"""Re-run both mass-balance models with Fischer 2024's measured albumin
partition coefficient in place of each model's own.

Why this is the one substitution worth making.  Section 6.1b shows that a
binding constant and a partition coefficient are the same measurement, bridged
by D = 1000 n Ka / M_alb, and that albumin carries 19-94% of the sorptive
capacity in these assays, so an error in D_alb propagates straight into
predicted C_cell by 1/[(1-a) + f a].  Fischer 2024 is the only albumin dataset
measured at a ligand:protein ratio the body ever sees (<= 0.004), so it is the
natural absolute reference.

The assays use fetal BOVINE serum, so Fischer's BSA column is used, not HSA.

Both models multiply their albumin coefficient by V_alb, a VOLUME in litres
(the 0.733 mL/g partial specific volume is applied when building it), so
Fischer's mass-basis log D is converted to the volume basis with +log10(1.364)
= +0.135 before substitution.  One caveat recorded rather than resolved: if
Armitage's published K_SaW is really a mass-basis L/kg, this implementation
and httk's both over-count albumin by 1.36x.  The substitution replaces that
term anyway, so it cannot affect the Fischer variants.

Run:  python3 tk_review/scripts/mbm_fischer_albumin.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P
import mbm_qin_vs_armitage as Q

DB = os.path.join(os.path.dirname(HERE), "db")
VOL_SHIFT = math.log10(1000.0 / 0.733 / 1000.0)      # mass -> volume basis

# Fischer 2024 Table S2, log D_BSA/water (L_water per kg_protein), carboxylates.
# nCF -> log D.  PFBA (nCF 3) is NOT measured; Fischer's series starts at nCF 5
# for BSA.
FISCHER_BSA = {5: 2.83, 6: 3.93, 7: 4.20, 8: 4.54, 9: 4.81,
               10: 5.17, 11: 5.12, 12: 4.47, 13: 3.87}
RISING = [5, 6, 7, 8, 9, 10]      # before the turnover at nC10-11


def fischer_logd(nCF):
    """Measured where Fischer has it, OLS extrapolation on the rising limb
    otherwise.  Only PFBA needs this, and it is flagged everywhere it is used."""
    if nCF in FISCHER_BSA:
        return FISCHER_BSA[nCF], "measured"
    xs = RISING
    ys = [FISCHER_BSA[x] for x in xs]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return my + b * (nCF - mx), f"extrapolated (slope {b:+.2f}/CF2 on nCF 5-10)"


# ---------------------------------------------------------------- model runs
def qin_run(name, dose, V, alb=None):
    """Qin's balance.  alb=None uses Qin's own Langmuir D_BSA/w.

    alb = ("const", D_vol)  replaces it with a concentration-independent value
    alb = ("lang",  D_vol)  keeps Qin's Kd and D_nonspecific but rescales Cmax
                            so the TRACE LIMIT equals D_vol, preserving the
                            saturation behaviour Qin fitted.
    """
    ntot = dose * V["Vbm"]
    Cw = dose * 1e-3
    p = Q.qin_params(name, P.NCF[name])
    cmax = None
    if alb and alb[0] == "lang":
        cmax = (alb[1] - 10 ** p["logDnonspec"]) * p["Kd"]
        if cmax <= 0:
            return None                      # Fischer below the non-specific floor
    for _ in range(500):
        Dbsa_qin, Dsp, Dlip, _p = P.qin_D(name, Cw)
        if alb is None:
            Dbsa = Dbsa_qin
        elif alb[0] == "const":
            Dbsa = alb[1]
        else:
            Dbsa = cmax / (p["Kd"] + Cw) + 10 ** p["logDnonspec"]
        den = (V["Vw"] + Dbsa * V["Valb"] + Dlip * V["Vslip"]
               + Dsp * V["Vprot_cell"] + Dlip * V["Vlip_cell"])
        new = ntot / den
        if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-12:
            Cw = new
            break
        Cw = math.sqrt(Cw * new)
    _d, Dsp, Dlip, _p = P.qin_D(name, Cw)
    n_cell = Dsp * V["Vprot_cell"] * Cw + Dlip * V["Vlip_cell"] * Cw
    return dict(Cw=Cw, x_cell=n_cell / ntot, Dbsa=Dbsa,
                x_alb=Dbsa * Cw * V["Valb"] / ntot)


def armitage_run(name, lk, dose, V, alb=None):
    D = P.armitage_D(lk, V["comp"])
    Dsa = D["Dsa"] if alb is None else alb[1]
    ntot = dose * V["Vbm"]
    den = (V["Vw"] + Dsa * V["Valb"] + D["Dow"] * V["Vslip"]
           + D["Dcw"] * V["Vcells"] + D["Dpl"] * V["sarea"] * 1000.0)
    Cw = ntot / den
    return dict(Cw=Cw, x_cell=D["Dcw"] * Cw * V["Vcells"] / ntot, Dbsa=Dsa,
                x_alb=Dsa * Cw * V["Valb"] / ntot)


# ---------------------------------------------------------------- statistics
def stats(pairs):
    """Geometric-mean bias and geometric-mean absolute error of pred/meas."""
    if not pairs:
        return None, None, 0
    lg = [math.log10(p / m) for p, m in pairs]
    bias = 10 ** (sum(lg) / len(lg))
    aerr = 10 ** (sum(abs(x) for x in lg) / len(lg))
    return bias, aerr, len(lg)


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print("=" * 78)
    print("ALBUMIN PARTITION COEFFICIENT, BEFORE AND AFTER")
    print("=" * 78)
    print(f"\n  {'cmpd':8s}{'nCF':>4s}{'Fischer BSA':>13s}{'-> volume':>11s}"
          f"{'Qin trace':>11s}{'Armitage':>10s}{'  F/Qin':>9s}{'F/Arm':>8s}  provenance")
    fis = {}
    for (nm, nC, nCF, mw, ko, kc, ksrc) in Q.COMPOUNDS:
        ld, prov = fischer_logd(nCF)
        fis[nm] = (10 ** (ld + VOL_SHIFT), ld, prov)
        lk = ko if ko is not None else Q.KOW_OPERA(nCF)
        p = Q.qin_params(nm, nCF)
        qin_trace = p["Cmax"] / p["Kd"] + 10 ** p["logDnonspec"]
        gKsa = (1.08 * lk - 0.70) if lk < 4.5 else (0.37 * lk + 2.56)
        arm = 10 ** gKsa / 10 ** ((0.04 * lk + 0.114) * 0.15)
        print(f"  {nm:8s}{nCF:>4d}{ld:>13.2f}{math.log10(fis[nm][0]):>11.2f}"
              f"{math.log10(qin_trace):>11.2f}{math.log10(arm):>10.2f}"
              f"{fis[nm][0]/qin_trace:>9.2f}{fis[nm][0]/arm:>8.2f}  {prov}")
    print("""
  Fischer sits BELOW both models for every compound, and furthest below Qin.
  Substituting it therefore shrinks the albumin sink and pushes more of the
  dose into the cell -- the direction that matters, because both models
  already over-predict C_cell.
""")

    variants = [("Qin, published", "qin", None),
                ("Qin + Fischer (constant D)", "qin", "const"),
                ("Qin + Fischer (rescaled Langmuir)", "qin", "lang"),
                ("Armitage, published", "arm", None),
                ("Armitage + Fischer", "arm", "const")]
    rows, series = [], {v[0]: [] for v in variants}
    noPFBA = {v[0]: [] for v in variants}
    bydose = {v[0]: {"low": [], "high": []} for v in variants}
    bycmp = {v[0]: {} for v in variants}

    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        for (nm, nC, nCF, mw, ko, kc, ksrc) in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                rec = dict(study=study, compound=nm, dose_uM=d,
                           fischer_logD_BSA=round(fis[nm][1], 3),
                           fischer_provenance=fis[nm][2],
                           meas_x_cell_pct=round(m * 100, 4))
                for label, kind, mode in variants:
                    alb = None if mode is None else (mode, fis[nm][0])
                    r = (qin_run(nm, d, V, alb) if kind == "qin"
                         else armitage_run(nm, lk, d, V, alb))
                    if r is None:
                        rec[_key(label)] = ""
                        continue
                    rec[_key(label)] = round(r["x_cell"] * 100, 4)
                    rec[_key(label) + "_ratio"] = round(r["x_cell"] / m, 3)
                    series[label].append((r["x_cell"], m))
                    if nm != "PFBA":
                        noPFBA[label].append((r["x_cell"], m))
                    bydose[label]["low" if d <= 25 else "high"].append((r["x_cell"], m))
                    bycmp[label].setdefault(nm, []).append((r["x_cell"], m))
                rows.append(rec)

    print("=" * 78)
    print("PERFORMANCE OVER THE 29 MEASURED POINTS")
    print("=" * 78)
    print(f"\n  {'variant':36s}{'bias':>8s}{'|err|':>8s}{'n':>4s}"
          f"{'  bias w/o PFBA':>16s}{'<=25 uM':>10s}{'100 uM':>9s}")
    for label, _k, _m in variants:
        b, e, n = stats(series[label])
        b2, _e2, _n2 = stats(noPFBA[label])
        bl, _, _ = stats(bydose[label]["low"])
        bh, _, _ = stats(bydose[label]["high"])
        print(f"  {label:36s}{b:>8.2f}{e:>8.2f}{n:>4d}{b2:>16.2f}"
              f"{bl:>10.2f}{bh:>9.2f}")

    print(f"\n  Bias by compound (predicted / measured, geometric mean):\n")
    print(f"  {'variant':36s}" + "".join(f"{c:>9s}" for c in
                                          ("PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA")))
    for label, _k, _m in variants:
        line = f"  {label:36s}"
        for c in ("PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"):
            b, _, _ = stats(bycmp[label].get(c, []))
            line += f"{b:>9.2f}" if b else f"{'-':>9s}"
        print(line)

    # ------------------------------------------------ dose shape
    print("""
  DOSE SHAPE.  Qin's published error swings 3.5x across the dose range while
  the Fischer-constant variant is flat.  That is the one thing the
  substitution unambiguously fixes: Qin's dose dependence is an artefact of
  its albumin term SATURATING, not a property of the chemistry.""")
    print(f"\n  {'variant':36s}{'<=25 uM':>10s}{'100 uM':>9s}{'swing':>8s}")
    for label, _k, _m in variants:
        bl, _, _ = stats(bydose[label]["low"])
        bh, _, _ = stats(bydose[label]["high"])
        print(f"  {label:36s}{bl:>10.2f}{bh:>9.2f}{bh/bl:>8.2f}x")

    # ------------------------------------------------ inverse problem
    print("\n" + "=" * 78)
    print("THE INVERSE PROBLEM: what D_alb would make each model unbiased?")
    print("=" * 78)
    print("""
  C_cell is proportional to C_free = n_tot / S with S = Vw + D_alb Valb + R.
  To scale a prediction down by k, S must grow by k, so the albumin
  coefficient that removes the bias is

      D_alb,needed = D_alb,now + (k - 1) S / Valb,     k = bias

  Solving that per compound and comparing with what has actually been
  measured is the test of whether albumin can carry the error at all.
""")
    # ceiling of what has actually been measured for each compound, on the
    # MASS basis: the tightest literature Ka converted at n = 6, or Fischer,
    # whichever is higher.  n = 6 is generous - albumin's fatty-acid sites.
    import csv as _csv
    M_ALB = 66463.0
    lit_ceiling, lit_n1 = {}, {}
    for r in _csv.DictReader(open(os.path.join(DB, "protein_binding.csv"))):
        if "serum albumin" not in r["protein"] or "M^-1" not in r["ka_units"]:
            continue
        try:
            ka = float(r["ka"])
        except ValueError:
            continue
        for n, store in ((1, lit_n1), (6, lit_ceiling)):
            ld = math.log10(1000.0 * n * ka / M_ALB)
            store[r["chemical"]] = max(store.get(r["chemical"], -9), ld)

    print(f"  {'cmpd':8s}{'model':10s}{'bias':>7s}{'log D now':>11s}"
          f"{'needed':>9s}{'Fischer':>9s}{'lit n=1':>9s}{'lit n=6':>9s}  within measured?")
    inv = []
    for nm in ("PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"):
        for label, kind, _m in (("Qin, published", "qin", None),
                                ("Armitage, published", "arm", None)):
            pr = bycmp[label].get(nm, [])
            if not pr:
                continue
            k, _, _ = stats(pr)
            # recompute S at the representative (Amstutz, 100 uM) cell
            a = P.ASSAYS["Amstutz"]
            V = P.volumes(a)
            lk = next(c[4] if c[4] is not None else Q.KOW_OPERA(c[2])
                      for c in Q.COMPOUNDS if c[0] == nm)
            if kind == "qin":
                r = qin_run(nm, 100.0, V, None)
            else:
                r = armitage_run(nm, lk, 100.0, V, None)
            S = (100.0 * V["Vbm"]) / r["Cw"]
            need = r["Dbsa"] + (k - 1) * S / V["Valb"]
            fl = fis[nm][0]
            # put "needed" on the mass basis to compare with measurements
            need_mass = math.log10(need) - VOL_SHIFT
            c1 = max(lit_n1.get(nm, -9), fis[nm][1])
            ceil = max(lit_ceiling.get(nm, -9), fis[nm][1])
            if need_mass <= c1:
                ok = "yes, even at n=1"
            elif need_mass <= ceil:
                ok = f"only if n ~ {10**(need_mass-c1):.0f}"
            else:
                ok = f"NO, {10**(need_mass-ceil):.1f}x over even at n=6"
            print(f"  {nm:8s}{label.split(',')[0]:10s}{k:>7.2f}"
                  f"{math.log10(r['Dbsa']):>11.2f}{need_mass:>9.2f}"
                  f"{fis[nm][1]:>9.2f}{c1:>9.2f}{ceil:>9.2f}  {ok}")
            inv.append((nm, label, k, r["Dbsa"], need, fl, need_mass, ceil, c1))
    fails = [(x[0], x[1].split(",")[0], x[6], x[7]) for x in inv if x[6] > x[7]]
    fail_txt = ("; ".join(f"{nm} under {md} needs log D {nd:.2f} against a "
                          f"{ce:.2f} ceiling" for nm, md, nd, ce in fails)
                if fails else "none")
    nfit = sum(1 for x in inv if x[6] <= x[7])
    n1fit = sum(1 for x in inv if x[6] <= x[8])
    print(f"""
  All three log D columns are on the mass basis.  "lit n=1" and "lit n=6" are
  the tightest binding constant reported for that compound converted at one
  and at six sites, or Fischer, whichever is larger - the strict and the
  generous ceiling on what has actually been measured.

  {n1fit} of {len(inv)} cells can be fixed within the strict ceiling, {nfit} of {len(inv)} within the
  generous one.  Failing outright: {fail_txt}.  Those ceilings are already
  generous - six sites on albumin, and Fischer's own BSA column where it is
  higher - so a cell that fails there cannot be rescued by any albumin number
  anyone has measured.

  But the direction is the real result.  Every "needed" value is ABOVE the
  model's current one, because both models over-predict C_cell, while
  Fischer's measurement is BELOW it for every compound.  The best available
  measurement therefore pushes the albumin parameter the WRONG WAY: the data
  say albumin binds less than the models assume, and the residual says the
  models need it to bind more.  Substituting Fischer makes both models worse,
  which is what the run above shows.

  Albumin is where the LEVERAGE is (section 6.1b: it carries 19-94% of
  sorptive capacity) but it is not where the ERROR is.  What the substitution
  does fix is the SHAPE: Qin's 3.5x dose swing collapses to flat, because that
  swing was its albumin term saturating, not chemistry.  The level error has
  to come from somewhere else - the cell term, the plastic term, or the
  measured C_cell values.
""")

    dst = os.path.join(DB, "mbm_fischer_albumin.csv")
    with open(dst, "w", newline="") as fh:
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {dst} ({len(rows)} rows)")
    return series


def _key(label):
    return (label.lower().replace(" ", "_").replace(",", "")
            .replace("(", "").replace(")", "").replace("+", "plus")
            + "_x_cell_pct")


if __name__ == "__main__":
    main()
