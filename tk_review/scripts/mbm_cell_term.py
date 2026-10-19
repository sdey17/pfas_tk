#!/usr/bin/env python3
"""Probing the cell term of both mass-balance models.

Section 6.1c showed albumin cannot carry the models' over-prediction: the
albumin coefficient that would fix it is above anything measured, while the
best measurement is below what the models already use.  That leaves the cell
term, which is what this script takes apart.

The cell term is where the two models disagree most about MECHANISM while
agreeing most about VALUE, which is the signature of compensating errors:

  Qin       cell = protein + lipid only, no water phase.  D_SP/w is a MEASURED
            non-specific protein affinity, D_lip/w a measured lipid partition.
  Armitage  cell = storage lipid + membrane lipid + 3.5% of protein treated as
            octanol + WATER as a sorptive phase in its own right.

Qin also published MEASURED log D_cell/w for four compounds (Table 2), so
unlike albumin there is an internal reference that does not depend on
importing another laboratory's dataset.

Run:  python3 tk_review/scripts/mbm_cell_term.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P
import mbm_qin_vs_armitage as Q
import mbm_fischer_albumin as F   # stats(), and the same assay loop

DB = os.path.join(os.path.dirname(HERE), "db")
LINE = "=" * 78

# The three cell compositions in play, by VOLUME fraction.
COMPS = {
    "Armitage AREc2 (used in both runs)":
        dict(prot=0.0514, lip=0.00046 + 0.00414, water=0.944),
    "Qin HEK293H (Table A5)":
        dict(prot=0.0814, lip=0.00432, water=None),
    "project sheet, HepG2":
        dict(prot=0.164, lip=0.00870, water=None),
}
CELL = P.CELL_COMP["AREc2"]
F_AREC2 = CELL["prot"] + CELL["sl"] + CELL["ml"]

ORDER = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]


def qin_measured_protlip(name, nCF):
    """Qin Table 2 log D_cell/w, per litre of cell protein+lipid."""
    if name in Q.QIN_T2:
        return Q.QIN_T2[name][1], "measured"
    # carboxylate anchors PFBA (nCF 3) and PFOA (nCF 7)
    a, b = Q.QIN_T2["PFBA"][1], Q.QIN_T2["PFOA"][1]
    slope = (b - a) / 4.0
    return a + slope * (nCF - 3), f"interpolated on the measured carboxylate slope ({slope:+.3f}/CF2)"


def model_protlip(name, nCF, lk):
    """Both models' D_cell/w expressed per litre of cell protein+lipid, the
    one basis that does not depend on how much water the cell carries."""
    p = Q.qin_params(name, nCF)
    Dsp, Dlip = 10 ** p["logDnonspec"], 10 ** p["logDlipw"]
    fp = CELL["prot"] / F_AREC2
    fl = (CELL["sl"] + CELL["ml"]) / F_AREC2
    qin = fp * Dsp + fl * Dlip
    D = P.armitage_D(lk, CELL)
    arm_org = (CELL["sl"] * D["Dow"] + CELL["ml"] * D["Dmw"]
               + P.P_NLOM * CELL["prot"] * D["Dow"]) / F_AREC2
    return qin, arm_org, D


def main():
    print(__doc__.split("Run:")[0].rstrip())

    # ------------------------------------------------------------------ 1
    print(LINE)
    print("1.  WHOSE CELL?  Three compositions, spanning 3.1x")
    print(LINE)
    print(f"\n  {'source':36s}{'protein':>9s}{'lipid':>9s}{'prot+lip':>10s}{'  water'}")
    for k, v in COMPS.items():
        f = v["prot"] + v["lip"]
        print(f"  {k:36s}{v['prot']:>9.4f}{v['lip']:>9.5f}{f:>10.4f}"
              f"{('  %.3f' % v['water']) if v['water'] else '  not used'}")
    fs = [v["prot"] + v["lip"] for v in COMPS.values()]
    print(f"""
  prot+lip volume fraction spans {min(fs):.3f} to {max(fs):.3f}, a factor of {max(fs)/min(fs):.1f}.  A
  D_cell/w quoted per litre of WHOLE CELL therefore means three different
  things depending on whose cell is meant, and Qin's Table 2 is explicitly
  "[Lw/Lcell,prot+lip]".  Everything below is on the protein+lipid basis for
  that reason: it is the only one that does not silently carry a cell-water
  assumption.
""")

    # ------------------------------------------------------------------ 2
    print(LINE)
    print("2.  WHAT EACH MODEL'S CELL TERM IS MADE OF")
    print(LINE)
    print(f"\n  Armitage, D_cell/w = sl*Dow + ml*Dmw + 0.035*prot*Dow + water")
    print(f"  {'cmpd':8s}{'logKow':>8s}{'stor.lip':>10s}{'memb.lip':>10s}"
          f"{'protein':>9s}{'WATER':>8s}{'  log D (whole cell)':>20s}")
    for nm in ORDER:
        nCF, lk = _ncf_lk(nm)
        D = P.armitage_D(lk, CELL)
        t = [CELL["sl"] * D["Dow"], CELL["ml"] * D["Dmw"],
             P.P_NLOM * CELL["prot"] * D["Dow"], CELL["water"]]
        s = sum(t)
        print(f"  {nm:8s}{lk:>8.2f}" + "".join(f"{100*x/s:>9.1f}%" for x in t[:3])
              + f"{100*t[3]/s:>7.1f}%{math.log10(s):>20.2f}")
    print(f"\n  Qin, D_cell/w = D_nonspecific*prot + D_lip/w*lipid, no water")
    print(f"  {'cmpd':8s}{'protein':>10s}{'lipid':>9s}{'  log D (prot+lip basis)':>26s}")
    for nm in ORDER:
        nCF, lk = _ncf_lk(nm)
        p = Q.qin_params(nm, nCF)
        tp = 10 ** p["logDnonspec"] * CELL["prot"]
        tl = 10 ** p["logDlipw"] * (CELL["sl"] + CELL["ml"])
        q, _a, _D = model_protlip(nm, nCF, lk)
        print(f"  {nm:8s}{100*tp/(tp+tl):>9.1f}%{100*tl/(tp+tl):>8.1f}%{math.log10(q):>26.2f}")
    print("""
  The two models could hardly disagree more about what a cell is.  Armitage's
  cell is 79% WATER at PFBA and 99.7% MEMBRANE LIPID at PFNA -- its protein
  contributes 0.1-0.2% throughout, because only 3.5% of protein is treated as
  sorptive and that fraction is modelled as octanol.  Qin's cell is 92-95%
  PROTEIN at every chain length and has no water phase at all.
""")

    # ------------------------------------------------------------------ 3
    print(LINE)
    print("3.  AGAINST QIN'S OWN MEASUREMENT (Table 2), protein+lipid basis")
    print(LINE)
    print(f"\n  {'cmpd':8s}{'nCF':>4s}{'measured':>10s}{'Qin model':>11s}{'Armitage':>10s}"
          f"{'Qin/meas':>10s}{'Arm/meas':>10s}  provenance")
    rows = []
    for nm in ORDER + ["PFHxS", "PFOS"]:
        nCF, lk = _ncf_lk(nm)
        if nCF is None:
            continue
        meas, prov = qin_measured_protlip(nm, nCF)
        q, a, _D = model_protlip(nm, nCF, lk)
        print(f"  {nm:8s}{nCF:>4d}{meas:>10.2f}{math.log10(q):>11.2f}"
              f"{math.log10(a):>10.2f}{q/10**meas:>10.2f}{a/10**meas:>10.2f}  {prov}")
        rows.append((nm, nCF, meas, q, a, prov))
    print("""
  Both models are LOW at the short chain and HIGH at the long one, and they
  cross Qin's measurement between PFHxA and PFHpA.  That is not a level
  offset a scaling factor could absorb -- the two disagree about SLOPE.
""")
    carb = [r for r in rows if r[0] in ORDER]
    ms = _slope([(r[1], r[2]) for r in carb])
    qs = _slope([(r[1], math.log10(r[3])) for r in carb])
    as_ = _slope([(r[1], math.log10(r[4])) for r in carb])
    print(f"  decades of log D_cell/w per CF2, carboxylates nCF 3-8:")
    print(f"    measured (Qin Table 2, 2 anchors) {ms:+.3f}")
    print(f"    Qin model                         {qs:+.3f}   ({qs/ms:.1f}x too steep)")
    print(f"    Armitage model                    {as_:+.3f}   ({as_/ms:.1f}x too steep)")
    print(f"""
  Qin's cell-homogenate measurement rises {ms:.2f} decades per CF2 where both
  models build in {qs:.2f}-{as_:.2f}.  Do NOT read that as the models being wrong yet:
  it says one measurement disagrees with them, and section 6 brings a second,
  independent measurement to adjudicate.  Worth holding alongside: Ryu 2024's
  whole-TISSUE f_unbound falls 0.46-0.53 decades per CF2 (section 6.3a of the
  main report), which is the models' slope, not Qin's Table 2 slope.
""")

    # ------------------------------------------------------------------ 4
    print(LINE)
    print("4.  COMPENSATING ERRORS: same answer, opposite mechanism")
    print(LINE)
    print(f"\n  {'cmpd':8s}{'Qin logD':>10s}{'Arm logD':>10s}{'Qin/Arm':>9s}"
          f"{'   x_cell Qin/Arm':>18s}")
    ratios, xr = [], []
    Vref = P.volumes(P.ASSAYS["Amstutz"])
    for nm in ORDER:
        nCF, lk = _ncf_lk(nm)
        q, a, _D = model_protlip(nm, nCF, lk)
        xq = _run("qin", nm, lk, 100.0, Vref)
        xa = _run("arm", nm, lk, 100.0, Vref)
        ratios.append(q / a)
        xr.append(xq / xa)
        print(f"  {nm:8s}{math.log10(q):>10.2f}{math.log10(a):>10.2f}{q/a:>9.2f}"
              f"{xq/xa:>18.2f}")
    print(f"""
  The two models disagree about the CELL PARTITION by {min(ratios):.2f}-{max(ratios):.1f}x, yet their
  predicted FRACTION OF DOSE IN THE CELL agrees to {min(xr):.2f}-{max(xr):.2f}x.  The albumin
  term dominates the denominator and damps the difference, so agreement on the
  output is not agreement on the mechanism.  Qin's protein term and Armitage's
  membrane-lipid term happen to have similar magnitude over this window while
  being built from entirely different phases; that is coincidence, not
  corroboration, and any "two independent models agree" argument for the cell
  term is worthless.  Note where the coincidence is worst: at PFBA, where
  Armitage's cell is 79% water and Qin's is 92% protein, the two differ {max(ratios):.0f}x.
""")

    # ------------------------------------------------------------------ 5
    print(LINE)
    print("5.  SUBSTITUTION: Qin's measured D_cell/w into both models")
    print(LINE)
    variants = [("Qin, published", "qin", None),
                ("Qin + measured D_cell/w", "qin", "cell"),
                ("Armitage, published", "arm", None),
                ("Armitage + measured D_cell/w", "arm", "cell")]
    series = {v[0]: [] for v in variants}
    bydose = {v[0]: {"low": [], "high": []} for v in variants}
    bycmp = {v[0]: {} for v in variants}
    out = []
    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        for nm in ORDER:
            nCF, lk = _ncf_lk(nm)
            meas, prov = qin_measured_protlip(nm, nCF)
            Dcell_meas = 10 ** meas * F_AREC2        # -> whole-cell basis
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                rec = dict(study=study, compound=nm, dose_uM=d,
                           qin_measured_logDcellw_protlip=round(meas, 3),
                           provenance=prov, meas_x_cell_pct=round(m * 100, 4))
                for label, kind, mode in variants:
                    x = _run(kind, nm, lk, d, V, Dcell_meas if mode else None)
                    rec[_k(label)] = round(x * 100, 4)
                    rec[_k(label) + "_ratio"] = round(x / m, 3)
                    series[label].append((x, m))
                    bydose[label]["low" if d <= 25 else "high"].append((x, m))
                    bycmp[label].setdefault(nm, []).append((x, m))
                out.append(rec)
    print(f"\n  {'variant':32s}{'bias':>8s}{'|err|':>8s}{'<=25 uM':>10s}{'100 uM':>9s}{'swing':>8s}")
    for label, _k_, _m in variants:
        b, e, n = F.stats(series[label])
        bl, _, _ = F.stats(bydose[label]["low"])
        bh, _, _ = F.stats(bydose[label]["high"])
        print(f"  {label:32s}{b:>8.2f}{e:>8.2f}{bl:>10.2f}{bh:>9.2f}{bh/bl:>8.2f}x")
    print(f"\n  bias by compound:\n")
    print(f"  {'variant':32s}" + "".join(f"{c:>9s}" for c in ORDER))
    for label, _k_, _m in variants:
        line = f"  {label:32s}"
        for c in ORDER:
            b, _, _ = F.stats(bycmp[label].get(c, []))
            line += f"{b:>9.2f}" if b else f"{'-':>9s}"
        print(line)

    # ------------------------------------------------------------------ 6
    print("\n" + LINE)
    print("6.  THE INVERSE PROBLEM: what D_cell/w do the uptake data imply?")
    print(LINE)
    print("""
  x_cell = D_cell Vcells / S with S = R + D_cell Vcells, so the cell partition
  that reproduces a measured uptake is

      D_cell,needed = [x_meas / (1 - x_meas)] R / Vcells

  with R re-solved each iteration because Qin's albumin term depends on Cw.
  This is a third, independent estimate of the same quantity, and it uses the
  cellular-uptake measurements rather than Qin's cell-homogenate experiment.
""")
    print(f"  {'cmpd':8s}{'nCF':>4s}{'uptake-implied':>16s}{'Qin meas':>10s}"
          f"{'Qin model':>11s}{'Armitage':>10s}{'  implied/measured':>20s}")
    implied = []
    for nm in ORDER:
        nCF, lk = _ncf_lk(nm)
        vals = []
        for study, a in P.ASSAYS.items():
            V = P.volumes(a)
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                vals.append(_needed_dcell(nm, lk, d, V, m))
        if not vals:
            continue
        gm = 10 ** (sum(math.log10(v) for v in vals) / len(vals))
        gm_pl = gm / F_AREC2                    # whole cell -> prot+lip basis
        meas, _prov = qin_measured_protlip(nm, nCF)
        q, arm, _D = model_protlip(nm, nCF, lk)
        print(f"  {nm:8s}{nCF:>4d}{math.log10(gm_pl):>16.2f}{meas:>10.2f}"
              f"{math.log10(q):>11.2f}{math.log10(arm):>10.2f}"
              f"{gm_pl/10**meas:>20.2f}")
        implied.append((nm, nCF, math.log10(gm_pl), meas, math.log10(q), math.log10(arm)))
    si = _slope([(r[1], r[2]) for r in implied])
    sm = _slope([(r[1], r[3]) for r in implied])
    print(f"""
  decades per CF2:  uptake-implied {si:+.3f}   Qin's cell-homogenate measurement {sm:+.3f}
                    Qin model {_slope([(r[1], r[4]) for r in implied]):+.3f}   Armitage model {_slope([(r[1], r[5]) for r in implied]):+.3f}

  The uptake data and Qin's cell-homogenate measurement are both measurements
  of cell sorption and they do not agree.  Three of the four slopes cluster:
  uptake-implied {si:+.3f}, Qin model {_slope([(r[1], r[4]) for r in implied]):+.3f}, Armitage {_slope([(r[1], r[5]) for r in implied]):+.3f}, and Ryu 2024's
  whole-tissue f_unbound gives 0.46-0.53 on the same axis.  Qin's Table 2 sits
  alone at {sm:+.3f}, and it is {10**(implied[0][3]-implied[0][2]):.0f}x above the uptake requirement at PFBA.

  SO SECTION 3 INVERTS.  The models' steep chain-length dependence is
  corroborated by two independent lines -- the cellular-uptake data themselves
  and whole-tissue dialysis -- and contradicted only by Qin's Table 2.  The
  slope is not the defect.
""")

    # ---------------------------------------------- robustness of section 6
    print(LINE)
    print("6b. IS THE IMPLIED VALUE AN ARTEFACT OF THE ALBUMIN ASSUMPTION?")
    print(LINE)
    print("""
  The inversion above holds everything except the cell term fixed, so it
  inherits whatever albumin parameterisation is in the denominator.  Redoing
  it under all three (Qin's Langmuir, Armitage's spLFER, Fischer's measured
  value) on one assay, Amstutz at 100 uM:
""")
    import mbm_fischer_albumin as FA
    Vr = P.volumes(P.ASSAYS["Amstutz"])
    print(f"  {'cmpd':8s}{'Qin alb':>10s}{'Armitage alb':>14s}{'Fischer alb':>13s}{'spread':>9s}")
    cols = {"Qin alb": [], "Armitage alb": [], "Fischer alb": []}
    for nm in ORDER:
        nCF, lk = _ncf_lk(nm)
        m = P.measured_x_cell("Amstutz", nm, 100.0, Vr, P.ASSAYS["Amstutz"])
        if m is None:
            continue
        ld, _pv = FA.fischer_logd(nCF)
        vals = []
        for Dalb in (None, P.armitage_D(lk, Vr["comp"])["Dsa"], 10 ** (ld + FA.VOL_SHIFT)):
            vals.append(math.log10(_needed_dcell(nm, lk, 100.0, Vr, m, Dalb) / F_AREC2))
        for key, v in zip(cols, vals):
            cols[key].append((nCF, v))
        print(f"  {nm:8s}{vals[0]:>10.2f}{vals[1]:>14.2f}{vals[2]:>13.2f}"
              f"{10**(max(vals)-min(vals)):>8.1f}x")
    sl = {k: _slope(v) for k, v in cols.items()}
    print(f"""
  LEVEL moves by several fold with the albumin choice, so the absolute
  numbers in section 7 are conditional, not settled.  SLOPE does not:
  {', '.join(f'{k} {v:+.2f}' for k, v in sl.items())} decades per CF2, against Qin's
  Table 2 at {sm:+.2f}.  The section 6 conclusion rests on the slope and survives;
  the section 7 magnitudes should be read as "several-fold too high", not as
  the specific factors printed.
""")

    # ------------------------------------------------------------------ 7
    print(LINE)
    print("7.  WHERE THE ERROR ACTUALLY IS")
    print(LINE)
    print(f"\n  {'cmpd':8s}{'implied':>9s}{'Qin model':>11s}{'Qin/implied':>13s}"
          f"{'Armitage':>10s}{'Arm/implied':>13s}")
    for nm, nCF, imp, meas, qm, am in implied:
        print(f"  {nm:8s}{imp:>9.2f}{qm:>11.2f}{10**(qm-imp):>13.2f}"
              f"{am:>10.2f}{10**(am-imp):>13.2f}")
    qo = [10 ** (r[4] - r[2]) for r in implied]
    ao = [10 ** (r[5] - r[2]) for r in implied]
    print(f"""
  With the slope exonerated, what is left is a LEVEL error in the cell term:
  Qin's cell partition runs {min(qo):.1f}-{max(qo):.1f}x above what the uptake data require,
  Armitage's {min(ao):.2f}-{max(ao):.1f}x, both peaking in the middle of the series
  (PFHxA-PFHpA) rather than at an end -- the same shape as the overall bias in
  section 5 and in the albumin run.  Section 6b is the caveat on those
  specific factors: they move several-fold with the albumin assumption, so
  read them as "several-fold too high, worst mid-chain" rather than as
  calibrated numbers.

  That closes the loop opened in section 6.1c of the main report.  Albumin
  could not carry the over-prediction: the coefficient needed was above
  anything measured and the best measurement pointed the other way.  The cell
  term can carry it, it is the right size to carry it, and it is wrong in the
  right direction.  What it is NOT is a chain-length-dependence problem, so
  re-fitting the lipid or octanol descriptor will not fix it.  The candidates
  that remain are the sorptive capacity assigned to cell protein -- Qin's
  D_nonspecific and Armitage's "3.5% of protein behaves as octanol" are the
  two very different guesses at the same quantity -- and the cell volume the
  whole thing is normalised by, which section 1 shows is uncertain to 3.1x.
""")

    dst = os.path.join(DB, "mbm_cell_term.csv")
    keys = []
    for r in out:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    print(f"\n  wrote {dst} ({len(out)} rows)")
    return series


def _ncf_lk(nm):
    for (n, nC, nCF, mw, ko, kc, src) in Q.COMPOUNDS:
        if n == nm:
            return nCF, (ko if ko is not None else Q.KOW_OPERA(nCF))
    extra = {"PFHxS": (6, 5.17), "PFOS": (8, 6.43)}   # Qin Table 1 sulfonates
    return extra.get(nm, (None, None))


def _slope(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            / sum((x - mx) ** 2 for x in xs))


def _run(kind, nm, lk, dose, V, Dcell=None):
    """x_cell with the cell partition optionally overridden (whole-cell basis)."""
    ntot = dose * V["Vbm"]
    if kind == "arm":
        D = P.armitage_D(lk, V["comp"])
        Dcw = D["Dcw"] if Dcell is None else Dcell
        den = (V["Vw"] + D["Dsa"] * V["Valb"] + D["Dow"] * V["Vslip"]
               + Dcw * V["Vcells"] + D["Dpl"] * V["sarea"] * 1000.0)
        Cw = ntot / den
        return Dcw * Cw * V["Vcells"] / ntot
    Cw = dose * 1e-3
    for _ in range(500):
        Dbsa, Dsp, Dlip, _p = P.qin_D(nm, Cw)
        cell_cap = (Dcell * V["Vcells"] if Dcell is not None
                    else Dsp * V["Vprot_cell"] + Dlip * V["Vlip_cell"])
        den = V["Vw"] + Dbsa * V["Valb"] + Dlip * V["Vslip"] + cell_cap
        new = ntot / den
        if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-12:
            Cw = new
            break
        Cw = math.sqrt(Cw * new)
    Dbsa, Dsp, Dlip, _p = P.qin_D(nm, Cw)
    cell_cap = (Dcell * V["Vcells"] if Dcell is not None
                else Dsp * V["Vprot_cell"] + Dlip * V["Vlip_cell"])
    return cell_cap * Cw / ntot


def _needed_dcell(nm, lk, dose, V, x_meas, Dalb=None):
    """Cell partition (whole-cell basis) that reproduces x_meas exactly.

    Dalb=None uses Qin's own Langmuir albumin term; pass a value to hold the
    albumin coefficient fixed instead (used by the robustness check)."""
    D = None
    guess = 1.0
    for _ in range(200):
        ntot = dose * V["Vbm"]
        Cw = dose * 1e-3
        for _j in range(200):
            Dbsa, Dsp, Dlip, _p = P.qin_D(nm, Cw)
            if Dalb is not None:
                Dbsa = Dalb
            den = (V["Vw"] + Dbsa * V["Valb"] + Dlip * V["Vslip"]
                   + guess * V["Vcells"])
            new = ntot / den
            if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-12:
                Cw = new
                break
            Cw = math.sqrt(Cw * new)
        Dbsa, Dsp, Dlip, _p = P.qin_D(nm, Cw)
        if Dalb is not None:
            Dbsa = Dalb
        R = V["Vw"] + Dbsa * V["Valb"] + Dlip * V["Vslip"]
        nd = (x_meas / (1 - x_meas)) * R / V["Vcells"]
        if D is not None and abs(math.log10(nd / D)) < 1e-10:
            return nd
        D, guess = nd, nd
    return D


def _k(label):
    return (label.lower().replace(" ", "_").replace(",", "").replace("+", "plus")
            .replace("/", "_") + "_x_cell_pct")


if __name__ == "__main__":
    main()
