#!/usr/bin/env python3
"""A third parameterisation: measured biopartitioning instead of spLFERs.

WHY.  The two models tested in the memo bracket the field in one respect --
predicted-from-logKow (Armitage) against measured-and-saturable (Qin) -- but
not in another.  The Escher-school balances (Fischer et al. 2017; Huchthausen
et al. 2020; Qin et al. 2023) share a different design decision: they take
membrane-water and protein-water partitioning from EXPERIMENT rather than from
a single-descriptor regression, while keeping a linear (non-saturable)
treatment.  Armitage et al. report that this substitution substantially
improves their own model's performance for ionisable organics.

That variant can be built here without implementing anyone else's code,
because the measured inputs now exist for PFAS:

    albumin        Fischer et al. 2024, Table S2, log D_BSA/w   (4 of 5 compounds)
    membrane lipid Chen et al. 2025, log K_MW on supported PC bilayers

Keeping Armitage's cell composition and plastic term, and substituting both
measured series, gives a "measured-parameter linear balance" -- structurally
the Escher-school design, parameterised from the best PFAS data available.

BASIS (resolved).  Chen's SI eq S1 defines K_MW as
(C_spike - C_w)*V_w / (C_w * M_PL), with a unit conversion applied so that the
reported coefficient is dimensionless.  The TRANSIL SSLM kit supplies the
membrane as a stated VOLUME (0.067-2.166 uL of bilayer per 500 uL tube), so
the lipid quantity is natively volumetric and the mass and volume bases
coincide at rho_lipid = 1 g/cm3 -- the same convention V_cells uses here
(10 ng/cell at 1 mg/uL).  Real phospholipid densities of 1.01-1.06 bound the
residual basis error at <=0.03 log units.  So log K_MW is directly usable as
Armitage's D_MW.  The one live caveat is pH: Chen reports at pH 7.0, the
balance runs at 7.4.  PFBA has no Fischer BSA
value (N.D.) and is excluded, which removes the compound with the most
leverage on any chain-length fit -- so slopes here are over nCF 5-8 only and
are directly comparable with the nCF 5-8 columns in the memo.

Run:  python3 tk_review/scripts/mbm_measured_parameter_variant.py
"""
import collections
import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P      # noqa: E402
import mbm_qin_vs_armitage as Q          # noqa: E402

DB = os.path.join(os.path.dirname(HERE), "db")
NCF = {"PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}

# Fischer et al. 2024 Table S2, log D_BSA/w (L_water per kg BSA), by nCF.
FISCHER_BSA = {5: 2.83, 6: 3.93, 7: 4.20, 8: 4.54}
VOL_SHIFT = math.log10(1.0 / 0.733)      # mass -> volume basis, as in the memo


def chen_kmw():
    out = {}
    for r in csv.DictReader(open(os.path.join(DB, "protein_binding.csv"))):
        if "Chen" in r["study"] and "log K_MW" in r["ka_units"]:
            out[r["chemical"]] = float(r["ka"])
    return out


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def run(variant, kmw):
    """variant: 'published' | 'measured_alb' | 'measured_mem' | 'measured_both'"""
    per = collections.defaultdict(list)
    rows = []
    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
            if nm not in NCF:
                continue
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            D = dict(P.armitage_D(lk, V["comp"]))
            if variant in ("measured_alb", "measured_both"):
                D["Dsa"] = 10 ** (FISCHER_BSA[nCF] + VOL_SHIFT)
            if variant in ("measured_mem", "measured_both"):
                D["Dmw"] = 10 ** kmw[nm]
                c = V["comp"]
                D["Dcw"] = (c["sl"] * D["Dow"] + c["ml"] * D["Dmw"]
                            + P.P_NLOM * c["prot"] * D["Dow"] + c["water"])
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                ntot = d * V["Vbm"]
                S = (V["Vw"] + D["Dsa"] * V["Valb"] + D["Dow"] * V["Vslip"]
                     + D["Dcw"] * V["Vcells"] + D["Dpl"] * V["sarea"] * 1000.0)
                x = D["Dcw"] * (ntot / S) * V["Vcells"] / ntot
                per[nm].append(x / m)
                rows.append(dict(variant=variant, study=study, compound=nm,
                                 dose_uM=d, meas_pct=round(m * 100, 4),
                                 pred_pct=round(x * 100, 4),
                                 ratio=round(x / m, 3)))
    allr = [v for vs in per.values() for v in vs]
    xs = [NCF[c] for c in NCF if per.get(c)]
    ys = [math.log10(gm(per[c])) for c in NCF if per.get(c)]
    slope = float(np.polyfit(xs, ys, 1)[0])
    return per, gm(allr), slope, rows


def main():
    print(__doc__.split("Run:")[0].rstrip())
    kmw = chen_kmw()

    print("\n" + "=" * 86)
    print("THE SUBSTITUTED PARAMETERS, AGAINST WHAT ARMITAGE PREDICTS")
    print("=" * 86)
    V = P.volumes(P.ASSAYS["Rosenmai"])
    print(f"\n  {'cmpd':8s}{'nCF':>4s} | {'D_alb spLFER':>13s}{'Fischer':>10s}"
          f"{'ratio':>8s} | {'D_mem spLFER':>13s}{'Chen':>9s}{'ratio':>8s}")
    for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
        if nm not in NCF:
            continue
        lk = ko if ko is not None else Q.KOW_OPERA(nCF)
        D = P.armitage_D(lk, V["comp"])
        fa, fm = FISCHER_BSA[nCF] + VOL_SHIFT, kmw[nm]
        print(f"  {nm:8s}{nCF:>4d} | {math.log10(D['Dsa']):>13.2f}{fa:>10.2f}"
              f"{10**(fa-math.log10(D['Dsa'])):>8.2f} |"
              f" {math.log10(D['Dmw']):>13.2f}{fm:>9.2f}"
              f"{10**(fm-math.log10(D['Dmw'])):>8.2f}")
    print("""
  Both measured series sit BELOW the spLFER predictions, albumin by 3-9x and
  membrane lipid by 25-900x.  A linear balance parameterised on measurement
  therefore has a smaller albumin sink AND a smaller cell sink than the
  predicted one; the two shifts act on the numerator and the denominator of
  x_cell in opposite directions, so the net effect is not obvious in advance.""")

    print("\n" + "=" * 86)
    print("PERFORMANCE, nCF 5-8 (PFBA excluded: no Fischer BSA value)")
    print("=" * 86)
    print(f"\n  {'variant':46s}{'bias':>9s}{'slope/CF2':>12s}   per-compound")
    out, allrows = {}, []
    for variant, label in (
            ("published", "Armitage as published (spLFER albumin + membrane)"),
            ("measured_alb", "+ measured albumin only (Fischer 2024)"),
            ("measured_mem", "+ measured membrane only (Chen 2025)"),
            ("measured_both", "+ both measured  = Escher-school design")):
        per, bias, slope, rows = run(variant, kmw)
        out[variant] = (bias, slope, per)
        allrows += rows
        cells = "  ".join(f"{c} {gm(per[c]):.1f}x" for c in NCF)
        print(f"  {label:46s}{bias:>8.2f}x{slope:>+12.3f}   {cells}")

    # Qin, same compounds, for reference
    perq = collections.defaultdict(list)
    for study, a in P.ASSAYS.items():
        V2 = P.volumes(a)
        for nm in NCF:
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V2, a)
                if m is None:
                    continue
                perq[nm].append(P.qin_run(nm, d, V2)["x_cell"] / m)
    xs = [NCF[c] for c in NCF if perq.get(c)]
    ys = [math.log10(gm(perq[c])) for c in NCF if perq.get(c)]
    sq = float(np.polyfit(xs, ys, 1)[0])
    allq = [v for vs in perq.values() for v in vs]
    print(f"  {'Qin as published (saturable albumin)':46s}{gm(allq):>8.2f}x"
          f"{sq:>+12.3f}   " + "  ".join(f"{c} {gm(perq[c]):.1f}x" for c in NCF))

    b0, s0, _ = out["published"]
    bb, sb, _ = out["measured_both"]
    print(f"""
  WHAT THIS SETTLES.

  Replacing both spLFERs with measured series moves the bias from {b0:.2f}x to
  {bb:.2f}x and the chain-length slope from {s0:+.3f} to {sb:+.3f} decades per CF2.
  The level improves; the chain-length drift does not go away.

  That is the informative result for the model landscape.  The design decision
  that separates the Escher-school balances from IV-MBM EQP v2.0 -- measured
  biopartitioning in place of single-descriptor regressions -- addresses the
  LEVEL and leaves the chain-length behaviour substantially intact, because the
  residual slope is set by the RATIO of the cell to the albumin term and both
  measured series shift in the same direction.  The design decision that
  separates Qin from both -- a saturable albumin term -- is the one that
  removes the chain-length drift ({sq:+.3f}) while leaving the level worst of the
  four.

  So the two failure modes documented in the memo map onto two independent
  design axes, and no published model currently occupies the corner that would
  fix both: measured, saturable biopartitioning with a protein-realistic cell.""")

    dst = os.path.join(DB, "mbm_measured_parameter_variant.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
        w.writeheader()
        w.writerows(allrows)
    print(f"\n  wrote {dst} ({len(allrows)} rows)")


if __name__ == "__main__":
    main()
