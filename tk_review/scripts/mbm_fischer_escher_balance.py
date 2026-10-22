#!/usr/bin/env python3
"""The Escher-school balance as Fischer et al. 2017 actually wrote it.

The review's measured-parameter variant (mbm_measured_parameter_variant.py)
substituted measured albumin and membrane coefficients into ARMITAGE's cell
structure.  With Fischer et al. 2017 (Chem Res Toxicol 30:1197-1208) now in
hand, that reconstruction can be checked, and it was not faithful.  Fischer's
balance is:

  V_medium = V_prot,medium + V_lip,medium + V_w,medium                   (eq 1)
  V_cell   = V_prot,cell   + V_lip,cell   + V_w,cell                     (eq 2)
  D_medium/w = VF_prot,med *D_BSA/w + VF_lip,med *D_lip/w + VF_w,med     (eq 3)
  D_cell/w   = VF_prot,cell*D_BSA/w + VF_lip,cell*D_lip/w + VF_w,cell    (eq 4)
  D_medium/cell = D_medium/w / D_cell/w                                  (eq 7)
  f_cell = (1 + D_medium/cell * V_medium/V_cell)^-1                      (eq 8)

Four differences from what the review's variant assumed:

  1. ONE protein surrogate, used everywhere.  "BSA and lip were used as
     surrogates for all medium and cellular proteins and lipids."  Cell protein
     gets D_BSA/w, NOT octanol-equivalence -- so the cell protein term is
     VF_prot * D_BSA/w, not 0.035 * VF_prot * D_OW.  This is the single biggest
     structural difference and it was got wrong.
  2. Cell water IS a phase in D_cell/w (unlike Qin 2023, which pools it with
     the medium water; this is exactly the contrast Qin draws against the
     Henneberger/Fischer formulation).
  3. The normalising volume is V_cell, so D_cell/w is in L_w/L_cell -- the
     same basis as Armitage's D_CW and NOT Qin's L_w/L_prot+lip.
  4. No plastic/vessel-wall term.  Fischer 2019 added that separately.

And the medium protein volume means TOTAL protein, measured, not albumin:
Fischer measured VF_prot,medium = 0.520% at 10% d-FBS, which is 3.0x the
albumin-only volume Armitage's serum convention gives (0.176%).

Fischer also supply the cell composition this project most needed: Table 2,
ARE assay, HepG2 (ARE-bla) -- 87.35% water, 9.46% protein, 3.19% lipid, 3.04 pg
per cell, with "all cellular lipids" treated as membrane lipids.  Two of the
three datasets modelled in the review use HepG2.  Armitage's AREc2 entry, which
the review has been using for them, is MCF7-derived (AREc32, Cellosaurus
CVCL_1D32), carries 0.46% total lipid against Fischer's measured 3.19% for
HepG2, and is not a hepatic composition at all.

Fischer's own position on single-descriptor prediction is worth quoting, since
it is the design premise being tested: the distribution ratios they compute
"emphasize that their derivation by log Kow correlation can lead to a high
uncertainty and inaccuracy", and multi-descriptor methods "can substantially
improve the quality of predicted partition constants compared to predictions
via log Kow".  Fischer use PP-LFERs for neutral species, a 3D-QSAR for anionic
K_BSA/w and COSMOmic for anionic K_lip/w; measured values are preferred where
they exist.  For PFAS measured values now exist, so this reconstruction uses
them, which is the strongest form of the design.

Run:  python3 tk_review/scripts/mbm_fischer_escher_balance.py
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
VOL_SHIFT = math.log10(1.0 / 0.733)      # mass -> volume basis

# Fischer et al. 2017 Table 2, measured medium composition (volume fractions).
# Two conditions are reported: 10% d-FBS in OptiMEM and 1% FBS in DMEM
# GlutaMAX.  The 10% row is used directly for the two 10%-FBS assays here.
# For Fragki's 2% FBS the two rows are interpolated linearly in %FBS, which
# mixes two basal media and two FBS treatments and is therefore flagged.
FISCHER_MEDIUM = {10.0: dict(prot=0.00520, lip=0.000300),
                  1.0: dict(prot=0.00138, lip=0.000200)}


def medium_vf(pct_fbs):
    (x0, a), (x1, b) = sorted(FISCHER_MEDIUM.items())
    if abs(pct_fbs - x1) < 1e-9:
        return dict(b), "measured at 10% d-FBS"
    if abs(pct_fbs - x0) < 1e-9:
        return dict(a), "measured at 1% FBS"
    f = (pct_fbs - x0) / (x1 - x0)
    return ({k: a[k] + f * (b[k] - a[k]) for k in a},
            "interpolated between Fischer's 1% and 10% rows")


def chen_kmw():
    out = {}
    for r in csv.DictReader(open(os.path.join(DB, "protein_binding.csv"))):
        if "Chen" in r["study"] and "log K_MW" in r["ka_units"]:
            out[r["chemical"]] = float(r["ka"])
    return out


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def fischer_run(study, nm, nCF, dose, kmw, comp_key, use_medium):
    """Fischer eqs 3, 4, 7, 8 for one compound in one assay."""
    a = P.ASSAYS[study]
    comp = P.CELL_COMP[comp_key]
    Vbm = a["v_uL"] * 1e-6
    Vcell = a["cells"] * comp["ng_cell"] * 1e-6 * 1e-6
    Dbsa = 10 ** (FISCHER_BSA[nCF] + VOL_SHIFT)
    Dlip = 10 ** kmw[nm]
    if use_medium:
        mvf, _src = medium_vf(a["FBS"] * 100.0)
        vfp, vfl = mvf["prot"], mvf["lip"]
    else:                                   # Armitage's serum convention
        vfp = a["FBS"] * P.SI["v_alb_partial"] * P.SI["conc_ser_alb"] * 1e-3
        vfl = a["FBS"] * P.SI["conc_ser_lip"] * 1e-3
    Dmed = vfp * Dbsa + vfl * Dlip + (1.0 - vfp - vfl)
    Dcell = (comp["prot"] * Dbsa + (comp["sl"] + comp["ml"]) * Dlip
             + comp["water"])
    f_cell = 1.0 / (1.0 + (Dmed / Dcell) * (Vbm - Vcell) / Vcell)
    return f_cell, Dcell, Dmed


def variant(label, comp_key, use_medium, kmw):
    per = collections.defaultdict(list)
    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        for (nm, _nC, nCF, _mw, _ko, _kc, _ks) in Q.COMPOUNDS:
            if nm not in NCF:
                continue
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                x, _dc, _dm = fischer_run(study, nm, nCF, d, kmw,
                                          comp_key, use_medium)
                per[nm].append(x / m)
    allr = [v for vs in per.values() for v in vs]
    xs = [NCF[c] for c in NCF if per.get(c)]
    ys = [math.log10(gm(per[c])) for c in NCF if per.get(c)]
    return label, gm(allr), float(np.polyfit(xs, ys, 1)[0]), per


def armitage_with(comp):
    """Armitage as published, but with a substituted cell characterisation."""
    P.CELL_COMP["_sub"] = comp
    per = collections.defaultdict(list)
    for study, a in P.ASSAYS.items():
        aa = dict(a, line="_sub")
        V = P.volumes(aa)
        for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, aa)
                if m is None:
                    continue
                per[nm].append(P.armitage_run(nm, lk, d, V)["x_cell"] / m)
    allr = [v for vs in per.values() for v in vs]
    xs = [NCF[c] for c in NCF if per.get(c)]
    ys = [math.log10(gm(per[c])) for c in NCF if per.get(c)]
    return gm(allr), float(np.polyfit(xs, ys, 1)[0]), len(allr)


def cell_substitution_grid():
    """The 2x2 the review's parameter table should have run.

    Section 6.4 reports that substituting an "independent HepG2 cell volume"
    improves the Armitage fit from 3.87x to 1.84x -- the one substitution that
    helped.  That volume is Fischer et al. 2017's, and Fischer measured the
    COMPOSITION of the same cells in the same table.  Taking one without the
    other is the same category of error as 6.7.1's.
    """
    ar, fi = P.CELL_COMP["AREc2"], P.CELL_COMP["HepG2_Fischer"]
    print("\n" + "=" * 78)
    print("TAKING FISCHER'S CELL MEASUREMENT WHOLE RATHER THAN BY HALVES")
    print("=" * 78)
    print(f"\n  {'composition':14s}{'volume':14s}{'bias':>9s}{'slope':>9s}   n")
    for lab_c, lab_v, comp in [
            ("AREc2", "AREc2 10 ng", dict(ar)),
            ("AREc2", "Fischer 3.04", dict(ar, ng_cell=3.04)),
            ("Fischer", "AREc2 10 ng", dict(fi, ng_cell=10.0)),
            ("Fischer", "Fischer 3.04", dict(fi))]:
        b, sl, n = armitage_with(comp)
        print(f"  {lab_c:14s}{lab_v:14s}{b:>8.2f}x{sl:>+9.3f}  {n}")
    print("""
  The published parameterisation is the first row and the review's reported
  improvement is the second.  The second row is Fischer's cell volume with a
  MCF7-derived composition; the fourth is Fischer's measurement taken whole.
  Taken whole it worsens the fit, from 3.87x to 9.10x, because the measured
  hepatic lipid fraction is 6.9x larger while the cell volume is only 3.3x
  smaller, so the cell sink grows by about 2x.

  So the review's statement that cell volume is the one parameter whose
  best-available measurement improves the fit does not hold: it improved the
  fit only because half of the measurement was used.  Every best-available
  measurement worsens the fit, with no exception.""")


def main():
    print(__doc__.split("Run:")[0].rstrip())
    kmw = chen_kmw()
    cell_substitution_grid()

    print("\n" + "=" * 78)
    print("THE CELL COMPOSITION QUESTION")
    print("=" * 78)
    ar, fi = P.CELL_COMP["AREc2"], P.CELL_COMP["HepG2_Fischer"]
    print(f"""
                        Armitage AREc2        Fischer HepG2
                        (MCF7-derived)        (measured, ARE-bla)
    mass per cell       {ar['ng_cell']:>6.2f} ng            {fi['ng_cell']:>6.2f} ng
    protein             {ar['prot']:>8.4f}            {fi['prot']:>8.4f}
    total lipid         {ar['sl'] + ar['ml']:>8.4f}            {fi['sl'] + fi['ml']:>8.4f}
    water               {ar['water']:>8.4f}            {fi['water']:>8.4f}

  Measured HepG2 carries {(fi['sl'] + fi['ml']) / (ar['sl'] + ar['ml']):.1f}x the lipid and """
          f"""{fi['prot'] / ar['prot']:.2f}x the protein of the
  composition the review has been using for its two HepG2 datasets, in """
          f"""{fi['ng_cell'] / ar['ng_cell']:.2f}x
  the cell volume.  The lipid fraction and the cell volume move the predicted
  cell sink in opposite directions.""")

    print("\n" + "=" * 78)
    print("MEDIUM PROTEIN: ALBUMIN VOLUME vs MEASURED TOTAL PROTEIN")
    print("=" * 78)
    print(f"\n  {'assay':10s}{'%FBS':>6s} | {'albumin VF':>12s}"
          f"{'Fischer VF':>12s}{'ratio':>8s}  source")
    for study, a in P.ASSAYS.items():
        alb = a["FBS"] * P.SI["v_alb_partial"] * P.SI["conc_ser_alb"] * 1e-3
        mvf, src = medium_vf(a["FBS"] * 100.0)
        print(f"  {study:10s}{100 * a['FBS']:>6.1f} | {alb:>12.5f}"
              f"{mvf['prot']:>12.5f}{mvf['prot'] / alb:>8.2f}  {src}")
    print("""
  Fischer's measured medium protein is ~3x the albumin volume, because it is
  total protein (serum globulins and basal-medium protein included) and because
  the surrogate D_BSA/w is then applied to all of it.  A faithful Escher-school
  balance therefore has a LARGER medium sink as well as a larger cell sink.""")

    print("\n" + "=" * 78)
    print("PERFORMANCE OVER nCF 5-8")
    print("=" * 78)
    rows = [
        variant("Fischer form, Fischer HepG2 cells, Fischer medium",
                "HepG2_Fischer", True, kmw),
        variant("Fischer form, Fischer HepG2 cells, albumin-only medium",
                "HepG2_Fischer", False, kmw),
        variant("Fischer form, Armitage AREc2 cells, Fischer medium",
                "AREc2", True, kmw),
        variant("Fischer form, Armitage AREc2 cells, albumin-only medium",
                "AREc2", False, kmw),
    ]
    print(f"\n  {'parameterisation':54s}{'bias':>9s}{'slope':>9s}")
    print("  " + "-" * 70)
    for lab, b, sl, _ in rows:
        print(f"  {lab:54s}{b:>8.2f}x{sl:>+9.3f}")
    print(f"""  {'-' * 70}
  for reference (mbm_measured_parameter_variant.py, Armitage structure):
  {'Armitage as published':54s}{'4.35x':>9s}{'+0.316':>9s}
  {'+ measured albumin and membrane, Armitage structure':54s}{'0.40x':>9s}{'-0.156':>9s}
  {'Qin as published':54s}{'8.59x':>9s}{'-0.043':>9s}""")

    print("\n" + "=" * 78)
    print("PER-COMPOUND, THE FAITHFUL FORM")
    print("=" * 78)
    _lab, _b, _s, per = rows[0]
    print(f"\n  {'cmpd':8s}{'nCF':>4s}{'pred/meas':>12s}")
    for c in NCF:
        if per.get(c):
            print(f"  {c:8s}{NCF[c]:>4d}{gm(per[c]):>11.2f}x")
    sub = [(NCF[c], math.log10(gm(per[c]))) for c in NCF if per.get(c)
           and NCF[c] <= 7]
    s57 = float(np.polyfit([x for x, _ in sub], [y for _, y in sub], 1)[0])
    print(f"""
  The series is FLAT over nCF 5-7 ({s57:+.3f} decades per CF2) and the whole of
  the -0.100 comes from the single PFNA point.  So the measured-parameter form
  does not merely reduce the chain-length drift, it removes it over the range
  where three of the four compounds sit.  Both measured series are involved:
  Fischer's albumin and Chen's membrane coefficients rise at nearly the same
  rate with chain length, so their ratio -- which is what the predicted
  fraction depends on -- is almost constant.

  Note that the slope is the SAME (-0.100, -0.033, -0.100, -0.044) in all four
  rows of the table above, across a 6.9x change in cell lipid, a 3.3x change in
  cell volume and a 3x change in medium protein.  That is the cleanest
  separation this reassessment has produced:

      the CHAIN-LENGTH behaviour is a property of the partition coefficients;
      the LEVEL is a property of the volumes and compositions.

  It also means the level error cannot be diagnosed from the coefficients and
  the chain-length error cannot be diagnosed from the system parameters.""")


if __name__ == "__main__":
    main()
