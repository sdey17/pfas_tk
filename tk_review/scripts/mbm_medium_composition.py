#!/usr/bin/env python3
"""The exposure medium is not the growth medium: Fragki dosed in 2% FBS.

WHAT THIS CHECKS.  Every term in both mass balances that lives outside the
cell is proportional to how much serum the medium carries:

    V_albumin      = V_medium x f_FBS x 0.733 mL/g x 42 g/L
    V_serum_lipid  = V_medium x f_FBS x 1.03 g/L

so f_FBS multiplies the two biggest sinks in the denominator.  The three
assays were all coded at f_FBS = 0.10, taken from each paper's GROWTH medium.
Two of the three state a different number for the medium the cells were
actually DOSED in.

    Amstutz et al. 2022, Toxicology 480:153312, doi 10.1016/j.tox.2022.153312
      2.2 "HepG2 cells were cultured ... in DMEM supplemented with 10% (v/v)
          FBS, 1% (v/v) L-glutamine, and 100 IU/mL of pen-strep."
      2.3 "exposed to the PFAS ... added from stocks in DMSO TO THE CULTURE
          MEDIUM for a final volume of 200 uL (final DMSO concentration 0.1%)"
      2.5 "Cells were exposed to 250 uL exposure medium, and at t = 0, a sample
          of 50 uL was taken for mass balance calculations ... exposed for 24 h
          to the remaining 200 uL"
      -> exposure in the full 10% FBS culture medium.  As coded.

    Fragki et al. 2022, Food Chem Toxicol 172:113559, doi 10.1016/j.fct.2022.113559
      2.2 "cultured in growth medium consisting of William's Medium E ...
          supplemented with 10% fetal bovine serum"
          "cells were cultured for 12 days in growth medium supplemented with
          1.7% DMSO (differentiation medium) for final differentiation"
          "Prior to biokinetic studies, differentiated HepaRG cells were
          incubated for 24 h in ASSAY MEDIUM (GROWTH MEDIUM CONTAINING 2% FBS)
          supplemented with 0.5% DMSO"
      2.3 "Test chemicals were diluted from 200-fold concentrated stock
          solutions IN ASSAY MEDIUM, providing exposure medium with a final
          DMSO concentration of 0.5%.  Cells were exposed to 550 uL exposure
          medium and at t = 0, a 50 uL medium sample was taken"
      -> exposure in 2% FBS, not 10%.  NOT as coded: 5x too much albumin.

    Rosenmai et al. 2018, J Appl Toxicol 38:219-226, doi 10.1002/jat.3515
      Full text paywalled here; 10% is the standard HepG2 condition and is
      what Amstutz followed ("based on the study of Rosenmai et al. (2018),
      with small modifications"), but it is NOT verified at source.  Flagged.

The two cell-associated datasets that the QIVIVE literature treats as
interchangeable were therefore measured five-fold apart in the one medium
property that the mass balance is most sensitive to.

Run:  python3 tk_review/scripts/mbm_medium_composition.py
"""
import copy
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P      # noqa: E402
import mbm_qin_vs_armitage as Q          # noqa: E402

DB = os.path.join(os.path.dirname(HERE), "db")

# Sourced exposure conditions.  fbs_exposure is the number this script argues
# for; fbs_growth is what the paper's culture section says and what was coded.
SOURCED = {
    "Rosenmai": dict(fbs_growth=0.10, fbs_exposure=0.10, dmso=None,
                     v_dosed_uL=3000.0, v_after_t0_uL=3000.0,
                     verified=False,
                     cite="Rosenmai 2018 J Appl Toxicol 38:219 (paywalled)"),
    "Amstutz":  dict(fbs_growth=0.10, fbs_exposure=0.10, dmso=0.001,
                     v_dosed_uL=250.0, v_after_t0_uL=200.0,
                     verified=True,
                     cite="Amstutz 2022 Toxicology 480:153312 S2.2, 2.3, 2.5"),
    "Fragki":   dict(fbs_growth=0.10, fbs_exposure=0.02, dmso=0.005,
                     v_dosed_uL=550.0, v_after_t0_uL=500.0,
                     verified=True,
                     cite="Fragki 2022 Food Chem Toxicol 172:113559 S2.2, 2.3"),
}

# Log-linear cosolvency: log K shifts by -sigma * f_cosolvent.  sigma for DMSO
# with hydrophobic organics is ~1.5-4 per unit volume fraction; 4 is the
# aggressive end, chosen so the estimate cannot be accused of being generous
# to the conclusion.
SIGMA_DMSO = 4.0


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def gmae(xs):
    return 10 ** (sum(abs(math.log10(x)) for x in xs) / len(xs))


def run(assays):
    """Both models over every measured point, grouped by study."""
    out, rows = {}, []
    for study, a in assays.items():
        V = P.volumes(a)
        for (nm, nC, nCF, mw, ko, kc, ksrc) in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                q = P.qin_run(nm, d, V)
                ar = P.armitage_run(nm, lk, d, V)
                out.setdefault(study, []).append((q["x_cell"] / m,
                                                  ar["x_cell"] / m))
                rows.append(dict(study=study, compound=nm, dose_uM=d,
                                 fbs=a["FBS"],
                                 meas_x_cell_pct=round(m * 100, 4),
                                 qin_x_cell_pct=round(q["x_cell"] * 100, 4),
                                 arm_x_cell_pct=round(ar["x_cell"] * 100, 4),
                                 qin_over_meas=round(q["x_cell"] / m, 3),
                                 arm_over_meas=round(ar["x_cell"] / m, 3)))
    return out, rows


def assays_with(fbs_map):
    A = copy.deepcopy(P.ASSAYS)
    for k, v in fbs_map.items():
        A[k]["FBS"] = v
    return A


def table(out, label):
    print(f"\n  {label}")
    print(f"  {'study':11s}{'n':>4s}{'Qin bias':>11s}{'Qin |err|':>11s}"
          f"{'Arm bias':>11s}{'Arm |err|':>11s}")
    for s in ("Rosenmai", "Amstutz", "Fragki"):
        q = [x[0] for x in out[s]]
        a = [x[1] for x in out[s]]
        print(f"  {s:11s}{len(q):>4d}{gm(q):>10.2f}x{gmae(q):>10.2f}x"
              f"{gm(a):>10.2f}x{gmae(a):>10.2f}x")
    q = [x[0] for v in out.values() for x in v]
    a = [x[1] for v in out.values() for x in v]
    print(f"  {'ALL':11s}{len(q):>4d}{gm(q):>10.2f}x{gmae(q):>10.2f}x"
          f"{gm(a):>10.2f}x{gmae(a):>10.2f}x")


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print("=" * 80)
    print("1.  WHAT EACH PAPER ACTUALLY SAYS")
    print("=" * 80)
    print(f"\n  {'study':11s}{'growth':>8s}{'EXPOSURE':>10s}{'DMSO':>8s}"
          f"{'dosed':>9s}{'after t0':>10s}  source")
    for s, d in SOURCED.items():
        dm = "-" if d["dmso"] is None else f"{d['dmso']:.1%}"
        mark = "" if d["verified"] else "   [NOT VERIFIED AT SOURCE]"
        print(f"  {s:11s}{d['fbs_growth']:>8.0%}{d['fbs_exposure']:>10.0%}"
              f"{dm:>8s}{d['v_dosed_uL']:>8.0f}uL{d['v_after_t0_uL']:>9.0f}uL"
              f"  {d['cite']}{mark}")

    # ------------------------------------------------------------------ 2
    print("\n" + "=" * 80)
    print("2.  WHAT 10% -> 2% DOES TO THE DENOMINATOR, FOR FRAGKI")
    print("=" * 80)
    print("""
  C_free = n_total / S, with S = Vw + D_alb Valb + D_lip Vslip + cell + plastic.
  Shrinking the serum shrinks S, which RAISES the free concentration and so
  raises every predicted uptake.  The share table says by how much.
""")
    for fbs in (0.10, 0.02):
        a = copy.deepcopy(P.ASSAYS["Fragki"])
        a["FBS"] = fbs
        V = P.volumes(a)
        print(f"  Fragki at {fbs:.0%} FBS:  Valb {V['Valb']*1e6:7.4f} uL   "
              f"Vslip {V['Vslip']*1e6:7.4f} uL   Vw {V['Vw']*1e6:7.2f} uL")
        print(f"  {'':14s}{'cmpd':8s}{'albumin':>10s}{'ser.lipid':>11s}"
              f"{'water':>9s}{'cells':>9s}{'plastic':>9s}{'S (L)':>12s}")
        for nm in ("PFBA", "PFOA", "PFNA"):
            lk = next(c[4] if c[4] is not None else Q.KOW_OPERA(c[2])
                      for c in Q.COMPOUNDS if c[0] == nm)
            D = P.armitage_D(lk, V["comp"])
            parts = {"albumin": D["Dsa"] * V["Valb"],
                     "ser.lipid": D["Dow"] * V["Vslip"],
                     "water": V["Vw"],
                     "cells": D["Dcw"] * V["Vcells"],
                     "plastic": D["Dpl"] * V["sarea"] * 1000.0}
            S = sum(parts.values())
            print(f"  {'':14s}{nm:8s}" + "".join(
                f"{parts[k]/S:>{w}.1%}" for k, w in
                (("albumin", 10), ("ser.lipid", 11), ("water", 9),
                 ("cells", 9), ("plastic", 9))) + f"{S:>12.3e}")
        print()

    # ------------------------------------------------------------------ 3
    print("=" * 80)
    print("3.  THE MODELS, BEFORE AND AFTER")
    print("=" * 80)
    before, _ = run(assays_with({}))
    after, rows = run(assays_with({"Fragki": 0.02}))
    table(before, "AS CODED -- all three at 10% FBS")
    table(after, "CORRECTED -- Fragki at 2% FBS, the medium it was dosed in")
    qb = gm([x[0] for x in before["Fragki"]])
    qa = gm([x[0] for x in after["Fragki"]])
    ab = gm([x[1] for x in before["Fragki"]])
    aa = gm([x[1] for x in after["Fragki"]])
    print(f"""
  Fragki's Qin bias moves {qb:.2f}x -> {qa:.2f}x ({qa/qb:.1f}x worse) and its Armitage bias
  {ab:.2f}x -> {aa:.2f}x ({aa/ab:.1f}x worse).  The dataset that looked like the models'
  best case was the one being handed five times more albumin than it had.

  Armitage moves by {aa/ab:.1f}x where the pure denominator argument predicts about
  4x, because Armitage's albumin term is linear.  Qin moves by {qa/qb:.1f}x, MORE than
  the denominator shrinks, because Qin's albumin is a saturating Langmuir: a
  higher free concentration drives D_BSA/w down, which shrinks the denominator
  again.  A model with a saturating sink amplifies an error in the size of
  that sink.""")

    # ------------------------------------------------------------------ 4
    print("\n" + "=" * 80)
    print("4.  IS THE 0.5% vs 0.1% DMSO WORTH ANYTHING?")
    print("=" * 80)
    shift = SIGMA_DMSO * (SOURCED["Fragki"]["dmso"] - SOURCED["Amstutz"]["dmso"])
    print(f"""
  Log-linear cosolvency: log K falls by sigma x f_cosolvent, with sigma ~ 1.5-4
  per unit volume fraction for a hydrophobic solute in aqueous DMSO.  Taking
  the aggressive sigma = {SIGMA_DMSO:.0f}:

    Fragki  0.5% DMSO  ->  {SIGMA_DMSO*SOURCED['Fragki']['dmso']:.3f} log units = {10**(SIGMA_DMSO*SOURCED['Fragki']['dmso']):.3f}x on any partition coefficient
    Amstutz 0.1% DMSO  ->  {SIGMA_DMSO*SOURCED['Amstutz']['dmso']:.3f} log units = {10**(SIGMA_DMSO*SOURCED['Amstutz']['dmso']):.3f}x
    difference between them: {shift:.3f} log units = {10**shift:.3f}x

  The serum difference is {math.log10(0.10/0.02):.2f} log units = 5.00x.  The solvent difference is
  {shift:.3f} log units -- {math.log10(0.10/0.02)/shift:.0f} times smaller, and inside the noise on any of these
  measurements.  So: the 2% vs 10% that matters is the SERUM, not the DMSO.

  One DMSO caveat that is NOT about partitioning.  HepaRG are differentiated in
  1.7% DMSO for twelve days and held there.  That is what makes them
  hepatocyte-like, and it also loads them with triglyceride -- the companion
  paper measures exactly that as an endpoint.  Both models run Fragki on the
  AREc2 composition (a HepG2-derived reporter line, membrane lipid 0.00414),
  because that is the only composition Armitage's SI publishes for a liver
  line.  A differentiated, DMSO-matured HepaRG is a different cell.  That is a
  separate error from the medium, and it is in the cell term, which section 14
  of the main memo already names as where the residual lives.""")

    # ------------------------------------------------------------------ 5
    print("\n" + "=" * 80)
    print("5.  THE VOLUME SUBTLETY, FOR COMPLETENESS")
    print("=" * 80)
    print("""
  Both papers take a 50 uL medium sample at t = 0 for mass balance, so the
  cells see less medium than was dosed: Fragki 500 of 550 uL, Amstutz 200 of
  250.  Both studies normalise the uptake fraction to the mass NOMINALLY
  APPLIED, so the coded volumes (550, 250) match the measurement's
  denominator.  Using the post-sample volume instead moves things like this:
""")
    A = assays_with({"Fragki": 0.02})
    for study in ("Amstutz", "Fragki"):
        a2 = copy.deepcopy(A)
        a2[study]["v_uL"] = SOURCED[study]["v_after_t0_uL"]
        o2, _ = run(a2)
        q1 = gm([x[0] for x in after[study]])
        q2 = gm([x[0] for x in o2[study]])
        a1 = gm([x[1] for x in after[study]])
        a2b = gm([x[1] for x in o2[study]])
        print(f"  {study:10s} {SOURCED[study]['v_dosed_uL']:.0f} -> "
              f"{SOURCED[study]['v_after_t0_uL']:.0f} uL:  Qin {q1:.2f}x -> "
              f"{q2:.2f}x,  Armitage {a1:.2f}x -> {a2b:.2f}x")
    print("""
  Under 1.25x either way, and in the direction that makes the models worse.
  Left as dosed, because that is the denominator the published uptake
  fractions use; recorded, not adopted.""")

    # ------------------------------------------------------------------ 6
    print("\n" + "=" * 80)
    print("6.  DOES THE MEASURED DATA SHOW THE SERUM EFFECT?")
    print("=" * 80)
    print("""
  If 2% FBS really leaves five times less albumin to compete with the cells,
  Fragki's MEASURED uptake fraction should sit above the two 10% FBS studies at
  the same nominal dose.  It does -- but the cell line differs too, so this
  corroborates the direction without isolating the cause.
""")
    print(f"  {'cmpd':7s}{'dose':>7s}{'Rosenmai 10%':>14s}{'Amstutz 10%':>13s}"
          f"{'Fragki 2%':>11s}{'Fragki / 10% studies':>22s}")
    for nm in ("PFOA", "PFNA"):
        for d in (6.25, 10.0, 100.0):
            vals = {}
            for study in ("Rosenmai", "Amstutz", "Fragki"):
                a = A[study]
                if d not in a["doses"]:
                    continue
                vals[study] = P.measured_x_cell(study, nm, d, P.volumes(a), a)
            if "Fragki" not in vals or vals["Fragki"] is None:
                continue
            tens = [v for k, v in vals.items() if k != "Fragki" and v]
            rat = (f"{vals['Fragki']/gm(tens):.1f}x" if tens else "-")
            cells = []
            for study in ("Rosenmai", "Amstutz", "Fragki"):
                v = vals.get(study)
                cells.append(f"{v*100:.3f}" if v else "-")
            print(f"  {nm:7s}{d:>7g}{cells[0]:>14s}{cells[1]:>13s}"
                  f"{cells[2]:>11s}{rat:>22s}")
    print("""
  PFNA at 100 uM is the one point Fragki themselves set aside: "given the
  observed low protein content this concentration level was excluded from the
  QIVIVE analysis, as this may indicate a possible cytotoxic effect".  It is
  still in the model comparison here, and it is the largest measured value in
  the set.""")

    dst = os.path.join(DB, "mbm_medium_composition.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {dst} ({len(rows)} rows, Fragki at the corrected 2% FBS)")


if __name__ == "__main__":
    main()
