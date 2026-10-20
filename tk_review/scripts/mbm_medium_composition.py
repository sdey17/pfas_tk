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
      2.2 "HepG2 cells were cultured in Dulbecco's modified Eagle medium ...
          and 10% fetal bovine serum"
      2.5 "cells were seeded in 12-well plates at a concentration of 18 x 10^4
          cells per well ... cells were exposed to PFASs in concentrations of
          10, 30 and 100 uM" -- no serum reduction stated, so 10%.  As coded.
      2.1 "PFAS stock solutions for in vitro investigations were prepared in
          ETHANOL in concentrations ranging from 4 to 350 mM"
      2.4 "Ethanol was used as control and the ethanol concentration was
          constant across all treatments and controls at 3% ETHANOL"
      -> 10% FBS as coded, but a 3% ethanol co-solvent: 6x Fragki's and 30x
         Amstutz's.  Stated for the PPARa assay; section 2.5 does not repeat
         it, though it doses the same stocks at the same concentrations.

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
    "Rosenmai": dict(fbs_growth=0.10, fbs_exposure=0.10,
                     solvent="ethanol", f_solvent=0.03,
                     v_dosed_uL=3000.0, v_after_t0_uL=3000.0,
                     verified=True,
                     cite="Rosenmai 2018 J Appl Toxicol 38:219 S2.1, 2.2, 2.4, 2.5"),
    "Amstutz":  dict(fbs_growth=0.10, fbs_exposure=0.10,
                     solvent="DMSO", f_solvent=0.001,
                     v_dosed_uL=250.0, v_after_t0_uL=200.0,
                     verified=True,
                     cite="Amstutz 2022 Toxicology 480:153312 S2.2, 2.3, 2.5"),
    "Fragki":   dict(fbs_growth=0.10, fbs_exposure=0.02,
                     solvent="DMSO", f_solvent=0.005,
                     v_dosed_uL=550.0, v_after_t0_uL=500.0,
                     verified=True,
                     cite="Fragki 2022 Food Chem Toxicol 172:113559 S2.2, 2.3"),
}

# Log-linear cosolvency: log K shifts by -sigma * f_cosolvent.  sigma for DMSO
# with hydrophobic organics is ~1.5-4 per unit volume fraction; 4 is the
# aggressive end, chosen so the estimate cannot be accused of being generous
# to the conclusion.
SIGMA = {"DMSO": 4.0, "ethanol": 4.5}


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
    print(f"\n  {'study':11s}{'growth':>8s}{'EXPOSURE':>10s}{'co-solvent':>14s}"
          f"{'dosed':>9s}{'after t0':>10s}  source")
    for s, d in SOURCED.items():
        sv = f"{d['f_solvent']:.1%} {d['solvent']}"
        mark = "" if d["verified"] else "   [NOT VERIFIED AT SOURCE]"
        print(f"  {s:11s}{d['fbs_growth']:>8.0%}{d['fbs_exposure']:>10.0%}"
              f"{sv:>14s}{d['v_dosed_uL']:>8.0f}uL{d['v_after_t0_uL']:>9.0f}uL"
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
    print("4.  THE CO-SOLVENTS: 3% ethanol vs 0.5% DMSO vs 0.1% DMSO")
    print("=" * 80)
    print("""
  Log-linear cosolvency: log K falls by sigma x f_cosolvent, with sigma ~ 1.5-4
  per unit volume fraction for a hydrophobic solute in aqueous DMSO and a
  little higher in ethanol.  Taking the aggressive end of each:
""")
    print(f"  {'study':11s}{'co-solvent':>15s}{'sigma':>7s}{'log shift':>11s}"
          f"{'as a factor':>13s}")
    for st, d in SOURCED.items():
        sig = SIGMA[d["solvent"]]
        sh = sig * d["f_solvent"]
        print(f"  {st:11s}{d['f_solvent']:>9.1%} {d['solvent']:<5s}{sig:>7.1f}"
              f"{sh:>11.3f}{10**sh:>12.2f}x")
    print(f"""
  Rosenmai's 3% ethanol is {SIGMA['ethanol']*0.03:.2f} log units, {SIGMA['ethanol']*0.03/(SIGMA['DMSO']*0.005):.0f}x the DMSO shift in Fragki
  and {SIGMA['ethanol']*0.03/(SIGMA['DMSO']*0.001):.0f}x the one in Amstutz -- so the co-solvent IS a between-study
  variable after all, just not the one the 2% vs 10% looked like.  The serum
  difference is still {math.log10(0.10/0.02):.2f} log units, {math.log10(0.10/0.02)/(SIGMA['ethanol']*0.03):.1f}x larger than even the ethanol.

  BUT -- and this is the part that matters -- a co-solvent lowers EVERY
  water-to-phase partition coefficient by the same factor, and the quantity
  being predicted is a RATIO of sorptive terms:

      x_cell = D_cell V_cell / (V_w + D_alb V_alb + D_lip V_lip + D_cell V_cell)

  Divide every D by k and the k cancels in numerator and denominator -- except
  against V_w, the one term with no D on it.  So the co-solvent only bites
  where water carries a real share of the denominator.  Run it and see:
""")
    print(f"  {'all D divided by':>18s}{'Armitage bias':>16s}{'PFBA':>9s}"
          f"{'PFHxA':>9s}{'PFOA':>9s}{'PFNA':>9s}   water share at PFBA")
    V = P.volumes(P.ASSAYS["Rosenmai"])
    for k in (1.0, 1.35, 1.5):
        per, allr = {}, []
        for (nm, nC, nCF, mw, ko, kc, ks) in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            D = P.armitage_D(lk, V["comp"])
            for d in P.ASSAYS["Rosenmai"]["doses"]:
                m = P.measured_x_cell("Rosenmai", nm, d, V, P.ASSAYS["Rosenmai"])
                if m is None:
                    continue
                ntot = d * V["Vbm"]
                den = (V["Vw"] + D["Dsa"] / k * V["Valb"]
                       + D["Dow"] / k * V["Vslip"] + D["Dcw"] / k * V["Vcells"]
                       + D["Dpl"] / k * V["sarea"] * 1000.0)
                x = D["Dcw"] / k * (ntot / den) * V["Vcells"] / ntot
                per.setdefault(nm, []).append(x / m)
                allr.append(x / m)
        lkb = next(c[4] if c[4] is not None else Q.KOW_OPERA(c[2])
                   for c in Q.COMPOUNDS if c[0] == "PFBA")
        Db = P.armitage_D(lkb, V["comp"])
        denb = (V["Vw"] + Db["Dsa"] / k * V["Valb"] + Db["Dow"] / k * V["Vslip"]
                + Db["Dcw"] / k * V["Vcells"] + Db["Dpl"] / k * V["sarea"] * 1000.0)
        print(f"  {k:>18.2f}{gm(allr):>15.2f}x" + "".join(
            f"{gm(per[c]):>8.2f}x" for c in ("PFBA", "PFHxA", "PFOA", "PFNA"))
            + f"{V['Vw']/denb:>21.0%}")
    print("""
  PFBA moves 1.42x -> 1.02x and PFNA does not move at all (9.33x -> 9.30x),
  because water is 80% of PFBA's denominator and 0.5% of PFNA's.  The overall
  bias barely stirs, 1.97x -> 1.82x.

  So Rosenmai's 3% ethanol is a real protocol difference, it is the largest
  co-solvent of the three, and it is worth almost nothing to this comparison
  -- for the specific reason that the fraction-in-cells axis is scale-free in
  the partition coefficients.  That is also the general lesson: an error that
  multiplies EVERY sorptive term cancels here, and only an error in the
  RELATIVE size of two sinks (such as how much albumin is present) survives.

  One DMSO point is not about partitioning at all.  HepaRG are differentiated
  in 1.7% DMSO for twelve days and held there.  That is what makes them
  hepatocyte-like, and it also loads them with triglyceride -- the companion
  paper measures exactly that as an endpoint.  Both models run Fragki on the
  AREc2 composition (a HepG2-derived reporter line, membrane lipid 0.00414),
  because that is the only liver-line composition Armitage's SI publishes.  A
  differentiated, DMSO-matured HepaRG is a different cell.  That error is in
  the cell term, which section 14 of the main memo already names as where the
  residual lives, and it does NOT cancel, because it changes one sink and not
  the others.""")

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
