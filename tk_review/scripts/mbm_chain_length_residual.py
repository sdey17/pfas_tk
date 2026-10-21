#!/usr/bin/env python3
"""What is left after the medium correction: a residual that rises with chain
length.  This localises it.

THE OBSERVATION.  On the corrected exposure media (section 15 of the memo),
Armitage's over-prediction climbs monotonically with perfluorocarbon chain
length -- +0.172 decades per CF2 over the pooled 29 points, r = +0.83.  A
level error would not do that.  Something in the model's chain-length
dependence is wrong, or something in the measurements is.

FIVE QUESTIONS, IN ORDER.  Each one can kill a candidate outright:

  1. Is the trend real, or an artefact of pooling three studies that measure
     different compound sets?  (Fragki contributes only PFOA and PFNA, and
     has the largest residuals, so pooling alone could manufacture a slope.)
  2. Is it a measurement artefact at the long-chain end -- cytotoxicity
     suppressing the measured uptake?
  3. Is the fitted slope itself an artefact of PFBA's leverage in a five-point
     regression?
  4. Which term in the balance carries it?  Substitute each candidate in turn
     and watch the SLOPE, not the level.
  5. What is the actual mechanism, and what would settle it?

Run:  python3 tk_review/scripts/mbm_chain_length_residual.py
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
NCF = {"PFBA": 3, "PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}


def gm(v):
    return 10 ** float(np.mean([math.log10(x) for x in v]))


def fit(per):
    """Slope and r of log10(residual) against nCF, over the compounds present."""
    cs = [c for c in NCF if per.get(c)]
    if len(cs) < 3:
        return None, None, {c: gm(per[c]) for c in cs}
    x = [NCF[c] for c in cs]
    y = [math.log10(gm(per[c])) for c in cs]
    return (float(np.polyfit(x, y, 1)[0]),
            float(np.corrcoef(x, y)[0, 1]),
            {c: gm(per[c]) for c in cs})


def row(label, per, width=44):
    sl, r, g = fit(per)
    cells = "".join(f"{g[c]:>8.2f}x" if c in g else f"{'-':>9s}" for c in NCF)
    tail = f"{sl:>+8.3f}{r:>+7.2f}" if sl is not None else f"{'(n<3)':>15s}"
    print(f"  {label:{width}s}{cells}{tail}")
    return sl


def residuals(study=None, dose=None, model="Armitage"):
    per = collections.defaultdict(list)
    for st, a in P.ASSAYS.items():
        if study and st != study:
            continue
        V = P.volumes(a)
        for (nm, nC, nCF, mw, ko, kc, ks) in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            for d in a["doses"]:
                if dose is not None and d != dose:
                    continue
                m = P.measured_x_cell(st, nm, d, V, a)
                if m is None:
                    continue
                r = (P.armitage_run(nm, lk, d, V) if model == "Armitage"
                     else P.qin_run(nm, d, V))
                per[nm].append(r["x_cell"] / m)
    return per


def from_csv(fn, col, study=None):
    per = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(DB, fn))):
        if study and r["study"] != study:
            continue
        per[r["compound"]].append(float(r[col]) / float(r["meas_x_cell_pct"]))
    return per


def header(w=44):
    print(f"  {'':{w}s}" + "".join(f"{c:>9s}" for c in NCF)
          + f"{'slope':>8s}{'r':>7s}")


def main():
    print(__doc__.split("Run:")[0].rstrip())

    # ------------------------------------------------------------------ 1
    print("\n" + "=" * 94)
    print("1.  IS IT REAL?  The same question asked inside one lab at a time")
    print("=" * 94)
    print("""
  Rosenmai is the only study that measures all five compounds, in one cell
  line, in one medium, at three doses.  If the trend is a pooling artefact it
  will vanish there.
""")
    header()
    for model in ("Armitage", "Qin"):
        for study in ("Rosenmai", "Amstutz", "Fragki"):
            row(f"{study} -- {model}", residuals(study, model=model))
        row(f"ALL 29 points -- {model}", residuals(model=model))
        print()
    print("""  VERDICT.  The trend is ARMITAGE'S, and it is real: Rosenmai alone gives
  +0.150 decades per CF2 at r = +0.70, with one lab, one cell line, one medium
  and one vehicle held fixed.  Qin has no trend inside a study (+0.022 at
  Rosenmai, and NEGATIVE at Amstutz), so Qin's pooled +0.092 is exactly the
  pooling artefact this test was built to catch -- Fragki contributes only
  PFOA and PFNA, at the largest residuals in the set.

  So from here the question is: why does Armitage's over-prediction grow with
  chain length where Qin's does not?""")

    # ------------------------------------------------------------------ 2
    print("\n" + "=" * 94)
    print("2.  IS IT CYTOTOXICITY AT THE LONG-CHAIN END?")
    print("=" * 94)
    print("""
  Long-chain PFAS are the cytotoxic ones.  Rosenmai report a significant drop
  in protein content at 100 uM PFDA and PFTeDA; Fragki excluded PFNA at 100 uM
  from their own QIVIVE for the same reason.  If dying cells were losing
  PFAS -- or losing the protein the uptake is normalised to -- the measured
  long-chain uptake would read low and the model would look like it
  over-predicts.  That hypothesis makes a sharp prediction: the trend should
  be worst at the HIGHEST dose and absent at the lowest.
""")
    header()
    for d in P.ASSAYS["Rosenmai"]["doses"]:
        row(f"Rosenmai, Armitage, {d:g} uM only", residuals("Rosenmai", d))
    print("""
  VERDICT: refuted, and backwards.  The trend is STRONGEST at 10 uM (+0.194)
  and weakest at 100 uM (+0.099) -- the opposite of a cytotoxicity artefact,
  which can only act where the cells are stressed.  It joins the wash
  hypothesis (section 11) on the list of measurement artefacts that do not
  survive their own direction test.""")

    # ------------------------------------------------------------------ 3
    print("\n" + "=" * 94)
    print("3.  PFBA HAS ENORMOUS LEVERAGE, AND IT IS THE WRONG KIND OF POINT")
    print("=" * 94)
    print("""
  PFBA sits at nCF 3, two carbons from the next compound in a five-point
  regression, so it swings the fitted slope hard.  It is also the one compound
  where WATER carries the denominator (80% of it, against 0.5% at PFNA), which
  makes it the compound least informative about the sorptive terms.  Every
  slope below is therefore reported twice: over all five, and over nCF 5-8.
""")
    header(38)
    for model in ("Armitage", "Qin"):
        for study in ("Rosenmai", None):
            per = residuals(study, model=model)
            sl5 = fit({c: v for c, v in per.items() if c != "PFBA"})[0]
            sl = row(f"{model} -- {study or 'all 29 points'}", per, 38)
            print(f"  {'':38s}{'':45s}  nCF 5-8 only: {sl5:+.3f}")
    print("""
  Armitage's rise over nCF 5-8 is +0.345 at Rosenmai and +0.316 pooled -- more
  than twice what the five-point fit reports, because PFBA sits ABOVE the line
  at the short end and drags the slope down.  Qin's is -0.074 and -0.043: no
  trend, or a slight fall.  The contrast between the two models is sharper, not
  weaker, once the leverage point is out.""")

    # ------------------------------------------------------------------ 4
    print("\n" + "=" * 94)
    print("4.  WHICH TERM CARRIES IT?  Substitute each, watch the SLOPE")
    print("=" * 94)
    print("""
  Sections 12-14 of the memo judged three substitutions on the LEVEL, where all
  three made the fit worse.  On the SLOPE they behave differently.
""")
    VAR = {"Armitage": [
        ("published", "mbm_fischer_albumin.csv", "armitage_published_x_cell_pct"),
        ("+ Fischer 2024 albumin", "mbm_fischer_albumin.csv",
         "armitage_plus_fischer_x_cell_pct"),
        ("+ Qin 2024 cell protein", "mbm_cell_protein.csv",
         "armitage_plus_qin_2024_cell_protein_x_cell_pct"),
        ("+ Qin 2023 measured D_cell/w", "mbm_cell_term.csv",
         "armitage_plus_measured_d_cell_w_x_cell_pct")],
        "Qin": [
        ("published", "mbm_fischer_albumin.csv", "qin_published_x_cell_pct"),
        ("+ Fischer 2024 albumin", "mbm_fischer_albumin.csv",
         "qin_plus_fischer_constant_d_x_cell_pct"),
        ("+ Qin 2024 cell protein", "mbm_cell_protein.csv",
         "qin_plus_qin_2024_cell_protein_x_cell_pct"),
        ("+ Qin 2023 measured D_cell/w", "mbm_cell_term.csv",
         "qin_plus_measured_d_cell_w_x_cell_pct")]}
    print(f"  {'':38s}{'all five':>12s}{'nCF 5-8':>11s}")
    for model in ("Armitage", "Qin"):
        for lab, fn, col in VAR[model]:
            per = from_csv(fn, col, "Rosenmai")
            a = fit(per)[0]
            b = fit({c: v for c, v in per.items() if c != "PFBA"})[0]
            print(f"  {model + ', ' + lab:38s}{a:>+12.3f}{b:>+11.3f}")
        print()
    print("""  VERDICT, on the nCF 5-8 column.  No single substitution zeroes the trend:

    + Fischer albumin        +0.345 -> +0.239   takes out about a third
    + Qin 2024 cell protein  +0.345 -> +0.248   takes out about a third
    + measured D_cell/w      +0.345 -> -0.243   over-corrects past zero

  (The five-compound column says albumin makes the slope STEEPER.  That is the
  PFBA leverage: Fischer's albumin is 10x below Armitage's at PFHxA but only
  2.9x below at PFNA, so the short chains gain most -- except PFBA, which is
  water-dominated and barely moves, so the five-point line tips the other way.
  An earlier draft of this section read the trap rather than the data.)""")

    # ------------------------------------------------------------------ 5
    print("\n" + "=" * 94)
    print("5.  THE ACTUAL MECHANISM: it is the RATIO, not either term alone")
    print("=" * 94)
    print("""
  Once albumin dominates the denominator -- 94-97% of it from PFOA up --

      x_cell  ~  D_cell . V_cell / (D_alb . V_alb)

  so the residual's chain-length slope is the chain-length slope of the RATIO
  D_cell/D_alb, and neither term alone predicts it.  Both models, as coded:
""")
    V = P.volumes(P.ASSAYS["Rosenmai"])
    print(f"  {'cmpd':8s}{'nCF':>4s} |{'Arm D_alb':>10s}{'Arm D_cell':>11s}"
          f"{'ratio':>8s} |{'Qin D_alb':>10s}{'Qin D_cell':>11s}{'ratio':>8s}"
          f" |{'Fischer BSA':>12s}")
    FIS = {5: 2.83, 6: 3.93, 7: 4.20, 8: 4.54}
    tab = []
    for (nm, nC, nCF, mw, ko, kc, ks) in Q.COMPOUNDS:
        lk = ko if ko is not None else Q.KOW_OPERA(nCF)
        D = P.armitage_D(lk, V["comp"])
        Db, Dsp, Dlip, _p = P.qin_D(nm, 1e-3)
        qcell = (Dsp * V["Vprot_cell"] + Dlip * V["Vlip_cell"]) / V["Vcells"]
        t = (nm, nCF, math.log10(D["Dsa"]), math.log10(D["Dcw"]),
             math.log10(Db), math.log10(qcell), FIS.get(nCF))
        tab.append(t)
        print(f"  {nm:8s}{nCF:>4d} |{t[2]:>10.2f}{t[3]:>11.2f}{t[3]-t[2]:>8.2f}"
              f" |{t[4]:>10.2f}{t[5]:>11.2f}{t[5]-t[4]:>8.2f}"
              f" |{('%.2f' % t[6]) if t[6] else '-':>12s}")

    def s58(idx):
        sub = [t for t in tab if t[1] >= 5]
        return float(np.polyfit([t[1] for t in sub],
                                [t[idx] if not callable(idx) else idx(t)
                                 for t in sub], 1)[0])
    arm_alb, arm_cell = s58(2), s58(3)
    qin_alb, qin_cell = s58(4), s58(5)
    fis = s58(6)
    print(f"""
  decades per CF2 over nCF 5-8:

    Armitage   albumin {arm_alb:+.3f}   cell {arm_cell:+.3f}   RATIO {arm_cell-arm_alb:+.3f}
    Qin        albumin {qin_alb:+.3f}   cell {qin_cell:+.3f}   RATIO {qin_cell-qin_alb:+.3f}
    Fischer 2024 measured BSA  {fis:+.3f}

  Compare those ratio slopes with the residual slopes over the same range:
  Armitage {arm_cell-arm_alb:+.3f} predicted against +0.345 observed; Qin {qin_cell-qin_alb:+.3f} against -0.074.
  The ratio IS the residual, to within the noise on a four-point fit.

  ARMITAGE GETS IT WRONG FROM BOTH ENDS, and the two errors stack:

    its albumin rises at {arm_alb:+.3f} where Fischer MEASURES {fis:+.3f}
        -> too shallow by {fis-arm_alb:.3f}, so the albumin sink grows too slowly
    its cell rises at {arm_cell:+.3f} where a flat ratio would need about {fis:+.3f}
        -> too steep by {arm_cell-fis:.3f}

  {(fis-arm_alb)/(arm_cell-arm_alb):.0%} of the gap is the albumin spLFER and {(arm_cell-fis)/(arm_cell-arm_alb):.0%} is the cell term.  Both have
  the same root: Armitage's piecewise albumin spLFER flattens to 0.37*logKow
  above logKow 4.5, and its cell is 99.7% membrane lipid at PFNA because
  "3.5% of cell protein behaves as octanol" puts the protein term 457-7,820x
  below the measured affinity (section 14).  A cell made of lipid climbs faster
  with chain length than a cell made of protein; an albumin spLFER that
  flattens climbs slower than albumin does.

  QIN GETS BOTH RIGHT.  Its albumin rises at {qin_alb:+.3f}, within {abs(qin_alb-fis):.3f} of Fischer's
  measurement, and its cell at {qin_cell:+.3f}, so the ratio is flat to {abs(qin_cell-qin_alb):.3f} per CF2.
  That is the headline: **Qin has the right chain-length physics and the wrong
  level; Armitage has the better level and the wrong chain-length physics.**
  The two failures are independent, which is why no single substitution fixes
  both and why their aggregate biases have been swapping places all along.

  WHAT WOULD SETTLE IT.  A measured D_cell/w for PFBA-PFNA on HepG2 or HepaRG
  rather than HEK293H homogenate, at a trace ligand:protein ratio, read
  alongside a measured albumin series on the same compounds.  The quantity to
  report is the RATIO and its slope, not either coefficient alone.""")

    dst = os.path.join(DB, "mbm_chain_length_residual.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cut", "model", "variant"] + list(NCF) + ["slope_per_CF2", "pearson_r"])
        for model in ("Armitage", "Qin"):
            for study in ("Rosenmai", "Amstutz", "Fragki", None):
                sl, r, g = fit(residuals(study, model=model))
                w.writerow([study or "all", model, "published"]
                           + [round(g[c], 4) if c in g else "" for c in NCF]
                           + [round(sl, 4) if sl is not None else "",
                              round(r, 3) if r is not None else ""])
            for d in P.ASSAYS["Rosenmai"]["doses"]:
                sl, r, g = fit(residuals("Rosenmai", d, model))
                w.writerow(["Rosenmai", model, f"published, {d:g} uM"]
                           + [round(g[c], 4) if c in g else "" for c in NCF]
                           + [round(sl, 4), round(r, 3)])
    print(f"\n  wrote {dst}")


if __name__ == "__main__":
    main()
