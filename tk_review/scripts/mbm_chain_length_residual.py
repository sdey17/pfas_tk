#!/usr/bin/env python3
"""What is left after the medium correction: a residual that rises with chain
length.  This localises it.

THE OBSERVATION.  On the corrected exposure media (section 15 of the memo),
Armitage's over-prediction climbs monotonically with perfluorocarbon chain
length -- +0.172 decades per CF2 over the pooled 29 points, r = +0.83.  A
level error would not do that.  Something in the model's chain-length
dependence is wrong, or something in the measurements is.

FOUR QUESTIONS, IN ORDER.  Each one can kill a candidate outright:

  1. Is the trend real, or an artefact of pooling three studies that measure
     different compound sets?  (Fragki contributes only PFOA and PFNA, and
     has the largest residuals, so pooling alone could manufacture a slope.)
  2. Is it a measurement artefact at the long-chain end -- cytotoxicity
     suppressing the measured uptake?
  3. Which term in the balance carries it?  Substitute each candidate in turn
     and watch the SLOPE, not the level.
  4. What would the right value have to be, and does anything measure it?

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
    print("3.  WHICH TERM CARRIES IT?  Substitute each, watch the SLOPE")
    print("=" * 94)
    print("""
  Section 12-14 substituted three measured parameters and judged each on the
  LEVEL, where all three made the fit worse.  Judged on the SLOPE they behave
  completely differently, and that is what localises the defect.
""")
    for study, tag in ((None, "all 29 points"), ("Rosenmai", "Rosenmai only")):
        print(f"\n  --- {tag} ---")
        header()
        for model, pub, fis, dcw, prot in (
            ("Armitage", "armitage_published_x_cell_pct",
             "armitage_plus_fischer_x_cell_pct",
             "armitage_plus_measured_d_cell_w_x_cell_pct",
             "armitage_plus_qin_2024_cell_protein_x_cell_pct"),
            ("Qin", "qin_published_x_cell_pct",
             "qin_plus_fischer_constant_d_x_cell_pct",
             "qin_plus_measured_d_cell_w_x_cell_pct",
             "qin_plus_qin_2024_cell_protein_x_cell_pct")):
            row(f"{model}, published",
                from_csv("mbm_fischer_albumin.csv", pub, study))
            row(f"{model} + Fischer 2024 albumin",
                from_csv("mbm_fischer_albumin.csv", fis, study))
            row(f"{model} + Qin 2024 cell protein",
                from_csv("mbm_cell_protein.csv", prot, study))
            row(f"{model} + Qin 2023 measured D_cell/w",
                from_csv("mbm_cell_term.csv", dcw, study))
            print()
    print("""  VERDICT.  Only one substitution touches the slope, and it does not nudge it,
  it INVERTS it:

    albumin        makes the slope STEEPER (+0.150 -> +0.212 at Rosenmai).
                   So the albumin term is not carrying the trend -- and
                   Fischer's measurement makes the chain-length story worse,
                   not better, which is its own small result.
    cell protein   barely moves it (+0.150 -> +0.067), and wrecks the level.
    D_cell/w       flips the sign outright: +0.150 -> -0.321, r = -0.97.

  A parameter that can swing the residual slope from +0.15 to -0.32 is the
  parameter the residual slope lives in.  The chain-length defect is in the
  CELL PARTITION COEFFICIENT, specifically in how fast it rises per CF2.""")

    # ------------------------------------------------------------------ 4
    print("\n" + "=" * 94)
    print("4.  THE SLOPE LEDGER: five estimates of one number")
    print("=" * 94)
    print("""
  d log D_cell/w per CF2, every independent estimate in hand:

    Qin 2023 Table 2, cell homogenate, measured      +0.177   (4 compounds)
    inverted from the uptake data themselves         +0.407   (this work)
    Ryu 2024 whole-tissue f_unbound                  +0.46 to +0.53
    Qin 2023 model, as parameterised                 +0.515
    Armitage 2021 model, as parameterised            +0.782

  Armitage sits 0.38 above the uptake requirement and 0.27 above the next
  highest estimate.  Qin sits 0.11 above it.  Their residual slopes order the
  same way: +0.150 and +0.022.  Qin's own Table 2 sits alone at the bottom,
  which is why substituting it overshoots into a negative slope rather than
  landing on zero -- the truth is between the model and that measurement.

  WHY ARMITAGE IS STEEP, mechanistically.  Its cell is 79% water at PFBA and
  99.7% MEMBRANE LIPID at PFNA, with protein contributing 0.1-0.2% throughout,
  because the convention "3.5% of cell protein behaves as octanol" puts the
  protein term 457-7,820x below the measured protein affinity (section 14).
  Define protein out of the cell and what is left is a lipid partition
  coefficient, and lipid affinity climbs with chain length faster than protein
  affinity does.  Qin's cell is 92-95% protein at every chain length and has
  essentially the right slope.

  So the same convention produces BOTH failures: the level error of section 14
  and the chain-length error here.  They are one defect seen twice.

  WHAT WOULD SETTLE IT.  A measured D_cell/w for PFBA through PFNA on HepG2 or
  HepaRG, not HEK293H homogenate, at a trace ligand:protein ratio.  The
  disagreement to resolve is +0.177 (Qin's homogenate) against +0.407 (what
  the uptake data require); everything else in the ledger clusters between
  +0.41 and +0.53.  Until that exists the slope is bounded, not known.""")

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
