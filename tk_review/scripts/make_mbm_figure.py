"""
Figure: Qin 2023 vs Armitage 2021 (IV-MBM EQP v2.0) in-vitro mass-balance
models for PFBA, PFHxA, PFHpA, PFOA, PFNA.

Four panels, each carrying one argument:
  A  measured vs predicted C_cell -- how far off both models are
  B  the Qin/Armitage ratio vs chain length -- WHERE they diverge and that
     they cross over between PFHpA and PFOA
  C  log D_cell/w vs chain length -- WHY: a slope disagreement, Qin measured
     0.18/CF2 against Armitage's constructed ~0.94/CF2
  D  Rosenmai's own cell-concentration axis read two ways

Palette follows the project's existing figure style; every series is directly
labelled as well as coloured.

Run:  python3 tk_review/scripts/make_mbm_figure.py
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_qin_vs_armitage as M

FIGS = os.path.join(os.path.dirname(HERE), "figures")
BLUE, ORANGE, AQUA, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#7b5ea7"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID, RULE = "#e6e5e1", "#d8d7d3"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": RULE, "axes.linewidth": 0.8, "axes.labelcolor": INK2,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
})

STUDY_MARK = {"Rosenmai2017": "o", "Amstutz2022": "s", "Fragki2023": "^"}
ORDER = [c[0] for c in M.COMPOUNDS]
NCF = {c[0]: c[2] for c in M.COMPOUNDS}
MW = {c[0]: c[3] for c in M.COMPOUNDS}


def collect():
    """measured + predicted C_cell per (compound, study, dose)."""
    pts = []
    for (name, nC, nCF, mw, ko, kc, ksrc) in M.COMPOUNDS:
        lk = ko if ko is not None else M.KOW_OPERA(nCF)
        qp = M.qin_params(name, nCF)
        for study, a in M.ASSAYS.items():
            V = M.volumes(a)
            for dose in a["doses"]:
                if study == "Rosenmai2017":
                    up = M.ROSENMAI_UPTAKE_PCT.get(name, {}).get(int(dose))
                    meas = M.ccell_from_fraction(up / 100.0, dose, V) if up else None
                elif study == "Amstutz2022":
                    v = M.AMSTUTZ_IF_PCT.get(name, (None, None))[0]
                    meas = M.ccell_from_fraction(v / 100.0, dose, V) if v else None
                else:
                    v = M.FRAGKI_NMOL_PER_MG.get(name, {}).get(dose)
                    meas = M.ccell_from_nmol_per_mg(v) if v else None
                if meas is None:
                    continue
                pts.append(dict(
                    compound=name, nCF=nCF, study=study, dose=dose, meas=meas,
                    qin=M.qin_predict(name, nCF, dose, V, "mechanistic")["C_cell"],
                    arm=M.armitage_predict(lk, dose, V)["C_cell"],
                    src=qp["source"]))
    return pts


def main():
    pts = collect()
    fig, axes = plt.subplots(1, 4, figsize=(13.6, 3.5))
    axA, axB, axC, axD = axes

    # ---- A: measured vs predicted -------------------------------------
    lo, hi = 3.0, 3e5
    axA.plot([lo, hi], [lo, hi], "-", color=INK2, lw=1.0, zorder=1)
    for f, lab in ((10, "10x"), (100, "100x")):
        axA.plot([lo, hi], [lo * f, hi * f], ":", color=RULE, lw=0.9, zorder=1)
        axA.text(6.0, 6.0 * f, lab, fontsize=6.5, color=MUTED)
    for p in pts:
        mk = STUDY_MARK[p["study"]]
        open_ = p["src"] != "measured"
        axA.plot(p["meas"], p["qin"], mk, ms=5.5, mew=1.1,
                 mfc="none" if open_ else ORANGE, mec=ORANGE, zorder=3)
        axA.plot(p["meas"], p["arm"], mk, ms=5.5, mew=1.1,
                 mfc="none" if open_ else BLUE, mec=BLUE, zorder=3)
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlim(lo, hi); axA.set_ylim(lo, hi)
    axA.set_xlabel("measured C$_{cell}$  (µM)")
    axA.set_ylabel("predicted C$_{cell}$  (µM)")
    axA.text(0.03, 0.95, "A", transform=axA.transAxes, fontsize=11, fontweight="bold")
    axA.text(0.97, 0.12, "Qin", transform=axA.transAxes, ha="right",
             color=ORANGE, fontsize=9, fontweight="bold")
    axA.text(0.97, 0.04, "Armitage", transform=axA.transAxes, ha="right",
             color=BLUE, fontsize=9, fontweight="bold")
    axA.text(0.03, 0.86, "open = extrapolated\nparameters", transform=axA.transAxes,
             fontsize=6.5, color=MUTED, va="top")
    axA.text(0.60, 0.33, "o Rosenmai\n□ Amstutz\n△ Fragki", transform=axA.transAxes,
             fontsize=6.5, color=MUTED, va="top")

    # ---- B: the crossover ---------------------------------------------
    axB.axhline(1.0, color=INK2, lw=1.0)
    axB.axhspan(1.0, 1e3, color=ORANGE, alpha=0.055, lw=0)
    axB.axhspan(1e-3, 1.0, color=BLUE, alpha=0.055, lw=0)
    for st, mk in STUDY_MARK.items():
        xs = [NCF[p["compound"]] for p in pts if p["study"] == st]
        ys = [p["qin"] / p["arm"] for p in pts if p["study"] == st]
        axB.plot(xs, ys, mk, ms=5.5, mfc="none", mec=PURPLE, mew=1.2, zorder=3)
    nn = np.array([NCF[c] for c in ORDER], float)
    rat = []
    for c in ORDER:
        qp = M.qin_params(c, NCF[c])
        lk = dict((x[0], x[4]) for x in M.COMPOUNDS)[c]
        lk = lk if lk is not None else M.KOW_OPERA(NCF[c])
        rat.append(10 ** qp["logDcellw"] * M.F_PROTLIP / M.armitage_drs(lk)["DR_kcw"])
    axB.plot(nn, rat, "-", color=PURPLE, lw=1.6, zorder=2)
    axB.set_yscale("log"); axB.set_ylim(2e-2, 4e2)
    axB.set_xticks(nn); axB.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axB.set_ylabel("Qin / Armitage")
    axB.set_xlabel("")
    axB.text(0.03, 0.95, "B", transform=axB.transAxes, fontsize=11, fontweight="bold")
    axB.text(0.96, 0.93, "Qin higher", transform=axB.transAxes, ha="right",
             color=ORANGE, fontsize=7.5, fontweight="bold")
    axB.text(0.04, 0.05, "Armitage higher", transform=axB.transAxes, ha="left",
             color=BLUE, fontsize=7.5, fontweight="bold")
    axB.text(0.5, 0.52, "cross over", transform=axB.transAxes, ha="center",
             color=INK2, fontsize=7)

    # ---- C: why -- a slope disagreement --------------------------------
    grid = np.linspace(2.5, 8.5, 60)
    arm_line = [math.log10(M.armitage_drs(M.KOW_OPERA(n))["DR_kcw"]) for n in grid]
    axC.plot(grid, arm_line, "-", color=BLUE, lw=1.8, label="Armitage")
    shift = math.log10(M.F_PROTLIP)   # Qin is per L cell prot+lip -> whole cell
    a, b = (M.QIN_T2["PFBA"][1], M.QIN_T2["PFOA"][1])
    qs = (b - a) / 4.0
    axC.plot(grid, [a + qs * (n - 3) + shift for n in grid], "-", color=ORANGE, lw=1.8)
    for c in ("PFBA", "PFOA"):
        axC.plot(NCF[c], M.QIN_T2[c][1] + shift, "o", ms=7, mfc=ORANGE,
                 mec=ORANGE, zorder=4)
    for c in ("PFHxA", "PFHpA", "PFNA"):
        axC.plot(NCF[c], M.qin_params(c, NCF[c])["logDcellw"] + shift, "o", ms=6.5,
                 mfc="none", mec=ORANGE, mew=1.3, zorder=4)
    for c in ORDER:
        lk = dict((x[0], x[4]) for x in M.COMPOUNDS)[c]
        lk = lk if lk is not None else M.KOW_OPERA(NCF[c])
        axC.plot(NCF[c], math.log10(M.armitage_drs(lk)["DR_kcw"]), "s", ms=6,
                 mfc=BLUE if c in ("PFBA", "PFOA") else "none", mec=BLUE, mew=1.3, zorder=4)
    axC.set_xticks(nn); axC.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axC.set_ylabel("log D$_{cell/w}$")
    axC.set_ylim(-0.3, 5.2)
    axC.axvline(5.24, color=MUTED, lw=0.9, ls=":")
    axC.text(0.03, 0.95, "C", transform=axC.transAxes, fontsize=11, fontweight="bold")
    axC.text(0.42, 0.93, "Armitage\n+0.75 / CF$_2$", transform=axC.transAxes,
             ha="center", va="top", color=BLUE, fontsize=7.5, fontweight="bold")
    axC.text(0.62, 0.30, "Qin, measured\n+0.18 / CF$_2$", transform=axC.transAxes,
             ha="center", color=ORANGE, fontsize=7.5, fontweight="bold")
    axC.text(0.97, 0.06, "whole-cell basis", transform=axC.transAxes, ha="right",
             fontsize=6.5, color=MUTED)

    # ---- D: Rosenmai's axis read two ways ------------------------------
    V = M.volumes(M.ASSAYS["Rosenmai2017"])
    mg = V["Vcells"] * 1000 * M.MG_PROTEIN_PER_ML_CELL
    w = 0.34
    for i, c in enumerate(ORDER):
        up = M.ROSENMAI_UPTAKE_PCT[c][100]
        nmol_mg = up / 100.0 * 100.0 * V["Vbm"] * 1000 / mg
        good = nmol_mg * M.MG_PROTEIN_PER_ML_CELL
        bad = good * MW[c]
        axD.bar(i - w / 2, good, w, color=AQUA, lw=0)
        axD.bar(i + w / 2, bad, w, color=ORANGE, lw=0)
        axD.text(i + w / 2, bad * 1.35, "%.0fx" % MW[c], ha="center",
                 fontsize=6.5, color=ORANGE, fontweight="bold")
    axD.set_yscale("log"); axD.set_ylim(50, 1e7)
    axD.axhline(1e6, color=INK2, lw=0.9, ls="--")
    axD.text(0.99, 0.815, "1 M — physically impossible", transform=axD.transAxes,
             ha="right", va="bottom", fontsize=6.5, color=INK2)
    axD.set_xticks(range(len(ORDER)))
    axD.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axD.set_ylabel("C$_{cell}$ at 100 µM  (µM)")
    axD.text(0.03, 0.95, "D", transform=axD.transAxes, fontsize=11, fontweight="bold")
    axD.text(0.11, 0.985, "\u25a0 read as ng/mg  (Abstract)", transform=axD.transAxes,
             color=AQUA, fontsize=7.5, fontweight="bold", va="top")
    axD.text(0.11, 0.905, "\u25a0 read as nmol/mg  (Results 3.3, Fig 2A)",
             transform=axD.transAxes, color=ORANGE, fontsize=7.5,
             fontweight="bold", va="top")


    for ax in axes:
        ax.grid(True, which="major", color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)

    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.235, top=0.90, wspace=0.33)
    out = os.path.join(FIGS, "fig10_mbm_qin_vs_armitage.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)
    # guard against the tight_layout/axes-collapse failure seen on fig 4
    for i, ax in enumerate(axes):
        h = ax.get_position().height
        assert h > 0.45, "panel %d collapsed to %.2f of figure height" % (i, h)
    print("panel heights OK:", ["%.2f" % ax.get_position().height for ax in axes])


if __name__ == "__main__":
    main()
