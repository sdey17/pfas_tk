#!/usr/bin/env python3
"""The ligand:protein regime gap between in vitro PFAS assays and human serum.

A: albumin site occupancy against nominal dose for the three uptake assays,
   with the physiological band marked
B: the apparent albumin distribution ratio collapsing from its trace limit as
   the sites fill, for the four compounds Qin measured

Run:  python3 scripts/make_mbm_regime_figure.py
"""
import copy
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P      # noqa: E402
import mbm_qin_vs_armitage as Q          # noqa: E402

OUT = os.path.join(HERE, "..", "figures")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PLUM, SAND = "#8a5cd1", "#c8a33a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
EMPTY = "#eceae6"
STUDY_COLOR = {"Rosenmai": BLUE, "Amstutz": AQUA, "Fragki": ORANGE}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
})


def albumin_uM(fbs):
    return P.SI["conc_ser_alb"] * fbs / Q.BSA_MW * 1e6


def panel_a(ax):
    """Nominal ligand:albumin ratio against dose, per assay."""
    doses = np.logspace(math.log10(0.0005), math.log10(200), 300)
    for study, a in P.ASSAYS.items():
        alb = albumin_uM(a["FBS"])
        # Rosenmai and Amstutz are at the same FBS, so their lines coincide;
        # dash the second so both remain visible.
        ls = "--" if study == "Amstutz" else "-"
        ax.plot(doses, doses / alb, color=STUDY_COLOR[study], lw=1.4, ls=ls,
                label=f"{study}  ({a['FBS']:.0%} FBS, {alb:.1f} µM albumin)")
        lo, hi = min(a["doses"]), max(a["doses"])
        ax.plot([lo, hi], [lo / alb, hi / alb], color=STUDY_COLOR[study],
                lw=5.5, alpha=0.30, solid_capstyle="butt", zorder=1)
        ax.scatter([lo, hi], [lo / alb, hi / alb], s=26,
                   c=STUDY_COLOR[study], edgecolor="white", linewidth=0.7,
                   zorder=4)
    ax.axhspan(1e-4, 4e-3, color=AQUA, alpha=0.13, zorder=0)
    ax.axhline(1.0, color=MUTED, lw=1.0, ls="--")
    ax.text(0.0009, 1.35, "one PFAS per albumin molecule  (n ≈ 1)",
            fontsize=7.5, color=INK2)
    ax.text(0.0009, 1.4e-3, "physiological range\n(Fischer 2024: ≤ 0.004;\nhuman serum ≈ 0.0005)",
            fontsize=7.5, color="#0e7a55", va="center")
    ax.annotate("", xy=(150, 2.0), xytext=(150, 0.0012),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.0))
    ax.text(120, 0.055, "3–4 orders\nof magnitude", fontsize=7.5, color=INK2,
            ha="right", va="center")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(5e-4, 200), ax.set_ylim(8e-5, 40)
    ax.set_xlabel("nominal PFAS concentration (µM)")
    ax.set_ylabel("PFAS : albumin  (mol/mol)")
    ax.set_title("A.  The assays sit far outside the regime they inform\n"
                 "bold segments = the doses actually used",
                 fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right",
              handletextpad=0.5)


def panel_b(ax):
    """D_BSA/w against free concentration, Qin's four measured compounds."""
    cw = np.logspace(-4, 3, 400)
    for nm, col in (("PFOS", PLUM), ("PFHxS", BLUE),
                    ("PFOA", ORANGE), ("PFBA", SAND)):
        Kd, Cmax, _n, _lk, lDns, _a, _b = Q.QIN_BSA_MEASURED[nm]
        d = Cmax / (Kd + cw) + 10 ** lDns
        trace = Cmax / Kd + 10 ** lDns
        ax.plot(cw, d / trace, color=col, lw=1.5, label=f"{nm}  (K$_d$ {Kd:.3g} µM)")
        ax.scatter([Kd], [(Cmax / (2 * Kd) + 10 ** lDns) / trace], s=28, c=col,
                   edgecolor="white", linewidth=0.7, zorder=4)
    # where the real wells sit
    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        xs, ys = [], []
        for d in a["doses"]:
            if P.measured_x_cell(study, "PFOA", d, V, a) is None:
                continue
            r = P.qin_run("PFOA", d, V)
            p = Q.qin_params("PFOA", 7)
            trace = p["Cmax"] / p["Kd"] + 10 ** p["logDnonspec"]
            xs.append(r["Cw"]), ys.append(r["Dbsa"] / trace)
        if xs:
            ax.scatter(xs, ys, s=52, marker="v", c=STUDY_COLOR[study],
                       edgecolor=INK, linewidth=0.8, zorder=5)
    ax.axvspan(1e-4, 2e-3, color=AQUA, alpha=0.13, zorder=0)
    ax.text(1.25e-4, 1.45, "free conc.\nin human serum",
            fontsize=7.5, color="#0e7a55", va="top")
    ax.text(0.9, 1.35, "▼ = PFOA in the three assays,\n     at the doses used",
            fontsize=7.5, color=INK2, va="top")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(1e-4, 1e3), ax.set_ylim(0.03, 2.2)
    ax.set_xlabel("freely dissolved concentration (µM)")
    ax.set_ylabel("D$_{BSA/w}$  ÷  its trace-limit value")
    ax.set_title("B.  Where the binding 'constant' stops being constant\n"
                 "markers on each curve = K$_d$, i.e. half the sites filled",
                 fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, loc="lower left", handletextpad=0.5)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6))
    panel_a(axes[0])
    panel_b(axes[1])
    fig.suptitle("Figure 1.  In vitro PFAS uptake assays operate three to four orders "
                 "of magnitude above the physiological ligand:protein ratio",
                 fontsize=11, color=INK, y=1.01)
    fig.tight_layout()
    path = os.path.join(OUT, "fig21_mbm_regime_gap.png")
    fig.savefig(path, dpi=190, bbox_inches="tight")
    print("wrote", os.path.abspath(path))


if __name__ == "__main__":
    main()
