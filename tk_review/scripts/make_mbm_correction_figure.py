#!/usr/bin/env python3
"""Four pictures of why the two mass-balance models disagreed with the data.

A: the three uptake datasets, which do not agree with each other either
B: predicted vs measured, before and after the exposure-medium correction
C: where each model puts the dose -- the denominator, compound by compound
D: every parameter substitution tried, and which way each moved the fit

Run:  python3 scripts/make_mbm_correction_figure.py
"""
import copy
import csv
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

DB = os.path.join(HERE, "..", "db")
OUT = os.path.join(HERE, "..", "figures")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PLUM, SAND = "#8a5cd1", "#c8a33a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
EMPTY = "#eceae6"

STUDY_COLOR = {"Rosenmai": BLUE, "Amstutz": AQUA, "Fragki": ORANGE}
CMPDS = ("PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA")
NCF = {"PFBA": 3, "PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
})


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def predictions(fbs_fragki):
    """(study, compound, dose) -> (measured, qin, armitage), all fractions."""
    A = copy.deepcopy(P.ASSAYS)
    A["Fragki"]["FBS"] = fbs_fragki
    out = {}
    for study, a in A.items():
        V = P.volumes(a)
        for (nm, nC, nCF, mw, ko, kc, ks) in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                out[(study, nm, d)] = (
                    m, P.qin_run(nm, d, V)["x_cell"],
                    P.armitage_run(nm, lk, d, V)["x_cell"])
    return out


# ------------------------------------------------------------------ panel A
def panel_a(ax, now):
    for study in ("Rosenmai", "Amstutz", "Fragki"):
        for nm in CMPDS:
            pts = sorted((d, v[0] * 100) for (s, c, d), v in now.items()
                         if s == study and c == nm)
            if not pts:
                continue
            xs, ys = zip(*pts)
            ax.plot(xs, ys, "-o", color=STUDY_COLOR[study], lw=1.1, ms=4,
                    markeredgecolor="white", markeredgewidth=0.6, alpha=0.9)
            dy = {"PFBA": 0, "PFHxA": 0, "PFHpA": 0, "PFOA": 9,
                  "PFNA": -9}[nm] if study == "Rosenmai" else 2
            dx = 5 if study != "Amstutz" else -6
            ax.annotate(nm, (xs[-1], ys[-1]), (dx, dy),
                        textcoords="offset points", fontsize=7,
                        ha="left" if study != "Amstutz" else "right",
                        color=STUDY_COLOR[study])
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlabel("nominal concentration (µM)")
    ax.set_ylabel("measured: % of dose found in the cells")
    ax.set_title("A.  The three datasets do not agree either\n"
                 "same axis, same compounds, 14× apart at PFNA",
                 fontsize=10, color=INK, pad=8)
    for i, (s, lab) in enumerate((("Rosenmai", "Rosenmai 2018  HepG2, 10% FBS"),
                                  ("Amstutz", "Amstutz 2022  HepG2, 10% FBS"),
                                  ("Fragki", "Fragki 2022  HepaRG, 2% FBS"))):
        ax.text(0.03, 0.97 - 0.07 * i, lab, transform=ax.transAxes,
                fontsize=7.5, color=STUDY_COLOR[s], va="top", fontweight="bold")


# ------------------------------------------------------------------ panel B
def panel_b(ax, before, now):
    lo, hi = 0.02, 60.0
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1.0, zorder=1)
    for f in (10.0, 0.1):
        ax.plot([lo, hi], [lo * f, hi * f], color=EMPTY, lw=1.0, zorder=1)
    for key, (m, q, a) in now.items():
        study = key[0]
        c = STUDY_COLOR[study]
        if study == "Fragki":
            m0, q0, a0 = before[key]
            for y0, y1 in ((a0, a), (q0, q)):
                ax.annotate("", xy=(m * 100, y1 * 100), xytext=(m * 100, y0 * 100),
                            arrowprops=dict(arrowstyle="-|>", color=c, lw=0.8,
                                            alpha=0.55, shrinkA=1, shrinkB=2))
            ax.scatter([m * 100, m * 100], [a0 * 100, q0 * 100], s=14,
                       facecolor="none", edgecolor=c, linewidth=0.8, zorder=2)
        ax.scatter(m * 100, a * 100, s=34, c=c, marker="o",
                   edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter(m * 100, q * 100, s=34, c=c, marker="^",
                   edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(lo, 10), ax.set_ylim(lo, hi)
    ax.set_xlabel("measured % of dose in cells")
    ax.set_ylabel("predicted % of dose in cells")
    ax.set_title("B.  Predicted vs measured, 29 points\n"
                 "arrows: Fragki moving from 10% to its real 2% FBS",
                 fontsize=10, color=INK, pad=8)
    ax.text(0.97, 0.14, "● Armitage   ▲ Qin\nopen = as coded at 10% FBS",
            transform=ax.transAxes, fontsize=7.5, color=INK2,
            ha="right", va="bottom")
    ax.text(0.055, 0.93, "10× over", fontsize=7, color=MUTED,
            transform=ax.transAxes)
    ax.text(0.40, 0.055, "1:1", fontsize=7, color=MUTED,
            transform=ax.transAxes)


# ------------------------------------------------------------------ panel C
def panel_c(ax):
    bars, labels, groups = [], [], []
    for study, fbs, tag in (("Rosenmai", 0.10, "Rosenmai\n10% FBS"),
                            ("Amstutz", 0.10, "Amstutz\n10% FBS"),
                            ("Fragki", 0.10, "Fragki\nas coded 10%"),
                            ("Fragki", 0.02, "Fragki\nas dosed 2%")):
        a = copy.deepcopy(P.ASSAYS[study])
        a["FBS"] = fbs
        V = P.volumes(a)
        for nm in ("PFBA", "PFOA", "PFNA"):
            lk = next(c[4] if c[4] is not None else Q.KOW_OPERA(c[2])
                      for c in Q.COMPOUNDS if c[0] == nm)
            D = P.armitage_D(lk, V["comp"])
            parts = np.array([D["Dsa"] * V["Valb"], D["Dow"] * V["Vslip"],
                              V["Vw"], D["Dcw"] * V["Vcells"],
                              D["Dpl"] * V["sarea"] * 1000.0])
            bars.append(parts / parts.sum())
            labels.append(nm)
            groups.append(tag)
    bars = np.array(bars)
    y = np.arange(len(bars))[::-1]
    cols = (BLUE, PLUM, EMPTY, ORANGE, SAND)
    names = ("albumin", "serum lipid", "water", "cells", "plastic")
    left = np.zeros(len(bars))
    for j, (cname, col) in enumerate(zip(names, cols)):
        ax.barh(y, bars[:, j], left=left, color=col, height=0.74, label=cname)
        left += bars[:, j]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    seen = set()
    for yy, g in zip(y, groups):
        if g not in seen:
            ax.text(-0.30, yy - 1.0, g, fontsize=7.5, color=INK,
                    ha="left", va="center", fontweight="bold")
            seen.add(g)
    ax.set_xlim(0, 1), ax.set_ylim(-0.6, len(bars) - 0.4)
    ax.set_xlabel("share of the balance's sorptive capacity")
    ax.set_title("C.  Where the balance says the dose sits\n"
                 "albumin is 94–97% of it at PFOA and PFNA",
                 fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, ncol=5, loc="lower center",
              bbox_to_anchor=(0.42, -0.21), handletextpad=0.4, columnspacing=1.0)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


# ------------------------------------------------------------------ panel D
def panel_d(ax):
    def ladder_from(fn, meas_col, cols):
        rows = list(csv.DictReader(open(os.path.join(DB, fn))))
        out = {}
        for key, col in cols.items():
            out[key] = gm([float(r[col]) / float(r[meas_col]) for r in rows])
        return out

    fa = ladder_from("mbm_fischer_albumin.csv", "meas_x_cell_pct", {
        "qin_pub": "qin_published_x_cell_pct",
        "qin_fis": "qin_plus_fischer_constant_d_x_cell_pct",
        "arm_pub": "armitage_published_x_cell_pct",
        "arm_fis": "armitage_plus_fischer_x_cell_pct"})
    ct = ladder_from("mbm_cell_term.csv", "meas_x_cell_pct", {
        "qin_dcw": "qin_plus_measured_d_cell_w_x_cell_pct",
        "arm_dcw": "armitage_plus_measured_d_cell_w_x_cell_pct"})
    cp = ladder_from("mbm_cell_protein.csv", "meas_x_cell_pct", {
        "qin_prot": "qin_plus_qin_2024_cell_protein_x_cell_pct",
        "arm_prot": "armitage_plus_qin_2024_cell_protein_x_cell_pct"})

    before = predictions(0.10)
    qin_b = gm([v[1] / v[0] for v in before.values()])
    arm_b = gm([v[2] / v[0] for v in before.values()])

    rows = [("Fragki coded at 10% FBS\n(the old baseline)", qin_b, arm_b, True),
            ("published models,\ncorrected medium", fa["qin_pub"], fa["arm_pub"], True),
            ("+ Fischer 2024 albumin", fa["qin_fis"], fa["arm_fis"], False),
            ("+ Qin 2023 measured D$_{cell/w}$", ct["qin_dcw"], ct["arm_dcw"], False),
            ("+ Qin 2024 cell protein", cp["qin_prot"], cp["arm_prot"], False),
            ("+ 3.04 pL cell (HepG2)", 3.71, 1.84, False)]
    y = np.arange(len(rows))[::-1]
    h = 0.34
    ax.axvline(1.0, color=MUTED, lw=1.0)
    for i, (lab, qv, av, base) in enumerate(rows):
        yy = y[i]
        ax.barh(yy + h / 2, av, height=h, color=INK2,
                alpha=0.40 if base else 1.0,
                label="Armitage" if i == 1 else None)
        ax.barh(yy - h / 2, qv, height=h, color=PLUM,
                alpha=0.40 if base else 1.0,
                label="Qin" if i == 1 else None)
        ax.text(av * 1.06, yy + h / 2, f"{av:.1f}×", fontsize=7,
                color=INK2, va="center")
        ax.text(qv * 1.06, yy - h / 2, f"{qv:.1f}×", fontsize=7,
                color=INK2, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlim(0.8, 40)
    ax.set_xlabel("predicted / measured, geometric mean over 29 points")
    ax.set_title("D.  Every substitution tried, and which way it moved\n"
                 "1× would be a perfect model", fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", handletextpad=0.4)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def main():
    before, now = predictions(0.10), predictions(0.02)
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 11.4))
    panel_a(axes[0][0], now)
    panel_b(axes[0][1], before, now)
    panel_c(axes[1][0])
    panel_d(axes[1][1])
    fig.suptitle("Two in vitro mass-balance models against three uptake "
                 "datasets: what disagreed, and what the medium correction fixed",
                 fontsize=12, color=INK, y=0.975)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    path = os.path.join(OUT, "fig20_mbm_three_datasets.png")
    fig.savefig(path, dpi=190, bbox_inches="tight")
    print("wrote", os.path.abspath(path))


if __name__ == "__main__":
    main()
