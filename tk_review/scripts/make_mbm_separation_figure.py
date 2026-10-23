#!/usr/bin/env python3
"""The two failures separate: coefficients set the slope, volumes set the level.

A: the lipid-water partition ratio against chain length -- the single-descriptor
   spLFER against three independent measured series.  Two laboratories, two
   lipid preparations and two analytical routes agree with each other and not
   with the prediction.
B: every parameterisation tested in this work, placed on (chain-length slope,
   geometric-mean bias).  Substituting measured COEFFICIENTS moves points
   horizontally onto slope ~ 0; changing cell and medium VOLUMES moves them
   vertically.  Nothing reaches the target at (0, 1).

All numbers are recomputed here over nCF 5-8 so the two axes are on one
footing; PFBA is excluded throughout because it has no measured albumin value.

Run:  python3 scripts/make_mbm_separation_figure.py
"""
import collections
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P              # noqa: E402
import mbm_qin_vs_armitage as Q                  # noqa: E402
import mbm_measured_parameter_variant as MV      # noqa: E402
import mbm_fischer_escher_balance as FE          # noqa: E402

OUT = os.path.join(HERE, "..", "figures")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PLUM, SAND = "#8a5cd1", "#c8a33a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e6e4e0"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
})

NCF = {"PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}
# Chen 2025 SI Table S5, perfluoroalkyl carboxylates, log K_MW at pH 7.0.
CHEN = {3: 1.63, 4: 2.02, 5: 2.42, 6: 2.85, 7: 3.28, 8: 3.75, 9: 4.18,
        10: 4.50, 11: 4.82, 12: 5.11, 13: 5.36, 15: 5.73}
# Qin et al. 2024 Table 1: MEASURED liposome-water ratios, pH 7.4, from Droge
# et al. and Ebert et al.  These replace the Qin 2023 Figure A1 regression
# this panel used to carry: the regression was a stand-in for exactly these
# points.  The authors' own fit over their wider Droge+Ebert+Endo set is SI
# eq S9 (+0.562, R2 0.969) for carboxylates and eq S10 (+0.535, R2 0.977) for
# sulfonates.
QIN24_PFCA = {3: 1.00, 5: 2.32, 6: 2.91, 7: 3.52, 8: 4.25, 10: 4.54}
QIN24_PFSA = {6: 4.13, 8: 4.89}
QIN24_S9 = (0.562, -0.494)


def gm(xs):
    return 10 ** (sum(math.log10(x) for x in xs) / len(xs))


def fit(xs, ys):
    return np.polyfit(xs, ys, 1)


# ----------------------------------------------------------------- panel A
def panel_a(ax):
    n = np.linspace(2.6, 15.4, 200)
    lk = np.array([Q.COMPOUNDS[0][4] or 0])            # touch, keep import used
    # spLFER: log D_MW for the anion, as the model evaluates it
    spl = []
    for v in n:
        kow = Q.KOW_OPERA(v)
        spl.append(math.log10(P.armitage_D(kow, P.CELL_COMP["AREc2"])["Dmw"]))
    ax.plot(n, spl, color=ORANGE, lw=2.0, solid_capstyle="round", zorder=4)

    cx = sorted(CHEN)
    cy = [CHEN[k] for k in cx]
    sc, ic = fit(cx, cy)
    ax.plot(n, sc * n + ic, color=BLUE, lw=2.0, solid_capstyle="round", zorder=5)
    ax.plot(cx, cy, "o", ms=5.5, mfc=BLUE, mec=SURFACE, mew=1.4, zorder=6)

    qx = sorted(QIN24_PFCA)
    qy = [QIN24_PFCA[k] for k in qx]
    sq, iq = fit(qx, qy)
    ax.plot(n, sq * n + iq, color=AQUA, lw=2.0, ls=(0, (5, 2.2)),
            solid_capstyle="round", zorder=5)
    ax.plot(qx, qy, "^", ms=6.5, mfc=AQUA, mec=SURFACE, mew=1.4, zorder=6)

    ax.annotate(f"IV-MBM EQP v2.0 spLFER\n1.01·log K$_{{OW}}$ + 0.12\n"
                f"slope {np.polyfit(n, spl, 1)[0]:+.2f} per CF$_2$",
                xy=(7.6, 7.6 * np.polyfit(n, spl, 1)[0] + np.polyfit(n, spl, 1)[1]),
                xytext=(3.1, 10.4),
                color=ORANGE, fontsize=8, ha="left", va="top", linespacing=1.45,
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8,
                                shrinkA=2, shrinkB=3))
    ax.annotate(f"Chen 2025, supported bilayers\nmeasured, n = {len(cx)}   "
                f"slope {sc:+.3f}",
                xy=(12.9, 12.9 * sc + ic), xytext=(9.0, 3.05),
                color=BLUE, fontsize=8, ha="left", va="top", linespacing=1.45,
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.8,
                                shrinkA=2, shrinkB=2))
    ax.annotate(f"Droge / Ebert liposomes, measured\n"
                f"Qin 2024 Table 1, n = {len(qx)}   slope {sq:+.3f}\n"
                f"(authors' own fit, wider set: {QIN24_S9[0]:+.3f})",
                xy=(8.4, 8.4 * sq + iq), xytext=(2.9, 7.9),
                color=AQUA, fontsize=8, ha="left", va="top", linespacing=1.45,
                arrowprops=dict(arrowstyle="-", color=AQUA, lw=0.8,
                                shrinkA=2, shrinkB=3))

    ax.set_xlabel("perfluorinated carbons, n$_{CF}$")
    ax.set_ylabel("log lipid–water partition ratio")
    ax.set_xlim(2.4, 15.6)
    ax.set_ylim(0.8, 11.4)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("A   The lipid term, predicted against measured",
                 loc="left", fontsize=9.5, color=INK, pad=9)


# ----------------------------------------------------------------- panel B
def armitage_variant(comp_key, measured_alb, measured_mem, kmw):
    per = collections.defaultdict(list)
    for study, a in P.ASSAYS.items():
        aa = dict(a, line=comp_key)
        V = P.volumes(aa)
        for (nm, _nC, nCF, _mw, ko, _kc, _ks) in Q.COMPOUNDS:
            if nm not in NCF:
                continue
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            D = dict(P.armitage_D(lk, V["comp"]))
            if measured_alb:
                D["Dsa"] = 10 ** (MV.FISCHER_BSA[nCF] + MV.VOL_SHIFT)
            if measured_mem:
                D["Dmw"] = 10 ** kmw[nm]
                c = V["comp"]
                D["Dcw"] = (c["sl"] * D["Dow"] + c["ml"] * D["Dmw"]
                            + P.P_NLOM * c["prot"] * D["Dow"] + c["water"])
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, aa)
                if m is None:
                    continue
                ntot = d * V["Vbm"]
                S = (V["Vw"] + D["Dsa"] * V["Valb"] + D["Dow"] * V["Vslip"]
                     + D["Dcw"] * V["Vcells"] + D["Dpl"] * V["sarea"] * 1000.0)
                per[nm].append((D["Dcw"] * (ntot / S) * V["Vcells"] / ntot) / m)
    return summarise(per)


def qin_variant():
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
                per[nm].append(P.qin_run(nm, d, V)["x_cell"] / m)
    return summarise(per)


def fischer_variant(comp_key, use_medium, kmw):
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
                x, _dc, _dm = FE.fischer_run(study, nm, nCF, d, kmw,
                                             comp_key, use_medium)
                per[nm].append(x / m)
    return summarise(per)


def summarise(per):
    allr = [v for vs in per.values() for v in vs]
    xs = [NCF[c] for c in NCF if per.get(c)]
    ys = [math.log10(gm(per[c])) for c in NCF if per.get(c)]
    return gm(allr), float(fit(xs, ys)[0])


def panel_b(ax, kmw):
    pts = [
        ("Armitage, as published", armitage_variant("AREc2", 0, 0, kmw),
         ORANGE, "o", "published"),
        ("+ measured albumin", armitage_variant("AREc2", 1, 0, kmw),
         ORANGE, "^", "coef"),
        ("+ measured membrane", armitage_variant("AREc2", 0, 1, kmw),
         ORANGE, "v", "coef"),
        ("+ both", armitage_variant("AREc2", 1, 1, kmw),
         ORANGE, "s", "coef"),
        ("Qin, as published", qin_variant(), PLUM, "o", "published"),
        ("Fischer form, measured cells\n+ measured medium",
         fischer_variant("HepG2_Fischer", True, kmw), BLUE, "o", "published"),
        ("albumin-only medium",
         fischer_variant("HepG2_Fischer", False, kmw), BLUE, "D", "vol"),
        ("AREc2 cells",
         fischer_variant("AREc2", True, kmw), BLUE, "P", "vol"),
        ("AREc2 cells +\nalbumin-only medium",
         fischer_variant("AREc2", False, kmw), BLUE, "X", "vol"),
    ]
    ax.axhline(1.0, color=INK2, lw=1.0, ls=(0, (4, 2.5)), zorder=2)
    ax.axvline(0.0, color=INK2, lw=1.0, ls=(0, (4, 2.5)), zorder=2)
    ax.plot([0], [1], "*", ms=15, mfc="none", mec=INK, mew=1.3, zorder=7)
    ax.annotate("perfect\nmodel", xy=(0, 1), xytext=(0.055, 0.60),
                fontsize=7.6, color=INK2, ha="left", va="center",
                linespacing=1.4)

    # the band the measured-coefficient variants collapse into
    ax.axvspan(-0.17, 0.0, color=BLUE, alpha=0.055, zorder=1, lw=0)
    ax.text(-0.085, 1.72, "every variant with a measured\nlipid coefficient",
            fontsize=7.6, color=BLUE, ha="center", va="center",
            style="italic", linespacing=1.45)

    for lab, (bias, slope), col, mk, _kind in pts:
        ax.plot([slope], [bias], mk, ms=9, mfc=col, mec=SURFACE, mew=1.5,
                zorder=6)
    off = {
        "Armitage, as published": (-0.012, 1.0, "right", "center"),
        "+ measured albumin": (0.012, 1.0, "left", "center"),
        "+ measured membrane": (0.010, 0.80, "left", "top"),
        "+ both": (0.012, 1.0, "left", "center"),
        "Qin, as published": (0.012, 1.0, "left", "center"),
        "Fischer form, measured cells\n+ measured medium":
            (0.014, 1.0, "left", "center"),
        "albumin-only medium": (0.013, 1.0, "left", "center"),
        "AREc2 cells": (-0.013, 1.0, "right", "center"),
        "AREc2 cells +\nalbumin-only medium": (-0.013, 1.0, "right", "center"),
    }
    for lab, (bias, slope), col, _mk, _kind in pts:
        dx, dy, ha, va = off[lab]
        ax.annotate(lab, xy=(slope + dx, bias * dy), fontsize=7.6, color=col,
                    ha=ha, va=va, linespacing=1.4, zorder=8)

    ax.set_yscale("log")
    ax.set_xlabel("chain-length slope of the residual, decades per CF$_2$")
    ax.set_ylabel("geometric-mean bias, predicted / measured")
    ax.set_xlim(-0.24, 0.40)
    ax.set_ylim(0.09, 46)
    ax.set_yticks([0.1, 0.3, 1, 3, 10, 30])
    ax.set_yticklabels(["0.1×", "0.3×", "1×", "3×", "10×", "30×"])
    ax.grid(color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("B   Every parameterisation tested, nCF 5–8",
                 loc="left", fontsize=9.5, color=INK, pad=9)
    return pts


def main():
    kmw = FE.chen_kmw()
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.0))
    panel_a(axes[0])
    pts = panel_b(axes[1], kmw)
    fig.text(0.008, 0.978,
             "Measured partition coefficients remove the chain-length error; "
             "nothing tested here removes the level error",
             fontsize=11.5, color=INK, va="top", fontweight="semibold")
    fig.text(0.008, 0.105,
             "A  log D$_{MW}$ for the anion at pH 7.4 as IV-MBM EQP v2.0 "
             "evaluates it, against two measured series: Chen 2025 SI Table "
             "S5 (solid-supported bilayers, pH 7.0) and the liposome-water "
             "ratios Qin 2024 Table 1 reports\n"
             "    from Droge et al. and Ebert et al. (pH 7.4). Lines are "
             "least-squares fits to the plotted points.\n"
             "B  bias and slope recomputed over nCF 5–8 (PFBA excluded: no "
             "measured albumin value) across the three uptake datasets, 24 "
             "observations. Horizontal moves come from coefficients, vertical "
             "moves from volumes and compositions.\n"
             "    Every variant whose lipid coefficient is measured — Qin 2023 "
             "included, which takes its own from Ebert 2020 — falls in the "
             "shaded band; the two that predict the lipid term from "
             "log K$_{OW}$ do not.",
             fontsize=7.6, color=MUTED, va="top", linespacing=1.6)
    fig.tight_layout(rect=(0.0, 0.155, 1.0, 0.94))
    path = os.path.join(OUT, "fig22_mbm_separation.png")
    fig.savefig(path, dpi=200)
    print("wrote", os.path.abspath(path))
    print(f"\n  {'parameterisation':46s}{'bias':>9s}{'slope':>9s}")
    for lab, (b, sl), _c, _m, _k in pts:
        print(f"  {lab.replace(chr(10), ' '):46s}{b:>8.2f}x{sl:>+9.3f}")


if __name__ == "__main__":
    main()
