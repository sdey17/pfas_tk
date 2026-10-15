"""
Figure: Qin 2023 vs Armitage 2021 IV-MBM EQP v2.0, from the PUBLISHED
equations and each paper's own supplementary parameters.

  A  measured vs predicted fraction of dose in cells -- the axis that does not
     depend on the 3.3x disagreement about cell VOLUME
  B  where each model puts the PFAS inside the cell: Qin in protein, Armitage
     in membrane lipid (or, at the short end, in cell water)
  C  dose exponent b in C_cell ~ C_nom^b.  Armitage is pinned at 1.00
  D  Rosenmai's cell-concentration axis read two ways

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
import mbm_published_equations as P
import mbm_qin_vs_armitage as Q

FIGS = os.path.join(os.path.dirname(HERE), "figures")
BLUE, ORANGE, AQUA, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#7b5ea7"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID, RULE = "#e6e5e1", "#d8d7d3"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": RULE, "axes.linewidth": 0.8, "axes.labelcolor": INK2,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False})

ORDER = [c[0] for c in Q.COMPOUNDS]
NCF = {c[0]: c[2] for c in Q.COMPOUNDS}
MW = {c[0]: c[3] for c in Q.COMPOUNDS}
MARK = {"Rosenmai": "o", "Amstutz": "s", "Fragki": "^"}


def lkow(c):
    v = dict((x[0], x[4]) for x in Q.COMPOUNDS)[c]
    return v if v is not None else Q.KOW_OPERA(NCF[c])


def slope(pairs):
    n = len(pairs)
    lx = [math.log10(x) for x, _ in pairs]
    ly = [math.log10(y) for _, y in pairs]
    sx, sy = sum(lx), sum(ly)
    return (n * sum(a * b for a, b in zip(lx, ly)) - sx * sy) / \
           (n * sum(v * v for v in lx) - sx ** 2)


def main():
    fig, (axA, axB, axC, axD) = plt.subplots(1, 4, figsize=(13.6, 3.6))

    # ---- A ----
    lo, hi = 2e-2, 2e1
    axA.plot([lo, hi], [lo, hi], "-", color=INK2, lw=1.0, zorder=2)
    for f, lab in ((10, "10x"), (0.1, "0.1x")):
        axA.plot([lo, hi], [lo * f, hi * f], ":", color=RULE, lw=0.9, zorder=1)
    for st, a in P.ASSAYS.items():
        V = P.volumes(a)
        for c in ORDER:
            for d in a["doses"]:
                m = P.measured_x_cell(st, c, d, V, a)
                if m is None:
                    continue
                op = Q.qin_params(c, NCF[c])["source"] != "measured"
                axA.plot(m * 100, P.qin_run(c, d, V)["x_cell"] * 100, MARK[st],
                         ms=5.5, mew=1.1, mfc="none" if op else ORANGE,
                         mec=ORANGE, zorder=3)
                axA.plot(m * 100, P.armitage_run(c, lkow(c), d, V)["x_cell"] * 100,
                         MARK[st], ms=5.5, mew=1.1,
                         mfc="none" if op else BLUE, mec=BLUE, zorder=3)
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlim(lo, hi); axA.set_ylim(lo, hi)
    axA.set_xlabel("measured  % of dose in cells")
    axA.set_ylabel("predicted  % of dose in cells")
    axA.text(0.03, 0.95, "A", transform=axA.transAxes, fontsize=11, fontweight="bold")
    axA.text(0.97, 0.13, "Qin", transform=axA.transAxes, ha="right",
             color=ORANGE, fontsize=9, fontweight="bold")
    axA.text(0.97, 0.05, "Armitage", transform=axA.transAxes, ha="right",
             color=BLUE, fontsize=9, fontweight="bold")
    axA.text(0.035, 0.90, "open = Qin params extrapolated",
             transform=axA.transAxes, fontsize=6.4, color=MUTED, va="top")
    axA.text(0.54, 0.30, "o Rosenmai\n□ Amstutz\n△ Fragki",
             transform=axA.transAxes, fontsize=6.4, color=MUTED, va="top")

    # ---- B ----
    x = np.arange(len(ORDER))
    w = 0.38
    V = P.volumes(P.ASSAYS["Rosenmai"]); comp = V["comp"]
    for i, c in enumerate(ORDER):
        a = P.armitage_D(lkow(c), comp)
        ml = comp["ml"] * a["Dmw"]; wa = comp["water"]
        ot = comp["sl"] * a["Dow"] + P.P_NLOM * comp["prot"] * a["Dow"]
        t = ml + wa + ot
        b = 0.0
        for frac, col in ((ml / t, BLUE), (wa / t, AQUA), (ot / t, MUTED)):
            axB.bar(i - w / 2, frac * 100, w, bottom=b * 100, color=col, lw=0); b += frac
        Db, Ds, Dl, p = P.qin_D(c, 1.0)
        qp = Ds * comp["prot"]; ql = Dl * (comp["sl"] + comp["ml"]); qt = qp + ql
        b = 0.0
        for frac, col in ((qp / qt, ORANGE), (ql / qt, PURPLE)):
            axB.bar(i + w / 2, frac * 100, w, bottom=b * 100, color=col, lw=0); b += frac
    axB.set_xticks(x); axB.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axB.set_ylabel("share of cell sorptive capacity (%)")
    axB.set_ylim(0, 168)
    axB.axhline(100, color=RULE, lw=0.8)
    axB.text(0.03, 0.74, "B", transform=axB.transAxes, fontsize=11, fontweight="bold")
    axB.text(0.5, 1.045, "left bar Armitage      right bar Qin", transform=axB.transAxes,
             ha="center", fontsize=7, color=INK2)
    for lab, col, xx, yy in (("membrane lipid", BLUE, 0.03, 0.965),
                             ("cell water", AQUA, 0.03, 0.885),
                             ("protein (D$_{non-spec}$)", ORANGE, 0.47, 0.965),
                             ("lipid", PURPLE, 0.47, 0.885)):
        axB.text(xx, yy, "■ " + lab, transform=axB.transAxes,
                 fontsize=6.6, color=col, fontweight="bold", va="top")

    # ---- C ----
    series = [("Rosenmai", "PFBA"), ("Rosenmai", "PFHxA"), ("Rosenmai", "PFOA"),
              ("Rosenmai", "PFNA"), ("Fragki", "PFOA"), ("Fragki", "PFNA")]
    xs = np.arange(len(series))
    axC.axhline(1.0, color=BLUE, lw=2.0, zorder=2)
    for i, (st, c) in enumerate(series):
        a = P.ASSAYS[st]; V = P.volumes(a); ds = a["doses"]
        m = [(d, P.measured_x_cell(st, c, d, V, a) * d) for d in ds
             if P.measured_x_cell(st, c, d, V, a)]
        q = [(d, P.qin_run(c, d, V)["C_cell"]) for d in ds]
        axC.plot(i, slope(m), "o", ms=8, mfc=INK, mec=INK, zorder=4)
        axC.plot(i, slope(q), "D", ms=7, mfc=ORANGE, mec=ORANGE, zorder=3)
    axC.set_xticks(xs)
    axC.set_xticklabels([s[0][:4] + "\n" + s[1] for s in series], fontsize=6.6)
    axC.set_ylabel("dose exponent  b   (C$_{cell}$ ~ C$_{nom}$$^b$)")
    axC.set_ylim(0.6, 2.05)
    axC.text(0.03, 0.95, "C", transform=axC.transAxes, fontsize=11, fontweight="bold")
    axC.text(0.97, 0.90, "Qin", transform=axC.transAxes, ha="right",
             color=ORANGE, fontsize=8, fontweight="bold")
    axC.text(0.97, 0.78, "measured", transform=axC.transAxes, ha="right",
             color=INK, fontsize=8, fontweight="bold")
    axC.text(0.97, 0.175, "Armitage: pinned at 1.00", transform=axC.transAxes,
             ha="right", color=BLUE, fontsize=7.5, fontweight="bold")

    # ---- D ----
    Vr = P.volumes(P.ASSAYS["Rosenmai"])
    mgp = Vr["Vcells"] * 1000 * P.mg_protein_per_mL_cell(P.ASSAYS["Rosenmai"])
    w = 0.34
    for i, c in enumerate(ORDER):
        f2a = P.ROSENMAI_FIG2A[c][100]
        good = f2a / MW[c] * P.mg_protein_per_mL_cell(P.ASSAYS["Rosenmai"])
        bad = good * MW[c]
        axD.bar(i - w / 2, good, w, color=AQUA, lw=0)
        axD.bar(i + w / 2, bad, w, color=ORANGE, lw=0)
        axD.text(i + w / 2, bad * 1.4, "%.0fx" % MW[c], ha="center",
                 fontsize=6.5, color=ORANGE, fontweight="bold")
    axD.set_yscale("log"); axD.set_ylim(5, 2e7)
    axD.axhline(1e6, color=INK2, lw=0.9, ls="--")
    axD.text(0.99, 0.80, "1 M — physically impossible", transform=axD.transAxes,
             ha="right", va="bottom", fontsize=6.5, color=INK2)
    axD.set_xticks(range(len(ORDER)))
    axD.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axD.set_ylabel("C$_{cell}$ at 100 µM  (µM)")
    axD.text(0.03, 0.95, "D", transform=axD.transAxes, fontsize=11, fontweight="bold")
    axD.text(0.11, 0.985, "■ Fig 2A as ng/mg (Abstract)", transform=axD.transAxes,
             color=AQUA, fontsize=7.2, fontweight="bold", va="top")
    axD.text(0.11, 0.905, "■ Fig 2A as nmol/mg (Results 3.3)",
             transform=axD.transAxes, color=ORANGE, fontsize=7.2,
             fontweight="bold", va="top")

    for ax in (axA, axB, axC, axD):
        ax.grid(True, which="major", color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.245, top=0.885, wspace=0.34)
    out = os.path.join(FIGS, "fig10_mbm_qin_vs_armitage.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)
    for i, ax in enumerate((axA, axB, axC, axD)):
        assert ax.get_position().height > 0.45, "panel %d collapsed" % i


if __name__ == "__main__":
    main()
