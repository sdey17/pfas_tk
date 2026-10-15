"""
Figure: Qin 2023 vs Armitage 2021 IV-MBM EQP v2.0 on the FRACTION-OF-DOSE axis.

Fraction of dose in cells is the axis to discriminate models on. Every source
agrees on protein per cell (Armitage AREc2 694 pg, project sheet 675 pg) and
disagrees 3.3x on cell VOLUME (10 pL vs 3.04 pL). All three datasets report per
mg protein, so the fraction is insensitive to that disagreement and C_cell
carries all of it.

  A  measured vs predicted fraction -- this repo's implementation of the
     published equations, and the project spreadsheet's, on one axis
  B  Rosenmai's internal contradiction, on the same unit-free axis: its Table 1
     uptake percentages against its Fig 2A values read both ways
  C  dose exponent b in C_cell ~ C_nom^b; Armitage is pinned at 1.00
  D  where each model puts the PFAS inside the cell

Run:  python3 tk_review/scripts/make_mbm_figure.py
"""
import csv
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
DB = os.path.join(os.path.dirname(HERE), "db")
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
PG_PROTEIN_PER_CELL = 675e-12 * 1e3        # mg; Armitage AREc2 gives 694 pg


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


def sheet_rows():
    """Project spreadsheet, nmol/mg protein -> fraction of nominal dose."""
    path = os.path.join(DB, "mbm_project_sheet.csv")
    out = []
    for r in csv.DictReader(open(path)):
        st, cp, d = r["study"], r["compound"], float(r["Ctest_uM"])
        a = P.ASSAYS[st]
        mg = a["cells"] * PG_PROTEIN_PER_CELL
        dose_nmol = d * a["v_uL"] * 1e-6 * 1000
        f = lambda v: float(v) * mg / dose_nmol
        out.append(dict(study=st, compound=cp, dose=d,
                        qin=f(r["Ccell_qin_nmolmg_sheet"]),
                        arm=f(r["Ccell_arm_nmolmg_sheet"])))
    return out


def main():
    fig, (axA, axB, axC, axD) = plt.subplots(1, 4, figsize=(13.8, 3.7))

    # ---------------- A: measured vs predicted fraction ----------------
    lo, hi = 2e-2, 3e1
    axA.plot([lo, hi], [lo, hi], "-", color=INK2, lw=1.0, zorder=2)
    for f in (10, 0.1):
        axA.plot([lo, hi], [lo * f, hi * f], ":", color=RULE, lw=0.9, zorder=1)
    meas = {}
    for st, a in P.ASSAYS.items():
        V = P.volumes(a)
        for c in ORDER:
            for d in a["doses"]:
                m = P.measured_x_cell(st, c, d, V, a)
                if m is None:
                    continue
                meas[(st, c, d)] = m
                op = Q.qin_params(c, NCF[c])["source"] != "measured"
                axA.plot(m * 100, P.qin_run(c, d, V)["x_cell"] * 100, MARK[st],
                         ms=5.5, mew=1.1, mfc="none" if op else ORANGE,
                         mec=ORANGE, zorder=4)
                axA.plot(m * 100, P.armitage_run(c, lkow(c), d, V)["x_cell"] * 100,
                         MARK[st], ms=5.5, mew=1.1, mfc="none" if op else BLUE,
                         mec=BLUE, zorder=4)
    for r in sheet_rows():
        m = meas.get((r["study"], r["compound"], r["dose"]))
        if m is None:
            continue
        axA.plot(m * 100, r["qin"] * 100, "x", ms=5, mew=1.1, color=ORANGE,
                 alpha=0.55, zorder=3)
        axA.plot(m * 100, r["arm"] * 100, "x", ms=5, mew=1.1, color=BLUE,
                 alpha=0.55, zorder=3)
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlim(lo, hi); axA.set_ylim(lo, hi)
    axA.set_xlabel("measured  % of dose in cells")
    axA.set_ylabel("predicted  % of dose in cells")
    axA.text(0.03, 0.96, "A", transform=axA.transAxes, fontsize=11, fontweight="bold")
    axA.text(0.97, 0.14, "Qin", transform=axA.transAxes, ha="right",
             color=ORANGE, fontsize=9, fontweight="bold")
    axA.text(0.97, 0.055, "Armitage", transform=axA.transAxes, ha="right",
             color=BLUE, fontsize=9, fontweight="bold")
    axA.text(0.985, 0.49,
             "● published eqns, measured params\n"
             "○ published eqns, extrapolated\n"
             "×  project spreadsheet",
             transform=axA.transAxes, ha="right", fontsize=6.3,
             color=MUTED, va="top", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.6, alpha=0.92))
    axA.text(0.985, 0.295, "o Rosenmai   □ Amstutz   △ Fragki",
             transform=axA.transAxes, ha="right", fontsize=6.4,
             color=MUTED, va="top", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.6, alpha=0.92))

    # ---------------- B: Rosenmai against itself ----------------
    mg = P.ASSAYS["Rosenmai"]["cells"] * PG_PROTEIN_PER_CELL
    Vb = 3e-3
    xs = np.arange(len(ORDER))
    w = 0.26
    for i, c in enumerate(ORDER):
        f2a = P.ROSENMAI_FIG2A[c][100]
        dose_nmol = 100 * Vb * 1000
        axB.bar(i - w, P.ROSENMAI_UPTAKE_PCT[c][100], w, color=INK2, lw=0)
        axB.bar(i, f2a / MW[c] * mg / dose_nmol * 100, w, color=AQUA, lw=0)
        axB.bar(i + w, f2a * mg / dose_nmol * 100, w, color=ORANGE, lw=0)
    axB.axhline(100, color=INK, lw=1.0, ls="--")
    axB.text(0.98, 0.845, "100% of dose", transform=axB.transAxes, ha="right",
             va="bottom", fontsize=6.5, color=INK)
    axB.set_yscale("log"); axB.set_ylim(5e-3, 1.2e3)
    axB.set_xticks(xs); axB.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axB.set_ylabel("% of dose in cells, at 100 µM")
    axB.text(0.03, 0.96, "B", transform=axB.transAxes, fontsize=11, fontweight="bold")
    for lab, col, yy in (("Table 1, uptake %", INK2, 0.145),
                         ("Fig 2A as ng/mg", AQUA, 0.085),
                         ("Fig 2A as nmol/mg", ORANGE, 0.025)):
        axB.text(0.035, yy, "■ " + lab, transform=axB.transAxes,
                 fontsize=6.6, color=col, fontweight="bold", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.6, alpha=0.92))
    axB.text(0.5, 0.955, "same measurement, 141–215× apart",
             transform=axB.transAxes, ha="center", fontsize=6.8, color=ORANGE,
             fontweight="bold")

    # ---------------- C: dose exponent ----------------
    series = [("Rosenmai", "PFBA"), ("Rosenmai", "PFHxA"), ("Rosenmai", "PFOA"),
              ("Rosenmai", "PFNA"), ("Fragki", "PFOA"), ("Fragki", "PFNA")]
    axC.axhline(1.0, color=BLUE, lw=2.0, zorder=2)
    for i, (st, c) in enumerate(series):
        a = P.ASSAYS[st]; V = P.volumes(a); ds = a["doses"]
        m = [(d, P.measured_x_cell(st, c, d, V, a) * d) for d in ds
             if P.measured_x_cell(st, c, d, V, a)]
        q = [(d, P.qin_run(c, d, V)["C_cell"]) for d in ds]
        axC.plot(i, slope(m), "o", ms=8, mfc=INK, mec=INK, zorder=4)
        axC.plot(i, slope(q), "D", ms=7, mfc=ORANGE, mec=ORANGE, zorder=3)
    axC.set_xticks(range(len(series)))
    axC.set_xticklabels([s[0][:4] + "\n" + s[1] for s in series], fontsize=6.6)
    axC.set_ylabel("dose exponent  b   (C$_{cell}$ ~ C$_{nom}$$^b$)")
    axC.set_ylim(0.6, 2.05)
    axC.text(0.03, 0.96, "C", transform=axC.transAxes, fontsize=11, fontweight="bold")
    axC.text(0.97, 0.90, "Qin", transform=axC.transAxes, ha="right",
             color=ORANGE, fontsize=8, fontweight="bold")
    axC.text(0.97, 0.78, "measured", transform=axC.transAxes, ha="right",
             color=INK, fontsize=8, fontweight="bold")
    axC.text(0.97, 0.175, "Armitage: pinned at 1.00", transform=axC.transAxes,
             ha="right", color=BLUE, fontsize=7.5, fontweight="bold")

    # ---------------- D: inside the cell ----------------
    V = P.volumes(P.ASSAYS["Rosenmai"]); comp = V["comp"]
    w = 0.38
    for i, c in enumerate(ORDER):
        a = P.armitage_D(lkow(c), comp)
        ml = comp["ml"] * a["Dmw"]; wa = comp["water"]
        ot = comp["sl"] * a["Dow"] + P.P_NLOM * comp["prot"] * a["Dow"]
        t = ml + wa + ot
        b = 0.0
        for frac, col in ((ml / t, BLUE), (wa / t, AQUA), (ot / t, MUTED)):
            axD.bar(i - w / 2, frac * 100, w, bottom=b * 100, color=col, lw=0); b += frac
        Db, Ds, Dl, p = P.qin_D(c, 1.0)
        qp = Ds * comp["prot"]; ql = Dl * (comp["sl"] + comp["ml"]); qt = qp + ql
        b = 0.0
        for frac, col in ((qp / qt, ORANGE), (ql / qt, PURPLE)):
            axD.bar(i + w / 2, frac * 100, w, bottom=b * 100, color=col, lw=0); b += frac
    axD.axhline(100, color=RULE, lw=0.8)
    axD.set_xticks(range(len(ORDER)))
    axD.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    axD.set_ylabel("share of cell sorptive capacity (%)")
    axD.set_ylim(0, 168)
    axD.text(0.03, 0.74, "D", transform=axD.transAxes, fontsize=11, fontweight="bold")
    axD.text(0.5, 1.045, "left bar Armitage      right bar Qin",
             transform=axD.transAxes, ha="center", fontsize=7, color=INK2)
    for lab, col, xx, yy in (("membrane lipid", BLUE, 0.03, 0.965),
                             ("cell water", AQUA, 0.03, 0.885),
                             ("protein (D$_{non-spec}$)", ORANGE, 0.47, 0.965),
                             ("lipid", PURPLE, 0.47, 0.885)):
        axD.text(xx, yy, "■ " + lab, transform=axD.transAxes,
                 fontsize=6.6, color=col, fontweight="bold", va="top")

    for ax in (axA, axB, axC, axD):
        ax.grid(True, which="major", color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.25, top=0.885, wspace=0.345)
    out = os.path.join(FIGS, "fig10_mbm_qin_vs_armitage.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)
    for i, ax in enumerate((axA, axB, axC, axD)):
        assert ax.get_position().height > 0.45, "panel %d collapsed" % i


if __name__ == "__main__":
    main()
