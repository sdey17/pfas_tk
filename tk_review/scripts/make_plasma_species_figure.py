#!/usr/bin/env python3
"""The human/trout plasma species difference, against what composition allows.

A: the measured human/trout plasma binding ratio against chain length, at a
   physiological and a bioassay concentration. The shaded band is everything
   a composition-only explanation can produce -- with one shared protein
   coefficient and one shared lipid coefficient, the ratio is confined
   between b_h/b_t and a_h/a_t by the two measured compositions alone, with
   no partition coefficient entering. Points outside it cannot be explained
   by composition at any choice of coefficients.
B: the shape, not the level. If both plasmas inherit albumin's saturation
   scaled by their protein fraction, every point lies on the solid line of
   slope f_prot,trout / f_prot,human = 0.76. Nine of thirteen sit on the
   x-axis instead: trout plasma does not detectably saturate where human
   plasma does.

Data: Qin et al. 2024, Environ Sci Technol 58(23):9954-9966,
doi:10.1021/acs.est.3c10824, Tables 1 and 3.

Run:  python3 tk_review/scripts/make_plasma_species_figure.py
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qin2024_plasma_species as S              # noqa: E402

OUT = os.path.join(HERE, "..", "figures")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID, BAND = "#e6e4e0", "#dfe8f4"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "axes.spines.top": False,
    "axes.spines.right": False,
})

HI = math.log10(S.HP_FPROT / S.TP_FPROT)              # +0.120
LO = math.log10((1 - S.HP_FPROT) / (1 - S.TP_FPROT))  # -0.517


def panelA(ax):
    carb = [n for n in S.PFAS if S.T1[n][1] == "carboxylate"]
    carb.sort(key=lambda n: S.T1[n][0])
    xs = [S.T1[n][0] for n in carb]
    ax.axhspan(LO, HI, color=BAND, zorder=0)
    ax.axhline(0.0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    for lab, cw, col, mk in (("1 nM  (physiological)", 1e-3, BLUE, "o"),
                             ("1 µM  (bioassay)", 1.0, ORANGE, "s")):
        ys = [S.logD(S.hp(n), cw) - S.logD(S.tp(n), cw) for n in carb]
        ax.plot(xs, ys, color=col, lw=2.0, zorder=3, label=lab,
                marker=mk, ms=8, mec=SURFACE, mew=1.6, clip_on=False)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{S.T1[n][0]}\n{n}" for n in carb], fontsize=7.6)
    ax.set_xlabel("perfluorinated carbons, and compound", labelpad=2)
    ax.set_ylabel("measured log$_{10}$ (D$_{human}$ / D$_{trout}$)")
    ax.set_ylim(-1.1, 2.8)
    yt = [-1, 0, 1, 2]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{10.0 ** v:g}×" if v else "1×"
                        for v in yt])
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)
    ax.text(3.12, (LO + HI) / 2, "everything a composition-only\n"
            "explanation can produce", color="#35567f", fontsize=7.4,
            va="center", ha="left", linespacing=1.25, zorder=4)
    ax.annotate("human binds harder", xy=(9.9, 0.17), fontsize=7.4,
                color=MUTED, ha="right", va="bottom")
    ax.annotate("trout binds harder", xy=(9.9, -0.78), fontsize=7.4,
                color=MUTED, ha="right", va="top")
    ax.legend(loc="upper right", fontsize=8, handletextpad=0.6)
    ax.set_title("A   the species ratio is not a species constant",
                 loc="left", fontsize=9.5, color=INK, pad=8)


def panelB(ax):
    pred = S.TP_FPROT / S.HP_FPROT
    lim = 0.78
    ax.plot([0, lim], [0, lim], color=MUTED, lw=0.8, ls=(0, (4, 3)),
            zorder=1)
    ax.plot([0, lim], [0, lim * pred], color=AQUA, lw=2.0, zorder=2)
    flat, curved = [], []
    for nm in S.T1:
        h = abs(S.hp(nm)[0] or 0.0)
        t = abs(S.tp(nm)[0] or 0.0)
        (flat if S.tp(nm)[0] is None else curved).append((h, t, nm))
    for pts, col, mk, lab in ((curved, BLUE, "o", "trout isotherm curved"),
                              (flat, ORANGE, "D", "trout isotherm flat")):
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ls="none",
                marker=mk, ms=8, color=col, mec=SURFACE, mew=1.6,
                zorder=3, label=lab, clip_on=False)
    # Label only the four points that carry information. The nine on the
    # x-axis are a group, not nine separate readings, so they get one
    # annotation between them rather than nine overlapping ones.
    for h, t, nm in curved:
        ax.annotate(nm, (h, t), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=7.2, color=INK2)
    ax.text(lim * 0.97, lim * pred * 0.985,
            "if composition\nexplained the shape", color="#157d57",
            fontsize=7.4, ha="right", va="top", linespacing=1.25)
    ax.annotate(f"{len(flat)} of {len(flat) + len(curved)} on this axis:\n"
                "no detectable saturation in trout\n"
                "plasma, at any human exponent",
                xy=(0.41, 0.006), xytext=(0.56, 0.045), fontsize=7.4,
                color=ORANGE, ha="center", va="bottom", linespacing=1.3,
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=1.0,
                                shrinkA=2, shrinkB=4))
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("|Freundlich exponent|, human plasma")
    ax.set_ylabel("|Freundlich exponent|, trout plasma")
    ax.grid(zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=8, handletextpad=0.6)
    ax.set_title("B   and the shape is not composition either",
                 loc="left", fontsize=9.5, color=INK, pad=8)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.4))
    panelA(axes[0])
    panelB(axes[1])
    fig.tight_layout(pad=1.6, w_pad=3.0)
    path = os.path.join(OUT, "fig23_plasma_species.png")
    fig.savefig(path, dpi=200)
    print(f"wrote {os.path.abspath(path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
