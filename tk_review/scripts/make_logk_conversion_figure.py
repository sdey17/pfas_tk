#!/usr/bin/env python3
"""Does a Kd-derived protein-water coefficient agree with a measured one?

A: the 33 pairs, converted against measured, with the 1:1 line and a band
   of +/-0.59 log around it -- the median disagreement between two DIRECT
   measurements of the same compound, protein and species in different
   laboratories. A point inside that band is as close as two measurements
   of this quantity normally get.
B: the same residuals by source, against the two baselines. Chen 2025 is
   the only well-powered source and the only one that assumes n = 1.
C: the question panel A cannot answer on its own -- is 0.80 log big? Four
   error scales on one axis, from the same database.

Colour carries where the stoichiometry came from, three slots, validated
all-pairs. Shape carries the head group. Wu et al. 2024 is drawn but
clipped with its value printed, because at -4.18 log it would otherwise
set the scale for everything else.

Data: db/logk_converted_vs_measured.csv, written by
      scripts/logk_converted_vs_measured.py.

Run:  python3 tk_review/scripts/make_logk_conversion_figure.py
"""
import collections
import csv
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "logk_converted_vs_measured.csv")
OUT = os.path.join(HERE, "..", "figures", "logk_converted_vs_measured.png")

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID, BAND, BAND2 = "#e6e4e0", "#dfe8f4", "#eef3fa"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "axes.spines.top": False,
    "axes.spines.right": False,
})

COLOUR = {"assumed": BLUE, "double-log": ORANGE, "isotherm": AQUA}
LABEL = {"assumed": "n assumed = 1", "double-log": "n from a double-log slope",
         "isotherm": "n fitted from an isotherm"}
SHAPE = {"carboxylate": "o", "sulfonate": "s", "other": "^"}

# From scripts/logk_converted_vs_measured.py, parts 2, 3 and 6.
WITHIN_LAB = 0.19
BETWEEN_LAB = 0.59
LOO_REGRESSION = 0.44


def load():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    for r in rows:
        r["converted_logK"] = float(r["converted_logK"])
        r["measured_logK"] = float(r["measured_logK"])
        r["residual"] = float(r["residual"])
        r["nCF"] = int(r["nCF"]) if r["nCF"] else 0
    return rows


def panelA(ax, rows):
    lo, hi = 1.9, 5.9
    ax.fill_between([lo, hi], [lo - BETWEEN_LAB, hi - BETWEEN_LAB],
                    [lo + BETWEEN_LAB, hi + BETWEEN_LAB],
                    color=BAND, zorder=0, lw=0)
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1.2, zorder=1)
    ax.text(hi - 0.12, hi - 0.30, "1:1", color=MUTED, fontsize=8,
            ha="right", va="top")
    ax.text(lo + 0.08, lo + BETWEEN_LAB - 0.05,
            f"±{BETWEEN_LAB:.2f} log: what two direct\nmeasurements "
            f"normally differ by",
            color="#4a6d9e", fontsize=7.2, ha="left", va="top", linespacing=1.3)

    clipped = []
    for r in rows:
        y = r["converted_logK"]
        if y < lo:
            clipped.append(r)
            continue
        ax.plot(r["measured_logK"], y, SHAPE.get(r["head"], "^"),
                ms=6.5, mfc=COLOUR[r["n_how"]], mec=SURFACE, mew=1.2,
                zorder=4, alpha=0.95)
    for r in clipped:
        ax.plot(r["measured_logK"], lo + 0.06, "v", ms=6.5,
                mfc=COLOUR[r["n_how"]], mec=SURFACE, mew=1.2, zorder=4)
    if clipped:
        xs = sorted(r["measured_logK"] for r in clipped)
        ax.annotate(f"Wu 2024: {len(clipped)} points off the bottom,\n"
                    f"log K {min(r['converted_logK'] for r in clipped):.1f} "
                    f"to {max(r['converted_logK'] for r in clipped):.1f}",
                    xy=(xs[-2], lo + 0.09),
                    xytext=(xs[-2] - 0.05, lo + 0.40),
                    fontsize=7.2, color=ORANGE, ha="center", va="bottom",
                    linespacing=1.35,
                    arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.9,
                                    shrinkA=3, shrinkB=3))

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("measured log K   (Fischer 2024, L$_{water}$/kg protein)")
    ax.set_ylabel("converted from Kd and n   (log K)")
    ax.set_title("A   Converted against measured, 33 pairs",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    ax.grid(True, lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    handles = [Line2D([], [], marker="o", ls="", ms=6.5, mfc=COLOUR[k],
                      mec=SURFACE, mew=1.2, label=LABEL[k])
               for k in ("assumed", "double-log", "isotherm")]
    handles += [Line2D([], [], marker=SHAPE[h], ls="", ms=6, mfc=MUTED,
                       mec=SURFACE, mew=1.0, label=h)
                for h in ("carboxylate", "sulfonate", "other")]
    ax.legend(handles=handles, loc="upper left", fontsize=7.4,
              labelspacing=0.35, handletextpad=0.4, borderpad=0.2)


def panelB(ax, rows):
    order = sorted({r["converted_source"] for r in rows},
                   key=lambda s: -len([r for r in rows
                                       if r["converted_source"] == s]))
    ax.axhspan(-WITHIN_LAB, WITHIN_LAB, color=BAND, zorder=0, lw=0)
    ax.axhspan(-BETWEEN_LAB, BETWEEN_LAB, color=BAND2, zorder=0, lw=0)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=1)
    lo = -2.4
    for i, src in enumerate(order):
        sub = [r for r in rows if r["converted_source"] == src]
        col = COLOUR[sub[0]["n_how"]]
        for j, r in enumerate(sub):
            x = i + (j - (len(sub) - 1) / 2) * min(0.055, 0.5 / len(sub))
            y = max(r["residual"], lo + 0.08)
            ax.plot(x, y, "v" if r["residual"] < lo else
                    SHAPE.get(r["head"], "^"),
                    ms=5.5, mfc=col, mec=SURFACE, mew=0.9, zorder=4,
                    alpha=0.95)
        med = statistics.median([r["residual"] for r in sub])
        ax.plot([i - 0.3, i + 0.3], [max(med, lo + 0.08)] * 2, color=INK,
                lw=2.0, zorder=5, solid_capstyle="butt")
        ax.text(i, lo + 0.30 if med < lo else med + 0.14,
                f"{med:+.2f}", ha="center", fontsize=7.4, color=INK,
                zorder=6)
    ax.set_xlim(-0.95, len(order) - 0.45)
    ax.text(-0.90, (BETWEEN_LAB + WITHIN_LAB) / 2 + 0.03, "between\nlabs",
            fontsize=6.9, color="#4a6d9e", va="center", ha="left",
            linespacing=1.25)
    ax.text(-0.90, 0.0, "within\none lab", fontsize=6.9, color="#2f5e9e",
            va="center", ha="left", linespacing=1.25)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{s}\nn={len([r for r in rows if r['converted_source'] == s])}"
                        for s in order], fontsize=7.6)
    ax.set_ylim(lo, 1.3)
    ax.set_ylabel("converted − measured   (log units)")
    ax.set_title("B   The residual by source, against the two baselines",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    ax.grid(True, axis="y", lw=0.6, zorder=0)
    ax.set_axisbelow(True)


def panelC(ax, rows):
    chen = abs(statistics.median([r["residual"] for r in rows
                                  if r["converted_source"] == "Chen 2025"]))
    bars = [
        ("one lab, two albumin batches", WITHIN_LAB, MUTED),
        ("chain-length line, leave-one-out", LOO_REGRESSION, AQUA),
        ("two labs, same measurement", BETWEEN_LAB, MUTED),
        ("Kd → K, n assumed = 1", chen, BLUE),
    ]
    ys = range(len(bars))
    for y, (lab, v, col) in zip(ys, bars):
        ax.barh(y, v, height=0.52, color=col, zorder=3)
        ax.text(v + 0.02, y, f"{v:.2f} log  ({10 ** v:.1f}×)", va="center",
                fontsize=7.8, color=INK)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([b[0] for b in bars], fontsize=7.6)
    ax.set_xlim(0, 1.18)
    ax.set_xlabel("typical disagreement   (log units)")
    ax.set_title("C   Is 0.80 log large? Four error scales from this database",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    ax.grid(True, axis="x", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def main():
    rows = load()
    fig = plt.figure(figsize=(11.6, 8.8))
    top = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.08], wspace=0.24,
                           left=0.065, right=0.975, top=0.845, bottom=0.40)
    bot = fig.add_gridspec(1, 1, left=0.215, right=0.86,
                           top=0.285, bottom=0.075)
    panelA(fig.add_subplot(top[0, 0]), rows)
    panelB(fig.add_subplot(top[0, 1]), rows)
    panelC(fig.add_subplot(bot[0, 0]), rows)

    fig.text(0.065, 0.962,
             "A dissociation constant converts to a protein-water "
             "coefficient; the two do not agree",
             ha="left", va="top", fontsize=13, color=INK)
    fig.text(0.065, 0.925,
             "33 pairs, same compound, protein and species. The gap is a "
             "missing site count, not non-specific binding: at the trace "
             "limit\nQin 2024's non-specific term is 1\u201314% of the total "
             "\u2014 0.07 log, where the gap is 0.80.",
             ha="left", va="top", fontsize=8.6, color=INK2, linespacing=1.5)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f"wrote {os.path.abspath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
