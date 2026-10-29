#!/usr/bin/env python3
"""Species differences in PFAS-albumin binding, against the noise floor.

A: Starnes et al. 2024, four species in one differential-scanning-
   fluorimetry run, seven compounds. Each compound's Kd is plotted
   RELATIVE TO ITS OWN GEOMETRIC MEAN across the four species, because
   the question is how far species moves a constant, not which compound
   binds tightest. Absolute Kd spans 0.33-3.06 mM and is in the label.
B: Fischer et al. 2024, human against bovine albumin, 27 compounds, one
   method, one laboratory. The band is +/-0.59 log, the median
   disagreement between two DIRECT measurements of the same compound,
   protein and species in different laboratories.
C: why A and B are small. Every scale is from this database; the two
   species effects are the coloured bars.

The point of the figure is C. Changing the mammal moves a PFAS-albumin
constant less than changing the laboratory does.

Run:  python3 tk_review/scripts/make_species_figure.py
"""
import collections
import csv
import importlib.util
import math
import os
import re
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
OUT = os.path.join(os.path.dirname(HERE), "figures", "species_differences.png")

_ns = importlib.util.spec_from_file_location(
    "pfas_names", os.path.join(HERE, "pfas_names.py"))
N = importlib.util.module_from_spec(_ns)
_ns.loader.exec_module(N)

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

SPECIES = ["porcine", "rat", "human", "bovine"]
PFCA = {"PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
        "PFUnDA", "PFDoDA", "PFTrDA", "PFTeDA"}
PFSA = {"PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS"}
NUM = re.compile(r"^\s*([\d.]+)\s*$")

# From scripts/logk_converted_vs_measured.py and
# scripts/binding_species_sex_resolution.py.
SCALES = [
    ("human vs bovine albumin,\none laboratory, 27 compounds", 0.09, BLUE),
    ("one laboratory,\ntwo albumin batches", 0.19, MUTED),
    ("four mammals,\none assay, 7 compounds", 0.36, BLUE),
    ("two laboratories,\nsame measured log K", 0.59, MUTED),
    ("two laboratories,\nthe same albumin Kd", 1.45, MUTED),
]


def head(c):
    return ("carboxylate" if c in PFCA else
            "sulfonate" if c in PFSA else "other")


COLOUR = {"carboxylate": BLUE, "sulfonate": ORANGE, "other": AQUA}


def load():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    ncf = {}
    for r in rows:
        v = r["n_perfluorinated_carbons"].strip()
        if v.isdigit():
            ncf.setdefault(N.canon(r["chemical"]), int(v))

    starnes = collections.defaultdict(dict)
    for r in rows:
        if r["study"].startswith("Starnes") and r["year"] == "2024":
            m = NUM.match(r["kd"])
            if m:
                starnes[N.canon(r["chemical"])][r["protein_species"]] = \
                    float(m.group(1))
    starnes = {c: v for c, v in starnes.items() if len(v) == 4}

    fis = collections.defaultdict(dict)
    for r in rows:
        if (r["study"] == "Fischer et al." and r["year"] == "2024"
                and r["ka"].strip()
                and "serum albumin" in r["protein"].lower()):
            try:
                fis[N.canon(r["chemical"])][r["protein_species"]] = \
                    float(r["ka"])
            except ValueError:
                pass
    fis = {c: v for c, v in fis.items() if len(v) == 2}
    return ncf, starnes, fis


def panelA(ax, ncf, st):
    xs = range(len(SPECIES))
    ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
    # Plotted as log10(mean / Kd), so UP IS TIGHTER BINDING. Plotting
    # the Kd ratio instead puts tighter binding at the bottom, which
    # reads backwards against every axis label one would write.
    ends = []
    for c in sorted(st, key=lambda c: ncf.get(c, 99)):
        v = st[c]
        gm = math.exp(statistics.fmean(math.log(v[s]) for s in SPECIES))
        ys = [math.log10(gm / v[s]) for s in SPECIES]
        col = COLOUR[head(c)]
        ax.plot(xs, ys, "-", color=col, lw=1.4, alpha=0.9, zorder=3)
        ax.plot(xs, ys, "o", ms=5, mfc=col, mec=SURFACE, mew=1.0, zorder=4)
        ends.append([ys[3], c, col])
    # Labels at the bovine end collide where two compounds land within
    # about 0.02 log of each other; nudge them apart in plot order.
    ends.sort(key=lambda e: -e[0])
    for i in range(1, len(ends)):
        gap = ends[i - 1][0] - ends[i][0]
        if gap < 0.030:
            ends[i][0] = ends[i - 1][0] - 0.030
    for y, c, col in ends:
        ax.annotate(c, xy=(3.06, y), fontsize=7.2, color=col,
                    va="center", ha="left")
    ax.set_xticks(list(xs))
    ax.set_xticklabels([s.capitalize() for s in SPECIES])
    ax.set_xlim(-0.25, 4.15)
    ax.set_ylabel("log$_{10}$ affinity relative to the compound's own\n"
                  "mean across the four species   (up = tighter)")
    ax.set_title("A   Four species, one assay (Starnes 2024)",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    ax.grid(True, axis="y", lw=0.6)
    ax.set_axisbelow(True)
    lo, hi = ax.get_ylim()
    ax.text(-0.18, hi - 0.01, "tighter binding", fontsize=6.9,
            color=MUTED, va="top", ha="left")
    ax.text(-0.18, lo + 0.01, "weaker binding", fontsize=6.9,
            color=MUTED, va="bottom", ha="left")


def panelB(ax, ncf, fis):
    lo, hi = 1.5, 6.0
    ax.fill_between([lo, hi], [lo - 0.59, hi - 0.59], [lo + 0.59, hi + 0.59],
                    color=BAND, zorder=0, lw=0)
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1.2, zorder=1)
    for c, v in fis.items():
        ax.plot(v["bovine"], v["human"], "o", ms=6,
                mfc=COLOUR[head(c)], mec=SURFACE, mew=1.1, zorder=4)
    d = [v["human"] - v["bovine"] for v in fis.values()]
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("bovine serum albumin,  log D")
    ax.set_ylabel("human serum albumin,  log D")
    ax.set_title(f"B   Human against bovine albumin, "
                 f"{len(fis)} compounds (Fischer 2024)",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    ax.grid(True, lw=0.6)
    ax.set_axisbelow(True)
    ax.text(lo + 0.12, hi - 0.12,
            f"median difference {statistics.median(d):+.2f} log\n"
            f"human higher in {sum(1 for x in d if x > 0)} of {len(d)}",
            fontsize=7.6, color=INK2, va="top", linespacing=1.4)
    ax.text(lo + 0.12, hi - 0.60,
            "shaded: ±0.59 log, what two laboratories\n"
            "differ by on the same measurement",
            fontsize=7, color="#4a6d9e", ha="left", va="top",
            linespacing=1.35)
    handles = [Line2D([], [], marker="o", ls="", ms=6, mfc=COLOUR[k],
                      mec=SURFACE, mew=1.0, label=k)
               for k in ("carboxylate", "sulfonate", "other")]
    ax.legend(handles=handles, loc="lower right", fontsize=7.4,
              labelspacing=0.3, handletextpad=0.4, borderpad=0.2)


def panelC(ax):
    ys = range(len(SCALES))
    for y, (lab, v, col) in zip(ys, SCALES):
        ax.barh(y, v, height=0.52, color=col, zorder=3)
        ax.text(v + 0.02, y, f"{v:.2f} log   ({10 ** v:.1f}×)",
                va="center", fontsize=7.8, color=INK)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([s[0] for s in SCALES], fontsize=7.4)
    ax.set_xlim(0, 1.85)
    ax.set_xlabel("typical disagreement   (log units)")
    ax.set_title("C   Species effects (blue) against laboratory noise "
                 "(grey), same database",
                 loc="left", fontsize=9.5, color=INK, pad=8)
    ax.grid(True, axis="x", lw=0.6)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.annotate("", xy=(0.36, 2.52), xytext=(0.59, 2.52),
                arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.9))
    ax.text(0.475, 2.52, "species < laboratory",
            bbox=dict(fc=SURFACE, ec="none", pad=1.5),
            ha="center", va="center", fontsize=7.4, color=INK2)


def main():
    ncf, st, fis = load()
    fig = plt.figure(figsize=(11.6, 8.9))
    top = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.92], wspace=0.30,
                           left=0.075, right=0.975, top=0.835, bottom=0.42)
    bot = fig.add_gridspec(1, 1, left=0.255, right=0.90,
                           top=0.315, bottom=0.075)
    panelA(fig.add_subplot(top[0, 0]), ncf, st)
    panelB(fig.add_subplot(top[0, 1]), ncf, fis)
    panelC(fig.add_subplot(bot[0, 0]))

    fig.text(0.075, 0.962,
             "Changing the mammal moves a PFAS-albumin constant less "
             "than changing the laboratory does",
             ha="left", va="top", fontsize=13, color=INK)
    fig.text(0.075, 0.925,
             "Both species comparisons are within a single laboratory, "
             "which is the only kind worth running: two laboratories "
             "differ by 1.45 log on the\nsame albumin Kd. Panel A is "
             "differential scanning fluorimetry — the ordering and "
             "spread are sound, the absolute Kd values are not "
             "physiological affinities.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.5)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f"wrote {os.path.abspath(OUT)}")
    print(f"  panel A {len(st)} compounds x 4 species; "
          f"panel B {len(fis)} compounds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
