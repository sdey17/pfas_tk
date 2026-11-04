"""Albumin log K against chain length, coloured by how each number arrived.

Notebook script. It reads one CSV from the working directory and plots
it. No paths, no arguments: keep the CSV beside the notebook.

WHY CHAIN LENGTH RATHER THAN LOADING. The loading axis was chosen to
test whether assay loading explained the between-laboratory spread. It
does not, and it cannot be tested with this data, so an axis that
existed only for that test has no further job. Chain length is known
for every row including Qin's, it carries the strongest real signal in
the data, and it reads straight onto the question actually being asked:
for this compound, how far apart are the laboratories and where does
the value being used sit among them.

FACETED BY SPECIES AND HEAD GROUP, not pooled. A chain-length trend
only means anything within one protein and one head group: a sulfonate
binds more strongly than the carboxylate of the same length, and human,
bovine and rat albumin are different proteins. Pooling them would draw
a trend across three species at once.

Points only, no fitted lines. The trend is not a straight line -- it
rises with chain length and turns over around C7 to C11, which Li 2021
on HSA, Zhang 2013 on L-FABP, Allendorf on BSA and Fan on rat L-FABP
all show independently -- so a regression line would be the wrong
model, not merely an extra detail.
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV = "albumin_chainlength_data.csv"
PNG = "albumin_chainlength.png"

FAMILIES = ["carboxylate", "sulfonate"]
SPECIES = ["human", "bovine", "rat"]

# Categorical slots 1, 2 and 3 of the validated palette, same
# assignment as the loading plot so the two read as one set. Checked
# all-pairs on the light surface: every check passes, worst pair
# dE 9.2 deutan / 24.0 normal. Aqua carries a contrast WARN at 2.74:1,
# which the legend and the shipped CSV answer.
STYLE = {
    "measured":  ("#eb6834", "full", "measured directly"),
    "converted": ("#2a78d6", "full", "converted from a Kd"),
    "quenching": ("#1baf7a", "none", "converted on a quenching exponent"),
    "qin":       ("#52514e", "full", "Qin at this study's operating point"),
}
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 9, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
})


def load(path=CSV):
    rows = []
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if r["family"] not in FAMILIES or r["species"] not in SPECIES:
            continue
        # nCF 0 is a non-fluorinated reference ligand -- butyric acid,
        # octanoic acid, the fatty acids these assays use as
        # competitors. Real data, but not PFAS, so not on a PFAS
        # chain-length axis.
        if int(r["n_cf"]) < 1:
            continue
        r["n_cf"] = int(r["n_cf"])
        r["log_k"] = float(r["log_k"])
        rows.append(r)
    return rows


def dodge(marks, step=0.17):
    """Lanes for marks sharing an x, over SORTED marks so the result
    does not depend on the order the file happened to be read in.

    Sorted by SOURCE, deliberately not by log_k. Ordering lanes by y
    makes every vertical stack come out as a rising diagonal, which
    reads as a trend inside one chain length and is an artefact of the
    lane assignment.
    """
    out = {}
    for x in {m["n_cf"] for m in marks}:
        same = sorted((m for m in marks if m["n_cf"] == x),
                      key=lambda m: (m["source"], m["category"]))
        start = -step * (len(same) - 1) / 2.0
        for j, m in enumerate(same):
            out[id(m)] = start + j * step
    return out


def plot(rows=None, save=PNG):
    rows = load() if rows is None else rows
    cells = defaultdict(list)
    for r in rows:
        cells[(r["family"], r["species"])].append(r)

    xs = [r["n_cf"] for r in rows]
    lo, hi = min(xs) - 0.8, max(xs) + 0.8
    ylo = min(r["log_k"] for r in rows) - 0.5
    yhi = max(r["log_k"] for r in rows) + 0.5

    fig, axes = plt.subplots(len(FAMILIES), len(SPECIES),
                             figsize=(12.6, 7.2), sharex=True,
                             sharey=True)
    for i, fam in enumerate(FAMILIES):
        for j, sp in enumerate(SPECIES):
            ax = axes[i][j]
            here = cells.get((fam, sp), [])
            lanes = dodge(here)
            for m in sorted(here, key=lambda r: (r["n_cf"], r["log_k"])):
                colour, fill, _ = STYLE[m["category"]]
                ax.plot(m["n_cf"] + lanes[id(m)], m["log_k"], "o",
                        ms=7.0,
                        mfc=colour if fill == "full" else "none",
                        mec=SURFACE if fill == "full" else colour,
                        mew=1.3 if fill == "full" else 1.7,
                        zorder=6 if m["category"] != "qin" else 7)
            ax.set_title(f"{sp} · {fam}s   (n = {len(here)})",
                         loc="left", fontsize=10.5, fontweight="bold",
                         color=INK, pad=7)
            ax.grid(True, color=GRID, lw=0.7, zorder=1)
            ax.set_axisbelow(True)
            ax.set_xlim(lo, hi)
            ax.set_ylim(ylo, yhi)
            if i == len(FAMILIES) - 1:
                ax.set_xlabel("perfluorinated carbons", fontsize=9.5,
                              labelpad=6)
            if j == 0:
                ax.set_ylabel("log K  (L water / kg protein)",
                              fontsize=9.5)

    handles = [
        Line2D([], [], marker="o", ls="none", ms=8,
               mfc=col if fill == "full" else "none",
               mec=SURFACE if fill == "full" else col,
               mew=1.3 if fill == "full" else 1.7, label=lab)
        for col, fill, lab in STYLE.values()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4,
               frameon=False, fontsize=9.5, handletextpad=0.8,
               columnspacing=2.4, bbox_to_anchor=(0.5, -0.004))

    fig.suptitle("Albumin partition coefficient against chain length",
                 x=0.04, y=0.985, ha="left", fontsize=14.5,
                 fontweight="bold", color=INK)
    fig.text(0.04, 0.928,
             "Coloured by how each number was arrived at.",
             ha="left", fontsize=10, color=INK2)
    fig.text(0.04, 0.899,
             "Faceted by species and head group: a chain-length trend "
             "only holds within one protein and one head group.",
             ha="left", fontsize=10, color=INK2)
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.868))
    if save:
        fig.savefig(save, dpi=200)
    return fig


fig = plot()
