"""One panel per PFAS: every albumin value, by who produced it.

Notebook script. It reads one CSV from the working directory and plots
it. No paths, no arguments: keep the CSV beside the notebook.

WHAT THIS ANSWERS that the pooled plots do not. With compounds on one
axis you can see the shape of the literature but you cannot follow a
single PFAS through it. Here each compound gets its own panel and each
laboratory its own row, on a y axis shared across panels, so a row
reads straight across: a source that sits low for one compound can be
checked against where it sits for the next.

The rows are grouped by species, because human, bovine and rat albumin
are different proteins and a value should never be read across them by
accident. Grouping keeps them adjacent and labelled rather than
pooling them into one cloud.

Qin is both a row and a dashed line, since it is the value being
compared against: the line is its number for that compound, so the
distance from any dot to the line is what the comparison is worth.
The dot and the line are the SAME NUMBER -- the dot keeps Qin in the
stack as one laboratory among many, the line lets every other row's
distance from it be read without tracing down the panel.

QIN'S VALUES ARE ALL BOVINE, and the line is drawn across the human
and rat rows too. That is the one place this figure crosses species,
so the legend says so. Grouping the rows by species keeps every other
comparison within one protein; the line is the exception and is
labelled as one.

EDIT `COMPOUNDS` to change which panels appear. Anything in the CSV
works; the default five are the ones under active comparison. The
layout adapts, but more than about six panels gets cramped.
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV = "albumin_chainlength_data.csv"
PNG = "albumin_by_compound.png"

COMPOUNDS = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]
SPECIES_ORDER = ["human", "bovine", "rat"]

# Categorical slots 1, 2 and 3 of the validated palette, the same
# assignment as the other two plots so all three read as one set.
# Checked all-pairs on the light surface: every check passes, worst
# pair dE 9.2 deutan / 24.0 normal. Aqua carries a contrast WARN at
# 2.74:1, which the legend and the shipped CSV answer.
STYLE = {
    "measured":  ("#eb6834", "full", "measured directly"),
    "converted": ("#2a78d6", "full", "converted from a Kd"),
    "quenching": ("#1baf7a", "none", "converted on a quenching exponent"),
    "qin":       ("#52514e", "full",
                  "Qin 2024 (BOVINE), at this study's operating point"),
}
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, BAND = "#e1e0d9", "#c3c2b7", "#f4f3ee"

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
        if r["pfas"] not in COMPOUNDS:
            continue
        r["log_k"] = float(r["log_k"])
        r["n_cf"] = int(r["n_cf"])
        rows.append(r)
    return rows


def plot(rows=None, save=PNG):
    rows = load() if rows is None else rows

    # One y slot per source-species pair, grouped by species so the
    # blocks are contiguous and a value is never read across species
    # by accident.
    keys = sorted({(r["species"], r["source"]) for r in rows},
                  key=lambda k: (SPECIES_ORDER.index(k[0])
                                 if k[0] in SPECIES_ORDER else 9, k[1]))
    slot = {k: i for i, k in enumerate(reversed(keys))}
    labels = [f"{s}  ·  {sp}" for sp, s in reversed(keys)]

    order = sorted({r["pfas"] for r in rows},
                   key=lambda c: next(r["n_cf"] for r in rows
                                      if r["pfas"] == c))
    ncf = {r["pfas"]: r["n_cf"] for r in rows}
    cells = defaultdict(list)
    for r in rows:
        cells[r["pfas"]].append(r)

    lo = min(r["log_k"] for r in rows) - 0.5
    hi = max(r["log_k"] for r in rows) + 0.5

    ncol = 3
    nrow = -(-(len(order) + 1) // ncol)
    fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                             figsize=(13.4, 1.0 + 0.30 * len(keys) * nrow))
    flat = axes.ravel() if hasattr(axes, "ravel") else [axes]

    for i, c in enumerate(order):
        ax = flat[i]
        here = cells[c]
        # Banded rows: a shared y axis is only readable if the eye can
        # follow one row across a wide panel.
        for y in range(len(keys)):
            if y % 2 == 0:
                ax.axhspan(y - 0.5, y + 0.5, color=BAND, zorder=0)

        qin = next((r["log_k"] for r in here if r["category"] == "qin"),
                   None)
        if qin is not None:
            ax.axvline(qin, color=MUTED, lw=1.1, ls=(0, (6, 3)), zorder=2)

        rowed = defaultdict(list)
        for r in here:
            rowed[(r["species"], r["source"])].append(r)
        for k, group in rowed.items():
            y0 = slot[k]
            group = sorted(group, key=lambda r: r["log_k"])
            step = 0.17
            start = -step * (len(group) - 1) / 2.0
            for j, r in enumerate(group):
                colour, fill, _ = STYLE[r["category"]]
                ax.plot(r["log_k"], y0 + start + j * step, "o", ms=6.8,
                        mfc=colour if fill == "full" else "none",
                        mec=SURFACE if fill == "full" else colour,
                        mew=1.3 if fill == "full" else 1.7, zorder=6)

        ax.set_title(f"{c}   (C{ncf[c]})", loc="left", fontsize=11.5,
                     fontweight="bold", color=INK, pad=7)
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.7, len(keys) - 0.3)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels(labels, fontsize=8.4)
        ax.tick_params(axis="y", length=0)
        ax.grid(True, axis="x", color=GRID, lw=0.7, zorder=1)
        ax.set_axisbelow(True)
        if i >= len(order) - ncol or len(order) - i <= ncol:
            ax.set_xlabel("log K  (L water / kg protein)", fontsize=9.5,
                          labelpad=6)
            ax.tick_params(labelbottom=True)

    for j in range(len(order), len(flat)):
        flat[j].axis("off")

    handles = [
        Line2D([], [], marker="o", ls="none", ms=8,
               mfc=col if fill == "full" else "none",
               mec=SURFACE if fill == "full" else col,
               mew=1.3 if fill == "full" else 1.7, label=lab)
        for col, fill, lab in STYLE.values()
    ] + [Line2D([], [], color=MUTED, lw=1.1, ls=(0, (6, 3)),
                label="the same value, drawn across every row —\n"
                      "reading a HUMAN or RAT row against it\n"
                      "crosses species")]
    flat[len(order)].legend(handles=handles, loc="upper left",
                            frameon=False, fontsize=9.5,
                            handletextpad=0.9, labelspacing=1.0,
                            bbox_to_anchor=(0.02, 0.92))

    fig.suptitle("Every albumin value, one panel per PFAS",
                 x=0.012, y=0.995, ha="left", fontsize=14.5,
                 fontweight="bold", color=INK)
    fig.text(0.012, 0.967,
             "Rows are shared across panels, so a laboratory can be "
             "followed from one compound to the next. Grouped by "
             "species; coloured by how each number was arrived at.",
             ha="left", fontsize=9.6, color=INK2)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.945))
    if save:
        fig.savefig(save, dpi=200)
    return fig


fig = plot()
