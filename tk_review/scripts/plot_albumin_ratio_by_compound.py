"""One panel per PFAS: the loading each laboratory worked at.

Notebook script. It reads one CSV from the working directory and plots
it. No paths, no arguments: keep the CSV beside the notebook.

The companion to plot_albumin_by_compound.py, which puts log K on the
same rows. Read the two side by side and a row answers both halves of
the question at once: what this laboratory got, and how hard it was
loading the protein when it got it.

MARKER SHAPE CARRIES THE PROVENANCE OF THE X VALUE, which on this plot
is not a detail. Only Allendorf and Fischer state a ligand:protein
ratio -- Allendorf per compound in SI Table SI-2.3, Fischer as a
design bound. Every other row states two concentrations and no ratio,
so its x is this project dividing one by the other. A circle is the
paper's number; a square is arithmetic done here. On the five
compounds under comparison that split is about 18 points to 78.

Qin has no row here. Its value is a Freundlich fit evaluated at a free
concentration this study's operating point implies, not a measurement
at a stated loading, so it has no honest x. Its absence from this
plot, against its presence on the log K one, is the point.

EDIT `COMPOUNDS` to change which panels appear.
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV = "albumin_chainlength_data.csv"
PNG = "albumin_ratio_by_compound.png"

COMPOUNDS = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]
SPECIES_ORDER = ["human", "bovine", "rat"]

# Categorical slots 1, 2 and 3 of the validated palette, the same
# assignment as the other plots so they read as one set. Checked
# all-pairs on the light surface: every check passes, worst pair
# dE 9.2 deutan / 24.0 normal. Aqua carries a contrast WARN at 2.74:1,
# which the legend and the shipped CSV answer.
STYLE = {
    "measured":  ("#eb6834", "full", "measured directly"),
    "converted": ("#2a78d6", "full", "converted from a Kd"),
    "quenching": ("#1baf7a", "none", "converted on a quenching exponent"),
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
        if not str(r["log_ratio"]).strip():
            continue           # no loading reported or derivable
        r["log_ratio"] = float(r["log_ratio"])
        r["n_cf"] = int(r["n_cf"])
        rows.append(r)
    return rows


def plot(rows=None, save=PNG):
    rows = load() if rows is None else rows

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

    lo = min(r["log_ratio"] for r in rows) - 0.6
    hi = max(r["log_ratio"] for r in rows) + 0.6

    ncol = 3
    nrow = -(-(len(order) + 1) // ncol)
    # HEAD is the title strip in INCHES, not in figure fractions. With
    # fractions the heading lands on the subtitle as soon as the panel
    # count, and with it the figure height, changes.
    head = 0.95
    h = head + 0.30 * len(keys) * nrow
    fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                             figsize=(13.4, h))
    flat = axes.ravel()

    for i, c in enumerate(order):
        ax = flat[i]
        here = cells[c]
        for y in range(len(keys)):
            if y % 2 == 0:
                ax.axhspan(y - 0.5, y + 0.5, color=BAND, zorder=0)
        # Unit loading: one ligand per protein, where occupancy stops
        # being a rounding error.
        ax.axvline(0.0, color=MUTED, lw=1.0, ls=(0, (2, 3)), zorder=2)

        rowed = defaultdict(list)
        for r in here:
            rowed[(r["species"], r["source"])].append(r)
        for k, group in rowed.items():
            y0 = slot[k]
            group = sorted(group, key=lambda r: r["log_ratio"])
            # Keep a stack inside its own band. Gao 2025 reports seven
            # PFHpA values at almost one loading, and at a fixed step
            # that column climbs out of the row and into its neighbour.
            n = len(group)
            step = min(0.17, 0.78 / (n - 1)) if n > 1 else 0.0
            shrink = 1.0 if n <= 4 else 0.82
            start = -step * (n - 1) / 2.0
            for j, r in enumerate(group):
                colour, fill, _ = STYLE[r["category"]]
                stated = r["ratio_reported"] == "yes"
                ax.plot(r["log_ratio"], y0 + start + j * step,
                        "o" if stated else "s",
                        ms=(6.8 if stated else 6.2) * shrink,
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
        if len(order) - i <= ncol:
            ax.set_xlabel("log₁₀  ligand : protein ratio",
                          fontsize=9.5, labelpad=6)
            ax.tick_params(labelbottom=True)

    for j in range(len(order), len(flat)):
        flat[j].axis("off")

    handles = [
        Line2D([], [], marker="o", ls="none", ms=8,
               mfc=col if fill == "full" else "none",
               mec=SURFACE if fill == "full" else col,
               mew=1.3 if fill == "full" else 1.7, label=lab)
        for col, fill, lab in STYLE.values()
    ] + [
        Line2D([], [], marker="o", ls="none", ms=8, mfc=MUTED,
               mec=SURFACE, mew=1.3, label="ratio stated by the paper"),
        Line2D([], [], marker="s", ls="none", ms=7.4, mfc=MUTED,
               mec=SURFACE, mew=1.3,
               label="ratio computed here, PFAS ÷ protein"),
        Line2D([], [], color=MUTED, lw=1.0, ls=(0, (2, 3)),
               label="one ligand per protein"),
    ]
    flat[len(order)].legend(handles=handles, loc="upper left",
                            frameon=False, fontsize=9.5,
                            handletextpad=0.9, labelspacing=1.0,
                            bbox_to_anchor=(0.02, 0.95))

    fig.suptitle("The loading each laboratory worked at, one panel "
                 "per PFAS", x=0.012, y=1 - 0.18 / h, ha="left",
                 va="top", fontsize=14.5, fontweight="bold", color=INK)
    fig.text(0.012, 1 - 0.50 / h,
             "Same rows as the log K panel, so the two read together. "
             "Shape carries where the x value came from: most of this "
             "axis is arithmetic, not a reported number.",
             ha="left", va="top", fontsize=9.6, color=INK2)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1 - head / h))
    if save:
        fig.savefig(save, dpi=200)
    return fig


fig = plot()
