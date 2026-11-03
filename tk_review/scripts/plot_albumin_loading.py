"""Albumin log K against loading, coloured by how each number arrived.

Notebook script. It reads one CSV from the working directory and plots
it. No paths, no arguments: keep the CSV beside the notebook.

Points only. There are no fitted lines here on purpose. A line through
these panels looks like a loading effect and is not one: every
converted point is a trace-limit K0 = 1e9*n/(M*Kd), and Kd is a fitted
constant, so K0 carries no concentration dependence to show. The
measured and converted clusters also sit 1.8 to 3.4 decades apart with
no overlap, so any line through both measures the step between the
groups rather than a trend inside either. The colours say which group
a point belongs to, which is the thing worth seeing.

Qin is a horizontal line because it has no loading of its own: its
number is the Table S6 fit evaluated at this study's operating point,
not a measurement at a stated ratio.
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV = "albumin_loading_vs_logk_data.csv"
PNG = "albumin_loading_vs_logk.png"

WANT = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]        # chain-length order
NCF = {"PFBA": 3, "PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}

# Categorical slots 1, 2 and 3 of the validated palette. Checked with
# the dataviz validator, all pairs, light surface: every check passes,
# worst pair dE 9.2 deutan / 24.0 normal. Aqua carries a contrast WARN
# at 2.74:1, which the legend and the shipped CSV answer.
STYLE = {
    "measured":  ("#eb6834", "full", "measured directly"),
    "converted": ("#2a78d6", "full", "converted from a Kd"),
    "quenching": ("#1baf7a", "none", "converted on a quenching exponent"),
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
        if r["pfas"] not in WANT:
            continue
        r["log_k"] = float(r["log_k"])
        r["log_ratio"] = (float(r["log_ratio"])
                          if str(r["log_ratio"]).strip() else None)
        rows.append(r)
    return rows


def dodge(marks, width=0.085, step=0.115):
    """Lanes for marks sharing an x, over SORTED marks so the result
    does not depend on the order the file happened to be read in."""
    out = {}
    for key in {round(m["log_ratio"] / width) for m in marks}:
        same = sorted((m for m in marks
                       if round(m["log_ratio"] / width) == key),
                      key=lambda m: (m["log_k"], m["source"]))
        start = -step * (len(same) - 1) / 2.0
        for j, m in enumerate(same):
            out[id(m)] = start + j * step
    return out


def plot(rows=None, save=PNG):
    rows = load() if rows is None else rows
    by = defaultdict(list)
    for r in rows:
        by[r["pfas"]].append(r)
    order = [c for c in WANT if c in by]

    xs = [r["log_ratio"] for r in rows if r["log_ratio"] is not None]
    lo, hi = min(xs) - 0.6, max(xs) + 0.6
    ylo = min(r["log_k"] for r in rows) - 0.5
    yhi = max(r["log_k"] for r in rows) + 0.5

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.4), sharex=True,
                             sharey=True)
    flat = axes.ravel()

    for i, c in enumerate(order):
        ax = flat[i]
        here = by[c]
        qin = next((r["log_k"] for r in here if r["category"] == "qin"),
                   None)
        if qin is not None:
            ax.axhline(qin, color=MUTED, lw=1.1, ls=(0, (6, 3)), zorder=2)
            ax.annotate(f"Qin  {qin:.2f}", xy=(hi, qin), xytext=(-3, 3),
                        textcoords="offset points", ha="right",
                        va="bottom", fontsize=7.5, color=MUTED)

        pts = [r for r in here if r["log_ratio"] is not None]
        lanes = dodge(pts)
        for m in sorted(pts, key=lambda r: (r["log_ratio"], r["log_k"])):
            colour, fill, _ = STYLE[m["category"]]
            ax.plot(m["log_ratio"] + lanes[id(m)], m["log_k"], "o",
                    ms=7.5, mfc=colour if fill == "full" else "none",
                    mec=SURFACE if fill == "full" else colour,
                    mew=1.4 if fill == "full" else 1.8, zorder=6)

        ax.set_title(f"{c}   (C{NCF[c]})", loc="left", fontsize=11.5,
                     fontweight="bold", color=INK, pad=8)
        ax.grid(True, color=GRID, lw=0.7, zorder=1)
        ax.set_axisbelow(True)
        ax.set_xlim(lo, hi)
        ax.set_ylim(ylo, yhi)

    cell = flat[len(order)]
    cell.axis("off")
    handles = [
        Line2D([], [], marker="o", ls="none", ms=8,
               mfc=col if fill == "full" else "none",
               mec=SURFACE if fill == "full" else col,
               mew=1.4 if fill == "full" else 1.8, label=lab)
        for col, fill, lab in STYLE.values()
    ] + [Line2D([], [], color=MUTED, lw=1.1, ls=(0, (6, 3)),
                label="Qin at this study's operating point")]
    cell.legend(handles=handles, loc="upper left", frameon=False,
                fontsize=9.5, handletextpad=0.9, labelspacing=0.9,
                bbox_to_anchor=(0.0, 0.95))

    for ax in (axes[1][0], axes[1][1], axes[0][2]):
        ax.set_xlabel("log₁₀  ligand : protein ratio",
                      fontsize=9.5, labelpad=6)
        ax.tick_params(labelbottom=True)
    for ax in (axes[0][0], axes[1][0]):
        ax.set_ylabel("log K  (L water / kg protein)", fontsize=9.5)

    fig.suptitle("Albumin partition coefficient against assay loading",
                 x=0.04, y=0.975, ha="left", fontsize=14.5,
                 fontweight="bold", color=INK)
    fig.text(0.04, 0.935, "Coloured by how each number was arrived at.",
             ha="left", fontsize=10, color=INK2)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.915))
    if save:
        fig.savefig(save, dpi=200)
    return fig


fig = plot()
