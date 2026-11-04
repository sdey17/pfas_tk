"""Every L-FABP value, one panel per PFAS, by who produced it.

Notebook script. It reads one CSV from the working directory and plots
it. No paths, no arguments: keep the CSV beside the notebook.

The companion to plot_albumin_by_compound.py, same layout and the same
three colours, so the two can be read side by side. What differs is
what the colour MEANS, and the difference is the whole point.

On albumin, colour carries how the number was arrived at, because some
albumin partition coefficients are measured directly and some are
converted from a Kd. On L-FABP there is nothing to choose between:
no paper reports a measured partition coefficient for any protein but
albumin, so every value here is converted. What varies between these
laboratories is the ASSAY, and that is what the colour carries:

  separation    equilibrium dialysis -- free and bound separated,
                then both measured
  direct        ITC or MST -- binding watched directly, no probe
                and no separation step
  displacement  a probe (ANS, DAUDA, a fatty acid) is competed off
                and the PFAS constant inferred from the signal lost

HALF THIS LITERATURE IS DISPLACEMENT, which is not the same quantity
as a direct Kd: it depends on the probe's own affinity and on the
model used to back the constant out. A row measured by displacement
and a row measured by ITC are not two attempts at one number.

MUTANT ROWS ARE SEPARATE ROWS. Sheng 2016 measured wild-type hL-FABP
and four site-directed point variants. Those are different proteins.
N111D is the one to look at: by ITC it takes PFOA from 6.4 to 284 uM,
a 44-fold weakening and 1.65 log, while R122G barely moves.

AND THE TWO ASSAYS DISAGREE ABOUT THOSE MUTANTS, which is the sharpest
thing on this figure. Fluorescence displacement reports NO BINDING for
N111D and R122G on both PFOA and PFNA; ITC on the same proteins returns
8.4 uM for R122G. Those no-binding rows cannot go on a log axis, so
they are in lfabp_by_compound_excluded.csv rather than silently absent.
A row that carries an orange point and no green one is where this
happened.

SPECIES ARE GROUPED AND LABELLED, never pooled -- human and rat L-FABP
are different proteins, as on the albumin plots.

EDIT `COMPOUNDS` to change which panels appear.
"""
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV = "lfabp_by_compound_data.csv"
PNG = "lfabp_by_compound.png"

COMPOUNDS = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]
SPECIES_ORDER = ["human", "rat"]

# Categorical slots 1, 2 and 3 of the validated palette, the same three
# the albumin plots use so the whole set reads as one family. Checked
# all-pairs on the light surface: every check passes, worst pair
# dE 9.2 deutan / 24.0 normal. Aqua carries a contrast WARN at 2.74:1,
# which the legend and the shipped CSV answer.
#
# Three colours, not four. The four assays would have wanted a fourth
# slot, and no fourth slot clears the all-pairs floor against these
# three -- yellow collides with orange at dE 13.7 normal, below the 15
# floor. So ITC and MST share one category, named for what they have in
# common: binding watched directly, with no probe and no separation.
STYLE = {
    "separation":   ("#2a78d6", "full",
                     "equilibrium dialysis"),
    "direct":       ("#eb6834", "full",
                     "ITC or MST, binding watched directly"),
    "displacement": ("#1baf7a", "none",
                     "probe displacement (ANS, DAUDA, a fatty acid)"),
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

    # One y slot per source-protein-species triple. The protein is in
    # the key, not just the label, so Sheng's mutants never share a row
    # with the wild type they are being compared against.
    keys = sorted({(r["species"], r["source"], r["protein_label"])
                   for r in rows},
                  key=lambda k: (SPECIES_ORDER.index(k[0])
                                 if k[0] in SPECIES_ORDER else 9,
                                 k[1], k[2]))
    slot = {k: i for i, k in enumerate(reversed(keys))}
    labels = [f"{s}  ·  {sp}" + ("" if p == "L-FABP" else f"  ·  {p[7:]}")
              for sp, s, p in reversed(keys)]

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
    # Title strip in INCHES, so the heading does not land on the
    # subtitle when the panel or row count changes.
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

        rowed = defaultdict(list)
        for r in here:
            rowed[(r["species"], r["source"], r["protein_label"])].append(r)
        for k, group in rowed.items():
            y0 = slot[k]
            group = sorted(group, key=lambda r: r["log_k"])
            n = len(group)
            step = min(0.17, 0.78 / (n - 1)) if n > 1 else 0.0
            shrink = 1.0 if n <= 4 else 0.82
            start = -step * (n - 1) / 2.0
            for j, r in enumerate(group):
                colour, fill, _ = STYLE[r["method"]]
                ax.plot(r["log_k"], y0 + start + j * step, "o",
                        ms=6.8 * shrink,
                        mfc=colour if fill == "full" else "none",
                        mec=SURFACE if fill == "full" else colour,
                        mew=1.3 if fill == "full" else 1.7, zorder=6)

        ax.set_title(f"{c}   (C{ncf[c]})   n = {len(here)}", loc="left",
                     fontsize=11.5, fontweight="bold", color=INK, pad=7)
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.7, len(keys) - 0.3)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels(labels, fontsize=8.4)
        ax.tick_params(axis="y", length=0)
        ax.grid(True, axis="x", color=GRID, lw=0.7, zorder=1)
        ax.set_axisbelow(True)
        if len(order) - i <= ncol:
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
    ]
    leg = flat[len(order)].legend(
        handles=handles, loc="upper left", frameon=False, fontsize=9.5,
        handletextpad=0.9, labelspacing=1.0, title="how it was measured",
        bbox_to_anchor=(0.02, 0.95))
    leg.get_title().set_fontsize(9.5)
    leg.get_title().set_color(INK2)
    leg._legend_box.align = "left"
    flat[len(order)].text(
        0.02, 0.30,
        "Every value here is converted from a Kd.\n"
        "No paper reports a measured partition\n"
        "coefficient for any protein but albumin.\n\n"
        "Sheng 2016's point variants are their own\n"
        "rows: they are different proteins. N111D is\n"
        "the knockout -- PFOA 6.4 to 284 uM by ITC.\n\n"
        "Where a row has an orange point and no green\n"
        "one, displacement reported NO BINDING on the\n"
        "protein ITC still bound. See the excluded CSV.",
        transform=flat[len(order)].transAxes, ha="left", va="top",
        fontsize=9.0, color=INK2, linespacing=1.5)

    fig.suptitle("Every L-FABP value, one panel per PFAS",
                 x=0.012, y=1 - 0.18 / h, ha="left", va="top",
                 fontsize=14.5, fontweight="bold", color=INK)
    fig.text(0.012, 1 - 0.50 / h,
             "Rows are shared across panels, so a laboratory can be "
             "followed from one compound to the next. Grouped by "
             "species; coloured by the assay, which is what differs "
             "between these laboratories.",
             ha="left", va="top", fontsize=9.6, color=INK2)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1 - head / h))
    if save:
        fig.savefig(save, dpi=200)
    print(f"{len(rows)} points, {len(keys)} rows, "
          f"{len(order)} panels")
    return fig


fig = plot()
