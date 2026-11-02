#!/usr/bin/env python3
"""Qin 2024 against the rest of the albumin literature, at one operating point.

THE QUESTION THE CHART ANSWERS. A study is using Qin's D_BSA/w at 10 uM
nominal in 10% FBS. Is that number an outlier against everything else
measured on albumin, and against what quality of comparator?

THE FORM. One row per compound, every albumin coefficient as a dot on a
log K axis. The reader's job is a position judgement -- where does the
anchor sit inside the field -- which is what a dot plot answers at a
glance and what a bar chart or a table of means destroys.

The measured range is drawn as a band so Qin's marker falling inside or
outside it is the whole reading. Converted values recede to a neutral,
because they are context rather than comparator: the argument is between
Qin and the other MEASUREMENTS, and the converted cloud says how much
noise the literature carries around both.

COLOUR. Two categorical hues only (#2a78d6 blue, #eb6834 orange), which
clears every all-pairs gate in the validator including contrast. The
converted values take the muted chrome grey rather than a third hue --
they are the background field, and a hue would claim they are a peer of
the other two. Quenching-derived conversions are the same grey as an
OPEN marker: a secondary channel, so the least trustworthy route is
flagged without spending a colour on it.

Writes figures/qin2024_albumin_comparison.png.

Run:  python3 tk_review/scripts/make_qin_comparison_figure.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")
OUT = os.path.join(HERE, "..", "figures", "qin2024_albumin_comparison.png")

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
QIN_C, MEAS_C = "#2a78d6", "#eb6834"
BAND = "#fbe4da"            # a pale wash of the measured hue

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False,
})


def read():
    from openpyxl import load_workbook
    ws = load_workbook(SRC)["Coefficients"]
    head = [c.value for c in ws[1]]
    ix = {h: i for i, h in enumerate(head)}
    rows = []
    for r in ws.iter_rows(min_row=2, values_only=True):
        rows.append({"pfas": r[ix["PFAS"]], "logK": r[ix["log K (L water / kg protein)"]],
                     "route": str(r[ix["Route"]]), "source": str(r[ix["Source"]])})
    return [r for r in rows if r["logK"] is not None]


def dodge(xs, width=0.085, step=0.115):
    """Vertical offsets so near-coincident dots stop hiding each other.

    A plain dot plot stacks every value within a few hundredths of a log
    on one pixel, which turns the densest rows -- the ones carrying the
    most evidence -- into a single blob. Points are swept left to right
    and anything inside `width` of its predecessor is pushed onto the
    next lane, alternating above and below the row.
    """
    out, lanes = {}, []
    for x in sorted(xs):
        lane = 0
        while lane < len(lanes) and abs(lanes[lane] - x) < width:
            lane += 1
        if lane == len(lanes):
            lanes.append(x)
        else:
            lanes[lane] = x
        out.setdefault(x, []).append(
            0.0 if lane == 0 else ((lane + 1) // 2) * step *
            (1 if lane % 2 else -1))
    return out


def main():
    rows = read()
    qin = {r["pfas"]: (float(r["logK"]), r["route"]) for r in rows
           if r["route"].startswith("measured; Qin")}
    order = sorted(qin, key=lambda c: qin[c][0])          # low to high, up the page

    fig, ax = plt.subplots(figsize=(10.4, 7.6))
    lo = min(float(r["logK"]) for r in rows if r["pfas"] in qin) - 0.4
    hi = max(float(r["logK"]) for r in rows if r["pfas"] in qin) + 0.25
    above = 0

    for y, c in enumerate(order):
        here = [r for r in rows if r["pfas"] == c]
        meas = [float(r["logK"]) for r in here
                if r["route"] == "measured directly"]
        conv = [(float(r["logK"]), "quenching" in r["route"]) for r in here
                if r["route"].startswith("converted")]
        vals = [v for v, _ in conv] + meas + [qin[c][0]]

        ax.plot([min(vals), max(vals)], [y, y], color=GRID, lw=1.4,
                zorder=1, solid_capstyle="round")
        if meas:
            ax.add_patch(Rectangle((min(meas), y - 0.26),
                                   max(meas) - min(meas) or 0.012, 0.52,
                                   facecolor=BAND, edgecolor="none", zorder=2))
        off = dodge([v for v, _ in conv] + meas)
        used = {}
        def lane(v):
            i = used.get(v, 0)
            used[v] = i + 1
            return off[v][min(i, len(off[v]) - 1)]
        for v, is_q in conv:
            ax.plot([v], [y + lane(v)], "o", ms=5.5, zorder=3,
                    mfc="none" if is_q else MUTED,
                    mec=MUTED, mew=1.2 if is_q else 0.9)
        for v in meas:
            ax.plot([v], [y + lane(v)], "o", ms=7.0, mfc=MEAS_C,
                    mec=SURFACE, mew=1.4, zorder=4)
        qv = qin[c][0]
        ax.plot([qv], [y], "D", ms=10.0, mfc=QIN_C, mec=SURFACE, mew=1.8,
                zorder=6)
        # The value sits in a column at the right edge, not beside the
        # marker: offset labels collided with the dots on seven of eleven
        # rows, and "4.62" overlaid by an orange dot read as "1.62".
        ax.annotate(f"{qv:.2f}", xy=(1.0, y),
                    xycoords=("axes fraction", "data"), xytext=(10, 0),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=9, color=QIN_C, fontweight="bold",
                    annotation_clip=False)
        if meas and qv > max(meas):
            above += 1

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=9.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("log K   (litres of water per kg of protein)")
    ax.grid(True, axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)

    handles = [
        Line2D([], [], marker="D", ls="", ms=9.5, mfc=QIN_C, mec=SURFACE,
               mew=1.6, label="Qin 2024, Table S6 equation at this "
                              "study's operating point"),
        Line2D([], [], marker="o", ls="", ms=7.0, mfc=MEAS_C, mec=SURFACE,
               mew=1.4, label="measured directly by another laboratory "
                              "(band = their range)"),
        Line2D([], [], marker="o", ls="", ms=5.5, mfc=MUTED, mec=MUTED,
               label="converted from a dissociation constant"),
        Line2D([], [], marker="o", ls="", ms=5.5, mfc="none", mec=MUTED,
               mew=1.2, label="converted on a fluorescence-quenching "
                              "exponent"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.085),
              ncol=2, fontsize=8.4, frameon=False, labelspacing=0.7,
              handletextpad=0.6, columnspacing=1.6, borderaxespad=0.0)

    fig.text(0.052, 0.975,
             "Qin 2024 against every other albumin coefficient",
             ha="left", va="top", fontsize=14.5, color=INK)
    fig.text(0.052, 0.940,
             f"At 10 µM nominal in 10% FBS, Qin sits above every "
             f"directly measured value for {above} of {len(order)} compounds. "
             "Free concentrations\n(0.011–3.8 µM) are 2–3 orders "
             "below where Qin measured, so these are the Table S6 fits "
             "extrapolated,\nnot read off data. The grey cloud is the "
             "literature's own spread, up to 4.6 log on one compound.",
             ha="left", va="top", fontsize=8.8, color=INK2, linespacing=1.6)
    ax.annotate("Qin", xy=(1.0, len(order) - 0.4),
                xycoords=("axes fraction", "data"), xytext=(10, 6),
                textcoords="offset points", va="bottom", ha="left",
                fontsize=8.4, color=MUTED, annotation_clip=False)
    fig.tight_layout(rect=(0.0, 0.085, 0.955, 0.878))
    fig.savefig(OUT, dpi=200)
    plt.close(fig)
    print(f"{len(order)} compounds; Qin above all measured on {above}")
    print(f"wrote {os.path.abspath(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
