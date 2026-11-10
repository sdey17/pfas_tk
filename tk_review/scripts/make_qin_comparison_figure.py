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

EDITING IT. The chart is drawn from a CSV, not from the workbook, so
every mark can be changed without touching this file:

    python3 .../make_qin_comparison_figure.py
        reads db/albumin_logK_shareable.xlsx, writes the CSV AND the PNG

    python3 .../make_qin_comparison_figure.py mine.csv
        reads mine.csv, writes the PNG only -- the workbook is not opened

So: run it once, edit
figures/qin2024_albumin_comparison_data.csv, run it again with that file
as the argument. Drop rows to remove marks, change `category` to recolour
one, change `row_order` to reorder or drop a compound from the chart
(every row of a compound shares its number). Colours, marker sizes and the
title live in the constants at the top of this file.

Writes figures/qin2024_albumin_comparison.png and
figures/qin2024_albumin_comparison_data.csv.

Run:  python3 tk_review/scripts/make_qin_comparison_figure.py [data.csv]
"""
import csv
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")
OUT = os.path.join(HERE, "..", "figures", "qin2024_albumin_comparison.png")
OUT_CSV = os.path.join(HERE, "..", "figures",
                       "qin2024_albumin_comparison_data.csv")

# category -> how it is drawn. Change `category` in the CSV to move a mark
# between these; add a key here to invent one.
CATEGORIES = ("qin", "measured", "converted", "converted_quenching")
FIELDS = ["row_order", "pfas", "category", "log_k", "source", "cw_uM",
          "route"]

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
QIN_C, MEAS_C = "#2a78d6", "#eb6834"

# No title or subtitle is drawn: the project adds those by hand.  main()
# prints the caption instead, so its numbers stay tied to the data.
QIN_LABEL = "Qin 2024, Table 1: the Freundlich fit read at IC10,free"
CAPTION = (
    "Qin 2024 sits above every directly measured value for {above} of {n} "
    "compounds. Each Qin point is its Table 1 log D_BSA/w: the Freundlich "
    "fit read at IC10,free, the free concentration causing 10% "
    "cytotoxicity, {lo:.3g} to {hi:.3g} uM. That is inside or above the "
    "range the isotherm was measured over, so these are interpolations "
    "rather than extrapolations -- but it is also a loading of 0.3 to 12 "
    "PFAS per albumin, far past the physiological trace limit. The grey "
    "cloud is the literature's own spread, up to 4.6 log on one compound.")
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


def category_of(route):
    if route.startswith("measured; Qin"):
        return "qin"
    if route == "measured directly":
        return "measured"
    return "converted_quenching" if "quenching" in route else "converted"


def from_workbook():
    """The plotted data, derived once from the shareable workbook."""
    from openpyxl import load_workbook
    ws = load_workbook(SRC)["Coefficients"]
    head = [c.value for c in ws[1]]
    ix = {h: i for i, h in enumerate(head)}
    rows = []
    for r in ws.iter_rows(min_row=2, values_only=True):
        v = r[ix["log K (L water / kg protein)"]]
        if v is None:
            continue
        route = str(r[ix["Route"]])
        # The anchor names the free concentration it was read at, as
        # "Cw = x uM" or "IC10,free = x uM" depending on which paper
        # qin2024_anchor.py has active.  Take either.
        cw = ""
        for tag in ("IC10,free = ", "Cw = "):
            if tag in route:
                cw = route.split(tag)[1].split(" uM")[0]
                break
        rows.append({"pfas": str(r[ix["PFAS"]]), "log_k": float(v),
                     "category": category_of(route),
                     "source": str(r[ix["Source"]]), "cw_uM": cw,
                     "route": route})
    # only the compounds the anchor covers, ordered by the anchor's value
    qin = {r["pfas"]: r["log_k"] for r in rows if r["category"] == "qin"}
    rows = [r for r in rows if r["pfas"] in qin]
    order = {c: i for i, c in enumerate(sorted(qin, key=qin.get))}
    for r in rows:
        r["row_order"] = order[r["pfas"]]
    return rows


def write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (-r["row_order"], r["category"],
                                             r["log_k"])):
            w.writerow({k: r.get(k, "") for k in FIELDS})


def read_csv(path):
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if not str(r.get("log_k", "")).strip():
            continue
        r["log_k"] = float(r["log_k"])
        r["row_order"] = int(float(r["row_order"]))
        if r["category"] not in CATEGORIES:
            raise SystemExit(f"unknown category {r['category']!r} on "
                             f"{r['pfas']}; expected one of "
                             f"{', '.join(CATEGORIES)}")
        out.append(r)
    return out


def dodge(marks, width=0.085, step=0.115):
    """Vertical offsets so near-coincident dots stop hiding each other.

    A plain dot plot stacks every value within a few hundredths of a log
    on one pixel, which turns the densest rows -- the ones carrying the
    most evidence -- into a single blob. Marks are swept left to right and
    anything inside `width` of the last mark in a lane is pushed onto the
    next one, alternating above and below the row.

    Lanes are assigned over the marks SORTED, not in the order they were
    read, so the same data draws the same picture whether it came from the
    workbook or from an edited CSV. Keying on the x value alone was not
    enough: two marks at the same x took lanes in call order, and the
    round trip through the CSV reorders them.
    """
    out = []
    lanes = []
    for i, m in sorted(enumerate(marks), key=lambda t: (t[1][0], t[1][1])):
        x = m[0]
        lane = 0
        while lane < len(lanes) and abs(lanes[lane] - x) < width:
            lane += 1
        if lane == len(lanes):
            lanes.append(x)
        else:
            lanes[lane] = x
        out.append((i, 0.0 if lane == 0 else
                    ((lane + 1) // 2) * step * (1 if lane % 2 else -1)))
    return [o for _, o in sorted(out)]


def main():
    if len(sys.argv) > 1:
        rows = read_csv(sys.argv[1])
        src = sys.argv[1]
    else:
        rows = from_workbook()
        write_csv(rows, OUT_CSV)
        src = os.path.relpath(SRC, os.path.join(HERE, ".."))

    order = sorted({r["row_order"] for r in rows})
    names = {}
    for r in rows:
        names.setdefault(r["row_order"], r["pfas"])
    qin = {r["row_order"]: r["log_k"] for r in rows
           if r["category"] == "qin"}

    fig, ax = plt.subplots(figsize=(10.4, 7.6))
    lo = min(r["log_k"] for r in rows) - 0.4
    hi = max(r["log_k"] for r in rows) + 0.25
    above = 0
    cws = [float(r["cw_uM"]) for r in rows
           if r["category"] == "qin" and str(r["cw_uM"]).strip()]

    for y, key in enumerate(order):
        here = [r for r in rows if r["row_order"] == key]
        meas = [r["log_k"] for r in here if r["category"] == "measured"]
        conv = [(r["log_k"], r["category"] == "converted_quenching")
                for r in here if r["category"].startswith("converted")]
        vals = [v for v, _ in conv] + meas + (
            [qin[key]] if key in qin else [])
        if not vals:
            continue

        ax.plot([min(vals), max(vals)], [y, y], color=GRID, lw=1.4,
                zorder=1, solid_capstyle="round")
        if meas:
            ax.add_patch(Rectangle((min(meas), y - 0.26),
                                   max(meas) - min(meas) or 0.012, 0.52,
                                   facecolor=BAND, edgecolor="none",
                                   zorder=2))
        # (x, sort-tiebreak, kind) -- the tiebreak keeps lane assignment
        # independent of the order the marks were read in
        # Sorted before drawing, not only before lane assignment: marks
        # share a zorder, so the one drawn last wins the overlap, and an
        # unsorted draw made the workbook and the CSV render differently
        # from identical data.
        marks = sorted([(v, "a_quench" if q else "b_conv",
                         "quench" if q else "conv") for v, q in conv]
                       + [(v, "c_meas", "meas") for v in meas])
        offs = dodge(marks)
        for (v, _t, kind), dy in zip(marks, offs):
            if kind == "meas":
                ax.plot([v], [y + dy], "o", ms=7.0, mfc=MEAS_C,
                        mec=SURFACE, mew=1.4, zorder=4)
            else:
                ax.plot([v], [y + dy], "o", ms=5.5, zorder=3,
                        mfc="none" if kind == "quench" else MUTED,
                        mec=MUTED, mew=1.2 if kind == "quench" else 0.9)
        if key in qin:
            qv = qin[key]
            ax.plot([qv], [y], "D", ms=10.0, mfc=QIN_C, mec=SURFACE,
                    mew=1.8, zorder=6)
            # The value sits in a column at the right edge, not beside the
            # marker: offset labels collided with the dots on seven of
            # eleven rows, and "4.62" overlaid by a dot read as "1.62".
            ax.annotate(f"{qv:.2f}", xy=(1.0, y),
                        xycoords=("axes fraction", "data"), xytext=(10, 0),
                        textcoords="offset points", va="center", ha="left",
                        fontsize=9, color=QIN_C, fontweight="bold",
                        annotation_clip=False)
            if meas and qv > max(meas):
                above += 1

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([names[k] for k in order], fontsize=9.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("log K   (litres of water per kg of protein)")
    ax.grid(True, axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)

    handles = [
        Line2D([], [], marker="D", ls="", ms=9.5, mfc=QIN_C, mec=SURFACE,
               mew=1.6, label=QIN_LABEL),
        Line2D([], [], marker="o", ls="", ms=7.0, mfc=MEAS_C, mec=SURFACE,
               mew=1.4, label="measured directly by another laboratory "
                              "(band = their range)"),
        Line2D([], [], marker="o", ls="", ms=5.5, mfc=MUTED, mec=MUTED,
               label="converted from a dissociation constant"),
        Line2D([], [], marker="o", ls="", ms=5.5, mfc="none", mec=MUTED,
               mew=1.2, label="converted on a fluorescence-quenching "
                              "exponent"),
    ]
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(0.0, -0.085), ncol=2, fontsize=8.4,
              frameon=False, labelspacing=0.7, handletextpad=0.6,
              columnspacing=1.6, borderaxespad=0.0)

    # No title or subtitle on the canvas: they are added by hand
    # downstream.  draw() returns the caption text it would have written,
    # so the numbers in it stay tied to the data rather than to memory.
    ax.annotate("Qin", xy=(1.0, len(order) - 0.4),
                xycoords=("axes fraction", "data"), xytext=(10, 6),
                textcoords="offset points", va="bottom", ha="left",
                fontsize=8.4, color=MUTED, annotation_clip=False)
    fig.tight_layout(rect=(0.0, 0.085, 0.955, 0.985))
    fig.savefig(OUT, dpi=200)
    plt.close(fig)
    print("\ncaption (not drawn on the figure):")
    for line in textwrap.wrap(CAPTION.format(
            above=above, n=len(order), lo=min(cws), hi=max(cws)), 74):
        print("  " + line)

    print(f"read {len(rows)} marks from {src}")
    print(f"  {len(order)} compounds; Qin above all measured on {above}")
    if len(sys.argv) == 1:
        print(f"wrote {os.path.abspath(OUT_CSV)}  <- edit this, then pass "
              f"it back as an argument")
    print(f"wrote {os.path.abspath(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
