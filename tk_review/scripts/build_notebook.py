#!/usr/bin/env python3
"""Write notebook/albumin_plots.ipynb -- the three figures, runnable.

WHY A GENERATOR.  A .ipynb is JSON with the source split into cells, which
is miserable to edit by hand and impossible to review in a diff.  The cells
live here as ordinary Python strings, so the notebook is generated, the
diff is readable, and `--py` dumps the same code as a flat script that can
be executed to prove the notebook runs before it ships.

NO PANDAS, NO PROJECT IMPORTS.  The notebook reads one CSV with the
standard library and draws with matplotlib.  It does not import anything
from tk_review/scripts, so it can be copied anywhere with its CSV beside
it.  An optional pandas cell at the end is for exploring, not for drawing.

Run:  python3 -I tk_review/scripts/build_notebook.py
      python3 -I tk_review/scripts/build_notebook.py --py  (flat script)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "notebook")
NB = os.path.join(OUT, "albumin_plots.ipynb")
PY = os.path.join(OUT, "albumin_plots.py")

MD_INTRO = """\
# PFAS albumin binding — the three figures

Everything here runs from one file, `albumin_logK_Kd_revised.csv`, with the
Python standard library and matplotlib. No pandas, no project imports.

**What the data is.** 285 protein–water partition coefficients for PFAS on
serum albumin, collected from the literature and put on one basis:

    log K = log10( litres of water per kg of protein )

Values measured as a dissociation constant were converted with

    log K = 4.1774 + log10(n) − log10(Kd in µM)

where 4.1774 = log10(1e9 / 66463) carries albumin's molar mass. The `Route`
column records how each row reached its number, from a direct measurement
down to a conversion resting on a stoichiometry the reader assumed. **Rows
in different routes are not the same kind of number** — that is what the
colour ramp in Figure 1 is for.

**Qin 2024** means the baseline-toxicity paper: Qin W, Henneberger L,
Glüge J, König M, Escher BI, *Environ Sci Technol* 2024;58(13):5727–5738,
doi 10.1021/acs.est.3c09950, Table 1, `log D_BSA/w`, converted from
L_water/L_protein to L_water/kg_protein by −log10(1.36) = −0.134.

Each Qin value is a Freundlich fit read at **IC₁₀,free**, the free
concentration causing 10% cytotoxicity. D_BSA/w falls about 0.3 log per
decade of free concentration, so these are values at one high loading —
0.3 to 12 PFAS per albumin — not trace-limit constants. Do not quote one
without its concentration.

A second 2024 paper by the same author (*Species Difference?…*, doi
10.1021/acs.est.3c10824) reports the same fits at a different
concentration, giving values up to 1.44 log higher. It is not in this CSV.

**Titles are deliberately absent** from all three figures — add them
yourself so the wording stays yours.
"""

SETUP = '''\
import csv
import math
import os
from collections import Counter, defaultdict

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV_PATH = "albumin_logK_Kd_revised.csv"   # <- edit if it lives elsewhere

# Arial where it exists, metric-compatible substitutes where it does not.
# Liberation Sans and Nimbus Sans have Arial's widths, so a layout tuned
# on one machine does not reflow on another.
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [
    "Arial", "Liberation Sans", "Nimbus Sans", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

# ---- colour -----------------------------------------------------------
# RAMP is ORDINAL, not categorical: Route is a quality ranking, so it gets
# one hue light-to-dark rather than five unrelated hues.  Dark = the most
# direct measurement.  ORANGE is the one categorical contrast in the whole
# set, reserved for Qin; the pair #2a78d6 / #eb6834 clears every check in
# the palette validator on this surface (worst CVD dE 24.7 protan).
RAMP = ["#0d366b", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef"]
ORANGE, GREY = "#eb6834", "#9a9a93"
INK, MUTED, SURFACE = "#1a1a17", "#6b6b64", "#fcfcfb"
GRID, SPINE = "#ececE6", "#d8d8d2"
PHYSIO = 0.01          # mol PFAS per mol albumin, the physiological limit

# ---- the Route ordering, most direct first ----------------------------
ROUTES = [
    ("measured directly", lambda r: r.startswith("measured")),
    ("converted, n counted by ITC", lambda r: "counted by ITC" in r),
    # Two data routes share this ramp step. The Excel keeps them
    # apart -- an n fitted from a binding isotherm is not an n = 1
    # fixed by a single-site model -- but a six-step single-hue
    # ordinal ramp cannot clear the validator's adjacent-lightness
    # gate on this surface, so the figure collapses the pair rather
    # than stepping the ramp past what it can carry. The same
    # collapse already applies to Kd inside vs beyond the range.
    ("converted, n from the paper's own fit",
     lambda r: "single-site model" in r or "fitted by the paper" in r),
    ("converted, n = the paper's quenching exponent",
     lambda r: "quenching-plot exponent" in r),
    ("converted, n = 1 assumed by the reader",
     lambda r: "assumed by the reader" in r),
]


def route_rank(route):
    for i, (_, test) in enumerate(ROUTES):
        if test(route):
            return i
    return len(ROUTES) - 1


def is_qin(source):
    return str(source).startswith("Qin 2024")
'''

LOAD = '''\
def num(v):
    """float, or None for a blank / unparseable cell."""
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def load(path=CSV_PATH):
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for d in csv.DictReader(fh):
            lk = num(d["log K (L water / kg protein)"])
            if lk is None:
                continue                       # no coefficient, no mark
            route = str(d["Route"])
            rows.append({
                "pfas": d["PFAS"], "species": d["Species"],
                "source": d["Source"], "log_k": lk, "route": route,
                "rank": route_rank(route), "qin": is_qin(d["Source"]),
                "ratio": num(d["ligand:protein ratio"]),
                "conc": num(d["PFAS conc (uM)"]),
                "prot": num(d["Protein conc (uM)"]),
            })
    return rows


rows = load()
print(f"{len(rows)} rows with a log K\\n")
print(f"{'route':50}{'n':>5}   log K range")
for i, (name, _) in enumerate(ROUTES):
    g = [r["log_k"] for r in rows if r["rank"] == i]
    if g:
        print(f"{name:50}{len(g):>5}   {min(g):5.2f} to {max(g):5.2f}")
print(f"\\nQin 2024 rows: {sum(r['qin'] for r in rows)}")
print("species:", dict(Counter(r["species"] for r in rows)))
'''

MD_F1 = """\
## Figure 1 — Qin 2024 against every other value, split by route

One row per compound, ordered by the Qin value. Each grey band is the full
span of the other values for that compound; each dot is one published
coefficient, **shaded by how it was obtained** and nudged off the centre
line by that same rank so the routes stay legible where they pile up.

Read it as a position judgement: does the orange diamond fall inside the
field, at its top, or outside it — and is the field made of direct
measurements or of conversions?
"""

FIG1 = '''\
def fig_route_split(rows, path="albumin_route_split.png"):
    qin_by = defaultdict(list)
    for r in rows:
        if r["qin"]:
            qin_by[r["pfas"]].append(r["log_k"])
    # only compounds that someone OTHER than Qin also reports
    order = [c for c in sorted(qin_by, key=lambda c: -max(qin_by[c]))
             if any(r["pfas"] == c and not r["qin"] for r in rows)]

    head = 0.92                       # strip left free for the legend
    fh = 0.62 * len(order) + head + 1.0
    fig, ax = plt.subplots(figsize=(10.6, fh))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    fig.subplots_adjust(top=1 - head / fh, bottom=1.0 / fh,
                        left=0.115, right=0.985)

    for i, cm in enumerate(order):
        y = len(order) - 1 - i
        others = [r for r in rows if r["pfas"] == cm and not r["qin"]]
        if others:
            lk = [r["log_k"] for r in others]
            ax.plot([min(lk), max(lk)], [y, y], color=GREY, lw=9,
                    alpha=0.16, solid_capstyle="round", zorder=1)
        for r in others:
            off = (r["rank"] - (len(ROUTES) - 1) / 2) * 0.115
            ax.scatter([r["log_k"]], [y + off], s=30, color=RAMP[r["rank"]],
                       edgecolor=SURFACE, linewidth=0.7, zorder=3)
        for v in qin_by[cm]:
            ax.scatter([v], [y], s=150, marker="D", color=ORANGE,
                       edgecolor=SURFACE, linewidth=2.1, zorder=5)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(list(reversed(order)), fontsize=10.5, color=INK)
    ax.set_xlabel("log K  (L water / kg protein)", fontsize=10.5,
                  color=MUTED)
    ax.tick_params(axis="x", colors=MUTED, labelsize=9.5)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(SPINE)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_ylim(-0.8, len(order) - 0.2)

    h = [Line2D([], [], marker="D", color="none", markerfacecolor=ORANGE,
                markeredgecolor=SURFACE, markersize=10, label="Qin 2024")]
    h += [Line2D([], [], marker="o", color="none", markerfacecolor=RAMP[i],
                 markeredgecolor=SURFACE, markersize=8, label=n)
          for i, (n, _) in enumerate(ROUTES)]
    ax.legend(handles=h, loc="lower left",
              bbox_to_anchor=(0.0, 1.005, 1.0, 0.14), ncol=2, frameon=False,
              fontsize=9.1, labelcolor=INK, handletextpad=0.45,
              columnspacing=1.5, borderaxespad=0)
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    print("wrote", path, f"({len(order)} compounds)")
    return fig


fig_route_split(rows);
'''

MD_F2 = """\
## Figure 2 — does loading explain the spread?

The hypothesis worth testing: albumin coefficients for one compound scatter
by two to four log across the literature, and binding is saturable, so a
coefficient measured at high ligand:protein is *genuinely* smaller than the
same coefficient at low loading. If the labs that disagree were simply
dosing differently, the spread would be an artefact of the operating point
rather than a disagreement about the chemistry.

x is the assay's ligand:protein ratio on a log scale, one panel per
compound for the five best-covered ones — pooling compounds would measure
the chain-length trend instead. The green band is physiological loading
(ν < 0.01): **almost nothing in the literature is measured there.**

Qin keeps its orange diamond from Figure 1 so both figures read with one
key.
"""

FIG2 = '''\
def fig_loading(rows, compounds=None, path="albumin_logk_vs_loading.png"):
    with_ratio = [r for r in rows if r["ratio"] and r["ratio"] > 0]
    if compounds is None:
        compounds = [c for c, _ in Counter(
            r["pfas"] for r in with_ratio).most_common(5)]

    n = len(compounds)
    fig, axes = plt.subplots(1, n, figsize=(2.55 * n, 4.75), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(top=0.78, bottom=0.155, left=0.075, right=0.985,
                        wspace=0.13)

    for ax, cm in zip(axes, compounds):
        ax.set_facecolor(SURFACE)
        g = [r for r in with_ratio if r["pfas"] == cm]
        ax.axvspan(1e-6, PHYSIO, color="#1baf7a", alpha=0.085, zorder=0)
        for r in g:                       # field first
            if not r["qin"]:
                ax.scatter([r["ratio"]], [r["log_k"]], s=34,
                           color=RAMP[r["rank"]], edgecolor=SURFACE,
                           linewidth=0.7, zorder=3)
        for r in g:                       # Qin on top, never hidden
            if r["qin"]:
                ax.scatter([r["ratio"]], [r["log_k"]], s=112, marker="D",
                           color=ORANGE, edgecolor=SURFACE, linewidth=1.8,
                           zorder=6)
        ax.set_xscale("log")
        ax.set_title(f"{cm}   (n={len(g)})", fontsize=10.5, color=INK,
                     pad=7)
        ax.set_xlabel("ligand : protein", fontsize=9.6, color=MUTED)
        ax.tick_params(colors=MUTED, labelsize=8.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("bottom", "left"):
            ax.spines[s].set_color(SPINE)
        ax.grid(color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        ax.set_xlim(1e-5, 1e5)
    axes[0].set_ylabel("log K  (L water / kg protein)", fontsize=10,
                       color=MUTED)

    h = [Line2D([], [], marker="s", color="none", markerfacecolor="#1baf7a",
                alpha=0.5, markeredgecolor="none", markersize=10,
                label="physiological loading (< 0.01)"),
         Line2D([], [], marker="D", color="none", markerfacecolor=ORANGE,
                markeredgecolor=SURFACE, markersize=9, label="Qin 2024")]
    h += [Line2D([], [], marker="o", color="none", markerfacecolor=RAMP[i],
                 markeredgecolor=SURFACE, markersize=8, label=nm)
          for i, (nm, _) in enumerate(ROUTES)]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.075, 0.998),
               ncol=3, frameon=False, fontsize=8.8, labelcolor=INK,
               handletextpad=0.4, columnspacing=1.4)
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    print("wrote", path, "-", ", ".join(compounds))
    return fig


fig_loading(rows);
# fig_loading(rows, ["PFOA", "PFOS", "PFHxS", "PFBA", "PFNA"]);  # pick your own
'''

MD_F3 = """\
## Figure 3 — the same comparison, split by measured vs converted

Figure 1 shades by all five routes. This one collapses them to the
distinction that matters most for trust — a directly measured partition
coefficient against one converted from a Kd — and draws the measured range
as a band so you can see whether Qin sits inside the directly measured
field or only inside the converted cloud.

`above` counts the compounds where Qin exceeds *every* directly measured
value. On the current data that is 1 of 11.
"""

FIG3 = '''\
def fig_qin_vs_all(rows, path="qin2024_albumin_comparison.png"):
    qin = {r["pfas"]: r["log_k"] for r in rows if r["qin"]}
    order = sorted((c for c in qin
                    if any(r["pfas"] == c and not r["qin"] for r in rows)),
                   key=lambda c: qin[c])

    fig, ax = plt.subplots(figsize=(10.4, 0.52 * len(order) + 2.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    above = 0

    for y, cm in enumerate(order):
        here = [r for r in rows if r["pfas"] == cm and not r["qin"]]
        meas = [r["log_k"] for r in here if r["rank"] == 0]
        conv = [r for r in here if r["rank"] > 0]
        allv = [r["log_k"] for r in here] + [qin[cm]]
        ax.plot([min(allv), max(allv)], [y, y], color=GRID, lw=1.4,
                zorder=1)
        if meas:
            ax.add_patch(plt.Rectangle(
                (min(meas), y - 0.30), max(meas) - min(meas) or 0.02, 0.60,
                facecolor=ORANGE, alpha=0.14, edgecolor="none", zorder=1))
            if qin[cm] > max(meas):
                above += 1
        for r in conv:
            hollow = "quenching" in r["route"]
            ax.scatter([r["log_k"]], [y], s=32,
                       facecolor="none" if hollow else MUTED,
                       edgecolor=MUTED, linewidth=1.2, zorder=3)
        for v in meas:
            ax.scatter([v], [y], s=46, color=ORANGE, edgecolor=SURFACE,
                       linewidth=1.3, zorder=4)
        ax.scatter([qin[cm]], [y], s=150, marker="D", color="#2a78d6",
                   edgecolor=SURFACE, linewidth=1.8, zorder=5)
        ax.annotate(f"{qin[cm]:.2f}", xy=(1.0, y),
                    xycoords=("axes fraction", "data"), xytext=(12, 0),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=9.6, color="#2a78d6", weight="bold",
                    annotation_clip=False)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=10.5, color=INK)
    ax.set_xlabel("log K  (L water / kg protein)", fontsize=10.5,
                  color=MUTED)
    ax.tick_params(axis="x", colors=MUTED, labelsize=9.5)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(SPINE)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_ylim(-0.7, len(order) - 0.3)

    h = [Line2D([], [], marker="D", ls="", ms=9.5, mfc="#2a78d6",
                mec=SURFACE, mew=1.6,
                label="Qin 2024, Table 1: the fit read at IC10,free"),
         Line2D([], [], marker="o", ls="", ms=7, mfc=ORANGE, mec=SURFACE,
                mew=1.4,
                label="measured directly (band = their range)"),
         Line2D([], [], marker="o", ls="", ms=5.5, mfc=MUTED, mec=MUTED,
                label="converted from a dissociation constant"),
         Line2D([], [], marker="o", ls="", ms=5.5, mfc="none", mec=MUTED,
                mew=1.2, label="converted on a quenching exponent")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(0.0, -0.10),
              ncol=2, fontsize=8.6, frameon=False, labelspacing=0.7,
              handletextpad=0.6, columnspacing=1.6, borderaxespad=0.0)
    fig.tight_layout(rect=(0.0, 0.07, 0.955, 1.0))
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    print("wrote", path,
          f"- Qin above every directly measured value on {above} "
          f"of {len(order)}")
    return fig


fig_qin_vs_all(rows);
'''

MD_END = """\
## Changing things

* **Another five compounds** — `fig_loading(rows, ["PFOA", "PFOS", ...])`.
* **Drop a route you don't trust** — filter before plotting, e.g.
  `rows = [r for r in rows if r["rank"] < 3]` removes the quenching
  conversions and the reader-assumed n = 1.
* **One species only** — `rows = [r for r in rows if r["species"] == "human"]`.
  Note Qin is *bovine*, so filtering to human removes the orange diamonds.
* **Titles** — none of the three draws one. Add with
  `fig.suptitle("...")` or in whatever you paste them into.

### Two things not to do with this table

1. **Do not take a mean down the log K column.** PFOA on human albumin
   appears at 1.64, 2.12, 2.67, 3.77, 4.48 and 4.70. That thousand-fold
   spread is the finding, not noise to average away — two labs measuring
   the same albumin Kd by the same method differ by 1.45 log.
2. **Do not compare across routes as if they were the same measurement.**
   A quenching-exponent conversion and a dialysis measurement are
   different quantities that happen to share a unit.

### Optional: the same data as a pandas DataFrame

Nothing above needs pandas. This is only for exploring.
"""

PANDAS = '''\
try:
    import pandas as pd
    df = pd.read_csv(CSV_PATH)
    df["rank"] = df["Route"].map(route_rank)
    df["qin"] = df["Source"].map(is_qin)
    display(df.groupby("Route")["log K (L water / kg protein)"]
              .agg(["count", "min", "median", "max"]).round(2))
except ImportError:
    print("pandas not installed - the figures above do not need it")
'''

CELLS = [("md", MD_INTRO), ("code", SETUP), ("code", LOAD),
         ("md", MD_F1), ("code", FIG1),
         ("md", MD_F2), ("code", FIG2),
         ("md", MD_F3), ("code", FIG3),
         ("md", MD_END), ("code", PANDAS)]


def lines(text):
    out = text.splitlines(keepends=True)
    if out and out[-1].endswith("\n"):
        out[-1] = out[-1][:-1]
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    if "--py" in sys.argv:
        body = "\n\n".join(
            s for k, s in CELLS if k == "code" and "import pandas" not in s)
        with open(PY, "w", encoding="utf-8") as fh:
            fh.write('"""Flat version of albumin_plots.ipynb, for testing."""\n'
                     + body)
        print("wrote", os.path.relpath(PY))
        return
    nb = {"cells": [], "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3"}},
        "nbformat": 4, "nbformat_minor": 5}
    for kind, src in CELLS:
        cell = {"cell_type": "markdown" if kind == "md" else "code",
                "metadata": {}, "source": lines(src)}
        if kind == "code":
            cell["execution_count"] = None
            cell["outputs"] = []
        nb["cells"].append(cell)
    with open(NB, "w", encoding="utf-8") as fh:
        json.dump(nb, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"wrote {os.path.relpath(NB)}  ({len(nb['cells'])} cells)")


if __name__ == "__main__":
    main()
