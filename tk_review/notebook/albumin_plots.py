"""Flat version of albumin_plots.ipynb, for testing."""
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
    ("converted, n from the paper's single-site model",
     lambda r: "single-site model" in r),
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
print(f"{len(rows)} rows with a log K\n")
print(f"{'route':50}{'n':>5}   log K range")
for i, (name, _) in enumerate(ROUTES):
    g = [r["log_k"] for r in rows if r["rank"] == i]
    if g:
        print(f"{name:50}{len(g):>5}   {min(g):5.2f} to {max(g):5.2f}")
print(f"\nQin 2024 rows: {sum(r['qin'] for r in rows)}")
print("species:", dict(Counter(r["species"] for r in rows)))


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
