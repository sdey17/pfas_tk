#!/usr/bin/env python3
"""Two views of the albumin field: by ROUTE, and against LOADING.

FIGURE A -- Qin 2024 against the literature, split by route.
  The earlier version of this comparison drew every albumin
  coefficient as one undifferentiated cloud.  That hides the thing a
  reader most needs: the cloud is not one population.  A value that
  was MEASURED and a value CONVERTED from a Kd with a stoichiometry
  the reader supplied are different kinds of number, and the spread
  between them is part of what the chart is being asked about.
  Splitting by route turns one blurred band into a stack of bands
  whose widths can be compared.

  ROUTE IS ORDINAL, NOT CATEGORICAL, and that decides the colour.
  The five routes run from a direct measurement to a conversion
  resting entirely on an assumption the paper never made.  That is an
  ordered scale, so it takes a single-hue ramp dark-to-light rather
  than five hues -- which is also what makes it drawable: five
  categorical hues cannot all clear the separation gate on a light
  surface, and a sequential ramp sidesteps the problem by encoding
  the order that is actually there.  Steps #86b6ef / #5598e7 /
  #2a78d6 / #1c5cab / #0d366b pass the skill validator's ordinal
  checks (monotone L, all adjacent gaps >= 0.06, light end 2.06:1 vs
  surface, hue spread 4 degrees).

FIGURE B -- log K against the loading it was measured at.
  One panel per compound for the five best-covered PFAS, log K on y,
  ligand:protein ratio on x, log scale.  The physiological band
  (nu < 0.01) is shaded, so the reader can see at a glance how much
  of the literature sits outside the regime it is being used to
  describe.

  This is the chart the whole loading argument rests on: if nu did
  not matter, these panels would be flat.

Writes figures/albumin_route_split.png, figures/albumin_logk_vs_loading.png
and the CSV each is drawn from.
"""
import collections
import csv
import math
import os
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
FIGS = os.path.join(HERE, "..", "figures")

# Ordinal: most direct first. The ramp is applied in this order.
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
RAMP = ["#0d366b", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef"]
ORANGE, GREY = "#eb6834", "#9a9a93"
INK, MUTED, SURFACE = "#1a1a17", "#6b6b64", "#fcfcfb"
PHYSIO = 0.01


def route_of(r):
    for i, (name, test) in enumerate(ROUTES):
        if test(r):
            return i, name
    return len(ROUTES) - 1, ROUTES[-1][0]


def load():
    ws = openpyxl.load_workbook(SRC, data_only=True)["Coefficients"]
    h = [c.value for c in ws[1]]
    out = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        d = dict(zip(h, row))
        if d["log K (L water / kg protein)"] in (None, ""):
            continue
        i, name = route_of(str(d["Route"]))
        try:
            ratio = float(str(d["ligand:protein ratio (max)"]))
        except (TypeError, ValueError):
            ratio = None
        out.append(dict(pfas=str(d["PFAS"]), species=str(d["Species"]),
                        source=str(d["Source"]),
                        log_k=float(d["log K (L water / kg protein)"]),
                        route=name, route_rank=i, ratio=ratio))
    return out


def fig_a(rows, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    qin = [r for r in rows if r["source"].startswith("Qin 2024")]
    qin_by = collections.defaultdict(list)
    for r in qin:
        qin_by[r["pfas"]].append(r["log_k"])
    order = [c for c in sorted(qin_by, key=lambda c: -max(qin_by[c]))
             if any(r["pfas"] == c and not r["source"].startswith("Qin")
                    for r in rows)]

    head = 1.55
    fh = 0.62 * len(order) + head + 1.0
    fig, ax = plt.subplots(figsize=(10.6, fh))
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    fig.subplots_adjust(top=1 - head / fh, bottom=1.0 / fh,
                        left=0.115, right=0.985)

    for i, cm in enumerate(order):
        y = len(order) - 1 - i
        others = [r for r in rows
                  if r["pfas"] == cm and not r["source"].startswith("Qin 2024")]
        if others:
            lk = [r["log_k"] for r in others]
            ax.plot([min(lk), max(lk)], [y, y], color=GREY, lw=9, alpha=0.16,
                    solid_capstyle="round", zorder=1)
        # jitter within the row by route rank so bands are legible
        for r in others:
            off = (r["route_rank"] - (len(ROUTES) - 1) / 2) * 0.115
            ax.scatter([r["log_k"]], [y + off], s=30,
                       color=RAMP[r["route_rank"]], edgecolor=SURFACE,
                       linewidth=0.7, zorder=3)
        for v in qin_by[cm]:
            ax.scatter([v], [y], s=150, marker="D", color=ORANGE,
                       edgecolor=SURFACE, linewidth=2.1, zorder=5)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(list(reversed(order)), fontsize=10.5, color=INK)
    ax.set_xlabel("log K  (L water / kg protein)", fontsize=10.5, color=MUTED)
    ax.tick_params(axis="x", colors=MUTED, labelsize=9.5)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#d8d8d2")
    ax.grid(axis="x", color="#ececE6", lw=0.8); ax.set_axisbelow(True)
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
    fig.text(0.115, 1 - 0.38 / fh, title, fontsize=13.2, color=INK, va="top")
    fig.text(0.115, 1 - 0.70 / fh,
             "Marks are shaded by how the number was obtained, darkest = "
             "measured. Rows ordered by Qin's value.",
             fontsize=9.6, color=MUTED, va="top")
    p = os.path.join(FIGS, "albumin_route_split.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE)
    print(f"  wrote {os.path.relpath(p, os.path.join(HERE,'..'))}")


def fig_b(rows, title, compounds):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    n = len(compounds)
    fig, axes = plt.subplots(1, n, figsize=(2.55 * n, 4.5), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(top=0.70, bottom=0.165, left=0.075, right=0.985,
                        wspace=0.13)

    for ax, cm in zip(axes, compounds):
        ax.set_facecolor(SURFACE)
        g = [r for r in rows if r["pfas"] == cm and r["ratio"]
             and r["ratio"] > 0]
        ax.axvspan(1e-6, PHYSIO, color="#1baf7a", alpha=0.085, zorder=0)
        for r in g:
            ax.scatter([r["ratio"]], [r["log_k"]], s=34,
                       color=RAMP[r["route_rank"]], edgecolor=SURFACE,
                       linewidth=0.7, zorder=3)
        ax.set_xscale("log")
        ax.set_title(f"{cm}   (n={len(g)})", fontsize=10.5, color=INK,
                     pad=7)
        ax.set_xlabel("ligand : protein", fontsize=9.6, color=MUTED)
        ax.tick_params(colors=MUTED, labelsize=8.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("bottom", "left"):
            ax.spines[s].set_color("#d8d8d2")
        ax.grid(color="#ececE6", lw=0.7); ax.set_axisbelow(True)
        ax.set_xlim(1e-5, 1e5)
    axes[0].set_ylabel("log K  (L water / kg protein)", fontsize=10,
                       color=MUTED)

    h = [Line2D([], [], marker="s", color="none", markerfacecolor="#1baf7a",
                alpha=0.5, markeredgecolor="none", markersize=10,
                label="physiological loading (< 0.01)")]
    h += [Line2D([], [], marker="o", color="none", markerfacecolor=RAMP[i],
                 markeredgecolor=SURFACE, markersize=8, label=nm)
          for i, (nm, _) in enumerate(ROUTES)]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.075, 0.875),
               ncol=3, frameon=False, fontsize=8.8, labelcolor=INK,
               handletextpad=0.4, columnspacing=1.4)
    fig.text(0.075, 0.965, title, fontsize=13.2, color=INK, va="top")
    fig.text(0.075, 0.918,
             "Every albumin coefficient for these five, against the molar "
             "loading it was measured at. If loading did not matter these "
             "panels would be flat.",
             fontsize=9.5, color=MUTED, va="top")
    p = os.path.join(FIGS, "albumin_logk_vs_loading.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE)
    print(f"  wrote {os.path.relpath(p, os.path.join(HERE,'..'))}")


def main():
    os.makedirs(FIGS, exist_ok=True)
    rows = load()
    with_ratio = [r for r in rows if r["ratio"]]
    top5 = [c for c, _ in collections.Counter(
        r["pfas"] for r in with_ratio).most_common(5)]
    print(f"  {len(rows)} rows, {len(with_ratio)} with a ratio")
    print(f"  five best-covered: {', '.join(top5)}\n")
    print(f"  {'route':48}{'n':>5}{'log K range':>18}")
    for i, (nm, _) in enumerate(ROUTES):
        g = [r["log_k"] for r in rows if r["route_rank"] == i]
        if g:
            print(f"  {nm:48}{len(g):>5}{f'{min(g):.2f} to {max(g):.2f}':>18}")

    for name, data in (("albumin_route_split_data.csv", rows),
                       ("albumin_logk_vs_loading_data.csv", with_ratio)):
        p = os.path.join(FIGS, name)
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
            w.writeheader(); w.writerows(data)

    fig_a(rows, sys.argv[1] if len(sys.argv) > 1 else
          "Where Qin 2024 sits, and what kind of number it is sitting among")
    fig_b(with_ratio, sys.argv[2] if len(sys.argv) > 2 else
          "A binding constant is a function of the loading it was measured at",
          top5)
    return 0


if __name__ == "__main__":
    sys.exit(main())
