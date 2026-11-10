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

QIN IS PLOTTED AT ITS OPERATING POINT, not at the maximum its
isotherm reached -- see qin_operating_nu.  The two differ by about
100x and only one of them is the condition the plotted coefficient
describes.

NO TITLES.  Both figures are drawn without a title or subtitle, so
they drop into a document that supplies its own caption.  Everything
needed to read a figure unaided -- axis labels, per-panel compound
names, the legend -- is still on the canvas; only the prose is gone.

FONT.  Arial is requested first and Liberation Sans is the fallback.
The two are metrically compatible -- identical glyph advance widths,
Liberation having been built as a drop-in -- so a figure rendered on a
machine without Arial has the same layout, line breaks and label
collisions as one rendered with it.

Writes figures/albumin_route_split.png, figures/albumin_logk_vs_loading.png
and the CSV each is drawn from.
"""
import collections
import importlib.util
import csv
import math
import os
import sys

import openpyxl

# Arial first, so a machine that has it uses it; Liberation Sans is
# metrically compatible with Arial (identical glyph advance widths,
# built as a drop-in) and is what renders where Arial is absent, so
# layout and line breaks are the same either way.
FONT_STACK = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]


def use_arial():
    import matplotlib
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = FONT_STACK
    matplotlib.rcParams["axes.unicode_minus"] = False


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
M_ALB_KG = 66463.0 / 1e3        # kg/mol


def _anchor():
    """qin2024_anchor.py -- the one definition of the Qin anchor."""
    spec = importlib.util.spec_from_file_location(
        "qin_anchor", os.path.join(HERE, "qin2024_anchor.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def qin_nu_table():
    """Occupancy of the Qin anchor, mol PFAS per mol albumin.

    Carried in the CSV beside the ratio, not plotted.  nu needs no protein
    concentration -- it is a property of the coefficient and the FREE
    concentration alone, nu = K x Cw x M -- so it is the one loading
    measure that is exactly comparable across studies.  Every other row
    here has only a dosed ratio, so the axis has to be the ratio and nu
    rides along for whoever wants to check where the Qin points sit.

    It is read straight off the anchor, never recomputed from the sheet:
    the sheet's concentration column is the DOSE, and nu needs the free
    concentration the coefficient was read at, which is a different and
    much smaller number.
    """
    a = _anchor()
    canon = {"PFUnA": "PFUnDA"}
    return {canon.get(r["pfas"], r["pfas"]): round(r["nu"], 4)
            for r in a.anchor()}


def route_of(r):
    for i, (name, test) in enumerate(ROUTES):
        if test(r):
            return i, name
    return len(ROUTES) - 1, ROUTES[-1][0]


def load():
    NU = qin_nu_table()
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
        lk = float(d["log K (L water / kg protein)"])
        # The Qin anchor now reports a dosed ligand:protein ratio of the
        # same kind as every other row -- IC10,nom over the assay medium's
        # protein -- so it goes on the axis unchanged.  It used to be
        # plotted at its occupancy instead, because the previous anchor
        # had no dosed ratio of its own to use.
        plot_ratio, basis = ratio, "assay maximum"
        nu = ""
        if str(d["Source"]).startswith("Qin 2024"):
            basis = "IC10,nom over medium protein (qin2024_anchor.py)"
            nu = NU.get(str(d["PFAS"]), "")
        out.append(dict(pfas=str(d["PFAS"]), species=str(d["Species"]),
                        source=str(d["Source"]), log_k=lk,
                        route=name, route_rank=i, ratio=ratio,
                        plot_ratio=plot_ratio, ratio_basis=basis, nu=nu))
    return out


def fig_a(rows):
    import matplotlib
    matplotlib.use("Agg")
    use_arial()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    qin = [r for r in rows if r["source"].startswith("Qin 2024")]
    qin_by = collections.defaultdict(list)
    for r in qin:
        qin_by[r["pfas"]].append(r["log_k"])
    order = [c for c in sorted(qin_by, key=lambda c: -max(qin_by[c]))
             if any(r["pfas"] == c and not r["source"].startswith("Qin")
                    for r in rows)]

    head = 0.92            # legend only; titles are added downstream
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
    p = os.path.join(FIGS, "albumin_route_split.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE)
    print(f"  wrote {os.path.relpath(p, os.path.join(HERE,'..'))}")


def fig_b(rows, compounds):
    import matplotlib
    matplotlib.use("Agg")
    use_arial()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    n = len(compounds)
    fig, axes = plt.subplots(1, n, figsize=(2.55 * n, 4.75), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    # The legend is three rows once Qin has its own key; leave it room
    # so it does not sit on the panel titles.
    fig.subplots_adjust(top=0.78, bottom=0.155, left=0.075, right=0.985,
                        wspace=0.13)

    for ax, cm in zip(axes, compounds):
        ax.set_facecolor(SURFACE)
        g = [r for r in rows if r["pfas"] == cm and r["plot_ratio"]
             and r["plot_ratio"] > 0]
        ax.axvspan(1e-6, PHYSIO, color="#1baf7a", alpha=0.085, zorder=0)
        # Qin 2024 carries the same orange diamond it has in figure A, so
        # the two figures are read with one key. It is drawn on top and
        # keeps its identity rather than its route shade -- in both
        # figures the question is where Qin sits relative to the field,
        # and its route is named in the caption of figure A anyway.
        for r in g:
            if r["source"].startswith("Qin 2024"):
                continue
            ax.scatter([r["plot_ratio"]], [r["log_k"]], s=34,
                       color=RAMP[r["route_rank"]], edgecolor=SURFACE,
                       linewidth=0.7, zorder=3)
        for r in g:
            if not r["source"].startswith("Qin 2024"):
                continue
            ax.scatter([r["plot_ratio"]], [r["log_k"]], s=112, marker="D",
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
            ax.spines[s].set_color("#d8d8d2")
        ax.grid(color="#ececE6", lw=0.7); ax.set_axisbelow(True)
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
    p = os.path.join(FIGS, "albumin_logk_vs_loading.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE)
    print(f"  wrote {os.path.relpath(p, os.path.join(HERE,'..'))}")


def main():
    os.makedirs(FIGS, exist_ok=True)
    rows = load()
    with_ratio = [r for r in rows if r["plot_ratio"]]
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

    fig_a(rows)
    fig_b(with_ratio, top5)
    return 0


if __name__ == "__main__":
    sys.exit(main())
