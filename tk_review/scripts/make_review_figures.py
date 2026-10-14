"""
Figures 2 and 3 for report/PFAS_TK_review.pdf.

Each panel carries an argument the prose otherwise has to assert. Palette and
mark specs follow the project's existing figure style, validated for
colour-vision deficiency: worst adjacent pair deltaE 9.2 (deutan), normal-vision
floor 27.1. The aqua fails the 3:1 contrast check against the surface, so every
series is directly labelled as well as coloured -- identity is never carried by
colour alone.

Run:  python3 scripts/make_review_figures.py
"""

import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS = os.path.join(HERE, "figures")

BLUE, ORANGE, AQUA, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#7b5ea7"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID, RULE = "#e6e5e1", "#d8d7d3"

# One colour per compound, shared across the figures that name compounds.
COMPOUND = {"PFOA": BLUE, "PFHxS": ORANGE, "PFOS": AQUA,
            "PFNA": PURPLE, "PFBA": INK2, "PFHxA": MUTED}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": RULE, "axes.linewidth": 0.8, "axes.labelcolor": INK2,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "figure.dpi": 150,
})


def rows(p):
    return list(csv.DictReader(open(os.path.join(HERE, p), newline="")))


# ---------------------------------------------------------------- figure 2
def figure2():
    """A: both ratios in the SAME direction, so the asymmetry is not an
    artefact of how each was written. B: the species gap split into its two
    factors."""
    d = [r for r in rows("db/sex_decomposition.csv")
         if r["vd_ratio_F_over_M"] and r["clearance_ratio_F_over_M"]]
    sp_colour = {"rat": BLUE, "mouse": ORANGE, "primate": PURPLE}

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.6, 4.5), gridspec_kw=dict(width_ratios=[1.45, 1]))

    # ---- panel A: strip plot of the two ratio distributions ----
    for i, (key, lab) in enumerate([("vd_ratio_F_over_M", "volume of distribution"),
                                    ("clearance_ratio_F_over_M", "clearance")]):
        vals = [float(r[key]) for r in d]
        cols = [sp_colour[r["species"]] for r in d]
        jitter = np.linspace(-0.16, 0.16, len(vals))
        axA.scatter(vals, [i] * len(vals) + jitter, s=60, c=cols, zorder=4,
                    edgecolor=SURFACE, linewidth=1.4)
        lo, hi = min(vals), max(vals)
        axA.plot([lo, hi], [i - 0.30, i - 0.30], color=MUTED, lw=1.4, zorder=2)
        axA.text(math.sqrt(lo * hi), i - 0.47,
                 f"spans {hi/lo:.0f}×" if hi / lo >= 10
                 else f"spans {hi/lo:.1f}×",
                 ha="center", va="top", fontsize=9, color=INK,
                 fontweight="bold")

    axA.axvline(1.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    axA.text(1.0, -0.70, "1.0 = no sex difference", fontsize=7.5,
             color=MUTED, ha="center", va="bottom")
    axA.set_xscale("log")
    axA.set_yticks([0, 1])
    axA.set_yticklabels(["volume of\ndistribution", "clearance"], fontsize=9.5)
    axA.set_ylim(-0.75, 1.6)
    axA.set_xlabel("female / male ratio   (log scale, both measures same direction)")
    axA.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    axA.set_axisbelow(True)
    axA.set_title("A · Only one of the two terms varies with sex",
                  loc="left", fontsize=10.5, color=INK, fontweight="bold",
                  pad=10)
    for sp, col in sp_colour.items():
        axA.scatter([], [], s=60, c=col, label=sp, edgecolor=SURFACE,
                    linewidth=1.4)
    axA.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=8,
               ncol=3, handletextpad=0.3, columnspacing=1.1)

    # ---- panel B: the 333x gap decomposed ----
    total, escape, gfr = 333.0, 50.0, 6.5
    le, lg = math.log10(escape), math.log10(gfr)
    axB.barh([0], [le], color=BLUE, height=0.42, zorder=3)
    axB.barh([0], [lg], left=[le + 0.012], color=ORANGE, height=0.42, zorder=3)
    axB.text(le / 2, 0, f"escape fraction\n{escape:.0f}×",
             ha="center", va="center", color="white", fontsize=9,
             fontweight="bold", zorder=5)
    axB.text(le + lg / 2 + 0.012, 0, f"GFR\n{gfr:.1f}×", ha="center",
             va="center", color="white", fontsize=9, fontweight="bold",
             zorder=5)
    axB.text(le / 2, 0.30, f"{100*le/(le+lg):.0f}% of the gap", ha="center",
             fontsize=8.5, color=INK2)
    axB.text(le + lg / 2, 0.30, f"{100*lg/(le+lg):.0f}%", ha="center",
             fontsize=8.5, color=INK2)
    axB.set_xlim(0, (le + lg) * 1.06)
    axB.set_ylim(-0.75, 0.75)
    axB.set_yticks([])
    axB.set_xlabel("log$_{10}$ of the fold-contribution")
    axB.set_title("B · The axis is two factors, not one",
                  loc="left", fontsize=10.5, color=INK, fontweight="bold",
                  pad=10)
    axB.text(0, -0.52,
             f"male mouse → human renal clearance gap: {total:.0f}× total.\n"
             "A third of it is filtration rate, which no transporter\n"
             "mechanism explains.",
             fontsize=8.3, color=INK2, va="top")
    axB.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    axB.set_axisbelow(True)

    fig.suptitle("Where the sex and species differences actually sit",
                 x=0.012, y=0.985, ha="left", fontsize=11.5, color=INK,
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    out = os.path.join(FIGS, "fig02_sex_and_species_decomposition.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    vd = [float(r["vd_ratio_F_over_M"]) for r in d]
    cl = [float(r["clearance_ratio_F_over_M"]) for r in d]
    print(f"  fig02: {len(d)} pairs; Vd spans {max(vd)/min(vd):.1f}x, "
          f"CL spans {max(cl)/min(cl):.0f}x")
    return max(vd) / min(vd), max(cl) / min(cl), len(d)


# ---------------------------------------------------------------- figure 3
def figure3():
    """A: the structure-activity endpoint by compound. B: how far the whole
    endpoint moves on one contested input."""
    q = [r for r in rows("db/qsar/qsar_endpoint_table.csv")
         if r["log10R_gfr_argoul"]]
    q.sort(key=lambda r: float(r["log10R_gfr_argoul"]))

    head_colour = {"carboxylate": BLUE, "sulfonate": ORANGE,
                   "ether-carboxylate": AQUA}
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.6, 4.6), gridspec_kw=dict(width_ratios=[1.5, 1]))

    y = np.arange(len(q))
    vals = [float(r["log10R_gfr_argoul"]) for r in q]
    cols = [head_colour[r["head_group"]] for r in q]
    axA.hlines(y, 0, vals, color=GRID, lw=1.6, zorder=2)
    axA.scatter(vals, y, s=85, c=cols, zorder=4, edgecolor=SURFACE,
                linewidth=1.6)
    axA.axvline(0, color=INK, lw=1.0, zorder=3)
    # Compound identity goes on the axis rather than into floating labels, so
    # nothing can collide as the values move.
    axA.set_yticks(y)
    axA.set_yticklabels([f"{r['chemical']}   C{r['n_fluorinated_c']}"
                         for r in q], fontsize=8.8)
    axA.tick_params(axis="y", length=0)
    axA.set_ylim(-0.95, len(q) - 0.1)
    axA.set_xlim(-3.15, 1.5)
    axA.set_xlabel(r"log$_{10}$ R   (R = CL$_{renal}$ / f$_u \cdot$ GFR)")
    axA.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    axA.set_axisbelow(True)
    axA.text(-0.08, -0.86, "← reabsorbed", fontsize=8, color=INK2,
             ha="right")
    axA.text(0.08, -0.86, "secreted →", fontsize=8, color=INK2)
    axA.set_title("A · Chain length alone does not order the endpoint",
                  loc="left", fontsize=10.5, color=INK, fontweight="bold",
                  pad=24)
    for h, col in head_colour.items():
        axA.scatter([], [], s=85, c=col, label=h.replace("-", " "),
                    edgecolor=SURFACE, linewidth=1.6)
    axA.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), fontsize=8.2,
               ncol=3, handletextpad=0.3, columnspacing=1.2)
    # bracket the three five-carbon compounds: same chain length, 1.4 apart
    c5 = [i for i, r in enumerate(q) if r["n_fluorinated_c"] == "5"]
    xb = 1.05
    axA.plot([xb, xb], [min(c5), max(c5)], color=MUTED, lw=1.1, zorder=3)
    for i in c5:
        axA.plot([vals[i] + 0.07, xb], [i, i], color=GRID, lw=0.9, zorder=1)
        axA.plot([xb - 0.05, xb], [i, i], color=MUTED, lw=1.1, zorder=3)
    axA.text(xb + 0.1, (min(c5) + max(c5)) / 2,
             "all C5,\n1.4 log units\napart", fontsize=7.8, color=INK2,
             va="center", ha="left")

    # ---- panel B: the fu sensitivity ----
    fus = [(0.10, "read as\n'~90% bound'", -3.93, 354),
           (0.02, "assumed\n(Han 2012)", -3.23, 71),
           (0.00061, "measured\n(Fischer)", -1.72, 2.2)]
    mouse = -1.384
    xs = np.arange(len(fus))
    axB.axhline(mouse, color=ORANGE, lw=2.0, zorder=3)
    axB.text(2.42, mouse + 0.10, "female mouse PFOA", fontsize=8.4,
             color=ORANGE, ha="right", fontweight="bold")
    for x, (fu, lab, lr, gap) in zip(xs, fus):
        axB.vlines(x, lr, mouse, color=GRID, lw=2.2, zorder=2)
        axB.scatter([x], [lr], s=95, c=BLUE, zorder=5, edgecolor=SURFACE,
                    linewidth=1.6)
        axB.text(x, lr - 0.17, f"{gap:g}×", ha="center", va="top",
                 fontsize=9.5, color=INK, fontweight="bold")
        axB.text(x, lr - 0.46, "gap to mouse", ha="center", va="top",
                 fontsize=7.4, color=MUTED)
    axB.set_xticks(xs)
    axB.set_xticklabels([f"f$_u$ = {f:g}\n{l}" for f, l, _, _ in fus],
                        fontsize=8.2)
    axB.set_xlim(-0.6, 2.6)
    axB.set_ylim(-4.6, -0.95)
    axB.set_ylabel("log$_{10}$ R, human PFOA")
    axB.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    axB.set_axisbelow(True)
    axB.set_title("B · One contested input moves the species gap 160-fold",
                  loc="left", fontsize=10.5, color=INK, fontweight="bold",
                  pad=10)

    fig.suptitle("The structure-activity endpoint, and how much it rests "
                 "on one number",
                 x=0.012, y=0.985, ha="left", fontsize=11.5, color=INK,
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    out = os.path.join(FIGS, "fig03_structure_activity_and_sensitivity.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig03: {len(q)} compounds; span "
          f"{max(vals)-min(vals):.2f} log units; fu sensitivity 2.21 log units")


# ---------------------------------------------------------------- figure 4
def figure4():
    """Argoul's cocktail against the single-compound literature. Both terms
    sit low; their ratio does not. Reads the comparison table written by
    scripts/argoul_vs_single_compound.py so no number is duplicated here."""
    q = rows("db/argoul_vs_single_compound.csv")
    gmean = lambda v: math.exp(sum(math.log(x) for x in v) / len(v))

    def pick(param):
        return [(r["chemical"], r["comparator"], float(r["ratio"]))
                for r in q if r["parameter"] == param]

    cl = pick("clearance")
    vd = pick("volume_of_distribution")
    hl = pick("half_life")

    def decade_axis(ax):
        """Log axis labelled only where we put a tick. Matplotlib's default
        minor labels (3x10^0, 4x10^0, ...) overprint each other here."""
        ax.set_xscale("log")
        ax.set_xlim(0.25, 2.6)
        ax.xaxis.set_minor_locator(mticker.NullLocator())
        ax.xaxis.set_major_locator(mticker.FixedLocator([0.3, 0.5, 1, 2]))
        ax.xaxis.set_major_formatter(
            mticker.FixedFormatter(["0.3", "0.5", "1.0", "2.0"]))
        ax.axvline(1.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
        ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)

    def gmark(ax, v, y0, y1, ytext, color=INK, note=False):
        """Geometric-mean tick, drawn clear of the data rather than on them.
        Panel A carries two of these, so its caption goes to the side once
        instead of under each."""
        g = gmean(v)
        ax.plot([g, g], [y0, y1], color=color, lw=2.4, zorder=5)
        ax.text(g, ytext, f"{g:.2f}\u00d7", ha="center", va="top",
                fontsize=9.5, color=color, fontweight="bold")
        if note:
            # Offset in points, not data units: a data-unit drop scales with
            # the panel's y range and lands differently on each.
            ax.annotate("geometric mean", xy=(g, ytext), xytext=(0, -13),
                        textcoords="offset points", ha="center", va="top",
                        fontsize=7.2, color=MUTED)
        return g

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.4, 4.0),
                                   gridspec_kw=dict(width_ratios=[1, 1]))

    # ---- A: one row per comparator study, clearance and volume side by
    #         side. Shape and colour both carry the parameter, and a joining
    #         line shows the two moving together rather than apart.
    cl_by = {(c, src): v for c, src, v in cl}
    vd_by = {(c, src): v for c, src, v in vd}
    order = [("PFHxS", "Sundstrom 2012, 20 mg/kg", "PFHxS \u00b7 Sundstr\u00f6m, 20 mg/kg"),
             ("PFHxS", "Sundstrom 2012, 1 mg/kg", "PFHxS \u00b7 Sundstr\u00f6m, 1 mg/kg"),
             ("PFOA", "Fujii 2015", "PFOA \u00b7 Fujii 2015")]
    for i, (c, src, _) in enumerate(order):
        a, b = cl_by.get((c, src)), vd_by.get((c, src))
        if a and b:
            axA.plot([a, b], [i, i], color=GRID, lw=2.4, zorder=2,
                     solid_capstyle="round")
        if a:
            axA.scatter([a], [i], s=90, c=BLUE, marker="o", zorder=4,
                        edgecolor=SURFACE, linewidth=1.5)
        if b:
            axA.scatter([b], [i], s=95, c=AQUA, marker="s", zorder=4,
                        edgecolor=SURFACE, linewidth=1.5)
    gmark(axA, [x[2] for x in cl], -0.72, -0.52, -0.76, BLUE)
    gmark(axA, [x[2] for x in vd], -1.42, -1.22, -1.46, AQUA)
    axA.text(0.262, -1.30, "geometric\nmeans", fontsize=7.2, color=MUTED,
             ha="left", va="center", linespacing=1.35)
    axA.text(1.07, 2.55, "1.0 = agrees with the\nsingle-compound study",
             fontsize=7.4, color=MUTED, ha="left", va="top", linespacing=1.35)
    for lab, col, mk in [("clearance", BLUE, "o"),
                         ("volume of distribution", AQUA, "s")]:
        axA.scatter([], [], s=90, c=col, marker=mk, label=lab,
                    edgecolor=SURFACE, linewidth=1.4)
    axA.legend(loc="center right", bbox_to_anchor=(1.0, 0.52), fontsize=8,
               handletextpad=0.35, labelspacing=0.6)
    axA.set_yticks(range(len(order)))
    axA.set_yticklabels([x[2] for x in order], fontsize=8.5)
    axA.set_ylim(-1.90, 2.75)
    decade_axis(axA)
    axA.set_xlabel("Argoul cocktail / single-compound study   (log scale)")
    axA.set_title("A \u00b7 Both terms sit two- to threefold low",
                  loc="left", fontsize=10.5, color=INK, fontweight="bold",
                  pad=20)

    # ---- B: half-life, the ratio of the two. Compound on the y-axis, as in
    #         figure 3 -- no label to collide, no colour to decode.
    order = ["PFBA", "PFHxS", "PFOS"]
    for i, c in enumerate(order):
        v = [x[2] for x in hl if x[0] == c]
        axB.scatter(v, [i] * len(v), s=85, c=COMPOUND[c], zorder=4,
                    edgecolor=SURFACE, linewidth=1.5)
    vals = [x[2] for x in hl]
    g = gmark(axB, vals, -0.86, -0.68, -0.90, note=True)
    axB.text(1.06, 2.70, f"{len(vals)} comparisons, {min(vals):.2f}\u00d7 to "
             f"{max(vals):.2f}\u00d7,\nscattered either side of 1.0",
             fontsize=7.4, color=MUTED, ha="left", va="top", linespacing=1.35)
    axB.set_yticks(range(len(order)))
    axB.set_yticklabels(order, fontsize=9.5)
    axB.set_ylim(-1.30, 2.85)
    decade_axis(axB)
    axB.set_xlabel("Argoul ln2\u00b7MRT / single-compound half-life")
    axB.set_title("B \u00b7 Their ratio does not",
                  loc="left", fontsize=10.5, color=INK, fontweight="bold",
                  pad=20)

    fig.suptitle("A shared scaling, not a disagreement about elimination: "
                 "the cocktail study against the single-compound literature",
                 x=0.012, y=0.96, ha="left", fontsize=11.2, color=INK,
                 fontweight="bold")
    fig.subplots_adjust(left=0.175, right=0.985, bottom=0.21, top=0.74,
                        wspace=0.30)
    out = os.path.join(FIGS, "fig04_argoul_vs_single_compound.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig04: CL {gmean([x[2] for x in cl]):.2f}x, "
          f"Vd {gmean([x[2] for x in vd]):.2f}x, "
          f"half-life {g:.2f}x over {len(vals)} comparisons")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    figure2()
    figure3()
    figure4()
    print("  wrote figures/fig02_*.png, fig03_*.png and fig04_*.png")
