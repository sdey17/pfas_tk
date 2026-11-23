#!/usr/bin/env python3
"""Figure: the tubular factor T against chain length, and what is undetermined.

    python3 -I scripts/stage6_figure_T_vs_chain.py

Writes report/fig_T_vs_chain_length.png. Idempotent.

THE DATA IS READ FROM db/tubular_factor_consolidated.csv, NOT TYPED IN HERE.
An earlier version of this figure carried its own hardcoded copy of the
numbers, and when the Fujii 2015 mouse rows were added to the pipeline the
chart silently kept drawing the old ones -- the mouse series simply stopped at
C10 while the data went to C12, and nothing complained. A figure that holds
its own copy of the data WILL drift from it. This one cannot.

TWO DIFFERENT SPREADS, AND THEY MUST NOT BE CONFLATED. An earlier version of
this figure drew one bar per cell and captioned it as the disagreement between
laboratories on the unbound fraction. That caption did not describe the data:
every cell in the consolidated file carries exactly ONE fu, so its T_low..T_high
spans that species' CLEARANCE sources, not its binding sources. The two spreads
are plotted separately here:

  COLOURED BAR   the range across that species' own clearance sources, both
                 sexes pooled, under a single unbound fraction.
  GREY WHISKER   human only: the range of T obtained by swapping in each of
                 the four independent human fu datasets (Smeltz, Ryu, Fischer,
                 Argoul) while holding the clearance fixed. Where it crosses
                 T = 1 the SIGN of the human tubular term is decided by which
                 binding laboratory one believes, which is the central open
                 question of the fu strand and must not disappear from the
                 figure merely because the consolidated file fixes one fu.

The dashed line at T = 1 is filtration alone. A range that lies wholly below it
means net reabsorption is established; one that crosses it means THE SIGN IS
NOT DETERMINED BY THE DATA, which is the chart's main message for the short
chains.

Palette: categorical slots 1-3 of the reference instance, validated with
scripts/validate_palette.js --pairs all --mode light. All checks pass: worst
all-pairs CVD dE 9.2 (deutan), worst normal-vision dE 24.0. Slot 3 (aqua)
draws a contrast WARN against the light surface at 2.74:1, which obliges
relief rather than dismissal: every species carries a direct text label in its
own colour, and the underlying table ships as the CSV this script reads. Light
mode only, deliberately -- this renders into a light-background PDF, so a
selected dark mode would have no surface to be validated against.
"""
import csv
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent.parent

HUMAN, MOUSE, RAT = "#2a78d6", "#eb6834", "#1baf7a"
INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"

COL = {"human": HUMAN, "mouse": MOUSE, "rat": RAT}
OFF = {"human": -0.17, "mouse": 0.0, "rat": 0.17}
YLO, YHI = 2e-3, 8.0


def load():
    """Collapse the consolidated rows to one range per species-compound.

    A cell's low end is the smallest T_low among rows that HAVE one (the
    direct measurements); its high end is the largest T_high. `bounded` is
    True only when no row in the cell carries a lower end at all, i.e. the
    whole cell rests on upper bounds and the true T lies at or below.
    """
    cells = {}
    f = HERE / "db" / "tubular_factor_consolidated.csv"
    for r in csv.DictReader(f.open()):
        sp = r["species"]
        if sp not in COL or not r["T_high"]:
            continue
        k = (r["family"], int(r["n_carbons"]), r["chemical"], sp)
        c = cells.setdefault(k, dict(lows=[], highs=[]))
        if r["T_low"]:
            c["lows"].append(float(r["T_low"]))
        c["highs"].append(float(r["T_high"]))
    out = []
    for (fam, nc, chem, sp), c in cells.items():
        hi = max(c["highs"])
        lo = min(c["lows"]) if c["lows"] else hi
        out.append((fam, nc, chem, sp, lo, hi, not c["lows"]))
    return sorted(out)


def load_fu_spread():
    """Human only: min/max T across the four independent fu datasets.

    Read from db/human_tubular_sign_by_fu_source.csv, which holds one row per
    (compound, fu source) with the clearance held fixed.
    """
    f = HERE / "db" / "human_tubular_sign_by_fu_source.csv"
    if not f.exists():
        return {}
    by = {}
    for r in csv.DictReader(f.open()):
        if not r["T"]:
            continue
        by.setdefault(r["chemical"], []).append(float(r["T"]))
    return {k: (min(v), max(v)) for k, v in by.items() if len(v) > 1}


def panel(ax, data, fam, title, xs, fuspread):
    ax.set_facecolor(SURFACE)
    ax.set_yscale("log")
    ax.axhline(1.0, color=AXIS, lw=1.6, ls=(0, (5, 3)), zorder=1)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    ax.grid(axis="x", visible=False)

    # the fu-source spread first, so the coloured marks sit on top of it
    for f, nc, lab, sp, lo, hi, bnd in data:
        if f != fam or sp != "human" or lab not in fuspread:
            continue
        flo, fhi = fuspread[lab]
        x = nc + OFF["human"]
        a, b = max(flo, YLO), min(fhi, YHI)
        if b <= a:
            continue
        ax.plot([x, x], [a, b], color=MUTED, lw=1.2, zorder=2,
                solid_capstyle="butt")
        for y, real in ((a, flo >= YLO), (b, fhi <= YHI)):
            if real:
                ax.plot([x - 0.07, x + 0.07], [y, y], color=MUTED, lw=1.2,
                        zorder=2)
        if fhi > YHI:
            ax.plot([x], [YHI], "^", ms=6, mfc="none", mec=MUTED, mew=1.2,
                    zorder=2)

    for f, nc, lab, sp, lo, hi, bnd in data:
        if f != fam:
            continue
        x = nc + OFF[sp]
        c = COL[sp]
        # a value off the top of the axis is drawn to the edge with an open
        # arrow, never silently clipped to look like a finite endpoint
        shown_hi, over = min(hi, YHI), hi > YHI
        if shown_hi > lo:
            ax.plot([x, x], [lo, shown_hi], color=c, lw=2.0,
                    solid_capstyle="butt", zorder=3)
        mk = "v" if bnd else "o"
        if over:
            ax.plot([x], [YHI], "^", ms=7, mfc="none", mec=c, mew=1.6,
                    zorder=4)
            ends = {lo} if shown_hi > lo else set()
        else:
            ends = {lo, shown_hi} if shown_hi > lo else {lo}
        for y in ends:
            ax.plot([x], [y], mk, ms=8, mfc=c, mec=SURFACE, mew=1.4,
                    zorder=4)

    ax.set_xticks(xs)
    ax.set_xticklabels([str(v) for v in xs])
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.set_ylim(YLO, YHI)
    ax.set_title(title, color=INK, fontsize=11.5, loc="left", pad=8,
                 fontweight="semibold")
    ax.set_xlabel("perfluorocarbon chain length (carbons)", color=SECOND,
                  fontsize=9.5, labelpad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)


def main():
    data = load()
    fuspread = load_fu_spread()
    skipped = sorted({(d[2], d[0]) for d in data
                      if d[0] not in ("carboxylate", "sulfonate")})
    fig, (a1, a2) = plt.subplots(
        1, 2, figsize=(11.2, 6.1), sharey=True,
        gridspec_kw=dict(width_ratios=[9, 5], wspace=0.07))
    fig.patch.set_facecolor(SURFACE)

    carb = [d for d in data if d[0] == "carboxylate"]
    sulf = [d for d in data if d[0] == "sulfonate"]
    xs1 = list(range(min(d[1] for d in carb), max(d[1] for d in carb) + 1))
    xs2 = sorted({d[1] for d in sulf})
    panel(a1, data, "carboxylate", "Carboxylates", xs1, fuspread)
    panel(a2, data, "sulfonate", "Sulfonates", xs2, fuspread)
    a2.set_ylabel("")

    a1.set_ylabel("tubular factor  T  =  CL$_R$ / (f$_u$ × GFR)",
                  color=SECOND, fontsize=10, labelpad=8)

    a1.annotate("T = 1   filtration alone", xy=(3.28, 1.0),
                xytext=(3.28, 1.16), color=SECOND, fontsize=9, va="bottom")

    # Direct series labels, so identity is never colour-alone -- and the
    # relief the aqua contrast WARN obliges. Each sits beside a REAL mark of
    # that species: human PFPeA (C5), mouse PFBA (C4), rat PFNA (C9).
    for sp, txt, xy in (("human", "human", (4.95, 0.74)),
                        ("mouse", "mouse", (4.13, 0.050)),
                        ("rat", "rat", (9.34, 2.45))):
        a1.annotate(txt, xy=xy, color=COL[sp], fontsize=10.5,
                    fontweight="bold")

    handles = [
        Line2D([], [], color=COL[s], lw=2, marker="o", ms=8, mec=SURFACE,
               mew=1.4, label=s) for s in ("human", "rat", "mouse")
    ] + [
        Line2D([], [], color=MUTED, lw=0, marker="o", ms=8, mec=SURFACE,
               mew=1.4, label="direct renal clearance"),
        Line2D([], [], color=MUTED, lw=0, marker="v", ms=8, mec=SURFACE,
               mew=1.4, label="upper bound (true T at or below)"),
        Line2D([], [], color=MUTED, lw=0, marker="^", ms=7, mfc="none",
               mew=1.6, label="value runs off the top of the axis"),
        Line2D([], [], color=MUTED, lw=1.2, marker="_", ms=9,
               label="human: spread across the four f$_u$ datasets"),
    ]
    # Placed below the panels, not inside them: in the sulfonate panel the
    # rat PFHxS bar ran straight through the legend text.
    leg = fig.legend(handles=handles, loc="lower center", frameon=False,
                     fontsize=9, ncol=3, handletextpad=0.6,
                     columnspacing=2.2, bbox_to_anchor=(0.52, 0.100))
    for t in leg.get_texts():
        t.set_color(SECOND)

    fig.suptitle("Net tubular handling of PFAS in three species: what the "
                 "clearances settle, and what the binding disagreement "
                 "leaves open",
                 color=INK, fontsize=12.5, x=0.012, ha="left", y=0.985,
                 fontweight="semibold")
    note = ("Below the dashed line the kidney returns filtered PFAS to the "
            "blood; above it, the kidney secretes. A coloured bar spans that "
            "species' own CLEARANCE sources under a single\nunbound "
            "fraction; the grey whisker spans the four human f$_u$ datasets "
            "with the clearance held fixed. Where either crosses T = 1 the "
            "sign is not determined by\nthe data. The renal clearances "
            "agree across laboratories; the binding measurements do not — "
            "which is why the grey whiskers are the wider ones."
            "\nDrawn from db/tubular_factor_consolidated.csv and "
            "db/human_tubular_sign_by_fu_source.csv.")
    if skipped:
        note += ("  Not shown (neither family): "
                 + ", ".join(c for c, _ in skipped) + ".")
    fig.text(0.012, 0.010, note, color=MUTED, fontsize=8.2, va="bottom")

    fig.subplots_adjust(left=0.085, right=0.985, top=0.895, bottom=0.355)
    out = HERE / "report" / "fig_T_vs_chain_length.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    n = len(data)
    print(f"  {n} species-compound cells drawn "
          f"({len(carb)} carboxylate, {len(sulf)} sulfonate)")
    for sp in ("human", "rat", "mouse"):
        d = [x for x in data if x[3] == sp]
        print(f"    {sp:7} {len(d):>2} cells, C{min(x[1] for x in d)}"
              f"-C{max(x[1] for x in d)}")
    print(f"wrote {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
