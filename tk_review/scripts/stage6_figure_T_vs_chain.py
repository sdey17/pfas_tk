#!/usr/bin/env python3
"""Figure: the tubular factor T against chain length, and what is undetermined.

    python3 -I scripts/stage6_figure_T_vs_chain.py

Writes report/fig_T_vs_chain_length.png. Idempotent.

WHY A RANGE AND NOT A POINT. For every species-compound cell there is more than
one defensible value of T, because more than one laboratory measured the
unbound fraction and T scales inversely with it. Plotting one number would
hide that. Each cell is therefore drawn as the RANGE across the available
within-species sources, with the marker style saying whether the endpoints are
direct measurements or one-sided bounds.

The dashed line at T = 1 is filtration alone. A range that lies wholly below it
means net reabsorption is established; one that crosses it means THE SIGN IS
NOT DETERMINED BY THE DATA, which is the chart's main message for the short
chains.

Palette: categorical slots 1 and 2 of the reference instance, validated with
scripts/validate_palette.js --pairs all --mode light (all six checks PASS:
worst-pair CVD dE 24.7, normal-vision 33.6). Light mode only, deliberately --
this renders into a light-background PDF, so a selected dark mode has no
surface to be validated against.
"""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent.parent

HUMAN, MOUSE = "#2a78d6", "#eb6834"
INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"

# (family, nC, label, species, T_low, T_high, bounded)
# T_low/high span the within-species sources; bounded=True means the endpoints
# are upper bounds, so the true T lies at or below them.
D = [
    # ---- carboxylates, human: Ryu fu vs Abraham's own fup ----
    ("carboxylate", 4,  "PFBA",   "human", 0.0300, 0.0946, False),
    ("carboxylate", 5,  "PFPeA",  "human", 0.6840, 1.7780, False),
    ("carboxylate", 6,  "PFHxA",  "human", 0.4820, 3.4805, False),
    ("carboxylate", 7,  "PFHpA",  "human", 0.0350, 0.5111, False),
    ("carboxylate", 8,  "PFOA",   "human", 0.0066, 0.0182, True),
    ("carboxylate", 9,  "PFNA",   "human", 0.0059, 0.0089, True),
    ("carboxylate", 10, "PFDA",   "human", 0.0106, 0.0365, True),
    ("carboxylate", 11, "PFUnDA", "human", 0.0039, 0.1060, True),
    ("carboxylate", 12, "PFDoDA", "human", 1.2230, 1.2230, True),
    # ---- carboxylates, mouse: Ryu fu + OEHHA GFR vs Argoul's own fu + GFR ----
    ("carboxylate", 4,  "PFBA",   "mouse", 0.0611, 0.1017, False),
    ("carboxylate", 6,  "PFHxA",  "mouse", 2.2307, 2.7120, False),
    ("carboxylate", 8,  "PFOA",   "mouse", 0.0151, 0.0413, False),
    ("carboxylate", 9,  "PFNA",   "mouse", 0.0082, 0.0171, False),
    ("carboxylate", 10, "PFDA",   "mouse", 0.0031, 0.0061, False),
    # ---- sulfonates ----
    ("sulfonate", 4, "PFBS",  "human", 0.0128, 0.1288, False),
    ("sulfonate", 6, "PFHxS", "human", 0.0295, 0.0300, True),
    ("sulfonate", 8, "PFOS",  "human", 0.0098, 0.0454, True),
    ("sulfonate", 6, "PFHxS", "mouse", 0.0027, 0.0045, False),
    ("sulfonate", 8, "PFOS",  "mouse", 0.0290, 0.0496, False),
]
COL = {"human": HUMAN, "mouse": MOUSE}
OFF = {"human": -0.11, "mouse": 0.11}


def panel(ax, fam, title, xs):
    ax.set_facecolor(SURFACE)
    ax.set_yscale("log")
    ax.axhline(1.0, color=AXIS, lw=1.6, ls=(0, (5, 3)), zorder=1)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    ax.grid(axis="x", visible=False)

    for f, nc, lab, sp, lo, hi, bnd in D:
        if f != fam:
            continue
        x = nc + OFF[sp]
        c = COL[sp]
        if hi > lo:
            ax.plot([x, x], [lo, hi], color=c, lw=2.0, solid_capstyle="butt",
                    zorder=3)
        # endpoints: filled circle = direct, downward caret = upper bound
        mk = "v" if bnd else "o"
        for y in ({lo, hi} if hi > lo else {lo}):
            ax.plot([x], [y], mk, ms=8, mfc=c, mec=SURFACE, mew=1.4,
                    zorder=4)

    ax.set_xticks(xs)
    ax.set_xticklabels([str(v) for v in xs])
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.set_ylim(2e-3, 8.0)
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
    fig, (a1, a2) = plt.subplots(
        1, 2, figsize=(11.2, 6.1), sharey=True,
        gridspec_kw=dict(width_ratios=[9, 5], wspace=0.07))
    fig.patch.set_facecolor(SURFACE)

    panel(a1, "carboxylate", "Carboxylates", list(range(4, 13)))
    panel(a2, "sulfonate", "Sulfonates", [4, 6, 8])
    a2.set_ylabel("")

    a1.set_ylabel("tubular factor  T  =  CL$_R$ / (f$_u$ × GFR)",
                  color=SECOND, fontsize=10, labelpad=8)

    # the reference line earns a label, since it is the decision boundary
    a1.annotate("T = 1   filtration alone", xy=(8.9, 1.0),
                xytext=(8.9, 1.22), color=SECOND, fontsize=9, va="bottom")

    # direct series labels, so identity is never colour-alone
    a1.annotate("human", xy=(7 - 0.11, 0.5111), xytext=(7.22, 0.42),
                color=HUMAN, fontsize=10.5, fontweight="bold")
    a1.annotate("mouse", xy=(6 + 0.11, 2.71), xytext=(6.42, 3.6),
                color=MOUSE, fontsize=10.5, fontweight="bold")

    handles = [
        Line2D([], [], color=HUMAN, lw=2, marker="o", ms=8, mec=SURFACE,
               mew=1.4, label="human"),
        Line2D([], [], color=MOUSE, lw=2, marker="o", ms=8, mec=SURFACE,
               mew=1.4, label="mouse"),
        Line2D([], [], color=MUTED, lw=0, marker="o", ms=8, mec=SURFACE,
               mew=1.4, label="direct renal clearance"),
        Line2D([], [], color=MUTED, lw=0, marker="v", ms=8, mec=SURFACE,
               mew=1.4, label="upper bound (true T at or below)"),
    ]
    leg = a2.legend(handles=handles, loc="upper left", frameon=False,
                    fontsize=9, labelcolor=SECOND, handletextpad=0.6,
                    borderaxespad=0.4, labelspacing=0.5)
    for t in leg.get_texts():
        t.set_color(SECOND)

    fig.suptitle("Net tubular handling of PFAS, human and mouse: a bar spans "
                 "the range across that species' own sources",
                 color=INK, fontsize=12.5, x=0.012, ha="left", y=0.985,
                 fontweight="semibold")
    fig.text(0.012, 0.015,
             "Below the dashed line the kidney returns filtered PFAS to "
             "the blood; above it, the kidney secretes. A bar wholly below "
             "establishes net reabsorption;\na bar crossing it means the "
             "sign is not determined — the spread is the disagreement "
             "between laboratories on the unbound fraction, which the "
             "filtered\nload scales with. The renal clearances agree across "
             "laboratories; the binding measurements do not.",
             color=MUTED, fontsize=8.6, va="bottom")

    fig.subplots_adjust(left=0.085, right=0.985, top=0.895, bottom=0.235)
    out = HERE / "report" / "fig_T_vs_chain_length.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
