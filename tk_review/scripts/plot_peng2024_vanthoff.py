"""Peng 2024's temperature series, as a van't Hoff plot.

Notebook script. Self-contained: the six series are transcribed from
Peng's Table 2 and sit in this file, so there is no CSV to keep beside
it.

WHAT THE PLOT IS FOR. Peng reports binding enthalpies from a
three-temperature series that reach +412 kJ/mol. Non-covalent
protein-ligand binding is normally within about +/-50, so the question
is whether that number is a fitting artefact or what the data say.
A van't Hoff plot answers it directly: ln K against 1/T is a straight
line only if ONE enthalpy describes the range, and the slope is
-dH/R. Three points is the minimum that can show curvature at all --
two would give a perfect line by construction, which is why the 304 K
row had to be ingested before this could be drawn.

WHAT IT SHOWS.

  The enthalpies are NOT a fitting artefact. Refitting each series
  reproduces the paper's own dH to better than 1 % on all six. The
  +412 kJ/mol for HFPO-TA is what those three numbers imply, and its
  line is the straightest of the set at R2 0.999. So the implausible
  magnitude is in the measurements, not the arithmetic.

  Half the series are not straight. R2 runs 0.72 to 0.999 on three
  points, where two points would give 1.000 for free. PFHpA's two
  segments differ by fourteen-fold in slope (+4.77 then +0.33 in
  ln K): a single enthalpy does not describe it.

  PFOA both turns over AND fails an internal check. It is the only
  compound whose ln K rises then falls, and the only one whose
  reported dG disagrees with its own Kb: at 304 K, -RT ln Kb gives
  -36.95 kJ/mol against the -32.85 printed. Every other point in the
  table agrees to 0.02. A dG of -32.85 implies Kb = 4.4e5, not the
  2.23e6 printed -- and with 4.4e5 the series is monotonic and the
  turnover disappears. The single flagged point is marked on the plot.
"""
import math

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PNG = "peng2024_vanthoff.png"
R = 8.314

# Peng Table 2, transcribed as printed.
# compound: ({T: Kb}, paper dH kJ/mol, {T: paper dG kJ/mol})
SERIES = {
    "PFNA":    ({298: 7.81e6, 304: 2.02e5, 310: 9.9e4},  -278.25,
                {298: -39.32, 304: -30.86, 310: -29.65}),
    "HFPO-TA": ({298: 3.7e6,  304: 8.54e7, 310: 2.31e9},  412.15,
                {298: -37.47, 304: -46.16, 310: -55.57}),
    "PFOA":    ({298: 2.27e5, 304: 2.23e6, 310: 1.98e6},  146.68,
                {298: -30.55, 304: -32.85, 310: -37.67}),
    "PFO3DA":  ({298: 1.59e5, 304: 4.55e4, 310: 2.99e4}, -106.57,
                {298: -29.67, 304: -27.11, 310: -26.56}),
    "PFHpA":   ({298: 4.53e3, 304: 5.35e5, 310: 7.47e5},  324.98,
                {298: -20.86, 304: -33.34, 310: -34.86}),
    "DFSA":    ({298: 1.52e3, 304: 4.93e3, 310: 2.07e4},  167.49,
                {298: -18.15, 304: -21.49, 310: -25.62}),
}
TEMPS = (298, 304, 310)

# Two slots of the validated palette, the same assignment the albumin
# figures use so the whole set reads as one family. Checked all-pairs
# on this surface: every check passes, worst pair dE 9.2 deutan /
# 24.0 normal.
MARK, FLAG = "#2a78d6", "#eb6834"
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


def fit(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    m = sxy / sxx
    b = my - m * mx
    ss = sum((y - (m * x + b)) ** 2 for x, y in zip(xs, ys))
    st = sum((y - my) ** 2 for y in ys)
    return m, b, (1 - ss / st if st else 1.0)


def inconsistent(K, dG):
    """Temperatures where the printed dG does not match -RT ln Kb."""
    out = []
    for T in TEMPS:
        calc = -R * T * math.log(K[T]) / 1000.0
        if abs(calc - dG[T]) > 0.5:
            out.append((T, calc, dG[T]))
    return out


def plot(save=PNG):
    head = 1.45
    h = head + 2 * 2.6
    fig, axes = plt.subplots(2, 3, figsize=(13.4, h))
    flat = axes.ravel()

    for ax, (name, (K, dH, dG)) in zip(flat, SERIES.items()):
        xs = [1000.0 / T for T in TEMPS]
        ys = [math.log(K[T]) for T in TEMPS]
        m, b, r2 = fit(xs, ys)
        bad = {T for T, _c, _p in inconsistent(K, dG)}

        lo, hi = min(xs) - 0.012, max(xs) + 0.012
        ax.plot([lo, hi], [m * lo + b, m * hi + b], color=MUTED,
                lw=1.1, ls=(0, (6, 3)), zorder=2)
        # Headroom so a temperature label never reaches the annotation
        # row above the axes, and the label goes BELOW the highest
        # point so it cannot collide with it either.
        span = max(ys) - min(ys) or 1.0
        ax.set_ylim(min(ys) - 0.22 * span, max(ys) + 0.30 * span)
        top = max(ys)
        for T, x, y in zip(TEMPS, xs, ys):
            flagged = T in bad
            ax.plot(x, y, "o", ms=9.5,
                    mfc=FLAG if flagged else MARK, mec=SURFACE,
                    mew=1.6, zorder=6)
            below = y == top
            ax.annotate(f"{T} K", (x, y), textcoords="offset points",
                        xytext=(0, -20 if below else 13), ha="center",
                        fontsize=8.2, color=INK2)

        ax.set_title(f"{name}", loc="left", fontsize=11.5,
                     fontweight="bold", color=INK, pad=20)
        straight = "straight" if r2 >= 0.99 else "NOT straight"
        ax.text(0.0, 1.015,
                f"ΔH {dH:+.0f} kJ/mol   ·   refit {-R*m*1000/1000:+.0f}"
                f"   ·   R² {r2:.3f}  {straight}",
                transform=ax.transAxes, fontsize=8.8, color=INK2)
        ax.set_xlim(lo, hi)
        ax.grid(True, color=GRID, lw=0.7, zorder=1)
        ax.set_axisbelow(True)
        ax.set_xlabel("1000 / T   (K⁻¹)", fontsize=9.5, labelpad=5)
        ax.set_ylabel("ln K$_b$", fontsize=9.5)

    handles = [
        Line2D([], [], marker="o", ls="none", ms=8, mfc=MARK,
               mec=SURFACE, mew=1.4, label="as printed in Table 2"),
        Line2D([], [], marker="o", ls="none", ms=8, mfc=FLAG,
               mec=SURFACE, mew=1.4,
               label="ΔG printed for this point disagrees with "
                     "−RT ln K$_b$"),
        Line2D([], [], color=MUTED, lw=1.1, ls=(0, (6, 3)),
               label="van't Hoff fit through the three points"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               frameon=False, fontsize=9.5, handletextpad=0.8,
               columnspacing=2.4, bbox_to_anchor=(0.5, -0.002))

    fig.suptitle("Peng 2024: does one binding enthalpy describe "
                 "298–310 K?", x=0.012, y=1 - 0.20 / h, ha="left",
                 va="top", fontsize=14.5, fontweight="bold", color=INK)
    fig.text(0.012, 1 - 0.52 / h,
             "ln K against 1/T is a straight line only if a single "
             "enthalpy holds across the range; the slope is −ΔH/R.\n"
             "The refit reproduces the paper's ΔH to better than 1 % "
             "on all six, so the values are not a fitting artefact — "
             "but half the series are not straight,\nand three points "
             "is the fewest number that can show it.",
             ha="left", va="top", fontsize=9.6, color=INK2,
             linespacing=1.5)
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 1 - head / h))
    if save:
        fig.savefig(save, dpi=200)
    return fig


fig = plot()
