#!/usr/bin/env python3
"""Does loading explain the albumin spread? No -- and this shows why not.

WHAT THIS FIGURE USED TO CLAIM, AND WHY THAT WAS WRONG. An earlier
version of this script fitted one line per compound through all the
points, found a negative slope close to Qin's own within-laboratory
slope, and read it as the saturation effect: load more PFAS per albumin
and the apparent coefficient falls.

That reading cannot be right. Every converted point on these panels is
K0 = 1e9*n/(M*Kd), the trace-limit coefficient. Kd is a fitted constant,
so K0 carries no concentration dependence at all -- a quantity with
saturation divided out of it by construction cannot exhibit saturation.
Qin's slope IS a saturation effect, because she measures D at each free
concentration; the pooled slope here is something else wearing the same
number.

WHAT IT ACTUALLY IS. Splitting each fit by route settles it. The pooled
slope is steep and tidy; the two within-group slopes are weak and
disagree with each other and with the pooled line on sign. The reason
is drawn on the panels as a shaded band: the measured and converted
clusters DO NOT OVERLAP on any compound, with 2.1 to 3.7 decades of
empty x between them. A line through two separated clusters measures
the step between the groups, not a trend within either. Route and
loading are collinear in this literature, so no amount of fitting
separates them.

That is the finding, and it is consistent with the arithmetic reached
another way: restricting to directly measured rows removes 61% of the
pooled spread, correcting for loading removes 14%.

WHY LOADING AND NOT CONCENTRATION on the x axis. Fischer reports a
ligand:protein ratio and no molar concentration, so a concentration
axis drops that cluster entirely. The ratio keeps both groups and is
dimensionless, so assays at different protein loads are comparable.

ALLENDORF'S X IS MEASURED, NOT ASSUMED. Its main text gives only
"molar ratio held below 0.1", but SI Table SI-2.3 reports the
analyte:BSA ratio actually reached per compound, 4e-4 to 3.8e-2. SI-2.3
covers batch A3803 only; the A7906 row inherits it, which is an
assumption the data CSV records.

QIN IS A LINE, NOT A POINT, and is drawn in neutral ink rather than a
series colour: it is a reference value, not one of the two groups. Its
number is the Table S6 fit evaluated at this study's operating point,
so it has no honest x position.

Writes figures/albumin_loading_vs_logk.png and its data CSV. Given a
CSV as an argument it re-renders from that instead of the workbook.

Run:  python3 tk_review/scripts/make_conc_vs_logk_figure.py [data.csv]
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
OUT = os.path.join(HERE, "..", "figures", "albumin_loading_vs_logk.png")
OUT_CSV = os.path.join(HERE, "..", "figures",
                       "albumin_loading_vs_logk_data.csv")

WANT = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]        # chain-length order
NCF = {"PFBA": 3, "PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}

# Categorical slots 1 and 2 of the validated palette. Checked with
# scripts/validate_palette.js "#2a78d6,#eb6834" --mode light: all six
# checks pass, worst adjacent pair dE 24.7 protan / 33.6 normal. Light
# surface only -- this is a PNG for a printed report, so the dark-mode
# step is out of scope and deliberately not claimed.
CONV_C, MEAS_C = "#2a78d6", "#eb6834"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, BAND = "#e1e0d9", "#c3c2b7", "#f0efe9"
FIELDS = ["pfas", "n_cf", "log_ratio", "log_k", "category", "source"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 9, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
})


def fit(xs, ys):
    n = len(xs)
    if n < 3 or len(set(xs)) < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    return b, a, (1 - ss_res / ss_tot if ss_tot else float("nan"))


def from_workbook():
    from openpyxl import load_workbook
    ws = load_workbook(SRC)["Coefficients"]
    head = [c.value for c in ws[1]]
    ix = {h: i for i, h in enumerate(head)}
    rows = []
    for r in ws.iter_rows(min_row=2, values_only=True):
        name = str(r[ix["PFAS"]])
        if name not in WANT:
            continue
        route = str(r[ix["Route"]])
        src = str(r[ix["Source"]])
        cat = ("qin" if src.startswith("Qin 2024") else
               "measured" if route == "measured directly" else
               "quenching" if "quenching" in route else "converted")
        try:
            rat = float(str(r[ix["ligand:protein ratio (max)"]]).strip())
        except (TypeError, ValueError):
            rat = None
        rows.append({
            "pfas": name, "n_cf": NCF[name],
            "log_ratio": "" if not rat or rat <= 0 else round(
                math.log10(rat), 4),
            "log_k": float(r[ix["log K (L water / kg protein)"]]),
            "category": cat, "source": src})
    return rows


def read_csv(path):
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8")):
        r["log_k"] = float(r["log_k"])
        r["n_cf"] = int(float(r["n_cf"]))
        r["log_ratio"] = (float(r["log_ratio"])
                          if str(r["log_ratio"]).strip() else "")
        out.append(r)
    return out


def dodge(marks, width=0.085, step=0.115):
    """Lanes for marks sharing an x, assigned over SORTED marks so the
    result does not depend on the order they were read in."""
    out = {}
    for key in {round(m["log_ratio"] / width) for m in marks}:
        same = sorted((m for m in marks
                       if round(m["log_ratio"] / width) == key),
                      key=lambda m: (m["log_k"], m["source"]))
        start = -step * (len(same) - 1) / 2.0
        for j, m in enumerate(same):
            out[id(m)] = start + j * step
    return out


def main():
    if len(sys.argv) > 1:
        rows = read_csv(sys.argv[1])
    else:
        rows = from_workbook()
        with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(sorted(rows, key=lambda r: (r["n_cf"],
                                                    r["category"],
                                                    r["log_k"])))

    order = [c for c in WANT if any(r["pfas"] == c for r in rows)]
    fig, axes = plt.subplots(2, 3, figsize=(12.8, 8.2), sharex=True,
                             sharey=True)
    flat = axes.ravel()
    xs_all = [r["log_ratio"] for r in rows if r["log_ratio"] != ""]
    lo, hi = min(xs_all) - 0.6, max(xs_all) + 0.6
    ylo = min(r["log_k"] for r in rows) - 0.5
    yhi = max(r["log_k"] for r in rows) + 0.5
    summary = []

    for i, c in enumerate(order):
        ax = flat[i]
        here = [r for r in rows if r["pfas"] == c]
        qin = next((r["log_k"] for r in here if r["category"] == "qin"), None)
        pts = [r for r in here if r["log_ratio"] != ""]
        meas = [r for r in pts if r["category"] == "measured"]
        conv = [r for r in pts if r["category"] in ("converted", "quenching")]

        # The gap: the empty x between the two clusters. This band IS
        # the argument, so it is drawn first and everything sits on it.
        gap = None
        if meas and conv:
            left, right = max(m["log_ratio"] for m in meas), \
                min(v["log_ratio"] for v in conv)
            if right > left:
                gap = right - left
                ax.axvspan(left, right, color=BAND, zorder=0)
                ax.annotate(f"{gap:.1f} decades\nwith no data",
                            xy=((left + right) / 2, ylo + 0.25),
                            ha="center", va="bottom", fontsize=7.5,
                            color=MUTED, linespacing=1.3)

        if qin is not None:
            ax.axhline(qin, color=MUTED, lw=1.2, ls=(0, (6, 3)), zorder=2)
            ax.annotate(f"Qin  {qin:.2f}", xy=(hi, qin), xytext=(-3, 3),
                        textcoords="offset points", ha="right",
                        va="bottom", fontsize=7.5, color=MUTED)

        # Pooled fit -- the misleading one, so it is dashed and thin.
        pooled = fit([r["log_ratio"] for r in pts],
                     [r["log_k"] for r in pts])
        if pooled:
            b, a, r2 = pooled
            ax.plot([lo, hi], [a + b * lo, a + b * hi], color=INK2,
                    lw=1.6, ls=(0, (4, 2.5)), zorder=4)

        # Within-group fits, each spanning only its own group's x range.
        within = {}
        for grp, col in ((meas, MEAS_C), (conv, CONV_C)):
            f = fit([r["log_ratio"] for r in grp],
                    [r["log_k"] for r in grp])
            within["m" if col == MEAS_C else "c"] = f
            if not f:
                continue
            b, a, _ = f
            x0 = min(r["log_ratio"] for r in grp) - 0.25
            x1 = max(r["log_ratio"] for r in grp) + 0.25
            ax.plot([x0, x1], [a + b * x0, a + b * x1], color=col,
                    lw=2.2, solid_capstyle="round", zorder=5)

        lanes = dodge(pts)
        for m in sorted(pts, key=lambda r: (r["log_ratio"], r["log_k"])):
            x = m["log_ratio"] + lanes[id(m)]
            if m["category"] == "measured":
                ax.plot(x, m["log_k"], "o", ms=7.5, mfc=MEAS_C,
                        mec=SURFACE, mew=1.4, zorder=7)
            elif m["category"] == "quenching":
                ax.plot(x, m["log_k"], "o", ms=7.5, mfc="none",
                        mec=CONV_C, mew=1.8, zorder=6)
            else:
                ax.plot(x, m["log_k"], "o", ms=7.5, mfc=CONV_C,
                        mec=SURFACE, mew=1.4, zorder=6)

        def lab(f):
            return f"{f[0]:+.2f}" if f else "n/a"
        summary.append((c, pooled, within["m"], within["c"], gap))
        ax.text(0.0, 1.012,
                f"pooled {lab(pooled)}"
                + (f" (r² {pooled[2]:.2f})" if pooled else "")
                + f"     measured {lab(within['m'])}"
                + f"     converted {lab(within['c'])}",
                transform=ax.transAxes, fontsize=7.8, color=MUTED,
                va="bottom")
        ax.set_title(f"{c}   (C{NCF[c]})", loc="left", fontsize=11.5,
                     fontweight="bold", color=INK, pad=20)
        ax.grid(True, color=GRID, lw=0.7, zorder=1)
        ax.set_axisbelow(True)
        ax.set_xlim(lo, hi)
        ax.set_ylim(ylo, yhi)

    # Sixth cell: legend and the reading note.
    cell = flat[len(order)]
    cell.axis("off")
    handles = [
        Line2D([], [], marker="o", ls="none", ms=8, mfc=MEAS_C,
               mec=SURFACE, mew=1.4, label="measured directly"),
        Line2D([], [], marker="o", ls="none", ms=8, mfc=CONV_C,
               mec=SURFACE, mew=1.4, label="converted from a Kd"),
        Line2D([], [], marker="o", ls="none", ms=8, mfc="none",
               mec=CONV_C, mew=1.8,
               label="converted on a quenching exponent"),
        Line2D([], [], color=MEAS_C, lw=2.2, label="fit, measured only"),
        Line2D([], [], color=CONV_C, lw=2.2, label="fit, converted only"),
        Line2D([], [], color=INK2, lw=1.6, ls=(0, (4, 2.5)),
               label="fit through both groups pooled"),
        Line2D([], [], color=MUTED, lw=1.2, ls=(0, (6, 3)),
               label="Qin at this study's operating point"),
        Patch(facecolor=BAND, edgecolor="none",
              label="x range where neither group has data"),
    ]
    cell.legend(handles=handles, loc="upper left", frameon=False,
                fontsize=9.0, handletextpad=0.9, labelspacing=0.60,
                bbox_to_anchor=(0.0, 1.0))
    cell.text(0.0, 0.015,
              "The pooled line is steep and tidy; the two within-group\n"
              "lines are weak and disagree with it on sign. A line through\n"
              "two clusters that never overlap measures the step between\n"
              "them, not a trend inside either.",
              transform=cell.transAxes, fontsize=8.6, color=INK2,
              va="bottom", linespacing=1.5)

    for ax in (axes[1][0], axes[1][1], axes[0][2]):
        ax.set_xlabel("log₁₀  ligand : protein ratio",
                      fontsize=9.5, labelpad=6)
        ax.tick_params(labelbottom=True)
    for ax in (axes[0][0], axes[1][0]):
        ax.set_ylabel("log K  (L water / kg protein)", fontsize=9.5)

    fig.suptitle("The loading slope is a gap between two groups, "
                 "not a trend within either", x=0.045, y=0.983,
                 ha="left", fontsize=15.5, fontweight="bold", color=INK)
    fig.text(0.045, 0.944,
             "Every converted point is a trace-limit K⁰, which has "
             "saturation divided out of it by construction, so it cannot "
             "carry a saturation effect.",
             ha="left", fontsize=10, color=INK2)
    # The gap range is computed, not asserted: the previous caption
    # went stale the moment Allendorf's measured ratios replaced its
    # design ceiling and every gap moved.
    gaps = [g for _, _, _, _, g in summary if g]
    fig.text(0.045, 0.921,
             f"Measured and converted rows share no x range on any "
             f"compound — {min(gaps):.1f} to {max(gaps):.1f} decades "
             f"separate them — so route and loading are collinear.",
             ha="left", fontsize=10, color=INK2)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.905))
    fig.savefig(OUT, dpi=200)
    print(f"wrote {OUT}")
    print(f"wrote {OUT_CSV}\n")
    print(f"{'PFAS':<8}{'pooled':>18}{'measured':>14}{'converted':>14}"
          f"{'gap':>10}")
    for c, p, m, v, g in summary:
        def s(f):
            return f"{f[0]:+.2f} r2 {f[2]:.2f}" if f else "n/a"
        print(f"{c:<8}{s(p):>18}{(f'{m[0]:+.2f}' if m else 'n/a'):>14}"
              f"{(f'{v[0]:+.2f}' if v else 'n/a'):>14}"
              f"{(f'{g:.1f} dec' if g else '-'):>10}")


if __name__ == "__main__":
    main()
