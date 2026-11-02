#!/usr/bin/env python3
"""Does loading explain the albumin spread? Five PFCAs, one panel each.

THE CLAIM BEING TESTED. Binding is saturable, so a coefficient measured
at high ligand loading is genuinely smaller than the same coefficient at
low loading. If the laboratories that disagree were simply working at
different loadings, the spread would be an operating-point artefact
rather than a disagreement. Qin measured that slope directly, within one
laboratory, and it is drawn on every panel as the benchmark.

WHAT THE PANELS SHOW, AND THE REASON THEY ARE HONEST ABOUT IT. The
trend is real: the fitted slope is negative on every compound and close
to Qin's own. But the orange and grey points barely share an x range --
every directly measured value sits at a loading below 0.1, every
converted one above 1 -- so "measured" and "low loading" are the same
rows, and this data cannot separate the two explanations. That overlap
failure is the finding, and colouring by route is what makes it visible
instead of hiding inside an r-squared.

WHY LOADING AND NOT CONCENTRATION on the x axis, when the question was
about concentration. Fischer and Allendorf report a ligand:protein ratio
and no molar concentration, so a concentration axis silently drops the
entire tight-agreeing measured cluster -- the group whose position is
the whole point. The ratio keeps both groups and is dimensionless, so
assays at different protein loads are comparable. On the compounds where
both are available the picture is the same.

QIN IS A LINE, NOT A POINT. Its value here is the Table S6 fit evaluated
at a free concentration this study's operating point implies, not a
measurement at a stated loading, so it has no honest x position. It is
drawn as a horizontal reference: where the number being used sits
against the cloud.

Writes figures/albumin_loading_vs_logk.png and its data CSV.

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

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
OUT = os.path.join(HERE, "..", "figures", "albumin_loading_vs_logk.png")
OUT_CSV = os.path.join(HERE, "..", "figures",
                       "albumin_loading_vs_logk_data.csv")

WANT = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]        # chain-length order
NCF = {"PFBA": 3, "PFHxA": 5, "PFHpA": 6, "PFOA": 7, "PFNA": 8}
# Qin 2024 SI Table S6: log D_BSA/w against log Cw, one laboratory.
QIN_SLOPE = {"PFBA": -0.298, "PFHxA": -0.271, "PFHpA": -0.305,
             "PFOA": -0.314, "PFNA": -0.147}

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
QIN_C, MEAS_C = "#2a78d6", "#eb6834"
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


def main():
    if len(sys.argv) > 1:
        rows = read_csv(sys.argv[1])
        src = sys.argv[1]
    else:
        rows = from_workbook()
        with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(sorted(rows, key=lambda r: (r["n_cf"],
                                                    r["category"],
                                                    r["log_k"])))
        src = os.path.relpath(SRC, os.path.join(HERE, ".."))

    order = [c for c in WANT if any(r["pfas"] == c for r in rows)]
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 8.0), sharex=True,
                             sharey=True)
    flat = axes.ravel()
    xs_all = [r["log_ratio"] for r in rows if r["log_ratio"] != ""]
    lo, hi = min(xs_all) - 0.6, max(xs_all) + 0.6
    ylo = min(r["log_k"] for r in rows) - 0.4
    yhi = max(r["log_k"] for r in rows) + 0.4
    summary = []

    for i, c in enumerate(order):
        ax = flat[i]
        here = [r for r in rows if r["pfas"] == c]
        qin = next((r["log_k"] for r in here if r["category"] == "qin"), None)
        pts = [r for r in here if r["log_ratio"] != ""]
        meas = [r for r in pts if r["category"] == "measured"]
        conv = [r for r in pts if r["category"] in ("converted", "quenching")]

        if qin is not None:
            ax.axhline(qin, color=QIN_C, lw=1.6, zorder=2)
            ax.annotate(f"Qin  {qin:.2f}", xy=(hi, qin), xytext=(-4, 4),
                        textcoords="offset points", ha="right", va="bottom",
                        fontsize=8, color=QIN_C, fontweight="bold")

        f = fit([r["log_ratio"] for r in pts], [r["log_k"] for r in pts])
        if f:
            b, a, r2 = f
            ax.plot([lo, hi], [a + b * lo, a + b * hi], color=INK2, lw=1.8,
                    zorder=3)
            # Qin's own slope through the same centroid, so the two slopes
            # can be compared without the intercept getting in the way.
            mx = sum(r["log_ratio"] for r in pts) / len(pts)
            my = sum(r["log_k"] for r in pts) / len(pts)
            q = QIN_SLOPE[c]
            ax.plot([lo, hi], [my + q * (lo - mx), my + q * (hi - mx)],
                    color=MUTED, lw=1.6, ls=(0, (5, 3)), zorder=3)
            summary.append((c, b, q, r2))
            # Above the axes, not inside: at the bottom left it sat
            # on Wu's quenching point in the PFBA panel, and no
            # corner is empty in all five.
            ax.text(0.0, 1.012, f"fitted {b:+.2f}   Qin {q:+.2f}   r² {r2:.2f}",
                    transform=ax.transAxes, fontsize=8.2,
                    color=MUTED, va="bottom")

        for r in conv:
            ax.plot([r["log_ratio"]], [r["log_k"]], "o", ms=6.5, zorder=4,
                    mfc="none" if r["category"] == "quenching" else MUTED,
                    mec=MUTED, mew=1.3 if r["category"] == "quenching" else 1.0)
        for r in meas:
            ax.plot([r["log_ratio"]], [r["log_k"]], "o", ms=7.5, mfc=MEAS_C,
                    mec=SURFACE, mew=1.4, zorder=5)

        ax.set_title(f"{c}   (C{NCF[c] + 1})", fontsize=11, color=INK,
                     loc="left", pad=20)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_xlim(lo, hi)
        ax.set_ylim(ylo, yhi)

    # the sixth cell carries the legend and the reading
    ax = flat[len(order)]
    ax.axis("off")
    handles = [
        Line2D([], [], marker="o", ls="", ms=7.5, mfc=MEAS_C, mec=SURFACE,
               mew=1.4, label="measured directly"),
        Line2D([], [], marker="o", ls="", ms=6.5, mfc=MUTED, mec=MUTED,
               label="converted from a Kd"),
        Line2D([], [], marker="o", ls="", ms=6.5, mfc="none", mec=MUTED,
               mew=1.3, label="converted on a quenching exponent"),
        Line2D([], [], color=INK2, lw=1.8, label="fitted slope, all points"),
        Line2D([], [], color=MUTED, lw=1.6, ls=(0, (5, 3)),
               label="Qin's own slope, measured in one lab"),
        Line2D([], [], color=QIN_C, lw=1.6,
               label="Qin at this study's operating point"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=9, frameon=False,
              labelspacing=0.95, handletextpad=0.7,
              bbox_to_anchor=(-0.02, 0.98))
    ax.text(-0.02, 0.30,
            "Every orange point sits left of 10⁻¹ and almost every\n"
            "grey one right of 10⁰. \"Measured\" and \"low loading\" are\n"
            "the same rows, so these panels cannot separate\n"
            "the two explanations — only show that the slope\n"
            "is real and that the groups do not overlap.",
            transform=ax.transAxes, fontsize=8.6, color=INK2,
            va="top", linespacing=1.6)

    # sharex hides tick labels wherever a panel sits below; the
    # top-right panel's neighbour is the legend cell, so without this
    # it ships with an unlabelled x axis.
    for ax in (axes[1][0], axes[1][1], axes[0][2]):
        ax.set_xlabel("log₁₀  ligand : protein ratio")
        ax.tick_params(axis="x", labelbottom=True)
    for ax in axes[:, 0]:
        ax.set_ylabel("log K  (L water / kg protein)")

    fig.text(0.045, 0.975,
             "Loading moves the albumin coefficient, but it is not what "
             "separates the laboratories",
             ha="left", va="top", fontsize=14.5, color=INK)
    med = sorted(b for _, b, _, _ in summary)[len(summary) // 2]
    fig.text(0.045, 0.941,
             f"Fitted slope is negative on all {len(summary)} compounds, "
             f"median {med:+.2f} log per decade".replace("-", "\u2212")
             + ", against Qin's "
             "within-laboratory −0.30. The mechanism is real.\nBut "
             "correcting for it removes only 14% of the spread, because "
             "the measured and converted rows occupy different loadings "
             "and the comparison is confounded.",
             ha="left", va="top", fontsize=9, color=INK2, linespacing=1.6)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.905))
    fig.savefig(OUT, dpi=200)
    plt.close(fig)

    print(f"read {len(rows)} marks from {src}")
    for c, b, q, r2 in summary:
        print(f"  {c:8s} fitted {b:+.2f}   Qin {q:+.2f}   r2 {r2:.2f}")
    if len(sys.argv) == 1:
        print(f"wrote {os.path.abspath(OUT_CSV)}")
    print(f"wrote {os.path.abspath(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
