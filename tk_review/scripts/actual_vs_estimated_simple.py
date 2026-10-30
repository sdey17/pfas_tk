#!/usr/bin/env python3
"""Measured against estimated protein-water coefficients, three sources only.

A deliberately small view of the comparison. The full one
(logk_converted_vs_measured.py) carries 29 admitted pairs, ten excluded
ones, five error scales and the reasoning behind each; this file answers
one question -- does an estimate land on its measurement -- for three
sources picked because each estimates differently:

  Chen 2025   CONVERTED from a Kd, with the stoichiometry ASSUMED. The
              paper fits a single-site Langmuir, which fixes n = 1 by
              construction, so every value moves one for one with that
              choice. The large, well-powered case.

  Moro 2022   CONVERTED from a Kd, with the stoichiometry MEASURED, by
              isothermal titration calorimetry. The same arithmetic on a
              better input, and the only such source in the database.

  Qin 2024    MEASURED, and the one that needs the CORRECTION FACTOR. Its
              albumin column is reported at half-saturation, as Eq. 15's
              K_spec = C_max/(2 K_d), so it sits 0.301 log (a factor of 2)
              below the Cw -> 0 plateau everyone else reports. It is also
              split into a specific and a non-specific term, neither of
              which is a total. Both are handled here: the specific term is
              lifted by log10(2) and the two are summed in linear space.
              Without that it would look like a disagreement; with it, it
              is a measurement.

TWO VALUES FOR ONE THING, TWICE, HANDLED OPPOSITE WAYS. Moro reports
each compound twice, and Qin reports albumin in two columns. They look
alike and are not. Moro's pair are REPLICATE determinations -- each row
is a complete fit, carrying both site classes where the compound has two
(PFOA: "0.4 (high-affinity set); 29.7 (low-affinity set)" in one row) --
so they average, in log space, and the count travels in `n_values`.
Qin's pair are COMPONENTS, a specific and a non-specific term, neither a
total on its own -- so they sum, in linear space, after the specific one
is lifted onto the plateau. Averaging components or summing replicates
would both be wrong, and nothing but reading the rows says which is
which.

Fischer 2024 is the yardstick throughout -- the only source reporting
albumin partition coefficients directly on a mass basis, for both human
and bovine albumin, which is what lets Chen and Moro (human) and Qin
(bovine) each be compared without mixing species.

Writes db/actual_vs_estimated_simple.csv and
figures/actual_vs_estimated_simple.png.

Run:  python3 tk_review/scripts/actual_vs_estimated_simple.py
"""
import collections
import csv
import importlib.util
import math
import os
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "..", "db", "actual_vs_estimated_simple.csv")
OUT_PNG = os.path.join(HERE, "..", "figures",
                       "actual_vs_estimated_simple.png")


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B = _load("baf", "build_binding_affinity_table.py")
L = _load("lkc", "logk_converted_vs_measured.py")

REFERENCE = "Fischer 2024"
# (source, how it was arrived at, colour, marker, the one-line caveat)
SERIES = [
    ("Chen 2025", "converted from Kd, n assumed = 1", "#2a78d6", "o",
     "single-site Langmuir fixes n at 1"),
    ("Moro 2022", "converted from Kd, n measured by ITC", "#1baf7a", "D",
     "n fitted from the titration isotherm"),
    ("Qin 2024 (plasma)", "measured, +0.301 correction applied", "#eb6834",
     "s", "Eq. 15 is at half-saturation; specific + non-specific summed"),
]

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
GRID, BAND = "#e6e4e0", "#dfe8f4"
BETWEEN_LAB = 0.59

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "axes.spines.top": False,
    "axes.spines.right": False,
})

FIELDS = ["pfas", "protein", "species", "source", "estimate_kind",
          "estimated_logK", "measured_logK", "residual", "n_values", "note"]


def _mean_log(vals):
    """Replicates of one quantity, averaged where they belong: in log space."""
    return sum(vals) / len(vals)


def collect():
    """One row per (compound, protein, species) a series shares with Fischer."""
    recs = B.build(filter_compounds=False)[1]
    B.stamp_basis(recs, B.basis_index(recs))
    # totals() folds Qin's two albumin columns into one comparable number,
    # applying the half-saturation correction on the way.
    # Lists, not single values: a source may determine the same compound
    # more than once, and quietly keeping whichever row came last would
    # pick one of Moro's two PFHxA fits at random.
    measured = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in L.totals(recs):
        measured[r["source"]][L.key(r)].append(float(r["logK"]))

    converted = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in recs:
        if L.route(r) == "converted" and r["logK"]:
            converted[r["source"]][L.key(r)].append(float(r["logK"]))

    ref = measured.get(REFERENCE, {})
    rows = []
    for source, kind, _, _, note in SERIES:
        pool = converted.get(source) or measured.get(source) or {}
        for k, vals in sorted(pool.items()):
            if k not in ref:
                continue
            est, meas = _mean_log(vals), _mean_log(ref[k])
            rows.append({
                "pfas": k[0], "protein": k[1], "species": k[2],
                "source": source, "estimate_kind": kind,
                "estimated_logK": round(est, 3),
                "measured_logK": round(meas, 3),
                "residual": round(est - meas, 3),
                "n_values": len(vals), "note": note})
    return rows


def plot(rows):
    fig, ax = plt.subplots(figsize=(7.0, 7.4))
    lo = math.floor(min(min(r["estimated_logK"], r["measured_logK"])
                        for r in rows) * 2) / 2 - 0.3
    hi = math.ceil(max(max(r["estimated_logK"], r["measured_logK"])
                       for r in rows) * 2) / 2 + 0.3

    ax.fill_between([lo, hi], [lo - BETWEEN_LAB, hi - BETWEEN_LAB],
                    [lo + BETWEEN_LAB, hi + BETWEEN_LAB],
                    color=BAND, lw=0, zorder=1)
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1.0, zorder=2)
    ax.text(hi - 0.08, hi - 0.22, "1:1", color=MUTED, ha="right", fontsize=8)
    ax.text(lo + 0.15, lo + 0.34,
            f"the pale band is ±{BETWEEN_LAB:.2f} log: what two direct\n"
            f"measurements of the same thing normally differ by",
            color="#3f6ea8", fontsize=7.6, va="bottom")

    handles = []
    for source, kind, colour, marker, _ in SERIES:
        sub = [r for r in rows if r["source"] == source]
        if not sub:
            continue
        ax.plot([r["measured_logK"] for r in sub],
                [r["estimated_logK"] for r in sub], marker, ms=6.2,
                mfc=colour, mec=SURFACE, mew=1.1, ls="", zorder=4)
        med = statistics.median(r["residual"] for r in sub)
        handles.append(Line2D([], [], marker=marker, ls="", ms=6.2,
                              mfc=colour, mec=SURFACE, mew=1.1,
                              label=f"{source} — {kind}\n"
                                    f"    {len(sub)} compound"
                                    f"{'s' if len(sub) != 1 else ''}, "
                                    f"median {med:+.2f} log"))

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.grid(True, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"measured log K  ({REFERENCE}, "
                  f"L$_{{water}}$/kg protein)")
    ax.set_ylabel("estimated log K   (same units)")
    med = {src: statistics.median(r["residual"] for r in rows
                                  if r["source"] == src)
           for src, *_ in SERIES if any(r["source"] == src for r in rows)}
    n_chen = sum(1 for r in rows if r["source"] == "Chen 2025")
    # An earlier title read "an estimate lands on its measurement only when
    # the stoichiometry is measured too" -- the right reading of the project
    # as a whole, and far more than these points can carry: Moro contributes
    # ONE compound. The title now says what is on the chart.
    fig.text(0.075, 0.972,
             "Three routes to a partition coefficient, one yardstick",
             ha="left", va="top", fontsize=13, color=INK)
    fig.text(0.075, 0.936,
             f"Chen's {n_chen} converted values sit "
             f"{med.get('Chen 2025', 0):.2f} log low, all resting on an "
             "assumed n = 1.\nMoro's single ITC-based conversion lands at "
             f"{med.get('Moro 2022', 0):+.2f}. Qin's measurement, once "
             "lifted\noff half-saturation, sits "
             f"{med.get('Qin 2024 (plasma)', 0):+.2f} away — the "
             "ordinary gap between two laboratories.",
             ha="left", va="top", fontsize=8.6, color=INK2, linespacing=1.5)
    ax.legend(handles=handles, loc="upper left", fontsize=7.4,
              labelspacing=1.0, handletextpad=0.5, borderpad=0.4)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)


def main():
    rows = collect()
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} comparisons against {REFERENCE}:")
    print(f"  {'source':20s}{'n':>4s}{'median':>9s}{'|median|':>10s}  "
          f"protein")
    for source, _, _, _, _ in SERIES:
        sub = [r for r in rows if r["source"] == source]
        if not sub:
            print(f"  {source:20s}   0        -         -  (no shared "
                  f"compound with {REFERENCE})")
            continue
        res = [r["residual"] for r in sub]
        prot = sorted({r["protein"] for r in sub})
        rep = sum(1 for r in sub if r["n_values"] > 1)
        print(f"  {source:20s}{len(sub):4d}"
              f"{statistics.median(res):+9.2f}"
              f"{statistics.median(abs(x) for x in res):10.2f}  "
              f"{', '.join(prot)}"
              + (f"   ({rep} averaged over replicates)" if rep else ""))
    print()
    print(f"wrote {os.path.abspath(OUT_CSV)}")
    print(f"wrote {os.path.abspath(OUT_PNG)}")
    plot(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
