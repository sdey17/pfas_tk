#!/usr/bin/env python3
"""Two plain pictures of what the binding database covers.

A: compound x quantity   -- which of Kd, K_protein/water, K_membrane/water
                            and f_unbound exists for each compound
B: compound x species    -- which species each compound was measured in

Run:  python3 scripts/make_binding_coverage_figure.py
"""
import collections
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_binding_affinity_table import build  # noqa: E402

OUT = os.path.join(HERE, "..", "figures")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
EMPTY = "#eceae6"

KIND = {"Kd": "Kd", "Ka": "Kd", "log D": "K prot/water",
        "log K_membrane": "K mem/water", "f_unbound": "f unbound"}
ORDER = ["Kd", "K prot/water", "K mem/water", "f unbound"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
})


def grid(ax, data, rows, cols, title, counts=True):
    cmap = ListedColormap([EMPTY, "#c5dbf4", "#7fb0e6", BLUE])
    m = [[min(data.get((r, c), 0), 3) for c in cols] for r in rows]
    ax.imshow(m, cmap=cmap, aspect="auto", vmin=0, vmax=3)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=7)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            n = data.get((r, c), 0)
            if n and counts:
                ax.text(j, i, str(n), ha="center", va="center", fontsize=6.5,
                        color="white" if n >= 3 else INK)
    for j in range(len(cols) + 1):
        ax.axvline(j - .5, color=SURFACE, lw=1.2)
    for i in range(len(rows) + 1):
        ax.axhline(i - .5, color=SURFACE, lw=1.2)
    ax.set_title(title, fontsize=10, color=INK, pad=8)
    ax.tick_params(length=0)


def main():
    _, out, _, _, _ = build()

    byq = collections.Counter()
    bysp = collections.Counter()
    for r in out:
        k = KIND.get(r["measure"])
        if not k:
            continue
        byq[(r["pfas"], k)] += 1
        sp = r["species"].split(" (")[0].split(";")[0].strip()
        if sp.startswith("not applicable"):
            sp = "synthetic membrane"
        bysp[(r["pfas"], sp)] += 1

    # order compounds by how many of the four quantities they carry, then rows
    score = collections.Counter()
    for (c, k), n in byq.items():
        score[c] += 100 + n
    chems = [c for c, _ in score.most_common(40)]
    species = [s for s, _ in collections.Counter(
        s for (_, s) in bysp).most_common()]
    species = [s for s in species if any(bysp.get((c, s)) for c in chems)]

    fig, (a1, a2) = plt.subplots(
        1, 2, figsize=(13.5, 9.5),
        gridspec_kw={"width_ratios": [1, 1.5], "wspace": 0.42})
    grid(a1, byq, chems, ORDER,
         "A. What quantity exists, per compound\n"
         "cell = number of measurements")
    grid(a2, bysp, chems, species,
         "B. Which species it was measured in\n"
         "cell = number of measurements")
    fig.suptitle("Binding database coverage: the 40 best-covered compounds "
                 f"of {len(score)}",
                 fontsize=12, color=INK, y=0.955)
    path = os.path.join(OUT, "fig18_binding_coverage.png")
    fig.savefig(path, dpi=190, bbox_inches="tight")
    print(f"wrote {os.path.abspath(path)}")

    # --- the completeness ranking ---------------------------------------
    print("\nMost complete single studies (compounds x quantities):")
    per = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in out:
        k = KIND.get(r["measure"])
        if k:
            per[r["source"]][k].add(r["pfas"])
    rank = []
    for s, d in per.items():
        allc = set().union(*d.values())
        rank.append((len(allc) * len(d), len(allc), len(d), s, d))
    for _, nc, nk, s, d in sorted(rank, reverse=True)[:8]:
        kinds = ", ".join(f"{k} ({len(v)})" for k, v in
                          sorted(d.items(), key=lambda kv: -len(kv[1])))
        print(f"  {s:18s} {nc:3d} compounds x {nk} quantit{'y' if nk==1 else 'ies'}"
              f"  -- {kinds}")


if __name__ == "__main__":
    main()
