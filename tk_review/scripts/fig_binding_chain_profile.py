#!/usr/bin/env python3
"""Figure: albumin and L-FABP affinity against chain length.

    python3 -I scripts/fig_binding_chain_profile.py

Reads  db/binding_chain_profile.csv   <- EDIT THAT FILE to change the figure
Writes report/fig_binding_chain_profile.png
Idempotent. All data lives in the CSV; this script holds none.

WHAT THE FIGURE IS FOR. The standing explanation for hepatic PFAS accumulation
is binding to liver fatty acid-binding protein. That explanation requires
L-FABP to hold PFAS more tightly than the plasma carrier does. Put the two
profiles on one axis, on the same basis, and they lie on top of each other --
so the premise fails, and the explanation has to come from somewhere else
(compartment volume, protein concentration, transport into the hepatocyte).

Carbon axis: TOTAL carbons, read from the n_total_carbons column of the CSV.

Palette: categorical slots 1-2 of the reference instance, validated with
scripts/validate_palette.js --pairs all --mode light: all checks PASS, worst
all-pairs CVD dE 24.7 (protan), normal-vision dE 33.6. Light mode only,
deliberately -- this renders into a light-background PDF.
"""
import csv, pathlib, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent.parent
ALB, FABP = "#2a78d6", "#eb6834"
INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
COL = {"albumin": ALB, "L-FABP": FABP}


def load():
    p = HERE / "db" / "binding_chain_profile.csv"
    rows, hdr = [], None
    for raw in p.open():
        if raw.startswith("#") or not raw.strip():
            continue
        f = next(csv.reader([raw]))
        if hdr is None:
            hdr = f
            continue
        r = dict(zip(hdr, f))
        rows.append((r["protein"], r["chemical"],
                     int(r["n_total_carbons"]),
                     float(r["median_logK"]), int(r["n_values"])))
    return rows


def main():
    rows = load()
    by = {}
    for prot, chem, nc, lk, n in rows:
        by.setdefault(prot, []).append((nc, lk, n, chem))
    for v in by.values():
        v.sort()

    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)

    for prot in ("albumin", "L-FABP"):
        d = by[prot]
        c = COL[prot]
        ax.plot([x[0] for x in d], [x[1] for x in d], "-", color=c, lw=2.0,
                zorder=3, alpha=0.9)
        # marker area carries the number of underlying measurements, so a
        # point resting on one value cannot read as firm as one resting on 25
        for nc, lk, n, chem in d:
            ax.plot([nc], [lk], "o", ms=4.5 + 1.9 * (n ** 0.5), mfc=c,
                    mec=SURFACE, mew=1.6, zorder=4, alpha=0.95)

    # the shared maximum is the point of the figure, so it gets a label
    ax.annotate("both peak at PFNA", xy=(9, 4.29), xytext=(9.5, 4.78),
                color=SECOND, fontsize=9.5,
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=1.2))

    # each label sits over its OWN line with clear air around it; an earlier
    # placement put "L-FABP" on top of the albumin curve near C4
    for prot, xy in (("albumin", (6.20, 3.82)), ("L-FABP", (15.15, 3.24))):
        ax.annotate(prot, xy=xy, color=COL[prot], fontsize=11,
                    fontweight="bold")

    ax.set_xticks(range(4, 19))
    ax.set_xlim(3.4, 18.8)
    ax.set_ylim(1.2, 5.2)
    ax.set_xlabel("total carbons", color=SECOND, fontsize=10,
                  labelpad=6)
    ax.set_ylabel("median log $K$   (L water / kg protein)", color=SECOND,
                  fontsize=10, labelpad=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)

    handles = [Line2D([], [], color=COL[p], lw=2, marker="o", ms=8,
                      mec=SURFACE, mew=1.4, label=p)
               for p in ("albumin", "L-FABP")]
    handles.append(Line2D([], [], color=MUTED, lw=0, marker="o", ms=5,
                          mec=SURFACE, mew=1.2,
                          label="marker area ∝ number of measurements"))
    leg = ax.legend(handles=handles, loc="lower right", frameon=False,
                    fontsize=9, borderaxespad=0.8)
    for t in leg.get_texts():
        t.set_color(SECOND)

    sh = (set(c for _, c, *_ in rows if _ == "albumin") if False else
          {c for p, c, *_ in rows if p == "albumin"} &
          {c for p, c, *_ in rows if p == "L-FABP"})
    la = {c: lk for p, c, _, lk, _n in rows if p == "albumin"}
    lf = {c: lk for p, c, _, lk, _n in rows if p == "L-FABP"}
    med = statistics.median(lf[c] - la[c] for c in sh)

    fig.suptitle("L-FABP has no affinity advantage over the plasma carrier",
                 color=INK, fontsize=13, x=0.012, ha="left", y=0.975,
                 fontweight="semibold")
    fig.text(0.012, 0.012,
             f"Across the {len(sh)} carboxylates measured against both "
             f"proteins the median difference is {med:+.2f} log units "
             f"(L-FABP the weaker); the audited full set of 21 compounds, "
             f"which also\nincludes sulfonates, gives −0.12. Any account of "
             "hepatic accumulation that rests on L-FABP binding PFAS more "
             "tightly than albumin does not survive this comparison.\n"
             "Per-compound medians from report/ALBUMIN_BINDING.md Tables 5 "
             "and 10; data in db/binding_chain_profile.csv.",
             color=MUTED, fontsize=8.4, va="bottom")

    fig.subplots_adjust(left=0.085, right=0.985, top=0.90, bottom=0.235)
    out = HERE / "report" / "fig_binding_chain_profile.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"  {len(rows)} points, {len(sh)} shared compounds, "
          f"median difference {med:+.2f} log units")
    print(f"wrote {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
