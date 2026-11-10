#!/usr/bin/env python3
"""Albumin binding against Qin 2024's cellular structural protein.

THE QUESTION. Qin et al. 2024, "Baseline Toxicity Model to Identify the
Specific and Nonspecific Effects of Per- and Polyfluoroalkyl Substances
in Cell-Based Bioassays" (Environ Sci Technol 58:5727, doi
10.1021/acs.est.3c09950), Table 1, reports a protein-water partition
coefficient for the STRUCTURAL PROTEIN of HEK293H cells -- a generic
intracellular protein pool, not a transport protein.

Serum albumin is the opposite kind of protein: a dedicated carrier with
evolved fatty-acid pockets.  Plotting them on one axis asks whether
PFAS binding to albumin is a special property of albumin, or whether
PFAS stick to protein generally and albumin is merely the protein they
meet first.

Answering it needs two corrections that are easy to miss.

UNIT BASIS.  Qin's D_SP/w is L_water per L_PROTEIN; the albumin values
in this project are L_water per KG_protein.  They differ by the protein
density, 1/0.733 kg/L, so Qin's values are shifted by
-log10(1/0.733) = -0.135 to put them on the albumin basis.  Plotting
them unshifted would overstate the structural protein by that amount.

WHICH QIN PAPER.  Two different 2024 Qin papers share the label "Qin et
al. 2024" in this database: the plasma-binding paper (3c10824) and this
baseline-toxicity one (3c09950).  Everything here is selected on DOI,
because the label alone has already caused one cross-contamination.

THE FORM.  One row per compound, albumin coefficients as a field of
dots, Qin's structural protein as a single distinct marker, and the
albumin range drawn as a band so the reading is a position judgement:
does the structural protein fall inside the albumin field or below it.
A dot plot answers that at a glance; a bar chart of means destroys it.

COLOUR.  Two categorical hues, #2a78d6 and #eb6834, which clear every
check in the skill's validator on the light surface (worst adjacent CVD
dE 24.7 protan, normal-vision 33.6, contrast >= 3:1).  The albumin
cloud takes a neutral chrome grey rather than a third hue, because it
is the background field the comparison is made against, not a third
peer.  Identity never rests on colour alone: both series are
direct-labelled and the markers differ in shape.

EDITING IT.  Drawn from a CSV, so marks can be changed without editing
this file:

    python3 .../make_qin_companion_figure.py
        reads the database, writes the CSV and the PNG
    python3 .../make_qin_companion_figure.py mine.csv
        reads mine.csv, writes the PNG only

Writes figures/qin_companion_structural_protein.png and
figures/qin_companion_structural_protein_data.csv
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
ALB = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
FIGS = os.path.join(HERE, "..", "figures")
PNG = os.path.join(FIGS, "qin_companion_structural_protein.png")
CSV = os.path.join(FIGS, "qin_companion_structural_protein_data.csv")

COMPANION_DOI = "10.1021/acs.est.3c09950"
V_BAR = 0.733                      # mL/g, albumin partial specific volume
VOL_TO_MASS = -math.log10(1.0 / V_BAR)   # -0.135

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#9a9a93"
INK, MUTED, SURFACE = "#1a1a17", "#6b6b64", "#fcfcfb"


def build():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    sp = {}
    for r in rows:
        if r["doi"] == COMPANION_DOI and "structural protein" in r["protein"]:
            sp[r["chemical"]] = float(r["ka"])
    cell = {}
    for r in rows:
        if r["doi"] == COMPANION_DOI and r["protein"].startswith("HEK293H"):
            cell[r["chemical"]] = float(r["ka"])

    import openpyxl
    ws = openpyxl.load_workbook(ALB, data_only=True)["Coefficients"]
    h = [c.value for c in ws[1]]
    alb = {}
    for r in ws.iter_rows(min_row=2, values_only=True):
        d = dict(zip(h, r))
        alb.setdefault(str(d["PFAS"]), []).append(
            (float(d["log K (L water / kg protein)"]),
             str(d["Source"]), str(d["Species"])))

    out = []
    for cm in sorted(sp, key=lambda c: -sp[c]):
        if cm not in alb:
            continue
        for lk, src, spc in alb[cm]:
            out.append(dict(pfas=cm, series="albumin", log_k=f"{lk:.3f}",
                            source=src, species=spc, basis="kg_protein"))
        out.append(dict(pfas=cm, series="structural protein (Qin 2024)",
                        log_k=f"{sp[cm] + VOL_TO_MASS:.3f}",
                        source="Qin 2024 companion Table 1, D_SP/w",
                        species="human cell line",
                        basis="converted from L_protein by %+.3f" % VOL_TO_MASS))
        if cm in cell:
            out.append(dict(pfas=cm, series="whole cell (Qin 2024)",
                            log_k=f"{cell[cm] + VOL_TO_MASS:.3f}",
                            source="Qin 2024 companion Table 1, D_cell/w",
                            species="human cell line",
                            basis="converted from L_protein by %+.3f"
                                  % VOL_TO_MASS))
    return out


def draw(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    order, seen = [], set()
    for r in rows:
        if r["pfas"] not in seen:
            seen.add(r["pfas"]); order.append(r["pfas"])

    head_in = 1.45
    fig_h = 0.52 * len(order) + head_in + 0.95
    fig, ax = plt.subplots(figsize=(10.0, fig_h))
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    fig.subplots_adjust(top=1 - head_in / fig_h, bottom=0.95 / fig_h,
                        left=0.115, right=0.985)

    for i, cm in enumerate(order):
        y = len(order) - 1 - i
        a = [float(r["log_k"]) for r in rows
             if r["pfas"] == cm and r["series"] == "albumin"]
        if a:
            ax.plot([min(a), max(a)], [y, y], color=GREY, lw=7, alpha=0.22,
                    solid_capstyle="round", zorder=1)
            ax.scatter(a, [y] * len(a), s=26, facecolor="none",
                       edgecolor=GREY, linewidth=1.1, zorder=2)
        for r in rows:
            if r["pfas"] != cm:
                continue
            if r["series"].startswith("structural"):
                ax.scatter([float(r["log_k"])], [y], s=108, marker="D",
                           color=ORANGE, edgecolor=SURFACE, linewidth=2,
                           zorder=4)
            elif r["series"].startswith("whole cell"):
                ax.scatter([float(r["log_k"])], [y], s=96, marker="s",
                           color=BLUE, edgecolor=SURFACE, linewidth=2,
                           zorder=3)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(list(reversed(order)), fontsize=10.5, color=INK)
    ax.set_xlabel("log K  (L water / kg protein)", fontsize=10.5, color=MUTED)
    ax.tick_params(axis="x", colors=MUTED, labelsize=9.5)
    ax.tick_params(axis="y", length=0)
    for sp_ in ("top", "right", "left"):
        ax.spines[sp_].set_visible(False)
    ax.spines["bottom"].set_color("#d8d8d2")
    ax.grid(axis="x", color="#ececE6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(-0.7, len(order) - 0.3)

    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor="none",
               markeredgecolor=GREY, markersize=7,
               label="serum albumin (this database)"),
        Line2D([], [], marker="D", color="none", markerfacecolor=ORANGE,
               markeredgecolor=SURFACE, markersize=9,
               label="HEK293H structural protein  (Qin 2024)"),
        Line2D([], [], marker="s", color="none", markerfacecolor=BLUE,
               markeredgecolor=SURFACE, markersize=9,
               label="HEK293H whole cell  (Qin 2024)"),
    ]
    ax.legend(handles=handles, loc="lower left",
              bbox_to_anchor=(0.0, 1.005, 1.0, 0.1), ncol=3, frameon=False,
              fontsize=9.5, labelcolor=INK, handletextpad=0.45,
              columnspacing=1.6, borderaxespad=0)
    fig.text(0.115, 1 - 0.40 / fig_h,
             "Is PFAS–albumin binding special, or do PFAS stick to protein "
             "generally?", fontsize=13.2, color=INK, va="top")
    # Two lines: one line overflowed the right edge at this figure width.
    fig.text(0.115, 1 - 0.70 / fig_h,
             "Qin's cellular structural protein against every albumin "
             "coefficient in this database, on one basis.",
             fontsize=9.6, color=MUTED, va="top")
    fig.text(0.115, 1 - 0.90 / fig_h,
             f"Qin's values shifted {VOL_TO_MASS:+.3f} log from "
             "L_water/L_protein to L_water/kg_protein (v-bar 0.733 mL/g).",
             fontsize=9.6, color=MUTED, va="top")
    fig.savefig(PNG, dpi=200, facecolor=SURFACE)
    print(f"  wrote {os.path.relpath(PNG, os.path.join(HERE, '..'))}")


def main():
    os.makedirs(FIGS, exist_ok=True)
    if len(sys.argv) > 1:
        rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8")))
        print(f"  {len(rows)} rows from {sys.argv[1]}")
    else:
        rows = build()
        with open(CSV, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"  wrote {os.path.relpath(CSV, os.path.join(HERE, '..'))} "
              f"({len(rows)} rows)")
    draw(rows)

    # The reading, printed so it is not left to the eye alone.
    import collections
    g = collections.defaultdict(dict)
    for r in rows:
        g[r["pfas"]].setdefault(r["series"], []).append(float(r["log_k"]))
    print(f"\n  {'PFAS':10}{'albumin min-max':>20}{'struct prot':>13}"
          f"{'whole cell':>12}  structural protein sits")
    below = inside = 0
    for cm, d in g.items():
        a = d.get("albumin", [])
        s = d.get("structural protein (Qin 2024)", [None])[0]
        c = d.get("whole cell (Qin 2024)", [None])
        if not a or s is None:
            continue
        where = ("INSIDE the albumin range" if min(a) <= s <= max(a)
                 else "BELOW all albumin values" if s < min(a)
                 else "ABOVE all albumin values")
        inside += where.startswith("INSIDE"); below += where.startswith("BELOW")
        print(f"  {cm:10}{f'{min(a):.2f} to {max(a):.2f}':>20}{s:>13.2f}"
              f"{(c[0] if c[0] is not None else float('nan')):>12.2f}  {where}")
    print(f"\n  {inside} inside the albumin field, {below} below it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
