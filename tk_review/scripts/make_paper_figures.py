#!/usr/bin/env python3
"""Build the paper's four figures (print, light surface).

    python3 -I scripts/make_paper_figures.py

Writes figures/paper_fig1..4.png. Idempotent: same inputs, same bytes.

COLOUR. Slots 1-3 of the validated categorical palette, which are the three
that clear the all-pairs CVD and normal-vision floors in light mode
(validator: worst all-pairs CVD dE 9.2, normal-vision dE 24.0). Aqua sits
below 3:1 contrast on the light surface, so the relief rule applies and every
series carries a direct label or legend entry as well as its colour; each
figure's data also appears in a table in the paper. Identity is never
colour-alone.

FORM. Figure 1 is the paper's central result, which until now had no figure.
Its job is comparing two quantities across compounds in two species, so it is
a grouped bar rather than a stack: the shares include negative values (the
distribution term), which a stack cannot show honestly.
"""
import csv
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent.parent
FIGS = HERE / "figures"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8982"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"   # validated slots 1-3
GRID = "#e4e3de"


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, length=3, width=0.8)
    ax.grid(axis="x", color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.patch.set_facecolor(SURFACE)
    p = FIGS / name
    fig.savefig(p, dpi=200, facecolor=SURFACE, bbox_inches="tight",
                metadata={"Software": None, "Creation Time": None})
    plt.close(fig)
    print(f"wrote figures/{name}")


def rows(path):
    return list(csv.DictReader((HERE / path).open()))


# ------------------------------------------------------------------ fig 1
def fig1():
    d = rows("db/species_decomposition.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharex=True)
    for ax, comp, lbl in zip(axes, ("human/rat", "human/mouse"),
                             ("Rat", "Mouse")):
        sub = sorted([r for r in d if r["comparison"] == comp],
                     key=lambda r: (int(r["n_total_carbons"]), r["chemical"]))
        y = range(len(sub))
        b = [float(r["pct_binding"]) for r in sub]
        tu = [float(r["pct_tubule"]) for r in sub]
        h = 0.36
        ax.barh([i + h / 2 + 0.01 for i in y], b, height=h, color=S1,
                zorder=3, label="binding")
        ax.barh([i - h / 2 - 0.01 for i in y], tu, height=h, color=S2,
                zorder=3, label="tubular")
        ax.axvline(0, color=MUTED, linewidth=0.9, zorder=2)
        ax.set_yticks(list(y))
        ax.set_yticklabels([f"{r['chemical']} (C{r['n_total_carbons']})"
                            for r in sub], fontsize=8, color=INK)
        ax.set_title(f"{lbl} vs human", fontsize=10, color=INK, pad=8,
                     loc="left")
        ax.set_xlabel("share of the log half-life ratio (%)", fontsize=8.5,
                      color=INK2)
        style(ax)
        ax.invert_yaxis()
    h_, l_ = axes[0].get_legend_handles_labels()
    fig.legend(h_, l_, frameon=False, fontsize=8.5, labelcolor=INK,
               loc="upper right", bbox_to_anchor=(0.995, 1.012), ncol=2)
    fig.suptitle("Binding carries the mouse difference and not the rat's",
                 fontsize=11.5, color=INK, x=0.055, ha="left", y=1.012)
    save(fig, "paper_fig1.png")


# ------------------------------------------------------------------ fig 2
def fig2():
    d = rows("db/lfabp_within_method.csv")
    xs = [int(r["total_carbons"]) for r in d]
    ys = [float(r["diff_logK_lfabp_minus_albumin"]) for r in d]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    icept = my - slope * mx
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.axhline(0, color=MUTED, linewidth=0.9, zorder=2)
    gx = [min(xs) - 0.4, max(xs) + 0.4]
    ax.plot(gx, [icept + slope * v for v in gx], color=S1, linewidth=2,
            zorder=3, alpha=0.55)
    ax.scatter(xs, ys, s=58, color=S1, zorder=4, edgecolor=SURFACE,
               linewidth=1.6)
    for r, x, y in zip(d, xs, ys):
        if r["chemical"] in ("PFBS", "PFHxS", "PFOS", "PFUnDA", "PFDoDA",
                             "PFTrDA", "PFBA"):
            ax.annotate(r["chemical"], (x, y), textcoords="offset points",
                        xytext=(6, -3), fontsize=7.5, color=INK2)
    ax.set_xlabel("total carbons", fontsize=8.5, color=INK2)
    # the direction lives on the axis label, not in a floating note that
    # collided with the PFBS and PFHxS point labels in the first render
    ax.set_ylabel("Δ log K  (L-FABP − albumin)\n← L-FABP binds more weakly",
                  fontsize=8.5, color=INK2)
    ax.text(0.98, 0.06, f"slope +{slope:.3f} per carbon\n"
                        f"crosses 0 at C{-icept/slope:.1f}",
            transform=ax.transAxes, ha="right", fontsize=8, color=INK2)
    ax.annotate("parity", xy=(max(xs) + 0.35, 0), fontsize=7.5, color=MUTED,
                va="center", ha="left")
    style(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_title("L-FABP's disadvantage closes as the chain lengthens",
                 fontsize=11, color=INK, loc="left", pad=8)
    save(fig, "paper_fig2.png")


# ------------------------------------------------------------------ fig 3
def fig3():
    d = rows("db/dose_sign_test.csv")
    d = sorted(d, key=lambda r: (r["species"], float(r["cl_ratio"])))
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.axvline(1, color=MUTED, linewidth=1.0, zorder=2)
    lab, ys = [], []
    for i, r in enumerate(d):
        c = S1 if r["species"] == "rat" else S2
        pfos_rat = r["chemical"] == "PFOS" and r["species"] == "rat"
        ax.scatter(float(r["cl_ratio"]), i, s=72 if pfos_rat else 46,
                   color=c, zorder=4, edgecolor=SURFACE,
                   linewidth=1.8 if pfos_rat else 1.2)
        lab.append(f"{r['chemical']} {r['sex'][:1].upper()}")
        ys.append(i)
    ax.set_yticks(ys)
    ax.set_yticklabels(lab, fontsize=6.6, color=INK)
    ax.set_xscale("log")
    ax.set_xticks([0.2, 0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"])
    ax.set_xlabel("clearance at high dose / clearance at low dose",
                  fontsize=8.5, color=INK2)
    style(ax)
    # ratio < 1 is a FALL (left); ratio > 1 is a RISE (right). These were
    # transposed in the first render, which inverted the paper's mechanism
    # claim; the sides are asserted against the axis direction below.
    assert ax.get_xlim()[0] < 1 < ax.get_xlim()[1]
    ax.text(0.02, 0.985, "← falls\nsecretion saturates",
            transform=ax.transAxes, ha="left", va="top", fontsize=8,
            color=MUTED, style="italic")
    ax.text(0.98, 0.985, "rises →\nreabsorption saturates",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            color=MUTED, style="italic")
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", color=S1,
                              label="rat", markersize=7),
                       Line2D([], [], marker="o", linestyle="", color=S2,
                              label="mouse", markersize=7)],
              frameon=False, fontsize=8.5, loc="lower right", labelcolor=INK)
    ax.set_title("No rat series falls with dose except PFOS",
                 fontsize=11, color=INK, loc="left", pad=8)
    save(fig, "paper_fig3.png")


# ------------------------------------------------------------------ fig 4
def fig4():
    d = rows("db/primary_2026/thakur2024_table1_transporters.csv")
    order = ["OAT1", "OAT3", "OAT2", "URAT1", "OAT4", "OATP1A1", "OATP1A3"]
    d = sorted(d, key=lambda r: order.index(r["protein"]))

    def val(s):
        return float(s.split("+/-")[0]) if "+/-" in s else 0.0

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    w, x = 0.26, range(len(d))
    for off, key, col, name in ((-w - 0.01, "human_pmol_mg", S1, "human"),
                                (0.0, "rat_pmol_mg", S2, "rat"),
                                (w + 0.01, "mouse_pmol_mg", S3, "mouse")):
        vs = [val(r[key]) for r in d]
        ax.bar([i + off for i in x], vs, width=w, color=col, zorder=3,
               label=name)
        for i, (r, v) in enumerate(zip(d, vs)):
            if v == 0:
                # the rodent Oatp1a family has no close human ortholog; the
                # other blanks are genes absent from that genome
                txt = ("no ortholog" if (name == "human"
                                         and r["protein"].startswith("OATP"))
                       else "absent")
                ax.text(i + off, 0.6, txt, rotation=90, fontsize=6.4,
                        ha="center", va="bottom", color=MUTED)
    ax.axvline(3.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)),
               zorder=2)
    ax.text(1.75, 34, "secretory / exchange", fontsize=8.5, color=INK2,
            ha="center")
    ax.text(5.0, 34, "REABSORPTIVE", fontsize=8.5, color=INK, ha="center",
            fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels([r["protein"] for r in d], fontsize=8, color=INK)
    ax.set_ylabel("pmol / mg membrane protein", fontsize=8.5, color=INK2)
    ax.set_ylim(0, 40)
    style(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)
    ax.set_title("The reabsorptive transporters are not shared",
                 fontsize=11, color=INK, loc="left", pad=8)
    save(fig, "paper_fig4.png")




# ------------------------------------------------------------------ fig 5
def fig5():
    """Human tubular factors: every direct value sits far below unity.

    T is recomputed here from the primary renal clearances and unbound
    fractions rather than copied from the paper, then asserted against the
    published Table 9 so the figure cannot drift from the text.
    """
    GFR_H = 2570.0
    fu = {}
    for r in rows("db/plasma_unbound_fraction.csv"):
        if r["experiment"] == "Ryu 2024" and r["fu"] \
                and r["species"].split(" (")[0] == "human":
            fu[r["chemical"]] = float(r["fu"])
    # CL_R, source, total carbons, and the former upper bound (Table 9)
    DATA = [("PFHxS", 0.027, "Zhang 2013", 6, 0.036),
            ("PFHpA", 0.674, "Fujii 2015", 7, 0.051),
            ("PFOA", 0.030, "Harada 2005", 8, 0.019),
            ("PFOA", 0.044, "Fujii 2015", 8, 0.019),
            ("PFOS", 0.015, "Harada 2005", 8, 0.046),
            ("PFNA", 0.038, "Fujii 2015", 9, 0.032),
            ("PFDA", 0.015, "Fujii 2015", 10, 0.094),
            ("PFUnDA", 0.005, "Fujii 2015", 11, 0.266),
            ("PFDoDA", 0.005, "Fujii 2015", 12, 3.08)]
    PUBLISHED_T = {("PFHxS", 0.027): 0.0151, ("PFHpA", 0.674): 0.0589,
                   ("PFOA", 0.030): 0.0048, ("PFOA", 0.044): 0.0070,
                   ("PFOS", 0.015): 0.0078, ("PFNA", 0.038): 0.0079,
                   ("PFDA", 0.015): 0.0094, ("PFUnDA", 0.005): 0.0046,
                   ("PFDoDA", 0.005): 0.0185}
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    ax.axhline(1, color=MUTED, linewidth=1.1, zorder=2)
    ax.text(12.45, 1.12, "filtration alone", fontsize=7.5, color=MUTED,
            ha="right")
    for chem, clr, src, nc, bound in DATA:
        T = clr / (fu[chem] * GFR_H)
        assert abs(T - PUBLISHED_T[(chem, clr)]) < 0.0003, (chem, T)
        ax.plot([nc, nc], [T, bound], color=GRID, linewidth=1.4, zorder=2)
        ax.scatter(nc, bound, s=40, facecolor="none", edgecolor=S2,
                   linewidth=1.6, zorder=3)
        ax.scatter(nc, T, s=58, color=S1, zorder=4, edgecolor=SURFACE,
                   linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("total carbons", fontsize=8.5, color=INK2)
    ax.set_ylabel("tubular factor T", fontsize=8.5, color=INK2)
    ax.set_xlim(5.3, 12.7)
    ax.set_xticks(range(6, 13))
    # PFHpA is the one point where the measured value exceeds its nominal
    # bound. The bound came from Abraham's TOTAL clearance (0.584) while the
    # measurement is Fujii's RENAL clearance (0.674): two studies, not a
    # violated inequality. Marked rather than hidden.
    ax.annotate("PFHpA: measured renal clearance (Fujii)\n"
                "exceeds the bound's total clearance\n(Abraham) by 15%",
                xy=(7, 0.0589), xytext=(7.45, 0.30), fontsize=7.2,
                color=INK2, ha="left",
                arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.8))
    ax.legend(handles=[
        Line2D([], [], marker="o", linestyle="", color=S1, markersize=7,
               label="measured, from direct renal clearance"),
        Line2D([], [], marker="o", linestyle="", markerfacecolor="none",
               markeredgecolor=S2, markersize=7,
               label="former upper bound, from total clearance")],
        frameon=False, fontsize=8, loc="upper center",
        bbox_to_anchor=(0.5, -0.17), ncol=2, labelcolor=INK)
    style(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_title("Human reabsorption is 94-99.5% complete across the chain",
                 fontsize=11, color=INK, loc="left", pad=8)
    save(fig, "paper_fig5.png")


# ------------------------------------------------------------------ fig 6
def fig6():
    """The cynomolgus reductio: the required unbound fraction is impossible."""
    # compound, required f_u, multiple of the human value (Table 7)
    DATA = [("PFBA", 17.7, 77), ("PFHxA", 1.96, 51), ("PFHxS", 0.0170, 25),
            ("PFOA", 0.729, 297), ("PFOS", 0.0095, 13)]
    DATA = sorted(DATA, key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    ax.axvspan(1, 60, color=S2, alpha=0.08, zorder=1)
    ax.axvline(1, color=S2, linewidth=1.6, zorder=3)
    y = range(len(DATA))
    ax.barh(list(y), [d[1] for d in DATA], height=0.5, color=S1, zorder=4)
    for i, (chem, req, mult) in enumerate(DATA):
        ax.annotate(f"{req:g}   ({mult}x the human value)",
                    (req, i), xytext=(7, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK2)
    ax.set_yticks(list(y))
    ax.set_yticklabels([d[0] for d in DATA], fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ax.set_xlim(0.006, 60)
    ax.set_xlabel("unbound fraction the monkey would need for binding alone "
                  "to explain its half-life", fontsize=8.5, color=INK2)
    ax.text(1.3, 0.28, "IMPOSSIBLE:\nf_u cannot exceed 1",
            fontsize=8.5, color=S2, va="center")
    style(ax)
    ax.set_title("Binding alone cannot explain the cynomolgus monkey",
                 fontsize=11, color=INK, loc="left", pad=8)
    save(fig, "paper_fig6.png")




# ------------------------------------------------------------------ fig 7
def fig7():
    """Two measurements of rat-vs-human binding that point opposite ways.

    Panel (a) is the directional disagreement; panel (b) is why neither
    settles it, the fluorimetric assay missing the plasma scale 100-fold.
    """
    d = [r for r in rows("db/starnes_vs_ryu_rat.csv")
         if r["fu_ratio_rat_over_human"]]
    d = sorted(d, key=lambda r: float(r["kd_ratio_rat_over_human"]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9),
                                 gridspec_kw={"width_ratios": [1.25, 1]})

    # (a) the two ratios, on either side of parity
    y = range(len(d))
    a1.axvline(1, color=MUTED, linewidth=1.1, zorder=2)
    for i, r in enumerate(d):
        k, f = (float(r["kd_ratio_rat_over_human"]),
                float(r["fu_ratio_rat_over_human"]))
        a1.plot([k, f], [i, i], color=GRID, linewidth=1.6, zorder=2)
        a1.scatter(k, i, s=58, color=S1, zorder=4, edgecolor=SURFACE,
                   linewidth=1.4)
        a1.scatter(f, i, s=58, color=S2, zorder=4, edgecolor=SURFACE,
                   linewidth=1.4)
    a1.set_yticks(list(y))
    a1.set_yticklabels([r["chemical"] for r in d], fontsize=8.5, color=INK)
    a1.set_xscale("log")
    # matplotlib adds log MINOR ticks that collided with these labels
    a1.set_xticks([0.5, 0.7, 1, 1.5, 2])
    a1.set_xticks([], minor=True)
    a1.set_xticklabels(["0.5", "0.7", "1", "1.5", "2"])
    a1.set_xlim(0.4, 2.3)
    a1.set_ylim(-1.75, len(d) - 0.25)
    a1.set_xlabel("rat / human ratio", fontsize=8.5, color=INK2)
    a1.text(0.62, -0.5, "rat binds tighter", fontsize=7.5, color=MUTED,
            ha="center", style="italic")
    a1.text(1.6, -0.5, "rat leaves more free", fontsize=7.5, color=MUTED,
            ha="center", style="italic")
    a1.legend(handles=[
        Line2D([], [], marker="o", linestyle="", color=S1, markersize=7,
               label="K$_d$, isolated albumin (Starnes, DSF)"),
        Line2D([], [], marker="o", linestyle="", color=S2, markersize=7,
               label="f$_u$, whole plasma (Ryu, dialysis)")],
        frameon=False, fontsize=7.8, loc="lower left", labelcolor=INK)
    style(a1)
    a1.set_title("(a) they disagree in direction", fontsize=9.5, color=INK,
                 loc="left", pad=6)

    # (b) why neither settles it: the gap to the plasma value, per compound
    gaps = []
    lim = [5e-4, 3]
    a2.plot(lim, lim, color=MUTED, linewidth=1.1, linestyle=(0, (4, 3)),
            zorder=2)
    # alternate the label offsets: adjacent points sit close enough on the
    # log axis that same-side labels overlapped
    for j, r in enumerate(sorted(d, key=lambda x: float(x["fu_rat"]))):
        m = float(r["fu_rat"])
        imp = float(r["fu_implied_by_kd_rat"])
        gaps.append(imp / m)
        a2.plot([m, m], [m, imp], color=S2, linewidth=1.5, zorder=3,
                alpha=0.75)
        a2.scatter(m, imp, s=58, color=S1, zorder=4, edgecolor=SURFACE,
                   linewidth=1.4)
        dy = 9 if j % 2 else -14
        a2.annotate(f"{imp/m:.0f}x", (m, imp), xytext=(-2, dy),
                    textcoords="offset points", fontsize=7.4, color=INK2,
                    ha="center")
    a2.set_xscale("log"); a2.set_yscale("log")
    a2.set_xlim(*lim); a2.set_ylim(*lim)
    a2.set_xlabel("measured f$_u$, rat plasma", fontsize=8.5, color=INK2)
    a2.set_ylabel("f$_u$ implied by the K$_d$", fontsize=8.5, color=INK2)
    a2.text(0.06, 0.88, "identity", transform=a2.transAxes, fontsize=7.5,
            color=MUTED, rotation=34)
    a2.text(0.97, 0.06, f"gap {min(gaps):.0f}x to {max(gaps):.0f}x,\n"
                        f"largest where binding\nis tightest",
            transform=a2.transAxes, ha="right", fontsize=7.8, color=INK2)
    style(a2)
    a2.grid(axis="y", color=GRID, linewidth=0.7)
    a2.set_title("(b) and the assay misses the plasma scale",
                 fontsize=9.5, color=INK, loc="left", pad=6)

    fig.suptitle("Neither rat binding measurement can adjudicate the other",
                 fontsize=11.5, color=INK, x=0.04, ha="left", y=1.03)
    save(fig, "paper_fig7.png")


# ------------------------------------------------------------------ fig 8
def fig8():
    """How big the biliary alternative is, against the renal route."""
    # Harada 2007, human, n = 4 paired serum/bile samples. Renal ranges are
    # the project's human values (Harada 2005, Zhang 2013, Fujii 2015).
    DATA = [("PFOA", 1.06, 0.89, 0.030, 0.044),
            ("PFOS", 2.98, 0.97, 0.012, 0.037)]
    # A DOT PLOT, NOT BARS. On a log axis a bar begins at the truncated
    # axis edge rather than at zero, so its length encodes nothing; and the
    # renal figure is a RANGE, which a bar from the origin misrepresents.
    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    for i, (chem, gross, reab, rlo, rhi) in enumerate(DATA):
        net = gross * (1 - reab)
        assert abs(net - (0.117 if chem == "PFOA" else 0.089)) < 0.001
        ax.plot([rlo, rhi], [i - 0.22, i - 0.22], color=S3, linewidth=5,
                solid_capstyle="butt", zorder=3,
                label="renal (range across sources)" if i == 0 else None)
        ax.plot([net, gross], [i, i], color=GRID, linewidth=1.4, zorder=2)
        ax.scatter(gross, i + 0.22, s=70, color=S1, zorder=4,
                   edgecolor=SURFACE, linewidth=1.4,
                   label="biliary, gross" if i == 0 else None)
        ax.scatter(net, i, s=70, color=S2, zorder=4, edgecolor=SURFACE,
                   linewidth=1.4,
                   label="biliary, net of reabsorption" if i == 0 else None)
        ax.annotate(f"{gross:.2f}", (gross, i + 0.22), xytext=(8, 0),
                    textcoords="offset points", va="center", fontsize=7.8,
                    color=INK2)
        ax.annotate(f"{net:.3f}  ({reab*100:.0f}% reabsorbed)", (net, i),
                    xytext=(8, 0), textcoords="offset points", va="center",
                    fontsize=7.8, color=INK2)
        ax.annotate(f"{rlo:.3f}-{rhi:.3f}", (rhi, i - 0.22), xytext=(8, 0),
                    textcoords="offset points", va="center", fontsize=7.8,
                    color=INK2)
        ax.annotate("", xy=(net, i - 0.1), xytext=(rhi, i - 0.1),
                    arrowprops=dict(arrowstyle="->", color=MUTED,
                                    linewidth=0.9))
        ax.text(math.sqrt(net * rhi), i - 0.09, f"{net/rhi:.1f}-{net/rlo:.1f}x",
                fontsize=7.2, color=MUTED, ha="center", va="bottom")
    ax.set_yticks(range(len(DATA)))
    ax.set_yticklabels([d[0] for d in DATA], fontsize=9.5, color=INK)
    ax.set_ylim(-0.65, len(DATA) - 0.35)
    ax.set_xscale("log")
    ax.set_xlim(0.008, 30)
    ax.set_xlabel("clearance (mL/kg/day), human", fontsize=8.5, color=INK2)
    ax.legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK)
    style(ax)
    ax.set_title("The non-renal route is the larger one, even after "
                 "reabsorption", fontsize=11, color=INK, loc="left", pad=8)
    # The caveat belongs ON the figure, not only in its caption: these are
    # HUMAN clearances, and using them to size the alternative behind the
    # RAT dose response is a cross-species step rather than a measurement.
    fig.text(0.0, -0.10,
             "Human values throughout. Their bearing on the rat dose "
             "response rests on Harada et al. (2007)\nreporting biliary "
             "excretion as comparable in rats and humans "
             "\u2014 a cross-species inference,\nnot a rat measurement. "
             "No rat biliary or renal dose series exists for PFOS.",
             fontsize=7.4, color=MUTED, ha="left", va="top", style="italic")
    save(fig, "paper_fig8.png")


if __name__ == "__main__":
    FIGS.mkdir(exist_ok=True)
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6(); fig7(); fig8()
