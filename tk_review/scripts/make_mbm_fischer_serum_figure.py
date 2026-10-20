#!/usr/bin/env python3
"""Four pictures of the Fischer 2024 serum mass balance.

A: predicted vs measured log D_serum/w, nothing fitted
B: what carries the binding -- albumin, globulin, water -- compound by compound
C: the residual against the one thing it tracks, D_glob/D_HSA
D: how far f_unbound moves on serum composition alone, and why that depends
   on whether albumin and globulin travel together

Run:  python3 scripts/make_mbm_fischer_serum_figure.py
"""
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_fischer_serum as M  # noqa: E402

OUT = os.path.join(HERE, "..", "figures")
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PLUM, SAND = "#8a5cd1", "#c8a33a"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
EMPTY = "#eceae6"

CLS_COLOR = {"carboxylate": ORANGE, "sulfonate": BLUE, "sulfonamide": SAND,
             "fluorotelomer": AQUA, "ether": PLUM}
CLS_MARK = {"carboxylate": "o", "sulfonate": "^", "sulfonamide": "s",
            "fluorotelomer": "D", "ether": "d"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 8.5, "axes.edgecolor": "#d8d7d3", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
})


def panel_a(ax, recs, meas, pred):
    lo, hi = 1.6, 4.3
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1.0, zorder=1)
    for f, lab in ((0.5, None), (2.0, None)):
        ax.plot([lo, hi], [lo + math.log10(f), hi + math.log10(f)],
                color=EMPTY, lw=1.0, zorder=1)
    seen = set()
    for r, m, p in zip(recs, meas, pred):
        c = CLS_COLOR[r["cls"]]
        ax.scatter(m, p, s=38, c=c, marker=CLS_MARK[r["cls"]],
                   edgecolor="white", linewidth=0.7, zorder=3,
                   label=r["cls"] if r["cls"] not in seen else None)
        seen.add(r["cls"])
    for nm, dx, dy, ha in (("PFHpS", 0.05, 0.08, "left"),
                           ("6:2 FtS", 0.05, -0.16, "left"),
                           ("7:3 FTCA", 0.05, -0.12, "left"),
                           ("5:3 FTCA", -0.05, -0.17, "right")):
        for r, m, p in zip(recs, meas, pred):
            if r["cmpd"] == nm:
                ax.annotate(nm, (m, p), (m + dx, p + dy), fontsize=7.5,
                            color=INK2, ha=ha)
    ax.text(lo + 0.08, lo + 0.08 + math.log10(2), "2x over", fontsize=7,
            color=MUTED, ha="left", va="bottom", rotation=38)
    ax.text(hi - 0.08, hi - 0.08 - math.log10(2), "2x under", fontsize=7,
            color=MUTED, ha="right", va="top", rotation=38)
    ax.set_xlim(lo, hi), ax.set_ylim(lo, hi)
    ax.set_xlabel("measured  log D$_{serum/w}$  (L$_w$ / L$_{serum}$)")
    ax.set_ylabel("predicted from HSA + globulin")
    ax.set_title("A.  Components predict the whole matrix\n"
                 "26 PFAS, no fitted parameter", fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", handletextpad=0.3)
    res = pred - meas
    ax.text(0.97, 0.06,
            "bias %.2fx   |err| %.2fx\nRMSE %.3f log   R$^2$ %.2f"
            % (10 ** res.mean(), 10 ** np.abs(res).mean(),
               float(np.sqrt((res ** 2).mean())),
               1 - float((res ** 2).sum()) / float(((meas - meas.mean()) ** 2).sum())),
            transform=ax.transAxes, fontsize=7.5, color=INK2,
            ha="right", va="bottom")


def panel_b(ax, recs):
    order = sorted(recs, key=lambda r: r["glob"] - r["HSA"])
    vw = M.vf_water(M.C_ALB_NHANES, M.C_GLOB_NHANES)
    y = np.arange(len(order))
    alb, glo, wat = [], [], []
    for r in order:
        a = 10 ** r["HSA"] * M.C_ALB_NHANES
        g = 10 ** r["glob"] * M.C_GLOB_NHANES
        t = vw + a + g
        alb.append(a / t), glo.append(g / t), wat.append(vw / t)
    alb, glo, wat = np.array(alb), np.array(glo), np.array(wat)
    ax.barh(y, alb, color=BLUE, height=0.72, label="albumin")
    ax.barh(y, glo, left=alb, color=ORANGE, height=0.72, label="globulin")
    ax.barh(y, wat, left=alb + glo, color=EMPTY, height=0.72, label="water")
    ax.set_yticks(y)
    ax.set_yticklabels([r["cmpd"] for r in order], fontsize=7.5)
    ax.set_xlim(0, 1), ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("share of D$_{serum/w}$")
    ax.set_title("B.  Albumin carries it until globulin catches up\n"
                 "sorted by D$_{glob}$/D$_{HSA}$", fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, -0.17), handletextpad=0.4)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def panel_c(ax, recs, res):
    x = np.array([r["glob"] - r["HSA"] for r in recs])
    ax.axhline(0, color=MUTED, lw=1.0)
    for f in (math.log10(2), -math.log10(2)):
        ax.axhline(f, color=EMPTY, lw=1.0)
    for r, xx, rr in zip(recs, x, res):
        ax.scatter(xx, rr, s=38, c=CLS_COLOR[r["cls"]],
                   marker=CLS_MARK[r["cls"]], edgecolor="white",
                   linewidth=0.7, zorder=3)
    for nm, dx, dy, ha in (("PFHpS", 0.08, 0.0, "left"),
                           ("6:2 FtS", 0.08, 0.02, "left"),
                           ("7:3 FTCA", 0.08, 0.02, "left"),
                           ("5:3 FTCA", -0.08, 0.02, "right")):
        for r, xx, rr in zip(recs, x, res):
            if r["cmpd"] == nm:
                ax.annotate(nm, (xx, rr), (xx + dx, rr + dy),
                            fontsize=7.5, color=INK2, ha=ha, va="center")
    ax.set_ylim(res.min() - 0.10, res.max() + 0.10)
    sl, ic = np.polyfit(x, res, 1)
    xs = np.linspace(x.min(), x.max(), 2)
    ax.plot(xs, sl * xs + ic, color=MUTED, lw=1.0, ls="--")
    ax.axvline(-1.0, color=MUTED, lw=0.8, ls=":")
    ax.text(-0.97, 0.62, "Fischer's\nD$_{glob}$/D$_{HSA}$ = 0.1", fontsize=7,
            color=MUTED, va="top")
    ax.set_xlabel("log$_{10}$ (D$_{glob/w}$ / D$_{HSA/w}$)")
    ax.set_ylabel("log$_{10}$ (predicted / measured)")
    ax.set_title("C.  The residual, and the four compounds Fischer\n"
                 "also named as misses", fontsize=10, color=INK, pad=8)
    ax.text(0.97, 0.95, "Pearson r = %+.2f" % float(np.corrcoef(x, res)[0, 1]),
            transform=ax.transAxes, fontsize=7.5, color=INK2,
            ha="right", va="top")


def panel_d(ax, varrows):
    order = sorted(varrows, key=lambda t: t[5] / t[4])
    y = np.arange(len(order))
    indep = [t[3] / t[1] for t in order]
    anti = [t[5] / t[4] for t in order]
    for i, (lo, hi) in enumerate(zip(anti, indep)):
        ax.plot([lo, hi], [i, i], color=EMPTY, lw=3.0, solid_capstyle="round",
                zorder=1)
    ax.scatter(anti, y, s=34, c=BLUE, edgecolor="white", linewidth=0.7,
               zorder=3, label="albumin and globulin anti-correlated")
    ax.scatter(indep, y, s=34, c=ORANGE, edgecolor="white", linewidth=0.7,
               zorder=3, label="albumin and globulin independent")
    ax.axvline(1.0, color=MUTED, lw=0.9)
    for nm, val, lab in (("PFHxS", 2.5, "Fischer"), ("10:2 FtS", 1.6, None)):
        for i, t in enumerate(order):
            if t[0]["cmpd"] == nm:
                ax.scatter([val], [i], s=52, facecolor="none",
                           edgecolor=INK, linewidth=1.1, zorder=4,
                           label="Fischer's reported range" if lab else None)
    ax.set_yticks(y)
    ax.set_yticklabels([t[0]["cmpd"] for t in order], fontsize=7.5)
    ax.set_xlim(0.9, 3.0), ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("fold range in f$_{unbound}$ across NHANES serum composition")
    ax.set_title("D.  Same chemical, same dose: what serum composition\n"
                 "alone does to the free fraction", fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", handletextpad=0.3)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def main():
    data = M.load()
    recs = sorted((r for r in data.values()
                   if "HSA" in r and "glob" in r and "serum" in r),
                  key=lambda r: (r["cls"], r["nCF"]))
    meas = np.array([r["serum"] for r in recs])
    pred = M.predict(recs, M.C_ALB_NHANES, M.C_GLOB_NHANES)
    ca, cg, _cl, _r = M.fit(recs)
    res = M.predict(recs, ca, cg) - meas

    tot = M.C_ALB_NHANES + M.C_GLOB_NHANES
    varrows = []
    for r in recs:
        indep = [M.vf_water(a, g) / 10 ** M.predict([r], a, g)[0]
                 for a in M.NHANES_ALB_RANGE for g in M.NHANES_GLOB_RANGE]
        anti = [M.vf_water(a, max(tot - a, 0.0))
                / 10 ** M.predict([r], a, max(tot - a, 0.0))[0]
                for a in M.NHANES_ALB_RANGE]
        varrows.append((r, min(indep), None, max(indep), min(anti), max(anti)))

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 11.0))
    panel_a(axes[0][0], recs, meas, pred)
    panel_b(axes[0][1], recs)
    panel_c(axes[1][0], recs, res)
    panel_d(axes[1][1], varrows)
    fig.suptitle("A serum mass balance fitted to one species and one lab: "
                 "Fischer et al. 2024, 26 PFAS",
                 fontsize=12, color=INK, y=0.975)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    path = os.path.join(OUT, "fig19_fischer_serum_mbm.png")
    fig.savefig(path, dpi=190, bbox_inches="tight")
    print("wrote", os.path.abspath(path))


if __name__ == "__main__":
    main()
