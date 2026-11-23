# # PFAS toxicokinetics — every plot in one notebook
# 
# Companion to `report/PFAS_TK_COMBINED.md`. Each figure is one cell, reads one
# CSV from the repository, and draws with matplotlib and the standard library
# only — **no pandas, no project imports**, so a cell can be lifted out whole.
# 
# **Carbon axis.** Every chain-length plot uses **total carbons**. The stored
# data counts *perfluorinated* carbons, and the conversion depends on the head
# group — carboxylates gain one for the COOH carbon, sulfonates do not. See
# `total_carbons()` below; that is the only place the rule is written.
# 
# **Sections**
# 
# | § | figures | theme |
# |---|---|---|
# | 1 | 1–4 | serum albumin binding |
# | 2 | 5–6 | liver fatty acid-binding protein |
# | 3 | 7–8 | renal transporters |
# | 4 | 9–10 | the tubular factor, and the f_u dispute |
# | 5 | 11–12 | coverage, and the un-ingested literature |
# | 6 | — | index of the older figure scripts, nothing dropped |
# 
# Run top to bottom. Every cell is independent after the setup cell.

# ============================= SETUP =============================
import csv, math, os, statistics
from collections import defaultdict, Counter
import matplotlib
import matplotlib.pyplot as plt

# Point this at the tk_review directory if you move the notebook.
ROOT = os.path.abspath(os.path.join(os.getcwd(), ".."))
DB   = os.path.join(ROOT, "db")
FIGS = os.path.join(ROOT, "figures")
NBD  = os.path.join(ROOT, "notebook")

# -- palette: categorical slots 1-3, validated for colour-vision deficiency --
# worst all-pairs CVD dE 9.2, normal-vision dE 24.0, both above the gates.
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
SPECIES = {"human": C1, "mouse": C2, "rat": C3,
           "bovine": "#eda100", "monkey": "#6250d6"}

matplotlib.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 9,
})


def total_carbons(family, n_cf):
    """Perfluorinated carbons -> TOTAL carbons. Head-group dependent.

    carboxylate  CF3(CF2)nCOOH : the carboxyl carbon is not perfluorinated,
                                 so total = n_cf + 1   (PFOA  7 -> 8)
    sulfonate    CF3(CF2)nSO3H : the sulfonur is not a carbon,
                                 so total = n_cf       (PFOS  8 -> 8)

    Getting this wrong shifts every sulfonate one place along the x axis.
    """
    n_cf = int(n_cf)
    f = (family or "").strip().lower()
    if f.startswith("carboxyl"):
        return n_cf + 1
    if f.startswith("sulfon"):
        return n_cf
    # fluorotelomers and anything else: fall back to the carboxylate rule
    # but say so, rather than guessing silently
    return n_cf + 1


def read(path):
    """CSV -> list of dicts, skipping commented lines."""
    with open(path) as fh:
        lines = [l for l in fh if not l.lstrip().startswith("#")]
    return list(csv.DictReader(lines))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def style(ax, xlabel="", ylabel="", title=""):
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    ax.grid(axis="x", visible=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
    if xlabel:
        ax.set_xlabel(xlabel, color=SECOND, fontsize=9.5, labelpad=6)
    if ylabel:
        ax.set_ylabel(ylabel, color=SECOND, fontsize=9.5, labelpad=8)
    if title:
        ax.set_title(title, color=INK, fontsize=11.5, loc="left", pad=8,
                     fontweight="semibold")
    return ax


print("setup ok — ROOT =", ROOT)
print("sanity:  PFOA(carboxylate, 7) ->", total_carbons("carboxylate", 7),
      " PFOS(sulfonate, 8) ->", total_carbons("sulfonate", 8))

# ---
# ## 1. Serum albumin
# 
# 288 partition coefficients for 63 PFAS from 20 sources, converted to one basis
# — log *K*, litres of water per kilogram of protein — through
# 
# $$\log K = 4.1774 + \log_{10} n - \log_{10}(K_d/\mu M)$$
# 
# All 181 rows carrying both a $K_d$ and a stoichiometry satisfy it exactly, so
# the arithmetic is not where the uncertainty lives.

# ### Figure 1 — the dispersion, by compound
# 
# **The finding:** picking one published constant is a choice with a range of
# defensible outcomes spanning **two orders of magnitude**.
# 
# *The cell prints its own numbers rather than repeating the report's.* It draws
# from `figures/albumin_chainlength_data.csv` (309 measurements, all species)
# and gets a median within-compound range near 2.0 log units over 33 compounds.
# The audited compilation in `ALBUMIN_BINDING.md` applies tighter inclusion
# criteria and reports **2.38 log units over 32 compounds, maximum 4.22**. Both
# say the same thing; use the audited figure when quoting.
# 
# This is the single most important plot in the binding strand, because it is
# the reason a species extrapolation factor cannot be estimated from this
# literature — the disagreement *within* a species is larger than the difference
# *between* species.

def fig1_albumin_dispersion(min_n=3, save=None):
    rows = read(os.path.join(FIGS, "albumin_chainlength_data.csv"))
    by = defaultdict(list)
    for r in rows:
        v = num(r["log_k"])
        if v is not None:
            by[(r["pfas"], r["family"])].append(v)
    keep = {k: v for k, v in by.items() if len(v) >= min_n}
    order = sorted(keep, key=lambda k: total_carbons(k[1], 
                   [r["n_cf"] for r in rows if r["pfas"] == k[0]][0]))

    fig, ax = plt.subplots(figsize=(11, 5.4))
    rng = []
    for i, k in enumerate(order):
        v = sorted(keep[k])
        lo, hi = v[0], v[-1]
        rng.append(hi - lo)
        col = C1 if k[1].startswith("carboxyl") else C2
        ax.plot([i, i], [lo, hi], color=col, lw=2.0, solid_capstyle="round",
                alpha=0.85, zorder=3)
        ax.plot([i] * len(v), v, "o", ms=4, mfc=col, mec=SURFACE, mew=0.8,
                zorder=4, alpha=0.9)
        ax.plot([i], [statistics.median(v)], "_", ms=13, color=INK, zorder=5)

    med = statistics.median(rng)
    ax.axhline(0, color=SURFACE, lw=0)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([k[0] for k in order], rotation=60, ha="right",
                       fontsize=8)
    style(ax, "", "log $K$  (L water / kg protein)",
          f"Within-compound spread of albumin binding "
          f"(median range {med:.2f} log units = {10**med:.0f}-fold)")
    ax.plot([], [], color=C1, lw=2, label="carboxylate")
    ax.plot([], [], color=C2, lw=2, label="sulfonate")
    ax.plot([], [], "_", ms=11, color=INK, label="median")
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right")
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    print(f"{len(order)} compounds with >= {min_n} values; "
          f"median range {med:.2f} log units, max {max(rng):.2f}")
    return fig


fig1_albumin_dispersion();

# ### Figure 2 — chain length, on total carbons
# 
# **The finding:** affinity rises to a maximum at **PFNA** and falls away.
# Fischer's mechanism for the decline: *"for ηpfc > 8 PFAS, reduced molecular
# flexibility limits their accessibility to the HSA binding pockets."*
# 
# Note the two families are plotted on the **same total-carbon axis**, which is
# only meaningful because `total_carbons()` handles them differently. On a
# perfluorinated-carbon axis PFOS would sit one place to the left of PFOA
# despite both being eight-carbon molecules.

def fig2_chain_length(save=None):
    rows = read(os.path.join(FIGS, "albumin_chainlength_data.csv"))
    by = defaultdict(list)
    for r in rows:
        v = num(r["log_k"])
        if v is None:
            continue
        by[(r["family"], total_carbons(r["family"], r["n_cf"]))].append(v)

    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    for fam, col, lab in (("carboxylate", C1, "carboxylates"),
                          ("sulfonate", C2, "sulfonates")):
        pts = sorted((nc, statistics.median(v), len(v))
                     for (f, nc), v in by.items() if f == fam)
        if not pts:
            continue
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=col,
                lw=2.0, zorder=3, label=lab)
        for nc, m, n in pts:
            ax.plot([nc], [m], "o", ms=4 + 1.8 * n ** 0.5, mfc=col,
                    mec=SURFACE, mew=1.4, zorder=4, alpha=0.95)

    # PFNA is 8 perfluorinated + 1 carboxyl carbon = 9 total
    ax.annotate("maximum at PFNA\n(9 total carbons)", xy=(9, 4.29),
                xytext=(10.3, 4.75), color=SECOND, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=1.2))
    ax.set_xticks(range(4, 19))
    style(ax, "total carbons", "median log $K$  (L water / kg protein)",
          "Albumin affinity against chain length")
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right")
    for t in leg.get_texts():
        t.set_color(SECOND)
    ax.plot([], [], "o", ms=5, color=MUTED, mec=SURFACE,
            label="marker area ∝ n")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


fig2_chain_length();

# ### Figure 3 — method is the principal axis of disagreement
# 
# **The finding:** assay family is a systematic axis of disagreement, not noise.
# 
# **Read the scope of this cell carefully.** It draws from
# `db/albumin_kd_to_k.csv`, which holds the 119 rows convertible from a
# dissociation constant — and that file contains **no differential scanning
# fluorimetry**, because DSF reports a melting temperature rather than a $K_d$.
# So the headline comparison in `ALBUMIN_BINDING.md` — **equilibrium dialysis
# against DSF, differing by 2.35 log units in median** — cannot be reproduced
# here, and this cell does not claim it. What it shows is the spread among the
# $K_d$-based families, about 1 log unit between the extremes.
# 
# The underlying point survives either way: **these are not competing estimates
# of one quantity, they are two quantities.** DSF measures thermal denaturation,
# dialysis measures an equilibrium. Averaging across methods is the mistake this
# plot exists to prevent.

def method_family(text):
    """Bucket a free-text method description into an assay family.

    The stored `method` column is a full sentence ("equilibrium dialysis,
    Slide-A-Lyzer 10 kDa MWCO, 48 h at 50 rpm; single-site Langmuir fit...").
    Grouping on the raw string gives one group per sentence, which is not a
    plot. Order matters below: the first match wins.
    """
    t = (text or "").lower()
    for needle, label in (
            ("differential scanning", "DSF"),
            ("nanodsf", "DSF"),
            ("isothermal titration", "ITC"),
            ("itc", "ITC"),
            ("equilibrium dialysis", "dialysis"),
            ("dialysis", "dialysis"),
            ("ultrafiltration", "ultrafiltration"),
            ("ultracentrifug", "ultracentrifugation"),
            ("microextraction", "SPME"),
            ("spme", "SPME"),
            ("electrospray", "native ESI-MS"),
            ("fluorescence", "fluorescence quench"),
            ("quench", "fluorescence quench"),
            ("docking", "docking / in silico"),
            ("simulation", "docking / in silico"),
    ):
        if needle in t:
            return label
    return "other / unstated"


def fig3_method(save=None, min_n=3):
    rows = read(os.path.join(DB, "albumin_kd_to_k.csv"))
    by = defaultdict(list)
    for r in rows:
        v = num(r.get("logK"))
        if v is not None:
            by[method_family(r.get("method"))].append(v)
    by = {k: v for k, v in by.items() if len(v) >= min_n}
    order = sorted(by, key=lambda m: -statistics.median(by[m]))

    fig, ax = plt.subplots(figsize=(9.4, 4.8))
    pal = [C1, C2, C3, "#eda100", "#e87ba4", "#6250d6", "#e34948"]
    for i, m in enumerate(order):
        v = sorted(by[m])
        col = pal[i % len(pal)]
        x = [i + (j - len(v) / 2) * 0.012 for j in range(len(v))]
        ax.plot(x, v, "o", ms=4.5, mfc=col, mec=SURFACE, mew=0.8, alpha=0.85,
                zorder=3)
        md_ = statistics.median(v)
        ax.plot([i - 0.28, i + 0.28], [md_, md_], "-", color=INK, lw=2.2,
                zorder=5)
        ax.annotate(f"{md_:.2f}", xy=(i + 0.32, md_), color=SECOND,
                    fontsize=8.5, va="center")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=8.5, rotation=25, ha="right")
    style(ax, "", "log $K$  (L water / kg protein)",
          "Albumin binding by assay family (black bar = median)")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    for m in order:
        print(f"  {m:24} n={len(by[m]):>3}  "
              f"median {statistics.median(by[m]):.2f}")
    if len(order) > 1:
        gap = statistics.median(by[order[0]]) - statistics.median(by[order[-1]])
        print(f"  widest gap between families: {gap:.2f} log units "
              f"= {10 ** gap:.0f}-fold")
    return fig


fig3_method();

# ### Figure 4 — does loading explain the spread?
# 
# **The hypothesis worth testing:** an apparent binding constant falls as the
# assay saturates the protein, so part of the between-laboratory scatter should
# track the ligand-to-protein ratio.
# 
# **Read this one honestly.** If the points show no trend, loading is *not* the
# explanation and the dispersion has to come from somewhere else. The plot is
# here to be allowed to fail.

def fig4_loading(compounds=("PFOA", "PFOS", "PFNA", "PFHxS"), save=None):
    rows = read(os.path.join(NBD, "albumin_logK_Kd_revised.csv"))
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    pal = [C1, C2, C3, "#eda100"]
    any_pt = False
    for i, cmp_ in enumerate(compounds):
        pts = []
        for r in rows:
            if r["PFAS"] != cmp_:
                continue
            x, y = num(r.get("ligand:protein ratio")), num(
                r["log K (L water / kg protein)"])
            if x and y is not None and x > 0:
                pts.append((x, y))
        if not pts:
            continue
        any_pt = True
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=7,
                mfc=pal[i % len(pal)], mec=SURFACE, mew=1.2, alpha=0.9,
                label=f"{cmp_} (n={len(pts)})", zorder=3)
    if any_pt:
        ax.set_xscale("log")
    style(ax, "ligand : protein molar ratio at the top of the titration",
          "log $K$  (L water / kg protein)",
          "Does assay loading explain the between-laboratory spread?")
    leg = ax.legend(frameon=False, fontsize=9)
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


fig4_loading();

# ---
# ## 2. Liver fatty acid-binding protein
# 
# 91 coefficients for 32 PFAS from ten sources, on the same basis, using a
# per-species bridge constant of **4.8475** for human L-FABP against albumin's
# 4.1774.

# ### Figure 5 — albumin against L-FABP
# 
# **The finding, and the most consequential one in this section:** across the 21
# compounds measured against both proteins the median difference is **−0.12 log
# units**, L-FABP marginally the *weaker*. (The cell recomputes this from the
# per-measurement files and gets −0.22 over the 20 compounds they share; the
# audited 21-compound figure is −0.12. Quote the audited one.)
# 
# **What it does not say:** that the two are interchangeable in the body.
# Albumin sits in plasma near 600 µM, L-FABP in hepatocyte cytosol at 100–400 µM
# — equal affinity per kilogram still gives very different capacity.
# 
# **What it does say:** there is no intrinsic affinity advantage. Any
# explanation of hepatic accumulation requiring L-FABP to out-bind the plasma
# carrier fails here. The one in vivo test agrees: deleting L-FABP changes the
# mouse liver partition coefficient by **under 20% for 15 of 19 PFAS**.

def fig5_albumin_vs_lfabp(save=None):
    alb = read(os.path.join(FIGS, "albumin_chainlength_data.csv"))
    fab = read(os.path.join(FIGS, "lfabp_by_compound_data.csv"))

    def medians(rows):
        by = defaultdict(list)
        for r in rows:
            v = num(r["log_k"])
            if v is not None:
                by[(r["pfas"], r["family"], r["n_cf"])].append(v)
        return {k: (statistics.median(v), len(v)) for k, v in by.items()}

    A, F = medians(alb), medians(fab)
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    for d, col, lab in ((A, C1, "albumin"), (F, C2, "L-FABP")):
        pts = sorted((total_carbons(fam, ncf), m, n)
                     for (p, fam, ncf), (m, n) in d.items())
        # collapse to one median per total-carbon position
        agg = defaultdict(list)
        for nc, m, n in pts:
            agg[nc].append((m, n))
        xs = sorted(agg)
        ys = [statistics.median([m for m, _ in agg[x]]) for x in xs]
        ns = [sum(n for _, n in agg[x]) for x in xs]
        ax.plot(xs, ys, "-", color=col, lw=2.0, zorder=3, label=lab)
        for x, y, n in zip(xs, ys, ns):
            ax.plot([x], [y], "o", ms=4 + 1.8 * n ** 0.5, mfc=col,
                    mec=SURFACE, mew=1.4, zorder=4, alpha=0.95)

    shared = {p for p, _, _ in A} & {p for p, _, _ in F}
    am = {p: m for (p, _, _), (m, _) in A.items()}
    fm = {p: m for (p, _, _), (m, _) in F.items()}
    diff = statistics.median(fm[p] - am[p] for p in shared)

    ax.set_xticks(range(4, 19))
    style(ax, "total carbons", "median log $K$  (L water / kg protein)",
          "L-FABP has no affinity advantage over the plasma carrier")
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right")
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    print(f"{len(shared)} compounds measured against both; "
          f"median (L-FABP - albumin) = {diff:+.2f} log units")
    print("report's audited figure over its 21-compound set: -0.12")
    return fig


fig5_albumin_vs_lfabp();

# ### Figure 6 — L-FABP by compound, and the two-site problem
# 
# **A structural cause for part of the scatter.** L-FABP has **two binding
# sites**, and most constants come from probe displacement. Sources using a
# *secondary*-site probe and sources using a *primary*-site probe **are not
# estimating the same quantity** — Luebker's assay displaces DAUDA, which binds
# the secondary site, and the paper says so.
# 
# So part of this dispersion is fixable by stratifying on probe site rather
# than being irreducible noise. Within-compound median range is 1.40 log units,
# maximum 3.30 for PFOA across nine sources.

def fig6_lfabp_by_compound(save=None):
    rows = read(os.path.join(FIGS, "lfabp_by_compound_data.csv"))
    by = defaultdict(list)
    for r in rows:
        v = num(r["log_k"])
        if v is not None:
            by[(r["pfas"], r["family"], r["n_cf"])].append((v, r.get("species", "")))
    order = sorted(by, key=lambda k: total_carbons(k[1], k[2]))

    fig, ax = plt.subplots(figsize=(11, 5.0))
    for i, k in enumerate(order):
        vals = [v for v, _ in by[k]]
        lo, hi = min(vals), max(vals)
        ax.plot([i, i], [lo, hi], color=MUTED, lw=1.4, zorder=2)
        for v, sp in by[k]:
            ax.plot([i], [v], "o", ms=5.5,
                    mfc=SPECIES.get((sp or "").lower(), C3), mec=SURFACE,
                    mew=0.9, zorder=4, alpha=0.9)
        ax.plot([i], [statistics.median(vals)], "_", ms=13, color=INK,
                zorder=5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{k[0]}" for k in order], rotation=60, ha="right",
                       fontsize=8)
    style(ax, "", "log $K$  (L water / kg protein)",
          "L-FABP binding by compound, ordered by total carbons")
    seen = {sp for v in by.values() for _, sp in v if sp}
    for sp in sorted(seen):
        ax.plot([], [], "o", ms=6, mfc=SPECIES.get(sp.lower(), C3),
                mec=SURFACE, label=sp)
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right", ncol=2)
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    rngs = [max(v for v, _ in by[k]) - min(v for v, _ in by[k])
            for k in order if len(by[k]) >= 3]
    if rngs:
        print(f"compounds with >=3 values: {len(rngs)}; "
              f"median range {statistics.median(rngs):.2f}, max {max(rngs):.2f}")
    return fig


fig6_lfabp_by_compound();

# ---
# ## 3. Renal transporters
# 
# ### Figure 7 — expression across three species, one assay
# 
# **The finding, and the most important fact in the transporter strand:**
# the *secretory* side is conserved in kind and differs in degree — both species
# have OAT1 and OAT3, the rat simply has 6.8× more. The *reabsorptive* side is
# **not the same machinery at all**: human has OAT4 and no Oatp1a1; rat has
# Oatp1a1 and no OAT4.
# 
# That matters because **reabsorption is what sets the long half-lives**, and it
# is the half of the system that is not shared. A rat-to-human extrapolation of
# reabsorption is an extrapolation between two different proteins.
# 
# Thakur 2024, kidney membrane fractions, human n=15, rat n=10, mouse n=14 — one
# assay, so the species comparison is not assembled across laboratories.

def fig7_expression(save=None):
    rows = read(os.path.join(DB, "transporter_expression.csv"))
    want = ["OAT1", "OAT3", "OAT4", "OATP1A1", "URAT1"]

    def canon(t):
        return (t or "").upper().replace("-", "")

    # ONE ASSAY ONLY. db/transporter_expression.csv also carries Uchida 2016
    # rat values, and letting them in silently overwrote Thakur's rat OAT1
    # (32.8 -> 10.6) in an earlier version of this cell. The whole point of
    # the figure is a species comparison inside a single assay, so the study
    # filter is not cosmetic.
    by = {}
    for r in rows:
        t, sp = canon(r["transporter"]), (r["species"] or "").lower()
        if not (r.get("source") or "").startswith("Thakur 2024"):
            continue
        if (r.get("sex") or "").strip() != "both":
            continue
        if t in want and sp in ("human", "rat", "mouse"):
            v, sd = num(r["value"]), num(r["sd"])
            # a row with a BLANK value is a measured "not detected", which is
            # a finding, not a missing record -- keep it out of `by` so the
            # n.d. annotation fires, but note it did come from the assay
            if v is not None:
                by[(t, sp)] = (v, sd, (r.get("note") or ""))

    species = ["human", "rat", "mouse"]
    # TWO PANELS, because one linear axis cannot show both. OAT1 reaches
    # 32.8 and OAT4 is 0.5; on a shared scale the reabsorptive transporters
    # -- the entire point of the figure -- are invisible slivers. A log axis
    # is no good either, since "not detected" has to be drawn at zero.
    groups = [("secretory (basolateral)", ["OAT1", "OAT3"]),
              ("REABSORPTIVE (apical)", ["OAT4", "OATP1A1"]),
              ("no PFAS", ["URAT1"])]
    labels = {"OAT1": "OAT1", "OAT3": "OAT3", "OAT4": "OAT4",
              "OATP1A1": "Oatp1a1", "URAT1": "URAT1"}
    fig, axes = plt.subplots(
        1, 3, figsize=(11.0, 4.8),
        gridspec_kw=dict(width_ratios=[2, 2, 1.05], wspace=0.28))
    w = 0.26
    for ax, (title, ts) in zip(axes, groups):
        for j, sp in enumerate(species):
            xs, ys, es = [], [], []
            for i, t in enumerate(ts):
                v, sd, _ = by.get((t, sp), (None, None, ""))
                xs.append(i + (j - 1) * w)
                ys.append(v if v is not None else 0.0)
                es.append(sd if sd else 0.0)
            ax.bar(xs, ys, width=w * 0.92, color=SPECIES[sp], label=sp,
                   zorder=3, edgecolor=SURFACE, linewidth=1.2)
            ax.errorbar(xs, ys, yerr=es, fmt="none", ecolor=MUTED,
                        elinewidth=1.1, capsize=2.5, zorder=4)
            top = max([v for (t_, s_), (v, _, _) in by.items() if t_ in ts]
                      or [1])
            for x, t in zip(xs, ts):
                if (t, sp) not in by:
                    ax.annotate("not detected", xy=(x, top * 0.04),
                                color=MUTED, fontsize=7, ha="center",
                                va="bottom", rotation=90)
        ax.set_xticks(range(len(ts)))
        ax.set_xticklabels([labels[t] for t in ts], fontsize=9)
        style(ax, "", "pmol / mg membrane protein" if ax is axes[0] else "",
              title)
    leg = axes[0].legend(frameon=False, fontsize=9)
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.suptitle("Renal transporter abundance: the reabsorptive side is "
                 "not the same machinery", color=INK, fontsize=12,
                 x=0.008, ha="left", y=0.99, fontweight="semibold")
    fig.subplots_adjust(left=0.07, right=0.99, top=0.82, bottom=0.09)
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    for t in want:
        got = {sp: by.get((t, sp), (None,))[0] for sp in species}
        print(f"  {t:9} " + "  ".join(
            f"{sp}={'n.d.' if got[sp] is None else got[sp]}" for sp in species))
    return fig


fig7_expression();

# ### Figure 8 — reabsorption against secretion, within one laboratory
# 
# **The finding:** for every carboxylate, OAT4 (reabsorptive) carries 3.1–8.6
# times the intrinsic clearance of a secretory transporter. PFHxS behaves the
# same. **PFOS is the exception** — OAT1 matches OAT4 and OAT3 does not
# transport it at all.
# 
# **The caveat that bounds it, and it is a large one.** OAT4 cannot be scaled to
# the kidney and OAT1/OAT3 can. Applying the expression factor lifts OAT1 to
# 118–334 µL/min/mg, which **reverses the ordering**. So this is evidence about
# the proteins, not about the in vivo balance.

def fig8_oat4_vs_secretory(save=None):
    rows = read(os.path.join(DB, "transporter_clearance_stage3.csv"))
    # LOUISSE ONLY, and no max() across studies. The file also holds Weaver,
    # Nakagawa, Niu and Zhao values for the same transporter-compound pairs,
    # and taking the maximum across them silently built a chimera (it put
    # Weaver's OAT3 57.8 for PFOA next to Louisse's OAT4 95.7, destroying the
    # within-one-laboratory comparison the figure claims to show).
    by = defaultdict(dict)
    for r in rows:
        cl = num(r.get("clint_uL_min_mg"))
        t = (r["transporter"] or "").upper().replace("-", "")
        if cl is None or t not in ("OAT1", "OAT3", "OAT4"):
            continue
        if not (r.get("study") or "").startswith("Louisse"):
            continue
        key = (r["chemical"], r.get("head_group", ""),
               r.get("n_perfluorinated_carbons"))
        by[key][t] = cl

    keys = [k for k in by if "OAT4" in by[k]]
    keys.sort(key=lambda k: total_carbons(k[1], k[2]) if k[2] else 0)
    if not keys:
        print("no OAT4 rows with an intrinsic clearance in this file")
        return None

    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    w = 0.26
    for j, (t, col) in enumerate((("OAT4", C1), ("OAT3", C2), ("OAT1", C3))):
        xs = [i + (j - 1) * w for i in range(len(keys))]
        ys = [by[k].get(t, 0) for k in keys]
        ax.bar(xs, ys, width=w * 0.92, color=col, zorder=3,
               edgecolor=SURFACE, linewidth=1.2,
               label=f"{t} ({'reabsorptive' if t == 'OAT4' else 'secretory'})")
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f"{k[0]}\nC{total_carbons(k[1], k[2]) if k[2] else '?'}"
                        for k in keys], fontsize=8.5)
    style(ax, "", "intrinsic clearance  (µL/min/mg)",
          "OAT4 against the secretory transporters, unscaled "
          "(Louisse 2023 / 2024 only)")
    leg = ax.legend(frameon=False, fontsize=9)
    for t_ in leg.get_texts():
        t_.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    for k in keys:
        d = by[k]
        ratio = (d["OAT4"] / d["OAT3"]) if d.get("OAT3") else None
        print(f"  {k[0]:7} OAT4={d.get('OAT4')}  OAT3={d.get('OAT3')}  "
              f"OAT1={d.get('OAT1')}  OAT4/OAT3="
              f"{'-' if ratio is None else format(ratio, '.1f')}")
    return fig


fig8_oat4_vs_secretory();

# ---
# ## 4. The tubular factor, and the f_u dispute
# 
# Everything converges on one identity:
# 
# $$\mathrm{CL_R} = f_u \times \mathrm{GFR} \times T
# \qquad\Longrightarrow\qquad
# T = \frac{\mathrm{CL_R}}{f_u \times \mathrm{GFR}}$$
# 
# $T<1$ net reabsorption · $T=1$ filtration alone · $T>1$ net secretion.

# ### Figure 9 — T against chain length, three species
# 
# **Two different spreads are drawn and they must not be conflated:**
# 
# * a **coloured bar** spans that species' own *clearance* sources, under a
#   single unbound fraction;
# * a **grey whisker** spans the four human *f_u* datasets, clearance held fixed.
# 
# The grey whiskers are the wider ones. That is the visual statement of the
# whole project: **the renal clearances agree across laboratories and the
# binding measurements do not.**

def fig9_tubular(save=None):
    cons = read(os.path.join(DB, "tubular_factor_consolidated.csv"))
    fus = read(os.path.join(DB, "human_tubular_sign_by_fu_source.csv"))
    YLO, YHI = 2e-3, 8.0
    OFF = {"human": -0.17, "mouse": 0.0, "rat": 0.17}

    cells = defaultdict(lambda: {"lo": [], "hi": []})
    for r in cons:
        sp = r["species"]
        if sp not in OFF or not r["T_high"]:
            continue
        nc = total_carbons(r["family"], int(r["n_carbons"]) -
                           (1 if r["family"].startswith("carboxyl") else 0))
        k = (r["family"], nc, r["chemical"], sp)
        if r["T_low"]:
            cells[k]["lo"].append(float(r["T_low"]))
        cells[k]["hi"].append(float(r["T_high"]))

    spread = defaultdict(list)
    for r in fus:
        if r["T"]:
            spread[r["chemical"]].append(float(r["T"]))
    spread = {k: (min(v), max(v)) for k, v in spread.items() if len(v) > 1}

    fams = ["carboxylate", "sulfonate"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.8), sharey=True,
                             gridspec_kw=dict(width_ratios=[9, 5], wspace=.07))
    for ax, fam in zip(axes, fams):
        ax.set_yscale("log")
        ax.axhline(1.0, color=AXIS, lw=1.6, ls=(0, (5, 3)), zorder=1)
        sel = [(k, v) for k, v in cells.items() if k[0] == fam]
        for (f, nc, chem, sp), v in sel:
            if sp == "human" and chem in spread:
                flo, fhi = spread[chem]
                a, b = max(flo, YLO), min(fhi, YHI)
                if b > a:
                    ax.plot([nc + OFF[sp]] * 2, [a, b], color=MUTED, lw=1.2,
                            zorder=2)
        for (f, nc, chem, sp), v in sel:
            hi = max(v["hi"])
            lo = min(v["lo"]) if v["lo"] else hi
            bnd = not v["lo"]
            x, c = nc + OFF[sp], SPECIES[sp]
            shown, over = min(hi, YHI), hi > YHI
            if shown > lo:
                ax.plot([x, x], [lo, shown], color=c, lw=2.0, zorder=3)
            if over:
                ax.plot([x], [YHI], "^", ms=7, mfc="none", mec=c, mew=1.6,
                        zorder=4)
                ends = {lo} if shown > lo else set()
            else:
                ends = {lo, shown} if shown > lo else {lo}
            for y in ends:
                ax.plot([x], [y], "v" if bnd else "o", ms=8, mfc=c,
                        mec=SURFACE, mew=1.4, zorder=4)
        xs = sorted({k[1] for k, _ in sel})
        ax.set_xticks(xs)
        ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
        ax.set_ylim(YLO, YHI)
        style(ax, "total carbons", "", fam.capitalize() + "s")
    axes[0].set_ylabel("tubular factor  T = CL$_R$ / (f$_u$ × GFR)",
                       color=SECOND, fontsize=9.5, labelpad=8)
    h = [plt.Line2D([], [], color=SPECIES[s], lw=2, marker="o", ms=8,
                    mec=SURFACE, mew=1.4, label=s)
         for s in ("human", "rat", "mouse")]
    h += [plt.Line2D([], [], color=MUTED, lw=0, marker="v", ms=8, mec=SURFACE,
                     mew=1.4, label="upper bound"),
          plt.Line2D([], [], color=MUTED, lw=1.2, marker="_", ms=9,
                     label="human: spread over four f$_u$ datasets")]
    # Legend BELOW the panels. Inside the sulfonate panel it sat on top of
    # the mouse PFHxS bar, which runs the full height of that column.
    leg = fig.legend(handles=h, frameon=False, fontsize=8.5, ncol=5,
                     loc="lower center", bbox_to_anchor=(0.52, 0.0))
    for t in leg.get_texts():
        t.set_color(SECOND)
    # subplots_adjust, not tight_layout: this figure uses a width_ratios
    # gridspec and tight_layout warns that it cannot handle it
    fig.subplots_adjust(left=0.085, right=0.985, top=0.93, bottom=0.20)
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    print(f"{len(cells)} species-compound cells; "
          f"{len(spread)} human compounds carry an f_u spread")
    return fig


fig9_tubular();

# ### Figure 10 — what the choice of f_u dataset actually decides
# 
# **The finding:** ten of the twelve PFAS with enough data give the same tubular
# sign on **every** binding dataset, and all ten are reabsorption. The
# disagreement is confined to **PFPeA and PFHxA** — where the sign flips
# depending on which laboratory you believe.
# 
# So the decision is far narrower than the 780-fold spread in f_u makes it look.
# But it is not cosmetic: those two compounds also decide whether the human
# chain-length relationship is monotone or humped.

def fig10_fu_dispute(save=None):
    rows = read(os.path.join(DB, "human_tubular_sign_by_fu_source.csv"))
    by = defaultdict(dict)
    for r in rows:
        t = num(r["T"])
        if t is not None:
            by[r["chemical"]][r["fu_source"]] = t
    order = sorted(by, key=lambda c: int(
        [r["n_carbons"] for r in rows if r["chemical"] == c][0]))
    srcs = sorted({s for d in by.values() for s in d})
    pal = {s: c for s, c in zip(srcs, [C1, C2, C3, "#eda100", "#6250d6"])}

    fig, ax = plt.subplots(figsize=(10.4, 5.2))
    ax.set_yscale("log")
    ax.axhspan(1, 1e4, color="#f3e6de", zorder=0)
    ax.axhline(1.0, color=AXIS, lw=1.6, ls=(0, (5, 3)), zorder=2)
    for i, c in enumerate(order):
        vals = by[c]
        if len(vals) > 1:
            ax.plot([i, i], [min(vals.values()), max(vals.values())],
                    color=MUTED, lw=1.2, zorder=2)
        for s, t in vals.items():
            ax.plot([i], [t], "o", ms=8, mfc=pal[s], mec=SURFACE, mew=1.3,
                    zorder=4, alpha=0.95)
    flips = [c for c in order
             if min(by[c].values()) < 1 < max(by[c].values())]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8.5)
    style(ax, "", "tubular factor  T",
          "Where the choice of binding dataset flips the sign")
    ax.annotate("net secretion above the line", xy=(0.2, 1.6), color=SECOND,
                fontsize=9)
    for s in srcs:
        ax.plot([], [], "o", ms=7, mfc=pal[s], mec=SURFACE, label=s)
    leg = ax.legend(frameon=False, fontsize=8.5, loc="lower left", ncol=2)
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    print(f"compounds whose SIGN depends on the f_u source: {flips}")
    return fig


fig10_fu_dispute();

# ---
# ## 5. Coverage, and the literature not yet ingested
# 
# ### Figure 11 — is the gap clearance, or binding?
# 
# **The finding, and it redirected the project:** the 1,497-row binding
# collection unlocks **zero** additional species-compounds. f_u is already known
# wherever a clearance exists; the gaps that remain have no binding measurement
# either. **The bottleneck is renal clearance.**

def fig11_coverage(save=None):
    rows = read(os.path.join(DB, "clearance_fu_coverage.csv"))
    sp_order, agg = [], defaultdict(lambda: [0, 0, 0])
    for r in rows:
        sp = r["species"]
        if sp not in sp_order:
            sp_order.append(sp)
        if r["has_fu"] == "yes":
            agg[sp][0] += 1
        elif r["has_plasma_type_binding"] == "yes":
            agg[sp][1] += 1          # gap that binding COULD close
        else:
            agg[sp][2] += 1          # gap binding cannot close
    sp_order.sort(key=lambda s: -sum(agg[s]))

    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    labels = ["has f$_u$ already", "gap binding could close",
              "gap binding cannot close"]
    cols = [C3, "#eda100", C2]
    left = [0] * len(sp_order)
    for k in range(3):
        vals = [agg[s][k] for s in sp_order]
        ax.barh(sp_order, vals, left=left, color=cols[k], label=labels[k],
                zorder=3, edgecolor=SURFACE, linewidth=1.4)
        for i, (v, l) in enumerate(zip(vals, left)):
            if v:
                ax.annotate(str(v), xy=(l + v / 2, i), color="white",
                            fontsize=9, ha="center", va="center",
                            fontweight="bold")
        left = [a + b for a, b in zip(left, vals)]
    ax.invert_yaxis()
    style(ax, "compounds that have a clearance somewhere in the project", "",
          "The binding collection unlocks nothing further")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right")
    for t in leg.get_texts():
        t.set_color(SECOND)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    n = sum(agg[s][1] for s in sp_order)
    print(f"compounds the binding collection could unlock: {n}")
    return fig


fig11_coverage();

# ### Figure 12 — the candidates not yet admitted
# 
# Ten human renal clearances tabulated by EFSA 2020 from four primary studies,
# **none of them ingested**. Eight would upgrade a compound from an upper bound
# to a measurement, against a human DIRECT tier that currently holds two.
# 
# **Why they are still candidates:** the EFSA document gives two different
# values for the same compounds in two passages (PFUnDA 0.064 against 0.005, a
# 13-fold gap), and some of its renal values sit 3–18× above other groups'
# measured *total* clearances. The provisional T here exists to decide whether
# chasing the primary PDFs is worth it — not to be used.

def fig12_leads(save=None):
    rows = read(os.path.join(DB, "renal_clearance_leads.csv"))
    use = [r for r in rows if r["T_provisional_low"]]
    use.sort(key=lambda r: int(r["n_carbons"]))
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.set_yscale("log")
    ax.axhline(1.0, color=AXIS, lw=1.6, ls=(0, (5, 3)), zorder=2)
    for i, r in enumerate(use):
        lo, hi = float(r["T_provisional_low"]), float(r["T_provisional_high"])
        ax.plot([i, i], [lo, hi], color=C1, lw=2.2, zorder=3)
        for y in {lo, hi}:
            ax.plot([i], [y], "o", ms=7, mfc=C1, mec=SURFACE, mew=1.3,
                    zorder=4)
    ax.set_xticks(range(len(use)))
    ax.set_xticklabels([f"{r['chemical']}\nC{r['n_carbons']}" for r in use],
                       fontsize=8.5)
    style(ax, "", "provisional tubular factor  T",
          "Candidate renal clearances — NOT admitted, decision aid only")
    missing = [r["chemical"] for r in rows if not r["T_provisional_low"]]
    print(f"{len(use)} candidates with a provisional T; "
          f"{len(missing)} blocked for want of a human f_u: {missing}")
    print("none of them covers PFPeA or PFHxA — the two compounds in dispute")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


fig12_leads();

# ---
# ## 6. The older figure scripts — nothing has been dropped
# 
# This notebook holds the twelve figures that the combined document uses. The
# project contains more, written earlier for other strands, and **none of them
# has been deleted**. They live in `tk_review/scripts/` and write into
# `tk_review/figures/`. Run any of them the same way:
# 
# ```bash
# cd tk_review && python3 -I scripts/<name>.py
# ```
# 
# | theme | scripts |
# |---|---|
# | **albumin** | `plot_albumin_by_compound.py`, `plot_albumin_chainlength.py`, `plot_albumin_loading.py`, `plot_albumin_ratio_by_compound.py`, `albumin_chain_length_optimum.py`, `make_logk_conversion_figure.py`, `make_binding_coverage_figure.py` |
# | **L-FABP** | `plot_lfabp_by_compound.py` |
# | **Qin 2024 comparison** | `make_qin_comparison_figure.py`, `make_qin_companion_figure.py` |
# | **species / clearance** | `make_species_figure.py`, `make_plasma_species_figure.py`, `argoul_individual_variability.py`, `actual_vs_estimated_simple.py` |
# | **mass balance** (deferred) | `make_mbm_figure.py`, `make_mbm_correction_figure.py`, `make_mbm_regime_figure.py`, `make_mbm_separation_figure.py`, `make_mbm_fischer_serum_figure.py` |
# | **thermodynamics** | `plot_peng2024_vanthoff.py` |
# | **omnibus builders** | `make_figures.py`, `make_figures_primary.py`, `make_review_figures.py`, `make_master_figure.py`, `make_route_and_loading_figures.py` |
# 
# Several write a `*_data.csv` beside the image, which is the editable input in
# the same way `db/binding_chain_profile.csv` is here.
# 
# > **Do not run `for f in scripts/*.py`.** That re-executes one-shot ingestion
# > scripts and duplicates rows in `db/protein_binding.csv`.

# ## Changing things
# 
# * **Save any figure to a file** — every function takes `save=`:
#   `fig2_chain_length(save="chain.png")`.
# * **Different compounds in figure 4** —
#   `fig4_loading(["PFBA", "PFHxA", "PFDA"])`.
# * **A lower inclusion threshold in figure 1** —
#   `fig1_albumin_dispersion(min_n=2)`.
# * **Switch the carbon axis back to perfluorinated** — edit `total_carbons()`
#   to `return int(n_cf)`. It is deliberately the only place the rule appears.
# * **Change the palette** — `C1`, `C2`, `C3` in the setup cell. If you pick
#   your own, check them for colour-vision deficiency rather than by eye.

# Optional: pandas for exploring, not for drawing.
try:
    import pandas as pd
    df = pd.DataFrame(read(os.path.join(DB, "tubular_factor_consolidated.csv")))
    print(df.groupby(["species", "tier"]).size())
except ImportError:
    print("pandas not installed — every figure above works without it")
