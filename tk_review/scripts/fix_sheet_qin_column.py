#!/usr/bin/env python3
"""
Rebuild the project spreadsheet's Qin column with (a) the composed cell term
and (b) C_free rather than C_nominal, and plot measured vs estimated C_cell in
the spreadsheet's own layout.

The two fixes
-------------
(a) COMPOSED CELL TERM.  Qin's cell is not one lumped D_cell/w applied to a
    cell volume.  His equation splits it, and his text says what goes in each
    slot: "D_cell/w was predicted by the MBM (Eq. 17) with D_non-specific
    (Table 1) as a surrogate for partitioning to structural proteins and
    D_lip/w (Table 2) for lipid", with no cell-water term at all.

        n_cell = C_free x ( D_SP/w x V_protein,cell + D_lip/w x V_lipid,cell )
        D_SP/w = D_non-specific          (Qin Table 1)

(b) C_FREE.  Qin's distribution ratios are defined against the freely
    dissolved concentration, which the full equation solves for:

        C_free = C_nom x V_tot / ( V_w
                                 + D_BSA/w x V_protein,medium
                                 + D_lip/w x V_lipid,medium
                                 + D_SP/w  x V_protein,cell
                                 + D_lip/w x V_lipid,cell )

    D_BSA/w is concentration-dependent (Eq. 16), so this is solved
    self-consistently rather than evaluated once.

Cell protein volume is anchored on the one quantity every source agrees on --
675 pg protein per cell (Armitage's AREc2 gives 694) -- at a protein density of
1.35 g/mL, so V_protein,cell = cells x 5.0e-13 L.  Cell lipid is set from Qin's
own measured volume ratio, lipid/protein = 0.00432/0.0814.

Run:  python3 tk_review/scripts/fix_sheet_qin_column.py
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P
import mbm_qin_vs_armitage as Q

DB = os.path.join(os.path.dirname(HERE), "db")
FIGS = os.path.join(os.path.dirname(HERE), "figures")
MW = {c[0]: c[3] for c in Q.COMPOUNDS}
NCF = {c[0]: c[2] for c in Q.COMPOUNDS}
LKOW = dict((c[0], c[4]) for c in Q.COMPOUNDS)

MG_PROTEIN_PER_CELL = 675e-9          # mg
PROTEIN_DENSITY = 1.35                # g/mL
L_PROTEIN_PER_CELL = 675e-12 / PROTEIN_DENSITY / 1000.0     # 5.0e-13 L
LIP_OVER_PROT = 0.00432 / 0.0814      # Qin Table A5, by volume


def lkow(c):
    return LKOW[c] if LKOW[c] is not None else Q.KOW_OPERA(NCF[c])


def cell_volumes(cells):
    vp = cells * L_PROTEIN_PER_CELL
    return vp, vp * LIP_OVER_PROT


def qin_fixed(compound, dose_uM, study):
    """Qin's equation with the composed cell term, solved for C_free.

    Returns C_cell in nmol per mg protein -- the spreadsheet's unit.
    """
    a = P.ASSAYS[study]
    V = P.volumes(a)
    Vp, Vl = cell_volumes(a["cells"])
    mg_protein = a["cells"] * MG_PROTEIN_PER_CELL
    ntot = dose_uM * V["Vbm"]                       # umol
    Cw = dose_uM * 1e-3                             # seed
    for _ in range(500):
        Dbsa, Dsp, Dlip, p = P.qin_D(compound, Cw)
        den = (V["Vw"]
               + Dbsa * V["Valb"]                   # D_BSA/w x V_protein,medium
               + Dlip * V["Vslip"]                  # D_lip/w x V_lipid,medium
               + Dsp * Vp                           # D_SP/w  x V_protein,cell
               + Dlip * Vl)                         # D_lip/w x V_lipid,cell
        new = ntot / den
        if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-12:
            Cw = new
            break
        Cw = math.sqrt(Cw * new)
    Dbsa, Dsp, Dlip, p = P.qin_D(compound, Cw)
    n_cell_umol = Cw * (Dsp * Vp + Dlip * Vl)
    return dict(nmol_per_mg=n_cell_umol * 1000.0 / mg_protein,
                Cfree_uM=Cw, Cfree_over_Cnom=Cw / dose_uM,
                Dbsa=Dbsa, Dsp=Dsp, Dlip=Dlip,
                x_cell=n_cell_umol / ntot,
                param_source=Q.qin_params(compound, NCF[compound])["source"])


def armitage_nmol_per_mg(compound, dose_uM, study):
    a = P.ASSAYS[study]
    V = P.volumes(a)
    mg_protein = a["cells"] * MG_PROTEIN_PER_CELL
    r = P.armitage_run(compound, lkow(compound), dose_uM, V)
    return r["C_cell"] * V["Vcells"] * 1000.0 / mg_protein


def build():
    rows = []
    for r in csv.DictReader(open(os.path.join(DB, "mbm_project_sheet.csv"))):
        st, cp, d = r["study"], r["compound"], float(r["Ctest_uM"])
        q = qin_fixed(cp, d, st)
        rows.append(dict(
            study=st, compound=cp, Ctest_uM=d,
            Ccell_measured_nmolmgprotein=float(r["Ccell_meas_nmolmg"]),
            Ccell_measured_MWcorrected=round(
                float(r["Ccell_meas_nmolmg"]) / (MW[cp] if st == "Rosenmai" else 1.0), 5),
            Qin_OLD_sheet=float(r["Ccell_qin_nmolmg_sheet"]),
            Qin_FIXED=round(q["nmol_per_mg"], 5),
            Qin_fix_factor=round(float(r["Ccell_qin_nmolmg_sheet"]) / q["nmol_per_mg"], 1),
            Armitage_sheet=float(r["Ccell_arm_nmolmg_sheet"]),
            Armitage_published_eqns=round(armitage_nmol_per_mg(cp, d, st), 5),
            Cfree_uM=round(q["Cfree_uM"], 4),
            Cfree_over_Cnom=round(q["Cfree_over_Cnom"], 4),
            D_BSAw=round(q["Dbsa"], 1), D_SPw=round(q["Dsp"], 1),
            D_lipw=round(q["Dlip"], 1),
            qin_param_source=q["param_source"]))
    dst = os.path.join(DB, "mbm_sheet_qin_fixed.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return rows, dst


# ggplot default hue palette, 5 levels, in the spreadsheet's alphabetical order
PAL = {"PFBA": "#F8766D", "PFHpA": "#A3A500", "PFHxA": "#00BF7D",
       "PFNA": "#00B0F6", "PFOA": "#E76BF3"}
SHAPE = {"Amstutz": "o", "Fragki": "s", "Rosenmai": "^"}
INK, INK2, MUTED, RULE, GRID = "#0b0b0b", "#52514e", "#8a8984", "#d8d7d3", "#e6e5e1"
SURFACE = "#fcfcfb"


def plot(rows):
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
        "font.size": 9, "axes.edgecolor": RULE, "axes.linewidth": 0.8,
        "axes.labelcolor": INK2, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2,
        "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    # axis limits match the project's original plot, for direct comparison
    XLO, XHI, YLO, YHI = 2e-2, 2e3, 1e-3, 4e2
    lo, hi = 1e-4, 1e4
    ax.plot([lo, hi], [lo, hi], "-", color=INK, lw=1.1, zorder=2)
    for f, ls in ((10, "--"), (0.1, "--"), (100, ":"), (0.01, ":")):
        ax.plot([lo, hi], [lo * f, hi * f], ls, color=MUTED, lw=0.8, zorder=1)
    for r in rows:
        x = r["Ccell_measured_MWcorrected"]
        col = PAL[r["compound"]]; mk = SHAPE[r["study"]]
        ax.plot(x, r["Armitage_published_eqns"], mk, ms=7, mfc=col, mec=INK,
                mew=0.8, zorder=4)                       # filled  = Armitage
        ax.plot(x, r["Qin_FIXED"], mk, ms=7, mfc="none", mec=col,
                mew=1.6, zorder=3)                       # open    = Qin
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(XLO, XHI); ax.set_ylim(YLO, YHI)
    ax.set_xlabel("Measured C$_{cell}$  (nmol/mg protein)")
    ax.set_ylabel("Estimated C$_{cell}$  (nmol/mg protein)")
    ax.set_title("Measured vs estimated C$_{cell}$\n"
                 "Qin with composed cell term and C$_{free}$; "
                 "Armitage from SI S1.  Rosenmai ÷ MW.",
                 fontsize=10, loc="left", color=INK)
    ax.grid(True, which="major", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    h = []
    h.append(plt.Line2D([], [], ls="", marker="o", ms=7, mfc=MUTED, mec=INK,
                        mew=0.8, label="Armitage (filled)"))
    h.append(plt.Line2D([], [], ls="", marker="o", ms=7, mfc="none", mec=MUTED,
                        mew=1.6, label="Qin (open)"))
    h.append(plt.Line2D([], [], ls="", label=""))
    for c, col in PAL.items():
        h.append(plt.Line2D([], [], ls="", marker="s", ms=7, color=col, label=c))
    h.append(plt.Line2D([], [], ls="", label=""))
    for s, mk in SHAPE.items():
        h.append(plt.Line2D([], [], ls="", marker=mk, ms=7, mfc="none",
                            mec=INK2, label=s))
    ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, fontsize=8, handletextpad=0.6, labelspacing=0.55)
    ax.text(0.98, 0.04, "solid 1:1   dashed 10×   dotted 100×",
            transform=ax.transAxes, ha="right", fontsize=7, color=MUTED)
    fig.subplots_adjust(left=0.095, right=0.775, bottom=0.115, top=0.845)
    out = os.path.join(FIGS, "fig11_measured_vs_estimated_ccell.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    ROWS, DST = build()
    print("wrote", DST)
    plot(ROWS)
