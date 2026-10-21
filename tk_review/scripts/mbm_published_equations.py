#!/usr/bin/env python3
"""
Qin 2023 and Armitage 2021 IV-MBM EQP v2.0, implemented from the PUBLISHED
equations and the papers' own supplementary parameters.

This supersedes the httk-derived parameterisation in mbm_qin_vs_armitage.py.
Armitage's Supplementary Section S1 is now in hand and differs from EPA's httk
implementation in four ways that matter:

  httk R/armitage.R                     Armitage 2021 SI S1 (used here)
  ------------------------------------  ------------------------------------
  log K_SaW = 0.71 logKow + 0.42        PIECEWISE is the default:
    (piecewise only if option.kbsa2)      < 4.5: 1.08 logKow - 0.70
                                          >= 4.5: 0.37 logKow + 2.56
  d_OW = 3.5, d_PlW = 3.5               d_OW = 3.1, d_PlW = 3.1
  one generic cell: 0.01/0.04/0.05      CELL-LINE-SPECIFIC (Tables S3-S6).
                                          The HepG2-derived line AREc2 is
                                          0.00046 / 0.00414 / 0.0514 / 0.944
  anionic phospholipids + lysosomes     neither appears in the published model

The cell-composition point is the big one: Armitage's PFAS cell term is ~99%
membrane lipid, and the hepatic line's membrane lipid is 0.00414, not 0.04.

Qin's equation, as given in the project notes and confirmed by his text:

  C_free = C_nom x V_tot / ( V_w
                           + D_BSA/w        x V_protein,medium
                           + D_lip/w        x V_lipid,medium
                           + D_SP/w         x V_protein,cell
                           + D_lip/w        x V_lipid,cell )

with D_SP/w = D_non-specific: "D_cell/w was predicted by the MBM (Eq. 17) with
D_non-specific (Table 1) as a surrogate for partitioning to structural proteins
and D_lip/w (Table 2) for lipid". There is NO cell-water term -- "we considered
only proteins and lipids of the cell homogenate as the main sorptive phase for
D_cell/w, which is different from Henneberger et al. where the cells were
recognized as an entity containing water, protein and lipid".

That single difference is most of the story: Armitage's cell contains water as
a sorptive phase and Qin's does not, and Armitage's cell protein runs on
octanol while Qin's runs on a measured non-specific protein affinity.

Run:  python3 tk_review/scripts/mbm_published_equations.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_qin_vs_armitage as Q   # compounds, Qin Tables 1-2, extrapolation rules

DB = os.path.join(os.path.dirname(HERE), "db")
MW = {c[0]: c[3] for c in Q.COMPOUNDS}
NCF = {c[0]: c[2] for c in Q.COMPOUNDS}

# ---------------------------------------------------------------- Armitage
SI = dict(
    d_OW=3.1, d_MW=1.0, d_SaW_acid=0.0, d_PlW=3.1,    # Table S1
    conc_ser_alb=24.0,          # g/L albumin in serum
    conc_ser_lip=1.9,           # g/L lipid in serum
    v_alb_partial=0.733,        # mL/g
    csalt=0.15, pH=7.4, cell_pH=7.4,
)
# Cell composition by line, Armitage SI Tables S3-S6.
# AREc2 is HepG2-derived and is the closest published match for HepG2/HepaRG.
CELL_COMP = {
    "AREc2":   dict(sl=0.00046, ml=0.00414, prot=0.0514, water=0.944, ng_cell=10.0),
    "PPARg":   dict(sl=0.001,   ml=0.009,   prot=0.08,   water=0.91,  ng_cell=7.3),
    "generic": dict(sl=0.01,    ml=0.04,    prot=0.055,  water=0.895, ng_cell=1.0),
}
# Octanol-equivalence of bulk cell protein.  This is httk's value.  The
# Armitage 2021 SI is in hand (literature/fulltext/) and confirms every OTHER
# parameter below verbatim -- both spLFERs, all four scaling factors (Table S1)
# and the AREc2/PPARg compositions (Table S5) -- but states nothing about how
# bulk protein sorbs.  The bulk-protein phase is new in v2.0 ("expanded to
# include ... bulk proteins") so the 2014 SI cannot hold it either.  It exists
# only in the Excel/VBA tool, by request.  Findings that turn on this constant
# are findings about httk's implementation.  See MBM_Qin_vs_Armitage.md section 14.
P_NLOM = 0.035


def armitage_K(logKow_n, piecewise=True):
    """spLFERs, Armitage 2021 SI S1 Eqs. 2-3 and main-text Eqs. 5-6."""
    gKmw = 1.01 * logKow_n + 0.12
    if piecewise:
        gKsa = (1.08 * logKow_n - 0.70) if logKow_n < 4.5 else (0.37 * logKow_n + 2.56)
    else:
        gKsa = 0.71 * logKow_n + 0.42          # httk default, for comparison
    gKpl = 0.97 * logKow_n - 6.94              # Kramer, Eq. 5
    ksalt = 0.04 * logKow_n + 0.114
    return gKmw, gKsa, gKpl, ksalt


def armitage_D(logKow_n, comp, piecewise=True, Fneutral=0.0):
    """Distribution ratios for an acidic IOC at pH 7.4. Eq. 2 terms."""
    Fn, Fi = Fneutral, 1.0 - Fneutral
    gKmw, gKsa, gKpl, ksalt = armitage_K(logKow_n, piecewise)
    salt = 10 ** (-ksalt * SI["csalt"])
    Dow = (Fn * 10 ** logKow_n + Fi * 10 ** (logKow_n - SI["d_OW"])) / salt
    Dmw = (Fn * 10 ** gKmw + Fi * 10 ** (gKmw - SI["d_MW"])) / salt
    Dsa = (Fn * 10 ** gKsa + Fi * 10 ** (gKsa - SI["d_SaW_acid"])) / salt
    Dpl = (Fn * 10 ** gKpl + Fi * 10 ** (gKpl - SI["d_PlW"])) / salt
    ioc = 10 ** (SI["cell_pH"] - SI["pH"])     # Eq. 3; 1.0 when both are 7.4
    Dcw = (comp["sl"] * Dow + comp["ml"] * Dmw
           + P_NLOM * comp["prot"] * Dow + comp["water"]) * ioc
    return dict(Dow=Dow, Dmw=Dmw, Dsa=Dsa, Dpl=Dpl, Dcw=Dcw,
                gKsa=gKsa, gKmw=gKmw)


# ---------------------------------------------------------------- Qin
def qin_D(name, Cw):
    """Qin's three distribution ratios at aqueous concentration Cw (uM)."""
    p = Q.qin_params(name, NCF[name])
    return (p["Cmax"] / (p["Kd"] + Cw) + 10 ** p["logDnonspec"],   # D_BSA/w, Eq.16
            10 ** p["logDnonspec"],                                # D_SP/w  (= D_ns)
            10 ** p["logDlipw"],                                   # D_lip/w
            p)


# ---------------------------------------------------------------- geometry
def well_sarea(diam_mm, v_uL, n=1):
    r = diam_mm / 20.0
    base = math.pi * r ** 2
    h = (v_uL / n / 1000.0) / base
    return n * (base + 2 * math.pi * r * h) * 1e-4      # m2


ASSAYS = {
    # Rosenmai 2017: HepG2, 12-well, 18e4 seeded, 3 wells pooled, 10% FBS.
    # Working volume not stated; 1 mL/well is the standard and the internal
    # ng/mg reconciliation implies ~0.7 mL/well.
    "Rosenmai": dict(line="AREc2", diam_mm=22.1, wells=3, v_uL=3000.0,
                     cells=5.434e5, FBS=0.10, doses=[10.0, 30.0, 100.0]),
    # Amstutz 2022 (Toxicology 480:153312) S2.2/2.3/2.5: HepG2, 96-well, 4e4
    # cells, 250 uL dosed, exposed in the GROWTH medium -- DMEM + 10% FBS,
    # 0.1% DMSO carrier.  FBS verified at source.
    "Amstutz": dict(line="AREc2", diam_mm=6.4, wells=1, v_uL=250.0,
                    cells=4.0e4, FBS=0.10, doses=[100.0]),
    # Fragki 2022 (Food Chem Toxicol 172:113559) S2.2/2.3: HepaRG, 24-well,
    # seeded 55,000/well in 500 uL, dosed in 550 uL.  The growth medium is
    # William's E + 10% FBS, but exposure is in "assay medium (growth medium
    # containing 2% FBS) supplemented with 0.5% DMSO" -- FBS = 0.02, NOT 0.10.
    # This was 0.10 until mbm_medium_composition.py traced it to the source;
    # the 5x matters because albumin carries 79-97% of the denominator here.
    # Differentiated HepaRG reach roughly 4x the seeding density.
    "Fragki": dict(line="AREc2", diam_mm=15.6, wells=1, v_uL=550.0,
                   cells=2.2e5, FBS=0.02,
                   doses=[6.25, 12.5, 25.0, 50.0, 100.0]),
}


def volumes(a):
    comp = CELL_COMP[a["line"]]
    Vbm = a["v_uL"] * 1e-6                                   # L
    Vcells = a["cells"] * comp["ng_cell"] * 1e-6 * 1e-6      # ng -> mg -> L at 1 mg/uL
    Valb = Vbm * a["FBS"] * SI["v_alb_partial"] * SI["conc_ser_alb"] * 1e-3
    Vslip = Vbm * a["FBS"] * SI["conc_ser_lip"] * 1e-3
    Vw = Vbm - Valb - Vslip
    return dict(Vbm=Vbm, Vcells=Vcells, Valb=Valb, Vslip=Vslip, Vw=Vw,
                Vprot_cell=Vcells * comp["prot"],
                Vlip_cell=Vcells * (comp["sl"] + comp["ml"]),
                sarea=well_sarea(a["diam_mm"], a["v_uL"], a["wells"]),
                comp=comp)


def mg_protein_per_mL_cell(a, protein_density=1.35):
    return CELL_COMP[a["line"]]["prot"] * protein_density * 1000.0


# ---------------------------------------------------------------- balances
def armitage_run(name, logKow_n, dose, V, piecewise=True):
    D = armitage_D(logKow_n, V["comp"], piecewise)
    ntot = dose * V["Vbm"]
    den = (V["Vw"] + D["Dsa"] * V["Valb"] + D["Dow"] * V["Vslip"]
           + D["Dcw"] * V["Vcells"] + D["Dpl"] * V["sarea"] * 1000.0)
    Cw = ntot / den
    Ccell = D["Dcw"] * Cw
    return dict(Cw=Cw, C_cell=Ccell, D=D,
                x_cell=Ccell * V["Vcells"] / ntot,
                x_alb=D["Dsa"] * Cw * V["Valb"] / ntot,
                x_plastic=D["Dpl"] * Cw * V["sarea"] * 1000.0 / ntot)


def qin_run(name, dose, V, linear=False):
    """Qin's equation, solved self-consistently for the non-linear D_BSA/w.

    C_cell is reported on a WHOLE-CELL basis (amount in cell protein + cell
    lipid, divided by total cell volume) so it is comparable with Armitage's
    and with a measured per-cell concentration.
    """
    ntot = dose * V["Vbm"]
    Cw = dose * 1e-3
    for _ in range(500):
        Dbsa, Dsp, Dlip, p = qin_D(name, Cw)
        if linear:
            Dbsa = 10 ** p["logKspec"]
        den = (V["Vw"] + Dbsa * V["Valb"] + Dlip * V["Vslip"]
               + Dsp * V["Vprot_cell"] + Dlip * V["Vlip_cell"])
        new = ntot / den
        if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-12:
            Cw = new
            break
        Cw = math.sqrt(Cw * new)
    Dbsa, Dsp, Dlip, p = qin_D(name, Cw)
    if linear:
        Dbsa = 10 ** p["logKspec"]
    n_cell = Dsp * V["Vprot_cell"] * Cw + Dlip * V["Vlip_cell"] * Cw
    return dict(Cw=Cw, C_cell=n_cell / V["Vcells"], x_cell=n_cell / ntot,
                x_alb=Dbsa * Cw * V["Valb"] / ntot,
                Dbsa=Dbsa, Dsp=Dsp, Dlip=Dlip, params=p,
                Dcw_effective=(Dsp * V["Vprot_cell"] + Dlip * V["Vlip_cell"])
                / V["Vcells"])


# ---------------------------------------------------------------- measured
ROSENMAI_UPTAKE_PCT = Q.ROSENMAI_UPTAKE_PCT          # Table 1, % of moles added
AMSTUTZ_IF_PCT = {"PFBA": None, "PFHxA": 0.039, "PFHpA": 0.044,
                  "PFOA": 0.245, "PFNA": 0.302}      # Table 5, 100 uM
FRAGKI_NMOLMG = {                                     # main paper dose series
    "PFOA": {6.25: 0.0562, 12.5: 0.169, 25.0: 0.506, 50.0: 1.24, 100.0: 2.81},
    "PFNA": {6.25: 0.281, 12.5: 0.73, 25.0: 1.63, 50.0: 4.49, 100.0: 15.8},
}
# Rosenmai Fig 2A, as plotted (axis says nmol/mg; the Abstract says ng/mg)
ROSENMAI_FIG2A = {
    "PFBA": {10: 18, 30: 45, 100: 126}, "PFHxA": {10: 26, 30: 69, 100: 164},
    "PFHpA": {10: 28, 30: 116, 100: 297}, "PFOA": {10: 39, 30: 127, 100: 525},
    "PFNA": {10: 24, 30: 86, 100: 458},
}


def measured_x_cell(study, cp, dose, V, a):
    """Fraction of nominal dose in cells -- the axis that does NOT depend on
    the assumed cell volume, only on protein per cell, which all sources agree
    on (~675-694 pg/cell)."""
    if study == "Rosenmai":
        v = ROSENMAI_UPTAKE_PCT.get(cp, {}).get(int(dose))
        return v / 100.0 if v else None
    if study == "Amstutz":
        v = AMSTUTZ_IF_PCT.get(cp)
        return v / 100.0 if v else None
    v = FRAGKI_NMOLMG.get(cp, {}).get(dose)
    if v is None:
        return None
    mg_prot = V["Vcells"] * 1000 * mg_protein_per_mL_cell(a)
    return v * mg_prot / (dose * V["Vbm"] * 1000)


def main():
    print("Armitage: SI S1 spLFERs, d_OW/d_PlW = 3.1, AREc2 cell composition")
    print("Qin: Eq. with D_SP/w = D_non-specific, cell = protein + lipid only")
    print()
    print("PROTEIN PER CELL -- the one thing every source agrees on")
    c = CELL_COMP["AREc2"]
    print("  Armitage AREc2: %.0f pg/cell (%.0f pL x %.0f mg/mL)"
          % (c["ng_cell"] * 1000 * c["prot"] * 1.35, c["ng_cell"],
             c["prot"] * 1350))
    print("  project sheet : 675 pg/cell (3.04 pL x 222 mg/mL)")
    print("  -> agree to 3% on protein/cell; disagree 3.3x on CELL VOLUME.")
    print("  -> fraction-of-dose is robust; C_cell carries that 3.3x.")
    print()
    hdr = ("%-9s %-6s %6s | %8s %8s %8s | %7s %7s"
           % ("study", "cmpd", "dose", "meas %", "Qin %", "Arm %", "Q/M", "A/M"))
    print(hdr); print("-" * len(hdr))
    rows = []
    for study, a in ASSAYS.items():
        V = volumes(a)
        for nm, nC, nCF, mw, ko, kc, ksrc in Q.COMPOUNDS:
            lk = ko if ko is not None else Q.KOW_OPERA(nCF)
            for d in a["doses"]:
                m = measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                q = qin_run(nm, d, V)
                ar = armitage_run(nm, lk, d, V)
                print("%-9s %-6s %6g | %8.3f %8.3f %8.3f | %7.1f %7.1f"
                      % (study, nm, d, m * 100, q["x_cell"] * 100,
                         ar["x_cell"] * 100, q["x_cell"] / m, ar["x_cell"] / m))
                rows.append(dict(
                    study=study, compound=nm, dose_uM=d, logKow=round(lk, 3),
                    qin_param_source=Q.qin_params(nm, nCF)["source"],
                    meas_x_cell_pct=round(m * 100, 4),
                    qin_x_cell_pct=round(q["x_cell"] * 100, 4),
                    arm_x_cell_pct=round(ar["x_cell"] * 100, 4),
                    qin_over_meas=round(q["x_cell"] / m, 2),
                    arm_over_meas=round(ar["x_cell"] / m, 2),
                    qin_Dcellw=round(q["Dcw_effective"], 2),
                    arm_Dcellw=round(ar["D"]["Dcw"], 2),
                    qin_Dbsaw=round(q["Dbsa"], 1),
                    arm_Dsaw=round(ar["D"]["Dsa"], 1),
                    qin_Cw_uM=round(q["Cw"], 4), arm_Cw_uM=round(ar["Cw"], 4),
                    qin_x_alb_pct=round(q["x_alb"] * 100, 2),
                    arm_x_alb_pct=round(ar["x_alb"] * 100, 2),
                    arm_x_plastic_pct=round(ar["x_plastic"] * 100, 3)))
    dst = os.path.join(DB, "mbm_published_equations.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print()
    print("wrote %s (%d rows)" % (dst, len(rows)))
    return rows


if __name__ == "__main__":
    main()
