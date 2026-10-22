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
                                          AREc2 (MCF7-derived) is
                                          0.00046 / 0.00414 / 0.0514 / 0.944
  anionic phospholipids + lysosomes     neither appears in the published model

The cell-composition point is the big one: Armitage's PFAS cell term is ~99%
membrane lipid, and AREc2's membrane lipid is 0.00414, not 0.04.  Note that
AREc2 is not a hepatic line -- see the CELL_COMP comment below -- and that
Fischer et al. 2017 measured 0.0319 total lipid for HepG2.

Qin's equations, now verified against the full text (Environ Int 173:107857):

  Eq. 17   D_i/w = D_BSA/w x V_prot/V_prot+lip  +  D_lip/w x V_lip/V_prot+lip
  Eq. 20   C_free = n_tot / [ (D_cell/w   x V_cell,prot+lip   + V_cell,w)
                            + (D_medium/w x V_medium,prot+lip + V_medium,w) ]

Three points of fidelity that the tabulated parameters alone did not settle:

  * The normalising volume in Eq. 17 is the SORPTIVE volume V_prot+lip, not
    the volume of the biomaterial.  Table 2 gives the units as
    [L_w / L_cell,prot+lip], confirming it.  Substituting Eq. 17 into Eq. 20
    therefore gives exactly D_SP/w x V_prot,cell + D_lip/w x V_lip,cell, which
    is the form used below -- the two are algebraically identical.
  * Cell water IS in Eq. 20, but as a non-sorbing phase at D = 1 pooled with
    the medium water, NOT inside D_cell/w.  The contrast the authors draw with
    Henneberger et al. 2019a is about the NUMERATOR: "we considered only
    proteins and lipids of the cell homogenate as the main sorptive phase for
    D_cell/w, which is different from Henneberger et al. where the cells were
    recognized as an entity containing water, protein and lipid and all of
    them contributed to cell-water distribution".  So cell water counts
    towards the denominator but not towards the amount "in the cell".
  * D_SP/w = D_non-specific: "D_cell/w was predicted by the MBM (Eq. 17) with
    D_non-specific (Table 1) as a surrogate for partitioning to structural
    proteins and D_lip/w (Table 2) for lipid" (Sec 3.5).

Qin's own cell composition (Table A5, HEK293H) is water 91.6%, protein 8.14%,
lipid 0.432% by volume -- within 7% of Armitage's AREc2 lipid but 1.6x his
protein.  The composition used below is Armitage's, for both models, so that
the two differ only in their partitioning treatment.

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

# Perfluoroalkane sulfonates.  Qin measured EVERY parameter for these two, so
# they need no interpolation at all -- they are the cleanest case in his
# dataset, cleaner than three of the five carboxylates.  They are deliberately
# absent from Q.COMPOUNDS, which is the carboxylate series that supplies the
# interpolation anchors; putting a sulfonate in there would let it anchor a
# carboxylate estimate, and sulfonates partition about a log unit more
# strongly into lipid at the same chain length.
#
# nCF counts perfluorinated carbons.  A sulfonate has no carboxyl carbon, so
# PFHxS (C6F13SO3-) is nCF 6 and PFOS (C8F17SO3-) is nCF 8, whereas the
# carboxylate PFOA (C7F15COO-) is nCF 7 at nC 8.
#
# log K_OW is from Qin Table 1.  Note the OPERA/COSMOtherm disagreement for
# PFHxS -- 3.69 against 5.29, a factor of 40 - which the Armitage limb rides
# on directly.  Both are carried so the spread can be reported.
SULFONATES = [
    # name     nC  nCF    MW       logKow_opera  logKow_cosmo
    ("PFHxS",   6,   6,  400.11,    3.69,        5.29),
    ("PFOS",    8,   8,  500.13,    5.77,        None),
]
NCF.update({c[0]: c[2] for c in SULFONATES})

# ---------------------------------------------------------------- Armitage
SI = dict(
    d_OW=3.1, d_MW=1.0, d_SaW_acid=0.0, d_PlW=3.1,    # Table S1
    conc_ser_alb=24.0,          # g/L albumin in serum
    conc_ser_lip=1.9,           # g/L lipid in serum
    v_alb_partial=0.733,        # mL/g
    csalt=0.15, pH=7.4, cell_pH=7.4,
)
# Cell composition by line, Armitage SI Tables S3-S6.
#
# CAUTION on AREc2.  Armitage's cell-line-specific entries are for the two
# lines Huchthausen et al. 2020 used: AREc2 (= AREc32) and a PPARgamma
# reporter.  AREc32 is MCF7-derived, not HepG2-derived (Wang, Hayes & Wolf
# 2006; Cellosaurus CVCL_1D32 "MCF7 AREc32"), so it is NOT a hepatic
# composition.  It is retained as the published Armitage parameterisation --
# the thing being tested -- but it is a breast-carcinoma line applied to
# HepG2 and HepaRG assays, which is a real mismatch, not a close match.
# Armitage's PPARgamma entry reproduces Fischer et al. 2017's measured HEK293
# composition (8.43% protein, 1.00% lipid, 90.57% water) to the digits shown,
# which is where it evidently comes from.
#
# HepG2_Fischer is the MEASURED composition of the actual cell line used by
# two of the three datasets here: Fischer et al. 2017 Table 2, ARE assay,
# HepG2 (ARE-bla) -- 87.35% water, 9.46% protein, 3.19% lipid, and 3.04 pg per
# cell from their own volumes (6.077 nL for 2000 cells).  Fischer assumes "all
# cellular lipids are membrane lipids", so the whole 3.19% is entered as ml.
CELL_COMP = {
    "AREc2":   dict(sl=0.00046, ml=0.00414, prot=0.0514, water=0.944, ng_cell=10.0),
    "PPARg":   dict(sl=0.001,   ml=0.009,   prot=0.08,   water=0.91,  ng_cell=7.3),
    "generic": dict(sl=0.01,    ml=0.04,    prot=0.055,  water=0.895, ng_cell=1.0),
    "HepG2_Fischer": dict(sl=0.0, ml=0.0319, prot=0.0946, water=0.8735,
                          ng_cell=3.04),
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
                Vw_cell=Vcells * comp["water"],      # Qin Eq. 20, V_cell,w
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
    and with a measured per-cell concentration.  Note that this is NOT the
    basis of Qin's own tabulated log D_cell/w, which is per litre of cell
    protein+lipid; Dcw_effective below is the whole-cell rebasing, and differs
    from Qin's Table 2 convention by the factor V_prot+lip/V_cell (0.060 for
    the AREc2 composition, i.e. 1.22 log units).  mbm_qin_internal_check.py
    does the comparison on Qin's own basis.

    The V_cell,w term of Eq. 20 is included in the denominator.  It is a
    non-sorbing phase at D = 1 and shifts C_free by under 0.2% at these cell
    densities, but it is in the published equation.
    """
    ntot = dose * V["Vbm"]
    Cw = dose * 1e-3
    for _ in range(500):
        Dbsa, Dsp, Dlip, p = qin_D(name, Cw)
        if linear:
            Dbsa = 10 ** p["logKspec"]
        den = (V["Vw"] + V["Vw_cell"] + Dbsa * V["Valb"] + Dlip * V["Vslip"]
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
# Fragki 2022 Table 3: "Cellular uptake (%) of the four PFASs ... Mean values
# of three individual experiments (three replicates each)".
#
# PREFER THIS over the nmol/mg series above.  Converting nmol/mg protein to a
# fraction of dose needs the protein per well, which the paper does not state
# and which this file was estimating as 2.2e5 cells x 694 pg (0.153 mg).  The
# paper's own two numbers pin it instead: protein = uptake_fraction x n_total
# / (nmol/mg) gives 0.30-0.34 mg per well from the PFOA series (CV 6%), i.e.
# about 2.1x the estimate, implying ~4.8e5 differentiated cells per 24-well
# rather than the 4x-seeding-density guess of 2.2e5.  Table 3 is the fraction
# of dose directly, so it needs no protein assumption at all and is the
# quantity both models predict.  Using it drops the Fragki bias from 10.80x
# to 6.66x (Armitage) and 14.91x to 9.20x (Qin).
#
# The PFNA nmol/mg series is not self-consistent with Table 3 to the same
# precision (it implies 0.14-0.22 mg per well), which is a further reason to
# take the directly reported percentages.
FRAGKI_UPTAKE_PCT = {
    "PFOA":  {6.25: 0.55, 12.5: 0.84, 25.0: 1.23, 50.0: 1.34, 100.0: 1.55},
    "PFNA":  {6.25: 1.12, 12.5: 2.34, 25.0: 2.65, 50.0: 3.46, 100.0: 5.00},
    "PFHxS": {6.25: 0.40, 12.5: 1.23, 25.0: 1.33, 50.0: 1.49, 100.0: 1.64},
    "PFOS":  {6.25: 3.87, 12.5: 4.52, 25.0: 5.91, 50.0: 11.22, 100.0: 10.94},
}
# Rosenmai Fig 2A, as plotted (axis says nmol/mg; the Abstract says ng/mg)
ROSENMAI_FIG2A = {
    "PFBA": {10: 18, 30: 45, 100: 126}, "PFHxA": {10: 26, 30: 69, 100: 164},
    "PFHpA": {10: 28, 30: 116, 100: 297}, "PFOA": {10: 39, 30: 127, 100: 525},
    "PFNA": {10: 24, 30: 86, 100: 458},
}


def measured_x_cell(study, cp, dose, V, a):
    """Fraction of nominal dose in cells.

    This axis does not depend on the assumed cell volume.  Where a study
    reports the fraction directly (Rosenmai, Amstutz, and Fragki Table 3) it
    is used as reported and no cell parameter enters at all.  Only the
    nmol/mg fallback needs protein per well, and the Fragki case shows why
    that route is worth avoiding: it was out by 2.1x.
    """
    if study == "Rosenmai":
        v = ROSENMAI_UPTAKE_PCT.get(cp, {}).get(int(dose))
        return v / 100.0 if v else None
    if study == "Amstutz":
        v = AMSTUTZ_IF_PCT.get(cp)
        return v / 100.0 if v else None
    v = FRAGKI_UPTAKE_PCT.get(cp, {}).get(dose)
    if v is not None:                      # the paper's own fraction of dose
        return v / 100.0
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
