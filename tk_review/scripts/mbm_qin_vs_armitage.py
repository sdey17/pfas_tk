#!/usr/bin/env python3
"""
Qin (2023) vs Armitage (2021, IV-MBM EQP v2.0) in-vitro mass-balance models
run head-to-head on PFBA, PFHxA, PFHpA, PFOA, PFNA, against measured
cell-associated PFAS from Rosenmai 2017, Amstutz 2022 and Fragki 2023.

Why this script exists
----------------------
The two models disagree, and the question is where. They differ in four
structural ways, each isolated below so its contribution can be read off:

  (1) ALBUMIN.   Qin *measures* D_BSA/w and it is NON-LINEAR (saturable sites
                 plus a non-specific term, Eq. 16). Armitage uses a linear
                 spLFER on the NEUTRAL log Kow and -- for acids -- applies no
                 ionisation penalty at all (httk SFbsa_acidic = 0).
  (2) CELL.      Qin *measures* a lumped D_cell/w. Armitage composes the cell
                 from storage lipid + membrane lipid + structural protein,
                 each driven by the ANION's Kow = 10^(logKow_n - 3.5).
  (3) IONISATION Qin uses measured anionic log Kow / log Dow. Armitage uses
                 generic constant offsets (SFkow 3.5, SFmw 1, SFbsa 0).
  (4) DOSE.      Qin's albumin sink shrinks as C_w rises, so the cell fraction
                 is dose-dependent. Armitage is linear everywhere: the cell
                 fraction is a dose-invariant constant.

Parameter provenance is tracked per compound. Qin measured only PFBA, PFOA,
PFHxS and PFOS; PFHxA / PFHpA / PFNA are flagged `extrapolated`.

Sources
-------
Qin et al. 2023, Environ Int, "Role of bioavailability and protein binding of
  four anionic PFAS in cell-based bioassays for QIVIVE" -- Tables 1, 2, Fig A1,
  Eqs. 12-17, 20.
Armitage et al. 2021, Toxics 9:315, "Update and Evaluation of a High-Throughput
  In Vitro Mass Balance Distribution Model: IV-MBM EQP v2.0" -- Eqs. 1-3, 5, 6.
  The paper puts the biopartitioning spLFERs in Supplementary Section S1, which
  is paywalled/unavailable here; the numerical parameterisation below is taken
  from EPA's own implementation of the same model, httk R/armitage.R
  (USEPA/CompTox-ExpoCast-httk), whose plastic spLFERs and acid pH-trapping
  term reproduce Armitage 2021 Eqs. 5, 6 and 3 exactly. See ARMITAGE_PROVENANCE.
Rosenmai et al. 2017, J Appl Toxicol -- Table 1 (uptake, % of moles added).
Amstutz et al. 2022, Toxicology 480:153312 -- Table 5 (intracellular fraction).
Fragki et al. 2023, suppl. Table 2 -- nmol/mg protein in HepaRG.

Run:  python3 tk_review/scripts/mbm_qin_vs_armitage.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DB = os.path.join(ROOT, "db")

# --------------------------------------------------------------------------
# Compounds.  nCF = number of fully fluorinated carbons (Qin's Fig A1 x-axis):
# PFBA = CF3-CF2-CF2-COOH -> 3.  All five are strong acids (pKa < 1), so in a
# pH 7.4 medium Fneutral ~ 0 and Fnegative ~ 1.
# --------------------------------------------------------------------------
COMPOUNDS = [
    # name    nC nCF    MW      logKow_opera logKow_cosmo  kow_src
    ("PFBA",   4,  3,  214.04,   2.60,        2.95,        "Qin Table 1"),
    ("PFHxA",  6,  5,  314.05,   None,        None,        "interpolated"),
    ("PFHpA",  7,  6,  364.06,   None,        None,        "interpolated"),
    ("PFOA",   8,  7,  414.07,   5.68,        5.73,        "Qin Table 1"),
    ("PFNA",   9,  8,  464.08,   None,        None,        "extrapolated"),
]
TARGETS = [c[0] for c in COMPOUNDS]

# log Kow is the single input the whole Armitage limb rides on, and Qin reports
# it for only PFBA and PFOA.  Fill the gaps with a per-CF2 regression through
# those two anchors, separately for each prediction method, and carry the flag.
def _kow_fill(col):
    anchors = [(nCF, c[col]) for (nm, nC, nCF, mw, ko, kc, src) in COMPOUNDS
               for c in [{4: ko, 5: kc}] if c[col] is not None]
    (x0, y0), (x1, y1) = anchors[0], anchors[-1]
    slope = (y1 - y0) / (x1 - x0)
    return lambda nCF: y0 + slope * (nCF - x0), slope

KOW_OPERA, SLOPE_OPERA = _kow_fill(4)   # OPERA average, comptox.epa.gov
KOW_COSMO, SLOPE_COSMO = _kow_fill(5)   # COSMOtherm

# --------------------------------------------------------------------------
# QIN 2023
# --------------------------------------------------------------------------
# Measured BSA binding constants, Qin Table 1.  PFHxS/PFOS are kept only as
# validation rows -- they are sulfonates and are NOT used for extrapolation.
QIN_BSA_MEASURED = {
    # name     Kd[uM]    Cmax[umol/L_BSA] Cmax[mol/mol] logKspec logDnonspec  nFr-1    logKFr
    "PFBA":  (18.3,      9.22e3,          0.450,        2.40,    1.63,        -0.391,  2.80),
    "PFOA":  (0.317,     2.54e4,          1.24,         4.60,    3.70,        -0.381,  4.32),
    "PFHxS": (7.44e-2,   1.88e4,          0.919,        5.10,    3.47,        -0.480,  4.32),
    "PFOS":  (6.24e-3,   8.01e3,          0.391,        5.81,    4.84,        -0.381,  4.91),
}
# Measured lumped partition ratios, Qin Table 2 (log10).
QIN_T2 = {   # name:   (logDlipw, logDcellw, logDplasmaw_exp, logDplasmaw_mbm)
    "PFBA":  (1.61, 2.56, 2.12, 2.36),   # logDlipw itself predicted from Fig A1
    "PFOA":  (3.52, 3.27, 4.60, 4.56),
    "PFHxS": (4.13, 3.23, 5.18, 5.06),
    "PFOS":  (4.89, 4.19, 5.96, 5.77),
}
# Qin Fig A1: log D_lip/w vs number of fluorinated carbons, R2 = 0.95.
QIN_FIGA1 = (0.476, 0.172)

# Qin et al. 2024 (Environ Sci Technol 58:5727, the baseline-toxicity paper by
# the same group) Table 1: MEASURED log D_lip/w, pH 7.4, from Droge et al. and
# Ebert et al.  These cover three compounds this file was interpolating
# (PFHxA, PFHpA, PFNA) and supply a measured PFBA value where Qin 2023 Table 2
# carried a predicted one.  They are NOT substituted into the published model
# by default -- "Qin as published" must stay as published -- but are available
# as a variant via qin_params(dlipw_2024=True).
#
# Their own regression over the wider Droge+Ebert+Endo set is SI eq S9,
# log D_lip/w = 0.562 N - 0.494 (R2 0.969) for carboxylates and eq S10,
# 0.535 N + 0.658 (R2 0.977) for sulfonates -- both steeper than the Fig A1
# regression above, which is fit to a narrower set.
QIN2024_DLIPW = {"PFBA": 1.00, "PFHxA": 2.32, "PFHpA": 2.91, "PFOA": 3.52,
                 "PFNA": 4.25, "PFUnA": 4.54, "PFHxS": 4.13, "PFOS": 4.89,
                 "HFPO-DA": 2.41, "6:2 FTSA": 3.87, "PFOSA": 4.94}
QIN2024_EQ_S9 = (0.562, -0.494)     # PFCA,  R2 0.969
QIN2024_EQ_S10 = (0.535, 0.658)     # PFSA,  R2 0.977

# --------------------------------------------------------------------------
# Qin et al. 2024, Environ Sci Technol 58(23):9954-9966, doi
# 10.1021/acs.est.3c10824 ("Species Difference? Bovine, Trout, and Human
# Plasma Protein Binding..."), PMC11171458.  This paper REFITTED the BSA
# isotherms of QIN_BSA_MEASURED above, pooling the 2023 single-fibre
# C18-SPME points with new BioSPME points, and extended the set from 4
# compounds to 11.  Its Table 1 is reprinted verbatim and uncited as SI
# Table S6 of the baseline-toxicity paper (10.1021/acs.est.3c09950).
#
# WHY THIS MATTERS FOR PFOS.  Qin 2023's PFOS row is the one parameter set in
# QIN_BSA_MEASURED that should not be used, and the authors say so in the
# refit: "5 mg/mL of BSA was used for the C18-SPME, resulting in a bound
# fraction of PFOS > 99% in the low concentration range.  The bound fraction
# reduced to 40-90% after adjusting BSA to 0.1 mg/mL in this study.  log
# D_BSA/w of PFOS were also fitted with data from two methods, but several
# points in the low concentrations were excluded."  At >99% bound the free
# concentration is below what the method resolves, so the SPECIFIC term --
# the one that lives at low concentration -- was fitted where there was
# least information.
#
# Two independent checks agree.  Qin 2023's PFOS Kd of 0.00624 uM is 111x
# tighter than Gao 2025's bovine value (0.693 uM) and 151x tighter than Chen
# 2025's human one (0.94 uM); and the plateau it implies, log(Cmax/Kd +
# D_nonspec) = 6.13, sits 1.2-2.0 log units above every independent
# sub-saturating measurement of BSA partitioning (Allendorf 2019, Bischel
# 2011, Fischer 2024: 4.10-4.89).  See qin2024_tableS6_provenance.py.
#
# The refit does NOT publish a new Kd/Cmax pair, so a corrected Kd is not
# recoverable from the tables; what it publishes is the specific and
# non-specific distribution ratios and the Freundlich parameters, below.
# Holding Cmax fixed would imply Kd = Cmax/(2*10^logKspec) = 0.0215 uM, still
# 32x tighter than Gao's, so even the refit does not reconcile with the
# independent affinities.  That assumption is not made anywhere here.
PFOS_KD_2023_IS_UNRELIABLE = True

# Table 2, BSA columns: (log D specific, log D non-specific).  None for
# specific means no plateau of K_spec appeared at nu < 1, so the paper could
# not resolve a specific term and reports an average instead.
# NOTE the basis: "specific" is Qin Eq. 15, the ratio AT HALF-SATURATION
# (Cmax/(2 Kd)), which is 0.301 log BELOW the Cw -> 0 plateau of the specific
# term.  Reading it as a plateau makes it 2x low.
QIN2024_BSA_SPEC_NONSPEC = {
    "PFBA":     (2.44, 1.97),
    "PFHxA":    (3.26, 2.69),
    "PFHpA":    (4.12, 3.47),
    "PFOA":     (4.58, 3.88),
    "PFNA":     (None, 4.61),
    "PFUnA":    (None, 4.75),
    "HFPO-DA":  (3.31, 2.17),
    "PFHxS":    (5.02, 3.58),
    "PFOS":     (5.27, 4.17),
    "6:2 FTSA": (None, 3.87),
    "PFOSA":    (None, 4.32),
}

# Table 1, BSA column: (nFr - 1, log K_Fr), i.e. slope and intercept of
# log D_BSA/w against log Cw with Cw in umol/L.  A slope of None is a value
# the paper reports as concentration-independent, n_Fr having been fixed to 1.
QIN2024_BSA_FREUNDLICH = {
    "PFBA":     (-0.298, 2.91),
    "PFHxA":    (-0.271, 3.43),
    "PFHpA":    (-0.305, 3.95),
    "PFOA":     (-0.314, 4.38),
    "PFNA":     (-0.147, 4.52),
    "PFUnA":    (None,   4.75),
    "HFPO-DA":  (-0.493, 3.44),
    "PFHxS":    (-0.472, 4.28),
    "PFOS":     (-0.379, 4.74),
    "6:2 FTSA": (-0.092, 3.86),
    "PFOSA":    (-0.105, 4.28),
}

# BSA volume basis.  Qin reports Cmax both per litre of BSA and per mole of
# BSA; the ratio pins the density of pure albumin and must agree with the
# 0.733 mL/g partial specific volume httk uses -- a shared-basis check that
# lets Qin's D_BSA/w and Armitage's K_bsa be compared at all.
BSA_MW = 66463.0          # g/mol
BSA_PARTIAL_SPEC_VOL = 0.733   # mL/g (httk armitage.R, Valb)

def _loglin(nCF, anchors):
    """Linear interpolation/extrapolation in nCF.  anchors = [(nCF, value), ...]"""
    (x0, y0), (x1, y1) = anchors
    return y0 + (y1 - y0) * (nCF - x0) / (x1 - x0)

# Carboxylate anchors for every Qin parameter that has no published rule.
_PFBA_nCF, _PFOA_nCF = 3, 7
def qin_params(name, nCF, dlipw_2024=False, bsa_2024=False):
    """Qin's parameters for one compound, with provenance.

    Measured -> straight from Table 1/2.
    Otherwise -> interpolated in nCF between the two CARBOXYLATE anchors
    (PFBA, PFOA).  Only log D_lip/w has a published rule (Fig A1); everything
    else is a two-point interpolation with no basis in the paper, so it is
    flagged.  PFNA (nCF 8) sits beyond both anchors = true extrapolation.
    """
    a, b = QIN_FIGA1
    d_lip = a * nCF + b                      # published rule, carboxylates only
    if name in QIN_BSA_MEASURED and name in QIN_T2:
        Kd, Cmax, CmaxM, lKs, lDns, nFr1, lKFr = QIN_BSA_MEASURED[name]
        d_cell = QIN_T2[name][1]
        # Use the MEASURED lipid-water ratio, not the Fig A1 regression.  The
        # regression is fit to perfluoroCARBOXYLATES, and for those the two
        # agree (PFBA 1.600 vs 1.61, PFOA 3.504 vs 3.52).  For the sulfonates
        # it is out by about a log unit (PFHxS 3.03 vs 4.13, PFOS 3.98 vs
        # 4.89), because sulfonates partition into lipid far more strongly at
        # the same chain length.
        d_lip = QIN_T2[name][0]
        src = "measured"
    else:
        pb = QIN_BSA_MEASURED["PFBA"]; po = QIN_BSA_MEASURED["PFOA"]
        Kd    = 10 ** _loglin(nCF, [(_PFBA_nCF, math.log10(pb[0])), (_PFOA_nCF, math.log10(po[0]))])
        Cmax  = 10 ** _loglin(nCF, [(_PFBA_nCF, math.log10(pb[1])), (_PFOA_nCF, math.log10(po[1]))])
        CmaxM = Cmax / (1e6 * BSA_PARTIAL_SPEC_VOL / BSA_MW * 1e-3 * 1e6)  # informational
        lKs   = _loglin(nCF, [(_PFBA_nCF, pb[3]), (_PFOA_nCF, po[3])])
        lDns  = _loglin(nCF, [(_PFBA_nCF, pb[4]), (_PFOA_nCF, po[4])])
        nFr1  = _loglin(nCF, [(_PFBA_nCF, pb[5]), (_PFOA_nCF, po[5])])
        lKFr  = _loglin(nCF, [(_PFBA_nCF, pb[6]), (_PFOA_nCF, po[6])])
        d_cell = _loglin(nCF, [(_PFBA_nCF, QIN_T2["PFBA"][1]), (_PFOA_nCF, QIN_T2["PFOA"][1])])
        src = "extrapolated" if not (_PFBA_nCF <= nCF <= _PFOA_nCF) else "interpolated"
    if dlipw_2024 and name in QIN2024_DLIPW:
        d_lip = QIN2024_DLIPW[name]
    if bsa_2024:
        # The 2024 refit, which supersedes Qin 2023 for the four overlapping
        # compounds and extends the set to eleven. Kd and Cmax are NOT
        # refitted by that paper, so they are left as they were and the
        # source is marked, rather than mixing a 2023 Kd into a 2024 fit
        # without saying so.
        if name in QIN2024_BSA_SPEC_NONSPEC:
            spec, nonspec = QIN2024_BSA_SPEC_NONSPEC[name]
            lDns = nonspec
            if spec is not None:
                lKs = spec
        if name in QIN2024_BSA_FREUNDLICH:
            sl, it = QIN2024_BSA_FREUNDLICH[name]
            nFr1, lKFr = (0.0 if sl is None else sl), it
        src += " + 2024 BSA refit (Kd/Cmax still 2023)"
    return dict(Kd=Kd, Cmax=Cmax, logKspec=lKs, logDnonspec=lDns,
                nFr_minus1=nFr1, logKFr=lKFr,
                logDlipw=d_lip, logDcellw=d_cell, source=src)

def qin_dbsaw(p, Cw, mode="mechanistic"):
    """Qin D_BSA/w at aqueous concentration Cw (uM).

    mechanistic (Eq. 16): saturable specific sites + linear non-specific term
        D = Cmax/(Kd + Cw) + D_non-specific
    freundlich (Eq. 12):  log D = (nFr-1) log Cw + log K_Fr
    linear:               the Henry-regime limit, K_specific (Eq. 15), i.e.
                          what Armitage's linear treatment is equivalent to.
    """
    if mode == "mechanistic":
        return p["Cmax"] / (p["Kd"] + Cw) + 10 ** p["logDnonspec"]
    if mode == "freundlich":
        return 10 ** (p["nFr_minus1"] * math.log10(max(Cw, 1e-12)) + p["logKFr"])
    if mode == "linear":
        return 10 ** p["logKspec"]
    raise ValueError(mode)

# --------------------------------------------------------------------------
# ARMITAGE 2021 / httk IV-MBM EQP
# --------------------------------------------------------------------------
ARMITAGE_PROVENANCE = (
    "Numerical parameterisation from EPA httk R/armitage.R "
    "(USEPA/CompTox-ExpoCast-httk, armitage_eval); its plastic spLFERs "
    "(0.97*logKow-6.94 Kramer; 0.56*logKow-4.635 Fischer) and acid ion-trapping "
    "term 10^(cell_pH - pH) reproduce Armitage 2021 Eqs. 5, 6 and 3."
)
ARM = dict(
    # cell composition (volume fractions) and octanol-proportionality constants
    pseudooct=0.01, memblip=0.04, nlom=0.05, P_cells=1.0, P_nlom=0.035,
    # ionisation scaling factors: log PC of the ion = log PC of the neutral - SF
    SFkow=3.5, SFmw=1.0, SFbsa_acidic=0.0, SFplw=3.5,
    # anionic phospholipid handling
    Anionic_VF=0.175, A_Prop_acid=0.05,
    # lysosomes (inert for acids: Lyso_Pump = 1)
    Lyso_VF=0.0068, Lyso_Diam=500.0, Lyso_MemThick=7.0,
    # medium / serum
    conc_ser_alb=24.0, conc_ser_lip=1.9,   # g/L in serum
    csalt=0.15, cell_pH=7.4, pH=7.4,
    cellmass=3.04, celldensity=1.0,        # ng/cell, kg/L -> 3.04 pL/cell
)

# Qin's Table 2 header reads "log Dcell/w [Lw/Lcell,prot+lip]": the ratio is per
# litre of cell PROTEIN+LIPID, not per litre of cell.  Armitage's DR_kcw is per
# litre of whole cell (its composition sums to 1.0 and includes 0.90 water), so
# Qin's cell term must be scaled by that fraction before the two share an axis.
#
# Qin measured it in the same HEK293H cells D_cell/w came from: "The volume
# fractions of water, protein and lipid of cells were 91.6%, 8.14% and 0.432%
# (Table A5)".  Those two numbers also fix the protein-per-cell-volume factor
# that converts a measured nmol/mg protein into a concentration, so they are
# kept as ONE pair -- setting them independently would assume two different
# cells in the same calculation.  The ratio predicted/measured is nearly
# invariant to the pair, because both sides scale with the protein fraction;
# see composition_invariance().
# The measured datasets are HepG2 and HepaRG, so the default is the project
# spreadsheet's HepG2 basis: 675 pg protein in a 3.04 pL cell = 222 mg/mL =
# a 0.164 protein volume fraction at 1.35 g/mL.  Qin's own HEK293H value is
# 0.0814 -- 2x lower -- and is carried as the alternative.  Transferring Qin's
# D_cell/w to HepG2 means applying it to HepG2's prot+lip volume, which is the
# point of normalising to prot+lip in the first place.
PROTEIN_DENSITY_G_ML = 1.35
CELL_VF_PROTEIN = 0.164        # HepG2, project spreadsheet
CELL_VF_LIPID = 0.00870        # scaled from Qin's protein:lipid volume ratio
CELL_VF_PROTEIN_QIN = 0.0814   # Qin Table A5, HEK293H
CELL_VF_LIPID_QIN = 0.00432    # Qin Table A5

F_PROTLIP = CELL_VF_PROTEIN + CELL_VF_LIPID

def armitage_pcs(logKow_n, kbsa2=False, kpl2=False):
    """spLFERs, httk armitage.R lines 663-694."""
    gkmw_n = 1.01 * logKow_n + 0.12                      # Endo 2011, membrane lipid
    if kbsa2:                                            # Endo & Goss 2011 Fig 2
        gkbsa_n = (1.08 * logKow_n - 0.7) if logKow_n < 4.5 else (0.37 * logKow_n + 2.56)
    else:
        gkbsa_n = 0.71 * logKow_n + 0.42                 # Endo & Goss 2011 (httk default)
    gkpl_n = (0.56 * logKow_n - 4.635) if kpl2 else (0.97 * logKow_n - 6.94)
    ksalt = 0.04 * logKow_n + 0.114                      # Ni & Yalkowsky 2003
    return gkmw_n, gkbsa_n, gkpl_n, ksalt

def armitage_drs(logKow_n, kbsa2=False, kpl2=False, Fneutral=0.0):
    """Distribution ratios for an acidic IOC.  httk armitage.R lines 736-790."""
    A = ARM
    Fn, Fneg = Fneutral, 1.0 - Fneutral
    gkmw_n, gkbsa_n, gkpl_n, ksalt = armitage_pcs(logKow_n, kbsa2, kpl2)
    kow_n, kow_i = 10 ** logKow_n, 10 ** (logKow_n - A["SFkow"])
    kmw_n, kmw_i = 10 ** gkmw_n, 10 ** (gkmw_n - A["SFmw"])
    kbsa_n = 10 ** gkbsa_n
    kbsa_i_acidic = 10 ** (gkbsa_n - A["SFbsa_acidic"])   # SF = 0 -> anion keeps it all
    kpl_n, kpl_i = 10 ** gkpl_n, 10 ** (gkpl_n - A["SFplw"])
    salt = 10 ** (-ksalt * A["csalt"])                    # divide by this

    # --- cell ---
    cell_DR_kow = (Fn * kow_n + Fneg * kow_i) / salt
    avf = A["Anionic_VF"]
    cell_DR_kmw = ((1 - avf) * (Fn * kmw_n + Fneg * kmw_i)
                   + avf * (Fn * kmw_n + Fneg * kmw_i * A["A_Prop_acid"])) / salt
    cellwat = 1 - (A["pseudooct"] + A["memblip"] + A["nlom"])
    IOC_mult = 10 ** (A["cell_pH"] - A["pH"])             # acid ion trapping, Eq. 3
    pre = (A["P_cells"] * A["pseudooct"] * cell_DR_kow
           + A["memblip"] * cell_DR_kmw
           + A["P_nlom"] * A["nlom"] * cell_DR_kow
           + cellwat) * IOC_mult
    r = A["Lyso_Diam"] / 2.0
    lyso_memvf = (r ** 3 - (r - A["Lyso_MemThick"]) ** 3) / r ** 3
    DR_kcw = (1 - A["Lyso_VF"]) * pre + A["Lyso_VF"] * (lyso_memvf * cell_DR_kmw * IOC_mult)

    # --- medium ---
    DR_kow = (Fn * kow_n + Fneg * kow_i) / salt
    DR_kbsa = (Fn * kbsa_n + Fneg * kbsa_i_acidic) / salt
    DR_kpl = (Fn * kpl_n + Fneg * kpl_i) / salt
    return dict(DR_kcw=DR_kcw, DR_kbsa=DR_kbsa, DR_kow=DR_kow, DR_kpl=DR_kpl,
                gkmw_n=gkmw_n, gkbsa_n=gkbsa_n, gkpl_n=gkpl_n,
                cell_DR_kmw=cell_DR_kmw, cell_DR_kow=cell_DR_kow)

# --------------------------------------------------------------------------
# Assay geometry.  Volumes in L, areas in m2.
# --------------------------------------------------------------------------
def well_sarea(diam_mm, v_working_uL, bottom=True):
    """Wetted plastic area: base + side wall up to the medium meniscus."""
    r = diam_mm / 20.0                       # cm
    base = math.pi * r ** 2                  # cm2
    h = (v_working_uL / 1000.0) / base       # cm
    side = 2 * math.pi * r * h
    return (base * bottom + side) * 1e-4     # cm2 -> m2

ASSAYS = {
    # Rosenmai 2017: HepG2, 12-well, 18e4 cells/well, DMEM + 10% FBS, 24 h.
    # Working volume is not stated in the paper; 1.0 mL is the standard 12-well
    # volume and is flagged as an assumption.
    # "lysates from three wells per treatment were pooled" -- so the unit of
    # measurement is 3 wells: 543,400 cells at end of exposure and 3x the
    # medium.  Uptake % is invariant to the pooling, but the protein mass and
    # the absolute nmol are not.
    "Rosenmai2017": dict(cells="HepG2", plate="12-well x3 pooled", diam_mm=22.1,
                         v_working_uL=3000.0, cell_yield=5.434e5, FBSf=0.10,
                         doses=[10.0, 30.0, 100.0], vol_assumed=True, wells=3,
                         basis="nominal"),
    # Amstutz 2022: HepG2, 96-well, 4e4 cells/well, 250 uL exposure medium,
    # DMEM + 10% FBS, 100 uM, 24 h.  IF is normalised to RECOVERED mass
    # (cells + medium at 24 h), not to nominal -- mass balance was 93-106%.
    "Amstutz2022": dict(cells="HepG2", plate="96-well", diam_mm=6.4,
                        v_working_uL=250.0, cell_yield=4.0e4, FBSf=0.10,
                        doses=[100.0], vol_assumed=False,
                        basis="recovered"),
    # Fragki 2022 (Food Chem Toxicol 172:113559).  HepaRG, 24-well, dosed in
    # 550 uL of "assay medium (growth medium containing 2% FBS)" -- the growth
    # medium's 10% is not what the cells were exposed in.  See
    # mbm_medium_composition.py.  Volume kept at the older 500 uL here because
    # this script is superseded by mbm_published_equations.py.
    "Fragki2023": dict(cells="HepaRG", plate="24-well(assumed)", diam_mm=15.6,
                       v_working_uL=500.0, cell_yield=2.5e5, FBSf=0.02,
                       doses=[6.25, 25.0, 100.0], vol_assumed=True,
                       basis="nominal"),
}

# Protein -> cell-volume factor.  Needed only to put "per mg protein" data on a
# concentration basis.  A cell of density 1 g/mL that is ~18% protein by mass
# holds ~180 mg protein per mL of cell volume; equivalently 0.54 ng protein per
# httk's 3 ng cell.  This is a COMMON-MODE factor: it shifts Rosenmai and
# Fragki together and does not touch Amstutz, so it cannot create scatter
# between studies -- only a uniform offset.
MG_PROTEIN_PER_ML_CELL = CELL_VF_PROTEIN * PROTEIN_DENSITY_G_ML * 1000.0

def volumes(a):
    """Compartment volumes in L.  httk armitage.R lines 820-830."""
    Vbm = a["v_working_uL"] * 1e-6
    Vcells = a["cell_yield"] * (ARM["cellmass"] * 1e-6) / ARM["celldensity"] * 1e-6
    Valb = Vbm * a["FBSf"] * BSA_PARTIAL_SPEC_VOL * ARM["conc_ser_alb"] * 1e-3
    Vslip = Vbm * a["FBSf"] * ARM["conc_ser_lip"] * 1e-3
    Vm = Vbm - Valb - Vslip
    n = a.get("wells", 1)
    sarea = n * well_sarea(a["diam_mm"], a["v_working_uL"] / n)
    return dict(Vbm=Vbm, Vcells=Vcells, Valb=Valb, Vslip=Vslip, Vm=Vm, sarea=sarea)

# --------------------------------------------------------------------------
# Mass balance.  Both models share this algebra; they differ only in the
# partition ratios fed into it.  Armitage 2021 Eq. 2 == Qin Eq. 20 in form.
#     n_tot = C_w * ( Vm + D_bsa*Valb + D_ow*Vslip + D_cw*Vcells + D_pl*sarea*1000 )
# Headspace is dropped: PFAS anions are non-volatile (K_aw ~ 0).
# --------------------------------------------------------------------------
def balance(nomconc_uM, V, D_cell, D_bsa, D_ow=0.0, D_pl=0.0, cell_frac=1.0):
    """cell_frac = fraction of cell volume D_cell is expressed against.
    Armitage: 1.0 (whole cell).  Qin: F_PROTLIP (cell protein+lipid only).
    C_cell is always returned on the WHOLE-CELL basis so the two are
    comparable with each other and with a measured per-cell concentration."""
    Vc_eff = V["Vcells"] * cell_frac
    ntot = nomconc_uM * V["Vbm"]                       # umol
    denom = (V["Vm"] + D_bsa * V["Valb"] + D_ow * V["Vslip"]
             + D_cell * Vc_eff + D_pl * V["sarea"] * 1000.0)
    Cw = ntot / denom
    C_cell = D_cell * Cw * cell_frac                   # whole-cell basis
    n_cell = C_cell * V["Vcells"]
    n_med = Cw * V["Vm"] + D_bsa * Cw * V["Valb"] + D_ow * Cw * V["Vslip"]
    n_pl = D_pl * Cw * V["sarea"] * 1000.0
    return dict(Cw=Cw, C_cell=C_cell,
                x_cell_nominal=n_cell / ntot,
                x_cell_recovered=n_cell / (n_cell + n_med),
                x_alb=D_bsa * Cw * V["Valb"] / ntot,
                x_plastic=n_pl / ntot, x_water=Cw * V["Vm"] / ntot)

def qin_predict(name, nCF, nomconc_uM, V, mode="mechanistic"):
    """Qin's model.  Non-linear albumin binding needs a self-consistent C_w."""
    p = qin_params(name, nCF)
    D_cell = 10 ** p["logDcellw"]                      # [Lw / L_cell,prot+lip]
    Cw = nomconc_uM * 1e-3                             # seed
    for _ in range(400):
        D_bsa = qin_dbsaw(p, Cw, mode)
        new = balance(nomconc_uM, V, D_cell, D_bsa, cell_frac=F_PROTLIP)["Cw"]
        if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-10:
            Cw = new
            break
        Cw = math.sqrt(Cw * new)                       # damped, in log space
    out = balance(nomconc_uM, V, D_cell, qin_dbsaw(p, Cw, mode), cell_frac=F_PROTLIP)
    out.update(D_cell=D_cell * F_PROTLIP, D_cell_protlip=D_cell,
               D_bsa=qin_dbsaw(p, Cw, mode), params=p)
    return out

def armitage_predict(logKow_n, nomconc_uM, V, kbsa2=False, kpl2=False):
    d = armitage_drs(logKow_n, kbsa2, kpl2)
    out = balance(nomconc_uM, V, d["DR_kcw"], d["DR_kbsa"], d["DR_kow"], d["DR_kpl"])
    out.update(D_cell=d["DR_kcw"], D_bsa=d["DR_kbsa"], drs=d)
    return out

# --------------------------------------------------------------------------
# Measured data
# --------------------------------------------------------------------------
# Rosenmai 2017 Table 1: cellular uptake as % of total moles added.
ROSENMAI_UPTAKE_PCT = {
    "PFBA":  {10: 0.14, 30: 0.13, 100: 0.10},
    "PFHxA": {10: 0.19, 30: 0.19, 100: 0.14},
    "PFHpA": {10: 0.16, 30: 0.22, 100: 0.21},
    "PFOA":  {10: 0.24, 30: 0.27, 100: 0.33},
    "PFNA":  {10: 0.15, 30: 0.22, 100: 0.31},
}
# Amstutz 2022 Table 5: intracellular fraction (%) at 100 uM, 24 h, mean+-SD.
# Normalised to RECOVERED mass.  Mass balance 93-106%.
AMSTUTZ_IF_PCT = {
    "PFBA":  (None, None),      # not detected in any sample
    "PFHxA": (0.039, 0.015),
    "PFHpA": (0.044, 0.016),
    "PFOA":  (0.245, 0.151),
    "PFNA":  (0.302, 0.213),
}
# Fragki 2023 suppl. Table 2: nmol/mg protein in HepaRG, 24 h.
FRAGKI_NMOL_PER_MG = {
    "PFOA": {6.25: 0.09226, 25.0: 0.52306, 100.0: 2.79309},
    "PFNA": {6.25: 0.22203, 100.0: 15.69467},
}
# Rosenmai's own absolute scale is internally contradictory: the Abstract says
# "4-2500 ng mg-1 protein", Results 3.3 says "18-2490 nmol mg-1 protein".
ROSENMAI_UNIT_CONFLICT = dict(abstract="ng/mg protein, 4-2500",
                              results="nmol/mg protein, 18-2490")

def ccell_from_fraction(x_cell, nomconc_uM, V):
    """Fraction of nominal dose in cells -> C_cell (uM)."""
    return x_cell * nomconc_uM * V["Vbm"] / V["Vcells"]

def ccell_from_nmol_per_mg(v, mg_per_mL=None):
    """nmol/mg protein -> uM.  nmol/mg * mg/mL = nmol/mL = umol/L.

    mg_per_mL is resolved at CALL time, not bound as a default, so the
    sensitivity runs below actually move the measured side.
    """
    if mg_per_mL is None:
        mg_per_mL = MG_PROTEIN_PER_ML_CELL
    return v * mg_per_mL

# --------------------------------------------------------------------------
# Validation against the papers' own published numbers
# --------------------------------------------------------------------------
def composition_invariance():
    """The cell-composition pair cancels out of predicted/measured for Qin.

    Measured C_cell scales with mg protein per mL of cell, i.e. with the
    protein volume fraction.  Qin's predicted C_cell scales with prot+lip.
    Both move together, so the ratio that the comparison turns on barely
    shifts.  Armitage's cell is on a whole-cell basis and does not scale, so
    its ratio moves with the full factor -- stated here rather than hidden.
    """
    global CELL_VF_PROTEIN, F_PROTLIP, MG_PROTEIN_PER_ML_CELL
    keep = (CELL_VF_PROTEIN, F_PROTLIP, MG_PROTEIN_PER_ML_CELL)
    out = []
    for vf, vl in ((CELL_VF_PROTEIN, CELL_VF_LIPID),
                   (CELL_VF_PROTEIN_QIN, CELL_VF_LIPID_QIN)):
        CELL_VF_PROTEIN = vf
        F_PROTLIP = vf + vl
        MG_PROTEIN_PER_ML_CELL = vf * PROTEIN_DENSITY_G_ML * 1000.0
        V = volumes(ASSAYS["Rosenmai2017"])
        meas = ccell_from_nmol_per_mg(126.0 / 214.04)      # Rosenmai PFBA 100 uM
        q = qin_predict("PFBA", 3, 100.0, V, "mechanistic")["C_cell"]
        a = armitage_predict(2.60, 100.0, V)["C_cell"]
        out.append((vf, F_PROTLIP, MG_PROTEIN_PER_ML_CELL, q / meas, a / meas))
    CELL_VF_PROTEIN, F_PROTLIP, MG_PROTEIN_PER_ML_CELL = keep
    return out


def validate():
    checks = []
    def ck(label, got, want, tol):
        ok = abs(got - want) <= tol * max(abs(want), 1e-12)
        checks.append((label, got, want, ok))
        return ok

    # Qin Table 1 is internally consistent two ways.
    for nm, (Kd, Cmax, CmaxM, lKs, lDns, _n, _k) in QIN_BSA_MEASURED.items():
        # (a) K_specific = Cmax/(2*Kd)   (Eq. 15)
        ck(f"Qin Eq15 logKspec {nm}", math.log10(Cmax / (2 * Kd)), lKs, 0.01)
        # (b) Cmax[umol/L_BSA] / Cmax[mol/mol] pins the density of pure albumin,
        #     which must equal 1/0.733 g/mL -- the basis httk's Valb uses.  This
        #     is what makes Qin's D_BSA/w and Armitage's K_bsa comparable.
        dens = (Cmax / CmaxM) * 1e-6 * BSA_MW          # g per L of BSA
        ck(f"BSA density from {nm}", dens, 1000.0 / BSA_PARTIAL_SPEC_VOL, 0.01)

    # The 2024 refit moved every compound measured twice by at most 0.34 log
    # EXCEPT PFOS, which moved 0.54 on the specific term and 0.67 on the
    # non-specific one. That isolation is the point: it is not a general
    # recalibration, it is one compound's experiment being corrected.
    moved = {}
    for nm, (Kd, Cmax, CmaxM, lKs, lDns, _n, _k) in QIN_BSA_MEASURED.items():
        spec, nonspec = QIN2024_BSA_SPEC_NONSPEC[nm]
        moved[nm] = max(abs(spec - lKs) if spec is not None else 0.0,
                        abs(nonspec - lDns))
    others = max(v for k, v in moved.items() if k != "PFOS")
    ck("2024 refit moves PFOS more than any other compound",
       1.0 if moved["PFOS"] > others else 0.0, 1.0, 0.01)
    ck("PFOS moves >0.5 log, every other compound <0.4",
       1.0 if moved["PFOS"] > 0.5 and others < 0.4 else 0.0, 1.0, 0.01)

    # Qin Fig A1 reproduces the measured carboxylate log D_lip/w ...
    a, b = QIN_FIGA1
    ck("FigA1 -> PFBA logDlipw", a * 3 + b, QIN_T2["PFBA"][0], 0.02)
    ck("FigA1 -> PFOA logDlipw", a * 7 + b, QIN_T2["PFOA"][0], 0.02)
    # ... but NOT the sulfonates: PFHxS (nCF 6) is 1.1 decades off, which is why
    # only PFBA and PFOA are used as extrapolation anchors.
    checks.append(("FigA1 fails on PFHxS (expected)",
                   a * 6 + b, QIN_T2["PFHxS"][0],
                   abs(a * 6 + b - QIN_T2["PFHxS"][0]) > 0.5))

    # Armitage: httk's plastic spLFERs must reproduce Armitage 2021 Eqs. 5/6.
    for lk in (2.0, 5.0, 7.0):
        ck(f"Kramer Eq5 logKpl @{lk}", armitage_pcs(lk)[2], 0.97 * lk - 6.94, 1e-9)
        ck(f"Fischer Eq6 logKpl @{lk}", armitage_pcs(lk, kpl2=True)[2],
           0.56 * lk - 4.64, 0.01)

    # Fragki's unit chain: nmol/mg * 126 mg protein/g liver = umol/kg liver.
    ck("Fragki PFOA 100uM umol/kg", 2.79309 * 126.0, 351.93, 0.01)
    ck("Fragki PFOA 100uM ug/L", 2.79309 * 126.0 * 414.07, 145699.0, 0.01)

    bad = [c for c in checks if not c[3]]
    print("VALIDATION  %d/%d checks pass" % (len(checks) - len(bad), len(checks)))
    for lab, got, want, ok in checks:
        if not ok:
            print("  FAIL %-34s got %.6g want %.6g" % (lab, got, want))
    assert not bad, "validation failed"
    return checks

# --------------------------------------------------------------------------
def main():
    validate()
    print()
    print("Cell composition (HepG2 basis): protein %.3f, lipid %.5f by volume"
          " -> F_protlip %.4f, %.0f mg protein/mL cell"
          % (CELL_VF_PROTEIN, CELL_VF_LIPID, F_PROTLIP, MG_PROTEIN_PER_ML_CELL))
    print("   sensitivity to that pair (Qin Table A5 gives 0.0814 for HEK293H):")
    for vf, fp, mg, qr, ar in composition_invariance():
        print("   vf_protein %.4f -> Qin/measured %6.2fx   Armitage/measured %6.2fx"
              % (vf, qr, ar))
    print("Armitage parameterisation:", ARMITAGE_PROVENANCE)
    print("log Kow infill: OPERA %+.3f per CF2, COSMOtherm %+.3f per CF2 "
          "(2 anchors: PFBA, PFOA)" % (SLOPE_OPERA, SLOPE_COSMO))
    print()

    rows = []
    for (name, nC, nCF, MW, ko, kc, kow_src) in COMPOUNDS:
        lkow_op = ko if ko is not None else KOW_OPERA(nCF)
        lkow_cs = kc if kc is not None else KOW_COSMO(nCF)
        qp = qin_params(name, nCF)
        for study, a in ASSAYS.items():
            V = volumes(a)
            for dose in a["doses"]:
                # measured
                meas_x, meas_basis, meas_note = None, a["basis"], ""
                if study == "Rosenmai2017":
                    v = ROSENMAI_UPTAKE_PCT.get(name, {}).get(int(dose))
                    if v is not None:
                        meas_x = v / 100.0
                elif study == "Amstutz2022":
                    v, sd = AMSTUTZ_IF_PCT.get(name, (None, None))
                    if v is not None:
                        meas_x = v / 100.0
                    else:
                        meas_note = "not detected"
                elif study == "Fragki2023":
                    v = FRAGKI_NMOL_PER_MG.get(name, {}).get(dose)
                    if v is not None:
                        meas_cc = ccell_from_nmol_per_mg(v)
                        meas_x = meas_cc * V["Vcells"] / (dose * V["Vbm"])
                        meas_basis = "nmol/mg protein"
                if meas_x is None and not meas_note:
                    continue
                meas_cc = (ccell_from_fraction(meas_x, dose, V)
                           if meas_x is not None else None)

                # predictions
                q_nl = qin_predict(name, nCF, dose, V, "mechanistic")
                q_fr = qin_predict(name, nCF, dose, V, "freundlich")
                q_li = qin_predict(name, nCF, dose, V, "linear")
                a_op = armitage_predict(lkow_op, dose, V)
                a_cs = armitage_predict(lkow_cs, dose, V)
                a_alt = armitage_predict(lkow_op, dose, V, kbsa2=True)

                key = "x_cell_recovered" if a["basis"] == "recovered" else "x_cell_nominal"
                rows.append(dict(
                    compound=name, nC=nC, nCF=nCF, MW=MW, study=study,
                    cells=a["cells"], plate=a["plate"], dose_uM=dose,
                    logKow_opera=round(lkow_op, 3), logKow_cosmo=round(lkow_cs, 3),
                    logKow_source=kow_src,
                    qin_param_source=qp["source"],
                    qin_logDcellw=round(qp["logDcellw"], 3),
                    qin_logDlipw=round(qp["logDlipw"], 3),
                    arm_logKbsa=round(a_op["drs"]["gkbsa_n"], 3),
                    arm_logDcellw=round(math.log10(a_op["D_cell"]), 3),
                    meas_basis=meas_basis, meas_note=meas_note,
                    meas_x_cell_pct=(round(meas_x * 100, 4) if meas_x is not None else ""),
                    meas_C_cell_uM=(round(meas_cc, 1) if meas_cc is not None else ""),
                    qin_nonlin_C_cell_uM=round(q_nl["C_cell"], 1),
                    qin_freund_C_cell_uM=round(q_fr["C_cell"], 1),
                    qin_linear_C_cell_uM=round(q_li["C_cell"], 1),
                    arm_opera_C_cell_uM=round(a_op["C_cell"], 1),
                    arm_cosmo_C_cell_uM=round(a_cs["C_cell"], 1),
                    arm_kbsa2_C_cell_uM=round(a_alt["C_cell"], 1),
                    qin_x_cell_pct=round(q_nl[key] * 100, 4),
                    arm_x_cell_pct=round(a_op[key] * 100, 4),
                    qin_Cw_uM=round(q_nl["Cw"], 4),
                    arm_Cw_uM=round(a_op["Cw"], 4),
                    qin_x_alb_pct=round(q_nl["x_alb"] * 100, 2),
                    arm_x_alb_pct=round(a_op["x_alb"] * 100, 2),
                    arm_x_plastic_pct=round(a_op["x_plastic"] * 100, 3),
                    qin_over_meas=(round(q_nl["C_cell"] / meas_cc, 2) if meas_cc else ""),
                    arm_over_meas=(round(a_op["C_cell"] / meas_cc, 2) if meas_cc else ""),
                    qin_over_arm=round(q_nl["C_cell"] / a_op["C_cell"], 2),
                ))

    os.makedirs(DB, exist_ok=True)
    out = os.path.join(DB, "mbm_qin_vs_armitage.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote %s  (%d rows x %d cols)" % (out, len(rows), len(rows[0])))
    return rows

if __name__ == "__main__":
    ROWS = main()
