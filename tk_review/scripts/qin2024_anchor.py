#!/usr/bin/env python3
"""The single definition of the "Qin 2024" albumin anchor.

WHY THIS FILE EXISTS.  Three downstream scripts need the Qin albumin
coefficient, the free concentration it belongs to, and the loading it
describes: albumin_shareable.py (the sheet), albumin_conc_vs_logk.py (the
concentration columns) and make_route_and_loading_figures.py (the figures).
Each used to carry its own hard-coded copy, from a different Qin paper and a
different concentration.  One switch now sets all three.

THE TWO PAPERS.  Weiping Qin published twice in Environ Sci Technol in 2024,
and BOTH main-text Table 1s carry a column headed log D_BSA/w holding
DIFFERENT numbers.  verify_qin2024_bsa_provenance.py settles which is which.

  BASELINE  "Baseline Toxicity Model to Identify the Specific and
            Nonspecific Effects of Per- and Polyfluoroalkyl Substances in
            Cell-Based Bioassays."  Qin W, Henneberger L, Gluge J, Konig M,
            Escher BI.  Environ Sci Technol 2024;58(13):5727-5738.
            doi 10.1021/acs.est.3c09950.
            Table 1 log D_BSA/w = the Freundlich fit READ AT IC10,free, the
            free concentration causing 10% cytotoxicity in its PPARgamma
            assay.  PFOA 3.98.

  PLASMA    "Species Difference?  Bovine, Trout, and Human Plasma Protein
            Binding of Per- and Polyfluoroalkyl Substances."  Qin W, Escher
            BI, Huchthausen J, Fu Q, Henneberger L.  Environ Sci Technol
            2024;58(23):9954-9966.  doi 10.1021/acs.est.3c10824.
            Table 1 log D_BSA/w = the FIT ITSELF, slope and intercept, from
            its own BSA isotherms (SI Table S5).  PFOA 4.38 at Cw = 1 uM.

They are not in conflict: the baseline paper's column is the plasma paper's
fit evaluated at a particular concentration, and reproduces it to within
0.008 log for all eleven compounds (`python3 -I qin2024_anchor.py` checks
this).  They differ because D_BSA/w is not a constant.

ACTIVE SETS WHICH ONE IS THE ANCHOR.  Changing it moves the sheet, the
concentration columns and both figures together.

WHAT THE CHOICE COSTS.  The two anchors sit at loadings three orders apart,
so this is not a cosmetic switch:

    nu = D [L_w/kg] x Cw [mol/L_w] x M_albumin [kg/mol]      mol PFAS / mol

  baseline   nu 0.34 to 11.5   -- at or past albumin's site count, because
             IC10,free is a cytotoxic dose.  The numbers are printed in a
             table and need no calculation of ours.
  plasma     nu 0.10 to 0.16   -- a tenth-occupied carrier, nearer (though
             still well above) the physiological nu < 0.01.  The numbers
             need our own mass balance to place them at a Cw.

Neither is the trace-limit constant this review argues a TK model wants.
The baseline anchor is the more conservative choice for "what does a
measured albumin coefficient look like", and the more honest one about
provenance, since nothing is recomputed on our side.

Run:  python3 -I tk_review/scripts/qin2024_anchor.py
      prints the anchor, the alternative, and the discrepancy between them.
"""
import math

ACTIVE = "baseline"

DOI = {"baseline": "10.1021/acs.est.3c09950",
       "plasma": "10.1021/acs.est.3c10824"}
LABEL = {"baseline": "Qin 2024", "plasma": "Qin 2024 (plasma)"}

M_ALB_KG = 66463.0 / 1e3          # kg/mol, bovine serum albumin
RHO_PROT = 1.36                   # kg/L, the density Qin states (base, p3)
VOL_TO_MASS = -math.log10(RHO_PROT)      # -0.134, L_prot basis -> kg basis

# Qin's PPARgamma assay medium, 98% OptiMEM + 2% csFBS: SI Table S7 gives
# 0.94 mL protein per L.  Expressed as albumin-equivalent molarity this is
# the protein the IC10 was measured against -- NOT the BSA solution the
# isotherm was measured in, which the baseline paper never ran.
MEDIUM_PROT_uM = 0.94e-3 * RHO_PROT / 66463.0 * 1e9      # 19.24 uM

# --- BASELINE paper, Table 1 -------------------------------------------
# (IC10,nom [mol/L], IC10,free [mol/L], log D_BSA/w [L_w/L_prot])
BASELINE_T1 = {
    "PFBA":     (3.19e-3, 1.78e-3, 1.94),
    "PFHxA":    (9.95e-4, 4.72e-4, 2.70),
    "PFHpA":    (3.75e-4, 6.76e-5, 3.40),
    "PFOA":     (1.30e-4, 1.84e-5, 3.98),
    "PFNA":     (9.07e-5, 7.72e-6, 4.39),
    "PFUnA":    (2.16e-5, 3.83e-7, 4.75),
    "HFPO-DA":  (6.54e-4, 2.20e-4, 2.28),
    "PFHxS":    (1.46e-4, 4.22e-5, 3.52),
    "PFOS":     (4.78e-5, 3.89e-6, 4.52),
    "6:2 FTSA": (3.07e-4, 4.25e-5, 3.71),
    "PFOSA":    (7.85e-6, 3.28e-7, 4.33),
}

# --- PLASMA paper, Table 1: (slope or None, intercept, R2) -------------
# Reprinted verbatim as the baseline paper's SI Table S6 BSA column.
PLASMA_T1 = {
    "PFBA": (-0.298, 2.91, 0.60), "PFHxA": (-0.271, 3.43, 0.88),
    "PFHpA": (-0.305, 3.95, 0.81), "PFOA": (-0.314, 4.38, 0.93),
    "PFNA": (-0.147, 4.52, 0.61), "PFUnA": (None, 4.75, 0.96),
    "HFPO-DA": (-0.493, 3.44, 0.86), "PFHxS": (-0.472, 4.28, 0.92),
    "PFOS": (-0.379, 4.74, 0.92), "6:2 FTSA": (-0.092, 3.86, 0.54),
    "PFOSA": (-0.105, 4.28, 0.50),
}

# Free Cw [uM] reached from 10 uM nominal in a generic 10% FBS medium,
# from Qin's eq 6 with the Table S7 / S8 compositions.  Ours, not his.
PLASMA_OP_CW = {
    "PFBA": 3.79, "PFHxA": 1.134, "PFHpA": 0.2344, "PFOA": 0.05577,
    "PFNA": 0.06658, "PFUnA": 0.05805, "HFPO-DA": 1.141, "PFHxS": 0.03635,
    "PFOS": 0.01085, "6:2 FTSA": 0.3965, "PFOSA": 0.1278,
}

# Plasma paper SI Table S5, BSA columns: (Ctot [uM], BSA [mg/L])
PLASMA_S5 = {
    "PFBA": (2712, 5000), "PFHxA": (731, 1500), "PFHpA": (78, 500),
    "PFOA": (100, 500), "PFNA": (100, 500), "PFUnA": (76, 500),
    "HFPO-DA": (203, 1500), "PFHxS": (54, 100), "PFOS": (55, 100),
    "6:2 FTSA": (106, 500), "PFOSA": (96, 500),
}

COMPOUNDS = list(BASELINE_T1)


def fit_at(name, cw_uM):
    """The plasma paper's Freundlich fit read at a free concentration,
    on the L_water/L_protein basis the paper states."""
    slope, b, _ = PLASMA_T1[name]
    return b if slope is None else slope * math.log10(cw_uM) + b


def nu(logk_mass, cw_uM):
    """Loading the coefficient describes, mol PFAS per mol albumin.

    Needs no protein concentration: it is a property of the coefficient
    and the free concentration alone.  nu = K x Cw x M."""
    return 10 ** logk_mass * cw_uM * M_ALB_KG / 1e6


def anchor(which=None):
    """One dict per compound for the active (or named) paper."""
    which = which or ACTIVE
    out = []
    for name in COMPOUNDS:
        if which == "baseline":
            nom, free, logd = BASELINE_T1[name]
            cw, lig = free * 1e6, nom * 1e6
            prot = MEDIUM_PROT_uM
            basis = "Table 1 as printed: the Freundlich fit read at IC10,free"
            where = (f"IC10,free = {cw:.3g} uM, read off the fit; dosed "
                     f"IC10,nom = {lig:.4g} uM over {prot:.3g} uM medium "
                     f"protein (PPARgamma assay)")
        else:
            cw = PLASMA_OP_CW[name]
            logd = fit_at(name, cw)
            ctot, mgL = PLASMA_S5[name]
            lig, prot = ctot, mgL / 66463.0 * 1e6
            basis = ("Table 1 fit read at the Cw our mass balance gives "
                     "for 10 uM nominal in 10% FBS; the dose and ratio are "
                     "SI Table S5's own BSA isotherm")
            where = f"free Cw = {cw:.3g} uM (10 uM nominal, 10% FBS)"
        logk = logd + VOL_TO_MASS
        out.append({"pfas": name, "logK": round(logk, 3), "cw_uM": cw,
                    "logD_vol": logd, "nu": nu(logk, cw),
                    "ligand_uM": lig, "protein_uM": prot,
                    "ratio": lig / prot, "where": where, "basis": basis,
                    "doi": DOI[which], "label": LABEL[which]})
    return out


def discrepancy():
    """Baseline anchor minus plasma anchor, per compound."""
    b = {r["pfas"]: r for r in anchor("baseline")}
    p = {r["pfas"]: r for r in anchor("plasma")}
    return [{"pfas": c, "baseline": b[c], "plasma": p[c],
             "dlogK": b[c]["logK"] - p[c]["logK"],
             "slope": PLASMA_T1[c][0]} for c in COMPOUNDS]


def main():
    print(f"ACTIVE = {ACTIVE!r}  ->  {LABEL[ACTIVE]}, doi {DOI[ACTIVE]}")
    print(f"medium protein (Table S7, 0.94 mL/L x {RHO_PROT} kg/L) = "
          f"{MEDIUM_PROT_uM:.2f} uM albumin-equivalent")
    print(f"volume -> mass basis: {VOL_TO_MASS:+.3f} log\n")

    print("CHECK  baseline Table 1 == plasma Table 1 fit read at IC10,free?")
    worst = 0.0
    for c in COMPOUNDS:
        _, free, printed = BASELINE_T1[c]
        got = fit_at(c, free * 1e6)
        worst = max(worst, abs(got - printed))
    print(f"       worst disagreement over 11 compounds: {worst:.3f} log"
          f"  -> {'same fit' if worst < 0.01 else 'DIFFERENT'}\n")

    print("%-9s %7s %7s %6s   %9s %9s %7s   %6s" % (
        "PFAS", "base", "plasma", "d log", "nu base", "nu plasma",
        "x more", "slope"))
    for d in sorted(discrepancy(), key=lambda d: d["dlogK"]):
        s = "flat" if d["slope"] is None else f"{d['slope']:.3f}"
        print("%-9s %7.3f %7.3f %+6.2f   %9.2f %9.3f %7.0f   %6s" % (
            d["pfas"], d["baseline"]["logK"], d["plasma"]["logK"],
            d["dlogK"], d["baseline"]["nu"], d["plasma"]["nu"],
            d["baseline"]["nu"] / d["plasma"]["nu"], s))
    ds = sorted(d["dlogK"] for d in discrepancy())
    print(f"\nmedian d log K = {ds[len(ds) // 2]:+.2f}, "
          f"range {ds[0]:+.2f} to {ds[-1]:+.2f}")
    print("The gap IS the Freundlich slope: a flat isotherm (PFOSA -0.105,")
    print("6:2 FTSA -0.092, PFUnA flat) barely moves, a steep one")
    print("(PFHxS -0.472, HFPO-DA -0.493) moves more than a log.")


if __name__ == "__main__":
    main()
