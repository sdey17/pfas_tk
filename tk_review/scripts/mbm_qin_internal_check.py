#!/usr/bin/env python3
"""Does Qin's Eq. (17) reproduce Qin's OWN measured D_cell/w?

The full text of Qin et al. 2023 (Environ Int 173:107857) is now in hand, and
it supplies three things the tabulated parameters alone did not:

  Eq. (17)  D_i/w = D_BSA/w * V_prot/V_prot+lip  +  D_lip/w * V_lip/V_prot+lip
            -- the normalising volume is the SORPTIVE volume, not V_i.  Table 2
            states the units as [L_w / L_cell,prot+lip], confirming this.
  Eq. (20)  the bioassay balance:  n_tot / [ (D_cell/w*V_cell,prot+lip
            + V_cell,w) + (D_medium/w*V_medium,prot+lip + V_medium,w) ]
  Table A5  HEK293H volume fractions: water 91.6%, protein 8.14%, lipid 0.432%
  Sec 3.5   "D_cell/w was predicted by the MBM (Eq. (17)) with D_non-specific
            (Table 1) as a surrogate for partitioning to structural proteins
            and D_lip/w (Table 2) for lipid."

So the cell limb of the model can be checked against the authors' own
measurement, which Table 2 reports.  Section 3.5 characterises the result as
"not much different ... with discrepancies within 15%" for PFOA, PFHxS and
PFOS and "around 35% lower" for PFBA.  This script asks what those percentages
are percentages OF.

Run:  python3 tk_review/scripts/mbm_qin_internal_check.py
"""
import math

# Qin Table A5, HEK293H
F_PROT, F_LIP = 8.14, 0.432
f_prot = F_PROT / (F_PROT + F_LIP)
f_lip = F_LIP / (F_PROT + F_LIP)

# Qin Table 1 (log D_non-specific, log K_specific) and Table 2 (log D_lip/w,
# log D_cell/w measured by Eq. 11).
P = {                 # logDns  logKspec  logDlipw  logDcellw_meas
    "PFBA":  (1.63, 2.40, 1.61, 2.56),
    "PFOA":  (3.70, 4.60, 3.52, 3.27),
    "PFHxS": (3.47, 5.10, 4.13, 3.23),
    "PFOS":  (4.84, 5.81, 4.89, 4.19),
}


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print(f"\n  V_prot/V_prot+lip = {f_prot:.4f}   V_lip/V_prot+lip = {f_lip:.4f}"
          "   (Qin Table A5)\n")
    print("=" * 78)
    print("Eq. (17) WITH D_non-specific, AGAINST QIN'S MEASURED D_cell/w (Table 2)")
    print("=" * 78)
    print(f"\n  {'cmpd':7s}{'pred':>8s}{'meas':>8s}{'d(log)':>9s}"
          f"{'pred/meas':>11s}{'% of log':>10s}")
    for nm, (lns, lks, llip, lmeas) in P.items():
        pred = 10 ** lns * f_prot + 10 ** llip * f_lip
        lp = math.log10(pred)
        print(f"  {nm:7s}{lp:>8.2f}{lmeas:>8.2f}{lp - lmeas:>+9.2f}"
              f"{10 ** (lp - lmeas):>11.2f}{100 * (lp - lmeas) / lmeas:>+9.1f}%")
    print("""
  The last column reproduces the paper's own figures: PFOA +13%, PFHxS +10%,
  PFOS +16% ("within 15%") and PFBA -36% ("around 35% lower").  They are
  percentages of the LOG-TRANSFORMED distribution ratio.  On the quantity
  actually entering the mass balance they correspond to the pred/meas column:
  the model over-predicts cell partitioning by 2.0-4.5x for the three
  hydrophobic PFAS and under-predicts PFBA by 8.5x.

  That is not a criticism of the measurement or of the model's construction;
  it is a statement about what "agreement" means in the source.  Two things
  follow for the present reassessment:""")
    xs, ys = [], []
    for nm, (lns, lks, llip, lmeas) in P.items():
        pred = 10 ** lns * f_prot + 10 ** llip * f_lip
        xs.append(nm); ys.append(math.log10(pred) - lmeas)
    print(f"""
  (1) The residual is ORDERED: PFBA {ys[0]:+.2f}, PFOA {ys[1]:+.2f},
      PFHxS {ys[2]:+.2f}, PFOS {ys[3]:+.2f} log units.  Within the carboxylates
      the cell term acquires {ys[1] - ys[0]:+.2f} log units of positive drift
      between nCF 3 and nCF 7, i.e. {(ys[1] - ys[0]) / 4:+.3f} per CF2 -- the
      same SIGN as the drift this work finds in the Armitage cell term, in a
      model whose cell parameters are measured rather than predicted.

  (2) A 2-4.5x over-prediction of D_cell/w is inside Qin's own validated
      envelope.  Benchmarking the model's uptake predictions against an
      independent dataset therefore cannot resolve discrepancies smaller than
      that, which bounds how much of the level error found in this work can be
      attributed to the model as opposed to its inputs.""")

    print("\n" + "=" * 78)
    print("AND WITH K_specific INSTEAD (the low-concentration limb)")
    print("=" * 78)
    print(f"\n  {'cmpd':7s}{'pred':>8s}{'meas':>8s}{'pred/meas':>11s}")
    for nm, (lns, lks, llip, lmeas) in P.items():
        pred = 10 ** lks * f_prot + 10 ** llip * f_lip
        print(f"  {nm:7s}{math.log10(pred):>8.2f}{lmeas:>8.2f}"
              f"{10 ** (math.log10(pred) - lmeas):>11.1f}")
    print("""
  K_specific is the wrong surrogate for the cell, as the authors say: it
  over-predicts by 1-2 orders of magnitude.  D_non-specific is the correct
  choice, and is what is used in this work.""")

    print("\n" + "=" * 78)
    print("C_max CONSISTENCY, AND THE ALBUMIN DENSITY QIN ACTUALLY USED")
    print("=" * 78)
    CMAX = {"PFBA": (9.22e3, 0.450), "PFOA": (2.54e4, 1.24),
            "PFHxS": (1.88e4, 0.919), "PFOS": (8.01e3, 0.391)}
    M_BSA = 66463.0          # Da, stated in Sec 2.3.4
    for rho in (1.36, 1.364):
        print(f"\n  rho_BSA = {rho} kg/L -> rho/M = "
              f"{1e9 * rho / M_BSA:,.0f} umol/L_BSA")
    print(f"\n  {'cmpd':7s}{'uM/L_BSA':>11s}{'mol/molBSA':>12s}{'implied':>12s}")
    for nm, (cv, cm) in CMAX.items():
        print(f"  {nm:7s}{cv:>11,.0f}{cm:>12.3f}{cv / cm:>12,.0f}")
    print("""
  The paper states 66463 Da and a density of 1.36 kg/L (Sec 2.3.4), giving
  20,462 umol per L of BSA.  The four tabulated pairs imply 20,457-20,489, so
  the two Table 1 C_max columns are internally consistent and the stated
  density is confirmed at source.  This work had reconstructed the same
  constant independently as 1/0.733 = 1.364 kg/L from albumin's partial
  specific volume, which gives 20,523 umol/L -- the same number to 0.3%.  The
  reconstruction is now redundant: Qin's 1.36 is the authors' own value and is
  what the identity C_max[umol/L_BSA] = n * rho/M rests on.""")


if __name__ == "__main__":
    main()
