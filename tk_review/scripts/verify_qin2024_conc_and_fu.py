#!/usr/bin/env python3
"""Where Qin 2024's concentration columns come from, and whether the
paper carries a fraction unbound.

THE PAPER.  Qin W, Henneberger L, Gluge J, Konig M, Escher BI.
"Baseline Toxicity Model to Identify the Specific and Nonspecific
Effects of Per- and Polyfluoroalkyl Substances in Cell-Based
Bioassays."  Environ Sci Technol 2024;58(13):5727-5738.
doi 10.1021/acs.est.3c09950.  This is the BASELINE paper, the active
anchor in qin2024_anchor.py -- not the plasma paper, whose Table 1
carries a different column of the same name.

THE THREE COLUMNS, AND WHERE EACH IS READ.

  PFAS conc (uM)      Table 1, column IC10,nom [mol/L], x 1e6.
                      The NOMINAL concentration causing 10%
                      cytotoxicity in the PPARgamma-GeneBLAzer assay.
                      It is a dosed concentration, not a measured
                      free one, which is what makes it the right
                      numerator for a loading.

  Protein conc (uM)   SI Table S7, PPARgamma row: the assay medium is
                      98% OptiMEM + 2% csFBS and its protein volume
                      fraction is 0.94 mL/L.  Main text (Cell Binding
                      Experiments) states the density of protein as
                      1.36 kg/L.  So

                        0.94e-3 L/L x 1.36 kg/L / 66463 g/mol
                                               = 19.23 umol/L

                      This is TOTAL MEDIUM PROTEIN expressed as an
                      albumin-equivalent molarity.  It is not a BSA
                      solution: the baseline paper ran no BSA
                      isotherm.  Its log D_BSA/w column is the plasma
                      paper's Freundlich fit (SI Table S6) read at
                      IC10,free, per Table 1 footnote b.  The one
                      number is the same for all eleven rows because
                      it is a property of the medium, not of the
                      compound.

  ligand:protein      IC10,nom / 19.23.  A quotient of ours, and the
                      only one of the three that is not read off a
                      page.

  (log K itself is Table 1's log D_BSA/w minus log10(1.36) = 0.134,
   the volume-to-mass basis change, which is why no K_d or n appears
   on these rows: nothing was converted.)

THE FRACTION UNBOUND.  The paper prints no such column, but it prints
both concentrations it is made of.  SI eq S3 defines

    C_free,medium = C_nom x f_free,medium x V_tot / V_w,medium

and V_tot/V_w = 1.0014 for this assay, so to within 0.15%

    f_free,medium = IC10,free / IC10,nom

which runs 1.8% to 55.8% free, i.e. 44% to 98% bound.  IC10,free was
MEASURED, by solid-phase microextraction on the free
concentration-cytotoxicity curve, so this is an experimental free
fraction and not a model output.

WHAT IT IS A FREE FRACTION OF.  Not albumin.  SI eq S4 defines
f_free,medium over the whole bioassay system, and eq 6 names its four
sorptive phases: medium protein, medium lipid, cell structural
protein, cell lipid.  Part B below decomposes the denominator and
finds medium protein carries 92.8% to 99.2% of all bound PFAS, median
97.6%, so for THIS assay the distinction costs at most 1.4 percentage
points -- the albumin-only free fraction and the system free fraction
are the same number to within the data's resolution.  That is a
property of a 2% serum medium with 0.94 mL/L protein and 0.0147 mL/L
lipid, and it would not survive a 10% FBS medium.

Part C closes the loop the way this project closes every other one:
two numbers from one paper that should reproduce each other.  Eq 3
with Table 1's own composite D_medium/w and D_cell/w returns the
measured f_free to a median of 0.11 log.  The deviations are
systematically NEGATIVE (9 of 11), so the measurement finds less free
PFAS than protein and lipid account for -- consistent with sorption to
the plate, which eq 6 does not carry a term for.

Run:  python3 -I tk_review/scripts/verify_qin2024_conc_and_fu.py
"""
import math

# --- Table 1, main text p5731, transcribed as printed ----------------
# (IC10_nom M, IC10_free M, logD_med/w, logD_cell/w,
#  logD_BSA/w, logD_SP/w, logD_lip/w)
T1 = {
    "PFBA":     (3.19e-3, 1.78e-3, 2.60, 2.40, 1.94, 2.11, 1.00),
    "PFHxA":    (9.95e-4, 4.72e-4, 2.78, 3.45, 2.70, 2.45, 2.32),
    "PFHpA":    (3.75e-4, 6.76e-5, 3.21, 3.69, 3.40, 2.65, 2.91),
    "PFOA":     (1.30e-4, 1.84e-5, 3.67, 3.94, 3.98, 3.14, 3.52),
    "PFNA":     (9.07e-5, 7.72e-6, 4.14, 3.74, 4.39, 3.39, 4.25),
    "PFUnA":    (2.16e-5, 3.83e-7, 4.74, 3.78, 4.75, 3.97, 4.54),
    "HFPO-DA":  (6.54e-4, 2.20e-4, 3.18, 3.09, 2.28, 2.77, 2.41),
    "PFHxS":    (1.46e-4, 4.22e-5, 3.14, 3.56, 3.52, 2.95, 4.13),
    "PFOS":     (4.78e-5, 3.89e-6, 4.04, 4.32, 4.52, 3.94, 4.89),
    "6:2 FTSA": (3.07e-4, 4.25e-5, 3.45, 3.90, 3.71, 3.39, 3.87),
    "PFOSA":    (7.85e-6, 3.28e-7, 4.18, 3.81, 4.33, 3.52, 4.94),
}

# --- SI Table S7, PPARgamma row --------------------------------------
VF_PROT_MED, VF_LIP_MED, VF_W_MED = 0.94e-3, 1.47e-5, 0.999
V_MED = 4.00e-5                      # L, 384-well medium volume
# --- SI Table S8, PPARgamma row (HEK293H) ----------------------------
VF_PROT_CELL, VF_LIP_CELL = 33.9e-3, 5.3e-3
V_CELL = 1.71e-8                     # L

RHO_PROT = 1.36                      # kg/L, main text
M_ALB = 66463.0                      # g/mol, BSA

Vw = VF_W_MED * V_MED
Vp_med, Vl_med = VF_PROT_MED * V_MED, VF_LIP_MED * V_MED
Vp_cell, Vl_cell = VF_PROT_CELL * V_CELL, VF_LIP_CELL * V_CELL
Vpl_med, Vpl_cell = Vp_med + Vl_med, Vp_cell + Vl_cell
Vtot = V_MED + V_CELL
PROT_uM = VF_PROT_MED * RHO_PROT / M_ALB * 1e9


def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2]


def rule(t):
    print("\n" + "=" * 70 + f"\n{t}\n" + "=" * 70)


def main():
    rule("A.  The three columns, and the page each is read from")
    print(f"  Protein conc, every row: {VF_PROT_MED*1e3:.2f} mL/L "
          f"(SI Table S7) x {RHO_PROT} kg/L (main text)")
    print(f"                           / {M_ALB:.0f} g/mol "
          f"= {PROT_uM:.2f} uM\n")
    print(f"  {'PFAS':9s} {'IC10,nom':>11s} {'= PFAS conc':>12s} "
          f"{'ratio':>8s}   {'logD_BSA/w':>10s} {'- 0.134':>8s}")
    print(f"  {'':9s} {'[mol/L]':>11s} {'[uM]':>12s} "
          f"{'[mol/mol]':>8s}   {'Table 1':>10s} {'= log K':>8s}")
    print("  " + "-" * 64)
    for c, v in T1.items():
        nom, lbsa = v[0], v[4]
        print(f"  {c:9s} {nom:11.3g} {nom*1e6:12.3g} "
              f"{nom*1e6/PROT_uM:8.2f}   {lbsa:10.2f} "
              f"{lbsa - math.log10(RHO_PROT):8.3f}")

    rule("B.  Is there a fraction unbound?  Yes -- and it is albumin's")
    print("  f_free,medium = IC10,free / IC10,nom   (SI eq S3; the")
    print(f"  V_tot/V_w factor is {Vtot/Vw:.4f}, so it is a 0.15% "
          "correction)\n")
    print(f"  {'PFAS':9s} {'f_free':>8s} {'f_bound':>8s} | "
          f"{'f_free, albumin':>15s} {'cost of':>8s} | {'albumin share':>13s}")
    print(f"  {'':9s} {'measured':>8s} {'':>8s} | "
          f"{'term alone':>15s} {'the swap':>8s} | {'of binding':>13s}")
    print("  " + "-" * 70)
    shares, costs, fb = [], [], {}
    for c, (nom, free, lmed, lcell, lbsa, lsp, llip) in T1.items():
        f_obs = free / nom
        t_alb, t_lipm = 10**lbsa * Vp_med, 10**llip * Vl_med
        t_psc, t_lipc = 10**lsp * Vp_cell, 10**llip * Vl_cell
        bound = t_alb + t_lipm + t_psc + t_lipc
        f_alb = Vtot / (Vw + t_alb)
        f_sys = Vtot / (Vw + bound)
        share, cost = t_alb / bound * 100, (f_alb - f_sys) * 100
        shares.append(share)
        costs.append(cost)
        fb[c] = (1 - f_obs) * 100
        print(f"  {c:9s} {f_obs*100:7.2f}% {(1-f_obs)*100:7.2f}% | "
              f"{f_alb*100:14.2f}% {cost:+7.2f}pp | {share:12.1f}%")
    print(f"\n  medium protein carries {min(shares):.1f}% to "
          f"{max(shares):.1f}% of all bound PFAS, median "
          f"{med(shares):.1f}%")
    print(f"  swapping the system f_free for the albumin-only one costs "
          f"at most {max(costs):.2f} pp")
    print("  -- so for THIS assay they are the same number.  That rests")
    print("     on 2% serum: 0.94 mL/L protein against 0.0147 mL/L lipid.")

    rule("C.  Do the paper's own two numbers reproduce each other?")
    print("  eq 3 : V_tot / (V_w + D_med/w.V_pl,med + D_cell/w.V_pl,cell)")
    print("  eq 6 : the same, split into its four named phases\n")
    print(f"  {'PFAS':9s} {'f_obs':>8s} | {'eq 3':>8s} {'dev':>7s} | "
          f"{'eq 6':>8s} {'dev':>7s}")
    print("  " + "-" * 56)
    d3, d6 = [], []
    for c, (nom, free, lmed, lcell, lbsa, lsp, llip) in T1.items():
        f_obs = free / nom
        f3 = Vtot / (Vw + 10**lmed * Vpl_med + 10**lcell * Vpl_cell)
        f6 = Vtot / (Vw + 10**lbsa * Vp_med + 10**llip * Vl_med
                     + 10**lsp * Vp_cell + 10**llip * Vl_cell)
        e3, e6 = math.log10(f_obs / f3), math.log10(f_obs / f6)
        d3.append(e3)
        d6.append(e6)
        print(f"  {c:9s} {f_obs*100:7.2f}% | {f3*100:7.2f}% {e3:+7.2f} | "
              f"{f6*100:7.2f}% {e6:+7.2f}")
    print(f"\n  eq 3  |dev| median {med([abs(x) for x in d3]):.2f}  "
          f"max {max(abs(x) for x in d3):.2f} log   "
          f"signed median {med(d3):+.2f}")
    print(f"  eq 6  |dev| median {med([abs(x) for x in d6]):.2f}  "
          f"max {max(abs(x) for x in d6):.2f} log   "
          f"signed median {med(d6):+.2f}")
    neg = sum(1 for x in d3 if x < 0)
    print(f"\n  {neg} of {len(d3)} eq-3 deviations are negative: the "
          "measurement finds")
    print("  LESS free PFAS than protein and lipid account for. eq 6 "
          "carries no")
    print("  term for sorption to the plate, which is the obvious "
          "candidate.")
    print("\n  Note this is a consistency check, not a validation: "
          "Table 1's")
    print("  D_BSA/w is itself the Table S6 fit READ AT IC10,free "
          "(footnote b),")
    print("  so the measured free concentration is an input to the "
          "prediction.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
