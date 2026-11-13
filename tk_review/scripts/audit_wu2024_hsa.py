#!/usr/bin/env python3
"""Why Wu et al. 2024 sits four logs below every other HSA study.

Wu Y, Bao J, Liu Y, Wang X, Lu X, Wang K.  In Vitro and In Silico
Analysis of the Bindings between Legacy and Novel Per- and
Polyfluoroalkyl Substances and Human Serum Albumin.  Toxics
2024;12(1):46.  doi 10.3390/toxics12010046, PMID 38251003.  Open
access.

The paper reports binding constants KA of 20 to 395 M^-1 for ten PFAS
on HSA.  Every other albumin study in this database lands between 10^3
and 10^6.  Converted onto the log K axis, Wu's ten compounds come out
at -0.84 to +0.71 against a literature median near 3.7: a uniform gap
of -2.7 to -3.8 logs with no compound-specific structure.  A uniform
gap is the signature of a method or an analysis, not of chemistry.

EVERY NUMBER BELOW IS ARITHMETIC ON THE PAPER'S OWN TABLES 1 TO 6.
Nothing is refitted and nothing is assumed.  The four checks are
independent and all four fail.

  1  f_a, the accessible fraction of fluorophore, comes out above 1
     for all ten compounds -- up to 4.07.  HSA has ONE tryptophan, as
     the paper itself states, so f_a cannot exceed 1.  The modified
     Stern-Volmer model does not describe these data, and the Kq that
     the static-quenching argument rests on is built on that model.

  2  In lg[(F0-F)/F] = lg KA + n lg C, KA carries units of M^-n, not
     M^-1.  Inverting it to a Kd in molar is only valid at n = 1, and
     Wu's n runs 0.47 to 0.87.  Re-expressing each fit as the 1:1
     constant it implies INSIDE the measured window reproduces the
     paper's own independent Stern-Volmer constant to within a factor
     of about 2 on 8 of 10 compounds, which is what says the
     re-expression is the right one.

  3  Re-expressed that way, all ten compounds collapse to 1.21 to 1.35
     on the log K axis -- a 0.14-log spread where the paper reports
     2.2 logs.  The paper's headline ordering (PFSAs above PFCAs,
     affinity rising with chain length) does not survive.

  4  lg KA regresses on n with R^2 = 0.992 and a slope of 3.23, within
     reach of the 3.50 that pure extrapolation geometry predicts.  The
     measured quenching at the middle concentration is 0.34 to 0.46
     for every compound -- a 1.35-fold spread -- while the reported KA
     spans 19-fold.  lg KA is a restatement of the fitted exponent,
     not a measurement of affinity.

THE ASSAY ITSELF.  Section 2.1: each PFAS was prepared as a 1 mM stock
in ANHYDROUS ETHANOL.  Section 2.2: 5 mL of 1 mM HSA was mixed with 1,
3 and 5 mL of that stock.  So the three tubes are 17, 38 and 50 %
ethanol v/v, against the under 1 % that albumin binding work normally
tolerates.

  READ THE RANGE AS A FLOOR, NOT AS THE EXTENT.  Section 2.2 is the
  three-tube screening step, run to pick a ratio for the 3-D spectra,
  and it yields no binding constant; the MIX below is its geometry,
  kept because the concentrations are the only ones the paper states
  numerically.  The constants in Tables 4-6 come from Section 2.3,
  whose series runs 5:0, 5:1, 5:5, 5:10, 5:15, 5:20, 5:30 -- so the
  top point of every binding curve is 86 % ethanol v/v and HSA falls
  5.8-fold, 833 to 143 uM.  Everything below therefore understates
  the cosolvent and the protein drift.

Two further consequences of that protocol:

  - the HSA concentration is NOT held constant.  It falls from 833 to
    500 uM across the series while the PFAS rises, so protein and
    ligand are varied together.  This is a dilution series of a fixed
    total, not a titration.  (A pure dilution artefact would give
    n = 1.47; the fitted n is 0.47 to 0.87, so dilution alone does not
    account for it -- but F0 is never defined against a
    concentration-matched blank, and the text does not say whether one
    was used.)

  - three concentrations, one replicate, no error bars.  A two-
    parameter line through three points reaches R^2 = 0.99 for free.

WHAT FOLLOWS.  The rows stay in db/protein_binding.csv: they are a real
paper, correctly transcribed, and the database is a record of what the
literature says.  They are excluded from the albumin figures and from
any pooled statistic, by name and with this file as the reason.  See
EXCLUDED in albumin_chainlength_data.py.

Run:  python3 tk_review/scripts/audit_wu2024_hsa.py
"""
import math

MIX = [(5, 1), (5, 5), (5, 10), (5, 15), (5, 20), (5, 30)]
# Section 2.3, the two-dimensional titration that produces the
# constants -- NOT section 2.2's 5:1/5:3/5:5, which is the
# screening step the paper runs to pick a ratio for its 3-D
# spectra and from which no binding constant is derived.
C = [p / (h + p) * 1e-3 for h, p in MIX]            # PFAS, mol/L
P = [h / (h + p) * 1e-3 for h, p in MIX]            # HSA,  mol/L
ETOH = [p / (h + p) for h, p in MIX]                # ethanol, v/v
CSTAR = math.exp(sum(math.log(c) for c in C) / len(C))   # geometric mean

M_ALB = 66463.0
BRIDGE = math.log10(1e9 / M_ALB)                    # 4.1774
TAU0 = 1e-8                                         # s, as the paper takes it

# Tables 1-3 (modified Stern-Volmer, 1/C on the x axis) and Tables 4-6
# (Lineweaver-Burk double logarithmic), transcribed as printed.
#  name          SV slope   SV int    Kq printed   n      lg KA    KA
T = [
    ("PFOA",         0.0005035, 0.6350, 1.26e11, 0.6811, 2.0143, 103.3475),
    ("PFOS",         0.0007969, 0.2456, 3.08e10, 0.8694, 2.5962, 394.6390),
    ("PFHxS",        0.0005641, 0.4741, 8.40e10, 0.7328, 2.1985, 157.9429),
    ("PFPrS",        0.0004813, 0.7054, 1.47e11, 0.6487, 1.9001,  79.4511),
    ("PFBS",         0.0005327, 0.6013, 1.13e11, 0.6986, 2.0663, 116.4930),
    ("PFBA",         0.0003302, 0.9979, 3.02e11, 0.4692, 1.3089,  20.3657),
    ("9Cl-PF3ONS",   0.0006538, 0.7081, 1.08e11, 0.7240, 2.0686, 117.1120),
    ("11Cl-PF3OUdS", 0.0005285, 0.6068, 1.15e11, 0.7102, 2.1059, 127.6145),
    ("HFPO-DA",      0.0005095, 0.8051, 1.58e11, 0.6214, 1.7708,  58.9929),
    ("NaDONA",       0.0005118, 0.6743, 1.32e11, 0.6718, 1.9689,  93.0894),
]


def assay():
    print("THE ASSAY, as mixed (Section 2.3, the titration)")
    for (h, p), c, pr, e in zip(MIX, C, P, ETOH):
        print(f"  {h} mL HSA : {p} mL PFAS   PFAS {c*1e6:6.1f} uM"
              f"   HSA {pr*1e6:6.1f} uM   L:P {c/pr:4.2f}"
              f"   ethanol {e*100:4.1f} % v/v")
    print(f"  geometric-mean PFAS across the {len(MIX)} points: "
          f"{CSTAR*1e6:.1f} uM")
    print(f"  HSA is not held constant: it falls "
          f"{P[0]/P[-1]:.2f}-fold across the series.\n")


def check1():
    print("CHECK 1  f_a = 1 / intercept of the modified Stern-Volmer plot.")
    print("         HSA has one tryptophan, so f_a cannot exceed 1.0.")
    print(f"  {'compound':14}{'f_a':>7}{'Kq re-derived':>16}"
          f"{'Kq printed':>13}  verdict")
    for name, slope, inter, kq, n, lgka, ka in T:
        fa = 1.0 / inter
        ksv = 1.0 / (slope * fa)
        print(f"  {name:14}{fa:7.2f}{ksv/TAU0:16.2e}{kq:13.2e}"
              f"  {'ok' if fa <= 1.0 else 'IMPOSSIBLE'}")
    print("  (Kq re-derived matches Kq printed, so the transcription is")
    print("   right and the slope/intercept really are what produced it.)\n")


def check2():
    print("CHECK 2  KA is in M^-n, so 1/KA is a Kd only at n = 1.")
    print("         Re-express inside the measured window:")
    print("             K_eff = KA * C*^(n-1)   [M^-1]")
    print("         and test against the paper's own Stern-Volmer K'_SV.")
    print(f"  {'compound':14}{'n':>6}{'KA (M^-n)':>11}{'K_eff (M^-1)':>14}"
          + "K'_SV".rjust(10) + "ratio".rjust(7))
    out = {}
    for name, slope, inter, kq, n, lgka, ka in T:
        ksv = 1.0 / (slope * (1.0 / inter))
        keff = ka * CSTAR ** (n - 1.0)
        out[name] = (n, ka, keff)
        print(f"  {name:14}{n:6.3f}{ka:11.1f}{keff:14.0f}{ksv:10.0f}"
              f"{keff/ksv:7.2f}")
    return out


def check3(eff):
    print("\nCHECK 3  the two routes on the log K axis of the figures.")
    print("         as ingested : log K = 4.1774 + log n - log Kd(uM),"
          "  Kd = 1/KA")
    print("         re-expressed: log K = 4.1774 - log Kd(uM),       "
          "  Kd = 1/K_eff")
    print(f"  {'compound':14}{'as ingested':>13}{'re-expressed':>14}"
          f"{'shift':>8}")
    fixed = []
    for name, (n, ka, keff) in eff.items():
        a = BRIDGE + math.log10(n) - math.log10(1e6 / ka)
        b = BRIDGE - math.log10(1e6 / keff)
        fixed.append(b)
        print(f"  {name:14}{a:13.2f}{b:14.2f}{b-a:+8.2f}")
    print(f"  re-expressed spread: {min(fixed):.2f} to {max(fixed):.2f}"
          f"  =  {max(fixed)-min(fixed):.2f} log, for ten different"
          f" compounds.\n")


def check4():
    print("CHECK 4  is lg KA anything but a restatement of n?")
    print("         If the quenching actually measured at C* is similar")
    print("         for every compound then lg KA = y* - n lg C*, so")
    print(f"         lg KA must be linear in n with slope"
          f" {-math.log10(CSTAR):.2f}.")
    ns = [r[4] for r in T]
    lgs = [r[5] for r in T]
    mn, ml = sum(ns) / len(ns), sum(lgs) / len(lgs)
    sxy = sum((a - mn) * (b - ml) for a, b in zip(ns, lgs))
    sxx = sum((a - mn) ** 2 for a in ns)
    syy = sum((b - ml) ** 2 for b in lgs)
    print(f"         observed: slope {sxy/sxx:.2f}, "
          f"R^2 {(sxy**2)/(sxx*syy):.4f}")
    print(f"\n  {'compound':14}{'n':>7}{'lg KA':>8}{'lg[(F0-F)/F] at C*':>21}"
          f"{'(F0-F)/F':>11}")
    q = []
    for name, slope, inter, kq, n, lgka, ka in T:
        y = lgka + n * math.log10(CSTAR)
        q.append(10 ** y)
        print(f"  {name:14}{n:7.3f}{lgka:8.3f}{y:21.3f}{10**y:11.3f}")
    ka_all = [r[6] for r in T]
    print(f"\n  measured quenching at C* spans {max(q)/min(q):.2f}-fold.")
    print(f"  reported KA spans              "
          f"{max(ka_all)/min(ka_all):.2f}-fold.")


def main():
    assay()
    check1()
    eff = check2()
    check3(eff)
    check4()


if __name__ == "__main__":
    main()
