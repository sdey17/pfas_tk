"""
Rosenmai 2017: deciding whether Fig 2A is nmol/mg or ng/mg, using only
Table 1 and arithmetic.

THE IDEA
--------
Table 1 and Fig 2A are two normalisations of ONE measurement:

  Table 1  uptake %   = moles in cells / moles added      (dimensionless)
  Fig 2A   C_cell     = moles in cells / mg cell protein   (unit in question)

Divide one by the other and the "moles in cells" cancels, leaving a pure
property of the experiment:

      C_cell                 moles in cells / m          1
  ------------------  =  ------------------------  =  -------
   uptake x C_nom          uptake x C_nom              m / V

  because  moles added = C_nom x V.  Rearranged:

      V / m  =  C_cell / (uptake_frac x C_nom)          [L per mg protein]

V/m is "millilitres of medium per milligram of cell protein". It is the same
number for every compound and every dose, and it is bounded by the plate: the
paper states 3 pooled wells of a 12-well plate at 18e4 cells/well seeded.

If Fig 2A is nmol/mg,  use C_cell as printed.
If Fig 2A is ng/mg,    divide by MW first.

Only one of the two gives a V/m that a 12-well plate can produce.
"""
MW = {"PFBA": 214.04, "PFHxA": 314.05, "PFHpA": 364.06,
      "PFOA": 414.07, "PFNA": 464.08}
# Rosenmai Table 1: cellular uptake, % of total moles added
UPTAKE = {"PFBA": (0.14, 0.13, 0.10), "PFHxA": (0.19, 0.19, 0.14),
          "PFHpA": (0.16, 0.22, 0.21), "PFOA": (0.24, 0.27, 0.33),
          "PFNA": (0.15, 0.22, 0.31)}
# Rosenmai Fig 2A, as plotted (axis labelled nmol/mg protein)
FIG2A = {"PFBA": (18, 45, 126), "PFHxA": (26, 69, 164),
         "PFHpA": (28, 116, 297), "PFOA": (39, 127, 525),
         "PFNA": (24, 86, 458)}
DOSES = (10.0, 30.0, 100.0)          # uM

def V_over_m(C_cell_nmol_per_mg, uptake_pct, dose_uM):
    """mL of medium per mg of cell protein."""
    n_per_mg = C_cell_nmol_per_mg                     # nmol / mg
    added_per_L = dose_uM * 1000.0                    # nmol / L
    return n_per_mg / added_per_L / (uptake_pct / 100.0) * 1000.0   # mL/mg

print("V/m implied by each row  (mL of medium per mg of cell protein)")
print("%-7s %6s | %14s %14s" % ("cmpd", "dose", "if nmol/mg", "if ng/mg"))
both = {"nmol": [], "ng": []}
for c in FIG2A:
    for d, up, f in zip(DOSES, UPTAKE[c], FIG2A[c]):
        a = V_over_m(f, up, d)              # taking the axis at face value
        b = V_over_m(f / MW[c], up, d)      # dividing by molecular weight
        both["nmol"].append(a); both["ng"].append(b)
        print("%-7s %6g | %11.0f mL %11.2f mL" % (c, d, a, b))

print()
print("                     range of V/m     median")
for k, v in both.items():
    v2 = sorted(v)
    print("  as %-5s %12.1f - %-8.1f %8.1f mL/mg"
          % (k + "/mg", min(v), max(v), v2[len(v2) // 2]))

print()
print("What can the plate actually give?")
cells = 3 * 18e4                     # 3 pooled wells, 18e4 cells/well
for pg, src in ((675, "project spreadsheet"), (694, "Armitage AREc2 Table S4")):
    m = cells * pg * 1e-9            # mg protein
    print("  %d cells x %d pg protein/cell = %.3f mg   (%s)" % (cells, pg, m, src))
    for vol, lab in ((0.5, "0.5 mL/well"), (1.0, "1 mL/well, the usual"),
                     (2.0, "2 mL/well, brim-ful")):
        print("      %-22s -> V/m = %6.1f mL/mg" % (lab, 3 * vol / m))
