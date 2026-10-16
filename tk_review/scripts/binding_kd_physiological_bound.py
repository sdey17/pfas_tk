#!/usr/bin/env python3
"""
Which albumin dissociation constants are physiologically possible?

Reproduces report section 6.1a. Derivation, for one binding site:

    L + P <-> PL            Kd = [L][P] / [PL]
    f_bound = [PL] / ([PL] + [L])

substitute [PL] = [L][P]/Kd and the free-ligand term cancels:

    f_bound = [P] / ([P] + Kd)        ->      Kd = [P] x (1/f_bound - 1)

n sites per albumin scales [P] to n x [albumin]. The cancellation assumes
trace ligand, [L] << [P], which holds in plasma and fails in most of the
assays that produced the constants being scored.

Kd values are from OEHHA 2024 PFOA/PFOS PHG Table A6.2 (albumin rows), plus
Starnes et al. 2024 Table 2 for the DSF comparison.

Run:  python3 tk_review/scripts/binding_kd_physiological_bound.py
"""
MW_ALB = 66500.0        # g/mol
MW_PFOA = 414.07        # g/mol
ALB_UM = 600.0          # uM, from 40 g/L

# (label, species, PFAS, Kd in M) -- OEHHA Table A6.2, serum-albumin rows
A62 = [
    ("Beesoon & Martin 2015", "human", "n-PFOS", 8e-8),
    ("Bischel et al. 2010 (low)", "human", "PFOA", 0.7e-6),
    ("Chen & Guo 2009", "human", "PFOA", 3.7e-6),
    ("Bischel et al. 2010 (high)", "human", "PFOA", 5e-6),
    ("Messina et al. 2005 (step 1)", "human", "PFOA", 6.9e-6),
    ("MacManus-Spencer 2010 (low)", "bovine", "PFOA", 0.7e-5),
    ("Qin et al. 2010", "bovine", "PFOA", 2.1e-5),
    ("Wu et al. 2009 (step 1)", "human", "PFOA", 2.6e-5),
    ("MacManus-Spencer 2010 (high)", "bovine", "PFOA", 3e-5),
    ("Messina et al. 2005 (step 2)", "human", "PFOA", 3.2e-5),
    ("Chen & Guo 2009", "human", "PFOS", 4.5e-5),
    ("Beesoon & Martin 2015", "human", "n-PFOA", 1e-4),
    ("Li et al. 2009", "human", "PFOS", 2e-4),
    ("Han et al. 2003", "rat/human", "PFOA", 4e-4),
    ("Starnes et al. 2024 (DSF)", "human", "PFOA", 0.79e-3),
    ("Wu et al. 2009 (step 2)", "human", "PFOA", 1.3e-3),
    ("Klevens & Ellenbogen 1954", "bovine", "PFOA", 3.1e-3),
]


def f_bound(kd_M, n=1, alb_uM=ALB_UM):
    return 1.0 / (1.0 + (kd_M * 1e6) / (n * alb_uM))


def kd_needed(f, n=1, alb_uM=ALB_UM):
    return n * alb_uM * (1.0 / f - 1.0)


def main():
    print("[albumin] from %g g/L at MW %g = %.0f uM" % (40, MW_ALB, 40 / MW_ALB * 1e6))
    print()
    print("Kd required, by target binding and stoichiometry")
    print("%-10s %12s %12s %12s" % ("f_bound", "n=1", "n=3", "n=6"))
    for f in (0.95, 0.99, 0.995, 0.999):
        print("%-10s %9.1f uM %9.1f uM %9.1f uM"
              % ("%.1f%%" % (100 * f), kd_needed(f, 1), kd_needed(f, 3), kd_needed(f, 6)))
    print()
    print("OEHHA Table A6.2 albumin rows scored against >99%% bound")
    print("%-30s %-9s %-7s %10s %9s %9s" %
          ("source", "species", "PFAS", "Kd (uM)", "n=1", "n=6"))
    for lab, sp, ch, kd in sorted(A62, key=lambda x: x[3]):
        f1, f6 = f_bound(kd, 1), f_bound(kd, 6)
        mark = "" if f6 >= 0.99 else ("   <-- too weak" if f6 < 0.95 else "   ~")
        print("%-30s %-9s %-7s %10.2f %8.1f%% %8.1f%%%s"
              % (lab, sp, ch, kd * 1e6, 100 * f1, 100 * f6, mark))
    print()
    print("Trace-ligand assumption, [L] << [P]")
    print("%-36s %12s %12s" % ("", "[PFOA] uM", "L / P"))
    for lab, ngml in (("human background, 4 ng/mL", 4.0),
                      ("US 95th percentile, 10 ng/mL", 10.0),
                      ("Ronneby high exposure, 200 ng/mL", 200.0),
                      ("fluorochemical worker, 5000 ng/mL", 5000.0)):
        L = ngml / MW_PFOA          # ng/mL / (g/mol) = umol/L
        print("  %-34s %10.4f %12.2g" % (lab, L, L / ALB_UM))
    for lab, L, P in (("Starnes 2024 DSF assay", 1000.0, ALB_UM),
                      ("Chen 2020 assay vs 2 uM HSA", 80.0, 2.0)):
        print("  %-34s %10.1f %12.2g" % (lab, L, L / P))


if __name__ == "__main__":
    main()
