"""What K_a,1, K_a,2, K_a,3 are, and why they fall."""
import math
BR = 4.1774

print("=" * 68)
print("1.  THE DEFINITION: they are SEQUENTIAL, not parallel")
print("=" * 68)
print("""
   P   + L  <=>  PL      K_a,1 = [PL]   / ([P]  [L])
   PL  + L  <=>  PL2     K_a,2 = [PL2]  / ([PL] [L])
   PL2 + L  <=>  PL3     K_a,3 = [PL3]  / ([PL2][L])

   K_a,2 is NOT "the affinity of the second site".  It is the
   affinity for adding a second ligand to a protein that ALREADY
   HAS ONE.  Every step starts from the previous complex.

   Chi's eq 3 is exactly this.  With R_j = [PL_j]/[P],

        K_a,j = R_j / (R_{j-1} [L])
              = ([PL_j]/[P]) / (([PL_{j-1}]/[P]) [L])
              = [PL_j] / ([PL_{j-1}] [L])        <- same thing

   Native ESI-MS can do this because it SEES PL1, PL2, PL3 ... as
   separate peaks and measures their relative abundances directly.
""")

print("=" * 68)
print("2.  THREE DIFFERENT 'CONSTANTS' PEOPLE WRITE")
print("=" * 68)
print("""
   stepwise    K_a,j    one more ligand onto PL_{j-1}   <- Chi's
   cumulative  beta_j   j ligands onto bare P
                        beta_j = K_a,1 x K_a,2 x ... x K_a,j
   intrinsic   k        one ligand onto ONE empty site,
                        if all n sites were identical
""")
KA = {"BSA-PFOA": (2.1e5, 1.1e5, 5.5e4), "BSA-PFOS": (7.6e5, 2.4e5, 1.2e5),
      "HSA-PFOA": (1.1e5, 6.1e4, 5.2e4), "HSA-PFOS": (1.1e5, 7.6e4, 5.4e4)}
print(f"   Chi's 1:2 column, and the cumulative constants it implies:")
print(f"   {'system':10s} {'K_a,1':>9s} {'K_a,2':>9s} {'K_a,3':>9s}"
      f" | {'beta_2':>10s} {'beta_3':>10s}")
for s, (k1, k2, k3) in KA.items():
    print(f"   {s:10s} {k1:9.3g} {k2:9.3g} {k3:9.3g} |"
          f" {k1*k2:10.3g} {k1*k2*k3:10.3g}")

print()
print("=" * 68)
print("3.  WHY THEY FALL -- and most of it is just COUNTING")
print("=" * 68)
print("""
   Take n IDENTICAL, INDEPENDENT sites, each with intrinsic k.
   No cooperativity at all.  The stepwise constants STILL fall:

        K_a,j = k x (n - j + 1) / j

   (n-j+1) empty sites to land on, j filled ones to fall off.
   With n = 8, which is what Chi's spectra reach:

        K_a,1 = 8k        K_a,2 = 3.5k      K_a,3 = 2k

   so the pattern 1 : 0.44 : 0.25 is what you get from pure
   statistics.  Any REAL cooperativity is the departure from it.
""")
n = 8
stat = [(n - j + 1) / j for j in (1, 2, 3)]
stat = [x / stat[0] for x in stat]
print(f"   statistical, n=8        1.000  {stat[1]:.3f}  {stat[2]:.3f}")
print()
print(f"   {'system':10s} {'K2/K1':>7s} {'K3/K1':>7s}   verdict")
for s, (k1, k2, k3) in KA.items():
    r2, r3 = k2 / k1, k3 / k1
    if r3 < stat[2] * 0.85:
        v = "falls FASTER -> negative cooperativity"
    elif r3 > stat[2] * 1.15:
        v = "falls SLOWER than counting alone"
    else:
        v = "matches counting -> no cooperativity needed"
    print(f"   {s:10s} {r2:7.3f} {r3:7.3f}   {v}")

print()
print("=" * 68)
print("4.  WHICH ONE DOES THE BRIDGE WANT?  (and a free lunch)")
print("=" * 68)
print("""
   log K is a TRACE-LIMIT quantity: vanishing ligand, so only the
   FIRST binding event matters.  That is K_a,1.  K_a,2 and K_a,3
   describe a protein that is already loaded, which is not the
   regime a partition coefficient is defined in.

   Now the free lunch.  You might worry that splitting K_a,1 into
   n identical sites changes the answer.  It does not:

        k = K_a,1 / n          intrinsic site constant
        Kd_site = 1/k = n / K_a,1

        log K = 4.1774 + log10(n) - log10(Kd_site)
              = 4.1774 + log10(n) - log10(n x 1e6 / K_a,1)
              = 4.1774 - log10(1e6 / K_a,1)       <- n cancels
""")
for s, (k1, _, _) in KA.items():
    direct = BR - math.log10(1e6 / k1)
    for nn in (1, 4, 8):
        kd = nn * 1e6 / k1
        assert abs(BR + math.log10(nn) - math.log10(kd) - direct) < 1e-12
    print(f"   {s:10s} K_a,1 {k1:7.3g} -> log K {direct:.3f}"
          f"   (identical for n = 1, 4 or 8)")
print("""
   So the n = 1 in the sheet is a LABEL, not an assumption that
   moves the number.  At trace the protein's capacity is n x k,
   and K_a,1 already is that product.
""")
