#!/usr/bin/env python3
"""Beesoon & Martin 2015 read from the figures, not the abstract.

THE PAPER.  Beesoon S, Martin JW.  "Isomer-Specific Binding Affinity
of Perfluorooctanesulfonate (PFOS) and Perfluorooctanoate (PFOA) to
Serum Proteins."  Environ Sci Technol 2015;49(9):5722-5731.
doi 10.1021/es505399w.  PMID 25826685.

WHAT IT COVERS.  PFOS and PFOA only -- no PFNA, no PFHxS, nothing
else.  What it has that no other source in this compilation has is
the ISOMER split: linear plus three branched (3m, 4m, 5m) for each,
by ultrafiltration against 6 uM HSA at 37 C, each isomer as a pure
standard.  Eight Scatchard panels, Figures 1 and 2.  Seven carry a
Kd; 4m-PFOA does not (p > 0.05, R2 = 0.53, "binding was very weak").

Our four existing rows came from the ABSTRACT (db source_table says
so) and collapse the branched isomers into "branched PFOS isomers
(range)" at an assumed n = 1.  The figures give each isomer its own
Kd AND its own n, and none of the branched n is 1: they run 0.26 to
0.74.

THREE THINGS THE FIGURES SAY THAT THE ABSTRACT DOES NOT.

(A) Six of the eight boxes follow exactly from their own printed
    regression, read as n = 1/slope and Kd = intercept/slope, which
    is what the Scatchard form the paper prints requires.  Two do
    not:

      5m-PFOS      the printed fit has a NEGATIVE intercept
                   (y = 2.6822x - 0.0004), so Kd = intercept/slope is
                   negative -- not a dissociation constant at all.
                   The box reports +9(+/-5)e-5 M.  Its slope gives
                   n = 0.37; the box says 0.6 +/- 0.03.  Neither
                   number in that box can be got from its own plot.

      linear PFOA  slope 0.6076 gives n = 1.65, the box says
                   2 +/- 0.06.  The Kd does reproduce.  The text
                   says only that n was "close to 2".

(B) THE AXIS UNITS AND THE REPORTED Kd DISAGREE BY 1000x.  Both axes
    of all eight panels are labelled (mM).  Kd = intercept/slope
    therefore carries mM.  Every Kd is reported in M.

    The mass balance decides it, and it is not close.  Methods: each
    isomer was spiked "from 0.002 mM to 0.012 mM" against 6 uM HSA,
    and the paper states that at the top linear-PFOS point "the ratio
    of PFOS to HSA did not exceed 2".  Taking r = LF/y off the two
    linear panels and closing total = bound + free:

      reading the axes as mM   linear PFOA totals 3.8, 6.3, 7.9,
                               10.1, 11.9 uM -- five points spanning
                               the stated 2-12 uM window, the top
                               within 1% of 12.
      reading the axes as M    the same totals are 54, 256, 657,
                               1309, 2459 uM, up to 200x the window
                               and far more ligand than the donated
                               stocks (107-203 uM) could supply.

    So the free concentrations are micromolar, the axes are mM as
    labelled, and the Kd values as printed are 1000x too large.  The
    branched panels make the same point harder still: 3m-PFOA reaches
    0.015 on its x axis, which as M would be 15 mM of free ligand.

    This does NOT touch the paper's conclusion.  Linear binds more
    tightly than branched by the same factor either way, because both
    sides of every comparison scale together.  It touches only the
    absolute constants, which is exactly what this compilation
    ingests.

(C) The whole-serum experiments (technical mixtures in calf and human
    serum) give relative enrichment of branched isomers in the
    ultrafiltrate -- ratios of 1.3-1.8 -- and no Kd.  Not usable as
    binding constants.  SI Table S1 is a literature compilation (Han
    2003, Hebert 2010, Wu 2009, Chen 2009), all already in our
    database, and is worth noting for one reason: the authors place
    their own linear PFOA at 1e-4 M beside those studies' 2.2e-4 to
    3.8e-4 M and the agreement looks good.  Under the axis-consistent
    reading it would be 1.3e-7 M and would not.  That agreement is
    the one argument on the other side, and it is a soft one against
    a mass balance that fails by three orders.

NOTHING IS REPAIRED HERE.  Consistent with the Peng 2024 and Wu 2024
decisions, the published values are what the sheet carries and the
defect is recorded.  Part C prints both readings so the cost of the
choice is visible: about +3 log per row.

Run:  python3 -I tk_review/scripts/audit_beesoon2015_scatchard.py
"""
import math
BR = 4.1774   # log10(1e9/66463)

# panel: (slope, intercept, boxed Kd in M, boxed n)   None = not reported
P = {
 "linear PFOS": (0.4888,  4e-8,   8e-8,  2.0),
 "3m-PFOS":     (1.9554,  7e-4,   4e-4,  0.5),
 "4m-PFOS":     (3.8316,  3e-4,   8e-5,  0.3),
 "5m-PFOS":     (2.6822, -4e-4,   9e-5,  0.6),
 "linear PFOA": (0.6076,  8e-5,   1e-4,  2.0),
 "3m-PFOA":     (1.9077,  8e-4,   4e-4,  0.5),
 "4m-PFOA":     (0.3822,  1.08e-2, None, None),
 "5m-PFOA":     (1.3534,  4e-4,   3e-4,  0.7),
}
print("A.  Does each box follow from its own printed regression?")
print("    Scatchard as the paper prints it:  LF/r = LF/n + Kd/n")
print("    so  n = 1/slope   and   Kd = intercept/slope\n")
print(f"  {'isomer':12s} {'slope':>8s} {'intercept':>10s} | {'n=1/s':>7s} "
      f"{'boxed n':>8s} | {'Kd=i/s':>10s} {'boxed Kd':>10s} {'ok?':>5s}")
print("  " + "-"*78)
for k,(s,i,kd,n) in P.items():
    nd, kdd = 1/s, i/s
    okn = "-" if n is None else ("yes" if abs(nd-n) < 0.06 else "NO")
    okk = "-" if kd is None else ("yes" if kd and 0.7 < kdd/kd < 1.4 else "NO")
    print(f"  {k:12s} {s:8.4f} {i:10.2e} | {nd:7.3f} "
          f"{('--' if n is None else f'{n:.1f}'):>8s} | {kdd:10.2e} "
          f"{('--' if kd is None else f'{kd:.0e}'):>10s} "
          f"{okn+'/'+okk:>5s}")

print("\nB.  What do the AXES say the numbers are?")
print("    Both axes are labelled (mM).  Methods: isomers spiked")
print("    'from 0.002 mM to 0.012 mM' = 2 to 12 uM, HSA = 6 uM.")
print("    r = LF/y, bound = r x 6 uM, total = bound + free.\n")
# points read off the two linear panels at 300 dpi
PTS = {
 "linear PFOA": [(5e-5,8e-5),(2.5e-4,2.5e-4),(6.5e-4,5.4e-4),
                 (1.3e-3,8.9e-4),(2.45e-3,1.55e-3)],
 "linear PFOS": [(3e-8,3.2e-8),(5e-8,5.2e-8),(7.5e-8,7.8e-8),
                 (3.4e-7,2.1e-7),(7.2e-7,3.8e-7)],
}
for k,pts in PTS.items():
    print(f"  {k}:")
    print(f"    {'free':>12s} {'r':>6s} {'bound uM':>9s} {'total uM':>9s}"
          f"   {'total if axes were M':>22s}")
    for x,y in pts:
        r = x/y
        b = r*6.0
        print(f"    {x:12.2e} {r:6.2f} {b:9.2f} {b+x*1e3:9.2f}"
              f"   {b+x*1e6:19.0f} uM")
    print()
print("    Under mM every total lands inside the stated 2-12 uM window,")
print("    the top point within 1% of 12.  Under M the totals are 1000x")
print("    the window and 1000x more ligand than the stocks held.")
print("    => the axes are mM, so Kd = intercept/slope is in mM too.\n")

print("C.  What the two readings do to log K")
print(f"  {'isomer':12s} {'Kd pub (uM)':>12s} {'logK pub':>9s} | "
      f"{'Kd axis (uM)':>13s} {'logK axis':>10s}")
print("  " + "-"*64)
for k,(s,i,kd,n) in P.items():
    if kd is None: continue
    nn = n
    kd_pub = kd*1e6                      # M -> uM
    kd_axis = (i/s)*1e3                  # mM -> uM
    lk = BR + math.log10(nn) - math.log10(kd_pub)
    if kd_axis <= 0:
        print(f"  {k:12s} {kd_pub:12.4g} {lk:9.3f} | "
              f"{kd_axis:13.4g} {'NEGATIVE':>10s}   <- unusable either way")
        continue
    la = BR + math.log10(1/s) - math.log10(kd_axis)
    print(f"  {k:12s} {kd_pub:12.4g} {lk:9.3f} | {kd_axis:13.4g} {la:10.3f}")
print()
print("  The mM/M factor alone is exactly +3 log.  The spread above")
print("  (+2.80 to +3.06) is that 3 log plus the difference between")
print("  each boxed value and its own regression -- biggest for")
print("  linear PFOA, where the box says n = 2 and the slope says 1.65.")
