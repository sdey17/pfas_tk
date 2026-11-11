#!/usr/bin/env python3
"""How far does one albumin dataset move when you refit it legitimately?

Bischel BJ, MacManus-Spencer LA, Luthy RG.  "Noncovalent interactions of
long-chain perfluoroalkyl acids with serum albumin."  Environ Sci Technol
2010;44(13):5263-5269.  doi 10.1021/es101334s.  SI Tables S2, S3, S4.

WHY THIS MATTERS HERE.  what_drives_albumin_spread.py ranks the sources of
the 1.08 log scatter in albumin log K and finds the LABORATORY the biggest
identifiable term at 45%, with ~0.59 log left inside a single laboratory on
a single compound.  That residue has no obvious owner.  One candidate is
the fit: a binding isotherm is reduced to (Ka, n) by a choice of model,
weighting and concentration window, and those are judgement calls.

Bischel's SI happens to measure exactly that.  It reports the SAME two
dialysis datasets -- PFOA and PFNA on 1 uM BSA -- refitted many ways:

  Table S2  PFNA, two-class, SEVEN weighting factors
  Table S3  PFOA and PFNA, two-class, two weightings
  Table S4  PFOA and PFNA, one-class, three model/window combinations

Every row is the same experiment.  None is wrong; the paper reports one
from each family in its Table 2.  So the spread across them is fitting
choice with the chemistry, the bench, the day and the operator all held
fixed -- a floor on what no amount of harmonising can remove.

NOTHING HERE IS INGESTED.  These are refits of data the database already
holds once, through Table 2.  Adding them would count one experiment
seventeen times and would inflate exactly the spread this project is
trying to explain.  The script reports; it does not write.

Run:  python3 -I tk_review/scripts/bischel2010_refit_spread.py
"""
import collections
import math
import statistics as st

CONV = 4.1774              # log10(1e9 / 66463), the Kd -> K bridge

# (table, compound, model, weighting, Ka1 [uM^-1], n1, R2, [L]:[P] window)
FITS = [
    ("S4", "PFOA", "one-class", "none",          1.1,  2.1, 0.988, "0.04-1"),
    ("S4", "PFOA", "one-class", "none",          0.20, 4.2, 0.767, "0.04-5"),
    ("S4", "PFOA", "one-class", "1/[free]",      1.7,  1.5, 0.894, "0.04-5"),
    ("S3", "PFOA", "two-class", "none",          3.0,  1.0, 0.917, "0.04-70"),
    ("S3", "PFOA", "two-class", "1/[free]",      3.1,  0.9, 0.891, "0.04-70"),
    ("S4", "PFNA", "one-class", "none",          1.3,  4.6, 0.900, "0.02-1"),
    ("S4", "PFNA", "one-class", "none",          1.2,  4.6, 0.948, "0.02-5"),
    ("S4", "PFNA", "one-class", "1/[free]",      1.6,  4.2, 0.930, "0.04-5"),
    ("S3", "PFNA", "two-class", "none",          1.7,  3.9, 0.968, "0.02-120"),
    ("S3", "PFNA", "two-class", "1/[free]",      1.7,  3.8, 0.940, "0.02-120"),
    ("S2", "PFNA", "two-class", "none",          3.3,  2.9, 0.953, "0.02-120"),
    ("S2", "PFNA", "two-class", "1/[free]",      1.9,  3.5, 0.938, "0.02-120"),
    ("S2", "PFNA", "two-class", "1/[total]",     1.7,  3.6, 0.955, "0.02-120"),
    ("S2", "PFNA", "two-class", "1/nu",          2.0,  3.3, 0.959, "0.02-120"),
    ("S2", "PFNA", "two-class", "1/[free]^2",    2.4,  3.0, 0.766, "0.02-120"),
    ("S2", "PFNA", "two-class", "1/[total]^2",   1.6,  3.6, 0.918, "0.02-120"),
    ("S2", "PFNA", "two-class", "1/nu^2",        1.6,  3.5, 0.946, "0.02-120"),
]

# What the main text selected, and what the database therefore holds.
TABLE2 = {("PFOA", "one-class"): (0.20, 4.3), ("PFOA", "two-class"): (1.4, 1.4),
          ("PFNA", "one-class"): (1.1, 4.6), ("PFNA", "two-class"): (3.3, 2.9)}


def logk(ka_uM, n):
    """log K [L_water/kg_protein] from Ka in uM^-1 and a site count."""
    return CONV + math.log10(n) + math.log10(ka_uM)


def main():
    print(__doc__.split("Run:")[0].strip().split("\n\n")[0])
    print("\n" + "=" * 72)
    print("EVERY FIT BISCHEL PUBLISHED OF THESE TWO DATASETS")
    print("=" * 72)
    print("%-4s %-6s %-10s %-12s %8s %5s %7s %10s %8s" % (
        "tbl", "PFAS", "model", "weighting", "Ka uM-1", "n", "R2",
        "window", "log K"))
    for t, c, model, w, ka, n, r2, win in FITS:
        star = " *" if TABLE2.get((c, model)) == (ka, n) else ""
        print("%-4s %-6s %-10s %-12s %8.2f %5.1f %7.3f %10s %8.3f%s" % (
            t, c, model, w, ka, n, r2, win, logk(ka, n), star))
    print("  * the fit the main text's Table 2 reports, and the database "
          "holds")

    print("\n" + "=" * 72)
    print("THE SPREAD, WITH EVERYTHING BUT THE FIT HELD FIXED")
    print("=" * 72)
    by = collections.defaultdict(list)
    for t, c, model, w, ka, n, r2, win in FITS:
        by[c].append(logk(ka, n))
        by[(c, model)].append(logk(ka, n))
    print("%-22s %4s %8s %8s %8s   %s" % (
        "held fixed", "n", "min", "max", "span", "what varies"))
    for c in ("PFOA", "PFNA"):
        for model in ("one-class", "two-class"):
            v = by[(c, model)]
            if len(v) > 1:
                print("%-22s %4d %8.2f %8.2f %8.2f   %s" % (
                    f"{c}, {model}", len(v), min(v), max(v),
                    max(v) - min(v), "weighting and window only"))
        v = by[c]
        print("%-22s %4d %8.2f %8.2f %8.2f   %s" % (
            f"{c}, either model", len(v), min(v), max(v), max(v) - min(v),
            "model, weighting and window"))

    spans = [max(by[c]) - min(by[c]) for c in ("PFOA", "PFNA")]
    print("\n" + "=" * 72)
    print("READING")
    print("=" * 72)
    print("  One dataset, one bench, one operator, one protein batch.")
    print("  Refitting alone moves log K by %.2f log (PFOA) and %.2f "
          "(PFNA)," % (spans[0], spans[1]))
    print("  a factor of %.0f and %.0f." % (10 ** spans[0], 10 ** spans[1]))
    print()
    print("  Against the project's own decomposition:")
    print("    1.08 log   total scatter inside a compound")
    print("    0.59 log   what survives after the laboratory is known")
    print("    %.2f log   what refitting ONE dataset can do, here"
          % max(spans))
    print()
    print("  So for PFOA the fit alone covers most of the within-laboratory")
    print("  residue.  The scatter the literature shows is not all")
    print("  measurement: a large part is how the isotherm was reduced to")
    print("  two numbers, and that part is invisible in any table that")
    print("  prints one (Ka, n) pair per compound.")
    print()
    print("  Note which way the model pushes.  For PFOA the two-class fits")
    print("  give the HIGHEST log K and the lowest n -- n = 0.9 to 1.0,")
    print("  against 1.5 to 4.2 for one-class on the same data.  n and Ka")
    print("  trade off against each other, and log K carries both, so the")
    print("  compensation is only partial.")
    print()
    print("  PFNA moves less (%.2f log) because its one-class and two-class"
          % min(spans))
    print("  fits happen to agree on n near 3-4.6.  That is luck, not a")
    print("  property of the method.")


if __name__ == "__main__":
    main()
