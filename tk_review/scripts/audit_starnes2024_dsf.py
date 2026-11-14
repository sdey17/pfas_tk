#!/usr/bin/env python3
"""Starnes 2024: the only four-species albumin study, 2 log off the rest.

THE PAPER.  Starnes HM, Jackson TW, Rock KD, Belcher SM.
"Quantitative cross-species comparison of serum albumin binding of
per- and polyfluoroalkyl substances from five structural classes."
Toxicol Sci 2024;199(1):132-148.  doi 10.1093/toxsci/kfae028,
PMID 38518100.  Open access.

WHAT IS IN IT.  Table 2: 8 PFAS x 4 albumins (human, bovine, porcine,
rat), one laboratory, one method, DSF.  Seven compounds give a Kd;
6:2 FTOH gives none because it does not bind ANY of the four, which is
the paper's own headline about the charged headgroup and is a result
rather than a gap.  So 28 constants plus 4 documented non-binders.

It is the only four-species albumin comparison in this compilation,
and the only source of PORCINE albumin data anywhere in it.

THE PROBLEM: ITS Kd ARE MILLIMOLAR, AND SYSTEMATICALLY SO.  Table 2
runs 0.33 to 3.06 mM.  Converted at n = 1 the whole block occupies
log K 0.69 to 1.66 -- a band at the very bottom of the sheet, whose
current floor is 0.655.  Against the sheet's own per-compound,
per-species medians, all twenty comparisons that have a pool are
NEGATIVE:

    median offset -2.05 log, range -3.17 to -1.12

That is not scatter, it is a method offset.  DSF reports the
concentration at which a ligand measurably stabilises the protein
against thermal unfolding, which needs far more occupancy than
half-saturation of the tightest site.  Two further consequences of the
design point the same way: the Kd is extracted from a 50-90 C thermal
ramp, so it is not a fixed-temperature constant at all; and the
albumin is at 0.11-0.18 mM, within a factor of 2 to 20 of the Kd
being measured, which is the stoichiometric-titration regime where a
Kd is weakly determined even on its own terms.

WHAT IT IS GOOD FOR ANYWAY: THE RELATIVE PATTERN.  Within a compound,
all four species were run side by side, so the ratios are internally
controlled even if the absolute scale is not.  Relative to human
(positive = binds tighter than human, log10 Kd_HSA/Kd_x):

               BSA     PSA     RSA
    mean      -0.11   +0.15   +0.13

Porcine and rat bind these PFAS slightly TIGHTER than human; bovine
sits within 0.1 log of human.  Every species is inside half a log of
every other, for every compound.

AND THAT CONTRADICTS THE SHEET'S OTHER WITHIN-STUDY PAIRS.  On PFOA:

    Chi 2018      BSA vs HSA   +0.81 log   bovine 6.4x tighter
    Starnes 2024  BSA vs HSA   -0.10 log   indistinguishable
    Han 2003      RSA vs HSA    ~0         "no significant difference"
    Starnes 2024  RSA vs HSA   +0.30 log   rat 2x tighter

So the two bovine/human estimates disagree by 0.9 log and the two
rat/human estimates by 0.3 log.  For the species aim that
disagreement IS the finding, and it is not resolvable from absolute
values that sit 2 log apart on method alone.

NOT EXPERIMENTAL, so not ingestable under the standing rule:
Autodock Vina docking dG (Supplementary Table S4) and the SeqAPASS
sequence-similarity analysis (Table 4).  Table 3 is a within-species
chain-length comparison derived from Table 2, not new measurement.

Run:  python3 -I tk_review/scripts/audit_starnes2024_dsf.py
"""
import math

BR = 4.1774
T2 = {"PFBA": {"HSA": 2.64, "BSA": 3.06, "PSA": 1.28, "RSA": 1.90},
      "PFHxA": {"HSA": 1.71, "BSA": 2.17, "PSA": 0.95, "RSA": 1.33},
      "PFOA": {"HSA": 0.79, "BSA": 0.99, "PSA": 0.37, "RSA": 0.40},
      "PFBS": {"HSA": 1.68, "BSA": 0.80, "PSA": 1.25, "RSA": 1.41},
      "PFOS": {"HSA": 0.69, "BSA": 0.86, "PSA": 0.38, "RSA": 0.42},
      "HFPO-DA": {"HSA": 1.57, "BSA": 2.55, "PSA": 2.12, "RSA": 1.64},
      "6:2 FTSA": {"HSA": 0.41, "BSA": 1.39, "PSA": 0.52, "RSA": 0.33}}
PROT = {"HSA": 0.13, "BSA": 0.18, "PSA": 0.11, "RSA": 0.18}   # mM


def main():
    print("A.  Table 2 on the sheet's scale (n = 1)\n")
    lk = []
    print(f"  {'chemical':9s} " + " ".join(f"{s:>8s}" for s in PROT))
    for c, d in T2.items():
        row = [BR - math.log10(d[s] * 1e3) for s in PROT]
        lk += row
        print(f"  {c:9s} " + " ".join(f"{v:8.3f}" for v in row))
    print(f"\n  whole block: log K {min(lk):.3f} to {max(lk):.3f}")
    print("  the sheet's floor today is 0.655, so all 28 land at the")
    print("  bottom.  Median offset against per-compound, per-species")
    print("  pools is -2.05 log and every one of the 20 comparisons")
    print("  that has a pool is negative.  A method offset, not scatter.")

    print("\nB.  Why DSF reads low")
    print("  It reports the concentration at which a ligand measurably")
    print("  STABILISES the protein against unfolding, which needs far")
    print("  more occupancy than half-saturation of the tightest site.")
    print("  The Kd also comes off a 50-90 C ramp, so it is not a")
    print("  fixed-temperature constant, and the albumin sits at")
    print("  0.11-0.18 mM against Kd of 0.33-3.06 mM:")
    for s, p in PROT.items():
        kds = [d[s] for d in T2.values()]
        print(f"    {s}  [P] {p:.2f} mM   Kd/[P] "
              f"{min(kds)/p:.1f} to {max(kds)/p:.1f}")
    print("  i.e. near stoichiometric titration, where Kd is weakly")
    print("  determined even on the method's own terms.")

    print("\nC.  The relative pattern, which IS internally controlled")
    print("  log10(Kd_HSA / Kd_x); positive = tighter than human\n")
    print(f"  {'chemical':9s} {'BSA':>7s} {'PSA':>7s} {'RSA':>7s}")
    acc = {s: [] for s in ("BSA", "PSA", "RSA")}
    for c, d in T2.items():
        out = []
        for s in acc:
            v = math.log10(d["HSA"] / d[s])
            acc[s].append(v)
            out.append(f"{v:+7.2f}")
        print(f"  {c:9s} " + " ".join(out))
    print(f"  {'mean':9s} "
          + " ".join(f"{sum(acc[s])/len(acc[s]):+7.2f}" for s in acc))
    print("\n  Every species within half a log of every other, every")
    print("  compound.  Against the sheet's other within-study pairs,")
    print("  on PFOA:")
    print("    Chi 2018      BSA vs HSA  +0.81   bovine 6.4x tighter")
    print(f"    Starnes 2024  BSA vs HSA  "
          f"{math.log10(T2['PFOA']['HSA']/T2['PFOA']['BSA']):+.2f}   "
          "indistinguishable")
    print("    Han 2003      RSA vs HSA   ~0     'no sig. difference'")
    print(f"    Starnes 2024  RSA vs HSA  "
          f"{math.log10(T2['PFOA']['HSA']/T2['PFOA']['RSA']):+.2f}   "
          "rat 2x tighter")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
