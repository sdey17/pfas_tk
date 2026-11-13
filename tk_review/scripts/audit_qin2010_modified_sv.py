#!/usr/bin/env python3
"""Qin 2010 read from the primary: six constants, and an f_a above 1.

THE PAPER.  Qin P, Liu R, Pan X, Fang X, Mou Y.  "Impact of Carbon
Chain Length on Binding of Perfluoroalkyl Acids to Bovine Serum
Albumin Determined by Spectroscopic Methods."  J Agric Food Chem
2010;58(9):5561-5567.  doi 10.1021/jf100412q, PMID 20397730.
PFPA, PFOA and PFDA against BSA, pH 7.40 phosphate, three
temperatures.

WHAT THE STUB SAID.  One row, source_table "abstract; Ka range from
Starnes & Belcher 2026 Table 1", chemical "PFOA and PFDA (PFPeA
tested, little effect)", ka "0.22e5-6.85e5 M^-1", Kd 8.146 uM taken as
the geometric midpoint, n = 1 assumed, log K 3.266.

That apparent range is TWO COMPOUNDS AND THREE TEMPERATURES collapsed:
0.22e5 is PFOA at 320 K (Table 2: 2.17e4) and 6.85e5 is PFDA at 300 K
(6.846e5).  Its midpoint describes no compound at any temperature.

WHAT TABLE 2 ACTUALLY HOLDS.  Six modified Stern-Volmer constants:

           300 K        310 K        320 K
  PFOA   4.36e4       2.93e4       2.17e4     M^-1
  PFDA   6.846e5      4.417e5      2.722e5    M^-1

All six reproduce exactly from their own Fig. 3 fits.  Eq 2 is

    F0/dF = 1/(f_a K_a [Q]) + 1/f_a

so a fit y = b + m/[Q] gives f_a = 1/b and K_a = b/m, and the six
printed fits return Table 2 to four figures.

THE DEFECT: f_a EXCEEDS 1 FOR PFOA AT EVERY TEMPERATURE.

    PFOA   300 K  f_a 1.190     PFDA   300 K  f_a 0.437
           310 K      1.449            310 K      0.551
           320 K      1.761            320 K      0.730

f_a is the FRACTION of accessible fluorescence -- the paper's own
words, "f_a is the fraction of accessible fluorescence".  A fraction
cannot exceed 1.  BSA has two tryptophans, Trp-134 and Trp-212, both
named in this paper.  PFDA's three values are physical; PFOA's three
are not, and they get worse with temperature.

This is the same test that held out Wu 2024, whose accessible fraction
exceeded 1 for all ten compounds.  The difference in degree matters:
Wu's constants were also 10^1-10^2 M^-1, three to four logs below
every other HSA study, whereas these are 10^4-10^5 M^-1, an ordinary
range.  So PFOA here fails one check, not two.

WHAT K_a IS.  Not a binding constant: the paper calls it "the
effective quenching constant for the accessible fluorophores".  It has
units of L/mol and is dimensionally a Ka, but it describes how
accessible the tryptophans are to quencher, not how many ligands sit
on the protein.

NO STOICHIOMETRY IS REPORTED.  Checked: the paper has a Stern-Volmer
plot (Table 1, K_sv) and a modified Stern-Volmer plot (Table 2, K_a)
and NO double-logarithmic plot, so there is no fitted n anywhere in
it.  Unlike Klevens 1954 and Sheng 2020, where the sheet's "the paper
reports none" turned out to be false, here it is TRUE.

PFPA (perfluoropentanoic acid, C5) was tested and gave no constant:
"PFPA had little effect on BSA".  That is a result, not a gap.

CONDITIONS, from the Fig. 1 caption and Materials and Methods:
BSA 1 uM throughout; PFOA dosed 0-20 uM (ratio 20), PFDA 0-7 uM
(ratio 7), PFPA 0-10 uM; 0.02 mol/L NaH2PO4-Na2HPO4 at pH 7.40;
300, 310 and 320 K.

WHERE THE NUMBERS WOULD LAND, at 310 K with n = 1 assumed:

    PFOA   Kd 34.13 uM   log K 2.644   -1.20 log vs the PFOA median
    PFDA   Kd  2.26 uM   log K 3.823   -0.14 log vs the PFDA median

Nothing is changed here.

Run:  python3 -I tk_review/scripts/audit_qin2010_modified_sv.py
"""
import math

BR = 4.1774
# (compound, T, Fig. 3 intercept, Fig. 3 slope, Table 2 Ka)
F3 = [("PFOA", 300, 0.840, 1.926e-5, 4.36e4),
      ("PFOA", 310, 0.690, 2.358e-5, 2.93e4),
      ("PFOA", 320, 0.568, 2.619e-5, 2.17e4),
      ("PFDA", 300, 2.289, 3.344e-6, 6.846e5),
      ("PFDA", 310, 1.814, 4.107e-6, 4.417e5),
      ("PFDA", 320, 1.370, 5.033e-6, 2.722e5)]
MEDIAN = {"PFOA": 3.846, "PFDA": 3.966}


def main():
    print("A.  Do Table 2's constants follow from Fig. 3's own fits?")
    print("    eq 2:  F0/dF = 1/(f_a K_a [Q]) + 1/f_a")
    print("    a fit y = b + m/[Q] gives f_a = 1/b and K_a = b/m\n")
    print(f"  {'cmpd':5s} {'T':>4s} {'b':>7s} {'m':>11s} "
          f"{'f_a=1/b':>8s} {'K_a=b/m':>11s} {'Table 2':>11s}")
    print("  " + "-" * 60)
    for c, T, b, m, ka in F3:
        fa, k = 1 / b, b / m
        bad = "   f_a > 1" if fa > 1 else ""
        print(f"  {c:5s} {T:4d} {b:7.3f} {m:11.3e} {fa:8.3f} "
              f"{k:11.4g} {ka:11.4g}{bad}")
    print("\n  All six return Table 2 to four figures.")

    print("\nB.  f_a is a FRACTION and cannot exceed 1")
    print("    The paper: 'f_a is the fraction of accessible "
          "fluorescence'.")
    print("    BSA has two tryptophans, Trp-134 and Trp-212, both named")
    print("    in this paper.  PFOA returns 1.190, 1.449 and 1.761, and")
    print("    the excess grows with temperature.  PFDA's 0.437, 0.551")
    print("    and 0.730 are physical.")
    print("\n    Same test that held out Wu 2024.  Degree differs: Wu's")
    print("    constants were ALSO 10^1-10^2 M^-1, three to four logs")
    print("    below every other HSA study.  These are 10^4-10^5 M^-1,")
    print("    an ordinary range, so PFOA fails one check here, not two.")

    print("\nC.  At 310 K, n = 1 assumed (the paper reports no n at all)")
    for c, T, b, m, ka in F3:
        if T != 310:
            continue
        kd = 1e6 / ka
        lk = BR - math.log10(kd)
        print(f"    {c}  K_a {ka:.4g} M^-1 -> Kd {kd:7.3f} uM -> "
              f"log K {lk:.3f}   ({lk - MEDIAN[c]:+.2f} vs its median)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
