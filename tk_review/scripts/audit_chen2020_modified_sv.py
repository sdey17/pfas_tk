#!/usr/bin/env python3
"""Chen 2020: Tables 1 and 2 do not cohere, by a factor of 2.2 to 2.4.

THE PAPER.  Chen Y, et al.  Int J Environ Res Public Health
2020;17(4):1319.  doi 10.3390/ijerph17041319, PMID 32085632.  Open
access.  PFBA, PFHxA, PFOA and PFDA against HSA by fluorescence
quenching at 296, 303 and 310 K, PBS pH 7.4.

THE STUB said chemical "PFOA and PFDA" with ka "0.45e4-2.68e4 M^-1"
from a secondary, and its own note admitted why: "Table 2 Ka values
are not in the saved flattened text, so the Ka here remains the
Starnes 2026 secondary range".  Read from the PDF, that apparent range
is -- for the third source running -- TWO COMPOUNDS AND THREE
TEMPERATURES collapsed: 0.45e4 is PFOA at 310 K (Table 2: 0.4463e4)
and 2.68e4 is PFDA at 296 K (2.6788e4).  Its geometric midpoint,
Kd 91.06 uM and log K 2.218, is no compound at any temperature.

TABLE 2, modified Stern-Volmer association constants K_a (x10^4
L/mol):

           296 K     303 K     310 K
  PFOA     0.6153    0.5332    0.4463
  PFDA     2.6788    2.0093    1.4514

TABLE 1, ordinary Stern-Volmer quenching constants K_SV (x10^4
L/mol):

           296 K     303 K     310 K
  PFOA     1.328     1.203     1.076
  PFDA     1.731     1.647     1.431

THE TWO TABLES CONTRADICT EACH OTHER.  The paper's eq 2 is the Lehrer
modified Stern-Volmer form,

    F0/dF = 1/(f_a K_a [Q]) + 1/f_a

in which f_a "represents the fraction of accessible fluorescence".
Expanding it at small [Q] gives F0/F = 1 + f_a K_a [Q], so the
ordinary Stern-Volmer slope IS f_a K_a and therefore

    f_a = K_SV / K_a

which needs nothing but the two published tables -- no plot reading,
no assumption of ours.  It gives

  PFOA   296 K  2.158     PFDA   296 K  0.646
         303 K  2.256            303 K  0.820
         310 K  2.411            310 K  0.986

PFOA requires an accessible fraction above 2 at every temperature,
rising with temperature.  A fraction cannot exceed 1, and HSA has a
SINGLE tryptophan, Trp-214, which is what the paper's 295 nm
excitation selects -- so there is not even a second fluorophore for
the excess to hide in.  PFDA's values are physical throughout.

Read without the Lehrer expansion the conclusion is the same in
weaker words: Table 1 and Table 2 disagree by 2.2 to 2.4 for PFOA and
cannot both be right.

RELATION TO THE OTHER TWO.  This is the same method and the same
shape of failure as Qin 2010 (scripts/audit_qin2010_modified_sv.py),
PFOA unphysical and PFDA physical in both, and the same accessible
fraction test that held out Wu 2024.  Chen's is the worst of the
three: Qin's PFOA f_a runs 1.19 to 1.76 on BSA, which has two
tryptophans; Chen's runs 2.16 to 2.41 on HSA, which has one.

WHERE THE NUMBERS WOULD LAND at 310 K with n = 1 assumed:

    PFOA   Kd 224.06 uM   log K 1.827   -1.98 vs the PFOA median,
                                        and BELOW the lowest PFOA row
                                        now in the sheet (2.117)
    PFDA   Kd  68.90 uM   log K 2.339   -1.48 vs the PFDA median

CONDITIONS were already right in the stub and the PDF confirms them.
Fig. 1's caption: "[PFCAs] (x10^-6 mol/L) 1-10: 0, 2, 4, 8, 12, 16,
20, 30, 40 and 80. [HSA] = 2 x 10^-6 mol/L; Buffer: PBS, pH = 7.4."
So 2 uM HSA, PFAS to 80 uM, loading 40.  PFBA and PFHxA were tested
and gave no constants -- Fig. 1E shows both essentially flat, and the
FT-IR section says they "have little influence on the HSA structure".

Nothing is changed here.

Run:  python3 -I tk_review/scripts/audit_chen2020_modified_sv.py
"""
import math

BR = 4.1774
D = [("PFOA", 296, 1.328, 0.6153), ("PFOA", 303, 1.203, 0.5332),
     ("PFOA", 310, 1.076, 0.4463), ("PFDA", 296, 1.731, 2.6788),
     ("PFDA", 303, 1.647, 2.0093), ("PFDA", 310, 1.431, 1.4514)]
MEDIAN = {"PFOA": 3.806, "PFDA": 3.823}


def main():
    print("A.  Table 1 (K_SV) against Table 2 (K_a)\n")
    print("    eq 2:  F0/dF = 1/(f_a K_a [Q]) + 1/f_a")
    print("    at small [Q] this is F0/F = 1 + f_a.K_a.[Q], so the")
    print("    Stern-Volmer slope IS f_a.K_a and f_a = K_SV / K_a,")
    print("    from the two published tables alone.\n")
    print(f"  {'cmpd':5s} {'T':>4s} {'K_SV':>8s} {'K_a':>8s} {'f_a':>8s}")
    print("  " + "-" * 40)
    for c, T, sv, ka in D:
        fa = sv / ka
        print(f"  {c:5s} {T:4d} {sv*1e4:8.0f} {ka*1e4:8.0f} {fa:8.3f}"
              + ("   > 1, impossible" if fa > 1 else ""))
    print("\n  HSA has a SINGLE tryptophan, Trp-214, which the paper's")
    print("  295 nm excitation selects, so there is not even a second")
    print("  fluorophore for the excess to hide in.  Without the Lehrer")
    print("  expansion the statement is only weaker, not different:")
    print("  the two tables disagree by 2.2-2.4 for PFOA.")

    print("\nB.  At 310 K, n = 1 assumed")
    for c, T, sv, ka in D:
        if T != 310:
            continue
        kd = 1e6 / (ka * 1e4)
        lk = BR - math.log10(kd)
        print(f"    {c}  Kd {kd:7.2f} uM -> log K {lk:.3f}"
              f"   ({lk - MEDIAN[c]:+.2f} vs its median)")
    print("\n    PFOA at 1.827 would be the lowest PFOA row in the sheet;")
    print("    the lowest now is 2.117.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
