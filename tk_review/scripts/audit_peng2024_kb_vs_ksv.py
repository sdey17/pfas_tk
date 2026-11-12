#!/usr/bin/env python3
"""Peng 2024's K_b is not a binding constant where n is not 1.

Peng M, Xu Y, Wu Y, Cai X, Zhang W, Zheng L, Du E, Fu J.  "Binding
Affinity and Mechanism of Six PFAS with Human Serum Albumin: Insights
from Multi-Spectroscopy, DFT and Molecular Dynamics Approaches."
Toxics 2024;12(1):43.  doi 10.3390/toxics12010043.  Table 2.

WHAT PROMPTED THIS.  Peng's HFPO-TA at 310 K converts to log K 7.86, the
highest value anywhere in this compilation -- about 200x above the
highest DIRECTLY MEASURED albumin coefficient (~5.5).  plot_peng2024_
vanthoff.py already established that the implausible enthalpies behind it
are in the measurements rather than the arithmetic.  This asks a
different question: is K_b the quantity the conversion assumes it is?

THE TEST, and it needs no data this table does not already contain.
Table 2 prints a Stern-Volmer constant K_SV and a binding constant K_b
side by side, from the same titration.  For static 1:1 quenching -- which
is what the paper argues for, and what its K_q near 1e12 against the 2e10
diffusion limit supports -- the two describe the same complex and should
agree.

WHAT IT FINDS.  K_b/K_SV runs from 0.3 to 116,667 over the eighteen rows,
and log10 of that ratio is a straight line in n (r2 ~0.95).  The whole
divergence is n.

WHY.  K_b comes from the double-logarithmic plot

    log[(F0 - F)/F] = log K_b + n log[Q]

where K_b is the INTERCEPT and n the SLOPE.  Fitting n > 1 raises the
intercept exponentially: that is arithmetic, not affinity.  It also means
K_b carries units of (L/mol)^n.  HFPO-TA at 310 K has n = 2.09, so its
"2.31e9" is (L/mol)^2.09.  It is not the same quantity as DFSA's 1.52e3
at n = 0.88, and it is not a binding constant.

WHAT FOLLOWS.  db/protein_binding.csv already refuses a quenching-derived
n for conversion (build_binding_affinity_table, n_provenance ->
"quenching fit").  albumin_logK_shareable.xlsx deliberately revives those
rows, on their own route, because a comparison across the whole
literature was wanted.  Nothing here changes that: the rows stay, the
route already says what they are, and this script says how far the
quantity has drifted from a binding constant so a reader can decide.

Run:  python3 -I tk_review/scripts/audit_peng2024_kb_vs_ksv.py
"""
import math

R = 8.314
CONV = 4.1774

# Table 2: compound -> [(T K, K_SV x1e4 L/mol, K_b, n, dG kJ/mol)]
T2 = {
 "PFNA":    [(298, 2.52, 7.81e6, 1.51, -39.32), (304, 2.35, 2.02e5, 1.18, -30.86),
             (310, 1.27, 9.9e4, 1.10, -29.65)],
 "HFPO-TA": [(298, 1.78, 3.7e6, 1.49, -37.47), (304, 1.95, 8.54e7, 1.78, -46.16),
             (310, 1.98, 2.31e9, 2.09, -55.57)],
 "PFOA":    [(298, 1.39, 2.27e5, 1.26, -30.55), (304, 1.46, 2.23e6, 1.31, -32.85),
             (310, 1.47, 1.98e6, 1.47, -37.67)],
 "PFO3DA":  [(298, 1.05, 1.59e5, 1.25, -29.67), (304, 1.03, 4.55e4, 1.33, -27.11),
             (310, 0.87, 2.99e4, 1.10, -26.56)],
 "PFHpA":   [(298, 0.69, 4.53e3, 0.96, -20.86), (304, 0.53, 5.35e5, 1.42, -33.34),
             (310, 0.32, 7.47e5, 1.50, -34.86)],
 "DFSA":    [(298, 0.58, 1.52e3, 0.88, -18.15), (304, 0.51, 4.93e3, 0.98, -21.49),
             (310, 0.45, 2.07e4, 1.13, -25.62)],
}
TAU0 = 1e-8          # s, the fluorescence lifetime assumed for a biopolymer
DIFFUSION = 2e10     # L/mol/s, the maximum diffusion-controlled rate


def fit(pts):
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    b = sum((x - mx) * (y - my) for x, y in pts) / sxx
    ssr = sum((y - (my + b * (x - mx))) ** 2 for x, y in pts)
    sst = sum((y - my) ** 2 for _, y in pts)
    return b, (1 - ssr / sst if sst else float("nan"))


def main():
    print("1.  Is the quenching static?  K_q against the diffusion limit\n")
    kq = [k * 1e4 / TAU0 for v in T2.values() for _, k, _, _, _ in v]
    print("    K_q = K_SV / tau0 runs %.2e to %.2e L/mol/s," % (min(kq), max(kq)))
    print("    against a diffusion maximum of %.0e. So yes, static --" % DIFFUSION)
    print("    a ground-state complex, which is what makes K_SV and K_b")
    print("    descriptions of the SAME thing and lets them be compared.\n")

    print("2.  K_b against K_SV, same row, same titration\n")
    print("    %-9s %4s %10s %11s %11s %6s %8s" % (
        "PFAS", "T", "K_SV", "K_b", "K_b/K_SV", "n", "log K"))
    pts, worst = [], (0, None)
    for c, rows in T2.items():
        for T, ksv4, kb, n, dg in rows:
            ksv = ksv4 * 1e4
            lk = CONV + math.log10(n) - math.log10(1e6 / kb)
            pts.append((n, math.log10(kb / ksv)))
            if kb / ksv > worst[0]:
                worst = (kb / ksv, (c, T, n, lk))
            print("    %-9s %4d %10.3g %11.3g %11.1f %6.2f %8.2f"
                  % (c, T, ksv, kb, kb / ksv, n, lk))
    b, r2 = fit(pts)
    print("\n    log10(K_b/K_SV) vs n:  slope %+.2f per unit n, r2 %.2f"
          % (b, r2))
    print("    The divergence is n, and almost nothing else.\n")

    print("3.  What that means for the units\n")
    c, T, n, lk = worst[1]
    print("    log[(F0-F)/F] = log K_b + n log[Q]  ->  K_b is the intercept")
    print("    and carries units of (L/mol)^n.\n")
    print("    worst case: %s at %d K, n = %.2f, K_b/K_SV = %.0f"
          % (c, T, n, worst[0]))
    print("    its K_b is in (L/mol)^%.2f, and converts to log K %.2f --"
          % (n, lk))
    print("    the highest value in this compilation, against a ceiling of")
    print("    about 5.5 for anything measured directly.\n")

    print("4.  Does dG at least follow K_b?  dG = -RT ln K_b\n")
    bad = []
    for c, rows in T2.items():
        for T, ksv4, kb, n, dg in rows:
            calc = -R * T * math.log(kb) / 1000
            if abs(calc - dg) > 1.0:
                bad.append((c, T, dg, calc))
    if bad:
        for c, T, dg, calc in bad:
            print("    %-9s %d K: printed %+.2f, from K_b %+.2f  (off %.2f)"
                  % (c, T, dg, calc, abs(calc - dg)))
        print("\n    Everything else reproduces within 1 kJ/mol, so K_b and")
        print("    dG are the same quantity throughout -- the units problem")
        print("    is carried consistently, not introduced by a typo.")
    else:
        print("    all 18 reproduce within 1 kJ/mol")

    print("\n5.  Standing\n")
    print("    The main database already refuses a quenching n for")
    print("    conversion. The shareable sheet revives these rows on their")
    print("    own route because a whole-literature comparison was asked")
    print("    for. They stay. This script is the reason to read that route")
    print("    label rather than the number beside it.")


if __name__ == "__main__":
    main()
