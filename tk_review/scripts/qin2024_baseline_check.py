#!/usr/bin/env python3
"""Qin et al. 2024 baseline-toxicity paper, checked against this project.

Qin W, Henneberger L, Gluge J, Konig M, Escher BI.  "Baseline Toxicity Model
to Identify the Specific and Nonspecific Effects of PFAS in Cell-Based
Bioassays."  Environ Sci Technol 2024;58(13):5727-5738.
doi:10.1021/acs.est.3c09950 (open access, PMC10993398, SI obtained).

It is the direct successor to Qin et al. 2023, the model this project
implements, by the same group and on the same assay family.  Four things in
it bear on the review.

Run:  python3 tk_review/scripts/qin2024_baseline_check.py
"""
import math
import statistics as st

import numpy as np

# ---- Table 1, all logs.  D_cell/w and D_SP/w measured near IC10,nom in
# HEK293H; D_lip/w from Ebert 2020 (e) or Droge (d) or predicted (f).
T1 = {   # name: (D_cell/w, D_SP/w, D_lip/w, D_BSA/w, nCF, class, lip_src)
    "PFBA":    (2.40, 2.11, 1.00, 1.94,  3, "PFCA", "Droge"),
    "PFHxA":   (3.45, 2.45, 2.32, 2.70,  5, "PFCA", "Ebert"),
    "PFHpA":   (3.69, 2.65, 2.91, 3.40,  6, "PFCA", "Ebert"),
    "PFOA":    (3.94, 3.14, 3.52, 3.98,  7, "PFCA", "Ebert"),
    "PFNA":    (3.74, 3.39, 4.25, 4.39,  8, "PFCA", "Ebert"),
    "PFUnA":   (3.78, 3.97, 4.54, 4.75, 10, "PFCA", "Ebert"),
    "HFPO-DA": (3.09, 2.77, 2.41, 2.28, None, "ether", "Ebert"),
    "PFHxS":   (3.56, 2.95, 4.13, 3.52,  6, "PFSA", "Ebert"),
    "PFOS":    (4.32, 3.94, 4.89, 4.52,  8, "PFSA", "Ebert"),
    "6:2 FTSA":(3.90, 3.39, 3.87, 3.71, None, "FTSA", "predicted"),
    "PFOSA":   (3.81, 3.52, 4.94, 4.33, None, "FOSA", "predicted"),
}
# SI Table S8, measured cell composition, mL per L of cell
CELLS = {"HEK293H (PPARg)": (33.9, 5.3), "MCF7 (AREc32)": (18.7, 4.1),
         "H4IIe (AhR)": (25.5, 4.1), "SH-SY5Y (neuro)": (34.0, 4.9),
         "generic": (30.0, 5.0)}


def main():
    print(__doc__.split("Run:")[0].rstrip())

    print("\n" + "=" * 74)
    print("1.  THE MEASURED LIPID SERIES, PER COMPOUND")
    print("=" * 74)
    print("""
  The review's section 6.7 used Qin 2023's Figure A1 REGRESSION (+0.476) to
  represent the Ebert liposome data.  This paper tabulates the per-compound
  values, so the regression can be replaced by the points it was fit to.""")
    pfca = {v[4]: (k, v[2]) for k, v in T1.items() if v[5] == "PFCA"}
    pfsa = {v[4]: (k, v[2]) for k, v in T1.items() if v[5] == "PFSA"}
    for lab, d in (("carboxylates", pfca), ("sulfonates", pfsa)):
        xs = sorted(d)
        print(f"\n  {lab}: " + "  ".join(f"{d[n][0]} {d[n][1]:.2f}" for n in xs))
    for lab, keys in (("all six (incl. Droge PFBA)", sorted(pfca)),
                      ("Ebert only, nCF 5-10", [k for k in sorted(pfca) if k >= 5])):
        s, _i = np.polyfit(keys, [pfca[k][1] for k in keys], 1)
        r = np.corrcoef(keys, [pfca[k][1] for k in keys])[0, 1] ** 2
        print(f"\n  carboxylate slope, {lab:28s} {s:+.3f}  R2 {r:.3f}")
    print(f"  sulfonate slope, two points                   "
          f"{(pfsa[8][1] - pfsa[6][1]) / 2:+.3f}")
    print("""
  Against the review's other sources:
      Chen 2025 SSLM, 12 carboxylates           +0.361
      Qin 2023 Fig A1 regression (used so far)  +0.476
      IV-MBM EQP v2.0 spLFER                    +0.78

  So the per-compound data land at +0.46 (Ebert range) to +0.53 (with the
  Droge PFBA point), bracketing the regression this project has been using
  and leaving the spLFER comparison intact.  The conclusion of section 6.7
  does not change; its weakest link, a regression standing in for data, is
  now replaced by the data.""")

    print("\n" + "=" * 74)
    print("2.  CHEN AGAINST EBERT, THE SAME COMPOUNDS BY DIFFERENT METHODS")
    print("=" * 74)
    CHEN = {"PFBA": 1.63, "PFHxA": 2.42, "PFHpA": 2.85, "PFOA": 3.28,
            "PFNA": 3.75, "PFHxS": 3.65, "PFOS": 4.50}
    print(f"\n  {'cmpd':9s}{'Chen SSLM':>11s}{'Ebert lip':>11s}{'diff':>8s}"
          f"   (pH 7.0 vs 7.4)")
    diffs = []
    for k, c in CHEN.items():
        e = T1[k][2]
        diffs.append(e - c)
        print(f"  {k:9s}{c:>11.2f}{e:>11.2f}{e - c:>+8.2f}")
    print(f"""
  mean offset {st.mean(diffs):+.2f} log.  They agree closely at nCF 5-7 and
  diverge at both ends, which is the compression of Chen's chain-length
  profile that ALBUMIN_BINDING section 11b identified on the albumin axis --
  here it shows on the membrane axis, from an independent method.""")

    print("\n" + "=" * 74)
    print("3.  DOES EQ 5 REPRODUCE THE PAPER'S OWN MEASURED D_cell/w?")
    print("=" * 74)
    fp, fl = CELLS["HEK293H (PPARg)"]
    f_p, f_l = fp / (fp + fl), fl / (fp + fl)
    print(f"""
  Eq 5 is  D_cell/w = D_SP/w x Vprot/Vprot+lip + D_lip/w x Vlip/Vprot+lip,
  a weighted mean with weights summing to one.  HEK293H is {fp} and {fl} mL/L
  (SI Table S8), so the weights are {f_p:.3f} and {f_l:.3f}.

  A weighted mean CANNOT fall outside the range of the two values it
  averages.  That makes the check composition-independent: wherever the
  measured D_cell/w sits outside [D_SP/w, D_lip/w], no cell composition can
  reconcile eq 5 with Table 1.
""")
    print(f"  {'cmpd':9s}{'D_SP/w':>8s}{'D_lip/w':>9s}{'meas':>8s}"
          f"{'eq 5':>8s}{'pred/meas':>11s}   within range?")
    out = 0
    for k, (dc, dsp, dl, _db, _n, _c, _s) in T1.items():
        pred = math.log10(10 ** dsp * f_p + 10 ** dl * f_l)
        lo, hi = min(dsp, dl), max(dsp, dl)
        ok = lo - 1e-9 <= dc <= hi + 1e-9
        out += 0 if ok else 1
        flag = "yes" if ok else ("ABOVE both" if dc > hi else "BELOW both")
        print(f"  {k:9s}{dsp:>8.2f}{dl:>9.2f}{dc:>8.2f}{pred:>8.2f}"
              f"{10 ** (pred - dc):>10.2f}x   {flag}")
    print(f"""
  {out} of {len(T1)} compounds sit outside the achievable range.  Six are above
  both phases (PFBA, PFHxA, PFHpA, PFOA, HFPO-DA, 6:2 FTSA) and one below
  both (PFUnA).  For PFOA the measured cell-water ratio is 3.94 while the
  more sorptive of its two phases is 3.52, so the cell takes up 2.6x more
  per litre of protein+lipid than either surrogate can deliver.

  The paper reports these as agreeing well (Figure 2c).  On a log-log scatter
  spanning three decades a factor of 2-10 does sit close to the 1:1 line, so
  the two statements are not in conflict -- but the arithmetic above is not
  visible on that plot.  What it suggests is that the two surrogates, chicken
  breast protein for structural protein and liposomes for membrane lipid,
  under-represent real cellular sorption for the shorter carboxylates.  This
  cannot be settled further from the published record.""")

    print("\n" + "=" * 74)
    print("4.  CELL COMPOSITION IS NOT A SETTLED NUMBER")
    print("=" * 74)
    print(f"\n  {'source':34s}{'line':16s}{'protein':>9s}{'lipid':>8s}{'P:L':>7s}")
    rows = [("Qin 2024 SI Table S8", "HEK293H", 3.39, 0.53),
            ("Qin 2023 Table A5", "HEK293H", 8.14, 0.432),
            ("Qin 2024 SI Table S8", "MCF7 (AREc32)", 1.87, 0.41),
            ("Armitage 2021 SI Table S5", "AREc2 (=AREc32)", 5.14, 0.46),
            ("Fischer 2017 Table 2", "HepG2", 9.46, 3.19)]
    for src, line, p, l in rows:
        print(f"  {src:34s}{line:16s}{p:>8.2f}%{l:>7.2f}%{p / l:>7.1f}")
    print("""
  The same group reports HEK293H protein at 8.14% in 2023 and 3.39% in 2024,
  a factor of 2.4 on the dominant cell sink, and their MCF7 protein is 2.7x
  below the value Armitage's model uses for the same line.  Two independent
  confirmations fall out of this table for the review:

    - AREc32 IS MCF7-derived.  Qin 2024 Table 2 labels it "AREc32 (MCF7)"
      outright, which settles section 6.7.2's correction from the Escher
      group themselves rather than from a cell-line registry.
    - The lipid fractions agree far better across labs (0.41-0.53% for the
      non-hepatic lines) than the protein fractions do.  The quantity the
      review found to carry the Armitage cell term is the better measured
      one; the quantity it found to be undocumented is the volatile one.""")


if __name__ == "__main__":
    main()
