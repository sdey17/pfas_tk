#!/usr/bin/env python3
"""Chen 2025's PFNA sits below his PFOA. Is that real, and does it matter?

Section 11a of ALBUMIN_BINDING.md left this as the largest residual between
Chen 2025 and the other two albumin datasets: his K_D puts PFNA (3.75 uM)
WEAKER than PFOA (2.57 uM), where Fischer 2024 and Qin 2024 both put PFNA at
or above PFOA.

The transcription has been checked against the SI page image: Table S5 really
does print "3.75 +/- 0.163" for log K_MW and "3.75 +/- 0.77" for K_D on the
PFNA row. The repeated 3.75 is a coincidence of two different quantities --
K_MW 3.75 is exactly where the +0.36/CF2 membrane trend puts it.

Chen reports standard errors, so the question can be asked properly.

Run:  python3 tk_review/scripts/chen_pfna_inversion.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
M_ALB = 66463.0
VOL = math.log10(1.0 / 0.733)
LINE = "=" * 78
LN10 = math.log(10.0)

# Chen 2025 Table S5, K_D uM +/- SE, in chain order
CHEN_PFCA = [("PFBA", 3, 84.9, 15.95), ("PFPeA", 4, 16.27, 3.58),
             ("PFHxA", 5, 13.98, 2.79), ("PFHpA", 6, 6.2, 1.13),
             ("PFOA", 7, 2.57, 0.45), ("PFNA", 8, 3.75, 0.77),
             ("PFDA", 9, 3.64, 0.85), ("PFUnDA", 10, 4.73, 1.08),
             ("PFDoDA", 11, 15.76, 3.34), ("PFTrDA", 12, 31.44, 5.09),
             ("PFTeDA", 13, 13.09, 3.01), ("PFHxDA", 15, 10.45, 1.71)]
CHEN_PFSA = [("PFBS", 4, 15.6, 2.68), ("PFPeS", 5, 17.15, 2.58),
             ("PFHxS", 6, 1.31, 0.25), ("PFHpS", 7, 1.68, 0.36),
             ("PFOS", 8, 0.94, 0.23), ("PFNS", 9, 6.56, 1.36),
             ("PFDS", 10, 1.88, 0.4), ("PFDoS", 12, 1.29, 0.21)]


def logD(kd_uM, n=1):
    return math.log10(1000.0 * n * (1.0 / (kd_uM * 1e-6)) / M_ALB)


def logD_se(kd, se):
    """Delta method: SE of log10(1/Kd) is SE_Kd / (Kd ln10)."""
    return se / (kd * LN10)


def fischer():
    out = {}
    for r in csv.DictReader(open(DB)):
        if (r["study"], r["year"]) == ("Fischer et al.", "2024") and \
                "human serum albumin" in r["protein"] and \
                r["ka_units"].startswith("log D"):
            out[r["chemical"]] = float(r["ka"])
    return out


def qin_nonspec():
    out = {}
    for r in csv.DictReader(open(DB)):
        if (r["study"], r["year"]) == ("Qin et al.", "2024") and \
                "bovine" in r["protein"] and "nonspecific" in r["ka_units"]:
            out[r["chemical"]] = float(r["ka"]) - VOL
    return out


def z_test(a, sa, b, sb):
    d = a - b
    s = math.sqrt(sa ** 2 + sb ** 2)
    z = d / s
    # two-sided normal p, via erfc
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return d, s, z, p


def main():
    print(__doc__.split("Run:")[0].rstrip())
    F, Q = fischer(), qin_nonspec()

    # ------------------------------------------------------------------ 1
    print(LINE)
    print("1.  IS THE INVERSION SIGNIFICANT IN CHEN'S OWN DATA?")
    print(LINE)
    print("""
  Comparing PFOA against each of its neighbours, on log D with the standard
  error propagated from Chen's reported SE on K_D.
""")
    d = {c: (logD(k), logD_se(k, s), n) for c, n, k, s in CHEN_PFCA}
    print(f"  {'pair':22s}{'diff (log)':>12s}{'SE':>8s}{'z':>7s}{'p':>9s}  verdict")
    for other in ("PFHpA", "PFNA", "PFDA", "PFUnDA"):
        a, sa, _ = d["PFOA"]
        b, sb, _ = d[other]
        diff, se, z, p = z_test(a, sa, b, sb)
        v = "significant" if p < 0.05 else "NOT significant"
        print(f"  {'PFOA vs ' + other:22s}{diff:>12.3f}{se:>8.3f}{z:>7.2f}"
              f"{p:>9.3f}  {v}")
    print("""
  PFOA, PFNA and PFDA are statistically indistinguishable in Chen's data.
  Calling it an "inversion" overstates what the numbers support: PFOA is
  nominally the tightest, and PFNA/PFDA/PFUnDA form a plateau just above it
  that the measurement cannot separate.  That alone retires most of the
  puzzle.
""")

    # ------------------------------------------------------------------ 2
    print(LINE)
    print("2.  THE REAL DIFFERENCE IS NOT PFNA, IT IS DYNAMIC RANGE")
    print(LINE)
    print(f"\n  {'cmpd':9s}{'nCF':>4s}{'Chen logD':>11s}{'+/-':>7s}"
          f"{'Fischer':>9s}{'Fischer - Chen':>16s}")
    gaps, fvals = [], []
    for c, n, k, s in CHEN_PFCA:
        f = F.get(c)
        cd, cs = logD(k), logD_se(k, s)
        if f is None:
            print(f"  {c:9s}{n:>4d}{cd:>11.2f}{cs:>7.2f}{'-':>9s}{'-':>16s}")
            continue
        print(f"  {c:9s}{n:>4d}{cd:>11.2f}{cs:>7.2f}{f:>9.2f}{f - cd:>16.2f}")
        gaps.append(f - cd)
        fvals.append(f)
    crange = max(logD(k) for _c, _n, k, _s in CHEN_PFCA) - \
        min(logD(k) for _c, _n, k, _s in CHEN_PFCA)
    frange = max(F[c] for c, _n, _k, _s in CHEN_PFCA if c in F) - \
        min(F[c] for c, _n, _k, _s in CHEN_PFCA if c in F)
    cpeak = max(CHEN_PFCA, key=lambda t: logD(t[2]))
    fpeak = max(((c, F[c]) for c, _n, _k, _s in CHEN_PFCA if c in F),
                key=lambda t: t[1])
    print(f"""
  Chen's PFCA series spans {crange:.2f} log units; Fischer's spans {frange:.2f} over the same
  compounds.  Chen's optimum sits at {cpeak[0]} (nCF {cpeak[1]}), Fischer's at {fpeak[0]}.

  So the PFNA cell is not an isolated oddity.  Chen's whole chain-length
  profile is COMPRESSED, by about a factor of two in range, and the
  compression pushes the apparent optimum to shorter chains.  PFNA sitting
  below PFOA is one symptom of that, not a separate problem.
""")

    # ------------------------------------------------------------------ 3
    print(LINE)
    print("3.  WHAT COMPRESSES IT: three candidate explanations, compared")
    print(LINE)
    print("""
  A compression has a predictor.  Three candidates, each regressed against the
  shortfall (Fischer - Chen), using Fischer as the stand-in for true affinity:

    AFFINITY        a tight binder's specific site is substantially occupied
                    even at Chen's lowest ligand concentration, so the
                    single-site fit is dragged toward the nonspecific value.
                    The tighter the compound, the more complete that, so the
                    gap should scale with true affinity.
    CHAIN LENGTH    residual sorption of long-chain PFAS to the apparatus,
                    which Chen explicitly fought (they abandoned the RED
                    system for exactly this and moved to Slide-A-Lyzer).
                    Makes a long compound look weaker regardless of affinity.
    IDENTIFIABILITY a single-site K_D is best determined when it is near the
                    protein concentration.  Chen held HSA at 1 uM, so
                    |log10(K_D / 1 uM)| measures how far outside that window
                    each compound sits.
""")
    sets = {"all 16": [c for c, _n, _k, _s in CHEN_PFCA + CHEN_PFSA if c in F],
            "PFCA only": [c for c, _n, _k, _s in CHEN_PFCA if c in F]}
    kd = {c: k for c, _n, k, _s in CHEN_PFCA + CHEN_PFSA}
    ncf = {c: n for c, n, _k, _s in CHEN_PFCA + CHEN_PFSA}

    def fit(xs, ys):
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        sxx = sum((x - mx) ** 2 for x in xs)
        syy = sum((y - my) ** 2 for y in ys)
        return sxy / math.sqrt(sxx * syy), sxy / sxx

    print(f"  {'predictor':34s}" +
          "".join(f"{k + ' r2':>14s}" for k in sets))
    res = {}
    for label, getx in (
            ("AFFINITY (Fischer log D)", lambda c: F[c]),
            ("CHAIN LENGTH (nCF)", lambda c: ncf[c]),
            ("IDENTIFIABILITY |log Kd/[P]|", lambda c: abs(math.log10(kd[c] / 1.0))),
    ):
        line = f"  {label:34s}"
        for k, names in sets.items():
            g = [F[c] - logD(kd[c]) for c in names]
            r, sl = fit([getx(c) for c in names], g)
            res[(label, k)] = (r, sl)
            line += f"{r*r:>14.2f}"
        print(line)
    rA, slA = res[("AFFINITY (Fischer log D)", "PFCA only")]
    print(f"""
  AFFINITY wins, and not narrowly: r2 {res[('AFFINITY (Fischer log D)','PFCA only')][0]**2:.2f} against {res[('CHAIN LENGTH (nCF)','PFCA only')][0]**2:.2f} for chain length and
  {res[('IDENTIFIABILITY |log Kd/[P]|','PFCA only')][0]**2:.2f} for identifiability on the carboxylates.  The identifiability
  idea is simply wrong -- it explains nothing, and it is reported here because
  it was the first guess.  Chain length is a real but weaker correlate, which
  is expected since affinity rises with chain length over most of the series;
  the two cannot be fully separated in this dataset, but affinity fits better
  at both ends while chain length does not.

  The fitted slope is {slA:.2f} on the carboxylates.  A slope of 1 would mean
  Chen's values are flat regardless of true affinity -- complete compression.
  {slA:.2f} means roughly four fifths of any true affinity difference is lost.
""")
    # pivot prediction
    names = sets["PFCA only"]
    g = [F[c] - logD(kd[c]) for c in names]
    mF = sum(F[c] for c in names) / len(names)
    mg = sum(g) / len(g)
    pivot = mF - mg / slA
    print(f"""  A compression has a PIVOT -- the affinity at which Chen and Fischer agree.
  Fitted on all {len(names)} carboxylates it sits at Fischer log D = {pivot:.2f}, so below
  that Chen should read TIGHTER than Fischer and the gap should change sign.

  To make that a real test rather than a restatement of the fit, each compound
  is predicted from a regression REFITTED WITH THAT COMPOUND LEFT OUT:
""")
    print(f"    {'cmpd':8s}{'Fischer':>9s}{'predicted':>11s}{'observed':>10s}{'miss':>7s}")
    miss = []
    for c in names:
        rest = [x for x in names if x != c]
        gr = [F[x] - logD(kd[x]) for x in rest]
        _r, sl = fit([F[x] for x in rest], gr)
        mF2 = sum(F[x] for x in rest) / len(rest)
        mg2 = sum(gr) / len(gr)
        piv = mF2 - mg2 / sl
        pred = sl * (F[c] - piv)
        obs = F[c] - logD(kd[c])
        miss.append(abs(pred - obs))
        print(f"    {c:8s}{F[c]:>9.2f}{pred:>11.2f}{obs:>10.2f}{abs(pred-obs):>7.2f}")
    print(f"""
  Mean leave-one-out miss {sum(miss)/len(miss):.2f} log units, worst {max(miss):.2f}.  PFPeA, the one
  compound weak enough to sit below the pivot, is predicted at the correct
  SIGN and within {miss[names.index('PFPeA')]:.2f} log units when it is excluded from the fit.

  So this is a genuine compression about a pivot, not a uniform offset and not
  a handful of bad cells: Chen reads weak binders too tight and tight binders
  too weak, by a factor that one line predicts across three log units.

  One caveat on the mechanism.  The regression identifies affinity as the
  predictor; it does not prove the saturation story behind it.  Container
  sorption would also scale with hydrophobicity and so partly with affinity.
  What the data do support is that the compression is systematic, predictable
  from affinity, and a property of the assay rather than of the chemistry.
""")

    # ------------------------------------------------------------------ 4
    print(LINE)
    print("4.  WHAT IT MEANS FOR USING CHEN")
    print(LINE)
    qc = [(c, logD(k), Q[c]) for c, _n, k, _s in CHEN_PFCA if c in Q]
    print(f"\n  {'cmpd':9s}{'Chen':>8s}{'Qin nonspec':>13s}{'diff':>8s}")
    for c, cd, q in qc:
        print(f"  {c:9s}{cd:>8.2f}{q:>13.2f}{cd - q:>8.2f}")
    print("""
  Chen still matches Qin's non-specific term where both exist and the chain is
  short -- PFHpA +0.05, PFOA +0.02 -- which is section 11a's result and is
  unaffected.  PFNA is where that agreement breaks (-0.87), and it breaks in
  the direction compression predicts.

  Practical reading:

    * Chen's RANKING is sound at the short end and compressed at the long end.
      Use it for compounds up to about nCF 8.
    * Do NOT read a chain-length OPTIMUM off Chen alone.  Fischer puts it near
      nCF 10, Chen near nCF 7, and the difference is assay compression rather
      than chemistry.  All three datasets agree there IS an optimum and that
      binding falls away after it; they disagree on where inside C7-C10.
    * The PFNA cell itself is within noise of PFOA and should not be treated
      as a finding in either direction.
""")


if __name__ == "__main__":
    main()
