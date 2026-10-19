#!/usr/bin/env python3
"""Why Chen 2025 and Fischer 2024 disagree about albumin, and why they don't.

Section 14 of report/ALBUMIN_BINDING.md listed this as an open question: the
two best albumin datasets in the file, both measured at defensible
ligand:protein ratios, differ by about 5x for PFOA.  Converting Chen's K_D of
2.57 uM at one site gives log D 3.77; Fischer measured 4.48.

Two explanations compete:

  REGIME    Fischer works at nu <= 0.004, where SPECIFIC binding dominates.
            Chen works at nu = 0.06-2, which straddles the transition to
            NONSPECIFIC binding.  They would then be measuring two different
            real constants, not disagreeing about one.

  STOICHIOMETRY  Chen's conversion assumed n = 1.  Since log D = log Ka +
            log n - 1.823, an n of about 5 would close the gap by itself.

The two are distinguishable, because Qin et al. 2024 measured the specific
AND nonspecific terms separately for the same compounds.  If REGIME is right,
Fischer lands on Qin's specific term and Chen on his nonspecific term.
STOICHIOMETRY would move Chen up onto Fischer, which means moving him OFF
Qin's nonspecific term.  Only one of those can be true, and section 1 shows
which.

Run:  python3 tk_review/scripts/chen_fischer_reconciliation.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
M_ALB = 66463.0
VOL_SHIFT = math.log10(1.0 / 0.733)      # mass -> volume basis, +0.135
LINE = "=" * 78

rows = list(csv.DictReader(open(DB)))


def load():
    chen, fis, qspec, qnon = {}, {}, {}, {}
    for r in rows:
        sy = (r["study"], r["year"])
        if sy == ("Chen et al.", "2025") and r["kd"]:
            ka = 1.0 / (float(r["kd"]) * 1e-6)
            chen[r["chemical"]] = math.log10(1000.0 * ka / M_ALB)   # n = 1, mass
        elif sy == ("Fischer et al.", "2024") and \
                "human serum albumin" in r["protein"] and \
                r["ka_units"].startswith("log D"):
            fis[r["chemical"]] = float(r["ka"])                      # mass
        elif sy == ("Qin et al.", "2024") and "bovine" in r["protein"]:
            v = float(r["ka"]) - VOL_SHIFT                           # -> mass
            if "nonspecific" in r["ka_units"]:
                qnon[r["chemical"]] = v
            elif "specific" in r["ka_units"]:
                qspec[r["chemical"]] = v
    return chen, fis, qspec, qnon


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy), sxy / sxx


def main():
    print(__doc__.split("Run:")[0].rstrip())
    chen, fis, qspec, qnon = load()

    print(LINE)
    print("1.  ALL FOUR ON ONE BASIS (mass, L_water per kg_protein)")
    print(LINE)
    print("""
  Qin reports on a volume basis and is shifted by -0.135 to match.  Chen is
  converted from K_D at n = 1.  Fischer is measured directly.
""")
    order = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFHxS", "PFOS"]
    print(f"  {'cmpd':8s}{'Fischer':>9s}{'Qin spec':>10s}{'F - Qs':>9s}"
          f"{'Chen':>8s}{'Qin nonsp':>11s}{'C - Qn':>9s}")
    d_fs, d_cn = [], []
    for c in order:
        f, s, ch, n = fis.get(c), qspec.get(c), chen.get(c), qnon.get(c)
        def fmt(x, w):
            return f"{x:>{w}.2f}" if x is not None else f"{'-':>{w}}"
        a = f"{f - s:+.2f}" if f is not None and s is not None else "-"
        b = f"{ch - n:+.2f}" if ch is not None and n is not None else "-"
        print(f"  {c:8s}{fmt(f,9)}{fmt(s,10)}{a:>9s}{fmt(ch,8)}{fmt(n,11)}{b:>9s}")
        if f is not None and s is not None:
            d_fs.append(abs(f - s))
        if ch is not None and n is not None:
            d_cn.append(abs(ch - n))
    print(f"""
  Fischer against Qin's SPECIFIC term:    mean |gap| {sum(d_fs)/len(d_fs):.2f} log units
  Chen    against Qin's NONSPECIFIC term: mean |gap| {sum(d_cn)/len(d_cn):.2f} log units

  For PFOA the two matches are 0.03 and 0.02 log units -- 7% and 5%.  Three
  laboratories, three methods, two regimes, and each pair lands on top of the
  other.
""")

    print(LINE)
    print("2.  THE DISCRIMINATING TEST")
    print(LINE)
    print("""
  A weaker, secondary check: does the Fischer-Chen gap track Qin's own
  specific-to-nonspecific separation compound by compound?  It should if the
  regime explanation holds, though with seven compounds this is suggestive at
  best.
""")
    xs, ys, labels = [], [], []
    print(f"  {'cmpd':8s}{'Fischer - Chen':>16s}{'Qin spec - nonsp':>19s}")
    for c in order:
        f, ch, s, n = fis.get(c), chen.get(c), qspec.get(c), qnon.get(c)
        if None in (f, ch, s, n):
            continue
        xs.append(s - n)
        ys.append(f - ch)
        labels.append(c)
        print(f"  {c:8s}{f - ch:>16.2f}{s - n:>19.2f}")
    r, slope = pearson(xs, ys)
    print(f"""
  r = {r:.2f} over {len(xs)} compounds -- which settles nothing on its own.  Five
  points cannot carry a correlation test, and this one is dragged by the
  sulfonates.  Splitting them out:
""")
    for group, names in (("carboxylates", ["PFHxA", "PFHpA", "PFOA"]),
                         ("sulfonates", ["PFHxS", "PFOS"])):
        gx = [xs[labels.index(c)] for c in names if c in labels]
        gy = [ys[labels.index(c)] for c in names if c in labels]
        print(f"    {group:14s} Fischer-Chen {sum(gy)/len(gy):+.2f}  "
              f"Qin spec-nonsp {sum(gx)/len(gx):+.2f}  "
              f"ratio {10**(sum(gy)/len(gy)):.1f}x vs {10**(sum(gx)/len(gx)):.1f}x")
    print("""
  For the carboxylates the two separations agree (3.4x against 4.4x).  For the
  sulfonates Qin's own separation is much larger than the Fischer-Chen gap, so
  something else is also moving there -- most likely that Chen's single-site
  Langmuir over nu = 0.06-2 cannot cleanly resolve a sulfonate's very tight
  specific site.

  THE ARGUMENT THAT ACTUALLY DISCRIMINATES is simpler than the correlation,
  and it is in section 1.  Chen currently sits on Qin's NONSPECIFIC term to
  0.02 log units for PFOA and 0.05 for PFHpA.  Raising Chen's n from 1 to the
  ~5 the stoichiometry explanation needs would lift him onto Fischer -- and in
  doing so would BREAK an agreement that presently holds to 5%.  You cannot
  have both.  Since Qin measured his two terms on the same instrument in the
  same paper, the agreement Chen already has is the one worth keeping, and the
  stoichiometry explanation is the one to drop.

  Put the other way: n = 1 is what makes all three datasets mutually
  consistent.  That is also what the backwards conversion in REPORT.md 6.1b
  found independently (n_implied = 1.01 for PFOA, 1.12 for PFHpA), from a
  different direction.
""")

    print(LINE)
    print("3.  WHAT REMAINS GENUINELY ODD")
    print(LINE)
    pfoa, pfna = chen.get("PFOA"), chen.get("PFNA")
    print(f"""
  Chen's PFNA breaks its own series.  His K_D has PFNA (3.75 uM) WEAKER than
  PFOA (2.57 uM), i.e. log D {pfna:.2f} against {pfoa:.2f}.  Both Fischer
  ({fis['PFNA']:.2f} vs {fis['PFOA']:.2f}) and Qin's nonspecific term
  ({qnon['PFNA']:.2f} vs {qnon['PFOA']:.2f}) put PFNA at or above PFOA.  Chen
  is the only one of the three with that inversion, and it is the single
  largest residual in section 1 ({chen['PFNA'] - qnon['PFNA']:+.2f} log units).

  Note this does not touch the chain-length story: all three agree that
  binding peaks somewhere around C7-C9 and falls away after.  They disagree
  about where inside that plateau the maximum sits.
""")

    print(LINE)
    print("4.  THE ANSWER, AND WHAT TO USE")
    print(LINE)
    print("""
  Chen and Fischer are NOT in conflict.  They measure the two different
  binding regimes of the same protein, and each agrees with the corresponding
  term of the one study that measured both.  The ~5x "discrepancy" is the
  specific-to-nonspecific separation, which is a real property of albumin and
  not an error in either dataset.

  Which to use depends entirely on the question:

    a PBTK or plasma-binding calculation  -> the SPECIFIC constant
       (Fischer, or Qin's specific term).  The body runs at nu < 0.01.

    an in vitro mass-balance model of a dosed assay -> the NONSPECIFIC
       constant (Chen, or Qin's nonspecific term).  The assay runs at nu > 1.

  Using the wrong one costs about 5x for carboxylates and more for
  sulfonates, in a direction that depends on which way you cross.  And it
  retires the open question in ALBUMIN_BINDING.md section 14, item 2: the
  stoichiometry explanation is not needed, because n = 1 plus the regime
  accounts for the gap.
""")


if __name__ == "__main__":
    main()
