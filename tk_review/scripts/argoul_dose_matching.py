"""
What is a comparable dose to Argoul 2026, and does the dose gap explain the
twofold low clearance?

scripts/argoul_vs_single_compound.py established that Argoul's clearances run
about 2.1x below the single-compound literature, and offered three candidate
causes. DOSE was called the most plausible: the cocktail was dosed far below
its comparators, and section 8's saturable-reabsorption account predicts
clearance rising with dose.

This script tests that directly, and the test does not go the way the earlier
write-up assumed. Two things have to be separated first:

  PER-COMPOUND DOSE   Argoul did not dose eleven PFAS at one level. Each is
                      tuned, roughly inversely to its persistence, so serum
                      stays measurable across a 5,254x clearance span. The
                      dose gap is therefore per compound, not a single number.

  COCKTAIL BURDEN     If the saturable step is a transporter shared across
                      PFAS, what matters is the TOTAL load, not each compound
                      alone. Summed, the cocktail is 4.88 mg/kg -- the same
                      order as a single-compound 1 mg/kg study. On that
                      reading the cocktail is not a low-dose experiment at all.

Run:  python3 scripts/argoul_dose_matching.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (compound, comparator, comparator dose mg/kg, comparator CL mL/kg/d, note)
COMPARATORS = [
    ("PFOA",  "Fujii 2015",              0.1284, 11.8,
     "single IV 0.31 umol/kg"),
    ("PFOA",  "Lou 2009",                 1.0,    5.99,
     "1 and 10 mg/kg oral, pooled fit; low end used"),
    ("PFHxS", "Sundstrom 2012, 1 mg/kg",  1.0,    2.68, "oral"),
    ("PFHxS", "Sundstrom 2012, 20 mg/kg", 20.0,   3.79, "oral"),
]


def main():
    p = os.path.join(HERE, "db", "primary_2026", "argoul2026_mouse_tk.csv")
    a = {r["chemical"]: r for r in csv.DictReader(open(p))}

    print(__doc__.split("Run:")[0].rstrip())
    print("=" * 78)
    print("PART 1.  Argoul's dose, compound by compound")
    print("=" * 78)
    print(f"\n  {'cmpd':8s}{'IV mg/kg':>10s}{'oral mg/kg':>12s}{'IV umol/kg':>12s}"
          f"   {'CL mL/kg/d':>11s}")
    ivs = []
    for r in a.values():
        mw = float(r["mw_g_mol"])
        iv = float(r["iv_dose_mgkg"]) if r["iv_dose_mgkg"] else None
        orl = float(r["oral_dose_mgkg"]) if r["oral_dose_mgkg"] else None
        cl = r["cl_mL_kg_day"] or "-"
        if iv:
            ivs.append(iv)
        print(f"  {r['chemical']:8s}{(f'{iv:.3f}' if iv else '-'):>10s}"
              f"{(f'{orl:.3f}' if orl else '-'):>12s}"
              f"{(f'{iv/mw*1000:.3f}' if iv else '-'):>12s}   {cl:>11s}")
    print(f"""
  Range 0.019 to 1.550 mg/kg, but that range is not a dose-response design --
  it tracks persistence. The three short-chain compounds (PFBA, PFHxA, PFBS)
  sit near 1.3-1.55 mg/kg; the long-chain ones (PFOS, PFDS, PFOA, PFDA) sit at
  0.019-0.080. So "the cocktail was dosed 10-100x low" is only true compound by
  compound, and the factor differs per comparison.

  Summed across the eleven, the IV burden is {sum(ivs):.2f} mg/kg.""")

    print()
    print("=" * 78)
    print("PART 2.  The dose gap against each comparator, and whether it predicts")
    print("=" * 78)
    print("""
  If dose drives the discrepancy, the clearance ratio should approach 1 as the
  dose gap shrinks. Ordering the four comparisons by dose gap tests that.
""")
    print(f"  {'cmpd':7s}{'comparator':27s}{'their dose':>11s}{'Argoul':>8s}"
          f"{'gap':>8s}{'CL ratio':>10s}")
    rowsout = []
    for c, name, dose, ccl, note in COMPARATORS:
        acl = float(a[c]["cl_mL_kg_day"])
        adose = float(a[c]["iv_dose_mgkg"])
        gap = dose / adose
        rowsout.append((gap, c, name, dose, adose, acl / ccl, note))
    for gap, c, name, dose, adose, ratio, note in sorted(rowsout):
        print(f"  {c:7s}{name:27s}{dose:>11.3f}{adose:>8.3f}"
              f"{gap:>7.0f}x{ratio:>9.2f}x")
    print("""
  It does not. The comparison with almost NO dose gap -- Fujii at 0.128 mg/kg
  against Argoul's 0.080, a factor of 1.6 -- has the WORST agreement of the
  four (0.38x). The comparison with the largest gap among the PFOA pair, Lou at
  1 mg/kg and above, has the BEST (0.75x). The ordering is not monotonic in
  dose, and the one nearly dose-matched comparison is the one that disagrees
  most.

  That is evidence against dose as the main cause, and it was missed earlier
  because the doses were quoted as a range rather than matched per compound.""")

    print()
    print("=" * 78)
    print("PART 3.  Out-of-sample test on the one real dose series")
    print("=" * 78)
    print("""
  Sundstrom dosed PFHxS at 1 and 20 mg/kg in the same experiment, so its own
  dose-clearance slope can be fitted WITHOUT Argoul and extrapolated down to
  Argoul's dose. That is a genuine prediction rather than a fit through all
  three points, which is what the earlier write-up reported.
""")
    # First: how much did Argoul's own point define the slope it was checked
    # against? With n=3 and Argoul at one end of the x range, a lot.
    pts = [(0.094, 1.3, "Argoul"), (1.0, 2.68, "Sundstrom 1"),
           (20.0, 3.79, "Sundstrom 20")]
    xs = [math.log10(d) for d, _, _ in pts]
    ys = [math.log10(c) for _, c, _ in pts]
    mx, my = sum(xs) / 3, sum(ys) / 3
    sxx = sum((x - mx) ** 2 for x in xs)
    three = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    print(f"  The slope the earlier write-up quoted, {three:+.3f}, was fitted through "
          f"ALL\n  THREE points, Argoul's included. Leverage on that fit:\n")
    for (d, _, nm), x in zip(pts, xs):
        h = 1 / 3 + (x - mx) ** 2 / sxx
        print(f"    {nm:14s} dose {d:>6.3f} mg/kg    leverage {h:.2f}")
    print("""
  Argoul's point carries 79% leverage, so "the cocktail sits at the bottom of
  that trend rather than off it" was close to guaranteed by construction: the
  point helped draw the line it was then checked against. The honest test drops
  it.
""")

    cl1, cl20 = 2.68, 3.79
    slope = math.log(cl20 / cl1) / math.log(20.0 / 1.0)
    adose, acl = float(a["PFHxS"]["iv_dose_mgkg"]), float(a["PFHxS"]["cl_mL_kg_day"])
    pred = cl1 * (adose / 1.0) ** slope
    print(f"  Sundstrom's own slope, from its two doses:   {slope:+.3f} on log-log")
    print(f"  extrapolated to Argoul's {adose:.3f} mg/kg:          "
          f"{pred:.2f} mL/kg/d predicted")
    print(f"  Argoul observed:                             {acl:.2f} mL/kg/d")
    print(f"  residual after dose correction:              {acl/pred:.2f}x")
    print(f"""
  So the dose trend closes part of the PFHxS gap -- 0.49x becomes {acl/pred:.2f}x --
  but about {1/(acl/pred):.1f}x survives it. Dose is a contributor, not the explanation.

  Caveat: this extrapolates a slope fitted over 1-20 mg/kg down to 0.094, a
  factor of ten below the lowest observation, and the two points carry no
  reported uncertainty, so the slope itself is a point estimate from n=2.""")

    print()
    print("=" * 78)
    print("PART 4.  So what IS a comparable dose?")
    print("=" * 78)
    print(f"""
  Three defensible answers, and they disagree, which is the useful part:

   PER COMPOUND     Match Argoul's own per-compound dose. For PFOA that is
                    0.080 mg/kg IV, and Fujii at 0.128 mg/kg already does it.
                    For PFHxS it is 0.094 mg/kg, and nothing in the female-mouse
                    literature comes close -- the nearest is 1 mg/kg, 10.6x up.

   TOTAL BURDEN     If the saturable step is shared across PFAS, the comparable
                    single-compound dose is the summed {sum(ivs):.2f} mg/kg, which
                    sits between Sundstrom's 1 and 20 mg/kg arms. On this
                    reading the cocktail is a MID-dose experiment, and the
                    "10-100x below its comparators" framing is wrong.

   MOLAR            Dose per kg understates the short-chain compounds, whose
                    molar dose is 20-80x that of the long-chain ones in the same
                    animal. A molar-matched comparator would differ per compound
                    again.

  The experiment that would settle it is not a new comparator: it is Argoul's
  own design run at a second total dose. Nothing in the current literature
  separates per-compound dose from cocktail burden, because no one has dosed
  the same cocktail at two levels.""")


if __name__ == "__main__":
    main()
