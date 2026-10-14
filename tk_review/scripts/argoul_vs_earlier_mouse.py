"""
How do Argoul 2026's mouse clearances and volumes compare with earlier values?

Argoul 2026 is the only study that measures CL and Vd for eleven PFAS in one
experiment, so it carries a lot of weight in this review (sections 3 and 10).
That makes it worth asking whether its numbers agree with what was already
published for the mouse -- and if not, whether the disagreement has a
structure.

Three design differences matter before any number is compared:

  SEX     Argoul dosed FEMALE mice. Several earlier mouse values are male, or
          pooled. The mouse sex difference is small for PFOA (1.2-1.7x) but
          not zero.
  DOSE    Argoul dosed 0.019-1.55 mg/kg, one to two orders of magnitude below
          the earlier studies. If elimination saturates, that alone predicts a
          lower clearance and a smaller volume.
  DESIGN  Argoul gave all eleven compounds as one cocktail and fitted them
          simultaneously by nonlinear mixed effects; the earlier studies dosed
          one compound at a time and fitted non-compartmentally.

Run:  python3 scripts/argoul_vs_earlier_mouse.py
"""

import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (compound, sex, source, dose mg/kg, CL mL/kg/d, Vd mL/kg, note)
# Argoul rows from db/primary_2026/argoul2026_mouse_tk.csv (Vss).
# Comparators from db/combined/tk_parameters.csv and db/primary_2026/.
DATA = [
    ("PFOA", "female", "Argoul 2026",      0.080,   4.5,   89,
     "cocktail, IV+oral, NLME, Vss"),
    ("PFOA", "female", "Lou 2009",         None,    5.99, 135,
     "serum; CL derived here as ke*Vd"),
    ("PFOA", "female", "Fujii (via EPA)",  None,   11.8,  150, "IV"),
    ("PFOA", "male",   "Fujii (via EPA)",  None,   14.2,  180, "IV"),
    ("PFOA", "male",   "Lou 2009",         None,    7.21, 226,
     "serum; CL derived here as ke*Vd"),

    ("PFHxS", "female", "Argoul 2026",     0.094,   1.3,   62,
     "cocktail, IV+oral, NLME, Vss"),
    ("PFHxS", "female", "Sundstrom 2012",  1.0,     2.68,  96,  "oral"),
    ("PFHxS", "female", "Sundstrom 2012", 20.0,     3.79, 147,  "oral"),
    ("PFHxS", "male",   "Sundstrom 2012",  1.0,     2.94, 129,  "oral"),
    ("PFHxS", "male",   "Sundstrom 2012", 20.0,     4.83, 195,  "oral"),

    ("PFHxA", "female", "Argoul 2026",     1.28, 6830,   4000,
     "cocktail, IV+oral, NLME, Vss"),
    ("PFHxA", "female", "US EPA PFHxA IRIS", None, None,  780,  "compilation"),
    ("PFHxA", "male",   "US EPA PFHxA IRIS", None, None,  750,  "compilation"),
]


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print("=" * 76)
    print("PART 1.  Side by side, matched on sex where possible")
    print("=" * 76)
    for comp in ("PFOA", "PFHxS", "PFHxA"):
        rows = [d for d in DATA if d[0] == comp]
        print(f"\n  {comp}")
        print(f"    {'sex':<7s} {'source':<20s} {'dose':>7s} "
              f"{'CL':>8s} {'Vd':>7s}   note")
        for _, sex, src, dose, cl, vd, note in rows:
            print(f"    {sex:<7s} {src:<20s} "
                  f"{(f'{dose:g}' if dose is not None else '-'):>7s} "
                  f"{(f'{cl:g}' if cl is not None else '-'):>8s} "
                  f"{(f'{vd:g}' if vd is not None else '-'):>7s}   {note}")

    print()
    print("=" * 76)
    print("PART 2.  Argoul against each comparator, same sex only")
    print("=" * 76)
    print(f"\n  {'compound':<8s} {'comparator':<30s} {'CL ratio':>9s} "
          f"{'Vd ratio':>9s}")
    for comp in ("PFOA", "PFHxS", "PFHxA"):
        rows = [d for d in DATA if d[0] == comp]
        arg = next(r for r in rows if r[2].startswith("Argoul"))
        for r in rows:
            if r is arg or r[1] != arg[1]:
                continue
            lab = r[2] + (f" {r[3]:g} mg/kg" if r[3] is not None else "")
            clr = f"{arg[4]/r[4]:.2f}x" if (arg[4] and r[4]) else "-"
            vdr = f"{arg[5]/r[5]:.2f}x" if (arg[5] and r[5]) else "-"
            print(f"  {comp:<8s} {lab:<30s} {clr:>9s} {vdr:>9s}")

    print("""
  Argoul is LOWER than every same-sex comparator for PFOA and PFHxS, on both
  terms, by factors of 1.3-2.9. For PFHxA it is far HIGHER on volume (5.1x).
  So the disagreement is not a constant offset, and it is not noise either.""")

    print()
    print("=" * 76)
    print("PART 3.  Is it dose? The one comparator with its own dose series")
    print("=" * 76)
    s = [d for d in DATA if d[0] == "PFHxS" and d[1] == "female"]
    s.sort(key=lambda r: r[3])
    print(f"\n  PFHxS, female mouse, three doses spanning {s[-1][3]/s[0][3]:.0f}x:\n")
    print(f"    {'dose mg/kg':>10s} {'CL':>7s} {'Vd':>6s}   source")
    for _, _, src, dose, cl, vd, _n in s:
        print(f"    {dose:>10g} {cl:>7g} {vd:>6g}   {src}")
    import math
    ld = [math.log10(r[3]) for r in s]
    lc = [math.log10(r[4]) for r in s]
    lv = [math.log10(r[5]) for r in s]

    def slope(x, y):
        n = len(x); mx = sum(x)/n; my = sum(y)/n
        return (sum((a-mx)*(b-my) for a, b in zip(x, y))
                / sum((a-mx)**2 for a in x))
    print(f"\n    slope of log CL on log dose: {slope(ld, lc):+.3f}")
    print(f"    slope of log Vd on log dose: {slope(ld, lv):+.3f}")
    print("""
    Both rise with dose, monotonically, across three points from two
    laboratories. SUPERSEDED, read with care: that fit includes Argoul's own
    point, which carries 79% leverage on the slope, so "Argoul sits at the
    bottom of the trend" is largely true by construction. See
    scripts/argoul_dose_matching.py for the out-of-sample test, which leaves
    about 1.6x unexplained after the dose correction.
    A positive clearance slope is what saturable REABSORPTION predicts
    (section 8): more escapes when the transporter is swamped. Dose therefore
    explains the direction of the PFOA and PFHxS disagreements, and with only
    three points it is a consistency check rather than a fitted relationship.

    It does not explain PFHxA, where Argoul's volume is 5x HIGHER. PFHxA is
    also the compound Argoul finds in net tubular secretion, so its kinetics
    differ in kind from the others and a single explanation should not be
    expected to cover it.""")


if __name__ == "__main__":
    main()
