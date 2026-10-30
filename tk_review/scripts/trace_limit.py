#!/usr/bin/env python3
"""Is the trace-limit conversion applicable to THIS row? Answered per row.

THE QUESTION. The conversion this project runs on is

    K = 1000 n / [M (K_d + C_w)]        exact, nothing dropped
    K0 = 1000 n / (M K_d)               what we actually report

and the second line is the first with C_w dropped beside K_d. So every
converted value carries an error of exactly

    log K0 - log K = log10(1 + C_w/K_d)

which is zero when C_w << K_d and is NOT otherwise. The ratio C_w/K_d is
the whole test, and it is a number, not a judgement -- so it should be a
column rather than a paragraph.

WHICH C_w, AND WHY IT IS NOT THE DOSED CONCENTRATION. C_w is the FREE
aqueous concentration at equilibrium, and a binding assay's stated PFAS
concentration is the TOTAL: free plus bound. Using the total overstates
the error, sometimes by a lot, because the protein is there precisely to
take ligand out of solution. With one class of n sites on a protein at
P, mass balance is

    L_total = C_w + n P C_w / (K_d + C_w)

one quadratic in C_w with one positive root, solved exactly below. No
iteration, no approximation.

    C_w^2 + C_w (K_d + n P - L_total) - K_d L_total = 0

WHAT THE ANSWER IS AND IS NOT. This is the single most important thing
about this module, and getting it backwards is easy.

A fitted K_d is a PARAMETER OF A MODEL, not a measurement taken at a
concentration. Every albumin row in this database that reaches a molar
PFAS concentration states a titration SERIES -- "0.0625-2 uM (0.0625,
0.125, 0.25, 0.5, 1, 2)", "0-1000 uM", "0-40 uM" -- so in every case the
series starts at or near zero and the fitted K_d is already the limiting
constant. K0 = 1000 n/(M K_d) is therefore the right number to report,
whatever the top of the series was. The trace limit is not an
approximation being made on top of the paper's work; it is the paper's
own model evaluated where that model is defined.

So this module does NOT compute an error in the conversion. It computes
WHERE THE ASSAY'S TOP POINT SAT ON THE ISOTHERM, which answers a
different and sharper question: how far the K0 we report sits from the
apparent coefficient at the loading the experiment actually used. That
is the number that matters the moment a converted K0 is compared against
someone else's directly measured coefficient, because that measurement
was made at ITS own loading and is not a K0 at all.

Both tails of it are flags, in opposite directions:

  a large value means the top point is deep into saturation, so the fit
  saw curvature and K_d is well constrained, but the coefficient at that
  loading is far below the reported K0;

  a value near zero means the series never reached K_d, the isotherm
  never curved, and the fitted K_d is an extrapolation beyond the data.
  `kd_over_top_conc` above 1 is the same fact stated directly.

FOUR ROUTES, and three of them are bounds rather than values. Where the
protein concentration is unknown the total stands in for the free; where
only a molar ligand:protein ratio is stated, C_w/K_d = ratio/n with the
protein cancelling out. That last one needs no K_d, which is what makes
it the only route open to a DIRECTLY MEASURED coefficient -- a log D has
no K_d to divide by, so without it there would be no way to ask whether
the measurement itself was taken at zero loading. Fischer's "PFAS:protein
molar ratio < 0.004" is enough, and the answer it gives is an upper
bound, so "within 0.002 log of zero loading" is a conclusion that holds
however the experiment was actually loaded.

Run:  python3 tk_review/scripts/trace_limit.py --check
"""
import math
import sys

# Below this the trace limit is exact to the second decimal of a log.
NEGLIGIBLE = 0.01
# Above this the reported K0 and the coefficient at the assay's own
# concentration are further apart than two labs measuring the same thing.
SERIOUS = 0.30


def free_uM(kd_uM, n, protein_uM, total_uM):
    """Free ligand at equilibrium, uM. One class of n sites.

    The positive root of C^2 + C (Kd + nP - Ltot) - Kd Ltot = 0. With
    Kd, Ltot > 0 the constant term is negative, so exactly one root is
    positive and no sign test is needed.
    """
    b = kd_uM + n * protein_uM - total_uM
    return (-b + math.sqrt(b * b + 4.0 * kd_uM * total_uM)) / 2.0


def cost_log(kd_uM, cw_uM):
    """log10(1 + Cw/Kd) -- what dropping Cw costs, in log units."""
    return math.log10(1.0 + cw_uM / kd_uM)


def verdict(cost, bounded):
    """Where the assay's TOP point sat on the isotherm. Both tails flag.

    The names say what the isotherm was doing and not how good the row
    is, because the two tails are opposite faults and neither is "good":

      saturating        the top point is well past K_d, so the fit saw
                        curvature and K_d is constrained -- but a
                        coefficient measured at that loading is far
                        below the K0 this conversion reports.
      curving           the useful middle.
      linear to the top the series never reached K_d at all, so the
                        isotherm never curved and the fitted K_d is an
                        extrapolation beyond the data. A small number in
                        this column is NOT reassurance.
    """
    if cost is None:
        return "unknown"
    tag = ("linear to the top" if cost < NEGLIGIBLE else
           "curving" if cost < SERIOUS else "saturating")
    return tag + " (upper bound)" if bounded else tag


# The low-loading shortcut below is good to ~10% of itself while
# C_w/K_d stays under this; past it the neglected K_d term matters and the
# protein concentration is genuinely needed.
RATIO_LIMIT = 0.1


def from_ratio(ratio, n):
    """C_w/K_d from a stated ligand:protein ratio alone, or None.

    Several papers state the loading only as a molar ratio -- Fischer's
    "PFAS:protein molar ratio < 0.004", Allendorf's "held below 0.1" --
    and give no PFAS concentration at all. That is enough, in the regime
    those papers work in, and the algebra says why.

    While C_w << K_d, mass balance gives C_w = L / (1 + nP/K_d), so

        C_w/K_d = L / (K_d + nP)  ->  ratio / n   when nP >> K_d

    which is independent of P: the protein cancels. The step needs
    C_w << K_d to have been true, and the result IS C_w/K_d, so the
    estimate checks its own premise -- a ratio of 0.004 returns 0.004 and
    is self-consistent, a ratio of 4 returns 4 and is not. So this
    refuses above RATIO_LIMIT rather than returning a number it has just
    shown to be invalid, which is also why it cannot be used to rescue
    the high-ratio assays (Bischel's 4:1, MacManus-Spencer's 5:1): there
    the sites saturate, C_w tends to L - nP, and nothing short of the
    protein concentration will do.

    WHAT IT RETURNS IS AN UPPER BOUND, provably and not approximately.
    Dropping K_d from L/(K_d + nP) can only make the quotient larger, so
    ratio/n is above the truth by exactly the factor 1 + K_d/(nP) -- 5.7%
    at Allendorf's lowest protein concentration, under 1% at Fischer's.
    It is labelled as a bound for that reason, and the conclusion drawn
    from it ("this measurement is within X log of zero loading") is the
    stronger statement for being one.
    """
    if not ratio or ratio <= 0:
        return None
    est = ratio / (n if n else 1.0)
    return est if est <= RATIO_LIMIT else None


def assess(kd_uM, n, protein_uM, total_uM, stated="", ratio=None):
    """(C_w uM, C_w/K_d, cost in log, basis, verdict) for one row.

    Four routes, in order of how much they assume. Returns None in the
    first three slots where no route reaches a number, because an
    unknown is not a zero.

      1. K_d, a molar total and a molar protein concentration -> exact
         mass balance. The only route that yields C_w itself.
      2. the same, with the loading stated as a ratio instead of a
         concentration -> ratio x P IS the total, so route 1 exactly.
      3. K_d and a molar total but no protein concentration -> the total
         stands in for the free. An upper bound.
      4. NO K_d at all, or no molar total, but a low stated ratio ->
         C_w/K_d = ratio/n with P cancelling. Also an upper bound, and
         the only route open to a directly measured coefficient, which
         reports no K_d to divide by.
    """
    nn = n if n else 1.0
    if "syringe" in (stated or "").lower():
        # Isothermal titration calorimetry. The syringe concentration is
        # what is INJECTED, 20 to 40x the cell concentration by design,
        # and the ligand in the cell ramps from zero over the titration.
        # Putting a syringe number where a cell number belongs would read
        # 4000 uM against a 110 uM protein and call every ITC row
        # saturated, which is an artefact of the wrong input.
        return (None, None, None,
                "an ITC syringe concentration, which is not the "
                "concentration in the cell", "unknown (ITC)")

    # Route 2 reduces to route 1 the moment the product is formed.
    from_ratio_total = False
    if (not total_uM or total_uM <= 0) and ratio and protein_uM:
        total_uM = ratio * protein_uM
        from_ratio_total = True

    # The ratio, however it happens to be known.
    r = ratio
    if r is None and total_uM and protein_uM and protein_uM > 0:
        r = total_uM / protein_uM

    if kd_uM and kd_uM > 0 and total_uM and total_uM > 0:
        if protein_uM and protein_uM > 0:
            cw = free_uM(kd_uM, nn, protein_uM, total_uM)
            basis = "free concentration by mass balance on the assay's own P"
            if from_ratio_total:
                basis += ", total ligand from the stated ratio x that P"
            bounded = False
        else:
            cw = total_uM
            basis = ("protein concentration not in molar units, so the "
                     "TOTAL PFAS is used: an upper bound on the error, "
                     "not the error")
            bounded = True
        cost = cost_log(kd_uM, cw)
        return cw, cw / kd_uM, cost, basis, verdict(cost, bounded)

    # Route 4. Needs no K_d, which is what makes it the only route open
    # to a directly measured coefficient.
    est = from_ratio(r, nn)
    if est is not None:
        return (None, est, cost_log(1.0, est),
                f"from the stated ligand:protein ratio alone ({r:g}): "
                f"Cw/Kd is at most ratio/n in the low-loading limit, "
                f"whatever the protein concentration was",
                verdict(cost_log(1.0, est), True) + " from the ratio")
    if r:
        return (None, None, None,
                f"stated only as a ligand:protein ratio of {r:g}, too "
                f"high for the low-loading shortcut; the protein "
                f"concentration in molar units is needed", "unknown")
    if "not reported" in (stated or "").lower():
        # Distinct from a unit this parser cannot read: the primary
        # source does not carry the number, so no extraction fixes it.
        return (None, None, None,
                "the paper does not report the concentrations it tested",
                "unknown (not reported)")
    if not kd_uM or kd_uM <= 0:
        return None, None, None, "no Kd and no stated ratio", "unknown"
    return (None, None, None,
            "the assay's PFAS concentration is not in molar units",
            "unknown")


def _check():
    """Properties that must hold, checked rather than asserted in prose."""
    ok = True

    def say(name, got, want, tol=1e-9):
        nonlocal ok
        good = abs(got - want) <= tol * max(1.0, abs(want))
        ok = ok and good
        print(f"  {'ok  ' if good else 'FAIL'} {name:52s} "
              f"{got:.6g}  (expected {want:.6g})")

    # 1. No protein: nothing binds, so free is total.
    say("no protein, free == total", free_uM(10.0, 1.0, 0.0, 5.0), 5.0)

    # 2. Mass balance, on a case where binding is substantial.
    kd, n, p, lt = 9.06, 1.0, 1.0, 2.0
    cw = free_uM(kd, n, p, lt)
    say("mass balance closes (Chen-like row)",
        cw + n * p * cw / (kd + cw), lt, 1e-12)

    # 3. The textbook approximation, which holds only when Cw << Kd.
    say("vs Ltot/(1+nP/Kd) in the trace regime",
        free_uM(100.0, 1.0, 1.0, 0.01), 0.01 / (1 + 1.0 / 100.0), 1e-3)

    # 4. Strong binding with sites to spare drives free well below total.
    cw = free_uM(0.01, 1.0, 10.0, 1.0)
    print(f"  ok   {'strong binding, excess sites: free < total':52s} "
          f"{cw:.6g}  (total 1, bound {1 - cw:.4f})")
    ok = ok and cw < 0.01

    # 5. The cost is monotone in Cw and hits the known log-2 point.
    say("cost at Cw == Kd is log10(2)", cost_log(5.0, 5.0), math.log10(2.0))
    say("cost at Cw == 0", cost_log(5.0, 0.0), 0.0)

    # 6. The two in vivo anchors quoted in the write-up. PFOA on albumin
    #    at Kd 2.57 uM (Chen 2025), serum albumin 600 uM, PFOA molar mass
    #    414.07 g/mol: 2 ng/mL is 4.83e-6 uM, 1000 ng/mL is 2.41e-3 uM.
    for ng_per_mL in (2.0, 1000.0):
        tot = ng_per_mL / 414.07 * 1e-3
        cw, ratio, cost, _b, v = assess(2.57, 1.0, 600.0, tot)
        print(f"  ok   in vivo, {ng_per_mL:6.0f} ng/mL PFOA, albumin 600 uM"
              f"      occupancy {100 * cw / (2.57 + cw):.5f}%  "
              f"cost {cost:.5f} log  -> {v}")
        ok = ok and cost < NEGLIGIBLE

    # 7. An in vitro case that is NOT in the trace limit, as a contrast.
    cw, ratio, cost, _b, v = assess(2.57, 1.0, 1.0, 100.0)
    ok = ok and cost > SERIOUS
    # An ITC syringe concentration must not be read as a cell concentration.
    _c, _r, _k, _bas, v_itc = assess(2.4, 1.6, 110.0, 4000.0,
                                     "4000 (syringe) uM")
    print(f"  {'ok  ' if v_itc == 'unknown (ITC)' else 'FAIL'} "
          f"{'ITC syringe concentration refused':52s} -> {v_itc}")
    ok = ok and v_itc == "unknown (ITC)"
    print(f"  ok   in vitro, 100 uM on 1 uM albumin"
          f"                     Cw/Kd {ratio:.1f}  "
          f"cost {cost:.2f} log  -> {v}")
    ok = ok and cost > SERIOUS

    # 8. The ratio-only route. Fischer states no concentration at all,
    #    only "PFAS:protein molar ratio < 0.004", and the answer should
    #    match the explicit mass-balance one across any plausible P.
    _c, r_only, c_only, _b, v = assess(2.57, 1.0, None, None,
                                       "molar ratio mol/mol", 0.004)
    #    It must BOUND that answer from above, by exactly 1 + Kd/(nP).
    for P in (45.0, 200.0, 600.0):
        _c2, r_exp, c_exp, _b2, _v2 = assess(2.57, 1.0, P, 0.004 * P)
        bound = 1.0 + 2.57 / P
        good = r_exp <= r_only <= r_exp * bound * 1.0001
        ok = ok and good
        print(f"  {'ok  ' if good else 'FAIL'} "
              f"ratio-only 0.004 bounds mass balance at P={P:5.0f} uM   "
              f"{r_only:.5f} >= {r_exp:.5f}, within {100 * (bound - 1):.1f}%")
    say("ratio-only cost is independent of P", c_only, math.log10(1.004),
        1e-6)

    # 9. Exact route: a ratio AND a molar protein concentration.
    _c3, r3, c3, b3, _v3 = assess(2.57, 1.0, 45.1379, None,
                                  "held below 0.1 mol/mol", 0.1)
    _c4, r4, _c5, _b4, _v4 = assess(2.57, 1.0, 45.1379, 0.1 * 45.1379)
    say("ratio x protein equals the explicit total", r3, r4)
    good = "total ligand from the stated ratio" in b3
    ok = ok and good
    print(f"  {'ok  ' if good else 'FAIL'} "
          f"{'and the basis says the total came from the ratio':52s}")

    # 10. A high ratio must be refused, not estimated.
    _c6, r6, _c7, b6, v6 = assess(10.0, 1.0, None, None, "4:1 mol/mol", 4.0)
    print(f"  {'ok  ' if r6 is None and v6 == 'unknown' else 'FAIL'} "
          f"{'a 4:1 ratio is refused, not estimated':52s} {b6[:44]}...")
    ok = ok and r6 is None and v6 == "unknown"

    print("\n" + ("all checks passed" if ok else "CHECKS FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    if "--check" in sys.argv:
        raise SystemExit(_check())
    print(__doc__)
