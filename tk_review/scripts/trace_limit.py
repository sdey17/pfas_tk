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

One further honesty: where the protein concentration is unknown the free
concentration cannot be had and the total is used instead, which is an
UPPER BOUND and is labelled as one.

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


def assess(kd_uM, n, protein_uM, total_uM, stated=""):
    """(C_w uM, C_w/K_d, cost in log, basis, verdict) for one row.

    Returns None in the first three slots where the inputs are missing,
    because an unknown is not a zero.
    """
    if not kd_uM or kd_uM <= 0:
        return None, None, None, "no Kd", "unknown"
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
    if not total_uM or total_uM <= 0:
        return (None, None, None,
                "the assay's PFAS concentration is not in molar units",
                "unknown")
    if protein_uM and protein_uM > 0:
        cw = free_uM(kd_uM, n if n else 1.0, protein_uM, total_uM)
        basis = "free concentration by mass balance on the assay's own P"
        bounded = False
    else:
        cw = total_uM
        basis = ("protein concentration not in molar units, so the TOTAL "
                 "PFAS is used: an upper bound on the error, not the error")
        bounded = True
    cost = cost_log(kd_uM, cw)
    return cw, cw / kd_uM, cost, basis, verdict(cost, bounded)


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

    print("\n" + ("all checks passed" if ok else "CHECKS FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    if "--check" in sys.argv:
        raise SystemExit(_check())
    print(__doc__)
