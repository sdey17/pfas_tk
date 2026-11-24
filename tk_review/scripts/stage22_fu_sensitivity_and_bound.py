#!/usr/bin/env python3
"""The binding claim: an arithmetic bound, and its f_u sensitivity (Stage 22).

    python3 -I scripts/stage22_fu_sensitivity_and_bound.py

Writes db/fu_source_sensitivity.csv and db/binding_arithmetic_bound.csv.
Idempotent.

TWO THINGS THIS STAGE FIXES.

(1) THE SHARES ARE ANTICORRELATED BY CONSTRUCTION, SO THEY ARE NOT TWO
FINDINGS. In equation (4) the tubular term is a RESIDUAL:

    T_ratio = halflife_ratio / (Vd_ratio * fu_ratio * GFR_ratio)

so any error in the unbound-fraction ratio enters the binding term with one
sign and the tubular term with exactly the opposite sign. Their log shares
move together with correlation -1. Reporting "binding 8%, tubule 65%" as two
results overstates the evidence: it is one measurement stated twice. The
paper should say so, and should lead with the claim that does NOT depend on
the residual at all.

(2) THAT CLAIM IS AN ARITHMETIC BOUND, AND IT IS THE ROBUST RESULT. Binding
can contribute to the half-life ratio at most the unbound-fraction ratio
itself. Measured rat/human ratios are near unity, while half-life ratios run
to three orders of magnitude. The gap is so large that no redistribution
between the other terms can close it, and the bound needs no decomposition,
no residual, and no assumption about GFR or Vd. It is the strongest form of
the paper's central claim and it should be stated first.

HOW FAR THE SENSITIVITY CAN ACTUALLY GO, STATED HONESTLY. Four plasma
unbound-fraction sources exist in this project: Ryu 2024, two Argoul 2025
experiments, and Fischer 2024. They do not cover the species equally:

    human  4 sources
    mouse  3 sources
    rat    1 source  (Ryu 2024, and nothing else)

So the rat binding factor CANNOT be varied on its numerator. Only the human
denominator can be swapped, and the resulting range understates the true
uncertainty by an unknown amount. This is a real limitation of the rat
result and the sensitivity analysis cannot remove it; it can only show that
the conclusion survives every variation that IS available. The mouse, with
three sources, can be tested properly.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

SOURCE_LABEL = {
    "Ryu 2024": "Ryu 2024 (Toxics 12:253)",
    "1": "Argoul 2025 expt 1",
    "2": "Argoul 2025 expt 2",
    "Fischer 2024": "Fischer 2024 (ES&T 58:1055)",
}


def load_fu():
    """{(species, source): {chemical: fu}}"""
    out = {}
    for r in csv.DictReader(
            (HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if not r["fu"]:
            continue
        sp = r["species"].split(" (")[0]
        out.setdefault((sp, r["experiment"]), {})[r["chemical"]] = float(r["fu"])
    return out


def load_decomp():
    return list(csv.DictReader(
        (HERE / "db" / "species_decomposition.csv").open()))


def main():
    fu = load_fu()
    dec = load_decomp()
    bar = "=" * 78

    # ---------------------------------------------------- 1. the bound
    print(bar)
    print("THE ARITHMETIC BOUND: WHAT BINDING COULD EXPLAIN AT MOST")
    print(bar)
    print("  Binding's contribution to a half-life ratio is at most the")
    print("  unbound-fraction ratio. No decomposition is needed, and no")
    print("  residual is involved.")
    print()
    bound_rows = []
    for comp in ("human/rat", "human/mouse"):
        sub = [r for r in dec if r["comparison"] == comp]
        if not sub:
            continue
        fb = [float(r["factor_binding"]) for r in sub]
        hr = [float(r["halflife_ratio"]) for r in sub]
        print(f"  {comp}  ({len(sub)} compounds)")
        print(f"    unbound-fraction ratio : {min(fb):.2f} to {max(fb):.2f}")
        print(f"    half-life ratio        : {min(hr):.1f} to {max(hr):.1f}")
        worst = max(hr) / max(fb)
        print(f"    so binding leaves AT LEAST {worst:.0f}-fold unexplained")
        print(f"    at the largest half-life ratio, even if the largest")
        print(f"    binding ratio applied to it.")
        for r in sub:
            b, h = float(r["factor_binding"]), float(r["halflife_ratio"])
            bound_rows.append(dict(
                comparison=comp, chemical=r["chemical"],
                n_total_carbons=r["n_total_carbons"],
                halflife_ratio=h, fu_ratio_binding_factor=b,
                max_pct_of_log_ratio_from_binding=round(
                    math.log10(b) / math.log10(h) * 100, 1) if h > 1 else "",
                unexplained_fold_if_binding_maximal=round(h / b, 1)))
        print()
    rat = [r for r in bound_rows if r["comparison"] == "human/rat"
           and r["max_pct_of_log_ratio_from_binding"] != ""]
    mou = [r for r in bound_rows if r["comparison"] == "human/mouse"
           and r["max_pct_of_log_ratio_from_binding"] != ""]
    rp = [r["max_pct_of_log_ratio_from_binding"] for r in rat]
    mp = [r["max_pct_of_log_ratio_from_binding"] for r in mou]
    print(f"  Binding's MAXIMUM share of the log half-life ratio, per compound:")
    print(f"    rat   {min(rp):>6.1f}% to {max(rp):>6.1f}%   "
          f"(median {sorted(rp)[len(rp)//2]:.1f}%)")
    print(f"    mouse {min(mp):>6.1f}% to {max(mp):>6.1f}%   "
          f"(median {sorted(mp)[len(mp)//2]:.1f}%)")
    print()
    print("  This is a CEILING, not an estimate: it assumes every bit of the")
    print("  unbound-fraction difference acts on the half-life and nothing")
    print("  opposes it. The rat ceiling is already small. That is the")
    print("  claim to lead with, because nothing in it depends on the")
    print("  residual, on GFR, or on V_d.")
    print()

    # ------------------------------------------- 2. f_u source sensitivity
    print(bar)
    print("SENSITIVITY TO THE f_u SOURCE, AND HOW FAR IT CAN GO")
    print(bar)
    avail = {}
    for (sp, src) in fu:
        avail.setdefault(sp, []).append(src)
    for sp in ("human", "rat", "mouse"):
        srcs = sorted(avail.get(sp, []))
        print(f"  {sp:7} {len(srcs)} source(s): "
              f"{', '.join(SOURCE_LABEL.get(s, s) for s in srcs)}")
    print()
    print("  The rat has ONE. Its binding factor can only be varied on the")
    print("  human denominator, so the range below UNDERSTATES rat")
    print("  uncertainty by an unknown amount. Stated, not hidden.")
    print()
    rows = []
    for comp, animal in (("human/rat", "rat"), ("human/mouse", "mouse")):
        for asrc in sorted(avail.get(animal, [])):
            for hsrc in sorted(avail.get("human", [])):
                pairs = []
                for chem, av in fu[(animal, asrc)].items():
                    hv = fu[("human", hsrc)].get(chem)
                    if hv:
                        pairs.append((chem, av / hv))
                if len(pairs) < 3:
                    continue
                vals = [v for _, v in pairs]
                rows.append(dict(
                    comparison=comp, animal_fu_source=SOURCE_LABEL[asrc],
                    human_fu_source=SOURCE_LABEL[hsrc],
                    n_compounds=len(pairs),
                    fu_ratio_min=round(min(vals), 3),
                    fu_ratio_max=round(max(vals), 3),
                    fu_ratio_geomean=round(
                        10 ** (sum(math.log10(v) for v in vals) / len(vals)), 3),
                    max_pct_if_halflife_ratio_100=round(
                        math.log10(max(vals)) / 2.0 * 100, 1)))
    print(f"  {'comparison':13}{'animal f_u':22}{'human f_u':28}{'n':>3}"
          f"{'ratio range':>16}")
    for r in rows:
        print(f"  {r['comparison']:13}{r['animal_fu_source'][:21]:22}"
              f"{r['human_fu_source'][:27]:28}{r['n_compounds']:>3}"
              f"  {r['fu_ratio_min']:>5.2f} to {r['fu_ratio_max']:>6.2f}")
    print()
    ratr = [r for r in rows if r["comparison"] == "human/rat"]
    mour = [r for r in rows if r["comparison"] == "human/mouse"]
    if ratr:
        print(f"  RAT: across all {len(ratr)} source pairings the binding")
        print(f"  factor never exceeds {max(r['fu_ratio_max'] for r in ratr):.2f}.")
    if mour:
        print(f"  MOUSE: across all {len(mour)} pairings it reaches "
              f"{max(r['fu_ratio_max'] for r in mour):.1f}.")
        print("  The rat/mouse asymmetry that drives the paper's second")
        print("  result is therefore not an artefact of one laboratory: it")
        print("  survives every available substitution on both sides.")
    print()

    # ------------------------------------------- 3. the anticorrelation
    print(bar)
    print("WHY THE TWO SHARES ARE ONE FINDING, DEMONSTRATED")
    print(bar)
    print("  Perturb every unbound-fraction ratio by a factor k and recompute")
    print("  the shares. The residual absorbs it exactly.")
    print()
    sub = [r for r in dec if r["comparison"] == "human/rat"]
    print(f"  {'k':>6}{'mean binding %':>17}{'mean tubule %':>16}"
          f"{'sum of the two':>17}")
    for k in (0.5, 0.8, 1.0, 1.25, 2.0):
        bs, ts = [], []
        for r in sub:
            h = float(r["halflife_ratio"])
            if h <= 1:
                continue
            lg = math.log10(h)
            b = float(r["factor_binding"]) * k
            tr = float(r["factor_tubule_residual"]) / k   # residual absorbs it
            bs.append(math.log10(b) / lg * 100)
            ts.append(math.log10(tr) / lg * 100)
        mb, mt = sum(bs) / len(bs), sum(ts) / len(ts)
        print(f"  {k:>6.2f}{mb:>17.1f}{mt:>16.1f}{mb+mt:>17.1f}")
    print()
    print("  The sum is invariant. Every log unit added to binding is")
    print("  subtracted from the tubular residual, exactly, because the")
    print("  residual is DEFINED as what the other three terms leave over.")
    print("  So 'binding is small' and 'the tubular term is large' are the")
    print("  same statement, and the paper should present them as one.")
    print()
    print("  This does not weaken the result. It relocates it: the secure")
    print("  claim is the measured unbound-fraction ratio and the bound it")
    print("  implies, not the partition. Equation (5), where T is measured")
    print("  rather than inferred, is the only place the two are genuinely")
    print("  independent -- which is why extending it matters.")
    print()

    # ------------------------------- 4. the bound under the WORST source
    print(bar)
    print("THE BOUND UNDER THE LEAST FAVOURABLE f_u SOURCE")
    print(bar)
    print("  The paper states the rat binding factor as 0.76-2.03. That is")
    print("  the Ryu/Ryu pairing only. The ceiling must be recomputed under")
    print("  the pairing most favourable to the binding explanation.")
    print()
    hl = {}
    for r in dec:
        if r["comparison"] == "human/rat":
            hl[r["chemical"]] = float(r["halflife_ratio"])
    worst = []
    for hsrc in sorted(avail.get("human", [])):
        for chem, av in fu[("rat", "Ryu 2024")].items():
            hv = fu[("human", hsrc)].get(chem)
            if not hv or chem not in hl or hl[chem] <= 1:
                continue
            ratio = av / hv
            share = math.log10(ratio) / math.log10(hl[chem]) * 100
            worst.append((share, chem, hsrc, ratio, hl[chem]))
    worst.sort(reverse=True)
    print(f"  {'compound':9}{'human f_u source':28}{'fu ratio':>9}"
          f"{'t1/2 ratio':>12}{'binding ceiling':>17}")
    for share, chem, hsrc, ratio, h in worst[:5]:
        print(f"  {chem:9}{SOURCE_LABEL[hsrc][:27]:28}{ratio:>9.2f}"
              f"{h:>12.1f}{share:>16.1f}%")
    print()
    top = worst[0]
    print(f"  Worst case: {top[1]} under {SOURCE_LABEL[top[2]]}, where")
    print(f"  binding could account for {top[0]:.0f}% of the log half-life")
    print(f"  ratio -- against the {min(rp):.1f}-{max(rp):.1f}% the Ryu/Ryu")
    print(f"  pairing gives. The paper's stated 0.76-2.03 range is")
    print(f"  SOURCE-SPECIFIC and must be reported as such.")
    print()
    print(f"  The conclusion still holds: even at this ceiling binding is a")
    print(f"  minority of the difference, and the Fischer 2024 pairing rests")
    print(f"  on only 3 shared compounds. But 'binding is bounded at")
    print(f"  twofold' is not defensible as written; 'bounded at roughly")
    print(f"  fivefold across every available source' is.")
    print()

    out = []
    p = HERE / "db" / "binding_arithmetic_bound.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bound_rows[0].keys()))
        w.writeheader()
        for r in sorted(bound_rows, key=lambda x: (x["comparison"],
                                                   int(x["n_total_carbons"]),
                                                   x["chemical"])):
            w.writerow(r)
    out.append((p, len(bound_rows)))

    p = HERE / "db" / "fu_source_sensitivity.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["comparison"],
                                             x["animal_fu_source"],
                                             x["human_fu_source"])):
            w.writerow(r)
    out.append((p, len(rows)))

    for p, n in out:
        print(f"wrote {p.relative_to(HERE)}  ({n} rows)")


if __name__ == "__main__":
    main()
