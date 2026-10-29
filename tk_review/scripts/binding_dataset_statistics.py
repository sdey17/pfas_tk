#!/usr/bin/env python3
"""How much data is there, on what, and where does the conversion work?

Descriptive statistics for db/protein_binding.csv, plus the one result
that turns the Kd -> K conversion from a caveat into a usable rule.

THE RULE, for the impatient: the conversion works below nine
perfluorinated carbons and fails above it. All eleven converted values
that land within the between-laboratory floor are nCF <= 8, and not one
of the six pairs at nCF >= 9 does. That is exactly what the flat
specific site predicts -- at short chain lengths the saturable site IS
the binding, so a single-site Kd describes the whole of it; by C9 the
non-specific term has overtaken it and a single-site Kd describes a
minority of the capacity.

Run:  python3 tk_review/scripts/binding_dataset_statistics.py
"""
import collections
import csv
import importlib.util
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
PAIRS = os.path.join(os.path.dirname(HERE), "db",
                     "logk_converted_vs_measured.csv")

BETWEEN_LAB, WITHIN_LAB = 0.59, 0.19


def rule(ch="=", n=78):
    print(ch * n)


def is_albumin(r):
    p = r["protein"].lower()
    return "serum albumin" in p and "globulin" not in p


def is_single(r):
    """One named compound, not a prose list of several."""
    c = r["chemical"]
    return "," not in c and len(c) < 30


def part1(rows):
    rule()
    print("1. HOW MANY PFAS BIND ALBUMIN")
    rule()
    alb = [r for r in rows if is_albumin(r)]
    single = [r for r in alb if is_single(r)]
    kd = {r["chemical"] for r in single if r["kd"].strip()}
    ka = {r["chemical"] for r in single if r["ka"].strip()}
    print(f"  {len(alb):4d} albumin rows, {len({r['chemical'] for r in alb})} "
          f"distinct chemical entries")
    print(f"  {len(single):4d} of those name ONE compound "
          f"({len({r['chemical'] for r in single})} distinct); the rest are "
          f"prose lists")
    print()
    print(f"  {len(kd | ka):4d} compounds carry a numeric affinity on albumin")
    print(f"       {len(kd):3d} with a dissociation constant")
    print(f"       {len(ka):3d} with an association constant")
    print(f"       {len(kd & ka):3d} with both")
    print()
    print("  'Binds albumin' is doing some work in that sentence. A few of")
    print("  these are explicit NON-binders measured and reported as such,")
    print("  which is a result rather than an absence: 6:2 FTOH is a")
    print("  non-binder in all four species Starnes et al. 2024 tested.")
    print()


def part2(rows):
    rule()
    print("2. THE COMPOUNDS EVERYONE MEASURES")
    rule()
    alb = [r for r in rows if is_albumin(r) and is_single(r)
           and (r["kd"].strip() or r["ka"].strip())]
    n = collections.Counter(r["chemical"] for r in alb)
    st = collections.defaultdict(set)
    for r in alb:
        st[r["chemical"]].add((r["study"], r["year"]))
    print(f"  {'compound':16s}{'rows':>6s}{'studies':>9s}")
    for c, v in n.most_common(12):
        print(f"  {c:16s}{v:>6d}{len(st[c]):>9d}")
    print()
    top = [c for c, _ in n.most_common(10)]
    print(f"  Those ten carry {sum(n[c] for c in top)} of {len(alb)} albumin "
          f"affinity rows ({100 * sum(n[c] for c in top) / len(alb):.0f}%).")
    once = sum(1 for c, v in n.items() if v == 1)
    print(f"  {once} compounds appear exactly once "
          f"({100 * once / len(n):.0f}% of the {len(n)} measured).")
    print()
    print("  The distribution is the point: PFOA and PFOS are measured 38")
    print("  and 30 times across 21 and 16 studies, and the long tail is")
    print("  measured once. Any statement about 'PFAS binding to albumin'")
    print("  rests on two compounds unless it says otherwise.")
    print()


def part3():
    rule()
    print("3. HOW MANY PARTITION COEFFICIENTS EXIST")
    rule()
    spec = importlib.util.spec_from_file_location(
        "cmp", os.path.join(HERE, "logk_converted_vs_measured.py"))
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)
    _, out, _, _, _ = M.B.build()
    meas = [r for r in out if M.route(r) == "measured"]
    conv = [r for r in out if M.route(r) == "converted"]
    bio = [r for r in out if M.route(r) == "biomimetic"]
    print(f"  {len(meas):4d} MEASURED protein-water coefficients, on a stated")
    print(f"       basis, covering {len({r['pfas'] for r in meas})} compounds "
          f"from {len({r['source'] for r in meas})} sources")
    print(f"  {len(bio):4d} biomimetic, basis never stated by the paper")
    print(f"  {len(conv):4d} converted here from an affinity and a "
          f"stoichiometry")
    print()
    print("  the measured ones, by what they were measured in:")
    for k, v in collections.Counter(r["protein"] for r in meas).most_common():
        print(f"    {v:4d}  {k}")
    print()
    aff = {r["pfas"] for r in out if r["measure"] in ("Kd", "Ka")}
    par = {r["pfas"] for r in meas}
    print("  COVERAGE OVERLAP, by compound:")
    print(f"    {len(aff - par):3d} have an affinity but NO measured "
          f"partition coefficient")
    print(f"    {len(par - aff):3d} have a partition coefficient but no "
          f"affinity")
    print(f"    {len(aff & par):3d} have both -- the only ones where the "
          f"conversion can be checked")
    print()
    print("  Forty compounds with a measured coefficient against 192 with")
    print("  albumin data of some kind. The affinity literature is five")
    print("  times the size of the partition literature, which is why the")
    print("  conversion is worth having at all.")
    print()


def part4():
    rule()
    print("4. WHERE THE CONVERSION ACTUALLY WORKS")
    rule()
    pairs = list(csv.DictReader(open(PAIRS, encoding="utf-8")))
    for p in pairs:
        p["res"] = float(p["residual"])
        p["n"] = int(p["nCF"]) if p["nCF"] else 0
    good = sorted(pairs, key=lambda p: abs(p["res"]))[:12]
    print("  THE BEST AGREEMENTS, converted against measured:")
    print(f"  {'compound':14s}{'nCF':>4s}{'converted':>11s}{'measured':>10s}"
          f"{'residual':>10s}   from")
    for p in good:
        print(f"  {p['pfas']:14s}{p['n']:>4d}"
              f"{float(p['converted_logK']):>11.2f}"
              f"{float(p['measured_logK']):>10.2f}{p['res']:>+10.2f}   "
              f"{p['converted_source']}")
    print()
    w = sum(1 for p in pairs if abs(p["res"]) <= BETWEEN_LAB)
    t = sum(1 for p in pairs if abs(p["res"]) <= WITHIN_LAB)
    print(f"  {w} of {len(pairs)} pairs agree within the {BETWEEN_LAB} log "
          f"between-laboratory floor;")
    print(f"  {t} agree within the {WITHIN_LAB} log within-laboratory floor.")
    print()
    print("  NOW SORT THEM BY CHAIN LENGTH, WHICH IS THE WHOLE RESULT:")
    print()
    print(f"  {'band':12s}{'pairs':>7s}{'median':>9s}"
          f"{'within 0.59':>14s}")
    for lo, hi, lab in ((0, 6, "nCF <= 6"), (7, 8, "nCF 7-8"),
                        (9, 99, "nCF >= 9")):
        s = [p["res"] for p in pairs if lo <= p["n"] <= hi]
        if s:
            k = sum(1 for x in s if abs(x) <= BETWEEN_LAB)
            print(f"  {lab:12s}{len(s):>7d}{statistics.median(s):>+9.2f}"
                  f"{k:>9d}/{len(s):<4d}")
    print()
    ok = sorted(p["n"] for p in pairs if abs(p["res"]) <= BETWEEN_LAB)
    print(f"  Every good match is nCF <= {max(ok)}: {ok}")
    print(f"  Not one of the {sum(1 for p in pairs if p['n'] >= 9)} pairs at "
          f"nCF >= 9 agrees within the floor.")
    print()
    print("  THE RULE. Convert below nine perfluorinated carbons; do not")
    print("  convert above it. And this is not a fitted cutoff -- it is")
    print("  what the flat specific site predicts. At short chain lengths")
    print("  the saturable site IS the binding, so a single-site Kd")
    print("  describes the whole of it. By C9 the non-specific term has")
    print("  overtaken it, and a single-site Kd describes a minority of")
    print("  the capacity. The slope table in")
    print("  binding_all_data_protein_specific.py is the mechanism; this")
    print("  is where it bites.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    part1(rows)
    part2(rows)
    part3()
    part4()
    return 0


if __name__ == "__main__":
    sys.exit(main())
