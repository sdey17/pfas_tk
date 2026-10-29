#!/usr/bin/env python3
"""What extending past the top ten sources buys, and what differs by protein.

TWO QUESTIONS, and the second is the interesting one.

1. WHAT DOES EXTENDING THE SCOPE ADD? Very little that is usable, and
   the reason is structural rather than a gap in extraction: of the rows
   measuring a NAMED SINGLE protein, roughly 95 per cent are albumin.
   The bulk of what lies outside the top ten is phospholipid membrane
   (476 rows) and whole plasma (98 rows) -- different sorbents answering
   different questions, not more protein-binding data.

2. DOES BINDING DIFFER BY PROTEIN? Yes, and not in the way the level
   suggests. The LEVELS differ by about a log, which anyone would
   expect. What matters is that the CHAIN-LENGTH SLOPE differs by
   two-fold, and albumin is the outlier:

     sorbent                        per CF2       behaves like
     albumin, specific site         +0.026        a saturating site
     albumin, whole partition       +0.193        a mixture
     liver FABP (n = 4, thin)       +0.347        hydrophobic partition
     phospholipid membrane          +0.386        hydrophobic partition
     gamma-globulin                 +0.396        hydrophobic partition

   (carboxylates; the sulfonate slopes are steeper throughout and tell
   the same story.) Everything that is a hydrophobic partition climbs at
   0.35 to 0.40 per CF2, and three unrelated sorbents agree on it.
   Albumin's saturable site is FLAT. Albumin's whole partition sits
   between the two because it IS the two: a flat specific site plus a
   climbing non-specific term.

   That is the same result part 3 of logk_converted_vs_measured reached
   from the other end -- the stoichiometry implied by reconciling a Kd
   with a partition coefficient rises with chain length -- and it is why.

THE CONSEQUENCE, which is actionable. Albumin out-binds gamma-globulin
by 1.96 log at four perfluorinated carbons and by 0.37 log at twelve.
Every model that treats plasma binding as albumin-only is therefore
making a chain-length-dependent error, small for the short compounds and
not small for the long ones.

Every slope here is fitted WITHIN one paper, one method and one head
group. Across papers the same albumin constant moves 1.45 log, which
would swamp all of it.

Run:  python3 tk_review/scripts/binding_all_data_protein_specific.py
"""
import collections
import csv
import math
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

TOP10 = {"Fischer 2024", "Chen 2025", "Muensterman 2026", "Gao 2025",
         "Qin 2024 (plasma)", "Starnes 2024", "Starnes 2025",
         "Allendorf 2019", "Jia 2022", "Bischel 2011"}

# Both spellings of the C11 and C13 acids occur across sources.
PFCA = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
        "PFUnDA", "PFUnA", "PFDoDA", "PFDoA", "PFTrDA", "PFTriDA",
        "PFTeDA"]
PFSA = ["PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS"]
NUM = re.compile(r"^\s*([\d.]+(?:[eE][-+]?\d+)?)")

# Plasma protein split, carried as a range for the same reason the albumin
# concentration is: it varies between people more than any constant here.
# Albumin 35-50 g/L, total protein 60-80 g/L, so globulin is 20-35 g/L.
GLOB_OVER_ALB = (20.0 / 50.0, 35.0 / 35.0)      # 0.40x to 1.00x by mass


def rule(ch="=", n=78):
    print(ch * n)


def load():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    ncf = {}
    for r in rows:
        v = r["n_perfluorinated_carbons"].strip()
        if v.isdigit():
            ncf.setdefault(r["chemical"], int(v))
    return rows, ncf


def kd_uM(r):
    m = NUM.match(r["kd"].strip())
    if not m:
        return None
    u = (r["kd_units"].strip().lower().split() or [""])[0]
    s = {"m": 1e6, "mm": 1e3, "um": 1.0, "nm": 1e-3}.get(u)
    return float(m.group(1)) * s if s else None


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if not sxx:
        return 0.0, my, 0.0
    m = sxy / sxx
    syy = sum((y - my) ** 2 for y in ys)
    return m, my - m * mx, (sxy / math.sqrt(sxx * syy) if syy else 0.0)


def pclass(p):
    q = p.lower()
    for key, lab in (("serum albumin", "albumin"), ("globulin", "globulin"),
                     ("transthyretin", "transthyretin"),
                     ("l-fabp", "L-FABP"), ("liver fatty acid", "L-FABP"),
                     ("i-fabp", "I-FABP"),
                     ("intestinal fatty acid", "I-FABP"),
                     ("peroxisome", "PPAR"), ("ppar", "PPAR"),
                     ("glycoprotein", "AGP"),
                     ("plasma", "whole plasma/serum"),
                     ("serum", "whole plasma/serum"),
                     ("membrane", "phospholipid membrane"),
                     ("bilayer", "phospholipid membrane"),
                     ("liposome", "phospholipid membrane"),
                     ("phospholipid", "phospholipid membrane"),
                     ("cell", "cells/tissue"), ("tissue", "cells/tissue"),
                     ("homogenate", "cells/tissue")):
        if key in q:
            return lab
    return "other"


def part1(rows):
    rule()
    print("1. WHAT IS ACTUALLY IN THE FULL DATABASE")
    rule()
    c = collections.Counter(pclass(r["protein"]) for r in rows)
    named = {"albumin", "globulin", "L-FABP", "I-FABP", "PPAR",
             "transthyretin", "AGP"}
    tot_named = sum(v for k, v in c.items() if k in named)
    print(f"  {len(rows)} rows.\n")
    print(f"  {'sorbent':26s}{'rows':>6s}{'cmpds':>7s}{'studies':>9s}"
          f"{'share of named':>16s}")
    for k, v in c.most_common():
        sub = [r for r in rows if pclass(r["protein"]) == k]
        share = (f"{100 * v / tot_named:.1f}%" if k in named else "-")
        print(f"  {k:26s}{v:>6d}"
              f"{len({r['chemical'] for r in sub}):>7d}"
              f"{len({(r['study'], r['year']) for r in sub}):>9d}{share:>16s}")
    print()
    print(f"  Of the {tot_named} rows against a NAMED SINGLE PROTEIN, "
          f"{100 * c['albumin'] / tot_named:.0f}% are albumin.")
    print("  The field is not short of PFAS-protein data; it is short of")
    print("  PFAS-protein data about any protein except albumin.")
    print()
    print("  The two largest blocks outside the top ten are not protein")
    print("  binding at all: phospholipid membrane is the LIPID limb of a")
    print("  mass balance, and whole plasma is the thing a mass balance")
    print("  has to predict. Both are useful; neither adds a protein.")
    print()


def part2():
    rule()
    print("2. WHAT EXTENDING THE SCOPE COSTS AND BUYS")
    rule()
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "cmp", os.path.join(HERE, "logk_converted_vs_measured.py"))
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)
    _, out, _, _, _ = M.B.build()
    conv = [r for r in out if M.route(r) == "converted"]
    meas = [r for r in out if M.route(r) == "measured"]
    ci = sum(1 for r in conv if r["source"] in TOP10)
    mi = sum(1 for r in meas if r["source"] in TOP10)
    print(f"  convertible rows   {len(conv):4d} all   {ci:4d} top-10   "
          f"{len(conv) - ci:+4d}")
    print(f"  measured log K     {len(meas):4d} all   {mi:4d} top-10   "
          f"{len(meas) - mi:+4d}")
    print()
    extra = collections.Counter(r["source"] for r in conv
                                if r["source"] not in TOP10)
    print("  The 25 extra convertible rows, by source:")
    for s, v in extra.most_common():
        note = {"Wu 2024": "affinities 3-4 orders below every other "
                           "measurement of the same compound",
                "Peng 2024": "n from a fluorescence double-log slope",
                "Chen 2009": "n from a fluorescence double-log slope",
                "Moro 2022": "ITC, n fitted -- the best of these",
                "Khazaee 2021": "L-FABP, the only non-albumin protein "
                                "with a usable n",
                "Woodcroft 2010": "rat L-FABP, multi-site"}.get(s, "")
        print(f"    {v:4d}  {s:16s}{note}")
    print()
    print("  So extending the scope adds 30% more convertible rows, and")
    print("  16 of the 25 come from double-log slopes or from a source")
    print("  whose affinities are three orders out. The honest summary:")
    print("  extending the SOURCE list adds almost no reliability. What")
    print("  it adds is the only non-albumin protein that converts at")
    print("  all, which is why Khazaee and Woodcroft were worth having.")
    print()


def series(rows, ncf, doi, want, value, label, species=None):
    """(nCF, y, compound) points for one paper, one sorbent, one quantity.

    `species` is not optional in spirit: albumin is measured in two
    species in the same paper, and pooling them puts a species
    difference inside a chain-length slope.
    """
    out = []
    for r in rows:
        if r["doi"] != doi or pclass(r["protein"]) != want:
            continue
        if species and r["protein_species"] != species:
            continue
        c = r["chemical"]
        if c not in ncf:
            continue
        y = value(r)
        if y is None:
            continue
        out.append((ncf[c], y, c))
    return out


def part3(rows, ncf):
    rule()
    print("3. THE CHAIN-LENGTH SLOPE IS PROTEIN-SPECIFIC")
    rule()
    print("  Every fit below is WITHIN one paper, one method, one head")
    print("  group. Across papers the same albumin constant moves 1.45")
    print("  log, which would swamp every difference here.")
    print()

    def logd(r):
        try:
            return float(r["ka"])
        except ValueError:
            return None

    def negkd(r):
        v = kd_uM(r)
        return -math.log10(v) if v else None

    FISCHER = "10.1021/acs.est.3c07415"
    CHEN = "10.1021/acs.est.4c06734"
    JIA = "10.1021/acs.est.1c08493"
    # All human, so no species difference hides inside any slope.
    sets = [
        ("albumin, specific site (-log Kd)", "Chen 2025", CHEN, "albumin",
         negkd, "human"),
        ("albumin, whole partition (log D)", "Fischer 2024", FISCHER,
         "albumin", logd, "human"),
        ("phospholipid membrane (log K_MW)", "Chen 2025", CHEN,
         "phospholipid membrane", logd, None),
        ("gamma-globulin (log D)", "Fischer 2024", FISCHER, "globulin",
         logd, "human"),
        ("liver FABP (-log Kd)", "Jia 2022", JIA, "L-FABP", negkd, "human"),
    ]
    print(f"  {'sorbent':34s}{'source':16s}{'head group':12s}"
          f"{'per CF2':>9s}{'r':>7s}{'n':>4s}")
    rowsout = []
    for label, src, doi, want, val, sp in sets:
        pts = series(rows, ncf, doi, want, val, label, sp)
        for cls, S in (("carboxylate", PFCA), ("sulfonate", PFSA)):
            p = [(x, y) for x, y, c in pts if c in S]
            if len(p) < 4:
                continue
            m, b, r = fit([x for x, _ in p], [y for _, y in p])
            rowsout.append((label, cls, m, r, len(p)))
            print(f"  {label:34s}{src:16s}{cls:12s}{m:>+9.3f}{r:>+7.2f}"
                  f"{len(p):>4d}")
    print()
    print("  READ THE COLUMN, NOT THE ROWS. On the carboxylates, liver")
    print("  FABP, membrane and globulin climb at +0.347, +0.386 and")
    print("  +0.396 per CF2 -- three unrelated sorbents agreeing on what")
    print("  hydrophobic partitioning looks like. Albumin's SATURABLE")
    print("  SITE is flat at +0.026. Albumin's whole partition sits in")
    print("  between, at +0.193, because it is the sum of the two.")
    print()
    print("  The sulfonates are steeper everywhere and preserve the same")
    print("  ordering. The liver FABP fit rests on four compounds and is")
    print("  the one line here that should be treated as indicative.")
    print()
    print("  The membrane series is the cleanest (r = 1.00) and still not")
    print("  a straight line: its per-CF2 steps rise from +0.39 at C4 to")
    print("  +0.47 at C8 and fall to +0.25 at C13. Even pure lipid")
    print("  partitioning saturates in chain length eventually.")
    print()
    return rowsout


def part4(rows, ncf):
    rule()
    print("4. THE CONSEQUENCE: ALBUMIN'S ADVANTAGE DIES WITH CHAIN LENGTH")
    rule()
    print("  Fischer et al. 2024 measured albumin and gamma-globulin on")
    print("  the same compounds in the same hands, so the comparison")
    print("  needs no model at all -- just subtraction.")
    print()
    FISCHER = "10.1021/acs.est.3c07415"
    alb, glob = {}, {}
    for r in rows:
        if r["doi"] != FISCHER or not r["ka"].strip():
            continue
        try:
            v = float(r["ka"])
        except ValueError:
            continue
        if pclass(r["protein"]) == "albumin" and r["protein_species"] == "human":
            alb[r["chemical"]] = v
        elif pclass(r["protein"]) == "globulin":
            glob[r["chemical"]] = v
    shared = [c for c in alb if c in glob and c in ncf]
    shared.sort(key=lambda c: (ncf[c], c))
    print(f"  {'compound':12s}{'nCF':>4s}{'albumin':>9s}{'globulin':>10s}"
          f"{'advantage':>11s}{'globulin share of plasma binding':>34s}")
    for c in shared:
        gap = alb[c] - glob[c]
        lo = GLOB_OVER_ALB[0] * 10 ** (-gap)
        hi = GLOB_OVER_ALB[1] * 10 ** (-gap)
        print(f"  {c:12s}{ncf[c]:>4d}{alb[c]:>9.2f}{glob[c]:>10.2f}"
              f"{gap:>+11.2f}"
              + f"{100 * lo / (1 + lo):>14.0f}% to"
              + f"{100 * hi / (1 + hi):>4.0f}%")
    print()
    print("  The last column is globulin's share of total plasma protein")
    print("  binding at trace, from these two measurements and the mass")
    print(f"  split alone ({GLOB_OVER_ALB[0]:.2f}x to "
          f"{GLOB_OVER_ALB[1]:.2f}x globulin:albumin by mass, from the")
    print("  35-50 and 60-80 g/L clinical intervals). No partition")
    print("  coefficient and no model enters it.")
    print()
    short = [c for c in shared if ncf[c] <= 6]
    long = [c for c in shared if ncf[c] >= 10]
    if short and long:
        gs = statistics.median(alb[c] - glob[c] for c in short)
        gl = statistics.median(alb[c] - glob[c] for c in long)
        print(f"  Median albumin advantage: {gs:+.2f} log at nCF <= 6, "
              f"{gl:+.2f} log at nCF >= 10.")
        print(f"  It falls by {gs - gl:.2f} log across the measured range.")
    print()
    print("  SO: treating plasma binding as albumin-only is a")
    print("  chain-length-dependent error, not a constant one. For the")
    print("  short compounds globulin is a rounding error. For the long")
    print("  ones it is not, and the correction goes the way that makes")
    print("  the free fraction SMALLER than an albumin-only model says.")
    print()


def part5(rows):
    rule()
    print("5. WHAT STILL CANNOT BE ANSWERED, AT ANY SCOPE")
    rule()
    counts = collections.Counter(pclass(r["protein"]) for r in rows)
    thin = [(k, counts[k]) for k in
            ("transthyretin", "I-FABP", "AGP", "PPAR", "L-FABP")
            if counts[k]]
    print("  Proteins with too little data to fit anything:")
    for k, v in sorted(thin, key=lambda kv: kv[1]):
        sub = [r for r in rows if pclass(r["protein"]) == k]
        n = len({r["chemical"] for r in sub})
        print(f"    {k:16s}{v:4d} rows, {n:3d} compounds, "
              f"{len({(r['study'], r['year']) for r in sub})} studies")
    print()
    print("  Transthyretin is the sharpest gap. It carries thyroxine, the")
    print("  displacement hypothesis for PFAS thyroid effects runs through")
    print("  it, and this database holds TWO rows for it.")
    print()
    print("  No protein other than albumin has a measured protein-water")
    print("  partition coefficient anywhere in the database. Every")
    print("  non-albumin number here is an affinity, which is why the")
    print("  conversion matters and why it cannot be checked.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    rows, ncf = load()
    part1(rows)
    part2()
    part3(rows, ncf)
    part4(rows, ncf)
    part5(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
