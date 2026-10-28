#!/usr/bin/env python3
"""Does a partition coefficient CONVERTED from a Kd agree with one that was
MEASURED directly?

build_binding_affinity_table.py now fills `logK` by two quite different
routes, and the column does not distinguish them:

    MEASURED   the paper reports a protein-water partition coefficient.
               187 rows (plus 56 biomimetic ones whose basis the paper
               never states, held apart throughout).
    CONVERTED  the paper reports an affinity and a stoichiometry, and
               log K = log10(1e9/M) + log n - log Kd(uM). 82 rows.

This file asks whether the two agree. They do not, and the answer is
more specific than "the conversion is unreliable".

THE ANSWER: 33 pairs are available, Chen et al. 2025 supplies 22 of
them, and that one well-powered comparison sits 0.80 log LOW with a
consistent sign, outside the 0.59 log between-laboratory baseline. The
size of that gap is the size of a missing site count: 10^0.80 = 6.3, so
closing it needs about six sites where Chen's single-site Langmuir fit
allowed one. Maso et al. 2021 co-crystallised four PFOA molecules per
human serum albumin and found one high-affinity and three low-affinity
sites in solution -- independent, structural, and of the right order.

    source        pairs   median   where its n came from
    Chen 2025        22    -0.80   single-site Langmuir, n fixed at 1
    Wu 2024           4    -4.18   fluorescence double-log slope
    Chen 2009         2    -0.53   fluorescence double-log slope
    Peng 2024         3    +0.13   fluorescence double-log slope
    Moro 2022         2    +0.28   ITC, n fitted from the isotherm

The remaining 11 pairs are too few to characterise and should not be
read as a trend. They do establish one negative: three papers that take
n off a fluorescence double-logarithmic slope land at -4.18, -0.53 and
+0.13, so that slope is not a site count. Wu 2024's problem is not even
its n -- its association constants sit three to four orders below every
other measurement of the same compound on albumin, which this database
already flags.

Note what CANNOT be separated here: Chen 2025 is the only assumed-n
source in the comparison, so "assumed n" and "Chen 2025" are the same
column of numbers. The 0.80 log is solid; attributing it to the
assumption rather than to the laboratory needs a second assumed-n paper
measured against the same set, and there is none.

I TESTED THE OBVIOUS ALTERNATIVE AND IT FAILS. The natural story is that
a Kd describes one saturable site while a measured partition coefficient
also counts non-specific surface, so the converted value should be low by
the non-specific share. Qin et al. 2024 Table 2 reports those two terms
SEPARATELY for bovine serum albumin -- the only such split in this
database -- and at the trace limit the non-specific term contributes
1 to 14 per cent of the total, which is 0.01 to 0.07 log. That is one to
two orders too small to be the explanation. The story is wrong; the
arithmetic that killed it is in part 4.

Three further things keep this honest:

  1. A BASELINE. Two DIRECT measurements of the same compound, protein and
     species in different laboratories disagree by a median 0.59 log over
     46 pairs, and two batches of albumin in ONE laboratory by 0.19 log.
     Any claim about the conversion has to clear those bars.
  2. AN INTERNAL CHECK on Qin's split: specific + non-specific should
     recover a direct total, and Fischer et al. 2024 measured one on the
     same protein and species.
  3. A PREDICTION BASELINE. A chain-length regression fitted inside one
     study, scored leave-one-out, says how well this quantity can be
     predicted at all when no conversion is involved.

Qin's specific column is Eq. 15, K_spec = Cmax/(2 Kd), the ratio AT
half-saturation. That is 0.301 log below the Cw -> 0 plateau that a
converted Kd gives, so +0.301 is applied wherever the two are compared.
It is an exact correction with the paper's own equation behind it, not an
adjustment.

NOTHING HERE POOLS ACROSS STUDIES. Every comparison names both studies and
both methods, and no value from one paper is ever used to correct a value
in another. The bovine decomposition says what SIZE of gap to expect in
the human comparison; it does not touch a human number.

Writes db/logk_converted_vs_measured.csv.

Run:  python3 tk_review/scripts/logk_converted_vs_measured.py
"""
import collections
import csv
import importlib.util
import itertools
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "logk_converted_vs_measured.csv")

_spec = importlib.util.spec_from_file_location(
    "baf", os.path.join(HERE, "build_binding_affinity_table.py"))
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)

HALF_SAT = math.log10(2.0)      # Qin Eq. 15 is at half-saturation

PFCA = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
        "PFUnDA", "PFDoDA", "PFTrDA", "PFTeDA"]
PFSA = ["PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS"]


def rule(ch="=", n=78):
    print(ch * n)


def ncf_map():
    db = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
    out = {}
    for r in csv.DictReader(open(db, encoding="utf-8")):
        v = r["n_perfluorinated_carbons"].strip()
        if v.isdigit():
            out.setdefault(r["chemical"], int(v))
    return out


def head(pfas):
    if pfas in PFCA:
        return "carboxylate"
    if pfas in PFSA:
        return "sulfonate"
    return "other"


def protein_class(r):
    p = r["protein"].lower()
    if "globulin" in p:
        return "globulin"
    if "serum albumin" in p:
        return "albumin"
    if "l-fabp" in p or "liver fatty acid" in p:
        return "L-FABP"
    return p[:26]


def route(r):
    b = r["logK_basis"]
    if b.startswith("from the paper's own Kd"):
        return "converted"
    if b.startswith("as reported (mass") or b.startswith("reported per L"):
        return "measured"
    if "biomimetic" in b:
        return "biomimetic"
    return None


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0:
        return 0.0, my, 0.0
    m = sxy / sxx
    syy = sum((y - my) ** 2 for y in ys)
    r = sxy / math.sqrt(sxx * syy) if syy else 0.0
    return m, my - m * mx, r


def describe(vals, label):
    v = sorted(vals)
    print(f"    {label:28s} n {len(v):3d}   median {statistics.median(v):+.2f}"
          f"   mean {statistics.fmean(v):+.2f}"
          f"   range {v[0]:+.2f} to {v[-1]:+.2f}")


def part1(out):
    rule()
    print("1. THE TWO ROUTES INTO ONE COLUMN")
    rule()
    byroute = collections.Counter(route(r) for r in out if r["logK"])
    for k in ("measured", "converted", "biomimetic"):
        print(f"    {k:12s}{byroute[k]:5d} rows")
    print()
    print("  CONVERTED, by source and protein:")
    agg = collections.Counter((r["source"], protein_class(r), r["species"])
                              for r in out if route(r) == "converted")
    for (s, p, sp), v in agg.most_common():
        print(f"    {v:4d}  {s:18s}{p:10s}{sp}")
    print()
    print("  MEASURED, by source and protein:")
    agg = collections.Counter((r["source"], protein_class(r), r["species"])
                              for r in out if route(r) == "measured")
    for (s, p, sp), v in agg.most_common():
        print(f"    {v:4d}  {s:18s}{p:12s}{sp}")
    print()


def key(r):
    return (r["pfas"], protein_class(r), r["species"])


def totals(out):
    """Measured rows, one TOTAL per source x compound x protein x species.

    Every source but one reports a single coefficient. Qin et al. 2024
    reports a specific and a non-specific term for albumin, and neither
    is comparable with another laboratory's single number: they are
    components. They are added in linear space here, with the specific
    term first lifted by log10(2) off Eq. 15's half-saturation basis
    onto the Cw -> 0 plateau, so the sum is a total on the same footing
    as everyone else's.
    """
    out_rows, qin = [], collections.defaultdict(dict)
    for r in out:
        if route(r) != "measured":
            continue
        if r["source"] == "Qin 2024 (plasma)" and "column" in r["location"]:
            which = ("nonspecific" if "nonspecific column" in r["location"]
                     else "specific")
            qin[key(r)][which] = float(r["logK"])
            qin[key(r)]["row"] = r
            continue
        out_rows.append(r)
    for k, d in qin.items():
        if "specific" not in d or "nonspecific" not in d:
            continue                       # a lone component is not a total
        tot = math.log10(10 ** (d["specific"] + HALF_SAT)
                         + 10 ** d["nonspecific"])
        out_rows.append(dict(d["row"], logK="%.3f" % tot,
                             logK_basis="specific (plateau) + non-specific, "
                                        "summed from the paper's two columns"))
    return out_rows


def part2(out):
    """How far apart are two DIRECT measurements? The bar to clear."""
    rule()
    print("2. THE BASELINE -- TWO DIRECT MEASUREMENTS OF THE SAME THING")
    rule()
    print("  Before asking whether a converted value agrees with a measured")
    print("  one, it is worth knowing how well two MEASURED ones agree. Same")
    print("  compound, same protein, same species, different laboratory.")
    print()
    print("  Qin et al. 2024 reports its albumin coefficient as TWO")
    print("  columns, a specific and a non-specific term. Neither is a")
    print("  total, so comparing either against another laboratory's")
    print("  single number would compare a part with a whole. They are")
    print("  combined into one total per compound here, and that total is")
    print("  what enters every pair below.")
    print()
    # A source can measure the same compound twice. Allendorf et al. 2019
    # ran two BSA batches, one essentially fatty-acid free and one
    # heat-shock. Those replicates are a floor on the whole comparison and
    # are reported as one, then averaged so that a lab with replicates
    # does not count twice in every between-lab pair.
    reps = collections.defaultdict(list)
    for r in totals(out):
        reps[(r["source"],) + key(r)].append(float(r["logK"]))
    within = [max(v) - min(v) for v in reps.values() if len(v) > 1]
    if within:
        src = {k[0] for k, v in reps.items() if len(v) > 1}
        print(f"  WITHIN one laboratory, same compound, repeated "
              f"({', '.join(sorted(src))}):")
        print(f"    {len(within)} replicate pairs, median spread "
              f"{statistics.median(within):.2f} log, worst "
              f"{max(within):.2f} log.")
        print("    That is the floor: the same people, the same protocol,")
        print("    two preparations of the protein.")
        print()
    meas = collections.defaultdict(list)
    for (src, *k), v in reps.items():
        meas[tuple(k)].append((src, statistics.fmean(v)))
    pairs = collections.defaultdict(list)
    for k, v in meas.items():
        for (sa, va), (sb, vb) in itertools.combinations(sorted(v), 2):
            if sa == sb:
                continue
            pairs[(sa, sb)].append((k[0], va - vb))
    allgaps = []
    print(f"  {'study pair':44s}{'n':>4s}{'median |gap|':>14s}{'worst':>8s}")
    for (a, b), v in sorted(pairs.items(), key=lambda kv: -len(kv[1])):
        g = [abs(d) for _, d in v]
        allgaps += g
        worst = max(v, key=lambda t: abs(t[1]))
        print(f"  {a + ' vs ' + b:44s}{len(v):>4d}{statistics.median(g):>14.2f}"
              f"{abs(worst[1]):>8.2f}  {worst[0]}")
    print()
    med = statistics.median(allgaps)
    print(f"  Over all {len(allgaps)} pairs: median |gap| {med:.2f} log "
          f"({10 ** med:.1f}x), worst {max(allgaps):.2f} log.")
    print(f"  Direct measurement of this quantity is reproducible between")
    print(f"  laboratories to about a factor of {10 ** med:.1f}. That is the bar")
    print("  anything else has to clear.")
    print()
    return statistics.median(allgaps)


def n_used(rec):
    """The stoichiometry that went into this row's conversion."""
    b = rec["logK_basis"]
    if "one class of site" in b:
        return 1.0
    if "n = " in b:
        return float(b.split("n = ")[1].split(",")[0].strip())
    if "summed over" in b:
        return float(b.split("summed over")[1].split()[0])
    return 1.0


# Where each paper's stoichiometry comes from. This is read off the
# `method` field of the source rows, not inferred, and it is kept as a
# literal table because it is the axis the residuals sort on.
N_PROVENANCE = {
    "Chen 2025":      ("assumed",
                       "single-site Langmuir; n fixed at 1"),
    "Moro 2022":      ("isotherm",
                       "ITC, n fitted from the isotherm"),
    "Woodcroft 2010": ("isotherm",
                       "ITC and displacement, several site classes"),
    "Peng 2024":      ("double-log",
                       "fluorescence double-log slope"),
    "Chen 2009":      ("double-log",
                       "fluorescence quenching double-log slope"),
    "Wu 2024":        ("double-log",
                       "fluorescence double-log, no inner-filter correction"),
    "Khazaee 2021":   ("assumed",
                       "single-site model assumed in the fit"),
}


def n_how(rec):
    return N_PROVENANCE.get(rec["source"], ("unstated", "unstated"))[0]


def n_how_long(rec):
    return N_PROVENANCE.get(rec["source"], ("unstated", "unstated"))[1]


def part3(out, ncf, bar):
    """Converted against measured, same compound, protein and species."""
    rule()
    print("3. CONVERTED AGAINST MEASURED")
    rule()
    conv = collections.defaultdict(list)
    meas = collections.defaultdict(list)
    for r in out:
        if route(r) == "converted":
            conv[key(r)].append(r)
        elif route(r) == "measured":
            meas[key(r)].append(r)
    rows = []
    for k in sorted(set(conv) & set(meas)):
        for c, m in itertools.product(conv[k], meas[k]):
            rows.append({
                "pfas": k[0], "protein": k[1], "species": k[2],
                "nCF": ncf.get(k[0], ""), "head": head(k[0]),
                "converted_source": c["source"], "converted_logK": float(c["logK"]),
                "converted_Kd_uM": c["kd_uM"] or "",
                "n_used": n_used(c), "n_how": n_how(c), "n_long": n_how_long(c),
                "measured_source": m["source"], "measured_logK": float(m["logK"]),
                "residual": float(c["logK"]) - float(m["logK"]),
            })
    if not rows:
        print("  no pairs")
        return rows
    print(f"  {len(rows)} pairs over "
          f"{len({(r['pfas'], r['protein'], r['species']) for r in rows})} "
          f"(compound, protein, species) keys.")
    print(f"  Converted from: "
          f"{', '.join(sorted({r['converted_source'] for r in rows}))}")
    print(f"  Measured by:    "
          f"{', '.join(sorted({r['measured_source'] for r in rows}))}")
    print()
    res = [r["residual"] for r in rows]
    neg = sum(1 for x in res if x < 0)
    describe(res, "converted - measured, all")
    print(f"    {neg} of {len(res)} pairs are negative, "
          f"{len(res) - neg} positive.")
    print()
    print("  Pooling those is the wrong move, and the by-source breakdown")
    print("  says why: the sources do not disagree by degree, they fall")
    print("  into three groups that behave differently.")
    print()
    print(f"  {'source':16s}{'n':>3s}{'median':>9s}{'where n came from':>24s}")
    for s2 in sorted({r["converted_source"] for r in rows}):
        sub = [r for r in rows if r["converted_source"] == s2]
        med = statistics.median([r["residual"] for r in sub])
        print(f"  {s2:16s}{len(sub):>3d}{med:>+9.2f}   {sub[0]['n_long']}")
    print()
    print("  GROUPED BY WHERE THE STOICHIOMETRY CAME FROM:")
    groups = {}
    for g in ("isotherm", "double-log", "assumed"):
        v = [r["residual"] for r in rows if r["n_how"] == g]
        if v:
            describe(v, g)
            groups[g] = v
    print()
    print("  Read that honestly. Only ONE source is well powered: Chen 2025")
    print(f"  supplies {sum(1 for r in rows if r['converted_source'] == 'Chen 2025')}"
          f" of the {len(rows)} pairs, and it is the only assumed-n source")
    print("  here, so 'assumed' and 'Chen 2025' are the same column of")
    print("  numbers and cannot be separated by this dataset.")
    print()
    if "double-log" in groups:
        dl = sorted({r["converted_source"] for r in rows
                     if r["n_how"] == "double-log"})
        meds = {d: statistics.median([r["residual"] for r in rows
                                      if r["converted_source"] == d])
                for d in dl}
        print("  The double-log group is not a group at all: "
              + ", ".join(f"{d} {meds[d]:+.2f}" for d in dl) + ".")
        print("  A fluorescence double-log slope is not a site count, and")
        print("  these three show it. Wu 2024's problem is not even its n:")
        print("  its association constants sit three to four orders below")
        print("  every other measurement of the same compound on albumin,")
        print("  which this database already flags.")
    print()
    chen = [r["residual"] for r in rows if r["converted_source"] == "Chen 2025"]
    if chen:
        print(f"  WHAT CAN BE SAID. Chen 2025, {len(chen)} pairs against one")
        print(f"  measured set, sits at {statistics.median(chen):+.2f} log, "
              f"outside the {bar:.2f} log")
        print("  between-laboratory baseline, with a consistent sign. The")
        print("  other four sources contribute 11 pairs between them and")
        print("  span +0.28 to -4.18; nothing can be concluded from them")
        print("  except that a double-log n is not to be trusted.")
    print()
    print("  THE GAP IS A SITE COUNT. Rearranging, the stoichiometry that")
    print("  would reconcile each pair is n_implied = n_used x 10^(-residual).")
    print()
    print(f"  {'source':16s}{'n used':>8s}{'n implied (median)':>21s}")
    for s2 in sorted({r["converted_source"] for r in rows}):
        sub = [r for r in rows if r["converted_source"] == s2]
        imp = statistics.median([r["n_used"] * 10 ** (-r["residual"])
                                 for r in sub])
        print(f"  {s2:16s}{statistics.median([r['n_used'] for r in sub]):>8.2f}"
              f"{imp:>21.1f}")
    print()
    print("  Chen 2025 needs about six sites where its Langmuir fit allowed")
    print("  one. Maso et al. 2021 co-crystallised FOUR PFOA per human")
    print("  serum albumin and reported one high-affinity plus three")
    print("  low-affinity sites in solution. That is independent structural")
    print("  evidence of the right order, from a third laboratory, and it")
    print("  is not used to correct anything here -- only to say that a")
    print("  single-site fit of albumin is known to be an underestimate.")
    print()
    xs = [r["nCF"] for r in rows if r["nCF"] != ""]
    ys = [r["residual"] for r in rows if r["nCF"] != ""]
    m, b0, rr = fit(xs, ys)
    print(f"  RESIDUAL AGAINST CHAIN LENGTH: {m:+.3f} log per CF2 "
          f"(r {rr:+.2f}, n = {len(xs)}) -- weak, and it does not hold")
    print("  within head group:")
    for h in ("carboxylate", "sulfonate"):
        sub = [(r["nCF"], r["residual"]) for r in rows if r["head"] == h]
        if len(sub) > 2:
            m2, _, r2 = fit([x for x, _ in sub], [y for _, y in sub])
            print(f"    {h:12s} {m2:+.3f} per CF2 (r {r2:+.2f}, n = {len(sub)})")
    print("  Carboxylates and sulfonates give opposite signs, so there is")
    print("  no chain-length story here to tell. An earlier draft of this")
    print("  file claimed one off the pooled slope; the split kills it.")
    print()
    return rows


def qin_split(out):
    """Qin 2024 Table 2 bovine albumin, specific and non-specific, per kg.

    The specific column is Eq. 15, the ratio at half-saturation, so +0.301
    puts it on the Cw -> 0 plateau that a converted Kd gives.
    """
    spec, nons = {}, {}
    for r in out:
        if r["source"] != "Qin 2024 (plasma)" or protein_class(r) != "albumin":
            continue
        if not r["logK"]:
            continue
        if "specific column" in r["location"] and "nonspecific" not in r["location"]:
            spec[r["pfas"]] = float(r["logK"]) + HALF_SAT
        elif "nonspecific column" in r["location"]:
            nons[r["pfas"]] = float(r["logK"])
    return spec, nons


def part4(out, ncf):
    """The missing term, read off the only paper that reports it."""
    rule()
    print("4. THE OBVIOUS EXPLANATION, TESTED AND REJECTED")
    rule()
    print("  The natural story for part 3 is non-specific binding. A Kd")
    print("  describes one saturable site; a measured partition")
    print("  coefficient also counts whatever sticks to the rest of the")
    print("  protein. On that story the converted value must come out low")
    print("  by the non-specific share.")
    print()
    print("  It is testable. Qin et al. 2024 Table 2 is the only place in")
    print("  this database where the two contributions are reported")
    print("  SEPARATELY for one protein: bovine serum albumin, a saturable")
    print("  specific term and a linear non-specific one, fitted to the")
    print("  same isotherm.")
    print()
    print("  The specific column is Eq. 15, K_spec = Cmax/(2 Kd), the ratio")
    print("  at half-saturation. A converted Kd gives the Cw -> 0 plateau,")
    print(f"  which is {HALF_SAT:.3f} log higher, so that is added here.")
    print()
    spec, nons = qin_split(out)
    shared = sorted(set(spec) & set(nons), key=lambda c: ncf.get(c, 99))
    print(f"  {'cmpd':15s}{'nCF':>4s}{'spec(plateau)':>14s}{'nonspec':>9s}"
          f"{'total':>8s}{'spec share':>12s}{'total-spec':>11s}")
    shares, gaps, xs = [], [], []
    for c in shared:
        s, nsp = spec[c], nons[c]
        tot = math.log10(10 ** s + 10 ** nsp)
        share = 10 ** s / (10 ** s + 10 ** nsp)
        print(f"  {c:15s}{ncf.get(c, ''):>4}{s:>14.2f}{nsp:>9.2f}{tot:>8.2f}"
              f"{100 * share:>11.0f}%{s - tot:>11.2f}")
        shares.append(share)
        gaps.append(s - tot)
        xs.append(ncf.get(c, 0))
    print()
    print("  The last column is the quantity part 3 would have to")
    print("  explain: how far below the total a specific-only number sits.")
    describe(gaps, "specific - total, here")
    m, _, r = fit(xs, gaps)
    print(f"    against chain length         {m:+.3f} log per CF2 (r {r:+.2f})")
    print()
    big = max(abs(g) for g in gaps)
    print("  THE STORY FAILS. At the trace limit the specific term carries")
    print(f"  {100 * min(shares):.0f} to {100 * max(shares):.0f} per cent "
          f"of the total, so dropping the non-specific")
    print(f"  term costs at most {big:.2f} log. Part 3's gap is 0.80. "
          f"Non-specific")
    print(f"  binding is {0.80 / big:.0f}x too small to be the explanation, and the")
    print("  sign of its chain-length trend is wrong too.")
    print()
    print("  That is not a surprise once stated: the specific term is")
    print("  saturable and therefore LARGEST at trace, which is exactly")
    print("  the limit a converted Kd describes. Non-specific partitioning")
    print("  takes over higher up the isotherm, where neither number")
    print("  applies. The missing capacity in part 3 is extra SITES, not")
    print("  extra surface.")
    print()
    return spec, nons


def part5(out, ncf):
    """Does specific + non-specific recover a direct total?"""
    rule()
    print("5. THE SPLIT CHECKED AGAINST A DIRECT TOTAL")
    rule()
    print("  If the decomposition is right, specific + non-specific should")
    print("  reproduce a directly measured total. Fischer et al. 2024")
    print("  measured one on the same protein and the same species, so the")
    print("  check is available and stays inside bovine albumin.")
    print()
    print("  Both laboratories use solid-phase microextraction. This is")
    print("  therefore not a comparison of techniques: it is two")
    print("  implementations of the SAME technique on the same protein.")
    print()
    spec, nons = qin_split(out)
    fis = {r["pfas"]: float(r["logK"]) for r in out
           if r["source"] == "Fischer 2024" and protein_class(r) == "albumin"
           and r["species"] == "bovine" and r["logK"]}
    shared = sorted(set(spec) & set(nons) & set(fis),
                    key=lambda c: ncf.get(c, 99))
    if not shared:
        print("  no shared compounds")
        return []
    print(f"  {'cmpd':15s}{'nCF':>4s}{'Qin spec+nonspec':>18s}"
          f"{'Fischer total':>15s}{'diff':>8s}{'spec alone':>12s}{'vs Fischer':>12s}")
    d_sum, d_spec, rows = [], [], []
    for c in shared:
        tot = math.log10(10 ** spec[c] + 10 ** nons[c])
        print(f"  {c:15s}{ncf.get(c, ''):>4}{tot:>18.2f}{fis[c]:>15.2f}"
              f"{tot - fis[c]:>8.2f}{spec[c]:>12.2f}{spec[c] - fis[c]:>12.2f}")
        d_sum.append(tot - fis[c])
        d_spec.append(spec[c] - fis[c])
        rows.append((c, ncf.get(c, ""), tot, fis[c], spec[c]))
    print()
    describe(d_sum, "sum vs Fischer total")
    describe(d_spec, "specific alone vs Fischer")
    a = statistics.median([abs(x) for x in d_sum])
    b = statistics.median([abs(x) for x in d_spec])
    print()
    print(f"  Adding the non-specific term moves the agreement from "
          f"{b:.2f} to {a:.2f} log.")
    print("  It changes nothing, which is part 4's result stated a second")
    print("  way: at trace the specific term IS the coefficient.")
    print()
    print(f"  What is left is a flat {statistics.median(d_sum):+.2f} log "
          f"offset between the two")
    print("  laboratories, the same sign on all five compounds and almost")
    print("  the same size on each. Two SPME implementations of the same")
    print("  measurement on the same protein differ by a factor of")
    print(f"  {10 ** abs(statistics.median(d_sum)):.1f}. That is the part-2 "
          f"baseline showing up in a single pair,")
    print("  and it is a floor under every comparison in this file: no")
    print("  conversion can be validated to better than this.")
    print()
    print("  It also means part 3's -0.80 is only modestly outside what")
    print("  two laboratories do to each other. The reason to believe it")
    print("  is a real missing term is not its size but its CONSISTENCY")
    print("  across 22 compounds and the independent structural evidence")
    print("  for four or more PFOA sites on albumin -- not this margin.")
    print()
    return rows


def part6(out, ncf):
    """How well can this quantity be predicted with no conversion at all?"""
    rule()
    print("6. A PREDICTION BASELINE WITH NO CONVERSION IN IT")
    rule()
    print("  Part 3 puts Chen's converted values 0.80 log below Fischer's")
    print("  measured ones. Is that bad? Only against an alternative. The")
    print("  alternative that needs no affinity data at all is a")
    print("  chain-length regression fitted INSIDE one study, one head")
    print("  group and one protein, scored leave-one-out so that no")
    print("  compound helps predict itself.")
    print()
    print(f"  {'study / protein / species / class':46s}{'n':>3s}"
          f"{'slope':>8s}{'r':>7s}{'LOO median |err|':>18s}")
    best = []
    groups = collections.defaultdict(list)
    for r in out:
        if route(r) != "measured":
            continue
        n = ncf.get(r["pfas"])
        h = head(r["pfas"])
        if n is None or h == "other":
            continue
        groups[(r["source"], protein_class(r), r["species"], h)].append(
            (n, float(r["logK"]), r["pfas"]))
    for g, v in sorted(groups.items()):
        if len(v) < 5:
            continue
        xs = [a for a, _, _ in v]
        ys = [b for _, b, _ in v]
        m, b0, rr = fit(xs, ys)
        errs = []
        for i in range(len(v)):
            ox = xs[:i] + xs[i + 1:]
            oy = ys[:i] + ys[i + 1:]
            mm, bb, _ = fit(ox, oy)
            errs.append(abs(mm * xs[i] + bb - ys[i]))
        med = statistics.median(errs)
        best.append(med)
        label = f"{g[0]} / {g[1]} / {g[2]} / {g[3]}"
        print(f"  {label:46s}{len(v):>3d}{m:>8.3f}{rr:>7.2f}{med:>18.2f}")
    print()
    if best:
        print(f"  Median across groups: {statistics.median(best):.2f} log "
              f"({10 ** statistics.median(best):.1f}x).")
    print("  Put that beside part 3 before reading too much into either.")
    print("  A within-study chain-length line interpolates to about the")
    print("  same accuracy as the between-laboratory spread, and the")
    print("  conversion's offset is under twice that -- close, not an")
    print("  order apart. Neither route is good; the regression is")
    print("  somewhat less bad, and it is unbiased by construction where")
    print("  the conversion is one-sided.")
    print()
    print("  They also answer different questions, which matters more")
    print("  than the margin. 'What is log K for a compound this study")
    print("  did not measure' is the regression's question and it needs")
    print("  no affinity data. 'What does the saturable site contribute,")
    print("  and at what free concentration does it stop contributing' is")
    print("  the conversion's, and no regression can answer it, because a")
    print("  line through log K against chain length contains no Kd and")
    print("  therefore knows nothing about saturation.")
    print()


def part7(out):
    rule()
    print("7. WHAT COULD NOT BE PUT ON THIS PLOT")
    rule()
    conv = {key(r) for r in out if route(r) == "converted"}
    meas = {key(r) for r in out if route(r) == "measured"}
    print(f"  {len(conv - meas)} converted keys have no measured partner:")
    for (p, sp), v in collections.Counter(
            (k[1], k[2]) for k in conv - meas).most_common():
        print(f"    {v:4d}  {p} / {sp}")
    print()
    print("  The L-FABP rows are the whole of the non-albumin conversion")
    print("  and there is no directly measured L-FABP partition")
    print("  coefficient anywhere in the database to test them against.")
    print("  That gap is the most useful thing this analysis found: the")
    print("  protein that matters most for liver accumulation is the one")
    print("  with no measured protein-water coefficient at all.")
    print()
    print("  Also absent, and for a reason:")
    print("    - 56 biomimetic rows (Muensterman 2026). The paper does not")
    print("      state the basis of its log K_HSA, so it cannot be placed")
    print("      on a per-kg axis without assuming one.")
    print("    - 7 PPAR rows. Affinity and stoichiometry both present, the")
    print("      construct mass absent.")
    print("    - 42 Gao 2025 rows, whose n comes from reading site counts")
    print("      off Hill coefficients in prose.")
    print("    - every Freundlich row, which is a ratio at one stated")
    print("      aqueous concentration and not a constant.")
    print()


def write(rows):
    if not rows:
        return
    fields = list(rows[0].keys())
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {os.path.abspath(OUT)}  ({len(rows)} pairs)")


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    rows_db, out, skipped, qual, dropped = B.build()
    ncf = ncf_map()
    part1(out)
    bar = part2(out)
    pairs = part3(out, ncf, bar)
    part4(out, ncf)
    part5(out, ncf)
    part6(out, ncf)
    part7(out)
    write(pairs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
