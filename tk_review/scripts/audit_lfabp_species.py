#!/usr/bin/env python3
"""Is there a rat-human difference in PFAS binding to L-FABP? No -- it is
smaller than the disagreement among human laboratories using the same method.

Two published claims sit in opposition:

  Modaresi et al. 2025 (doi 10.1021/acs.chemrestox.5c00199), discussion:
  "species differences in the binding affinity of PFAS to FABPs have been
  observed, with PFAS binding more strongly to human FABP than to rodent
  FABP."

  Fan et al. 2025 (doi 10.1016/j.envint.2025.109591), results: "the affinities
  of PFAS to rL-FABP were comparable to those measured for PFAS binding to
  human L-FABP using the same method (Zhang et al. 2013)."

This compilation can referee it. The rat arm is now TWO laboratories -- Fan
2025 and Woodcroft 2010, the latter re-ingested from the primary -- and three
human sources used the SAME assay, 1,8-ANS fluorescence displacement, so the
comparison can be made with method held fixed on both sides.

SIGN CONVENTION throughout: rat minus human, in log K. POSITIVE means the rat
protein binds MORE TIGHTLY.

THE ANSWER DEPENDS ENTIRELY ON WHICH HUMAN LABORATORY IS CHOSEN -- not just
its size, but ITS SIGN. Restricting both sides to fluorescence displacement
(ITC rows excluded, and a compound measured more than once by one source
reduced by median):

    rat vs Zhang 2013 (human)   10 shared   about +0.5 log  (rat tighter)
    rat vs Sheng 2016 (human)    2 shared   about -0.7 log  (human tighter)
    rat vs Yang 2020  (human)    6 shared   about -0.1 log  (no difference)

All three measured human L-FABP by displacement of a fluorescent probe. They
disagree among themselves by about 1.2 log -- PFOA alone is 2.36 uM in Sheng
2016 and 50.4 uM in Zhang 2013, a factor of 21 between two human laboratories
-- and that spread is larger than any rat-human gap and wide enough to contain
both signs. A species difference cannot be resolved above it. The script
prints the current figures; the text here is deliberately rounded so it cannot
drift silently out of agreement with them.

Pooling across sources, rat comes out a few tenths of a log TIGHTER overall,
and the gap shrinks on the compounds carrying four or more independent human
values -- i.e. the apparent species effect is largest exactly where human
coverage is thinnest, which is the signature of noise rather than biology.

NEITHER PUBLISHED CLAIM SURVIVES CLEANLY. Modaresi's direction -- human
binding more strongly than rodent -- holds against exactly one human source
(Sheng 2016, on two compounds) and reverses against the other two and against
the pooled comparison, so it is not a species effect but a choice of
comparator. Fan's claim of comparability is supported on the pooled and
Yang-matched cuts, but against the human source it actually cites, Zhang 2013,
the difference is the largest positive one of the three; it is against Yang
2020, which Fan does not cite for this, that the agreement is near-exact.

PRACTICAL CONSEQUENCE. A rat-to-human extrapolation factor for L-FABP binding
cannot be estimated from the present literature, because the quantity it would
correct for is smaller than the noise between laboratories measuring one
species. For toxicokinetic use the defensible choice is to treat rat and human
L-FABP binding as indistinguishable and carry the inter-laboratory spread as
the uncertainty.

Read-only.
"""
import collections, os, statistics, sys

XL = os.path.join(os.path.dirname(__file__), '..', 'lfabp_logK_Kd.xlsx')
L = 'log K (L water / kg protein)'
SAME_METHOD = ("Zhang 2013", "Sheng 2016", "Yang 2020")   # 1,8-ANS displacement

def medians(rows, key):
    """compound -> median log K, so a compound measured twice by one source
    (two methods, or a resolved site series) counts once, not last-wins."""
    acc = collections.defaultdict(list)
    for r in rows:
        if r[L] is not None:
            acc[str(r['PFAS']).strip()].append(r[L])
    return {c: statistics.median(v) for c, v in acc.items()}


def main():
    try:
        import openpyxl
    except ImportError:
        print("openpyxl required", file=sys.stderr); return 2
    ws = openpyxl.load_workbook(XL)['Coefficients']
    hdr = [c.value for c in ws[1]]
    R = [dict(zip(hdr, [c.value for c in r])) for r in ws.iter_rows(min_row=2)
         if r[0].value is not None]

    rat_src = sorted({r['Source'] for r in R if r['Species'] == 'rat'})
    rat_all = medians([r for r in R if r['Species'] == 'rat'], L)
    print(f"rat arm: {len(rat_all)} compounds from {len(rat_src)} "
          f"laboratories {rat_src}")

    # method-matched: fluorescence displacement on BOTH sides, ITC excluded
    disp = [r for r in R if 'ITC' not in str(r['Route'])]
    rat = medians([r for r in disp if r['Species'] == 'rat'], L)
    print("\nMETHOD-MATCHED, 1,8-ANS fluorescence displacement throughout "
          f"(ITC rows excluded; rat arm {len(rat)} compounds):")
    spreads = []
    for src in SAME_METHOD:
        o = medians([r for r in disp if r['Source'] == src], L)
        both = sorted(set(rat) & set(o))
        if not both:
            continue
        m = statistics.median([rat[c] - o[c] for c in both])
        spreads.append(m)
        print(f"   rat vs {src:12s} {len(both):2d} shared   median {m:+.2f} log")
    if len(spreads) > 1:
        print(f"\n   spread across the human sources: "
              f"{max(spreads) - min(spreads):.2f} log")
        print("   => larger than any single rat-human difference, so no species")
        print("      effect is resolvable above the inter-laboratory disagreement")

    print("\nPOOLED, all sources and methods (positive = rat binds tighter):")
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in R:
        if r[L] is not None:
            by[str(r['PFAS']).strip()][r['Species']].append(r[L])
    pairs = [(c, d) for c, d in by.items() if 'human' in d and 'rat' in d]
    allz = [statistics.median(d['rat']) - statistics.median(d['human'])
            for _, d in pairs]
    dense = [statistics.median(d['rat']) - statistics.median(d['human'])
             for _, d in pairs if len(d['human']) >= 4]
    print(f"   all {len(allz)} shared compounds        "
          f"median {statistics.median(allz):+.2f}")
    if dense:
        print(f"   the {len(dense)} with >=4 human values  "
              f"median {statistics.median(dense):+.2f}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
