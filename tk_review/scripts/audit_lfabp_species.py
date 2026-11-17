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

This compilation can referee it, because Fan 2025 is the whole rat arm and
three human sources used the SAME assay -- 1,8-ANS fluorescence displacement --
so the comparison can be made with method held fixed.

THE ANSWER DEPENDS ENTIRELY ON WHICH HUMAN LABORATORY IS CHOSEN:

    Fan 2025 (rat) vs Zhang 2013 (human)   10 shared compounds   +0.50 log
    Fan 2025 (rat) vs Sheng 2016 (human)    3 shared compounds   +0.24 log
    Fan 2025 (rat) vs Yang 2020  (human)    6 shared compounds   -0.02 log

All three human sources measured human L-FABP by the same displacement assay
with the same probe. The spread between them (0.52 log) exceeds any of the
individual rat-human differences. So the between-laboratory difference among
human studies is larger than the between-species difference, and no species
effect can be resolved above it.

That also holds pooling across sources. Across all 12 compounds measured in
both species the median rat-human difference is +0.47 log; restricted to the
five compounds carrying four or more independent human values it is -0.09.
The apparent species effect is an artefact of thin human coverage on the
compounds where it appears.

NEITHER PUBLISHED CLAIM SURVIVES CLEANLY. Modaresi's direction (human stronger
than rodent) is not supported: where the data are densest the difference is
-0.09, i.e. rat marginally stronger, and against Yang 2020 it is -0.02. Fan's
claim of comparability is supported, but against the human source it cites --
Zhang 2013 -- the difference is the LARGEST of the three at +0.50; it is
against Yang 2020, which Fan does not cite for this, that the agreement is
near-exact.

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

def main():
    try:
        import openpyxl
    except ImportError:
        print("openpyxl required", file=sys.stderr); return 2
    ws = openpyxl.load_workbook(XL)['Coefficients']
    hdr = [c.value for c in ws[1]]
    R = [dict(zip(hdr, [c.value for c in r])) for r in ws.iter_rows(min_row=2)]

    rat = {str(r['PFAS']).strip(): r[L] for r in R if r['Species'] == 'rat'}
    print(f"rat arm: {len(rat)} compounds, all from "
          f"{sorted({r['Source'] for r in R if r['Species']=='rat'})}")

    print("\nMETHOD-MATCHED, 1,8-ANS fluorescence displacement throughout:")
    spreads = []
    for src in SAME_METHOD:
        o = {str(r['PFAS']).strip(): r[L] for r in R if r['Source'] == src}
        both = sorted(set(rat) & set(o))
        if not both:
            continue
        d = [rat[c] - o[c] for c in both]
        m = statistics.median(d)
        spreads.append(m)
        print(f"   rat vs {src:12s} {len(both):2d} shared   median {m:+.2f} log")
    if len(spreads) > 1:
        print(f"\n   spread across the three human sources: "
              f"{max(spreads) - min(spreads):.2f} log")
        print("   => larger than any single rat-human difference, so no species")
        print("      effect is resolvable above the inter-laboratory disagreement")

    print("\nPOOLED, all sources:")
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in R:
        by[str(r['PFAS']).strip()][r['Species']].append(r[L])
    pairs = [(c, d) for c, d in by.items() if 'human' in d and 'rat' in d]
    allz = [statistics.median(d['rat']) - statistics.median(d['human'])
            for c, d in pairs]
    dense = [statistics.median(d['rat']) - statistics.median(d['human'])
             for c, d in pairs if len(d['human']) >= 4]
    print(f"   all {len(allz)} shared compounds        median {statistics.median(allz):+.2f}")
    print(f"   the {len(dense)} with >=4 human values  median {statistics.median(dense):+.2f}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
