#!/usr/bin/env python3
"""Verify the Li 2021 constants came from the EXPERIMENTAL column, not a docking one.

Li W, Hu Y, Bischel HN. Toxics 2021;9(3):63. doi:10.3390/toxics9030063
SI workbook toxics-1134374-SI.xlsx, sheet 'S7.2 Ka_docking and exp ' (the
sheet name ends in a space; openpyxl will not find it without one).

WHY THIS CHECK EXISTS. That sheet holds THREE log Ka columns side by side:

    B, C   1E7G                 GM / GSD log Ka   -- MOLECULAR DOCKING
    D, E   1AO6                 GM / GSD log Ka   -- MOLECULAR DOCKING
    F, G   Experimental result  log Ka / SD       -- the measurement

Only column F is a measurement. 1E7G and 1AO6 are PDB structures of albumin
and those columns are AutoDock output, which this project does not ingest as
data.

Taking the wrong column would not announce itself. For PFBA the docking value
is 4.099 against an experimental 4.097 -- a difference of 0.002, invisible in
any spot check. The columns separate only as the chain grows, because docking
rises monotonically with chain length while the measurement peaks around C7
and falls again (the paper's "inverted-V"). By PFTeDA it is 7.502 docked
against 3.962 measured: 3.5 log units, or a factor of 3000 in Kd.

So a spot check on a short-chain compound proves nothing here, and this script
checks all of them.

NOTE ON Kd. Li reports no dissociation constant anywhere -- "Kd" and
"dissociation constant" occur zero times in the paper. The Kd column in the
working sheet is derived: Kd/uM = 10^(6 - log Ka). This script does not
validate Kd against the paper, because there is nothing there to validate it
against.

Usage:  python3 -I scripts/verify_li2021_against_si.py <toxics-1134374-SI.xlsx>
"""
import csv, os, re, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.3390/toxics9030063'
SHEET = 'S7.2 Ka_docking and exp '      # trailing space is real
COL_EXP, COL_1E7G, COL_1AO6 = 5, 1, 3   # zero-based

def key(s):
    s = str(s).lower().strip()
    s = re.sub(r'\s*\((linear|branched)\)', '', s)
    s = re.sub(r'^(l-|br-|n-)', '', s).replace(' ', '')
    return {'pfuda': 'pfunda', 'pfdoa': 'pfdoda', 'fosa': 'pfosa'}.get(s, s)

def main(argv):
    if not argv:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    try:
        import openpyxl
    except ImportError:
        print("openpyxl required", file=sys.stderr)
        return 2
    ws = openpyxl.load_workbook(argv[0], data_only=True)[SHEET]

    si = {}
    for r in ws.iter_rows(min_row=3, values_only=True):
        if not r or not r[0]:
            continue
        si.setdefault(key(r[0]), []).append(
            (str(r[0]).strip(), r[COL_EXP], r[COL_1E7G], r[COL_1AO6]))

    db = [r for r in csv.DictReader(open(DB, encoding='utf-8'))
          if (r.get('doi') or '') == DOI]
    ok = bad = 0
    worst = (0.0, None)
    for r in sorted(db, key=lambda x: x['chemical']):
        k = key(r['chemical'])
        rows = si.get(k)
        if not rows:
            print(f"  NOT IN SI   {r['chemical']}"); bad += 1; continue
        # a key can hold both isomers; match on the value we stored
        ours = float(r['ka'])
        hit = next((x for x in rows if x[1] not in (None, 'NA')
                    and abs(float(x[1]) - ours) < 0.002), None)
        if hit is None:
            vals = [x[1] for x in rows]
            print(f"  MISMATCH    {r['chemical']:20s} ours {ours:.3f} "
                  f"vs SI exp {vals}"); bad += 1; continue
        ok += 1
        if hit[2] is not None and hit[2] != 'NA':
            gap = abs(float(hit[2]) - ours)
            if gap > worst[0]:
                worst = (gap, f"{r['chemical']} (docked {float(hit[2]):.3f} "
                              f"vs measured {ours:.3f})")

    print(f"matched the experimental column exactly: {ok}")
    print(f"mismatched or missing:                   {bad}")
    if worst[1]:
        print(f"largest docking-vs-measured gap: {worst[0]:.3f} log -- {worst[1]}")

    carried = {key(r['chemical']) for r in db}
    for k, rows in si.items():
        if k not in carried:
            for nm, e, _, _ in rows:
                print(f"in the SI, not carried: {nm} (exp log Ka = {e!r})")
    return 1 if bad else 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
