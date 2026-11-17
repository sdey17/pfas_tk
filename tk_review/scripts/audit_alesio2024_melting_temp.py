#!/usr/bin/env python3
"""Why Alesio 2024 is not ingested: its constants are measured at BSA's melting
temperature, and the same lab's own 37 C data shows what that costs.

Alesio J, Bothun GD. "Differential scanning fluorimetry to assess PFAS binding
to bovine serum albumin protein." Sci Rep 2024;14:6501.
doi:10.1038/s41598-024-57140-9

THE DISQUALIFYING FACT, stated by the paper itself. The abstract: "binding
constants were determined by (a) an equilibrium-based model AT THE MELTING
TEMPERATURE OF BSA". Table 2's caption: "Association constants at 61 C". The
text: BSA's native Tm is 62 C, rising to 69-71 C with PFAS bound, and 61 C "was
selected because there were measurable differences between" the conditions.

So every Ka is evaluated where BSA is at its unfolding transition -- roughly
half-denatured. That is a defensible choice for the rapid screening the paper
is designed for, and the authors are entirely transparent about it. It is not
the quantity this sheet holds, which is binding to folded albumin near 37 C.
The sheet has already made that choice three times: Chen 2020, Qin 2010 and
Peng 2024 each measured at several temperatures and only the 310 K rows were
taken. Admitting a 61 C constant would contradict it.

WHAT MAKES THIS UNUSUALLY CLEAN. Alesio 2022 (doi 10.1021/acs.est.2c00500) is
the SAME LAB on the SAME SEVEN PFAS and the SAME PROTEIN, by equilibrium
dialysis at 37 C, and those 7 rows are already in the sheet. So the comparison
below is not cross-laboratory: it isolates the temperature and method change
with everything else held fixed.

    PFAS       2022 dialysis 37 C   2024 DSF 61 C      gap
    PFDA             4.967              1.696        -3.27
    PFOS             4.433              2.614        -1.82
    PFOA             3.754              2.565        -1.19
    PFHxS            3.989              3.151        -0.84
    PFNA             4.821              4.162        -0.66
    PFBS             3.146              2.495        -0.65
    HFPO-DA          2.990              2.474        -0.52

Every compound falls, by 0.5 to 3.3 log units, and the chain-length ordering
breaks: at 61 C PFNA (4.162) exceeds PFDA (1.696), reversing the C9 < C10
trend that the same lab's 37 C data gets right (4.821 < 4.967). A factor of
1860 on PFDA is not an offset that could be calibrated out.

Four of the 11 rows are not PFAS anyway -- octanoic, nonanoic and decanoic
acid as reference ligands, and warfarin as a Sudlow site I control -- so the
NON_PFAS rule would have taken those regardless.

Reviewer dropped the paper. This records the reason so it is not revisited.

Read-only.
"""
import csv, math, os, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
XL = os.path.join(os.path.dirname(__file__), '..', 'albumin_logK_Kd_revised.xlsx')
BRIDGE = math.log10(1e9 / 66463.0)
DOI24 = '10.1038/s41598-024-57140-9'

def main():
    rows = [r for r in csv.DictReader(open(DB, encoding='utf-8'))
            if DOI24 in (r.get('doi') or '')]
    print("Alesio 2024 -- why it stays out")
    print("=" * 66)
    print(f"  {len(rows)} database rows, species "
          f"{sorted({r.get('protein_species') for r in rows})}")
    print(f"  temperature: {(rows[0].get('temperature_c') or '')[:70]}")
    non = [r['chemical'] for r in rows
           if 'reference ligand' in r['chemical'] or 'warfarin' in r['chemical']]
    print(f"  of which not PFAS: {len(non)} -> {non}")

    try:
        import openpyxl
    except ImportError:
        print("\n(openpyxl unavailable; comparison skipped)"); return 0
    ws = openpyxl.load_workbook(XL)['Coefficients']
    hdr = [c.value for c in ws[1]]
    L = 'log K (L water / kg protein)'
    a22 = {}
    for r in ws.iter_rows(min_row=2):
        d = dict(zip(hdr, [c.value for c in r]))
        if d.get('Source') == 'Alesio 2022':
            a22[str(d['PFAS']).strip()] = d[L]

    print("\n  SAME LAB, SAME PROTEIN, SAME COMPOUNDS -- only T and method differ")
    print(f"  {'PFAS':10s}{'2022 dialysis 37C':>19s}{'2024 DSF 61C':>14s}{'gap':>8s}")
    gaps = []
    for r in rows:
        c = r['chemical'].split(' (')[0]
        if c not in a22 or not r.get('ka'):
            continue
        lk = math.log10(float(r['ka'])) - 6.0 + BRIDGE
        gaps.append(lk - a22[c])
        print(f"  {c:10s}{a22[c]:19.3f}{lk:14.3f}{lk - a22[c]:+8.2f}")
    if gaps:
        print(f"\n  every compound falls; gap {min(gaps):+.2f} to {max(gaps):+.2f} log")
        print("  => a 61 C constant at the unfolding transition is a different")
        print("     quantity, not an offset that calibrates out")
    return 0

if __name__ == '__main__':
    sys.exit(main())
