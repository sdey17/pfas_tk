#!/usr/bin/env python3
"""Correct the Muensterman 2026 basis note and record the published Correction.

Two things were wrong or missing in db/protein_binding.csv:

  1. The note read "basis not stated in the paper, so not converted". The SI
     DOES state a basis -- the values are described as dimensionless partition
     coefficients. What is missing is not the basis but the phase ratio and the
     immobilised-albumin mass, without which a dimensionless ratio on an
     albumin-on-silica stationary phase cannot be put on a per-kg-albumin
     footing. The note now says that, because "basis not stated" would invite a
     future reader to go looking for a basis statement that is already there.

  2. The rows carried no reference to the published Correction
     (Environ Sci Technol 2026;60(23):17051, doi:10.1021/acs.est.6c07668),
     which fixes eq 4's intercept from "42" to "0.42". Verified in
     scripts/audit_muensterman2026_biomimetic.py: the published values were
     computed with 0.42 throughout, so no reported number changes, and eq 4
     feeds the IAM chain only -- the HSA rows never touch it.

Idempotent: re-running makes no further change.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.1021/acs.est.5c11716'
CORR = 'Correction: Environ Sci Technol 2026;60(23):17051, doi:10.1021/acs.est.6c07668 (eq 4 intercept 42 -> 0.42; affects the IAM chain only, verified not to change any published value)'

NEW_HSA = ('log K_HSA (albumin-water, biomimetic); retention-derived index from '
           'eq 6 log K_HSA = e^(log k_HSA) on an immobilised-HSA column. The SI '
           'states the values are dimensionless partition coefficients, but no '
           'phase ratio or immobilised-albumin mass is reported anywhere in the '
           'paper or SI, so there is no route to L water / kg albumin: the '
           'albumin quantity is absorbed into the calibration constants. '
           'Calibration rests on equilibrium-dialysis percent-bound data, which '
           'cannot resolve above 99.5% bound; the authors state that absolute '
           'values may not be confirmable and that the elution order is what is '
           'reproducible. Not converted. ' + CORR)

NEW_IAM = ('log K_IAM (phospholipid membrane-water, biomimetic); retention-derived '
           'index from eq 4 (corrected intercept 0.42) feeding eq 5 '
           'log K_IAM = 0.29 e^(log k_IAM) + 0.70 on an IAM.PC.DD2 column. A '
           'phospholipid stationary phase, not a protein, so out of scope for '
           'the albumin coefficient sheet regardless of basis; and as for HSA, no '
           'phase ratio is reported. Not converted. ' + CORR)

OLD_PREFIXES = ('log K_HSA (albumin-water, biomimetic); basis not stated',
                'log K_IAM (phospholipid membrane-water, biomimetic); basis not stated')

def main():
    with open(DB, newline='', encoding='utf-8') as fh:
        rdr = csv.DictReader(fh)
        fields, rows = rdr.fieldnames, list(rdr)

    touched = already = 0
    for r in rows:
        if (r.get('doi') or '') != DOI:
            continue
        is_hsa = 'HSA' in (r.get('protein') or '')
        want = NEW_HSA if is_hsa else NEW_IAM
        cur = r.get('ka_units') or ''
        if cur == want:
            already += 1
        elif cur.startswith(OLD_PREFIXES):
            r['ka_units'] = want
            touched += 1
        elif cur:
            print(f"  ! unexpected note on {r.get('chemical')!r}, left alone:\n"
                  f"    {cur[:90]}", file=sys.stderr)

    if not touched:
        print(f"no change needed ({already} rows already carry the corrected note)")
        return 0

    shutil.copy2(DB, DB + '.bak')
    with open(DB, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"updated {touched} Muensterman rows ({already} were already correct)")
    print("  - replaced 'basis not stated' with the real obstruction "
          "(no phase ratio / albumin mass)")
    print("  - recorded the Correction and its verified scope")
    return 0

if __name__ == '__main__':
    sys.exit(main())
