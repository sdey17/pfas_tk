#!/usr/bin/env python3
"""Correct the Li 2021 conditions from the primary, and mark its n = 1 as the
paper's own.

Li W, Hu Y, Bischel HN. "In-Vitro and In-Silico Assessment of Per- and
Polyfluoroalkyl Substances (PFAS) in Aqueous Film-Forming Foam (AFFF) Binding
to Human Serum Albumin." Toxics 2021;9(3):63. doi:10.3390/toxics9030063

Three things were wrong or missing on the 26 rows:

  1. protein_conc was annotated "stated for the AFFF bound-fraction
     experiment", which reads as though the 600 uM does not cover the 26
     standards. Section 2.2 says otherwise: "One side of each dialysis cell
     was dosed with 7.97 mg HSA ... to a final concentration of 600 uM", and
     "The other side of each cell was then dosed with 200 uL of an AFFF
     dilution (4000 to 16,000-fold in PBS) OR an in-house mixture of 26
     PFAS". One protein side serves both arms, so 600 uM is the protein
     concentration for the standards too.

  2. ligand_conc_range said "not stated for the 26 standards". The paper does
     not print a concentration directly, but it states both numbers that give
     it: 600 uM HSA (Section 2.2) and "we tested PFAS-HSA association
     constants from 0.001 to 0.1 PFAS:HSA" (Discussion). 0.001-0.1 x 600 uM
     is 0.6-60 uM. That is a derivation from two stated quantities, and it is
     recorded as derived, not as reported -- the distinction that the Starnes
     2024 1000 uM error turned on.

  3. The route. The sheet printed "n = 1 assumed by the reader, the paper
     reports none", which is false here. Equations (1) and (2) ARE a
     single-site model: "Noncovalent binding affinities, measured as
     association constants (KA), were calculated assuming a one site specific
     binding". The n = 1 is the paper's, so the row belongs in the
     PAPER_OK/PAPER_EX class, not the reader class. Handled by
     PAPER_SINGLE_SITE in albumin_shareable.py; this script records why.

Also sets temperature (37 C, Section 2.2: incubated at 37 C, 30 RPM, 108 h)
and buffer (PBS pH 7.4), both stated and both previously blank.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.3390/toxics9030063'

WANT = {
    'protein_conc': '600',
    'protein_conc_units': 'uM',
    'ligand_conc_range': ('0.6-60 (DERIVED, not printed as a concentration: '
                          '0.001-0.1 PFAS:HSA, stated in the Discussion, times '
                          'the 600 uM HSA of Section 2.2)'),
    'ligand_conc_units': 'uM',
    'ligand_over_protein_max': '0.1',
    'temperature_c': '37',
    'buffer': ('PBS pH 7.4, 10 kDa regenerated cellulose membrane, 96-well '
               'Harvard Apparatus, 37 C, 30 RPM, 108 h; six dilution levels '
               'x 4 replicate trials'),
}

def main():
    with open(DB, newline='', encoding='utf-8') as fh:
        rdr = csv.DictReader(fh); fields, rows = rdr.fieldnames, list(rdr)
    for k in WANT:
        if k not in fields:
            print(f"column {k!r} absent from the database", file=sys.stderr)
            return 1
    changed = 0
    for r in rows:
        if (r.get('doi') or '') != DOI:
            continue
        for k, v in WANT.items():
            if (r.get(k) or '') != v:
                r[k] = v; changed += 1
    if not changed:
        print("no change needed (Li 2021 conditions already correct)")
        return 0
    shutil.copy2(DB, DB + '.bak')
    with open(DB, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"updated {changed} cell(s) across the Li 2021 rows")
    for k, v in WANT.items():
        print(f"  {k:24s} = {v[:68]}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
