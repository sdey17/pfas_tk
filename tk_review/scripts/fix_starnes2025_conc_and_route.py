#!/usr/bin/env python3
"""Remove a dosing concentration Starnes 2025 never states, and record that its
n = 1 is the paper's own.

Starnes HM, Green AJ, Reif DM, Belcher SM. "An in vitro and machine learning
framework for quantifying serum albumin binding of per- and polyfluoroalkyl
substances." Toxicol Sci 2025;203(1):67-78. doi:10.1093/toxsci/kfae124

PROBLEM 1 -- a concentration the paper does not give.
ligand_conc_range held prose that ENDED in parseable numbers:

  "not reported by the paper; stocks serially diluted to unstated working
   concentrations, and Table S3 gives the stock concentration and solvent
   only; positive controls octanoic acid 1-3 mM and PFHxA 0.3-3 mM"

The note said "not reported" in words, but the trailing "0.3-3 mM" parsed,
so every one of the 49 rows reached the sheet with PFAS conc 3000 uM and a
loading of 23.1 over the 130 uM HSA. Those are the POSITIVE CONTROLS, run on
every plate alongside the test compound (Methods, Differential scanning
fluorimetry): "Controls on each plate included a no protein control,
solvent-matched vehicle control (no ligand), and 3 concentrations of octanoic
acid (1, 2, 3 mM) and/or PFHxA (0.3, 1, 3 mM) as positive controls."

What the paper says about the test compounds is only: "PFAS and fatty acid
stock solutions were serially diluted to achieve necessary working stock
concentrations for DSF concentration-response analysis" -- stocks of 5, 10 or
20 mM, and no working series anywhere. So the cell must be blank.

This is the same failure already found in Starnes 2024, where 1000 uM came
from one figure caption for PFOS on BSA alone. A loading of 23 would also have
been the highest in the sheet by a wide margin, on rows whose dosing is
unknown. The fix puts the control concentrations in assay_volume_note, where
nothing parses them as a dosing range.

PROBLEM 2 -- the route misattributed the stoichiometry.
The sheet printed "n = 1 assumed by the reader, the paper reports none". The
paper reports no n column, but it does fix n by its own model (Methods, DSF
data analysis): "a single-site ligand binding model was used to calculate
dissociation constants". That is the same basis as Li 2021, so the rows belong
in the PAPER_OK/PAPER_EX class. Handled by PAPER_SINGLE_SITE in
albumin_shareable.py.

The Ka values, the 130 uM HSA (stated: "Each well contained a final
concentration of 0.13 mM HSA") and everything else are untouched.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.1093/toxsci/kfae124'

NOTE = ('20 uL reaction in a sealed optical 96-well plate; >=2 independent '
        'plates per protein-ligand pair, 3 replicates per concentration per '
        'plate; 0.13 mM HSA per well, the concentration giving optimal '
        'signal-to-noise. Plate controls (NOT a dosing series for the test '
        'compounds): no-protein control, solvent-matched vehicle control, '
        'octanoic acid at 1, 2, 3 mM and/or PFHxA at 0.3, 1, 3 mM as positive '
        'controls. Test-compound working concentrations are not reported: '
        'stocks of 5, 10 or 20 mM (Table S3) were "serially diluted to achieve '
        'necessary working stock concentrations".')
NEW_RANGE = ''

def main():
    with open(DB, newline='', encoding='utf-8') as fh:
        rdr = csv.DictReader(fh); fields, rows = rdr.fieldnames, list(rdr)
    changed = 0
    for r in rows:
        if (r.get('doi') or '') != DOI:
            continue
        for k, v in (('ligand_conc_range', NEW_RANGE),
                     ('ligand_over_protein_max', ''),
                     ('assay_volume_note', NOTE)):
            if (r.get(k) or '') != v:
                r[k] = v; changed += 1
    if not changed:
        print("no change needed (Starnes 2025 already corrected)")
        return 0
    shutil.copy2(DB, DB + '.bak')
    with open(DB, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"updated {changed} cell(s) across the Starnes 2025 rows")
    print("  cleared ligand_conc_range (was prose ending in '0.3-3 mM', which")
    print("          parsed to 3000 uM on all 49 rows -- the positive control)")
    print("  cleared ligand_over_protein_max (was 23.08)")
    print("  moved the control concentrations into assay_volume_note")
    return 0

if __name__ == '__main__':
    sys.exit(main())
