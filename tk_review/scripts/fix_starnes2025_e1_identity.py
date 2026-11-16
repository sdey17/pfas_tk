#!/usr/bin/env python3
"""Record what Starnes 2025 calls "E1", and correct a co-solvent coverage claim.

E1 is the paper's own abbreviation (Table S4, class PFAE, CAS 3330-15-2) for
1,1,1,2,2,3,3-heptafluoro-3-(1,2,2,2-tetrafluoroethoxy)propane, known as
FLUOROETHER E-1 or heptafluoropropyl 1,2,2,2-tetrafluoroethyl ether.
C5HF11O, MW 286.04, PubChem CID 94258, DTXSID8052017.

It matters for three reasons, none of which "E1" conveys:

  1. It is NEUTRAL. No carboxylate or sulfonate -- an ether with a hydrogen
     endcap. Almost everything else in the sheet is anionic at pH 7.4, and
     albumin binds PFAS largely through electrostatic contact with the head
     group, so this row is a different binding mode, not a different chain
     length. The database already had head_group "ether (neutral)".

  2. It belongs to GenX chemistry. Fluoroether E-1 is used in the GenX
     process and can arise from degradation of GenX chemicals including
     FRD-903. Starnes 2025 measures HFPO-DA (GenX itself) in the same run, and
     the two land almost on top of each other: log K 0.963 for E1 against
     0.977 for HFPO-DA.

  3. It boils at 40-42 C and is practically insoluble in water. DSF reads a
     shift in albumin's melting temperature, and albumin melts far above 42 C
     (this study reports dTm up to 20.6 C). So the compound's boiling point
     sits below the temperature at which its effect is measured. The plate is
     sealed, which contains it, but the solution-phase concentration during
     the transition is not the dosed concentration. This is recorded as a
     caveat on the row, not acted on: the Ka and its R2 of 0.93 are what the
     paper reports, and the sheet has no other volatile-ligand row to
     calibrate against.

CORRECTION TO AN EARLIER CLAIM. It was previously stated that "3 of the
admitted Starnes 2025 rows carry a co-solvent". That overstated what is known.
The extracted Table S3 holds 41 CAS numbers while the paper tested 68
chemicals, and it covers only 14 of the 36 admitted sheet rows. So the right
statement is: 3 of the 14 rows whose solvent is recorded were prepared in 50%
DMSO, and the solvent for the other 22 -- PFOA, PFOS, PFBA and the other
legacy compounds among them, and E1 -- is not known from what has been
extracted. Those are mostly water-soluble anionic surfactants that would not
need a co-solvent, but that is an inference, not a reading. Whether the gap is
Table S3's own scope or the extraction's is unresolved: the 135 MB
supplementary workbook is not kept in the repository.

Sources for the identity: PubChem CID 94258; Wikipedia "Fluoroether E-1",
citing Zhang et al., Environ Sci Technol 2022;56(10):6103-6112
(doi 10.1021/acs.est.1c03979) and Wickersham, J Air Waste Manag Assoc
2023;73(7):533-552 (doi 10.1080/10962247.2023.2192009).

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.1093/toxsci/kfae124'

NOTE = ('Fluoroether E-1 = heptafluoropropyl 1,2,2,2-tetrafluoroethyl ether, '
        'CAS 3330-15-2, C5HF11O, MW 286.04, PubChem CID 94258, '
        'DTXSID8052017. "E1" is the paper\'s abbreviation (Table S4, class '
        'PFAE). NEUTRAL ether with a hydrogen endcap, not an acid, so not '
        'anionic at pH 7.4 like most of this sheet. Used in the GenX process '
        'and can arise from degradation of GenX chemicals including FRD-903. '
        'CAVEAT: boiling point 40-42 C and practically insoluble in water, '
        'while DSF reads a shift in albumin\'s melting temperature, which is '
        'far above 42 C -- so the solution-phase concentration during the '
        'measured transition is not the dosed concentration. Sealed plate. '
        'Not in the extracted Table S3, so its solvent is unrecorded.')

def main():
    with open(DB, newline='', encoding='utf-8') as fh:
        rdr = csv.DictReader(fh); fields, rows = rdr.fieldnames, list(rdr)
    col = 'binding_site_identified'
    if col not in fields:
        print(f"no {col!r} column", file=sys.stderr); return 1
    hit = [r for r in rows
           if (r.get('doi') or '') == DOI and r.get('chemical', '').strip() == 'E1']
    if not hit:
        print("no E1 row found", file=sys.stderr); return 1
    if all((r.get(col) or '') == NOTE for r in hit):
        print("no change needed (E1 identity already recorded)"); return 0
    for r in hit:
        r[col] = NOTE
    shutil.copy2(DB, DB + '.bak')
    with open(DB, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"recorded the E1 identity on {len(hit)} row(s)")
    print("  Fluoroether E-1, CAS 3330-15-2, neutral ether, GenX-process related")
    print("  caveat recorded: bp 40-42 C, below the melting transition measured")
    return 0

if __name__ == '__main__':
    sys.exit(main())
