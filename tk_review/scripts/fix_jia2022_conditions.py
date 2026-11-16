#!/usr/bin/env python3
"""Record the Jia 2022 HSA assay conditions, which are not the sheet's norm.

Jia Y, Zhu Y, Xu D, Feng X, Yu X, Shan G, Zhu L. "Insights into the
Competitive Mechanisms of Per- and Polyfluoroalkyl Substances Partition in
Liver and Blood." Environ Sci Technol 2022;56(9):6192-6200.
doi:10.1021/acs.est.1c08493

WHAT IS INGESTABLE. Of 17 database rows, 11 are hL-FABP (a liver protein, out
of scope for an albumin sheet) and 6 are HSA. One of those 6 is a collapsed
stub naming five compounds in one cell with no Kd. So five HSA rows carry a
per-compound constant: HFPO-DA 1156, HFPO-TA 84, HFPO-TeA 179, PFEESA 874 and
8:2 Cl-PFESA 31 uM.

WHY ONLY FIVE, WHEN TABLE 1 PRINTS ELEVEN. The HSA column of Table 1 carries a
superscript 19 on six of its values, and reference 19 is Gao et al., Environ
Sci Technol 2019;53:6529-6538 -- a paper this project already holds as a
primary. Those six are CITED, not measured here: PFHxA 195, PFOA 115, PFNA 73,
PFDoA 103, PFOS 38 and 6:2 Cl-PFESA 67 uM. Verified by rebuilding them through
the bridge against our own Gao 2019 rows: PFHxA 1.887, PFOA 2.117, PFNA 2.314,
PFOS 2.598 and 6:2 Cl-PFESA 2.351 all match to 0.000, and PFDoA 2.165 matches
the Gao 2019 row carried under the name PFDoDA. Ingesting them would have
double-counted Gao 2019 under Jia's name. In a layout where the citation
marker runs straight into the digits ("3819" for 38 with superscript 19) that
is an easy error to make and a hard one to see.

CONDITIONS, which differ from the sheet's norm on three counts (SI Text, MST):

  pH 8.0, not 7.4. The running buffer is 20 mM Tris, 150 mM NaCl, 10 mM
  MgCl2, pH 8.0. Albumin's PFAS affinity is pH-sensitive through
  conformational change -- Bischel 2011 attributes an observed pH dependence
  to exactly that -- so this is not a free parameter.

  0.05% Tween-20 in the buffer. A non-ionic surfactant present throughout an
  assay whose ligands are themselves surfactants. Above its CMC Tween-20 forms
  micelles that can sequester PFAS and lower the free concentration the fit
  assumes.

  The HSA is covalently labeled. RED-NHS is amine-reactive, so it modifies
  lysines, and lysine contacts are part of how albumin binds an anionic PFAS
  head group.

Also: the paper states no binding model for the HSA fits ("the Kd was finally
fitted by MO Affinity Analysis software"), so n = 1 is the READER's assumption
here, not the paper's. That is weaker than Li 2021 or Starnes 2025, each of
which names a single-site model in its own text, and the route label must say
so rather than borrowing their wording.

Temperature: Table 1's title gives 25 C for hL-FABP and does not state one for
the HSA column; the SI says measurements ran "at different temperatures from 21
to 40 oC". 25 C is recorded as the database already had it, with the caveat
that the paper does not pin the HSA column to it.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
DOI = '10.1021/acs.est.1c08493'

BUF = ('20 mM Tris, 150 mM NaCl, 10 mM MgCl2, pH 8.0, 0.05% Tween-20 (MST '
       'running buffer). NOTE pH 8.0, not the 7.4 of most of this sheet; '
       'albumin PFAS affinity is pH-sensitive via conformational change. NOTE '
       '0.05% Tween-20, a non-ionic surfactant present throughout an assay '
       'whose ligands are surfactants. HSA covalently labeled with RED-NHS '
       '(amine-reactive, modifies lysines, which participate in head-group '
       'binding). Final assay protein concentration not stated: 10 uM protein '
       'was used for labeling, then 10 uL labeled protein was mixed 1:1 with '
       '10 uL ligand.')
TEMP = ('25 (as the database had it; Table 1 states 25 C for the hL-FABP '
        'column and gives no temperature for the HSA column, and the SI says '
        'runs were made "at different temperatures from 21 to 40 oC")')

def main():
    with open(DB, newline='', encoding='utf-8') as fh:
        rdr = csv.DictReader(fh); fields, rows = rdr.fieldnames, list(rdr)
    changed = 0
    for r in rows:
        if (r.get('doi') or '') != DOI or 'HSA' not in (r.get('protein') or ''):
            continue
        for k, v in (('buffer', BUF), ('temperature_c', TEMP)):
            if (r.get(k) or '') != v:
                r[k] = v; changed += 1
    if not changed:
        print("no change needed (Jia 2022 HSA conditions already recorded)")
        return 0
    shutil.copy2(DB, DB + '.bak')
    with open(DB, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"updated {changed} cell(s) on the Jia 2022 HSA rows")
    print("  recorded pH 8.0, 0.05% Tween-20, RED-NHS labeling, and the")
    print("  temperature caveat for the HSA column")
    return 0

if __name__ == '__main__':
    sys.exit(main())
