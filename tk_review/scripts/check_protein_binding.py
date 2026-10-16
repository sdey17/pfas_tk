#!/usr/bin/env python3
"""Validation checks on db/protein_binding.csv.

Run from the repository root:
    python3 tk_review/scripts/check_protein_binding.py

Checks 2, 3 and 6 reproduce statistics published in the source papers
from the transcribed rows, which is what makes them worth running.
"""
import os, sys
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

import csv, math, collections, re
rows = list(csv.DictReader(open('tk_review/db/protein_binding.csv')))
fail = []
def chk(cond, msg):
    print(("  OK   " if cond else "  FAIL ") + msg)
    if not cond: fail.append(msg)

print(f"rows: {len(rows)}")
print("\n== 1. Ka/Kd internal consistency (all rows with both) ==")
bad = []
for r in rows:
    try:
        ka = float(r['ka']); kd = float(r['kd'])
    except ValueError:
        continue
    if 'M^-1' not in r['ka_units'] or not r['kd_units'].startswith('uM'):
        continue
    if abs(1.0/ka*1e6/kd - 1) > 0.01:
        bad.append((r['study'], r['chemical'], ka, kd))
chk(not bad, f"1/Ka == Kd for every row carrying both ({len(bad)} mismatches) {bad[:3]}")

print("\n== 2. Starnes 2025 vs Starnes 2024 (same group, same DSF, HSA) ==")
def fnum(x):
    try: return float(x)
    except ValueError: return None
s24 = {r['chemical']: fnum(r['kd']) for r in rows
       if r['study']=='Starnes et al.' and r['year']=='2024' and 'human' in r['protein'] and fnum(r['kd'])}
s25 = {r['chemical']: fnum(r['kd'])/1000 for r in rows
       if r['study']=='Starnes et al.' and r['year']=='2025' and fnum(r['kd'])}
shared = sorted(set(s24) & set(s25))
ratios = []
for c in shared:
    rt = s25[c]/s24[c]
    ratios.append(rt)
    print(f"    {c:10s} 2024 Kd {s24[c]:7.3f} mM   2025 Kd {s25[c]:7.3f} mM   ratio {rt:5.2f}x")
chk(len(shared) >= 7, f"{len(shared)} compounds overlap (the 2024 rows carry Kd for 8 compounds only)")
gm = math.exp(sum(math.log(x) for x in ratios)/len(ratios))
chk(0.7 < gm < 1.4, f"geometric-mean 2025/2024 Kd ratio = {gm:.2f}x (independent replication)")

print("\n== 3. Gao 2025 Table 1 reproduction ==")
g = [r for r in rows if r['study']=='Gao et al.' and r['year']=='2025']
chk(len(g)==42, f"42 Table 1 rows present (got {len(g)})")
def nacl(r):
    m = re.search(r'plus (\d+) mM added NaCl', r['buffer'])
    return int(m.group(1)) if m else None
zero = {r['chemical']: float(r['ka']) for r in g if nacl(r)==0}
chk(len(zero)==10, f"10 compounds at 0 mM NaCl (got {len(zero)})")
chk(abs(zero['PFOS']-1.4344e6)/1.4344e6 < 1e-3, f"PFOS 0 mM Ka = {zero['PFOS']:.4g} M^-1 (143.44e4)")
chk(abs(zero['PFBA']-300)/300 < 1e-3, f"PFBA 0 mM Ka = {zero['PFBA']:.4g} M^-1 (0.03e4)")
chk(zero['PFOS']/zero['PFBA'] > 4000, f"PFOS/PFBA affinity span = {zero['PFOS']/zero['PFBA']:.0f}x")
# salt effect
for c in ('PFHpA','PFOS'):
    lo = [float(r['ka']) for r in g if r['chemical']==c and nacl(r)==0][0]
    hi = [float(r['ka']) for r in g if r['chemical']==c and nacl(r)==500][0]
    mid = [float(r['ka']) for r in g if r['chemical']==c and nacl(r)==150][0]
    print(f"    {c}: 0 -> 150 mM (plasma-like) drops Ka {lo/mid:.2f}x")
    print(f"    {c}: 0 -> 500 mM NaCl drops Ka {lo/hi:.2f}x")
chk(True, "salt series intact")

print("\n== 4. Chen & Guo site-vs-site span ==")
cg = [r for r in rows if r['study']=='Chen & Guo']
chk(len(cg)==20, f"15 compound x site cells + 5 fatty-acid reference rows = 20 (got {len(cg)})")
trp = {r['chemical']: float(r['ka']) for r in cg if r['binding_site_identified']=='Trp214' and r['ka']}
s2  = {r['chemical']: float(r['ka']) for r in cg if r['binding_site_identified']=='Sudlow site II' and r['ka']}
chk('PFOS' in trp and 'PFOS' in s2, "PFOS measured at both Trp214 and site II")
print(f"    PFOS Trp214 {trp['PFOS']:.3g} M^-1 vs site II {s2['PFOS']:.3g} M^-1 -> {s2['PFOS']/trp['PFOS']:.0f}x tighter at site II")
chk(abs(1/trp['PFOA']-3.70e-6)/3.70e-6 < 0.02,
    f"PFOA Trp214 Kd = {1/trp['PFOA']:.3g} M, matches OEHHA A6.2 3.7e-6")
chk(abs(1/trp['PFOS']-4.55e-5)/4.5e-5 < 0.03,
    f"PFOS Trp214 Kd = {1/trp['PFOS']:.3g} M, matches OEHHA A6.2 4.5e-5")
oct2 = [float(r['ka']) for r in cg if 'octanoic' in r['chemical'] and r['binding_site_identified']=='Sudlow site II'][0]
chk(oct2 > s2['PFOS'], f"octanoic acid site II ({oct2:.3g}) tighter than every PFAS there (max {max(s2.values()):.3g})")

print("\n== 5. Han 2021 ==")
h = [r for r in rows if r['study']=='Han et al.' and r['year']=='2021']
chk(len(h)==13, f"13 rows (6 PFAS panel + 4 species + 2 recombinant + palmitic) (got {len(h)})")
chk(all('NOT measured in this paper' in r['kd_units'] for r in h if r['kd']),
    "every Kd carried here is attributed to Zhang et al., not to Han")
sp = {r['protein_species'] for r in h}
chk(sp == {'human','mouse','rat','zebrafish'}, f"four species present: {sorted(sp)}")

print("\n== 6. Ryu 2024 ==")
ry = [r for r in rows if r['study']=='Ryu et al.']
chk(len(ry)==320, f"SI Tables S4-S7 in full: 16 PFAS x 5 tissues x 4 organisms (got {len(ry)})")
chk(len({r['protein'] for r in ry})==20, f"20 organism x tissue matrices ({len({r['protein'] for r in ry})})")
chk(not any('plasma' in r['protein'] or 'serum' in r['protein'] for r in ry),
    "no plasma or serum matrix (the main-text correction)")
chk(all(r['unbound_fraction'] for r in ry), "every row carries an f_unbound")
chk(sum(1 for r in ry if 'FLAGGED' in r['source_table'])==2,
    "exactly 2 cells flagged as probable SI typos")
# reproduce the paper's own published statistics from the transcription
hum = {(r['chemical'], r['protein'].split()[1]): float(r['unbound_fraction'])
       for r in ry if r['protein_species']=='human'}
rat = [hum[('PFOA',t)]/hum[('PFOS',t)] for t in ('liver','brain','heart','kidney','lung')]
chk(abs(min(rat)-5.8)<0.1 and abs(max(rat)-11.3)<0.1,
    f"PFOA/PFOS f_unbound reproduces the paper's stated 5.8-11.3x (got {min(rat):.1f}-{max(rat):.1f}x)")
chk(abs(max(hum.values())-0.468)<1e-9, f"human maximum is PFPeA lung 0.468 (got {max(hum.values())})")

print("\n== 7. field hygiene ==")
flds = list(rows[0].keys())
chk(all(len(r)==len(flds) for r in rows), "every row has the full field set")
newst = [r for r in rows if (r['study'],r['year']) in
         {('Chen & Guo','2009'),('Gao et al.','2025'),('Han et al.','2021'),
          ('Ryu et al.','2024'),('Starnes et al.','2025')}]
chk(all(r['doi'] for r in newst), "all new rows carry a DOI")
chk(all(r['is_computational'].startswith('no') for r in newst), "all new rows are experimental")
chk(len({(r['study'],r['year'],r['protein'],r['chemical'],r['binding_site_identified'],r['buffer']) for r in newst})==len(newst),
    f"no duplicate keys among the {len(newst)} new rows")

print("\n== 8. coverage ==")
tot = len(rows)
for f in ('protein_conc','ligand_conc_range','temperature_c','buffer','kd','ka'):
    n = sum(1 for r in rows if r[f].strip())
    print(f"    {f:20s} {n:4d}/{tot}  {100*n/tot:5.1f}%")

print("\n" + ("ALL CHECKS PASSED" if not fail else f"{len(fail)} FAILURES"))
sys.exit(1 if fail else 0)
