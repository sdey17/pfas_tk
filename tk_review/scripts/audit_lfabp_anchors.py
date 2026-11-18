#!/usr/bin/env python3
"""Two published L-FABP claims, checked against this compilation's own rows.

Both are claims made BY papers ABOUT the binding literature, and both can be
tested from the convertible L-FABP rows without any new measurement.

(1) ZHAO 2023's QSAR ANCHOR IS 16 COMPOUNDS, NOT 1599.
    Zhao et al., Ecotoxicol Environ Saf 2023;262:115310
    (doi 10.1016/j.ecoenv.2023.115310) docks a 1599-PFAS library against PDB
    2LKK, trains five regression algorithms on the docking scores, and
    predicts L-FABP affinity across that library. Its only tie to measurement
    is a Pearson r = 0.79 against the Kd values of Yang et al. 2020.

    Yang 2020 is now in this database: 25 L-FABP rows, 19 with a numeric Kd.
    But three of those 19 are salt forms of compounds already in the list, and
    for two of them the Kd is IDENTICAL to the free acid (PFOA and PFOA
    ammonium salt both 18.74 uM; PFOS and PFOS potassium salt both 8.1 uM) --
    the same measurement entered twice under two names. So the anchor is 16
    distinct compounds. r = 0.79 is r^2 = 0.62, i.e. the docking leaves 38% of
    the variance in those 16 unexplained, and that is what carries predictions
    for 1599 library compounds in Zhao and, in later work, for tens of
    thousands.

(2) MODARESI 2025's DOSE DESIGN RESTS ON THE MOST EXTREME Kd IN THE
    LITERATURE.
    Modaresi et al., Chem Res Toxicol 2025 (doi 10.1021/acs.chemrestox.5c00199)
    reports that L-FABP and I-FABP knockout does not change PFOS
    toxicokinetics in mice. Its dose selection argues the doses span
    sub-saturating to saturating conditions, using "the reported FABP
    dissociation constant (Kd = 0.18 uM)", with free liver PFOS predicted at
    about Kd at the top dose of 5 mg/kg.

    0.18 uM is Khazaee 2021's PFOS value. It is the lowest PFOS-L-FABP Kd in
    this compilation by a factor of 45 against the next lowest, and 103x below
    the median of the other five sources. The same paper's own introduction
    quotes the literature range as 38 to 879 uM.

    The consequence is quantitative: if free liver PFOS at the top dose is
    0.18 uM, then on a Kd of 18.5 uM (the median of the other five) L-FABP is
    about 1% occupied there, and under 0.5% on the paper's own quoted range.
    The arm of the design intended to saturate L-FABP almost certainly did
    not.

    This does NOT overturn the knockout result. A null at low occupancy is
    still a null, the SPME binding assay independently found 15 of 19 PFAS
    within 20% between genotypes, and no compensatory FABP upregulation was
    found. What it removes is the claim to have tested saturating conditions.

Read-only.
"""
import math
import csv, os, re, statistics, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')

def num(s):
    m = re.match(r'\s*([0-9.]+(?:[eE][-+]?\d+)?)', str(s or ''))
    return float(m.group(1)) if m else None

def lfabp_kd_uM():
    out = []
    for r in csv.DictReader(open(DB, encoding='utf-8')):
        p = (r.get('protein') or '').lower()
        if 'liver fatty acid' not in p or 'albumin' in p:
            continue
        src = f"{(r.get('study') or '').replace('et al.', '').strip().rstrip(',')} {r.get('year') or ''}"
        if src.strip() == 'Han 2021':
            continue                      # quoted values, not measured
        kd, ku = r.get('kd'), (r.get('kd_units') or '').strip()
        if not kd or ',' in str(kd) or num(kd) is None:
            continue
        f = {'uM': 1.0, 'nM': 1e-3, 'mM': 1e3, 'M': 1e6}.get(ku.split()[0] if ku else '')
        if f is None:
            continue
        out.append((src, r.get('protein_species'), r['chemical'], num(kd) * f))
    return out

def main():
    rows = lfabp_kd_uM()
    print("=" * 74)
    print("(1) Zhao 2023's experimental anchor")
    y = [r for r in rows if r[0].startswith('Yang 2020')]
    base = {}
    for _, _, c, k in y:
        b = re.sub(r',?\s*(potassium|ammonium|sodium)\s*salt', '', c, flags=re.I).strip()
        base.setdefault(b, []).append((c, k))
    dups = {b: v for b, v in base.items() if len(v) > 1}
    print(f"  Yang 2020 rows with a numeric Kd : {len(y)}")
    print(f"  distinct compounds behind them   : {len(base)}")
    for b, v in sorted(dups.items()):
        same = "IDENTICAL Kd" if len({k for _, k in v}) == 1 else "different Kd"
        print(f"    {b:10s} appears {len(v)}x, {same}: {[k for _, k in v]}")
    print(f"  => r = 0.79 rests on n = {len(base)}; r^2 = {0.79**2:.2f}")

    print("\n" + "=" * 74)
    print("(2) Modaresi 2025's 0.18 uM against the compilation")
    pf = sorted((k, s) for s, sp, c, k in rows
                if c.strip().upper().startswith('PFOS') and 'SULFONAMIDE' not in c.upper()
                and 'FOSA' not in c.upper())
    for k, s in pf:
        print(f"    {s:16s} Kd {k:8.2f} uM")
    if len(pf) > 1:
        others = [k for k, _ in pf[1:]]
        med = statistics.median(others)
        print(f"  lowest value {pf[0][0]} uM ({pf[0][1]}) is "
              f"{med / pf[0][0]:.0f}x below the median of the rest ({med:.1f} uM)")
        free_top = 0.18          # the paper's own predicted free liver conc at 5 mg/kg
        for lbl, kd in (("paper's Kd", 0.18), ("median of the rest", med),
                        ("paper's quoted low end", 38.0),
                        ("paper's quoted high end", 879.0)):
            print(f"    occupancy at the top dose if Kd = {kd:7.1f} uM "
                  f"({lbl:24s}): {100 * free_top / (free_top + kd):5.1f}%")
    allk = [k for _, _, _, k in rows]
    print(f"\n  whole compilation: {len(allk)} Kd values, median {statistics.median(allk):.1f} uM, "
          f"{sum(1 for k in allk if k <= 0.18)} at or below 0.18 uM")

    print("\n" + "=" * 74)
    print("(3) WHY Khazaee 2021 is the outlier: its capacity, from SI Figure S5")
    print("    The ordinate is specific binding in umol PFAS/umol protein, so each")
    print("    curve carries Bmax as well as Kd. The paper fits a one-site model in")
    print("    GraphPad Prism and reports Kd only. Two plotted points determine both;")
    print("    the recovered Kd checks against the published value.\n")
    B = math.log10(1e9 / 14208.0)
    for nm, (l1, b1, l2, b2), pub in (
            ("PFOA  (S5A)", (0.33, 0.027, 3.75, 0.034), 0.099),
            ("PFHxS (S5D)", (0.25, 0.12, 3.1, 0.65), 1.7)):
        r = b1 / b2
        kd = (l1 * l2 * (1 - r)) / (r * l2 - l1)
        bmax = b1 * (kd + l1) / l1
        print(f"    {nm}  Kd recovered {kd:6.3f} vs published {pub:6.3f} uM"
              f"   ->  Bmax {bmax:6.3f} mol/mol")
        print(f"                log K with n = 1 is {B - math.log10(pub):.2f}; "
              f"with n = Bmax it is {B + math.log10(bmax) - math.log10(pub):.2f}")
    print("\n    PFHxS saturates at one site, as it should. PFOA saturates at 0.035,")
    print("    about one PFOA per thirty L-FABP molecules, so its 0.099 uM describes")
    print("    a small competent sub-population rather than the bulk protein. That is")
    print("    the whole of its 97x gap from every other source. NOT applied to the")
    print("    sheet: n stays 1, because the albumin rule is to use a reported n and")
    print("    Khazaee reports no Bmax. Recorded so the outlier is explained.")
    return 0

if __name__ == '__main__':
    sys.exit(main())
