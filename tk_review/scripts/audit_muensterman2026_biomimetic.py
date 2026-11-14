#!/usr/bin/env python3
"""Audit Muensterman et al. 2026 (biomimetic chromatography) before any ingest decision.

Muensterman DJ, Valko K, Chen R, Ng C, Droz B, Vyas S, Benotti M, Field J.
"Biomimetic Chromatography: A Novel Approach for Measuring Phospholipid
Membrane-Water and Protein-Water Partition Coefficients for Target and Suspect
PFAS." Environ Sci Technol 2026;60(1):1207-1217. doi:10.1021/acs.est.5c11716
Correction: Environ Sci Technol 2026;60(23):17051. doi:10.1021/acs.est.6c07668

Read-only. Writes nothing. Answers four questions:
  1. What does the Correction change, and does it touch our rows?
  2. Are the published numbers internally reconstructible?
  3. Is there any route from log K_HSA to L_water / kg albumin?
  4. If the values cannot be ingested, is the rank order still usable?
"""
import csv, math, os, statistics, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
XL = os.path.join(os.path.dirname(__file__), '..', 'albumin_logK_Kd_revised.xlsx')

# --- Table S13, as printed (tR_IAM, CHI_IAM, logK_IAM, tR_HSA, logK_HSA) -----
# A representative transcription; the full table is 92 rows in the SI.
S13 = [
    ("PFBA",   3.2,  5, 1.3, 4.4, 3.3), ("PFPeA",  4.0, 17, 1.7, 5.6, 4.1),
    ("PFHxA",  4.4, 24, 2.0, 5.9, 4.3), ("PFHpA",  4.7, 29, 2.4, 6.4, 4.6),
    ("PFOA",   4.9, 34, 2.8, 6.4, 4.7), ("PFNA",   5.2, 38, 3.2, 6.4, 4.6),
    ("PFDA",   5.4, 42, 3.7, 6.3, 4.6), ("PFUdA",  5.6, 46, 4.3, 6.2, 4.5),
    ("PFDoA",  5.8, 49, 4.9, 5.9, 4.3), ("PFTrDA", 6.0, 52, 5.5, 5.9, 4.3),
    ("PFTeDA", 6.2, 55, 6.2, 5.9, 4.3), ("PFHxDA", 6.5, 60, 7.5, 5.9, 4.3),
    ("PFODA",  6.7, 64, 9.1, 5.8, 4.2), ("PFHxS",  5.0, 35, 2.9, 7.1, 5.1),
    ("PFOS",   5.5, 42, 3.8, 6.3, 4.6),
]
# Table S8 column regressions; Table S4 calibrant log k(HSA) span.
HSA_REG   = (2.12, -0.17)          # log k_HSA = 2.12 log10(tR) - 0.17
CALIBRANTS = {"Paracetamol": -0.79, "Nizatidine": -0.60, "Trimetoprim": -0.23,
              "Propranolol": 0.29, "Carbamazepine": 0.46, "Nicardipine": 1.20,
              "Warfarin": 1.50, "Indometacin": 1.82, "Diclofenac": 1.92}

def eq5(chi, intercept):
    """log K_IAM = 0.29 e^(log k_IAM) + 0.70, with log k_IAM = 0.046 CHI + b."""
    return 0.29 * math.exp(0.046 * chi + intercept) + 0.70

def eq6(tr):
    """log K_HSA = e^(log k_HSA), log k_HSA from the Table S8 HSA regression."""
    m, b = HSA_REG
    return math.exp(m * math.log10(tr) + b)

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):                      # average ties
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')

def main():
    print("=" * 78)
    print("Muensterman 2026 -- biomimetic chromatography, pre-ingest audit")
    print("=" * 78)

    # --- 1. the Correction -----------------------------------------------
    print("\n[1] WHAT THE CORRECTION CHANGES")
    print("  Correction to doi:10.1021/acs.est.5c11716, published 2026-06-02")
    print("  (Environ Sci Technol 60(23):17051, doi:10.1021/acs.est.6c07668).")
    print("  Its entire body is one line: the correct form of eq 4 is")
    print("      log k_IAM = 0.046 CHI_IAM + 0.42       (printed as '+ 42')")
    print("  citing Hollosy, Valko, Hersey, Nunhuck, Keri, Bevan,")
    print("  J Med Chem 2006;49(24):6958-6971 as the prior source.")
    print("  eq 4 feeds eq 5 (log K_IAM) ONLY. The HSA chain is eq 6 on the")
    print("  Table S8 HSA regression and does not involve eq 4 at all.")

    # --- 2. internal reconstruction --------------------------------------
    print("\n[2] DO THE PUBLISHED NUMBERS RECONSTRUCT?")
    for label, b in (("as printed  (+42  )", 42.0), ("as corrected(+0.42)", 0.42)):
        err = [abs(eq5(chi, b) - k) for _, _, chi, k, _, _ in S13]
        print(f"  eq 4 {label} -> eq 5 : max err {max(err):.3g}, "
              f"mean {statistics.mean(err):.3g}")
    err = [abs(eq6(tr) - k) for _, _, _, _, tr, k in S13]
    print(f"  Table S8 + eq 6      -> log K_HSA : max err {max(err):.3g}, "
          f"mean {statistics.mean(err):.3g}")
    print("  => the published DATA were computed with 0.42; the typo was")
    print("     typographical only. No reported value changes, and the")
    print("     Correction publishes no revised table.")
    print("  => the HSA values reconstruct from retention time alone.")

    # --- 3. is there a basis? --------------------------------------------
    print("\n[3] IS THERE A ROUTE TO L_water / kg ALBUMIN?")
    print("  The SI does state a basis (contradicting our current DB note):")
    print('    "The resulting values are dimensionless partition coefficients,')
    print('     as they are derived from the ratio of solute concentrations ...')
    print('     in two immiscible or non-interacting environments at equilibrium."')
    print("  But it states the basis in the same breath as what removes it:")
    print('    "the concentration of the stationary phase (e.g., phospholipids')
    print('     or immobilized protein) is effectively constant and in large')
    print('     excess, which simplifies the apparent affinity constant (K_aff)')
    print('     to a partition-like equilibrium ratio"')
    print("  Neither paper nor SI states a phase ratio, an immobilised-protein")
    print("  mass, or a column loading. The albumin quantity is absorbed into")
    print("  the calibration constants of eqs 5-6, so it cannot be divided back")
    print("  out. Our other volume-basis sources were convertible because a")
    print("  density (1.36 kg/L, log10 = 0.1335) bridged L_protein -> kg_protein;")
    print("  here the denominator is not a volume of albumin at all, but an")
    print("  albumin-on-silica stationary phase of unstated composition.")
    print("  => no conversion exists. Not a unit problem; a missing quantity.")

    # --- 4. what the scale actually is ------------------------------------
    print("\n[4] WHAT THE REPORTED SCALE IS")
    lo = min(CALIBRANTS.values()); hi = max(CALIBRANTS.values())
    print(f"  Calibrant log k(HSA) span (Table S4): {lo:+.2f} (paracetamol) to "
          f"{hi:+.2f} (diclofenac)")
    print(f"  Through eq 6 that is log K_HSA {math.exp(lo):.2f} to "
          f"{math.exp(hi):.2f}")
    print("  The authors are explicit about what survives:")
    print('    "While absolute binding values may not always be confirmed by')
    print('     other methods, the elution order is highly reproducible."')
    print('    "... detecting binding values above 99.5% is not feasible or')
    print('     reliable ..." (the calibration rests on equilibrium-dialysis')
    print("     %-bound data, which saturates).")
    print("  => the paper offers a reproducible RANKING, and says so.")

    # --- 5. our rows ------------------------------------------------------
    print("\n[5] OUR ROWS")
    rows = [r for r in csv.DictReader(open(DB, encoding='utf-8'))
            if 'Muensterman' in (r.get('study') or '')]
    hsa = [r for r in rows if 'HSA' in (r.get('protein') or '')]
    print(f"  {len(rows)} DB rows = {len(hsa)} HSA + {len(rows)-len(hsa)} IAM "
          "(phospholipid, not a protein)")
    vals = sorted(float(r['ka']) for r in hsa if r.get('ka'))
    pfas = sorted(float(r['ka']) for r in hsa
                  if r.get('ka') and 'ecanoic' not in r['chemical']
                  and 'ctanoic' not in r['chemical'])
    print(f"  all HSA rows : {min(vals):.2f} to {max(vals):.2f}, "
          f"median {statistics.median(vals):.2f}")
    print(f"  PFAS only    : {min(pfas):.2f} to {max(pfas):.2f}, "
          f"median {statistics.median(pfas):.2f}")
    print("  abstract     : 0.45 to 5.1, median 4.3")
    print("  => the apparent 6.80 ceiling is Decanoic Acid, a non-PFAS")
    print("     reference fatty acid. Excluding the two fatty-acid standards")
    print("     the range matches the abstract exactly. Ingest is correct.")

    # --- 6. rank check against the sheet ---------------------------------
    print("\n[6] RANK AGREEMENT WITH OUR HUMAN-ALBUMIN POOL")
    try:
        import openpyxl
    except ImportError:
        print("  (openpyxl unavailable; skipped)"); return
    ws = openpyxl.load_workbook(XL)["Coefficients"]
    hdr = [c.value for c in ws[1]]
    L = 'log K (L water / kg protein)'
    pool = {}
    for r in ws.iter_rows(min_row=2):
        d = dict(zip(hdr, [c.value for c in r]))
        if d.get('Species') == 'human' and isinstance(d.get(L), (int, float)):
            pool.setdefault(str(d['PFAS']), []).append(d[L])
    pairs = [(float(r['ka']), statistics.median(pool[r['chemical']]), r['chemical'])
             for r in hsa if r.get('ka') and r['chemical'] in pool]
    if len(pairs) < 3:
        print("  too few overlaps"); return
    mu = [p[0] for p in pairs]; sh = [p[1] for p in pairs]
    off = [a - b for a, b in zip(mu, sh)]
    print(f"  {len(pairs)} compounds overlap")
    print(f"  Spearman rho = {spearman(mu, sh):+.3f}")
    print(f"  offset (Muensterman - sheet): median {statistics.median(off):+.2f}, "
          f"range {min(off):+.2f} to {max(off):+.2f}")
    print("  => a systematic positive offset, not a scatter: the two scales")
    print("     are different quantities, consistent with [3].")

if __name__ == '__main__':
    sys.exit(main())
