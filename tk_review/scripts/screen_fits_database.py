#!/usr/bin/env python3
"""Screen a source for whether it can EVER reach log K [L water / kg protein].

Run this BEFORE reading a paper. The working sheet is built on one identity:

    log K = 4.1774 + log10(n) - log10(Kd / uM)

So a source earns a sheet row only if it reports, for albumin, in a named
species, a genuine dissociation constant and a site count (or a site count the
reader is willing to set to 1). Nothing else about the paper matters -- not
how careful the method is, not whether its arithmetic checks out, not whether
it has a published correction. If the two quantities are absent, there is no
row, and reading further is wasted effort.

This exists because Muensterman 2026 was read end to end before that test was
applied. Its rows carry a populated `ka` column, which is why a naive "does it
have a number?" check passed it: the number is a chromatographic retention
index (log K_HSA = e^(log k_HSA)) with no phase ratio or albumin mass behind
it. Hence NOT_A_CONSTANT below -- a populated cell is not a binding constant.

Usage:
    python3 -I scripts/screen_fits_database.py              # all albumin sources
    python3 -I scripts/screen_fits_database.py Jia Li       # just these
"""
import collections, csv, os, sys

DB = os.path.join(os.path.dirname(__file__), '..', 'db', 'protein_binding.csv')
XL = os.path.join(os.path.dirname(__file__), '..', 'albumin_logK_Kd_revised.xlsx')
SHEET_SPECIES = {'human', 'bovine'}

# Sources the reviewer has already ruled on. They pass or fail the screen on
# their own merits, but they are not worklist: the decision is made. Keeping
# them here stops the screen from re-proposing settled papers every run.
DECIDED = {
    'Wu 2024':        'named holdout (reviewer removed; barred from figures)',
    'Beesoon 2015':   'scope drop (branched/linear isomer mixture)',
    'Wu 2009':        'scope drop (Temkin dropped, ITC merged)',
    'Klevens 1954':   'scope drop (reviewer: not worth it)',
    'Peng 2025':      'scope drop (secondary stub: 7 compounds x 2 species collapsed)',
    'Fedorenko 2021': 'reviewer skipped (19F NMR, 6 values stated in words)',
    'Fan 2025':       'reviewer dropped (rat only)',
    'Han 2003':       'reviewer dropped (rat only)',
    'Muensterman 2026': 'cannot fit (biomimetic retention index, no phase ratio)',
}

# A populated ka/kd cell that is NOT a dissociation or association constant.
NOT_A_CONSTANT = (
    'biomimetic', 'retention', 'not a dissociation constant', 'binding energy',
    'kcal/mol', 'log k_hsa', 'log k_iam', 'log k_agp', 'docking', 'temkin',
    'gas phase', 'predicted', 'dimensionless partition',
)
# Methods that cannot yield a solution-phase Kd at all.
NO_SOLUTION_KD = (
    'docking', 'crystallograph', 'qsar', 'pbpk', 'bioconcentration model',
    'literature review', 'systematic literature review', 'read-across',
    'mechanistic protein-binding', 'physiologically based', 'qualitative',
    'size separation', 'two competing mechanistic',
)

def txt(r, k):
    return (r.get(k) or '').strip()

def blocked(s):
    s = s.lower()
    return next((w for w in NOT_A_CONSTANT if w in s), None)

def verdict(rows):
    """Return (fits: bool, reason: str) for one source's albumin rows."""
    sp = {txt(r, 'protein_species').lower() for r in rows}
    if not (sp & SHEET_SPECIES) and not any(
            any(x in s for x in SHEET_SPECIES) for s in sp):
        return False, f"species not carried by the sheet ({'; '.join(sorted(sp))})"

    if all(txt(r, 'is_computational').lower() in ('true', 'yes', '1')
           for r in rows):
        return False, "entirely computational"

    usable = []
    bad_units = set()
    bad_method = set()
    for r in rows:
        meth = txt(r, 'method').lower()
        hit = next((w for w in NO_SOLUTION_KD if w in meth), None)
        if hit:
            bad_method.add(hit); continue
        kd, ka = txt(r, 'kd'), txt(r, 'ka')
        if not kd and not ka:
            continue
        b = blocked(txt(r, 'kd_units') if kd else txt(r, 'ka_units'))
        if b:
            bad_units.add(b); continue
        usable.append(r)

    if not usable:
        why = []
        if bad_units:
            why.append("the reported constant is not a dissociation constant "
                       f"({', '.join(sorted(bad_units))})")
        if bad_method:
            why.append(f"method yields no solution Kd ({', '.join(sorted(bad_method))})")
        if not why:
            why.append("no Kd and no Ka reported")
        return False, "; ".join(why)

    with_n = [r for r in usable if txt(r, 'n_binding_sites')]
    note = (f"{len(usable)} row(s) with a constant, "
            f"{len(with_n)} of them with n")
    if not with_n:
        return True, note + " -- CONVERTIBLE ONLY IF n = 1 IS ASSUMED"
    return True, note

def in_sheet():
    try:
        import openpyxl
    except ImportError:
        return None
    ws = openpyxl.load_workbook(XL)['Coefficients']
    hdr = [c.value for c in ws[1]]
    i = hdr.index('Source')
    return collections.Counter(
        str(r[i].value) for r in ws.iter_rows(min_row=2) if r[i].value)

def key(study, year):
    first = (study or '').replace('et al.', '').strip().split()
    return f"{first[0].rstrip(',') if first else '?'} {year}"

def main(argv):
    rows = list(csv.DictReader(open(DB, encoding='utf-8')))
    alb = [r for r in rows
           if 'albumin' in txt(r, 'protein').lower()
           or txt(r, 'protein').upper().startswith(('HSA', 'BSA'))]
    by = collections.defaultdict(list)
    for r in alb:
        by[key(txt(r, 'study'), txt(r, 'year'))].append(r)

    if argv:
        want = [a.lower() for a in argv]
        by = {k: v for k, v in by.items() if any(w in k.lower() for w in want)}

    sheet = in_sheet() or {}
    fits_new, fits_in, fails, decided = [], [], [], []
    for k, rs in by.items():
        ok, why = verdict(rs)
        n_sheet = sheet.get(k, 0)
        rec = (k, len(rs), n_sheet, why)
        if k in DECIDED:
            decided.append((k, len(rs), n_sheet, DECIDED[k]))
        elif ok and n_sheet:
            fits_in.append(rec)
        elif ok:
            fits_new.append(rec)
        else:
            fails.append(rec)

    def show(title, recs, note=''):
        print(f"\n{title}  ({len(recs)})")
        if note: print(f"  {note}")
        for k, n, ns, why in sorted(recs, key=lambda x: -x[1]):
            print(f"  {k:28s} {n:4d} db rows  {('sheet ' + str(ns)) if ns else '':10s} {why}")

    print("=" * 78)
    print("SCREEN: can this source reach log K [L water / kg protein]?")
    print("=" * 78)
    show("FITS, AND NOT YET IN THE SHEET -- the worklist", fits_new,
         "read these; each one can produce rows")
    show("CANNOT FIT -- do not read for the sheet", fails,
         "no route to log K; reading further adds nothing")
    show("FITS, ALREADY REPRESENTED", fits_in)
    show("ALREADY DECIDED -- not worklist", decided)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
