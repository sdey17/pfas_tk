#!/usr/bin/env python3
"""Does each reported Kd lie inside the concentration range tested?

A dissociation constant is a property of the curve where half the
protein is occupied.  If the titration never went near that
concentration, the reported Kd is an extrapolation, and how far an
extrapolation it is was not previously recorded anywhere in this
project.  This check records it for every row that carries both a
constant and a ligand range.

IT WAS WRITTEN AFTER a specific case.  Gao 2019 and Khazaee 2021 both
measured PFOA on human L-FABP by equilibrium dialysis and differ
1990-fold.  Gao's titration ran 2-12 uM against a reported Kd of
197 uM: it never came within twenty-fold of the constant it reports.
Khazaee's ran 0.05-5 uM against 0.099 uM, which brackets it.  Nobody
had to read either paper to see that -- it follows from the
conditions already in the database -- so it should be computed rather
than discovered one paper at a time.  See
audit_khazaee2021_lfabp.py for what it does and does not settle.

TWO SEPARATE THINGS ARE REPORTED, and they are not the same.

  range_flag   purely geometric: where the Kd sits relative to the
               ligand concentrations used.  Needs only the ligand
               range.
                 inside        the titration brackets the Kd
                 edge          within 3x of an edge
                 outside       3-30x beyond an edge
                 far outside   more than 30x beyond an edge

  basis        which identifiability test applies, decided by the
               method. A free-concentration shift is only meaningful
               for an assay that measures a free concentration, so
               separation methods get that, ITC gets its own Wiseman
               c value, and everything else gets the geometric flag
               alone rather than a number that does not apply to it.

  shift_pct    for a separation method, what the assay could see: the largest change
               in free ligand concentration, anywhere in the
               titration, that a ten-fold tighter Kd would produce.
               Needs the protein concentration too, and assumes
               one site.  A Kd that moves the observable by less than
               a typical method precision of about 15 % is not
               measurable by that titration whatever the fit says.
                 identifiable      shift > 36 %
                 marginal          shift 18-36 %
                 not identifiable  shift < 18 %

A row can be geometrically "outside" and still identifiable, or
"inside" and not, which is why both are kept.

WHAT THIS IS NOT.  It is not a verdict on a paper.  A Kd outside the
tested range can still be right; it is a statement about what that
experiment could constrain on its own.  Nothing is excluded on the
strength of this flag -- it is written into the shipped workbooks as
a column so a reader can weigh it.

Run:  python3 tk_review/scripts/kd_within_tested_range.py
Writes db/kd_range_check.csv
"""
import collections
import csv
import importlib.util
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "db", "kd_range_check.csv")

NUM = re.compile(r"^\s*([\d.]+(?:[eE][-+]?\d+)?)")
PRECISION = 18.0          # per cent; a typical whole-method RSD ceiling
FACTOR = 10.0             # the Kd perturbation the shift is measured against

FIELDS = ["study", "year", "chemical", "protein", "species", "method",
          "kd_uM", "kd_source", "ligand_lo_uM", "ligand_hi_uM",
          "protein_uM", "range_flag", "fold_outside", "basis",
          "shift_pct", "wiseman_c", "identifiability"]


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def kd_uM(r):
    """Kd in micromolar, from kd or from Ka, or (None, '')."""
    if r["kd"].strip():
        m = NUM.match(r["kd"])
        u = r["kd_units"].strip()
        if m:
            v = float(m.group(1))
            if u.startswith("uM"):
                return v, "Kd"
            if u.startswith("nM"):
                return v / 1e3, "Kd"
            if u.startswith("mM"):
                return v * 1e3, "Kd"
            if u.startswith("M") and not u.startswith("M^"):
                return v * 1e6, "Kd"
    if r["ka"].strip():
        m = NUM.match(r["ka"])
        u = r["ka_units"].strip()
        if m:
            v = float(m.group(1))
            if u.startswith("log Ka") or u.startswith("log10 M"):
                return 1e6 / (10 ** v), "1/Ka"
            if u.startswith("M^-1"):
                return 1e6 / v, "1/Ka"
    return None, ""


# conc_range_uM takes "2 to 10" but not "0.002, 0.004, 0.006" -- and a
# comma-separated list of the actual concentrations used is the most
# informative form a paper can give, so it should not be the one that
# falls through. This is unit arithmetic on an explicit list, the same
# arithmetic conc_range_uM does, and it is used only as a fallback.
LISTED = re.compile(r"\d*\.?\d+(?:[eE][-+]?\d+)?")


def listed_range(value, units, mw, B, MM):
    """(lo, hi) in uM from a comma-separated list of concentrations."""
    v, u = value.strip(), units.strip()
    if "," not in v or not u:
        return None, None
    head = u.split()[0]
    factor = B.TO_UM.get(head)
    if factor is None and mw:
        per = MM.MASS_UNITS.get(head)
        factor = None if per is None else per / float(mw)
    if factor is None:
        return None, None
    # Stop at the first run of text; "0.002, 0.004, 0.010 mM in PBS"
    # should not pick up a 7 out of a trailing sentence.
    head_txt = re.split(r"[;(]|\bin\b|\bper\b", v)[0]
    vals = [float(x) * factor for x in LISTED.findall(head_txt)]
    vals = [x for x in vals if x > 0]
    if len(vals) < 2:
        return None, None
    return min(vals), max(vals)


# How the method senses binding decides which identifiability test
# applies. The free-concentration shift is only meaningful where the
# assay actually measures a free concentration.
SEPARATION = re.compile(
    r"dialysis|ultrafiltration|ultracentrifug|centrifug|SPME|"
    r"solid.phase micro", re.I)
ITC = re.compile(r"isothermal titration|\bITC\b", re.I)


def free(L, P, Kd):
    """Free ligand from L = f + P f/(Kd+f), one site."""
    b = Kd + P - L
    return (-b + math.sqrt(b * b + 4.0 * Kd * L)) / 2.0


def classify(kd, lo, hi):
    if lo <= kd <= hi:
        return "inside", 1.0
    fold = kd / hi if kd > hi else lo / kd
    flag = ("edge" if fold <= 3 else
            "outside" if fold <= 30 else "far outside")
    return flag, fold


def shift(kd, lo, hi, P):
    """Largest free-concentration change a 10x tighter Kd would make,
    over the titration. Returns per cent, or None without a protein
    concentration."""
    if P is None or P <= 0:
        return None
    worst = 0.0
    n = 7
    for i in range(n):
        L = lo * (hi / lo) ** (i / (n - 1.0)) if hi > lo else lo
        f0, ft = free(L, P, kd), free(L, P, kd / FACTOR)
        if f0 > 0:
            worst = max(worst, abs(f0 - ft) / f0 * 100.0)
    return worst


def main():
    B = _load("B", "build_binding_affinity_table.py")
    MM = B.MM
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    out, no_kd, no_range = [], 0, 0
    for r in rows:
        kd, src = kd_uM(r)
        if kd is None or kd <= 0:
            no_kd += 1
            continue
        mw = MM.molar_mass(r["chemical"], r["n_perfluorinated_carbons"],
                           r["head_group"])[0]
        lo, hi = B.conc_range_uM(r["ligand_conc_range"],
                                 r["ligand_conc_units"], mw)
        if lo is None or lo <= 0:
            lo, hi = listed_range(r["ligand_conc_range"],
                                  r["ligand_conc_units"], mw, B, MM)
        if lo is None or lo <= 0:
            no_range += 1
            continue
        plo, phi = B.conc_range_uM(r["protein_conc"],
                                   r["protein_conc_units"])
        P = phi if phi else None
        flag, fold = classify(kd, lo, hi)

        # The free-concentration shift answers "could this assay see
        # the Kd", but only for an assay that measures a free
        # concentration. ITC measures heat, and its own criterion is
        # the Wiseman c value; everything else gets the geometric
        # flag alone rather than a number that does not apply to it.
        meth = r["method"]
        sh, cval, ident, basis = None, None, "", ""
        if SEPARATION.search(meth):
            basis = "free-conc shift"
            sh = shift(kd, lo, hi, P)
            ident = ("" if sh is None else
                     "identifiable" if sh > 2 * PRECISION else
                     "marginal" if sh > PRECISION else
                     "not identifiable")
        elif ITC.search(meth):
            basis = "Wiseman c = [P]/Kd"
            if P:
                c_val = P / kd
                ident = ("identifiable" if 1 <= c_val <= 1000 else
                         "marginal" if 0.3 <= c_val <= 3000 else
                         "not identifiable")
                cval = c_val
        out.append(dict(
            study=r["study"], year=r["year"], chemical=r["chemical"],
            protein=r["protein"][:46], species=r["protein_species"],
            method=r["method"][:60], kd_uM=f"{kd:.6g}", kd_source=src,
            ligand_lo_uM=f"{lo:.6g}", ligand_hi_uM=f"{hi:.6g}",
            protein_uM="" if P is None else f"{P:.6g}",
            basis=basis, range_flag=flag,
            fold_outside="" if flag == "inside" else f"{fold:.3g}",
            shift_pct="" if sh is None else f"{sh:.1f}",
            wiseman_c="" if cval is None else f"{cval:.3g}",
            identifiability=ident))

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: (r["study"], r["year"],
                                               r["chemical"])))

    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print(f"  {len(rows)} rows in the database")
    print(f"  {no_kd} with no usable constant, {no_range} with a "
          f"constant but no parseable ligand range")
    print(f"  {len(out)} CHECKED\n")

    c = collections.Counter(r["range_flag"] for r in out)
    print("  where the Kd sits relative to the concentrations tested")
    for k in ("inside", "edge", "outside", "far outside"):
        if c[k]:
            print(f"    {c[k]:5d}  {k}")
    bad = sum(v for k, v in c.items() if k != "inside")
    print(f"    {100.0*bad/len(out):.0f} % of checked rows report a Kd "
          f"the titration did not reach")

    i = collections.Counter(r["identifiability"] for r in out
                            if r["identifiability"])
    print(f"\n  of the {sum(i.values())} that also state a protein "
          f"concentration")
    for k in ("identifiable", "marginal", "not identifiable"):
        if i[k]:
            print(f"    {i[k]:5d}  {k}")

    worst = sorted((r for r in out if r["shift_pct"]),
                   key=lambda r: float(r["shift_pct"]))[:10]
    print("\n  SEPARATION METHODS, least identifiable. 'shift' is how")
    print("  much a ten-fold tighter Kd would move the free")
    print("  concentration; under about 18 % the assay cannot see it.")
    print(f"    {'study':22}{'chemical':12}{'Kd uM':>10}"
          f"{'ligand uM':>16}{'shift':>9}")
    for r in worst:
        rng = f"{r['ligand_lo_uM']}-{r['ligand_hi_uM']}"
        print(f"    {r['study'][:16]+' '+r['year']:22}{r['chemical'][:11]:12}"
              f"{r['kd_uM']:>10}{rng:>16}{r['shift_pct']+' %':>9}")

    itc = sorted((r for r in out if r["wiseman_c"]),
                 key=lambda r: float(r["wiseman_c"]))
    if itc:
        print("\n  ITC, by Wiseman c = [protein]/Kd. The isotherm has a")
        print("  usable inflection for roughly 1 < c < 1000; outside")
        print("  that its shape does not determine the Kd.")
        print(f"    {'study':22}{'chemical':12}{'Kd uM':>10}"
              f"{'[P] uM':>9}{'c':>10}  verdict")
        for r in itc[:6] + (itc[-3:] if len(itc) > 9 else []):
            print(f"    {r['study'][:16]+' '+r['year']:22}"
                  f"{r['chemical'][:11]:12}{r['kd_uM']:>10}"
                  f"{r['protein_uM']:>9}{r['wiseman_c']:>10}"
                  f"  {r['identifiability']}")


if __name__ == "__main__":
    main()
