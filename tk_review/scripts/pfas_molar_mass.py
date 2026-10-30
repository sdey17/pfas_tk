#!/usr/bin/env python3
"""Molar masses for the PFAS in this database, computed from formulae.

WHY. A concentration stated in ug/L cannot be turned into micromolar
without a molar mass, so build_binding_affinity_table.py refused all of
them -- 128 of the 1,241 binding rows, Fischer 2024's 111 and Khazaee
2021's 17.

HOW MANY THAT ACTUALLY FREES: 17, not 128. The molar mass was never the
only thing in the way, and for Fischer it was not in the way at all.
Fischer states its ligand concentration as "0.2-42 (C18-water step);
32-PFAS mixture at sub-saturating levels in protein steps" -- the number
belongs to the C18-water calibration step, and the protein measurements
are given no concentration at all. Those 111 rows stay blank, correctly,
and would have stayed blank with every molar mass in the world. Khazaee's
17 are the ones this file frees.

NOT QUOTED FROM MEMORY. Every mass here is computed from a molecular
formula, and every formula is derived from something the database or the
source paper already states:

  * the regular families follow a RULE applied to `n_perfluorinated_carbons`
    and `head_group`, two columns db/protein_binding.csv already carries
    for each row, so nothing new is asserted about the chemistry;

  * the irregular ones are read off the IUPAC names Fischer et al. 2024
    prints in its own compound table, quoted in EXCEPTIONS below. That
    table is also what settles what "FDSA" means -- Perfluoro-1-
    decanesulfonamide, so the sulfonamide rule covers it and it needs no
    exception -- and what separates the three sulfonamidoacetic acids from
    the plain sulfonamides this database files them with.

THE RULES ARE CHECKED, NOT TRUSTED. --check recomputes the 17 rule-covered
compounds for which this repository already holds a molar mass from an
independent source (db/qsar/qsar_endpoint_table.csv, db/ryu2024_funbound.csv
and db/mbm_qin_vs_armitage.csv) and prints the agreement. Sixteen land
within 0.01 g/mol, which is what a different set of atomic weights and a
two-decimal rounding give; the seventeenth is out by 0.18 and is left
standing as a disagreement rather than reconciled: see FTS_NOTE.

Run:  python3 tk_review/scripts/pfas_molar_mass.py --check
"""
import argparse
import csv
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db")

# IUPAC conventional atomic weights.
ATOMIC = {"C": 12.011, "H": 1.008, "F": 18.998403163,
          "N": 14.007, "O": 15.999, "S": 32.06}

# Mass-per-volume units that appear in db/protein_binding.csv. ug/L and
# ng/mL are the same number, which is why both convert with the same
# arithmetic: uM = (ug/L) / (g/mol).
MASS_UNITS = {"ug/L": 1.0, "ng/mL": 1.0, "mg/L": 1e3, "ug/mL": 1e3,
              "ng/L": 1e-3, "mg/mL": 1e6}

# n is n_perfluorinated_carbons, the column the database already carries.
# Each rule returns the FREE ACID, which is the form every molar mass in
# this repository is quoted on (PFOA 414.07, not the 413.06 anion).
RULES = {
    # CF3(CF2)n-1COOH -- the carboxyl carbon is not perfluorinated, so C is
    # n+1 while F stays 2n+1.
    "carboxylate":
        lambda n: {"C": n + 1, "H": 1, "F": 2 * n + 1, "O": 2},
    # CF3(CF2)n-1SO3H
    "sulfonate":
        lambda n: {"C": n, "H": 1, "F": 2 * n + 1, "S": 1, "O": 3},
    # CF3(CF2)n-1SO2NH2
    "sulfonamide":
        lambda n: {"C": n, "H": 2, "F": 2 * n + 1, "N": 1, "S": 1, "O": 2},
    # F(CF2)nCH2CH2SO3H -- the two telomer carbons carry four H, the acid one
    "fluorotelomer sulfonate":
        lambda n: {"C": n + 2, "H": 5, "F": 2 * n + 1, "S": 1, "O": 3},
    # F(CF2)nCH2CH2COOH
    "fluorotelomer carboxylate":
        lambda n: {"C": n + 3, "H": 5, "F": 2 * n + 1, "O": 2},
}

# Compounds no rule reaches, each with the name it was derived from.
EXCEPTIONS = {
    # Fischer 2024 files these under "Perfluoroalkyl Sulfonamidoacetic
    # Acids", a head group db/protein_binding.csv does not distinguish from
    # the plain sulfonamides -- so the rule would return the sulfonamide and
    # be wrong by the whole acetic acid.
    "FOSAA": ({"C": 10, "H": 4, "F": 17, "N": 1, "O": 4, "S": 1},
              "Perfluoro-1-octanesulfonamidoacetic acid (Fischer 2024 "
              "compound table): C8F17-SO2-NH-CH2-COOH"),
    "N-MeFOSAA": ({"C": 11, "H": 6, "F": 17, "N": 1, "O": 4, "S": 1},
                  "N-methyl perfluoro-1-octanesulfonamidoacetic acid "
                  "(Fischer 2024): FOSAA with CH3 on the amide nitrogen"),
    "N-EtFOSAA": ({"C": 12, "H": 8, "F": 17, "N": 1, "O": 4, "S": 1},
                  "N-ethyl perfluoro-1-octanesulfonamidoacetic acid "
                  "(Fischer 2024): FOSAA with C2H5 on the amide nitrogen"),
    # Fischer 2024 names it "Dodecafluoro-3H-4,8-dioxanonanoate". Read
    # straight off that: a nonanoate chain, so nine positions, of which 4
    # and 8 are oxygen -- seven carbons and two ether oxygens -- twelve
    # fluorines, one hydrogen at position 3, two carboxyl oxygens, and the
    # acidic hydrogen the free acid adds to the named anion.
    "ADONA": ({"C": 7, "H": 2, "F": 12, "O": 4},
              "Dodecafluoro-3H-4,8-dioxanonanoate (Fischer 2024), as the "
              "free acid"),
}

# How close two molar masses have to be to count as the same number. The
# held values are printed to two decimals and were computed with some
# other set of atomic weights, so a disagreement below this is bookkeeping
# rather than chemistry: a formula error in any of these families would
# show up as a whole CF2 (50 g/mol), a hydrogen (1), or an oxygen (16) --
# never as 0.01.
SAME = 0.015

# The one disagreement, left standing rather than reconciled.
FTS_NOTE = (
    "6:2 FTSA computes to 428.16 from C8H5F13O3S; db/ryu2024_funbound.csv "
    "carries 427.98 for the same compound. The gap is 0.18 g/mol, 0.04%, "
    "which is immaterial for a ug/L to uM conversion -- a hundredth of the "
    "width of the narrowest effect anywhere in this project -- but it is a "
    "real disagreement and not a rounding artifact, so it is recorded here "
    "rather than papered over. The formula is the one Fischer's own naming "
    "implies and is used; Ryu's value is not overwritten in its own file."
)

_ns = importlib.util.spec_from_file_location(
    "pfas_names", os.path.join(HERE, "pfas_names.py"))
N = importlib.util.module_from_spec(_ns)
_ns.loader.exec_module(N)


def mass_of(counts):
    return sum(ATOMIC[el] * k for el, k in counts.items())


def formula_string(counts):
    parts = []
    for el in ("C", "H", "F", "N", "O", "S"):
        k = counts.get(el, 0)
        if k:
            parts.append(f"{el}{k if k > 1 else ''}")
    return "".join(parts)


def molar_mass(name, n_cf="", head=""):
    """(g/mol, formula, how it was reached) or (None, None, why not)."""
    canon = N.canon(name)
    if canon in EXCEPTIONS:
        counts, why = EXCEPTIONS[canon]
        return mass_of(counts), formula_string(counts), why
    rule = RULES.get((head or "").strip().lower())
    if rule is None:
        return None, None, (f"no rule for head group {head!r}, and no "
                            f"formula recorded for {canon}")
    try:
        n = int(float(str(n_cf).strip()))
    except (TypeError, ValueError):
        return None, None, (f"{canon} has head group {head!r} but no "
                            f"perfluorinated-carbon count to apply it to")
    counts = rule(n)
    return (mass_of(counts), formula_string(counts),
            f"rule for {head} at n_perfluorinated_carbons = {n}")


# ------------------------------------------------------------------ check
KNOWN = [("db/qsar/qsar_endpoint_table.csv", "chemical", "mw_g_mol"),
         ("db/ryu2024_funbound.csv", "chemical", "mw_g_mol"),
         ("db/mbm_qin_vs_armitage.csv", "compound", "MW")]


def known_masses():
    out = {}
    for path, namecol, mwcol in KNOWN:
        with open(os.path.join(DB, os.path.relpath(path, "db"))) as fh:
            for r in csv.DictReader(fh):
                v = (r.get(mwcol) or "").strip()
                if not v:
                    continue
                try:
                    out.setdefault(N.canon(r[namecol]), []).append(
                        (float(v), os.path.basename(path)))
                except ValueError:
                    pass
    return out


def descriptors():
    """n_perfluorinated_carbons and head_group, as the database records them."""
    out = {}
    with open(os.path.join(DB, "protein_binding.csv")) as fh:
        for r in csv.DictReader(fh):
            c = N.canon(r["chemical"])
            if c not in out and r["head_group"].strip():
                out[c] = (r["n_perfluorinated_carbons"], r["head_group"])
    return out


def check():
    desc, known = descriptors(), known_masses()
    print("PFAS MOLAR MASSES")
    print("=" * 74)
    print("Every mass computed from a formula; every formula from the "
          "database's own\nn_perfluorinated_carbons and head_group, or from "
          "a name the source paper prints.")
    print()
    print("Against the molar masses this repository already holds, from "
          "three files\nwritten independently of this one:")
    print(f"  {'compound':12s}{'computed':>10s}{'held':>10s}{'diff':>9s}  "
          f"formula")
    agree = disagree = 0
    for c in sorted(known):
        if c not in desc:
            continue
        mw, formula, _ = molar_mass(c, *desc[c])
        if mw is None:
            continue
        held = known[c][0][0]
        d = mw - held
        flag = "" if abs(d) < SAME else "   <-- DISAGREES"
        agree, disagree = (agree + 1, disagree) if not flag else (agree,
                                                                  disagree + 1)
        print(f"  {c:12s}{mw:10.2f}{held:10.2f}{d:+9.3f}  {formula}{flag}")
    print()
    print(f"  {agree} agree within {SAME} g/mol -- the width of a different "
          f"atomic-weight set,")
    print(f"  not of a formula error, which would show as a whole CF2, H or "
          f"O. {disagree} do not.")
    if disagree:
        print()
        for line in (FTS_NOTE[i:i + 70] for i in range(0, len(FTS_NOTE), 70)):
            print("  " + line)
    print()
    # Coverage is only interesting where a molar mass is actually needed:
    # a row whose ligand concentration is stated in mass units. The database
    # holds plenty of betaines, diols and compound MIXTURES that no rule
    # reaches and that no row asks about.
    need, rows_needing = {}, 0
    with open(os.path.join(DB, "protein_binding.csv")) as fh:
        for r in csv.DictReader(fh):
            u = (r["ligand_conc_units"] or "").strip().split()
            if u and u[0] in MASS_UNITS:
                rows_needing += 1
                need.setdefault(N.canon(r["chemical"]),
                                (r["n_perfluorinated_carbons"],
                                 r["head_group"]))
    reached = {c: molar_mass(c, *d) for c, d in need.items()}
    got = [c for c, (m, _, _) in reached.items() if m is not None]
    miss = [(c, w) for c, (m, _, w) in reached.items() if m is None]
    print(f"WHERE A MOLAR MASS IS ACTUALLY NEEDED -- {rows_needing} rows "
          f"state the ligand")
    print(f"concentration in mass units, over {len(need)} compounds:")
    print(f"  {len(got)} compounds reached, {len(miss)} not.")
    for c, why in sorted(miss):
        print(f"    {c:16s} {why}")
    # Having a molar mass is necessary and not sufficient: the value string
    # still has to parse. Report both, because the gap between them is the
    # whole Fischer story and quoting the first number alone would overstate
    # what this file achieved.
    import importlib.util as _iu
    _b = _iu.spec_from_file_location(
        "baf", os.path.join(HERE, "build_binding_affinity_table.py"))
    B = _iu.module_from_spec(_b)
    _b.loader.exec_module(B)
    covered = converts = 0
    blocked = {}
    with open(os.path.join(DB, "protein_binding.csv")) as fh:
        for r in csv.DictReader(fh):
            u = (r["ligand_conc_units"] or "").strip().split()
            if not (u and u[0] in MASS_UNITS):
                continue
            c = N.canon(r["chemical"])
            if c not in got:
                continue
            covered += 1
            mw = molar_mass(c, *need[c])[0]
            if B.conc_range_uM(r["ligand_conc_range"],
                               r["ligand_conc_units"], mw)[1] is not None:
                converts += 1
            else:
                blocked.setdefault(r["study"], r["ligand_conc_range"])
    print(f"  {covered} of {rows_needing} rows now have a molar mass, and "
          f"{converts} of them convert.")
    if blocked:
        print()
        print("  The rest are blocked by the concentration STRING, not by "
              "the molar mass, and")
        print("  a molar mass was never going to free them:")
        for study, raw in sorted(blocked.items()):
            print(f"    {study}: {raw[:96]}")
        print("    There the number belongs to a different step of the "
              "experiment; the binding")
        print("    measurements themselves are given no concentration at "
              "all, so a blank is the")
        print("    only honest cell.")
    print()
    print("The database also holds betaines, diols, acrylates, neutral "
          "ethers and rows naming\na MIXTURE of compounds. No rule reaches "
          "those and none of them states a mass\nconcentration, so they are "
          "not a gap this file has to close.")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    raise SystemExit(check())
