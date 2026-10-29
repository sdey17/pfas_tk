#!/usr/bin/env python3
"""One name per compound, for the whole database.

The same PFAS is written differently by different groups. Left alone it
inflates every compound count and, worse, invents coverage gaps: a paper
reporting a partition coefficient for "PFUnA" looks like a compound with
no affinity data when PFUnDA has plenty. That error made the
reverse-conversion gap look like ten compounds when it is three.

This module is the single place that decides compound identity. It does
NOT rewrite db/protein_binding.csv: the raw `chemical` string stays as
the source published it, because that is provenance. Reconciliation is a
layer on top, auditable and reversible.

EVERY GROUP BELOW WAS CHECKED against the database's own
n_perfluorinated_carbons, head_group and source fields. Three rules were
used, in this order:

  1. If ONE source uses BOTH names, they are different compounds. This
     is what separates N-TAmP-FHxSA from PFHxSaAm: Chen et al. 2025
     reports both, at 33.68 and 6.19 uM.
  2. If the names share nCF and head group and no source uses both,
     they are candidates -- but only merged when the chemistry is
     unambiguous (an N- prefix, a trade name in brackets, a case
     difference, a known abbreviation pair).
  3. Anything left over goes in UNRESOLVED and is NOT merged.

Run:  python3 tk_review/scripts/pfas_names.py        # audit the database
"""
import collections
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

# canonical -> the spellings that mean it. The canonical form is the one
# this project already uses most, so existing output does not shift.
ALIASES = {
    # --- trade name or systematic name carried in brackets -------------
    "HFPO-DA":      ("HFPO-DA", "HFPO-DA (GenX)"),
    "ADONA":        ("ADONA", "NaDONA (ADONA)"),
    "6:2 Cl-PFESA": ("6:2 Cl-PFESA", "9Cl-PF3ONS", "9Cl-PF3ONS (F-53B)",
                     "9Cl-PF3ONS (F-53B, 6:2 Cl-PFESA)"),
    "8:2 Cl-PFESA": ("8:2 Cl-PFESA", "11Cl-PF3OUdS",
                     "11Cl-PF3OUdS (8:2 Cl-PFESA)"),
    # --- two accepted abbreviations for one perfluoroalkyl acid --------
    # nCF and head group agree in every case; no source uses both.
    "PFUnDA":       ("PFUnA", "PFUnDA"),            # nCF 10, carboxylate
    "PFDoDA":       ("PFDoA", "PFDoDA"),            # nCF 11, carboxylate
    "PFTrDA":       ("PFTrDA", "PFTriDA"),          # nCF 12, carboxylate
    "PFOSA":        ("FOSA", "PFOSA"),              # nCF 8, sulfonamide
    "PFECHS":       ("PFECHS", "PFEtCHxS"),         # nCF 8, cyclic sulfonate
    # --- the N- prefix, present in one source and absent in another ----
    # Chen et al. 2025 writes the N- form with nCF and head group;
    # Muensterman et al. 2026 writes the short form with neither.
    "N-EtFOSA":     ("EtFOSA", "NEtFOSA"),
    "N-MeFOSA":     ("MeFOSA", "NMeFOSA"),
    "N-EtFOSAA":    ("EtFOSAA", "N-EtFOSAA"),
    "N-MeFOSAA":    ("MeFOSAA", "N-MeFOSAA"),
    "N-EtFOSE":     ("EtFOSE", "NEtFOSE"),
    "N-MeFOSE":     ("MeFOSE", "NMeFOSE"),
    # --- transposed letters in the unsaturated telomer acids -----------
    "6:2 FTUCA":    ("6:2 FTUCA", "6:2 UFTCA"),
    "8:2 FTUCA":    ("8:2 FTUCA", "8:2 UFTCA"),
    # --- capitalisation and a parenthetical chain length ---------------
    "Octanoic acid": ("Octanoic Acid", "Octanoic acid", "octanoic acid (C8)"),
    "Decanoic acid": ("Decanoic Acid", "Decanoic acid"),
}
_CANON = {v: k for k, vs in ALIASES.items() for v in vs}

# The fluorotelomer sulfonates need a pattern, not a list: FTS, FtS and
# FTSA are one compound but 4:2, 6:2, 8:2 and 10:2 are four of them, so
# the chain length is kept and only the spelling is dropped.
_FTS = re.compile(r"^(\d+):2[ ]?F[Tt]S(A)?$")

# NOT merged, and each for a stated reason.
ISOMER_VARIANTS = {
    # A linear isomer is not the technical mixture and neither is the
    # branched fraction. Merging them would average three different
    # measurements of three different things.
    "PFOA": ("linear PFOA (n-PFOA)", "branched PFOA isomers (range)"),
    "PFOS": ("linear PFOS (n-PFOS)", "branched PFOS isomers (range)"),
}
UNRESOLVED = {
    # Muensterman's "Cl-PFOS" carries no nCF and no head group, and no
    # source uses both spellings, so it cannot be shown to be Chen's
    # 8Cl-PFOS rather than another chlorinated isomer. Left apart.
    ("Cl-PFOS", "8Cl-PFOS"): "Muensterman gives no nCF or head group; "
                             "the chlorine position cannot be confirmed",
    # Chen et al. 2025 reports BOTH, at 33.68 and 6.19 uM. Same nCF and
    # head group, different compounds. Rule 1.
    ("N-TAmP-FHxSA", "PFHxSaAm"): "one source reports both, with "
                                  "different constants",
}


def canon(name):
    """The canonical spelling of one compound name."""
    n = (name or "").strip()
    m = _FTS.match(n)
    if m:
        return f"{m.group(1)}:2 FTSA"
    return _CANON.get(n, n)


def is_single_compound(name):
    """False for prose lists, ranges and narrative entries."""
    n = (name or "").strip()
    return bool(n) and "," not in n and " and " not in n and len(n) < 32


# ---------------------------------------------------------------- audit
def candidates(rows):
    """Alias candidates the table does not yet cover.

    Two generators: names that collapse on case and bracket stripping,
    and names sharing nCF with head group. Rule 1 removes any group a
    single source uses more than one member of.
    """
    singles = [r for r in rows if is_single_compound(r["chemical"])]
    byname = collections.defaultdict(set)
    bykey = collections.defaultdict(set)
    srcs = collections.defaultdict(set)
    for r in singles:
        n = r["chemical"].strip()
        srcs[n].add(f"{r['study']} {r['year']}")
        flat = re.sub(r"\s*\(.*?\)\s*", "", n).upper().replace(" ", "")
        byname[flat].add(n)
        v, h = r["n_perfluorinated_carbons"].strip(), r["head_group"].strip()
        if v.isdigit() and h:
            bykey[(int(v), h.lower())].add(n)

    out = []
    for how, groups in (("spelling", byname), ("nCF+head", bykey)):
        for k, v in groups.items():
            if len(v) < 2:
                continue
            if len({canon(x) for x in v}) == 1:
                continue                       # already reconciled
            pair = tuple(sorted(v))
            if any(set(p) <= set(pair) for p in UNRESOLVED):
                continue
            shared = any(len({x for x in v if s in srcs[x]}) > 1
                         for s in set().union(*(srcs[x] for x in v)))
            out.append((how, k, sorted(v), shared))
    return out


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    names = {r["chemical"].strip() for r in rows}
    singles = {n for n in names if is_single_compound(n)}
    print(f"{len(rows)} rows, {len(names)} distinct chemical strings, "
          f"{len(singles)} single compounds")
    print(f"after reconciliation: {len({canon(n) for n in singles})} "
          f"compounds "
          f"({len(singles) - len({canon(n) for n in singles})} merged)")
    print()
    print(f"{len(ALIASES)} explicit groups + the n:2 FTS pattern:")
    for k, v in sorted(ALIASES.items()):
        used = [x for x in v if x in names]
        if len(used) > 1:
            print(f"  {k:16s} <- {', '.join(used)}")
    fts = sorted(n for n in singles if _FTS.match(n))
    if fts:
        print(f"  {'n:2 FTSA':16s} <- {', '.join(fts)}")
    print()
    print("DELIBERATELY NOT MERGED:")
    for base, v in ISOMER_VARIANTS.items():
        print(f"  {base} vs {', '.join(v)}")
        print("      an isomer fraction is not the technical mixture")
    for pair, why in UNRESOLVED.items():
        print(f"  {' vs '.join(pair)}")
        print(f"      {why}")
    print()
    cand = candidates(rows)
    if not cand:
        print("AUDIT: no unreconciled alias candidates.")
        return 0
    print(f"AUDIT: {len(cand)} candidate group(s) not covered by the table.")
    print("  'one source uses several' means rule 1 says they are")
    print("  different compounds, which is the usual reason a group is")
    print("  listed here and left alone.")
    for how, k, v, shared in sorted(cand, key=lambda c: str(c[1])):
        tag = "  <- one source uses several" if shared else ""
        print(f"    [{how:8s}] {k if how == 'nCF+head' else ''} "
              f"{v}{tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
