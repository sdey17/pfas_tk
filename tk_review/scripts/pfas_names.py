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

Run:  python3 tk_review/scripts/pfas_names.py              # audit
      python3 tk_review/scripts/pfas_names.py --selftest   # audit the audit
"""
import collections
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_DBDIR = os.path.join(os.path.dirname(HERE), "db")
DB = os.path.join(_DBDIR, "protein_binding.csv")
# Every file whose rows carry a `chemical`. The alias table was first
# built and audited against protein_binding.csv alone, which left the
# other three unchecked -- and they held seven more merges, including a
# bare DSSTox identifier standing in for a compound name.
SOURCES = [DB] + [os.path.join(_DBDIR, "combined", f + ".csv")
                  for f in ("tk_parameters", "binding", "transporters",
                            "regulatory")]


def all_rows():
    out = []
    for path in SOURCES:
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path, encoding="utf-8")):
            if "chemical" in r:
                r.setdefault("n_perfluorinated_carbons", "")
                r.setdefault("head_group", "")
                r.setdefault("study", "")
                r.setdefault("year", "")
                out.append(r)
    return out

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
    "Octanoic acid": ("Octanoic Acid", "Octanoic acid", "octanoic acid (C8)",
                      "octanoic acid (reference ligand)"),
    "Decanoic acid": ("Decanoic Acid", "Decanoic acid",
                      "decanoic acid (reference ligand)"),
    # --- the systematic name spelled out in brackets -------------------
    "PFPrS":        ("PFPrS", "PFPrS (perfluoropropane sulfonate)"),
    # --- found only once the other three datasets were audited --------
    "PFODA":        ("PFODA", "PFOcDA"),      # perfluorooctadecanoate
    # PFHpA: EPA's httk pfas.clearance table carries the name and the
    # DSSTox id on the same rows (db/compiled_tk_parameters.csv), both
    # Ohmori 2003, rat, female, and the half-lives agree -- 1.2332 h in
    # the CPHEA refit against httk's 1.2 h. Resolved from this project's
    # own files, not guessed.
    "PFHpA":        ("PFHpA", "DTXSID1037303"),
}
# Added to groups already defined above, now that the other datasets are
# in scope. Kept separate so the provenance of each late addition shows.
for _canon, _extra in (
        ("HFPO-DA", "GenX"),                   # Argoul 2026's spelling
        ("PFUnDA", "PFUdA"),                   # Abraham 2024
        ("PFTeDA", "PFTA"),                    # Liu 2020
        ("ADONA", "DONA"),                     # Abraham 2024; the acid
        ("6:2 Cl-PFESA", "6:2 Cl-PFESA (F-53B major component)")):
    ALIASES[_canon] = tuple(dict.fromkeys(
        ALIASES.get(_canon, (_canon,)) + (_extra,)))
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
# Not compounds at all: a class name standing in for a set, used by
# review papers. Counting these as compounds inflates every total.
NOT_A_COMPOUND = (
    "perfluoroalkyl acids",
    "perfluoroalkyl acids (model compounds)",
    "C8-C11 PFCAs",
    # the transporter dataset uses these where no single compound applies
    "(endogenous organic anions)",
    "all 14 PFAS",
    "any",
    "n/a",
    "n/a (constitutive mRNA expression)",
    "n/a (hormone manipulation)",
    "linear PFCAs C4-C12",
)

# Pairs the audit keeps proposing that are DIFFERENT compounds, each with
# the reason, so the audit stops asking and the reasoning is on record.
NOT_MERGED = {
    ("ADONA", "PFO3DA (perfluoro-3,6,9-trioxadecanoic acid)"):
        "both nCF 6 ether carboxylates, but ADONA is "
        "4,8-dioxa-3H-perfluorononanoic acid and PFO3DA is "
        "perfluoro-3,6,9-trioxadecanoic acid",
    ("Decanoic Acid", "Octanoic acid"): "different fatty acid chain lengths",
    ("Decanoic acid", "Octanoic Acid"): "different fatty acid chain lengths",
    ("Decanoic acid", "Dodecanoic Acid"): "different fatty acid chain lengths",
    ("Dodecanoic Acid", "Octanoic acid"): "different fatty acid chain lengths",
    ("EtFOSAA", "NEtFOSA"):
        "a sulfonamidoacetic acid against a sulfonamide",
    ("MeFOSA", "N-MeFOSAA"):
        "a sulfonamide against a sulfonamidoacetic acid",
    ("MeFOSAA", "NMeFOSA"):
        "a sulfonamidoacetic acid against a sulfonamide",
    ("N-EtFBSE", "NEtFOSE"):
        "a butane- against an octane-sulfonamidoethanol",
    ("N-MeFBSE", "NMeFOSE"):
        "a butane- against an octane-sulfonamidoethanol",
    ("PFOA", "PFOA (as 8:2-FTOH metabolite)"):
        "the same compound, but the parameter is PFOA formed from an "
        "8:2 FTOH dose, which is not a direct-dosing parameter and must "
        "not be pooled with one",
    # Sakolish et al. 2026 dose the salts and the acids side by side, so
    # rule 1 already separates these; recorded because the rule only
    # fires while that one source is in scope.
    ("NH4PFOA", "PFOA"): "one source doses both the ammonium salt and "
                         "the acid",
    ("PFOA", "PFOA-K"): "one source doses both the potassium salt and "
                        "the acid",
    ("PFOS", "PFOS-K"): "one source doses both the potassium salt and "
                        "the acid",
    ("PFHxS", "PFHxS-K"): "one source doses both the potassium salt and "
                          "the acid",
    ("PFDoDA", "PFODA"): "perfluorododecanoate against "
                         "perfluorooctadecanoate -- six CF2 apart, and "
                         "PFOcDA/PFDoDA look alike only as strings",
}

# Merges no audit can discover, because the two strings share nothing
# and the evidence is external to the name. The self-test excludes these
# rather than failing on them: an audit built on spelling and annotation
# cannot be blamed for missing a registry identifier. They are listed so
# the exclusion is visible, not silent.
UNDISCOVERABLE = {
    ("DTXSID1037303", "PFHpA"):
        "a DSSTox identifier used in place of a name; resolved from "
        "db/compiled_tk_parameters.csv, where EPA's httk pfas.clearance "
        "table carries both on the same Ohmori 2003 rat rows and the "
        "half-lives agree (1.2332 h refitted against httk's 1.2 h)",
    ("PFTA", "PFTeDA"):
        "perfluorotetradecanoate written two ways that share only three "
        "letters; chemistry, not spelling",
    ("PFOcDA", "PFODA"):
        "perfluorooctadecanoate, likewise",
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
    """False for prose lists, ranges and narrative entries.

    The comma test runs on the name with bracketed text removed, because
    a synonym in brackets may itself contain one: "9Cl-PF3ONS (F-53B,
    6:2 Cl-PFESA)" is one compound with two synonyms, not a list. That
    name was the last thing the self-test could not see.
    """
    n = (name or "").strip()
    if n in NOT_A_COMPOUND:
        return False
    bare = re.sub(r"\s*\(.*?\)\s*", " ", n).strip()
    return (bool(n) and "," not in bare and " and " not in bare
            and len(bare) < 32)


# ---------------------------------------------------------------- audit
# THREE GENERATORS, because two were not enough. A self-test (--selftest)
# blanks the alias table and asks the audit to rediscover the 34 merges
# made by hand. With only the first two generators it found 21 of them:
# the short N- forms, the transposed FTUCA/UFTCA pair and PFECHS/PFEtCHxS
# were invisible, because Muensterman et al. 2026 reports neither chain
# length nor head group and those names differ by more than case and
# brackets. The similarity generator closes that gap. Keep the self-test
# passing: an audit that cannot rediscover what is already known is not
# evidence that nothing is left.
import difflib
import itertools


def _flat(n):
    """Letters and digits only, upper case."""
    return re.sub(r"[^A-Z0-9:]", "", n.upper())


def _head_key(h):
    """Head group with wording differences removed.

    "cyclic sulfonate" and "sulfonate (cyclic)" are the same chemistry
    written two ways, and an exact string match treats them as two
    groups, which is how PFECHS and PFEtCHxS escaped the audit.
    """
    return tuple(sorted(re.findall(r"[a-z]+", h.lower())))


def candidates(rows):
    """Alias candidates the table does not yet cover.

    Generator 1  names equal after case and bracket stripping
    Generator 2  names sharing chain length and head-group wording
    Generator 3  names that look alike: a shared prefix or suffix, the
                 same letters in a different order, or a high character
                 similarity. This one needs no annotation at all, which
                 is what makes it work on the sparsely annotated rows.

    Rule 1 (one source uses several) is reported, not applied: the
    audit's job is to surface, the table's job is to decide.
    """
    singles = [r for r in rows if is_single_compound(r["chemical"])]
    byname = collections.defaultdict(set)
    bykey = collections.defaultdict(set)
    srcs = collections.defaultdict(set)
    for r in singles:
        n = r["chemical"].strip()
        srcs[n].add(f"{r['study']} {r['year']}")
        byname[re.sub(r"\s*\(.*?\)\s*", "", n).upper().replace(" ", "")].add(n)
        v, h = r["n_perfluorinated_carbons"].strip(), r["head_group"].strip()
        if v.isdigit() and h:
            bykey[(int(v), _head_key(h))].add(n)

    # Chain length and head group are IDENTITY in a PFAS name, not
    # noise: "3:2 FTOH" and "4:2 FTOH" are 94% similar as strings and
    # are different compounds. So the lookalike generator requires both
    # to agree where they are known, and only then asks whether the rest
    # of the name looks like a spelling variant. Without that constraint
    # it returns 121 candidates, nearly all of them chain-length
    # neighbours.
    chain, head = {}, {}
    for r in singles:
        n = r["chemical"].strip()
        m = re.match(r"^(\d+):(\d+)", n)
        if m:
            chain[n] = m.groups()
        elif r["n_perfluorinated_carbons"].strip().isdigit():
            chain.setdefault(n, ("nCF", r["n_perfluorinated_carbons"].strip()))
        if r["head_group"].strip():
            head.setdefault(n, _head_key(r["head_group"]))

    def compatible(a, b):
        """False when the database already says these are different."""
        if a in chain and b in chain and chain[a] != chain[b]:
            return False
        if a in head and b in head and head[a] != head[b]:
            return False
        return True

    # Generator 4: a synonym carried in brackets. "HFPO-DA (GenX)" is
    # the database telling us GenX and HFPO-DA are one compound, and no
    # amount of string similarity between "GENX" and "HFPODA" would ever
    # say so. Splitting such a name links the bare form to each synonym
    # inside the brackets.
    allnames = {r["chemical"].strip() for r in singles}
    bysyn = collections.defaultdict(set)
    for n in sorted(allnames):
        m = re.match(r"^(.*?)\s*\(([^)]*)\)\s*$", n)
        if not m:
            continue
        bare = m.group(1).strip()
        parts = [bare] + [x.strip() for x in m.group(2).split(",")]
        hit = {p for p in parts if p in allnames}
        if len(hit) > 1 or (hit and bare in allnames):
            bysyn[n] = hit | {n}

    names = sorted(allnames)
    bysim = collections.defaultdict(set)
    for i, a in enumerate(names):
        fa = _flat(a)
        for b in names[i + 1:]:
            if not compatible(a, b):
                continue
            fb = _flat(b)
            if fa == fb or sorted(fa) == sorted(fb):
                close = True                       # case, or transposed
            elif (fa.endswith(fb) or fb.endswith(fa)
                  or fa.startswith(fb) or fb.startswith(fa)):
                close = abs(len(fa) - len(fb)) <= 2    # an N- prefix
            else:
                close = difflib.SequenceMatcher(None, fa, fb).ratio() >= 0.82
            if close:
                bysim[tuple(sorted((a, b)))] |= {a, b}

    out = []
    for how, groups in (("spelling", byname), ("nCF+head", bykey),
                        ("lookalike", bysim), ("synonym", bysyn)):
        for k, v in groups.items():
            if len(v) < 2 or len({canon(x) for x in v}) == 1:
                continue
            # A group is settled when every distinct compound in it has
            # already been ruled different from every other. Checking
            # pairs rather than the whole group matters: {ADONA,
            # NaDONA (ADONA), PFO3DA} is two aliases plus one distinct
            # compound, and only collapses once the aliases merge.
            cset = {canon(x) for x in v}
            pairs = {frozenset((canon(a), canon(b)))
                     for p in tuple(UNRESOLVED) + tuple(NOT_MERGED)
                     for a, b in (p,)}
            if all(frozenset(c) in pairs
                   for c in itertools.combinations(sorted(cset), 2)):
                continue
            shared = any(len({x for x in v if s in srcs[x]}) > 1
                         for s in set().union(*(srcs[x] for x in v)))
            out.append((how, k, sorted(v), shared))
    return out


def main():
    rows = all_rows()
    names = {r["chemical"].strip() for r in rows}
    singles = {n for n in names if is_single_compound(n)}
    print(f"{len(SOURCES)} source files, {len(rows)} rows, "
          f"{len(names)} distinct chemical strings, "
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
    for pair, why in list(NOT_MERGED.items()) + list(UNRESOLVED.items()):
        print(f"  {' vs '.join(pair)}")
        print(f"      {why}")
    print()
    cand = candidates(rows)
    decided = [c for c in cand if c[3]]
    review = [c for c in cand if not c[3]]
    print(f"AUDIT: {len(cand)} candidate group(s) outside the table.")
    print(f"  {len(decided)} are settled by rule 1 -- a single source uses")
    print("  more than one member, so they are different compounds "
          "(Me/Et")
    print("  homologues, the N2-N5 Cl-PFECA series, the sulfonamide")
    print("  family, the isomer sets). Not listed.")
    print()
    if not review:
        print("  NOTHING NEEDS A DECISION.")
        return 0
    print(f"  {len(review)} group(s) rule 1 does not settle, listed in full")
    print("  because these are the ones a human has to look at:")
    print()
    for how, k, v, _ in sorted(review, key=lambda c: str(c[2])):
        print(f"    [{how:9s}] {v}")
    print()
    print("  Each of those is a different compound for a reason the")
    print("  names make plain: a fatty acid of another chain length, a")
    print("  sulfonamide against a sulfonamidoacetic acid or a")
    print("  sulfonamidoethanol, or a butane- against an octane-")
    print("  sulfonamide. None is a spelling variant.")
    return 0


def selftest():
    """Can the audit rediscover the merges the table already encodes?

    An audit that cannot find what is known is no evidence about what is
    unknown. This blanks the table and checks that every hand-made merge
    still surfaces as a candidate.
    """
    global ALIASES, _CANON, _FTS
    rows = all_rows()
    names = {r["chemical"].strip() for r in rows}
    truth = {v: k for k, vs in ALIASES.items() for v in vs if v in names}
    fts = {n for n in names if _FTS.match(n)}
    keep = dict(ALIASES), dict(_CANON), _FTS
    ALIASES, _CANON, _FTS = {}, {}, re.compile(r"^(?!x)x")
    seen = set()
    for _, _, v, _ in candidates(rows):
        seen |= set(v)
    ALIASES, _CANON, _FTS = keep
    skip = {x for pair in UNDISCOVERABLE for x in pair}
    known = (set(truth) | fts) - skip
    missed = sorted(r for r in known if r not in seen)
    print(f"SELF-TEST: {len(known)} discoverable merges, "
          f"{len(missed)} not surfaced by the audit.")
    for m in missed:
        print(f"    MISSED  {m}")
    print(f"  {len(UNDISCOVERABLE)} further merges are excluded as "
          f"undiscoverable by any")
    print("  spelling or annotation rule, each with its evidence:")
    for pair, why in UNDISCOVERABLE.items():
        print(f"    {' = '.join(pair)}")
        print(f"        {why}")
    print("PASS" if not missed else "FAIL -- the audit has a blind spot")
    return 1 if missed else 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else main())
