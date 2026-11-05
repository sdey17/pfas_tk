#!/usr/bin/env python3
"""Where a paper reports both a constant and its thermodynamics, check them.

A binding constant and a free energy are the same measurement twice:

    dG = -RT ln K                                              (1)

and where an enthalpy and entropy are reported too, they close the
triangle:

    dG = dH - T dS                                             (2)

Neither is a model. Both are definitions, so a paper that prints all
of these has given a free consistency test on its own table, and a
failure localises to one printed cell.

WHAT IT FOUND.  Sheng 2016's Table 3 and SI Table S3 are internally
inconsistent on 12 of their 13 checkable rows.  An ITC free energy
MUST equal -RT ln Ka, because Ka comes from the same fit, yet the
printed dG for PFNA wild type is -5.275 kJ/mol where its own
Ka of 3.18e5 requires -31.41.  Only PFOA wild type is consistent on
all three routes.  The transcription in this database was verified
character by character against the published table before this was
reported.

It does NOT touch the binding constants.  Ka and Kd are reciprocals
and agree to better than 0.2 % on all 15 rows, so the log K values
this project builds from Sheng's Kd are unaffected.  What fails is
the thermodynamic block, and two distinct patterns separate within
it -- see the verdicts.

THIS WAS WRITTEN AFTER TWO HITS FOUND BY HAND.  Peng 2024's PFOA at
304 K is the only point in its Table 2 whose printed dG disagrees
with -RT ln Kb -- and it is also the only compound whose series turns
over, which the disagreement explains (see plot_peng2024_vanthoff.py).
Gao's SI Table S9 misattributes one of Zhang's constants.  Two single
bad cells, both found by reading. This does the reading.

SCOPE, STATED PLAINLY.  The check needs a constant, a temperature and
a free energy IN THE SAME ROW. This database transcribed
thermodynamics for only a few studies, so the check covers far fewer
rows than carry a constant, and the gap is reported rather than
hidden: UNTRANSCRIBED below lists studies known to report
thermodynamics whose numbers were never extracted, which is a
to-do, not a pass.

A docking score in kcal/mol is NOT a free energy of binding for this
purpose and is excluded: it is a scoring function's output, not a
measured dG.

Run:  python3 tk_review/scripts/thermo_consistency.py
Writes db/thermo_consistency.csv
"""
import csv
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "db", "thermo_consistency.csv")

R = 8.314                       # J/mol/K
TOL = 1.0                       # kJ/mol; below this, call it rounding

NUM = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
RE_DG = re.compile(r"dG\s*" + NUM + r"\s*(?:\+/-\s*[\d.]+\s*)?kJ", re.I)
RE_DH = re.compile(r"dH\s*" + NUM + r"\s*(?:\+/-\s*[\d.]+\s*)?kJ", re.I)
RE_DS = re.compile(r"dS\s*" + NUM + r"\s*(?:\+/-\s*[\d.]+\s*)?J", re.I)
RE_KA = re.compile(r"Ka\s*\(?" + NUM + r"(?:\s*\+/-\s*[\d.]+)?\)?e"
                   r"([-+]?\d+)", re.I)
RE_KD = re.compile(r"Kd\s*" + NUM + r"\s*uM", re.I)
RE_T = re.compile(NUM)

# Studies that report thermodynamics whose values are not in this
# database. Listed so the scope of the check is honest.
UNTRANSCRIBED = [
    ("Jia et al.", "2022", "Table 1 gives dS, dH and dG from a van't "
     "Hoff fit over 21-40 C; only the 25 C Kd was transcribed"),
    ("Fedorenko et al.", "2021", "reports binding as dominantly "
     "entropic; no numbers transcribed"),
]

FIELDS = ["study", "year", "chemical", "protein", "variant", "T_K",
          "K_source", "K_M^-1", "dG_reported", "dG_from_K",
          "gap_from_K", "dH", "dS", "dG_from_dH_dS",
          "gap_from_dH_dS", "verdict"]


def temp_K(text):
    m = re.search(r"(\d+)\s*\(?(\d{3})\s*K", text or "")
    if m:
        return float(m.group(2))
    m = RE_T.search(text or "")
    return float(m.group(1)) + 273.15 if m else None


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    out = []
    no_dg = no_K = no_T = 0
    for r in rows:
        st = r["source_table"]
        if "kcal/mol" in (r["kd_units"] + r["ka_units"] + r["method"]):
            continue                       # docking score, not a dG
        mg = RE_DG.search(st)
        if not mg:
            if re.search(r"\bdG\b|ΔG", st):
                no_dg += 1
            continue
        dG = float(mg.group(1))
        T = temp_K(r["temperature_c"])
        if T is None:
            no_T += 1
            continue
        mk, md = RE_KA.search(st), RE_KD.search(st)
        if mk:
            K, ksrc = float(mk.group(1)) * 10 ** float(mk.group(2)), "Ka"
        elif md:
            K, ksrc = 1e6 / float(md.group(1)), "1/Kd"
        else:
            no_K += 1
            continue

        dG_K = -R * T * math.log(K) / 1000.0
        gap_K = dG - dG_K
        mh, ms = RE_DH.search(st), RE_DS.search(st)
        dH = float(mh.group(1)) if mh else None
        dS = float(ms.group(1)) if ms else None
        dG_HS = (dH - T * dS / 1000.0) if (dH is not None
                                           and dS is not None) else None
        gap_HS = None if dG_HS is None else dG - dG_HS

        bad_K = abs(gap_K) > TOL
        bad_HS = gap_HS is not None and abs(gap_HS) > TOL
        # Two failures that mean different things, so they are named
        # separately rather than lumped as "inconsistent".
        #
        #   dG is the odd one out -- it disagrees with BOTH the
        #   constant and dH - TdS, which agree with each other. The
        #   fault localises to the printed dG cell.
        #
        #   the dG/dH/dS block is self-consistent but disagrees with
        #   the constant. Nothing localises: either the constant or
        #   the whole thermodynamic block is from a different fit.
        agree_KH = (gap_HS is not None
                    and abs(dG_K - dG_HS) <= 2 * TOL)
        verdict = ("consistent" if not bad_K and not bad_HS else
                   "dG is the odd one out; K and dH-TdS agree"
                   if bad_K and bad_HS and agree_KH else
                   "dG and dH-TdS agree; the CONSTANT disagrees"
                   if bad_K and not bad_HS else
                   "dG disagrees with K only" if bad_K else
                   "dG disagrees with dH-TdS only")
        variant = ""
        mv = re.search(r"Table S?\d+,\s*([A-Za-z0-9]+)", st)
        if mv:
            variant = mv.group(1)
        out.append(dict(
            study=r["study"], year=r["year"], chemical=r["chemical"],
            protein=r["protein"][:40], variant=variant, T_K=f"{T:.2f}",
            K_source=ksrc, **{"K_M^-1": f"{K:.4g}"},
            dG_reported=f"{dG:.2f}", dG_from_K=f"{dG_K:.2f}",
            gap_from_K=f"{gap_K:+.2f}",
            dH="" if dH is None else f"{dH:.2f}",
            dS="" if dS is None else f"{dS:.1f}",
            dG_from_dH_dS="" if dG_HS is None else f"{dG_HS:.2f}",
            gap_from_dH_dS="" if gap_HS is None else f"{gap_HS:+.2f}",
            verdict=verdict))

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)

    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print(f"  {len(rows)} rows in the database")
    print(f"  {len(out)} carry a constant, a temperature AND a free "
          f"energy -- the only ones checkable")
    if no_dg or no_K or no_T:
        print(f"  skipped: {no_dg} mention dG without a parseable "
              f"value, {no_K} have a dG but no constant, {no_T} no "
              f"temperature")
    print()
    bad = [r for r in out if r["verdict"] != "consistent"]
    print(f"  {len(out) - len(bad)} consistent, {len(bad)} NOT, at a "
          f"{TOL:.0f} kJ/mol tolerance\n")
    if bad:
        import collections
        for v, n in collections.Counter(
                r["verdict"] for r in bad).most_common():
            print(f"    {n:3d}  {v}")
        print()
        print(f"  {'study':17}{'cmpd':8}{'var':7}{'dG':>8}"
              f"{'-RTlnK':>9}{'dH-TdS':>9}   pattern")
        for r in sorted(bad, key=lambda x: (x["verdict"],
                                            -abs(float(x["gap_from_K"])))):
            short = ("dG odd" if "odd one out" in r["verdict"] else
                     "K odd" if "CONSTANT" in r["verdict"] else "mixed")
            print(f"  {r['study'][:12]+' '+r['year']:17}"
                  f"{r['chemical'][:7]:8}{r['variant'][:6]:7}"
                  f"{r['dG_reported']:>8}{r['dG_from_K']:>9}"
                  f"{r['dG_from_dH_dS'] or '-':>9}   {short}")
    print()
    # The constants themselves are checked separately, because the
    # question a reader will ask next is whether an inconsistent
    # thermodynamic column means the Kd is unsafe too.
    print("  DOES IT TOUCH THE CONSTANTS? Ka and Kd are reciprocals,")
    print("  so a table that prints both checks itself:")
    recip = 0
    tot = 0
    for r in rows:
        st = r["source_table"]
        mk, md = RE_KA.search(st), RE_KD.search(st)
        if not (mk and md):
            continue
        tot += 1
        ka = float(mk.group(1)) * 10 ** float(mk.group(2))
        if abs((1e6 / ka) / float(md.group(1)) - 1.0) > 0.02:
            recip += 1
            print(f"    MISMATCH {r['study']} {r['year']} "
                  f"{r['chemical']}")
    print(f"    {tot} rows print both; {recip} disagree by more "
          f"than 2 %.")
    if not recip:
        print("    So the binding constants are internally sound even")
        print("    where the thermodynamic columns are not, and the")
        print("    log K values this project builds from Kd are")
        print("    unaffected.")
    print()
    print("  UNTRANSCRIBED -- reports thermodynamics, values not in "
          "this database:")
    for s, y, what in UNTRANSCRIBED:
        print(f"    {s} {y}: {what}")


if __name__ == "__main__":
    main()
