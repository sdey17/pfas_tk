#!/usr/bin/env python3
"""Is rat albumin like human albumin, the way bovine albumin is?

THE TRAP THIS AVOIDS.  The obvious way to answer is to pool every rat
value and every human value and compare the two piles.  That answer
is worthless here, because the piles come from different
laboratories, and the laboratory term in this literature is far
larger than the species term -- Fischer 2024's own within-lab
human/bovine pairs differ by a median 0.00 log while the same
compound across labs spans whole log units.  A pooled comparison
would be measuring who ran the assay, not which species.

So only WITHIN-LABORATORY pairs count: the same bench, the same
method, the same compound, two species.  There are exactly two such
sources for rat albumin in this database, and one of them is the
only study that ran four species side by side.

SIGN CONVENTION.  Everything is reported as a difference in log K,
binding strength, positive meaning the other species binds TIGHTER
than human.  Starnes reports Kd, where smaller is tighter, so the
conversion flips the sign: dlogK = log10(Kd_human / Kd_species).

WHAT STARNES IS AND IS NOT.  Starnes 2024 is the four-species study,
by differential scanning fluorimetry, and it is HELD OUT of this
project's plotted figures: a DSF midpoint in the millimolar range is
not a physiological binding constant, and it fails the bound in
scripts/binding_kd_physiological_bound.py.  That holdout is about the
ABSOLUTE values.  It does not disqualify the study from answering a
RELATIVE question, because every species here was run on one
instrument, in one buffer, against one protein preparation, so the
systematic offset that makes the absolute number unusable cancels in
the ratio.  This is the single legitimate use of a held-out study in
this project, and it is confined to this script.

Run:  python3 tk_review/scripts/species_albumin_comparison.py
Reads only; writes nothing.
"""
import collections
import csv
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
LIT = os.path.join(HERE, "..", "figures",
                   "albumin_chainlength_data.csv")


def num(s):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return None


def starnes():
    """Four species, one bench. Returns {species: {compound: Kd_mM}}."""
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if not r["study"].startswith("Starnes"):
            continue
        if r["kd_units"].strip() != "mM":      # the four-species set
            continue
        v = num(r["kd"])
        if v is None:                          # "nonbinding"
            continue
        out[r["protein_species"]][r["chemical"].split(" (")[0]] = v
    return out


def within_lab_pairs():
    """Every source in the plotted CSV reporting two species."""
    by = collections.defaultdict(lambda: collections.defaultdict(dict))
    for r in csv.DictReader(open(LIT, encoding="utf-8")):
        by[r["source"]][r["species"]][r["pfas"]] = float(r["log_k"])
    return by


def main():
    print("DOES RAT ALBUMIN BEHAVE LIKE HUMAN ALBUMIN?")
    print("Within-laboratory pairs only. Positive = binds TIGHTER "
          "than human.\n")

    # ---- the four-species study -----------------------------------
    S = starnes()
    human = S.get("human", {})
    print("1. STARNES 2024 -- four species, one laboratory, DSF")
    print("   (held out of the figures for its ABSOLUTE values; used "
          "here only\n    for ratios, which cancel the offset)\n")
    print(f"   {'compound':12}" + "".join(f"{s:>11}"
          for s in ("human", "bovine", "porcine", "rat"))
          + "    (Kd, mM -- lower binds tighter)")
    for c in sorted(human):
        row = "".join(f"{S[s].get(c, float('nan')):>11.2f}"
                      for s in ("human", "bovine", "porcine", "rat"))
        print(f"   {c:12}{row}")

    print(f"\n   {'vs human':12}{'n':>4}{'median dlogK':>15}"
          f"{'sd':>8}{'fold':>9}")
    summary = {}
    for sp in ("bovine", "porcine", "rat"):
        common = sorted(set(human) & set(S.get(sp, {})))
        d = [math.log10(human[c] / S[sp][c]) for c in common]
        if not d:
            continue
        m, s = st.median(d), st.pstdev(d)
        summary[sp] = (len(d), m, s)
        print(f"   {sp:12}{len(d):>4}{m:>+15.3f}{s:>8.3f}"
              f"{10**abs(m):>8.2f}x")

    if {"bovine", "rat"} <= set(summary):
        nb, mb, sb = summary["bovine"]
        nr, mr, sr = summary["rat"]
        print(f"\n   Rat sits {10**abs(mr):.2f}x from human; bovine "
              f"sits {10**abs(mb):.2f}x.")
        print(f"   They are the same size of difference, and rat is "
              f"the more CONSISTENT")
        print(f"   of the two across compounds (sd {sr:.3f} vs "
              f"{sb:.3f}).")

    # ---- every other within-lab pair in the plotted data ----------
    print("\n2. EVERY OTHER WITHIN-LAB CROSS-SPECIES PAIR IN THE "
          "PLOTTED DATA\n")
    for src, sp in sorted(within_lab_pairs().items()):
        if len(sp) < 2:
            continue
        for a in sorted(sp):
            for b in sorted(sp):
                if a >= b:
                    continue
                common = sorted(set(sp[a]) & set(sp[b]))
                if not common:
                    continue
                d = [sp[b][c] - sp[a][c] for c in common]
                tail = (f"  sd {st.pstdev(d):.3f}" if len(d) > 1
                        else "")
                print(f"   {src:14} {b:7} - {a:7}  n={len(d):2d}  "
                      f"median {st.median(d):+.3f} log{tail}")
    # Han's zero is not two measurements that agreed. Worth saying,
    # because +0.000 log otherwise reads as a precise result.
    print("\n   Han 2003's exact zero is NOT an independent "
          "agreement: the paper")
    print("   reports the SAME range, 0.3-0.4 mM, for rat and for "
          "human PFOA, so")
    print("   it is one shared number. Read it as 'no difference "
          "resolved at this")
    print("   resolution', not as a measured equality.")

    # ---- the disagreement inside the one big rat study ------------
    print("\n3. THE CAVEAT: FAN 2025, the largest rat albumin set")
    fan = collections.defaultdict(dict)
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if not (r["study"].startswith("Fan")
                and "albumin" in r["protein"].lower()):
            continue
        v = num(r["kd"])
        if v is None or r["kd_units"].strip() != "uM":
            continue
        meth = ("ITC" if "calorimetry" in r["method"].lower()
                else "ultrafiltration"
                if "ultrafilt" in r["method"].lower() else None)
        if meth:
            fan[meth][r["chemical"].split(" (")[0]] = v
    both = sorted(set(fan.get("ITC", {})) & set(fan.get(
        "ultrafiltration", {})))
    if both:
        d = [math.log10(fan["ITC"][c] / fan["ultrafiltration"][c])
             for c in both]
        print(f"   Fan measured rat albumin TWICE, by two methods, on "
              f"{len(both)} shared compounds.")
        print(f"   Median log10(Kd_ITC / Kd_ultrafiltration) = "
              f"{st.median(d):+.2f}  -> the two")
        print(f"   disagree by about {10**abs(st.median(d)):.0f}x "
              f"WITHIN ONE PAPER.")
        print("   So the rat number you get depends more on the "
              "method than on the")
        print("   species. That is the honest limit on any rat "
              "albumin value here.")

    fold = (f"~{10**abs(st.median(d)):.0f}x" if both else "several-fold")
    print("\nBOTTOM LINE")
    print("  Rat albumin is as close to human as bovine is -- same "
          "order of")
    print("  difference, slightly more consistent across compounds. "
          "What it is")
    print("  NOT is well measured: 35 rows, 3 studies, and the "
          "largest of them")
    print(f"  disagrees with itself by {fold} depending on method.")
    print("  Bovine's comparability rests on 25 within-lab pairs "
          "(Fischer 2024);")
    print("  rat's rests on 7 (Starnes, one held-out study) plus a "
          "Han 2003 pair")
    print("  that is really one shared number. Treat the equivalence "
          "as")
    print("  provisional, and never substitute a rat value for a "
          "human one")
    print("  without saying which study it came from.")


if __name__ == "__main__":
    main()
