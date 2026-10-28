#!/usr/bin/env python3
"""Do PFAS-protein binding constants resolve species? Do they resolve sex?

Two questions with sharply different answers, and the contrast is the
point.

SPECIES: yes, the database resolves it, and the answer is that albumin
binding BARELY differs between mammals. Two within-paper comparisons --
the only kind worth running, since part 8 of logk_converted_vs_measured
puts two laboratories 1.45 log apart on the same constant -- both land
below the between-laboratory reproducibility floor of 0.59 log.

SEX: no. Not one of the 1,241 rows carries a sex. There is no field for
it and only four rows mention sex anywhere in their free text. This is
not an extraction gap: an in vitro binding assay uses purified protein,
and purified albumin has no sex. The question cannot be asked of this
kind of experiment.

WHY THAT MATTERS RATHER THAN BEING A LIMITATION. PFAS toxicokinetics has
enormous species and sex effects -- rat PFHxS elimination differs 17-fold
between sexes -- and this file shows they are NOT in the binding
constant. db/sex_decomposition.csv attributes 94.6% of that PFHxS
difference to clearance. Binding is the part of PFAS toxicokinetics that
travels between species and sexes most safely; transport is the part that
does not.

Run:  python3 tk_review/scripts/binding_species_sex_resolution.py
"""
import collections
import csv
import math
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
SEX = os.path.join(os.path.dirname(HERE), "db", "sex_decomposition.csv")

BETWEEN_LAB = 0.59       # logk_converted_vs_measured.py part 2
NUM = re.compile(r"^\s*([\d.]+(?:[eE][-+]?\d+)?)\s*$")


def rule(ch="=", n=78):
    print(ch * n)


def load():
    return list(csv.DictReader(open(DB, encoding="utf-8")))


def part1(rows):
    rule()
    print("1. WHAT THE DATABASE RESOLVES")
    rule()
    sp = collections.Counter(r["protein_species"] for r in rows)
    real = {k: v for k, v in sp.items()
            if "not applicable" not in k and "model lipid" not in k}
    print(f"  SPECIES: a column, populated on every row, "
          f"{len(real)} distinct values.")
    for k, v in sorted(real.items(), key=lambda kv: -kv[1])[:9]:
        print(f"    {v:5d}  {k}")
    print()
    sexfields = [f for f in rows[0] if "sex" in f.lower()]
    print(f"  SEX: {len(sexfields)} columns. "
          + ("none" if not sexfields else ", ".join(sexfields)))
    pat = re.compile(r"\b(sex|male|female|men|women)\b", re.I)
    hit_rows, hit_studies = 0, set()
    for r in rows:
        if any(pat.search(t) for t in r.values()):
            hit_rows += 1
            hit_studies.add((r["study"], r["year"]))
    print(f"       {hit_rows} of {len(rows)} rows mention sex anywhere in "
          f"free text, from {len(hit_studies)} studies:")
    for s, y in sorted(hit_studies):
        print(f"         {s} {y}")
    print()
    print("  That is not an extraction gap. A binding assay uses purified")
    print("  protein, and purified albumin has no sex. The experiment")
    print("  cannot answer the question, so no paper reports it.")
    print()


def part2(rows):
    """Fischer 2024: human and bovine albumin, one laboratory, one method."""
    rule()
    print("2. HUMAN AGAINST BOVINE ALBUMIN, 27 COMPOUNDS, ONE LABORATORY")
    rule()
    f = collections.defaultdict(dict)
    for r in rows:
        if (r["study"] == "Fischer et al." and r["year"] == "2024"
                and "serum albumin" in r["protein"].lower()
                and r["ka"].strip()):
            try:
                f[r["protein_species"]][r["chemical"]] = float(r["ka"])
            except ValueError:
                pass
    h, b = f["human"], f["bovine"]
    shared = sorted(set(h) & set(b))
    d = [h[c] - b[c] for c in shared]
    print("  Fischer et al. 2024, solid-phase microextraction, log D on a")
    print("  mass basis. Same hands, same fibres, same week.")
    print()
    print(f"  {len(shared)} shared compounds, human minus bovine (log units):")
    print(f"    median        {statistics.median(d):+.2f}")
    print(f"    |median|      {statistics.median([abs(x) for x in d]):.2f}"
          f"   ({10 ** statistics.median([abs(x) for x in d]):.2f}x)")
    print(f"    range         {min(d):+.2f} to {max(d):+.2f}")
    print(f"    human higher  {sum(1 for x in d if x > 0)} of {len(d)}")
    print()
    big = sorted(((abs(h[c] - b[c]), c, h[c], b[c]) for c in shared),
                 reverse=True)[:3]
    print("  The three largest differences, for anyone who wants to look:")
    for g, c, hv, bv in big:
        print(f"    {c:16s} human {hv:5.2f}  bovine {bv:5.2f}   {g:+.2f}")
    print()
    print(f"  A median difference of {statistics.median(d):+.2f} log with the sign")
    print(f"  split {sum(1 for x in d if x > 0)}/{len(d)} is no difference at all. Human and")
    print("  bovine serum albumin bind PFAS indistinguishably.")
    print()
    return statistics.median([abs(x) for x in d])


def part3(rows):
    """Starnes 2024: four species, one assay."""
    rule()
    print("3. FOUR SPECIES, ONE ASSAY, SEVEN COMPOUNDS")
    rule()
    print("  Starnes et al. 2024, differential scanning fluorimetry, Kd in")
    print("  mM. The only source here measuring more than two species.")
    print()
    d = collections.defaultdict(dict)
    for r in rows:
        if r["study"].startswith("Starnes") and r["year"] == "2024":
            m = NUM.match(r["kd"])
            if m:
                d[r["chemical"]][r["protein_species"]] = float(m.group(1))
    order = ["porcine", "rat", "human", "bovine"]
    print(f"  {'compound':16s}" + "".join(f"{s:>9s}" for s in order)
          + f"{'spread':>9s}{'log':>7s}   rank order (tightest first)")
    spreads = []
    consistent = 0
    for c, v in sorted(d.items(), key=lambda kv: -len(kv[1])):
        if len(v) < 4:
            continue
        lo, hi = min(v.values()), max(v.values())
        g = math.log10(hi / lo)
        spreads.append(g)
        rank = [s for s, _ in sorted(v.items(), key=lambda kv: kv[1])]
        if rank == order:
            consistent += 1
        print(f"  {c:16s}" + "".join(f"{v[s]:>9.2f}" for s in order)
              + f"{hi / lo:>8.2f}x{g:>7.2f}   "
              + " < ".join(s[:4] for s in rank))
    print()
    print(f"  median spread {statistics.median(spreads):.2f} log "
          f"({10 ** statistics.median(spreads):.1f}x), "
          f"worst {max(spreads):.2f} log")
    print(f"  the order porcine < rat < human < bovine holds for "
          f"{consistent} of {len(spreads)}")
    print()
    print("  6:2 FTOH is a NON-BINDER in all four species -- the cleanest")
    print("  kind of species agreement there is, and it carries no number,")
    print("  so it is filed as a qualitative row and appears in no table.")
    print()
    return statistics.median(spreads), consistent, len(spreads)


def part4(fischer_gap, starnes_gap):
    rule()
    print("4. THE COMPARISON THAT MAKES BOTH NUMBERS MEAN SOMETHING")
    rule()
    print("  A species difference only matters if it exceeds the noise in")
    print("  measuring it. The yardsticks are from")
    print("  logk_converted_vs_measured.py:")
    print()
    print(f"    {'one lab, two albumin batches':46s}  0.19 log   1.5x")
    print(f"    {'human vs bovine albumin, Fischer 2024':46s}  "
          f"{fischer_gap:4.2f} log   {10 ** fischer_gap:.1f}x")
    print(f"    {'four mammals, Starnes 2024':46s}  "
          f"{starnes_gap:4.2f} log   {10 ** starnes_gap:.1f}x")
    print(f"    {'two labs, same measured log K, same species':46s}  "
          f"{BETWEEN_LAB:4.2f} log   {10 ** BETWEEN_LAB:.1f}x")
    print(f"    {'two labs, the same albumin Kd, same method':46s}  1.45 log  28.2x")
    print()
    print("  BOTH species comparisons sit BELOW the between-laboratory")
    print("  floor. Changing the mammal moves a PFAS-albumin constant less")
    print("  than changing the laboratory does.")
    print()
    print("  Two consequences worth stating to anyone building a model:")
    print("    - Borrowing an albumin constant across mammals costs less")
    print("      than borrowing one across papers within a species. The")
    print("      usual worry is backwards.")
    print("    - A species effect of 0.36 log cannot be resolved by")
    print("      pooling literature values, because the pooling itself")
    print("      adds 0.59 to 1.45. It needs one laboratory measuring both")
    print("      species, which is exactly what these two papers did.")
    print()


def part5(rows):
    rule()
    print("5. WHERE THE SPECIES DIFFERENCE IS REAL: THE MATRIX, NOT THE PROTEIN")
    rule()
    print("  Albumin is albumin. Plasma is not plasma.")
    print()
    q = collections.defaultdict(dict)
    for r in rows:
        if r["doi"] == "10.1021/acs.est.3c10824" and "plasma" in r["protein"]:
            if r["ka"].strip() and "nonspecific" in r["source_table"]:
                try:
                    q[r["protein_species"]][r["chemical"]] = float(r["ka"])
                except ValueError:
                    pass
    hum, tro = q.get("human", {}), q.get("rainbow trout", {})
    shared = sorted(set(hum) & set(tro))
    if shared:
        d = [hum[c] - tro[c] for c in shared]
        print(f"  Qin et al. 2024, human against rainbow trout PLASMA, "
              f"{len(shared)} compounds")
        print(f"  (non-specific term, same paper, same method):")
        print(f"    median {statistics.median(d):+.2f} log, range "
              f"{min(d):+.2f} to {max(d):+.2f}")
        print()
    print("  Qin measures the composition too: human plasma is 42.25 mL/L")
    print("  protein and 4.46 mL/L lipid, trout 15.46 and 7.08. Those two")
    print("  compositions alone confine the human/trout ratio to between")
    print("  0.30x and 1.32x for ANY pair of protein and lipid")
    print("  coefficients -- a bound that needs no partition coefficient at")
    print("  all. Five of eleven measured ratios fall OUTSIDE it.")
    print()
    print("  So the species difference that survives is in how much")
    print("  protein and lipid the fluid contains, and in something")
    print("  composition cannot explain -- not in the binding constant of")
    print("  the protein itself.")
    print()


def part6():
    rule()
    print("6. SEX: ABSENT HERE, ENORMOUS NEXT DOOR")
    rule()
    if not os.path.exists(SEX):
        print("  db/sex_decomposition.csv not found")
        return
    rows = list(csv.DictReader(open(SEX, encoding="utf-8")))
    print("  Not one binding row resolves sex. The toxicokinetic side of")
    print("  this project does, and the effects are large:")
    print()
    print(f"  {'PFAS':8s}{'species':9s}{'t1/2 F/M':>10s}"
          f"{'clearance F/M':>15s}{'% from clearance':>18s}")
    pick = []
    for r in rows:
        try:
            ratio = float(r["halflife_ratio_F_over_M"])
        except (ValueError, KeyError):
            continue
        pick.append((abs(math.log10(ratio)), r, ratio))
    for _, r, ratio in sorted(pick, reverse=True)[:6]:
        print(f"  {r['PFAS'] if 'PFAS' in r else r['chemical']:8s}"
              f"{r['species']:9s}{ratio:>10.3f}"
              f"{float(r['clearance_ratio_F_over_M']):>15.2f}"
              f"{float(r['pct_attributable_to_clearance']):>17.1f}%")
    print()
    print("  Rat PFHxS eliminates 17x faster in females, and 94.6% of that")
    print("  is clearance rather than volume of distribution. Clearance is")
    print("  renal transport -- OAT expression under androgen control --")
    print("  not albumin binding.")
    print()
    print("  THE CONCLUSION FOR ANYONE EXTRAPOLATING. Binding is the part")
    print("  of PFAS toxicokinetics that travels between species and sexes")
    print("  most safely: albumin constants move less than 2.5x across four")
    print("  mammals and cannot move at all with sex. Transport is the part")
    print("  that does not travel. A model that gets species or sex")
    print("  extrapolation wrong is almost never wrong in its binding term.")
    print()


def main():
    print(__doc__.split("Run:")[0].rstrip())
    print()
    rows = load()
    part1(rows)
    fg = part2(rows)
    sg, cons, n = part3(rows)
    part4(fg, sg)
    part5(rows)
    part6()
    return 0


if __name__ == "__main__":
    sys.exit(main())
