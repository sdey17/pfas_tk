#!/usr/bin/env python3
"""What is actually in hand for a mass balance, and where it came from.

Three quantities matter for fitting a serum mass-balance model, and they are
not interchangeable:

  Kd            a dissociation constant for one binding site
  K_prot/water  a protein-water partition coefficient -- a whole phase's
                sorptive capacity, which is what the model needs
  K_mem/water   the same for the phospholipid membrane, the lipid limb

plus f_unbound, which is not an input at all but the quantity a fitted model
is tested against.

This counts them, measures the overlap compound by compound, and lists the
references each came from.

Run:  python3 scripts/binding_data_inventory.py
"""
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_binding_affinity_table import build, conditions  # noqa: E402

SRC = os.path.join(HERE, "..", "db", "protein_binding.csv")

KIND = {
    "Kd": "Kd",
    "Ka": "Kd",                      # same quantity, reciprocal
    "log D": "K_prot/water",
    "log K_membrane": "K_mem/water",
    "f_unbound": "f_unbound",
}


def banner(t):
    print()
    print(t)
    print("-" * len(t))


def main():
    rows, out, skipped, qual, dropped = build()

    # --- 1. counts -------------------------------------------------------
    banner("1. How many of each quantity")
    kinds = collections.Counter()
    usable = collections.Counter()
    for r in out:
        k = KIND.get(r["measure"])
        if not k:
            continue
        kinds[k] += 1
        if k == "Kd":
            usable["Kd"] += 1 if r["kd_uM"] else 0
        elif k == "K_prot/water":
            usable[k] += 1 if r["logK"] else 0
        elif k == "K_mem/water":
            # the membrane coefficient lives in `value`; logK is reserved for
            # the protein limb so the two cannot be summed by accident
            usable[k] += 1 if r["value"] else 0
        else:
            usable[k] += 1
    print(f"{'quantity':16s} {'rows':>6s} {'with a usable number':>21s}")
    for k in ("Kd", "K_prot/water", "K_mem/water", "f_unbound"):
        print(f"{k:16s} {kinds[k]:6d} {usable[k]:21d}")
    print(f"{'(total)':16s} {sum(kinds.values()):6d}")

    # --- 2. overlap ------------------------------------------------------
    banner("2. Overlap, by compound")
    have = collections.defaultdict(set)
    for r in out:
        k = KIND.get(r["measure"])
        if not k:
            continue
        if k == "Kd" and not r["kd_uM"]:
            continue
        if k == "K_prot/water" and not r["logK"]:
            continue
        have[r["pfas"]].add(k)
    combos = collections.Counter(frozenset(v) for v in have.values())
    print(f"{len(have)} compounds carry at least one of the four.\n")
    want = ("Kd", "K_prot/water", "K_mem/water")
    all3 = sorted(c for c, v in have.items() if set(want) <= v)
    print(f"  all three (Kd + both partition coefficients): {len(all3)}")
    print("    " + ", ".join(all3))
    for a, b in (("Kd", "K_prot/water"), ("Kd", "K_mem/water"),
                 ("K_prot/water", "K_mem/water")):
        n = sorted(c for c, v in have.items() if {a, b} <= v)
        print(f"\n  {a} + {b}: {len(n)}")
        print("    " + ", ".join(n))
    print()
    print("  the fittable set -- both limbs AND a measured unbound fraction "
          "to test against:")
    fit = sorted(c for c, v in have.items()
                 if {"K_prot/water", "K_mem/water", "f_unbound"} <= v)
    print(f"    {len(fit)}: " + ", ".join(fit))
    print()
    print(f"{'carries':44s} {'compounds':>9s}")
    for combo, n in combos.most_common():
        print(f"{' + '.join(sorted(combo)) or '(none)':44s} {n:9d}")

    # --- 3. references ---------------------------------------------------
    banner("3. Where each quantity came from")
    with open(SRC) as fh:
        raw = list(csv.DictReader(fh))
    ident = {}
    for r in raw:
        key = (r["study"], r["year"])
        if key not in ident and (r["doi"] or r["pmid"]):
            ident[key] = r["doi"] or ("PMID " + r["pmid"])
    bysrc = collections.defaultdict(collections.Counter)
    for r in out:
        k = KIND.get(r["measure"])
        if k:
            bysrc[r["source"]][k] += 1
    for quantity in ("K_prot/water", "K_mem/water", "Kd", "f_unbound"):
        print(f"\n  {quantity}")
        ranked = sorted(((c[quantity], s) for s, c in bysrc.items()
                         if c[quantity]), reverse=True)
        for n, s in ranked:
            surname, _, year = s.rpartition(" ")
            doi = next((v for (st, yr), v in ident.items()
                        if st.startswith(surname) and yr == year), "")
            print(f"    {n:4d}  {s:18s} {doi}")


if __name__ == "__main__":
    main()
