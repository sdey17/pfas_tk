#!/usr/bin/env python3
"""Write db/pfas_name_map.csv: every raw name, its canonical form, evidence.

The map is the auditable record of what scripts/pfas_names.py decided.
It is written, not hand-maintained, so it cannot drift from the module.

Run:  python3 tk_review/scripts/build_pfas_name_map.py
"""
import collections
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")
OUT = os.path.join(os.path.dirname(HERE), "db", "pfas_name_map.csv")

_s = importlib.util.spec_from_file_location(
    "pfas_names", os.path.join(HERE, "pfas_names.py"))
N = importlib.util.module_from_spec(_s)
_s.loader.exec_module(N)

FIELDS = ["chemical_raw", "chemical_canonical", "merged", "rows",
          "n_perfluorinated_carbons", "head_group", "sources", "note"]


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    agg = collections.defaultdict(
        lambda: {"rows": 0, "ncf": set(), "head": set(), "src": set()})
    for r in rows:
        n = r["chemical"].strip()
        a = agg[n]
        a["rows"] += 1
        if r["n_perfluorinated_carbons"].strip():
            a["ncf"].add(r["n_perfluorinated_carbons"].strip())
        if r["head_group"].strip():
            a["head"].add(r["head_group"].strip())
        a["src"].add(f"{r['study']} {r['year']}")

    unres = {}
    for pair, why in N.UNRESOLVED.items():
        for p in pair:
            unres[p] = f"NOT merged with {', '.join(x for x in pair if x != p)}: {why}"
    iso = {}
    for base, vs in N.ISOMER_VARIANTS.items():
        for v in vs:
            iso[v] = (f"NOT merged with {base}: an isomer fraction is not "
                      f"the technical mixture")

    out = []
    for n in sorted(agg):
        c = N.canon(n)
        a = agg[n]
        note = unres.get(n, "") or iso.get(n, "")
        if not N.is_single_compound(n):
            note = note or "prose list or range, not one compound"
        out.append({
            "chemical_raw": n, "chemical_canonical": c,
            "merged": "yes" if c != n else "",
            "rows": a["rows"],
            "n_perfluorinated_carbons": "; ".join(sorted(a["ncf"])),
            "head_group": "; ".join(sorted(a["head"]))[:80],
            "sources": "; ".join(sorted(a["src"]))[:110],
            "note": note})
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)

    merged = [r for r in out if r["merged"]]
    print(f"wrote {os.path.abspath(OUT)}  ({len(out)} raw names)")
    print(f"  {len(merged)} merged into another spelling")
    print(f"  {len({r['chemical_canonical'] for r in out if N.is_single_compound(r['chemical_raw'])})}"
          f" distinct compounds after reconciliation")
    print(f"  {sum(1 for r in out if r['note'])} carry a note")
    print()
    print("  merged names and where their rows go:")
    for r in sorted(merged, key=lambda r: r["chemical_canonical"]):
        print(f"    {r['chemical_raw']:32s} -> {r['chemical_canonical']:16s}"
              f"{r['rows']:>4d} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
