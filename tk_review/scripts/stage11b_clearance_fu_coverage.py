#!/usr/bin/env python3
"""Where is the gap -- clearance or binding? (Stage 11b)

    python3 -I scripts/stage11b_clearance_fu_coverage.py

Writes db/clearance_fu_coverage.csv. Idempotent.

THE POINT. Before spending effort on the 1497-row protein binding collection
in the hope of extending the tubular endpoint, establish which of the two
terms in T = CL_R / (fu x GFR) is actually missing. This script counts, per
species, the compounds that have a clearance ANYWHERE in the project and asks
how many of them already have an unbound fraction.

THE UNION IS STATED EXPLICITLY rather than remembered, because the answer
depends entirely on what counts as "has a clearance":

    db/vd_clearance.csv        the curated kinetic record
    db/efsa_appendix_c.csv     EFSA's own tabulation (Fujii and others)
    db/tubular_factor_evidence.csv   every clearance wired into stage 6

A compound counts for a species if any of those three carries a non-empty
clearance for it. Species labels are normalised ("mice" -> mouse,
"monkey (cynomolgus)" -> monkey) and chemical aliases folded, or the counts
would be wrong in the flattering direction.
"""
import csv, collections, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
ALIAS = {"PFTrDA": "PFTriDA", "PFTrdA": "PFTriDA", "PFDoA": "PFDoDA",
         "PFUnA": "PFUnDA", "PFTeA": "PFTeDA"}


def species(s):
    s = (s or "").strip().lower()
    if s.startswith("mouse") or s.startswith("mice"):
        return "mouse"
    if s.startswith("monkey") or "macaca" in s or "cynomolgus" in s:
        return "monkey"
    for k in ("human", "rat"):
        if s.startswith(k):
            return k
    return s.split(" (")[0]


def chem(c):
    c = (c or "").strip()
    return ALIAS.get(c, c)


def main():
    clear = collections.defaultdict(lambda: collections.defaultdict(set))
    def add(sp, c, src):
        if sp and c:
            clear[sp][c].add(src)

    for r in csv.DictReader((HERE / "db" / "vd_clearance.csv").open()):
        if r.get("clearance"):
            add(species(r["species"]), chem(r["chemical"]), "vd_clearance")
    for r in csv.DictReader((HERE / "db" / "efsa_appendix_c.csv").open()):
        if r.get("clearance_total_mL_kg_day"):
            add(species(r["species"]), chem(r["chemical"]), "efsa_appendix_c")
    f = HERE / "db" / "tubular_factor_evidence.csv"
    if f.exists():
        for r in csv.DictReader(f.open()):
            add(species(r["species"]), chem(r["chemical"]), "stage6_evidence")

    fu = collections.defaultdict(set)
    for r in csv.DictReader(
            (HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r.get("fu"):
            fu[species(r["species"])].add(chem(r["chemical"]))

    # does the binding collection hold a PLASMA-RELEVANT measurement for a
    # species-compound? albumin or plasma or serum, in that species.
    bind = collections.defaultdict(set)
    for r in csv.DictReader((HERE / "db" / "protein_binding.csv").open()):
        prot = (r.get("protein") or "").lower()
        if any(w in prot for w in ("albumin", "plasma", "serum")):
            bind[species(r["protein_species"])].add(chem(r["chemical"]))

    bar = "=" * 78
    print(bar)
    print("IS THE GAP CLEARANCE OR BINDING?")
    print(bar)
    print(f"  {'species':9}{'w/ clearance':>14}{'have fu':>9}{'gap':>5}"
          f"   gap compounds (plasma-type binding in db?)")
    rows = []
    for sp in sorted(clear, key=lambda s: -len(clear[s])):
        have = set(clear[sp])
        withfu = have & fu.get(sp, set())
        gap = sorted(have - withfu)
        notes = []
        for c in gap:
            notes.append(f"{c}[{'binding YES' if c in bind.get(sp, set()) else 'no binding'}]")
        print(f"  {sp:9}{len(have):>14}{len(withfu):>9}{len(gap):>5}   "
              f"{', '.join(notes) if notes else '-'}")
        for c in sorted(have):
            rows.append([sp, c, "yes" if c in withfu else "no",
                         "yes" if c in bind.get(sp, set()) else "no",
                         ";".join(sorted(clear[sp][c]))])
    print()
    unlockable = [r for r in rows if r[2] == "no" and r[3] == "yes"]
    print(f"  COMPOUNDS THE BINDING COLLECTION COULD UNLOCK: {len(unlockable)}")
    for r in unlockable:
        print(f"    {r[0]} {r[1]}")
    if not unlockable:
        print("    none -- every gap is a compound the collection does not")
        print("    cover in that species either. The bottleneck is CLEARANCE.")
    print()
    p = HERE / "db" / "clearance_fu_coverage.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "chemical", "has_fu",
                    "has_plasma_type_binding", "clearance_sources"])
        w.writerows(rows)
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
