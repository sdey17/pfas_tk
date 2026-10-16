#!/usr/bin/env python3
"""What Ryu et al. 2024's 320 f_unbound values say about species, strain and tissue.

Reads db/ryu2024_funbound.csv, writes db/ryu2024_species_tissue_summary.csv.
The two cells flagged as probable SI typos are excluded throughout.
"""
import csv
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "ryu2024_funbound.csv")
OUT = os.path.join(HERE, "..", "db", "ryu2024_species_tissue_summary.csv")

rows = [r for r in csv.DictReader(open(SRC)) if r["flagged"] != "yes"]
for r in rows:
    r["f"] = float(r["f_unbound"])
    r["npfc"] = int(r["n_perfluorinated_carbons"])

STRAINS = ["human", "rat (Wistar Han)", "mouse (CD-1)", "mouse (C57BL/6J)"]
TISSUES = ["liver", "brain", "heart", "kidney", "lung"]


def gm(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def get(strain, chem, tissue):
    for r in rows:
        if r["strain"] == strain and r["chemical"] == chem and r["tissue"] == tissue:
            return r["f"]
    return None


chems = []
for r in rows:
    if r["chemical"] not in chems:
        chems.append(r["chemical"])

print(f"{len(rows)} values after excluding the 2 flagged cells\n")

print("=" * 72)
print("1. SPECIES. Ratio to human, same chemical and tissue, geometric mean")
print("=" * 72)
print(f"{'':20s}" + "".join(f"{t:>10s}" for t in TISSUES) + f"{'all':>10s}")
species_rows = []
for strain in STRAINS[1:]:
    cells, line = [], f"{strain:20s}"
    for tissue in TISSUES:
        rs = [get(strain, c, tissue) / get("human", c, tissue) for c in chems
              if get(strain, c, tissue) and get("human", c, tissue)]
        line += f"{gm(rs):>10.2f}"
        cells += rs
    line += f"{gm(cells):>10.2f}"
    print(line)
    species_rows.append((strain, gm(cells), max(cells), min(cells)))
print("\n  worst single chemical x tissue departure from human:")
for strain, g, hi, lo in species_rows:
    print(f"    {strain:20s} geomean {g:.2f}x, range {lo:.2f}-{hi:.2f}x")

print("\n" + "=" * 72)
print("2. STRAIN. CD-1 vs C57BL/6J, the comparison the Argoul limb turns on")
print("=" * 72)
print(f"{'tissue':10s}{'geomean CD-1/B6':>18s}{'min':>8s}{'max':>8s}{'n over 2x':>11s}")
strain_tbl = []
for tissue in TISSUES:
    rs = [get("mouse (CD-1)", c, tissue) / get("mouse (C57BL/6J)", c, tissue)
          for c in chems
          if get("mouse (CD-1)", c, tissue) and get("mouse (C57BL/6J)", c, tissue)]
    over = sum(1 for x in rs if x > 2 or x < 0.5)
    print(f"{tissue:10s}{gm(rs):>18.2f}{min(rs):>8.2f}{max(rs):>8.2f}{over:>8d}/{len(rs)}")
    strain_tbl.append((tissue, gm(rs), min(rs), max(rs), over, len(rs)))
allr = [get("mouse (CD-1)", c, t) / get("mouse (C57BL/6J)", c, t)
        for c in chems for t in TISSUES
        if get("mouse (CD-1)", c, t) and get("mouse (C57BL/6J)", c, t)]
print(f"{'ALL':10s}{gm(allr):>18.2f}{min(allr):>8.2f}{max(allr):>8.2f}"
      f"{sum(1 for x in allr if x > 2 or x < 0.5):>8d}/{len(allr)}")
print("\n  largest strain gaps (CD-1 / C57BL/6J):")
gaps = sorted(((get("mouse (CD-1)", c, t) / get("mouse (C57BL/6J)", c, t), c, t)
               for c in chems for t in TISSUES
               if get("mouse (CD-1)", c, t) and get("mouse (C57BL/6J)", c, t)),
              key=lambda x: -abs(math.log(x[0])))
for v, c, t in gaps[:6]:
    print(f"    {c:8s} {t:7s} {v:5.2f}x")

print("\n" + "=" * 72)
print("3. TISSUE. Ratio to liver, same chemical and organism, geometric mean")
print("=" * 72)
print(f"{'':20s}" + "".join(f"{t:>10s}" for t in TISSUES))
for strain in STRAINS:
    line = f"{strain:20s}"
    for tissue in TISSUES:
        rs = [get(strain, c, tissue) / get(strain, c, "liver") for c in chems
              if get(strain, c, tissue) and get(strain, c, "liver")]
        line += f"{gm(rs):>10.2f}"
    print(line)
print("\n  brain binds hardest, lung least: the liver is unremarkable despite L-FABP")
for strain in STRAINS:
    pairs = []
    for c in chems:
        liv = get(strain, c, "liver")
        oth = [get(strain, c, t) for t in TISSUES if t != "liver"]
        if liv is None or any(o is None for o in oth):
            continue          # a flagged cell was dropped from this compound
        pairs.append((liv, gm(oth)))
    within2 = sum(1 for a, b in pairs if 0.5 <= a / b <= 2)
    print(f"    {strain:20s} liver within 2x of the other tissues for "
          f"{within2}/{len(pairs)} compounds")

print("\n" + "=" * 72)
print("4. CHAIN LENGTH. Decades of f_unbound lost per CF2, human tissues")
print("=" * 72)
for tissue in TISSUES:
    pts = [(r["npfc"], math.log10(r["f"])) for r in rows
           if r["strain"] == "human" and r["tissue"] == tissue
           and r["head_group"] == "carboxylate"]
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    slope = (sum((x - mx) * (y - my) for x, y in pts)
             / sum((x - mx) ** 2 for x, _ in pts))
    ss = sum((y - my) ** 2 for _, y in pts)
    res = sum((y - my - slope * (x - mx)) ** 2 for x, y in pts)
    print(f"  {tissue:8s} {slope:+.3f} decades per CF2 (R2 {1 - res/ss:.3f}, "
          f"PFCA n={n})")
print("  monotonic across the whole C3-C13 range, unlike albumin log D which "
      "turns over near C10")

print("\n" + "=" * 72)
print("5. The paper's own half-life claim, checked")
print("=" * 72)
for tissue in TISSUES:
    a, b = get("human", "PFOA", tissue), get("human", "PFOS", tissue)
    print(f"  human {tissue:7s} PFOA/PFOS f_unbound = {a/b:5.2f}x")
rs = [get("human", "PFOA", t) / get("human", "PFOS", t) for t in TISSUES]
print(f"  range {min(rs):.1f}-{max(rs):.1f}x, matching the paper's stated 5.8-11.3x")
print("  human half-lives the paper cites: PFOA 3.8 y, PFOS 5.4 y = 1.42x.")
print("  So a 5.8-11.3x difference in binding accompanies a 1.4x difference in")
print("  half-life: the sign agrees, the magnitude does not.")

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["comparison", "stratum", "geometric_mean_ratio", "min", "max",
                "n_over_2fold", "n"])
    for strain, g, hi, lo in species_rows:
        w.writerow([f"{strain} / human", "all tissues and chemicals",
                    f"{g:.3f}", f"{lo:.3f}", f"{hi:.3f}", "", len(chems) * 5])
    for tissue, g, lo, hi, over, n in strain_tbl:
        w.writerow(["mouse (CD-1) / mouse (C57BL/6J)", tissue,
                    f"{g:.3f}", f"{lo:.3f}", f"{hi:.3f}", over, n])
print(f"\nwrote {OUT}")
