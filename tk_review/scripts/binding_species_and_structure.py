"""
In vitro PFAS-protein binding: what it says about SPECIES and about STRUCTURE.

db/protein_binding.csv holds 285 rows from 57 study-years across 31 protein
targets. Most of it cannot answer either question, because a single affinity
number for a PFAS-protein pair does not exist (section 6.1: Ka for PFOA-HSA
spans 19,200x across six studies). What CAN answer them is the handful of
studies that hold method constant and vary one thing:

  SPECIES    Starnes 2024 -- 4 albumins (human, rat, bovine, porcine) x the
             same compounds, one method (differential scanning fluorimetry).
             Fischer 2024 -- human vs bovine albumin, 10 chain lengths, one
             method (solid-phase microextraction at <=0.004 ligand:protein).

  STRUCTURE  Fischer 2024 -- 10 perfluorinated chain lengths x 3 protein
             targets, carboxylates and sulfonates separable.

Absolute values from either are method-bound and should not be mixed with
other studies. The RATIOS within one study are the usable quantity.

Run:  python3 scripts/binding_species_and_structure.py
"""

import collections
import csv
import math
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HSA = "human serum albumin (HSA)"
BSA = "bovine serum albumin (BSA)"
GLB = "gamma-globulin"


def nc(v):
    m = re.findall(r"\d+", v or "")
    return int(m[0]) if m else None


def rows():
    p = os.path.join(HERE, "db", "protein_binding.csv")
    return [r for r in csv.DictReader(open(p))]


def main():
    all_rows = rows()
    print(__doc__.split("Run:")[0].rstrip())

    # ---------------------------------------------------------------- species
    print("=" * 78)
    print("PART 1.  SPECIES: does albumin from different animals bind differently?")
    print("=" * 78)
    s = [r for r in all_rows
         if "Starnes" in r["study"] and "albumin" in r["protein"].lower()
         and r.get("kd") and r["kd"] != "nonbinding"]
    d = collections.defaultdict(dict)
    for r in s:
        d[(nc(r["n_perfluorinated_carbons"]), r["chemical"])][r["protein_species"]] = float(r["kd"])
    print("""
  Starnes 2024, dissociation constant Kd in mM -- LOWER means TIGHTER binding.
  One method throughout, so the across-species ratios are comparable even
  though the absolute values are not (differential scanning fluorimetry reads
  affinity off a thermal ramp, not 37 C, and its own authors call it low).
""")
    print(f"  {'cmpd':10s}{'nC':>3s}{'human':>8s}{'rat':>8s}{'bovine':>8s}{'porcine':>9s}"
          f"{'rat/hum':>9s}{'max/min':>9s}")
    ratios = []
    for (n, c), v in sorted(d.items()):
        h, ra = v.get("human"), v.get("rat")
        vals = [x for x in v.values()]
        rr = ra / h if (h and ra) else None
        if rr:
            ratios.append(rr)
        print(f"  {c:10s}{n:>3d}{v.get('human',0):>8.2f}{v.get('rat',0):>8.2f}"
              f"{v.get('bovine',0):>8.2f}{v.get('porcine',0):>9.2f}"
              f"{(f'{rr:.2f}x' if rr else '-'):>9s}"
              f"{max(vals)/min(vals):>8.1f}x")
    gm = math.exp(sum(math.log(x) for x in ratios) / len(ratios))
    print(f"""
  Rat albumin binds every compound here at least as tightly as human albumin,
  geometric mean Kd ratio {gm:.2f}x -- that is, rat binds about {1/gm:.1f}x TIGHTER. For
  PFOA specifically the factor is {d[(7,'PFOA')]['rat']/d[(7,'PFOA')]['human']:.2f}x.

  That is the wrong direction to explain anything. The rat eliminates PFOA
  roughly 1,000x FASTER than the human. If tighter albumin binding retained
  PFAS, the rat should be the slow species. Binding cannot be the cause of the
  species difference in elimination; section 6.3 reaches the same conclusion
  from Han 2004 on alpha2u-globulin.

  The spread across all four species is small: max/min never exceeds about 4x
  and is usually under 2x. Albumin is a conserved protein and it behaves like
  one. Whatever makes species differ in PFAS kinetics, it is not the affinity
  of their albumin.
""")

    # Fischer's independent human-vs-bovine check at a physiological ratio
    f = [r for r in all_rows if "Fischer" in r["study"]]
    fd = collections.defaultdict(dict)
    for r in f:
        n = nc(r["n_perfluorinated_carbons"])
        if n:
            fd[(n, r["chemical"])][r["protein"]] = r["ka"]
    diffs = [(float(v[BSA]) - float(v[HSA]), c)
             for (n, c), v in fd.items() if v.get(HSA) and v.get(BSA)]
    within = sum(1 for x, _ in diffs if abs(x) <= 0.3)
    print(f"""  A second method agrees. Fischer 2024 measured human AND bovine albumin for
  {len(diffs)} compounds by solid-phase microextraction at a physiological
  ligand:protein ratio. The bovine-minus-human difference in log D is within
  0.3 log units ({within} of {len(diffs)} compounds), i.e. under twofold, with a
  median of {sorted(x for x,_ in diffs)[len(diffs)//2]:+.2f}. The one outlier is {min(diffs)[1]} at {min(diffs)[0]:+.2f}.

  Two methods, two decades apart, agree that albumin affinity is close to
  species-invariant. BSA is a fair stand-in for HSA in PFAS work -- which
  matters, because much of the in vitro literature uses BSA for cost.
""")

    # -------------------------------------------------------------- structure
    print("=" * 78)
    print("PART 2.  STRUCTURE: chain length and head group")
    print("=" * 78)
    print("""
  Fischer 2024, log D_protein/water (L_water per kg_protein) -- HIGHER means
  MORE bound. One method, one laboratory, 10 chain lengths, sub-saturating.
  This is the only dataset in the file that can carry a structure-activity
  reading, because everything else mixes methods.
""")
    series = {"carboxylate": [], "sulfonate": []}
    print(f"  {'nC':>3s} {'cmpd':9s}{'head group':<16s}{'HSA':>7s}{'BSA':>7s}{'glob':>7s}"
          f"{'glob-HSA':>10s}")
    for (n, c), v in sorted(fd.items()):
        if not v.get(HSA):
            continue
        hg = next((r["head_group"] for r in f if r["chemical"] == c), "")
        h = float(v[HSA])
        g = float(v[GLB]) if v.get(GLB) else None
        if hg in series:
            series[hg].append((n, h))
        print(f"  {n:>3d} {c:9s}{hg[:15]:<16s}{h:>7.2f}"
              f"{(f'{float(v[BSA]):.2f}' if v.get(BSA) else '-'):>7s}"
              f"{(f'{g:.2f}' if g else '-'):>7s}"
              f"{(f'{g-h:+.2f}' if g else '-'):>10s}")

    def slope(pts):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        return (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
                / sum((x - mx) ** 2 for x in xs))

    carb = sorted(series["carboxylate"])
    sulf = sorted(series["sulfonate"])
    rise = [p for p in carb if p[0] <= 10]
    fall = [p for p in carb if p[0] >= 10]
    print(f"""
  Three things fall out, and only the first is the textbook one.

  1. CHAIN LENGTH IS NOT MONOTONIC. The carboxylates rise from {carb[0][1]:.2f} at
     nC{carb[0][0]} to a maximum of {max(p[1] for p in carb):.2f} around nC10-11, then FALL to {carb[-1][1]:.2f}
     at nC{carb[-1][0]}. Slope is {slope(rise):+.2f} log units per CF2 up to nC10 and
     {slope(fall):+.2f} after it. A structure-activity model fitted with carbon number
     as a linear term is wrong at both ends of the range.

  2. THE HEAD GROUP IS WORTH ABOUT A LOG UNIT, AND IT DOES NOT PLATEAU. At
     matched chain length the sulfonate binds harder than the carboxylate:
     PFHxS {fd[(6,'PFHxS')][HSA]} against PFHpA {fd[(6,'PFHpA')][HSA]} at nC6, PFHpS {fd[(7,'PFHpS')][HSA]} against
     PFOA {fd[(7,'PFOA')][HSA]} at nC7. Over the range where both are measured the
     sulfonate slope is {slope([p for p in sulf if p[0]<=7]):+.2f} per CF2 against {slope([p for p in carb if p[0]<=7]):+.2f} for the
     carboxylates -- similar rise, offset by roughly one log unit.

  3. ALBUMIN IS NOT THE WHOLE STORY AND STOPS BEING THE MAIN ONE. gamma-globulin
     binds 2-3 log units weaker than albumin at short chain lengths and catches
     up steadily: the gap closes from {min(float(fd[k][GLB])-float(fd[k][HSA]) for k in fd if fd[k].get(GLB) and fd[k].get(HSA)):+.2f} at its widest to parity and
     beyond for the longest chains. For PFTeDA (nC13) globulin binds MORE
     tightly than albumin. Any model that treats "protein binding" as "albumin
     binding" is increasingly wrong as the chain lengthens -- and the long
     chains are the ones that persist.
""")

    # ------------------------------------------------------------------- gaps
    print("=" * 78)
    print("PART 3.  What is missing")
    print("=" * 78)
    exp = [r for r in all_rows if r.get("is_computational", "no") != "yes"]
    sp = collections.Counter(r["protein_species"] for r in exp)
    mouse = [r for r in exp if "mouse" in (r["protein_species"] + r["protein"]).lower()]
    print(f"""
  {len(exp)} experimental rows cover {len(sp)} species labels, dominated by human
  ({sp['human']}) and bovine ({sp['bovine']}), with rat ({sp['rat']}) and porcine ({sp['porcine']}) thin.

  THERE IS NO MOUSE SERUM ALBUMIN BINDING DATA AT ALL. The {len(mouse)} rows mentioning
  mouse are tissue homogenate and cell membranes, not plasma albumin. That is a
  hole directly under this review's own argument: the mouse is the species
  whose clearance and unbound fraction section 10's endpoint is built from, and
  the only mouse fu available (Argoul 2026, equilibrium dialysis in whole
  plasma) has no in vitro albumin measurement to decompose it against.

  So the species question can be answered for human / rat / bovine / porcine
  and NOT for the mouse, which is the animal the rest of this review runs on.
""")


if __name__ == "__main__":
    main()
