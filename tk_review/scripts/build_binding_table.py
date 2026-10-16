#!/usr/bin/env python3
"""Summarise the in vitro PFAS-protein binding database into readable tables.

Two quantities live in the source file and must NOT be pooled:

  * an association constant Ka, in M^-1 (or L/mol), from fluorescence quenching,
    ITC, equilibrium dialysis and the like;
  * a protein/water distribution ratio log D_protein/w, in L_water per
    kg_protein, from solid-phase microextraction (Fischer et al. 2024).

They answer different questions and are numerically incomparable. The source
file labels its units explicitly; this script honours that label and emits one
table per quantity.

Within the Ka table the spread is still large, and the summary reports a range
rather than a central value, because the reported value for one chemical-protein
pair is driven by the ligand:protein ratio and the fitting model more than by the
protein.
"""
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DB = HERE.parent / "db"

SCI = re.compile(r"([0-9]*\.?[0-9]+)\s*(?:[x×*]\s*10\s*\^?|[eE])\s*([-+]?[0-9]+)")
BARE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*$")

KA_UNITS = {"m^-1", "l/mol (= m^-1)", "m-1", "l/mol"}
LOGKA_UNITS = {"log m^-1"}
LOGD_UNITS_PREFIX = "log d_protein/w"
KD_SCALE = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "µm": 1e-6, "nm": 1e-9, "pm": 1e-12}


def parse(v):
    if not v:
        return None
    v = v.strip().replace("−", "-")
    m = SCI.search(v)
    if m:
        try:
            return float(m.group(1)) * 10 ** float(m.group(2))
        except ValueError:
            return None
    m = BARE.match(v)
    return float(m.group(1)) if m else None


def ka_of(row):
    """Return Ka in M^-1 for a row, or None if the row is not an association
    constant (or is computational, or is a distribution ratio)."""
    if row.get("is_computational", "").strip().lower() == "yes":
        return None
    units = row.get("ka_units", "").strip().lower()
    ka = parse(row.get("ka", ""))
    if ka is not None:
        if units in LOGKA_UNITS:
            return 10 ** ka
        # Accept a unit string that STARTS with M^-1 even if it carries a note.
        if units in KA_UNITS or units.startswith("m^-1"):
            return ka
        return None                      # log D and anything unlabelled
    kd = parse(row.get("kd", ""))
    if kd and kd > 0:
        u = row.get("kd_units", "").strip().lower().replace("μ", "µ")
        scale = KD_SCALE.get(u)
        if scale:
            return 1.0 / (kd * scale)
    return None


def summarise(groups, out_path, value_name, fmt="{:.3g}"):
    out = []
    for (chem, prot), vals in sorted(groups.items()):
        xs = [v[0] for v in vals]
        out.append({
            "chemical": chem, "protein": prot,
            "n_measurements": len(vals),
            "n_methods": len(sorted({v[1][:46] for v in vals})),
            f"{value_name}_min": fmt.format(min(xs)),
            f"{value_name}_max": fmt.format(max(xs)),
            f"{value_name}_geomean": fmt.format(
                math.exp(sum(math.log(x) for x in xs) / len(xs))) if min(xs) > 0
                else "",
            "fold_spread": f"{max(xs)/min(xs):.3g}" if min(xs) > 0 else "",
            "methods": " | ".join(sorted({v[1][:46] for v in vals})),
            "studies": " | ".join(sorted({f"{v[2]} {v[3]}" for v in vals})),
        })
    out.sort(key=lambda r: -(float(r["fold_spread"]) if r["fold_spread"] else 0))
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    return out


def main() -> None:
    rows = list(csv.DictReader((DB / "protein_binding.csv").open()))

    ka_groups, logd_groups = defaultdict(list), defaultdict(list)
    skipped, funbound = 0, 0
    for r in rows:
        key = (r["chemical"], r["protein"])
        ka = ka_of(r)
        if ka and ka > 0:
            ka_groups[key].append((ka, r["method"], r["study"], r["year"]))
            continue
        units = r.get("ka_units", "").strip().lower()
        v = parse(r.get("ka", ""))
        if v is not None and units.startswith(LOGD_UNITS_PREFIX):
            logd_groups[key].append((v, r["method"], r["study"], r["year"]))
        elif parse(r.get("unbound_fraction", "")) is not None:
            # An unbound fraction in whole tissue or plasma is a real
            # measurement, not a failed affinity. Counting it as "skipped"
            # made 320 Ryu 2024 rows look like junk.
            funbound += 1
        else:
            skipped += 1

    ka_out = summarise(ka_groups, DB / "binding_summary_ka.csv", "ka_M-1")
    logd_out = summarise(logd_groups, DB / "binding_summary_logD.csv",
                         "logD", fmt="{:.2f}")

    print(f"association constants : {sum(len(v) for v in ka_groups.values())} "
          f"values over {len(ka_out)} chemical-protein pairs "
          f"-> {DB/'binding_summary_ka.csv'}")
    print(f"protein/water log D   : {sum(len(v) for v in logd_groups.values())} "
          f"values over {len(logd_out)} pairs "
          f"-> {DB/'binding_summary_logD.csv'}")
    print(f"unbound fraction      : {funbound} rows (whole tissue or plasma; "
          f"an f_u, not an affinity constant)")
    print(f"not a usable affinity : {skipped} rows (computational, relative "
          f"potency, IC50, rank order, or no units)")
    big = Counter(f"{r['study']} {r['year']}" for r in rows).most_common(3)
    print("  file is not uniform - largest contributors: "
          + ", ".join(f"{k} ({v})" for k, v in big))
    print("  group by study before taking any average over it.\n")

    multi = [r for r in ka_out if r["n_measurements"] > 1]
    print("ASSOCIATION CONSTANTS, pairs measured more than once "
          "(worst disagreement first)")
    print(f"  {'chemical':<9} {'protein':<28} {'n':>2} "
          f"{'Ka range (M^-1)':>24} {'spread':>8}")
    print("  " + "-" * 78)
    for r in multi[:12]:
        rng = f"{r['ka_M-1_min']} - {r['ka_M-1_max']}"
        print(f"  {r['chemical'][:9]:<9} {r['protein'][:28]:<28} "
              f"{r['n_measurements']:>2} {rng:>24} "
              f"{float(r['fold_spread']):>7.0f}x")
    print(f"\n  {len(multi)} of {len(ka_out)} pairs have more than one "
          "measurement.")
    print("  Read the range, not a central value: Alesio & Bothun 2022 obtained")
    print("  Ka from 0.07 to 6.16 x 10^6 M^-1 from ONE dataset using three")
    print("  different binding models. Physiological PFAS:albumin molar ratios")
    print("  are <= 0.0005; most in vitro work runs orders of magnitude above.")

    print("\nPROTEIN/WATER DISTRIBUTION, Fischer et al. 2024 (log D, L_water/kg)")
    print("  measured at a PFAS:protein molar ratio <= 0.004, the only wide-")
    print("  coverage dataset at a physiologically realistic ratio\n")
    print(f"  {'chemical':<10} {'protein':<34} {'log D':>7}")
    print("  " + "-" * 54)
    want = ("human serum albumin (HSA)", "gamma-globulin")
    for chem in ("PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFHxS", "PFOS"):
        for prot in want:
            vals = logd_groups.get((chem, prot))
            if vals:
                print(f"  {chem:<10} {prot:<34} {vals[0][0]:>7.2f}")


if __name__ == "__main__":
    main()
