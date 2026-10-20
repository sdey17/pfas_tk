#!/usr/bin/env python3
"""Chen et al. 2025: the 60 phospholipid membrane-water partition coefficients.

These values were already transcribed from SI Table S5 -- they have been
sitting inside the source_table note on every Chen HSA row, where nothing
could use them. This promotes them to rows of their own, because a serum
mass balance needs a lipid limb as well as a protein limb, and K_MW is the
one fixed partition coefficient Chen reports.

Why K_MW and not D_HSA/W. Chen measures two things and derives a third:

  K_MW      membrane-water, by solid-supported lipid membrane (SSLM).
            A true partition coefficient -- the paper is explicit that it
            "does not change with PFAS concentrations in the water phase
            and thus the log K_MW is a fixed intercept".
  K_D       HSA binding, by equilibrium dialysis. Already in the file.
  D_HSA/W   protein-water, DERIVED from the binding curve by eq 4, and
            concentration-dependent by construction: "The D_HSA/W for each
            PFAS decreases when the concentration of PFAS in water
            increases, which is a feature of saturable specific binding."

So there is no single Chen protein-water number to transcribe. Eq 4 needs
B_max and a chosen aqueous concentration, and B_max is in SI Table S5,
which is not in hand. Only K_MW is a constant, and only K_MW is added here.

The basis of K_MW (per kg lipid or per L lipid) is defined in SI eq S1,
which is also not in hand, so the unit is recorded as the paper labels it
rather than converted.

Run:  python3 scripts/add_chen2025_membrane_kmw.py [--dry-run]
"""
import argparse
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

KMW = re.compile(r"log K_MW \(phospholipid membrane-water, pH 7\.0\) "
                 r"([-\d.]+)(?: \+/- ([\d.]+))?")

METHOD = ("solid-supported lipid membrane (SSLM; TRANSIL Membrane Affinity "
          "Kit, phosphatidylcholine beads), 4 h at 50 rpm; aqueous phase by "
          "LC-MS/MS and membrane concentration by difference")
UNITS = ("log K_MW (phospholipid membrane-water, pH 7.0); basis as published "
         "- SI eq S1 defines it and is not in hand, so not converted")
TK = ("the lipid limb of a serum mass balance: the sorptive capacity of "
      "membrane phospholipid, measured on the same 60 compounds and in the "
      "same study as the HSA binding constants")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, list(reader)

    existing = {(r["study"], r["year"], r["chemical"], r["protein"])
                for r in rows}
    new = []
    for r in rows:
        if r["study"] != "Chen et al." or r["year"] != "2025":
            continue
        if "membrane" in r["protein"]:
            continue
        m = KMW.search(r["source_table"])
        if not m:
            print(f"  no K_MW on the {r['chemical']} row -- skipped")
            continue
        val, se = m.group(1), m.group(2)
        protein = "phospholipid bilayer membrane (solid-supported, SSLM)"
        if (r["study"], r["year"], r["chemical"], protein) in existing:
            continue
        row = dict.fromkeys(fields, "")
        row.update(
            study=r["study"], year=r["year"], pmid=r["pmid"], doi=r["doi"],
            chemical=r["chemical"],
            n_perfluorinated_carbons=r["n_perfluorinated_carbons"],
            head_group=r["head_group"],
            protein=protein,
            protein_species="not applicable (synthetic phosphatidylcholine)",
            method=METHOD,
            ka=val, ka_units=UNITS,
            n_binding_sites="not applicable (partitioning, not site binding)",
            protein_conc="0.067-2.166 uL lipid per vial (6-point series)",
            protein_conc_units="",
            ligand_conc_range="one spike level per PFAS group",
            ligand_conc_units="",
            temperature_c="room temperature (not stated more precisely)",
            buffer="10 mM ammonium acetate, pH 7.0",
            is_computational="no",
            tk_relevance_noted=TK,
            source_table=("SI Table S5 | log K_MW " + val
                          + (f" +/- {se}" if se else "")
                          + " | measured at pH 6 and pH 7 in triplicate; the "
                            "pH 7.0 value is taken"))
        new.append(row)

    print(f"{len(new)} membrane-water rows to add")
    if args.dry_run:
        for r in new[:5]:
            print(f"   {r['chemical']:12s} log K_MW {r['ka']}")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows + new)
    print(f"wrote {os.path.abspath(CSV)}  ({len(rows) + len(new)} rows total)")


if __name__ == "__main__":
    main()
