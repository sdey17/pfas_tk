#!/usr/bin/env python3
"""Allendorf 2019 batch A7906: wrong conditions, and four rows missing.

WHAT WENT WRONG.  add_partition_bischel_allendorf.py ingested SI Table
SI-2.3, which carries the dialysis concentrations for batch A3803, and
noted "All of SI-2.3 is batch A3803".  True, but it led to treating
A7906 as though the paper reported no conditions for it.  It does:
Table SI-2.4.2 is the exact counterpart of SI-2.3 for A7906, with its
own c_ref, c_free, c_BSA and v.  It was never read.

Three consequences, all verified against the SI PDF fetched from
pubs.rsc.org (c9em00290a1_suppl):

  PROTEIN CONCENTRATION wrong on 5 of 10 A7906 rows, because the
  A3803 c_BSA was reused.  A7906 used a different protein level for
  six compounds -- PFBA 10 g/L not 3, PFHpA and DONA 0.5 not 0.1,
  PFBS 2 not 1, PFNA 0.1 not 0.05.

  FRACTION BOUND wrong on all 10, also copied from A3803.  SI-2.4.1
  reports its own, and they differ a lot: PFHxA 35% not 58%, DONA
  61% not 30%, PFHpA 62% not 34%.

  LIGAND:PROTEIN RATIO wrong on all 10.  They carried the design
  ceiling 0.1 because no per-compound value was thought to exist.
  SI-2.4.2 gives one for every compound, 1.5e-4 to 2.7e-2.

  FOUR COMPOUNDS MISSING.  SI-2.4.1 reports A7906 values for PFHxS,
  PFOS, 9Cl-PF3ONS and PFECHS; the ingestion left those slots empty.
  Allendorf is 28 rows, not 24.

The A3803 half was checked against SI-2.3 at the same time and is
correct in every cell, including PFECHS's v = 1.5e-3 -- which does
not reconcile with c_ref/c_BSA by a factor of 8, but is what the
paper prints.  That anomaly is Allendorf's, not this project's, and
is left as published with a note.

Run:  python3 tk_review/scripts/fix_allendorf_a7906.py [--dry-run]
Idempotent: corrects rows in place and adds the missing four only if
they are absent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

BATCH_A7906 = ("bovine serum albumin (BSA), batch A7906 (heat-shock "
               "fraction, NOT defatted; SI Table SI-2.4.1). SI-2.4 "
               "states the paper could not determine whether the "
               "batch difference comes from residual fatty acids or "
               "from conformational change during extraction")
METHOD = ("equilibrium dialysis in a physiologically representative "
          "system at 37 C, compound:BSA molar ratio held below 0.1; "
          "LC-MS/MS with a mass-balance calculation")

# SI Table SI-2.4.2 + SI-2.4.1, batch A7906.
# chemical: (c_ref, SD, c_free, SD, c_BSA g/L, v, logK, SD, Ka,
#            f_bound %, SD, nCF, head_group)
A7906 = {
 "PFBA":      (8.42, 0.42, 6.00, 0.29, "1e1",   "1.5e-4", "1.90", "0.07",
               "5.4e3", "29", "4", "3",  "carboxylate"),
 "PFHxA":     (8.49, 0.54, 5.54, 0.30, "1e0",   "1.2e-3", "3.02", "0.07",
               "7.0e4", "35", "4", "5",  "carboxylate"),
 "PFHpA":     (7.16, 0.35, 2.74, 0.19, "5e-1",  "3.1e-3", "3.79", "0.05",
               "4.1e5", "62", "3", "6",  "carboxylate"),
 "PFOA":      (12.10, 0.70, 6.93, 0.38, "1e-1", "1.6e-2", "4.17", "0.05",
               "1.0e6", "43", "3", "7",  "carboxylate"),
 "PFNA":      (6.11, 0.36, 3.34, 0.17, "1e-1",  "1.5e-2", "4.48", "0.05",
               "2.1e6", "45", "3", "8",  "carboxylate"),
 "PFDA":      (4.11, 0.26, 1.51, 0.12, "5e-2",  "1.4e-2", "4.84", "0.05",
               "4.7e6", "63", "3", "9",  "carboxylate"),
 "PFUnDA":    (2.99, 0.44, 1.67, 0.26, "2.5e-2", "1.0e-2", "4.72", "0.1",
               "3.6e6", "44", "6", "10", "carboxylate"),
 "HFPO-DA (GenX)": (9.30, 0.44, 5.24, 0.23, "3e0", "5.4e-4", "2.71",
               "0.04", "3.4e4", "44", "2", "",
               "ether carboxylate (polyfluoroether)"),
 "ADONA":     (9.70, 1.02, 3.79, 0.20, "5e-1",  "4.2e-3", "3.80", "0.04",
               "4.2e5", "61", "2", "",
               "ether carboxylate (polyfluoroether)"),
 "PFBS":      (8.42, 0.42, 6.00, 0.29, "2e0",   "7.0e-4", "2.65", "0.12",
               "3.1e4", "32", "6", "4",  "sulfonate"),
 # --- the four the ingestion left out ---
 "PFHxS":     (9.20, 0.49, 5.18, 0.75, "5e-2",  "2.7e-2", "4.49", "0.13",
               "2.2e6", "44", "7", "6",  "sulfonate"),
 "PFOS":      (5.71, 0.41, 1.32, 0.24, "1e-1",  "1.2e-2", "4.83", "0.11",
               "4.7e6", "77", "4", "8",  "sulfonate"),
 "9Cl-PF3ONS": (4.59, 0.10, 1.94, 0.11, "2.5e-2", "2.5e-2", "5.02",
               "0.03", "7.2e6", "58", "2", "8", "ether sulfonate"),
 "PFECHS":    (6.69, 1.09, 3.33, 0.11, "1e-1",  "9.6e-4", "4.30", "0.03",
               "1.3e6", "50", "2", "8",  "cyclic sulfonate"),
}


def fields(name, d):
    cr, crs, cf, cfs, cbsa, v, lk, lks, ka, fb, fbs, ncf, head = d
    return dict(
        protein_conc=cbsa,
        ligand_conc_range=(f"c_ref {cr:.2f} +/- {crs:.2f} ng/mL, "
                           f"c_free {cf:.2f} +/- {cfs:.2f} ng/mL at "
                           f"equilibrium (SI Table SI-2.4.2)"),
        ligand_conc_units="ng/mL",
        ligand_over_protein_max=v,
        source_table=(f"SI Table SI-2.4.1 | log K_albumin/water {lk} "
                      f"+/- {lks} | Ka {ka} M^-1, converted from the "
                      f"partition coefficient at one molar ratio, NOT "
                      f"fitted from a binding isotherm | fraction "
                      f"bound {fb} +/- {fbs}%"),
        ka=lk)


def new_row(name, d, template):
    r = dict(template)
    r.update(chemical=name, n_perfluorinated_carbons=d[11],
             head_group=d[12], protein=BATCH_A7906)
    r.update(fields(name, d))
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    a7 = [r for r in rows if r["study"].startswith("Allendorf")
          and "SI-2.4.1" in r["source_table"]]
    template = dict(a7[0])
    have = {r["chemical"] for r in a7}

    fixed = added = 0
    for r in a7:
        d = A7906.get(r["chemical"])
        if not d:
            print(f"  ?? no SI-2.4.2 entry for {r['chemical']}")
            continue
        upd = fields(r["chemical"], d)
        changed = [k for k, v in upd.items() if r[k] != v]
        if changed:
            fixed += 1
            print(f"  fix   {r['chemical']:16} {', '.join(changed)}")
            for k, v in upd.items():
                r[k] = v

    tail = max(i for i, r in enumerate(rows)
               if r["study"].startswith("Allendorf"))
    for name in reversed(list(A7906)):
        if name in have:
            continue
        rows.insert(tail + 1, new_row(name, A7906[name], template))
        added += 1
        print(f"  ADD   {name:16} A7906, log K {A7906[name][6]}")

    print(f"\n  {fixed} rows corrected, {added} rows added")
    n = sum(1 for r in rows if r["study"].startswith("Allendorf"))
    print(f"  Allendorf is now {n} rows")
    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
