#!/usr/bin/env python3
"""Close the K_MW basis caveat in db/protein_binding.csv.

The 60 Chen et al. 2025 membrane-water rows were annotated "SI eq S1 defines
it and is not in hand, so not converted".  The SI was obtained on 2026-10-22
and eq S1 reads

    K_MW = (C_spike - C_PFAS,W) * V_W / (C_PFAS,W * M_PL)

with the authors noting that a unit conversion is applied so the reported
coefficient is dimensionless.  The measurement uses a TRANSIL solid-supported
lipid membrane kit, which supplies the membrane as a stated VOLUME -- 0.067 to
2.166 uL of bilayer per 500 uL tube -- so the lipid quantity is natively
volumetric and the mass and volume bases coincide at a phospholipid density of
1 g/cm3.  The reported log K_MW is therefore a dimensionless concentration
ratio on a volume-of-membrane to volume-of-water basis.  Real phospholipid
densities of 1.01-1.06 g/cm3 bound the residual basis error at <= 0.03 log
units.

Two checks support the reading: the tabulated series reproduces the authors'
own unweighted regression of log K_MW on chain length (+0.361, R2 0.983 for
the carboxylates; +0.374, R2 0.985 for the sulfonates, against their stated
+0.361 +/- 0.014 and +0.374 +/- 0.019), and the absolute values sit in the
published range for PFAS membrane-water partitioning by the same technique.

The one caveat that remains is pH: Chen reports at 7.0, and anything
evaluating a balance at 7.4 is extrapolating 0.4 units.

Run:  python3 tk_review/scripts/update_chen_kmw_basis.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

OLD = ("log K_MW (phospholipid membrane-water, pH 7.0); basis as published - "
       "SI eq S1 defines it and is not in hand, so not converted")
NEW = ("log K_MW (phospholipid membrane-water, pH 7.0); dimensionless "
       "volume-of-membrane / volume-of-water concentration ratio, resolved "
       "from SI eq S1 + the TRANSIL kit's stated bilayer volumes at "
       "rho_lipid = 1 g/cm3 (residual basis error <= 0.03 log units); "
       "directly usable as a volumetric D_MW; pH 7.0 not 7.4")


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    fields = list(rows[0])
    n = 0
    for r in rows:
        if r["ka_units"] == OLD:
            r["ka_units"] = NEW
            n += 1
    if not n:
        print("nothing to update (already applied?)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"updated {n} Chen 2025 K_MW rows in {DB}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
