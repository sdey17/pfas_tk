#!/usr/bin/env python3
"""
Fill the assay-condition columns of db/protein_binding.csv for Starnes 2024,
and add a derived ligand:protein molar ratio column.

Why the ratio matters: a dissociation constant is only well determined when it
is comparable to, or larger than, the protein concentration in the assay. If
Kd << [protein] the assay is in the tight-binding regime and reports [protein]/2
rather than the true Kd -- the measurement has a FLOOR. Recording [protein] and
the ligand range alongside every Ka/Kd makes that floor visible.

Starnes et al. 2024 (Toxicol Sci 199:132-149), Methods, verbatim:
  "Final protein concentrations for each assay were 0.13 mM for HSA, 0.18 mM
   for BSA, 0.11 mM for PSA, and 0.18 mM for RSA."
  Assay volume 20 uL in sealed 96-well plates; PFAS serially diluted in HBS;
  stocks 10-20 mM (Table 1).

Run:  python3 tk_review/scripts/add_binding_assay_conditions.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

# Starnes 2024 Methods: final assay protein concentration, by species
STARNES2024_PROTEIN_MM = {"human": 0.13, "bovine": 0.18,
                          "porcine": 0.11, "rat": 0.18}
STARNES2024_NOTE = ("20 uL assay volume, sealed 96-well plate; albumin "
                    "concentration optimised per species for signal-to-noise")


def to_uM(value, units):
    f = {"mm": 1e3, "um": 1.0, "nm": 1e-3, "m": 1e6, "mg/ml": None}
    return None if not value else (lambda k: None if k is None else value * k)(
        f.get((units or "").strip().lower()))


def main():
    rows = list(csv.DictReader(open(DB)))
    fields = list(rows[0].keys())
    for c in ("assay_volume_note", "ligand_over_protein_max"):
        if c not in fields:
            fields.append(c)

    filled = 0
    for r in rows:
        if r["study"].startswith("Starnes") and r["year"] == "2024":
            p = STARNES2024_PROTEIN_MM.get(r["protein_species"])
            if p and not r["protein_conc"].strip():
                r["protein_conc"] = "%.2f" % p
                r["protein_conc_units"] = "mM"
                r["assay_volume_note"] = STARNES2024_NOTE
                filled += 1

    # derived: highest ligand concentration / protein concentration, molar
    for r in rows:
        r.setdefault("assay_volume_note", "")
        r["ligand_over_protein_max"] = ""
        try:
            p_uM = to_uM(float(r["protein_conc"]), r["protein_conc_units"])
        except (ValueError, TypeError):
            p_uM = None
        hi = None
        rng = (r.get("ligand_conc_range") or "").replace("~", "").strip()
        if rng:
            parts = [x for x in rng.replace("to", "-").split("-") if x.strip()]
            try:
                hi = to_uM(float(parts[-1]), r.get("ligand_conc_units"))
            except (ValueError, TypeError, IndexError):
                hi = None
        if p_uM and hi:
            r["ligand_over_protein_max"] = "%.2f" % (hi / p_uM)

    with open(DB, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    print("filled protein_conc on %d Starnes 2024 rows" % filled)
    for c in ("protein_conc", "ligand_conc_range", "temperature_c",
              "buffer", "ligand_over_protein_max"):
        k = sum(1 for r in rows if (r.get(c) or "").strip())
        print("  %-24s %3d / %d  (%.0f%%)" % (c, k, n, 100 * k / n))


if __name__ == "__main__":
    main()
