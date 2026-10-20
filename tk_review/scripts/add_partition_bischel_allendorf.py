#!/usr/bin/env python3
"""Per-compound protein-water partition coefficients: Bischel 2011, Allendorf 2019.

Both were in the file as a single summary row quoting only the range from
their abstracts ("log K_protein-water = 3.3-4.3", "log K_BSA/w = 2.8-4.8").
The tables are now in hand, so those summaries are replaced by the rows.

Bischel et al. 2011, Environ Toxicol Chem 30(11):2423-2430, Table 1.
  Equilibrium dialysis against 200 uM BSA, reservoir PFAA 0.009-0.68 uM, so
  a ligand:protein ratio of at most 0.0034 -- the specific-binding regime.
  K_PW is defined in the paper as
      K_PW = f_bound / (C_W * rho_albumin * [P] * (1 - f_bound))
  with rho_albumin the partial specific volume, 0.733 mL/g, giving units of
  "(g bound PFAA/mL BSA)/(g free PFAA/mL water)". That is a VOLUME basis,
  L_water per L_albumin, the same basis Qin 2024 uses and 0.135 log units
  above the mass basis Fischer and Allendorf report.

Allendorf et al. 2019, Environ Sci Process Impacts 21:1852-1863, Table 1.
  Dialysis at 37 C with the compound:BSA molar ratio held below 0.1, so also
  the specific regime. Reported directly as log K_albumin/water in
  L_water per kg_albumin -- a MASS basis, no conversion needed.

  Allendorf measured two BSA batches and they do not agree. Table 1 uses
  A3803, essentially fatty-acid-free, which the authors chose for the
  headline values. SI Table SI-2.4.1 uses A7906, heat-shock fraction, and
  runs up to 0.62 log units lower at short chain while agreeing to within
  0.11 at nCF 9-10. Both are ingested, flagged by batch, because a
  batch-to-batch spread that large on the same method in the same lab is
  itself a datum about how far these numbers can be trusted.

  The Ka column in both of Allendorf's tables is "converted from
  K_albumin/water at one molar ratio and not determined by measuring the
  binding isotherms", so it is recorded in the source note rather than as a
  measured association constant.

Run:  python3 scripts/add_partition_bischel_allendorf.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

# Bischel Table 1: (analyte, nCF, head, log K_PW, SE, fraction bound %, +/-,
#                   reservoir range uM)
BISCHEL = [
    ("PFBS",   4,  "sulfonate",   "3.86", "0.07", "99.0", "0.5", "0.010 to 0.57"),
    ("PFHxS",  6,  "sulfonate",   "4.3",  "0.1",  "99.2", "0.5", "0.012 to 0.27"),
    ("PFOS",   8,  "sulfonate",   "4.1",  "0.1",  "99.1", "0.4", "0.013 to 0.26"),
    ("PFPeA",  4,  "carboxylate", "3.40", "0.02", "96.6", "0.8", "0.024 to 0.56"),
    ("PFHxA",  5,  "carboxylate", "4.05", "0.02", "99.2", "0.3", "0.009 to 0.43"),
    ("PFHpA",  6,  "carboxylate", "4.23", "0.08", "99.3", "0.1", "0.054 to 0.30"),
    ("PFOA",   7,  "carboxylate", "4.14", "0.04", "99.1", "0.2", "0.023 to 0.32"),
    ("PFNA",   8,  "carboxylate", "4.05", "0.08", "98.9", "0.3", "0.064 to 0.46"),
    ("PFDA",   9,  "carboxylate", "3.86", "0.08", "98",   "1",   "0.016 to 0.47"),
    ("PFUnDA", 10, "carboxylate", "3.7",  "0.2",  "95",   "3",   "0.026 to 0.38"),
    ("PFDoDA", 11, "carboxylate", "3.3",  "0.1",  "80",   "10",  "0.031 to 0.68"),
]

# Allendorf: (analyte, nCF or "", head, logK_A3803, SD, Ka_A3803, fbound%, SD,
#             cBSA g/L, logK_A7906, SD_A7906, Ka_A7906)
ALLENDORF = [
    ("PFBA",        3,  "carboxylate", "2.52", "0.11", "2.3e4", "34", "5",  "3e0",
     "1.90", "0.07", "5.4e3"),
    ("PFHxA",       5,  "carboxylate", "3.43", "0.04", "1.8e5", "58", "2",  "1e0",
     "3.02", "0.07", "7.0e4"),
    ("PFHpA",       6,  "carboxylate", "4.02", "0.08", "7.2e5", "34", "4",  "1e-1",
     "3.79", "0.05", "4.1e5"),
    ("PFOA",        7,  "carboxylate", "4.20", "0.05", "1.1e6", "44", "3",  "1e-1",
     "4.17", "0.05", "1.0e6"),
    ("PFNA",        8,  "carboxylate", "4.32", "0.07", "1.4e6", "37", "4",  "5e-2",
     "4.48", "0.05", "2.1e6"),
    ("PFDA",        9,  "carboxylate", "4.73", "0.04", "3.6e6", "57", "2",  "5e-2",
     "4.84", "0.05", "4.7e6"),
    ("PFUnDA",      10, "carboxylate", "4.60", "0.30", "3.1e6", "35", "13", "2.5e-2",
     "4.72", "0.1",  "3.6e6"),
    ("HFPO-DA (GenX)", "", "ether carboxylate (polyfluoroether)",
     "3.06", "0.03", "2.3e4", "64", "1", "3e0", "2.71", "0.04", "3.4e4"),
    ("ADONA",       "", "ether carboxylate (polyfluoroether)",
     "3.93", "0.03", "1.8e5", "30", "2", "1e-1", "3.80", "0.04", "4.2e5"),
    ("PFBS",        4,  "sulfonate",   "3.20", "0.11", "1.1e5", "44", "6",  "1e0",
     "2.65", "0.12", "3.1e4"),
    ("PFHxS",       6,  "sulfonate",   "4.81", "0.05", "4.4e6", "62", "2",  "5e-2",
     "", "", ""),
    ("PFOS",        8,  "sulfonate",   "4.67", "0.07", "3.2e6", "37", "4",  "2.5e-2",
     "", "", ""),
    ("9Cl-PF3ONS",  8,  "ether sulfonate", "5.01", "0.04", "6.9e6", "56", "2", "2.5e-2",
     "", "", ""),
    ("PFECHS",      8,  "cyclic sulfonate", "4.55", "0.03", "2.4e6", "66", "1", "1e-1",
     "", "", ""),
]

B_METHOD = ("equilibrium dialysis against BSA; K_PW from the linear "
            "regression of bound-phase against aqueous concentration")
B_UNITS = ("log K_PW (protein-water distribution coefficient), L_water per "
           "L_albumin -- VOLUME basis, from the 0.733 mL/g partial specific "
           "volume; 0.135 log units above the mass basis")
A_METHOD = ("equilibrium dialysis in a physiologically representative system "
            "at 37 C, compound:BSA molar ratio held below 0.1; LC-MS/MS with "
            "a mass-balance calculation")
A_UNITS = ("log K_albumin/water, L_water per kg_albumin -- MASS basis, as "
           "reported")


def blank(fields):
    return dict.fromkeys(fields, "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, list(reader)

    before = len(rows)
    rows = [r for r in rows
            if not (r["study"] == "Bischel et al." and r["year"] == "2011")
            and not (r["study"] == "Allendorf et al." and r["year"] == "2019")]
    print(f"removed {before - len(rows)} summary rows")

    new = []
    for (chem, ncf, head, lk, se, fb, fbsd, rng) in BISCHEL:
        r = blank(fields)
        r.update(
            study="Bischel et al.", year="2011", pmid="21842491",
            doi="10.1002/etc.647", chemical=chem,
            n_perfluorinated_carbons=str(ncf), head_group=head,
            protein="bovine serum albumin (BSA)", protein_species="bovine",
            method=B_METHOD, ka=lk, ka_units=B_UNITS,
            unbound_fraction="%.4g" % (1 - float(fb) / 100.0),
            n_binding_sites="not applicable (partitioning, not site binding)",
            protein_conc="200", protein_conc_units="uM",
            ligand_conc_range=rng, ligand_conc_units="uM",
            temperature_c="not stated for the dialysis arm",
            buffer="50 mM sodium phosphate, pH 7.0",
            is_computational="no",
            tk_relevance_noted=("protein-water partition coefficient measured "
                                "at a ligand:protein ratio of at most 0.0034, "
                                "the specific-binding regime"),
            source_table=(f"Table 1 | log K_PW {lk} +/- {se} | fraction bound "
                          f"{fb} +/- {fbsd}% to 200 uM BSA"),
            assay_volume_note="", ligand_over_protein_max="0.0034")
        new.append(r)

    for rec in ALLENDORF:
        (chem, ncf, head, lk, sd, ka, fb, fbsd, cbsa, lk7, sd7, ka7) = rec
        for batch, val, s, kav, src in (
                ("A3803 (essentially fatty-acid-free; the paper's Table 1)",
                 lk, sd, ka, "Table 1"),
                ("A7906 (heat-shock fraction; SI Table SI-2.4.1)",
                 lk7, sd7, ka7, "SI Table SI-2.4.1")):
            if not val:
                continue
            r = blank(fields)
            r.update(
                study="Allendorf et al.", year="2019", pmid="31475719",
                doi="10.1039/c9em00290a", chemical=chem,
                n_perfluorinated_carbons=str(ncf), head_group=head,
                protein=f"bovine serum albumin (BSA), batch {batch}",
                protein_species="bovine",
                method=A_METHOD, ka=val, ka_units=A_UNITS,
                unbound_fraction="",
                n_binding_sites="not applicable (partitioning, not site binding)",
                protein_conc=cbsa, protein_conc_units="g/L",
                ligand_conc_range="compound:BSA molar ratio held below 0.1",
                ligand_conc_units="mol/mol",
                temperature_c="37",
                buffer="physiologically representative dialysis buffer",
                is_computational="no",
                tk_relevance_noted=("protein-water partition coefficient on a "
                                    "mass basis, measured below a 0.1 "
                                    "ligand:protein ratio"),
                source_table=(f"{src} | log K_albumin/water {val} +/- {s} | "
                              f"Ka {kav} M^-1, converted from the partition "
                              f"coefficient at one molar ratio, NOT fitted "
                              f"from a binding isotherm | fraction bound "
                              f"{fb} +/- {fbsd}%"),
                ligand_over_protein_max="0.1")
            new.append(r)

    print(f"{len(new)} per-compound rows to add "
          f"({sum(1 for r in new if r['study'].startswith('Bischel'))} Bischel, "
          f"{sum(1 for r in new if r['study'].startswith('Allendorf'))} Allendorf)")
    if args.dry_run:
        for r in new[:4]:
            print(f"   {r['study']:18s} {r['chemical']:14s} {r['ka']}")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows + new)
    print(f"wrote {os.path.abspath(CSV)}  ({len(rows) + len(new)} rows)")


if __name__ == "__main__":
    main()
