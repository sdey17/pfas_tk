#!/usr/bin/env python3
"""Li et al. 2021 SI: per-congener HSA association constants for 26 PFAS.

Li Y, Yao J, Zhang J, Pan Y, Dai J, Ji C, Tang J.  (AFFF PFAS binding to
human serum albumin by equilibrium dialysis, with companion docking.)
Toxics 2021;9(3):63.  doi 10.3390/toxics9030063, PMID 33803062.
Open access.  Values from SI spreadsheet toxics-1134374-SI.xlsx, sheet
"S7.2 Ka_docking and exp", column "log Ka (M-1)" with its SD.

This database held ONE aggregate row for the paper -- "26 PFAS
standards from a 1999 3M AFFF, range 0.10e5 to 3.16e5 M^-1" -- which is
the abstract's range and nothing per compound.  The 26 individual
constants are now in.

ONLY THE EXPERIMENTAL COLUMN IS INGESTED.  The same sheet carries two
docking columns, geometric-mean log Ka against PDB 1E7G and 1AO6.
Those are computed, they run up to 7.5 where the measurements top out
near 5.5, and they rise monotonically with chain length where the
measurements turn over.  They are not binding data and are left out.

THE MEASUREMENTS TURN OVER, which is the result worth having.  Among
the carboxylates log Ka rises from PFBA 4.10 to PFHpA 5.53 and then
falls away to PFTeDA 3.96 -- the same non-monotonic shape Zhang 2013
found on L-FABP, Allendorf found on albumin and Fan found on rat
L-FABP.  The docking columns in this very spreadsheet miss it, which is
a concrete limit on structure-only prediction and relevant to the QSAR
aim.

HFPO-DA is docked but carries NA in the experimental column, so it gets
no row here.

CONDITIONS, partial.  The setup is equilibrium dialysis with four
trials and six cells per trial, "adapted and adjusted from the
published supplementary information by Allendorf et al".  HSA at
600 uM is stated for the AFFF bound-fraction work (Figure S1.4).  The
ligand concentrations for the 26 standards are not stated in the main
text or in the SI sections read, so that column stays empty.

Run:  python3 scripts/add_li2021_hsa_si.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Li et al.", "2021"
DOI, PMID = "10.3390/toxics9030063", "33803062"
PROTEIN = "human serum albumin (HSA)"
METHOD = ("equilibrium dialysis, n = 4 trials with six cells each "
          "(blank and control cells for contamination, membrane "
          "fouling, equilibrium and surface adsorption); setup adapted "
          "from the Allendorf et al. 2019 supplementary information; "
          "LC-QTOF-MS")
PROT = ("600 uM HSA (stated for the AFFF bound-fraction experiment, "
        "SI Figure S1.4)", "uM")
LIG = ("not stated for the 26 standards in the main text or in the SI "
       "sections read", "")

# SI sheet S7.2, experimental columns only.
# (chemical, nCF or None, head group, log Ka, SD, fitting model)
S72 = [
    ("PFBA",    3,  "carboxylate", "4.097", "0.490", "Langmuir isotherm"),
    ("PFPeA",   4,  "carboxylate", "4.114", "0.277", "1st order association"),
    ("PFHxA",   5,  "carboxylate", "4.596", "0.246", "1st order association"),
    ("PFHpA",   6,  "carboxylate", "5.528", "0.385", "1st order association"),
    ("PFOA",    7,  "carboxylate", "5.046", "0.142", "1st order association"),
    ("PFNA",    8,  "carboxylate", "5.365", "0.202", "1st order association"),
    ("PFDA",    9,  "carboxylate", "4.925", "0.078", "1st order association"),
    ("PFUnDA",  10, "carboxylate", "4.625", "0.235", "1st order association"),
    ("PFDoDA",  11, "carboxylate", "4.536", "0.486", "1st order association"),
    ("PFTrDA",  12, "carboxylate", "4.349", "0.505", "1st order association"),
    ("PFTeDA",  13, "carboxylate", "3.962", "0.865", "1st order association"),
    ("PFOS (linear)", 8, "sulfonate", "4.744", "0.537", "Langmuir isotherm"),
    ("PFOS (branched)", 8, "sulfonate", "4.528", "0.266",
     "Langmuir isotherm"),
    ("PFBS",    4,  "sulfonate", "4.603", "0.133", "1st order association"),
    ("PFPeS",   5,  "sulfonate", "4.752", "0.526", "1st order association"),
    ("PFHxS",   6,  "sulfonate", "4.988", "0.438", "1st order association"),
    ("PFHpS",   7,  "sulfonate", "4.632", "0.339", "1st order association"),
    ("PFPrS (linear)", 3, "sulfonate", "4.068", "0.178",
     "1st order association"),
    ("PFNS",    9,  "sulfonate", "5.093", "0.153", "1st order association"),
    ("PFDS",    10, "sulfonate", "4.599", "0.245", "1st order association"),
    ("FOSA",    8,  "sulfonamide", "4.769", "0.075",
     "1st order association"),
    ("MeFOSAA", 8,  "sulfonamidoacetate", "5.148", "0.158",
     "1st order association"),
    ("EtFOSAA", 8,  "sulfonamidoacetate", "5.148", "0.155",
     "1st order association"),
    ("4:2 FTS", None, "fluorotelomer sulfonate", "4.920", "0.124",
     "1st order association"),
    ("6:2 FTS", None, "fluorotelomer sulfonate", "5.166", "0.113",
     "1st order association"),
    ("8:2 FTS", None, "fluorotelomer sulfonate", "4.634", "0.075",
     "1st order association"),
]

TK = ("yes - 26 HSA association constants measured in one laboratory "
      "by one method, the widest single-laboratory albumin set in this "
      "database. The measured values turn over above PFHpA while the "
      "paper's own docking columns rise monotonically, which bounds "
      "what structure-only prediction can do")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []
    for chem, ncf, head, ka, sd, model in S72:
        r = blank(fields)
        r.update(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            n_perfluorinated_carbons="" if ncf is None else str(ncf),
            head_group=head, protein=PROTEIN, protein_species="human",
            method=f"{METHOD}; fitted as a {model}",
            ka=ka, ka_units="log Ka, log M^-1",
            n_binding_sites="", protein_conc=PROT[0],
            protein_conc_units=PROT[1], ligand_conc_range=LIG[0],
            ligand_conc_units=LIG[1], is_computational="no",
            tk_relevance_noted=TK,
            source_table=(
                f"SI spreadsheet toxics-1134374-SI.xlsx, sheet "
                f"'S7.2 Ka_docking and exp', experimental column | "
                f"log Ka {ka} +/- {sd} | fitted as a {model} | the two "
                f"docking columns in the same sheet are computed and "
                f"are NOT ingested"))
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, existing = reader.fieldnames, list(reader)

    before = len(existing)
    # Keep the aggregate row only if it is NOT this paper's; the
    # abstract-range row for this doi is superseded.
    existing = [r for r in existing
                if not (r["study"] == STUDY and r["year"] == YEAR
                        and r["doi"] == DOI
                        and "serum albumin" in r["protein"].lower())]
    print(f"removed {before - len(existing)} existing HSA rows for {DOI}")

    new = rows(fields)
    print(f"adding  {len(new)} per-congener rows "
          f"-> {len(existing) + len(new)} total\n")
    cx = [r for r in new if r["head_group"] == "carboxylate"]
    print("   carboxylates, showing the turn:")
    for r in cx:
        print(f"      C{r['n_perfluorinated_carbons']:<3}"
              f"{r['chemical']:<10}log Ka {r['ka']}")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(existing + new)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    main()
