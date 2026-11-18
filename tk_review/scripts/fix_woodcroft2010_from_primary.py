#!/usr/bin/env python3
"""Woodcroft 2010 re-ingested from the primary, replacing a secondary stub.

Woodcroft MW, Ellis DA, Rafferty SP, Burns DC, March RE, Stock NL,
Trumpour KS, Yee J, Munro K. "Experimental characterization of the mechanism
of perfluorocarboxylic acids' liver protein bioaccumulation: the key role of
the neutral species." Environ Toxicol Chem 2010;29(8):1669-77.
doi:10.1002/etc.199  PMID 20821618.

WHAT WAS THERE BEFORE. Two rows, both PFOA, both carrying several Kd values in
one cell ("13.1e-6, 23.8e-6" and "3.1e-6, 26.2e-6, 52.6e-6"), taken from a
SECONDARY source (OEHHA 2024 PFOA/PFOS PHG Table A6.2) with source_table
saying in as many words "primary paper not retrieved, no assay conditions
available". They were held out of the L-FABP sheet for exactly that reason.

WHAT THE PRIMARY ADDS.

  1. FOUR MORE COMPOUNDS. Table 1 reports ANS-displacement inhibition
     constants for five PFCAs, not one. The secondary table carried only PFOA.

  2. THE SITES ARE RESOLVED. Each composite cell is a multi-site fit, and the
     paper prints the sites separately, so each becomes its own row with
     binding_site_identified naming which site it is. The albumin sheet
     already works this way -- Moro 2022 sits in it as "ITC, sequential site 1
     of 2, 37 C" -- so the high-affinity site is the one that reaches the
     working table and the rest stay in the database as the record.

  3. THE SPECIES IS RAT, and it was already stored as rat, which the primary
     confirms: "the rat L-FABP coding sequence was PCR amplified from 2 ng of
     the pTRC99a vector", expressed as an N-terminal His6 fusion in E. coli
     BL21 (DE3) Rosetta. This matters because it makes Woodcroft the SECOND
     rat laboratory alongside Fan 2025, which is the minimum needed before a
     rat-human difference can be told apart from a laboratory effect.

  4. THE CONDITIONS, which the secondary had none of -- and one of them is a
     problem worth seeing (point 5).

  Fluorescence displacement (Table 1): 1 uM L-FABP pre-bound with 70 uM ANS in
  2.5 ml of Buffer C (20 mM potassium phosphate pH 7.4, 30 mM NaCl) at 20 C,
  titrated with competitor. PFNA and PFOA fit a TWO-site competitive model;
  PFHpA, PFHxA and PFPA fit a ONE-site model.

  ITC (Table 2): 0.162 mM L-FABP in the sample cell, 4.05 mM PFOA in the
  syringe, 29 injections of 10 ul, 30 C, fitted to a 3:1 sequential model.

  5. THE ITC BUFFER IS 1 M AMMONIUM SULFATE. The paper states the ITC was run
     in Buffer B -- "50 mM Tris-HCl pH 8.0, 0.5 mM EDTA, 1.0 mM DTT, 1 M
     (NH4)2SO4" -- while the fluorescence work was in Buffer C at 20 mM
     phosphate plus 30 mM NaCl. That is roughly a thirtyfold difference in
     ionic strength between the two methods in the same paper, and the two
     methods disagree by 4.2x on PFOA:

         fluorescence, Buffer C   Ki1 = 13.14 uM
         ITC, Buffer B (1 M salt) Kd1 =  3.13 uM

     The direction is what a salt effect would give for an anionic ligand, so
     the method difference and the ionic-strength difference are confounded
     and neither can be assigned the gap. Recorded in the buffer field of each
     row so it cannot be read off the sheet without seeing it. Gao 2025's salt
     series forced the same decision on the albumin side.

VALUES, Table 1 (Ki, uM, +/- SE) and Table 2 (Kd, uM, +/- SE), transcribed
from the primary and cross-checked against the Results text, which quotes
"Kd1 = 3.1 +/- 0.3 uM", "Kd2 = 26 +/- 2 uM" and "Kd3 = 53 +/- 2 uM" -- the
rounded forms of Table 2's 3.13, 26.25 and 52.56. The two previously stored
composite cells agree with the primary to the digits they carried.

NAMING. The paper writes PFPA for perfluoropentanoic acid. This database
calls that compound PFPeA (Zhang 2013's rows use PFPeA), so it is stored as
PFPeA with the paper's own abbreviation kept in source_table.

ANS itself is stored too, as the reference ligand whose constants the
competitive fit depends on. It is not an analyte and the working table drops
it as a probe.

Octanoic acid was titrated (Fig. 3) but the paper states that its solubility
limit meant "it was not possible to obtain accurate binding information by
nonlinear regression", so no constant is stored for it.

Idempotent: the Woodcroft rows are deleted and rewritten from the table below,
so re-running changes nothing.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "protein_binding.csv")
DOI = "10.1002/etc.199"

BUF_FLUOR = ("Buffer C, 20 mM potassium phosphate pH 7.4, 30 mM NaCl "
             "(fluorescence displacement)")
BUF_ITC = ("Buffer B, 50 mM Tris-HCl pH 8.0, 0.5 mM EDTA, 1.0 mM DTT, "
           "1 M (NH4)2SO4 -- NOTE the ITC was run at ~30x the ionic strength "
           "of the fluorescence work, so the 4.2x PFOA gap between the two "
           "methods is confounded with a salt effect")

PROT_FLUOR = "1 uM L-FABP pre-bound with 70 uM ANS in a final volume of 2.5 ml"
LIG_FLUOR = ("competitor titration; Fig. 3 plots log[competitor] from -6 to -2 "
             "(1 uM to 10 mM) but the per-compound endpoints are not tabulated")
PROT_ITC = "0.162 mM (162 uM) L-FABP in the ITC sample cell"
LIG_ITC = "4.05 mM PFOA in the injection syringe, 29 injections of 10 ul"

M_FLUOR = ("competitive displacement of 1-anilino-8-naphthalenesulfonic acid "
           "(ANS) from recombinant rat His6-L-FABP, fluorimetric, nonlinear "
           "regression to a one- or two-site competitive binding model")
M_ITC = ("isothermal titration calorimetry (VP-ITC), fitted to a 3:1 "
         "(protein:ligand) sequential binding model")
M_ANS = ("direct fluorimetric titration of L-FABP with ANS, two-site binding "
         "model, Scatchard analysis for stoichiometry")

SRC = "primary: Woodcroft et al. 2010, Environ Toxicol Chem 29:1669-1677, "

#  chemical, nC, head, kd_uM, n_sites, site label, method, prot, lig, T, buffer, table
ROWS = [
    # --- Table 1: ANS displacement, 20 C, Buffer C ---
    ("PFNA", "8", "carboxylate", "6.07", "2 (two-site competitive fit)",
     "site 1 of 2 (high-affinity site)", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki1 = 6.07 +/- 2.70 uM"),
    ("PFNA", "8", "carboxylate", "9.84", "2 (two-site competitive fit)",
     "site 2 of 2 (low-affinity site)", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki2 = 9.84 +/- 3.95 uM"),
    ("PFOA", "7", "carboxylate", "13.14", "2 (two-site competitive fit)",
     "site 1 of 2 (high-affinity site)", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki1 = 13.14 +/- 0.27 uM"),
    ("PFOA", "7", "carboxylate", "23.81", "2 (two-site competitive fit)",
     "site 2 of 2 (low-affinity site)", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki2 = 23.81 +/- 0.89 uM"),
    ("PFHpA", "6", "carboxylate", "61.17", "1 (one-site competitive fit)",
     "", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki1 = 61.17 +/- 1.95 uM, one-site model"),
    ("PFHxA", "5", "carboxylate", "261.70", "1 (one-site competitive fit)",
     "", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki1 = 261.70 +/- 7.00 uM, one-site model"),
    ("PFPeA", "4", "carboxylate", "608.40", "1 (one-site competitive fit)",
     "", M_FLUOR, PROT_FLUOR, LIG_FLUOR, "20",
     BUF_FLUOR, SRC + "Table 1, Ki1 = 608.40 +/- 16.06 uM, one-site model; "
     "the paper abbreviates perfluoropentanoic acid as PFPA"),
    # --- Table 2: ITC, 30 C, Buffer B ---
    ("PFOA", "7", "carboxylate", "3.13", "3 (3:1 sequential ITC fit)",
     "site 1 of 3 (high-affinity site)", M_ITC, PROT_ITC, LIG_ITC, "30",
     BUF_ITC, SRC + "Table 2, Kd1 = 3.13 +/- 0.26 uM"),
    ("PFOA", "7", "carboxylate", "26.25", "3 (3:1 sequential ITC fit)",
     "site 2 of 3 (low-affinity site)", M_ITC, PROT_ITC, LIG_ITC, "30",
     BUF_ITC, SRC + "Table 2, Kd2 = 26.25 +/- 2.08 uM"),
    ("PFOA", "7", "carboxylate", "52.56", "3 (3:1 sequential ITC fit)",
     "site 3 of 3 (low-affinity site)", M_ITC, PROT_ITC, LIG_ITC, "30",
     BUF_ITC, SRC + "Table 2, Kd3 = 52.56 +/- 2.31 uM"),
    # --- the reference ligand, for the record ---
    ("ANS (1-anilino-8-naphthalenesulfonic acid, reference ligand)", "", "",
     "1.0", "2 (two-site, 2:1 ANS:protein by Scatchard)",
     "site 1 of 2 (high-affinity site)", M_ANS,
     "1 uM L-FABP in 2.5 ml of Buffer C",
     "ANS titration series, 1-105 uM", "20", BUF_FLUOR,
     SRC + "Results text, Kd1 = 1.0 +/- 0.2 uM for the high-affinity site"),
    ("ANS (1-anilino-8-naphthalenesulfonic acid, reference ligand)", "", "",
     "10.5", "2 (two-site, 2:1 ANS:protein by Scatchard)",
     "site 2 of 2 (low-affinity site)", M_ANS,
     "1 uM L-FABP in 2.5 ml of Buffer C",
     "ANS titration series, 1-105 uM", "20", BUF_FLUOR,
     SRC + "Results text, Kd2 = 10.5 +/- 0.9 uM for the low-affinity site"),
]

PROTEIN = ("liver fatty acid-binding protein (rL-FABP), recombinant "
           "N-terminal His6 fusion, expressed in E. coli BL21(DE3) Rosetta, "
           "delipidated on Phenyl Sepharose")


def build(fields):
    out = []
    for (ch, nc, head, kd, nsites, site, meth, prot, lig, T, buf, tbl) in ROWS:
        r = {f: "" for f in fields}
        r.update({
            "study": "Woodcroft et al.", "year": "2010", "pmid": "20821618",
            "doi": DOI, "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "protein": PROTEIN, "protein_species": "rat",
            "method": meth, "kd": kd, "kd_units": "uM",
            "n_binding_sites": nsites, "binding_site_identified": site,
            "protein_conc": prot, "protein_conc_units": "uM",
            "ligand_conc_range": lig, "temperature_c": T, "buffer": buf,
            "is_computational": "no", "source_table": tbl,
        })
        out.append(r)
    return out


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        fields, rows = rdr.fieldnames, list(rdr)

    old = [r for r in rows if (r.get("doi") or "") == DOI]
    new = build(fields)
    if old and len(old) == len(new) and all(
            a.get(k) == b.get(k) for a, b in zip(old, new) for k in fields):
        print("no change needed (Woodcroft 2010 already ingested from the primary)")
        return 0

    merged, placed = [], False
    for r in rows:
        if (r.get("doi") or "") == DOI:
            if not placed:
                merged.extend(new); placed = True
            continue          # drop every old Woodcroft row
        merged.append(r)
    if not placed:
        merged.extend(new)

    shutil.copy2(DB, DB + ".bak")
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(merged)

    print(f"Woodcroft 2010: replaced {len(old)} secondary-source row(s) "
          f"with {len(new)} rows from the primary")
    print("   5 PFCAs by ANS displacement (PFNA, PFOA, PFHpA, PFHxA, PFPeA), "
          "sites resolved")
    print("   PFOA by ITC, 3 sequential sites resolved")
    print("   ANS reference ligand, 2 sites")
    print("   conditions recorded; the ITC buffer carries 1 M (NH4)2SO4")
    return 0


if __name__ == "__main__":
    sys.exit(main())
