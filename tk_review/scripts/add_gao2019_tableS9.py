#!/usr/bin/env python3
"""Gao et al. 2019 SI Table S9, read at source: three protein matrices.

Gao K, Zhuang T, Liu X, Fu J, Zhang J, Fu J, Wang L, Zhang A, Liang Y,
Song M, Jiang G.  "Prenatal Exposure to Per- and Polyfluoroalkyl
Substances (PFASs) and Association between the Placental Transfer
Efficiencies and Dissociation Constant of Serum Proteins-PFAS
Complexes."  Environ Sci Technol 2019;53(11):6529-6538.
doi 10.1021/acs.est.9b00715, PMID 31099564.

This replaces the 26 rows transcribed from Fan et al. 2025 Table 1.
Fan's two columns reproduce Table S9 EXACTLY on all 26 values, so the
secondary transcription was sound. The primary adds three things.

A THIRD MATRIX. Table S9 has Kd-SP as well as Kd-HP and Kd-LP: the
dissociation constant against POOLED HUMAN SERUM PROTEINS, 14 values
that no secondary source carried. It sits between the two purified
proteins on almost every compound, which is what a mixture of them
should do, and it is the matrix a toxicokinetic model actually faces.

A COMPOUND Fan dropped: 6:2 Cl-PFESA, on all three matrices.

TWO NAMING CORRECTIONS, both inherited from Fan:

  Fan labels the C13 row PFTeDA. Gao calls it PFTrA, perfluoro-
  TRIdecanoic acid, 12 perfluorinated carbons. PFTeDA is the C14. The
  rows are renamed and n_perfluorinated_carbons corrected from 13 to
  12. Zhang et al. 2013 separately reports a real PFTeDA, so the two
  were being merged under one label.

  Fan labels Gao's PFUdA as PFUnDA. Those are the same compound and
  the project's canonical name is PFUnDA, so that one stands.

AND AN ERROR IN GAO'S OWN TABLE, recorded but not propagated: its
"Zhang et al. 2013" column puts 1034 +/- 110 on the PFHxA row. Zhang's
own Table 1 gives PFHxA as ND and 1034 +/- 110 as PFBS. Nine other
values in that column match Zhang exactly, so it is a single
misplacement. This database holds Zhang as primary and is unaffected.

CONDITIONS, from SI page S5 and the main text: protein at 0.006 mM
(6 uM) for all three matrices, PFAS at 0.002 to 0.012 mM (2, 4, 6, 8,
10, 12 uM), incubated 1 h at 37 C, then a rapid equilibrium dialysis
device with an 8 kDa cutoff for 4 h. The ligand:protein ratio therefore
reaches 2 at the top of the series.

Run:  python3 scripts/add_gao2019_tableS9.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Gao et al.", "2019"
DOI, PMID = "10.1021/acs.est.9b00715", "31099564"

METHOD = ("equilibrium dialysis: protein or serum spiked with PFAS, "
          "1 h at 37 C in a circulating water bath, then 200 uL into a "
          "rapid equilibrium dialysis device with an 8 kDa molecular "
          "weight cutoff for 4 h; protein precipitated with 90/10 "
          "acetonitrile/water + 0.1% formic acid, then LC-MS/MS. "
          "Incubation after Beesoon and Martin")
PROT = ("0.006 mM (6 uM)", "mM")
LIG = ("0.002, 0.004, 0.006, 0.008, 0.010 and 0.012 mM "
       "(2 to 12 uM, a six-point series)", "mM")
RATIO = "2"      # 12 uM PFAS over 6 uM protein at the top of the series

MATRIX = {
    "HP": ("human serum albumin (HSA)", "human", "Kd-HP"),
    "SP": ("pooled human serum proteins (whole serum)", "human",
           "Kd-SP"),
    "LP": ("liver fatty acid-binding protein (L-FABP)", "human",
           "Kd-LP"),
}

# SI Table S9. (chemical, nCF or None, head, Kd-HP, Kd-SP, Kd-LP)
S9 = [
    ("PFBA",   3,  "carboxylate", "419 +/- 37", "482 +/- 45", "879 +/- 51"),
    ("PFPeA",  4,  "carboxylate", "286 +/- 17", "387 +/- 11", "685 +/- 23"),
    ("PFHxA",  5,  "carboxylate", "195 +/- 14", "243 +/- 22", "443 +/- 25"),
    ("PFHpA",  6,  "carboxylate", "172 +/- 23", "209 +/- 17", "358 +/- 16"),
    ("PFOA",   7,  "carboxylate", "115 +/- 16", "166 +/- 10", "197 +/- 13"),
    ("PFNA",   8,  "carboxylate", "73 +/- 21",  "84 +/- 17",  "90 +/- 8"),
    ("PFDA",   9,  "carboxylate", "65 +/- 8",   "70 +/- 13",  "77 +/- 14"),
    ("PFUnDA", 10, "carboxylate", "79 +/- 13",  "75 +/- 10",  "54 +/- 5"),
    ("PFDoDA", 11, "carboxylate", "103 +/- 12", "100 +/- 13", "95 +/- 11"),
    ("PFTrDA", 12, "carboxylate", "265 +/- 21", "284 +/- 27", "317 +/- 20"),
    ("PFBS",   4,  "sulfonate",   "169 +/- 15", "239 +/- 18", "436 +/- 33"),
    ("PFHxS",  6,  "sulfonate",   "79 +/- 12",  "158 +/- 23", "210 +/- 14"),
    ("PFOS",   8,  "sulfonate",   "38 +/- 5",   "49 +/- 8",   "81 +/- 7"),
    ("6:2 Cl-PFESA", None, "ether sulfonate",
     "67 +/- 6", "82 +/- 5", "102 +/- 10"),
]

# Fan et al. 2025 Table 1, the secondary transcription this replaces.
# Asserted before anything is written.
FAN = {
    "PFBA": (419, 879), "PFPeA": (286, 685), "PFHxA": (195, 443),
    "PFHpA": (172, 358), "PFOA": (115, 197), "PFNA": (73, 90),
    "PFDA": (65, 77), "PFUnDA": (79, 54), "PFDoDA": (103, 95),
    "PFTrDA": (265, 317), "PFBS": (169, 436), "PFHxS": (79, 210),
    "PFOS": (38, 81),
}

TK = ("yes - Gao measured HSA, whole serum and L-FABP side by side and "
      "related the Kd to placental transfer efficiency, which is the "
      "same liver-versus-blood competition Fan et al. 2025 later "
      "measured in rats")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []
    for chem, ncf, head, hp, sp, lp in S9:
        for key, val in (("HP", hp), ("SP", sp), ("LP", lp)):
            protein, species, label = MATRIX[key]
            r = blank(fields)
            r.update(
                study=STUDY, year=YEAR, pmid=PMID, doi=DOI,
                chemical=chem,
                n_perfluorinated_carbons="" if ncf is None else str(ncf),
                head_group=head, protein=protein,
                protein_species=species, method=METHOD,
                kd=val.split("+/-")[0].strip(), kd_units="uM",
                n_binding_sites="",
                protein_conc=PROT[0], protein_conc_units=PROT[1],
                ligand_conc_range=LIG[0], ligand_conc_units=LIG[1],
                ligand_over_protein_max=RATIO,
                temperature_c="37", buffer="0.1 M PBS, pH 7.4",
                is_computational="no", tk_relevance_noted=TK,
                source_table=(f"SI Table S9, {label} = {val} uM | read "
                              f"from the Supporting Information, primary"))
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    prim = {c: (float(hp.split("+/-")[0]), float(lp.split("+/-")[0]))
            for c, _, _, hp, _, lp in S9}
    for chem, (fan_hp, fan_lp) in FAN.items():
        got = prim.get(chem)
        assert got and abs(got[0] - fan_hp) < 1e-9 \
            and abs(got[1] - fan_lp) < 1e-9, (
            f"{chem}: primary {got}, Fan said ({fan_hp}, {fan_lp}) -- "
            f"the two readings disagree, not proceeding")
    print(f"cross-check passed: the primary reproduces all "
          f"{2 * len(FAN)} values transcribed from Fan's Table 1")

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, existing = reader.fieldnames, list(reader)

    before = len(existing)
    existing = [r for r in existing
                if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"removed {before - len(existing)} secondary {STUDY} {YEAR} rows")

    new = rows(fields)
    print(f"adding  {len(new)} primary rows "
          f"-> {len(existing) + len(new)} total")
    for key in ("HP", "SP", "LP"):
        got = [r for r in new if MATRIX[key][0] == r["protein"]]
        print(f"   {MATRIX[key][2]:<6} {MATRIX[key][0][:44]:<46}"
              f"{len(got):>3} rows")

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
