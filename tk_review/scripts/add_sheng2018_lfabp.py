#!/usr/bin/env python3
"""Sheng 2018 from the primary: a rank order becomes eight dissociation constants.

Sheng N, Cui R, Wang J, Guo Y, Wang J, Dai J. "Cytotoxicity of novel
fluorinated alternatives to long-chain perfluoroalkyl substances to human
liver cell line and their binding capacity to human liver fatty acid binding
protein." Arch Toxicol 2018;92(1):359-369. doi:10.1007/s00204-017-2055-1
PMID 28864880.

WHAT WAS THERE BEFORE. One row, chemical "6:2 FTCA, 6:2 FTSA, HFPO-DA (GenX),
PFOA, PFOS, 6:2 Cl-PFESA (F-53B), HFPO-TA", kd "affinity rank order: 6:2 FTCA
< 6:2 FTSA < HFPO-DA < PFOA < PFOS = 6:2 Cl-PFESA = HFPO-TA (numeric values
not in abstract)", source_table "abstract". The working table dropped it under
"several compounds in one cell", correctly.

THE PRIMARY HAS TABLE 2, with IC50 and Kd for all eight compounds and the
probe. Note the rank order in the abstract omits HFPO-TeA, which Table 2 has
and which turns out to be the TIGHTEST binder of the eight.

Table 2 verbatim -- "R2, number of binding sites (N), dissociation constant
(Kd) of 1,8-ANS, IC50 value, and dissociation constant (Kd) of eight compounds
with hL-FABP determined fluorometrically by displacement of 1,8-ANS":

                  IC50 (uM)        Kd (uM)           n            R2
    1,8-ANS       --               7.47 +/- 0.97     0.93 +/- 0.03  0.996
    PFOA          2.15 +/- 0.54    8.03 +/- 2.10     --             0.999
    PFOS          1.34 +/- 0.21    4.99 +/- 0.97     --             0.997
    6:2 FTCA    116.88 +/- 0.52  436.55 +/- 1.92     --             0.999
    6:2 FTSA     78.97 +/- 0.56  345.54 +/- 2.18     --             0.996
    6:2 Cl-PFESA  1.14 +/- 0.12    4.05 +/- 0.42     --             0.999
    HFPO-DA       4.76 +/- 0.31   15.36 +/- 1.09     --             0.997
    HFPO-TA       1.21 +/- 0.37    4.36 +/- 1.17     --             0.998
    HFPO-TeA      0.98 +/- 0.18    3.35 +/- 0.46     --             0.998

INDEPENDENTLY CONFIRMED. Khazaee 2021's SI Table S2 quotes five of these from
Sheng 2018 as a literature comparison: PFOA 8.03 (IC50 2.15), PFOS 4.99
(1.34), 6:2 FTCA 436.55 (116.88), 6:2 FTSA 345.54 (78.97), 6:2 Cl-PFESA 4.05
(1.14). Every digit agrees with the table read here. That is a transcription
check against a copy made by people with no stake in this compilation.

N IS THE PROBE'S, NOT THE PFAA's -- the same trap as Sheng 2016. The 0.93 sits
in the 1,8-ANS row and the PFAS rows carry a dash. So these are
"n = 1 by single-site model", not a fitted PFAA stoichiometry.

CONDITIONS. "With 1,8-ANS as a fluorometric probe, fluorescence displacement
experiments were performed as described in earlier studies (Sheng et al.
2016a; Velkov et al. 2005), though with an adjustment in protein concentration
(200 nM)." So 0.2 uM hL-FABP, recombinant, expressed in E. coli BL21(DE3).
The paper does NOT restate the 1,8-ANS concentration or the PFAS series; the
binding parameters came from "titrating 200 nM hL-FABP with increasing
concentrations of ANS". Buffer and temperature are by reference to Sheng
2016a. Recorded as not stated rather than borrowed from the 2016 paper.

THE SAME LABORATORY DISAGREES WITH ITSELF. Sheng 2016 and Sheng 2018 both
measured PFOA against hL-FABP by 1,8-ANS displacement, two years apart:

    Sheng 2016   0.25 uM protein, ANS Kd 3.56 uM   ->  PFOA Kd 2.36 uM
    Sheng 2018   0.20 uM protein, ANS Kd 7.47 uM   ->  PFOA Kd 8.03 uM

a factor of 3.4 on PFOA, tracking a factor of 2.1 on the probe's own constant.
Both rows are carried, as the albumin sheet carries Hebert 2010's two protein
concentrations. It puts a floor under how much of the 1.2 log inter-laboratory
spread in this compilation is method rather than laboratory.

TWO COMPOUNDS ARE NEW TO THE L-FABP SHEET: 6:2 FTCA and 6:2 FTSA. The other
six already appear (PFOA, PFOS, 6:2 Cl-PFESA, HFPO-DA, HFPO-TA, HFPO-TeA).

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "protein_binding.csv")
DOI = "10.1007/s00204-017-2055-1"

PROTEIN = ("liver fatty acid-binding protein (hL-FABP), recombinant, "
           "expressed in E. coli BL21(DE3)")
METHOD = ("1,8-ANS fluorescence displacement, after Sheng et al. 2016 and "
          "Velkov et al. 2005 with protein adjusted to 200 nM; IC50 and Kd "
          "fitted in GraphPad Prism V5.0")
PROT_C = ("0.2 uM (200 nM) hL-FABP, stated as an adjustment to the Sheng 2016 "
          "protocol; the 1,8-ANS working concentration is not restated")
LIG_C = ("not stated; the paper gives IC50 values but not the PFAS "
         "concentration series behind them")
BUF = ("not stated; by reference to Sheng et al. 2016a, which used buffer C, "
       "50 mM Tris-HCl pH 8.0")

#  chemical, nC, head, IC50, Kd, n, R2
ROWS = [
    ("1,8-ANS (reference probe)", "", "fluorescent probe",
     "", "7.47", "0.93 +/- 0.03 (the PROBE's stoichiometry, not a PFAA's)",
     "0.996"),
    ("PFOA", "7", "carboxylate", "2.15", "8.03", "", "0.999"),
    ("PFOS", "8", "sulfonate", "1.34", "4.99", "", "0.997"),
    ("6:2 FTCA", "6", "carboxylate (fluorotelomer)",
     "116.88", "436.55", "", "0.999"),
    ("6:2 FTSA", "6", "sulfonate (fluorotelomer)",
     "78.97", "345.54", "", "0.996"),
    ("6:2 Cl-PFESA", "8", "chlorinated ether sulfonate",
     "1.14", "4.05", "", "0.999"),
    ("HFPO-DA", "5", "ether carboxylate", "4.76", "15.36", "", "0.997"),
    ("HFPO-TA", "8", "ether carboxylate", "1.21", "4.36", "", "0.998"),
    ("HFPO-TeA", "11", "ether carboxylate", "0.98", "3.35", "", "0.998"),
]


def build(fields):
    out = []
    for ch, nc, head, ic50, kd, n, r2 in ROWS:
        r = {f: "" for f in fields}
        src = f"Table 2 | Kd {kd} uM"
        if ic50:
            src += f" | IC50 {ic50} uM"
        src += f" | R2 {r2}"
        r.update({
            "study": "Sheng et al.", "year": "2018", "pmid": "28864880",
            "doi": DOI, "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "protein": PROTEIN,
            "protein_species": "human", "method": METHOD,
            "kd": kd, "kd_units": "uM", "n_binding_sites": n,
            "protein_conc": PROT_C, "protein_conc_units": "uM",
            "ligand_conc_range": LIG_C, "buffer": BUF,
            "is_computational": "no",
            "tk_relevance_noted": "yes - safety assessment of replacement "
                                  "chemicals",
            "source_table": src,
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
        print("no change needed (Sheng 2018 already ingested from the primary)")
        return 0

    merged, placed = [], False
    for r in rows:
        if (r.get("doi") or "") == DOI:
            if not placed:
                merged.extend(new); placed = True
            continue
        merged.append(r)
    if not placed:
        merged.extend(new)

    shutil.copy2(DB, DB + ".bak")
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(merged)
    print(f"Sheng 2018: replaced {len(old)} abstract-only rank-order row(s) "
          f"with {len(new)} rows from Table 2")
    print("   8 compounds with Kd, plus the 1,8-ANS probe")
    print("   new to the L-FABP sheet: 6:2 FTCA, 6:2 FTSA")
    print("   five values cross-checked against Khazaee 2021 SI Table S2: all agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
