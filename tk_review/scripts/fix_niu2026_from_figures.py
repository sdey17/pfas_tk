#!/usr/bin/env python3
"""Niu 2026: the per-congener Km/Vmax from the figures, and a swap in the text.

Niu S, Ma Y, Tewari A, Ng C. "Transport of perfluoroalkyl substances (PFAS)
by three renal transporters: implications for PFAS bioaccumulation
mechanisms." Environ Health (Wash) 2026;4(9):2150-2157.
doi:10.1021/envhealth.6c00088  PMID 42775203, PMC13595387.

The record held eight rows with no constants, because the Km and Vmax values
are printed inside Figures 2f and 3c,f rather than in a table. The reviewer
supplied the full text and the figures are legible.

    OAT1  (HEK293-transfected cells, net of wild-type)
        PFOA     Km 17.0 uM   Vmax 523 pmol/mg protein/min   r2 0.90
        PFOS     Km 12.6      Vmax 674                        r2 0.84
        HFPO-TA  Km 30.7      Vmax 428                        r2 0.77
    P-gp  (inside-out membrane vesicles, net of control vesicles)
        HFPO-TA  Km 2.78      Vmax 631                        r2 0.69
        PFOS     Km 2.74      Vmax 382                        r2 0.72
        PFOA     Km 1.41      Vmax 299                        r2 0.82
    BCRP  (inside-out membrane vesicles, net of control vesicles)
        HFPO-TA  Km 2.67      Vmax 621                        r2 0.75
        PFOA     Km 1.29      Vmax 360                        r2 0.80
        PFOS     Km 7.92      Vmax 192                        r2 0.73

THE PAPER'S TEXT AND ABSTRACT SWAP P-gp AND BCRP IN THE Vmax/Km PARAGRAPH.
The text reads "The highest ratio of Vmax/Km of P-gp was observed for PFOA
(279 ...), followed by HFPO-TA (234) and PFOS (24)", and for BCRP
"HFPO-TA (227) > PFOA (212) > PFOS (139)". Dividing the figures' own fitted
constants gives the opposite assignment, on all six values:

    from Fig 3c (P-gp)   HFPO-TA 227.0   PFOA 212.1   PFOS 139.4
    from Fig 3f (BCRP)   PFOA    279.1   HFPO-TA 232.6  PFOS 24.2

Each set matches the text's figures for the OTHER transporter, to the
rounding. Three independent checks say the FIGURES are right and the text
paragraph is transposed:

  1. The abstract's own Km ranges track the figures correctly -- "Km ranging
     from 1.29 to 7.92 umol/L and 1.41 to 2.78" for BCRP and P-gp, which are
     exactly the Km sets in Figures 3f and 3c respectively.
  2. The same arithmetic on OAT1 (Fig 2f) reproduces the text exactly
     (523/17.0 = 30.8, 674/12.6 = 53.5, 428/30.7 = 13.9), so the method is
     right and only the P-gp/BCRP paragraph is wrong.
  3. The text also cites the wrong panels, "Figure 3b,d" for constants that
     are in 3c,f; b and d are the time-course and control-vesicle panels.

Consequence: the abstract's headline, "The highest Vmax/Km ratio was observed
for PFOA transport by P-gp (279 pmol/mg protein/min/(umol/L))", should read
BCRP. This database stores the figure values and the ratios computed from
them, and records the discrepancy rather than the abstract.

THE AUTHORS DOWNGRADE THEIR OWN CONSTANTS: "The fitted kinetic parameters are
intended to describe the overall uptake trend rather than exact values." Four
of the nine fits have r2 below 0.80. Carried as a quality flag.

OAT1 IS NOT COMPARABLE WITH P-gp AND BCRP HERE, and the paper says so: "the
specific rates of PFAS uptake cannot be directly compared between these
transporters" because the concentrations differed, and OAT1 was measured in
whole cells while P-gp and BCRP used vesicles whose background uptake was
1000-2000 pmol/mg against 100-400 in the cells.

A SCALING NUMBER WORTH KEEPING FOR STAGE 3. The paper measures transporter
protein in human kidney cortex: OAT1 5.3 +/- 1.9 pmol/mg protein, P-gp
2.1 +/- 0.8, BCRP below detection (n = 20). That is the expression side of
an IVIVE scaling factor, measured in the tissue rather than assumed.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "transporter_kinetics.csv")
DOI = "10.1021/envhealth.6c00088"
STUDY = "Niu S, Ma Y, Tewari A, Ng C"
SWAP = (
    "THE PAPER'S TEXT AND ABSTRACT SWAP P-gp AND BCRP in the Vmax/Km "
    "paragraph. Dividing the figures' own constants gives {here} uL/min/mg "
    "for this row, which matches what the text attributes to the OTHER "
    "transporter. The abstract's Km ranges (BCRP 1.29-7.92, P-gp 1.41-2.78) "
    "track the figures correctly, and the identical arithmetic on OAT1 "
    "reproduces the text exactly, so the figures are right and the "
    "P-gp/BCRP paragraph is transposed. The abstract's headline claim that "
    "the highest Vmax/Km is PFOA by P-gp should read BCRP")
TREND = ("The authors state the fitted constants are \"intended to describe "
         "the overall uptake trend rather than exact values\" (r2 {r2})")
NOCOMP = ("Measured in {sys}. The paper states that specific uptake rates "
          "cannot be compared between OAT1 and the two efflux transporters, "
          "because the concentrations differed and the vesicle background "
          "(1000-2000 pmol/mg) far exceeds the cell background (100-400)")
EXPR = ("Transporter protein in human kidney cortex, measured in this paper "
        "(n = 20): OAT1 5.3 +/- 1.9 pmol/mg protein, P-gp 2.1 +/- 0.8, BCRP "
        "below detection -- the expression side of an IVIVE scaling factor")
CELLS = "human OAT1-transfected HEK293 cells, net of wild-type"
VES = "inside-out human membrane vesicles, net of control vesicles"

#  transporter, chemical, nC, head, km, vmax, r2, figure
ROWS = [
    ("OAT1", "PFOA", "7", "carboxylate", "17.0", "523", "0.90", "Fig 2f"),
    ("OAT1", "PFOS", "8", "sulfonate", "12.6", "674", "0.84", "Fig 2f"),
    ("OAT1", "HFPO-TA", "8", "ether carboxylate", "30.7", "428", "0.77",
     "Fig 2f"),
    ("P-gp", "HFPO-TA", "8", "ether carboxylate", "2.78", "631", "0.69",
     "Fig 3c"),
    ("P-gp", "PFOS", "8", "sulfonate", "2.74", "382", "0.72", "Fig 3c"),
    ("P-gp", "PFOA", "7", "carboxylate", "1.41", "299", "0.82", "Fig 3c"),
    ("BCRP", "HFPO-TA", "8", "ether carboxylate", "2.67", "621", "0.75",
     "Fig 3f"),
    ("BCRP", "PFOA", "7", "carboxylate", "1.29", "360", "0.80", "Fig 3f"),
    ("BCRP", "PFOS", "8", "sulfonate", "7.92", "192", "0.73", "Fig 3f"),
]
#  transporter, chemical, nC, head, why no constant
NEG = [
    ("OAT1", "HFPO-DA", "5", "ether carboxylate",
     "not transported by OAT1 (and neither P-gp nor BCRP showed transport "
     "capability for HFPO-DA, data not shown)"),
    ("OAT1", "9Cl-PF3ONS", "8", "ether sulfonate",
     "transported, but excluded from the concentration-dependent assay for "
     "low solubility, so no constant. Equilibrium uptake 896 pmol/mg protein "
     "by OAT1 against 161 by P-gp and 137 by BCRP -- the paper reads this as "
     "a high bioaccumulation potential, consistent with a reported half-life "
     "of 18.5 years"),
    ("OAT1", "11Cl-PF3OUdS", "10", "ether sulfonate",
     "transported, but excluded from the concentration-dependent assay for "
     "low solubility. Uptake order here is P-gp > BCRP > OAT1, the reverse "
     "of 9Cl-PF3ONS"),
]
GENE = {"OAT1": "SLC22A6", "P-gp": "ABCB1", "BCRP": "ABCG2"}
DIRN = {"OAT1": "secretion (basolateral uptake into the proximal tubule cell)",
        "P-gp": "secretion (apical efflux into the tubular lumen)",
        "BCRP": "secretion (apical efflux into the tubular lumen)"}


def build(fields):
    out = []
    for tr, ch, nc, head, km, vm, r2, fig in ROWS:
        cl = round(float(vm) / float(km), 1)
        note = TREND.format(r2=r2) + ". " + NOCOMP.format(
            sys=CELLS if tr == "OAT1" else VES) + ". " + EXPR
        if tr in ("P-gp", "BCRP"):
            note = SWAP.format(here=cl) + ". " + note
        r = {f: "" for f in fields}
        r.update({
            "study": STUDY, "year": "2026", "pmid": "42775203", "doi": DOI,
            "transporter": tr, "gene": GENE[tr], "species_orthologue": "human",
            "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head,
            "expression_system": CELLS if tr == "OAT1" else VES,
            "measurement_type": "Km/Vmax (Michaelis-Menten)",
            "km": km, "km_units": "uM",
            "vmax": vm, "vmax_units": "pmol/mg protein/min",
            "substrate_or_inhibitor": "substrate",
            "direction_reabsorption_or_secretion": DIRN[tr],
            "sex_or_hormone_effect": "not assessed",
            "assay_volume_note": note,
            "source_table": (f"{fig} (constants printed on the panel) | "
                             f"Vmax/Km = {cl} uL/min/mg protein, computed "
                             f"here from the figure"),
        })
        out.append(r)
    for tr, ch, nc, head, why in NEG:
        r = {f: "" for f in fields}
        r.update({
            "study": STUDY, "year": "2026", "pmid": "42775203", "doi": DOI,
            "transporter": tr, "gene": GENE[tr], "species_orthologue": "human",
            "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "expression_system": CELLS,
            "measurement_type": "uptake measured; no Michaelis-Menten fit",
            "substrate_or_inhibitor": (
                "not a substrate" if "not transported" in why else "substrate"),
            "direction_reabsorption_or_secretion": DIRN[tr],
            "sex_or_hormone_effect": "not assessed",
            "assay_volume_note": why + ". " + EXPR,
            "source_table": "Results and Conclusions (full text)",
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
        print("no change needed (Niu 2026 figures already ingested)")
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
    print(f"Niu 2026: {len(old)} row(s) -> {len(new)}")
    print("   9 rows gain Km + Vmax read from Figures 2f and 3c,f")
    print("   P-gp/BCRP text-vs-figure swap recorded on all 6 efflux rows")
    print("   3 further rows: HFPO-DA not transported, two Cl-ethers without "
          "constants")
    return 0


if __name__ == "__main__":
    sys.exit(main())
