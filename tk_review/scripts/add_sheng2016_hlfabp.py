#!/usr/bin/env python3
"""Sheng et al. 2016: PFOA and PFNA on human L-FABP, two methods.

Sheng N, Li J, Liu H, Zhang A, Dai J.  "Interaction of perfluoroalkyl
acids with human liver fatty acid-binding protein."  Arch Toxicol
2016;90:217-227.  doi 10.1007/s00204-014-1391-7, PMID 25370009.
Read from the PDF; values are Tables 2 and 3.

This database previously held the paper as three abstract-only rows with
no numbers ("moderate affinity (numeric Kd not in abstract)").  The
tables are now in hand.

WHAT MAKES IT WORTH THE ROWS: the same laboratory measured the same two
compounds on the same protein by TWO methods, which is the within-
laboratory, between-method comparison this project otherwise has to
infer across papers.

    WT hL-FABP       fluorescence      ITC        difference
    PFOA             2.36 uM           6.49 uM    0.44 log
    PFNA             1.32 uM           3.14 uM    0.38 log

Both methods agree on the ORDER (PFNA binds tighter than PFOA) and on
the stoichiometry (1,8-ANS n = 1.01, ITC N = 0.999 for PFOA), and they
disagree on the LEVEL by a factor of about 2.5.

THE FOUR VARIANTS are engineered point mutants, not wild type, and are
ingested with that stated in the protein string so they can never be
read as a human measurement.  They carry the mechanism: N111D and R122G
lose binding entirely by fluorescence, and R122G's ITC stoichiometry
falls to about 0.55, which the paper reads as one of the two sites
being lost.

NOT RECOVERABLE FROM THE MAIN TEXT: the ITC protein and ligand
concentrations are in SI Table S1, and the PFHxS values are in SI Table
S3.  The main text says only that PFHxS binding was weak "with Kd
values ranging from 26.95 to 61.35 uM for the variants", which does not
make clear whether wild type is inside that range.  Recorded as stated,
without a number assigned to wild type.

Run:  python3 scripts/add_sheng2016_hlfabp.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Sheng et al.", "2016"
DOI, PMID = "10.1007/s00204-014-1391-7", "25370009"

WT = "liver fatty acid-binding protein (hL-FABP), recombinant wild type"
# The variant string must NOT contain the words "wild type", even in a
# negation: any downstream substring filter for wild type would then
# match the variants too. Same trap as "ratio" inside "titration".
VAR = ("liver fatty acid-binding protein (hL-FABP), "
       "ENGINEERED POINT VARIANT {v} (site-directed mutant, "
       "not the native protein)")

FLUOR = ("1,8-ANS fluorescence displacement, Horiba Fluoromax-4, buffer "
         "C (50 mM Tris-HCl pH 8.0), room temperature; apparent Kd from "
         "the IC50 of the displacement curve")
ITC = ("isothermal titration calorimetry (iTC200, MicroCal/GE), 25 C, "
       "buffer B, 20 injections (first 0.5 uL then 19 x 2 uL), 1000 rpm, "
       "120 s between injections, one-site binding model, triplicate")

FLUOR_PROT = ("0.25 uM hL-FABP (250 nM) in the displacement assay; "
              "1,8-ANS at 12.25 uM", "")
FLUOR_LIG = ("titrated against 0.25 uM hL-FABP; see IC50 in Table 2", "")
ITC_PROT = ("not stated in the main text -- SI Table S1", "")
ITC_LIG = ("not stated in the main text -- SI Table S1", "")

# Table 2: fluorescence displacement.  (variant, ANS Kd, ANS n,
#          PFOA IC50, PFOA Kd, PFNA IC50, PFNA Kd)  "" = no binding
T2 = [
    ("WT",    "3.56 +/- 0.48", "1.01 +/- 0.08",
     "8.14 +/- 1.17", "2.36 +/- 0.34", "4.55 +/- 0.31", "1.32 +/- 0.20"),
    ("S39G",  "6.64 +/- 0.96", "1.19 +/- 0.07",
     "7.16 +/- 0.46", "3.80 +/- 0.18", "4.06 +/- 0.53", "2.20 +/- 0.26"),
    ("M74G",  "4.79 +/- 0.79", "1.07 +/- 0.08",
     "6.68 +/- 0.92", "2.62 +/- 0.36", "4.93 +/- 0.45", "1.84 +/- 0.26"),
    ("N111D", "", "", "", "", "", ""),
    ("R122G", "", "", "", "", "", ""),
]

# Table 3: ITC at 25 C.  (compound, variant, N, Ka, Kd uM, dH, dS, dG)
# "ND" in the paper where the fit did not resolve.
T3 = [
    ("PFOA", "WT",    "0.999 +/- 0.016", "(1.54 +/- 0.14)e5", "6.49",
     "-29.49 +/- 0.68", "0.4", "-29.50 +/- 0.672"),
    ("PFOA", "S39G",  "1.08 +/- 0.024", "(1.21 +/- 0.13)e5", "8.26",
     "-24.47 +/- 0.72", "15", "-17.67 +/- 0.35"),
    ("PFOA", "M74G",  "0.984 +/- 0.013", "(1.12 +/- 0.09)e5", "8.93",
     "-5.57 +/- 0.11", "4.0", "-5.675 +/- 0.006"),
    ("PFOA", "N111D", "ND", "(3.52 +/- 1.89)e3", "284.09", "ND", "ND",
     "ND"),
    ("PFOA", "R122G", "0.554 +/- 0.042", "(1.19 +/- 0.32)e5", "8.40",
     "-17.30 +/- 2.71", "27", "-17.97 +/- 2.04"),
    ("PFNA", "WT",    "1.20 +/- 0.016", "(3.18 +/- 0.35)e5", "3.14",
     "-5.07 +/- 0.09", "8.2", "-5.275 +/- 0.113"),
    ("PFNA", "S39G",  "1.08 +/- 0.012", "(3.00 +/- 0.26)e5", "3.33",
     "-5.10 +/- 0.08", "7.9", "-5.301 +/- 0.122"),
    ("PFNA", "M74G",  "1.09 +/- 0.026", "(1.08 +/- 0.19)e5", "9.26",
     "-5.39 +/- 0.20", "5.0", "-5.513 +/- 0.073"),
    ("PFNA", "N111D", "ND", "(8.84 +/- 1.77)e3", "113.12", "ND", "ND",
     "ND"),
    ("PFNA", "R122G", "0.562 +/- 0.010", "(1.55 +/- 0.16)e5", "6.45",
     "-8.49 +/- 0.21", "-4.7", "-8.609 +/- 0.086"),
]

NCF = {"PFOA": "7", "PFNA": "8", "PFHxA": "5", "PFHxS": "6",
       "1,8-ANS": ""}
HEAD = {"PFOA": "carboxylate", "PFNA": "carboxylate",
        "PFHxA": "carboxylate", "PFHxS": "sulfonate",
        "1,8-ANS": "reference probe, not a PFAS"}

TK = ("yes - hL-FABP is the liver counterpart to albumin in the "
      "blood, and the paper resolves two binding sites with the "
      "residues that control them (Asn111 at the outer site, Arg122 "
      "at the interface, Ser39 in the inner site)")


def blank(fields):
    return dict.fromkeys(fields, "")


def num(x):
    return x.split("+/-")[0].strip() if x else ""


def rows(fields):
    out = []

    def base(chem):
        return dict(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            n_perfluorinated_carbons=NCF.get(chem, ""),
            head_group=HEAD.get(chem, ""), protein_species="human",
            buffer="buffer C, 50 mM Tris-HCl pH 8.0 (fluorescence)",
            is_computational="no", tk_relevance_noted=TK)

    # --- Table 2, fluorescence displacement -------------------------
    for (v, ans_kd, ans_n, pfoa_ic, pfoa_kd, pfna_ic, pfna_kd) in T2:
        prot = WT if v == "WT" else VAR.format(v=v)
        for chem, ic, kd in (("PFOA", pfoa_ic, pfoa_kd),
                             ("PFNA", pfna_ic, pfna_kd)):
            r = blank(fields)
            r.update(base(chem))
            r.update(
                protein=prot, method=FLUOR, temperature_c="",
                kd=num(kd), kd_units="uM" if kd else
                "no binding detected (NB in Table 2)",
                n_binding_sites=num(ans_n),
                binding_site_identified=(
                    "two sites resolved; Asn111 at the outer site and "
                    "Arg122 at the interface are required"),
                protein_conc=FLUOR_PROT[0], protein_conc_units="",
                ligand_conc_range=FLUOR_LIG[0], ligand_conc_units="",
                source_table=(
                    f"Table 2, {v} | " + (
                        f"Kd {kd} uM | IC50 {ic} uM | 1,8-ANS Kd "
                        f"{ans_kd} uM, n {ans_n}" if kd else
                        "NB, no binding detected")))
            out.append(r)
        # the reference probe itself
        if ans_kd:
            r = blank(fields)
            r.update(base("1,8-ANS"))
            r.update(
                protein=prot, method=FLUOR, temperature_c="",
                kd=num(ans_kd), kd_units="uM",
                n_binding_sites=num(ans_n),
                protein_conc=FLUOR_PROT[0], protein_conc_units="",
                ligand_conc_range="1,8-ANS titration series",
                ligand_conc_units="uM",
                tk_relevance_noted=(
                    "no - the fluorescent probe the PFAS displace, "
                    "recorded so the displacement Kd can be checked"),
                source_table=(f"Table 2, {v} | 1,8-ANS Kd {ans_kd} uM, "
                              f"n {ans_n}, R2 given in the table"))
            out.append(r)

    # --- Table 3, ITC -----------------------------------------------
    for (chem, v, n, ka, kd, dh, ds, dg) in T3:
        prot = WT if v == "WT" else VAR.format(v=v)
        r = blank(fields)
        r.update(base(chem))
        r.update(
            protein=prot, method=ITC, temperature_c="25",
            buffer="buffer B (ITC)",
            kd=kd, kd_units="uM",
            ka=ka, ka_units="M^-1, as reported in Table 3",
            n_binding_sites="" if n == "ND" else num(n),
            protein_conc=ITC_PROT[0], protein_conc_units="",
            ligand_conc_range=ITC_LIG[0], ligand_conc_units="",
            source_table=(
                f"Table 3, {v} | Kd {kd} uM | Ka {ka} M^-1 | N {n} | "
                f"dH {dh} kJ/mol | dS {ds} J/mol/deg | dG {dg} kJ/mol"))
        out.append(r)

    # --- the two compounds reported only in words -------------------
    for chem, note in (
        ("PFHxA",
         "no binding detected to WT hL-FABP or any variant, by either "
         "ITC or fluorescence displacement (Results text and Fig. S5)"),
        ("PFHxS",
         "weak binding; the Results text gives Kd 26.95 to 61.35 uM "
         "'for the variants' (SI Table S3) and does not make clear "
         "whether wild type falls inside that range, so no number is "
         "assigned here"),
    ):
        r = blank(fields)
        r.update(base(chem))
        r.update(protein=WT, method=f"{ITC}; and {FLUOR}",
                 temperature_c="25", kd="", kd_units=note,
                 protein_conc=ITC_PROT[0],
                 ligand_conc_range=ITC_LIG[0],
                 source_table="Results text; SI Table S3 not retrieved")
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
    existing = [r for r in existing
                if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"removed {before - len(existing)} existing {STUDY} {YEAR} rows")

    new = rows(fields)
    wt = sum(1 for r in new if "VARIANT" not in r["protein"])
    print(f"adding  {len(new)} rows ({wt} wild type, {len(new) - wt} "
          f"engineered variants) -> {len(existing) + len(new)} total")
    for r in new:
        if "VARIANT" not in r["protein"] and r["kd"]:
            kind = "ITC" if "calorimetry" in r["method"] else "fluor"
            print(f"   {r['chemical']:<9}{kind:<7}Kd {r['kd']:>7} uM  "
                  f"n {r['n_binding_sites'] or '-'}")

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
