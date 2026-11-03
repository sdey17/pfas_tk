#!/usr/bin/env python3
"""Fan et al. 2025: 14 PFAS against rat L-FABP and rat serum albumin.

Fan X, Li X, Li T, Shao B, Niu S, Fan W, Wang Q, Zhong L, Wang X, Ma F,
Zheng G, Peng H, Chen L, Dong Z.  "Binding of per- and polyfluoroalkyl
substances with liver and serum proteins in rats: implications for
physiologically based pharmacokinetic modelling."  Environ Int
2025;201:109591.  doi 10.1016/j.envint.2025.109591, PMID 40494107.
Open access, CC BY-NC.  Read from the PDF; values are Table 1.

Why this paper matters to the review, beyond adding rows:

  * It is the ONLY source that measures serum albumin and L-FABP for the
    SAME species in ONE laboratory.  Every cross-protein comparison in
    this project up to now has had to cross a laboratory boundary too.
  * It reports the STOICHIOMETRY.  n is the parameter the Kd-to-
    partition bridge needs and almost nobody states; every albumin
    conversion here runs on an assumed n = 1.  Fan gives n for RSA by
    two methods.
  * Its Table 1 is a cross-species compilation that carries the
    per-congener hL-FABP values of Zhang et al. 2013, which this
    database has held as "abstract only, per-congener constants NOT
    retrieved" -- and which Han et al. 2021 quoted without this project
    being able to check them.

THREE SETS ARE INGESTED, all Fan's own measurements:

  rL-FABP   1,8-ANS fluorescence displacement, single binding site.
  RSA UF    ultrafiltration, n and Kd.
  RSA ITC   isothermal titration calorimetry, n and Kd.

The paper states ITC is the set to prefer: 20-95% of PFAS adsorbs to the
ultrafiltration membrane at low concentration, correlated with log Kow
(p < 0.001), so the ultrafiltration affinities are biased high.  Both are
stored, with that caution on the ultrafiltration rows.

A TRANSCRIPTION CAUTION, recorded because it will otherwise trip up a
later reader.  The Results text quotes "PFOA (Kd = 26.77 uM, n = 3.11)"
and "GenX (Kd = 17.87 uM, n = 2.47)", each pairing an ULTRAFILTRATION Kd
with an ITC n.  Table 1 keeps them apart: PFOA is n 4.24 / Kd 26.77 by
ultrafiltration and n 3.11 / Kd 94.54 by ITC.  The table is internally
consistent with the ranges the text gives for each method (ITC n ~1-4,
ultrafiltration n ~3-5), so the table is followed and the text's paired
quotes are not.

Run:  python3 scripts/add_fan2025_rat_proteins.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Fan et al.", "2025"
DOI, PMID = "10.1016/j.envint.2025.109591", "40494107"

# Table 1.  (chemical, nCF, head group,
#            rL-FABP Kd uM, +/-,
#            RSA ultrafiltration n, +/-, Kd uM, +/-,
#            RSA ITC n, +/-, Kd uM, +/-)
# "" means the table leaves the cell empty; PFBA's rL-FABP cell is "---",
# which the Results text explains: binding was undetectable.
T1 = [
    ("PFBA",    3,  "carboxylate", "", "",
     "4.46", "0.21", "16.21", "2.95",   "2.84", "0.35", "168.40", "33.08"),
    ("PFPeA",   4,  "carboxylate", "27.57", "3.68",
     "3.92", "0.16", "13.9", "1.86",    "3.33", "0.21", "303.08", "111.38"),
    ("PFHxA",   5,  "carboxylate", "132.16", "6.16",
     "3.36", "0.14", "31.36", "2.86",   "1.93", "0.65", "166.70", "33.24"),
    ("PFHpA",   6,  "carboxylate", "44.54", "1.12",
     "3.02", "0.62", "10.15", "1.37",   "2.79", "0.15", "79.96", "8.83"),
    ("PFOA",    7,  "carboxylate", "12.30", "3.60",
     "4.24", "0.52", "26.77", "2.77",   "3.11", "0.33", "94.54", "8.16"),
    ("PFNA",    8,  "carboxylate", "6.55", "0.19",
     "4.21", "0.49", "13.27", "1.20",   "1.06", "0.05", "31.18", "2.45"),
    ("PFDA",    9,  "carboxylate", "5.36", "0.73",
     "3.70", "0.32", "10.33", "0.85",   "1.91", "0.49", "161.20", "9.41"),
    ("PFUnDA",  10, "carboxylate", "6.95", "0.19",
     "", "", "", "",                    "1.42", "0.3", "59.18", "9.79"),
    ("PFDoDA",  11, "carboxylate", "13.46", "4.18",
     "", "", "", "",                    "1.17", "0.95", "32.78", "21.41"),
    ("PFTeDA",  13, "carboxylate", "6.32", "1.09",
     "", "", "", "",                    "1.02", "0.09", "195.30", "8.54"),
    ("PFBS",    4,  "sulfonate",   "110.82", "4.80",
     "3.68", "0.14", "2.78", "0.58",    "1.59", "0.59", "176.30", "39.11"),
    ("PFHxS",   6,  "sulfonate",   "13.02", "0.11",
     "4.24", "0.53", "78.87", "2.74",   "2.15", "0.18", "26.87", "1.18"),
    ("PFOS",    8,  "sulfonate",   "20.73", "0.90",
     "4.61", "0.32", "13.77", "2.65",   "2.79", "0.40", "129.80", "17.54"),
    ("HFPO-DA (GenX)", None, "ether carboxylate (polyfluoroether)",
     "108.23", "6.47",
     "3.26", "0.23", "17.87", "0.52",   "2.47", "0.4", "67.57", "20.34"),
]

FABP_METHOD = (
    "1,8-ANS fluorescence displacement against recombinant rat L-FABP "
    "expressed in a prokaryotic system, after Yang et al. 2020; oleic "
    "acid as positive control; all experiments in triplicate")
FABP_PROT = ("1 uM", "rL-FABP in buffer, 100 uL per well")
FABP_LIG = ("1 uL of PFAS at different concentrations added to 100 uL of "
            "the 1,8-ANS-rL-FABP complex with 1,8-ANS at saturation", "")

UF_METHOD = (
    "ultrafiltration (30 kDa Millipore, 1500 g, 30 min) then UPLC-MS/MS, "
    "adapted from Gao et al. 2019 and Sheng et al. 2016; 1 h at 37 C")
UF_PROT = ("1", "uM")
UF_LIG = ("1 uL of a 0.2-15 mM PFAS stock added to 100 uL of 1 uM RSA", "mM")

ITC_METHOD = (
    "isothermal titration calorimetry (Nano-ITC, TA Instruments), 37 C, "
    "350 rpm, 25 injections of 1.96 uL at 150 s intervals, independent "
    "model; 1% (v/v) DMSO in both cell and syringe")
ITC_PROT = ("0.02 mM RSA in the sample cell (effective volume 170 uL)",
            "mM")
ITC_LIG = ("0.6 mM PFAS in the injection syringe (effective volume 50 uL)",
           "mM")

BUFFER = "0.1 M PBS pH 7.4"
UF_CAUTION = (
    "the paper reports that 20-95% of PFAS adsorbs to the ultrafiltration "
    "membrane at the lower concentration (2 uM), significantly correlated "
    "with log Kow (p < 0.001), so these affinities are biased high; the "
    "paper uses the ITC values in its subsequent analysis")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []
    for (chem, ncf, head, fk, fksd,
         un, unsd, uk, uksd, inn, insd, ik, iksd) in T1:
        base = dict(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            n_perfluorinated_carbons="" if ncf is None else str(ncf),
            head_group=head, temperature_c="37", buffer=BUFFER,
            is_computational="no",
            tk_relevance_noted=(
                "yes - the paper's PBPK model uses these constants as fixed "
                "parameters and shows the PFAS liver-blood partition tracks "
                "the Ka(rL-FABP)/Ka(RSA) ratio, R2 0.93 male and 0.67 "
                "female rats, extended to humans at R2 0.67"))

        r = blank(fields)
        r.update(base)
        r.update(
            protein="liver fatty acid-binding protein (rL-FABP), "
                    "recombinant",
            protein_species="rat", method=FABP_METHOD,
            kd=fk, kd_units="uM" if fk else "",
            n_binding_sites="1 (the paper states L-FABP has a single "
                            "binding cavity for the fluorescent probe)",
            protein_conc=FABP_PROT[0], protein_conc_units=FABP_PROT[1],
            ligand_conc_range=FABP_LIG[0], ligand_conc_units=FABP_LIG[1],
            source_table=(
                "Table 1, rL-FABP Kd column (this study)"
                + (f" | Kd {fk} +/- {fksd} uM" if fk else
                   " | binding UNDETECTABLE; the Results text attributes "
                   "this to the fluorescence method's insensitivity or to "
                   "PFBA's low affinity, and notes the same result in "
                   "Yang et al. 2020 on human L-FABP")))
        out.append(r)

        if uk:
            r = blank(fields)
            r.update(base)
            r.update(
                protein="rat serum albumin (RSA), fatty-acid-free, "
                        "purified from rat serum",
                protein_species="rat", method=UF_METHOD,
                kd=uk, kd_units="uM", n_binding_sites=un,
                protein_conc=UF_PROT[0], protein_conc_units=UF_PROT[1],
                ligand_conc_range=UF_LIG[0], ligand_conc_units=UF_LIG[1],
                assay_volume_note=UF_CAUTION,
                source_table=(
                    "Table 1, RSA ultrafiltration columns | "
                    f"Kd {uk} +/- {uksd} uM | n {un} +/- {unsd}"))
            out.append(r)

        if ik:
            r = blank(fields)
            r.update(base)
            r.update(
                protein="rat serum albumin (RSA), fatty-acid-free, "
                        "purified from rat serum",
                protein_species="rat", method=ITC_METHOD,
                kd=ik, kd_units="uM", n_binding_sites=inn,
                protein_conc=ITC_PROT[0], protein_conc_units=ITC_PROT[1],
                ligand_conc_range=ITC_LIG[0], ligand_conc_units=ITC_LIG[1],
                source_table=(
                    "Table 1, RSA ITC columns | "
                    f"Kd {ik} +/- {iksd} uM | n {inn} +/- {insd} | the "
                    "paper's preferred set, ultrafiltration being biased "
                    "by membrane adsorption"))
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
    removed = before - len(existing)

    new = rows(fields)
    print(f"removed {removed} existing {STUDY} {YEAR} rows")
    print(f"adding  {len(new)} rows -> {len(existing) + len(new)} total")
    by_protein = {}
    for r in new:
        by_protein.setdefault(
            ("rL-FABP" if "FABP" in r["protein"] else "RSA")
            + ("/ITC" if "calorimetry" in r["method"] else
               "/UF" if "ultrafiltration" in r["method"] else ""), []
        ).append(r)
    for k, v in sorted(by_protein.items()):
        got = sum(1 for r in v if r["kd"])
        print(f"   {k:<12} {len(v):>3} rows, {got} with a Kd")

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
