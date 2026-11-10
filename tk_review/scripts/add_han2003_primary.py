#!/usr/bin/env python3
"""Han 2003 read at source, replacing four abstract-level rows.

Han X, Snow TA, Kemper RA, Jepson GW.  "Binding of perfluorooctanoic
acid to rat and human plasma proteins."  Chem Res Toxicol
2003;16(6):775-781.  doi 10.1021/tx034005w, PMID 12807361.
Closed access; the full text was supplied directly.

WHAT WAS THERE.  Four rows, every one sourced "abstract", carrying
the abstract's range "0.3-0.4" mM with NO stoichiometry and no
conditions.  Downstream that became a single Kd of 346.4 uM -- the
geometric mean of the range -- converted with n = 1.

THE STOICHIOMETRY IS THE CORRECTION.  Table 1 reports n, and it is
not 1:

                        RSA              HSA
    Kd (mM)  NMR      0.29 +/- 0.10       -
             micro-SEC 0.36 +/- 0.08   0.38 +/- 0.04
    n        micro-SEC 7.8 +/- 1.5     7.2 +/- 1.3

Seven to eight sites, not one.  Carrying n = 1 understates log K by
log10(7.2) = 0.86, which is most of the gap between this study and
the rest of the file.  The abstract says "6-9"; Table 1 gives the
fitted values per protein, and those are what go in.

THE TWO METHODS ARE SEPARATE ROWS, not an averaged range.  NMR and
micro-SEC are different measurements of the same quantity, and the
paper itself checks them against each other: "the difference of Kd
values determined by NMR and micro-SEC is statistically
insignificant" by unpaired t-test.  Keeping them apart preserves
that agreement as evidence; collapsing them to "0.3-0.4" threw it
away.  The NMR row carries no n, because the NMR experiment does not
measure one -- only micro-SEC does.

THE SPECIES COMPARISON IS THE POINT OF THE PAPER, and it survives
only if RSA and HSA stay separate rows: Kd 0.36 vs 0.38 mM and n 7.8
vs 7.2, i.e. "no significant difference was found between PFOA
binding to RSA and PFOA binding to HSA".  This is one of very few
within-study, same-method rat/human albumin comparisons in the file.

CONCENTRATIONS, all stated:

  micro-SEC   50 or 60 uM protein titrated with 0.1-3 mM PFOA in PBS,
              final volume 70 uL, with 0.5 uM 14C-PFOA for counting;
              30 min at room temperature.  Ratio up to 3000/50 = 60
  19F NMR     RSA at 12.5 or 25 uM against a wide PFOA range; the two
              protein levels give Kd 0.31 +/- 0.15 and 0.27 +/- 0.05
              mM, whose average is the 0.29 in Table 1.  Isotherms at
              23 C, PBS pH 7.4
  ESI MS      25 uM RSA with 200 uM PFOA in 50 uM ammonium
              bicarbonate pH 7.0, ratio 8 -- qualitative, gas-phase
              confirmation of the complex, no constant

A NOTE ON ULTRAFILTRATION, worth keeping because other studies in
this file rely on it: Han tried it and it FAILED, because "PFOA
completely nonspecifically bound to the membrane of the
ultrafiltration unit".  That is a direct caution against comparing
ultrafiltration-derived constants without checking membrane losses.

Run:  python3 tk_review/scripts/add_han2003_primary.py [--dry-run]
Idempotent: removes any Han 2003 rows, then writes its own.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Han et al.", "2003"
PMID, DOI = "12807361", "10.1021/tx034005w"

MICROSEC = ("microdesalting column separation (Bio-Spin), Scatchard "
            "variant [L]F/r vs [L]F: slope 1/n, intercept Kd/n")
NMR = "19F NMR titration, chemical-shift change of the PFOA CF2 resonance"
ESIMS = "electrospray ionization MS, native conditions (gas phase)"

COND_SEC = ("n is stored bare for the converter; Table 1 gives it as "
            "7.8 +/- 1.5 (RSA) and 7.2 +/- 1.3 (HSA), from three "
            "independent experiments. protein_conc is stored as 50, the "
            "LOWER of the two levels used, so that with the top ligand it "
            "gives the maximum ratio the ratio column is defined as. "
            "micro-SEC: 50 or 60 uM protein titrated with 0.1-3 mM PFOA in "
            "PBS, final volume 70 uL, 0.5 uM 14C-PFOA included for "
            "counting, 30 min at room temperature; protein recovery "
            "(103 +/- 6)%. Ligand:protein reaches 3000/50 = 60")
COND_NMR = ("protein_conc stored as 12.5, the lower of the two levels. "
            "19F NMR: RSA at 12.5 or 25 uM titrated over a wide PFOA range, "
            "isotherms at 23 C in PBS pH 7.4. The two protein levels give "
            "Kd 0.31 +/- 0.15 and 0.27 +/- 0.05 mM and Table 1 reports "
            "their average, 0.29 +/- 0.10; the paper notes the y-intercept "
            "is -Kd and is independent of protein concentration")
COND_MS = ("ESI MS: 25 uM fatty-acid-free RSA with 200 uM PFOA in 50 uM "
           "ammonium bicarbonate pH 7.0, ligand:protein 8. Qualitative "
           "confirmation of a noncovalent complex in the gas phase; no "
           "constant is derived from it")
UF_WARNING = ("ULTRAFILTRATION FAILED in this study and the authors say why: "
              "'PFOA completely nonspecifically bound to the membrane of the "
              "ultrafiltration unit' (Amicon Ultrafree, Biomax). They used "
              "microdesalting columns instead. A direct caution against "
              "pooling ultrafiltration-derived constants without checking "
              "membrane losses")
AGREE = ("NMR and micro-SEC are kept as SEPARATE rows rather than averaged "
         "into the abstract's '0.3-0.4 mM': the paper tests them against "
         "each other and reports that 'the difference of Kd values "
         "determined by NMR and micro-SEC is statistically insignificant' "
         "(unpaired t-test, 95%). That agreement is evidence and is lost "
         "if the two are collapsed")
SPECIES = ("SPECIES COMPARISON, same method, same study: RSA Kd 0.36 mM with "
           "n 7.8 against HSA Kd 0.38 mM with n 7.2 -- 'no significant "
           "difference was found between PFOA binding to RSA and PFOA "
           "binding to HSA'. One of very few within-study rat/human albumin "
           "comparisons in this database")

ROWS = [
 dict(protein="rat serum albumin (RSA)", species="rat", method=MICROSEC,
      kd="360", n="7.8", cond=COND_SEC,
      src="Table 1, micro-SEC column, RSA",
      extra=SPECIES + " | " + AGREE + " | " + UF_WARNING,
      prot="50", lig="100 to 3000", ratio="60"),
 dict(protein="human serum albumin (HSA)", species="human", method=MICROSEC,
      kd="380", n="7.2", cond=COND_SEC,
      src="Table 1, micro-SEC column, HSA",
      extra=SPECIES + " | " + AGREE + " | " + UF_WARNING,
      prot="50", lig="100 to 3000", ratio="60"),
 dict(protein="rat plasma proteins (whole plasma, SEC + ligand blotting)",
      species="rat", method="size-exclusion chromatography + ligand blotting",
      kd="", n="", cond=("whole plasma from 14C-PFOA dosed rats, SEC on a "
                         "Bio-Sil SEC 125-5 column eluted with PBS, "
                         "fractions counted; bound PFOA identified on "
                         "membrane by autoradiography"),
      src="Results, SEC and ligand-blotting sections",
      extra=("QUALITATIVE, no constant: establishes that the PRIMARY PFOA "
             "binding protein in plasma is serum albumin, which is the "
             "premise the albumin constants rest on. Kept as a row because "
             "that identification is a finding, not a condition"),
      prot="", lig="", ratio=""),
 dict(protein="rat serum albumin (RSA)", species="rat", method=NMR,
      kd="290", n="", cond=COND_NMR,
      src="Table 1, NMR row, RSA",
      extra=AGREE + " | n is EMPTY because the NMR experiment does not "
                    "measure a stoichiometry; only micro-SEC does",
      prot="12.5", lig="", ratio=""),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    old = [r for r in rows if r["study"] == STUDY and r["year"] == YEAR]
    keep = [r for r in rows if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"  removing {len(old)} abstract-level Han 2003 row(s):")
    for r in old:
        print(f"    {r['protein'][:34]:36} kd={r['kd'] or '-':10} "
              f"n={r['n_binding_sites'] or '-':6} src={r['source_table']}")

    new = []
    for d in ROWS:
        new.append(dict(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI,
            chemical="PFOA", n_perfluorinated_carbons="7",
            head_group="carboxylate", protein=d["protein"],
            protein_species=d["species"], method=d["method"],
            ka=(f"{1e6/float(d['kd']):.4g}" if d["kd"] else ""),
            ka_units=("M^-1, computed as 1/K_D" if d["kd"] else ""),
            kd=d["kd"], kd_units=("uM" if d["kd"] else ""),
            n_binding_sites=d["n"], binding_site_identified="",
            unbound_fraction="",
            protein_conc=d["prot"], protein_conc_units="uM",
            ligand_conc_range=d["lig"], ligand_conc_units="uM",
            temperature_c="23 (NMR isotherms); room temperature (micro-SEC)",
            buffer="PBS pH 7.4", is_computational="no",
            tk_relevance_noted=(
                "yes - the paper's own conclusion is that 'greater than 90% "
                "of PFOA would be bound to serum albumin in both rat and "
                "human blood', from these constants and plasma albumin"),
            source_table=f"{d['src']} | read from the full text, primary",
            assay_volume_note=d["cond"] + " | " + d["extra"] + " | " + COND_MS,
            ligand_over_protein_max=d["ratio"]))

    print(f"\n  adding {len(new)} rows:")
    print(f"    {'protein':30}{'method':14}{'Kd uM':>8}{'n':>14}{'ratio':>7}")
    for r in new:
        m = ("micro-SEC" if "microdesalting" in r["method"]
             else "19F NMR" if "NMR" in r["method"] else "SEC+blot")
        print(f"    {r['protein'][:29]:30}{m:14}{r['kd']:>8}"
              f"{r['n_binding_sites'] or '-':>14}{r['ligand_over_protein_max'] or '-':>7}")

    if args.dry_run:
        print("\n  (dry run, nothing written)")
        return 0
    out = keep + new
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader(); w.writerows(out)
    print(f"\n  {len(rows)} -> {len(out)} rows; wrote "
          f"{os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
