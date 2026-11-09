#!/usr/bin/env python3
"""Hebert & MacManus-Spencer 2010, SI Tables S1-S3, read at source.

Hebert PC, MacManus-Spencer LA.  "Development of a Fluorescence Model
for the Binding of Medium- to Long-Chain Perfluoroalkyl Acids to Human
Serum Albumin Through a Mechanistic Evaluation of Spectroscopic
Evidence."  Anal Chem 2010;82(15):6463-6471.
doi 10.1021/ac100721e, PMID 20590160.

WHY THIS EXISTS.  The database held ONE row for this paper, a lumped
qualitative entry listing nine compounds in its `chemical` cell with
"Ka ~1e4" -- an abstract-level summary, not data.  The paper reports
per-compound constants.  The gap surfaced while diffing Gao 2025's SI
Table S1 against this database: Gao quotes nine Hebert values we did
not have.  Those were NOT imported from Gao, which is a secondary
compilation already shown to misstate Allendorf's protein
concentrations by nine orders of magnitude.  They were read from
Hebert's own SI.

THE CONSTANT IS SOUND, which had to be checked before using it.
K_Hill is reported in M^-1 while n_Hill runs 1.3-3.3, and that
combination is exactly the defect that got Wu 2024 held out of this
project -- an association constant in M^-n cannot be inverted to a
molar Kd.  Hebert's is not that case.  Testing the two Hill forms
against the paper's own saturation data for PFOS (Figure S23,
K_Hill 8.89e3, n_Hill 2.08, half-saturation observed near 9e-5 M):

    Y = (K[L])^n / (1 + (K[L])^n)   -> half at 1/K      = 1.1e-4 M
    Y =  K[L]^n  / (1 +  K[L]^n)    -> half at K^(-1/n) = 1.3e-2 M

The first lands on the data and the second is 140x off, so K_Hill is
a genuine M^-1 association constant and 1/K_Hill is a valid Kd.

TABLE S2, NOT TABLE S1.  Hebert fits every curve twice: Table S1
against TOTAL PFAA concentration, Table S2 against CALCULATED FREE
PFAA.  A binding constant's denominator is the free concentration, so
S2 is the right set -- the same reasoning that makes c_ref the
condition and c_free the denominator of K for Allendorf.  The S1
value is kept in the source note, because it is what Gao 2025 and
other secondary compilations quote, and the two differ by up to 30 %
for the long chains (PFDA 47.7 -> 62.3 x1e3 M^-1).

REPLICATES ARE AVERAGED.  Hebert reports 2-4 independent fits per
compound per HSA level.  One row per compound x [HSA], carrying the
mean, with every individual fit in the note.  That gives 9 rows and
preserves the one comparison this paper uniquely offers: PFOA and
PFNA measured at BOTH 4 and 10 uM HSA, where the constant falls by
40 % and 58 % as the protein is raised 2.5x.  A binding constant
should not depend on protein concentration; that it does is a datum
about the method, and Hebert's Table S3 standardises for it.

WHAT IS LEFT BLANK, AND WHY.  The ligand concentration range is not
recorded.  The SI gives it only as per-figure data tables behind each
titration, and parsing them did not produce trustworthy values -- two
compounds came out at 90 mM, which is far above their solubility, so
the extraction was picking up saturation values as concentrations.
Rather than store a figure mis-parse, the field is empty and the
ligand:protein ratio with it.  Buffer and temperature are in the main
text, which is paywalled and was not retrieved; Gao 2025's Table S1
says 50 mM phosphate at pH 7, but that is secondary and is recorded
as such rather than as fact.

Run:  python3 tk_review/scripts/add_hebert2010_si.py [--dry-run]
Idempotent: replaces the lumped row, then rewrites its own 9 rows.
"""
import argparse
import csv
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Hebert & MacManus-Spencer", "2010"
PMID, DOI = "20590160", "10.1021/ac100721e"

# SI Table S2 (free PFAA) and Table S1 (total PFAA), verbatim.
# compound: nCF, head, {[HSA] uM: ([S2 K x1e3], [S2 n], [S1 K x1e3])}
DATA = {
 "PFHpA": (6, "carboxylate", {
    4.0:  ([10.25, 9.20], [1.55, 1.325], [9.87, 8.90])}),
 "PFOA":  (7, "carboxylate", {
    4.0:  ([26.5, 22.6], [1.61, 2.08], [24.1, 20.1]),
    10.0: ([17.74, 13.51, 13.46], [1.72, 2.14, 2.53],
           [15.17, 11.33, 11.28])}),
 "PFNA":  (8, "carboxylate", {
    4.0:  ([26.6, 33.6], [2.63, 1.86], [23.10, 29.7]),
    10.0: ([13.10, 12.33, 12.16], [3.25, 2.67, 2.46],
           [10.50, 10.47, 10.36])}),
 "PFDA":  (9, "carboxylate", {
    4.0:  ([57.9, 72.9, 63.3, 55.1], [1.83, 2.49, 2.18, 1.93],
           [47.4, 51.4, 46.7, 45.4])}),
 "PFUnDA": (10, "carboxylate", {        # the paper writes PFUnA
    4.0:  ([53.0, 54.5], [1.98, 2.57], [44.1, 41.6])}),
 "PFHxS": (6, "sulfonate", {
    4.0:  ([14.48, 11.75], [1.76, 1.490], [13.72, 11.25])}),
 "PFOS":  (8, "sulfonate", {
    4.0:  ([9.27, 9.35], [2.27, 2.06], [8.81, 8.89])}),
}

METHOD = ("intrinsic Trp fluorescence quenching of HSA, degree of "
          "saturation Y fitted to a Hill model (eq 6); K_Hill is a "
          "true M^-1 constant, verified against the paper's own "
          "saturation data")
PROTEIN = "human serum albumin (HSA)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    old = [r for r in rows if r["study"].startswith("Hebert")]
    keep = [r for r in rows if not r["study"].startswith("Hebert")]
    print(f"  removing {len(old)} existing Hebert row(s):")
    for r in old:
        print(f"    {r['chemical'][:68]}...  Ka={r['ka']}")

    new = []
    for cm, (ncf, head, levels) in DATA.items():
        for hsa, (k2, n2, k1) in sorted(levels.items()):
            Kfree = st.mean(k2) * 1e3            # M^-1
            Ktot = st.mean(k1) * 1e3
            kd = 1e6 / Kfree                     # uM
            nh = st.mean(n2)
            new.append(dict(
                study=STUDY, year=YEAR, pmid=PMID, doi=DOI,
                chemical=cm, n_perfluorinated_carbons=str(ncf),
                head_group=head, protein=PROTEIN,
                protein_species="human", method=METHOD,
                ka=f"{Kfree:.4g}", ka_units="M^-1 (K_Hill, free-PFAA basis)",
                kd=f"{kd:.4g}", kd_units="uM",
                n_binding_sites=(f"n_Hill {nh:.2f} (mean of "
                                 f"{len(n2)}); a Hill coefficient, NOT "
                                 f"a site count"),
                binding_site_identified="",
                unbound_fraction="",
                protein_conc=f"{hsa:g}", protein_conc_units="uM",
                ligand_conc_range="", ligand_conc_units="",
                temperature_c="not retrieved (main text, paywalled)",
                buffer=("not retrieved (main text, paywalled); Gao 2025 "
                        "SI Table S1 cites 50 mM phosphate pH 7 for this "
                        "study, which is SECONDARY"),
                is_computational="no",
                tk_relevance_noted=(
                    "per-compound HSA constants by Trp quenching. PFOA "
                    "and PFNA measured at BOTH 4 and 10 uM HSA, the only "
                    "within-study test in this database of whether a "
                    "binding constant depends on protein concentration "
                    "-- it does, falling 40 % and 58 % as HSA rises 2.5x"),
                source_table=(
                    f"SI Table S2 (free-PFAA basis), K_Hill "
                    f"{'/'.join(f'{x:g}' for x in k2)} x1e3 M^-1, mean "
                    f"{st.mean(k2):.4g}; n_Hill "
                    f"{'/'.join(f'{x:g}' for x in n2)} | SI Table S1 "
                    f"(total-PFAA basis) {'/'.join(f'{x:g}' for x in k1)} "
                    f"x1e3, mean {st.mean(k1):.4g} -- the value secondary "
                    f"compilations quote | read from the Supporting "
                    f"Information, primary"),
                assay_volume_note=(
                    "ligand range not recorded: the SI gives it only as "
                    "per-figure data tables and the extraction was not "
                    "trustworthy (two compounds parsed to 90 mM, above "
                    "their solubility), so it is left blank rather than "
                    "stored as a figure mis-parse"),
                ligand_over_protein_max=""))

    print(f"\n  adding {len(new)} rows:")
    print(f"    {'PFAS':9}{'[HSA] uM':>10}{'K_free M^-1':>14}{'Kd uM':>9}"
          f"{'n_Hill':>9}{'K_tot':>12}")
    for r in new:
        k1m = r["source_table"].split("mean ")[-1].split(" --")[0]
        print(f"    {r['chemical']:9}{r['protein_conc']:>10}"
              f"{float(r['ka']):>14,.0f}{float(r['kd']):>9.1f}"
              f"{r['n_binding_sites'].split()[1]:>9}"
              f"{float(k1m)*1e3:>12,.0f}")

    if args.dry_run:
        print("\n  (dry run, nothing written)")
        return 0
    out = keep + new
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(out)
    print(f"\n  {len(rows)} -> {len(out)} rows; wrote "
          f"{os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
