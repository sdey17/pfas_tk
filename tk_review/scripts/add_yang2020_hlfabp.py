#!/usr/bin/env python3
"""Yang et al. 2020: 74 PFAS screened against human L-FABP.

Yang D, Han J, Hall DR, Sun J, Fu J, Kutarna S, Houck KA, LaLone CA,
Doering JA, Ng CA, Peng H.  "Nontarget Screening of Per- and
Polyfluoroalkyl Substances Binding to Human Liver Fatty Acid Binding
Protein."  Environ Sci Technol 2020;54(9):5676-5686.
doi 10.1021/acs.est.0c00049, PMID 32249562.  Read from the PDF.

WHY IT MATTERS OUT OF PROPORTION TO ITS ROWS: this is the measurement
two other papers in this database rest on.  Fan et al. 2025 follows its
1,8-ANS displacement protocol for rat L-FABP, and Zhao et al. 2023
validates an entire machine-learning QSAR over 1599 compounds against
its Kd values at Pearson r = 0.79.  Hui Peng is this paper's
corresponding author and a co-author on Fan, so the method lineage is
one group, not an independent replication.

WHAT THE MAIN TEXT GIVES, and what it does not.  The per-congener Kd
table is SI Table S2, which was not supplied, so PFOA and PFNA -- named
in the text as high-affinity ligands -- have no individual number here.
What the main text does give is three affinity BANDS over all 74
compounds, and a handful of named values, and those are ingested.  A
band is stored with an empty kd and the range in kd_units, so nothing
downstream can mistake it for a point estimate.

THE ASSAY CONCENTRATIONS ARE FULLY REPORTED, which is unusual in this
collection: 1 uM hL-FABP, 1,8-ANS at 40 uM, and PFAS fixed at 300 uM
for the initial screen.

Chemical IDs 7, 28, 34 and 57 are the paper's own numbering from SI
Table S1.  The text names 46 = PFOS, 57 = PFHxS, 9 = NEt-FOSE and
16 = NMe-FOSA; 7, 28 and 34 are identified only as ether PFAS, so they
are recorded under the paper's ID with no name invented for them.

Run:  python3 scripts/add_yang2020_hlfabp.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Yang et al.", "2020"
DOI, PMID = "10.1021/acs.est.0c00049", "32249562"

PROTEIN = ("liver fatty acid-binding protein (hL-FABP), His-tagged, "
           "overexpressed in E. coli BL21(DE3)")
METHOD = ("1,8-ANS fluorescence displacement: 74 PFAS spiked into a "
          "mixture of hL-FABP and 1,8-ANS, the fall in fluorescence "
          "measuring displacement of the probe from the hydrophobic "
          "pocket; the 30 positive hits taken on to dose-response "
          "curves for Kd")
PROT_CONC = ("1 uM hL-FABP, with 1,8-ANS at 40 uM as the displaced "
             "probe", "uM")
LIG_CONC = ("300 uM for the initial screen of all 74; dose-response "
            "series for the 30 positive hits", "uM")

# (chemical, kd, kd_units, note).  kd empty where the paper gives a
# band or a bound rather than a value.
ROWS = [
    ("PFOS", "8.1", "uM",
     "chemical 46 | Kd 8.1 uM | named in the text and in Figure 3A"),
    ("PFHxS", "125.25", "uM",
     "chemical 57 | Kd 125.25 uM | IC10 44.86 uM | Figure 1B"),
    ("NEt-FOSE", "", "> 1000 uM, no binding detected",
     "chemical 9 | Kd > 1000 uM | Figure 3A; also the one high-affinity "
     "ligand by fluorescence that SECC failed to capture"),
    ("NMe-FOSA", "", "approx 1000 uM (weak ligand)",
     "chemical 16 | Kd approx 1000 uM | detected under ESI- at high "
     "intensity but a weak ligand, which the paper uses to show the "
     "LC-MS criteria are a prioritisation filter, not a binding test"),
    ("ether PFAS, chemical 34 (SI Table S1)", "8.22", "uM",
     "chemical 34 | Kd 8.22 uM | IC50 26.49 uM | an ether PFAS with "
     "high affinity; the paper notes 6:2 Cl-PFESA and HFPO-TA bind "
     "strongly too"),
    ("ether PFAS, chemical 28 (SI Table S1)", "814.92", "uM",
     "chemical 28 | Kd 814.92 uM | IC10 291.88 uM"),
    ("ether PFAS, chemical 7 (SI Table S1)", "", "> 1000 uM",
     "chemical 7 | Kd > 1000 uM | C18 retention 1.50 min"),
    ("9 high-affinity PFAS of 74 (band, incl. PFOS, PFOA, PFNA)", "",
     "4.3 to 18.7 uM (band over 9 compounds, not a point estimate)",
     "Results text and Figure 2A | the individual values are in SI "
     "Table S2, NOT retrieved, so PFOA and PFNA carry no number here"),
    ("11 low-affinity PFAS of 74 (band)", "",
     "116 to 916 uM (band over 11 compounds, not a point estimate)",
     "Figure 2A"),
    ("54 inactive PFAS of 74 (band)", "",
     "> 1000 uM (band over 54 compounds; 39 were also undetectable "
     "by LC-MS under ESI- and none of those bound)",
     "Figure 2A | 44 of 74 showed no fluorescence displacement at all"),
]

TK = ("yes - the experimental anchor under two other entries in this "
      "database: Fan et al. 2025 follows this protocol for rat L-FABP, "
      "and Zhao et al. 2023 validates a 1599-compound QSAR against "
      "these Kd at r = 0.79. Same group (Hui Peng), so not an "
      "independent replication")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []
    for chem, kd, units, note in ROWS:
        r = blank(fields)
        r.update(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            protein=PROTEIN, protein_species="human", method=METHOD,
            kd=kd, kd_units=units,
            n_binding_sites=(
                "not fitted here; the paper cites the palmitic acid "
                "crystal structure showing two sites"),
            binding_site_identified=(
                "Arg-122, Ser-39 and Ser-124 for the acid head group; "
                "Phe-50 and Ile-52 for the alkyl chain"),
            protein_conc=PROT_CONC[0], protein_conc_units=PROT_CONC[1],
            ligand_conc_range=LIG_CONC[0], ligand_conc_units=LIG_CONC[1],
            buffer=("buffer A, 50 mM Tris-HCl, 150 mM NaCl, 1 mM "
                    "dithiothreitol, pH 8.0 (SECC arm); 10 mM PBS as "
                    "SEC mobile phase"),
            is_computational="no", tk_relevance_noted=TK,
            source_table=note)
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
    pts = sum(1 for r in new if r["kd"])
    print(f"adding  {len(new)} rows ({pts} with a point Kd, "
          f"{len(new) - pts} bands or bounds) -> "
          f"{len(existing) + len(new)} total")
    for r in new:
        print(f"   {r['chemical'][:46]:<48}"
              f"{(r['kd'] + ' uM') if r['kd'] else r['kd_units'][:34]}")

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
