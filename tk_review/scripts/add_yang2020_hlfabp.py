#!/usr/bin/env python3
"""Yang et al. 2020: 74 PFAS screened against human L-FABP, SI Table S2.

Yang D, Han J, Hall DR, Sun J, Fu J, Kutarna S, Houck KA, LaLone CA,
Doering JA, Ng CA, Peng H.  "Nontarget Screening of Per- and
Polyfluoroalkyl Substances Binding to Human Liver Fatty Acid Binding
Protein."  Environ Sci Technol 2020;54(9):5676-5686.
doi 10.1021/acs.est.0c00049, PMID 32249562.
Read from the article PDF and its Supporting Information.

WHY IT MATTERS OUT OF PROPORTION TO ITS ROWS: this is the measurement
two other entries in this database rest on.  Fan et al. 2025 follows
its 1,8-ANS displacement protocol for rat L-FABP, and Zhao et al. 2023
validates a machine-learning QSAR over 1599 compounds against these Kd
at Pearson r = 0.79.  Hui Peng is this paper's corresponding author and
a co-author on Fan, so the lineage is one group, not an independent
replication.

A LABELLING ERROR IN THE PUBLISHED MAIN TEXT.  The Figure 1 caption
reads "PFAS 46, 57, and 9 refer to PFOS, PFHxS, and NEt-FOSE as
described in Table S1."  Table S1 does not agree: chemical 57 is
POTASSIUM PERFLUOROBUTANESULFONATE (PFBS, CAS 29420-49-3) and chemical
26 is potassium perfluorohexanesulfonate (PFHxS, CAS 3871-99-6).

Table S1 is right and the caption is wrong, on three independent
checks.  The library carries several compounds as both the acid and a
salt, and in Table S2 each pair lands on the same Kd:

    PFOS   acid (46) 8.10 uM    potassium salt (6)  8.10 uM   1.00x
    PFOA   acid (74) 18.74 uM   ammonium salt (21)  18.74 uM  1.00x
    PFBS   acid (73) 116.40 uM  potassium salt (57) 125.25 uM 1.08x

Chemical 57 pairs with 73, which Table S1 names perfluorobutanesulfonic
acid; a PFHxS would not.  The Kd ordering also fits: PFNA 4.30 < PFOS
8.10 < PFHxS 11.19 < PFOA 18.74 < PFBS ~120 < PFHxA 140.76, which is
the chain-length pattern every other L-FABP source shows.  The caption
is followed by nothing here.

WHAT TABLE S2 GIVES: Kd for 19 of the 74.  The caption states the other
55 are inactive and that NMe-FOSA (16) is excluded despite a high peak
intensity because its Kd exceeds 1000 uM, which is 19 + 55 = 74 and
20 ligands counting NMe-FOSA.  The 55 inactive are not listed
individually, so a compound absent from Table S2 is recorded as
"> 1000 uM" only where Table S1 names it and this review needs it --
PFBA (chemical 12) is the case that matters.

ASSAY CONCENTRATIONS ARE FULLY REPORTED, unusually for this
collection: 1 uM hL-FABP, 1,8-ANS at 40 uM, PFAS at 300 uM for the
initial screen, from 30 mM DMSO stocks.

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
           "overexpressed in E. coli BL21(DE3) from plasmid pNIC28-Bsa4")
METHOD = ("1,8-ANS fluorescence displacement: each PFAS spiked into a "
          "mixture of hL-FABP and 1,8-ANS, the fall in fluorescence "
          "measuring displacement of the probe from the hydrophobic "
          "pocket; the 30 positive hits taken on to dose-response "
          "curves, Kd derived from IC50 or IC10")
PROT_CONC = ("1 uM hL-FABP, with 1,8-ANS at 40 uM as the displaced "
             "probe", "uM")
LIG_CONC = ("300 uM for the initial screen of all 74, from 30 mM DMSO "
            "stocks; dose-response series for the 30 positive hits",
            "uM")
BUFFER = ("buffer A, 50 mM Tris-HCl, 150 mM NaCl, 1 mM dithiothreitol, "
          "pH 8.0 (SECC arm); 10 mM PBS as the SEC mobile phase")

# SI Table S2, every row, joined to SI Table S1 by chemical number.
# (No, name as this review records it, CAS, nCF, head group,
#  IC50 uM or "", IC10 uM or "", Kd uM, peak intensity, RT min)
S2 = [
    (17, "PFNA", "375-95-1", 8, "carboxylate",
     "13.88", "", "4.30", "2.54E+08", "4.16"),
    (60, "PFOSA (FOSA, perfluorooctanesulfonamide)", "754-91-6", 8,
     "sulfonamide", "18.92", "", "5.87", "3.85E+08", "4.81"),
    (6, "PFOS, potassium salt", "2795-39-3", 8, "sulfonate",
     "26.12", "", "8.10", "5.25E+08", "4.58"),
    (46, "PFOS", "1763-23-1", 8, "sulfonate",
     "26.12", "", "8.10", "5.25E+08", "4.58"),
    (34, "perfluoro-3,6,9-trioxatridecanoic acid", "330562-41-9", None,
     "ether carboxylate", "26.49", "", "8.22", "2.88E+08", "4.87"),
    (26, "PFHxS, potassium salt", "3871-99-6", 6, "sulfonate",
     "36.06", "", "11.19", "3.96E+08", "3.97"),
    (24, "N-ethylperfluorooctanesulfonamide (NEtFOSA)", "4151-50-2", 8,
     "sulfonamide", "53.46", "", "16.58", "2.14E+08", "5.52"),
    (21, "PFOA, ammonium salt", "3825-26-1", 7, "carboxylate",
     "60.39", "", "18.74", "3.76E+08", "3.84"),
    (74, "PFOA", "335-67-1", 7, "carboxylate",
     "60.39", "", "18.74", "3.76E+08", "3.84"),
    (73, "PFBS", "375-73-5", 4, "sulfonate",
     "", "41.69", "116.40", "3.74E+08", "3.29"),
    (57, "PFBS, potassium salt", "29420-49-3", 4, "sulfonate",
     "", "44.86", "125.25", "3.74E+08", "3.29"),
    (61, "3H-perfluoro-4-hydroxy-3-penten-2-one", "1694-30-0", None,
     "alcohol/ketone", "", "48.44", "135.24", "1.67E+07", "2.39"),
    (20, "PFHxA", "307-24-4", 5, "carboxylate",
     "", "50.42", "140.76", "7.43E+07", "3.22"),
    (14, "2H,2H,3H,3H-perfluorooctanoic acid", "914637-49-3", None,
     "carboxylate", "", "56.60", "158.03", "9.54E+07", "4.26"),
    (54, "perfluoro-3,6-dioxaheptanoic acid", "151772-58-6", None,
     "ether carboxylate", "", "90.37", "252.31", "3.64E+08", "3.16"),
    (4, "4H-perfluorobutanoic acid", "679-12-9", None, "carboxylate",
     "", "125.90", "351.53", "2.72E+07", "1.23"),
    (52, "difluoromethyl 1H,1H-perfluoropropyl ether", "56860-81-2",
     None, "ether", "", "256.14", "715.15", "7.36E+07", "3.09"),
    (28, "perfluoro(4-methoxybutanoic) acid", "863090-89-5", None,
     "ether carboxylate", "", "291.88", "814.93", "4.31E+07", "2.99"),
    (48, "1H,1H,5H-perfluoropentanol", "355-80-6", None, "alcohol",
     "", "328.10", "916.06", "6.70E+06", "3.12"),
]

# Named in Table S1, absent from Table S2, therefore among the 55 the
# caption calls inactive. Only the ones this review needs are recorded.
INACTIVE = [
    (12, "PFBA", "375-22-4", 3, "carboxylate"),
    (16, "N-methylperfluorooctanesulfonamide (NMeFOSA)", "31506-32-8",
     8, "sulfonamide"),
    (9, "N-ethyl-N-(2-hydroxyethyl)perfluorooctanesulfonamide "
        "(NEt-FOSE)", "1691-99-2", 8, "sulfonamide"),
]

BANDS = [
    ("9 high-affinity PFAS of 74 (band)", "4.3 to 18.7 uM"),
    ("11 low-affinity PFAS of 74 (band)", "116 to 916 uM"),
    ("54 inactive PFAS of 74 (band)", "> 1000 uM"),
]

TK = ("yes - the experimental anchor under two other entries here: Fan "
      "et al. 2025 follows this protocol for rat L-FABP and Zhao et al. "
      "2023 validates a 1599-compound QSAR against these Kd at r = "
      "0.79. Same group (Hui Peng), so not an independent replication")
CAPTION_NOTE = (
    " | NOTE: the published Figure 1 caption calls chemical 57 'PFHxS'; "
    "SI Table S1 names it potassium perfluorobutanesulfonate and names "
    "chemical 26 as potassium perfluorohexanesulfonate. Table S1 is "
    "followed here -- see the script docstring for the three checks")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []

    def base(chem, cas, ncf, head):
        return dict(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            n_perfluorinated_carbons="" if ncf is None else str(ncf),
            head_group=head, protein=PROTEIN, protein_species="human",
            method=METHOD, protein_conc=PROT_CONC[0],
            protein_conc_units=PROT_CONC[1],
            ligand_conc_range=LIG_CONC[0],
            ligand_conc_units=LIG_CONC[1], buffer=BUFFER,
            is_computational="no", tk_relevance_noted=TK,
            n_binding_sites=("not fitted here; the paper cites the "
                             "palmitic acid crystal structure showing "
                             "two sites"),
            binding_site_identified=(
                "Arg-122, Ser-39 and Ser-124 for the acid head group; "
                "Phe-50 and Ile-52 for the alkyl chain"))

    for (no, name, cas, ncf, head, ic50, ic10, kd, inten, rt) in S2:
        r = blank(fields)
        r.update(base(name, cas, ncf, head))
        note = (f"SI Table S2, chemical {no} (CAS {cas}) | Kd {kd} uM | "
                + (f"IC50 {ic50} uM" if ic50 else f"IC10 {ic10} uM")
                + f" | LC-MS peak intensity {inten} at 1 uM | C18 "
                  f"retention {rt} min")
        if no in (57, 26):
            note += CAPTION_NOTE
        r.update(kd=kd, kd_units="uM", source_table=note)
        out.append(r)

    for (no, name, cas, ncf, head) in INACTIVE:
        r = blank(fields)
        r.update(base(name, cas, ncf, head))
        r.update(
            kd="", kd_units="> 1000 uM, no binding detected",
            source_table=(
                f"SI Table S1 chemical {no} (CAS {cas}), ABSENT from "
                f"Table S2, whose caption states the other 55 of the 74 "
                f"are inactive"
                + (" and that NMe-FOSA (16) is excluded despite a very "
                   "high peak intensity because its Kd exceeds 1000 uM"
                   if no == 16 else "")))
        out.append(r)

    for name, band in BANDS:
        r = blank(fields)
        r.update(base(name, "", None, ""))
        r.update(
            kd="",
            kd_units=f"{band} (band over the group, not a point value)",
            source_table="Results text and Figure 2A")
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # The acid/salt pairs are the evidence for the Table S1 mapping.
    # If a future edit breaks one, stop rather than publish the error.
    kd = {no: float(k) for no, _, _, _, _, _, _, k, _, _ in S2}
    for acid, salt, nm in ((46, 6, "PFOS"), (74, 21, "PFOA"),
                           (73, 57, "PFBS")):
        q = kd[salt] / kd[acid]
        assert 0.9 < q < 1.15, f"{nm}: acid/salt Kd differ by {q:.2f}x"
    print("acid/salt check passed: PFOS 1.00x, PFOA 1.00x, PFBS 1.08x "
          "-- Table S1's identities are self-consistent")

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, existing = reader.fieldnames, list(reader)

    before = len(existing)
    existing = [r for r in existing
                if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"removed {before - len(existing)} existing {STUDY} {YEAR} rows")

    new = rows(fields)
    pts = sum(1 for r in new if r["kd"])
    print(f"adding  {len(new)} rows ({pts} with a Kd, "
          f"{len(new) - pts} bounds or bands) -> "
          f"{len(existing) + len(new)} total\n")
    for r in sorted((r for r in new if r["kd"]),
                    key=lambda x: float(x["kd"])):
        print(f"   {float(r['kd']):>8.2f} uM   {r['chemical'][:52]}")

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
