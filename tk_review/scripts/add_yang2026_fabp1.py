#!/usr/bin/env python3
"""Yang 2026 (PNAS): FABP1 constants, four measured here and five carried over.

Yang D, Wang X, Liu J, Gong Y, Nair P, Sun J, Qian X, Cui J, Zeng H, Dong A,
Harding RJ, Burgess-Brown N, Beyett TS, Song D, Krause HM, Diamond ML,
Bolhuis DL, Brown NG, Arrowsmith CH, Edwards AM, Halabelian L, Peng H.
"High-throughput protein target mapping enables accelerated bioactivity
discovery for ToxCast and PFAS compounds." Proc Natl Acad Sci USA
2026;123(25):e2515940123. doi:10.1073/pnas.2515940123  PMID 42247526.

This is the source the L-FABP section listed as unassessable: its PMC copy is
embargoed to December 2026. The reviewer supplied the paper and both SI parts.

WHAT IT IS. Affinity-selection mass spectrometry against 50 human proteins and
481 chemicals (446 ToxCast, 35 PFAS), 24,050 pairs, finding FABPs the most
ligandable family. AS-MS itself returns enrichment, not affinity; the numbers
used here come from the ORTHOGONAL VALIDATION, a 1,8-ANS fluorescence
displacement assay on 1 uM purified His-tagged human FABP1 -- the same assay
family as Yang 2020, Sheng 2016/2018 and Fan 2025. Sixteen FABP1 ligands were
validated with Kd from 0.045 to 31.0 uM.

FIVE OF ITS SIX LEGACY-PFAS VALUES ARE YANG 2020's, NOT NEW MEASUREMENTS, and
the paper says so for two of them: "PFOS and perfluorooctanoic acid (PFOA) were
found to bind to FABP1 with KD of 8.10 and 18.8 uM, respectively, IN THE
PREVIOUS STUDY (26)", reference 26 being Yang et al. 2020. The arithmetic
confirms it for the rest -- Table S3 against the Yang 2020 rows already in this
database:

    PFOSA     5.87  vs 5.87     identical
    PFOS      8.10  vs 8.10     identical
    PFHxS    11.2   vs 11.19    identical to the digits printed
    N-EtFOSA 16.6   vs 16.58    identical to the digits printed
    PFOA     18.7   vs 18.74    identical to the digits printed
    PFNA      5.49  vs 4.30     DIFFERENT

Five independent re-measurements agreeing to three significant figures is not
credible; one that moved is. So the five are stored with no constant and the
reason recorded, which is the Han 2021 treatment -- a value quoted from another
paper is that paper's value, and entering it twice would turn one measurement
into two and tighten the apparent consensus for free.

FOUR ARE MEASURED HERE and are stored with their constants. PFNA and PFO3TDA
have their dose-response curves printed in Fig. 3D with the concentrations
tested; HFPO-TA's is in SI Fig. S5 with its own series; HFPO-TeA is the fourth.

    PFNA       5.49 uM   tested at 1.78-100 uM (Fig. 3D)
    PFO3TDA    0.37 uM   tested at 0.37-50 uM  (Fig. 3D)
    HFPO-TA    0.54 uM   tested at 0.41-100 uM (SI Fig. S5)
    HFPO-TeA   0.05 uM   series not separately stated

THE SAME LABORATORY DISAGREES WITH ITSELF, AGAIN, AND HARDER THAN SHENG DID.
Where Yang 2026 re-measured rather than carried over, it does not reproduce
Yang 2020:

    PFNA      5.49 uM here      4.30 uM in Yang 2020      1.3x
    PFO3TDA   0.37 uM here      8.22 uM in Yang 2020       22x
    PFO2HpA   no binding here   252.31 uM in Yang 2020
    PFMOBA    no binding here   814.93 uM in Yang 2020

Both Yang rows are kept for PFNA and PFO3TDA, as Sheng 2016 and Sheng 2018 are
both kept. A 22-fold within-laboratory gap on one compound, in the same assay
with the same probe, is a floor on how much of this compilation's spread can be
attributed to anything finer than the measurement itself.

HFPO-TeA AT 0.05 uM IS NOW THE TIGHTEST CONSTANT IN THE WHOLE L-FABP
COMPILATION, ahead of Khazaee 2021's PFOA at 0.099 uM; it takes the sheet's
ceiling from log K 5.85 to 6.15. Sheng 2018 puts the same compound
at 3.35 uM and Jia 2022 at 1.22 uM, so this is 24-67x tighter than either. It
is admitted because it is measured, reported and bracketed, but it is the kind
of value section 5.4 exists to warn about.

NO BINDING DETECTED, recorded without a constant: HFPO-DA/GenX, PFO2HpA,
PFMOBA (short-chain perfluoroether carboxylic acids). N-MeFOSA was tested and
not validated, the paper attributing that to low purity or degradation of the
standard.

Oleic acid at 0.045 uM is the assay's endogenous benchmark, stored and dropped
by the working table as a natural ligand rather than an analyte.

THREE COMPOUNDS GAIN AN ACRONYM. Yang 2020 entered three perfluoroether
carboxylic acids under IUPAC names because no abbreviation was given; Yang 2026
supplies them, matching by CAS. perfluoro-3,6,9-trioxatridecanoic acid is
PFO3TDA (330562-41-9), perfluoro-3,6-dioxaheptanoic acid is PFO2HpA
(151772-58-6), perfluoro(4-methoxybutanoic) acid is PFMOBA (640-688-5). The
working table maps the long names onto the acronyms so one compound is not
counted twice.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "protein_binding.csv")
DOI = "10.1073/pnas.2515940123"

PROTEIN = ("liver fatty acid-binding protein (hL-FABP, FABP1), recombinant "
           "His-tagged, expressed in E. coli, purity by SDS-PAGE and anti-His "
           "western blot (SI Fig. S6)")
METHOD = ("1,8-ANS fluorescence displacement, used as the orthogonal "
          "validation of an affinity-selection mass spectrometry screen "
          "(AS-MS itself returns enrichment, not affinity)")
PROT_C = "1 uM FABP1 (SI Fig. S7)"
CARRIED = ("NOT MEASURED HERE. Table S3 prints this value but it is Yang et "
           "al. 2020's: the main text attributes PFOS and PFOA to \"the "
           "previous study\" explicitly, and all five agree with the Yang 2020 "
           "rows in this database to the digits printed. Stored without a "
           "constant so one measurement is not counted twice")

#  chemical, nC, head, kd, ligand series, note
ROWS = [
    # --- measured here ---
    ("PFNA", "8", "carboxylate", "5.49", "1.78-100 uM, seven concentrations (Fig. 3D)",
     "measured here; Yang et al. 2020, the same laboratory and assay, "
     "reported 4.30 uM for this compound"),
    ("PFO3TDA", "", "ether carboxylate", "0.37",
     "0.37-50 uM, six concentrations (Fig. 3D)",
     "measured here; Yang et al. 2020, the same laboratory and assay, "
     "reported 8.22 uM for this compound as perfluoro-3,6,9-trioxatridecanoic "
     "acid -- a 22-fold within-laboratory difference"),
    ("HFPO-TA", "8", "ether carboxylate", "0.54",
     "0.41-100 uM, six concentrations (SI Fig. S5)",
     "measured here, selected for the structure-activity analysis"),
    ("HFPO-TeA", "11", "ether carboxylate", "0.05", "not separately stated",
     "measured here, selected for the structure-activity analysis. This is the "
     "tightest constant in the L-FABP compilation; Sheng 2018 reports "
     "3.35 uM and Jia 2022 1.22 uM for the same compound"),
    # --- carried over from Yang 2020, stored without a constant ---
    ("PFOA", "7", "carboxylate", "", "", CARRIED + " (18.7 here, 18.74 there)"),
    ("PFOS", "8", "sulfonate", "", "", CARRIED + " (8.10 here, 8.10 there)"),
    ("PFHxS", "6", "sulfonate", "", "", CARRIED + " (11.2 here, 11.19 there)"),
    ("PFOSA", "8", "sulfonamide", "", "", CARRIED + " (5.87 here, 5.87 there)"),
    ("N-EtFOSA", "8", "sulfonamide", "", "", CARRIED + " (16.6 here, 16.58 there)"),
    # --- tested, no constant obtained ---
    ("HFPO-DA", "5", "ether carboxylate", "", "0.41-100 uM, six concentrations",
     "no binding observed (Table S3)"),
    ("PFO2HpA", "", "ether carboxylate", "", "0.37-50 uM, six concentrations",
     "no binding observed (Table S3); Yang et al. 2020 reported 252.31 uM for "
     "this compound as perfluoro-3,6-dioxaheptanoic acid"),
    ("PFMOBA", "", "ether carboxylate", "", "0.37-50 uM, six concentrations",
     "no binding observed (Table S3); Yang et al. 2020 reported 814.93 uM for "
     "this compound as perfluoro(4-methoxybutanoic) acid"),
    ("N-MeFOSA", "8", "sulfonamide", "", "",
     "tested and NOT validated; the paper attributes this to low purity or "
     "degradation of the standard"),
    # --- the assay's benchmark ligand ---
    ("oleic acid", "", "fatty acid (endogenous ligand)", "0.045", "",
     "benchmark for the 1,8-ANS assay, not an analyte"),
]


def build(fields):
    out = []
    for ch, nc, head, kd, lig, note in ROWS:
        r = {f: "" for f in fields}
        r.update({
            "study": "Yang et al.", "year": "2026", "pmid": "42247526",
            "doi": DOI, "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "protein": PROTEIN,
            "protein_species": "human", "method": METHOD,
            "kd": kd, "kd_units": "uM" if kd else "",
            "protein_conc": PROT_C, "protein_conc_units": "uM",
            "ligand_conc_range": lig, "ligand_conc_units": "uM" if lig else "",
            "is_computational": "no",
            "tk_relevance_noted": "yes - FABPs identified as the most "
                                  "ligandable protein family for PFAS",
            "assay_volume_note": note,
            "source_table": ("SI Appendix Table S3, FABP1 block"
                             + (f" | Kd {kd} uM" if kd else " | no constant")),
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
        print("no change needed (Yang 2026 already ingested)")
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
    print(f"Yang 2026: {len(new)} FABP1 rows ingested")
    print("   4 measured here (PFNA, PFO3TDA, HFPO-TA, HFPO-TeA)")
    print("   5 carried over from Yang 2020, stored WITHOUT a constant")
    print("   3 no binding, 1 not validated, 1 benchmark ligand")
    return 0


if __name__ == "__main__":
    sys.exit(main())
