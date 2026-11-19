#!/usr/bin/env python3
"""Nakagawa 2008: add the Vmax values from Table 1, and record a units defect.

Nakagawa H, Hirata T, Terada T, Jutabha P, Miura D, Harada KH, Inoue K,
Anzai N, Endou H, Inui K, Kanai Y, Koizumi A. "Roles of organic anion
transporters in the renal excretion of perfluorooctanoic acid."
Basic Clin Pharmacol Toxicol 2008;103(1):1-8.
doi:10.1111/j.1742-7843.2007.00155.x  PMID 18373647.

The record held the four Km values from the abstract and no Vmax. The
reviewer supplied the full text; Table 1 has both.

                 Km (uM)        Vmax (pmol/mg protein/min)
    hOAT1     48.0 +/- 6.4         479.3 +/- 19.9
    rOAT1     51.0 +/- 12.0        421.3 +/- 31.0
    hOAT3     49.1 +/- 21.4        399.1 +/- 54.0
    rOAT3     80.2 +/- 17.8        365.6 +/- 28.7

Transient transfection in HEK293, 2-min uptake of [14C]PFOA, 5-400 uM, so
every Km sits inside the tested range. Uptake is the difference between
transporter-expressing and mock-transfected cells. hOAT2 and rOAT2
transported nothing.

TABLE 1's THIRD COLUMN IS MISLABELLED, AND THE ERROR INVERTS THE RANKING.
It is headed "Km/Vmax (ul/mg protein/min.)" and prints 0.10, 0.12, 0.12,
0.22. The numbers are indeed Km divided by Vmax (48.0/479.3 = 0.100,
51.0/421.3 = 0.121, 49.1/399.1 = 0.123, 80.2/365.6 = 0.219), but ul/mg
protein/min are the units of Vmax/Km, a clearance; Km/Vmax carries the
reciprocal units. The intrinsic clearance is Vmax/Km:

    hOAT1  9.99     rOAT1  8.26     hOAT3  8.13     rOAT3  4.56  uL/min/mg

so a reader taking the printed 0.10 as a clearance is wrong by a factor of
100, and the ranking reverses: by the printed column rOAT3 is highest at
0.22, by actual clearance it is the LOWEST of the four.

WHY THIS PAPER MATTERS FOR AIM 1. It is one of the few studies to measure
human and rat orthologues in the SAME expression system, which is the only
way to compare species without mixing laboratories. Its own conclusion is
that "the species differences in its renal elimination are not attributable
to affinity differences in these OATs between human beings and rats", and
the clearances above bear that out: human/rat is 1.21-fold on OAT1 and
1.78-fold on OAT3, against the roughly 50-fold difference in half-life.
Argoul 2026 reaches the same conclusion by a different route.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "transporter_kinetics.csv")
DOI = "10.1111/j.1742-7843.2007.00155.x"
STUDY = ("Nakagawa H, Hirata T, Terada T, Jutabha P, Miura D, Harada KH, "
         "Inoue K, Anzai N, Endou H, Inui K, Kanai Y, Koizumi A")
SYSTEM = ("HEK293 cells transiently transfected with the OAT cDNA vs "
          "mock-transfected controls; 2-min uptake of [14C]PFOA")
MEAS = ("Km/Vmax (Michaelis-Menten, non-linear least squares), "
        "[14C]PFOA varied 5-400 uM")
DEFECT = (
    "TABLE 1 UNITS DEFECT: the third column is headed \"Km/Vmax (ul/mg "
    "protein/min.)\" and prints {kv}. The number is Km/Vmax, but ul/mg "
    "protein/min are the units of Vmax/Km. The intrinsic clearance here is "
    "Vmax/Km = {cl} uL/min/mg protein. Taking the printed value as a "
    "clearance is wrong by a factor of ~100 and reverses the ranking: by the "
    "printed column rOAT3 looks highest, by clearance it is the lowest")
SPECIES_NOTE = (
    "Human and rat orthologues measured in the SAME expression system, which "
    "is the only way to compare species without also changing laboratory. "
    "Human/rat clearance ratio is 1.21 on OAT1 and 1.78 on OAT3, against a "
    "roughly 50-fold species difference in half-life -- the paper's own "
    "conclusion is that species differences in renal elimination are NOT "
    "attributable to OAT affinity")

#  transporter, species, km, km_sd, vmax, vmax_sd, printed Km/Vmax, CL_int
ROWS = [
    ("OAT1", "human", "48.0", "6.4", "479.3", "19.9", "0.10", "9.99"),
    ("OAT1", "rat", "51.0", "12.0", "421.3", "31.0", "0.12", "8.26"),
    ("OAT3", "human", "49.1", "21.4", "399.1", "54.0", "0.12", "8.13"),
    ("OAT3", "rat", "80.2", "17.8", "365.6", "28.7", "0.22", "4.56"),
]
GENE = {"OAT1": "SLC22A6 / Slc22a6", "OAT3": "SLC22A8 / Slc22a8",
        "OAT2": "SLC22A7 / Slc22a7"}


def build(fields):
    out = []
    for tr, sp, km, ksd, vm, vsd, kv, cl in ROWS:
        r = {f: "" for f in fields}
        r.update({
            "study": STUDY, "year": "2008", "pmid": "18373647", "doi": DOI,
            "transporter": tr, "gene": GENE[tr], "species_orthologue": sp,
            "chemical": "PFOA", "n_perfluorinated_carbons": "7",
            "head_group": "carboxylate", "expression_system": SYSTEM,
            "measurement_type": MEAS,
            "km": f"{km} +/- {ksd}", "km_units": "uM",
            "vmax": f"{vm} +/- {vsd}", "vmax_units": "pmol/mg protein/min",
            "substrate_or_inhibitor": "substrate",
            "direction_reabsorption_or_secretion":
                "secretion (basolateral uptake into the proximal tubule cell)",
            "sex_or_hormone_effect": "not assessed",
            "assay_volume_note": DEFECT.format(kv=kv, cl=cl) + ". "
                                 + SPECIES_NOTE,
            "source_table": (f"Table 1 (full text) | Vmax/Km = {cl} "
                             f"uL/min/mg protein, computed here"),
        })
        out.append(r)
    # the two OAT2 negatives, refreshed from the full text
    for sp in ("human", "rat"):
        r = {f: "" for f in fields}
        r.update({
            "study": STUDY, "year": "2008", "pmid": "18373647", "doi": DOI,
            "transporter": "OAT2", "gene": GENE["OAT2"],
            "species_orthologue": sp, "chemical": "PFOA",
            "n_perfluorinated_carbons": "7", "head_group": "carboxylate",
            "expression_system": SYSTEM,
            "measurement_type": "tested for uptake; no transport detected",
            "substrate_or_inhibitor": "not a substrate",
            "direction_reabsorption_or_secretion":
                "basolateral (predominantly hepatic in human)",
            "sex_or_hormone_effect": "not assessed",
            "assay_volume_note": "neither hOAT2 nor rOAT2 transported PFOA",
            "source_table": "Results and Table 1 (full text)",
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
        print("no change needed (Nakagawa 2008 Table 1 already ingested)")
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
    print(f"Nakagawa 2008: {len(old)} row(s) -> {len(new)}")
    print("   4 rows gain a Vmax from Table 1 (pmol/mg protein/min)")
    print("   Table 1 units defect recorded; CL_int 9.99 / 8.26 / 8.13 / 4.56")
    print("   human/rat clearance ratio 1.21 (OAT1) and 1.78 (OAT3)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
