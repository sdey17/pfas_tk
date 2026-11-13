#!/usr/bin/env python3
"""Replace the Chen 2020 secondary stub with its six measured rows.

THE PAPER.  Chen Y, et al.  "Interaction of Perfluoroalkyl Carboxylic
Acids of Different Carbon Chain Length with Human Serum Albumin."
Int J Environ Res Public Health 2020;17(4):1319.
doi 10.3390/ijerph17041319, PMID 32085632, open access.

THE STUB carried ka "0.45e4-2.68e4 M^-1" from a secondary, and said
why in its own note: "Table 2 Ka values are not in the saved flattened
text, so the Ka here remains the Starnes 2026 secondary range".  Read
from the PDF, that apparent range is two compounds and three
temperatures collapsed -- 0.45e4 is PFOA at 310 K (Table 2: 0.4463e4),
2.68e4 is PFDA at 296 K (2.6788e4) -- so its geometric midpoint,
Kd 91.06 uM and log K 2.218, was no compound at any temperature.

WHAT IS INGESTED.  Table 2's six modified Stern-Volmer constants.
Only the 310 K pair reaches the sheet, via BY_TEMPERATURE in
build_revised_working_table.py, as Peng 2024's and Qin 2010's do.

CONDITIONS, from Fig. 1's caption, which states the series verbatim:
"[PFCAs] (x10^-6 mol/L) 1-10: 0, 2, 4, 8, 12, 16, 20, 30, 40 and 80.
[HSA] = 2 x 10^-6 mol/L; Buffer: PBS, pH = 7.4."  So HSA 2 uM, PFCAs
to 80 uM, maximum loading 40, PBS pH 7.4, excitation 295 nm.  The stub
already had these right; the PDF confirms them.

PFBA and PFHxA were dosed over the same series and gave NO constants:
Fig. 1E shows both essentially flat, and the FT-IR section reports
that they "have little influence on the HSA structure".  A result, not
a gap, so no rows -- the same situation as Qin 2010's PFPA.

NO STOICHIOMETRY IS REPORTED.  Checked against the full text: the
paper has a Stern-Volmer plot (Table 1) and a modified Stern-Volmer
plot (Table 2) and no double-logarithmic plot, so there is no fitted n
anywhere in it.  The sheet's "n = 1 assumed by the reader, the paper
reports none" is therefore true here, as it is for Qin 2010 and unlike
Klevens 1954 and Sheng 2020.

THE DEFECT, carried on the PFOA rows.  The paper's Table 1 and Table 2
do not cohere.  Eq 2 is the Lehrer modified Stern-Volmer form, F0/dF =
1/(f_a K_a [Q]) + 1/f_a, in which f_a "represents the fraction of
accessible fluorescence"; expanding at small [Q] gives F0/F =
1 + f_a K_a [Q], so the ordinary Stern-Volmer slope IS f_a K_a and
f_a = K_SV/K_a follows from the two tables with no plot reading:

  PFOA   296 K  2.158     PFDA   296 K  0.646
         303 K  2.256            303 K  0.820
         310 K  2.411            310 K  0.986

PFOA needs an accessible fraction above 2 at every temperature.  A
fraction cannot exceed 1, and HSA has a SINGLE tryptophan, Trp-214,
which the 295 nm excitation selects, so there is no second fluorophore
for the excess to hide in.  Stated without the expansion: the two
tables disagree by 2.2 to 2.4 for PFOA and cannot both hold.  Same
method and same shape of failure as Qin 2010, and worse -- Qin's PFOA
f_a runs 1.19 to 1.76 on BSA, which has two tryptophans.

Reported as published with the defect recorded, as Peng 2024 and
Qin 2010 are.  Full working:
scripts/audit_chen2020_modified_sv.py

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/add_chen2020_primary.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "Table 2 | read from the full text, primary 2026-11-13"

METHOD = ("fluorescence quenching, modified Stern-Volmer plot "
          "F0/dF vs 1/[Q]: intercept 1/f_a, slope 1/(f_a K_a)")

CONC = ("CONCENTRATIONS, from Fig. 1's caption verbatim: '[PFCAs] "
        "(x10^-6 mol/L) 1-10: 0, 2, 4, 8, 12, 16, 20, 30, 40 and 80. "
        "[HSA] = 2 x 10^-6 mol/L; Buffer: PBS, pH = 7.4.' Excitation "
        "295 nm, 1.0 cm quartz cell, 5.0 nm slits. Maximum loading 40. "
        "PFBA and PFHxA were dosed over the same series and gave NO "
        "constants: Fig. 1E shows both essentially flat and the FT-IR "
        "section reports they 'have little influence on the HSA "
        "structure' -- a result, not a gap, so they get no rows.")

KA = ("K_a is the modified Stern-Volmer 'effective quenching constant "
      "of accessible fluorophores', not a site binding constant. "
      "Table 1's K_SV is the ordinary Stern-Volmer quenching constant "
      "and is not a binding constant either. The paper reports NO "
      "stoichiometry anywhere -- Stern-Volmer and modified "
      "Stern-Volmer plots only, no double-logarithmic plot -- so n = 1 "
      "is the reader's and 'the paper reports none' is true.")

FA = ("TABLES 1 AND 2 DO NOT COHERE. Eq 2 is the Lehrer form F0/dF = "
      "1/(f_a K_a [Q]) + 1/f_a with f_a 'the fraction of accessible "
      "fluorescence'; expanding at small [Q] gives F0/F = 1 + "
      "f_a.K_a.[Q], so the Stern-Volmer slope IS f_a.K_a and f_a = "
      "K_SV/K_a follows from the two published tables with no plot "
      "reading. PFOA gives 2.158, 2.256 and 2.411 at 296, 303 and "
      "310 K; PFDA gives 0.646, 0.820 and 0.986. A fraction cannot "
      "exceed 1, and HSA has a SINGLE tryptophan, Trp-214, which the "
      "295 nm excitation selects, so there is no second fluorophore "
      "for the excess to hide in. Without the expansion the statement "
      "is only weaker: the two tables disagree by 2.2-2.4 for PFOA and "
      "cannot both hold. Same method and same failure shape as Qin "
      "2010, and worse -- Qin's PFOA f_a runs 1.19-1.76 on BSA, which "
      "has two tryptophans. Full working: "
      "scripts/audit_chen2020_modified_sv.py")

COMMON = {
    "study": "Chen et al.", "year": "2020", "pmid": "32085632",
    "doi": "10.3390/ijerph17041319",
    "protein": "human serum albumin (HSA)", "protein_species": "human",
    "method": METHOD, "ka_units": "M^-1 (modified Stern-Volmer K_a)",
    "kd_units": "uM", "n_binding_sites": "",
    "protein_conc": "2.0", "protein_conc_units": "uM",
    "ligand_conc_range": "0 to 80", "ligand_conc_units": "uM",
    "ligand_over_protein_max": "40", "head_group": "carboxylate",
    "buffer": "PBS, pH 7.4", "is_computational": (
        "no (molecular modelling into sub-domain IIA is reported "
        "alongside in Section 3.4, as a separate computational strand)"),
    "tk_relevance_noted": ("partly - the paper notes PFOA/PFDA 'can be "
                           "stored and transported in vivo by HSA' but "
                           "draws no kinetic conclusion"),
    "binding_site_identified": ("sub-domain IIA, Sudlow site I, by "
                                "molecular modelling; FT-IR shows "
                                "beta-sheet and alpha-helix both fall"),
}

SPEC = [("PFOA", "7", [(296, 0.6153e4), (303, 0.5332e4),
                       (310, 0.4463e4)], True),
        ("PFDA", "9", [(296, 2.6788e4), (303, 2.0093e4),
                       (310, 1.4514e4)], False)]


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    if any(MARK in (r.get("source_table") or "")
           and "Chen" in (r.get("study") or "") for r in rows):
        print("  already ingested; nothing to do")
        return 0

    before = len(rows)
    rows = [r for r in rows
            if not ("Chen" in (r.get("study") or "")
                    and r.get("year") == "2020")]
    removed = before - len(rows)

    added = []
    for chem, ncarb, series, flag in SPEC:
        for T, ka in series:
            note = CONC + " | " + KA + (" | " + FA if flag else "")
            added.append(dict(
                COMMON, chemical=chem, n_perfluorinated_carbons=ncarb,
                ka="%g" % ka, kd="%.4g" % (1e6 / ka),
                temperature_c=f"{T} K",
                source_table=(f"Table 2, {chem} at {T} K | " + MARK),
                assay_volume_note=note))
    for spec in added:
        rows.append({f: spec.get(f, "") for f in fields})

    with open(DB, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  removed {removed} stub row, added {len(added)}:")
    for s in added:
        print(f"    {s['chemical']:5s} {s['temperature_c']:6s} "
              f"K_a {s['ka']:>8s} M^-1  Kd {s['kd']:>7s} uM")
    print("  only the 310 K pair reaches the sheet (BY_TEMPERATURE)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
