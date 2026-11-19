#!/usr/bin/env python3
"""IVIVE scaling factors for renal OAT transporters, and the albumin problem.

    python3 -I scripts/add_ivive_scaling_factors.py

Writes db/ivive_scaling_factors.csv. This is what Stage 3 needs to turn a
CL_int measured in transfected cells into a clearance in a kidney.

TWO WAYS TO SCALE, AND THEY DISAGREE.

  REF, relative EXPRESSION factor: the ratio of transporter protein abundance
       between the in vitro system and the tissue, by quantitative proteomics.
  RAF, relative ACTIVITY factor: the ratio of transport rates for a selective
       PROBE substrate between the same two, which needs no proteomics but
       assumes the probe behaves like the test compound.

    Kumar AR, Prasad B, Bhatt DK, Mathialagan S, Varma MVS, Unadkat JD.
    "In Vivo-to-In Vitro Extrapolation of Transporter-Mediated Renal
    Clearance: Relative Expression Factor Versus Relative Activity Factor
    Approach." Drug Metab Dispos 2021;49(6):470-478.
    doi:10.1124/dmd.121.000367  PMID 33824168, PMC8232064.

compared the two over 26 OAT substrates and found RAF the better predictor:
62% and 73% of intrinsic secretory clearances within 2- and 3-fold against
46% and 62% for REF; for plasma secretory clearance, 65%/81% against 50%/69%.

                 REF (HEK293 / kidney)    RAF (probe)
    OAT1               0.16                  0.64     (tenofovir)
    OAT2               0.15                  7.3      (acyclovir, ganciclovir)
    OAT3               0.37                  4.1      (oseltamivir acid,
                                                       benzylpenicillin)

THE DIRECTION OF REF IS NOT A CONVENTION YOU CAN ASSUME, and getting it
backwards is a factor of 40 on OAT1. Kumar defines it as "transporter
abundance in transporter-expressing cells / transporter abundance in the
human kidney", so REF < 1 means the cells UNDER-express and an in vitro
clearance is scaled UP by dividing by REF. Tan 2024 (below) reports REF 8.0
for OAT1 and 32.0 for OAT3 in a different system while stating that its
expression was "considerably lower compared to in vivo" -- which only holds
if their REF is the reciprocal, kidney over cells. Both conventions appear in
the same literature. Every row here carries its direction explicitly.

REF IS ALSO CELL-LINE SPECIFIC, NOT A CONSTANT. Antonescu 2020 measured OAT1
at 1.32 fmol/ug total protein in HEK293-Flp-In-OAT1 against 5.33 fmol/ug in
human renal cortex, a ratio of 0.25, where Kumar's HEK293 line gives 0.16.
Same transporter, same tissue reference, different transfected line, 1.5-fold
apart. A REF borrowed from another laboratory's cells is an estimate.

THE BIGGER PROBLEM FOR THIS PROJECT IS ALBUMIN.

    Tan SPF, Tillmann A, Murby SJ, Rostami-Hodjegan A, Scotcher D,
    Galetin A. "Albumin-Mediated Drug Uptake by Organic Anion Transporter 1/3
    Is Real: Implications for the Prediction of Active Renal Secretion
    Clearance." Mol Pharm 2024;21(9):4603-4617.
    doi:10.1021/acs.molpharmaceut.4c00504  PMID 39166754, PMC11372837.

With 4% human serum albumin present and the UNBOUND substrate concentration
held constant, the unbound intrinsic active uptake clearance rose 3.5- to
122-fold for OAT1 and up to 28-fold for OAT3, for six highly protein-bound
substrates. The effect was absent for weakly bound ones. The proposed
mechanism is facilitated dissociation: albumin interacts with the cell
surface, changes conformation and releases its ligand at the membrane, so the
transporter sees more than the free fraction in bulk solution would predict.
The authors recommend including HSA in routine OAT1/3 assays.

WHY THAT MATTERS HERE MORE THAN ANYWHERE. The substrates showing the effect
had fu,inc between 0.002 and 0.18. PFOA's plasma unbound fraction is 0.0017
(Argoul 2026 Table 1) -- at or below the bottom of that range, and PFOS,
PFNA, PFHxS are lower still. PFAS are the extreme case of the compound class
this effect was demonstrated on.

NO PFAS TRANSPORT ASSAY IN THIS DATABASE IS RECORDED AS INCLUDING PROTEIN.
Of the 52 tier-A rows, across Kimura, Louisse, Nakagawa, Niu, Weaver and
Zhao, not one mentions albumin, serum or BSA in its expression system,
method or notes. Where the buffer is stated outright it is protein-free:
Argoul 2026 used Hank's balanced salt solution, and Krebs-Henseleit for rat
Oat1. For the rest this is an absence of any mention rather than a confirmed
composition, which is the normal practice for transfected-cell uptake assays
but has not been verified paper by paper. If
facilitated dissociation applies to PFAS, every CL_int in
db/transporter_kinetics_working.csv understates OAT1/OAT3-mediated uptake,
possibly by one to two orders of magnitude, and understates it MOST for the
longest-chain, most tightly bound compounds -- which is the same direction as
the chain-length trend the compilation is trying to explain.

This is not a reason to discard the table. It is a reason to state, in the
write-up, that these clearances are lower bounds under protein-free
conditions, and it is a testable prediction: repeat one PFAS uptake assay
with and without HSA. It also connects the strands. The albumin association
constants this project has spent its effort on are exactly the input a
facilitated-dissociation correction needs.

Idempotent.
"""
import csv, io, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                   "ivive_scaling_factors.csv")
KUMAR = ("Kumar 2021, Drug Metab Dispos 49:470, doi:10.1124/dmd.121.000367, "
         "PMC8232064")
TAN = ("Tan 2024, Mol Pharm 21:4603, "
       "doi:10.1021/acs.molpharmaceut.4c00504, PMC11372837")
ANTO = ("Antonescu 2020, Pharmaceutics 12:390, doi:10.3390/pharmaceutics"
        "12040390, PMC7238232")
MATH = ("Mathialagan 2017, Drug Metab Dispos 45:409, "
        "doi:10.1124/dmd.116.074294")

FIELDS = ["factor_type", "transporter", "in_vitro_system", "value",
          "direction", "how_to_apply", "probe_substrate", "note", "source"]
IV_OVER_VIVO = "in vitro / in vivo (cells over kidney)"
VIVO_OVER_IV = "in vivo / in vitro (kidney over cells)"
DIVIDE = "multiply the in vitro CL_int by 1/value to reach the kidney"
MULT = "multiply the in vitro CL_int by value to reach the kidney"

ROWS = [
    ("REF", "OAT1", "OAT1-overexpressing HEK293", "0.16", IV_OVER_VIVO,
     DIVIDE, "",
     "cells under-express ~6.3-fold relative to renal cortex", KUMAR),
    ("REF", "OAT2", "OAT2-overexpressing HEK293", "0.15", IV_OVER_VIVO,
     DIVIDE, "", "OAT2 transports no PFAS (Nakagawa 2008)", KUMAR),
    ("REF", "OAT3", "OAT3-overexpressing HEK293", "0.37", IV_OVER_VIVO,
     DIVIDE, "", "cells under-express ~2.7-fold", KUMAR),
    ("RAF", "OAT1", "OAT1-overexpressing HEK293", "0.64", IV_OVER_VIVO,
     DIVIDE, "tenofovir",
     "4-fold apart from the REF for the same transporter and system", MATH),
    ("RAF", "OAT2", "OAT2-overexpressing HEK293", "7.3", IV_OVER_VIVO,
     DIVIDE, "acyclovir, ganciclovir", "", MATH),
    ("RAF", "OAT3", "OAT3-overexpressing HEK293", "4.1", IV_OVER_VIVO,
     DIVIDE, "oseltamivir acid, benzylpenicillin",
     "11-fold apart from the REF for the same transporter and system: the "
     "choice of method changes a scaled OAT3 clearance by an order of "
     "magnitude", MATH),
    ("REF", "OAT1", "ciPTEC transduced with OAT1", "8.0", VIVO_OVER_IV,
     MULT, "",
     "OPPOSITE CONVENTION to Kumar: the paper states expression in ciPTEC "
     "was considerably lower than in vivo, which only holds if this ratio "
     "is kidney over cells", TAN),
    ("REF", "OAT3", "ciPTEC transduced with OAT3", "32.0", VIVO_OVER_IV,
     MULT, "", "OPPOSITE CONVENTION to Kumar", TAN),
    ("RAF", "OAT1", "ciPTEC transduced with OAT1", "2.0", VIVO_OVER_IV,
     MULT, "clinical OAT1 probe", "", TAN),
    ("RAF", "OAT3", "ciPTEC transduced with OAT3", "161.2", VIVO_OVER_IV,
     MULT, "clinical OAT3 probe",
     "5-fold apart from the same paper's REF of 32.0", TAN),
    ("abundance", "OAT1", "HEK293-Flp-In-OAT1", "1.32", "fmol/ug protein",
     "divide the renal value by this to get a REF", "",
     "against 5.33 fmol/ug in human renal cortex, a ratio of 0.25 where "
     "Kumar's HEK293 line gives 0.16 -- REF is cell-line specific", ANTO),
    ("abundance", "OAT1", "human renal cortex", "5.33", "fmol/ug protein",
     "the in vivo reference value", "",
     "agrees with Thakur 2024 (4.8 pmol/mg) and Niu 2026 (5.3 pmol/mg) "
     "within the measurement spread", ANTO),
]


def main():
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader()
    for r in ROWS:
        w.writerow(dict(zip(FIELDS, r)))
    text = s.getvalue()
    if os.path.exists(OUT) and open(OUT, encoding="utf-8").read() == text:
        print("no change needed (scaling factors already current)")
        return 0
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"wrote {len(ROWS)} rows to db/ivive_scaling_factors.csv")
    print("   REF (HEK293/kidney): OAT1 0.16, OAT2 0.15, OAT3 0.37")
    print("   RAF (probe):         OAT1 0.64, OAT2 7.3,  OAT3 4.1")
    print("   RAF outpredicts REF (65%/81% vs 50%/69% within 2-/3-fold)")
    print("   REF direction is NOT a safe assumption: Kumar and Tan use "
          "opposite conventions")
    print()
    print("   ALBUMIN: 4% HSA raised OAT1 unbound CL_int 3.5-122x and OAT3 "
          "up to 28x")
    print("   for substrates with fu 0.002-0.18. PFOA's plasma fu is 0.0017.")
    print("   Every PFAS assay in this database was run protein-free.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
