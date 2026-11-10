#!/usr/bin/env python3
"""Which Qin 2024 do the eleven albumin rows come from, and is every
number on them from that paper?

WHY THIS SCRIPT EXISTS.  Two papers by Weiping Qin appeared in Environ Sci
Technol in 2024, both carry a main-text Table 1 WITH A COLUMN HEADED
log D_BSA/w, and the two columns hold different numbers.  Reading the label
"Qin et al. 2024 Table 1, BSA column" is therefore not enough to know what a
row holds; this script decides it from the numbers.

    PLASMA PAPER   Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L.
    "Species Difference?  Bovine, Trout, and Human Plasma Protein Binding of
    Per- and Polyfluoroalkyl Substances."  Environ Sci Technol
    2024;58(23):9954-9966.  doi 10.1021/acs.est.3c10824.  PMC11171458.

        Table 1    BSA Freundlich FITS: log D_BSA/w = m log Cw + b, with R2,
                   in L_water / L_protein, Cw in umol/L.   <- THE VALUE
        SI Tab S5  Ctot, mass concentration of BSA, bound fraction, for the
                   isotherms those fits were made on.      <- THE CONDITIONS

    BASELINE PAPER  Qin W, Henneberger L, Gluge J, Konig M, Escher BI.
    "Baseline Toxicity Model to Identify the Specific and Nonspecific Effects
    of Per- and Polyfluoroalkyl Substances in Cell-Based Bioassays."  Environ
    Sci Technol 2024;58(13):5727-5738.  doi 10.1021/acs.est.3c09950.

        SI Tab S6  the SAME eleven BSA fits, reprinted and uncited -- slope,
                   intercept and R2 agree with the plasma paper's Table 1 to
                   the last digit for all eleven.
        Table 1    those fits EVALUATED AT IC10,free (footnote b), which is a
                   cytotoxic concentration.  Different numbers: PFBA 1.94 not
                   2.91, PFOA 3.98 not 4.38.  Not a binding constant at a
                   stated loading, and not what this database stores.

HOW IT DECIDES.  The plasma paper's own BSA assay is the primary source of
the constants, so the eleven rows should carry its DOI, its Table 1 fit as
the value, and its Table S5 concentrations and bound fractions as the
conditions.  Eight checks per row test exactly that, including the negative
one: the stored value must NOT be the baseline paper's Table 1 number.

    PFOSA fails that last check on a tolerance, not on substance.  Its slope
    is -0.105, so the fit is nearly flat and its intercept (4.28) and its
    value at IC10,free (4.33) are 0.05 apart by construction.  PFUnA is
    exempt outright: n_Fr = 1 makes it concentration-independent, so the two
    tables must agree at 4.75 and there is nothing to tell apart.

Run:  python3 tk_review/scripts/verify_qin2024_bsa_provenance.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

# --- plasma Table 1, BSA column: (slope or None, intercept, R2) ----------
T1_BSA = {
    "PFBA": (-0.298, 2.91, 0.60), "PFHxA": (-0.271, 3.43, 0.88),
    "PFHpA": (-0.305, 3.95, 0.81), "PFOA": (-0.314, 4.38, 0.93),
    "PFNA": (-0.147, 4.52, 0.61), "PFUnA": (None, 4.75, 0.96),
    "HFPO-DA": (-0.493, 3.44, 0.86), "PFHxS": (-0.472, 4.28, 0.92),
    "PFOS": (-0.379, 4.74, 0.92), "6:2 FTSA": (-0.092, 3.86, 0.54),
    "PFOSA": (-0.105, 4.28, 0.50),
}
# --- plasma SI Table S5: (Ctot uM, BSA mg/L, BSA fbound %) ---------------
S5 = {
    "PFBA": (2712, 5000, 46), "PFHxA": (731, 1500, 59),
    "PFHpA": (78, 500, 65), "PFOA": (100, 500, 87),
    "PFNA": (100, 500, 86), "PFUnA": (76, 500, 87),
    "HFPO-DA": (203, 1500, 58), "PFHxS": (54, 100, 57),
    "PFOS": (55, 100, 72), "6:2 FTSA": (106, 500, 65),
    "PFOSA": (96, 500, 78),
}
# --- baseline Table 1, log D_BSA/w column: NOT what the db should hold ---
BASE_T1 = {"PFBA": 1.94, "PFHxA": 2.70, "PFHpA": 3.40, "PFOA": 3.98,
           "PFNA": 4.39, "PFUnA": 4.75, "HFPO-DA": 2.28, "PFHxS": 3.52,
           "PFOS": 4.52, "6:2 FTSA": 3.71, "PFOSA": 4.33}

MW = 66463.0            # g/mol, BSA
PLASMA = "10.1021/acs.est.3c10824"
BASELINE = "10.1021/acs.est.3c09950"

rows = [r for r in csv.DictReader(open(DB, encoding="utf-8"))
        if r["study"].startswith("Qin et al.") and r["year"] == "2024"
        and "Table 1, BSA column" in r["source_table"]]

bad = 0
print(f"{len(rows)} rows cite 'Table 1, BSA column'\n")
hdr = ("chem", "doi ok", "log(ka)", "T1 b", "prot uM", "S5 uM", "ratio",
       "S5 ratio", "fub", "S5 fub", "= base T1?")
print("%-9s %-6s %8s %6s %9s %8s %8s %9s %5s %7s %10s" % hdr)
for r in sorted(rows, key=lambda r: r["chemical"]):
    c = r["chemical"]
    slope, b, r2 = T1_BSA[c]
    ctot, bsa_mgL, fb = S5[c]
    bsa_uM = bsa_mgL / MW * 1e3
    lka = math.log10(float(r["ka"]))
    prot = float(r["protein_conc"])
    ratio = float(r["ligand_over_protein_max"])
    fub = int(r["unbound_fraction"].split("%")[0])
    checks = [
        ("doi", r["doi"] == PLASMA),
        ("value=T1 intercept", abs(lka - b) < 0.006),
        ("equation in cite", (f"{b:.2f}" in r["source_table"]
                              and (slope is None
                                   or f"{slope:.3f}" in r["source_table"]))),
        ("R2 in cite", f"R2 {r2:g}" in r["source_table"]),
        ("prot_conc=S5", abs(prot - bsa_uM) < 0.01),
        ("ratio=S5 Ctot/BSA", abs(ratio - ctot / bsa_uM) < 0.02),
        ("fub=100-S5 fbound", fub == 100 - fb),
        ("NOT baseline T1", abs(lka - BASE_T1[c]) > 0.05 or c == "PFUnA"),
    ]
    fails = [n for n, ok in checks if not ok]
    bad += len(fails)
    print("%-9s %-6s %8.3f %6.2f %9.4g %8.4g %8.4g %9.4g %4d%% %6d%% %10s%s"
          % (c, "ok" if r["doi"] == PLASMA else "WRONG", lka, b, prot,
             bsa_uM, ratio, ctot / bsa_uM, fub, 100 - fb,
             f"{BASE_T1[c]:.2f}", "  <-- " + "; ".join(fails) if fails else ""))
print(f"\n{bad} failed checks over {len(rows)} rows x 8 checks")
print("PFUnA is the one compound where the plasma Table 1 value (4.75) and"
      "\nthe baseline Table 1 value (4.75) coincide, because n_Fr = 1 makes"
      "\nthe equation concentration-independent -- so there is nothing to"
      "\nevaluate and the two tables must agree.")
