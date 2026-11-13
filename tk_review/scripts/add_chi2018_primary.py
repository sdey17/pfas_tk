#!/usr/bin/env python3
"""Replace the Chi 2018 stub with its four trace-limit constants.

THE PAPER.  Chi Q, Li Z, Huang J, Ma J, Wang X.  "Interactions of
perfluorooctanoic acid and perfluorooctanesulfonic acid with serum
albumins by native mass spectrometry, fluorescence and molecular
docking."  Chemosphere 2018;198:442-449.
doi 10.1016/j.chemosphere.2018.01.152, PMID 29425944.

THE STUB was one row -- chemical "PFOA and PFOS", protein_species
"human; bovine", ka "0.52e5-9.90e5 M^-1" from a secondary, geometric
midpoint Kd 4.407 uM, log K 3.533.  Table 1 holds TWENTY-FOUR
constants and those two numbers are exactly their minimum and maximum.

WHAT IS INGESTED, AND WHY ONLY FOUR.  Chi's constants are STEPWISE:

    P   + L <=> PL     K_a,1 = [PL]  /([P]  [L])
    PL  + L <=> PL2    K_a,2 = [PL2] /([PL] [L])
    PL2 + L <=> PL3    K_a,3 = [PL3] /([PL2][L])

so K_a,2 is not a second site's affinity, it is the affinity for
adding a second ligand to a protein that already holds one.  log K is
a TRACE-LIMIT quantity, where only the first binding event exists, so
K_a,1 is the only one of the three that answers the question this
compilation asks.  Ingesting all 24 would put 20 rows into
db/binding_affinity.xlsx that the converter would turn into log K as
though they were binding constants, which they are not.  The full
table is therefore recorded verbatim in the note instead, which is
the record without the false rows.

Of the three mole ratios, 1:0.5 is taken: it is the lowest loading
and so the closest to the trace limit, and K_a,1 drifts with loading
(BSA-PFOA 8.3e5 -> 2.6e5 -> 2.1e5 across 1:0.5, 1:1, 1:2).  That
drift is recorded on the rows as a caveat; a mass-action constant
should not depend on the dose used to measure it.

  BSA-PFOA  K_a,1 8.3e5 M^-1  Kd 1.205 uM  log K 4.096
  BSA-PFOS  K_a,1 9.9e5       Kd 1.010     log K 4.173
  HSA-PFOA  K_a,1 1.3e5       Kd 7.692     log K 3.291
  HSA-PFOS  K_a,1 2.0e5       Kd 5.000     log K 3.478

n = 1 IS A LABEL HERE, NOT AN ASSUMPTION THAT MOVES THE NUMBER.
Splitting K_a,1 into n identical sites cancels out of the bridge:
log10(n) - log10(n x 1e6/K_a,1) = -log10(1e6/K_a,1).  Verified in
notebook/chi2018_stepwise_constants.py for n = 1, 4 and 8.  That is a
better footing than the n = 1 on Qin 2010 or Chen 2020, where the
number genuinely depends on an assumption nobody made.

CONDITIONS.  Albumin 5 uM, fatty-acid-free and IgG-free; at the 1:0.5
ratio the PFAS is 2.5 uM, so the LOADING IS 0.5 -- among the lowest in
this sheet, which is what makes these rows trace-like.  45 min at
37 C.  BUFFER IS 10 mM AMMONIUM ACETATE at pH 7.4, not phosphate:
ESI-MS needs a volatile buffer, so ionic strength is about a
fifteenth of blood's and the electrostatics are not physiological.
Recorded as a condition, not a defect.

WITHIN-STUDY SPECIES COMPARISON, one experiment, one method: BSA binds
6.4x tighter than HSA for PFOA (+0.81 log) and 5.0x for PFOS (+0.69).
Only three sources in the sheet compare species within a study
(Bischel 2010, Fischer 2024, Han 2003) and none measures both
compounds on both albumins in one run.

NOT INGESTED: the molecular docking (computational), and the
fluorescence strand, whose only number is a Stern-Volmer slope for
BSA-PFOA (F0/F = 1 + 7.781e4 [Q]) that the paper itself offers as
corroboration rather than a constant, with the PFOS plot explicitly
nonlinear.

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/add_chi2018_primary.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "Table 1 | read from the full text, primary 2026-11-13"

METHOD = ("native electrospray ionization mass spectrometry "
          "(ESI-Q-TOF), stepwise association constants from "
          "deconvolved PL_j peak abundances")

TABLE1 = (
    "TABLE 1 IN FULL, so the four ingested rows can be read against "
    "the twenty they come from (K_a, L/mol; columns BSA-PFOA, "
    "BSA-PFOS, HSA-PFOA, HSA-PFOS). 1:0.5 K_a,1 = 8.3e5, 9.9e5, "
    "1.3e5, 2.0e5. 1:1 K_a,1 = 2.6e5, 9.1e5, 1.2e5, 1.4e5; K_a,2 = "
    "1.5e5, 2.7e5, 6.3e4, 6.9e4. 1:2 K_a,1 = 2.1e5, 7.6e5, 1.1e5, "
    "1.1e5; K_a,2 = 1.1e5, 2.4e5, 6.1e4, 7.6e4; K_a,3 = 5.5e4, 1.2e5, "
    "5.2e4, 5.4e4. The stub's 'range' 0.52e5-9.90e5 is exactly the "
    "minimum and maximum of these 24.")

STEPWISE = (
    "ONLY K_a,1 IS INGESTED. The constants are STEPWISE: K_a,j = "
    "[PL_j]/([PL_{j-1}][L]), so K_a,2 is the affinity for adding a "
    "SECOND ligand to a protein that already holds one, not a second "
    "site's affinity. log K is a trace-limit quantity in which only "
    "the first binding event exists, so K_a,1 is the only one that "
    "answers it; K_a,2 and K_a,3 describe an already-loaded protein. "
    "Their decrease is also largely statistical rather than "
    "cooperative: for n identical independent sites K_a,j = "
    "k(n-j+1)/j, which for n = 8 already gives 1 : 0.44 : 0.25. "
    "See notebook/chi2018_stepwise_constants.py.")

DRIFT = (
    "CAVEAT, K_a,1 DRIFTS WITH LOADING: across the 1:0.5, 1:1 and 1:2 "
    "ratios it moves 8.3e5 -> 2.6e5 -> 2.1e5 for BSA-PFOA (4.0x), "
    "9.9e5 -> 9.1e5 -> 7.6e5 for BSA-PFOS (1.3x), 1.3e5 -> 1.2e5 -> "
    "1.1e5 for HSA-PFOA (1.2x) and 2.0e5 -> 1.4e5 -> 1.1e5 for "
    "HSA-PFOS (1.8x). A mass-action constant should not depend on the "
    "dose used to measure it. The paper remarks on the stepwise "
    "decrease WITHIN a ratio but not on K_a,1 moving BETWEEN ratios. "
    "The 1:0.5 column is taken here because it is the lowest loading "
    "and so the closest to the trace limit.")

NCANCEL = (
    "n = 1 IS A LABEL, NOT A LOAD-BEARING ASSUMPTION. Splitting K_a,1 "
    "into n identical sites cancels out of the bridge: log10(n) - "
    "log10(n x 1e6/K_a,1) = -log10(1e6/K_a,1). Asserted for n = 1, 4 "
    "and 8 in notebook/chi2018_stepwise_constants.py. At trace the "
    "protein's capacity is n x k and K_a,1 already is that product.")

COND = (
    "Albumin 5 uM (BSA 66,430 Da and HSA 66,550 Da as the paper "
    "states), fatty-acid-free and IgG-free, Sigma-Aldrich. PFOA/PFOS "
    "dosed at 0, 2.5, 5, 10, 20, 30 and 40 uM over the whole study, "
    "i.e. loadings 0 to 8; THESE ROWS are the 1:0.5 spectrum, where "
    "the PFAS is 2.5 uM and the loading is 0.5. Incubated 45 min at "
    "37 C, injected into an Agilent 6520 ESI-Q-TOF. BUFFER IS 10 mM "
    "AMMONIUM ACETATE at pH 7.4, NOT phosphate: ESI-MS needs a "
    "volatile buffer, so the ionic strength is about a fifteenth of "
    "blood's and the electrostatics are not physiological. Up to 8 "
    "ligands are resolved per albumin at the highest ratios.")

COMMON = {
    "study": "Chi et al.", "year": "2018", "pmid": "29425944",
    "doi": "10.1016/j.chemosphere.2018.01.152",
    "method": METHOD, "ka_units": "M^-1 (stepwise K_a,1, Table 1)",
    "kd_units": "uM", "protein_conc": "5", "protein_conc_units": "uM",
    "ligand_conc_range": "2.5 (the 1:0.5 spectrum; 0 to 40 overall)",
    "ligand_conc_units": "uM", "ligand_over_protein_max": "0.5",
    "temperature_c": "37", "buffer": "10 mM ammonium acetate, pH 7.4",
    "is_computational": ("no (AutoDock molecular docking is reported "
                         "alongside in Section 2.5, as a separate "
                         "computational strand)"),
    "tk_relevance_noted": ("partly - the paper frames its aim as the "
                           "toxic mechanism of PFCAs in vivo but draws "
                           "no kinetic conclusion"),
    "n_binding_sites": ("1 (the first stepwise constant K_a,1 from "
                        "native ESI-MS; up to 8 ligands bind overall, "
                        "but K_a,1 is the first binding event)"),
    "binding_site_identified": ("hydrophobic pocket proximate to "
                                "Trp214 in HSA is one dominant site, "
                                "from warfarin/ibuprofen displacement "
                                "and docking; more than one pocket is "
                                "occupied"),
}

ROWS = [("PFOA", "bovine", "bovine serum albumin (BSA)", "7",
         "carboxylate", 8.3e5),
        ("PFOS", "bovine", "bovine serum albumin (BSA)", "8",
         "sulfonate", 9.9e5),
        ("PFOA", "human", "human serum albumin (HSA)", "7",
         "carboxylate", 1.3e5),
        ("PFOS", "human", "human serum albumin (HSA)", "8",
         "sulfonate", 2.0e5)]


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    if any(MARK in (r.get("source_table") or "")
           and "Chi" in (r.get("study") or "") for r in rows):
        print("  already ingested; nothing to do")
        return 0

    before = len(rows)
    rows = [r for r in rows
            if not ("Chi" in (r.get("study") or "")
                    and r.get("year") == "2018")]
    removed = before - len(rows)

    note = " | ".join([COND, STEPWISE, DRIFT, NCANCEL, TABLE1])
    added = []
    for chem, sp, prot, ncarb, head, ka in ROWS:
        added.append(dict(
            COMMON, chemical=chem, protein_species=sp, protein=prot,
            n_perfluorinated_carbons=ncarb, head_group=head,
            ka="%g" % ka, kd="%.4g" % (1e6 / ka),
            source_table=(f"Table 1, {prot[:3]}-{chem} K_a,1 at the "
                          f"1:0.5 mole ratio | " + MARK),
            assay_volume_note=note))
    for spec in added:
        rows.append({f: spec.get(f, "") for f in fields})

    with open(DB, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  removed {removed} stub row, added {len(added)}:")
    for s in added:
        print(f"    {s['protein'][:3]}-{s['chemical']:5s} "
              f"K_a,1 {s['ka']:>8s} M^-1  Kd {s['kd']:>7s} uM  "
              f"loading 0.5")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
