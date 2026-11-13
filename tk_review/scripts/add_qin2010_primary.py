#!/usr/bin/env python3
"""Replace the Qin 2010 secondary stub with its six measured rows.

THE PAPER.  Qin P, Liu R, Pan X, Fang X, Mou Y.  "Impact of Carbon
Chain Length on Binding of Perfluoroalkyl Acids to Bovine Serum
Albumin Determined by Spectroscopic Methods."  J Agric Food Chem
2010;58(9):5561-5567.  doi 10.1021/jf100412q, PMID 20397730.

THE STUB said chemical "PFOA and PFDA (PFPeA tested, little effect)"
with ka "0.22e5-6.85e5 M^-1" from a secondary, took the geometric
midpoint as Kd 8.146 uM and assumed n = 1, giving log K 3.266.  That
apparent range is TWO COMPOUNDS AND THREE TEMPERATURES collapsed:
0.22e5 is PFOA at 320 K (Table 2: 2.17e4), 6.85e5 is PFDA at 300 K
(6.846e5).  The midpoint belongs to no compound at any temperature.

WHAT IS INGESTED.  Table 2's six modified Stern-Volmer constants, all
verified against their own Fig. 3 fits in
scripts/audit_qin2010_modified_sv.py.  Only the 310 K pair reaches the
sheet, via BY_TEMPERATURE in build_revised_working_table.py; the other
four stay here as the record, as Peng 2024's non-310 K rows do.

WHERE THE CONCENTRATIONS COME FROM, and they cross-check:

  BSA 1 uM    Fig. 1's caption states it outright, "The concentration
              of BSA: 1 x 10^-6 mol/L", and Figs. 2 and 3 repeat it.
              It also falls out of the method: Reagents prepares BSA
              at 1.0 x 10^-5 mol/L, and Fluorescence Quenching
              Measurements adds "1 mL of BSA stock solution ... to a
              10 mL colorimetric tube and made up to the mark", so
              1.0e-5 x 1/10 = 1.0e-6 mol/L.  Two independent routes to
              the same number.

  PFOA        Fig. 1's caption lists the dosed series verbatim:
  0-20 uM     "PFOA (x10^-6 mol/L) 1-10: 0, 1, 2, 3, 4, 5, 6, 7, 10,
              and 20", so the top point is 20 uM and the maximum
              loading is 20.

  PFDA        same caption: "PFDA (x10^-6 mol/L) 1-10: 0, 0.25, 0.5,
  0-7 uM      0.75, 1, 1.5, 2, 3, 5, and 7", top 7 uM, loading 7.

  buffer      Reagents and Methods: 0.02 mol/L NaH2PO4-Na2HPO4 at
              pH 7.40, which is also in the captions of Figs. 1-3.

  PFPA        "PFPA (x10^-6 mol/L) 1-6: 0, 1, 3, 5, 7, and 10" was
              tested to 10 uM and gave NO constant -- "PFPA had little
              effect on BSA".  A result, not a gap, so no row.

NOTE ON THE PROTEIN.  BSA is "electrophoretic grade reagent, 90%
purity", so the nominal 1 uM is about 0.9 uM of albumin.  Recorded
rather than corrected.

THE DEFECT, carried on the PFOA rows.  Eq 2 is F0/dF =
1/(f_a K_a [Q]) + 1/f_a, so Fig. 3's intercept is 1/f_a.  PFOA's three
intercepts give f_a = 1.190, 1.449 and 1.761, and f_a is "the fraction
of accessible fluorescence" in the paper's own words -- a fraction
cannot exceed 1, and BSA has two tryptophans, Trp-134 and Trp-212,
both named in this paper.  PFDA's 0.437, 0.551 and 0.730 are physical.
This is the test that held out Wu 2024, but the degree differs: Wu's
constants were ALSO 10^1-10^2 M^-1, three to four logs below every
other albumin study, whereas these are 10^4-10^5 M^-1.  Reported as
published with the defect recorded, as Peng 2024 is.

K_a is not a binding constant in any case: the paper calls it "the
effective quenching constant for the accessible fluorophores".  And
the paper reports NO stoichiometry anywhere -- Stern-Volmer in Table 1
and modified Stern-Volmer in Table 2, no double-logarithmic plot -- so
the sheet's "n = 1 assumed by the reader, the paper reports none" is,
unusually, true here.

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/add_qin2010_primary.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "Table 2 | read from the full text, primary 2026-11-13"

METHOD = ("fluorescence quenching, modified Stern-Volmer plot "
          "F0/dF vs 1/[Q]: intercept 1/f_a, slope 1/(f_a K_a)")

CONC = ("CONCENTRATIONS, from Fig. 1's caption, which states the dosed "
        "series verbatim: 'PFOA (x10^-6 mol/L) 1-10: 0, 1, 2, 3, 4, 5, "
        "6, 7, 10, and 20. PFDA (x10^-6 mol/L) 1-10: 0, 0.25, 0.5, "
        "0.75, 1, 1.5, 2, 3, 5, and 7. The concentration of BSA: 1 x "
        "10^-6 mol/L'. Figs. 2 and 3 repeat the BSA concentration. It "
        "cross-checks against the method: Reagents prepares BSA at 1.0 "
        "x 10^-5 mol/L and Fluorescence Quenching Measurements adds '1 "
        "mL of BSA stock solution ... to a 10 mL colorimetric tube and "
        "made up to the mark', giving 1.0e-5 x 1/10 = 1.0e-6 mol/L -- "
        "two independent routes to the same 1 uM. BSA is "
        "'electrophoretic grade reagent, 90% purity', so the nominal 1 "
        "uM is about 0.9 uM of albumin; recorded, not corrected. PFPA "
        "(C5) was dosed to 10 uM and gave NO constant: 'PFPA had "
        "little effect on BSA', a result rather than a gap.")

FA = ("ACCESSIBLE-FLUOROPHORE FRACTION ABOVE 1. Eq 2 is F0/dF = "
      "1/(f_a K_a [Q]) + 1/f_a, so Fig. 3's intercept is 1/f_a. PFOA's "
      "three intercepts (0.840, 0.690, 0.568 at 300/310/320 K) give "
      "f_a = 1.190, 1.449 and 1.761. f_a is 'the fraction of "
      "accessible fluorescence' in the paper's own words and cannot "
      "exceed 1; BSA has two tryptophans, Trp-134 and Trp-212, both "
      "named here. The excess grows with temperature. PFDA's "
      "intercepts (2.289, 1.814, 1.370) give 0.437, 0.551 and 0.730, "
      "which are physical. Same test that held out Wu 2024, but Wu's "
      "constants were ALSO 10^1-10^2 M^-1; these are 10^4-10^5 M^-1, "
      "an ordinary range, so PFOA fails one check here and not two. "
      "Full working: scripts/audit_qin2010_modified_sv.py")

KA = ("K_a is the paper's 'effective quenching constant for the "
      "accessible fluorophores' from the MODIFIED Stern-Volmer "
      "equation, not a site binding constant. Table 1's K_sv is the "
      "ordinary Stern-Volmer quenching constant and is not a binding "
      "constant either. The paper reports NO stoichiometry anywhere: "
      "there is no double-logarithmic plot in it, so n = 1 is the "
      "reader's and 'the paper reports none' is true.")

COMMON = {
    "study": "Qin et al.", "year": "2010", "pmid": "20397730",
    "doi": "10.1021/jf100412q",
    "protein": "bovine serum albumin (BSA)", "protein_species": "bovine",
    "method": METHOD, "ka_units": "M^-1 (modified Stern-Volmer K_a)",
    "kd_units": "uM", "n_binding_sites": "",
    "protein_conc": "1", "protein_conc_units": "uM",
    "ligand_conc_units": "uM", "head_group": "carboxylate",
    "buffer": "0.02 mol/L NaH2PO4-Na2HPO4, pH 7.40",
    "is_computational": "no",
    "tk_relevance_noted": ("partly - the paper's framing is toxicity at "
                           "the protein level and chain-length trend, "
                           "not toxicokinetics"),
    "binding_site_identified": ("alpha-helix content reduced (CD); van "
                                "der Waals forces and hydrogen bonds "
                                "dominant from dH < 0 and dS < 0"),
}

SPEC = [("PFOA", "7", "0 to 20", "20", [(300, 4.36e4), (310, 2.93e4),
                                        (320, 2.17e4)], True),
        ("PFDA", "9", "0 to 7", "7", [(300, 6.846e5), (310, 4.417e5),
                                      (320, 2.722e5)], False)]


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())
    if any(MARK in (r.get("source_table") or "") for r in rows):
        print("  already ingested; nothing to do")
        return 0

    before = len(rows)
    rows = [r for r in rows
            if not ("Qin" in (r.get("study") or "")
                    and r.get("year") == "2010")]
    removed = before - len(rows)

    added = []
    for chem, ncarb, lrange, lmax, series, flag in SPEC:
        for T, ka in series:
            note = CONC + " | " + KA + (" | " + FA if flag else "")
            added.append(dict(
                COMMON, chemical=chem, n_perfluorinated_carbons=ncarb,
                ka="%g" % ka, kd="%.4g" % (1e6 / ka),
                ligand_conc_range=lrange, ligand_over_protein_max=lmax,
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
              f"K_a {s['ka']:>9s} M^-1  Kd {s['kd']:>7s} uM")
    print("  only the 310 K pair reaches the sheet (BY_TEMPERATURE)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
