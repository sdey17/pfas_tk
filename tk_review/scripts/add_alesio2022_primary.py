#!/usr/bin/env python3
"""Replace the Alesio 2022 narrative stub with its Table 1, from full text.

Alesio JL, Slitt A, Bothun GD.  "Critical new insights into the binding of
poly- and perfluoroalkyl substances (PFAS) to albumin protein."
Chemosphere 2021;287:131979.  doi 10.1016/j.chemosphere.2021.131979.
Table 1, "Equilibrium dialysis results for PFAS-BSA binding at 37 C".

WHAT WAS THERE.  One narrative row carrying a Ka RANGE, 0.07e6-6.16e6
M^-1, taken from Starnes & Belcher 2026's compilation, with no compound
resolution, no conditions, and `method` describing the paper's
fluorescence-quenching work.  That range is in fact Table 1's EQUILIBRIUM
DIALYSIS Ka column rounded, so the stub was both secondary and mislabelled.
Seven primary rows replace it.

IT IS BOVINE ALBUMIN.  The paper's title says "albumin protein" and the
abstract, Methods and Table 1 all say bovine serum albumin, lyophilised,
99% fatty-acid free, Sigma-Aldrich.  There is no HSA in it.  Worth stating
because the title does not.

WHAT IS INGESTED, AND WHAT IS NOT.  Table 1 has four numeric columns:

  f        bound fraction            -> INGESTED as unbound_fraction
  K_PW     partition coefficient     -> NOT ingested, see below
  Ka       association constant      -> INGESTED, and Kd = 1/Ka
  Ka (b)   "Reported by Allendorf et al. (2019)" -- a DIFFERENT STUDY,
           reprinted for comparison.  NOT ingested: Allendorf 2019
           already holds 28 albumin rows here and this would double-count
           it under Alesio's name.

The K_PW column is left out because the table does not reproduce itself.
The paper's eq 1 derives K_PW from f and the BSA concentration, so f, K_PW
and Ka are three views of one experiment and must agree.  Two of them do:
going from Ka to a partition coefficient by this project's conversion, and
going from f through eq 1 with the stated 10 uM BSA, land within a median
0.085 log of each other.  The printed K_PW column sits a consistent 10.3x
(9.5-10.7x over seven compounds) above both -- a factor-of-ten slip, not a
scientific disagreement.  So Ka and f are ingested, the derived column is
not, and the note on every row records the discrepancy.

(The two compounds where the two routes agree least, PFBS at 0.45 log and
GenX at 0.26, are exactly the two run at a 2:1 PFAS:BSA ratio instead of
0.1:1.  Eq 2 assumes 1:1 binding, which strains hardest at the highest
loading.  The pattern is coherent, which is part of why the 10x is read as
a units slip in one column rather than noise.)

n = 1 IS THE PAPER'S OWN MODEL, not a reader's assumption: eq 2 states
that ([P] - [PFAS]bound) "represents the residual unbound BSA protein
assuming 1:1 binding".  The rows therefore carry n = 1 with that basis,
and the downstream classifier splits them on whether Kd falls inside the
concentration actually tested -- for PFOA and PFHxS it does not.

THE LOADING IS LOW BUT NOT THE LOWEST.  Five of seven compounds were run
at 1 uM PFAS against 10 uM BSA, a ligand:protein ratio of 0.1, and the
other two at 2.0.  That is well below most of this compilation, but
Allendorf 2019 and Bischel 2010/2011 reach 0.0004 to 0.0015 -- two to
three orders lower again, and the only rows anywhere near the
physiological nu < 0.01.

Idempotent: re-running changes nothing.

Run:  python3 -I tk_review/scripts/add_alesio2022_primary.py
"""
import csv
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Alesio et al.", "2022"
DOI = "10.1016/j.chemosphere.2021.131979"
PMID = "34450356"

METHOD = ("equilibrium dialysis (Thermo RED rapid-equilibrium device, "
          "4 h at 37 C under constant shaking, 1 h pre-incubation) with "
          "LC/MS of the buffer side; Ka from the paper's eq 2, which "
          "assumes 1:1 binding")
BUFFER = "phosphate-buffered saline (PBS), pH 7.4"
NOTE = ("Ka is Table 1's equilibrium-dialysis column, NOT the paper's "
        "fluorescence-quenching work (Stern-Volmer / modified "
        "Stern-Volmer / Hill, reported separately). Table 1's printed "
        "K_PW column is NOT used: recomputing it from the paper's own eq "
        "1 with the stated 10 uM BSA puts it 10.3x (9.5-10.7x over the "
        "seven compounds) below what is printed, while the Ka route and "
        "the eq 1 route agree with each other to a median 0.085 log. "
        "Table 1's fourth column is Allendorf et al. 2019's Ka, reprinted "
        "for comparison, and is NOT this study's measurement. BSA was "
        "lyophilised, 99% fatty-acid free (Sigma-Aldrich)")

# (chemical, nC, head group, f, Ka x1e5 M^-1, PFAS:BSA ratio)
T1 = [("PFOA",           7, "carboxylate", 0.755,  3.77, 0.1),
      ("PFNA",           8, "carboxylate", 0.975, 44.05, 0.1),
      ("PFDA",           9, "carboxylate", 0.983, 61.62, 0.1),
      ("PFBS",           4, "sulfonate",   0.249,  0.93, 2.0),
      ("PFHxS",          6, "sulfonate",   0.842,  6.48, 0.1),
      ("PFOS",           8, "sulfonate",   0.942, 17.99, 0.1),
      ("HFPO-DA (GenX)", 5, "ether carboxylate", 0.265, 0.65, 2.0)]
SE = {"PFOA": 0.80, "PFNA": 8.37, "PFDA": 15.72, "PFBS": 0.55,
      "PFHxS": 0.39, "PFOS": 1.44, "HFPO-DA (GenX)": 0.04}
F_SE = {"PFOA": 0.038, "PFNA": 0.005, "PFDA": 0.004, "PFBS": 0.064,
        "PFHxS": 0.008, "PFOS": 0.004, "HFPO-DA (GenX)": 0.006}
P_uM = 10.0


def sig(x, n=4):
    s = f"{x:.{n}g}"
    return s


def build(fields):
    out = []
    for chem, nc, head, f, ka5, ratio in T1:
        ka = ka5 * 1e5
        kd_uM = 1e6 / ka
        lig = P_uM * ratio
        r = {k: "" for k in fields}
        r.update({
            "study": STUDY, "year": YEAR, "pmid": PMID, "doi": DOI,
            "chemical": chem, "n_perfluorinated_carbons": str(nc),
            "head_group": head,
            "protein": "bovine serum albumin (BSA)",
            "protein_species": "bovine",
            "method": METHOD,
            "ka": f"{ka:.3g}",
            "ka_units": f"M^-1 (SE {SE[chem] * 1e5:.3g})",
            "kd": sig(kd_uM), "kd_units": "uM, computed as 1/Ka",
            "n_binding_sites": "1",
            "binding_site_identified": "",
            "unbound_fraction": (f"{(1 - f) * 100:.1f}% unbound "
                                 f"({f * 100:.1f}% bound, f = {f:.3f} "
                                 f"+/- {F_SE[chem]:.3f} SE, Table 1)"),
            "protein_conc": f"{P_uM:g}", "protein_conc_units": "uM",
            "ligand_conc_range": f"{lig:g}", "ligand_conc_units": "uM",
            "temperature_c": "37", "buffer": BUFFER,
            "is_computational": "no", "tk_relevance_noted": "yes",
            "source_table": (f"Table 1, equilibrium dialysis, Ka column | "
                             f"n = 1 is the paper's own eq 2 assumption, "
                             f"not a reader's | PFAS:BSA molar ratio "
                             f"{ratio:g}:1"),
            "assay_volume_note": NOTE,
            "ligand_over_protein_max": f"{ratio:g}",
        })
        out.append(r)
    return out


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        fields, rows = rd.fieldnames, list(rd)

    new = build(fields)
    keep = [r for r in rows
            if not (r["study"] == STUDY and r["year"] == YEAR)]
    dropped = len(rows) - len(keep)
    if dropped:
        print(f"  removing {dropped} existing '{STUDY} {YEAR}' row(s)")
        for r in rows:
            if r["study"] == STUDY and r["year"] == YEAR:
                print(f"    was: {r['chemical'][:58]}")
                print(f"         {r['source_table'][:68]}")

    have = {(r["study"], r["year"], r["chemical"]) for r in keep}
    if all((STUDY, YEAR, r["chemical"]) in have for r in new) and not dropped:
        print("  already present and unchanged; nothing to do")
        return

    shutil.copy(DB, DB + ".bak")
    at = max((i for i, r in enumerate(keep)
              if r["study"].startswith("Alesio")), default=len(keep) - 1) + 1
    out = keep[:at] + new + keep[at:]
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out)

    print(f"  added {len(new)} primary rows; {len(rows)} -> {len(out)}\n")
    print("  %-16s %7s %10s %9s %8s %7s" % (
        "chemical", "f", "Ka M^-1", "Kd uM", "lig uM", "ratio"))
    for r in new:
        print("  %-16s %7s %10s %9s %8s %7s" % (
            r["chemical"], r["unbound_fraction"].split("(")[1].split("%")[0],
            r["ka"], r["kd"], r["ligand_conc_range"],
            r["ligand_over_protein_max"]))
    print("\n  species: bovine (BSA) -- the paper has no human albumin")
    print("  five of seven at a ligand:protein ratio of 0.1, two at 2.0;")
    print("  low for this compilation, but Allendorf 2019 and Bischel")
    print("  2010/2011 reach 0.0004-0.0015, which is where the")
    print("  physiological limit actually sits")


if __name__ == "__main__":
    main()
