#!/usr/bin/env python3
"""Remove a dosing concentration Starnes 2024 never states.

The Starnes 2024 rows carried ligand_conc_range "0-1000" uM and a
ligand_over_protein_max derived from it, giving every row a PFAS
concentration of 1000 uM and a loading of 5.56 (BSA) or 7.69 (HSA).

1000 uM appears EXACTLY ONCE in the paper, in the caption of Figure 4:
"(C-H) Melt curves for BSA at individual concentrations of PFOS from 0
to 1000 uM."  That is one panel, one compound, one albumin.  It was
generalised to all eight compounds and all four species.

The paper states no general dosing series.  What it does state:

  Table 1   stock concentrations, 10 or 20 mM depending on solubility
  Methods   "PFAS stock solutions were serially diluted into HBS"
  Methods   positive controls at OA 1, 2, 3 mM and PFHxA 0.3, 1, 3 mM

and nothing else.  The series is compound-specific and only ever shown
graphically: Figure 3's x axes for PFBS reach about log[Q] = -2.2,
i.e. roughly 6 mM, and Figure 2's octanoic acid panel is captioned
0-6 mM.

The arithmetic also refuses 1000 uM.  Table 2's Kd run to 3.06 mM, and
a Kd of 3.06 mM cannot be fitted from a series topping out at 1 mM --
the curve would not reach half saturation.  So the stored value is not
merely unsourced, it is inconsistent with the constants on the same
rows.

WHAT THIS SCRIPT DOES.  Clears ligand_conc_range and
ligand_over_protein_max so the sheet's PFAS-conc and ratio cells go
blank, and records in the note what the paper does and does not say,
including the one stated figure for PFOS on BSA.  The Kd, the protein
concentrations (0.13 mM HSA, 0.18 mM BSA, both stated in Methods) and
everything else are untouched.

A blank is the right cell here.  The alternative, keeping a number the
paper never gives, is the failure this project has now found in four
other sources: a value present in the sheet that cannot be pointed at
in the primary.

Idempotent: keyed on MARK, safe to re-run.

Run:  python3 -I tk_review/scripts/fix_starnes2024_ligand_conc.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
MARK = "ligand series corrected 2026-11-14"

NOTE = (
    "LIGAND CONCENTRATIONS ARE NOT STATED NUMERICALLY for this assay, "
    "so PFAS conc and the loading ratio are left blank (" + MARK + "). "
    "The paper gives stock concentrations only -- Table 1, 10 or 20 mM "
    "depending on solubility -- and says 'PFAS stock solutions were "
    "serially diluted into HBS' without giving the series. The one "
    "numeric dosing range anywhere in the paper is Figure 4's caption, "
    "'Melt curves for BSA at individual concentrations of PFOS from 0 "
    "to 1000 uM', which is one panel, one compound and one albumin; an "
    "earlier version of these rows generalised that 1000 uM to all "
    "eight compounds and all four species. The series is otherwise "
    "shown only graphically: Figure 3's PFBS x axes reach about "
    "log[Q] = -2.2, roughly 6 mM, and Figure 2's octanoic acid panel is "
    "captioned 0-6 mM. The arithmetic refuses 1000 uM independently: "
    "Table 2's Kd run to 3.06 mM and a Kd of 3.06 mM cannot be fitted "
    "from a series topping out at 1 mM. Protein concentrations ARE "
    "stated and are kept: 0.13 mM HSA, 0.18 mM BSA, 0.11 mM PSA, "
    "0.18 mM RSA, each optimised per species for signal-to-noise.")


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fields = list(rows[0].keys())

    hit = 0
    for r in rows:
        if "Starnes" not in (r.get("study") or "") or r.get("year") != "2024":
            continue
        if MARK in (r.get("assay_volume_note") or ""):
            print("  already corrected; nothing to do")
            return 0
        r["ligand_conc_range"] = ""
        r["ligand_over_protein_max"] = ""
        base = (r.get("assay_volume_note") or "").strip(" |")
        r["assay_volume_note"] = (base + " | " + NOTE) if base else NOTE
        hit += 1

    if not hit:
        print("  no Starnes 2024 rows found")
        return 1
    with open(DB, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  cleared the unsourced ligand series on {hit} rows")
    print("  PFAS conc and ratio now blank; Kd and protein conc kept")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
