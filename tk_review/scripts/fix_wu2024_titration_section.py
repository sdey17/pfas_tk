#!/usr/bin/env python3
"""Wu 2024: the conditions came from the wrong section of the paper.

Wu Y, Bao J, Liu Y, Wang X, Lu X, Wang K.  "In Vitro and In Silico
Analysis of the Bindings between Legacy and Novel Per- and
Polyfluoroalkyl Substances and Human Serum Albumin."  Toxics
2024;12(1):46.  doi 10.3390/toxics12010046, PMID 38251003, PMC10818824.
Open access; read at source.

THE PAPER HAS TWO MIXING SCHEMES AND THEY ARE NOT THE SAME EXPERIMENT.

  Section 2.2, three-dimensional fluorescence, mixes 5 mL of HSA with
  1, 3 or 5 mL of PFAS -- and says why: "to determine the optimal
  concentration for the instrumental analysis".  It is a SCREENING
  step that picks a ratio for the 3-D spectra.  No binding constant
  comes out of it.

  Section 2.3, two-dimensional fluorescence, is the titration that
  produces every constant in Tables 4-6: "mixtures were configured
  with HSA:PFASs concentration ratios of 5:0, 5:1, 5:5, 5:10, 5:15,
  5:20, and 5:30".

The database recorded Section 2.2.  Its top point, 5:5, gives a
ligand:protein ratio of 1.00, and that is the 1.00 every Wu 2024 row
was carrying.  The titration that actually produced the constants
reaches 5:30 -- a ratio of SIX.

  v_HSA:v_PFAS   [HSA] uM   [PFAS] uM   molar ratio   ethanol v/v
       5:0          1000           0          0.0            0 %
       5:1           833         167          0.2           17 %
       5:5           500         500          1.0           50 %
      5:10           333         667          2.0           67 %
      5:15           250         750          3.0           75 %
      5:20           200         800          4.0           80 %
      5:30           143         857          6.0           86 %

THE RATIO OF 6 DOES NOT DEPEND ON THE STOCK CONCENTRATION, because
both stocks are made up at the same molarity and the ratio is a ratio
of volumes.  The absolute concentrations do, and the PMC text
extraction loses the exponent in "1 x mol/L" for both stocks; 1 mM is
what the database already assumed and what the printed micromolar
figures here rest on.  The ratio is safe; treat the absolute
concentrations as resting on that reading.

THE ETHANOL IS THE HEADLINE.  Every PFAS stock is in ANHYDROUS
ETHANOL, so the ethanol fraction is the PFAS volume fraction.  The
database said 17-50 %, from Section 2.2.  The real titration ends at
**86 % v/v ethanol** -- the top point of every binding curve in the
paper is albumin in more ethanol than water.  Albumin work normally
keeps cosolvent below 1 %.

That matters because audit_wu2024_hsa.py already shows Wu's constants
sitting 2.7-3.8 logs below every other HSA study with no
compound-specific structure, which is the signature of a method
rather than of chemistry.  An 86 % ethanol endpoint is a mechanism
for exactly that, and it was invisible while the conditions were
being read off the screening section.

HSA IS NOT HELD CONSTANT either: it falls 5.8-fold across the series,
833 down to 143 uM, so the "quenching versus PFAS concentration" plots
vary both axes at once.

This script fixes the 20 database rows and the MIX constant in
audit_wu2024_hsa.py, whose C* moves from 315 uM (the Section 2.2
geometric mean) to 553 uM.  No binding constant is changed: the
constants are the paper's and they stay as reported.

Run:  python3 tk_review/scripts/fix_wu2024_titration_section.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
AUDIT = os.path.join(HERE, "audit_wu2024_hsa.py")

STOCK_UM = 1000.0                                  # see docstring
TITRATION = [(5, 0), (5, 1), (5, 5), (5, 10), (5, 15), (5, 20), (5, 30)]
WITH_LIGAND = [(h, p) for h, p in TITRATION if p]

HSA = [STOCK_UM * h / (h + p) for h, p in WITH_LIGAND]
PFAS = [STOCK_UM * p / (h + p) for h, p in WITH_LIGAND]
RATIO = max(PFAS) / min(HSA)
ETOH = max(p / (h + p) for h, p in WITH_LIGAND)

NOTE = (
    f"CONDITIONS ARE FROM SECTION 2.3, the two-dimensional fluorescence "
    f"titration that produces the constants in Tables 4-6: HSA:PFAS "
    f"concentration ratios 5:0, 5:1, 5:5, 5:10, 5:15, 5:20, 5:30. An "
    f"earlier version of this row recorded SECTION 2.2 instead, which is "
    f"the three-dimensional screening step the paper runs 'to determine "
    f"the optimal concentration for the instrumental analysis' and which "
    f"yields no binding constant; its top point 5:5 is where the old "
    f"ligand:protein of 1.00 came from. The real titration reaches "
    f"{RATIO:.0f}. Both stocks are at the same molarity so the ratio is a "
    f"volume ratio and does not depend on it; the absolute micromolar "
    f"values assume the 1 mM stock reading, whose exponent the PMC text "
    f"loses. ETHANOL: every PFAS stock is in ANHYDROUS ETHANOL, so the "
    f"cosolvent fraction equals the PFAS volume fraction and the top "
    f"point of every binding curve is {ETOH:.0%} v/v ethanol -- albumin "
    f"work normally stays below 1%. HSA is NOT held constant: it falls "
    f"{max(HSA)/min(HSA):.1f}-fold across the series, {max(HSA):.0f} to "
    f"{min(HSA):.0f} uM, so both axes of the quenching plots vary at "
    f"once. protein_conc stores {min(HSA):.0f}, the HSA concentration at "
    f"the 5:30 point where the ratio is maximal, so that the two "
    f"concentration columns divide to the ratio column. See audit_wu2024_hsa.py: these constants sit 2.7-3.8 logs "
    f"below every other HSA study with no compound-specific structure")
OLD_MARKS = ("HSA 1 mM stock mixed 5:1, 5:3, 5:5",
             "CONDITIONS ARE FROM SECTION 2.3")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"].startswith("Wu") and r["year"] == "2024"]
    print(f"  {len(hit)} Wu 2024 rows\n")
    print(f"  {'field':24}{'was':>26}{'now':>26}")
    was = (hit[0]["protein_conc"], hit[0]["ligand_conc_range"],
           hit[0]["ligand_over_protein_max"])
    # Store the HSA concentration at the point where the ratio is
    # maximal (the 5:30 end), not the range, so the two concentration
    # columns divide to the ratio column. Storing "143-833" made the
    # workbook show 833 beside a ratio of 6, which reads as broken.
    # The full range is in the note.
    new_prot = f"{min(HSA):.0f}"
    new_lig = f"0-{max(PFAS):.0f}"
    new_ratio = f"{RATIO:.2f}"
    print(f"  {'protein_conc':24}{was[0] + ' mM':>26}"
          f"{new_prot + ' uM (at the 5:30 point)':>26}")
    print(f"  {'ligand_conc_range':24}{was[1] + ' mM':>26}{new_lig + ' uM':>26}")
    print(f"  {'ligand_over_protein_max':24}{was[2]:>26}{new_ratio:>26}")
    print(f"  {'ethanol at top point':24}{'50% (section 2.2)':>26}"
          f"{f'{ETOH:.0%} (section 2.3)':>26}")

    n = 0
    for r in hit:
        r["protein_conc"] = new_prot
        r["protein_conc_units"] = "uM"
        r["ligand_conc_range"] = new_lig
        r["ligand_conc_units"] = "uM"
        r["ligand_over_protein_max"] = new_ratio
        parts = [x for x in r["assay_volume_note"].split(" | ")
                 if not any(m in x for m in OLD_MARKS)]
        r["assay_volume_note"] = " | ".join(parts + [NOTE])
        n += 1
    print(f"\n  {n} rows updated; no binding constant changed")

    # --- the audit script reads the same wrong section -------------
    src = open(AUDIT, encoding="utf-8").read()
    old_mix = "MIX = [(5, 1), (5, 3), (5, 5)]                      # mL HSA : mL PFAS"
    new_mix = ("MIX = [(5, 1), (5, 5), (5, 10), (5, 15), (5, 20), (5, 30)]\n"
               "# Section 2.3, the two-dimensional titration that produces the\n"
               "# constants -- NOT section 2.2's 5:1/5:3/5:5, which is the\n"
               "# screening step the paper runs to pick a ratio for its 3-D\n"
               "# spectra and from which no binding constant is derived.")
    cstar_old = math.exp(sum(math.log(p / (h + p) * 1e-3)
                             for h, p in [(5, 1), (5, 3), (5, 5)]) / 3)
    cstar_new = math.exp(sum(math.log(p / (h + p) * 1e-3)
                             for h, p in WITH_LIGAND) / len(WITH_LIGAND))
    if old_mix in src:
        src = src.replace(old_mix, new_mix)
        src = src.replace('print("THE ASSAY, as mixed (Section 2.2)")',
                          'print("THE ASSAY, as mixed (Section 2.3, the '
                          'titration)")')
        print(f"  audit_wu2024_hsa.py: MIX moved to section 2.3; "
              f"C* {cstar_old*1e6:.0f} -> {cstar_new*1e6:.0f} uM")
    else:
        print(f"  audit_wu2024_hsa.py: already on section 2.3 "
              f"(C* {cstar_new*1e6:.0f} uM)")

    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader(); w.writerows(rows)
    open(AUDIT, "w", encoding="utf-8").write(src)
    print(f"  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))} "
          f"and {os.path.relpath(AUDIT, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
