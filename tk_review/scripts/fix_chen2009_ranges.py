#!/usr/bin/env python3
"""Chen & Guo 2009: ligand ranges that the stock solutions make impossible.

Chen YM, Guo LH.  "Fluorescence study on site-specific binding of
perfluoroalkyl acids to human serum albumin."  Arch Toxicol
2009;83(3):255-261.  doi 10.1007/s00204-008-0359-x, PMID 18854981.

WHAT IS WRONG.  Six of the fifteen PFAS rows carry a ligand range
read off a figure's x-axis rather than off the compound's own series.
Figures 2, 4 and 5 each plot all five PFAAs on a SHARED axis running
0-6 mM, so "0-6" is a property of the panel, not of the compound --
and for two compounds it is not even reachable.

THE CEILING NOBODY CAN EXCEED.  The Experimental section gives a
stock concentration per compound:

    PFBA 20 mM   PFBS 5 mM   PFOA 2.5 mM   PFDoA 10 mM
    PFOS 200 mM in DMSO, diluted to 20 mM in buffer

A titration adds stock to a protein solution, so the cuvette
concentration approaches the stock from below and can never pass it.
PFBS cannot reach 6 mM from a 5 mM stock.  PFOA cannot reach 6 mM
from a 2.5 mM stock -- it is 2.4x over the ceiling.  Those are not
readings to be refined; they are arithmetically impossible, and the
ligand:protein ratio built on them (60) is impossible with them.

WHAT REPLACES THEM, in order of preference:

  1. a concentration the TEXT states          -> use it
  2. a figure caption's explicit POINT LIST   -> use its top
  3. the SHARED 0-6 mM axis                   -> upper bound only,
     and only where the stock allows it
  4. the STOCK concentration                  -> where (3) does not

That gives, per compound:

    PFOS   4 mM    Fig 1 lists 0/0.4/0.7/1.5/2/4; Fig 3 lists
                   0/0.2/0.7/2/4 and the text says "from 0.7 to
                   4 mM".  The Trp row said 0-3, which is the
                   PLATEAU ("steady at higher concentrations"), not
                   the top point.  This column is a maximum.
    PFOA   0.8 mM  text, "in the concentration range of 0-0.8 mM";
                   <=2.5 mM at the two displacement sites, stock
    PFDoA  2 mM    text, "addition of up to 2 mM", solubility-limited
    PFBS   <=5 mM  stock; the 0-6 axis is above it
    PFBA   <=6 mM  shared axis, and the 20 mM stock allows it

THE SECOND LIGAND, recorded because the ratio column cannot show it.
The twelve displacement rows are not single-ligand experiments.
"Competitive fluorescence measurements were carried out using 0.1 mM
HSA and 0.1 mM DA or DP" -- the probe is EQUIMOLAR with the albumin,
one dansyl ligand per protein already bound before any PFAS is added.
Every displacement constant in Table 1 is therefore a competition
constant against an occupied site, which is the point of the method
but is not what a ligand:protein ratio of 50 describes.

WHAT IS NOT TOUCHED.  Every constant in Table 1 was checked against
the PDF and all fifteen cells agree, qualitative cells included
(PFDoA is `*` at Site I, no binding, and `**` at Trp, binds but no
constant -- the database has both the right way round).  Protein
0.1 mM, buffer 100 mM phosphate pH 7.4, and the absence of a stated
temperature are all confirmed.  The five butyric/octanoic acid rows
are from the Discussion, not this paper's own experiments, and have
no range to correct.

Run:  python3 tk_review/scripts/fix_chen2009_ranges.py [--dry-run]
Idempotent.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Chen & Guo", "2009"
HSA_MM = 0.1                      # Experimental; same 0.1 mM in every caption

# Experimental section, stock concentrations in mM.  A titration cannot
# exceed its stock, so this is a hard ceiling per compound.
STOCK_MM = {"PFBA": 20.0, "PFBS": 5.0, "PFOA": 2.5, "PFDoA": 10.0,
            "PFOS": 20.0}

# (chemical, site) -> (top mM, provenance).  Sites spelled as the database
# spells them.
TOP = {
 ("PFOS", "Trp214"):
    (4.0, "Fig 1 caption lists 0/0.4/0.7/1.5/2/4 mM; the text's 3 mM is "
          "where the signal plateaus, not the top point"),
 ("PFOS", "Sudlow site I"):
    (4.0, "Fig 3 caption lists 0/0.2/0.7/2/4 mM; text 'from 0.7 to 4 mM'"),
 ("PFOS", "Sudlow site II"):
    (6.0, "Fig 5 shared x-axis, 0-6 mM, an upper bound not a series; the "
          "20 mM stock allows it. Fig S7 holds the actual series"),
 ("PFOA", "Trp214"):
    (0.8, "text, 'in the concentration range of 0-0.8 mM'"),
 ("PFOA", "Sudlow site I"):
    (2.5, "STOCK CEILING. Fig 4's shared 0-6 mM axis is 2.4x above the "
          "2.5 mM PFOA stock, so 6 mM was never in the cuvette; the "
          "series (Fig S4) is not given numerically"),
 ("PFOA", "Sudlow site II"):
    (2.5, "STOCK CEILING. Fig 5's shared 0-6 mM axis is 2.4x above the "
          "2.5 mM PFOA stock; the series (Fig S4) is not given "
          "numerically"),
 ("PFDoA", "Trp214"):
    (2.0, "text, 'addition of up to 2 mM PFDoA'; solubility-limited, "
          "'could not be increased further to reach a plateau'"),
 ("PFDoA", "Sudlow site I"):
    (2.0, "text, 2 mM; 'in the limited concentration range examined'"),
 ("PFDoA", "Sudlow site II"):
    (2.0, "text, 2 mM, solubility-limited"),
 ("PFBS", "Trp214"):
    (5.0, "STOCK CEILING. The shared 0-6 mM axis is above the 5 mM PFBS "
          "stock, so 6 mM was never reached"),
 ("PFBS", "Sudlow site I"):
    (5.0, "STOCK CEILING. The shared 0-6 mM axis is above the 5 mM PFBS "
          "stock"),
 ("PFBS", "Sudlow site II"):
    (5.0, "STOCK CEILING. The shared 0-6 mM axis is above the 5 mM PFBS "
          "stock"),
 ("PFBA", "Trp214"):
    (6.0, "Figs 2/4/5 shared x-axis, 0-6 mM, an upper bound not a series; "
          "the 20 mM stock allows it"),
 ("PFBA", "Sudlow site I"):
    (6.0, "Figs 2/4/5 shared x-axis, 0-6 mM; the 20 mM stock allows it"),
 ("PFBA", "Sudlow site II"):
    (6.0, "Figs 2/4/5 shared x-axis, 0-6 mM; the 20 mM stock allows it"),
}

PROBE = ("the probe is EQUIMOLAR with the protein: 'competitive "
         "fluorescence measurements were carried out using 0.1 mM HSA and "
         "0.1 mM DA or DP', so one dansyl ligand per albumin is already "
         "bound before any PFAS is added and this constant is a "
         "competition constant, not a ratio of PFAS to free protein")
STOCKNOTE = ("Experimental stock concentrations, the hard ceiling on any "
             "titration point: PFBA 20 mM, PFBS 5 mM, PFOA 2.5 mM, "
             "PFDoA 10 mM, PFOS 20 mM (from 200 mM in DMSO)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    hit = [r for r in rows if r["study"] == STUDY and r["year"] == YEAR]
    print(f"  {len(hit)} Chen & Guo 2009 rows "
          f"({sum(1 for r in hit if r['n_perfluorinated_carbons'] == '0')} "
          f"non-fluorinated reference, no range to correct)\n")

    print(f"  {'PFAS':7}{'site':16}{'was':>9}{'now':>9}{'ratio':>14}  basis")
    nr = nn = 0
    for r in sorted(hit, key=lambda r: (r["chemical"],
                                        r["binding_site_identified"])):
        k = (r["chemical"], r["binding_site_identified"])
        if k not in TOP:
            continue
        top, why = TOP[k]
        cap = STOCK_MM[r["chemical"]]
        assert top <= cap, f"{k} top {top} over stock {cap}"
        # The cell stays a plain numeric range so the range parser in
        # build_binding_affinity_table.py still reads it; whether the top
        # is a measured point or a ceiling is stated in the note.
        rng = f"0-{top:g}"
        ratio = f"{top / HSA_MM:g}"

        was_r, was_v = r["ligand_conc_range"], r["ligand_over_protein_max"]
        if was_r != rng:
            nr += 1
        if was_v != ratio:
            nn += 1
        r["ligand_conc_range"] = rng
        r["ligand_over_protein_max"] = ratio

        # Provenance, and the probe for the two displacement sites.
        note = r["assay_volume_note"]
        for tag in (f"ligand range: {why}", STOCKNOTE,
                    PROBE if "Sudlow" in k[1] else None):
            if tag and tag not in note:
                note += " | " + tag
        r["assay_volume_note"] = note

        flag = "  <-- IMPOSSIBLE BEFORE" if was_r == "0-6" and top < 6 else ""
        print(f"  {k[0]:7}{k[1]:16}{was_r:>9}{rng:>9}"
              f"{was_v + ' -> ' + ratio:>14}{flag}")

    print(f"\n  {nr} ranges and {nn} ratios changed")
    print(f"  ratio now spans {min(float(r['ligand_over_protein_max']) for r in hit if r['ligand_over_protein_max']):g} "
          f"to {max(float(r['ligand_over_protein_max']) for r in hit if r['ligand_over_protein_max']):g} "
          f"ligands per albumin at 0.1 mM HSA")

    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
