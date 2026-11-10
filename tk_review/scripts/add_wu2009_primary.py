#!/usr/bin/env python3
"""Wu 2009 read at source, replacing a secondary one-line stub.

Wu LL, Gao HW, Gao NY, Chen FF, Chen L.  "Interaction of
perfluorooctanoic acid with human serum albumin."  BMC Struct Biol
2009;9:31.  doi 10.1186/1472-6807-9-31, PMID 19442292, PMC2688008.
Open access.

WHAT WAS THERE.  One row, Ka 3.12e4, method "equilibrium dialysis",
source "secondary (Starnes & Belcher 2026 Table 1); primary not
retrieved", and no concentrations at all.  The primary is open access
and carries four methods, a second constant the database did not
have, and every concentration stated explicitly.

THE UNITS HAD TO BE RECOVERED.  Table 1's association-constant column
header is "(10 M)" with the exponent lost, so the constants are
ambiguous as printed -- 24.7 and 0.55.  The printed Gibbs energies
are not ambiguous, and K = exp(-dG/RT) at 310 K gives 23,002 and 542
against the printed 24.7 and 0.55, i.e. a factor of 931 and 986.
Both land on 1e3, so the column is 10^3 M^-1:

    step 1  Ka 2.47e4 M^-1  ->  Kd   40.5 uM
    step 2  Ka 5.5e2  M^-1  ->  Kd 1818   uM

That is corroborated from outside the table: the equilibrium-dialysis
Temkin constant, 3.12e4, is 1.3x the ITC step-1 value, and the two
methods agree exactly on the total binding number (13).

THE TEMKIN CONSTANT IS NOT A SITE CONSTANT, which matters because it
is the one value the database already had.  It is the adsorption
constant of a Temkin monolayer ADSORPTION isotherm fitted to the
dialysis data, not a 1:1 equilibrium binding constant, so inverting
it to a molar Kd is not the same operation as inverting the ITC
constant.  Recorded as reported, with Kd left empty and the reason
in the note; the ITC rows carry the invertible constants.

STEP 2 IS SURFACE COATING, NOT A POCKET.  The paper is explicit:
"only one PFOA molecule entered the hydrophobic intracavity in the
first step and 12 PFOA molecules binding on the hydrophilic outer
surface in the second step", with Ka 550 M^-1, i.e. Kd 1.8 mM.  That
is the same shape as Moro 2022's HFPO-DA second event, and the same
conclusion follows: at the trace occupancy of blood only step 1 is
occupied, and a log K built on step 2 describes a protein coated
with PFOA, not a binding site.

CONCENTRATIONS, all stated, all at 37 C, pH 7.40 Britton-Robinson
buffer, 0.15 M NaCl:

  equilibrium dialysis  0.033 mM HSA in a 12.0 mL bag against 36.0 mL,
                        5 kDa MWCO, 4 h; saturation at 13 bound per
                        HSA once PFOA:HSA exceeds 20
  ITC                   10.0 uM HSA in a 1.4685 mL cell -- the cell
                        volume is STATED, in the Figure 3 caption --
                        titrated with 2.50 mM PFOA, ~40 x 6.0 uL, so
                        the run reaches ~41 ligands per albumin (44
                        with the overflow correction).  Table 1 fits
                        step 1 below PFOA:HSA 1.5 and step 2 from
                        1.5 to 20
  fluorescence          PFOA:HSA 0, 1, 1.5, 2, 3, 5, 10, 20, 30
  circular dichroism    0.010 mM HSA with 1.00 mM PFOA, ratio 100

Run:  python3 tk_review/scripts/add_wu2009_primary.py [--dry-run]
Idempotent: removes any Wu 2009 rows, then writes its own 3.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Wu et al.", "2009"
PMID, DOI = "19442292", "10.1186/1472-6807-9-31"
BUFFER = ("Britton-Robinson buffer pH 7.40 with 0.15 M NaCl, chosen to "
          "match blood: 'normal physiological condition, pH 7.40 and "
          "0.15 M electrolyte and 37 C'")
TEMP = "37 (310 K)"
COMMON = ("All four methods at 37 C, pH 7.40 BR buffer, 0.15 M NaCl. "
          "Equilibrium dialysis 0.033 mM HSA in a 12.0 mL bag against "
          "36.0 mL, 5 kDa MWCO, 4 h. ITC 10.0 uM HSA in a 1.4685 mL cell "
          "(the cell volume is STATED, Figure 3 caption) titrated with "
          "2.50 mM PFOA, ~40 x 6.0 uL at 3 min intervals, so the run "
          "reaches ~41 ligands per albumin, 44 with the overflow "
          "correction. Fluorescence at PFOA:HSA 0/1/1.5/2/3/5/10/20/30. "
          "CD at 0.010 mM HSA with 1.00 mM PFOA, ratio 100")
UNITS = ("Table 1's constant column header is '(10 M)' with the exponent "
         "lost, so 24.7 and 0.55 are ambiguous as printed. The printed "
         "Gibbs energies are not: K = exp(-dG/RT) at 310 K gives 23,002 "
         "and 542, a factor of 931 and 986 on the printed values, so the "
         "column is 10^3 M^-1")

ROWS = [
 dict(site="Temkin monolayer adsorption (whole protein)",
      method=("equilibrium dialysis, Temkin isothermal adsorption model, "
              "5 kDa MWCO membrane"),
      ka="3.12e4", ka_units="M^-1 (Temkin adsorption constant, NOT a 1:1 "
                            "site constant)",
      kd="", kd_units="",
      n="13 (saturation binding number, PFOA bound per HSA)",
      prot="33", lig="series up to >20x the HSA concentration",
      ratio="20",
      note=("Temkin ADSORPTION constant from a monolayer adsorption fit "
            "to the dialysis data, with saturation adsorption energy "
            "-3.89 kJ/mol. It is not a 1:1 equilibrium binding constant, "
            "so Kd is left EMPTY rather than inverting it -- the ITC rows "
            "carry the invertible constants. This 3.12e4 is the single "
            "value the database previously held, via Starnes & Belcher "
            "2026 Table 1; the primary confirms it. Ratio 20 is where "
            "the paper observes saturation ('approaches a constant "
            "maximum at 13 when PFOA:HSA is more than 20'), so the "
            "series ran at least that far; no explicit top is given"),
      src="Results, equilibrium dialysis section, and Figure 2"),
 dict(site="step 1, hydrophobic intracavity (subdomains IIA/IIB, across "
           "Trp214)",
      method=("isothermal titration calorimetry (ITC), MicroCal MCS-ITC, "
              "two-step sequential binding model"),
      ka="2.47e4", ka_units="M^-1",
      kd="40.5", kd_units="uM",
      n="1",
      prot="10", lig="2500 (syringe); fitted below PFOA:HSA 1.5",
      ratio="1.5",
      note=("THE PHYSIOLOGICALLY RELEVANT ROW. One PFOA in the "
            "hydrophobic intracavity, the step that is occupied at the "
            "trace loading of blood. Table 1: Ka 24.7 +/- 6.5 x10^3 "
            "M^-1, dH -10.0 kcal/mol, dS -12.3 cal/mol/K, dG -6.19 "
            "kcal/mol, fitted over PFOA:HSA < 1.5. " + UNITS),
      src="Table 1, step 1"),
 dict(site="step 2, hydrophilic OUTER SURFACE (12 sites)",
      method=("isothermal titration calorimetry (ITC), MicroCal MCS-ITC, "
              "two-step sequential binding model"),
      ka="5.5e2", ka_units="M^-1",
      kd="1818", kd_units="uM",
      n="12",
      prot="10", lig="2500 (syringe); fitted over PFOA:HSA 1.5 to 20",
      ratio="20",
      note=("SURFACE COATING, NOT A POCKET -- do not pool with site "
            "constants. The paper is explicit: 'only one PFOA molecule "
            "entered the hydrophobic intracavity in the first step and "
            "12 PFOA molecules binding on the hydrophilic outer surface "
            "in the second step'. Kd 1.8 mM is ~4500x weaker than step 1 "
            "and this regime starts above PFOA:HSA 1.5, roughly 1e5 "
            "times physiological occupancy. Same shape as Moro 2022's "
            "HFPO-DA second event. Table 1: Ka 0.55 +/- 0.14 x10^3 "
            "M^-1, dH -13.8 kcal/mol, dS -32.0 cal/mol/K, dG -3.88 "
            "kcal/mol. " + UNITS),
      src="Table 1, step 2"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    old = [r for r in rows if r["study"] == STUDY and r["year"] == YEAR]
    keep = [r for r in rows if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"  removing {len(old)} existing Wu 2009 row(s):")
    for r in old:
        print(f"    ka={r['ka']:10} {r['source_table'][:72]}")

    new = []
    for d in ROWS:
        new.append(dict(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI,
            chemical="PFOA", n_perfluorinated_carbons="7",
            head_group="carboxylate",
            protein="human serum albumin (HSA)", protein_species="human",
            method=d["method"], ka=d["ka"], ka_units=d["ka_units"],
            kd=d["kd"], kd_units=d["kd_units"], n_binding_sites=d["n"],
            binding_site_identified=d["site"], unbound_fraction="",
            protein_conc=d["prot"], protein_conc_units="uM",
            ligand_conc_range=d["lig"], ligand_conc_units="uM",
            temperature_c=TEMP, buffer=BUFFER, is_computational="no",
            tk_relevance_noted=(
                "yes - run deliberately at blood conditions (37 C, "
                "pH 7.40, 0.15 M NaCl). Two methods agree that HSA "
                "carries 13 PFOA at saturation, but only ONE of those "
                "is a pocket; the other 12 coat the surface and need "
                "millimolar PFOA"),
            source_table=f"{d['src']} | read from PMC2688008, open "
                         f"access, primary",
            assay_volume_note=COMMON + " | " + d["note"],
            ligand_over_protein_max=d["ratio"]))

    print(f"\n  adding {len(new)} rows:")
    print(f"    {'site':46}{'Ka M^-1':>10}{'Kd uM':>9}{'n':>4}"
          f"{'[HSA] uM':>10}{'ratio':>7}")
    for r in new:
        print(f"    {r['binding_site_identified'][:45]:46}{r['ka']:>10}"
              f"{r['kd'] or '-':>9}{r['n_binding_sites'].split()[0]:>4}"
              f"{r['protein_conc']:>10}{r['ligand_over_protein_max']:>7}")

    if args.dry_run:
        print("\n  (dry run, nothing written)")
        return 0
    out = keep + new
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(out)
    print(f"\n  {len(rows)} -> {len(out)} rows; wrote "
          f"{os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
