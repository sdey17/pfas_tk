#!/usr/bin/env python3
"""Moro 2022: one row per site class per temperature, read at source.

Moro G, Liberi S, Vascon F, Linciano S, De Felice S, Fasolato S,
Foresta C, De Toni L, Di Nisio A, Cendron L, Angelini A.
"Investigation of the Interaction between Human Serum Albumin and
Branched Short-Chain Perfluoroalkyl Compounds."  Chem Res Toxicol
2022;35(11):2049-2058.  doi 10.1021/acs.chemrestox.2c00211,
PMID 36148994, PMC9682524 (open access).

WHY THIS EXISTS.  The database held 8 rows, two per compound, and
PFOA and HFPO-DA were unusable in both of them:

  kd = "19.0 (KD1); 84.2 (KD2)"
  kd = "0.4 (high-affinity set); 29.7 (low-affinity set)"
  n_binding_sites = "sequential two-site model"
  temperature_c = "25, 37 (298, 310 K)"        <- on EVERY row

Two constants in one cell is not a constant, so albumin_kd_to_k
marked all four rows "two classes of site in one cell, not one
constant" and they reach no workbook, no figure, nothing.  That is
why PFOA and HFPO-DA look absent from this study.  Worse, the
temperature cell names BOTH temperatures on every row, so even C6O4
and PFHxA -- which did convert -- could not say which of their two
values was 25 C and which was 37 C.

Table 1 resolves all of it.  Its columns are ordered 298 K then
310 K within each compound, and the existing rows were already in
that order, so the pairing is recoverable rather than guessed.
12 rows out, one per (compound, site class, temperature).

A DEFECT IN THE PAPER, recorded and NOT silently corrected.
Table 1 gives PFOA's high-affinity Kd as 0.4 uM at 298 K and
39.3 uM at 310 K.  The 39.3 cannot be right:

  dG = RT ln(Kd) reproduces every other cell in Table 1 to within
  0.0-0.8 kJ/mol.  For PFOA high-affinity at 310 K the printed
  dG of -38.10 kJ/mol implies Kd = 0.38 uM, while the printed
  39.3 uM implies -26.1 -- a 12.0 kJ/mol contradiction, the only
  one of its size in the table, and exactly a factor of 100.

  The paper's own text says "no significant variations were
  observed in the stoichiometry nor in the binding affinities when
  ITC studies were performed at temperatures of 298 or 310 K".
  0.4 -> 39.3 uM is a 98-fold change; 0.4 -> 0.39 is none.  The
  text and the thermodynamics agree with each other and against
  the printed Kd.

  Figure 1e and the Discussion both quote 0.4 uM as PFOA's
  affinity with no temperature caveat.

The printed value is kept, because the database records what a
paper reports.  The row carries the contradiction, the implied
0.38 uM, and a warning not to pool it.

THE RATIO WAS NOT A RATIO.  The database cell was BLANK on all 8
rows, and the 36.36 that reaches the workbook and the loading
ranking is computed downstream, in build_binding_affinity_table's
ratio_max, as 4000 uM / 110 uM -- the SYRINGE concentration over
the CELL concentration.  In ITC those are never in the same vessel:
the syringe titrates into the cell.
The loading that matters is the molar ratio reached in the cell at
the end of the titration, and from the paper's stated quantities
(110 uM dhSA in the cell, 4 mM ligand in the syringe, an initial
0.4 uL injection plus 25 x 1.5 uL) with the PEAQ-ITC's 200 uL cell:

  simple dilution          4000 * 37.9/200 / 110        = 6.9
  with overflow correction [L]/[P] as the cell displaces = 7.6

so about 7 ligands per albumin, not 36.  This matters beyond this
study: 36 put Moro at the top of the loading ranking, above Peng's
18, when it belongs near Chen & Guo's 8.  The cell volume is an
instrument constant the paper does not state, so it is named in the
note as an assumption.

Filling the database cell is NOT enough on its own: _ratio_cell
prefers its own computed quotient unless the stated value is more
than 10x larger, so 36.36 would still win.  The companion change in
build_binding_affinity_table stops it computing a quotient for ITC
at all, because a syringe concentration over a cell concentration is
not a molar ratio.  That change moves this study only; every other
ITC study either already carries a stated value or has no
concentrations to divide.

NAME.  The paper writes "HPFO-DA" throughout, but its own Materials
section says "Hexafluoropropylene oxide-dimer acid (HFPO-DA or
GenX)".  HFPO-DA is correct and is what the database uses; HPFO-DA
is a transposition the paper repeats.  Noted so nobody "fixes" it
the wrong way.

STOICHIOMETRY.  Table 1 reports no n for HFPO-DA, because a
sequential-binding fit does not have one: the model is two ordered
sites of one ligand each, so n = 1 per site is fixed by the model,
not fitted -- the same standing as a paper's own single-site model.
The crystal structure (PDB 7Z57) independently finds exactly two
HFPO-DA molecules, in FA3 and FA4.

Run:  python3 tk_review/scripts/split_moro2022_sites.py [--dry-run]
Idempotent: removes its own 12 rows, or the original 8, and rewrites.
"""
import argparse
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Moro et al.", "2022"
PMID, DOI = "36148994", "10.1021/acs.chemrestox.2c00211"
BRIDGE = 4.1774                      # log10(1e9 / 66463)
R = 8.314e-3                         # kJ/mol/K

CELL_UL, SYR_UM, PROT_UM = 200.0, 4000.0, 110.0
V_INJ = 0.4 + 25 * 1.5               # uL actually entering the cell
# Standard ITC overflow treatment: both species dilute as the cell
# displaces, so [L]/[P] tends to (1-e^-x)/e^-x with x = V_inj/V_cell.
_x = V_INJ / CELL_UL
RATIO = f"{(SYR_UM * (1 - math.exp(-_x))) / (PROT_UM * math.exp(-_x)):.1f}"

PROTEIN = "defatted recombinant human serum albumin (dhSA)"
BUFFER = "50 mM sodium phosphate, 100 mM NaCl, pH 7.4"
METHOD = "isothermal titration calorimetry (ITC), Microcal PEAQ-ITC"

# compound: (nCF, head_group, model, [(site label, T_K, Kd_uM, n, n_basis)])
DATA = {
 "C6O4": (3, "branched ether carboxylate", "one set of sites", [
    ("single set of sites", 298, 2.4, "1.6", "fitted"),
    ("single set of sites", 310, 5.4, "1.5", "fitted")]),
 "PFHxA": (5, "carboxylate", "one set of sites", [
    ("single set of sites", 298, 4.5, "1.02", "fitted"),
    ("single set of sites", 310, 5.9, "1.03", "fitted")]),
 "HFPO-DA (GenX)": (4, "ether carboxylate", "sequential binding", [
    ("sequential site 1 (KD1, higher affinity)", 298, 19.0, "1", "model"),
    ("sequential site 2 (KD2, lower affinity)", 298, 84.2, "1", "model"),
    ("sequential site 1 (KD1, higher affinity)", 310, 9.9, "1", "model"),
    ("sequential site 2 (KD2, lower affinity)", 310, 27.0, "1", "model")]),
 "PFOA": (7, "carboxylate", "two sets of sites", [
    ("high-affinity set", 298, 0.4, "0.89", "fitted"),
    ("high-affinity set", 310, 39.3, "0.95", "fitted"),
    ("low-affinity set", 298, 29.7, "4.5", "fitted"),
    ("low-affinity set", 310, 46.0, "9.5", "fitted")]),
}

# Printed Gibbs free energies, for the cross-check that found the defect.
DG = {("C6O4", 298): -32.10, ("C6O4", 310): -31.70,
      ("PFHxA", 298): -30.60, ("PFHxA", 310): -31.70,
      ("PFOA high-affinity set", 298): -36.10,
      ("PFOA high-affinity set", 310): -38.10,
      ("PFOA low-affinity set", 298): -26.60,
      ("PFOA low-affinity set", 310): -29.10}

SITE_LOC = {
 "C6O4": "Sudlow site II region, FA4 and/or FA6 in subdomain IIIA; "
         "assigned by ibuprofen competition, not crystallography",
 "PFHxA": "Sudlow site II region, FA4 and/or FA6 in subdomain IIIA; "
          "assigned by ibuprofen competition, not crystallography",
 "HFPO-DA (GenX)": "Sudlow site II, subdomain IIIA; the crystal structure "
                   "(PDB 7Z57) places one molecule in FA3 and one in FA4. "
                   "The paper's mapping of ITC site 1 to FA3 and site 2 to "
                   "FA4 is its own speculation ('could probably'), so the "
                   "site class here is the ITC model's, not the structure's",
 "PFOA": "FA4 in subdomain IIIA by the earlier structure 7AAI; this paper "
         "measured PFOA by ITC only, as a long-chain control",
}

RATIO_NOTE = (
    f"LOADING: ligand_over_protein_max is the molar ratio reached in the "
    f"ITC CELL at the end of the titration, ~{RATIO}. It is NOT the "
    f"4000/110 = 36.4 that was reaching the workbook: 4 mM is the SYRINGE "
    f"concentration and 110 uM the CELL concentration, and in ITC those "
    f"are never in the same vessel. Derived from the paper's stated "
    f"110 uM dhSA, 4 mM ligand, 0.4 uL + 25 x 1.5 uL injections, with the "
    f"PEAQ-ITC's {CELL_UL:.0f} uL cell -- an instrument constant the paper "
    f"does NOT state, so this rests on that one assumption. Simple "
    f"dilution gives 6.9, the standard overflow correction {RATIO}")

DEFAT = (
    "DEFATTED albumin (activated charcoal, pH 3, then SEC). Endogenous "
    "fatty acids compete for these same subdomain IIIA sites, and this "
    "paper's own crystal structure shows HFPO-DA displacing myristate "
    "from FA3, so constants measured on defatted albumin are an UPPER "
    "bound relative to native albumin")

SUSPECT = (
    "SUSPECT VALUE, kept as printed. Table 1 gives this high-affinity Kd "
    "as 39.3 +/- 0.4 uM at 310 K against 0.4 +/- 0.01 uM at 298 K, but "
    "dG = RT ln(Kd) reproduces every other cell in Table 1 to within "
    "0.0-0.8 kJ/mol and here contradicts it by 12.0: the printed dG of "
    "-38.10 kJ/mol implies Kd = 0.38 uM, exactly 1/100 of the printed "
    "value. The paper's own text says 'no significant variations were "
    "observed in the stoichiometry nor in the binding affinities' between "
    "298 and 310 K, which holds at 0.39 uM and fails at 39.3; Figure 1e "
    "and the Discussion both quote 0.4 uM with no temperature caveat. "
    "DO NOT POOL this row or read a temperature dependence from it "
    "without checking the published table")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(DB, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        flds, rows = list(rd.fieldnames), list(rd)

    old = [r for r in rows if r["study"] == STUDY and r["year"] == YEAR]
    keep = [r for r in rows if not (r["study"] == STUDY and r["year"] == YEAR)]
    print(f"  removing {len(old)} existing Moro 2022 row(s):")
    for r in old:
        print(f"    {r['chemical']:16} kd={r['kd'][:42]:44} "
              f"T={r['temperature_c']}")

    new = []
    for cm, (ncf, head, model, sites) in DATA.items():
        for label, T, kd, n, nbasis in sites:
            nf = float(n)
            logk = BRIDGE + math.log10(nf) - math.log10(kd)
            dgk = DG.get((f"{cm} {label}", T), DG.get((cm, T)))
            cross = ""
            if dgk is not None:
                implied = math.exp(dgk / (R * T)) * 1e6
                cross = (f" | printed dG {dgk} kJ/mol implies Kd "
                         f"{implied:.3g} uM vs the printed {kd} uM")
            nb = ("fitted" if nbasis == "fitted" else
                  "FIXED AT 1 BY THE SEQUENTIAL MODEL, not fitted")
            note = [f"ITC on {PROTEIN}, '{model}' fit at {T} K "
                    f"({T - 273:.0f} C)", DEFAT, RATIO_NOTE]
            if cm == "PFOA" and label == "high-affinity set" and T == 310:
                note.insert(1, SUSPECT)
            new.append(dict(
                study=STUDY, year=YEAR, pmid=PMID, doi=DOI,
                chemical=cm, n_perfluorinated_carbons=str(ncf),
                head_group=head, protein=PROTEIN,
                protein_species="human", method=METHOD,
                # The paper reports Kd, so Kd is the measurement. The
                # reciprocal is stored beside it under the database's
                # existing marker so measure_of reads Kd from the kd
                # column: a 4-significant-figure Ka round-tripped back
                # would turn the paper's 9.9 uM into 9.901 and 84.2 into
                # 84.18.
                ka=f"{1e6 / kd:.6g}",
                ka_units="M^-1, computed as 1/K_D",
                kd=f"{kd:g}", kd_units="uM",
                n_binding_sites=n,
                binding_site_identified=f"{label}; {SITE_LOC[cm]}",
                unbound_fraction="",
                protein_conc=f"{PROT_UM:g}", protein_conc_units="uM",
                ligand_conc_range=f"{SYR_UM:g} (syringe, not the cell)",
                ligand_conc_units="uM",
                temperature_c=f"{T - 273:.0f} ({T} K)",
                buffer=BUFFER, is_computational="no",
                tk_relevance_noted=(
                    "partly - discussed for absorption/distribution of "
                    "short-chain alternatives. Both new-generation "
                    "compounds bind WEAKER than PFOA, which the paper "
                    "reads as a possibly shorter serum half-life"),
                source_table=(
                    f"Table 1, '{model}' block, {T} K column; "
                    f"n {nb}{cross} | read from PMC9682524, open "
                    f"access, primary"),
                assay_volume_note=" | ".join(note),
                ligand_over_protein_max=RATIO))

    new.sort(key=lambda r: (r["chemical"], r["binding_site_identified"],
                            r["temperature_c"]))
    print(f"\n  adding {len(new)} rows, one per (compound, site, T):\n")
    print(f"    {'PFAS':16}{'T':>6}{'site':34}{'Kd uM':>8}{'n':>7}"
          f"{'log K':>8}")
    for r in new:
        flag = "  <-- SUSPECT" if "SUSPECT" in r["assay_volume_note"] else ""
        lk = BRIDGE + math.log10(float(r["n_binding_sites"])) \
            - math.log10(float(r["kd"]))
        print(f"    {r['chemical']:16}{r['temperature_c'].split()[0]:>6}"
              f" {r['binding_site_identified'].split(';')[0][:32]:33}"
              f"{r['kd']:>8}{r['n_binding_sites']:>7}{lk:>8.3f}{flag}")

    print(f"\n  ligand:protein -> {RATIO} on all rows (the DB cell was "
          f"blank; 36.36 was computed downstream as syringe/cell "
          f"concentration, which is not a molar ratio)")

    if args.dry_run:
        print("  (dry run, nothing written)")
        return 0
    out = keep + new
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=flds)
        w.writeheader()
        w.writerows(out)
    print(f"  {len(rows)} -> {len(out)} rows; wrote "
          f"{os.path.relpath(DB, os.path.join(HERE, '..'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
