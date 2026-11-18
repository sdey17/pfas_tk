#!/usr/bin/env python3
"""Louisse 2024: the OAT1/OAT2/OAT3 Km table, replacing three qualitative rows.

Louisse J, Pedroni L, van den Heuvel JJMW, Rijkers D, Bokkers BGH, Dorne JCM,
Dellafiora L, Punt A, Russel FGM, Koenderink JB. "Transport of perfluoroalkyl
acids by human organic anion transporters OAT1, OAT2, OAT3 and OAT4."
Toxicology 2024;509:153961. doi:10.1016/j.tox.2024.153961  PMID 39343156.

WHAT WAS WRONG. The main table carried three Louisse 2024 rows, one per
transporter, each holding a substrate list in the `chemical` field and no
constant, with source_table reading "full Km table not retrieved (paywalled)".
Table 1 has since been obtained (db/primary_2026/louisse2024_human_oat_km.csv).
This script replaces those three rows with sixteen per-compound rows.

SIX ARE GENUINELY NEW -- two OAT1 and four OAT3, each with Km, Vmax and the
authors' own transport efficiency:

    OAT1  PFHpA   Km 159 +/- 98 uM    Vmax 3.4   eff 21 uL/min/mg
    OAT1  PFOA    Km 185 +/- 97       Vmax 3.5   eff 19
    OAT3  PFHpA   Km 367 +/- 199      Vmax 3.2   eff  9
    OAT3  PFOA    Km  90 +/- 63       Vmax 1.5   eff 17
    OAT3  PFNA    Km  42 +/- 38       Vmax 1.0   eff 24
    OAT3  PFHxS   Km 778 +/- 8        Vmax 5.8   eff  8

THE EXTRACTION VERIFIES AGAINST ITSELF. Efficiency is Vmax/Km, so the printed
third column is predicted by the first two: 3.4/159 = 21.4, 3.5/185 = 18.9,
3.2/367 = 8.7, 1.5/90 = 16.7, 1.0/42 = 23.8, 5.8/778 = 7.5 uL/min/mg. All six
round to the printed value. The carried-over OAT4 rows do too (4.5/47 = 95.7
against a printed 96), so the whole table is internally consistent.

FOUR ARE NOT NEW, and are stored WITHOUT constants. Table 1 also tabulates
OAT4, but those are Louisse et al. 2023's own earlier values, already in this
database from the primary. They are digit-identical:

    PFHpA    60 / 4.5   vs   60 (15) / 4.5 (0.4)
    PFOA     47 / 4.5   vs   47 (15) / 4.5 (0.5)
    PFNA     58 / 8.5   vs   58 (15) / 8.5 (0.7)
    PFHxS    92 / 7.3   vs   92 (33) / 7.3 (1.0)

and the re-tabulation even loses precision: it prints the standard error as 15
for all four, where the 2023 primary gives 33 for PFHxS. Entering them again
would turn one measurement into two and tighten the apparent consensus for
free. This is the treatment Han 2021 and Yang 2026 get on the binding side.

FIVE ARE NON-SUBSTRATES, stored with no constant and the finding recorded:
PFNA, PFBS, PFHxS and PFOS on OAT1; PFBS on OAT3. OAT2 transported nothing at
all, which is its own row.

THE TWO WEAKEST Km VALUES SHOULD NOT BE READ AS NUMBERS. OAT3/PFHxS is
778 +/- 8 uM, a standard error of 1% where the same table's OAT3/PFOA carries
70%; and OAT1/PFHpA is 159 +/- 98. A Km far above the highest concentration
tested is a lower bound on the half-saturation point, not an estimate of it,
and the efficiency column is the quantity the authors themselves rely on.
Flagged in assay_volume_note, carried, and left for the inclusion rule.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "transporter_kinetics.csv")
DOI = "10.1016/j.tox.2024.153961"
STUDY = ("Louisse J, Pedroni L, van den Heuvel JJMW, Rijkers D, Bokkers BGH, "
         "Dorne JCM, Dellafiora L, Punt A, Russel FGM, Koenderink JB")
SYSTEM = "OAT-transduced HEK293 cells vs control-transduced, 4-min uptake"
NCF = {"PFHpA": ("6", "carboxylate"), "PFOA": ("7", "carboxylate"),
       "PFNA": ("8", "carboxylate"), "PFBS": ("4", "sulfonate"),
       "PFHxS": ("6", "sulfonate"), "PFOS": ("8", "sulfonate")}
GENE = {"OAT1": "SLC22A6", "OAT2": "SLC22A7", "OAT3": "SLC22A8",
        "OAT4": "SLC22A11"}
DIRN = {"OAT1": "secretion (basolateral uptake into the proximal tubule cell)",
        "OAT2": "secretion (basolateral; predominantly hepatic in human)",
        "OAT3": "secretion (basolateral uptake into the proximal tubule cell)",
        "OAT4": "reabsorption (apical, lumen -> cell)"}
CARRIED = (
    "NOT MEASURED HERE. Table 1 tabulates this value but it is Louisse et al. "
    "2023's, already in this database from that primary, and digit-identical "
    "to it; the re-tabulation also prints the standard error as 15 for all "
    "four OAT4 compounds where the 2023 primary gives 33 for PFHxS. Stored "
    "without a constant so one measurement is not counted twice")
WEAK = ("Km POORLY DETERMINED: standard error {} of the estimate. A Km of this "
        "size is a lower bound on the half-saturation point, not an estimate "
        "of it; the transport efficiency Vmax/Km is the quantity the authors "
        "rely on")

#  transporter, chemical, km, km_se, vmax, efficiency, note
ROWS = [
    # --- measured here ---
    ("OAT1", "PFHpA", "159", "98", "3.4", "21", WEAK.format("62%")),
    ("OAT1", "PFOA", "185", "97", "3.5", "19", WEAK.format("52%")),
    ("OAT3", "PFHpA", "367", "199", "3.2", "9", WEAK.format("54%")),
    ("OAT3", "PFOA", "90", "63", "1.5", "17", WEAK.format("70%")),
    ("OAT3", "PFNA", "42", "38", "1.0", "24", WEAK.format("90%")),
    ("OAT3", "PFHxS", "778", "8", "5.8", "8",
     WEAK.format("1%, implausibly tight beside the 52-90% carried by every "
                 "other row of the same table")),
    # --- carried over from Louisse 2023, stored WITHOUT a constant ---
    ("OAT4", "PFHpA", "", "", "", "75", CARRIED + " (60 / 4.5 both places)"),
    ("OAT4", "PFOA", "", "", "", "96", CARRIED + " (47 / 4.5 both places)"),
    ("OAT4", "PFNA", "", "", "", "147", CARRIED + " (58 / 8.5 both places)"),
    ("OAT4", "PFHxS", "", "", "", "79", CARRIED + " (92 / 7.3 both places)"),
    # --- tested, not transported ---
    ("OAT1", "PFNA", "", "", "", "", "not transported (Table 1)"),
    ("OAT1", "PFBS", "", "", "", "", "not transported (Table 1)"),
    ("OAT1", "PFHxS", "", "", "", "", "not transported (Table 1)"),
    ("OAT1", "PFOS", "", "", "", "", "not transported (Table 1)"),
    ("OAT3", "PFBS", "", "", "", "", "not transported (Table 1)"),
    ("OAT2", "PFHpA, PFOA, PFNA, PFBS, PFHxS, PFOS", "", "", "", "",
     "OAT2 transported NONE of the six compounds tested (Table 1)"),
]


def build(fields):
    out = []
    for tr, ch, km, se, vmax, eff, note in ROWS:
        nc, head = NCF.get(ch, ("4-8", "carboxylate and sulfonate"))
        carried = note.startswith("NOT MEASURED HERE")
        r = {f: "" for f in fields}
        r.update({
            "study": STUDY, "year": "2024", "pmid": "39343156", "doi": DOI,
            "transporter": tr, "gene": GENE[tr], "species_orthologue": "human",
            "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "expression_system": SYSTEM,
            "measurement_type": (
                "Km/Vmax (Michaelis-Menten) and transport efficiency Vmax/Km"
                if km else
                "Km/Vmax tabulated from Louisse 2023, not measured here"
                if carried else
                "tested for saturable uptake; not transported"),
            "km": f"{km} +/- {se}" if km else "",
            "km_units": "uM" if km else "",
            "vmax": vmax, "vmax_units": "nmol/min/mg protein" if vmax else "",
            "uptake_ratio": "",
            "substrate_or_inhibitor": (
                "substrate" if (km or carried) else "not a substrate"),
            "direction_reabsorption_or_secretion": DIRN[tr],
            "sex_or_hormone_effect": "not assessed",
            "assay_volume_note": note,
            "source_table": (
                "Table 1" + (f" | transport efficiency {eff} uL/min/mg protein"
                             if eff else "")),
        })
        out.append(r)
    return out


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        fields, rows = rdr.fieldnames, list(rdr)
    if "assay_volume_note" not in fields:
        fields = list(fields) + ["assay_volume_note"]
        for r in rows:
            r.setdefault("assay_volume_note", "")
    old = [r for r in rows if (r.get("doi") or "") == DOI]
    new = build(fields)
    if old and len(old) == len(new) and all(
            a.get(k) == b.get(k) for a, b in zip(old, new) for k in fields):
        print("no change needed (Louisse 2024 Table 1 already ingested)")
        return 0
    merged, placed = [], False
    for r in rows:
        if (r.get("doi") or "") == DOI:
            if not placed:
                merged.extend(new); placed = True
            continue
        merged.append(r)
    if not placed:
        merged.extend(new)
    shutil.copy2(DB, DB + ".bak")
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(merged)
    print(f"Louisse 2024: {len(old)} qualitative row(s) -> {len(new)} "
          f"per-compound rows")
    print("    6 measured here with Km + Vmax (2 OAT1, 4 OAT3)")
    print("    4 carried over from Louisse 2023, stored WITHOUT a constant")
    print("    5 non-substrates + 1 OAT2 row (nothing transported)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
