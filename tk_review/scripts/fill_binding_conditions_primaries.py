#!/usr/bin/env python3
"""
Fill assay-condition columns in db/protein_binding.csv from the primary texts
saved under tk_review/papers/, for Peng 2024, Wu 2024 and Chen 2020.

All values below are quoted from the saved full texts; nothing is inferred.

Peng et al. 2024, Toxics 12:43 (PMID 38250999)
  "PBS buffer ... was used to prepare HSA stock solution (1 x 10-6 mol/L) to
   ensure the maintenance of appropriate ionic strength (pH = 7.4)"
  "3 mL of HSA stock solution was put into a 10 mm square quartz cuvette.
   Following that, PFAS stock solution was gradually added to achieve various
   final concentrations, namely 0, 3, 6, 9, 12, 15, 18 x 10-6 mol/L"
  "the solution temperature ... was set at 298, 304, and 310 K"

Wu et al. 2024, Toxics 12:46 (PMID 38251003)
  "Each PFAS was prepared in anhydrous ethanol as 1 x 10-3 mol/L stock
   solution ... HSA ... was configured as 1 x 10-3 mol/L stock solution with
   Tris-HCl buffer (pH 7.4, 1 x 10-3 mol/L)"
  "5 mL of HSA stock solution was mixed with 1 mL, 3 mL, and 5 mL of PFASs
   stock solution individually to obtain the concentration ratios of 5:1, 5:3,
   and 5:5 ... at room temperature"
  -> final HSA 0.83 / 0.63 / 0.50 mM; final PFAS 0.17 / 0.38 / 0.50 mM.
  NOTE the co-solvent: the PFAS stock is anhydrous ethanol, so those same
  mixtures are 17%, 38% and 50% ethanol v/v. That is far above the <1% usually
  tolerated in a protein binding assay and is recorded as a caveat.

Chen et al. 2020, Int J Environ Res Public Health 17:1319 (PMID 32085632)
  Figure 1 caption: "[HSA] = 2 x 10-6 mol/L; Buffer: PBS, pH = 7.4", PFCA
   series "(10-6 mol/L) 1-10: 0, 2, 4, 8, 12, 16, 20, 30, 40 and 80"
  "fluorescence emission spectra at different temperatures (296, 303, and 310 K)"

Not filled, and why:
  Chen & Guo 2009 (8 rows) and Chen et al. 2015 (2 rows) - abstract only, no
  full text retrieved. Fischer et al. 2024 (111 rows) already complete.

Run:  python3 tk_review/scripts/fill_binding_conditions_primaries.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

FILL = {
    ("Peng et al.", "2024"): dict(
        protein_conc="1.0", protein_conc_units="uM",
        ligand_conc_range="0-18", ligand_conc_units="uM",
        temperature_c="25, 31, 37 (298, 304, 310 K)",
        buffer="PBS, pH 7.4",
        assay_volume_note="3 mL in a 10 mm quartz cuvette; fluorescence "
                          "quenching titration"),
    ("Wu et al.", "2024"): dict(
        protein_conc="0.50-0.83", protein_conc_units="mM",
        ligand_conc_range="0.17-0.50", ligand_conc_units="mM",
        temperature_c="room temperature (not stated numerically)",
        buffer="Tris-HCl 1 mM, pH 7.4",
        assay_volume_note="HSA 1 mM stock mixed 5:1, 5:3, 5:5 v/v with 1 mM "
                          "PFAS in ANHYDROUS ETHANOL, so the assay is 17-50% "
                          "ethanol v/v - far above the <1% normally tolerated "
                          "in albumin binding work; treat these constants with "
                          "caution"),
    ("Chen et al.", "2020"): dict(
        protein_conc="2.0", protein_conc_units="uM",
        ligand_conc_range="0-80", ligand_conc_units="uM",
        temperature_c="23, 30, 37 (296, 303, 310 K)",
        buffer="PBS, pH 7.4",
        assay_volume_note="1.0 cm quartz cell, lambda_ex 295 nm, 5.0 nm slits; "
                          "Table 2 Ka values are not in the saved flattened "
                          "text, so the Ka here remains the Starnes 2026 "
                          "secondary range"),
}


def _uM(value, units):
    f = {"mm": 1e3, "um": 1.0, "nm": 1e-3, "m": 1e6}.get((units or "").strip().lower())
    return None if f is None else value * f


def _ends(s):
    """lowest and highest number in a '0-18' / '0.17-0.50' style field"""
    s = (s or "").replace("~", "").replace("to", "-").strip()
    vals = []
    for p in s.split("-"):
        try:
            vals.append(float(p.strip()))
        except ValueError:
            pass
    return (min(vals), max(vals)) if vals else (None, None)


def main():
    rows = list(csv.DictReader(open(DB)))
    fields = list(rows[0].keys())
    n_filled = 0
    for r in rows:
        key = (r["study"], r["year"])
        if key in FILL and not (r.get("protein_conc") or "").strip():
            r.update(FILL[key])
            n_filled += 1

    # recompute the derived ratio for every row that now has both
    for r in rows:
        r["ligand_over_protein_max"] = ""
        p_lo, _ = _ends(r.get("protein_conc"))
        _, l_hi = _ends(r.get("ligand_conc_range"))
        p = _uM(p_lo, r.get("protein_conc_units")) if p_lo is not None else None
        l = _uM(l_hi, r.get("ligand_conc_units")) if l_hi is not None else None
        if p and l:
            r["ligand_over_protein_max"] = "%.2f" % (l / p)

    with open(DB, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

    print("filled %d rows from primary texts" % n_filled)
    n = len(rows)
    for c in ("protein_conc", "ligand_conc_range", "temperature_c", "buffer",
              "ligand_over_protein_max"):
        k = sum(1 for r in rows if (r.get(c) or "").strip())
        print("  %-24s %3d / %d  (%.0f%%)" % (c, k, n, 100 * k / n))


if __name__ == "__main__":
    main()
