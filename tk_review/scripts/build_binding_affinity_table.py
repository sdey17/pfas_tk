#!/usr/bin/env python3
"""The plain binding-affinity table: one measurement per row, eight columns.

db/protein_binding.csv is the working file -- 28 columns, assay conditions,
transcription notes, provenance.  This builds the reading copy from it:

    pfas, species, protein, measure, value, units, source, location

with units canonicalised (the working file carries the paper's own wording,
commentary and all) and the source written as "First Author Year".

It REPLACES the binding limb of build_combined_datasets.py, which read the
two aggregate summary files and silently produced 311 of 327 rows with
species, value, units, study and year all blank -- those files carry
ka_M-1_geomean and studies, not value and study, so every .get() missed.

Run:  python3 scripts/build_binding_affinity_table.py
      python3 scripts/build_binding_affinity_table.py --check
"""
import argparse
import collections
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "db", "binding_affinity.csv")

FIELDS = ["pfas", "species", "protein", "measure", "value", "units",
          "source", "location"]

# Surnames that are two words, so the "first author" is not just the first token.
MULTIWORD = {"Sanchez Garcia"}


def first_author(study):
    """'Fischer et al.' -> 'Fischer';  'Chen & Guo' -> 'Chen';  'EFSA' -> 'EFSA'."""
    s = study.replace(" et al.", "").strip()
    if " & " in s:
        s = s.split(" & ")[0].strip()
    return s


def source_label(study, year):
    return f"{first_author(study)} {year}".strip()


# --- units -----------------------------------------------------------------
# The working file stores the paper's own unit wording, often with the
# derivation or a confidence interval appended.  These strip that back to the
# unit itself; anything unmatched is reported by --check rather than guessed.

def canon_units(raw, measure):
    u = " ".join(raw.split())
    if measure == "log D":
        if "kg_protein" in u or u.startswith("log K_PW"):
            return "log10 (L_water/kg_protein)"
        if "L_protein" in u:
            return "log10 (L_water/L_protein)"
        if "kg_albumin" in u:
            return "log10 (L_water/kg_albumin)"
        if "kg_membrane" in u:
            return "L_water/kg_membrane"
        if u == "log M^-1":
            return "log10 M^-1"
        return u
    if measure == "Ka":
        if u.startswith("M^-1") or u.startswith("L/mol") or u.startswith("no binding"):
            return "M^-1"
        return u
    if measure == "Kd":
        for unit in ("uM", "mM", "nM", "M"):
            if u == unit or u.startswith(unit + " ("):
                return unit
        if u.startswith("kcal/mol"):
            return "kcal/mol"
        return u
    if measure.startswith("docking"):
        return "kcal/mol"
    if measure == "f_unbound":
        return "fraction (0-1)"
    return u


def measure_of(r):
    """Which quantity this row actually reports, and from which column.

    The working file stores a reciprocal alongside whichever of Ka/Kd the
    paper printed, so the unit string says which one is the measurement:
      ka_units 'computed as 1/K_D'  -> the paper reported Kd
      kd_units '1/Ka'               -> the paper reported Ka
    """
    ka, kd = r["ka"].strip(), r["kd"].strip()
    kau, kdu = r["ka_units"], r["kd_units"]
    if ka or "no binding detected" in kau:
        if any(t in kau for t in ("log D", "log K_PW", "log L_water",
                                  "L_water/kg_membrane", "log M^-1")):
            return "log D", "ka"
        if "computed as 1/K_D" in kau and kd:
            return "Kd", "kd"
        return "Ka", "ka"
    if kd:
        if "kcal/mol" in kdu:
            return "docking dG (computational)", "kd"
        return "Kd", "kd"
    return None, None


# A value is a measurement when it starts with a number, optionally behind a
# comparator -- that admits ranges ("0.52e5-9.90e5"), approximations ("~1e4"),
# bounds ("<0.10") and multi-site sets ("19.0 (KD1); 84.2 (KD2)"), all of which
# are values the paper reported. Anything else is prose and is reported by
# --check instead of being written with an empty unit beside it.
QUANT = re.compile(r"^[~<>=]?\s*-?\s*\.?\d")
LOG_PREFIX = re.compile(r"^log\s+\S+\s*=\s*(.+)$", re.I)
WORD_VALUES = {"nonbinding", "not detected", "no binding detected"}


def normalise_value(v):
    """Strip a 'log K_x = ' prefix; the unit column already records the log."""
    m = LOG_PREFIX.match(v.strip())
    return m.group(1).strip() if m else v.strip()


def is_quantitative(v):
    v = v.strip()
    return bool(QUANT.match(v)) or v.lower() in WORD_VALUES


def location(src):
    """Where in the paper. The working file appends the transcribed value and
    any transcription note after a pipe; the location is what precedes it."""
    return src.split("|")[0].strip().rstrip(",").strip()


def build():
    with open(SRC) as fh:
        rows = list(csv.DictReader(fh))
    out, skipped, qualitative = [], [], []
    for r in rows:
        src = source_label(r["study"], r["year"])
        loc = location(r["source_table"])
        base = dict(pfas=r["chemical"], species=r["protein_species"],
                    protein=r["protein"], source=src, location=loc,
                    _method=r["method"], _id=(r["doi"] or r["pmid"]),
                    _study=r["study"], _year=r["year"])
        meas, col = measure_of(r)
        if meas:
            raw_u = r["ka_units"] if col == "ka" else r["kd_units"]
            val = r[col].strip()
            if not val and "no binding detected" in r["ka_units"]:
                val = "not detected"
            val = normalise_value(val)
            units = canon_units(raw_u, meas)
            if not units and val.lower() in WORD_VALUES:
                # a non-detection has no unit to carry
                val, units = "not detected", "not applicable"
            rec = dict(base, measure=meas, value=val, units=units)
            (out if is_quantitative(val) else qualitative).append(rec)
        if r["unbound_fraction"].strip():
            out.append(dict(base, measure="f_unbound",
                            value=r["unbound_fraction"].strip(),
                            units=canon_units("", "f_unbound")))
        if not meas and not r["unbound_fraction"].strip():
            skipped.append(r)
    for rec in out:
        if rec["species"] == "bovine; human":
            rec["species"] = "human; bovine"
    return rows, out, skipped, qualitative


def write(out):
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)


def check(rows, out, skipped, qualitative):
    print("=" * 74)
    print("BINDING AFFINITY TABLE")
    print("=" * 74)
    print(f"  source rows (db/protein_binding.csv)   {len(rows):5d}")
    print(f"  measurement rows written               {len(out):5d}")
    print(f"  source rows carrying no value          {len(skipped):5d}")
    print(f"  rows whose 'value' is prose, not a number {len(qualitative):3d}")
    print()
    print("Completeness of the eight columns:")
    for f in FIELDS:
        blank = sum(1 for r in out if not str(r[f]).strip())
        flag = "" if blank == 0 else "   <-- GAP"
        print(f"  {f:10s} blank {blank:4d} / {len(out)}{flag}")
    print()
    print("Rows by measure:")
    for m, n in collections.Counter(r["measure"] for r in out).most_common():
        print(f"  {n:5d}  {m}")
    print()
    print("Units, after canonicalisation:")
    for u, n in collections.Counter(r["units"] for r in out).most_common():
        print(f"  {n:5d}  {u}")
    print()
    print("Species:")
    for s, n in collections.Counter(r["species"] for r in out).most_common():
        print(f"  {n:5d}  {s}")
    print()
    print(f"Sources: {len(set(r['source'] for r in out))} distinct")
    print("  " + ", ".join(sorted(set(r["source"] for r in out))))
    print()
    if qualitative:
        print("Reported qualitatively, so carrying no value+unit pair. These "
              "stay in\ndb/protein_binding.csv and are listed here so the "
              "exclusion is visible:")
        for r in qualitative:
            print(f"  {r['source']:16s} {r['pfas'][:24]:24s} "
                  f"{r['protein'][:30]:30s} {r['value'][:52]}")
        print()
    if skipped:
        print("Source rows with no value, by study (qualitative or narrative "
              "entries;\nthey stay in db/protein_binding.csv and are not "
              "measurements):")
        for k, n in collections.Counter(
                (r["study"], r["year"]) for r in skipped).most_common():
            print(f"  {n:4d}  {k[0]} {k[1]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report only, do not write")
    args = ap.parse_args()
    rows, out, skipped, qualitative = build()
    if not args.check:
        write(out)
        print(f"wrote {os.path.abspath(OUT)}  ({len(out)} rows)")
        print()
    check(rows, out, skipped, qualitative)


if __name__ == "__main__":
    main()
