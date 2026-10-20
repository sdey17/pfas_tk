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
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "db", "binding_affinity.csv")

FIELDS = ["pfas", "species", "protein", "measure", "value", "units",
          "kd_uM", "kd_basis", "pfas_conc", "protein_conc",
          "source", "location"]

# --- putting everything on one scale -----------------------------------
# Kd in micromolar is the common currency: it is what most of this
# literature reports, and the rest converts into it exactly except for the
# partition coefficients, which need a stoichiometry and a molar mass.
#
# The bridge is the one derived in REPORT.md 6.1b. From Langmuir occupancy,
# dividing by sorbent mass makes the protein concentration cancel, and in
# the trace limit
#
#     log D (L_water/kg_protein) = log Ka + log n - log10(M_protein / 1000)
#
# so  log Ka = log D - log n + log10(M_protein / 1000).  n is the number of
# binding sites, which most papers do not fit, so n = 1 is assumed and said
# so in kd_basis on every row it touches.
MW_PROTEIN = {          # g/mol
    "albumin": 66500.0,     # the value the 6.1b derivation uses
    "globulin": 150000.0,   # IgG
}
# D on a volume basis (L_water per L_protein) is the mass basis times the
# protein density, 1.364 kg/L from a 0.733 mL/g partial specific volume.
VOL_TO_MASS = math.log10(1.364)
RT_KCAL = 0.5925        # kcal/mol at 298 K, for docking dG -> Ka

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
        if u == "log M^-1":
            return "log10 M^-1"
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
        if kau.strip() == "log M^-1":
            return "Ka", "ka"          # a log association constant, not a log D
        if any(t in kau for t in ("log D", "log K_PW", "log L_water",
                                  "L_water/kg_membrane")):
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



# --- the conversion ---------------------------------------------------

NUM = re.compile(r"-?\d*\.?\d+(?:[eE][-+]?\d+)?")


def parse_value(v):
    """One number from the printed value, plus how it was obtained.

    Most values are a bare number. Ranges ("0.52e5-9.90e5", "8e-5 to 4e-4")
    are reduced to their geometric midpoint and said so; bounds ("<0.10") and
    approximations ("~1e4") are taken at face value and said so. A value that
    names two sites ("19.0 (KD1); 84.2 (KD2)") is not reduced, because the two
    are different sites rather than a spread.
    """
    v = v.strip()
    if v.lower() in WORD_VALUES:
        return None, "no binding detected"
    if ";" in v or "class 1" in v or "high-affinity set" in v:
        return None, "two distinct sites reported, not one value"
    nums = [float(x) for x in NUM.findall(v)]
    if not nums:
        return None, "no number"
    bare = re.fullmatch(r"-?\d*\.?\d+(?:[eE][-+]?\d+)?", v)
    if bare:
        return nums[0], ""
    if v.startswith("~"):
        return nums[0], "approximate as printed"
    if v.startswith("<") or v.startswith(">"):
        return nums[0], "a bound, not a point estimate"
    # a strict number on each side, so the '-' inside an exponent ("1.2e-4")
    # is not mistaken for the range separator
    n = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    rng = re.fullmatch(rf"\s*({n})\s*(?:-|to|\u2013)\s*({n})\s*", v)
    if rng:
        try:
            a, b = float(rng.group(1)), float(rng.group(2))
            if a > 0 and b > 0:
                return math.sqrt(a * b), "geometric midpoint of the reported range"
        except ValueError:
            pass
    # a number with a parenthetical (uncertainty, derivation) after it
    if re.match(r"^-?[\d.eE+-]+\s*[(+]", v):
        return nums[0], "leading value; qualifier dropped"
    return None, "value is not reducible to one number"


def protein_mw(protein):
    """Molar mass for the partition-coefficient bridge, or None when the
    sorbent has no single one (whole serum, plasma, a membrane)."""
    p = protein.lower()
    if "globulin" in p and "corticosteroid" not in p and "thyroxine" not in p:
        return MW_PROTEIN["globulin"], "gamma-globulin, 150 kDa"
    if "albumin" in p:
        return MW_PROTEIN["albumin"], "albumin, 66.5 kDa"
    return None, None


def to_kd_uM(rec):
    """Standardise to a dissociation constant in micromolar.

    Returns (value_or_blank, basis). The basis is never empty: a converted
    number always says what it came from, and a blank value always says why.
    """
    meas, units, raw = rec["measure"], rec["units"], rec["value"]
    num, note = parse_value(raw)
    if num is None:
        return "", note

    def out(kd, basis):
        if kd <= 0 or not math.isfinite(kd):
            return "", "conversion out of range"
        return f"{kd:.4g}", basis + (f"; {note}" if note else "")

    if meas == "Kd":
        factor = {"uM": 1.0, "mM": 1e3, "nM": 1e-3, "M": 1e6}.get(units)
        if factor is None:
            return "", f"no conversion defined for {units}"
        return out(num * factor, "as reported" if factor == 1
                   else f"unit conversion from {units}")

    if meas == "Ka":
        if units == "M^-1":
            return out(1e6 / num, "1/Ka")
        if units == "log10 M^-1":
            return out(1e6 / 10 ** num, "1/Ka, from log10 Ka")
        return "", f"no conversion defined for {units}"

    if meas == "log D":
        mw, label = protein_mw(rec["protein"])
        if mw is None:
            return "", ("no single molar mass for this sorbent "
                        "(whole serum, plasma or membrane)")
        logd = num
        basis = f"from log D via REPORT 6.1b, n=1 assumed, {label}"
        if units == "log10 (L_water/L_protein)":
            logd -= VOL_TO_MASS
            basis += "; volume basis converted to mass basis"
        elif units not in ("log10 (L_water/kg_protein)",
                           "log10 (L_water/kg_albumin)"):
            return "", f"no conversion defined for {units}"
        return out(1e6 / 10 ** (logd + math.log10(mw / 1000.0)), basis)

    if meas.startswith("docking"):
        return out(1e6 / math.exp(-num / RT_KCAL),
                   "from docking dG at 298 K (COMPUTATIONAL, not measured)")

    if meas == "f_unbound":
        return "", ("an unbound fraction, not an affinity; converting needs "
                    "the binding-site molarity, which these papers do not give")
    return "", "no conversion defined"


def conc(value, units):
    """Concentration as one readable string: the units here are heterogeneous
    (uM, mM, a dilution factor, a mol/mol ratio) and do not normalise."""
    v, u = value.strip(), units.strip()
    if not v:
        return ""
    if not u or u.lower() in v.lower():
        return v
    return f"{v} {u}"


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
                    _study=r["study"], _year=r["year"],
                    pfas_conc=conc(r["ligand_conc_range"], r["ligand_conc_units"]),
                    protein_conc=conc(r["protein_conc"], r["protein_conc_units"]))
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
            rec["kd_uM"], rec["kd_basis"] = to_kd_uM(rec)
            (out if is_quantitative(val) else qualitative).append(rec)
        if r["unbound_fraction"].strip():
            rec = dict(base, measure="f_unbound",
                       value=r["unbound_fraction"].strip(),
                       units=canon_units("", "f_unbound"))
            rec["kd_uM"], rec["kd_basis"] = to_kd_uM(rec)
            out.append(rec)
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
    # Columns that must never be blank, against the three where a blank is a
    # legitimate state the row explains for itself.
    REQUIRED = ["pfas", "species", "protein", "measure", "value", "units",
                "kd_basis", "source", "location"]
    print(f"Completeness of the {len(FIELDS)} columns:")
    for f in FIELDS:
        blank = sum(1 for r in out if not str(r[f]).strip())
        if f in REQUIRED:
            flag = "" if blank == 0 else "   <-- MUST NOT BE BLANK"
        else:
            flag = f"   ({len(out) - blank} populated; every blank is " \
                   f"explained in kd_basis)" if f == "kd_uM" else \
                   f"   ({len(out) - blank} populated; not reported by the paper)"
        print(f"  {f:13s} blank {blank:4d} / {len(out)}{flag}")
    print()
    print("Why kd_uM is blank, where it is:")
    for b, n in collections.Counter(
            r["kd_basis"] for r in out if not r["kd_uM"]).most_common():
        print(f"  {n:5d}  {b[:86]}")
    print()
    print("How the standardised Kd was reached, where it exists:")
    for b, n in collections.Counter(
            r["kd_basis"].split(";")[0] for r in out if r["kd_uM"]).most_common():
        print(f"  {n:5d}  {b[:86]}")
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
    print("Conversion cross-check -- PFOA on serum albumin, standardised.")
    print("These four reach a Kd by four different routes; agreement between")
    print("pairs that REPORT.md 6.1b/11a predicts should agree is the test")
    print("that the log D bridge is being applied correctly:")
    print()
    want = {("Fischer 2024", "human", "log D"): "specific site, SPME depletion",
            ("Qin 2024", "bovine", "log D"): "specific term, BioSPME isotherm",
            ("Chen 2025", "human", "Kd"): "equilibrium dialysis",
            ("Starnes 2025", "human", "Ka"): "DSF, at the melting transition"}
    got = {}
    for r in out:
        k = (r["source"], r["species"], r["measure"])
        if r["pfas"] == "PFOA" and "albumin" in r["protein"].lower() \
                and k in want and k not in got and r["kd_uM"]:
            got[k] = float(r["kd_uM"])
            print(f"  {r['source']:14s} {r['value']:>8s} "
                  f"{r['units'][:26]:26s} -> {r['kd_uM']:>8s} uM   {want[k]}")
    f = got.get(("Fischer 2024", "human", "log D"))
    q = got.get(("Qin 2024", "bovine", "log D"))
    c = got.get(("Chen 2025", "human", "Kd"))
    print()
    if f and q:
        print(f"  Fischer / Qin specific  {max(f, q) / min(f, q):5.2f}x apart "
              f"-- 11a predicts these two match")
    if c and q:
        print(f"  Chen / Qin specific     {max(c, q) / min(c, q):5.2f}x apart "
              f"-- 11a predicts Chen sits with the NONSPECIFIC term instead")
    print("  (Starnes is three orders weaker, which is the known DSF offset:")
    print("   affinity measured at the protein's melting transition, not 37 C)")
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
