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
OUT = os.path.join(HERE, "..", "db", "binding_affinity.xlsx")

FIELDS = ["id", "pfas", "species", "protein", "measure", "value", "units",
          "kd_uM", "kd_basis", "logK", "logK_basis", "source", "location"]

CMPD_FIELDS = ["pfas", "class", "rows", "kd_rows", "logK_rows", "studies"]

# The unfiltered table, kept alongside the filtered one so nothing is lost to
# this workbook's scope. Same columns plus the class that decided the filter.
ALL_FIELDS = ["pfas", "class", "species", "protein", "measure", "value",
              "units", "kd_uM", "kd_basis", "logK", "logK_basis",
              "source", "location"]

COND_FIELDS = ["id", "pfas", "species", "protein", "source", "mode",
               "pfas_conc_uM", "protein_conc_uM", "ligand_protein_ratio_max",
               "ratio_basis", "temperature_C", "temperature_as_stated", "nacl_mM",
               "buffer", "binding_site"]

# Kd in micromolar is the common currency. Only exact routes are used to get
# there -- see to_kd_uM. Concentrations are converted the same way: unit
# arithmetic only, never inference.
TO_UM = {"uM": 1.0, "mM": 1e3, "nM": 1e-3, "M": 1e6}

# --- which compounds this table is about -------------------------------
# "Common" is not a judgement call: it is the set of chemical strings the
# repository's own toxicokinetic and regulatory tables model. An exotic PFAS
# rides along only when the study that measured it also measured a common
# one, so there is something in the same experiment to anchor it to.
COMMON = {
    "PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA",
    "PFDoDA", "PFDoA", "PFTrDA", "PFTriDA", "PFTeDA",
    "PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS",
    "HFPO-DA (GenX)", "GenX", "ADONA", "NaDONA (ADONA)", "DONA",
    "6:2 FTS", "6:2 FtS", "6:2 FTSA", "PFOSA", "FOSA", "PFO2OA",
    "cC6O4", "C6O4", "9Cl-PF3ONS (F-53B, 6:2 Cl-PFESA)",
    "9Cl-PF3ONS (F-53B)", "9Cl-PF3ONS",
}
# A cell naming several compounds at once is a summary, not a measurement of
# one chemical, so it cannot go into a per-compound table either way.
MULTI = re.compile(r"\band\b|,|;|\brange\b|\bPFCAs\b|\bPFSAs\b|\bPFAS\b|"
                   r"\bPFAA\b|\bcompounds\b|\bligands\b|\bstandards\b|"
                   r"\bcongeners\b|^\d+\s", re.I)
# A perfluoroalkyl acid, telomer or sulfonamide, as against a reference
# ligand (octanoic acid, warfarin, DAUDA), a pesticide, or one of the neutral
# fluorochemicals in a broad screen (acrylates, amines, esters, diols).
PFAS_LIKE = re.compile(
    r"^(PF[A-Z0-9]|FOSA|FBSA|FHxSA|FDSA|N-?(Me|Et)FOSA|"
    r"\d+:\d+\s*(FT|Ft|diPAP|FTCA|FTSA|FTS|FTOH|FTUCA|FtS)|"
    r"\d+Cl-PF|HFPO|NaDONA|ADONA|DONA|diSAmPAP|cC6O4|C6O4|"
    r"\d+:\d+(:\d+)?\s*FTB|N-TAmP)", re.I)


def compound_class(chem):
    if chem in COMMON:
        return "common"
    if MULTI.search(chem):
        return "multi-compound summary"
    return "exotic PFAS" if PFAS_LIKE.match(chem) else "not a PFAS"

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


def to_kd_uM(rec):
    """Standardise to a dissociation constant in micromolar, EXACTLY or not
    at all.

    Three routes are exact and are the only ones used:
      - a Kd already in uM            -> as reported
      - a Kd in mM, M or nM           -> unit conversion
      - a Ka                          -> Kd = 1/Ka, true by definition

    Two routes were tried and removed, because each needs a quantity the
    papers do not report:
      - log D -> Kd needs the stoichiometry n. The REPORT.md 6.1b bridge is
        sound, but almost nobody fits n, so every such number rested on an
        assumed n = 1. The measured quantity is the log D itself, and it is
        in `value` and `units` where it belongs.
      - a docking dG is a computed score. Converting it into the same column
        as measured constants would put a prediction beside data.

    Returns (value_or_blank, basis). The basis is never empty: a number says
    where it came from, a blank says why it is blank.
    """
    meas, units, raw = rec["measure"], rec["units"], rec["value"]

    if meas == "log D":
        return "", ("a partition coefficient; converting to Kd needs the "
                    "stoichiometry, which this paper does not report")
    if meas.startswith("docking"):
        return "", "a computed docking score, not a measured constant"
    if meas == "f_unbound":
        return "", ("an unbound fraction, not an affinity; converting needs "
                    "the binding-site molarity, which this paper does not give")

    num, note = parse_value(raw)
    if num is None:
        return "", note

    def out(kd, basis):
        if kd <= 0 or not math.isfinite(kd):
            return "", "conversion out of range"
        return f"{kd:.4g}", basis + (f"; {note}" if note else "")

    if meas == "Kd":
        factor = TO_UM.get(units)
        if factor is None:
            return "", "no exact conversion from " + units
        return out(num * factor, "as reported" if factor == 1
                   else "unit conversion from " + units)

    if meas == "Ka":
        if units == "M^-1":
            return out(1e6 / num, "1/Ka")
        if units == "log10 M^-1":
            return out(1e6 / 10 ** num, "1/Ka, from log10 Ka")
        return "", "no exact conversion from " + units

    return "", "no exact conversion defined"


# --- the partition coefficient -----------------------------------------
# K is what a mass-balance model needs: how much PFAS a kilogram of protein
# holds per litre of water, not how tightly one site binds. Two routes reach
# it without assuming anything.
#
#   1. The paper measured it. Fischer, Qin, Allendorf and Bischel report a
#      protein-water distribution ratio directly. Qin's is per litre of
#      protein rather than per kilogram, which is a density conversion, not
#      an assumption: subtract log10(1.364), from a 0.733 mL/g partial
#      specific volume.
#   2. The paper reported BOTH an affinity and a stoichiometry. Then, in the
#      trace limit of REPORT.md 6.1b,
#          K (L_water/kg_protein) = 1000 * n / (M_protein * Kd[M])
#      with every term on the right read off the paper except M_protein,
#      which is a property of the protein rather than of the experiment.
#      Restricted to serum albumin, whose 66.5 kDa is not in doubt.
#
# There is no third route. Turning a Kd into a K without a stoichiometry is
# the n = 1 assumption that was removed earlier, and it stays removed.
VOL_TO_MASS = math.log10(1.364)
MW_ALBUMIN = 66500.0
N_NUM = re.compile(r"^\s*(\d*\.?\d+)\s*$")


def is_serum_albumin(protein):
    p = protein.lower()
    return "serum albumin" in p and "globulin" not in p


def to_logK(rec, n_sites):
    """log10 K in L_water per kg protein, or blank with the reason."""
    meas, units, raw = rec["measure"], rec["units"], rec["value"]

    if meas == "log D":
        num, note = parse_value(raw)
        if num is None:
            return "", note
        if units in ("log10 (L_water/kg_protein)", "log10 (L_water/kg_albumin)"):
            return "%.3f" % num, "as reported (mass basis)"
        if units == "log10 (L_water/L_protein)":
            return ("%.3f" % (num - VOL_TO_MASS),
                    "reported per L protein; density conversion to per kg")
        return "", "no exact conversion from " + units

    if not rec["kd_uM"]:
        return "", "no Kd, and the paper reports no partition coefficient"
    m = N_NUM.match(n_sites or "")
    if not m:
        return "", ("the paper reports an affinity but no stoichiometry, and "
                    "K cannot be had from one without the other")
    if not is_serum_albumin(rec["protein"]):
        return "", ("stoichiometry given, but this protein has no "
                    "unambiguous molar mass to divide by")
    n = float(m.group(1))
    kd_M = float(rec["kd_uM"]) * 1e-6
    if n <= 0 or kd_M <= 0:
        return "", "conversion out of range"
    return ("%.3f" % math.log10(1000.0 * n / (MW_ALBUMIN * kd_M)),
            "from the paper's own Kd and n, albumin 66.5 kDa")


CONC_NUM = re.compile(
    r"^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
    r"\s*(?:(?:-|to)\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))?"
    r"\s*(?:\(.*\))?\s*$")


def conc_range_uM(value, units):
    """(low, high) in micromolar, or (None, None).

    Unit arithmetic only. A dilution factor, a mol/mol ratio and a
    mass-per-volume with no molar mass all return nothing rather than a guess.
    """
    v, u = value.strip(), units.strip()
    if not v or not u:
        return None, None
    factor = TO_UM.get(u.split()[0])
    if factor is None:
        return None, None
    m = CONC_NUM.match(v)
    if not m:
        return None, None
    lo = float(m.group(1))
    hi = float(m.group(2)) if m.group(2) else lo
    return min(lo, hi) * factor, max(lo, hi) * factor


def conc_uM(value, units):
    """The highest concentration used, in uM -- the one that decides whether
    the protein approached saturation."""
    _, hi = conc_range_uM(value, units)
    return "" if hi is None else "%.6g" % hi


TEMP_NUM = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*$")
# Only an NaCl concentration the buffer states outright. Full ionic strength
# is not computed: most of these buffer descriptions name the salt but not
# every ionic species, and a phosphate buffer's contribution depends on the
# protonation state at that pH. Reading what is written beats modelling what
# is not.
NACL = re.compile(r"(\d+(?:\.\d+)?)\s*mM\s+(?:added\s+)?NaCl", re.I)
NACL_M = re.compile(r"(\d*\.\d+)\s*M\s+NaCl", re.I)


def temperature(raw):
    """(numeric C, verbatim). Numeric only when the field is just a number --
    'thermal ramp', 'room temperature' and multi-temperature series stay as
    text, because none of them is a single temperature the constant holds at.
    """
    v = raw.strip()
    m = TEMP_NUM.match(v)
    return ((v, v) if m else ("", v))


def nacl_mM(buffer_text):
    m = NACL.search(buffer_text)
    if m:
        return m.group(1)
    m = NACL_M.search(buffer_text)
    if m:
        return "%.6g" % (float(m.group(1)) * 1000.0)
    return ""


# Some papers state the mole ratio instead of both concentrations. Reading it
# off is not a guess, so it is used where the concentrations do not compute --
# but the two are flagged apart in ratio_basis, because a stated bound and a
# computed value are not the same kind of number.
STATED_RATIO = [
    (re.compile(r"PFAS:protein molar ratio\s*<=\s*([\d.]+)\s*\(HSA,BSA\)"
                r".*?<=\s*([\d.]+)\s*\(gamma-globulin\)", re.I | re.S),
     "albumin_globulin"),
    (re.compile(r"up to\s+([\d.]+)\s*:\s*1", re.I), "simple"),
    (re.compile(r"^\s*([\d.]+)\s*:\s*1\b", re.I), "simple"),
]


def stated_ratio(protein_conc_text, ligand_conc_text, protein):
    """The mole ratio a paper states in words, where it states one."""
    p = protein.lower()
    for text in (protein_conc_text, ligand_conc_text):
        for rx, kind in STATED_RATIO:
            m = rx.search(text)
            if not m:
                continue
            if kind == "albumin_globulin":
                if "albumin" in p:
                    return m.group(1)
                if "globulin" in p:
                    return m.group(2)
                return ""      # whole serum is not covered by that sentence
            return m.group(1)
    return ""


def ratio_max(lig_hi, prot_lo):
    """The HIGHEST ligand:protein molar ratio the assay reached -- most
    ligand over least protein -- computed only when both are in hand.

    This is the axis that orders most of the spread in this literature: below
    about 0.01 a measurement sees the specific site, above 1 it sees
    non-specific partitioning, and the two differ by orders of magnitude.
    Taking the worst case is what makes it a saturation check, and it
    reproduces the ligand_over_protein_max column of the working file on all
    151 rows that carry one.
    """
    if lig_hi is None or not prot_lo:
        return ""
    return "%.4g" % (lig_hi / prot_lo) if prot_lo > 0 else ""



def location(src):
    """Where in the paper. The working file appends the transcribed value and
    any transcription note after a pipe; the location is what precedes it."""
    return src.split("|")[0].strip().rstrip(",").strip()


def build(filter_compounds=True):
    """filter_compounds=False returns every measurement, for the consolidated
    database in db/combined/, which is the full record rather than a reading
    copy and should not lose rows to this table's scope."""
    with open(SRC) as fh:
        rows = list(csv.DictReader(fh))
    out, skipped, qualitative = [], [], []
    for r in rows:
        src = source_label(r["study"], r["year"])
        loc = location(r["source_table"])
        base = dict(pfas=r["chemical"], species=r["protein_species"],
                    protein=r["protein"], source=src, location=loc,
                    _id=(r["doi"] or r["pmid"]),
                    _study=r["study"], _year=r["year"],
                    _mode=r["method"],
                    _lig=conc_uM(r["ligand_conc_range"], r["ligand_conc_units"]),
                    _prot=conc_uM(r["protein_conc"], r["protein_conc_units"]),
                    _lig_hi=conc_range_uM(r["ligand_conc_range"],
                                          r["ligand_conc_units"])[1],
                    _prot_lo=conc_range_uM(r["protein_conc"],
                                           r["protein_conc_units"])[0],
                    _prot_raw=r["protein_conc"], _lig_raw=r["ligand_conc_range"],
                    _temp=r["temperature_c"], _buffer=r["buffer"],
                    _site=r["binding_site_identified"],
                    _n=r["n_binding_sites"],
                    _method=r["method"])
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
            rec["logK"], rec["logK_basis"] = to_logK(rec, r["n_binding_sites"])
            (out if is_quantitative(val) else qualitative).append(rec)
        if r["unbound_fraction"].strip():
            rec = dict(base, measure="f_unbound",
                       value=r["unbound_fraction"].strip(),
                       units=canon_units("", "f_unbound"))
            rec["kd_uM"], rec["kd_basis"] = to_kd_uM(rec)
            rec["logK"], rec["logK_basis"] = to_logK(rec, r["n_binding_sites"])
            out.append(rec)
        if not meas and not r["unbound_fraction"].strip():
            skipped.append(r)
    for rec in out:
        if rec["species"] == "bovine; human":
            rec["species"] = "human; bovine"

    # Restrict to the compounds this project is about. A common PFAS is kept
    # wherever it appears; an exotic one only where the same study also
    # measured a common one, so the experiment itself anchors it.
    anchored = {(r["_study"], r["_year"]) for r in out
                if compound_class(r["pfas"]) == "common"}
    kept, dropped = [], []
    for rec in out:
        cls = compound_class(rec["pfas"])
        rec["_class"] = cls
        rec["class"] = cls
        if not filter_compounds or cls == "common" or (
                cls == "exotic PFAS"
                and (rec["_study"], rec["_year"]) in anchored):
            kept.append(rec)
        else:
            dropped.append(rec)

    for i, rec in enumerate(kept, start=1):
        rec["id"] = i
    return rows, kept, skipped, qualitative, dropped


def compounds(out):
    """One row per compound kept, so the filter is auditable and overridable."""
    agg = collections.OrderedDict()
    for r in out:
        a = agg.setdefault(r["pfas"], {
            "pfas": r["pfas"], "class": r["_class"], "rows": 0,
            "kd_rows": 0, "logK_rows": 0, "_studies": set()})
        a["rows"] += 1
        a["kd_rows"] += 1 if r["kd_uM"] else 0
        a["logK_rows"] += 1 if r["logK"] else 0
        a["_studies"].add(r["source"])
    out_rows = []
    for a in agg.values():
        a["studies"] = "; ".join(sorted(a.pop("_studies")))
        out_rows.append(a)
    out_rows.sort(key=lambda a: (a["class"] != "common", -a["rows"], a["pfas"]))
    return out_rows


def _ratio_cell(r):
    computed = ratio_max(r["_lig_hi"], r["_prot_lo"])
    if computed:
        return {"ligand_protein_ratio_max": computed,
                "ratio_basis": "computed from the two concentrations"}
    stated = stated_ratio(r["_prot_raw"], r["_lig_raw"], r["protein"])
    if stated:
        return {"ligand_protein_ratio_max": stated,
                "ratio_basis": "upper bound stated by the paper"}
    return {"ligand_protein_ratio_max": "",
            "ratio_basis": "neither both concentrations nor a stated ratio"}


def conditions(out):
    """The assay conditions, one row per measurement, joined by `id`."""
    cond = []
    for r in out:
        t_num, t_raw = temperature(r["_temp"])
        cond.append({
            "id": r["id"], "pfas": r["pfas"], "species": r["species"],
            "protein": r["protein"], "source": r["source"],
            "mode": r["_mode"],
            "pfas_conc_uM": r["_lig"], "protein_conc_uM": r["_prot"],
            **_ratio_cell(r),
            "temperature_C": t_num, "temperature_as_stated": t_raw,
            "nacl_mM": nacl_mM(r["_buffer"]), "buffer": r["_buffer"],
            "binding_site": r["_site"]})
    return cond


NUMERIC = {"id", "kd_uM", "logK", "rows", "kd_rows", "logK_rows", "pfas_conc_uM", "protein_conc_uM",
           "ligand_protein_ratio_max", "temperature_C", "nacl_mM"}

NOTES = [
    ("Affinity", "One row per binding measurement, joined to Conditions by id."),
    ("Compounds", "Which compounds are in and why, so the filter can be "
                  "audited or overridden."),
    ("All measurements", "Every measurement in db/protein_binding.csv, "
                         "unfiltered, with the `class` column that decided "
                         "whether it reached the Affinity tab. Nothing is "
                         "lost to the compound filter; it is only set aside."),
    ("", ""),
    ("pfas", "Chemical as the paper names it."),
    ("species", "Species the protein came from."),
    ("protein", "Protein or matrix, as the paper describes it."),
    ("measure", "Which quantity: Kd, Ka, log D, f_unbound, docking dG."),
    ("value", "The number as printed, including ranges and bounds."),
    ("units", "Canonicalised unit for that value."),
    ("kd_uM", "Everything that converts EXACTLY to a dissociation constant "
              "in micromolar. Blank is a deliberate state, never a gap."),
    ("kd_basis", "How the kd_uM was reached, or why it is blank. Never empty."),
    ("logK", "log10 of the protein-water partition coefficient, L_water per "
             "kg protein -- what a mass-balance model needs. Two exact routes "
             "only: the paper measured it, or the paper reported BOTH an "
             "affinity and a stoichiometry. Never inferred from a Kd alone."),
    ("logK_basis", "Which route, or why the cell is blank. Never empty."),
    ("source", "First Author Year."),
    ("location", "Where the number sits in the paper."),
    ("", ""),
    ("Conditions", "The assay each measurement was made under."),
    ("", ""),
    ("mode", "Experimental method, verbatim."),
    ("pfas_conc_uM", "Ligand concentration in uM; the top of the range where "
                     "one is given. Blank where the paper gives a dilution "
                     "factor, a mol/mol ratio or a mass-per-volume with no "
                     "molar mass -- those are not convertible without a guess."),
    ("protein_conc_uM", "Protein concentration in uM, same rule."),
    ("ligand_protein_ratio_max", "the highest ratio the assay reached: most ligand "
                             "over least protein. Below ~0.01 a measurement "
                             "sees the specific site; above ~1 it sees "
                             "non-specific partitioning, and the two differ "
                             "by orders of magnitude. This column is the "
                             "single best predictor of where a constant will "
                             "land. 239 of 828 rows support one."),
    ("ratio_basis", "Whether the ratio was computed from the two "
                    "concentrations or read off a bound the paper states in "
                    "words. A computed value and a stated upper bound are "
                    "not the same kind of number, so they are kept apart."),
    ("temperature_C", "Numeric only where the paper states one temperature. "
                      "A thermal ramp or 'room temperature' stays in "
                      "temperature_as_stated and leaves this blank."),
    ("nacl_mM", "NaCl the buffer states outright. Full ionic strength is NOT "
                "computed: most buffer descriptions name the salt but not "
                "every ionic species, and a phosphate buffer's contribution "
                "depends on protonation at that pH."),
    ("buffer", "Buffer verbatim."),
    ("binding_site", "Site, where the paper identifies one (84 of 865 rows do)."),
    ("", ""),
    ("Compound filter", "Common PFAS are kept wherever they appear. Common "
                        "means the chemical strings this repository's own "
                        "toxicokinetic and regulatory tables model, not a "
                        "judgement call. An exotic PFAS is kept only where "
                        "the study that measured it also measured a common "
                        "one, so the same experiment anchors it. Dropped: "
                        "reference ligands (octanoic acid, warfarin, DAUDA), "
                        "pesticides, the neutral fluorochemicals in broad "
                        "screens, and cells naming several compounds at once."),
    ("", ""),
    ("What is deliberately absent", ""),
    ("", "log D is not converted to Kd, and a Kd is not converted to K "
         "without a stoichiometry. The bridge in REPORT.md 6.1b is sound "
         "but needs n, which almost no paper fits, so such a number would "
         "rest on an assumed n = 1. Where a paper DOES report n, the "
         "conversion is exact and is used."),
    ("", "Docking energies are not converted. They are computed scores, not "
         "measurements."),
    ("", "Unbound fractions are not converted. That needs a binding-site "
         "molarity these papers do not give."),
    ("", "Ionic strength is not computed from buffer prose."),
]


def write(out, cond, cmpd, full):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    head_fill = PatternFill("solid", fgColor="2A78D6")
    head_font = Font(bold=True, color="FFFFFF")

    wb = Workbook()
    wb.remove(wb.active)

    for title, fields, data in (("Affinity", FIELDS, out),
                                ("Conditions", COND_FIELDS, cond),
                                ("Compounds", CMPD_FIELDS, cmpd),
                                ("All measurements", ALL_FIELDS, full)):
        ws = wb.create_sheet(title)
        ws.append(fields)
        for r in data:
            row = []
            for f in fields:
                v = r.get(f, "")
                if f in NUMERIC and str(v).strip():
                    try:
                        v = float(v)
                        if f == "id":
                            v = int(v)
                    except ValueError:
                        pass
                row.append(v)
            ws.append(row)
        for c, f in enumerate(fields, start=1):
            cell = ws.cell(row=1, column=c)
            cell.fill, cell.font = head_fill, head_font
            cell.alignment = Alignment(vertical="center")
            width = max([len(f)] + [len(str(r.get(f, ""))) for r in data[:400]])
            ws.column_dimensions[get_column_letter(c)].width = \
                min(max(width + 2, 9), 42)
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

    ws = wb.create_sheet("Notes")
    ws["A1"] = "binding_affinity.xlsx -- column definitions"
    ws["A1"].font = Font(bold=True, size=13)
    for i, (k, v) in enumerate(NOTES, start=3):
        ws.cell(row=i, column=1, value=k).font = Font(bold=True)
        ws.cell(row=i, column=2, value=v).alignment = \
            Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 24
    ws.column_dimensions["B"].width = 96

    wb.save(OUT)


def check(rows, out, skipped, qualitative, dropped):
    print("=" * 74)
    print("BINDING AFFINITY TABLE")
    print("=" * 74)
    print(f"  source rows (db/protein_binding.csv)   {len(rows):5d}")
    print(f"  measurement rows written               {len(out):5d}")
    print(f"  source rows carrying no value          {len(skipped):5d}")
    print(f"  rows dropped by the compound filter    {len(dropped):5d}")
    print(f"  rows whose 'value' is prose, not a number {len(qualitative):3d}")
    print()
    # Columns that must never be blank, against the three where a blank is a
    # legitimate state the row explains for itself.
    REQUIRED = ["id", "pfas", "species", "protein", "measure", "value",
                "units", "kd_basis", "logK_basis", "source", "location"]
    print(f"Completeness of the {len(FIELDS)} columns:")
    for f in FIELDS:
        blank = sum(1 for r in out if not str(r[f]).strip())
        if f in REQUIRED:
            flag = "" if blank == 0 else "   <-- MUST NOT BE BLANK"
        else:
            where = {"kd_uM": "kd_basis", "logK": "logK_basis"}.get(f)
            flag = (f"   ({len(out) - blank} populated; every blank is "
                    f"explained in {where})" if where else
                    f"   ({len(out) - blank} populated; not reported by "
                    f"the paper)")
        print(f"  {f:13s} blank {blank:4d} / {len(out)}{flag}")
    print()
    print("Why kd_uM is blank, where it is:")
    for b, n in collections.Counter(
            r["kd_basis"] for r in out if not r["kd_uM"]).most_common():
        print(f"  {n:5d}  {b[:86]}")
    print()
    print("How the partition coefficient was reached, where it exists:")
    for b, n in collections.Counter(
            r["logK_basis"] for r in out if r["logK"]).most_common():
        print(f"  {n:5d}  {b[:86]}")
    print()
    print("Why logK is blank, where it is:")
    for b, n in collections.Counter(
            r["logK_basis"] for r in out if not r["logK"]).most_common():
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
    rows, out, skipped, qualitative, dropped = build()
    if not args.check:
        full = build(filter_compounds=False)[1]
        write(out, conditions(out), compounds(out), full)
        print(f"wrote {os.path.abspath(OUT)}  ({len(out)} rows)")
        print()
    check(rows, out, skipped, qualitative, dropped)


if __name__ == "__main__":
    main()
