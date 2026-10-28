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

# --- the major-study subset -------------------------------------------
# "Major" is a count, not a judgement: a study is major here when it supplies
# at least this many standardised measurements. Five do, and they are the
# five the review's quantitative sections actually rest on -- Qin and Fischer
# for the mass-balance limb, Chen and Muensterman for the partition
# coefficients, Ryu for the free fractions. The threshold is written down so
# the selection can be re-derived rather than taken on trust, and so it moves
# on its own as the database grows.
MAJOR_MIN_ROWS = 100


def major_sources(out):
    """Sources carrying at least MAJOR_MIN_ROWS rows, largest first."""
    n = collections.Counter(r["source"] for r in out)
    return [s for s, c in n.most_common() if c >= MAJOR_MIN_ROWS]


# The unfiltered table, kept alongside the filtered one so nothing is lost to
# this workbook's scope. Same columns plus the class that decided the filter.
ALL_FIELDS = ["pfas", "class", "species", "protein", "measure", "value",
              "units", "kd_uM", "kd_basis", "logK", "logK_basis",
              "source", "location"]

COND_FIELDS = ["id", "pfas", "species", "protein", "source", "mode",
               "pfas_conc_uM", "pfas_conc_as_stated",
               "protein_conc_uM", "protein_conc_as_stated",
               "ligand_protein_ratio_max",
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


# One (study, year) pair in the database spans two different papers: Qin et
# al. 2024 published the plasma-binding paper (ES&T 58:9954, 3c10824) and the
# baseline-toxicity paper (ES&T 58:5727, 3c09950) in the same volume. A
# reading copy whose `source` column maps both to "Qin 2024" is ambiguous
# exactly where this project needs it not to be, since one measured albumin
# and plasma and the other measured cells and assay medium.
DISAMBIGUATE = {
    "10.1021/acs.est.3c10824": "Qin 2024 (plasma)",
    "10.1021/acs.est.3c09950": "Qin 2024 (bioassay)",
}


def source_label(study, year, doi=""):
    if doi in DISAMBIGUATE:
        return DISAMBIGUATE[doi]
    return f"{first_author(study)} {year}".strip()


# --- units -----------------------------------------------------------------
# The working file stores the paper's own unit wording, often with the
# derivation or a confidence interval appended.  These strip that back to the
# unit itself; anything unmatched is reported by --check rather than guessed.

def canon_units(raw, measure):
    u = " ".join(raw.split())
    if measure == "log K_membrane":
        return ("log10 K_IAM (basis as published)" if u.startswith("log K_IAM")
                else "log10 K_MW (basis as published)")
    if measure == "log D":
        # A Freundlich fit has no single constant: the stored value is the
        # ratio at one aqueous concentration. That is a different kind of
        # number from a concentration-independent ratio and must not be
        # canonicalised into the same string as one.
        if "Freundlich-type fit" in u:
            phase = ("L_prot+lip" if "L_prot+lip" in u else
                     "L_lipid" if "L_lipid" in u else "L_protein")
            return f"log10 (L_water/{phase}) at Cw = 1 uM, Freundlich fit"
        if "L_prot+lip" in u:
            return "log10 (L_water/L_prot+lip)"
        if "L_lipid" in u or "L_lip)" in u:
            return "log10 (L_water/L_lipid)"
        if u.startswith("log K_PW (protein-water distribution coefficient), "
                        "L_water per L_albumin"):
            return "log10 (L_water/L_protein)"
        if u.startswith("log K_HSA"):
            return "log10 K_HSA (basis as published)"
        if "L_serum" in u:
            return "log10 (L_water/L_serum, whole matrix)"
        if "kg_albumin" in u or "kg_protein" in u or u.startswith("log K_PW"):
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
        if kau.startswith("log K_MW") or kau.startswith("log K_IAM"):
            return "log K_membrane", "ka"
        if kau.startswith("log K_HSA"):
            return "log D", "ka"
        if kau.strip() == "log M^-1":
            return "Ka", "ka"          # a log association constant, not a log D
        if any(t in kau for t in ("log D", "log K_PW", "log L_water",
                                  "log K_albumin/water",
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

    if meas in ("log D", "log K_membrane"):
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
MW_ALBUMIN = 66463.0      # was 66500 here and 66463 everywhere else in the
                          # project; 0.00024 log apart, but one number only
N_NUM = re.compile(r"^\s*(\d*\.?\d+)\s*(?:\(.*\))?\s*$")

# The bridge constant is log10(1e9/M) and it is NOT universal. Albumin and
# L-FABP are 0.670 log apart -- a factor of 4.7 -- so dividing an L-FABP Kd
# by albumin's mass is a silent, one-directional error. Masses are the
# mature chain computed from the UniProt sequence (average residue masses,
# water added, 2 H per disulfide removed), which is what a purified
# preparation weighs.
PROTEIN_MW = [
    # (matcher, species or None for any, M g/mol, label for the basis cell)
    (lambda p: "serum albumin" in p and "globulin" not in p, None,
     MW_ALBUMIN, "albumin 66.5 kDa"),
    (lambda p: "l-fabp" in p or "liver fatty acid" in p, "human",
     14208.0, "human L-FABP 14.2 kDa, UniProt P07148"),
    (lambda p: "l-fabp" in p or "liver fatty acid" in p, "rat",
     14272.0, "rat L-FABP 14.3 kDa, UniProt P02692"),
    (lambda p: "i-fabp" in p or "intestinal fatty acid" in p, "human",
     15106.0, "human I-FABP 15.1 kDa, UniProt P12104 mature chain"),
]

# Species difference inside albumin: HSA 66,437 and BSA 66,398 by sequence,
# 0.00025 log apart. One mass is used for both; the distinction is below
# every other uncertainty in this table by three orders.

PPAR_REFUSAL = (
    "stoichiometry given, but the construct mass is not: Khazaee et al. "
    "2021 names only the supplier and never says whether the material is "
    "the full-length receptor (52.2, 57.6 and 49.9 kDa for alpha, gamma "
    "and delta) or the ligand-binding domain alone, about a third of "
    "that. A factor of 3 in M is 0.48 log in K")


def is_serum_albumin(protein):
    p = protein.lower()
    return "serum albumin" in p and "globulin" not in p


def protein_mw(protein, species):
    """(M, label) for a protein this table may divide by, or (None, reason).

    A refusal here is a statement that the molar mass of the material in
    the assay is not known, which is a different failure from a missing
    stoichiometry and is reported as such.
    """
    p, sp = protein.lower(), (species or "").strip().lower()
    if "peroxisome" in p or "ppar" in p:
        return None, PPAR_REFUSAL
    for match, want_sp, m, label in PROTEIN_MW:
        if match(p) and (want_sp is None or want_sp == sp):
            return m, label
    return None, ("stoichiometry given, but this protein has no "
                  "unambiguous molar mass to divide by")


SITE_PAREN = re.compile(r"^\s*(\d*\.?\d+)\s*(?:\(.*\))?\s*$")


def site_count(n_sites):
    """(n, how it was stated) from the n_binding_sites cell, or (None, None).

    Three wordings occur. A bare number is a fit. A number with the method
    in brackets -- Woodcroft's "2 (by displacement)" -- is the same thing.
    A single-site model statement fixes n at 1 as the paper's own
    assumption, which is weaker than a fit and is labelled so.
    """
    v = (n_sites or "").strip()
    m = SITE_PAREN.match(v)
    if m:
        return float(m.group(1)), "fitted"
    if "one site assumed" in v or "single-site model assumed" in v:
        return 1.0, "single-site model assumed by the paper"
    return None, None


def kd_list_uM(rec, row=None):
    """Every dissociation constant on the row, in uM.

    Usually one, already parsed into kd_uM. A paper that fits several
    classes of site writes them into one cell -- Woodcroft et al. 2010
    gives "3.1e-6, 26.2e-6, 52.6e-6" in molar -- which parse_value
    rightly refuses as a single number and which is read here instead.
    """
    if rec["kd_uM"]:
        return [float(rec["kd_uM"])]
    if row is None or rec["measure"] != "Kd":
        return []
    txt = row.get("kd", "").strip()
    if "," not in txt:
        return []
    out = []
    for part in txt.split(","):
        tok = part.split()[0] if part.split() else ""
        try:
            out.append(float(tok))
        except ValueError:
            return []
    u = row.get("kd_units", "").strip().lower().split()
    scale = {"m": 1e6, "mm": 1e3, "um": 1.0, "nm": 1e-3}.get(u[0] if u else "")
    return [v * scale for v in out] if scale else []


def to_logK(rec, n_sites, row=None):
    """log10 K in L_water per kg protein, or blank with the reason."""
    meas, units, raw = rec["measure"], rec["units"], rec["value"]

    if meas == "log K_membrane":
        return "", ("a membrane-water coefficient, not protein-water; it is "
                    "the lipid limb and belongs in its own column")
    if meas == "log D":
        num, note = parse_value(raw)
        if num is None:
            return "", note
        if units == "log10 K_HSA (basis as published)":
            return ("%.3f" % num,
                    "as reported, biomimetic surrogate; basis not stated by "
                    "the paper, so not placed on the mass basis")
        if units == "log10 (L_water/L_serum, whole matrix)":
            return "", ("a whole-matrix coefficient, per L of serum rather than "
                        "per kg of one protein; it is what a mass balance has to "
                        "PREDICT, not an input to one")
        if units in ("log10 (L_water/kg_protein)", "log10 (L_water/kg_albumin)"):
            return "%.3f" % num, "as reported (mass basis)"
        if units == "log10 (L_water/L_protein)":
            return ("%.3f" % (num - VOL_TO_MASS),
                    "reported per L protein; density conversion to per kg")
        return "", "no exact conversion from " + units

    kds = kd_list_uM(rec, row)
    if not kds:
        return "", "no Kd, and the paper reports no partition coefficient"
    n, how = site_count(n_sites)
    if n is None:
        return "", ("the paper reports an affinity but no stoichiometry, and "
                    "K cannot be had from one without the other")
    mw, label = protein_mw(rec["protein"], rec["species"])
    if mw is None:
        return "", label
    if n <= 0 or min(kds) <= 0:
        return "", "conversion out of range"

    # Several classes of site add: K0 = (1e9/M) x SUM_i n_i/Kd_i. A paper
    # that fits k classes and reports k constants has already told us
    # n_i = 1 for each of them, which is what "n = 3 (by ITC)" means --
    # three classes, one site apiece -- so the count must match the list.
    if len(kds) > 1:
        if abs(n - len(kds)) > 1e-9:
            return "", (f"{len(kds)} dissociation constants but n = {n:g}; "
                        "the site count and the constants do not correspond, "
                        "so the sum cannot be formed")
        total = sum(1.0 / k for k in kds)
        detail = f"summed over {len(kds)} classes of site"
    else:
        total = n / kds[0]
        detail = "one class of site" if n == 1 else f"n = {n:g}"
    return ("%.3f" % math.log10(1e9 * total / mw),
            f"from the paper's own Kd and n, {label}, {detail}")


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


def as_stated(value, units):
    """The concentration verbatim, value and units joined.

    conc_uM refuses anything it cannot convert by unit arithmetic alone,
    which is correct but leaves a blank that cannot be told from a study
    reporting nothing. This keeps the source text beside it, exactly as
    temperature_as_stated does for temperature.
    """
    v, u = (value or "").strip(), (units or "").strip()
    if not v:
        return ""
    return f"{v} {u}".strip()


def conc_uM(value, units):
    """The highest concentration used, in uM -- the one that decides whether
    the protein approached saturation."""
    _, hi = conc_range_uM(value, units)
    return "" if hi is None else "%.6g" % hi


# Albumin is the one protein here whose molar mass is a published constant
# rather than a property of the preparation, so a mass concentration of it
# converts to a molarity exactly. Serum albumin is 66.5 kDa (BSA 66,463
# g/mol, HSA 66,472 -- the difference is 0.01% and far below the precision
# of any concentration stated in these papers). This is unit arithmetic with
# one constant, not an inference about the sample, and it is applied ONLY to
# albumin: gamma-globulin is a heterogeneous fraction with no single mass,
# serum and plasma are mixtures, and an immobilised stationary phase has no
# concentration at all. Those stay blank.
ALBUMIN_G_PER_MOL = 66463.0
# Numbers in this column are written in several styles, scientific notation
# included ("3e0 g/L", "2.5e-2 g/L"), so the number pattern has to admit it.
_NUM = r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"
MASS_CONC = re.compile(
    rf"({_NUM})\s*(?:-|to|\u2013)\s*({_NUM})\s*(mg/mL|g/L)"
    rf"|({_NUM})\s*(mg/mL|g/L)", re.I)


def albumin_conc_uM(protein, raw):
    """Highest albumin concentration in uM, from a mass concentration.

    mg/mL and g/L are the same number, both being grams per litre. `raw` is
    the value and its unit already joined, because the database keeps them
    in separate columns and neither alone is interpretable.
    """
    p = protein.lower()
    if "albumin" not in p or "stationary phase" in p:
        return ""
    m = MASS_CONC.search(raw)
    if not m:
        return ""
    hi = float(m.group(2) or m.group(4))          # g/L
    return "%.6g" % (hi / ALBUMIN_G_PER_MOL * 1e6)


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


def mark_computed(rec, r):
    """Say in the `measure` cell when a value was computed, not measured.

    The working file carries that in is_computational, which this table has
    no column for, and the reason usually sits behind a pipe in source_table
    where `location` cuts it off. A predicted value reaching a reading copy
    unlabelled is the one error this table must not make, so the label goes
    where it cannot be missed. Applied after the conversions, which switch on
    the bare measure name.
    """
    if not r["is_computational"].strip().lower().startswith("no") and \
            "computational" not in rec["measure"]:
        rec["measure"] += " (computed)"


def build(filter_compounds=True):
    """filter_compounds=False returns every measurement, for the consolidated
    database in db/combined/, which is the full record rather than a reading
    copy and should not lose rows to this table's scope."""
    with open(SRC) as fh:
        rows = list(csv.DictReader(fh))
    out, skipped, qualitative = [], [], []
    for r in rows:
        src = source_label(r["study"], r["year"], r["doi"])
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
                    _prot_units=r["protein_conc_units"],
                    _lig_units=r["ligand_conc_units"],
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
            rec["logK"], rec["logK_basis"] = to_logK(rec, r["n_binding_sites"], r)
            mark_computed(rec, r)
            (out if is_quantitative(val) else qualitative).append(rec)
        if r["unbound_fraction"].strip():
            rec = dict(base, measure="f_unbound",
                       value=r["unbound_fraction"].strip(),
                       units=canon_units("", "f_unbound"))
            rec["kd_uM"], rec["kd_basis"] = to_kd_uM(rec)
            rec["logK"], rec["logK_basis"] = to_logK(rec, r["n_binding_sites"], r)
            mark_computed(rec, r)
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
            "pfas_conc_uM": r["_lig"],
            "pfas_conc_as_stated": as_stated(r["_lig_raw"], r["_lig_units"]),
            "protein_conc_uM": (r["_prot"] or albumin_conc_uM(
                r["protein"], as_stated(r["_prot_raw"], r["_prot_units"]))),
            "protein_conc_as_stated": as_stated(r["_prot_raw"], r["_prot_units"]),
            **_ratio_cell(r),
            "temperature_C": t_num, "temperature_as_stated": t_raw,
            "nacl_mM": nacl_mM(r["_buffer"]), "buffer": r["_buffer"],
            "binding_site": r["_site"]})
    return cond


NUMERIC = {"id", "kd_uM", "logK", "rows", "kd_rows", "logK_rows", "pfas_conc_uM", "protein_conc_uM",
           "ligand_protein_ratio_max", "temperature_C", "nacl_mM"}

NOTES = [
    ("Affinity", "One row per binding measurement, joined to Conditions by id."),
    ("Major studies", "The same rows as Affinity, restricted to the "
                      "studies that supply at least 100 standardised "
                      "measurements each. The threshold is a count, not a "
                      "judgement, and is written into the script "
                      "(MAJOR_MIN_ROWS) so the selection can be re-derived "
                      "and moves on its own as the database grows; the "
                      "--check output lists who qualifies and with what. "
                      "Note that counting rows is not the same as ranking "
                      "usefulness: Ryu 2024 is the largest contributor here "
                      "and reports unbound fractions in whole tissue "
                      "homogenate, which says nothing about any individual "
                      "protein. db/binding_major_sources.xlsx is the "
                      "companion selection made on that stricter basis, "
                      "over measured named-protein rows. Nothing in "
                      "either is absent from Affinity; both are views, not "
                      "different datasets."),
    ("", "`source` disambiguates the two Qin et al. 2024 papers, which share "
         "a study name and year but are different work: 'Qin 2024 (plasma)' "
         "is the albumin and plasma binding paper (doi 10.1021/acs.est."
         "3c10824) and 'Qin 2024 (bioassay)' the baseline-toxicity paper "
         "(10.1021/acs.est.3c09950), which measured cells and assay medium."),
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
    ("measure", "Which quantity: Kd, Ka, log D (protein-water), "
                "log K_membrane (phospholipid membrane-water), f_unbound, "
                "docking dG. A '(computed)' or '(computational)' suffix "
                "means the source reports that value as predicted rather "
                "than measured -- QSAR, COSMO-RS, a chain-length regression "
                "or docking. Four rows besides the ten docking ones are "
                "marked; everything unmarked is an experimental value."),
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
    ("protein_conc_uM", "Where the paper states a molarity, that. Where it "
                        "states a MASS concentration of albumin, converted "
                        "at 66,463 g/mol -- unit arithmetic with one "
                        "published constant, applied to albumin only. "
                        "Gamma-globulin (no single molar mass), serum, "
                        "plasma and immobilised stationary phases stay "
                        "blank, and the verbatim string sits beside it."),
    ("pfas_conc_as_stated", "The ligand concentration exactly as the paper "
                            "states it, with its own units. Present whenever "
                            "the paper says anything at all, so a blank in "
                            "pfas_conc_uM can be told apart from a study that "
                            "reported nothing."),
    ("protein_conc_as_stated", "The protein concentration as stated, same "
                               "rule. 93% of rows carry one, against the 27% "
                               "that are molar and so convertible: most of "
                               "the rest are dilution factors, mol/mol "
                               "ratios, an immobilised stationary phase or a "
                               "percentage of serum."),
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


def write(out, major, cond, cmpd, full):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    head_fill = PatternFill("solid", fgColor="2A78D6")
    head_font = Font(bold=True, color="FFFFFF")

    wb = Workbook()
    wb.remove(wb.active)

    for title, fields, data in (("Affinity", FIELDS, out),
                                ("Major studies", FIELDS, major),
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
    print("Cross-check -- every PFOA measurement on serum albumin.")
    print("An earlier version of this block claimed to compare four routes")
    print("to a Kd. It could not: two of the four report a partition")
    print("coefficient, and this table deliberately refuses to convert one")
    print("to a Kd without a stoichiometry the paper states. Those two rows")
    print("never printed and the two ratio lines beneath them never ran.")
    print("What follows instead is every PFOA-on-albumin row, in whichever")
    print("currency its source reports, so the split is visible rather than")
    print("silently filtered:")
    print()
    pfoa = [r for r in out if r["pfas"] == "PFOA"
            and "albumin" in r["protein"].lower()
            and "computational" not in r["measure"]]
    print(f"  {'source':20s}{'species':9s}{'measure':11s}{'value':>10s}  "
          f"{'kd_uM':>9s}{'logK':>8s}")
    for r in sorted(pfoa, key=lambda x: (x["measure"], x["source"])):
        print(f"  {r['source']:20s}{r['species']:9s}{r['measure']:11s}"
              f"{r['value'][:10]:>10s}  {r['kd_uM'] or '-':>9s}"
              f"{r['logK'] or '-':>8s}")
    print()
    kd_rows = [r for r in pfoa if r["kd_uM"]]
    lk_rows = [r for r in pfoa if r["logK"]]
    both = [r for r in pfoa if r["kd_uM"] and r["logK"]]
    print(f"  {len(kd_rows)} rows reach a Kd and {len(lk_rows)} reach a logK, "
          f"of {len(pfoa)} in all.")
    print(f"  Exactly {len(both)} reach BOTH, and only by the one exact route "
          f"there is:")
    print("  the paper reported an affinity AND a stoichiometry, so the")
    print("  conversion needs nothing assumed --")
    print("    " + ", ".join(f"{r['source']} ({r['species']})"
                             for r in both) + ".")
    print("  Every other row sits in one currency or the other. That is the")
    print("  constraint this project keeps running into: almost no study")
    print("  measures both on the same material, so the two halves of the")
    print("  table cannot be compared directly. The bridge in")
    print("  qin2024_tableS6_provenance.py goes around it, through the")
    print("  saturation slope, rather than through a conversion.")
    print()
    # Wu 2024's "binding constant KA" of 103 M^-1 is three to four orders
    # below every other affinity for PFOA on albumin and comes from a
    # docking-led study; it is reported, not excluded, but it would set the
    # range on its own, so the spread is given both ways.
    if lk_rows:
        vals = sorted((float(r["logK"]), r["source"], r["species"])
                      for r in lk_rows)
        for v, src, sp in vals:
            print(f"    {v:5.2f}   {src} ({sp})")
        main = [v for v in vals if v[0] > 1.0]
        print(f"  Spread over all {len(vals)}: {vals[0][0]:.2f} to "
              f"{vals[-1][0]:.2f} log units "
              f"({10 ** (vals[-1][0] - vals[0][0]):.0f}x), but that is set by")
        print(f"  {vals[0][1]} alone. Over the other {len(main)}, which agree "
              f"on the order of magnitude:")
        print(f"  {main[0][0]:.2f} to {main[-1][0]:.2f}, a spread of "
              f"{10 ** (main[-1][0] - main[0][0]):.1f}x across "
              f"{len({m[1] for m in main})} laboratories.")
    print()
    if kd_rows:
        vals = sorted(float(r["kd_uM"]) for r in kd_rows)
        print(f"  The {len(kd_rows)} dissociation constants span {vals[0]:.3g} "
              f"to {vals[-1]:.4g} uM ({vals[-1] / vals[0]:.0f}x),")
        print("  which is dominated by method rather than by chemistry --")
        print("  DSF at the melting transition against equilibrium dialysis")
        print("  at 37 C. Those are not replicates and must not be pooled.")
    print()
    majors = major_sources(out)
    nmaj = sum(1 for r in out if r["source"] in majors)
    print(f"Major-study subset -- sources with at least {MAJOR_MIN_ROWS} "
          f"standardised measurements:")
    for srcname in majors:
        rs = [r for r in out if r["source"] == srcname]
        meas = ", ".join(f"{m} {n}" for m, n in
                         collections.Counter(r["measure"] for r in rs).most_common())
        print(f"  {len(rs):5d}  {srcname:17s} "
              f"{len({r['pfas'] for r in rs}):3d} compounds  {meas}")
    print(f"  {nmaj:5d}  total, {100 * nmaj / len(out):.0f}% of the "
          f"{len(out)} rows, from {len(majors)} of "
          f"{len(set(r['source'] for r in out))} studies")
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
        majors = major_sources(out)
        major = [r for r in out if r["source"] in majors]
        write(out, major, conditions(out), compounds(out), full)
        print(f"wrote {os.path.abspath(OUT)}  ({len(out)} rows)")
        print()
    check(rows, out, skipped, qualitative, dropped)


if __name__ == "__main__":
    main()
