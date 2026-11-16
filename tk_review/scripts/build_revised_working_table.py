#!/usr/bin/env python3
"""Rebuild the reviewer's working table from the database.

The working copy had drifted: 228 rows against the workbook's 320,
one duplicated measurement, a study split across two n policies, and
nine sources absent.  Rather than patch it cell by cell, this rebuilds
it from db/albumin_logK_with_conc.xlsx in the working copy's own
columns and column order, and keeps the two reductions the reviewer
made ON PURPOSE:

  GAO 2025 is a salt-effect study with 42 rows, four salt levels per
  compound.  The working copy kept ONE row per compound and that is
  the right call for a cross-study table -- the salt series belongs in
  the salt analysis, not stamped four times into a compound list.  The
  no-added-salt row is kept, since that is the condition the other
  studies were run at.

  MORO 2022 has 12 rows, site class x temperature.  The working copy
  kept the four representative rows and that is also right: the
  weak/second site classes carry nothing at physiological loading.

What it does NOT preserve is the drift:

  PENG 2024.  The working copy had 12 rows -- six using the paper's
  fitted quenching exponent and six using n = 1 on a different
  temperature's Kd.  Peng reports a fitted n for EVERY row, so there
  is no reason to assume 1 on half of them, and mixing the policies
  inside one study makes the two halves non-comparable.  All 18 rows
  are written with the paper's own n.

  GAO 2019 PFDoDA appeared twice at Kd 103 / log K 2.165: the primary
  from Gao's SI Table S9 and a secondary transcription of the same
  measurement out of Jia 2022's comparison table, invisible because
  the two spell the compound PFDoDA and PFDoA.  The database already
  dropped the secondary, so the rebuild carries one.

Run:  python3 tk_review/scripts/build_revised_working_table.py [-o PATH]
"""
import argparse
import collections
import csv
import os
import sys

import subprocess

import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
OUT = os.path.join(HERE, "..", "albumin_logK_Kd_revised.xlsx")

COLS = ["PFAS", "Species", "Kd (uM)", "n used",
        "log K (L water / kg protein)", "Route", "Source",
        "PFAS conc (uM)", "Protein conc (uM)", "ligand:protein ratio",
        "f_bound (%)"]
FROM = {"ligand:protein ratio": "ligand:protein ratio (max)"}

# Reductions the reviewer made deliberately, preserved.
ONE_PER_COMPOUND = {
    "Gao 2025": "salt-effect study: one row per compound, at 150 mM added "
                "NaCl (Gao's physiological level) or the highest level "
                "measured below it, so the salt series does not enter a "
                "cross-study compound list several times",
}


# SPECIES THIS SHEET CARRIES, set 2026-11-14 by decision.
#
# The reviewer's working copy had already removed Fan 2025 (39 rows,
# rat serum albumin and rat L-FABP, rat only) and Han 2003's rat row,
# keeping the human one.  Stating it as a species filter makes that a
# rule rather than two deletions, so a re-ingestion cannot quietly
# undo it, and it applies to every source at once -- Starnes 2024's
# porcine and rat columns are held back by the same line.
#
# The rows are NOT lost: db/protein_binding.csv keeps every species,
# and the cross-species comparisons the review's first aim rests on
# are computed there, not here (species_albumin_comparison.py and
# albumin_cannot_explain_tk.py both read the database).
#
# To widen the sheet, add to this set; nothing else changes.
SHEET_SPECIES = {"human", "bovine"}
SPECIES_DROP_LABEL = "(species not carried)"


# SOURCES DROPPED FROM THIS SHEET FOR SCOPE, 2026-11-13, by decision.
#
# This is NOT the named-holdout category.  A named holdout (below) is a
# source whose numbers are not what they appear to be, and it is barred
# from figures too, by the leak scan.  A scope drop is a source whose
# numbers are reported as given but which does not earn a row in a
# cross-study working table.  Gao 2025's salt series, Moro 2022's weak
# site classes and Peng 2024's non-310 K rows are the same category,
# handled below; this dict is for dropping a source outright rather
# than reducing it.
#
# A scope drop does NOT touch the figures.  Beesoon's four rows stay in
# figures/albumin_chainlength_data.csv, where an isomer-resolved point
# is informative, and in db/protein_binding.csv.
DROP_SCOPE = {
    "Wu 2009":
        "doi 10.1186/1472-6807-9-31, PFOA on HSA at pH 7.40 and 37 C. "
        "Its Temkin equilibrium-dialysis constant is already refused "
        "by the converter on its own terms (see the Temkin check in "
        "build_binding_affinity_table.py and "
        "scripts/audit_wu2009_temkin.py), leaving only the ITC. That "
        "ITC row is a different KIND of measurement from the rest of "
        "the sheet: a two-step sequential fit whose number is a sum "
        "over two classes of site, n = 1 buried plus n = 12 on the "
        "surface, so it carries no single dissociation constant and "
        "its Kd cell is necessarily blank. Every other converted row "
        "here is one Kd with one stoichiometry. The merge that "
        "produced it (scripts/merge_wu2009_itc.py) and the per-class "
        "summing it needed stay in place, so re-admitting is one line "
        "here and nothing else; the rows remain in "
        "db/protein_binding.csv and in the figures.",
    "Klevens 1954":
        "doi 10.1039/df9541800277, the earliest PFAS-albumin constant "
        "in the literature. Dropped on reading the primary. It is a "
        "FLOCCULATION study -- its own abstract is about "
        "protein-protein aggregation, 'a marked loss of potential "
        "binding sites with increase in protein concentration' -- run "
        "in acetate at pH 4.3-5.8, straddling BSA's isoelectric point "
        "near 4.7, not at blood pH. Its stoichiometry is not a site "
        "count: n falls from about 150 at 0.5% BSA to about 45 at 4% "
        "because the protein associates, and the 85 in Fig. 3 belongs "
        "to the 1% solutions (147 uM) while the row stores 588 uM "
        "(4%), whose own n is about 45 -- so the protein column and "
        "the constant come from different conditions. Its ligand axis "
        "is a BOUND ratio, 'Moles C8 Acid Bound per Mole BSA', not a "
        "dose, so the row can carry no PFAS concentration. For scale, "
        "Han 2003 counts 6 to 9 sites on the same protein at pH 7.4 "
        "and Beesoon 2015 counts 2. Kept in db/protein_binding.csv "
        "with the primary now read and its pH, temperature and "
        "provenance corrected (fix_klevens1954_from_primary.py), "
        "which is where its historical value belongs.",
    "Peng 2025":
        "doi 10.1016/j.jes.2024.02.011. The row MISDESCRIBES its own "
        "paper. It was ingested from a secondary (Starnes & Belcher "
        "2026 Table 1) with source_table 'primary not retrieved', and "
        "records 7 compounds (PFHxA, PFOA, PFNA, PFDA, PFHxS, PFOS, "
        "HFPO-DA) on 'human; bovine'. The paper -- Peng S-Y, Yang "
        "Y-D, Tian R, Lu N, J Environ Sci (China) 2025;149:88-98, "
        "PMID 39181681 -- is 'Critical new insights into the "
        "interactions of hexafluoropropylene oxide-dimer acid (GenX "
        "or HFPO-DA) with albumin', and its abstract explores 'the "
        "binding characteristics of bovine serum albumin (BSA) to "
        "GenX' against PFOA. Two compounds, BSA only; its MeSH terms "
        "carry Serum Albumin, Bovine and no human albumin term. So "
        "both the compound list and the species are wrong, and the "
        "Ka range attributed to the row is unverified -- it probably "
        "belongs to a different row of the secondary's table. Not a "
        "conflation with Peng 2024, which is Peng M et al., Toxics, "
        "doi 10.3390/toxics12010043, a different group. Three "
        "further reasons even setting that aside: the stored Kd "
        "15.64 uM is the GEOMETRIC MIDPOINT of a 21-fold Ka range "
        "(Kd 3.42-71.4 uM, log K 2.324-3.643), so it is a midpoint "
        "and not a measurement; the 7-item PFAS string matches no "
        "compound bucket, so the row entered no per-compound "
        "statistic; and it carries no PFAS concentration, protein "
        "concentration, ratio, temperature or buffer at all. Method "
        "is fluorescence quenching, the family whose constants "
        "failed three separate audits in this compilation, and its "
        "'the paper reports none' claim about n is probably false "
        "for the same reason it was false for Klevens. The row stays "
        "in db/protein_binding.csv as the record of what the "
        "secondary said.",
    "Beesoon 2015":
        "PFOS and PFOA isomers by ultrafiltration, the only "
        "isomer-resolved source here. Dropped from the sheet for four "
        "reasons that compound: (1) its labels, 'linear PFOA (n-PFOA)' "
        "and 'linear PFOS (n-PFOS)', are kept out of the PFOA and PFOS "
        "pools by db/pfas_name_map.csv, so the rows enter no "
        "per-compound statistic and sit in buckets of one; (2) the "
        "Scatchard axis units and the reported Kd disagree by 1000x -- "
        "three independent checks (mass balance against the stated "
        "2-12 uM spike, the donated stock concentrations, and Kd versus "
        "the free concentrations actually sampled) all say the printed "
        "constants are 1000x too large, see "
        "audit_beesoon2015_scatchard.py; (3) the buffer is HPLC-grade "
        "water at pH 6.4 with 2-10% residual methanol, not the pH 7.4 "
        "phosphate the rest of the sheet uses; (4) linear PFOA's n does "
        "not reproduce -- slope 0.6076 gives 1.65, the box says 2. "
        "The isomer-mixture finding that made linear the right "
        "representative is preserved in "
        "scripts/beesoon2015_isomer_mixture.py and does not depend on "
        "these rows being in the sheet.",
}


# NAMED HOLDOUTS dropped from this sheet, 2026-11-13, by decision.
#
# scripts/excluded_sources.py already carried Wu 2024 as a named
# holdout, which kept it out of every FIGURE: figures/
# albumin_chainlength_data.csv and figures/lfabp_by_compound_data.csv
# have none of its rows, and the leak scan enforces that.  What the
# rule deliberately ALLOWED was its presence in the reading and
# working copies, "because these are the record".  This sheet is the
# one a reviewer copies into Excel and pools by hand, so the record
# argument does not protect it: the ten rows were the sheet's entire
# negative-log-K population and set its lower bound single-handed.
#
# They stay in db/protein_binding.csv with the full audit, which is
# where the record belongs.  The reason is imported rather than
# retyped so that this sheet and the leak scan cannot drift apart.
#
# To re-admit a source, take it out of this set; nothing else changes.
DROP_HOLDOUTS = {"Wu 2024"}

# Non-fluorinated compounds. Several DSF and displacement studies run
# straight-chain fatty acids as positive controls or competitors, because
# albumin's fatty-acid sites are what PFAS are thought to occupy. They are
# real measurements on the same assay, but they are not PFAS, and a PFAS
# coefficient sheet that carries them will have them swept into any
# chain-length fit or per-compound median. Dropped at the COMPOUND level, so
# the rest of each source stays.
#
# Reviewer asked for the three Starnes 2025 controls to go. The same rule
# reaches five pre-existing Chen 2009 rows (octanoic acid x3, butyric acid x2)
# that had been in the sheet all along -- flagged rather than silently kept,
# since leaving them would make the rule inconsistent.
NON_PFAS = {
    "octanoic acid", "decanoic acid", "hexadecanoic acid",
    "butyric acid (c4)",
}
NON_PFAS_DROP_LABEL = "(not a PFAS)"

# A cell naming SEVERAL compounds at once is a summary, not a measurement of
# any one of them, and its value is a midpoint describing none. This project
# has now hit that pattern in Peng 2025, Qin 2010, Chen 2020, Chi 2018,
# Fedorenko 2021, MacManus-Spencer 2010, Crisalli 2023, Jackson 2021 -- and in
# Jia 2022, whose row "HFPO-DA (GenX), HFPO-TA, HFPO-TeA, PFEESA, 8:2
# Cl-PFESA" reached the sheet at a Kd of 185.5 uM while the same five
# compounds also appear individually from the same paper. A comma outside
# brackets is the test; a chemical name's own commas sit inside them.
MULTI_DROP_LABEL = "(several compounds in one cell)"


def names_several(name, known):
    """True if the cell lists SEVERAL compounds.

    A bare comma test does not work: an IUPAC name's commas are locant
    separators, so "2,2,3,3,4,4-hexafluorobutyric acid" and
    "Perfluoro-3,6-dioxaheptanoic acid" are single compounds. Testing commas
    alone flagged ten Starnes 2025 rows that are perfectly good single
    measurements.

    So the test is data-driven instead of a regex guess: split on commas
    outside brackets and ask how many parts are themselves names that appear
    in this dataset as standalone compounds. Two or more means a list.
    "HFPO-DA (GenX), HFPO-TA, HFPO-TeA, PFEESA, 8:2 Cl-PFESA" yields four such
    names; "Perfluoro-3,6-dioxaheptanoic acid" yields none, because
    "Perfluoro-3" and "6-dioxaheptanoic acid" are not compounds anyone reports.
    """
    import re
    txt = str(name or "")
    if "," not in re.sub(r"\([^)]*\)", "", txt):
        return False
    parts, depth, cur = [], 0, ""
    for ch in txt:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            parts.append(cur); cur = ""
        else:
            cur += ch
    parts.append(cur)
    hits = sum(1 for q in parts
               if re.sub(r"\s*\([^)]*\)", "", q).strip().lower() in known)
    return hits >= 2

# Rows whose constant was measured in a co-solvent high enough to act on the
# quantity being measured. Keyed on (source, compound), NOT on compound alone:
# 6:2 FTCA and N-EtFOSA are also measured by Chen 2025 by equilibrium dialysis
# at log K 3.62 and 3.06, and those rows are sound. A compound-wide rule would
# have deleted them.
#
# Starnes 2025 Table S3 gives the solvent per chemical. Three of its admitted
# rows were prepared in 50% HBS / 50% DMSO because they would not dissolve in
# buffer. DSF reads out albumin's melting temperature, and DMSO both
# destabilises albumin and competes for its hydrophobic sites, so at 50% the
# co-solvent acts on the readout itself rather than sitting beside it. Note
# how far they fall below the same compounds measured in buffer by Chen 2025:
# 1.42 against 3.62, and 1.43 against 3.06.
#
# Same class of objection that made Wu 2024 a named holdout over 86% ethanol.
# Reviewer's decision, taken after the three were flagged.
COSOLVENT_HOLDOUT = {
    ("Starnes 2025", "9h-perfluorononanoic acid"),
    ("Starnes 2025", "6:2 ftca"),
    ("Starnes 2025", "n-etfosa"),
}
COSOLVENT_DROP_LABEL = "(50% DMSO co-solvent)"

# Minimum R2 for a fitted single-site constant, by source. A constant is only
# as good as the curve it was read off, and Starnes 2025 Table S4 publishes the
# R2 for every fit, so the weak ones can be identified rather than guessed at.
# Below 0.80 the concentration-response is not well described by the model that
# produced the Kd, which makes the Kd a number without a curve behind it.
#
# Applied ONLY to sources that publish a per-row R2. A source that reports no
# R2 is not silently admitted under this rule -- it is simply not covered by
# it, and r2_floor_drops() reports how many rows each source had to check so a
# source whose R2 went missing cannot pass unnoticed.
R2_FLOOR = {"Starnes 2025": 0.80}
R2_DROP_LABEL = "(R2 below the floor)"

# Sources admitted on condition that they add DEPTH, not BREADTH. A row from
# one of these enters only if the sheet already holds that compound from some
# other source, so these studies strengthen comparisons that already exist
# instead of adding single-source compounds nothing can be checked against.
#
# Reviewer's call: "keep only the PFAS for which we already have data and not
# add any unique/new PFAS at this stage". Reversible -- empty the set and the
# dropped rows return.
#
# An isomer suffix is stripped before the test, because "PFOS (linear)" and
# "PFOS (branched)" are PFOS, which the sheet holds from 14 sources; they are
# Li's isomer resolution of a compound already present, not new compounds.
# "PFPrS (linear)" is NOT spared by that: the sheet holds no PFPrS at all, so
# it is genuinely new and goes.
DEPTH_ONLY = {"Li 2021", "Starnes 2025", "Jia 2022"}
DEPTH_DROP_LABEL = "(compound new to the sheet)"

# Isomer-resolved rows. Li 2021 is the only source that splits a compound into
# linear and branched; everyone else reports a single value for technical
# material, which is 98.5-99.99% linear. So:
#
#   BRANCHED is dropped. It is a different measurement that does not pool with
#   a linear-dominated value, which is where the Beesoon 2015 discussion
#   landed. Reviewer deleted these by hand; this makes the builder agree.
#
#   LINEAR keeps its number and loses its label, because "PFOS (linear)" and
#   every other source's "PFOS" are the same quantity and must share a row key
#   or they never compare. Reviewer renamed these by hand.
#
# Stripping the label runs BEFORE the depth gate, so a renamed row is tested
# against the sheet under the name it will actually carry.
DROP_BRANCHED = True
BRANCHED_DROP_LABEL = "(branched isomer)"


def strip_isomer_label(name):
    """'PFOS (linear)' -> 'PFOS'. Leaves bracketed synonyms and chemical names
    alone -- only the isomer words go."""
    import re
    return re.sub(r"\s*\((linear|n-)\)\s*$", "", str(name or "").strip(),
                  flags=re.I).strip()


def is_branched(name):
    import re
    return bool(re.search(r"\((branched|iso)\)\s*$", str(name or "").strip(),
                          flags=re.I))


def isomer_base(name):
    """'PFOS (linear)' -> 'pfos'. Only strips a bracketed isomer word, so a
    bracketed synonym or a chemical name with brackets is left alone."""
    import re
    n = str(name or "").strip()
    n = re.sub(r"\s*\((linear|branched|n-|iso)\)\s*$", "", n, flags=re.I)
    return n.strip().lower()


GAO_SALT_PREFERRED = 150          # mM added NaCl; Gao marks this physiological


def gao_salt_pick():
    """{compound: log K at the preferred salt level}, read from the database.

    Gao 2025 reports a Hill association constant per compound per salt level
    (20 mM sodium phosphate pH 7 plus 0/10/30/80/150/300/500 mM NaCl), so
    log K = log10(Ka) + log10(n=1) + BRIDGE - 6. Only four compounds were run
    at 150 mM; the rest fall back to the highest level measured, which is
    80 mM. Returns the target value, and the caller picks the row matching it,
    so this never invents a number the pipeline did not already derive.
    """
    import math
    import re as _re
    db = os.path.join(HERE, "..", "db", "protein_binding.csv")
    best = {}
    with open(db, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if "Gao" not in (r.get("study") or "") or (r.get("year") or "") != "2025":
                continue
            if not r.get("ka"):
                continue
            m = _re.search(r"plus (\d+) mM added NaCl", r.get("buffer") or "")
            if not m:
                continue
            salt = int(m.group(1))
            if salt > GAO_SALT_PREFERRED:
                continue              # 300 and 500 mM are past physiological
            c = (r.get("chemical") or "").strip()
            lk = math.log10(float(r["ka"])) - 6.0 + math.log10(1e9 / 66463.0)
            if c not in best or salt > best[c][0]:
                best[c] = (salt, lk)
    # Key on BOTH the database spelling and the canonical one. The database
    # says "HFPO-DA (GenX)" where the sheet says "HFPO-DA", and keying only on
    # the raw name silently left that compound on the no-added-salt row while
    # the other nine moved -- a miss that looks like data, not a bug.
    out = {}
    for c, (salt, lk) in best.items():
        out[c] = lk
        bare = _re.sub(r"\s*\([^)]*\)\s*$", "", c).strip()
        if bare and bare != c:
            out.setdefault(bare, lk)
    return out


def r2_floor_drops():
    """(source, compound-lowercased) pairs failing their source's R2 floor,
    read from the database's source_table field, plus a per-source tally."""
    import re
    db = os.path.join(HERE, "..", "db", "protein_binding.csv")
    drop, seen, no_r2 = set(), collections.Counter(), collections.Counter()
    with open(db, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            src = f"{(r.get('study') or '').replace('et al.', '').strip()} {r.get('year') or ''}".strip()
            src = src.rstrip(',')
            floor = R2_FLOOR.get(src)
            if floor is None or not (r.get("ka") or r.get("kd")):
                continue
            seen[src] += 1
            m = re.search(r"R2 ([\d.]+)", r.get("source_table") or "")
            if not m:
                no_r2[src] += 1
                continue
            if float(m.group(1)) < floor:
                drop.add((src, (r.get("chemical") or "").strip().lower()))
    return drop, seen, no_r2


def holdout_reason(src):
    """The one recorded reason, read from scripts/excluded_sources.py."""
    try:
        sys.path.insert(0, HERE)
        from excluded_sources import SOURCES
        return " ".join(SOURCES[src].split())
    except Exception:
        return "named holdout; see scripts/excluded_sources.py"


# Studies whose rows differ only by assay temperature.  The review is not
# about temperature dependence, so one row per compound reaches the sheet
# and the rest stay in db/protein_binding.csv.  The temperature named here
# is the one kept.
BY_TEMPERATURE = {
    "10.3390/ijerph17041319": ("310 K", "Chen 2020",
                               "PFOA and PFDA each measured at 296, "
                               "303 and 310 K; 310 K is body "
                               "temperature"),
    "10.1021/jf100412q": ("310 K", "Qin 2010",
                          "PFOA and PFDA each measured at 300, 310 and "
                          "320 K; 310 K is body temperature"),
    "10.3390/toxics12010043": ("37 (310 K)", "Peng 2024",
                               "six PFAS measured at 298, 304 and 310 K; "
                               "310 K is body temperature"),
}


def temperature_picks():
    """(compound, n) of the rows to keep, for each study reduced by temperature.

    Keyed on the STOICHIOMETRY, not on Kd.  Kd is derived, and a derivation
    can change -- the double-log nth-root correction in
    build_binding_affinity_table moved every Peng Kd by orders of
    magnitude and silently emptied this list when it was keyed on Kd.  The
    n cell is transcribed verbatim from the paper, reaches the sheet
    unaltered, and is unique within a compound here: Peng's three
    temperatures give n = 1.51 / 1.18 / 1.10 for PFNA and so on for the
    rest.  Read from the database so a re-ingestion cannot put the sheet
    and this list out of step.
    """
    out = {}
    db = os.path.join(HERE, "..", "db", "protein_binding.csv")
    with open(db, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            spec = BY_TEMPERATURE.get(r["doi"])
            if not spec or r["temperature_c"] != spec[0]:
                continue
            out.setdefault(spec[1], set()).add(
                (r["chemical"],) + row_key(r))
    return out


def row_key(r):
    """A key that separates one temperature's row from another's.

    Prefer the STOICHIOMETRY: it is transcribed verbatim from the paper
    and reaches the sheet unaltered.  Where a paper reports none --
    Qin 2010 has a Stern-Volmer and a modified Stern-Volmer plot and no
    double-logarithmic one, so there is no fitted n anywhere in it --
    every row would carry the reader's assumed n = 1 and n could not
    tell the temperatures apart.  Fall back to Kd there.

    The docstring above warns against keying on Kd, and that warning
    was earned: the nth-root correction moved every Peng Kd and emptied
    this list.  It applies to a Kd that a DERIVATION can move.  Qin's
    is 1/Ka on a constant transcribed from Table 2, and nothing in the
    pipeline transforms it, so the objection does not reach it.  The
    fallback is used only when n is absent, so Peng keeps the n key.
    """
    try:
        return ("n", "%g" % float(r["n_binding_sites"]))
    except (TypeError, ValueError):
        pass
    try:
        return ("kd", "%g" % float(r["kd"]))
    except (TypeError, ValueError):
        return ("none", "")


def moro_picks():
    """The (compound, Kd) pairs moro2022_representative.py selects.

    Moro's reduction is a real decision -- tightest site class, 37 C
    where the cell is usable -- and that decision already lives in one
    script. Re-deriving it here would give two places to keep in step,
    so call that script and read its answer.
    """
    r = subprocess.run([sys.executable,
                        os.path.join(HERE, "moro2022_representative.py"),
                        "-o", "-"], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"moro2022_representative.py failed:\n{r.stderr}")
    lines = r.stdout.strip().split("\n")
    hdr = lines[0].split("\t")
    i_kd = hdr.index("Kd (uM)")
    # Match on Kd alone: it is unique inside this study (0.4, 5.4, 5.9,
    # 9.9), and the compound NAME is not a safe key here -- the
    # representative table says "HFPO-DA (GenX)" where the workbook says
    # "HFPO-DA", so matching on the name silently dropped that compound.
    return {c.split("\t")[i_kd] for c in lines[1:]}


def val(d, name):
    return d.get(FROM.get(name, name))


def setval(d, name, v):
    d[FROM.get(name, name)] = v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=OUT)
    args = ap.parse_args()

    ws = openpyxl.load_workbook(SRC, data_only=True)["Coefficients"]
    h = [c.value for c in ws[1]]
    rows = [dict(zip(h, r)) for r in ws.iter_rows(min_row=2, values_only=True)
            if r[0]]
    print(f"  source: {len(rows)} rows")

    picks = moro_picks()
    print(f"  Moro 2022 representative picks: {len(picks)}")
    temps = temperature_picks()
    for src, keep in sorted(temps.items()):
        spec = [v for v in BY_TEMPERATURE.values() if v[1] == src][0]
        print(f"  {src}: keeping the {spec[0]} rows only "
              f"({len(keep)} compounds) -- {spec[2]}")
    kept, dropped, renamed = [], collections.Counter(), collections.Counter()
    seen = collections.defaultdict(list)
    r2_failed, r2_seen, r2_missing = r2_floor_drops()
    # Compound names that appear in this dataset WITHOUT a comma, i.e. as
    # standalone measurements. names_several() tests a comma-bearing cell
    # against these, so the judgement comes from the data rather than a regex.
    single_names = {str(val(d, "PFAS")).strip().lower() for d in rows
                    if val(d, "PFAS") and "," not in str(val(d, "PFAS"))}
    # Compounds the sheet holds from a source OUTSIDE the depth-only set.
    #
    # This must be built from rows that SURVIVE the other exclusions, not from
    # raw input. Built naively it counted PFPrS as established because Wu 2024
    # carries it -- and Wu 2024 is a named holdout that never reaches the
    # sheet. That would have admitted Li's PFPrS (linear) as "depth" on a
    # compound the sheet does not actually have, which is the opposite of what
    # the gate is for. A compound is established only if a row for it is
    # actually kept.
    def excluded_otherwise(d):
        src_ = str(val(d, "Source"))
        sp_ = str(val(d, "Species") or "").strip()
        c_ = str(val(d, "PFAS") or "").strip().lower()
        return (bool(sp_) and sp_ not in SHEET_SPECIES
                or c_ in NON_PFAS
                or (src_, c_) in COSOLVENT_HOLDOUT
                or (src_, c_) in r2_failed
                or src_ in DROP_HOLDOUTS
                or src_ in DROP_SCOPE)

    established = {isomer_base(val(d, "PFAS")) for d in rows
                   if str(val(d, "Source")) not in DEPTH_ONLY
                   and not excluded_otherwise(d)}
    if DEPTH_ONLY:
        print(f"  depth-only sources {sorted(DEPTH_ONLY)}: admitted against "
              f"{len(established)} compounds already in the sheet")
    for src_name, floor in sorted(R2_FLOOR.items()):
        n = sum(1 for s_, _ in r2_failed if s_ == src_name)
        print(f"  R2 floor {floor} on {src_name}: checked "
              f"{r2_seen.get(src_name, 0)} db row(s), {n} below the floor"
              + (f", {r2_missing[src_name]} with NO R2 recorded"
                 if r2_missing.get(src_name) else ""))
    for d in rows:
        src = str(val(d, "Source"))
        sp = str(val(d, "Species") or "").strip()
        if sp and sp not in SHEET_SPECIES:
            dropped[SPECIES_DROP_LABEL] += 1
            continue
        if DROP_BRANCHED and is_branched(val(d, "PFAS")):
            dropped[BRANCHED_DROP_LABEL] += 1
            continue
        # Rename before anything else reads the compound: the depth gate, the
        # duplicate key and the output column must all see the final name.
        stripped = strip_isomer_label(val(d, "PFAS"))
        if stripped != str(val(d, "PFAS") or "").strip():
            setval(d, "PFAS", stripped)
            renamed[str(val(d, "Source"))] += 1
        if names_several(val(d, "PFAS"), single_names):
            dropped[MULTI_DROP_LABEL] += 1
            continue
        cmpd = str(val(d, "PFAS") or "").strip().lower()
        if cmpd in NON_PFAS:
            dropped[NON_PFAS_DROP_LABEL] += 1
            continue
        if (src, cmpd) in COSOLVENT_HOLDOUT:
            dropped[COSOLVENT_DROP_LABEL] += 1
            continue
        if (src, cmpd) in r2_failed:
            dropped[R2_DROP_LABEL] += 1
            continue
        if src in DEPTH_ONLY and isomer_base(cmpd) not in established:
            dropped[DEPTH_DROP_LABEL] += 1
            continue
        if src in DROP_HOLDOUTS:
            dropped[src] += 1
            continue
        if src in DROP_SCOPE:
            dropped[src] += 1
            continue
        if src == "Moro 2022":
            key = (f"{float(val(d, 'Kd (uM)')):g}"
                   if val(d, "Kd (uM)") is not None else "")
            if key in picks:
                kept.append(d)
            else:
                dropped["Moro 2022"] += 1
            continue
        if src in temps:
            keys = set()
            for kind, col in (("n", "n used"), ("kd", "Kd (uM)")):
                try:
                    keys.add((str(val(d, "PFAS")), kind,
                              "%g" % float(val(d, col))))
                except (TypeError, ValueError):
                    pass
            if keys & temps[src]:
                kept.append(d)
            else:
                dropped[src] += 1
            continue
        if src in ONE_PER_COMPOUND:
            k = (str(val(d, "PFAS")), str(val(d, "Species")), src)
            seen[k].append(d)
            continue
        kept.append(d)
    # For the reduced sources, keep the row with the LOWEST ratio.
    #
    # GAO 2025 IS THE EXCEPTION, and the default was wrong for it. Lowest
    # ratio selected the NO-ADDED-SALT condition: 20 mM phosphate and nothing
    # else. Gao's own finding is that salt SUPPRESSES PFAS-albumin binding, so
    # the saltless row is the one that most overstates affinity, and it is the
    # least like plasma. Reviewer's call: take 150 mM added NaCl, which Gao
    # marks as physiological, and fall back to the highest level measured for
    # compounds Gao did not run at 150. That shifts ten compounds down by 0.18
    # to 0.34 log. See gao_salt_pick().
    gao_target = gao_salt_pick()
    for k, g in seen.items():
        def rat(d):
            try:
                return float(str(val(d, "ligand:protein ratio")))
            except (TypeError, ValueError):
                return float("inf")
        if k[2] == "Gao 2025" and str(val(g[0], "PFAS")) in gao_target:
            want = gao_target[str(val(g[0], "PFAS"))]
            g2 = sorted(g, key=lambda d: abs(float(val(d, "log K (L water / kg protein)")) - want))
        else:
            g2 = sorted(g, key=rat)
        kept.append(g2[0])
        dropped[k[2]] += len(g2) - 1

    for s, n in sorted(dropped.items()):
        bytemp = [v for v in BY_TEMPERATURE.values() if v[1] == s]
        if s == MULTI_DROP_LABEL:
            why = ("the cell names several compounds, so its value is a "
                   "midpoint describing none of them. The rows stay in the "
                   "database")
        elif s == BRANCHED_DROP_LABEL:
            why = ("a branched-isomer row. Every other source reports one "
                   "value for technical material, which is overwhelmingly "
                   "linear, so a branched constant does not pool with it. "
                   "The rows stay in the database")
        elif s == DEPTH_DROP_LABEL:
            why = ("the compound is not in the sheet from any other source, "
                   "and " + repr(sorted(DEPTH_ONLY)) + " are admitted for "
                   "depth on existing compounds only. The rows stay in the "
                   "database")
        elif s == R2_DROP_LABEL:
            why = ("the paper's own R2 for the single-site fit is below "
                   + repr(R2_FLOOR) + ", so the Kd has no well-described "
                   "curve behind it. The rows stay in the database")
        elif s == COSOLVENT_DROP_LABEL:
            why = ("measured in 50% DMSO (Starnes 2025 Table S3), which acts "
                   "on the melting temperature DSF reads out. The rows stay "
                   "in the database; the same compounds measured in buffer by "
                   "other studies are unaffected")
        elif s == NON_PFAS_DROP_LABEL:
            why = ("not a PFAS: a non-fluorinated fatty acid run as a "
                   "positive control or competitor. The rows stay in the "
                   "database; " + repr(sorted(NON_PFAS)))
        elif s == SPECIES_DROP_LABEL:
            why = ("species outside SHEET_SPECIES "
                   + repr(sorted(SHEET_SPECIES)) + "; the rows stay in "
                   "db/protein_binding.csv, which is what the "
                   "cross-species scripts read")
        elif s in DROP_SCOPE:
            why = ("DROPPED FOR SCOPE (not a named holdout; still in the "
                   "figures and the database): " + DROP_SCOPE[s])
        elif s in DROP_HOLDOUTS:
            why = ("NAMED HOLDOUT, dropped from this sheet: "
                   + holdout_reason(s))
        elif bytemp:
            why = (f"kept only the {bytemp[0][0]} rows -- {bytemp[0][2]}; "
                   "the other temperatures stay in db/protein_binding.csv")
        else:
            why = ONE_PER_COMPOUND.get(
                s, "representative row per compound from "
                   "moro2022_representative.py: tightest site class, 37 C "
                   "where the cell is usable")
        print(f"  {s}: dropped {n} ({why})")

    out = []
    for d in kept:
        out.append([("" if val(d, c) is None else val(d, c)) for c in COLS])
    out.sort(key=lambda r: (str(r[6]), str(r[0]), str(r[1])))

    wb = openpyxl.Workbook()
    sh = wb.active
    sh.title = "Coefficients"
    sh.append(COLS)
    for r in out:
        sh.append(r)
    hdr_fill = PatternFill("solid", fgColor="DDE5EF")
    for c in sh[1]:
        c.font = Font(bold=True)
        c.fill = hdr_fill
        c.alignment = Alignment(wrap_text=True, vertical="top")
    sh.freeze_panes = "A2"
    for i, w in enumerate([26, 9, 11, 9, 16, 46, 14, 17, 15, 14, 24], 1):
        sh.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w

    nt = wb.create_sheet("Notes")
    for line in [
        ("log K = log10(1e9 / M) + log10(n) - log10(Kd in uM), "
         "M = 66,463 g/mol for albumin"),
        "",
        ("Read the Route column. Rows reach their log K by different "
         "routes, from a direct measurement to a conversion resting on an "
         "assumed stoichiometry. Two values in this column are not "
         "necessarily the same kind of number."),
        "",
        ("ligand:protein ratio is the molar loading the measurement was "
         "made at. Below ~0.01 the assay sees the specific site; above ~1 "
         "it increasingly sees non-specific partitioning. It is NOT "
         "always PFAS conc / protein conc. For Allendorf 2019 and "
         "Klevens 1954 the column is the paper's own reported BOUND "
         "ratio, moles of PFAS bound per mole of protein, not a dosing "
         "quotient -- Klevens plots exactly that quantity, 'Moles C8 "
         "Acid Bound per Mole BSA (MR)', and its 29 is where the "
         "binding curve breaks, not a concentration anyone dosed. That "
         "is why those rows leave PFAS conc blank: multiplying the "
         "ratio by the protein concentration would manufacture a dose "
         "out of a bound stoichiometry. For ITC studies the two "
         "concentrations are in different vessels."),
        "",
        ("Wu 2024 is NOT in this sheet. Its ten HSA fluorescence-"
         "quenching rows are a named holdout in "
         "scripts/excluded_sources.py and were already absent from "
         "every figure; as of 2026-11-13 they are out of this sheet "
         "too. They were its only negative log K values and set its "
         "lower bound alone. The rows, and the full audit, remain in "
         "db/protein_binding.csv."),
        "",
        ("This sheet carries HUMAN and BOVINE albumin only. Rat and "
         "porcine rows exist in db/protein_binding.csv and are used by "
         "the cross-species scripts, but are not pooled here."),
        "",
        ("Beesoon 2015 is NOT in this sheet. Its labels are held out "
         "of the PFOA and PFOS pools by db/pfas_name_map.csv so the "
         "rows pooled with nothing; its Scatchard axis units and "
         "reported Kd disagree by 1000x; and its buffer is pH 6.4 "
         "water with residual methanol, not pH 7.4 phosphate. This is "
         "a SCOPE drop, not a named holdout: the rows remain in "
         "db/protein_binding.csv and in the chain-length figure."),
        "",
        ("Gao 2025 is reduced to one row per compound, the no-added-salt "
         "condition; its full four-level salt series is in "
         "db/albumin_logK_with_conc.xlsx."),
        ("Moro 2022 is reduced to the representative row per compound: "
         "tightest site class, 37 C where usable."),
        "",
        ("Rebuilt by scripts/build_revised_working_table.py from "
         "db/albumin_logK_with_conc.xlsx."),
    ]:
        nt.append([line])
    nt.column_dimensions["A"].width = 110
    for c in nt["A"]:
        c.alignment = Alignment(wrap_text=True, vertical="top")

    wb.save(args.out)
    cnt = collections.Counter(str(r[6]) for r in out)
    print(f"\n  wrote {len(out)} rows / {len(cnt)} sources to "
          f"{os.path.relpath(args.out, os.path.join(HERE, '..'))}\n")
    for s, n in sorted(cnt.items()):
        print(f"    {s:22}{n:>5}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
