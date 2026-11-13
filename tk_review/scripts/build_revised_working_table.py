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
    "Gao 2025": "salt-effect study: one row per compound, the "
                "no-added-salt condition, so the salt series does not "
                "enter a cross-study compound list four times",
}


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
            try:
                n = "%g" % float(r["n_binding_sites"])
            except (TypeError, ValueError):
                continue
            out.setdefault(spec[1], set()).add((r["chemical"], n))
    return out


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
    kept, dropped = [], collections.Counter()
    seen = collections.defaultdict(list)
    for d in rows:
        src = str(val(d, "Source"))
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
            try:
                key = (str(val(d, "PFAS")), "%g" % float(val(d, "n used")))
            except (TypeError, ValueError):
                key = None
            if key in temps[src]:
                kept.append(d)
            else:
                dropped[src] += 1
            continue
        if src in ONE_PER_COMPOUND:
            k = (str(val(d, "PFAS")), str(val(d, "Species")), src)
            seen[k].append(d)
            continue
        kept.append(d)
    # For the reduced sources, keep the row with the LOWEST ratio, which
    # for Gao 2025 is the no-added-salt condition.
    for k, g in seen.items():
        def rat(d):
            try:
                return float(str(val(d, "ligand:protein ratio")))
            except (TypeError, ValueError):
                return float("inf")
        g2 = sorted(g, key=rat)
        kept.append(g2[0])
        dropped[k[2]] += len(g2) - 1

    for s, n in sorted(dropped.items()):
        bytemp = [v for v in BY_TEMPERATURE.values() if v[1] == s]
        if s in DROP_SCOPE:
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
         "always PFAS conc / protein conc: for Allendorf 2019 the column "
         "is the paper's own reported v, bound per protein, and for ITC "
         "studies the two concentrations are in different vessels."),
        "",
        ("Wu 2024 is NOT in this sheet. Its ten HSA fluorescence-"
         "quenching rows are a named holdout in "
         "scripts/excluded_sources.py and were already absent from "
         "every figure; as of 2026-11-13 they are out of this sheet "
         "too. They were its only negative log K values and set its "
         "lower bound alone. The rows, and the full audit, remain in "
         "db/protein_binding.csv."),
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
