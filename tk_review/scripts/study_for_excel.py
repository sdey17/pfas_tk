#!/usr/bin/env python3
"""One study's binding rows as a tab-separated table to paste into Excel.

Every study audit in this project ends the same way: the columns the
reviewer keeps in a spreadsheet, for one study, ready to paste.  Those
tables were being built by hand each time, which is how a dropped tab
once made the Allendorf columns look shifted.  This builds them from
db/protein_binding.csv so the table always matches the database --
including after a ratio or range correction.

COLUMNS, and where each comes from:

  PFAS, Species, Site      the database columns as they stand
  Kd (uM)                  `kd`, only rows that have one
  n used                   `n_binding_sites` when it parses as a bare
                           number, else 1.  A cell that merely
                           DESCRIBES n ("n_Hill 2.08 (mean of 2); a
                           Hill coefficient, NOT a site count") is not
                           a number and does not become one here.
  log K                    BRIDGE + log10(n) - log10(Kd uM), the
                           project bridge, albumin M = 66,463
  Route                    which of those two n cases applied
  PFAS conc (uM)           top of `ligand_conc_range`, unit-converted
  Protein conc (uM)        `protein_conc`, unit-converted
  ligand:protein ratio     `ligand_over_protein_max` as stored
  f_bound (%)              from `unbound_fraction` when present

It computes nothing the database does not already support and writes
no file inside db/ or figures/, so it cannot feed a figure or a
shipped workbook.  It is a reading copy.

Run:  python3 tk_review/scripts/study_for_excel.py "Chen & Guo" 2009
      python3 tk_review/scripts/study_for_excel.py Bischel 2011 -o -
`-o -` prints to stdout instead of writing a file.
"""
import argparse
import csv
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..")

# The project already decides how a study label becomes a figure/source
# label ("Chen & Guo" -> "Chen", and the two-word surnames it knows).
# Use that, so this table's Source column matches every other artifact.
_spec = importlib.util.spec_from_file_location(
    "_baf", os.path.join(HERE, "build_binding_affinity_table.py"))
_baf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_baf)
first_author = _baf.first_author
# WHAT COUNTS AS A PFAS HERE, and why it is not a name regex.
# build_binding_affinity_table's own classifier looked like the answer
# and is not.  Its multi-compound test fires on any comma, so the probe
# "1,8-ANS" reads as an aggregate rather than a non-PFAS, while its
# PFAS_LIKE pattern does not recognise EtFOSA, EtFBSA, 9Cl-PF3ONS,
# NFDHA, M3-Cl-PFECA or Nafion byproduct 2 -- so filtering on it DROPS
# REAL PFAS: 17 names in Starnes 2025 and 53 in Muensterman 2026.
# Hand-written fluorine patterns fail the same way in the other
# direction, and both failures are silent.
#
# So this does not try to recognise a PFAS from its name.  It asks only
# what the database already records about the row:
#
#   n_perfluorinated_carbons == 0   the curator counted zero
#   head_group or chemical says     "reference compound/ligand/probe",
#                                   "not a PFAS", "endogenous",
#                                   "pesticide", "... control"
#   the name is a plain fatty acid  octanoic acid, palmitic acid, ...
#                                   ANCHORED, so "Perfluorooctanoic
#                                   acid" and "2H,2H,3H,3H-perfluoro-
#                                   octanoic acid" cannot match
#
# Across the whole file that drops 19 distinct names -- fatty acids,
# two fluorescent probes, two pesticides and warfarin -- and no PFAS.
# Every dropped row is named on stderr, so a wrong call is visible
# rather than silent.
SAYS = re.compile(r"reference (compound|ligand|probe)|not a PFAS|endogenous|"
                  r"pesticide|\bcontrol\b", re.I)
FATTY = re.compile(r"^(n-)?(formic|acetic|propionic|butyric|valeric|caproic|"
                   r"enanthic|caprylic|pelargonic|capric|lauric|myristic|"
                   r"palmitic|stearic|oleic|linoleic|arachidonic|"
                   r"butanoic|pentanoic|hexanoic|heptanoic|octanoic|nonanoic|"
                   r"decanoic|undecanoic|dodecanoic|tridecanoic|"
                   r"tetradecanoic|hexadecanoic|octadecanoic)\s+acid", re.I)


def is_pfas(r):
    chem = (r["chemical"] or "").strip()
    return not (r["n_perfluorinated_carbons"].strip() == "0"
                or SAYS.search(r["head_group"] or "")
                or SAYS.search(chem)
                or FATTY.match(chem))


BRIDGE = 4.1774                       # log10(1e9 / 66463), albumin
TO_UM = {"mM": 1e3, "um": 1.0, "uM": 1.0, "nM": 1e-3, "M": 1e6}

HDR = ["PFAS", "Species", "Site", "Kd (uM)", "n used",
       "log K (L water / kg protein)", "Route", "Source",
       "PFAS conc (uM)", "Protein conc (uM)", "ligand:protein ratio",
       "f_bound (%)"]


def to_uM(value, units):
    """A concentration and its unit string, in uM, or '' if not usable."""
    if not value or not str(value).strip():
        return ""
    u = (units or "").strip().split()[0] if (units or "").strip() else ""
    f = TO_UM.get(u, TO_UM.get(u.lower()))
    try:
        return float(value) * f if f else ""
    except (TypeError, ValueError):
        return ""


def top_of(text, units):
    """Highest concentration in a range cell, in uM."""
    if not text or "not stated" in text.lower():
        return ""
    head = text.split("(")[0].strip()
    head = re.sub(r"\b(mM|uM|nM|M)\b", "", head).strip()
    m = re.match(r"\s*([\d.eE+-]+)\s*(?:-|to|–)\s*([\d.eE+-]+)\s*$", head)
    if m:
        try:
            hi = max(float(m.group(1)), float(m.group(2)))
        except ValueError:
            return ""
    else:
        m2 = re.match(r"\s*([\d.eE+-]+)\s*$", head)
        if not m2:
            return ""
        hi = float(m2.group(1))
    return to_uM(hi, units)


def n_used(cell):
    """(n, route).  A cell that describes n rather than stating it is not
    a number: the Hill-coefficient cells say so in words, and a Hill
    exponent is not a site count."""
    s = (cell or "").strip()
    if re.fullmatch(r"[\d.]+", s):
        return float(s), "n as reported by the paper"
    return 1.0, "n = 1 assumed; the paper reports no usable stoichiometry"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study", help="study label, or a prefix of one")
    ap.add_argument("year")
    ap.add_argument("-o", "--out", default="",
                    help="output path; '-' for stdout")
    ap.add_argument("--pfas-only", action="store_true",
                    help="drop non-PFAS reference ligands (octanoic acid, "
                         "warfarin, DAUDA and the rest)")
    ap.add_argument("--cols", default="",
                    help="comma-separated subset of the column headers, "
                         "in the order wanted")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    hit = [r for r in rows
           if r["study"].startswith(args.study) and r["year"] == args.year]
    if not hit:
        print(f"  no rows for {args.study} {args.year}", file=sys.stderr)
        return 1
    labels = sorted({r["study"] for r in hit})
    if len(labels) > 1:
        print(f"  {len(labels)} study labels match: {labels}", file=sys.stderr)

    src = f"{first_author(labels[0])} {args.year}"
    out, skipped, nonpfas = [], 0, []
    for r in hit:
        if args.pfas_only and not is_pfas(r):
            nonpfas.append(r["chemical"])
            continue
        if not r["kd"]:
            skipped += 1
            continue
        try:
            kd = float(r["kd"])
        except ValueError:
            skipped += 1
            continue
        n, route = n_used(r["n_binding_sites"])
        logk = BRIDGE + math.log10(n) - math.log10(kd) if kd > 0 else ""
        pf = top_of(r["ligand_conc_range"], r["ligand_conc_units"])
        pr = to_uM(r["protein_conc"], r["protein_conc_units"])
        fb = r["unbound_fraction"].strip()
        out.append([
            r["chemical"], r["protein_species"],
            r["binding_site_identified"] or "",
            f"{kd:.4g}", f"{n:g}",
            f"{logk:.3f}" if logk != "" else "", route, src,
            f"{pf:g}" if pf != "" else "",
            f"{pr:g}" if pr != "" else "",
            r["ligand_over_protein_max"],
            fb or "not reported"])
    out.sort(key=lambda x: (x[0], x[2]))

    hdr, pick = HDR, list(range(len(HDR)))
    if args.cols:
        want = [c.strip() for c in args.cols.split(",") if c.strip()]
        miss = [c for c in want if c not in HDR]
        if miss:
            print(f"  no such column(s): {miss}\n  have: {HDR}",
                  file=sys.stderr)
            return 1
        pick = [HDR.index(c) for c in want]
        hdr = want
    body = "\n".join(["\t".join(hdr)]
                     + ["\t".join(c[i] for i in pick) for c in out]) + "\n"
    print(f"  {len(hit)} rows for {labels[0]} {args.year}; "
          f"{len(out)} carry a Kd, {skipped} do not", file=sys.stderr)
    if nonpfas:
        u = sorted(set(nonpfas))
        print(f"  {len(nonpfas)} non-PFAS reference rows dropped "
              f"(--pfas-only): {', '.join(u)}", file=sys.stderr)

    if args.out == "-":
        sys.stdout.write(body)
        return 0
    name = args.out or (re.sub(r"[^a-z0-9]+", "",
                               first_author(labels[0]).lower())
                        + args.year + "_for_excel.tsv")
    path = name if os.path.isabs(name) else os.path.join(OUT, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    print(f"  wrote {os.path.relpath(path, OUT)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
