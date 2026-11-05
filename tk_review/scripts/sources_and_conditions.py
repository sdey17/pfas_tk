#!/usr/bin/env python3
"""Every source in the database: full citation, and where its
concentrations come from.

TWO QUESTIONS IN ONE FILE.

  WHAT IS THE PAPER?  A full citation for every study the database
  draws on, built from PubMed records in db/citations.csv.

  WHERE DID THE CONCENTRATIONS COME FROM?  For each study, the PFAS
  and protein concentrations exactly as the paper words them, and the
  table, figure or section they were read out of.

WHY THE SECOND MATTERS.  A concentration decides whether a binding
constant means anything -- it is what the range check and the loading
axis are built on -- and this database deliberately stores the
paper's own words rather than a parsed number, because "0.33 to
153.5 ng/mL (ligand:protein mole ratio 0.05-5)", "see Text S7" and
"not quantified in molar terms" are three different situations that a
numeric column would flatten into one blank.  This file keeps them
apart and says, per study, which is which.

STATED means the cell contains a NUMBER.  A cell can carry the
paper's own words and still report no quantity: "see Text S7", "not
stated for the 26 standards", "a range of concentrations per
compound". Those count as NOT stated, and the location column says
where the paper was looked at so the gap can be closed by reading
rather than by guessing.

Citations: db/citations.csv, from PubMed via the project's MCP
client.  63 of 71 studies have a PubMed record; the 8 without are
listed with whatever identifier exists, and are reports, pre-1960
papers or journals PubMed does not index.

Run:  python3 tk_review/scripts/sources_and_conditions.py
Writes db/sources_and_conditions.xlsx
"""
import collections
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
CITES = os.path.join(HERE, "..", "db", "citations.csv")
OUT = os.path.join(HERE, "..", "db", "sources_and_conditions.xlsx")

# A number, not a word. Matches 1, 0.5, 1e-3, 1.5-3, 10 to 20.
HAS_NUM = re.compile(r"\d")
# Phrases that carry the paper's wording but no quantity.
NO_NUMBER = re.compile(
    r"not (stated|reported|quantified|given|specified)|"
    r"see (Table|Text|Figure|SI|Supporting)|not applicable|"
    r"a range of concentrations|unknown", re.I)


def stated(cell):
    """True when the cell actually reports a quantity."""
    c = (cell or "").strip()
    if not c:
        return False
    if NO_NUMBER.search(c) and not re.match(r"^\s*[\d.]", c):
        return False
    return bool(HAS_NUM.search(c))


def location(source_table):
    """The table / figure / section the row was read from."""
    s = (source_table or "").split("|")[0].strip()
    return s[:120] or "(not recorded)"


def main():
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    cites = {(c["study"], c["year"]): c
             for c in csv.DictReader(open(CITES, encoding="utf-8"))}

    wb = Workbook()

    # ---- sheet 1: one row per study -------------------------------
    ws = wb.active
    ws.title = "Sources"
    hdr = ["Study", "Year", "Full citation", "DOI", "PMID",
           "Rows in database", "Proteins", "PFAS conc stated",
           "Protein conc stated", "Both stated", "Neither stated"]
    ws.append(hdr)
    by = collections.defaultdict(list)
    for r in rows:
        by[(r["study"], r["year"])].append(r)
    for key in sorted(by, key=lambda k: (k[0], k[1])):
        rs = by[key]
        c = cites.get(key, {})
        lig = sum(1 for r in rs if stated(r["ligand_conc_range"]))
        pro = sum(1 for r in rs if stated(r["protein_conc"]))
        both = sum(1 for r in rs if stated(r["ligand_conc_range"])
                   and stated(r["protein_conc"]))
        neither = sum(1 for r in rs if not stated(r["ligand_conc_range"])
                      and not stated(r["protein_conc"]))
        prots = sorted({r["protein"][:40] for r in rs})
        ws.append([key[0], key[1], c.get("citation", ""),
                   c.get("doi", ""), c.get("pmid", ""), len(rs),
                   "; ".join(prots)[:300],
                   f"{lig}/{len(rs)}", f"{pro}/{len(rs)}",
                   f"{both}/{len(rs)}", f"{neither}/{len(rs)}"])

    # ---- sheet 2: where each concentration came from --------------
    cs = wb.create_sheet("Where the concentrations are")
    hdr2 = ["Study", "Year", "Protein", "PFAS conc (the paper's "
            "words)", "stated?", "Protein conc (the paper's words)",
            "stated?", "ligand:protein ratio (as reported)",
            "Read from", "Rows"]
    cs.append(hdr2)
    seen = collections.Counter()
    first = {}
    for r in rows:
        k = (r["study"], r["year"], r["protein"][:40],
             r["ligand_conc_range"][:150], r["protein_conc"][:80],
             r["ligand_over_protein_max"], location(r["source_table"]))
        seen[k] += 1
        first.setdefault(k, r)
    for k in sorted(seen, key=lambda k: (k[0], k[1], k[2])):
        r = first[k]
        u_l = r["ligand_conc_units"].strip()
        u_p = r["protein_conc_units"].strip()
        cs.append([
            k[0], k[1], k[2],
            (k[3] + (" " + u_l if u_l and k[3] else "")) or "(empty)",
            "yes" if stated(k[3]) else "no",
            (k[4] + (" " + u_p if u_p and k[4] else "")) or "(empty)",
            "yes" if stated(k[4]) else "no",
            k[5] or "", k[6], seen[k]])

    bold = Font(bold=True)
    wrap = Alignment(vertical="top", wrap_text=True)
    for sh, head in ((ws, hdr), (cs, hdr2)):
        for j, name in enumerate(head):
            sh.cell(row=1, column=j + 1).font = bold
            wide = ("citation" in name.lower() or "words" in name
                    or name in ("Proteins", "Read from"))
            sh.column_dimensions[get_column_letter(j + 1)].width = (
                70 if "citation" in name.lower() else
                42 if wide else 16)
            if wide:
                for row in sh.iter_rows(min_row=2):
                    row[j].alignment = wrap
        sh.freeze_panes = "A2"

    wb.save(OUT)
    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print(f"  Sources: {len(by)} studies")
    nocite = [k for k in by if not cites.get(k, {}).get("citation")]
    print(f"    {len(by)-len(nocite)} with a full citation, "
          f"{len(nocite)} without:")
    for k in sorted(nocite):
        print(f"      {k[0]} {k[1]}")
    print(f"  Where the concentrations are: {len(seen)} distinct "
          f"study / protein / condition combinations")
    lig = sum(1 for r in rows if stated(r["ligand_conc_range"]))
    pro = sum(1 for r in rows if stated(r["protein_conc"]))
    both = sum(1 for r in rows if stated(r["ligand_conc_range"])
               and stated(r["protein_conc"]))
    print(f"\n  over all {len(rows)} rows, a NUMBER is given for:")
    print(f"    PFAS concentration     {lig:5d}  ({100*lig/len(rows):.0f} %)")
    print(f"    protein concentration  {pro:5d}  ({100*pro/len(rows):.0f} %)")
    print(f"    both                   {both:5d}  ({100*both/len(rows):.0f} %)")


if __name__ == "__main__":
    main()
