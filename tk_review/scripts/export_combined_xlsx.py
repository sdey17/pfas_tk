#!/usr/bin/env python3
"""Write the four combined datasets into one workbook, one sheet per table.

The CSVs in db/combined/ are the authoritative form; this is a convenience
export for opening alongside the original PFAS_TK.xlsx.

Run:  python3 scripts/export_combined_xlsx.py
"""
import csv
import os

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "combined")
DB = os.path.join(HERE, "..", "db")
DEST = os.path.join(HERE, "..", "db", "combined", "PFAS_TK_combined.xlsx")

# (file, sheet title, directory). Everything lives in db/combined/ except
# the converted-against-measured comparison, which is an analysis output of
# scripts/logk_converted_vs_measured.py and sits in db/.
SHEETS = [
    ("tk_parameters.csv", "TK parameters", SRC),
    ("binding.csv", "Protein binding", SRC),
    ("binding_conditions.csv", "Binding conditions", SRC),
    ("binding_basis.csv", "Conversion basis", SRC),
    ("logk_converted_vs_measured.csv", "Converted vs measured", DB),
    ("transporters.csv", "Transporters", SRC),
    ("regulatory.csv", "Regulatory", SRC),
]
HEAD_FILL = PatternFill("solid", fgColor="2A78D6")


def _pairs():
    with open(os.path.join(DB, "logk_converted_vs_measured.csv")) as fh:
        p = list(csv.DictReader(fh))
    ex = sum(1 for r in p if r.get("excluded", "").strip())
    return {"pairs": len(p) - ex, "pairs_excluded": ex}


def binding_stats():
    """Counts for the README, read off the files rather than typed in.

    A hand-written count in a generated file is a count that drifts; these
    three sentences are the ones a reader is most likely to check.
    """
    with open(os.path.join(SRC, "binding.csv")) as fh:
        b = list(csv.DictReader(fh))
    with open(os.path.join(SRC, "binding_conditions.csv")) as fh:
        c = list(csv.DictReader(fh))
    with open(os.path.join(SRC, "binding_basis.csv")) as fh:
        bas = list(csv.DictReader(fh))
    keyed = [r for r in b if r["condition_id"].strip()]
    return {
        "rows": len(b), "keyed": len(keyed), "keyless": len(b) - len(keyed),
        "assays": len(c), "routes": len(bas),
        "ratio": sum(1 for r in b if r["ligand_protein_ratio_max"].strip()),
        "nused": sum(1 for r in b if r["n_used"].strip()),
        "kd": sum(1 for r in b if r["kd_uM"].strip()),
        "logk": sum(1 for r in b if r["logK"].strip()),
        **_pairs(),
    }


def main():
    st = binding_stats()
    wb = Workbook()
    wb.remove(wb.active)

    readme = wb.create_sheet("README")
    for i, line in enumerate([
        ["PFAS toxicokinetics review - combined datasets"],
        [],
        ["Built by scripts/build_combined_datasets.py from the per-pass tables in"],
        ["db/ and the per-paper extractions in db/primary_2026/."],
        [],
        ["Every row carries a `provenance` column naming the file it came from, and"],
        ["a study / pmid_or_doi / source_table trio naming the primary source."],
        ["Nothing here is quoted from memory."],
        [],
        ["COMPOUND NAMES. `chemical` is the source's own spelling, left"],
        ["untouched because it is provenance. `chemical_canonical` is one"],
        ["name per compound across all four sheets -- GROUP AND COUNT ON"],
        ["THAT ONE. The same PFAS is written differently by different"],
        ["groups (GenX and HFPO-DA; PFUdA, PFUnA and PFUnDA; 6:2 FtS and"],
        ["6:2 FTSA; PFTA and PFTeDA), and counting raw names inflates"],
        ["every total. 256 raw names across these files reduce to 137"],
        ["compounds. The full mapping, the evidence for each merge and"],
        ["the groups deliberately left apart -- linear and branched"],
        ["isomers against the technical mixture, the potassium and"],
        ["ammonium salts that one study doses alongside the acids -- are"],
        ["in db/pfas_name_map.csv, built by scripts/pfas_names.py."],
        [],
        ["CONVERTED BINDING VALUES. `Protein binding` reports each paper's"],
        ["own currency in value/units, and beside it the two standardised"],
        ["ones this project runs on: `kd_uM`, a dissociation constant in"],
        ["micromolar, and `logK`, log10 of the protein-water partition"],
        [f"coefficient in L_water per kg protein. Of {st['rows']} rows, "
         f"{st['kd']} reach a"],
        [f"kd_uM and {st['logk']} reach a logK. Only exact routes are"],
        ["used. A Kd becomes a K only where the SAME paper also reports a"],
        ["stoichiometry, through log K = log10(1e9/M_protein) + log n -"],
        ["log Kd(uM); it is never derived from an affinity alone, because"],
        ["K fixes only the ratio n/Kd and an assumed n = 1 would be a"],
        [f"guess dressed as a number. `n_used` ({st['nused']} rows) is the"],
        ["stoichiometry each conversion actually used, as a number rather"],
        ["than a digit buried in a sentence."],
        [],
        ["HOW EACH ONE WAS REACHED. `basis_id` keys into `Conversion"],
        [f"basis`, {st['routes']} rows covering every measurement in the"],
        ["sheet. Neither route there is ever blank: where a converted value"],
        ["is blank, the route says why. That is the only way to tell a"],
        ["quantity the paper did not report from one this"],
        [f"project refuses to infer. It is keyed on the ROUTE and not on"],
        ["the study, because the basis is not a study constant: only 28 of"],
        ["the 50 sources use a single route pair, since the route follows"],
        ["what each individual ROW reports, not what the paper is."],
        [],
        [f"WHICH REGIME. `ligand_protein_ratio_max` ({st['ratio']} rows) is"],
        ["the highest ratio the assay reached, most ligand over least"],
        ["protein. Below about 0.01 a measurement sees the specific site;"],
        ["above about 1 it sees non-specific partitioning, and the two"],
        ["differ by orders of magnitude, so two constants measured at"],
        ["different ratios are not comparable even for the same compound"],
        ["and protein. It is a property of the assay, so it is also on"],
        ["`Binding conditions`; it is repeated on every measurement row"],
        ["because it is the column you filter on, and a lookup to learn"],
        ["which regime a number came from is a lookup too many."],
        [],
        [f"`Converted vs measured` holds the {st['pairs']} admitted"],
        ["cases where a converted value and an independently measured one"],
        ["exist for the same compound, protein and species, with the"],
        [f"residual between them, plus {st['pairs_excluded']} marked in an"],
        ["`excluded` column. Those are conversions this database refuses:"],
        ["their stoichiometry is the slope of a fluorescence"],
        ["double-logarithmic or Scatchard plot, not a counted site, and"],
        ["three such papers scatter across -4.01 to +0.13 log. They are"],
        ["computed and kept because they are the evidence for the rule,"],
        ["and they are in no statistic and in no figure."],
        [],
        ["ASSAY CONDITIONS. The concentrations, temperature and buffer are"],
        ["properties of the experiment, not of the individual number, so"],
        [f"they sit once in `Binding conditions` -- {st['assays']} assays "
         f"behind the"],
        [f"{st['keyed']} measurements that carry a `condition_id` key "
         f"(the other"],
        [f"{st['keyless']} come from per-paper extractions with their own"],
        ["schemas). Join on that key, as a lookup: many measurements to one"],
        ["assay. For 37 of the 50 sources the setting is a single study"],
        ["constant; for 13 it is not, and for Qin 2024, Bischel 2011 and"],
        ["Gao 2025 varying it IS the experiment, which is why the key"],
        ["points at an assay rather than at a study."],
        [],
        ["Sheet", "Rows", "What it holds"],
    ], 1):
        readme.append(line)
    readme["A1"].font = Font(bold=True, size=14)
    readme["A10"].font = readme["B10"].font = readme["C10"].font = Font(bold=True)

    for fn, title, where in SHEETS:
        path = os.path.join(where, fn)
        with open(path) as fh:
            rows = list(csv.reader(fh))
        ws = wb.create_sheet(title)
        for r in rows:
            ws.append(r)
        for c in range(1, len(rows[0]) + 1):
            cell = ws.cell(row=1, column=c)
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = HEAD_FILL
            cell.alignment = Alignment(vertical="center", wrap_text=True)
            width = max((len(str(rw[c - 1])) for rw in rows[:400]
                         if c - 1 < len(rw)), default=10)
            ws.column_dimensions[get_column_letter(c)].width = min(max(width + 2, 10), 46)
        ws.freeze_panes = "B2" if title == "Binding conditions" else "A2"
        ws.auto_filter.ref = ws.dimensions
        readme.append([title, len(rows) - 1, fn])

    wb.save(DEST)
    print("wrote", os.path.relpath(DEST, os.path.join(HERE, "..")))


if __name__ == "__main__":
    main()
