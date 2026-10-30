#!/usr/bin/env python3
"""Write the four combined datasets into one workbook, one sheet per table.

The CSVs in db/combined/ are the authoritative form; this is a convenience
export for opening alongside the original PFAS_TK.xlsx.

Run:  python3 scripts/export_combined_xlsx.py
"""
import collections
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
# the two analysis outputs -- the converted-against-measured comparison
# (scripts/logk_converted_vs_measured.py) and the flat albumin list
# (scripts/albumin_kd_to_k.py) -- which sit in db/.
SHEETS = [
    ("tk_parameters.csv", "TK parameters", SRC),
    ("binding.csv", "Protein binding", SRC),
    ("binding_conditions.csv", "Binding conditions", SRC),
    ("binding_basis.csv", "Conversion basis", SRC),
    ("albumin_kd_to_k.csv", "Albumin Kd to K", DB),
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


def _albumin():
    with open(os.path.join(DB, "albumin_kd_to_k.csv")) as fh:
        a = list(csv.DictReader(fh))
    c = collections.Counter(r["convertible"] for r in a)
    return {"alb": len(a), "alb_src": len({r["source"] for r in a}),
            "alb_yes": c["yes"], "alb_assume": c["assume n = 1"],
            "alb_no": c["no"],
            "alb_stated": sum(1 for r in a
                              if r["n_sites_as_stated"].strip())}


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
        **_pairs(), **_albumin(),
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
        ["That sheet also carries the LOADING question for each pair, which"],
        ["is what settles whether a residual is a real disagreement. A"],
        ["converted value is already at zero loading, because a fitted K_d"],
        ["defines it there; a measured one was determined AT its own"],
        ["loading and so sits below its own zero-loading limit by"],
        ["log10(1 + C_w/K_d). `measured_loading_max_log` bounds that from"],
        ["above, `measured_loading_how` says by which route, and"],
        ["`residual_at_zero_loading` is the residual with both sides put"],
        ["on the same footing. All 29 admitted pairs are measured by"],
        ["Fischer 2024, which held its PFAS:protein molar ratio below"],
        ["0.004, so the bound is 0.0017 log and the median residual moves"],
        ["from -0.810 to -0.812. Loading accounts for 0.2% of that gap:"],
        ["it is ruled out as the explanation, not merely unexamined."],
        ["`converted_isotherm_at_top` is the matching diagnostic on the"],
        ["other side -- not a correction, since nothing shifts a K0, but"],
        ["a statement of whether the K_d was fitted inside the data."],
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
        ["ALBUMIN, FLAT. `Albumin Kd to K` is the one sheet here that"],
        ["needs no lookup, and it exists because the normalisation above"],
        ["costs two of them to answer the question people actually arrive"],
        [f"with: can I convert THIS number. All {st['alb']} serum-albumin"],
        [f"affinities from {st['alb_src']} sources sit on one flat row"],
        ["each, with the protein concentration, the PFAS concentration,"],
        ["the ligand:protein ratio and the stoichiometry repeated inline."],
        ["It is not a subset of `Protein binding`: it is restricted to"],
        ["albumin -- HSA, BSA and the other mammalian albumins -- so the"],
        ["FABPs, the PPAR receptors, globulin, plasma, serum and whole"],
        ["cells are absent, and it spans every albumin source rather than"],
        ["the top ten."],
        [],
        [f"Its verdict column reads three ways. {st['alb_yes']} rows are"],
        ["`yes`: the paper gives a dissociation constant AND a"],
        [f"stoichiometry, and `logK` is filled. {st['alb_assume']} read"],
        ["`assume n = 1`: a constant and no usable site count, so one"],
        ["assumption short. Those carry the number in a SEPARATE column,"],
        ["`logK_if_n_assumed_1`, and `logK` stays blank, so the two kinds"],
        ["cannot be averaged together by accident. This does not reverse"],
        ["the refusal stated above -- `Protein binding` is unchanged and"],
        ["still derives no K from an affinity alone. The column measures"],
        [f"what the refusal costs and where. The last {st['alb_no']} read"],
        ["`no`: no single constant to start from, and 23 of those are"],
        ["negative results, where the paper looked and found no measurable"],
        ["binding, rather than gaps in reporting."],
        [],
        ["WHERE EACH ASSAY SAT ON THE ISOTHERM. A converted value is the"],
        ["coefficient at ZERO loading, because that is what a fitted K_d"],
        ["defines. Every row with a molar PFAS concentration states a"],
        ["titration series starting at or near zero, so the K_d is a model"],
        ["parameter and that extrapolation is the paper's own, not an"],
        ["extra approximation laid on top of it. What is worth knowing is"],
        ["where the experiment actually operated, since a coefficient"],
        ["measured at real loading is log10(1 + C_w/K_d) BELOW the"],
        ["reported value -- which is exactly the gap when a converted"],
        ["number meets someone else's measured one. Five columns say so,"],
        ["all at the top of the titration: `kd_over_top_conc`,"],
        ["`cw_free_at_top_uM` (free ligand by mass balance on the paper's"],
        ["own concentrations, not the dosed total), `cw_over_kd_at_top`,"],
        ["`saturation_at_top_log` and `isotherm_at_top`. BOTH TAILS ARE"],
        ["FLAGS: `saturating` means the fit saw curvature so K_d is"],
        ["constrained, but the coefficient at that loading sits far below"],
        ["the reported one; `linear to the top` means the series never"],
        ["reached K_d and the fit is an extrapolation beyond the data."],
        ["84 of the 133 computable rows have a K_d above the highest"],
        ["concentration their own series tested."],
        [],
        [f"`n_sites_as_stated` ({st['alb_stated']} rows) is the paper's own"],
        ["wording where a number could not be read from it -- `6-9`, `1-3"],
        ["primary sites`, `up to 8 occupied sites at 4:1`, `n_Hill"],
        ["0.69-2.5`. None of those is one number, so none converts; they"],
        ["are carried because a reader choosing what n to assume should"],
        ["see what the authors claimed."],
        [],
        ["Sheet", "Rows", "What it holds"],
    ], 1):
        readme.append(line)
    readme["A1"].font = Font(bold=True, size=14)
    # The index header is the last line appended, wherever that lands.
    head = readme.max_row
    for col in "ABC":
        readme[f"{col}{head}"].font = Font(bold=True)

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
