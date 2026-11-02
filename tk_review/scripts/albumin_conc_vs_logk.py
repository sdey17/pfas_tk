#!/usr/bin/env python3
"""Does the PFAS concentration explain the spread in albumin log K?

THE HYPOTHESIS, which is a good one. Albumin coefficients for the same
compound scatter by two to four log across the literature. Binding is
saturable, so a coefficient measured at a high ligand loading is
genuinely smaller than the same coefficient measured at a low one --
Qin measured that directly and put it at roughly -0.3 log per decade of
free concentration. If the labs that disagree were simply dosing at
different concentrations, the spread would be an artefact of the
operating point rather than a disagreement about the chemistry.

HOW IT IS TESTED. Compound by compound, because compounds differ in
affinity by orders of magnitude and pooling them would measure the
chain-length trend instead. Within each compound, log K is regressed on
log10 of the assay's PFAS concentration across the sources that report
one. The slopes are then compared with Qin's own measured slopes, which
are the only direct measurement of this effect on albumin and the
benchmark any cross-laboratory slope has to match if concentration is
really the mechanism.

WHAT IT ADDS TO THE TABLE. Four columns on every row: the paper's own
wording of the PFAS concentration, that concentration in micromolar
where unit arithmetic alone reaches one, the protein concentration, and
the highest ligand:protein ratio the assay reached. The wording column
is the fuller of the two concentration columns and is the one to read
when the number is blank -- a mass concentration or a bare molar ratio
is not converted, because converting it would need a molecular weight,
a volume or a dilution the paper did not state.

Reads an albumin log K workbook (the shareable extract, or an edited
copy), writes it back with the four columns added, and prints the test.

Run:  python3 -I tk_review/scripts/albumin_conc_vs_logk.py <workbook.xlsx>
"""
import collections
import importlib.util
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")

NEW = ["PFAS conc (as stated)", "PFAS conc (uM)", "Protein conc (uM)",
       "ligand:protein ratio (max)"]

# Qin 2024 SI Table S6: the slope of log D_BSA/w against log Cw, measured
# within one laboratory on one isotherm. This is the benchmark.
QIN_SLOPES = {"PFBA": -0.298, "PFHxA": -0.271, "PFHpA": -0.305,
              "PFOA": -0.314, "PFNA": -0.147, "HFPO-DA": -0.493,
              "PFHxS": -0.472, "PFOS": -0.379, "6:2 FTSA": -0.092,
              "PFOSA": -0.105}
# Free concentration each Qin row in the extract was evaluated at.
QIN_CW = {"PFBA": 3.79, "PFHxA": 1.134, "PFHpA": 0.2344, "PFOA": 0.05577,
          "PFNA": 0.06658, "PFUnDA": 0.05805, "HFPO-DA": 1.141,
          "PFHxS": 0.03635, "PFOS": 0.01085, "6:2 FTSA": 0.3965,
          "PFOSA": 0.1278}


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def is_serum_albumin(protein):
    p = protein.lower()
    if "serum albumin" not in p:
        return False
    return not any(w in p for w in ("globulin", "plasma", "whole serum",
                                    "serum)", "immobilis", "immobiliz"))


def conditions():
    """(compound, species, source) -> the assay's concentrations.

    Keyed on the assay rather than on the individual value, because a
    concentration is a property of the experiment. Fifteen of 346 groups
    carry more than one setting; the one that reaches a molar number
    wins, so the quantitative column is filled wherever the source
    allows it.
    """
    B = _load("baf", "build_binding_affinity_table.py")
    N = _load("pfas_names", "pfas_names.py")
    recs = B.build(filter_compounds=False)[1]
    blocks = B.condition_index(recs)
    B.stamp_conditions(recs, blocks)
    out = {}
    for r in recs:
        if not is_serum_albumin(r["protein"]):
            continue
        if r["measure"] not in ("Kd", "Ka", "log D"):
            continue
        c = blocks[B.block_key(B.condition_row(r))]
        k = (N.canon(r["pfas"]), r["species"], r["source"])
        cand = (c["pfas_conc_as_stated"], str(c["pfas_conc_uM"]),
                str(c["protein_conc_uM"]),
                str(c["ligand_protein_ratio_max"]))
        if k not in out or (not out[k][1].strip() and cand[1].strip()):
            out[k] = cand
    return out


def fit(xs, ys):
    """Least squares slope, intercept, r2. None if it cannot be fitted."""
    n = len(xs)
    if n < 3 or len(set(xs)) < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0:
        return None
    b = sxy / sxx
    a = my - b * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    return b, a, r2


def main():
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    src = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    wb = load_workbook(src)
    ws = wb["Coefficients"]
    head = [c.value for c in ws[1]]
    ix = {h: i for i, h in enumerate(head)}
    cond = conditions()

    rows, filled, qin_rows = [], 0, 0
    for r in ws.iter_rows(min_row=2, values_only=True):
        d = {h: r[i] for i, h in enumerate(head)}
        key = (str(d["PFAS"]), str(d["Species"]), str(d["Source"]))
        got = cond.get(key) or cond.get(
            (str(d["PFAS"]), str(d["Species"]),
             str(d["Source"]) + " (plasma)"))
        stated, uM, prot, ratio = got or ("", "", "", "")
        # Qin's rows are the Table S6 fit at a chosen free concentration,
        # not a dosed one, so its own Cw is the honest entry.
        if str(d["Source"]).startswith("Qin 2024") \
                and str(d["PFAS"]) in QIN_CW:
            cw = QIN_CW[str(d["PFAS"])]
            stated = (f"free Cw {cw:g} uM, from 10 uM nominal in 10% FBS "
                      "(not a dosed concentration)")
            uM = str(cw)
            qin_rows += 1
        d[NEW[0]], d[NEW[1]] = stated, uM
        d[NEW[2]], d[NEW[3]] = prot, ratio
        if str(uM).strip():
            filled += 1
        rows.append(d)

    # ---- write the enriched workbook --------------------------------
    for j, name in enumerate(NEW):
        ws.cell(row=1, column=len(head) + 1 + j, value=name)
    for i, d in enumerate(rows, start=2):
        for j, name in enumerate(NEW):
            v = d[name]
            if name != NEW[0] and str(v).strip():
                try:
                    v = float(v)
                except ValueError:
                    pass
            ws.cell(row=i, column=len(head) + 1 + j, value=v)
    for j, name in enumerate(NEW):
        c = ws.cell(row=1, column=len(head) + 1 + j)
        c.fill = PatternFill("solid", fgColor="2A78D6")
        c.font = Font(bold=True, color="FFFFFF")
        c.alignment = Alignment(vertical="center", wrap_text=True)
        w = max([len(name)] + [len(str(d[name])) for d in rows])
        ws.column_dimensions[get_column_letter(len(head) + 1 + j)].width = \
            min(max(w + 2, 12), 52)
    ws.auto_filter.ref = ws.dimensions
    wb.save(OUT)

    print(f"{len(rows)} rows; a molar PFAS concentration on {filled} "
          f"({qin_rows} of them Qin's own free Cw)")
    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print()

    # ---- the test ----------------------------------------------------
    by = collections.defaultdict(list)
    for d in rows:
        u = str(d[NEW[1]]).strip()
        if not u:
            continue
        try:
            c = float(u)
        except ValueError:
            continue
        if c <= 0:
            continue
        by[str(d["PFAS"])].append(
            (math.log10(c), float(d["log K (L water / kg protein)"]),
             str(d["Source"])))

    print("DOES CONCENTRATION EXPLAIN THE SPREAD?")
    print("Within each compound, log K regressed on log10(PFAS conc).")
    print("Negative and near Qin's slope would mean yes.")
    print()
    print(f"  {'PFAS':10s}{'n':>3s}{'sources':>8s}{'conc span':>11s}"
          f"{'logK span':>11s}{'slope':>9s}{'r2':>7s}   Qin's own slope")
    slopes = []
    for c in sorted(by, key=lambda k: -len(by[k])):
        pts = by[c]
        if len(pts) < 4:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        f = fit(xs, ys)
        if not f:
            continue
        b, _a, r2 = f
        slopes.append(b)
        q = QIN_SLOPES.get(c)
        print(f"  {c:10s}{len(pts):3d}{len({p[2] for p in pts}):>8d}"
              f"{max(xs) - min(xs):>10.1f}d{max(ys) - min(ys):>11.2f}"
              f"{b:>9.2f}{r2:>7.2f}   "
              f"{('%.2f' % q) if q is not None else '-'}")
    print()
    if slopes:
        print(f"  median cross-laboratory slope {statistics.median(slopes):+.2f} "
              f"log per decade, over {len(slopes)} compounds")
        qs = list(QIN_SLOPES.values())
        print(f"  median Qin within-laboratory slope "
              f"{statistics.median(qs):+.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
