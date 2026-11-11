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
import csv
import collections
import importlib.util
import math
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")
OUT = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

NEW = ["PFAS conc (as stated)", "PFAS conc (uM)", "Protein conc (uM)",
       "ligand:protein ratio (max)", "f_bound (%)"]

# The bound fraction travels with the VALUE, not with the assay block --
# two compounds measured side by side in one experiment share a
# concentration but not an f -- so it is joined per row rather than taken
# off the condition block the other four come from.
BOUND = re.compile(r"([\d.]+)\s*%\s*bound")
UNBOUND = re.compile(r"([\d.]+)\s*%\s*unbound")
# "99.95 +/- 0.01% bound" -- the regex would otherwise capture the standard
# error, because it is the number sitting against "% bound".
PM = re.compile(r"\s*(?:\+/-|\u00b1)\s*[\d.]+")

# Qin 2024 SI Table S6: the slope of log D_BSA/w against log Cw, measured
# within one laboratory on one isotherm. This is the benchmark.
QIN_SLOPES = {"PFBA": -0.298, "PFHxA": -0.271, "PFHpA": -0.305,
              "PFOA": -0.314, "PFNA": -0.147, "HFPO-DA": -0.493,
              "PFHxS": -0.472, "PFOS": -0.379, "6:2 FTSA": -0.092,
              "PFOSA": -0.105}
# Free concentration each Qin row in the extract was evaluated at.
# The Qin anchor's concentrations come from qin2024_anchor.py, so the sheet,
# these columns and both figures cannot drift apart.  Keyed by the canonical
# PFAS name the sheet uses, which renames PFUnA to PFUnDA.
ANCHOR_MOD = "qin2024_anchor.py"
CANON = {"PFUnA": "PFUnDA"}


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
    fb = bound_index(N)
    blocks = B.condition_index(recs)
    B.stamp_conditions(recs, blocks)
    out = {}
    for r in recs:
        if not is_serum_albumin(r["protein"]):
            continue
        if r["measure"] not in ("Kd", "Ka", "log D"):
            continue
        c = blocks[B.block_key(B.condition_row(r))]
        cand = (c["pfas_conc_as_stated"], str(c["pfas_conc_uM"]),
                str(c["protein_conc_uM"]),
                str(c["ligand_protein_ratio_max"]),
                fb.get((r["_id"], N.canon(r["pfas"]), r["species"],
                        str(r.get("value", "") or "").strip()))
                or fb.get((r["_id"], N.canon(r["pfas"]), "")) or "")
        # (compound, species, source) is NOT unique: one study can measure
        # the same compound on the same protein under several conditions,
        # and keying on those three alone hands them all the conditions of
        # whichever row was seen first. Bischel 2010 measures PFNA on BSA
        # four times -- at ligand:protein 6e-4, 1e-3, 5 and 120 -- and the
        # three-key collapse stamped 6e-4 on all four, which is wrong by
        # five orders of magnitude on two of them. Kd distinguishes them
        # and is present on both sides of the join, so key on it too and
        # fall back to the old key when a row has no Kd, leaving every
        # single-condition study exactly as it was.
        # Kd distinguishes most multi-condition studies. Bischel 2010 has
        # four PFNA conditions with NO Kd -- two albumins, a salt row and a
        # defatted control -- so log K is the tiebreaker there, and it is
        # on both sides of the join. Without it all four collapsed onto
        # whichever was read first and the salt row showed the wrong
        # ligand:protein ratio, 0.0006 against its real 0.001.
        kd = str(r.get("kd_uM", "") or "").strip()
        lk = str(r.get("logK", "") or "").strip()
        for k in ([(N.canon(r["pfas"]), r["species"], r["source"], kd)]
                  if kd else []) + \
                 ([(N.canon(r["pfas"]), r["species"], r["source"],
                    "logK:" + lk)] if lk and not kd else []) + \
                 [(N.canon(r["pfas"]), r["species"], r["source"])]:
            if k not in out or (not out[k][1].strip() and cand[1].strip()):
                out[k] = cand
    return out


def bound_index(N):
    """(doi, compound, Kd) -> bound percentage, read from the database.

    The affinity records do not carry unbound_fraction, so it is joined
    from db/protein_binding.csv.

    Kd is in the key because one study can report several bound fractions
    for one compound, and they are the thing that differs: Bischel 2010
    measures PFNA on BSA three times, at 99.95, 99.92 and 99.89% bound.
    A row whose Kd cell is empty cannot be told apart that way, so it
    never writes a keyed entry.

    The loose (doi, compound) key is added ONLY where the study reports a
    single value for that compound.  Where it reports several and nothing
    distinguishes them, the column is left empty rather than stamped with
    whichever row was read first -- the failure this project already hit
    once, on Qin 2024.

    Only an explicit "N% bound" or "N% unbound" is read.  Bischel 2011
    stores a bare "0.01" whose unit the field never states, so those
    eleven rows are skipped rather than guessed at.
    """
    seen = collections.defaultdict(set)
    out = {}
    with open(DB, newline="", encoding="utf-8") as fh:
        for d in csv.DictReader(fh):
            pct = bound_pct(d.get("unbound_fraction") or "")
            if not pct:
                continue
            cmp_ = N.canon(d["chemical"])
            # The raw cell the row reports, Kd or Ka/log K_PW, is what the
            # affinity records carry as `value`, so it joins without this
            # file having to redo any unit conversion.  Bischel 2010's four
            # PFNA conditions have no Kd at all and are told apart only by
            # it: 4.93, 4.74, 4.56, 4.80.
            for raw in {str(d.get("kd", "") or "").strip(),
                        str(d.get("ka", "") or "").strip()} - {""}:
                out.setdefault((d["doi"], cmp_, d["protein_species"], raw),
                               pct)
            seen[(d["doi"], cmp_)].add(pct)
    for (doi, cmp_), vals in seen.items():
        if len(vals) == 1:
            out.setdefault((doi, cmp_, ""), next(iter(vals)))
    return out


def bound_pct(text):
    """The bound percentage as a number, from whichever way it is written.

    Sources state it either way round -- "46% bound" or "54% unbound" --
    and some give both. Prefer an explicit bound figure; otherwise take
    the complement of the unbound one. Returns "" when neither is there,
    which is most rows: only 57 albumin rows report it at all.
    """
    text = PM.sub("", text)
    m = BOUND.search(text)
    if m:
        return f"{float(m.group(1)):g}"
    m = UNBOUND.search(text)
    if m:
        return f"{100 - float(m.group(1)):g}"
    return ""


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
    _a = _load("qin_anchor", ANCHOR_MOD)
    QIN = {CANON.get(a["pfas"], a["pfas"]): a for a in _a.anchor()}

    rows, filled, qin_rows = [], 0, 0
    for r in ws.iter_rows(min_row=2, values_only=True):
        d = {h: r[i] for i, h in enumerate(head)}
        kd = str(d.get("Kd (uM)", "") or "").strip()
        key = (str(d["PFAS"]), str(d["Species"]), str(d["Source"]))
        # Prefer the Kd-qualified key so multi-condition studies get their
        # own row's conditions; fall back to the unqualified one.
        # No "(plasma)" fallback here.  It used to let a "Qin 2024" row
        # borrow the plasma paper's Table S5 conditions; the anchor module
        # supplies Qin's conditions below, from the active paper only.
        lk = str(d.get("log K (L water / kg protein)", "") or "").strip()
        got = (cond.get(key + (kd,)) if kd else None) \
            or (cond.get(key + ("logK:" + lk,)) if lk and not kd else None) \
            or cond.get(key)
        stated, uM, prot, ratio, fb = got or ("", "", "", "", "")
        # A Qin row is a fit read at a chosen free concentration, not a
        # dosed one, so it takes its whole set of conditions from the
        # anchor module -- including the ligand, protein and ratio, which
        # must come from the SAME paper as the coefficient.
        if str(d["Source"]).startswith("Qin 2024") \
                and str(d["PFAS"]) in QIN:
            a = QIN[str(d["PFAS"])]
            stated = f"{a['where']}; {a['basis']}"
            # The conc column is the DOSED ligand everywhere else in the
            # sheet, so it is the dosed one here too and the ratio divides
            # into it.  The free concentration the fit was read at is a
            # different number and lives in the text, not this column.
            uM, prot, ratio = a["ligand_uM"], a["protein_uM"], a["ratio"]
            fb = ""      # the baseline paper reports no bound fraction
            qin_rows += 1
        d[NEW[0]], d[NEW[1]] = stated, uM
        d[NEW[2]], d[NEW[3]], d[NEW[4]] = prot, ratio, fb
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
