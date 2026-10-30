#!/usr/bin/env python3
"""Every serum-albumin affinity in the database, on one flat row each.

WHY THIS EXISTS. The main workbook normalises: the assay conditions live
in one sheet, the conversion basis in another, and a measurement row
carries keys into both. That is right for a database and wrong for the
question people actually arrive with, which is "can I convert this one,
and on what". Answering it from the normalised sheets takes two lookups
and a paragraph of prose.

This file answers it on one line. Every serum-albumin row that reports an
affinity is here, with the protein concentration, the PFAS concentration,
the stoichiometry used and a verdict, all in the same row. Nothing has to
be joined and nothing has to be read.

SCOPE is serum albumin only -- HSA and BSA, and the other mammalian
albumins where a paper measured them. Gamma-globulin, plasma, serum,
liver and intestinal FABP, the PPAR receptors, membranes and whole cells
are all out, because the question was about albumin.

THREE VERDICTS, NOT TWO, because the albumin literature does not divide
in two. `convertible` reads:

  yes            the paper supplies a dissociation constant AND a
                 stoichiometry that is a counted site or a model
                 assumption it states itself. `logK` is filled.

  assume n = 1   the paper supplies a dissociation constant and no
                 usable site count. One assumption away from a number,
                 and it is the SAME assumption Chen's single-site
                 Langmuir makes by construction -- which is why Chen
                 reads `yes` and these read `assume n = 1`: there the
                 paper owns the assumption, here the reader would.
                 `logK_if_n_assumed_1` is filled and `logK` is not, so
                 the two can never be summed together by accident.

  no             there is no single dissociation constant to start from.
                 No assumption rescues these. Most are not omissions:
                 nineteen of them are negative results -- the paper
                 looked and found no measurable binding -- and the rest
                 are ranges or two-class fits, which is a value that is
                 not one number rather than a value that is missing.

`why_not_as_reported` names, in a few words, what the row is short of.

TWO STOICHIOMETRY COLUMNS. `n_sites_used` is the number that entered a
conversion, or blank. `n_sites_as_stated` is the paper's own wording,
which is populated far more often -- "6-9", "1-3 primary sites", "up to 8
occupied sites at 4:1", "n_Hill 0.69-2.5" -- and which converts nowhere,
because none of those is one number. It is here because a reader choosing
what to assume should see what the authors actually claimed.

Writes db/albumin_kd_to_k.csv.

Run:  python3 tk_review/scripts/albumin_kd_to_k.py
"""
import collections
import csv
import importlib.util
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "albumin_kd_to_k.csv")


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B = _load("baf", "build_binding_affinity_table.py")
N = _load("pfas_names", "pfas_names.py")

# log10(1e9 / 66463) -- albumin's bridge constant, per kg of protein.
BRIDGE = math.log10(1e9 / B.MW_ALBUMIN)

FIELDS = ["pfas", "pfas_canonical", "protein", "species", "source",
          "method", "measure", "value_as_reported", "units_as_reported",
          "kd_uM", "n_sites_used", "n_sites_as_stated", "n_basis",
          "convertible",
          "logK", "logK_if_n_assumed_1", "why_not_as_reported",
          "protein_conc_uM", "protein_conc_as_stated",
          "pfas_conc_uM", "pfas_conc_as_stated",
          "ligand_protein_ratio_max", "location"]

# Why a row is short of a conversion, in a few words. The full basis
# strings in the main table say it at length, which is what that table is
# for; a reader scanning 300 rows for the usable ones wants a tag.
SHORT = [
    ("the paper reports an affinity but no stoichiometry",
     "no stoichiometry reported"),
    ("is the fitted slope of an optical quenching plot",
     "only n is a quenching-plot slope, which is not a site count"),
    ("a partition coefficient; converting to Kd needs",
     "reports a coefficient, not an affinity"),
    ("an unbound fraction, not an affinity", "unbound fraction only"),
    ("a computed docking score", "docking score, not a measurement"),
    ("dissociation constants but n =",
     "site count and the list of constants do not correspond"),
    ("conversion out of range", "conversion out of range"),
]


def no_kd_reason(rec):
    """Why a row reached no single Kd. Read off the value, not the prose.

    The main table's basis string says "no Kd" for all of these, which is
    true and unhelpful: a range and a two-class fit are different
    problems with different fixes.
    """
    v = (rec["value_as_reported"] or "").strip()
    low = v.lower()
    # A negative result is not a gap in reporting. Starnes' "not detected"
    # and "nonbinding" are the measurement, and reading them as missing
    # data would turn nineteen findings into nineteen holes.
    if any(w in low for w in ("not detected", "nonbinding", "no binding",
                              "not measurable", "n.d.", "below detection")):
        return "no measurable binding -- a negative result, not a gap"
    if any(w in low for w in ("kd1", "kd2", "high-affinity set",
                              "low-affinity set", "class 1", "class 2")):
        return "two classes of site in one cell, not one constant"
    if " to " in low or low.startswith("~") or "range" in low:
        return "reported as a range or an approximation, not one value"
    if ";" in v or "," in v:
        return "several values in one cell, not one constant"
    return "no single dissociation constant reported"


def short_reason(basis, rec):
    for needle, label in SHORT:
        if needle in basis:
            return label
    if "no Kd" in basis:
        return no_kd_reason(rec)
    return basis[:70]


def is_serum_albumin(protein):
    """Albumin as the sorbent, not a matrix that happens to contain it."""
    p = protein.lower()
    if "serum albumin" not in p:
        return False
    return not any(w in p for w in ("globulin", "plasma", "whole serum",
                                    "serum)", "immobilised", "immobilized"))


def collect():
    recs = B.build(filter_compounds=False)[1]
    blocks = B.condition_index(recs)
    B.stamp_conditions(recs, blocks)
    B.stamp_basis(recs, B.basis_index(recs))

    rows = []
    for r in recs:
        if not is_serum_albumin(r["protein"]):
            continue
        if r["measure"] not in ("Kd", "Ka"):
            continue
        cond = blocks[B.block_key(B.condition_row(r))]
        row = {
            "pfas": r["pfas"], "pfas_canonical": N.canon(r["pfas"]),
            "protein": r["protein"], "species": r["species"],
            "source": r["source"], "method": r["_method"],
            "measure": r["measure"], "value_as_reported": r["value"],
            "units_as_reported": r["units"],
            "kd_uM": r["kd_uM"],
            "n_sites_used": r["n_used"], "n_basis": r["n_basis"],
            # The paper's own wording, which is fuller than the number far
            # more often than not. Of the rows that reach no n_sites_used,
            # many still SAY something: Han 2003's "6-9", Bischel 2010's
            # "1-3 primary sites", Chi 2018's "up to 8 at high mole
            # ratios". None of those is one number, so none converts -- but
            # a reader deciding what n to assume should see them.
            "n_sites_as_stated": r["_n"],
            "logK": "", "logK_if_n_assumed_1": "",
            "protein_conc_uM": cond["protein_conc_uM"],
            "protein_conc_as_stated": cond["protein_conc_as_stated"],
            "pfas_conc_uM": cond["pfas_conc_uM"],
            "pfas_conc_as_stated": cond["pfas_conc_as_stated"],
            "ligand_protein_ratio_max": cond["ligand_protein_ratio_max"],
            "location": r["location"]}

        if r["logK"] and r["logK_basis"].startswith(
                "from the paper's own Kd and n"):
            row["convertible"] = "yes"
            row["logK"] = r["logK"]
            row["why_not_as_reported"] = ""
        elif r["kd_uM"]:
            # A Kd and no site count the paper stands behind. n = 1 is the
            # one assumption that closes it, kept in its own column.
            row["convertible"] = "assume n = 1"
            row["logK_if_n_assumed_1"] = "%.3f" % (
                BRIDGE - math.log10(float(r["kd_uM"])))
            row["why_not_as_reported"] = short_reason(r["logK_basis"], row)
        else:
            row["convertible"] = "no"
            row["why_not_as_reported"] = short_reason(r["logK_basis"], row)
        rows.append(row)

    order = {"yes": 0, "assume n = 1": 1, "no": 2}
    rows.sort(key=lambda r: (order[r["convertible"]], r["source"],
                             r["pfas_canonical"]))
    return rows


def main():
    rows = collect()
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    counts = collections.Counter(r["convertible"] for r in rows)
    print(f"{len(rows)} serum-albumin affinity rows, "
          f"{len({r['source'] for r in rows})} sources.")
    print(f"  albumin bridge constant log10(1e9/M) = {BRIDGE:.4f} "
          f"per kg, M = {B.MW_ALBUMIN:.0f} g/mol")
    print()
    for verdict in ("yes", "assume n = 1", "no"):
        sub = [r for r in rows if r["convertible"] == verdict]
        print(f"{verdict.upper():14s} {counts[verdict]:4d} rows")
        if verdict == "yes":
            for s, n in collections.Counter(
                    r["source"] for r in sub).most_common():
                part = [r for r in sub if r["source"] == s]
                ns = sorted({r["n_sites_used"] for r in part})
                basis = sorted({r["n_basis"] for r in part})
                print(f"     {n:4d}  {s:16s} n = {', '.join(ns[:4])}"
                      f"{' ...' if len(ns) > 4 else ''}"
                      f"   ({'/'.join(basis)})")
        else:
            for why, n in collections.Counter(
                    r["why_not_as_reported"] for r in sub).most_common():
                srcs = sorted({r["source"] for r in sub
                               if r["why_not_as_reported"] == why})
                print(f"     {n:4d}  {why}")
                print(f"           {', '.join(srcs)[:92]}")
        print()

    usable = [r for r in rows if r["convertible"] != "no"]
    print("CONDITIONS, on the rows that reach a number "
          f"({len(usable)} of {len(rows)}):")
    for f in ("protein_conc_uM", "protein_conc_as_stated", "pfas_conc_uM",
              "pfas_conc_as_stated", "ligand_protein_ratio_max"):
        got = sum(1 for r in usable if str(r[f]).strip())
        print(f"  {f:26s} {got:4d} / {len(usable)}")
    print()
    print(f"wrote {os.path.abspath(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
