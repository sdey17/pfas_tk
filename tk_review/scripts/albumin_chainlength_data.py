#!/usr/bin/env python3
"""Build the CSV the chain-length plot reads.

Same split as albumin_loading_data.py: this half knows where the
project keeps things, so the plotting half can live in a notebook with
no paths in it.

One row per value: pfas, n_cf, family, species, log_k, category,
source. `category` is how the number was arrived at, which is what the
plot colours by; `species` and `family` are what it facets by, because
a chain-length trend only means anything within one protein and one
head group.

Rows with no perfluorinated-carbon count are dropped. Those are the
aggregate rows that name several compounds at once, plus the
fluorotelomers and sulfonamides the database leaves uncounted -- 17 of
293.

Run:  python3 tk_review/scripts/albumin_chainlength_data.py
"""
import collections
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
SRC = os.path.join(HERE, "..", "db", "albumin_logK_shareable.xlsx")
# The same extract with the assay conditions joined on, which is where
# the ligand:protein ratio lives.
COND = os.path.join(HERE, "..", "db", "albumin_logK_with_conc.xlsx")
OUT = os.path.join(HERE, "..", "figures", "albumin_chainlength_data.csv")
FIELDS = ["pfas", "n_cf", "family", "species", "log_k", "category",
          "source", "log_ratio", "ratio_reported"]

# Which sources STATE a ligand:protein ratio rather than having one
# computed here as PFAS concentration over protein concentration.
# Allendorf gives it per compound in SI Table SI-2.3 and Fischer gives
# it as a design bound; everyone else states two concentrations and no
# ratio, so their ratio is this project's division and is marked as
# such. This matters because the ratio is the whole x axis of the
# loading plots.
RATIO_REPORTED = {"Allendorf 2019", "Fischer 2024"}


def lookup():
    ncf, head = {}, {}
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        c = r["chemical"]
        if r["n_perfluorinated_carbons"].strip().isdigit():
            ncf.setdefault(c, int(r["n_perfluorinated_carbons"]))
        if r["head_group"].strip():
            head.setdefault(c, r["head_group"])
    return ncf, head


def family(head):
    h = (head or "").lower()
    if "carboxyl" in h:
        return "carboxylate"
    if "sulfonate" in h:
        return "sulfonate"
    return "other"


def ratios():
    """(pfas, species, source) -> log10 ligand:protein ratio."""
    from openpyxl import load_workbook
    import math
    ws = load_workbook(COND)["Coefficients"]
    ix = {h: i for i, h in enumerate(c.value for c in ws[1])}
    out = {}
    for r in ws.iter_rows(min_row=2, values_only=True):
        try:
            v = float(str(r[ix["ligand:protein ratio (max)"]]).strip())
        except (TypeError, ValueError):
            continue
        if v <= 0:
            continue
        key = (str(r[ix["PFAS"]]), str(r[ix["Species"]]),
               str(r[ix["Source"]]).split(" Cw")[0].strip())
        out.setdefault(key, round(math.log10(v), 4))
    return out


def main():
    from openpyxl import load_workbook
    ncf, head = lookup()
    rat = ratios()
    ws = load_workbook(SRC)["Coefficients"]
    ix = {h: i for i, h in enumerate(c.value for c in ws[1])}

    out, dropped = [], 0
    for r in ws.iter_rows(min_row=2, values_only=True):
        name = str(r[ix["PFAS"]])
        if name not in ncf:
            dropped += 1
            continue
        route = str(r[ix["Route"]])
        src = str(r[ix["Source"]])
        cat = ("qin" if src.startswith("Qin 2024") else
               "measured" if route.startswith("measured") else
               "quenching" if "quenching" in route else "converted")
        source = src.split(" Cw")[0].strip()
        species = str(r[ix["Species"]])
        lr = rat.get((name, species, source))
        out.append({
            "pfas": name, "n_cf": ncf[name],
            "family": family(head.get(name)),
            "species": species,
            "log_k": float(r[ix["log K (L water / kg protein)"]]),
            "category": cat, "source": source,
            "log_ratio": "" if lr is None else lr,
            "ratio_reported": ("" if lr is None else
                               "yes" if source in RATIO_REPORTED
                               else "no")})

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: (r["family"], r["species"],
                                               r["n_cf"], r["log_k"])))
    print(f"wrote {OUT}")
    print(f"  {len(out)} rows, {dropped} dropped for having no nCF")
    got = [r for r in out if r["log_ratio"] != ""]
    rep = sum(1 for r in got if r["ratio_reported"] == "yes")
    print(f"  ratio    {len(got)} rows have one; {rep} stated by the "
          f"paper, {len(got) - rep} computed here")
    for k in ("family", "species", "category"):
        print(f"  {k:<9}{dict(collections.Counter(r[k] for r in out))}")


if __name__ == "__main__":
    main()
