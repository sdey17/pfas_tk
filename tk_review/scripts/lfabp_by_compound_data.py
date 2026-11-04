#!/usr/bin/env python3
"""Build the CSV the L-FABP per-compound plot reads.

Same split as albumin_chainlength_data.py: this half knows where the
project keeps things, so the plotting half can live in a notebook with
no paths in it.

ONLY LIVER FABP.  The filter is on "liver fatty acid", not on "fabp",
because the database also holds INTESTINAL FABP -- two Khazaee 2021
rows.  A substring match on "fabp" sweeps I-FABP into an L-FABP
comparison, which is a different protein; that mistake was made once in
this project already and produced a PFNA "inversion" that did not
exist.

EVERY L-FABP NUMBER HERE IS CONVERTED FROM A Kd.  No paper reports a
measured protein-water partition coefficient for any protein except
albumin, so there is nothing to colour by provenance the way the
albumin plots do.  What varies between these laboratories, and what
the plot colours by, is the ASSAY:

    separation   equilibrium dialysis -- free and bound physically
                 separated, then both measured
    direct       ITC or MST -- binding observed directly, with no
                 probe and no separation step
    displacement a probe (ANS, DAUDA, a fatty acid) is competed off
                 the protein and the PFAS constant is inferred from
                 how much probe signal is lost

Those are not the same quantity.  A displacement constant depends on
the probe's own affinity and on the model used to back it out; a
direct Kd does not.  Three quarters of the L-FABP literature is
displacement, which is why it gets its own colour rather than a
footnote.

MUTANTS GET THEIR OWN ROWS.  Sheng 2016 measured wild-type hL-FABP and
four site-directed point variants (S39G, M74G, N111D, R122G).  Those
are different proteins and pooling them into one "Sheng 2016" row
would hide the paper's actual result: by ITC, N111D raises the PFOA Kd
from 6.4 to 284 uM, a 44-fold weakening and 1.65 log, while R122G
barely moves it.

NO-BINDING ROWS GO TO THE SIDECAR, not into the drop count.  "No
binding detected" is a result, not a missing value, and on these
mutants it is the result that disagrees with the other assay: the
fluorescence-displacement assay reports NB for both N111D and R122G on
PFOA and PFNA, while ITC on the same protein gives R122G a Kd of
8.4 uM.  It cannot be plotted on a log axis, so it is written out
instead of vanishing into a count.

Run:  python3 tk_review/scripts/lfabp_by_compound_data.py
"""
import collections
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
SRC = os.path.join(HERE, "..", "db", "albumin_lfabp_reported_conc.xlsx")
OUT = os.path.join(HERE, "..", "figures", "lfabp_by_compound_data.csv")
OUT_EXCLUDED = os.path.join(HERE, "..", "figures",
                            "lfabp_by_compound_excluded.csv")
FIELDS = ["pfas", "n_cf", "family", "species", "protein_label", "log_k",
          "method", "source", "variant"]

# Held out, with the reason carried into the sidecar CSV. Nothing is
# excluded for disagreeing with the rest; each is excluded because the
# number is not an independent measurement of a dissociation constant.
EXCLUDED = {
    "no binding detected": (
        "reported as NB or ND -- a result, not a missing value, and "
        "not plottable on a log axis. Worth reading against the same "
        "protein's other assay: fluorescence displacement reports NB "
        "for Sheng 2016's N111D and R122G on PFOA and PFNA, where ITC "
        "gives R122G a Kd of 8.4 uM"),
    "Han 2021": (
        "its hL-FABP constants are explicitly quoted from other "
        "papers, not measured here -- the notes on those rows say so, "
        "and the PFBS value of 1034 uM is Zhang 2013's. Plotting it as "
        "a laboratory row would double-count whatever source it came "
        "from. Its own measurements are thermal-shift dTm, which is "
        "not a dissociation constant and does not belong on a log K "
        "axis"),
}

METHOD = [
    (re.compile(r"equilibrium dialysis", re.I), "separation"),
    (re.compile(r"\bITC\b|\bMST\b", re.I), "direct"),
    (re.compile(r"displacement", re.I), "displacement"),
]

VARIANT = re.compile(r"ENGINEERED POINT VARIANT (\S+)", re.I)


def lookup():
    """chemical -> (nCF, head group), from the database."""
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


def method_of(route):
    for pat, name in METHOD:
        if pat.search(route or ""):
            return name
    return None


def main():
    from openpyxl import load_workbook
    ncf, head = lookup()
    ws = load_workbook(SRC)["Coefficients"]
    ix = {h: i for i, h in enumerate(c.value for c in ws[1])}

    out, held = [], []
    no_logk = no_ncf = not_lfabp = no_method = 0
    for r in ws.iter_rows(min_row=2, values_only=True):
        prot = str(r[ix["Protein"]] or "")
        # "liver fatty acid", NOT "fabp" -- see the module docstring.
        if "liver fatty acid" not in prot.lower():
            not_lfabp += 1
            continue
        name = str(r[ix["PFAS"]])
        note = str(r[ix["Note"]] or "")
        nb = "no binding detected" in note.lower()
        if r[ix["log K (L water / kg protein)"]] is None and not nb:
            no_logk += 1
            continue
        if name not in ncf:
            no_ncf += 1
            continue
        # The holdout is decided BEFORE the assay filter. Han 2021's
        # route is thermal shift, which the filter would drop silently,
        # and a silently dropped row looks the same as one that was
        # never there. It belongs in the sidecar with its reason.
        source = str(r[ix["Source"]]).split(" Cw")[0].strip()
        hold = "no binding detected" if nb else source
        meth = method_of(str(r[ix["Route"]] or ""))
        if meth is None and hold not in EXCLUDED:
            no_method += 1
            continue
        m = VARIANT.search(prot)
        variant = m.group(1) if m else "native"
        row = {
            "pfas": name, "n_cf": ncf[name],
            "family": family(head.get(name)),
            "species": str(r[ix["Species"]]),
            "protein_label": ("L-FABP" if variant == "native"
                              else f"L-FABP {variant}"),
            "log_k": ("" if nb else
                      round(float(r[ix["log K (L water / kg protein)"]]), 4)),
            "method": meth or "", "source": source,
            "variant": variant}
        if hold in EXCLUDED:
            row["_why"] = EXCLUDED[hold]
            held.append(row)
        else:
            out.append(row)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: (r["species"], r["source"],
                                               r["n_cf"], r["log_k"])))
    if held:
        with open(OUT_EXCLUDED, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS + ["why"])
            w.writeheader()
            for r in sorted(held, key=lambda r: (r["source"], r["n_cf"])):
                why = r.pop("_why")
                w.writerow({**r, "why": why})

    print(f"wrote {OUT}")
    print(f"  {len(out)} L-FABP rows")
    print(f"  skipped  {not_lfabp:4d} not liver FABP (albumin, I-FABP, "
          f"reviews)")
    print(f"           {no_logk:4d} no log K (rank orders, dTm, "
          f"aggregates)")
    print(f"           {no_ncf:4d} no nCF (aggregate rows, reference "
          f"ligands)")
    print(f"           {no_method:4d} no recognised assay")
    nb_n = sum(1 for r in held if r["log_k"] == "")
    if nb_n:
        print(f"  held out {nb_n} 'no binding detected' rows "
              f"-> {OUT_EXCLUDED}")
    for src in sorted({r["source"] for r in held if r["log_k"] != ""}):
        n = sum(1 for r in held
                if r["source"] == src and r["log_k"] != "")
        print(f"  held out {n} rows from {src} -> {OUT_EXCLUDED}")
    for k in ("species", "method", "variant"):
        print(f"  {k:<9}{dict(collections.Counter(r[k] for r in out))}")
    print("  sources  " + ", ".join(
        f"{s} ({n})" for s, n in
        collections.Counter(r["source"] for r in out).most_common()))


if __name__ == "__main__":
    main()
