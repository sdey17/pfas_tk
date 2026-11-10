#!/usr/bin/env python3
"""Moro 2022 reduced to ONE row per compound, for the working table.

The database carries all 12 Moro rows -- 4 compounds x site class x
temperature -- because they are all real measurements and the 298/310 K
pairs are the only within-study temperature test in this project for
ITC.  A plot or a summary table wants one row per compound.  This
applies the selection rule and prints that table.

THE RULE, in two steps, and why it is in this order.

STEP 1, the SITE CLASS: take the class that dominates at PHYSIOLOGICAL
loading, not at the loading the assay reached.  Serum PFAS is a few nM
against ~600 uM albumin, so ligand:protein is ~1e-5 -- the trace
regime, the far left of the titration.  For two independent sets of
sites the occupancy from class i is

    nu_i = n_i [L]/(Kd_i + [L])  ->  n_i [L]/Kd_i   as [L] -> 0

so the classes split in proportion to n_i/Kd_i, independent of
concentration -- and n/Kd is exactly what log K encodes, so the class
with the highest log K is the one that carries the ligand.  For PFOA
that is the high-affinity set: n/Kd of 2.225 against 0.152, a 14.7x
split, so the tight set carries 93.6% of bound PFOA despite holding 5x
LESS capacity.  Adding both sets moves log K by +0.03.

For a sequential fit the comparison is not even close: the site-2
terms in

    nu = (K1[L] + 2 K1 K2 [L]^2) / (1 + K1[L] + K1 K2 [L]^2)

are SECOND order in [L] and vanish faster than site 1, so at trace
site 2 contributes nothing at all -- it only appears above ~10 uM free,
1e6 times physiological.  Moro says as much: HFPO-DA's second event
"could instead be plausibly attributed to unspecific hydrophobic
interactions occurring at high HPFO-DA concentrations", and its
enthalpy is POSITIVE (+28.2 kJ/mol), entropically driven, the
signature of desolvation rather than of a pocket.  n is 1 for both
steps by the model, so highest n/Kd picks site 1 either way.

STEP 2, the TEMPERATURE: prefer 310 K, which is body temperature and
the condition a TK model is about.  Fall back to 298 K only where the
310 K cell is flagged unusable -- which is PFOA's high-affinity set,
whose printed 39.3 uM contradicts its own printed dG by 12 kJ/mol and
implies 0.38 uM instead.  PFOA's tight site is genuinely temperature-
insensitive, so 298 K stands in for it honestly; that is NOT true of
the other compounds, where the two temperatures differ 2-3x in
OPPOSITE directions (C6O4 2.3x weaker at 310 K, HFPO-DA 1.9x tighter)
and the printed error bars do not overlap, despite the paper's claim
of "no significant variations".

The class is chosen before the temperature because doing it the other
way round picks PFOA's LOW-affinity set at 310 K -- the high-affinity
row there being unusable -- which is the wrong site by 1.2 log.

Nothing here is hard-coded per compound: the rule reads the site
classes, the n/Kd values and the SUSPECT flag out of the database, so
a correction to any of them changes this table.

Run:  python3 tk_review/scripts/moro2022_representative.py [-o PATH]
"""
import argparse
import collections
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "moro2022_representative.tsv")

BRIDGE = 4.1774                  # log10(1e9 / 66463), albumin
PREFER_K = 310                   # body temperature

# The concentration columns must describe ONE MOMENT, or they do not
# divide to the ratio and the table reads as broken.  There are three
# states in an ITC run and only two of them are the cell:
#
#   start of titration   PFAS   0 uM,  protein 110 uM
#   the syringe          PFAS 4000 uM, protein none     <- NOT the cell
#   end of titration     PFAS  ~680 uM, protein  91 uM
#
# The row used to pair the syringe with the start, which is no moment
# at all -- hence 4000/110 = 36 against a real ratio near 7.  It now
# reports the END of the titration: the highest PFAS the protein
# actually saw, beside the protein concentration at that same point,
# so [PFAS]/[protein] reproduces the measured ratio exactly.
#
# Protein falls from 110 to 91 because an ITC cell is sealed and full,
# so each injection displaces its own volume and slowly washes albumin
# out: 110 * exp(-V_inj/V_cell).  The 200 uL cell this uses is no
# longer an assumption -- it predicts ratios of 6.9-7.6 and Figure 1
# measures 6.7-7.6 for three of the four compounds.
V_CELL, V_INJ, P_PREP = 200.0, 0.4 + 25 * 1.5, 110.0
ORDER = ["PFOA", "C6O4", "PFHxA", "HFPO-DA (GenX)"]

HDR = ["PFAS", "Species", "Kd (uM)", "n used",
       "log K (L water / kg protein)", "Route", "Source",
       "PFAS conc (uM)", "Protein conc (uM)", "ligand:protein ratio",
       "f_bound (%)"]


def rows_for_study():
    return [r for r in csv.DictReader(open(DB, encoding="utf-8"))
            if r["study"] == "Moro et al." and r["year"] == "2022"]


def info(r):
    kd = float(r["kd"])
    n = float(r["n_binding_sites"])
    return dict(
        chem=r["chemical"],
        site=r["binding_site_identified"].split(";")[0].strip(),
        T=int(r["temperature_c"].split()[0]) + 273,
        kd=kd, n=n,
        trace=n / kd,                           # what sets the split
        logK=BRIDGE + math.log10(n) - math.log10(kd),
        suspect="SUSPECT VALUE" in r["assay_volume_note"],
        model=("sequential" if "sequential" in r["binding_site_identified"]
               else "two sets of sites" if "affinity set" in
               r["binding_site_identified"] else "one set of sites"),
        prot=r["protein_conc"], lig=r["ligand_conc_range"],
        ratio=r["ligand_over_protein_max"],
        fb=r["unbound_fraction"].strip() or "not reported",
        species=r["protein_species"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="", help="output path; '-' stdout")
    args = ap.parse_args()

    recs = [info(r) for r in rows_for_study()]
    if not recs:
        print("  no Moro 2022 rows", file=sys.stderr)
        return 1

    by_chem = collections.defaultdict(list)
    for x in recs:
        by_chem[x["chem"]].append(x)

    out, log = [], []
    for chem in sorted(by_chem, key=lambda c: ORDER.index(c)
                       if c in ORDER else 99):
        group = by_chem[chem]
        classes = collections.defaultdict(list)
        for x in group:
            classes[x["site"]].append(x)

        # STEP 1: the class that carries the ligand at trace. Score each
        # class by its best n/Kd among rows that are usable at all.
        def score(rs):
            ok = [x["trace"] for x in rs if not x["suspect"]]
            return max(ok) if ok else max(x["trace"] for x in rs)
        site = max(classes, key=lambda s: score(classes[s]))

        # STEP 2: body temperature unless that cell is unusable.
        cand = classes[site]
        pick = next((x for x in cand
                     if x["T"] == PREFER_K and not x["suspect"]), None)
        fell_back = pick is None
        if pick is None:
            pick = next((x for x in cand if not x["suspect"]), None)
        if pick is None:
            log.append(f"  {chem}: every row flagged, SKIPPED")
            continue

        other = sorted(s for s in classes if s != site)
        # Keep Route to one short readable clause: the model names the
        # site already, so do not print both.
        if pick["model"] == "sequential":
            what = f"sequential site 1 of {len(classes)}"
        elif pick["model"] == "two sets of sites":
            what = f"{site} of {len(classes)} sets"
        else:
            what = "one set of sites"
        route = (f"ITC, {what}, {pick['T'] - 273} C; concentrations are "
                 f"end of titration in the cell, from {P_PREP:.0f} uM dhSA "
                 f"prepared and 4 mM in the syringe")
        if fell_back:
            route += "; 37 C cell unusable, site is T-insensitive"
        # One moment: the end of the titration.
        p_end = P_PREP * math.exp(-V_INJ / V_CELL)
        pfas_end = float(pick["ratio"]) * p_end
        out.append([chem, pick["species"], f"{pick['kd']:g}",
                    f"{pick['n']:g}", f"{pick['logK']:.3f}", route,
                    "Moro 2022", f"{pfas_end:.0f}", f"{p_end:.0f}",
                    pick["ratio"], pick["fb"]])
        log.append(
            f"  {chem:16} kept {site} at {pick['T']} K"
            + (f"; dropped {', '.join(other)}" if other else "")
            + (" (fell back from 310 K)" if fell_back else ""))

    print(f"  {len(recs)} database rows -> {len(out)} representative rows",
          file=sys.stderr)
    for line in log:
        print(line, file=sys.stderr)

    body = "\n".join(["\t".join(HDR)] + ["\t".join(r) for r in out]) + "\n"
    if args.out == "-":
        sys.stdout.write(body)
        return 0
    path = args.out or OUT
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    print(f"  wrote {os.path.relpath(path, os.path.join(HERE, '..'))}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
