#!/usr/bin/env python3
"""Docking scores on the measurement axis: the conversion, and its limits.

CAN A DOCKING SCORE BE CONVERTED?  Mechanically, yes, and the
arithmetic is the same definition used everywhere else here:

    dG = -RT ln K     ->     K = exp(-dG / RT)
    Kd(uM) = 1e6 / K
    log K = log10(1e9/M) - log10(Kd(uM))

A docking score IS a predicted standard free energy of binding, in
kcal/mol, so it carries the right units for this and needs only the
4.184 kJ/kcal factor.  Nothing about the conversion is unsound.

WHAT THE CONVERSION DOES NOT FIX, and the reason these stay a
separate series rather than joining the measured ones:

  A scoring function is not an assay. Agreement is checked below
  against the measured human-albumin literature, and over the six
  compounds where both exist the median gap is -0.03 log -- better
  than most pairs of laboratories manage. That is a real result and
  it is why plotting them is worth doing.

  But the agreement is in the middle of the range only. Docking
  here spans 1.54 log across ten compounds while the measured
  literature spans five, so the scores are compressed: they cannot
  reproduce either a very weak or a very strong binder.

  And they are monotonic in chain length BY CONSTRUCTION. Wu 2024's
  own docking rises PFBA -> PFBS -> PFHxS -> PFOA -> PFOS without a
  turn, and Li 2021's docking columns do the same where that paper's
  own measurements turn over above PFHpA. The compounds Wu docked
  stop around C8, which is where the turnover begins, so this
  comparison does not test the regime where docking is known to
  fail. Agreement here is not a licence at the ends.

So: converted, plotted, visibly distinct, and never admitted to a
pooled statistic or to the coefficient-selection rule. The standing
rule in scripts/excluded_sources.py keeps computed values out of
every reading copy; the per-compound FIGURE is the one place they
now appear, marked as predictions, because seeing them beside the
measurements is the point.

Run:  python3 tk_review/scripts/docking_to_logk.py
Writes db/docking_derived_logk.csv
"""
import collections
import csv
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
LIT = os.path.join(HERE, "..", "figures", "albumin_chainlength_data.csv")
OUT = os.path.join(HERE, "..", "db", "docking_derived_logk.csv")

R, T = 8.314, 298.15          # J/mol/K, K
KCAL = 4184.0                 # J per kcal
M_ALB = 66463.0
BRIDGE = math.log10(1e9 / M_ALB)

FIELDS = ["study", "year", "chemical", "protein", "species",
          "dG_kcal_per_mol", "K_M^-1", "Kd_uM", "log_K",
          "lit_median_log_K", "lit_n", "gap_vs_lit", "pdb", "note"]


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    dock = [r for r in rows
            if "kcal/mol" in (r["kd_units"] + r["ka_units"])
            and "docking" in r["method"].lower()]

    lit = collections.defaultdict(list)
    for r in csv.DictReader(open(LIT, encoding="utf-8")):
        if r["species"] == "human" and r["source"] != "Wu 2024":
            lit[r["pfas"]].append(float(r["log_k"]))

    out, gaps = [], []
    for r in dock:
        val = (r["kd"] or r["ka"]).strip()
        try:
            dG = float(val)
        except ValueError:
            continue
        K = math.exp(-dG * KCAL / (R * T))
        kd = 1e6 / K
        lk = BRIDGE - math.log10(kd)
        key = r["chemical"].split(" (")[0]
        L = lit.get(key) or lit.get(r["chemical"]) or []
        med = st.median(L) if L else None
        if med is not None:
            gaps.append(lk - med)
        pdb = ""
        for tok in ("4E99", "1E7G", "1AO6", "7AAI"):
            if tok in r["method"] + r["source_table"]:
                pdb = tok
        out.append(dict(
            study=r["study"], year=r["year"], chemical=r["chemical"],
            protein=r["protein"][:40], species=r["protein_species"],
            dG_kcal_per_mol=f"{dG:.2f}", **{"K_M^-1": f"{K:.4g}"},
            Kd_uM=f"{kd:.3f}", log_K=f"{lk:.3f}",
            lit_median_log_K="" if med is None else f"{med:.2f}",
            lit_n=len(L),
            gap_vs_lit="" if med is None else f"{lk-med:+.2f}",
            pdb=pdb,
            note="PREDICTED from a docking score, not measured"))

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: -float(r["log_K"])))

    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print(f"  {len(out)} docking scores converted\n")
    print(f"  {'PFAS':17}{'dG kcal':>9}{'log K':>8}{'lit med':>9}"
          f"{'n':>4}{'gap':>8}")
    for r in out:
        print(f"  {r['chemical'][:16]:17}{r['dG_kcal_per_mol']:>9}"
              f"{r['log_K']:>8}{r['lit_median_log_K'] or '-':>9}"
              f"{r['lit_n']:>4}{r['gap_vs_lit'] or '-':>8}")
    if gaps:
        lk = [float(r["log_K"]) for r in out]
        print(f"\n  median docking - literature {st.median(gaps):+.2f} "
              f"log over {len(gaps)} compounds")
        print(f"  docking spans {max(lk)-min(lk):.2f} log across "
              f"{len(lk)} compounds;")
        allv = [v for vs in lit.values() for v in vs]
        print(f"  the measured human literature spans "
              f"{max(allv)-min(allv):.2f} log. The scores are "
              f"compressed.")


if __name__ == "__main__":
    main()
