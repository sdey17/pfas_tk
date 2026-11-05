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
import re
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
LIT = os.path.join(HERE, "..", "figures", "albumin_chainlength_data.csv")
OUT = os.path.join(HERE, "..", "db", "docking_derived_logk.csv")
# A SECOND, plot-ready file, in the schema the per-compound figure
# already reads, so the figure can overlay these without the measured
# CSV ever carrying a predicted value. Kept separate on purpose: the
# standing rule in excluded_sources.py is that no computed value
# reaches a plotted artifact, and this file is that rule's one named
# exception rather than a hole in it.
OUT_PLOT = os.path.join(HERE, "..", "figures",
                        "albumin_docking_data.csv")
PLOT_FIELDS = ["pfas", "n_cf", "family", "species", "log_k",
               "category", "source", "log_ratio", "ratio_reported"]

# The same compound under two names. The database records each paper's
# own spelling; the plotted CSV uses the short form the rest of the
# literature uses, and a panel only lines up if they are made to
# agree. Checked one at a time against the structure, not guessed from
# the string: the chlorinated ether sulfonates are the F-53B
# components, named by perfluorinated-chain length in the n:2
# convention, and NaDONA is the sodium salt of the ADONA acid.
ALIAS = {
    "9Cl-PF3ONS": "6:2 Cl-PFESA",    # F-53B major, C8 chain
    "11Cl-PF3OUdS": "8:2 Cl-PFESA",  # F-53B minor, C10 chain
    "NaDONA": "ADONA",               # sodium salt / free acid
    "HFPO-DA": "HFPO-DA",            # GenX; same name both sides
}

R, T = 8.314, 298.15          # J/mol/K, K
KCAL = 4184.0                 # J per kcal
M_ALB = 66463.0
BRIDGE = math.log10(1e9 / M_ALB)

FIELDS = ["study", "year", "chemical", "protein", "species",
          "dG_kcal_per_mol", "K_M^-1", "Kd_uM", "log_K",
          "lit_median_log_K", "lit_n", "gap_vs_lit", "pdb", "note"]


def canon(study, year):
    """DB study+year -> the short label the figures use."""
    short = re.sub(r"\s+et al\.?", "", study).strip()
    return f"{short} {year}"


def family(head):
    h = (head or "").lower()
    if "carboxyl" in h:
        return "carboxylate"
    if "sulfonate" in h:
        return "sulfonate"
    return "other"


def main():
    rows = list(csv.DictReader(open(DB, encoding="utf-8")))
    dock = [r for r in rows
            if "kcal/mol" in (r["kd_units"] + r["ka_units"])
            and "docking" in r["method"].lower()]

    ncf, head = {}, {}
    for r in rows:
        c = r["chemical"]
        if r["n_perfluorinated_carbons"].strip().isdigit():
            ncf.setdefault(c, int(r["n_perfluorinated_carbons"]))
        if r["head_group"].strip():
            head.setdefault(c, r["head_group"])

    lit = collections.defaultdict(list)
    for r in csv.DictReader(open(LIT, encoding="utf-8")):
        if r["species"] == "human" and r["source"] != "Wu 2024":
            lit[r["pfas"]].append(float(r["log_k"]))

    out, plot, gaps = [], [], []
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
        key = ALIAS.get(key, key)
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

        n = ncf.get(r["chemical"])
        if n is None:
            continue                 # no nCF: no panel to sit in
        plot.append(dict(
            pfas=key, n_cf=n, family=family(head.get(r["chemical"])),
            species=r["protein_species"], log_k=f"{lk:.4f}",
            category="docking",
            # The source label carries the word on the figure's row,
            # so the row cannot be read as a laboratory by accident.
            source=f"{canon(r['study'], r['year'])} DOCKING",
            # No loading: a score is not an assay, so it has no
            # ligand:protein ratio. Empty here means the companion
            # ratio figure drops these rows on its own, with no
            # special case needed there.
            log_ratio="", ratio_reported=""))

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: -float(r["log_K"])))

    with open(OUT_PLOT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=PLOT_FIELDS)
        w.writeheader()
        w.writerows(sorted(plot, key=lambda r: (r["n_cf"],
                                                r["pfas"])))

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
