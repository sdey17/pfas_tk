#!/usr/bin/env python3
"""Argoul 2025: measured plasma unbound fractions, human and mouse, 15 PFAS.

    python3 -I scripts/add_argoul2025_plasma_fu.py

Writes db/plasma_unbound_fraction.csv.

    Argoul CML, Toutain PL, Berard S, Picard-Hagen N, Gayrard V, Lacroix MZ.
    "Structure-related differences in plasma protein binding of per- and
    polyfluoroalkyl substances in mice and humans." Environ Toxicol Pharmacol
    2025;120:104879. doi:10.1016/j.etap.2025.104879  PMID 41314320.

This is the primary source for the unbound fractions Argoul 2026 Table 1
quotes. Argoul 2026's own reference list does NOT contain it: that list gives
the same LC-MS/MS method paper twice, as 2024 and 2025, with identical titles
and the identical DOI 10.1016/j.jpba.2024.116596, and that paper is about
quantifying PFAS in mouse plasma rather than about binding. The citation in
the 2026 text resolves to the wrong reference, which is why this file records
the real one.

METHOD. Equilibrium dialysis on the Dianorm system, chosen after an
equilibrium-time study over a 96-well dialyser (up to 120 h to equilibrium
for long-chain PFAS) and over RED. Human plasma pooled from 77 donors;
C57BL/6N mouse plasma. Two independent experiments, reported separately in
Table 2 and stored separately here -- no averaging where the paper gives a
range. PFDS is absent because equilibrium was not reached.

WHY IT MATTERS TO THIS PROJECT. OEHHA PHG Table A6.4, and Han 2012 before it,
computes net renal tubular reabsorption as fu x GFR - CLR with fu ASSUMED to
be 0.02 in every species. These are the measured values, and they are not
0.02 and not species-constant.

THE CONVERGENCE IS COMPOUND-SPECIFIC, and PFOA is not one of the converged
ones. The paper states that mouse values "were at least 2 times higher than
in humans for all PFAS excepted PFOS and PFDA":

    PFOS   1.4x apart     PFDA   1.2x apart      <- converged
    PFNA   2.4x
    PFOA   6.7x                                  <- NOT converged
    PFBS    16x           PFHxS   17x

So a rat fu for PFOA cannot be inferred from the human value, nor from the
mouse value, and the roadmap's earlier argument that long-chain fu converges
across species was too broad.

Idempotent.
"""
import csv, io, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                   "plasma_unbound_fraction.csv")
SRC = ("Argoul 2025, Environ Toxicol Pharmacol 120:104879, "
       "doi:10.1016/j.etap.2025.104879, Table 2")
METHOD = ("equilibrium dialysis, Dianorm system; human plasma pooled from 77 "
          "donors, C57BL/6N mouse plasma")
FIELDS = ["chemical", "mw_g_mol", "species", "experiment", "fu_percent",
          "fu", "range_low_percent", "range_high_percent", "n_replicates",
          "method", "source"]

#  chemical, MW, [(mouse exp1, lo, hi), (mouse exp2, lo, hi),
#                 (human exp1, lo, hi), (human exp2, lo, hi)]
T2 = [
    ("PFBS", 300.10, (29, 27, 32), (25, 24, 28), (1.7, 1.6, 1.7), (1.7, 1.5, 1.8)),
    ("PFPeS", 350.11, None, (6.5, 6.1, 7.2), None, (0.26, 0.24, 0.28)),
    ("PFHxS", 400.12, (1.3, 1.2, 1.5), (1.1, 0.87, 1.2), (0.08, 0.06, 0.10), (0.06, 0.05, 0.06)),
    ("PFHpS", 450.12, None, (0.40, 0.35, 0.43), None, (0.02, 0.01, 0.02)),
    ("PFOS", 500.13, (0.25, 0.23, 0.27), (0.21, 0.14, 0.31), (0.16, 0.10, 0.28), (0.16, 0.14, 0.19)),
    ("PFNS", 550.14, None, (1.57, 0.71, 2.5), None, (1.1, 0.8, 1.4)),
    ("PFPeA", 264.05, None, (50, 49, 51), None, (8.5, 7.9, 9)),
    ("PFHxA", 314.05, (25, 24, 26), (17, 17, 18), (2, 2, 2.1), (2, 1.9, 2.1)),
    ("PFHpA", 364.06, None, (3.9, 3.5, 4.1), None, (0.29, 0.27, 0.30)),
    ("PFOA", 414.07, (0.87, 0.82, 0.96), (1.2, 1.2, 1.3), (0.17, 0.13, 0.19), (0.13, 0.12, 0.14)),
    ("PFNA", 464.08, (0.35, 0.34, 0.38), (0.36, 0.34, 0.38), (0.17, 0.10, 0.24), (0.13, 0.12, 0.15)),
    ("PFDA", 514.08, (0.42, 0.36, 0.51), (0.36, 0.28, 0.51), (0.38, 0.16, 0.75), (0.25, 0.20, 0.30)),
    ("PFUdA", 564.09, None, (1.21, 0.94, 1.4), None, (1.7, 1.0, 2.3)),
    ("PFO2OA", 346.05, (10, 9, 11), None, (0.47, 0.42, 0.51), None),
    ("HFPO-DA (GenX)", 330.05, (26, 24, 27), None, (4.2, 3.9, 4.5), None),
]
NREP = {("mouse", 1): "3", ("mouse", 2): "3",
        ("human", 1): "3", ("human", 2): "2"}


def main():
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader()
    n = 0
    for chem, mw, m1, m2, h1, h2 in T2:
        for sp, exp, v in (("mouse", 1, m1), ("mouse", 2, m2),
                           ("human", 1, h1), ("human", 2, h2)):
            if v is None:
                continue
            fu, lo, hi = v
            w.writerow({
                "chemical": chem, "mw_g_mol": mw, "species": sp,
                "experiment": exp, "fu_percent": fu, "fu": round(fu / 100, 6),
                "range_low_percent": lo, "range_high_percent": hi,
                "n_replicates": NREP[(sp, exp)], "method": METHOD,
                "source": SRC,
            })
            n += 1
    text = s.getvalue()
    if os.path.exists(OUT) and open(OUT, encoding="utf-8").read() == text:
        print("no change needed (Argoul 2025 Table 2 already ingested)")
        return 0
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"Argoul 2025 Table 2: {n} rows -> db/plasma_unbound_fraction.csv")
    print(f"   {len(T2)} compounds, 2 species, up to 2 experiments each")
    print("   PFOA  mouse 0.87-1.2%   human 0.13-0.17%   6.7x apart")
    print("   PFOS  mouse 0.21-0.25%  human 0.16%        1.4x apart")
    print("   OEHHA PHG Table A6.4 assumes 2% for every species and compound")
    return 0


if __name__ == "__main__":
    sys.exit(main())
