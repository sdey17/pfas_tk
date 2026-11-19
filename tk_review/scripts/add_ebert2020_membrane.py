#!/usr/bin/env python3
"""Ebert 2020: membrane/water partitioning, and what the paper does NOT show.

Ebert A, Allendorf F, Berger U, Goss KU, Ulrich N. "Membrane/water
partitioning and permeabilities of perfluoroalkyl acids and four of their
alternatives and the effects on toxicokinetic behavior."
Environ Sci Technol 2020;54(8):5051-5061. doi:10.1021/acs.est.0c00175
PMID 32212724.

The roadmap listed this as request 2 and framed it as potentially
strand-reframing: "argues measured membrane permeabilities explain cellular
uptake with no need to postulate active uptake. If right, it reframes the
whole strand." The reviewer supplied the full text. IT DOES NOT REFRAME THE
STRAND, and the framing came from reading the abstract rather than the body.

WHAT THE ABSTRACT CLAIMS. "Experimental Pion values were high enough to
explain observed cellular uptake by passive diffusion with no need to
postulate the existence of active uptake processes."

WHAT THE BODY ACTUALLY DOES. It compares the authors' bilayer permeabilities
against effective permeabilities back-calculated from the CONTROL cells of
transporter studies -- the empty-vector or non-expressing cells, which is
where those papers put their own passive baseline. The paper says so
plainly: the cited data "reflected uptake rates by passive diffusion". So
the quantity explained is the passive component that the transporter papers
had already subtracted off, not the transporter-mediated component they
reported. The net uptake in expressing-minus-control cells is untouched by
this analysis.

The body is also much more hedged than the abstract:

  * PFSAs: effective anionic permeabilities match the literature values
    within a factor of 3, and anionic permeation dominates.
  * PFCAs: "the comparison to the literature seems less distinctive" --
    anionic and neutral effective permeabilities are the same order as the
    literature values and either could contribute.
  * "a definite conclusion on the species dominance in the transmembrane
    transport is not possible."
  * The comparison needs the artificial-bilayer permeabilities shifted by
    1.8 log units, a factor of 63, "to compensate for possible differences
    between artificial bilayers and biomembranes".

So the correct reading is that passive permeation is quantitatively
sufficient to account for the uptake seen in control cells, which is a real
and useful result, and that it says nothing about whether the saturable
component on top of that is transporter-mediated. It does not conflict with
any Km in this database.

The paper does supply an independent check on one of ours: it notes that in
Yang et al. PFOA uptake falls to 77% of its pH 6 rate at pH 7 and 67% at
pH 8, which is the same pH dependence behind the two Yang 2010 OAT4 Km
values already in the record (172.3 uM at pH 6, 310.3 uM at pH 7.4).

WHAT IS STORED HERE. Only values the running text states outright. The
per-compound log Kmem/w and log Pion tables are in SI Tables S3.1.1 and
S3.2.1, which were not supplied; Figures 1 and 2 plot them without printing
numbers, and digitising a plotted point to two decimals would invent
precision. The SI is therefore a live request. What the text states:

    log Kmem/w  HFPO-DA 2.4   DONA 3.0   9Cl-PF3ONS 5.1   PFECHS 4.9
                PFOS 4.5      full range 2.3-5.1 across all compounds
                PFBA below the measurable limit (< 1.7)

and two relationships that are more useful than any single value:

    +0.63 log Kmem/w per perfluorinated carbon (global fit; 0.63 for PFCAs
         and 0.61 for PFSAs fitted separately)
    PFSAs sorb about 1.2 log units more strongly than PFCAs with the same
         number of perfluorinated carbons

That slope is directly comparable to the per-carbon slopes this project
already computes for albumin and L-FABP binding, and the 1.2-log head-group
offset is the membrane counterpart of the head-group effect in the binding
sheets.

Idempotent. Writes db/membrane_partitioning.csv; touches no other file.
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "membrane_partitioning.csv")
DOI = "10.1021/acs.est.0c00175"
STUDY = "Ebert A, Allendorf F, Berger U, Goss KU, Ulrich N"
SYSTEM = ("POPC liposomes, equilibrium dialysis, 12-14 kDa cellulose "
          "membrane; sorbed:phospholipid molar ratio held below 0.08 to stay "
          "in the linear range of the isotherm")
QUAL = ("stated in the running text; the per-compound table is SI Table "
        "S3.1.1, not supplied. Experimental uncertainty about +/- 0.2 log")

FIELDS = ["study", "year", "pmid", "doi", "chemical",
          "n_perfluorinated_carbons", "head_group", "quantity", "value",
          "value_units", "method", "note"]

ROWS = [
    ("HFPO-DA (GenX)", "5", "ether carboxylate", "log Kmem/w", "2.4",
     "log (L water / kg membrane)",
     QUAL + ". Comparable to PFHxA; the paper notes neither the ether oxygen "
     "nor the side chain measurably changes membrane sorption"),
    ("DONA", "", "ether carboxylate", "log Kmem/w", "3.0",
     "log (L water / kg membrane)",
     QUAL + ". Comparable to PFHpA"),
    ("PFOS", "8", "sulfonate", "log Kmem/w", "4.5",
     "log (L water / kg membrane)", QUAL),
    ("9Cl-PF3ONS", "8", "ether sulfonate", "log Kmem/w", "5.1",
     "log (L water / kg membrane)",
     QUAL + ". Eight perfluorinated carbons, like PFOS"),
    ("PFECHS", "8", "cyclic sulfonate", "log Kmem/w", "4.9",
     "log (L water / kg membrane)",
     QUAL + ". Eight perfluorinated carbons, like PFOS"),
    ("PFBA", "3", "carboxylate", "log Kmem/w", "< 1.7",
     "log (L water / kg membrane)",
     "below what the setup could measure: the liposome concentration could "
     "not be raised past 10 g/L, so the bound fraction stayed under 20%"),
    ("PFCAs and PFSAs (homologous series)", "", "both",
     "slope of log Kmem/w per perfluorinated carbon", "0.63",
     "log per CF2", "global fit across both series; fitted separately, 0.63 "
     "for PFCAs and 0.61 for PFSAs. Directly comparable to the per-carbon "
     "slopes this project computes for albumin and L-FABP binding"),
    ("PFSAs vs PFCAs at equal perfluorinated carbon number", "", "both",
     "head-group offset in log Kmem/w", "1.2", "log",
     "PFSAs sorb about 1.2 log units more strongly to the membrane than "
     "PFCAs with the same number of perfluorinated carbons, from the "
     "intercepts of the two fits"),
    ("all measured compounds", "", "both", "log Kmem/w range", "2.3 to 5.1",
     "log (L water / kg membrane)",
     "thirteen compounds: six PFCAs, three PFSAs and four alternatives. "
     "PFUnDA is off the linear trend and the authors cannot say whether "
     "that is a real cut-off effect or its high sorption to glassware"),
]


def main():
    rows = []
    for ch, nc, head, q, v, u, note in ROWS:
        rows.append({
            "study": STUDY, "year": "2020", "pmid": "32212724", "doi": DOI,
            "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "quantity": q, "value": v,
            "value_units": u, "method": SYSTEM, "note": note,
        })
    new = ["﻿"] if False else []
    buf = []
    import io
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader(); w.writerows(rows)
    text = s.getvalue()
    if os.path.exists(OUT) and open(OUT, encoding="utf-8").read() == text:
        print("no change needed (Ebert 2020 already ingested)")
        return 0
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"Ebert 2020: {len(rows)} rows -> "
          f"{os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print("   5 per-compound log Kmem/w values stated in the text")
    print("   +0.63 log per perfluorinated carbon; PFSAs +1.2 log over PFCAs")
    print("   per-compound SI tables S3.1.1 / S3.2.1 NOT supplied - a live "
          "request")
    return 0


if __name__ == "__main__":
    sys.exit(main())
