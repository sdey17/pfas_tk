#!/usr/bin/env python3
"""Add the rat-clearance citations that Stage 18 and Stage 6 rely on.

    python3 -I scripts/add_stage18_citations.py

Appends to db/citations.csv, keyed on PMID. Idempotent: re-running changes
nothing.

WHY THESE WERE MISSING. All seven already carry weight in the project --
six of them supply rat or mouse total clearances in the TOTAL list of
scripts/stage6_consolidate_tubular.py, and Kim 2016 additionally supplies the
rat PFHxS renal clearance that is one of only three rat inputs to equation
(5) -- yet none had a row in db/citations.csv. Every field below comes from
PubMed metadata, not from memory.
"""
import csv
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

# Fields in db/citations.csv order. Metadata retrieved from PubMed.
NEW = [
    dict(study="Chang et al.", year="2008", pmid="18353799",
         doi="10.1093/toxsci/kfn057",
         authors="Chang SC; Das K; Ehresman DJ; Ellefson ME; Gorman GS; "
                 "Hart JA; Noker PE; Tan YM; Lieder PH; Lau C; Olsen GW; "
                 "Butenhoff JL",
         title="Comparative pharmacokinetics of perfluorobutyrate in rats, "
               "mice, monkeys, and humans and relevance to human exposure "
               "via drinking water",
         journal="Toxicological Sciences", pub_year="2008", volume="104",
         issue="1", pages="40-53", authority="PubMed"),
    dict(study="Olsen et al.", year="2009", pmid="19059455",
         doi="10.1016/j.tox.2008.11.008",
         authors="Olsen GW; Chang SC; Noker PE; Gorman GS; Ehresman DJ; "
                 "Lieder PH; Butenhoff JL",
         title="A comparison of the pharmacokinetics of "
               "perfluorobutanesulfonate (PFBS) in rats, monkeys, and humans",
         journal="Toxicology", pub_year="2009", volume="256", issue="1-2",
         pages="65-74", authority="PubMed"),
    dict(study="Chengelis et al.", year="2009", pmid="19429410",
         doi="10.1016/j.reprotox.2009.01.013",
         authors="Chengelis CP; Kirkpatrick JB; Myers NR; Shinohara M; "
                 "Stetson PL; Sved DW",
         title="Comparison of the toxicokinetic behavior of "
               "perfluorohexanoic acid (PFHxA) and nonafluorobutane-1-"
               "sulfonic acid (PFBS) in cynomolgus monkeys and rats",
         journal="Reproductive Toxicology", pub_year="2009", volume="27",
         issue="3-4", pages="400-406", authority="PubMed"),
    dict(study="Chang et al.", year="2012", pmid="21889587",
         doi="10.1016/j.reprotox.2011.07.002",
         authors="Chang SC; Noker PE; Gorman GS; Gibson SJ; Hart JA; "
                 "Ehresman DJ; Butenhoff JL",
         title="Comparative pharmacokinetics of perfluorooctanesulfonate "
               "(PFOS) in rats, mice, and monkeys",
         journal="Reproductive Toxicology", pub_year="2012", volume="33",
         issue="4", pages="428-440", authority="PubMed"),
    dict(study="Kim et al.", year="2016", pmid="27637925",
         doi="10.1016/j.fct.2016.09.017",
         authors="Kim SJ; Heo SH; Lee DS; Hwang IG; Lee YB; Cho HY",
         title="Gender differences in pharmacokinetics and tissue "
               "distribution of 3 perfluoroalkyl and polyfluoroalkyl "
               "substances in rats",
         journal="Food and Chemical Toxicology", pub_year="2016",
         volume="97", issue="", pages="243-255", authority="PubMed"),
    dict(study="Huang et al.", year="2019", pmid="31334035",
         doi="10.1016/j.toxrep.2019.06.016",
         authors="Huang MC; Dzierlenga AL; Robinson VG; Waidyanatha S; "
                 "DeVito MJ; Eifrid MA; Granville CA; Gibbs ST; Blystone CR",
         title="Toxicokinetics of perfluorobutane sulfonate (PFBS), "
               "perfluorohexane-1-sulphonic acid (PFHxS), and "
               "perfluorooctane sulfonic acid (PFOS) in male and female "
               "Hsd:Sprague Dawley SD rats after intravenous and gavage "
               "administration",
         journal="Toxicology Reports", pub_year="2019", volume="6",
         issue="", pages="645-655", authority="PubMed"),
    dict(study="Fan et al.", year="2026", pmid="42190831",
         doi="10.1016/j.envpol.2026.128427",
         authors="Fan X; Ben Y; Dong Z; Deng W; Li H; Wang Q",
         title="The pharmacokinetics of perfluorobutanoic acid in female rats",
         journal="Environmental Pollution", pub_year="2026", volume="403",
         issue="", pages="128427", authority="PubMed"),
]


def main():
    p = HERE / "db" / "citations.csv"
    rows = list(csv.DictReader(p.open()))
    cols = list(rows[0].keys())
    have = {r["pmid"] for r in rows if r["pmid"]}

    added = []
    for d in NEW:
        if d["pmid"] in have:
            continue
        d = dict(d)
        d["citation"] = (f"{d['authors']}. {d['title']}. {d['journal']} "
                         f"{d['pub_year']};{d['volume']}"
                         + (f"({d['issue']})" if d["issue"] else "")
                         + f":{d['pages']}.")
        rows.append({c: d.get(c, "") for c in cols})
        added.append(f"{d['study']} {d['year']}")

    if added:
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    print(f"db/citations.csv: {len(rows)} rows "
          f"({len(added)} added this run)")
    for a in added:
        print(f"  + {a}")
    if not added:
        print("  (all seven already present -- nothing to do)")

    missing = [r["study"] + " " + r["year"] for r in rows
               if not r["pmid"] and not r["doi"]]
    print(f"\nrows still lacking BOTH pmid and doi: {len(missing)}")
    for m in sorted(set(missing)):
        print(f"  ! {m}")


if __name__ == "__main__":
    main()
