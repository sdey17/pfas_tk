#!/usr/bin/env python3
"""Argoul 2026: ten PFAS against eight human and four rat renal transporters.

Argoul CML, Toutain PL, Picard-Hagen N, Klukovits A, Koncsos G, Gayrard V,
Lacroix MZ. "Comparative in vitro study of structure related uptake of ten
perfluoroalkyl substances by human and rat renal transporters."
Arch Toxicol 2026;100(11):5151-5158. doi:10.1007/s00204-026-04508-7
PMID 42446666.

The record held nothing from this paper; the roadmap listed it as request 1,
the species-difference backbone. The reviewer supplied the full text.

WHAT IT MEASURES, AND WHAT IT DOES NOT. Uptake RATIOS, not constants. A
mixture of eleven PFAS was applied at 0.2 uM and 2 uM each (0.2 and 0.02 uM
for PFDS), for 10 min, and a compound counts as a likely substrate if the
ratio of uptake in transporter-expressing to control cells exceeded 2 at
either concentration. So this paper adds substrate calls and direction, not
Km or Vmax, and nothing here enters the clearance column.

THE RESULT THAT MATTERS FOR AIM 1. The same molecular-weight threshold,
about 500 g/mol, governs transporter interaction in both species: PFCAs up
to nine carbons and PFSAs up to six are substrates of human OAT3/OAT4 and
rat Oat1/Oatp1a1, and the long-chain PFDA, PFOS and PFDS are substrates of
NONE of the twelve transporters in either species. The authors conclude that
"transporter selectivity may be largely conserved, and that additional
mechanisms may contribute to interspecies differences in PFAS renal
clearance". That is the same conclusion Nakagawa 2008 reached from matched
human/rat Km values, by a different route, and it argues against transporter
substrate specificity as the explanation for the ~50-fold species gap in
half-life. The authors name transporter EXPRESSION LEVEL as the missing
measurement.

THIS IS NOT "AGAINST EVERYONE", WHICH IS HOW THE ROADMAP HAD IT. It agrees
with Ryu 2024 on the chain-length threshold, with Louisse 2023 that URAT1
transports no PFAS, and with Louisse 2023/2024 on PFOA, PFNA and PFHxS. The
real disagreements are narrower and specific:

    OAT4/PFBS    substrate here; NOT a substrate in Louisse 2023
    OAT4/PFDA    not a substrate here; substrate in Louisse 2023 (Km 39)
    OAT1/PFOA    not a substrate here; substrate in Louisse 2024 (Km 185)
                 -- "the reasons for these discrepancies remain unclear"
    OAT3/PFBS    substrate here; not detected by Louisse 2024
    rat Oatp1a1/PFOS   no uptake here; substantial transport in Zhao 2017
    OSTa/b, ASBT       no uptake here for PFBS/PFHxS/PFOS; Zhao 2015
                       reported OSTa/b uptake of all three and PFOS by ASBT
    rat Asbt/PFBS, PFHxS   substrates here; Zhao 2015 found no transport

NEW: the first report of human OAT3-mediated uptake of GenX and PFO2OA.

TWO CAVEATS THE PAPER STATES. PFBA was tested but could not be quantified,
because solvent contamination prevented its separation from background, so
the eleven-compound mixture yields ten analysed compounds. And the URAT1
negative "should be interpreted with caution", because PFBS, PFHxS, PFOS and
GenX were below the limit of quantification in control or expressing cells
at both concentrations -- an absence of signal, not a measured absence of
transport. Both are carried on the affected rows.

A SIGN OF SATURATION AT 2 uM. OAT3 ratios for PFHxA, PFOA, PFBS and PFHxS,
and OAT4 ratios for PFHxA and PFBS, were LOWER at 2 uM than at 0.2 uM. The
authors read this as partial transporter saturation, which matters because
2 uM is only about tenfold above the highest human plasma levels.

TABLE 1 IS NOT TRANSPORTER DATA and is stored separately, at
db/primary_2026/argoul2026_table1_human_reabsorption.csv. It assembles human
unbound fractions with clearances from Abraham 2024 to estimate net tubular
reabsorption as a percentage of glomerular filtration clearance: 97-99% for
PFBS, PFHxS, PFOA, PFNA, PFOS and PFDA, 77% for GenX, and only 9.6% for
PFHxA. That is the in vivo counterpart of everything in this file.

Idempotent.
"""
import csv, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "transporter_kinetics.csv")
T1 = os.path.join(HERE, "..", "db", "primary_2026",
                  "argoul2026_table1_human_reabsorption.csv")
DOI = "10.1007/s00204-026-04508-7"
STUDY = ("Argoul CML, Toutain PL, Picard-Hagen N, Klukovits A, Koncsos G, "
         "Gayrard V, Lacroix MZ")
MEAS = ("transporter uptake ratio vs control cells (>= 2 at either "
        "concentration = likely substrate); 10-min uptake of an eleven-PFAS "
        "mixture at 0.2 and 2 uM each")

#  compound -> (n perfluorinated carbons, head group)
CPD = {
    "PFHxA": ("5", "carboxylate"), "PFOA": ("7", "carboxylate"),
    "PFNA": ("8", "carboxylate"), "PFDA": ("9", "carboxylate"),
    "PFBS": ("4", "sulfonate"), "PFHxS": ("6", "sulfonate"),
    "PFOS": ("8", "sulfonate"), "PFDS": ("10", "sulfonate"),
    "HFPO-DA": ("5", "ether carboxylate"), "PFO2OA": ("", "ether carboxylate"),
}
ORDER = ["PFHxA", "PFOA", "PFNA", "PFDA", "PFBS", "PFHxS", "PFOS", "PFDS",
         "HFPO-DA", "PFO2OA"]

#  transporter -> (gene, species, system, direction, set of substrates)
ACTIVE = {
    "OAT1": ("SLC22A6", "human", "HEK293",
             "secretion (basolateral uptake into the proximal tubule cell)",
             {"PFHxA"}),
    "OAT3": ("SLC22A8", "human", "HEK293",
             "secretion (basolateral uptake into the proximal tubule cell)",
             {"PFHxA", "PFOA", "PFNA", "PFBS", "PFHxS", "HFPO-DA", "PFO2OA"}),
    "OAT4": ("SLC22A11", "human", "HEK293",
             "reabsorption (apical, lumen -> cell)",
             {"PFHxA", "PFOA", "PFNA", "PFBS", "PFHxS"}),
    "Oat1": ("Slc22a6", "rat", "CHO-K1",
             "secretion (basolateral uptake into the proximal tubule cell)",
             {"PFHxA", "PFOA", "PFBS", "PFHxS", "HFPO-DA", "PFO2OA"}),
    "Oatp1a1": ("Slco1a1", "rat", "HEK293",
                "reabsorption (apical, lumen -> cell)",
                {"PFHxA", "PFOA", "PFNA", "PFBS", "PFHxS", "HFPO-DA",
                 "PFO2OA"}),
    "Asbt": ("Slc10a2", "rat", "HEK293",
             "reabsorption (apical, bile-acid transporter)",
             {"PFBS", "PFHxS"}),
}
#  transporter -> (gene, species, system, direction, note)
INACTIVE = {
    "ASBT": ("SLC10A2", "human", "HEK293", "reabsorption (apical)", ""),
    "OATP1A2": ("SLCO1A2", "human", "HEK293", "reabsorption (apical)",
                "included as the closest human homologue of rat Oatp1a1, the "
                "transporter proposed to drive PFOA reabsorption in rats; it "
                "transported none of the ten PFAS, which is the same "
                "negative Yang 2010 reported"),
    "OATP2B1": ("SLCO2B1", "human", "HEK293", "reabsorption (apical)", ""),
    "OSTalpha/beta": ("SLC51A/SLC51B", "human", "HEK293",
                      "basolateral efflux / return to the circulation",
                      "Zhao 2015 reported OSTalpha/beta-mediated uptake of "
                      "PFBS, PFHxS and PFOS; none was detected here"),
    "URAT1": ("SLC22A12", "human", "MDCKII",
              "reabsorption (apical, urate exchanger)",
              "INTERPRET WITH CAUTION: the paper states PFBS, PFHxS, PFOS "
              "and GenX were below the limit of quantification in control or "
              "expressing cell lysates at both concentrations, so this is an "
              "absence of signal rather than a measured absence of "
              "transport. Louisse 2023 also found no URAT1 transport"),
    "Ostalpha/beta": ("Slc51a/Slc51b", "rat", "HEK293",
                      "basolateral efflux / return to the circulation", ""),
}
SAT = ("Uptake ratio was LOWER at 2 uM than at 0.2 uM, which the authors read "
       "as partial transporter saturation; 2 uM is about tenfold above the "
       "highest human plasma concentrations reported")
SATURATING = {("OAT3", "PFHxA"), ("OAT3", "PFOA"), ("OAT3", "PFBS"),
              ("OAT3", "PFHxS"), ("OAT4", "PFHxA"), ("OAT4", "PFBS")}
LLOQ = {("OAT1", "PFHxA"), ("OAT1", "HFPO-DA"), ("Oat1", "PFHxA"),
        ("Oat1", "HFPO-DA"), ("Oatp1a1", "HFPO-DA"), ("Asbt", "PFBS"),
        ("OAT3", "HFPO-DA")}
LLOQ_NOTE = ("At the lower concentration the measurement was limited by "
             "values below the limit of quantification; the substrate call "
             "rests on the ratio exceeding 2 at the other concentration")
LONG = ("One of the three long-chain PFAS (PFDA, PFOS, PFDS, all above "
        "500 g/mol) that were substrates of NONE of the twelve transporters "
        "tested, in either species -- the paper's central structural result")
MIX = ("Measured in an eleven-PFAS mixture rather than singly; the authors "
       "argue the mixture did not compromise substrate identification "
       "because the calls agree with Ryu 2024, who tested single compounds "
       "at 1 uM")

TABLE1 = """pfas,mw_g_mol,fu_percent,cl_mL_kg_day,clfg_mL_kg_day,clr_mL_kg_day,net_reabsorption_mL_kg_day,net_reabsorption_percent_clfg
PFBS,300.1,1.7,1.9,44,1.2,43,97.3
PFHxA,314.06,2,57,52,47,5.0,9.6
HFPO-DA (GenX),330.05,4.2,43,109,25,84,77.1
PFO2OA,346.05,0.47,,12,,,
PFHxS,400.12,0.08,0.053,2.1,,2.0,97.5
PFOA,414.07,0.17,0.042,4.4,,4.4,99.0
PFNA,464.08,0.17,0.067,4.4,,4.4,98.5
PFOS,500.13,0.16,0.088,4.2,,4.1,97.9
PFDA,514.08,0.38,0.15,9.9,,9.7,98.5
"""


def build(fields):
    out = []
    for tr, (gene, sp, sys_, dirn, subs) in ACTIVE.items():
        for ch in ORDER:
            nc, head = CPD[ch]
            is_sub = ch in subs
            notes = [MIX]
            if not is_sub and ch in ("PFDA", "PFOS", "PFDS"):
                notes.insert(0, LONG)
            if is_sub and (tr, ch) in SATURATING:
                notes.insert(0, SAT)
            if is_sub and (tr, ch) in LLOQ:
                notes.insert(0, LLOQ_NOTE)
            r = {f: "" for f in fields}
            r.update({
                "study": STUDY, "year": "2026", "pmid": "42446666",
                "doi": DOI, "transporter": tr, "gene": gene,
                "species_orthologue": sp, "chemical": ch,
                "n_perfluorinated_carbons": nc, "head_group": head,
                "expression_system": f"{sys_} cells stably expressing {tr} "
                                     f"vs mock-transfected controls",
                "measurement_type": MEAS,
                "uptake_ratio": ">= 2 at 0.2 and/or 2 uM" if is_sub
                                else "< 2 at both concentrations",
                "substrate_or_inhibitor": ("likely substrate" if is_sub
                                           else "not a substrate"),
                "direction_reabsorption_or_secretion": dirn,
                "sex_or_hormone_effect": "not assessed",
                "assay_volume_note": ". ".join(notes),
                "source_table": ("Fig 1 (human) / Fig 2 (rat); per-compound "
                                 "uptake values in Supplemental Tables A2-A9"),
            })
            out.append(r)
    for tr, (gene, sp, sys_, dirn, note) in INACTIVE.items():
        r = {f: "" for f in fields}
        r.update({
            "study": STUDY, "year": "2026", "pmid": "42446666", "doi": DOI,
            "transporter": tr, "gene": gene, "species_orthologue": sp,
            "chemical": ", ".join(ORDER),
            "n_perfluorinated_carbons": "4-10", "head_group":
                "carboxylate, sulfonate and ether carboxylate",
            "expression_system": f"{sys_} cells stably expressing {tr} "
                                 f"vs mock-transfected controls",
            "measurement_type": MEAS,
            "uptake_ratio": "< 2 for all ten PFAS at both concentrations",
            "substrate_or_inhibitor": "no PFAS transported",
            "direction_reabsorption_or_secretion": dirn,
            "sex_or_hormone_effect": "not assessed",
            "assay_volume_note": ((note + ". ") if note else "") + MIX,
            "source_table": ("Supplemental Tables A2 (human) / A6 (rat); not "
                             "plotted, as no PFAS reached a ratio of 2"),
        })
        out.append(r)
    return out


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        fields, rows = rdr.fieldnames, list(rdr)
    old = [r for r in rows if (r.get("doi") or "") == DOI]
    new = build(fields)
    if not os.path.exists(T1) or open(T1, encoding="utf-8").read() != TABLE1:
        with open(T1, "w", encoding="utf-8") as fh:
            fh.write(TABLE1)
        print(f"wrote {os.path.relpath(T1, os.path.join(HERE, '..'))}")
    if old and len(old) == len(new) and all(
            a.get(k) == b.get(k) for a, b in zip(old, new) for k in fields):
        print("no change needed (Argoul 2026 already ingested)")
        return 0
    merged, placed = [], False
    for r in rows:
        if (r.get("doi") or "") == DOI:
            if not placed:
                merged.extend(new); placed = True
            continue
        merged.append(r)
    if not placed:
        merged.extend(new)
    shutil.copy2(DB, DB + ".bak")
    with open(DB, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(merged)
    subs = sum(1 for r in new if r["substrate_or_inhibitor"] == "likely substrate")
    print(f"Argoul 2026: {len(old)} row(s) -> {len(new)}")
    print(f"   {len(ACTIVE)} transporters x 10 PFAS = 60 per-compound calls, "
          f"{subs} positive")
    print(f"   {len(INACTIVE)} transporters transported nothing (1 row each)")
    print("   no Km or Vmax: this paper reports uptake ratios only")
    return 0


if __name__ == "__main__":
    sys.exit(main())
