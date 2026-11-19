#!/usr/bin/env python3
"""Stage 2: where each transporter sits, and which way it moves PFAS.

    python3 -I scripts/build_transporter_localisation.py

Writes db/transporter_localisation.csv, the authority table the kinetics and
direction tables are joined against. One row per (transporter, species,
tissue), because DIRECTION IS NOT A PROPERTY OF THE TRANSPORTER ALONE:

  * OATP2B1 is apical in human kidney, where it reabsorbs, and sinusoidal in
    human liver, where it takes PFAS out of blood into the hepatocyte.
  * Oatp1a1 is apical in rat kidney and sinusoidal in rat liver.
  * OAT2 is BASOLATERAL in human kidney and APICAL in rat kidney -- the same
    gene, opposite sides, so the same transport would secrete in one species
    and reabsorb in the other.

Nor is it a property of the compound. A transporter's membrane and polarity
are fixed; what varies by compound is only whether it is a substrate.

PRIMARY SOURCE FOR THE RENAL ROWS is the state-of-the-science review:

    Niu S, Cao Y, Chen R, Bedi M, Sanders AP, Ducatman A, Ng C. "A
    State-of-the-Science Review of Interactions of Per- and Polyfluoroalkyl
    Substances (PFAS) with Renal Transporters in Health and Disease."
    Environ Health Perspect 2023;131(7):076002. doi:10.1289/EHP11885
    PMID 37418334.  (This is the review Argoul 2026 cites as "Niu et al.
    2023". NOTE: papers/ holds a second copy of this same paper misnamed
    "Dzierlenga 2023 state-of-science ..."; same PMID, same DOI. Nothing
    cites it, and the Dzierlenga entries in db/ are that author's own
    half-life papers from 2019/2020, not this one.)

It states the renal polarity directly: "OAT1 and OAT3 reside in the
basolateral membrane of the proximal tubular cells and therefore participate
in PFAS renal secretion. In contrast, OAT4, URAT1, and OATP2B1 are involved
in PFAS renal reabsorption because they are expressed in the apical membrane
of proximal tubular cells." And for rat: "Oat1 and Oat3 are expressed at the
basolateral membrane ... whereas Oatp1a1 resides in the apical membrane."

THE ROW THAT MATTERS MOST FOR AIM 1 IS AN ABSENCE. OAT4 is human-specific.
Niu 2023: "OAT4 is a human-specific transporter that was previously
demonstrated to transport PFOA; thus, it may be one of the transporters
responsible for species-specific differences in PFAS toxicokinetics." EFSA
puts it the same way: "OAT4 (not expressed in rats) and urate transporter 1
(URAT1), but not Oatp1a2 (the closest orthologous of rat Oatp1a1), were shown
to mediate PFOA". So the rat has no counterpart to the apical reabsorptive
transporter that carries six of this compilation's measured human clearances.
That is a presence/absence difference between species, not an affinity one,
and it is a different kind of answer from the one section 8 reached from
Nakagawa and Argoul -- those two showed affinity and substrate specificity
are conserved WHERE A TRANSPORTER EXISTS IN BOTH. Both can be true.

TISSUE ROLES beyond the kidney follow OEHHA's PHG (March 2024, Appendix 6 and
section 4.3), which names three: renal secretion and reabsorption; hepatic
re-uptake of PFAS secreted in bile, via bile-acid transporters, giving
enterohepatic circulation; and placental OAT4 moving PFAS toward the fetus.

EFFECT ON HALF-LIFE is derived, not quoted, and is recorded as the direction
of the step rather than as a claim about the whole organism. Apical uptake
from the tubular lumen returns PFAS to the body and lengthens; basolateral
uptake from blood is the first step of secretion and shortens; biliary efflux
shortens only if the PFAS is not then reabsorbed, which is exactly what ASBT
and OST alpha/beta do, so the hepatobiliary loop is marked as a cycle.

Idempotent.
"""
import csv, io, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                   "transporter_localisation.csv")

NIU = ("Niu 2023 review, Environ Health Perspect 131:076002, "
       "doi:10.1289/EHP11885")
EFSA = "EFSA CONTAM opinion (OAT4 not expressed in rats)"
OEHHA = "OEHHA PFOA/PFOS PHG, March 2024, section 4.3 and Appendix 6 Table A6.3"
LJUB = "Ljubojevic 2007, via db/transporter_mechanism.csv"

FIELDS = ["transporter", "gene", "species", "tissue", "membrane",
          "transport_step", "net_process", "effect_on_half_life",
          "transports_pfas", "note", "source"]

#  transporter, gene, species, tissue, membrane, step, process, effect,
#  transports, note, source
ROWS = [
    # ---------------- kidney, human ----------------
    ("OAT1", "SLC22A6", "human", "kidney proximal tubule", "basolateral",
     "blood -> cell", "renal secretion (step 1 of 2)", "shortens", "yes",
     "", NIU),
    ("OAT3", "SLC22A8", "human", "kidney proximal tubule", "basolateral",
     "blood -> cell", "renal secretion (step 1 of 2)", "shortens", "yes",
     "", NIU),
    ("OAT2", "SLC22A7", "human", "kidney proximal tubule", "basolateral",
     "blood -> cell", "renal secretion (step 1 of 2)", "shortens", "no",
     "Nakagawa 2008 found neither human OAT2 nor rat Oat2 transports PFOA. "
     "Predominantly hepatic in human", NIU),
    ("OAT4", "SLC22A11", "human", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "yes",
     "HUMAN-SPECIFIC: no rat orthologue. Niu 2023 calls it \"a human-specific "
     "transporter ... it may be one of the transporters responsible for "
     "species-specific differences in PFAS toxicokinetics\". It carries 6 of "
     "the measured human clearances in this compilation (Louisse 2023)",
     NIU + "; " + EFSA),
    ("URAT1", "SLC22A12", "human", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens",
     "disputed",
     "Yang 2010 reports Km 64.1 uM but only with ZERO extracellular chloride, "
     "a non-physiological condition. Louisse 2023 and Argoul 2026 found no "
     "URAT1 transport; Argoul's negative is weakened by values below the "
     "limit of quantification", NIU),
    ("OATP2B1", "SLCO2B1", "human", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "yes",
     "Same gene is sinusoidal in liver, where it does the opposite job",
     NIU),
    ("OATP1A2", "SLCO1A2", "human", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "no",
     "The closest human homologue of rat Oatp1a1, and the obvious candidate "
     "for a human counterpart to rat apical reabsorption. Yang 2010 found it "
     "does NOT mediate saturable PFOA uptake; Argoul 2026 found no uptake of "
     "any of ten PFAS", NIU + "; " + EFSA),
    ("OSTalpha/beta", "SLC51A/SLC51B", "human", "kidney proximal tubule",
     "basolateral", "cell -> blood (efflux; bidirectional)",
     "renal reabsorption (step 2 of 2)", "lengthens", "disputed",
     "Niu 2023: the only human renal EFFLUX transporter tested, and shown to "
     "transport PFBS, PFHxS and PFOS, \"indicating it was involved in renal "
     "reabsorption of these PFAS\", and to \"transport bidirectionally, "
     "depending on the concentration gradient\". Zhao 2015 positive; Argoul "
     "2026 found no uptake", NIU),
    ("P-gp", "ABCB1", "human", "kidney proximal tubule", "apical",
     "cell -> lumen", "renal secretion (step 2 of 2)", "shortens", "yes",
     "Niu 2026 measured Km 1.41-2.78 uM", NIU),
    ("BCRP", "ABCG2", "human", "kidney proximal tubule", "apical",
     "cell -> lumen", "renal secretion (step 2 of 2)", "shortens", "yes",
     "Niu 2026 measured Km 1.29-7.92 uM; its kidney-cortex protein was below "
     "detection in the same paper", NIU),
    ("MRP2", "ABCC2", "human", "kidney proximal tubule", "apical",
     "cell -> lumen", "renal secretion (step 2 of 2)", "shortens", "partial",
     "Zhao: PFHxS and PFOS inhibit, PFBS does not", NIU),
    ("OCT2", "SLC22A2", "human", "kidney proximal tubule", "basolateral",
     "blood -> cell", "organic CATION secretion", "not applicable", "no",
     "A cation transporter; PFAS are anions at physiological pH. Tested by "
     "Ryu 2024 and found not to transport them", NIU),
    # ---------------- kidney, rat ----------------
    ("Oat1", "Slc22a6", "rat", "kidney proximal tubule", "basolateral",
     "blood -> cell", "renal secretion (step 1 of 2)", "shortens", "yes",
     "", NIU),
    ("Oat3", "Slc22a8", "rat", "kidney proximal tubule", "basolateral",
     "blood -> cell", "renal secretion (step 1 of 2)", "shortens", "yes",
     "", NIU),
    ("Oat2", "Slc22a7", "rat", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "no",
     "OPPOSITE MEMBRANE TO THE HUMAN ORTHOLOGUE. Ljubojevic 2007 localises "
     "rat Oat2 to the brush border of the S3 segment, female-predominant; "
     "Kudo 2002 and Weaver 2010 agree it is apical. It would therefore be "
     "reabsorptive IF it transported PFAS -- but Weaver 2010 and Nakagawa "
     "2008 both found it does not, which undercuts the Kudo 2002 inference "
     "that Oat2 drives the rat sex difference", LJUB + "; " + NIU),
    ("Oatp1a1", "Slco1a1", "rat", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "yes",
     "The central rat reabsorption story. Its closest human homologue, "
     "OATP1A2, does not transport PFOA", NIU),
    ("Oatp1a5", "Slco1a5", "rat", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "yes",
     "Zhao 2017 measured PFBS, PFHxS and PFOS", NIU),
    ("Urat1", "Slc22a12", "rat", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "no",
     "Weaver 2010: did not transport PFHpA, PFOA, PFNA or PFDA", NIU),
    ("Npt2 (NaPi-IIa)", "Slc34a1", "rat", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens", "no",
     "Katakura 2007: not involved in PFOA transport", NIU),
    ("Asbt", "Slc10a2", "rat", "kidney proximal tubule", "apical",
     "lumen -> cell", "renal reabsorption (step 1 of 2)", "lengthens",
     "disputed",
     "Zhao 2015 found no transport of PFBS, PFHxS or PFOS; Argoul 2026 found "
     "PFBS and PFHxS ARE substrates. Unresolved, and it matters because this "
     "is a reabsorptive position", NIU),
    ("Ostalpha/beta", "Slc51a/Slc51b", "rat", "kidney proximal tubule",
     "basolateral", "cell -> blood (efflux; bidirectional)",
     "renal reabsorption (step 2 of 2)", "lengthens", "disputed",
     "ADDED AFTER THE JOIN FLAGGED IT: Argoul 2026 tested rat Ostalpha/beta "
     "as a RENAL transporter and found no uptake of any of ten PFAS, so the "
     "authority needed a rat kidney row to mirror the human one. Zhao 2015 "
     "reported rat uptake of PFBS, PFHxS and PFOS", NIU),
    ("Mrp2", "Abcc2", "rat", "kidney proximal tubule", "apical",
     "cell -> lumen", "renal secretion (step 2 of 2)", "shortens", "no",
     "Katakura 2007: Mrp2-deficient rats had renal PFOA clearance similar to "
     "wild type", NIU + "; " + OEHHA),
    ("OAT4", "SLC22A11", "rat", "kidney proximal tubule", "ABSENT",
     "not applicable", "not applicable", "not applicable", "gene absent",
     "THE ASYMMETRY. The rat has no orthologue of the apical reabsorptive "
     "transporter that carries six measured human clearances here. A "
     "presence/absence species difference, distinct from the affinity "
     "question section 8 settles", EFSA + "; " + NIU),
    # ---------------- liver ----------------
    ("NTCP", "SLC10A1", "human", "liver", "sinusoidal (basolateral)",
     "blood -> hepatocyte", "hepatic uptake", "lengthens (sequestration and "
     "enterohepatic recycling)", "yes",
     "Carries the highest measured clearance in this compilation, PFOS at "
     "236 uL/min/mg", OEHHA),
    ("Ntcp", "Slc10a1", "rat", "liver", "sinusoidal (basolateral)",
     "blood -> hepatocyte", "hepatic uptake", "lengthens (sequestration and "
     "enterohepatic recycling)", "yes", "", OEHHA),
    ("OATP1B1", "SLCO1B1", "human", "liver", "sinusoidal (basolateral)",
     "blood -> hepatocyte", "hepatic uptake", "lengthens (sequestration and "
     "enterohepatic recycling)", "yes", "", OEHHA),
    ("OATP1B3", "SLCO1B3", "human", "liver", "sinusoidal (basolateral)",
     "blood -> hepatocyte", "hepatic uptake", "lengthens (sequestration and "
     "enterohepatic recycling)", "yes", "", OEHHA),
    ("OATP2B1", "SLCO2B1", "human", "liver", "sinusoidal (basolateral)",
     "blood -> hepatocyte", "hepatic uptake", "lengthens (sequestration and "
     "enterohepatic recycling)", "yes",
     "Same gene is APICAL in kidney, where it reabsorbs instead", OEHHA),
    ("BSEP", "ABCB11", "human", "liver", "canalicular (apical)",
     "hepatocyte -> bile", "biliary efflux", "cycle (shortens only if not "
     "reabsorbed downstream)", "tested",
     "Zhao: PFOS inhibits taurocholate transport at 100 uM", OEHHA),
    ("MRP2", "ABCC2", "human", "liver", "canalicular (apical)",
     "hepatocyte -> bile", "biliary efflux", "cycle (shortens only if not "
     "reabsorbed downstream)", "partial", "", OEHHA),
    ("OSTalpha/beta", "SLC51A/SLC51B", "human", "liver", "basolateral",
     "hepatocyte -> blood (efflux)",
     "hepatic efflux back to the circulation",
     "lengthens (returns PFAS to blood instead of letting it leave in bile)",
     "yes",
     "ADDED AFTER THE JOIN FLAGGED IT. Zhao 2015's rows describe "
     "\"basolateral efflux out of enterocyte/hepatocyte into blood\", so the "
     "authority needed a liver row as well as an ileal one; without it those "
     "rows resolved to nothing. Same gene, same membrane polarity, different "
     "organ", OEHHA),
    ("Ostalpha/beta", "Slc51a/Slc51b", "rat", "liver", "basolateral",
     "hepatocyte -> blood (efflux)",
     "hepatic efflux back to the circulation",
     "lengthens (returns PFAS to blood instead of letting it leave in bile)",
     "disputed",
     "Zhao 2015 positive; Argoul 2026 detected no uptake by rat Ostalpha/beta",
     OEHHA),
    # ---------------- intestine, the return leg ----------------
    ("ASBT", "SLC10A2", "human", "ileum", "apical",
     "gut lumen -> enterocyte", "intestinal reabsorption (enterohepatic "
     "step 1 of 2)", "lengthens", "yes",
     "OEHHA: bile-acid transporters \"appear to mediate re-uptake of PFAS "
     "following secretion in bile, likely contributing to enterohepatic "
     "circulation\". This is the step that turns biliary efflux into a cycle "
     "rather than an elimination route", OEHHA),
    ("OSTalpha/beta", "SLC51A/SLC51B", "human", "ileum", "basolateral",
     "enterocyte -> portal blood", "intestinal reabsorption (enterohepatic "
     "step 2 of 2)", "lengthens", "yes", "", OEHHA),
    ("Ostalpha/beta", "Slc51a/Slc51b", "rat", "ileum", "basolateral",
     "enterocyte -> portal blood", "intestinal reabsorption (enterohepatic "
     "step 2 of 2)", "lengthens", "disputed",
     "Zhao 2015 reported uptake of PFBS, PFHxS and PFOS; Argoul 2026 "
     "detected none", OEHHA),
    # ---------------- placenta ----------------
    ("OAT4", "SLC22A11", "human", "placenta", "basal (fetal-facing)",
     "transfer toward the fetus", "placental transfer", "not applicable",
     "yes",
     "Kummu 2015: OAT4 protein correlated NEGATIVELY with the transfer index "
     "of PFOA and PFOS, n = 4. OEHHA lists placental OAT4 as one of three "
     "proposed physiological roles", OEHHA),
]


def main():
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader()
    for r in ROWS:
        w.writerow(dict(zip(FIELDS, r)))
    text = s.getvalue()
    if os.path.exists(OUT) and open(OUT, encoding="utf-8").read() == text:
        print("no change needed (localisation table already current)")
        return 0
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    import collections
    by_sp = collections.Counter(r[2] for r in ROWS)
    by_pr = collections.Counter(r[7] for r in ROWS)
    print(f"wrote {len(ROWS)} rows to "
          f"{os.path.relpath(OUT, os.path.dirname(os.path.dirname(OUT)))}")
    print(f"   species: {dict(by_sp)}")
    print(f"   tissues: {sorted({r[3] for r in ROWS})}")
    for e, n in by_pr.most_common():
        print(f"   {n:3d}  {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
