#!/usr/bin/env python3
"""Measured renal transporter protein abundance, human vs rat vs mouse.

    python3 -I scripts/add_transporter_expression.py

Writes db/transporter_expression.csv. This is the expression side of IVIVE --
the quantity sections 8 and 10 named as the missing half of the species
question, and the scaling factor Stage 3 needs.

SOURCES, all three supplied by the reviewer or already held:

  Thakur A, Yue G, Ahire D, Mettu VS, Al Maghribi A, Ford K, Peixoto L,
  Leeder JS, Prasad B. "Sex and the Kidney Drug-Metabolizing Enzymes and
  Transporters: Are Preclinical Drug Disposition Data Translatable to
  Humans?" Clin Pharmacol Ther 2024;116(1):235-246. doi:10.1002/cpt.3277
  PMID 38711199.  78 DMETs, kidney membrane fractions, human (n=15, 11 male
  / 4 female), rat (n=10) and mouse (n=14), one assay, both sexes.

  Uchida Y, Toyohara T, Ohtsuki S, Moriyama Y, Abe T, Terasaki T. J Pharm
  Sci 2016;105(2):1011-1016. doi:10.1002/jps.24645  PMID 26367854.  Rat
  only, but SPLIT INTO BRUSH-BORDER AND BASOLATERAL FRACTIONS, which is the
  same apical/basolateral split db/transporter_localisation.csv uses.

  Niu S, Ma Y, Tewari A, Ng C. Environ Health (Wash) 2026;4(9):2150-2157.
  doi:10.1021/envhealth.6c00088.  Human kidney cortex, n=20, already in the
  record as a note on the kinetics rows.

A UNITS TRAP IN UCHIDA, worth recording because it nearly cost a factor of
1000. pdftotext silently drops the micro sign from this PDF, turning every
"fmol/ug protein" into "fmol/mg protein" -- in the table header, the caption
and the abstract alike, so the extracted text is internally consistent and
gives no hint of the loss. Reading the page as an image shows the header is
"Protein Expression Level (fmol/ug protein)". The PubMed abstract had it
right. An audit of all 68 PDFs in this session found Uchida is the only one
affected: it is the only file with zero micro signs AND unit strings of the
"fmol/mg" form. Lesson for the project: a quantity extracted by pdftotext
alone is not safe until the micro signs are accounted for.

WHAT THE NUMBERS SAY.

1. THE SECRETORY SIDE IS CONSERVED IN KIND AND DIFFERS IN DEGREE. Rat Oat1
   is 6.8-fold the human (32.8 vs 4.8 pmol/mg membrane protein) and rat Oat3
   1.5-fold (6.3 vs 4.3). Both species have both transporters; the rat simply
   has more.

2. THE REABSORPTIVE SIDE IS NOT THE SAME MACHINERY, AND THE PROTEOMICS SAYS
   SO DIRECTLY RATHER THAN BY CITATION. Human OAT4 is 0.5 pmol/mg and is NOT
   DETECTED in rat or mouse. Rat Oatp1a1 is 1.9 pmol/mg and is NOT DETECTED
   in human. Each species reabsorbs with a transporter the other does not
   have. Section 10 reached this from Niu 2023 and EFSA; here it is measured.

3. RAT Oatp1a1 IS 114-FOLD HIGHER IN MALES THAN FEMALES (p < 0.05), and
   Oatp1a3 7-fold. This is the mechanism for the rat sex difference, and it
   explains the shape of it rather than only its size: OEHHA PHG Table A6.4
   has male rats reabsorbing 93.7% of filtered PFOA while FEMALE RATS ARE IN
   NET TUBULAR SECRETION. A transporter that is effectively absent in females
   changes the SIGN of net tubular handling, which is exactly what Table A6.4
   shows and what a mere difference in abundance could not produce. It also
   fits Cheng 2006 ("in kidney, gender-divergent Oatp expression is
   exclusively caused by stimulation by androgens") and Kudo 2002 (castration
   raised renal clearance 14-fold).

4. NO HUMAN RENAL DMET PROTEIN SHOWED A SEX DIFFERENCE (male n=11, female
   n=4), against 28 in rat and 40 in mouse. Humans show no PFAS half-life sex
   difference either.

5. TWO CONFLICTS IN SECTION 4.3 ARE RESOLVED BY EXPRESSION. Asbt and
   OSTalpha/beta were NOT DETECTED in human or rat kidney membrane fractions
   (mouse has them). Zhao 2015 found transport by both in transfected cells;
   Argoul 2026 found none in cells expressing them. Both can be right: the
   proteins can carry PFAS when expressed, and are not present in the kidney
   to do it. A transfected-cell positive is a statement about capability, not
   about renal physiology.

6. HUMAN RENAL BCRP IS NOT DETECTED, by two independent laboratories --
   Thakur here and Niu 2026 ("below the detection limit"). Niu 2026's BCRP
   carries the HIGHEST CL_int in this whole compilation, 279 uL/min/mg for
   PFOA. Scaled by renal expression it contributes essentially nothing in
   human kidney. That is the clearest demonstration available of why Stage 3
   must not rank transporters by CL_int alone.

7. THE THREE SOURCES AGREE WHERE THEY OVERLAP. Human OAT1: 4.8 (Thakur) vs
   5.3 pmol/mg (Niu 2026). Human P-gp: 1.2 vs 2.1. Human BCRP: not detected
   in both. Rat Oat1: 32.8 in a total membrane fraction (Thakur) against 10.6
   fmol/ug in a purified basolateral fraction (Uchida) -- different
   preparations, same order.

Idempotent.
"""
import csv, io, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                   "transporter_expression.csv")
THAKUR = ("Thakur 2024, Clin Pharmacol Ther 116:235, doi:10.1002/cpt.3277, "
          "Table 1; kidney membrane fractions, human n=15, rat n=10, "
          "mouse n=14")
UCHIDA = ("Uchida 2016, J Pharm Sci 105:1011, doi:10.1002/jps.24645, "
          "Table 1; whole rat kidney, purified membrane fractions")
NIU = ("Niu 2026, Environ Health 4:2150, doi:10.1021/envhealth.6c00088; "
       "human kidney cortex, n=20")

FIELDS = ["transporter", "species", "sex", "value", "sd", "units",
          "preparation", "net_process", "note", "source"]
MEMB = "kidney membrane fraction (total)"
PM = "pmol/mg membrane protein"
SECR = "renal secretion (step 1 of 2)"
REAB = "renal reabsorption (step 1 of 2)"
EFF = "renal secretion (step 2 of 2)"

#  transporter, species, sex, value, sd, units, prep, process, note, source
ROWS = [
    # --- basolateral uptake: the secretory side ---
    ("OAT1", "human", "both", "4.8", "1.2", PM, MEMB, SECR,
     "Niu 2026 independently gives 5.3 +/- 1.9 pmol/mg in human kidney "
     "cortex (n=20)", THAKUR),
    ("Oat1", "rat", "both", "32.8", "13.7", PM, MEMB, SECR,
     "6.8-fold the human value in the same assay; 2-fold higher in male "
     "rats than female", THAKUR),
    ("Oat1", "mouse", "both", "19.1", "7.1", PM, MEMB, SECR,
     "2-fold higher in male mice than female", THAKUR),
    ("OAT3", "human", "both", "4.3", "1.1", PM, MEMB, SECR, "", THAKUR),
    ("Oat3", "rat", "both", "6.3", "1.1", PM, MEMB, SECR,
     "1.5-fold the human value; no significant sex difference in rat",
     THAKUR),
    ("Oat3", "mouse", "both", "6.5", "3.6", PM, MEMB, SECR,
     "3-fold higher in female mice than male", THAKUR),
    ("OAT2", "human", "both", "0.3", "0.1", PM, MEMB, SECR,
     "basolateral in human; transports no PFAS (Nakagawa 2008)", THAKUR),
    ("Oat2", "rat", "both", "3.1", "3.0", PM, MEMB, REAB,
     "APICAL in rat, so it would be reabsorptive -- but transports no PFAS. "
     "20-fold higher in female rats than male", THAKUR),
    # --- apical uptake: the reabsorptive side, where the species diverge ---
    ("OAT4", "human", "both", "0.5", "0.2", PM, MEMB, REAB,
     "NOT DETECTED in rat or mouse. Proteomic confirmation of the "
     "human-specific status that section 10 took from Niu 2023 and EFSA",
     THAKUR),
    ("OAT4", "rat", "both", "", "", PM, MEMB, REAB,
     "NOT DETECTED. The paper lists Oat4 among proteins 'absent in rats'",
     THAKUR),
    ("OAT4", "mouse", "both", "", "", PM, MEMB, REAB, "NOT DETECTED",
     THAKUR),
    ("Oatp1a1", "rat", "both", "1.9", "2.1", PM, MEMB, REAB,
     "NOT DETECTED IN HUMAN. The standard deviation exceeds the mean "
     "because the sexes differ 114-fold and both are pooled here", THAKUR),
    ("Oatp1a1", "rat", "male", "", "", PM, MEMB, REAB,
     "114-FOLD HIGHER THAN FEMALE RATS (p < 0.05). The mechanism for the "
     "rat sex difference: OEHHA PHG Table A6.4 has male rats reabsorbing "
     "93.7% of filtered PFOA while female rats are in NET TUBULAR "
     "SECRETION, and a transporter effectively absent in females changes "
     "the sign rather than the size", THAKUR),
    ("Oatp1a1", "mouse", "both", "3.1", "2.8", PM, MEMB, REAB,
     "ABSENT IN FEMALE MICE entirely", THAKUR),
    ("Oatp1a1", "human", "both", "", "", PM, MEMB, REAB,
     "NOT DETECTED. Thakur compares rodent Oatp1a1 with human OATP1B1 and "
     "finds only ~40% sequence similarity; Yang 2010 separately showed "
     "OATP1A2, the closest renal homologue, does not transport PFOA",
     THAKUR),
    ("Oatp1a3", "rat", "both", "0.7", "0.6", PM, MEMB, REAB,
     "7-fold higher in male rats than female; not detected in human or "
     "mouse", THAKUR),
    ("URAT1", "human", "both", "3.1", "1.2", PM, MEMB, REAB,
     "present in both species, but Weaver 2010 found rat Urat1 transports "
     "no PFAS and Louisse 2023 and Argoul 2026 found no human URAT1 "
     "transport either", THAKUR),
    ("Urat1", "rat", "both", "2.2", "0.7", PM, MEMB, REAB, "", THAKUR),
    ("Urat1", "mouse", "both", "10.8", "1.1", PM, MEMB, REAB,
     "3-fold higher than human", THAKUR),
    ("OATP4C1", "human", "both", "0.7", "0.2", PM, MEMB, REAB,
     "comparable across all three species (fold range 0.5-1.6)", THAKUR),
    ("Oatp4c1", "rat", "both", "1.0", "0.3", PM, MEMB, REAB, "", THAKUR),
    # --- apical efflux: the second step of secretion ---
    ("P-gp", "human", "both", "1.2", "0.3", PM, MEMB, EFF,
     "Abcb1a. Niu 2026 independently gives 2.1 +/- 0.8 pmol/mg in human "
     "cortex", THAKUR),
    ("P-gp", "rat", "both", "0.9", "0.3", PM, MEMB, EFF,
     "2-fold higher in male rats", THAKUR),
    ("BCRP", "human", "both", "", "", PM, MEMB, EFF,
     "NOT DETECTED, and Niu 2026 independently reports human renal BCRP "
     "below its detection limit. Niu 2026 measured the highest CL_int in "
     "this compilation for BCRP/PFOA, 279 uL/min/mg; scaled by renal "
     "expression it contributes essentially nothing in human kidney",
     THAKUR),
    ("Bcrp", "rat", "both", "7.1", "3.3", PM, MEMB, EFF,
     "2.5-fold higher in male rats", THAKUR),
    ("Bcrp", "mouse", "both", "46.1", "5.0", PM, MEMB, EFF, "", THAKUR),
    ("MRP2", "human", "both", "0.2", "0.1", PM, MEMB, EFF, "", THAKUR),
    ("Mrp2", "rat", "both", "0.1", "0.0", PM, MEMB, EFF, "", THAKUR),
    ("MRP4", "human", "both", "0.2", "0.1", PM, MEMB, EFF,
     "80% lower than in the rodents", THAKUR),
    ("Mrp4", "rat", "both", "0.8", "0.1", PM, MEMB, EFF, "", THAKUR),
    # --- the enterohepatic pair, absent from the kidney in both species ---
    ("ASBT", "human", "both", "", "", PM, MEMB, REAB,
     "NOT DETECTED in human kidney. Detected only in mouse (3.2 pmol/mg)",
     THAKUR),
    ("Asbt", "rat", "both", "", "", PM, MEMB, REAB,
     "NOT DETECTED in rat kidney. This reconciles section 4.3: Zhao 2015 "
     "found no transport by rat Asbt and Argoul 2026 found PFBS and PFHxS "
     "ARE substrates, both in transfected cells -- a statement about what "
     "the protein can do, not about whether it is there to do it", THAKUR),
    ("OSTalpha/beta", "human", "both", "", "", PM, MEMB, REAB,
     "NEITHER SUBUNIT DETECTED in human kidney; detected only in mouse "
     "(Ostalpha 0.2, Ostbeta 3.0 pmol/mg). Same reconciliation as Asbt",
     THAKUR),
    ("Ostalpha/beta", "rat", "both", "", "", PM, MEMB, REAB,
     "NEITHER SUBUNIT DETECTED in rat kidney", THAKUR),
    # --- Uchida's membrane-resolved rat values ---
    ("Oat1", "rat", "male", "10.6", "2.6", "fmol/ug protein",
     "purified BASOLATERAL membrane fraction", SECR,
     "Membrane-resolved, so directly comparable with the "
     "apical/basolateral split in db/transporter_localisation.csv. Same "
     "order as Thakur's 32.8 pmol/mg in a total membrane fraction",
     UCHIDA),
    ("Oat3", "rat", "male", "10.2", "1.0", "fmol/ug protein",
     "purified BASOLATERAL membrane fraction", SECR,
     "not significantly different from Oat1 in the same fraction", UCHIDA),
    ("Urat1", "rat", "male", "22.2", "1.9", "fmol/ug protein",
     "purified BRUSH-BORDER membrane fraction", REAB,
     "abundant on the apical membrane, yet Weaver 2010 found it transports "
     "no PFAS -- abundance is not relevance", UCHIDA),
    ("BCRP", "rat", "male", "40.3", "3.0", "fmol/ug protein",
     "purified BRUSH-BORDER membrane fraction", EFF,
     "the most abundant transporter in the rat brush border", UCHIDA),
    ("P-gp", "rat", "male", "3.22", "0.56", "fmol/ug protein",
     "purified BRUSH-BORDER membrane fraction", EFF, "", UCHIDA),
    ("Oatp1a1", "rat", "male", "", "", "fmol/ug protein",
     "purified membrane fractions", REAB,
     "NOT IN UCHIDA'S PANEL at all: of 28 transporters measured, Oatp1a1 "
     "appears neither among those quantified nor among those below the "
     "limit of quantification. Thakur 2024 is the source that has it",
     UCHIDA),
]


def main():
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader()
    for r in ROWS:
        w.writerow(dict(zip(FIELDS, r)))
    text = s.getvalue()
    if os.path.exists(OUT) and open(OUT, encoding="utf-8").read() == text:
        print("no change needed (expression table already current)")
        return 0
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"wrote {len(ROWS)} rows to db/transporter_expression.csv")
    print("   Oat1 rat/human 6.8x; Oat3 1.5x -- secretory side conserved")
    print("   OAT4 human only; Oatp1a1 rat/mouse only -- reabsorptive side "
          "is NOT the same machinery")
    print("   rat Oatp1a1 114-FOLD higher in males: the sex-difference "
          "mechanism, and it changes the SIGN")
    print("   Asbt and OSTa/b not detected in human OR rat kidney: resolves "
          "two section 4.3 conflicts")
    print("   human renal BCRP not detected by two laboratories, yet carries "
          "the compilation's highest CL_int")
    return 0


if __name__ == "__main__":
    sys.exit(main())
