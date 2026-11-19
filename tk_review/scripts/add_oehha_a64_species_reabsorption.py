#!/usr/bin/env python3
"""OEHHA PHG Table A6.4: net renal tubular reabsorption of PFOA by species.

    python3 -I scripts/add_oehha_a64_species_reabsorption.py

Source: OEHHA, "Public Health Goals for PFOA and PFOS in Drinking Water",
March 2024, Appendix 6 Table A6.4, adapted there from Han et al. (2012).

STATUS OF THIS SOURCE. The reviewer has ruled that the PHG may be relied on
directly: it went through dedicated external peer review, so unlike Han 2012
Table 6 or OEHHA Table A6.3 -- both of which are compilations whose values
this database holds under the papers that measured them -- Table A6.4 is a
DERIVED ANALYSIS and is stored as OEHHA's own. It is not double-counting,
because no other row in this database holds a net-reabsorption percentage.

WHAT IT IS. For each species and sex, net tubular reabsorption computed as

    net reabsorption = fu x GFR - CLR          percent reabsorbed =
                                               net reabsorption / (fu x GFR)

with fu, the plasma unbound fraction, ASSUMED TO BE 0.02 FOR EVERY SPECIES.
Where fu x GFR < CLR the animal is in net tubular SECRETION and no
reabsorption percentage is defined.

    species / sex        CLR (mL/kg-day)   % reabsorbed   % excreted
    human                      0.06            99.8          0.2
    rat male (Kemper)         18.2             93.7          6.3
    rat male (Ohmori)         21.2             92.7          7.3
    rat female               666-1009       NET SECRETION      -
    mouse male                10              97.0           3.0
    mouse female              16              95.2           4.8
    macaque male              15              91.2           8.8
    macaque female            32              81.2          18.8
    dog male                  43              59             41
    dog female                50.8            52             48
    rabbit male/female       640-670        NET SECRETION      -

TWO THINGS TO READ OFF IT, NEITHER OF WHICH IS THE PERCENTAGE ITSELF.

1. The percentages compress the difference and the FRACTION ESCAPING does
   not. Human reabsorbs 99.8% against the male rat's 93.7%, which looks like
   six points; but the fraction that escapes is 0.2% against 6.3%, a 31-fold
   difference, and the renal clearances themselves differ 303-fold (0.06
   against 18.2 mL/kg-day). The right denominator is what escapes.

2. The rat sex difference is QUALITATIVE, not quantitative. Male rats
   reabsorb ~93%; female rats are in net tubular SECRETION, both in Kemper
   (2003) and Ohmori (2003), so there is no reabsorption percentage to
   report. Whatever drives the ~70-fold male/female half-life difference in
   rats changes the SIGN of net tubular handling, not merely its size. Humans
   show no comparable sex difference.

THE fu ASSUMPTION IS THE WEAK POINT, AND IT IS WORTH BEING EXPLICIT. OEHHA
assumes fu = 0.02 for every species. Argoul 2026 Table 1 measures human PFOA
fu as 0.17%, which is 0.0017 -- twelve-fold lower -- and the filtered load it
computes (Clfg 4.4 mL/kg-day) is correspondingly twelve-fold below OEHHA's
51.3. Yet BOTH reach ~99-100% reabsorption for human PFOA, because renal
clearance is tiny either way (0.06 vs 0.042 mL/kg-day). So the qualitative
conclusion is robust to a twelve-fold disagreement about fu; the absolute net
reabsorption in mL/kg-day is not, and should not be quoted without the
assumption attached.

Idempotent. Writes a primary extract only; touches no other table.
"""
import csv, io, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                   "primary_2026", "oehha2024_a64_species_reabsorption.csv")
SRC = ("OEHHA PFOA/PFOS PHG, March 2024, Appendix 6 Table A6.4 "
       "(adapted there from Han et al. 2012)")
FIELDS = ["species", "sex", "reference_in_table", "gfr_mL_day_kg",
          "clr_mL_kg_day", "net_reabsorption_mL_kg_day",
          "percent_reabsorbed", "percent_excreted", "note", "source"]
FU = "fu assumed 0.02 for every species; net reabsorption = fu x GFR - CLR"
SECR = ("NET TUBULAR SECRETION (fu x GFR < CLR), so no reabsorption "
        "percentage is defined")

ROWS = [
    ("human", "not separated", "this document (PHG Table 4.5.1)", "2570",
     "0.06", "51.3", "99.8", "0.2",
     "Argoul 2026 measures human PFOA fu as 0.0017, twelve-fold below the "
     "0.02 assumed here, which moves the filtered load from 51.3 to about "
     "4.4 mL/kg-day but leaves the percentage at ~99-100% either way", SRC),
    ("rat", "male", "Kemper 2003 (unpublished study)", "14400", "18.2", "270",
     "93.7", "6.3", FU, SRC),
    ("rat", "female", "Kemper 2003 (unpublished study)", "14400", "666", "",
     "", "", SECR, SRC),
    ("rat", "male", "Ohmori 2003", "14400", "21.2", "267", "92.7", "7.3",
     FU + "; CLR digitized from a graph", SRC),
    ("rat", "female", "Ohmori 2003", "14400", "1009", "", "", "",
     SECR + "; CLR digitized from a graph", SRC),
    ("mouse", "male", "Kudo and Kawashima 2003", "16700", "10", "324", "97",
     "3", FU, SRC),
    ("mouse", "female", "Kudo and Kawashima 2003", "16700", "16", "318",
     "95.2", "4.8", FU, SRC),
    ("Japanese macaque", "male", "Harada 2005a", "8500", "15", "155", "91.2",
     "8.8", FU, SRC),
    ("Japanese macaque", "female", "Harada 2005a", "8500", "32", "138",
     "81.2", "18.8", FU, SRC),
    ("dog", "male", "Kudo and Kawashima 2003", "5300", "43", "63", "59", "41",
     FU, SRC),
    ("dog", "female", "Kudo and Kawashima 2003", "5300", "50.8", "55", "52",
     "48", FU, SRC),
    ("rabbit", "male", "Kudo and Kawashima 2003", "4000", "640", "", "", "",
     SECR, SRC),
    ("rabbit", "female", "Kudo and Kawashima 2003", "4000", "670", "", "", "",
     SECR, SRC),
]


def main():
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader()
    for r in ROWS:
        w.writerow(dict(zip(FIELDS, r)))
    text = s.getvalue()
    if os.path.exists(OUT) and open(OUT, encoding="utf-8").read() == text:
        print("no change needed (Table A6.4 already ingested)")
        return 0
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"OEHHA PHG Table A6.4: {len(ROWS)} rows")
    print("   human 99.8% reabsorbed vs male rat 93.7% -- but 0.2% vs 6.3% "
          "escaping, a 31-fold difference")
    print("   renal clearance 0.06 vs 18.2 mL/kg-day, a 303-fold difference")
    print("   FEMALE RATS AND RABBITS ARE IN NET TUBULAR SECRETION: the rat "
          "sex difference changes the sign, not the size")
    return 0


if __name__ == "__main__":
    sys.exit(main())
