#!/usr/bin/env python3
"""Thakur 2024 verified at source (Stage 23).

    python3 -I scripts/stage23_thakur_verified.py

Writes db/primary_2026/thakur2024_table1_transporters.csv. Idempotent.

WHY. Thakur 2024 is the sole source for the transporter-abundance section,
including the claim that the reabsorptive machinery is not shared between
human and rodent kidney. Until now it had been used at second hand. The two
other load-bearing sources checked at first hand in this project (Ohmori's
figure, Kim's table) each yielded a correction, so the base rate justified
checking this one too.

THE NUMBERS ARE CORRECT. All five transporter rows in the paper reproduce
Thakur's Table 1 exactly, as do the sample sizes (human n=15, rat n=10,
mouse n=14) and the 114-fold male:female Oatp1a1 ratio quoted in the text.
No correction is needed to any value.

BUT THE PAPER UNDERSTATES ITS OWN CASE ON ONE POINT, AND OVERSTATES ANOTHER.

UNDERSTATES: the paper says human OAT4 is "undetectable" in rat and mouse and
rodent Oatp1a1 "undetectable" in human. That reads as a limit-of-detection
statement. Thakur's Figure 1 legend is stronger: OAT4 is ABSENT from the rat
and mouse genomes, and the rodent Oatp1a1/Oatp1a3 family has no close human
ortholog, being only ~40% similar to human OATP1B1. These are not proteins
that fell below a detection threshold; the genes are not there. The claim
that the reabsorptive machinery is not shared is therefore firmer than the
paper states, and should be worded as orthology rather than detection.

OVERSTATES: Thakur reports no sex-dependent difference in any human renal
DMET protein, which the paper has used as a contrast with the rodents. Thakur
itself attributes this to possible lack of power -- the human sample is 11
male and 4 female -- so it is an absence of evidence, not evidence of
absence, and must not be cited as a species contrast without that caveat.

AND IT ADDS AN INDEPENDENT MECHANISM FOR THE SEX FINDING. Thakur reports Oat2
20-fold higher in female than male rats and states explicitly that higher
female Oat2 drives higher PFOA excretion in female rats. Two transporters
therefore push the same way: Oatp1a1, apical and reabsorptive, 114-fold
higher in males, and Oat2, 20-fold higher in females. That is independent
molecular support for the sign change in the tubular factor between sexes
found from Kim's clearances in Stage 20, derived from a different laboratory
by a different method.

ONE METHOD CAVEAT. Abundances come from label-free global proteomics (data
independent acquisition, total protein approach), not targeted assays with
stable-isotope standards. Thakur flags that proteins under 0.5 pmol/mg need
confirmation in enriched fractions. Human OAT4 sits at exactly 0.5.
"""
import csv
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

# Thakur 2024 Table 1, transporters column, verified character by character
# against the PMC full text (PMC11218045). ND = not detected.
TABLE1 = [
    # protein, membrane, net process, human, rat, mouse
    ("OAT1", "basolateral", "secretory uptake", "4.8+/-1.2", "32.8+/-13.7",
     "19.1+/-7.1", ""),
    ("OAT3", "basolateral", "secretory uptake", "4.3+/-1.1", "6.3+/-1.1",
     "6.5+/-3.6", ""),
    ("OAT2", "basolateral", "secretory uptake", "0.3+/-0.1", "3.1+/-3",
     "2.3+/-0.6", "20-fold higher in FEMALE rats; Thakur links this to "
                  "higher PFOA excretion in females"),
    ("OAT4", "apical", "reabsorptive", "0.5+/-0.2", "ND", "ND",
     "GENE ABSENT in rat and mouse (Thakur Fig. 1 legend), not merely "
     "undetected; human value sits at the 0.5 pmol/mg confirmation threshold"),
    ("OATP1A1", "apical", "reabsorptive", "ND", "1.9+/-2.1", "3.1+/-2.8",
     "no close human ortholog: ~40% similar to human OATP1B1. 114-fold "
     "higher in MALE than female rats; absent in female mice"),
    ("OATP1A3", "apical", "reabsorptive", "ND", "0.7+/-0.6", "ND",
     "rat only; absent from the mouse genome; 7-fold higher in male rats"),
    ("URAT1", "apical", "urate exchange", "3.1+/-1.2", "2.2+/-0.7",
     "10.8+/-1.1", ""),
]

SAMPLE_SIZES = {"human": 15, "rat": 10, "mouse": 14}
SEX_SIZES = {"human": "11 M / 4 F", "rat": "5 M / 5 F", "mouse": "8 M / 6 F"}

PAPER_CLAIMS = [
    ("OAT1 human 4.8+/-1.2, rat 32.8+/-13.7, mouse 19.1+/-7.1", True),
    ("OAT3 human 4.3+/-1.1, rat 6.3+/-1.1, mouse 6.5+/-3.6", True),
    ("OAT4 human 0.5+/-0.2, rat ND, mouse ND", True),
    ("OATP1A1 human ND, rat 1.9+/-2.1, mouse 3.1+/-2.8", True),
    ("URAT1 human 3.1+/-1.2, rat 2.2+/-0.7, mouse 10.8+/-1.1", True),
    ("sample sizes human 15, rat 10, mouse 14", True),
    ("Oatp1a1 114-fold higher in male than female rats", True),
]


def main():
    bar = "=" * 78
    print(bar)
    print("THAKUR 2024 TABLE 1 VERIFIED AT SOURCE (PMC11218045)")
    print(bar)
    print(f"  {'protein':9}{'membrane':13}{'process':19}{'human':>11}"
          f"{'rat':>13}{'mouse':>12}")
    for pr, mem, proc, h, r, m, note in TABLE1:
        print(f"  {pr:9}{mem:13}{proc:19}{h:>11}{r:>13}{m:>12}")
    print()
    print(f"  n = {SAMPLE_SIZES}")
    print(f"  sex split: {SEX_SIZES}")
    print()
    print("  EVERY VALUE IN THE PAPER REPRODUCES THIS TABLE EXACTLY:")
    for claim, ok in PAPER_CLAIMS:
        print(f"    {'OK ' if ok else 'NO '} {claim}")
    print()

    print(bar)
    print("THREE REFINEMENTS THE SOURCE REQUIRES")
    print(bar)
    print("  (1) UNDERSTATEMENT. The paper calls OAT4 'undetectable' in rodent")
    print("      and Oatp1a1 'undetectable' in human. Thakur's Figure 1")
    print("      legend is stronger: OAT4 is ABSENT from the rat and mouse")
    print("      genomes, and rodent Oatp1a1/1a3 have no close human ortholog")
    print("      (~40% similarity to OATP1B1). The genes are not there. State")
    print("      this as orthology, not detection -- it is a firmer claim.")
    print()
    print("  (2) OVERSTATEMENT. 'No sex difference in any human renal DMET")
    print("      protein' is underpowered: the human sample is 11 male and 4")
    print("      female, and Thakur itself gives low female n as a possible")
    print("      cause. Absence of evidence. Must carry that caveat wherever")
    print("      it is used as a species contrast.")
    print()
    print("  (3) ADDITION. Oat2 is 20-fold higher in female than male rats,")
    print("      and Thakur states that higher female Oat2 drives higher PFOA")
    print("      excretion in female rats. With Oatp1a1 114-fold higher in")
    print("      males, TWO transporters push the same way:")
    print()
    print("        males   more apical reabsorption (Oatp1a1)")
    print("        females more basolateral uptake for secretion (Oat2)")
    print()
    print("      That is independent molecular support, from another")
    print("      laboratory by another method, for the sign change in the")
    print("      tubular factor between sexes that Stage 20 found in Kim's")
    print("      renal clearances: T = 0.099 male against 7.541 female for")
    print("      PFOA. Two methods, two laboratories, one conclusion.")
    print()
    print(bar)
    print("METHOD CAVEAT")
    print(bar)
    print("  Abundances come from label-free global proteomics (DIA, total")
    print("  protein approach), not targeted assays with stable-isotope")
    print("  standards, so absolute values carry more uncertainty than the")
    print("  quoted SDs suggest. Thakur flags proteins under 0.5 pmol/mg as")
    print("  needing confirmation in enriched fractions; human OAT4 sits at")
    print("  exactly 0.5. The ABSENCE of rodent OAT4 does not depend on this,")
    print("  being a genomic fact rather than a measurement.")
    print()

    p = HERE / "db" / "primary_2026" / "thakur2024_table1_transporters.csv"
    cols = ["protein", "membrane", "net_process", "human_pmol_mg",
            "rat_pmol_mg", "mouse_pmol_mg", "note", "n_human", "n_rat",
            "n_mouse", "source"]
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for pr, mem, proc, h, r, m, note in sorted(TABLE1):
            w.writerow(dict(
                protein=pr, membrane=mem, net_process=proc,
                human_pmol_mg=h, rat_pmol_mg=r, mouse_pmol_mg=m, note=note,
                n_human=15, n_rat=10, n_mouse=14,
                source="Thakur 2024 Table 1, Clin Pharmacol Ther "
                       "116:235-246, doi:10.1002/cpt.3277, verified at "
                       "source via PMC11218045"))
    print(f"wrote {p.relative_to(HERE)}  ({len(TABLE1)} rows)")


if __name__ == "__main__":
    main()
