#!/usr/bin/env python3
"""Literature search against the four not-established items (Stage 7).

    python3 -I scripts/stage7_open_items_literature.py

Writes db/open_items_literature.csv. Idempotent.

The report's section 10 listed four things the measurements do not settle. This
searched PubMed, Europe PMC and Consensus against each. THE RESULT CHANGES ONE
OF THEM SUBSTANTIALLY: the human short-chain unbound fraction is not a gap in
the literature, it is a gap in what this review had found. Two further human
datasets exist, plus a mechanistic model that reframes the question.

Earlier in this project I reported, after a narrower search, that the OAT4 and
albumin items "do not exist". That was too strong for OAT4: half of it exists
and is confirmed twice over. The correction is recorded here rather than
quietly folded in.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

FINDINGS = [
    # item, status, source, what it supplies, doi, pmid
    ("11 human short-chain fu", "NEW DATASET FOUND",
     "Fischer FC, Ludtke S, Thackray C, Pickard HM, Haque F, Dassuncao C, "
     "Endo S, Schaider L, Sunderland EM (2024) Environ Sci Technol 58:1055-1063",
     "HSA-water AND globulin-water distribution coefficients for PFAS with "
     "4-13 perfluorinated carbons by C18 SPME. f_unbound in human serum "
     "0.01-1.2%: PFOA 0.061%, PFHxS 0.041%, PFOS 0.042%, PFHpS 0.035%. "
     "Predicts f_unbound from an individual's albumin and globulin "
     "concentrations, R2 = 0.99 against measured serum distribution.",
     "10.1021/acs.est.3c07415", "38166384"),
    ("11 human short-chain fu", "NEW DATASET FOUND",
     "Smeltz M, Wambaugh JF, Wetmore BA (2023) Chem Res Toxicol 36:870-881",
     "67 PFAS human plasma f_up by ultracentrifugation with LC-MS/MS and "
     "Bayesian uncertainty bounds, from the US EPA Center for Computational "
     "Toxicology. Reports RED failure for PFAS of high mass or logKow > 4 "
     "through failure to reach equilibrium.",
     "10.1021/acs.chemrestox.3c00003", "37184865"),
    ("11 human short-chain fu", "MECHANISM FOUND",
     "Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L (2024) "
     "Environ Sci Technol 58:9954-9966",
     "Binding of anionic PFAS to BSA and human plasma is HIGHLY "
     "CONCENTRATION-DEPENDENT. Specific binding dominates at a PFAS:protein "
     "molar ratio below 0.1; physiological is below 0.01; in vitro assays "
     "often dose in excess above 1, where nonspecific binding dominates.",
     "10.1021/acs.est.3c10824", "38804966"),
    ("12 OAT4 scaling factor", "KIDNEY SIDE CONFIRMED TWICE",
     "Prasad B, Johnson K, Billington S, Lee C, Chung GW, Brown CDA, "
     "Kelly EJ, Himmelfarb J, Unadkat JD (2016) Drug Metab Dispos 44:1920-1924",
     "OAT4 at 0.5 +/- 0.2 pmol/mg total membrane protein in human kidney "
     "cortex, n = 41, by quantitative targeted proteomics, alongside OAT1 "
     "5.3, OAT2 0.9, OAT3 3.5. INDEPENDENTLY IDENTICAL to Thakur 2024's "
     "0.5 +/- 0.2 from n = 15. The in vivo denominator of the REF exists and "
     "is replicated; the in vitro numerator still does not.",
     "10.1124/dmd.116.072066", "27621205"),
    ("13 transport with albumin", "STILL OPEN",
     "searched; nothing PFAS-specific found",
     "Tan 2024 remains the only demonstration of albumin-mediated uptake at "
     "OAT1/OAT3 and it used drugs. Fischer 2024 proposes the experiment in "
     "its closing paragraph, naming renal tubular reabsorption, so it is "
     "recognised as open by others too.",
     "", ""),
    ("14 monkey unbound fraction", "CONFIRMED ABSENT; CLEARANCE IS RICH",
     "Butenhoff JL, Kennedy GL, Hinderliter PM, Lieder PH, Jung R, "
     "Hansen KJ, Gorman GS, Noker PE, Thomford PJ (2004) Toxicol Sci "
     "82:394-406",
     "Cynomolgus PFOA, six-month oral dosing at 3/10/20 mg/kg plus a 10 "
     "mg/kg IV arm, with SERUM, URINE AND FAECES measured at steady state "
     "and a terminal half-life of 14-42 days. Urine is the primary route. So "
     "the monkey has the clearance AND the mass balance; it has no plasma "
     "protein binding measurement, and searches for one return nothing in "
     "any species of macaque.",
     "10.1093/toxsci/kfh302", "15470233"),
]


def main():
    p = HERE / "db" / "open_items_literature.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["open_item", "status", "source", "what_it_supplies",
                    "doi", "pmid"])
        for row in FINDINGS:
            w.writerow(row)

    bar = "=" * 78
    print(bar)
    print("ITEM 11  HUMAN SHORT-CHAIN fu -- THE VERDICT CHANGES")
    print(bar)
    print("  The report said the sign of T for C5-C7 is undetermined because")
    print("  three laboratories disagree. That stands. What does NOT stand is")
    print("  the implication that the evidence is exhausted.")
    print()
    print("  TWO MORE HUMAN DATASETS EXIST.")
    print()
    print("  Fischer 2024 (Harvard), Environ Sci Technol 58:1055. Measures")
    print("  albumin-water AND globulin-water distribution separately, for")
    print("  4 to 13 perfluorinated carbons. Human serum f_unbound:")
    print("      PFOA 0.061%   PFHxS 0.041%   PFOS 0.042%   PFHpS 0.035%")
    print()
    print("  Smeltz 2023 (US EPA), Chem Res Toxicol 36:870. SIXTY-SEVEN PFAS")
    print("  by ultracentrifugation, with Bayesian uncertainty bounds, from")
    print("  EPA's own computational toxicology centre. For a regulatory")
    print("  conversation this is the obvious reference dataset and this")
    print("  review did not have it.")
    print()
    print("  SO THERE ARE NOT THREE HUMAN fu SERIES BUT FIVE. For PFOA:")
    fu = [("Fischer 2024 (SPME)", 0.00061), ("Abraham 2024 (in vivo)", 0.001),
          ("Argoul 2025 (dialysis)", 0.0017), ("Ryu 2024 (dialysis)", 0.00245)]
    for n, v in fu:
        print(f"      {n:28}{v:.5f}")
    lo, hi = min(v for _, v in fu), max(v for _, v in fu)
    print(f"      spread {hi/lo:.1f}x")
    print()
    print("  AND THE DISAGREEMENT NOW HAS A CANDIDATE MECHANISM WITH NUMBERS.")
    print("  Qin 2024 shows PFAS binding to human plasma is concentration-")
    print("  dependent, with specific binding dominant below a PFAS:protein")
    print("  molar ratio of 0.1 and nonspecific binding taking over above 1.")
    print("  Fischer 2024 independently holds its ratios below 0.004-0.009,")
    print("  citing a physiological ratio of 0.0005 and non-linearity of BSA")
    print("  binding from 0.39 upward.")
    print()
    print("  Ryu dosed 5 uM into 4% albumin, about 600 uM: a ratio near")
    print("  0.008. Abraham's subject carried sub-nanomolar plasma, a ratio")
    print("  near 2e-6. BOTH SIT BELOW THE 0.39 THRESHOLD, so on Fischer's")
    print("  own criterion both are in the linear range and CONCENTRATION")
    print("  DOES NOT OBVIOUSLY EXPLAIN THE 5.6-FOLD PFHxA GAP. The")
    print("  mechanism is real and quantified; it just does not appear to be")
    print("  this gap's cause. Stating that is the point of checking.")
    print()
    print("  THE MORE USEFUL FINDING IS THAT fu IS NOT A CONSTANT AT ALL.")
    print("  Fischer predicts f_unbound from an individual's own albumin and")
    print("  globulin concentrations (R2 = 0.99) and finds it varies up to")
    print("  2.5-fold across 1919 NHANES participants, most for compounds")
    print("  below 500 g/mol -- which is the short-chain range in dispute.")
    print("  Compounds under 500 g/mol are over 98% albumin-bound; globulins")
    print("  matter only above 600.")
    print()
    print("  SO OEHHA'S SINGLE ASSUMED fu IS WRONG IN KIND AS WELL AS IN")
    print("  MAGNITUDE. It is not that 0.02 is the wrong constant; it is that")
    print("  there is no constant to be right about, and the variability is")
    print("  largest exactly where this review found the sign undetermined.")

    print()
    print(bar)
    print("ITEM 12  OAT4 SCALING -- HALF OF IT EXISTS, AND IS REPLICATED")
    print(bar)
    print("  Earlier in this project I reported that no absolute")
    print("  quantification of OAT4 exists. That was too strong.")
    print()
    print("    Prasad 2016, human kidney cortex, n = 41:  0.5 +/- 0.2 pmol/mg")
    print("    Thakur 2024, human kidney, n = 15:         0.5 +/- 0.2 pmol/mg")
    print()
    print("  Two laboratories, two cohorts, nine years apart, identical to")
    print("  the reported precision. The KIDNEY side of the relative")
    print("  expression factor is settled. What is missing is only the")
    print("  IN VITRO side: OAT4 abundance in the cell line that supplied the")
    print("  kinetics. That is one proteomics run on an existing cell line,")
    print("  not a new study.")

    print()
    print(bar)
    print("ITEM 13  TRANSPORT WITH ALBUMIN PRESENT -- STILL OPEN")
    print(bar)
    print("  Nothing PFAS-specific found. Tan 2024 remains the only")
    print("  demonstration of albumin-mediated uptake at OAT1 and OAT3 and it")
    print("  used drugs. Fischer 2024 proposes the same experiment in its")
    print("  closing paragraph and names renal tubular reabsorption as the")
    print("  target, so the gap is recognised outside this review as well.")

    print()
    print(bar)
    print("ITEM 14  MONKEY -- CONFIRMED ABSENT, AND THE PRIZE IS LARGER")
    print(bar)
    print("  No plasma protein binding measurement exists for any macaque,")
    print("  for any PFAS. Searches by species, by compound and by method")
    print("  return nothing.")
    print()
    print("  But the clearance side is richer than the report credited.")
    print("  Butenhoff 2004 measured cynomolgus PFOA in SERUM, URINE AND")
    print("  FAECES at steady state across three oral doses plus an IV arm,")
    print("  with a 14-42 day half-life and urine as the primary route. That")
    print("  is the same mass-balance design as Abraham 2024 in the human.")
    print()
    print("  SO ONE EQUILIBRIUM DIALYSIS RUN ON MACAQUE PLASMA WOULD NOT ADD")
    print("  A SPECIES WITH BOUNDED T -- IT WOULD ADD ONE WITH THE FULL")
    print("  THREE-TIER TREATMENT, INCLUDING THE FAECAL TERM. The report")
    print("  called this the highest-value single experiment available. This")
    print("  search raises rather than lowers that assessment.")

    print()
    print(bar)
    print("WHAT TO DO WITH THIS")
    print(bar)
    print("  Ingesting Fischer 2024 and Smeltz 2023 into")
    print("  db/plasma_unbound_fraction.csv would take the human series from")
    print("  three sources to five and let the monotonicity test of section 7")
    print("  run on all five rather than three. That is the next concrete")
    print("  step, and it is a reading task rather than a search.")
    print()
    print("  It does not pre-empt the OEHHA question. If anything it sharpens")
    print("  it: the question is no longer 'which of three series' but")
    print("  'which of five, and should a single value be used at all given")
    print("  that Fischer resolves fu to the individual'.")
    print(f"\nwrote {p.relative_to(HERE)}  ({len(FINDINGS)} rows)")


if __name__ == "__main__":
    main()
