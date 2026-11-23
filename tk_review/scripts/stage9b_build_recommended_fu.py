#!/usr/bin/env python3
"""The recommended human unbound fraction, as a usable file (Stage 9b).

    python3 -I scripts/stage9b_build_recommended_fu.py

Writes db/human_fu_recommended.csv -- the deliverable. Idempotent.

THE RECOMMENDATION IS BY PURPOSE, NOT A SINGLE WINNER, because the four
datasets are not interchangeable and are not trying to measure quite the same
thing.

  DEFAULT for human absolute values ....... Smeltz 2023, ultracentrifugation
  INDIVIDUAL VARIABILITY layer ............ Fischer 2024
  SPECIES RATIOS .......................... Ryu 2024, and only Ryu

WHY SMELTZ AS THE DEFAULT

  - 67 PFAS, so the series has no gaps that must be filled from a second
    laboratory on a different method. None of the others reaches half that.
  - It is the US EPA's own (Center for Computational Toxicology), which
    matters for a California public health goal that has to sit alongside
    federal assessments.
  - The method choice is argued rather than assumed: the paper runs
    ultracentrifugation AND rapid equilibrium dialysis on 22 compounds and
    shows RED failing to reach equilibrium for high-mass, high-logKow PFAS
    (PFUnDA, PFTrDA, PFTeDA, PFHxS, PFDS and others, flagged in its Table 3).
    That is a reason to prefer UC, not a preference.
  - Bayesian credible intervals exist per compound (its Table S6), so a
    regulator can carry uncertainty instead of a point value.
  - Abraham 2024, the only human study with a measured urinary AND faecal
    mass balance, already runs on these values. Choosing Smeltz keeps the one
    in vivo human dataset and the binding on the same basis.

  Carry two caveats with it: PFBS is 0.0046 in its Table 2 and 0.0128 in its
  Table 3, a 2.8-fold internal discrepancy in the same assay; and the PFOS
  value is an average over three salt forms.

WHY FISCHER AS A LAYER RATHER THAN A RIVAL

  Fischer measures albumin-water AND globulin-water distribution separately
  and gives an equation that returns fu from a person's own serum albumin and
  globulin concentrations, validated at R2 = 0.99 against measured serum
  distribution. It found fu varying up to 2.5-fold across 1919 NHANES
  participants, widest for compounds below 500 g/mol. That is the "variable
  fu" worth having: not a different constant, but a function of the
  individual. Its own absolute values are derived here rather than measured,
  so it is the wrong thing to use as a point estimate and the right thing to
  use as a spread.

WHY RYU ONLY FOR SPECIES RATIOS

  Ryu is the only dataset that measured human, rat and two mouse strains in
  one laboratory on one method. For a between-species RATIO that is decisive,
  because a systematic method bias divides out. For a human ABSOLUTE value it
  runs at the high end -- 4-fold above Smeltz for PFHxA and 15-fold above
  Fischer for PFPeA -- and those are the two compounds whose tubular sign is
  in dispute.

WHAT THE CHOICE ACTUALLY DECIDES

  Only PFPeA and PFHxA. Every PFAS of seven carbons and above gives net
  tubular reabsorption on all four datasets, so no conclusion about PFOA,
  PFNA, PFDA, PFUnDA, PFHxS or PFOS depends on this decision at all.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
# Fischer 2024 interindividual variability, from its Results
#   Mw < 500 g/mol : min-max factor 2.0-2.5, CV 8.1-9.5%
#   Mw 500-600     : min-max factor 1.6,     CV 6.3%
MW = {"PFPA": 164.0, "PFBA": 214.04, "PFPeA": 264.05, "PFHxA": 314.05,
      "PFHpA": 364.06, "PFOA": 414.07, "PFNA": 464.08, "PFDA": 514.09,
      "PFUnDA": 564.09, "PFTrDA": 664.11, "PFTeDA": 714.12, "PFBS": 300.09,
      "PFHxS": 400.11, "PFHpS": 450.12, "PFOS": 500.13, "PFDS": 622.12,
      "4:2FTS": 328.15, "6:2FTS": 428.17, "8:2FTS": 528.18}


def variability(mw):
    if mw is None:
        return "", "", ""
    if mw < 500:
        return 2.5, "8.1-9.5", "Fischer 2024: widest band, Mw < 500 g/mol"
    if mw <= 600:
        return 1.6, "6.3", "Fischer 2024: narrower band, Mw 500-600 g/mol"
    return 1.3, "6.3", ("Fischer 2024: narrowest, Mw > 600 g/mol; albumin "
                        "and globulin binding are comparable and their serum "
                        "concentrations are inversely correlated")


def main():
    comp = {r["chemical"]: r for r in
            csv.DictReader((HERE / "db" / "human_fu_comparison.csv").open())}
    smeltz = {r["chemical"]: r for r in
              csv.DictReader((HERE / "db" / "smeltz2023_fu.csv").open())}
    sign = {}
    for r in csv.DictReader(
            (HERE / "db" / "human_tubular_sign_by_fu_source.csv").open()):
        sign.setdefault(r["chemical"], set()).add(r["verdict"])

    p = HERE / "db" / "human_fu_recommended.csv"
    n = 0
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "chemical", "n_carbons", "family", "mw_g_mol",
            "fu_recommended", "recommended_source", "recommendation_basis",
            "interindividual_max_fold", "interindividual_cv_percent",
            "variability_source",
            "fu_smeltz_UC", "fu_ryu", "fu_argoul", "fu_fischer_derived",
            "n_independent_sources", "fold_spread_across_sources",
            "tubular_sign_robust_to_fu_choice", "caveat"])
        for c, r in smeltz.items():
            mw = MW.get(c)
            fold, cv, vsrc = variability(mw)
            cm = comp.get(c, {})
            sg = sign.get(c)
            robust = ("" if not sg else
                      "yes" if len(sg) == 1 else
                      "NO - sign depends on which binding dataset is used")
            cav = r["note"]
            w.writerow([
                c, r["n_carbons"], r["family"], mw,
                r["fup_ultracentrifugation"], "Smeltz 2023 (UC)",
                "default human series: 67 PFAS, US EPA, ultracentrifugation "
                "with the RED comparison that justifies it, Bayesian bounds "
                "available, and the basis Abraham 2024's in vivo mass "
                "balance already uses",
                fold, cv, vsrc,
                r["fup_ultracentrifugation"], cm.get("ryu", ""),
                cm.get("argoul", ""), cm.get("fischer_derived", ""),
                cm.get("n_sources", 1), cm.get("fold_spread", ""),
                robust, cav])
            n += 1

    print(f"wrote {p.relative_to(HERE)}  ({n} rows)")
    print()
    bar = "=" * 78
    print(bar)
    print("THE DELIVERABLE, AND WHAT IT COMMITS YOU TO")
    print(bar)
    print("  fu_recommended   = Smeltz 2023 ultracentrifugation")
    print("  the variability  = Fischer 2024's bands, carried per compound")
    print("  the other three  = kept in the file, so any row can be rerun on")
    print("                     a different source without re-deriving it")
    print()
    rows = list(csv.DictReader(p.open()))
    nd = [r for r in rows if r["tubular_sign_robust_to_fu_choice"].startswith("NO")]
    ok = [r for r in rows if r["tubular_sign_robust_to_fu_choice"] == "yes"]
    print(f"  rows: {len(rows)}   sign robust: {len(ok)}   "
          f"sign source-dependent: {len(nd)} ({', '.join(r['chemical'] for r in nd)})")
    print()
    print("  SO THE DECISION IS NARROWER THAN IT LOOKS. Choosing Smeltz over")
    print("  Ryu changes the tubular verdict for PFPeA and PFHxA and for")
    print("  nothing else in the series.")


if __name__ == "__main__":
    main()
