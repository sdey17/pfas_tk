#!/usr/bin/env python3
"""
Add the seven protein-binding studies that appear in OEHHA 2024 PHG
Table A6.2 but were missing from db/protein_binding.csv.

All values are transcribed from OEHHA 2024, "Public Health Goals for PFOA and
PFOS in Drinking Water", March 2024, Appendix 6, Table A6.2, as saved in
tk_review/papers/OEHHA 2024 PFOA PFOS PHG.txt. None of the PRIMARY papers has
been retrieved, so every row is marked secondary and carries no assay
conditions. They are added so the compilation is complete against a regulatory
source, not because the numbers are independently verified here.

OEHHA's own footnotes, carried through:
  a  Kerstner-Wood et al. 2003 is "as reported in US EPA (2016b)"
  b  Jones et al. 2003 SHBP values are relative to estrogen
  c  Jones et al. 2003 CBG values are relative to cortisol
  d  Wu 2009 and Qin 2010 dissociation constants are calculated from reported
     Gibbs energies, not measured directly

Run:  python3 tk_review/scripts/add_oehha_a62_rows.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

SRC = ("secondary (OEHHA 2024 PFOA/PFOS PHG Table A6.2); primary paper not "
       "retrieved, no assay conditions available")

NEW = [
    # Jones et al. 2003 - binding relative to the endogenous ligand
    dict(study="Jones et al.", year="2003", chemical="PFOA",
         protein="sex-hormone binding protein (SHBP)", protein_species="fish",
         method="competitive binding relative to estrogen", kd="1.7e-5",
         kd_units="M", tk_relevance_noted="OEHHA footnote b: relative to estrogen"),
    dict(study="Jones et al.", year="2003", chemical="PFOS",
         protein="sex-hormone binding protein (SHBP)", protein_species="fish",
         method="competitive binding relative to estrogen", kd="2.1e-5",
         kd_units="M", tk_relevance_noted="OEHHA footnote b: relative to estrogen"),
    dict(study="Jones et al.", year="2003", chemical="PFOA",
         protein="corticosteroid-binding globulin (CBG)", protein_species="bird",
         method="competitive binding relative to cortisol", kd="1.2e-4-3.1e-4",
         kd_units="M", tk_relevance_noted="OEHHA footnote c: relative to cortisol"),
    dict(study="Jones et al.", year="2003", chemical="PFOS",
         protein="corticosteroid-binding globulin (CBG)", protein_species="bird",
         method="competitive binding relative to cortisol", kd="2.5e-4-2.7e-4",
         kd_units="M", tk_relevance_noted="OEHHA footnote c: relative to cortisol"),
    # Kerstner-Wood et al. 2003 - percent bound only, no constant
    dict(study="Kerstner-Wood et al.", year="2003", chemical="PFOA",
         protein="protein fraction of plasma",
         protein_species="rat; monkey; human",
         method="plasma protein binding (percent bound)",
         unbound_fraction="0-0.03 (97-100% bound)",
         tk_relevance_noted="OEHHA footnote a: as reported in US EPA (2016b)"),
    # Messina et al. 2005 - two-step binding
    dict(study="Messina et al.", year="2005", chemical="PFOA",
         protein="human serum albumin (HSA)", protein_species="human",
         method="two-step binding model", kd="6.9e-6", kd_units="M",
         n_binding_sites="two-step reaction, step 1",
         tk_relevance_noted="first of two dissociation constants"),
    dict(study="Messina et al.", year="2005", chemical="PFOA",
         protein="human serum albumin (HSA)", protein_species="human",
         method="two-step binding model", kd="3.2e-5", kd_units="M",
         n_binding_sites="two-step reaction, step 2",
         tk_relevance_noted="second of two dissociation constants"),
    # Li et al. 2009
    dict(study="Li et al.", year="2009", chemical="PFOS",
         protein="human serum albumin (HSA)", protein_species="human",
         method="not stated in OEHHA Table A6.2", kd="2e-4", kd_units="M"),
    # Salvalaglio et al. 2010 - in silico, percent only
    dict(study="Salvalaglio et al.", year="2010", chemical="PFOA",
         protein="human serum albumin (HSA, in silico)", protein_species="human",
         method="molecular modelling; predicted binding rate in blood",
         unbound_fraction="0.05 (95% bound, predicted)", is_computational="yes",
         tk_relevance_noted="predicted over a range of concentrations"),
    dict(study="Salvalaglio et al.", year="2010", chemical="PFOS",
         protein="human serum albumin (HSA, in silico)", protein_species="human",
         method="molecular modelling; predicted binding rate in blood",
         unbound_fraction="0.02 (98% bound, predicted)", is_computational="yes",
         tk_relevance_noted="predicted over a range of concentrations"),
    # Woodcroft et al. 2010 - L-FABP, two models
    dict(study="Woodcroft et al.", year="2010", chemical="PFOA",
         protein="liver fatty acid-binding protein (L-FABP)",
         protein_species="rat", method="displacement assay",
         kd="13.1e-6, 23.8e-6", kd_units="M",
         n_binding_sites="2 (by displacement)"),
    dict(study="Woodcroft et al.", year="2010", chemical="PFOA",
         protein="liver fatty acid-binding protein (L-FABP)",
         protein_species="rat", method="isothermal titration calorimetry (ITC)",
         kd="3.1e-6, 26.2e-6, 52.6e-6", kd_units="M",
         n_binding_sites="3 (by ITC)"),
    # EFSA 2018 - qualitative, referencing an unpublished study
    dict(study="EFSA", year="2018", chemical="PFOS",
         protein="albumin; low-density lipoproteins", protein_species="human",
         method="qualitative, referenced unpublished study",
         tk_relevance_noted="OEHHA A6.2 records this as 'strong binding'; no "
                            "constant given"),
    dict(study="EFSA", year="2018", chemical="PFOS",
         protein="alpha-globulin; gamma-globulin", protein_species="human",
         method="qualitative, referenced unpublished study",
         tk_relevance_noted="OEHHA A6.2 records this as 'medium to weak "
                            "binding'; no constant given"),
]


def main():
    rows = list(csv.DictReader(open(DB)))
    fields = list(rows[0].keys())
    before = len(rows)
    for rec in NEW:
        row = {k: "" for k in fields}
        row.update({k: v for k, v in rec.items() if k in fields})
        row["source_table"] = SRC
        row.setdefault("is_computational", "")
        if not row["is_computational"]:
            row["is_computational"] = "no"
        rows.append(row)
    with open(DB, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print("added %d rows (%d -> %d)" % (len(NEW), before, len(rows)))
    for rec in NEW:
        print("  %-24s %-6s %-7s %s" % (rec["study"], rec["year"],
                                        rec.get("chemical", ""),
                                        rec.get("kd", rec.get("unbound_fraction", ""))))


if __name__ == "__main__":
    main()
