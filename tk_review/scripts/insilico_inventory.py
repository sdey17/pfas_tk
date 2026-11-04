#!/usr/bin/env python3
"""The in silico PFAS literature, tagged, with a double-counting check.

WHY THIS IS A SEPARATE FILE.  db/protein_binding.csv holds
measurements and db/combined/ holds measurements.  A model is not a
measurement and does not belong in a Kd column beside one -- but it is
not nothing either, since the QSAR aim of this review is exactly about
what structure can and cannot predict.  So the models live here, with
the fields that decide whether a model means anything and which a
binding table has nowhere to put.

EVERY ENTRY IS TAGGED on four axes, so the file can be filtered rather
than read:

  kind          what sort of model it is
                  qsar-ml              statistical model, ML algorithms
                  qsar-regression      statistical model, regression
                  read-across          analogue-based inference
                  docking              structure-based scoring, not fitted
                  homology-model       a modelled protein structure
                  pbtk                 physiologically based TK/PK model
                  mechanistic-model    process model fitted to data
                  empirical-regression regression on one physical variable
                  quantum-chemical     first-principles property prediction

  data_origin   where the numbers it learned from came from
                  own        the paper's own measurements
                  borrowed   somebody else's measurements
                  none       nothing is fitted (docking, homology)
                  unknown    not established

  held_out      what was actually held back
                  external-test  a test set the model never saw
                  cv-only        cross-validation folds only
                  none           no held-out evaluation reported
                  n/a            nothing fitted, so nothing to hold out

  detail        how deeply this entry is characterised
                  read        the paper was read for these fields
                  db-record   tagged from the method recorded in
                              db/protein_binding.csv; the paper has not
                              been read for descriptors and validation

THE OVERLAP CHECK is the second half of this file.  A model trained on
somebody else's measurements is not independent evidence about those
measurements, and if this database also holds them, then counting both
is counting once twice.  That is not hypothetical: Gallagher 2023's
QSAR is trained on the Starnes DSF dataset, which this database holds
as experimental rows.  `overlap()` resolves each model's declared
training sources against the database and says what the exposure is,
rather than leaving it to whoever remembers.

Sources: each paper as cited.  Gallagher 2023 full text via PubMed
Central (PMC10383382).  Ryu 2024, Starnes 2025, Dawson 2023 and
Gao 2025 from the copies in papers/.  The db-record entries from the
method and source_table columns of db/protein_binding.csv.

Run:  python3 tk_review/scripts/insilico_inventory.py
Writes db/insilico_models.csv
"""
import csv
import os
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
HALFLIFE = [os.path.join(HERE, "..", "db", f) for f in
            ("animal_halflife_measured.csv", "human_halflife_extended.csv")]
OUT = os.path.join(HERE, "..", "db", "insilico_models.csv")

FIELDS = ["study", "year", "doi", "pmid", "kind", "data_origin",
          "held_out", "detail", "endpoint", "model", "trained_on",
          "n_train", "n_compounds", "descriptors", "n_feat",
          "n_features", "validation", "performance", "applied_to",
          "overlaps_db", "caveat"]


def m(**kw):
    kw.setdefault("pmid", "")
    for f in FIELDS:
        kw.setdefault(f, "")
    kw.setdefault("overlaps", [])
    return kw


M = [
    m(study="Starnes et al.", year="2025", doi="10.1093/toxsci/kfae124",
      pmid="39298512", kind="qsar-ml", data_origin="own",
      held_out="cv-only", detail="read",
      endpoint="HSA binding affinity (Ka) and binary bind/no-bind, "
               "from its own differential scanning fluorimetry assay",
      model="7 classifiers and 9 regressors (python/sklearn); best "
            "classifier logistic regression, best regressor support "
            "vector regression",
      trained_on="ITS OWN measurements: 65 PFAS from 20 chemical "
                 "classes plus 3 fatty acids as non-fluorinated "
                 "analogues; 57 chemicals reached the training set, "
                 "one congener dropped for missing descriptors",
      n_train="57", n_compounds="57 (54 PFAS + 3 fatty acids)",
      descriptors="Saagar (834) + ToxPrints (729) + OPERA (12) = 1575 "
                  "generated; invariant features removed and min-max "
                  "scaled to 255 used",
      n_feat="255", n_features="1575 generated, 255 used",
      validation="3-by-3-fold cross-validation. No external test set.",
      performance="classifier AUROC 0.936; regressor R2 0.854",
      applied_to="EPAPFASINV inventory, 430 PFAS",
      caveat="255 features on 57 chemicals, about 4.5 per compound, "
             "and the only held-out data is the cross-validation "
             "fold. The DSF endpoint itself fails the physiological "
             "bound in binding_kd_physiological_bound.py (81 % bound "
             "at six sites, where plasma needs >99 %)."),

    m(study="Gallagher et al.", year="2023",
      doi="10.3390/molecules28145375", pmid="37513249",
      kind="qsar-regression", data_origin="borrowed",
      held_out="external-test", detail="read",
      endpoint="HSA binding affinity (EC50, mM) and a high/low class "
               "at a 1.45 mM threshold",
      model="GA-LDA classification (DTC Lab "
            "ClassificationBasedQSAR v1.0.0); GA-MLR and PLS "
            "regression (SmallDataModeler v1.0); GenRA read-across "
            "on ToxPrint chemotypes; Glide XP docking on PDB 7AAI",
      trained_on="NOT ITS OWN. 24 PFAS EC50 values collected from the "
                 "literature; the paper's description (24 PFAS, mM "
                 "EC50 by 4-parameter variable slope, 4:2 and 6:2 "
                 "FTOH non-binding, H/L threshold 1.45 mM) matches "
                 "the Starnes 2024 cross-species DSF dataset",
      n_train="16 train / 8 test",
      n_compounds="24 classification / 22 regression",
      descriptors="alvaDesc 2.0.16: 753 2D descriptors offered for "
                  "classification, 736 for regression. Classification "
                  "retains 3 (average connectivity index of order 4, "
                  "an augmented edge-adjacency eigenvalue, an "
                  "eccentric topological index); regression retains 4 "
                  "(packing density index, Moran autocorrelation lag "
                  "8 mass-weighted, Geary autocorrelation lag 8 "
                  "vdW-volume-weighted, QED drug-likeness)",
      n_feat="753 -> 3-4", n_features="753 offered, 3-4 retained",
      validation="classification split 2:1 at random, training n = 16 "
                 "and test n = 8; regression used the whole set under "
                 "exhaustive double cross-validation, no split. "
                 "Docking validated by redocking PFOA, RMSD 0.092 A.",
      performance="ROC AUC 0.97 training, 1.00 on the 8-compound test "
                  "set",
      applied_to="not applied to an external inventory",
      overlaps=[("Starnes et al.", "2024",
                 "its training EC50 values are this dataset")],
      caveat="753 descriptors searched against 16 training compounds, "
             "and a test AUC of 1.00 on eight. The retained "
             "descriptors are abstract topological indices with no "
             "stated mechanism. The authors say themselves that the "
             "GenRA predictions are 'not that reliable': every AUC "
             "was 0 and p ran 0.73 to 1."),

    m(study="Ryu et al.", year="2024", doi="10.1021/acs.est.4c04050",
      pmid="39161261", kind="qsar-regression", data_origin="own",
      held_out="cv-only", detail="read",
      endpoint="unbound fraction (f_unbound) in liver, brain, heart, "
               "kidney and lung",
      model="partial least squares regression, four separate models "
            "(human, C57, mouse, rat)",
      trained_on="ITS OWN measurements: 16 PFAS across four "
                 "organisms, the 320 experimental rows this database "
                 "holds",
      n_train="16", n_compounds="16",
      descriptors="20 physicochemical properties from MOE 2018.01 and "
                  "Dragon 6.0. Molecular weight and log D are the ones "
                  "the paper reads a trend off.",
      n_feat="20", n_features="20",
      validation="leave-one-out cross-validation over the 16 PFAS; "
                 "grid search for the component count (4 for human, 3 "
                 "for C57, mouse and rat). No external test set.",
      performance="not extracted here",
      applied_to="not applied to an external inventory",
      caveat="20 descriptors on 16 compounds. LOOCV on 16 is a weak "
             "test, but the descriptors are physicochemical rather "
             "than abstract and the training data are the paper's "
             "own, which makes this the cleanest of the binding "
             "models here."),

    m(study="Dawson et al.", year="2023", doi="10.3390/toxics11020098",
      pmid="36850973", kind="qsar-ml", data_origin="borrowed",
      held_out="cv-only", detail="read",
      endpoint="in vivo elimination half-life, as 4 bins (<12 h, "
               "<1 week, <2 months, >2 months)",
      model="random forest classifier",
      trained_on="NOT ITS OWN. In vivo half-lives curated by Lau et "
                 "al. in four reviews (2007, 2012, 2015, 2021), "
                 "across human, monkey, rat and mouse. The eleven "
                 "PFAS are not enumerated in the sections read here.",
      n_train="11", n_compounds="11 PFAS, 4 species",
      descriptors="22 physico-chemical (18 predicted by OPERA v2.7, a "
                  "ToxPrint binary for a backbone ether bond, average "
                  "molecular mass, 2 chain-length) + 88 "
                  "transport/re-uptake analogue descriptors "
                  "(structural similarity of defluorinated PFAS to "
                  "endogenous chemicals, as a surrogate for "
                  "transporter affinity) + physiological descriptors "
                  "including kidney geometry, as a surrogate for "
                  "transporter expression",
      n_feat="110+", n_features="110+",
      validation="y-randomised null model as comparator; "
                 "applicability domain by the method of Roy, Kar and "
                 "Ambure. No external test set.",
      performance="86.1 % accuracy against 32.2 % for the y-randomised "
                  "null",
      applied_to="CCD PFAS master list: 6603 of 8163 had enough "
                 "information, 4136 inside the domain across four "
                 "species, 3890 for human",
      overlaps=[("corpus", "in vivo half-life",
                 "trains on the same primary half-life literature "
                 "this project's half-life tables draw on; the "
                 "overlap is at the level of the underlying "
                 "measurements, not a single named study")],
      caveat="ELEVEN compounds in training, predicting 3890. The "
             "descriptor design is the most mechanistically motivated "
             "of any model here -- it builds explicit surrogates for "
             "renal transporter expression and affinity rather than "
             "searching a descriptor pool -- and the y-randomised "
             "comparator is the right control. The training set is "
             "still eleven compounds."),

    m(study="Gao et al.", year="2025", doi="10.1021/acs.est.5c08779",
      pmid="41344664", kind="mechanistic-model", data_origin="own",
      held_out="none", detail="read",
      endpoint="apparent BSA binding constant as a function of ionic "
               "strength and salt identity",
      model="mechanistic competitive-salt binding model fitted to the "
            "paper's own fluorescence data",
      trained_on="ITS OWN measurements: Hill-equation Ka across a "
                 "range of salt conditions",
      n_compounds="not extracted here",
      descriptors="n/a -- mechanistic parameters (Kas1, Kas2, KA1, "
                  "KA2), not structural descriptors",
      n_feat="4 fitted", n_features="4 fitted parameters",
      validation="parameters fitted simultaneously across multiple "
                 "salt conditions because a single condition does not "
                 "give a unique fit; uncertainties from the "
                 "covariance matrix",
      performance="not extracted here",
      applied_to="prediction of Ka at unmeasured salt conditions",
      caveat="Not a QSAR -- no structure enters it. Listed because "
             "its fitted constants are model output, and because the "
             "paper states plainly that a single salt condition does "
             "not identify the parameters, which is the kind of "
             "admission most of the models here do not make."),

    m(study="Wu et al.", year="2024", doi="10.3390/toxics12010046",
      pmid="38251003", kind="docking", data_origin="none",
      held_out="n/a", detail="read",
      endpoint="HSA binding energy, kcal/mol",
      model="molecular docking, AutoDock Vina, PDB 4E99",
      trained_on="nothing -- docking is not fitted to data",
      n_compounds="10", descriptors="n/a (structure and scoring "
                                    "function)",
      validation="none reported",
      caveat="Rises monotonically with chain length (PFPrS -6.2, PFBS "
             "-6.4, PFHxS -7.5, PFOS -8.3), the same pattern as Li "
             "2021's docking columns, where that paper's own "
             "measurements turn over above PFHpA. Two independent "
             "instances of structure-only scoring missing the turn. "
             "The paper's experimental constants are held out "
             "separately -- see audit_wu2024_hsa.py."),

    m(study="Li et al.", year="2021", doi="10.3390/toxics9030063",
      pmid="33803062", kind="docking", data_origin="none",
      held_out="n/a", detail="db-record",
      endpoint="HSA binding, geometric-mean log Ka from docking",
      model="molecular docking against PDB 1E7G and 1AO6",
      trained_on="nothing -- docking is not fitted to data",
      n_compounds="26 docked",
      descriptors="n/a (structure and scoring function)",
      caveat="The two docking columns of SI sheet S7.2 are "
             "DELIBERATELY NOT INGESTED as values -- see "
             "add_li2021_hsa_si.py. They run up to 7.5 where the "
             "same sheet's measurements top out near 5.5, and they "
             "rise monotonically where the measurements turn over "
             "above PFHpA. That contrast is the reason the paper is "
             "worth having."),

    m(study="Han et al.", year="2021", doi="10.1021/acs.est.1c00509",
      pmid="34133149", kind="docking", data_origin="none",
      held_out="n/a", detail="db-record",
      endpoint="hL-FABP binding energy",
      model="AutoDock4 binding energies, reported alongside the "
            "paper's own thermal-shift (QITSA) measurements",
      trained_on="nothing -- docking is not fitted to data",
      n_compounds="13 rows",
      descriptors="n/a (structure and scoring function)",
      overlaps=[("Zhang et al.", "2013",
                 "Han's quoted hL-FABP constants are other papers' "
                 "values -- its PFBS 1034 uM is Zhang 2013's")],
      caveat="Separate from the docking: Han's QUOTED hL-FABP "
             "constants are held out of the L-FABP figure, because "
             "its own notes say they were taken from other papers "
             "rather than measured. Its own measurements are "
             "thermal-shift dTm, not a dissociation constant."),

    m(study="Zhao et al.", year="2023",
      doi="10.1016/j.ecoenv.2023.115310", kind="qsar-ml",
      data_origin="unknown", held_out="", detail="db-record",
      endpoint="L-FABP binding",
      model="molecular docking (Glide, Schrodinger Maestro) plus a "
            "machine-learning QSAR built on the docking scores",
      trained_on="NOT ESTABLISHED -- full text not reachable. The "
                 "method string says the ML is built ON THE DOCKING "
                 "SCORES, which would make it a model of a model "
                 "rather than of any measurement, but that is not "
                 "confirmed against the paper",
      n_compounds="not established",
      descriptors="not established", validation="not established",
      caveat="No PMC record and the publisher gates the page, so "
             "neither full text nor abstract was reachable. NOTE "
             "there are TWO Zhao 2023 papers in this database and "
             "they are easy to confuse: this one, in Ecotoxicology "
             "and Environmental Safety, carries the L-FABP docking "
             "and the ML QSAR; 10.1016/j.envint.2023.107951 in "
             "Environment International is a review with no primary "
             "data. Needs institutional access."),

    m(study="Zhao et al.", year="2023",
      doi="10.1016/j.envint.2023.107951", pmid="37126916",
      kind="read-across", data_origin="borrowed", held_out="n/a",
      detail="db-record",
      endpoint="narrative review of PFAS binding across serum "
               "albumin, liver FABP and organic anion transporters",
      model="literature review; no primary measurements and no fitted "
            "model",
      trained_on="the published literature, not re-analysed",
      caveat="A review, listed so it is not mistaken for the other "
             "Zhao 2023 above, and so its single database row is not "
             "read as a measurement."),

    m(study="Salvalaglio et al.", year="2010", doi="10.1021/jp106584b",
      pmid="21028884", kind="docking", data_origin="none",
      held_out="n/a", detail="db-record",
      endpoint="predicted HSA binding and binding rate in blood",
      model="molecular modelling",
      trained_on="nothing fitted", n_compounds="2 (PFOA, PFOS)",
      descriptors="n/a",
      caveat="Held only as a SECONDARY transcription from OEHHA 2024 "
             "PFOA/PFOS PHG Table A6.2; the primary paper has not "
             "been retrieved and no assay or modelling conditions are "
             "available. Its predicted f_unbound values (0.05, 0.02) "
             "are predictions, not measurements."),

    m(study="Qin et al.", year="2024", doi="10.1021/acs.est.3c09950",
      pmid="38394616", kind="empirical-regression", data_origin="own",
      held_out="none", detail="db-record",
      endpoint="liposome-water partition coefficient, log D_lip/w",
      model="one-variable linear regression: log D_lip/w = 0.535 N + "
            "0.658 on the number of fluorinated carbons (eq S10)",
      trained_on="ITS OWN measured log D_lip/w values",
      n_feat="1", n_features="1 (number of fluorinated carbons)",
      descriptors="the number of fluorinated carbons -- one physical "
                  "variable, not a descriptor pool",
      performance="R2 0.977",
      caveat="Only the 2 PREDICTED rows are computed; Qin's other 94 "
             "rows in this database are measurements and are not "
             "affected. One mechanistically plain variable and a "
             "stated R2 -- the opposite end of the spectrum from a "
             "753-descriptor search, and worth holding up as such."),

    m(study="Ng & Hungerbuhler", year="2013", doi="10.1021/es400981a",
      pmid="23734664", kind="mechanistic-model", data_origin="borrowed",
      held_out="none", detail="db-record",
      endpoint="bioconcentration, via coupled protein binding",
      model="mechanistic protein-binding bioconcentration model "
            "coupling albumin in the circulatory and extracellular "
            "space, FABP in liver, and OAT-facilitated renal "
            "elimination and reabsorption",
      trained_on="parameterised from published in vitro binding data",
      caveat="A process model, not a QSAR: no structural descriptor "
             "enters it. Its parameters come from the same in vitro "
             "binding literature this database holds, so it is not "
             "independent of those measurements."),

    m(study="Ng & Hungerbuhler", year="2014", doi="10.1021/es404008g",
      pmid="24762048", kind="mechanistic-model", data_origin="borrowed",
      held_out="none", detail="db-record",
      endpoint="comparison of two competing mechanisms for tissue "
               "distribution",
      model="review and comparison of phospholipid partitioning "
            "against protein interaction as competing models",
      trained_on="the published literature",
      caveat="Concludes protein interactions are needed to explain "
             "tissue distribution, the chain-length relationship and "
             "species- and sex-specific half-lives, and that the two "
             "mechanisms are not mutually exclusive."),

    m(study="Khazaee & Ng", year="2018", doi="10.1039/c7em00474e",
      pmid="29265128", kind="pbtk", data_origin="borrowed",
      held_out="none", detail="db-record",
      endpoint="zebrafish PBPK, including hepatobiliary circulation",
      model="PBPK model adaptation plus a systematic literature "
            "review of which binding parameters are available at all",
      trained_on="published binding parameters",
      caveat="The review half is the useful half here: it documents "
             "which protein-binding parameters are and are not "
             "available for PBPK, which is the same gap this project "
             "keeps running into."),

    m(study="Fischer et al.", year="2025",
      doi="10.1021/acs.est.5c05473", kind="pbtk",
      data_origin="borrowed", held_out="none", detail="db-record",
      endpoint="PBTK distribution and elimination",
      model="physiologically based toxicokinetic model parameterised "
            "from in vitro binding, permeability and transporter data",
      trained_on="published in vitro binding, permeability and "
                 "transporter data",
      performance="sensitivity analysis: permeability and "
                  "phospholipid binding dominate elimination and "
                  "distribution of long-chain PFAA (nPFC >= 7), while "
                  "elimination of short-chain PFAA (nPFC <= 6) is more "
                  "sensitive to other terms",
      caveat="Full text on disk at literature/fulltext/Fischer "
             "2025.txt; the sensitivity result above is from the "
             "abstract and sensitivity analysis, and the model's "
             "parameter provenance has not been read in detail."),

    m(study="Ebert et al.", year="2020", doi="10.1021/acs.est.0c00175",
      pmid="32212724", kind="quantum-chemical", data_origin="none",
      held_out="n/a", detail="db-record",
      endpoint="membrane-water partitioning and ionic permeability",
      model="COSMO-RS prediction alongside dialysis (K_mem/w) and "
            "planar lipid bilayer permeability (P_ion) measurements",
      trained_on="nothing -- COSMO-RS is first-principles",
      caveat="The measurements are experimental and are kept; only "
             "the COSMO-RS strand is computed. A first-principles "
             "prediction with no training set is a different kind of "
             "claim from a 753-descriptor fit, and is tagged "
             "separately for that reason."),

    m(study="Ren et al.", year="2015", doi="10.1007/s00204-014-1258-y",
      pmid="24819616", kind="docking", data_origin="none",
      held_out="n/a", detail="db-record",
      endpoint="thyroid hormone receptor binding",
      model="molecular docking alongside a fluorescence competitive "
            "binding assay and TH-responsive cell proliferation",
      trained_on="nothing fitted",
      caveat="Experimental rows are kept; the docking component is "
             "the computed strand."),

    m(study="Ishibashi et al.", year="2019",
      doi="10.1021/acs.est.8b07273", pmid="30649875",
      kind="homology-model", data_origin="none", held_out="n/a",
      detail="db-record",
      endpoint="PPARalpha ligand-binding domain interaction",
      model="in silico homology modelling alongside an in vitro "
            "competitive binding assay with recombinant LBD protein",
      trained_on="nothing fitted",
      caveat="A modelled protein STRUCTURE rather than a predicted "
             "affinity. Experimental rows are kept."),

    m(study="Wilson et al.", year="2026",
      doi="10.1016/j.envres.2026.124118", pmid="41748002",
      kind="docking", data_origin="none", held_out="n/a",
      detail="db-record",
      endpoint="beta-lactoglobulin binding",
      model="molecular docking and molecular dynamics alongside "
            "spectroscopy and retinol-binding kinetics",
      trained_on="nothing fitted",
      caveat="Experimental rows are kept; docking and MD are the "
             "computed strand. Beta-lactoglobulin is a milk protein, "
             "not a plasma transport protein."),
]


def load_db():
    with open(DB, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def overlap(rows):
    """Resolve each model's declared training sources against the
    database, so a model and the measurements it learned from are
    never counted as two pieces of evidence."""
    print("=" * 78)
    print("DOUBLE-COUNTING CHECK")
    print("=" * 78)
    print(wrap("A model trained on somebody else's measurements is not "
               "independent evidence about those measurements. Where "
               "this database also holds them, counting both is "
               "counting once twice.", 2))
    hits = 0
    for mm in M:
        if not mm["overlaps"]:
            continue
        for study, year, why in mm["overlaps"]:
            print(f"\n  {mm['study']} {mm['year']}  [{mm['kind']}, "
                  f"{mm['data_origin']}]")
            if study == "corpus":
                n = 0
                for path in HALFLIFE:
                    if os.path.exists(path):
                        with open(path, encoding="utf-8") as fh:
                            n += sum(1 for _ in csv.DictReader(fh))
                print(f"      overlaps a CORPUS, not one study: {year}")
                print(wrap(why, 6))
                print(f"      this project holds {n} rows of that kind "
                      f"across {len(HALFLIFE)} half-life tables")
                hits += 1
                continue
            match = [r for r in rows
                     if r["study"] == study and r["year"] == year]
            exp = [r for r in match
                   if r["is_computational"].strip().lower()
                   .startswith("no")]
            cmp_ = sorted({r["chemical"] for r in exp})
            print(f"      trained on: {study} {year}")
            print(wrap(why, 6))
            if not match:
                print("      NOT IN THIS DATABASE -- no exposure here")
                continue
            hits += 1
            print(f"      this database holds {len(exp)} experimental "
                  f"rows from it, {len(cmp_)} compounds")
            print(wrap("compounds: " + ", ".join(cmp_), 6))
            print("      => NOT INDEPENDENT. Do not count the model "
                  "and these rows")
            print("         as separate support for the same claim.")
    print()
    print(f"  {hits} declared overlap(s) resolved.")
    print(wrap("This check only sees overlaps a model entry DECLARES. "
               "Adding a model means stating what it trained on; an "
               "undeclared training set is invisible here, which is "
               "the limit of the method, not a clean bill.", 2))


def wrap(s, i):
    return textwrap.fill(s, 76, initial_indent=" " * i,
                         subsequent_indent=" " * i)


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for mm in M:
            w.writerow({**mm, "overlaps_db": "; ".join(
                f"{s} {y}" for s, y, _ in mm["overlaps"])})
    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}  "
          f"({len(M)} entries)\n")

    print("=" * 78)
    print("TAGS")
    print("=" * 78)
    for axis in ("kind", "data_origin", "held_out", "detail"):
        tally = {}
        for mm in M:
            tally[mm[axis] or "(unset)"] = tally.get(
                mm[axis] or "(unset)", 0) + 1
        line = ", ".join(f"{k} {v}" for k, v in
                         sorted(tally.items(), key=lambda kv: -kv[1]))
        print(f"\n  {axis}")
        print(wrap(line, 6))

    print()
    print("=" * 78)
    print("THE STATISTICAL MODELS, SIDE BY SIDE")
    print("=" * 78)
    print(f"\n  {'study':22}{'origin':10}{'compounds':>19}"
          f"{'features':>13}  held out")
    for mm in M:
        if not mm["kind"].startswith("qsar") or not mm["n_feat"]:
            continue
        print(f"  {mm['study'].replace(' et al.', '') + ' ' + mm['year']:22}"
              f"{mm['data_origin']:10}{mm['n_train']:>19}"
              f"{mm['n_feat']:>13}  {mm['held_out']}")

    print()
    overlap(load_db())

    print()
    print("=" * 78)
    print("EVERY ENTRY")
    print("=" * 78)
    for mm in sorted(M, key=lambda x: (x["kind"], x["study"])):
        print(f"\n  [{mm['kind']} | {mm['data_origin']} | "
              f"{mm['detail']}]  {mm['study']} {mm['year']}")
        print(wrap(mm["caveat"], 6))


if __name__ == "__main__":
    main()
