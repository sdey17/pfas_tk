#!/usr/bin/env python3
"""The in silico PFAS-binding literature, tracked on its own terms.

WHY THIS IS A SEPARATE FILE.  db/protein_binding.csv holds
measurements, and db/combined/ holds measurements.  A model is not a
measurement and does not belong in a Kd column beside one, but it is
not nothing either -- the QSAR aim of this review is precisely about
what structure can and cannot predict.  So the models are kept here,
with the fields that actually decide whether a model means anything
and which a binding table has nowhere to put:

    WHAT IT WAS TRAINED ON.  A model trained on somebody else's
    measurements is not independent evidence about those
    measurements, and two models trained on the same dataset are one
    piece of evidence, not two.  This is not hypothetical here:
    Gallagher 2023's QSAR is trained on the Starnes DSF dataset,
    which this database already holds as experimental rows.

    WHAT DESCRIPTORS IT USED, and from which software.  "Chain
    length" and "Moran autocorrelation of lag 8 weighted by mass" are
    both descriptors; only one of them tells you anything about a
    mechanism.

    HOW MANY COMPOUNDS AGAINST HOW MANY FEATURES.  Every model below
    has more features available than compounds, most by two orders of
    magnitude.  That is not a reason to dismiss them, but it is the
    first thing to know about them.

    WHETHER ANYTHING WAS HELD OUT.  Cross-validation on 16 compounds
    and an external test set are not the same claim.

Rows marked `partial` are papers whose computational strand is real
but which have not been read for these fields yet.  They are listed
rather than omitted, because a silently missing study looks exactly
like a study that does not exist.

Sources: each paper as cited.  Gallagher 2023 full text via PubMed
Central (PMC10383382).  Ryu 2024, Starnes 2025, Dawson 2023 and
Gao 2025 from the copies in papers/.

Run:  python3 tk_review/scripts/insilico_inventory.py
Writes db/insilico_models.csv
"""
import csv
import os
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "db", "insilico_models.csv")

FIELDS = ["study", "year", "doi", "pmid", "status", "endpoint",
          "model", "trained_on", "n_train", "n_compounds",
          "descriptors", "n_feat", "n_features", "validation",
          "performance", "applied_to", "caveat"]

# n_train and n_feat are the bare numbers for the summary table;
# n_compounds and n_features carry the sentence that qualifies them.
COMPACT = {
    "Starnes et al.2025": ("57", "255"),
    "Gallagher et al.2023": ("16 train / 8 test", "753 -> 3-4"),
    "Ryu et al.2024": ("16", "20"),
    "Dawson et al.2023": ("11", "110+"),
    "Gao et al.2025": ("", "4 fitted"),
    "Wu et al.2024": ("", ""),
    "Zhao et al.2023": ("", ""),
}

M = [
    dict(
        study="Starnes et al.", year="2025",
        doi="10.1093/toxsci/kfae124", pmid="39298512", status="read",
        endpoint="HSA binding affinity (Ka) and binary bind/no-bind, "
                 "from its own differential scanning fluorimetry assay",
        model="7 classifiers and 9 regressors (python/sklearn); best "
              "classifier logistic regression, best regressor support "
              "vector regression",
        trained_on="ITS OWN measurements: 65 PFAS from 20 chemical "
                   "classes plus 3 fatty acids as non-fluorinated "
                   "analogues; 57 chemicals reached the training set "
                   "(54 PFAS + 3 fatty acids), one congener dropped "
                   "for missing descriptors",
        n_compounds="57",
        descriptors="Saagar (834) + ToxPrints (729) + OPERA (12) = "
                    "1575 generated; invariant features removed and "
                    "min-max scaled to 255 used",
        n_features="255",
        validation="3-by-3-fold cross-validation. No external test "
                   "set.",
        performance="classifier AUROC 0.936; regressor R2 0.854",
        applied_to="EPAPFASINV inventory, 430 PFAS",
        caveat="255 features on 57 chemicals, about 4.5 features per "
               "compound, and the only held-out data is the "
               "cross-validation fold. The DSF endpoint itself fails "
               "the physiological bound in "
               "binding_kd_physiological_bound.py (81 % bound at six "
               "sites, where plasma needs >99 %)."),
    dict(
        study="Gallagher et al.", year="2023",
        doi="10.3390/molecules28145375", pmid="37513249", status="read",
        endpoint="HSA binding affinity (EC50, mM) and a high/low "
                 "class at a 1.45 mM threshold",
        model="GA-LDA classification (DTC Lab "
              "ClassificationBasedQSAR v1.0.0); GA-MLR and PLS "
              "regression (SmallDataModeler v1.0); GenRA read-across "
              "on ToxPrint chemotypes; Glide XP docking on PDB 7AAI",
        trained_on="NOT ITS OWN. 24 PFAS EC50 values collected from "
                   "the literature; they match the Starnes DSF "
                   "dataset already held in db/protein_binding.csv, "
                   "so this is a model OF that dataset, not "
                   "independent evidence about it",
        n_compounds="24 classification / 22 regression (two "
                    "fluorotelomer alcohols did not bind)",
        descriptors="alvaDesc 2.0.16, 753 2D descriptors for "
                    "classification and 736 for regression. Final "
                    "classification model keeps 3: average "
                    "connectivity index of order 4, an augmented "
                    "edge-adjacency eigenvalue, an eccentric "
                    "topological index. Final regression model keeps "
                    "4: packing density index, Moran autocorrelation "
                    "lag 8 (mass), Geary autocorrelation lag 8 (vdW "
                    "volume), QED drug-likeness",
        n_features="753 offered, 3-4 retained",
        validation="classification split 2:1 at random, training "
                   "n = 16 and test n = 8; regression used the whole "
                   "set under exhaustive double cross-validation, no "
                   "split. Docking protocol validated by redocking "
                   "PFOA, RMSD 0.092 A.",
        performance="ROC AUC 0.97 training, 1.00 on the 8-compound "
                    "test set",
        applied_to="not applied to an external inventory",
        caveat="753 descriptors searched against 16 training "
               "compounds, and a test AUC of 1.00 on eight. The "
               "retained descriptors are abstract topological indices "
               "with no stated mechanism. The authors say themselves "
               "that the GenRA predictions are 'not that reliable': "
               "every AUC was 0 and p ran 0.73 to 1."),
    dict(
        study="Ryu et al.", year="2024",
        doi="10.1021/acs.est.4c04050", pmid="39161261", status="read",
        endpoint="unbound fraction (f_unbound) in liver, brain, "
                 "heart, kidney and lung",
        model="Partial least squares regression, four separate "
              "models (human, C57, mouse, rat)",
        trained_on="ITS OWN measurements: 16 PFAS across four "
                   "organisms, the 320 experimental rows this "
                   "database holds",
        n_compounds="16",
        descriptors="20 physicochemical properties from MOE 2018.01 "
                    "and Dragon 6.0. Molecular weight and log D are "
                    "the ones the paper reads a trend off.",
        n_features="20",
        validation="leave-one-out cross-validation over the 16 PFAS, "
                   "grid search for the component count (4 for human, "
                   "3 for C57, mouse and rat). No external test set.",
        performance="not extracted here",
        applied_to="not applied to an external inventory",
        caveat="20 descriptors on 16 compounds. LOOCV on 16 is a weak "
               "test, but the descriptors are physicochemical rather "
               "than abstract, and the training data are the paper's "
               "own, which makes this the cleanest of the binding "
               "models here."),
    dict(
        study="Dawson et al.", year="2023",
        doi="10.3390/toxics11020098", pmid="36850973", status="read",
        endpoint="in vivo elimination half-life, as 4 bins "
                 "(<12 h, <1 week, <2 months, >2 months)",
        model="random forest classifier",
        trained_on="NOT ITS OWN. In vivo half-lives curated by Lau "
                   "et al. in four reviews (2007, 2012, 2015, 2021), "
                   "across human, monkey, rat and mouse",
        n_compounds="11 PFAS",
        descriptors="22 physico-chemical (18 predicted by OPERA "
                    "v2.7, a ToxPrint binary for an ether bond in "
                    "the backbone, average molecular mass, 2 chain "
                    "length) + 88 transport/re-uptake analogue "
                    "descriptors (structural similarity of "
                    "defluorinated PFAS to endogenous chemicals, as "
                    "a surrogate for transporter affinity) + "
                    "physiological descriptors including kidney "
                    "geometry, as a surrogate for transporter "
                    "expression",
        n_features="110+",
        validation="y-randomised null model as the comparator; "
                   "applicability domain by the method of Roy, Kar "
                   "and Ambure",
        performance="86.1 % accuracy against 32.2 % for the "
                    "y-randomised null",
        applied_to="CCD PFAS master list: 6603 of 8163 had enough "
                   "information, 4136 fell inside the domain across "
                   "four species, 3890 for human",
        caveat="ELEVEN compounds in training, predicting 3890. The "
               "descriptor design is the most mechanistically "
               "motivated of any model here -- it builds explicit "
               "surrogates for renal transporter expression and "
               "affinity rather than searching a descriptor pool -- "
               "and the y-randomised comparator is the right control. "
               "The training set is still eleven compounds."),
    dict(
        study="Wu et al.", year="2024",
        doi="10.3390/toxics12010046", pmid="38251003", status="read",
        endpoint="HSA binding energy, kcal/mol",
        model="molecular docking, AutoDock Vina, PDB 4E99",
        trained_on="nothing -- docking is not fitted to data",
        n_compounds="10",
        descriptors="n/a (structure and scoring function)",
        n_features="n/a",
        validation="none reported",
        performance="n/a",
        applied_to="n/a",
        caveat="Rises monotonically with chain length (PFPrS -6.2, "
               "PFBS -6.4, PFHxS -7.5, PFOS -8.3), the same pattern "
               "as Li 2021's docking columns, where that paper's own "
               "measurements turn over above PFHpA. Two independent "
               "instances of structure-only scoring missing the "
               "turnover. The paper's experimental constants are "
               "held out separately -- see audit_wu2024_hsa.py."),
    dict(
        study="Gao et al.", year="2025",
        doi="10.1021/acs.est.5c08779", pmid="41344664", status="read",
        endpoint="apparent BSA binding constant as a function of "
                 "ionic strength and salt identity",
        model="mechanistic competitive-salt binding model, fitted "
              "to the paper's own fluorescence data",
        trained_on="ITS OWN measurements: Hill-equation Ka across a "
                   "range of salt conditions",
        n_compounds="not extracted here",
        descriptors="n/a -- mechanistic parameters (Kas1, Kas2, KA1, "
                    "KA2), not structural descriptors",
        n_features="4 fitted parameters",
        validation="parameters fitted simultaneously across multiple "
                   "salt conditions because a single condition does "
                   "not give a unique fit; uncertainties from the "
                   "covariance matrix",
        performance="not extracted here",
        applied_to="prediction of Ka at unmeasured salt conditions",
        caveat="Not a QSAR -- no structure enters it. Listed because "
               "its fitted constants are model output, and because "
               "the paper states plainly that a single salt "
               "condition does not identify the parameters, which is "
               "the kind of admission most of the models here do not "
               "make."),
    dict(
        study="Zhao et al.", year="2023",
        doi="10.1016/j.ecoenv.2023.115310", pmid="", status="UNRESOLVED",
        endpoint="L-FABP binding; a review with a docking and "
                 "machine-learning strand",
        model="molecular docking (Glide, Schrodinger Maestro) plus a "
              "machine-learning QSAR built on the docking scores",
        trained_on="NOT ESTABLISHED -- full text not reachable. The "
                   "method string says the ML is built ON THE "
                   "DOCKING SCORES, which would make it a model of a "
                   "model rather than of any measurement, but that "
                   "has not been confirmed against the paper",
        n_compounds="not established",
        descriptors="not established",
        n_features="not established",
        validation="not established",
        performance="not established",
        applied_to="not established",
        caveat="No PMC record and the publisher gates the page, so "
               "neither full text nor abstract was reachable. Note "
               "there are TWO Zhao 2023 papers in this database and "
               "they are easy to confuse: this one, in Ecotoxicology "
               "and Environmental Safety, carries the L-FABP docking "
               "and the ML QSAR; 10.1016/j.envint.2023.107951 in "
               "Environment International is a review with no "
               "primary data. Needs institutional access."),
]

# Papers with a real computational strand that have NOT been read for
# these fields. Listed so that a gap looks like a gap.
PARTIAL = [
    ("Han et al.", "2021", "AutoDock4 binding energies alongside "
     "experimental work; 13 rows"),
    ("Li et al.", "2021", "two docking columns in SI sheet S7.2, "
     "geometric-mean log Ka against PDB 1E7G and 1AO6; deliberately "
     "NOT ingested, see add_li2021_hsa_si.py. They rise monotonically "
     "where the same sheet's measurements turn over"),
    ("Ng & Hungerbuhler", "2013", "mechanistic protein-binding "
     "bioconcentration model, parameterised from in vitro data"),
    ("Ng & Hungerbuhler", "2014", "comparison of phospholipid "
     "partitioning against protein interaction as competing models"),
    ("Khazaee & Ng", "2018", "PBPK model plus a systematic review of "
     "which binding parameters are available at all"),
    ("Fischer et al.", "2025", "PBTK model parameterised from in "
     "vitro binding, permeability and transporter data"),
    ("Qin et al.", "2024", "one-descriptor regression: log D_lip/w = "
     "0.535 N + 0.658 on the number of fluorinated carbons, "
     "R2 0.977 (eq S10). Two predicted membrane values"),
    ("Salvalaglio et al.", "2010", "molecular modelling, predicted "
     "binding rate in blood; held only as a secondary transcription "
     "from OEHHA 2024 Table A6.2, primary not retrieved"),
    ("Ebert et al.", "2020", "COSMO-RS prediction alongside "
     "experiment"),
    ("Ren et al.", "2015", "docking component"),
    ("Ishibashi et al.", "2019", "homology modelling component"),
    ("Wilson et al.", "2026", "docking and MD component"),
]


def main():
    for m in M:
        n, f = COMPACT.get(m["study"] + m["year"], ("", ""))
        m.setdefault("pmid", "")
        m["n_train"], m["n_feat"] = n, f
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(M)
    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}  "
          f"({len(M)} characterised, {len(PARTIAL)} listed as partial)\n")

    wrap = lambda s, i: textwrap.fill(s, 74, initial_indent=" " * i,
                                      subsequent_indent=" " * i)
    print("=" * 78)
    print("WHAT THE MODELS WERE TRAINED ON")
    print("=" * 78)
    for m in M:
        if m["model"].startswith("molecular docking"):
            continue
        own = "OWN DATA " if m["trained_on"].startswith("ITS OWN") else \
              ("UNKNOWN  " if m["status"] == "UNRESOLVED" else "BORROWED ")
        print(f"\n{own}{m['study']} {m['year']}"
              + (f"  [n = {m['n_train']}]" if m["n_train"] else ""))
        print(wrap(m["trained_on"], 4))
    print()
    print("=" * 78)
    print("DESCRIPTORS, AND HOW MANY PER COMPOUND")
    print("=" * 78)
    print(f"\n  {'study':22}{'compounds':>19}{'features':>13}   software")
    for m in M:
        if not m["n_feat"]:
            continue
        soft = ("alvaDesc" if "alvaDesc" in m["descriptors"] else
                "MOE + Dragon" if "MOE" in m["descriptors"] else
                "Saagar + ToxPrint + OPERA"
                if "Saagar" in m["descriptors"] else
                "OPERA + ToxPrint" if "OPERA" in m["descriptors"] else
                "none (mechanistic)")
        print(f"  {m['study'].replace(' et al.','') + ' ' + m['year']:22}"
              f"{m['n_train']:>19}{m['n_feat']:>13}   {soft}")
    print()
    print("=" * 78)
    print("THE SHAPE OF THE EVIDENCE")
    print("=" * 78)
    for m in M:
        print(f"\n{m['study']} {m['year']}")
        print(wrap(m["caveat"], 4))
    print()
    print("=" * 78)
    print(f"COMPUTATIONAL STRAND PRESENT, NOT CHARACTERISED HERE "
          f"({len(PARTIAL)})")
    print("=" * 78)
    for s, y, what in PARTIAL:
        print(f"\n  {s} {y}")
        print(wrap(what, 6))


if __name__ == "__main__":
    main()
