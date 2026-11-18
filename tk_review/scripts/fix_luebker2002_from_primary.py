#!/usr/bin/env python3
"""Luebker 2002 from the primary. It reports IC50 only, and its abstract has a
thousandfold unit error.

Luebker DJ, Hansen KJ, Bass NM, Butenhoff JL, Seacat AM. "Interactions of
flurochemicals with rat liver fatty acid-binding protein." Toxicology
2002;176(3):175-185. doi:10.1016/S0300-483X(02)00081-1  PMID 12093614.

WHAT WAS THERE BEFORE. Two rows, both from the ABSTRACT: the DAUDA probe at
"0.47 nM", and one row naming five compounds at once with a rank order and no
numbers.

FINDING 1: THE ABSTRACT'S PROBE CONSTANT IS WRONG BY A FACTOR OF 1000, and we
had transcribed it. The abstract says "The dissociation constant (Kd) of the
DAUDA-L-FABP complex was 0.47 nM". Results section 3.3 says "The average Kd
for DAUDA-L-FABP binding was 0.35 +/- 0.02 uM (n = 3)". The Discussion then
anchors that against the literature: "equivalent to the Kds of 0.46 +/- 0.04
uM (Thumser et al., 1994a) and 0.38 +/- 0.02 uM (Thumser et al., 1996)
reported previously for recombinant wild type rat L-FABP".

So the body says 0.35 uM three times over and the abstract says 0.47 nM once.
MICROMOLAR is right. Stored as 0.35 uM with the discrepancy recorded. It is a
probe, so it never reached the working table, but a probe constant 1000x too
tight would corrupt any competitive conversion built on it.

FINDING 2: THE PAPER REPORTS NO DISSOCIATION CONSTANT FOR ANY PFAS. Table 1 is
IC50 values and nothing else:

    Oleic acid    0.01 uM        (endogenous ligand, not an analyte)
    PFOS          4.9  uM
    N-EtFOSA      9.7  uM
    WY-14643     10    uM        (non-fluorochemical comparator, not a PFAS)
    N-EtFOSE    > 10   uM
    PFOA        > 10   uM

Two of the four fluorochemicals have no finite IC50 at all. Under this
compilation's rule -- a source earns a row by reporting K, or Kd with n, not
by reporting that binding happened -- Luebker contributes NO rows to the
working table. The IC50 values are stored as the record, with kd_units marked
so the builder drops them for the stated reason rather than silently.

Converting these to Ki with Cheng-Prusoff would be the reader's arithmetic,
not the paper's, and would import the DAUDA constant into a quantity the
authors declined to publish. Not done. If it ever is, the inputs are here:
[DAUDA] = 1 uM and Kd(DAUDA) = 0.35 uM, giving Ki = IC50 / 3.857.

FINDING 3: THIS ASSAY PROBES THE SECONDARY SITE ONLY, and the authors say so.
"DAUDA is believed to bind exclusively to the secondary binding site on
L-FABP." And on their own competition data: "The shifting of wavelength and
decrease in FI seen in the competition assays are indicative of DAUDA being
displaced from the non-polar secondary binding site. Thus, secondary site
binding of the fluorochemicals and WY is supported by these data. However, to
investigate binding of fluorochemicals to the primary site and draw
conclusions about in vivo binding, a probe known to bind both sites, i.e.
radiolabeled oleic acid, would need to be used."

Woodcroft 2010 cites this paper for exactly that limitation and switches to
ANS for it. The probe therefore partly determines WHICH SITE a source is
measuring, which is a structural contributor to the spread in this
compilation and not simply noise:

    DAUDA            Luebker 2002                    secondary site only
    1,8-ANS / ANS    Woodcroft 2010 (2 sites, 2:1 by Scatchard),
                     Sheng 2016 (n 1.01), Sheng 2018 (n 0.93),
                     Yang 2020, Fan 2025
    NBD-lauric acid  Zhang 2013                      primary site

FINDING 4: THE SOLVENT CORRECTION IS LARGE. Competitors were added from 1 mM
stocks in 100% methanol, 0.4-20 ul into 2 ml, so up to 1% methanol at the top
point, and Fig. 5's caption states the methanol effect "ranged from 2%
inhibition at 1 uM to 21% inhibition at 10 uM". A fifth of the signal at the
top concentration is solvent, subtracted by correction.

PROTEIN. Native rat L-FABP purified from the livers of untreated male rats --
not recombinant, which makes it the only such source here. Measured mass
14,315 Da by LC/ESMSMS, with a presumed dimer at 28 kDa also present.

CONDITIONS. 1 uM L-FABP with 1 uM DAUDA in 2 ml of 50 mM KH2PO4 pH 7.2, room
temperature; competitor 0.2-10 uM final.

Idempotent.
"""
import csv, os, shutil, sys

DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "db",
                  "protein_binding.csv")
DOI = "10.1016/s0300-483x(02)00081-1"

PROTEIN = ("liver fatty acid-binding protein (rL-FABP), NATIVE, purified from "
           "the livers of untreated male Crl:CD rats; measured mass 14,315 Da "
           "by LC/ESMSMS, with a presumed dimer at 28 kDa also present")
M_COMP = ("competitive displacement of DAUDA (11-(5-dimethylaminonaphthalene"
          "sulphonyl)-undecanoic acid) monitored by fluorescence emission; "
          "IC50 only, no dissociation constant fitted for the competitors. "
          "DAUDA binds the SECONDARY site of L-FABP, so this measures "
          "secondary-site competition, as the paper states")
M_DIR = ("direct fluorimetric titration of 1 uM L-FABP with DAUDA, specific "
         "binding fitted to a Michaelis-Menten hyperbola")
PROT_C = "1 uM L-FABP with 1 uM DAUDA in 2 ml"
LIG_C = ("0.2-10 uM final, from 0.4-20 ul aliquots of a 1 mM competitor stock "
         "in 100% methanol into 2 ml; methanol alone inhibited 2% at 1 uM "
         "rising to 21% at 10 uM and was corrected for")
BUF = "50 mM KH2PO4, pH 7.2, room temperature"

DAUDA_NOTE = ("Results section 3.3 gives 0.35 +/- 0.02 uM (n = 3) and the "
              "Discussion calls it equivalent to Thumser's 0.46 +/- 0.04 uM "
              "and 0.38 +/- 0.02 uM. The ABSTRACT says 0.47 nM, which is "
              "wrong by a factor of 1000 and is what this database previously "
              "held. One site; DAUDA binds the secondary site only")
IC50_NOTE = ("Table 1 reports IC50 only. The paper fits no dissociation "
             "constant for any competitor, so this row carries no K or Kd and "
             "does not enter the working table. Cheng-Prusoff inputs if ever "
             "needed: [DAUDA] 1 uM, Kd(DAUDA) 0.35 uM, Ki = IC50 / 3.857")

#  chemical, nC, head, value, units, method, note
ROWS = [
    ("DAUDA (11-(dansylamino)undecanoic acid, reference ligand)", "",
     "fluorescent fatty acid probe", "0.35", "uM", M_DIR, DAUDA_NOTE,
     "Results 3.3 | Kd 0.35 +/- 0.02 uM, n = 3, one site"),
    ("oleic acid", "", "fatty acid (endogenous ligand)",
     "0.01", "IC50 uM", M_COMP, IC50_NOTE,
     "Table 1 | IC50 0.01 uM | 91% inhibition at 2 uM"),
    ("PFOS", "8", "sulfonate", "4.9", "IC50 uM", M_COMP, IC50_NOTE,
     "Table 1 | IC50 4.9 uM | 69% inhibition at 10 uM"),
    ("N-EtFOSA", "8", "sulfonamide", "9.7", "IC50 uM", M_COMP, IC50_NOTE,
     "Table 1 | IC50 9.7 uM | 51% inhibition at 10 uM"),
    ("Wyeth-14643 (WY, non-fluorochemical peroxisome proliferator)", "",
     "non-fluorochemical comparator", "10", "IC50 uM", M_COMP, IC50_NOTE,
     "Table 1 | IC50 10 uM | 50% inhibition at 10 uM"),
    ("N-EtFOSE", "8", "sulfonamidoethanol", "", "IC50 uM", M_COMP, IC50_NOTE,
     "Table 1 | IC50 > 10 uM, no finite value | 43% inhibition at 10 uM"),
    ("PFOA", "7", "carboxylate", "", "IC50 uM", M_COMP, IC50_NOTE,
     "Table 1 | IC50 > 10 uM, no finite value | 43% inhibition at 10 uM"),
]


def build(fields):
    out = []
    for ch, nc, head, val, units, meth, note, src in ROWS:
        r = {f: "" for f in fields}
        r.update({
            "study": "Luebker et al.", "year": "2002", "pmid": "12093614",
            "doi": DOI, "chemical": ch, "n_perfluorinated_carbons": nc,
            "head_group": head, "protein": PROTEIN, "protein_species": "rat",
            "method": meth, "kd": val, "kd_units": units,
            "n_binding_sites": "1 (DAUDA; one site)" if units == "uM" else "",
            "binding_site_identified": "secondary (non-polar) binding site",
            "protein_conc": PROT_C, "protein_conc_units": "uM",
            "ligand_conc_range": LIG_C, "ligand_conc_units": "uM",
            "ligand_over_protein_max": "10", "buffer": BUF,
            "is_computational": "no",
            "tk_relevance_noted": "partly - displacement of endogenous fatty "
                                  "acids from L-FABP as a toxicity mechanism",
            "assay_volume_note": note, "source_table": src,
        })
        out.append(r)
    return out


def main():
    with open(DB, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        fields, rows = rdr.fieldnames, list(rdr)

    old = [r for r in rows if (r.get("doi") or "").lower() == DOI]
    new = build(fields)
    if old and len(old) == len(new) and all(
            a.get(k) == b.get(k) for a, b in zip(old, new) for k in fields):
        print("no change needed (Luebker 2002 already ingested from the primary)")
        return 0

    merged, placed = [], False
    for r in rows:
        if (r.get("doi") or "").lower() == DOI:
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
    print(f"Luebker 2002: replaced {len(old)} abstract-only row(s) with "
          f"{len(new)} rows from the primary")
    print("   DAUDA probe Kd 0.47 nM -> 0.35 uM (the abstract is wrong by 1000x)")
    print("   Table 1's six IC50 values stored; no PFAS Kd exists in the paper,")
    print("   so Luebker still contributes NO rows to the working table")
    return 0


if __name__ == "__main__":
    sys.exit(main())
