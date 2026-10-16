#!/usr/bin/env python3
"""Ingest Ryu et al. 2024 SI Tables S1 and S4-S7 into db/protein_binding.csv.

Rapid equilibrium dialysis, f_unbound of 16 PFAS in five whole-tissue
homogenates from four organisms (human, Wistar Han rat, CD-1 mouse,
C57BL/6 mouse), n = 4, 37 C, 5 uM, 18 h.  320 values.

doi:10.1021/acs.est.4c04050, PMID 39161261, SI es4c04050_si_001.pdf.

Two printed cells are flagged rather than corrected; see SUSPECT below.
Replaces the single study-level row written by
add_binding_studies_2026_10b.py.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")
TIDY = os.path.join(HERE, "..", "db", "ryu2024_funbound.csv")

TISSUES = ["liver", "brain", "heart", "kidney", "lung"]

# Table S1: acronym -> (total carbons, MW g/mol, h_logD, pKa, class)
# n_perfluorinated_carbons = total C - 1 for the carboxylates, = total C for
# the sulfonates, the sulfonamide and the fluorotelomer sulfonate.
S1 = {
    "PFBA":   (4,  214.04, 2.23,  0.08,  "carboxylate", 3),
    "PFPeA":  (5,  264.05, 3.02,  0.40,  "carboxylate", 4),
    "PFHxA":  (6,  314.05, 3.81,  0.70,  "carboxylate", 5),
    "PFHpA":  (7,  364.06, 4.59, -0.15,  "carboxylate", 6),
    "PFOA":   (8,  414.07, 5.38,  2.50,  "carboxylate", 7),
    "PFNA":   (9,  464.08, 6.17, -0.17,  "carboxylate", 8),
    "PFDA":   (10, 514.08, 6.96, -0.17,  "carboxylate", 9),
    "PFUnDA": (11, 564.09, 7.75, -0.17,  "carboxylate", 10),
    "PFDoDA": (12, 614.10, 8.54, -0.17,  "carboxylate", 11),
    "PFTrDA": (13, 664.10, 9.33, None,   "carboxylate", 12),
    "PFTeDA": (14, 714.11, 10.11, -0.21, "carboxylate", 13),
    "PFBS":   (4,  300.10, 3.73,  0.14,  "sulfonate", 4),
    "PFHxS":  (6,  400.12, 5.30,  0.14,  "sulfonate", 6),
    "PFOS":   (8,  500.13, 6.88,  0.14,  "sulfonate", 8),
    "PFOSA":  (8,  499.14, 5.52,  6.24,  "sulfonamide", 8),
    "6:2 FtS": (8, 427.98, 4.44,  0.36,  "fluorotelomer sulfonate", 6),
}
# The SI writes PFPA, PFUDA and PFTrDA; mapped to the names already used in
# protein_binding.csv.
ORDER = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA",
         "PFDoDA", "PFTrDA", "PFTeDA", "PFBS", "PFHxS", "PFOS", "PFOSA",
         "6:2 FtS"]

# (f_unbound, %CV) in TISSUES order, transcribed from the SI tables.
S4_HUMAN = {
    "PFBA":   [(0.301, 15), (0.323, 26), (0.261, 4), (0.393, 9), (0.393, 13)],
    "PFPeA":  [(0.360, 14), (0.233, 8), (0.309, 7), (0.306, 5), (0.468, 12)],
    "PFHxA":  [(0.329, 6), (0.12, 10), (0.238, 7), (0.241, 11), (0.257, 4)],
    "PFHpA":  [(0.1073, 9), (0.046, 20), (0.12, 6), (0.147, 11), (0.153, 10)],
    "PFOA":   [(0.0360, 7), (0.0118, 7), (0.0438, 7), (0.0565, 7), (0.0488, 14)],
    "PFNA":   [(0.00753, 7), (0.00295, 12), (0.0108, 7), (0.0145, 8), (0.0105, 3)],
    "PFDA":   [(0.00228, 4), (0.000781, 4), (0.00415, 4), (0.00332, 6), (0.00334, 5)],
    "PFUnDA": [(0.000463, 6), (0.000194, 5), (0.0007, 5), (0.000765, 2), (0.00133, 12)],
    "PFDoDA": [(0.000212, 4), (0.0000311, 15), (0.000232, 20), (0.000204, 10), (0.000393, 6)],
    "PFTrDA": [(0.0000319, 13), (0.0000122, 5), (0.0000181, 16), (0.0000216, 13), (0.0000425, 13)],
    "PFTeDA": [(0.00000068, 24), (0.00000454, 15), (0.0000130, 29), (0.00000618, 9), (0.0000294, 14)],
    "PFBS":   [(0.204, 4), (0.11, 7), (0.219, 7), (0.21, 4), (0.211, 13)],
    "PFHxS":  [(0.0528, 5), (0.0282, 6), (0.05, 6), (0.068, 12), (0.0606, 6)],
    "PFOS":   [(0.00368, 5), (0.00202, 6), (0.00387, 5), (0.0074, 20), (0.00619, 6)],
    "PFOSA":  [(0.000411, 11), (0.00048, 14), (0.00103, 12), (0.00120, 8), (0.000832, 10)],
    "6:2 FtS": [(0.0365, 6), (0.0263, 8), (0.0606, 7), (0.0801, 13), (0.0595, 7)],
}
S5_RAT = {
    "PFBA":   [(0.354, 19), (0.242, 22), (0.489, 10), (0.393, 16), (0.433, 15)],
    "PFPeA":  [(0.336, 10), (0.215, 8), (0.832, 19), (0.379, 10), (0.504, 17)],
    "PFHxA":  [(0.232, 8), (0.251, 14), (0.28, 19), (0.297, 5), (0.795, 11)],
    "PFHpA":  [(0.0909, 8), (0.0577, 6), (0.0936, 12), (0.138, 6), (0.196, 6)],
    "PFOA":   [(0.0289, 7), (0.0179, 3), (0.042, 6), (0.0357, 8), (0.0459, 6)],
    "PFNA":   [(0.0069, 1), (0.00369, 5), (0.008, 8), (0.0122, 6), (0.0128, 8)],
    "PFDA":   [(0.00188, 11), (0.000848, 15), (0.00214, 10), (0.0041, 5), (0.00348, 4)],
    "PFUnDA": [(0.00032, 9), (0.000154, 6), (0.000512, 21), (0.000644, 6), (0.000876, 9)],
    "PFDoDA": [(0.000106, 10), (0.0000271, 14), (0.000102, 10), (0.0000675, 11), (0.000143, 10)],
    "PFTrDA": [(0.0000166, 18), (0.0000114, 9), (0.0000269, 10), (0.0000236, 9), (0.0000274, 17)],
    "PFTeDA": [(0.0000107, 24), (0.00000546, 14), (0.00000731, 10), (0.00001813, 8), (0.0000455, 25)],
    "PFBS":   [(0.242, 4), (0.159, 4), (0.259, 7), (0.303, 11), (0.650, 8)],
    "PFHxS":  [(0.0292, 6), (0.0208, 6), (0.0468, 3), (0.0554, 7), (0.0668, 9)],
    "PFOS":   [(0.00245, 4), (0.00129, 3), (0.00247, 17), (0.00315, 1), (0.00729, 16)],
    "PFOSA":  [(0.000501, 22), (0.000617, 13), (0.000404, 4), (0.000404, 10), (0.000712, 10)],
    "6:2 FtS": [(0.0343, 6), (0.0276, 6), (0.0376, 6), (0.0444, 7), (0.0822, 22)],
}
S6_CD1 = {
    "PFBA":   [(0.312, 9), (0.187, 5), (0.304, 10), (0.335, 5), (0.385, 8)],
    "PFPeA":  [(0.425, 5), (0.167, 10), (0.303, 14), (0.383, 14), (0.445, 12)],
    "PFHxA":  [(0.265, 10), (0.153, 4), (0.317, 11), (0.237, 8), (0.248, 7)],
    "PFHpA":  [(0.137, 13), (0.0756, 8), (0.212, 11), (0.124, 11), (0.139, 7)],
    "PFOA":   [(0.0368, 10), (0.0206, 4), (0.0706, 4), (0.0358, 12), (0.0346, 7)],
    "PFNA":   [(0.0045, 4), (0.0044, 4), (0.00636, 8), (0.00859, 5), (0.00872, 5)],
    "PFDA":   [(0.00253, 9), (0.00107, 6), (0.00193, 2), (0.00231, 3), (0.0025, 7)],
    "PFUnDA": [(0.000544, 5), (0.000182, 18), (0.000519, 25), (0.000593, 7), (0.000526, 24)],
    "PFDoDA": [(0.000126, 20), (0.0000309, 17), (0.0000816, 11), (0.000205, 12), (0.000352, 18)],
    "PFTrDA": [(0.0000267, 21), (0.00000609, 14), (0.0000133, 14), (0.000019, 6), (0.0000342, 22)],
    "PFTeDA": [(0.0000103, 11), (0.00000925, 9), (0.000011, 23), (0.000019, 11), (0.0000266, 11)],
    "PFBS":   [(0.286, 22), (0.107, 5), (0.261, 11), (0.201, 13), (0.262, 6)],
    "PFHxS":  [(0.0432, 5), (0.0282, 5), (0.0365, 2), (0.0499, 9), (0.0536, 5)],
    "PFOS":   [(0.00548, 11), (0.00207, 13), (0.00333, 9), (0.00554, 6), (0.00432, 12)],
    "PFOSA":  [(0.000428, 15), (0.000531, 9), (0.000558, 7), (0.000604, 12), (0.00113, 5)],
    "6:2 FtS": [(0.0496, 15), (0.0227, 3), (0.0383, 9), (0.0544, 19), (0.0519, 9)],
}
S7_B6 = {
    "PFBA":   [(0.39, 12), (0.225, 22), (0.301, 6), (0.548, 11), (0.328, 8)],
    "PFPeA":  [(0.337, 7), (0.192, 16), (0.341, 5), (0.275, 9), (0.432, 6)],
    "PFHxA":  [(0.245, 4), (0.151, 6), (0.209, 6), (0.177, 4), (0.301, 16)],
    "PFHpA":  [(0.0807, 5), (0.0831, 6), (0.0811, 5), (0.0782, 14), (0.134, 12)],
    "PFOA":   [(0.025, 6), (0.0214, 9), (0.0309, 5), (0.0236, 6), (0.0387, 6)],
    "PFNA":   [(0.00614, 7), (0.00400, 9), (0.00588, 3), (0.00499, 4), (0.00885, 5)],
    "PFDA":   [(0.00172, 10), (0.000972, 2), (0.0018, 5), (0.00131, 2), (0.00285, 6)],
    "PFUnDA": [(0.00487, 5), (0.000135, 9), (0.000773, 6), (0.000377, 7), (0.00101, 8)],
    "PFDoDA": [(0.000121, 10), (0.0000224, 11), (0.000168, 5), (0.0000625, 22), (0.00026, 9)],
    "PFTrDA": [(0.0000177, 12), (0.0000114, 9), (0.000017, 5), (0.0000165, 8), (0.0000384, 22)],
    "PFTeDA": [(0.00000918, 8), (0.00000467, 23), (0.0000184, 8), (0.0000128, 18), (0.0000174, 7)],
    "PFBS":   [(0.191, 3), (0.13, 9), (0.195, 9), (0.172, 8), (0.276, 7)],
    "PFHxS":  [(0.0421, 13), (0.0262, 4), (0.038, 9), (0.0284, 0), (0.0538, 3)],
    "PFOS":   [(0.00312, 10), (0.0018, 15), (0.00426, 8), (0.00176, 9), (0.00515, 9)],
    "PFOSA":  [(0.000678, 9), (0.000351, 2), (0.000625, 6), (0.000387, 5), (0.000687, 10)],
    "6:2 FtS": [(0.0298, 7), (0.0246, 3), (0.0408, 4), (0.0273, 5), (0.0448, 4)],
}

ORGANISMS = [
    ("human", "human", "pooled donors; brain pooled from 2 male + 2 female (NIH NeuroBioBank), other tissues from BioIVT", "Table S4", S4_HUMAN),
    ("rat", "rat (Wistar Han)", "pooled by sex; BioIVT and Pfizer DSRD tissue bank", "Table S5", S5_RAT),
    ("mouse", "mouse (CD-1)", "pooled by sex; BioIVT and Pfizer DSRD tissue bank", "Table S6", S6_CD1),
    ("mouse", "mouse (C57BL/6J)", "pooled by sex; BioIVT and Pfizer DSRD tissue bank", "Table S7", S7_B6),
]

# Cells recorded as printed but inconsistent with the paper's own text or with
# the pattern of every other row.  Keyed (strain, chemical, tissue).
SUSPECT = {
    ("human", "PFTeDA", "liver"): (
        "PRINTED AS 6.8e-7 AND RECORDED AS PRINTED. It is 15.7x below the same "
        "cell in the other three organisms and puts liver 6.7x BELOW brain, "
        "where liver exceeds brain for all 15 other compounds here; the main "
        "text also gives the study minimum as 6.8e-6. But the paper's own "
        "'liver within 2-fold of other tissues for 14 of 16 PFAS' reproduces "
        "only with 6.8e-7 (6.8e-6 would give 15 of 16), so this value, not the "
        "main text's, is what the paper analysed. Either the main text's stated "
        "minimum is the typo, or an error in the underlying datum propagated "
        "into the published statistic. Excluded from the derived ratios as an "
        "outlier either way"),
    ("mouse (C57BL/6J)", "PFUnDA", "liver"): (
        "PRINTED AS 4.87e-3 AND RECORDED AS PRINTED, but that is 9-15x above the "
        "same cell in the other three organisms (4.6e-4, 3.2e-4, 5.4e-4), 5-36x "
        "above the other four tissues in its own row, and higher than the "
        "shorter-chain PFDA in the same column, which breaks the monotonic "
        "chain-length decline seen everywhere else. Points to a missing zero, "
        "i.e. 4.87e-4. Excluded from the derived species ratios"),
}

# Non-monotonic steps that are real features of the data, not transcription
# errors.  Each was checked against the same cell in the other three organisms.
NOISE_OK = {
    ("rat (Wistar Han)", "PFPeA", "heart"): (
        "f_unbound in heart RISES from PFBA to PFPeA in three of the four "
        "organisms (human 1.18x, rat 1.70x, C57BL/6 1.13x; CD-1 flat). At the "
        "weakly-bound short-chain end the chain-length trend is flat and the "
        "scatter is wide, so the rat step is the same feature, larger. CV 19%. "
        "This 0.832 is the highest f_unbound in the whole study"),
    ("rat (Wistar Han)", "PFTeDA", "lung"): (
        "all four organisms put lung PFTeDA between 1.7e-5 and 4.6e-5, at the "
        "bottom of the assay's range; rat is the only one where it exceeds "
        "PFTrDA, and it carries the highest CV in that row (25%). Noise at the "
        "assay floor, not a transcription error"),
}

COMMON = dict(
    study="Ryu et al.",
    year="2024",
    pmid="39161261",
    doi="10.1021/acs.est.4c04050",
    method=("rapid equilibrium dialysis (HTD96, 12-14 kDa cellulose), 18 h at 37 C, "
            "200 rpm, 5% CO2, 75% RH, n = 4; LC-MS/MS; diluted f_unbound corrected "
            "to undiluted by the standard dilution-factor equation"),
    ka="", ka_units="", kd="", kd_units="",
    n_binding_sites="", binding_site_identified="",
    protein_conc_units="dilution factor",
    ligand_conc_range="5",
    ligand_conc_units="uM (single concentration, 1% DMSO; chosen to stay below albumin-type saturation)",
    temperature_c="37",
    buffer="Dulbecco's PBS, pH 7.4",
    is_computational="no (experimental; the PLSR descriptor model is a separate strand)",
    assay_volume_note=("150 uL spiked homogenate donor / 150 uL DPBS receiver; membranes "
                       "presaturated twice; matrix-matched samples for LC-MS/MS; "
                       "tolbutamide and carbamazepine internal standards"),
    ligand_over_protein_max="",
)

TK = ("yes - f_unbound and K_tissue/w offered directly as PBTK inputs. log f_unbound "
      "correlates across every species x tissue pair (R2 >= 0.97, Spearman 0.96-1), so "
      "whole-tissue binding is not a species difference; rat liver is the best single "
      "surrogate for human tissues (R2 >= 0.98). Liver binding was within 2-fold of "
      "other tissues for 14/16 PFAS despite L-FABP enrichment")


def rows():
    out = []
    for species, strain, provenance, table, data in ORGANISMS:
        for chem in ORDER:
            ncarb, mw, logd, pka, head, npfc = S1[chem]
            for tissue, (fu, cv) in zip(TISSUES, data[chem]):
                src = f"{table} ({strain} {tissue}); Table S1 for MW/logD/pKa"
                if (strain, chem, tissue) in SUSPECT:
                    src += " | FLAGGED: " + SUSPECT[(strain, chem, tissue)]
                elif (strain, chem, tissue) in NOISE_OK:
                    src += " | NOTE: " + NOISE_OK[(strain, chem, tissue)]
                out.append(dict(
                    COMMON,
                    chemical=chem,
                    n_perfluorinated_carbons=str(npfc),
                    head_group=head,
                    protein=f"{strain} {tissue} homogenate (whole tissue: total protein "
                            f"plus phospholipid, not an isolated protein)",
                    protein_species=species,
                    unbound_fraction=repr(fu),
                    protein_conc="1:10 in DPBS" if tissue == "brain" else "1:5 in DPBS",
                    tk_relevance_noted=TK,
                    source_table=src,
                ))
    return out


def tidy(rs):
    """Analysis-ready long form, one row per chemical x organism x tissue."""
    with open(TIDY, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_perfluorinated_carbons", "head_group", "mw_g_mol",
                    "h_logD", "pKa", "species", "strain", "tissue", "f_unbound",
                    "cv_percent", "flagged"])
        for species, strain, _prov, _tbl, data in ORGANISMS:
            for chem in ORDER:
                ncarb, mw, logd, pka, head, npfc = S1[chem]
                for tissue, (fu, cv) in zip(TISSUES, data[chem]):
                    w.writerow([chem, npfc, head, mw, logd,
                                "" if pka is None else pka,
                                species, strain, tissue, repr(fu), cv,
                                "yes" if (strain, chem, tissue) in SUSPECT else ""])
    print(f"wrote {TIDY}")


def checks(rs):
    """Transcription checks: anything that fails is either a typo of mine or a
    real feature of the paper, and must be explained before the data is used."""
    import collections
    bad = []
    byk = {}
    for r in rs:
        strain = r["protein"].split(" ")[0] if False else None
    lookup = {}
    for species, strain, _p, _t, data in ORGANISMS:
        for chem in ORDER:
            for tissue, (fu, cv) in zip(TISSUES, data[chem]):
                lookup[(strain, chem, tissue)] = fu

    print("\nCHECK 1  f_unbound falls with chain length inside each series")
    for series in (["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
                    "PFUnDA", "PFDoDA", "PFTrDA", "PFTeDA"],
                   ["PFBS", "PFHxS", "PFOS"]):
        for _s, strain, _p, _t, _d in ORGANISMS:
            for tissue in TISSUES:
                vals = [lookup[(strain, c, tissue)] for c in series]
                for a, b, ca, cb in zip(vals, vals[1:], series, series[1:]):
                    if b > a * 1.6:      # tolerate the known short-chain noise
                        known = ((strain, cb, tissue) in SUSPECT
                                 or (strain, ca, tissue) in SUSPECT
                                 or (strain, cb, tissue) in NOISE_OK)
                        bad.append((strain, tissue, ca, a, cb, b, known))
    unflagged = [b for b in bad if not b[-1]]
    for b in bad:
        mark = "  (documented)" if b[-1] else "  <-- UNEXPLAINED"
        print(f"   {b[0]:18s} {b[1]:7s} {b[2]}={b[3]:<10.3g} -> {b[4]}={b[5]:<10.3g}{mark}")
    print(f"   {len(bad)} non-monotonic steps, {len(unflagged)} unexplained")

    print("\nCHECK 2  each cell within 10x of the same cell in the other organisms")
    out = []
    for chem in ORDER:
        for tissue in TISSUES:
            vals = {strain: lookup[(strain, chem, tissue)]
                    for _s, strain, _p, _t, _d in ORGANISMS}
            hi, lo = max(vals.values()), min(vals.values())
            if hi / lo > 10:
                flagged = any((s, chem, tissue) in SUSPECT for s in vals)
                out.append((chem, tissue, hi / lo, vals, flagged))
    for chem, tissue, sp, vals, flagged in out:
        mark = "  (already flagged)" if flagged else "  <-- UNEXPLAINED"
        print(f"   {chem:8s} {tissue:7s} spread {sp:6.1f}x  "
              + " ".join(f"{k.split(' ')[0][:4]}={v:.3g}" for k, v in vals.items()) + mark)
    print(f"   {len(out)} cells over 10x spread, "
          f"{sum(1 for o in out if not o[-1])} unexplained")

    print("\nCHECK 3  study range matches the main text (6.8e-6 to 0.3)")
    allv = sorted(lookup.values())
    print(f"   min {allv[0]:.3g}  max {allv[-1]:.3g}  "
          f"({allv[-1]/allv[0]:,.0f}x, {len(allv)} values)")
    nosus = sorted(v for k, v in lookup.items() if k not in SUSPECT)
    print(f"   excluding the two flagged cells: min {nosus[0]:.3g}  max {nosus[-1]:.3g}")
    return len(unflagged) + sum(1 for o in out if not o[-1])


def main():
    with open(CSV) as fh:
        rdr = csv.DictReader(fh)
        fields = rdr.fieldnames
        existing = list(rdr)
    before = len(existing)
    dropped = sum(1 for r in existing if r["study"] == "Ryu et al.")
    existing = [r for r in existing if r["study"] != "Ryu et al."]
    print(f"dropped {dropped} existing Ryu row(s)")

    new = rows()
    unexplained = checks(new)

    for r in new:
        assert set(r) <= set(fields), set(r) - set(fields)
        existing.append({f: r.get(f, "") for f in fields})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(existing)
    print(f"\n{before} -> {len(existing)} rows ({len(new)} Ryu values added)")
    tidy(new)
    if unexplained:
        raise SystemExit(f"\n{unexplained} unexplained check failures - resolve before use")


if __name__ == "__main__":
    main()
