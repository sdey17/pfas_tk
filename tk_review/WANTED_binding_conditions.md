# Papers needed to finish the binding assay-condition columns

`db/protein_binding.csv` now has protein concentration on 196 of 285 rows
(69%), ligand range on 73%, temperature on 69%, buffer on 69%. The 89
remaining rows need the papers below.

Reproduce the fill: `python3 tk_review/scripts/add_binding_assay_conditions.py`
then `python3 tk_review/scripts/fill_binding_conditions_primaries.py`.

## Tier 1 — most rows gained per paper

| rows | study | DOI | PMID | why it matters |
|---|---|---|---|---|
| 10 | Gao et al. 2025, Environ Sci Technol | 10.1021/acs.est.5c08779 | 41344664 | paywalled, not in PMC. Only study we hold on how ionic strength, cation valency and pH change PFAS–albumin affinity. Abstract-only at present, with no numeric affinities at all |
| 8 | Chen & Guo 2009, Arch Toxicol | 10.1007/s00204-008-0359-x | 18854981 | dansyl-L-proline displacement, site-resolved. One of the few papers assigning PFAS to Sudlow I vs II |
| 6 | Bischel et al. 2010, Environ Sci Technol | 10.1021/es101334s | 20540534 | equilibrium dialysis + nanoESI-MS on BSA; a methods-comparison anchor |
| 5 | Beesoon & Martin 2015, Environ Sci Technol | 10.1021/es505399w | 25826685 | branched vs linear isomers — the only isomer-resolved binding data we have |
| 4 | Han et al. 2003, Chem Res Toxicol | 10.1021/tx034005w | 12807361 | 19F NMR, human and rat side by side. Already used elsewhere in the review for the free-fraction adjudication |

## Tier 2 — three rows each

| study | DOI | PMID |
|---|---|---|
| Bischel et al. 2011, Environ Toxicol Chem | 10.1002/etc.647 | 21842491 |
| MacManus-Spencer et al. 2010, Anal Chem | 10.1021/ac902238u | 20039637 |
| Li et al. 2021, Toxics | 10.3390/toxics9030063 | 33803062 |
| Sheng et al. 2016, Arch Toxicol | 10.1007/s00204-014-1391-7 | 25370009 |

Li et al. 2021: body text is already saved; what is missing is the
**Supporting Information** holding the 26-congener Ka table. That one is an SI
retrieval, not a new paper.

## Tier 3 — two rows each

| study | DOI | PMID |
|---|---|---|
| Han et al. 2004, Drug Chem Toxicol | 10.1081/dct-200039725 | 15573471 |
| Chen et al. 2015, Chemosphere | 10.1016/j.chemosphere.2014.11.040 | 25497588 |
| Luebker et al. 2002, Toxicology | 10.1016/s0300-483x(02)00081-1 | 12093614 |
| Ren et al. 2016, Toxicology | 10.1016/j.tox.2016.08.011 | 27528273 |

## Tier 4 — one row each (34 studies)

Not worth chasing individually unless a specific one matters to an argument.
The ones most likely to be wanted:

- **Han et al. 2021**, 10.1021/acs.est.1c00509, PMID 34133149 — quantitative
  interspecies thermal shift; the direct competitor to Starnes 2024 on
  cross-species comparison
- **Ryu et al. 2024**, 10.1021/acs.est.4c04050, PMID 39161261 — rapid
  equilibrium dialysis across human, C57BL/6 and CD-1 mouse, Wistar Han rat.
  The only row we hold with mouse strain resolution, which bears on the
  Argoul strain question in §3.2a of the main report
- **Jackson et al. 2021**, 10.1021/acs.est.1c01200, PMID 34495656 — the DSF
  method paper Starnes 2024 and 2025 both build on
- **Starnes et al. 2025**, 10.1093/toxsci/kfae124, PMID 39298512 — per-congener
  Kd for 65 PFAS are in the SI, which we did not retrieve
- **Klevens & Ellenbogen 1954**, 10.1039/df9541800277 — not indexed in PubMed
  or Europe PMC; would need a library copy of Discuss Faraday Soc 18:277

## Already complete, for reference

- Fischer et al. 2024 — 111 rows, complete from Table S2
- Starnes et al. 2024 — 32 rows, filled from Methods
- Wu et al. 2024 — 20 rows (NOTE: 17–50% ethanol v/v, flagged on every row)
- Peng et al. 2024 — 12 rows
- Alesio & Bothun 2024 — 11 rows
- Moro et al. 2022 — 8 rows (NOTE: defatted hSA, so an upper bound on affinity)
- Chen et al. 2020 — 1 row, conditions filled; its Table 2 Ka values are not in
  the saved flattened text, so the Ka stays a Starnes 2026 secondary range
