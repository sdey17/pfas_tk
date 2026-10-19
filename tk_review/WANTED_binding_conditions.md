# Papers needed to finish the binding assay-condition columns

`db/protein_binding.csv` now has protein concentration on 659 of 741 rows
(89%), ligand range on 90%, temperature on 89%, buffer on 89%. The 82
remaining rows need the papers below.

Reproduce the fill: `python3 tk_review/scripts/add_binding_assay_conditions.py`,
then `python3 tk_review/scripts/fill_binding_conditions_primaries.py`, then
`python3 tk_review/scripts/add_binding_studies_2026_10b.py` (which reads
`db/starnes2025_tableS4.csv`, the extracted SI table, now vendored in the repo),
then `python3 tk_review/scripts/add_ryu2024_si_funbound.py`.

**Cleared 2026-10-16** (five primaries plus the Ryu SI supplied): Gao 2025,
Chen & Guo 2009, Han 2021, Ryu 2024 with its Supporting Information, and the
Starnes 2025 SI are all extracted and in the database. See the "Already
complete" list at the bottom.

Validate the file at any time with
`python3 tk_review/scripts/check_protein_binding.py`, which reproduces
published statistics from three of the source papers out of the transcribed
rows.

## Tier 1 — most rows gained per paper

| rows | study | DOI | PMID | why it matters |
|---|---|---|---|---|
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

- **Jackson et al. 2021**, 10.1021/acs.est.1c01200, PMID 34495656 — the DSF
  method paper Starnes 2024 and 2025 both build on
- **Klevens & Ellenbogen 1954**, 10.1039/df9541800277 — not indexed in PubMed
  or Europe PMC; would need a library copy of Discuss Faraday Soc 18:277

## Still outstanding after the 2026-10-16 batch

- **Chen et al. 2025 Supporting Information**, 10.1021/acs.est.4c06734 — the
  per-congener K_MW (membrane-water) and K_A (HSA) for all 60 PFAS are in the
  SI; only the regression slopes and the group ordering were extractable from
  the article. This is the largest single membrane-water dataset in the field
  and would let Armitage's membrane-lipid spLFER be tested directly

- **Zhang et al.** (Han 2021 ref 34) — the fluorescence-displacement source of
  every L-FABP K_d that Han 2021 plots. Han measures ΔT_m, not K_d, so the
  affinities currently in the database for L-FABP are second-hand
- **Gao et al. 2025 Supporting Information** — Tables S3 (van't Hoff ΔH, ΔS,
  ΔG at four temperatures) and S12 (fitted n_Hill per compound per salt level).
  Only the Table 1 affinities were extracted
- **Starnes et al. 2025 Table S3** — per-chemical top concentrations, needed to
  compute ligand:protein for those 68 rows

## Already complete, for reference

- Fischer et al. 2024 — 111 rows, complete from Table S2
- Starnes et al. 2024 — 32 rows, filled from Methods
- Wu et al. 2024 — 20 rows (NOTE: 17–50% ethanol v/v, flagged on every row)
- Peng et al. 2024 — 12 rows
- Alesio & Bothun 2024 — 11 rows
- Moro et al. 2022 — 8 rows (NOTE: defatted hSA, so an upper bound on affinity)
- Chen et al. 2020 — 1 row, conditions filled; its Table 2 Ka values are not in
  the saved flattened text, so the Ka stays a Starnes 2026 secondary range
- Chen & Guo 2009 — 20 rows, the full 5 × 3 compound-by-site grid from Table 1
  plus five non-fluorinated fatty-acid reference constants from the discussion
- Gao et al. 2025 — 42 rows, the whole of Table 1 including the NaCl series
- Han et al. 2021 — 13 rows, ΔT_m panel + four species + two recombinant +
  palmitic acid. K_d values flagged as quoted from Zhang et al., not measured
- Ryu et al. 2024 — 320 rows, SI Tables S4-S7 complete: 16 PFAS x 5 tissues x
  4 organisms including CD-1 vs C57BL/6J. Two cells flagged as probable SI
  typos and excluded from derived ratios; see
  `scripts/add_ryu2024_si_funbound.py`. Tidy form in `db/ryu2024_funbound.csv`
- Starnes et al. 2025 — 68 rows from SI Table S4, 49 binders and 19 non-binders
- Qin et al. 2024 — 64 rows, Table 2 complete: 16 PFAS x BSA / human plasma /
  trout plasma, specific and nonspecific terms separated. The only dataset in
  the file that reports BOTH regimes for the same protein and compound
- Chen et al. 2025 — 1 study-level row; per-congener values need the SI (above)
