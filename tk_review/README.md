# PFAS toxicokinetics review

Why PFAS serum half-lives differ between humans, rats and mice; what each
reported half-life assumed; and whether exposure level is related to half-life
or to volume of distribution.

Built on top of the existing work in this repository (`../literature`,
`../species_dose`, `../pfas_dose`), which it extends rather than repeats.
See the [root README](../README.md) for how this fits with the other strands.

## Where things are

| Path | What it holds |
|---|---|
| `report/REPORT.md` | the written review — start here |
| `report/Argoul_comparison_memo.pdf` | 2-page standalone memo: how the Argoul cocktail experiment was run, which comparisons it supports, why six were dropped; rebuild with `scripts/make_argoul_memo.py` |
| `report/PFAS_TK_review.pdf` | the printable review — 16 pages, 14 sections, 4 figures, 57 numbered references; rebuild with `scripts/make_summary_pdf.py` (needs `reportlab`) |
| `EXECUTIVE_SUMMARY.md` | the same argument in two pages |
| `db/` | the structured extractions — one CSV per evidence type |
| `db/primary_2026/` | 24 per-paper extractions, read from the papers directly |
| `db/combined/` | the consolidated database — 1,438 rows across four tables, plus an `.xlsx` |
| `db/qsar/` | per-compound structure descriptors paired with the renal handling ratio (§11) |
| `papers/` | full-text extractions (`.txt`) and the retrieval logs (`SOURCES_*.md`) |
| `figures/` | 19 analysis figures; `fig00_master.png` is the whole argument in one panel pair |
| `scripts/` | everything that generates the above, runnable in order |
| `research_notes/` | the per-subtopic research notes behind the extractions |
| `WANTED.md` | papers still wanted, and why each one would settle something |

## The databases

| File | Rows | What it answers |
|---|---|---|
| `source_catalogue.csv` | 81 | the user-supplied citation list, each resolved to a verified PMID and classified by evidence type |
| `master_exposure_halflife.csv` | 46 | dose, measured serum concentration, half-life, Vd and clearance on one line per chemical × species × sex |
| `regulatory_values.csv` | 173 | what each agency adopted, from which primary study, under which assumptions |
| `agency_clearance_comparison.csv` | 16 | the human clearance factor each agency uses, normalised to mL/kg-day |
| `protein_binding.csv` | 273 | in vitro PFAS binding to albumin and other proteins, with the method for each |
| `transporter_kinetics.csv` | 217 | transporter Km, Vmax and IC50, with the reabsorption-or-secretion direction |
| `transporter_mechanism.csv` | — | transporter expression, hormonal regulation and species comparison |
| `animal_halflife_measured.csv` | — | measured (not modelled) animal half-lives by chemical, species, strain, sex, dose, route |
| `human_halflife_extended.csv` | — | human half-life estimates with their design and modelling assumptions |
| `vd_clearance.csv` | — | volume of distribution by species, flagged measured / fitted / assumed, with provenance |
| `dose_dependence.csv` | — | multi-dose studies: half-life at each dose level |
| `mouse_rat_decomposition.csv` | 14 | the mouse/rat half-life ratio split into its Vd and clearance terms |
| `sex_decomposition.csv` | 21 | the same split applied to the female/male difference within each species |
| `reabsorption_axis.csv` | 5 | renal reabsorption of PFOA by species, and the half-life it predicts |
| `cphea_fitted_halflives.csv` | 186 | terminal slopes refitted here from the EPA CPHEA raw curves — **read the `tail_selection_flag` column and §8 item 9 before using these** |

### `db/combined/` — the consolidated deliverable

| File | Rows | What it holds |
|---|---|---|
| `tk_parameters.csv` | 736 | half-life, clearance, Vd, MRT, bioavailability, GFR and reabsorption by chemical × species × sex × source, units normalised, every row carrying its provenance |
| `binding.csv` | 192 | protein binding constants and unbound fractions with the method for each |
| `transporters.csv` | 260 | transporter Km, PBPK Tm/KT, mRNA sex ratios and direction |
| `regulatory.csv` | 215 | what each agency adopted, from which study, under which assumptions |
| `PFAS_TK_combined.xlsx` | — | all four as one filterable workbook |

Rebuild with `python3 scripts/build_combined_datasets.py` then
`python3 scripts/export_combined_xlsx.py`.

### `db/primary_2026/` — the eleven full texts supplied 2026-10-02

| File | Rows | What it holds |
|---|---|---|
| `argoul2026_mouse_tk.csv` | 11 | Argoul 2026 Tables 1/2/4: clearance, Vss, MRT, bioavailability, renal/faecal split, unbound fraction and renal clearance for 11 PFAS in female mice |
| `argoul2026_allometry_human.csv` | 11 | mouse-to-human allometric scaling against measured human clearances |
| `kudo2002_rat_tk.csv` | 6 | Kudo 2002 Tables 2–3: Vd, clearance and half-life in both sexes, plus the probenecid arm |
| `kudo2002_transporter_mrna.csv` | 6 | the six renal transporters, their sex ratios, and whether each actually transports PFOA |
| `lou2009_mouse_tk.csv` | 9 | Lou 2009 Table 2: sex-resolved mouse PFOA Vd, ke and half-life in serum, liver and kidney |
| `lou2009_saturable_resorption.csv` | 10 | Lou 2009 Table 4: the saturable-resorption model parameters, Tm and KT |
| `thompson2010_vd_calibration.csv` | 2 | Thompson 2010 supplementary Table S1 — the two communities the 170 mL/kg human Vd was calculated from |
| `han2012_table4_reabsorption_axis.csv` | 11 | the reabsorption axis at its source, with OEHHA's adapted value alongside each |
| `han2012_table6_transporter_km.csv` | 13 | rat and human renal transporter Km values for PFCAs |
| `han2012_table7_pbpk_tm_kt.csv` | 8 | the Tm and KT constants published PFOA PBPK models run on |
| `yang2010_human_apical_transporters.csv` | 4 | human OATP1A2, OAT4 and URAT1 PFOA uptake, with the assay conditions |
| `kudo2001_chain_length_elimination.csv` | 8 | rat urinary and faecal recovery by chain length and sex |
| `tatumgibbs2011_pfna_rat_mouse.csv` | 4 | the only strain-matched rat-vs-mouse PFNA experiment |
| `sundstrom2012_pfhxs_three_species.csv` | 10 | PFHxS in rat, mouse and monkey, both sexes in all three, one laboratory |
| `cheng2005_mouse_oatp_sex.csv` | 7 | mouse renal and hepatic Oatp sex predominance and androgen dependence |

## Reproducing it

```bash
pip install openpyxl matplotlib
cd tk_review

python3 scripts/parse_xlsx_list.py       # classify the supplied citation list
python3 scripts/resolve_pmids.py         # resolve each to a PMID via NCBI
python3 scripts/verify_pmids.py          # verify by title match, retry misses
python3 scripts/finalize_source_list.py  # apply manual corrections

python3 scripts/build_master_table.py          # dose / serum / half-life table
python3 scripts/rat_mouse_decomposition.py     # where the species gap comes from
python3 scripts/vd_vs_dose.py                  # does Vd track exposure?
python3 scripts/agency_clearance_comparison.py # why agencies disagree
python3 scripts/reabsorption_axis.py           # the single mechanistic axis
python3 scripts/make_figures.py                # all five figures

# the five primary full texts supplied 2026-10-02
python3 scripts/thompson_vd_circularity.py            # is the human Vd a measurement?
python3 scripts/primary_sex_species_decomposition.py  # sex vs species, from primaries
python3 scripts/argoul_reabsorption_check.py          # the axis, independently recomputed
python3 scripts/validate_cphea_fits.py                # our own fits vs published values
python3 scripts/han2012_axis_at_source.py             # the axis at its origin
python3 scripts/saturation_margin_km_vs_kt.py         # in vitro Km vs PBPK KT

# the deliverables
python3 scripts/build_combined_datasets.py            # -> db/combined/*.csv
python3 scripts/export_combined_xlsx.py               # -> db/combined/*.xlsx
python3 scripts/make_figures_primary.py               # -> figures/fig10..fig15

# the QSAR endpoint, and how far its one dataset can be trusted
python3 scripts/qsar_endpoint_table.py                # -> db/qsar/qsar_endpoint_table.csv
python3 scripts/argoul_vs_single_compound.py          # -> db/argoul_vs_single_compound{,_sources}.csv
python3 scripts/argoul_dose_matching.py               # is the dose gap the cause? (no)
python3 scripts/make_review_figures.py                # -> figures/fig02..fig04
python3 scripts/make_summary_pdf.py                   # -> report/PFAS_TK_review.pdf
```

`EXECUTIVE_SUMMARY.md` is the two-page version: the five findings, what this
work corrects in the published record (including in its own earlier output), and
the one question still genuinely open.

`SUMMARY.md` is the illustrated standalone summary with full citations.

`WANTED.md` lists only what is still missing — the 11 full texts supplied during
this work have been removed from it.

`scripts/fetch_papers.sh` re-downloads the agency PDFs listed in
`papers/DOWNLOAD_MANIFEST.csv`.

## Provenance

Every number in `db/` carries the study, year, PMID or DOI, and the table or
page it came from. Nothing is quoted from memory. Where a value was taken from
another paper's citation rather than the original, the row says so.

**The PDFs are not in git.** They are ~110 MB of freely available government
documents at stable URLs. Their text extractions (`papers/*.txt`) *are*
committed, and those are what every number was read from. `papers/SOURCES_*.md`
logs every retrieval attempt, successes and failures alike, with URLs and dates.

## Two things to know before reading the databases

**The supplied citation list is a hazard corpus, not a toxicokinetic one.**
Of its 81 entries, 42 cannot contribute a half-life, Vd, clearance or binding
constant — they are immunotoxicity, steatosis and gene-expression studies. A
further 21 are repeat-dose toxicity studies that measured internal
concentration without reporting kinetics; those are flagged
`supplies_internal_dose` because a measured serum level at a known administered
dose is still a usable datum.

**Units are mixed in the upstream source.** `../species_dose/species_exposure.csv`
reports Chiu's human half-lives in *years* and the EPA animal fits in *days*,
with no units column. `master_exposure_halflife.csv` carries an explicit
`native_time_unit` column and normalised `halflife_d` / `halflife_y` columns.
