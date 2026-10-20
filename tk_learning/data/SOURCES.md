# Where the lesson data comes from

Every number in `tk_learning/data/*.csv` traces to a published figure, table
or appendix in a named study. This file is that trace, end to end.

**First, a correction worth making explicitly: there is no mouse data in these
five files.** They cover two species — rat and cynomolgus monkey (labelled
`primate`). Mouse PFAS kinetics are well covered in the upstream database, just
not in the subset exported for the lessons; see [§4](#4-the-mouse-data-and-everything-else-not-in-the-lessons).

## How to cite it

**EPA compiled these data; it did not generate them.** Citing only the EPA
paper under-credits the nine laboratories that ran the experiments, so cite
the compilation and cite the primary study for any specific result you lean
on:

> Animal serum time-course data were obtained from the US EPA CPHEA PFAS
> pharmacokinetic compilation (Zurlinden et al. 2025, *Toxicol Appl
> Pharmacol* 499:117336;
> https://github.com/USEPA/CPHEA-Animal-PFAS-PK), which digitised and
> transcribed them from the primary studies listed in §2.

The compilation paper is Zurlinden TJ, Dzierlenga MW, Kapraun DF, Ring C,
Bernstein AS, Schlosser PM, Morozov V (2025), *Estimation of species- and
sex-specific PFAS pharmacokinetics in mice, rats, and non-human primates
using a Bayesian hierarchical methodology*,
[doi:10.1016/j.taap.2025.117336](https://doi.org/10.1016/j.taap.2025.117336),
open access at
[PMC12172007](https://pmc.ncbi.nlm.nih.gov/articles/PMC12172007/). The
repository is MIT licensed, © U.S. Federal Government. Its methods confirm
the digitisation account below — "Graphical extraction of concentration vs.
time data was conducted using WebPlot Digitizer" — and its data-availability
statement points at the repository.

Nothing in `data/` is synthetic and nothing originates with this repository.
Every row is a serum or plasma concentration measured in a published animal
study, extracted by EPA, and re-exported here in a flat format so the
lessons run without a database.

## The chain

```
original study                           9 studies, §2
   └─ published figure / table / supplement
        └─ digitised or transcribed by EPA CPHEA      <- the lossy step, quantified below
             └─ USEPA/CPHEA-Animal-PFAS-PK  (extracted_data/*.csv, one per study)
                  └─ local copy in this repo
                     tk_review/papers/CPHEA-Animal-PFAS-PK extracted_data/
                       └─ subset + unit harmonisation  <- verified exact, 0 of 1,574 rows differ
                            └─ tk_learning/data/*.csv   (5 files, 1,574 rows)
```

Only one link in that chain loses information, and it is not ours: the last
step is verified value-for-value by `data/verify_against_epa.py`, so all the
uncertainty sits in EPA's extraction, where the table in §2 says per file
exactly how much of it there is.

Links for the two upstream layers:

| layer | where |
|---|---|
| EPA database, source code and extracted data | https://github.com/USEPA/CPHEA-Animal-PFAS-PK |
| EPA HERO (the reference IDs in the `study` column) | `https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/<ID>` |

The `study` column of every lesson CSV holds an **EPA HERO reference ID**, and
the `author` column holds the first author. Those two together identify the
original paper uniquely, and §2 resolves all nine.

## 1. What is in each lesson file

Reproduce this table with:

```bash
python3 - <<'PY'
import csv, collections, os
for f in sorted(os.listdir("data")):
    if not f.endswith(".csv"): continue
    rows = list(csv.DictReader(open(os.path.join("data", f))))
    print(f"=== {f} ({len(rows)} rows) ===")
    by = collections.Counter((r["study"], r["author"]) for r in rows)
    for (s, a), n in by.most_common():
        print(f"   HERO {s:<9s} {a:<12s} {n:>4d} rows")
PY
```

| lesson file | rows | contributing studies (HERO ID) |
|---|---|---|
| `PFOA_Male_primate.csv` | 43 | Butenhoff 3749227 (43) |
| `PFOS_Male_primate.csv` | 48 | Chang 1289832 (48) |
| `PFOA_Male_rat.csv` | 860 | Kemper 6302380 (682), Dzierlenga 5916078 (126), Kim 3749289 (18), Iwabuchi 3859701 (13), Kudo 2990271 (11), Ohmori 3858670 (10) |
| `PFOA_Female_rat.csv` | 377 | Kemper 6302380 (238), Dzierlenga 5916078 (110), Kim 3749289 (18), Kudo 2990271 (7), Ohmori 3858670 (4) |
| `PFOS_Male_rat.csv` | 246 | Huang 5387170 (108), Chang 1289832 (98), Kim 3749289 (29), Iwabuchi 3859701 (11) |

EPA ships these as plain CSVs, one per study, in
[`extracted_data/`](https://github.com/USEPA/CPHEA-Animal-PFAS-PK/tree/main/extracted_data),
named `Author_HEROID.csv`, with column definitions in that folder's
`README.txt`. The SQLite database the EPA notebooks use (`PFAS.db`) is built
*from* those CSVs by `auxiliary_notebooks/create_db.ipynb`; it is not itself
checked in, so nothing needs building to read the upstream values.

So the monkey curves come from **one study each**, and the rat files are
**pooled across six and four studies respectively**. That matters for lesson 10:
a species or sex contrast drawn from the rat files is a contrast across a study
mixture, not within one experiment.

## 2. The nine original studies

Resolved against PubMed; DOIs link to the publisher of record.

| # | HERO | citation | species / sex | what the lessons use |
|---|---|---|---|---|
| 1 | [3749227](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/3749227) | Butenhoff JL, Kennedy GL, Hinderliter PM, et al. (2004) Pharmacokinetics of perfluorooctanoate in cynomolgus monkeys. *Toxicol Sci* 82(2):394–406. [DOI](https://doi.org/10.1093/toxsci/kfh302) | monkey, M+F | the 10 mg/kg IV arm, 3 males, 123 d — **Table 5** |
| 2 | [1289832](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/1289832) | Chang SC, Noker PE, Gorman GS, et al. (2012) Comparative pharmacokinetics of perfluorooctanesulfonate (PFOS) in rats, mice, and monkeys. *Reprod Toxicol* 33(4):428–440. [DOI](https://doi.org/10.1016/j.reprotox.2011.07.002) | rat, mouse, monkey | monkey 2 mg/kg IV (**Fig 11**); rat gavage 2/4.2/15 mg/kg (**Figs 2, 3, 7, 8**) |
| 3 | [6302380](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/6302380) | **Kemper RA (2003) Perfluorooctanoic acid: toxicokinetics in the rat. Unpublished report, DuPont-7473; US EPA public docket AR-226-1499.** Haskell Laboratory, E.I. du Pont de Nemours. | rat, M+F | 0.1–25 mg/kg IV and gavage — **Appendices B, H, I, J, K, L, M** |
| 4 | [5916078](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/5916078) | Dzierlenga AL, Robinson VG, Waidyanatha S, et al. (2020) Toxicokinetics of PFHxA, PFOA and PFDA in male and female Hsd:Sprague Dawley SD rats following intravenous or gavage administration. *Xenobiotica* 50(6):722–732. [DOI](https://doi.org/10.1080/00498254.2019.1683776) | rat, M+F | PFOA IV and gavage arms — **supplementary data files** |
| 5 | [5387170](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/5387170) | Huang MC, Dzierlenga AL, Robinson VG, et al. (2019) Toxicokinetics of PFBS, PFHxS and PFOS in male and female Hsd:Sprague Dawley SD rats after intravenous and gavage administration. *Toxicol Rep* 6:645–655. [DOI](https://doi.org/10.1016/j.toxrep.2019.06.016) — **and its corrigendum**, *Toxicol Rep* 8:365 (2021). [DOI](https://doi.org/10.1016/j.toxrep.2021.02.001) | rat, M+F | PFOS IV and single-dose gavage — **supplementary data files** |
| 6 | [3749289](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/3749289) | Kim SJ, Heo SH, Lee DS, Hwang IG, Lee YB, Cho HY (2016) Gender differences in pharmacokinetics and tissue distribution of 3 perfluoroalkyl and polyfluoroalkyl substances in rats. *Food Chem Toxicol* 97:243–255. [DOI](https://doi.org/10.1016/j.fct.2016.09.017) | rat, M+F | PFOA and PFOS, 1–2 mg/kg — **Figs 4, 5, 6** |
| 7 | [2990271](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/2990271) | Kudo N, Katakura M, Sato Y, Kawashima Y (2002) Sex hormone-regulated renal transport of perfluorooctanoic acid. *Chem Biol Interact* 139(3):301–316. [DOI](https://doi.org/10.1016/s0009-2797(02)00006-6) | rat, M+F | 48.63 µmol/kg (= 20.14 mg/kg) IV — **Figure 1** |
| 8 | [3858670](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/3858670) | Ohmori K, Kudo N, Katayama K, Kawashima Y (2003) Comparison of the toxicokinetics between perfluorocarboxylic acids with different carbon chain length. *Toxicology* 184(2–3):135–140. [DOI](https://doi.org/10.1016/s0300-483x(02)00573-5) | rat, M+F | PFOA arm of the C7–C10 series — **Figs 1A, 1B** |
| 9 | [3859701](https://hero.epa.gov/hero/index.cfm/reference/details/reference_id/3859701) | Iwabuchi K, Senzaki N, Mazawa D, et al. (2017) Tissue toxicokinetics of perfluoro compounds with single and chronic low doses in male rats. *J Toxicol Sci* 42(3):301–317. [DOI](https://doi.org/10.2131/jts.42.301) | rat, M only | 0.1 mg/kg gavage, PFOA (C8A) and PFOS (C8S) — **Figs 1, 2** |

Worth noticing: **Kudo 2002 is titled "Sex hormone-regulated renal transport
of perfluorooctanoic acid."** The mechanism lesson 10 appeals to for the
28-fold sex difference is the subject of one of the papers the data come
from. That makes the observation and the explanation non-independent —
reassuring rather than circular, but worth knowing you are not getting two
separate pieces of evidence.

### Three things to know before trusting a row

**Entry 3 is the largest single source and it is not a paper.** Kemper 2003
supplies 682 of 860 rows in `PFOA_Male_rat.csv` and 238 of 377 in
`PFOA_Female_rat.csv` — 74% and 63%. It is an unpublished DuPont study
submitted to EPA under TSCA; there is no DOI and no journal version. The only
public route to the primary document is the EPA docket record (AR-226-1499) via
HERO. EPA's digitisation of its appendices is, in practice, the accessible form
of that data — which is why it is listed as unobtainable in
`tk_review/WANTED.md`.

**Group means vs individual animals.** The `animal_id` column is `-1` where the
source published only group means. Of the nine studies, Kudo, Ohmori, Kim and
Iwabuchi are **means only**; Butenhoff, Chang, Kemper, Dzierlenga and Huang give
individual animals. Fitting a hierarchical model (lesson 06) across a mixture of
the two is a modelling decision, not a neutral default.

**The digitisation layer is real and lossy — but it is quantified.** EPA
records a `source` field per *record*, naming the table, appendix or figure
each value came from, so the question is answerable row by row rather than
study by study. Resolved per lesson file:

| lesson file | from tables / appendices | digitised from figures |
|---|---|---|
| `PFOA_Male_primate.csv` | **43 (100%)** — Butenhoff Table 5 | 0 |
| `PFOA_Male_rat.csv` | 808 (94%) | 52 (6%) |
| `PFOA_Female_rat.csv` | 348 (92%) | 29 (8%) |
| `PFOS_Male_rat.csv` | 108 (44%) | 138 (56%) |
| `PFOS_Male_primate.csv` | 0 | **48 (100%)** — Chang figures |

This matters more than it looks, and it points the opposite way to a blanket
caveat. The monkey PFOA curve that lessons 02, 03, 04, 06, 07, 08 and 09 are
all built on is **transcribed from Butenhoff's Table 5**, not read off a
plot — so the dataset carrying most of the teaching is the most reliable one
in the set. `PFOS_Male_primate.csv`, by contrast, is wholly digitised, so
treat its third decimal place as decorative.

None of the lessons' conclusions turn on a few per cent: the effects they
rest on are 2-fold, 11-fold and 28-fold. If one ever did, go to the original
paper rather than the CSV.

### Verified against EPA, row by row

`data/verify_against_epa.py` re-derives the table above and checks every
value against the upstream files, normalising each study's own dose and
concentration units on the way:

```bash
# against the copy already in this repository
python3 data/verify_against_epa.py \
    "../tk_review/papers/CPHEA-Animal-PFAS-PK extracted_data"

# or against EPA's repository directly
git clone --depth 1 https://github.com/USEPA/CPHEA-Animal-PFAS-PK /tmp/epa
python3 data/verify_against_epa.py /tmp/epa/extracted_data
```

```
1574 values compared, 0 unmatched
largest disagreement anywhere: 1.421e-14 mg/L
```

That residual is floating-point noise from the unit conversions, not a data
difference. So the export changed no values, and the uncertainty lives
entirely in EPA's extraction step, where the table above says exactly how
much of it there is.

## 3. What was changed on the way in

Two transformations stand between the upstream CSVs and `data/`, and both are
documented in code rather than only here.

**Unit harmonisation.** Upstream carries each study's own units
(`nmol/ml`, `µg/ml`, `ng/ml`, `µg/L`; `umol/kg` or `mg/kg`; minutes, hours or
days). The lesson files are uniformly **days, mg/L, mg/kg**.

**One repaired bug.** Kudo (2990271) and Ohmori (3858670) report dose in
**µmol/kg**, and the original build copied that number into a column named
`dose_mgkg` without converting — so those rows read 48.63 where the dose is
20.14 mg/kg, a factor of 1000/414.07 = 2.415 (PFOA's molar mass). Concentration
and time were converted correctly, and the `dataset` label was built from a
correctly converted dose, so the label served as the authority for the repair.
This mattered because lesson 10 computes `CL = dose/AUC`, so the affected
datasets came out 2.415× too high in exactly the lesson about species and sex
differences. See `../fix_dose_units.py`, which reports before it writes:

```bash
python3 ../fix_dose_units.py          # report only
python3 ../fix_dose_units.py --write  # apply
```

Not affected: `pfas_dose/`, which reads EPA's own processed data through
`get_processed_data(dose_label="dose_mg")`; and
`tk_review/scripts/fit_cphea_halflives.py`, which filters on
`dose_units == "mg/kg"` explicitly.

## 4. The mouse data, and everything else not in the lessons

The local copy of the upstream extraction holds **20 studies, 7,962 rows**
across three species — far more than the five lesson files use:

| species | rows |
|---|---|
| rat | 5,071 |
| mouse | 2,445 |
| monkey (`primate`) | 446 |

```bash
ls "../tk_review/papers/CPHEA-Animal-PFAS-PK extracted_data/"
```

The eight studies carrying **mouse** rows, all gavage, both sexes except where
noted:

| file | mouse rows | compound |
|---|---|---|
| `Lou_2919359` | 1,152 | PFOA |
| `Tatum-Gibbs_2919268` | 480 | PFNA |
| `Lau_6579272` | 429 | PFBS |
| `Sundstrom_1289834` | 135 | PFHxS |
| `Chang_1289832` | 96 | PFOS |
| `Chang_2325359` | 57 | PFBA |
| `Daikin_6822782` | 53 | PFHxA (female only) |
| `Gannon_2850314` | 43 | PFHxA |

Compound identity is carried as a DTXSID. Dzierlenga and Huang label their
`source` column with the compound name as well, which gives a direct mapping
for six of the nine:

| DTXSID | compound | how established |
|---|---|---|
| `DTXSID8031865` | PFOA | Dzierlenga `source` labels |
| `DTXSID3031862` | PFHxA | Dzierlenga `source` labels |
| `DTXSID3031860` | PFDA | Dzierlenga `source` labels |
| `DTXSID3031864` | PFOS | Huang `source` labels |
| `DTXSID7040150` | PFHxS | Huang `source` labels |
| `DTXSID5030030` | PFBS | Huang `source` labels |
| `DTXSID8031863` | PFNA | inferred: the only compound common to Tatum-Gibbs (a PFNA study), Iwabuchi's C9A arm and Ohmori's series; cross-checked against the fitted male rat half-life of 17.8 d |
| `DTXSID1037303` | PFHpA | inferred: the one compound left in Ohmori's C7–C10 series; cross-checked against the fitted female half-life of 0.051 d vs Ohmori's published 0.05 d |
| `DTXSID4059916` | PFBA | inferred from the study (Chang 2325359 = rat + mouse + monkey, 3–300 mg/kg, which is the PFBA design); **not** read from a mapping file — verify before relying on it |

To build a mouse lesson file, the loader needs no change — only a subset
written in the same column layout. The sex contrast is the interesting one,
because its *sign* differs from the rat: see `tk_review/REPORT.md` §3.2a
and `db/primary_2026/lou2009_mouse_tk.csv`.

## 5. Related data elsewhere in this repository

| path | what it is | source |
|---|---|---|
| `../tk_review/db/primary_2026/*.csv` | 24 per-paper extractions made by reading the papers directly, not via EPA | each row carries `source_table` and a PMID or DOI |
| `../tk_review/db/cphea_fitted_halflives.csv` | half-lives fitted in this project from the upstream curves | **read `REPORT.md` §8 item 9 first** — the terminal-window selection is biased on biphasic curves |
| `../tk_review/db/combined/` | the consolidated 2,351-row database | every row carries `provenance` |
| `../../chiu_replication/data/*.R` | Chiu et al. 2022 human population model inputs | the Chiu paper and its SI are in `chiu_replication/paper/` |
