# PFAS toxicokinetics: why the half-lives disagree

Published PFAS half-lives disagree by orders of magnitude — between humans
and animals, between rat and mouse, between male and female, and between
papers describing the same experiment. This repository works out why, one
quantity at a time, by re-deriving the numbers from the tables they came
from rather than from the sentences about them.

It began as a Python re-fit of one paper (Chiu et al. 2022) and grew into six
connected strands. The through-line is a single identity:

```
t½ = ln2 · Vd / CL
```

Half-life is not an independent quantity. Any two of the three terms fix the
third, so nearly every disagreement in this literature turns out to be about
**which two were measured and which one was assumed**. Finding the assumed
constant is most of the work.

## Start here

| if you want… | read |
|---|---|
| what the whole project found | [`SUMMARY.md`](SUMMARY.md) — six phases, the short version |
| the species question answered in full | [`tk_review/report/REPORT.md`](tk_review/report/REPORT.md) (~2,000 lines) |
| that review in two pages | [`tk_review/EXECUTIVE_SUMMARY.md`](tk_review/EXECUTIVE_SUMMARY.md) |
| something to print or hand to someone | [`tk_review/report/PFAS_TK_review.pdf`](tk_review/report/PFAS_TK_review.pdf) — 16 pages, A4, 4 figures, 57 references |
| the whole argument in one figure | [`tk_review/figures/fig00_master.png`](tk_review/figures/fig00_master.png) |
| to learn toxicokinetic modelling | [`tk_learning/`](tk_learning/) — twelve runnable lessons |
| the data behind the lessons | [`tk_learning/data/SOURCES.md`](tk_learning/data/SOURCES.md) — nine studies, traced to figure and table |

## The six strands

```
chiu_replication/   Chiu et al. 2022's human Bayesian model, re-fitted in PyMC
tk_learning/        twelve lessons: C(t)=C0·exp(-kt) through to a working PBPK model
pfas_dose/          EPA animal PK pipeline; does clearance depend on dose within a species?
species_dose/       does dose explain the human/rodent gap? (no — ~13% at best)
literature/         23 human studies appraised by design rather than reputation
tk_review/          the species question at the source: 186 full texts, 1,438 database rows
```

`tk_learning/` and `tk_review/` are the pair to know: the first teaches the
methods on real animal data, the second applies them to the published
literature. Each folder has its own README with the detail. In brief:

### `chiu_replication/` — one paper, reproduced and debugged

Re-fits the Bayesian population model of Chiu et al. 2022 (*Environ Health
Perspect* 130(12):127001) in PyMC, against Chiu's own MCSim input files.
Reproduces all four published half-lives:

| PFAS | this code | paper, Table 3 |
|---|---|---|
| PFOA | 3.16 (2.65–3.77) | 3.14 (2.69–3.73) |
| PFOS | 3.34 (2.45–4.45) | 3.36 (2.52–4.42) |
| PFNA | 2.28 (1.53–3.19) | 2.35 (1.65–3.16) |
| PFHxS | 8.52 (5.66–13.01) | 8.30 (5.38–13.5) |

An earlier version gave PFOA 4.30 y. Four bugs explain the gap, and the
folder README works through all four — including why fixing them one at a
time makes each look harmless.

### `tk_learning/` — the methods, from scratch

Twelve runnable lessons on real animal data, 29 figures, 12 executed
notebooks, 53 worked questions and 35 regression checks. Builds from a
one-compartment IV bolus to a four-compartment PBPK model, and spends its
time on the places where the modelling actually goes wrong: four defensible
terminal-slope fits to one dataset giving half-lives 4.4–10.2 d, two
parameter sets 520× apart producing an identical curve, and four volumes of
distribution spanning 17× in the same two-compartment fit.

### `pfas_dose/` and `species_dose/` — testing the dose hypothesis

If PFAS elimination saturates, high-dose animal studies would read as
long half-lives and the species gap would partly be a dose artefact.
`pfas_dose/` fits the EPA animal curves (1- vs 2-compartment, chosen by
LOO) and finds a real but weak within-species dose effect for PFOA
(slope ≈ 0.11) and none for PFHxA. `species_dose/` matches human and rodent
studies on serum concentration and finds the dose explanation accounts for
~13% of the gap at best. Female rats clear PFOA **44× faster than males at
identical doses** — so the variance is in biology, not dose.

### `literature/` — trusting designs, not authors

Published human PFOA half-lives span 17-fold (0.5–8.5 y). The appraisal
tiers studies by design and locates the spread in six assumptions. The
structural finding: serum-decay designs read the rate off a slope and never
need a volume of distribution; mass-balance designs need both a Vd *and* a
complete excretion accounting, each uncertain 2–6×, multiplying straight
into the answer.

### `tk_review/` — the species question, answered

186 full texts read, 24 per-paper extractions, 1,438 rows across four
consolidated tables, 16 figures, 34 runnable scripts. Headline results:

- **It is clearance, not distribution.** Vd spans 1.33–2.18× between sexes;
  clearance spans 0.78–44.3× and changes sign between species.
- **One mechanistic axis orders every species** — fractional renal
  reabsorption, from 99.94% in humans to net secretion in the female rat —
  and predicts every rodent within 1.7× over a 37× span with no fitted
  parameter.
- **Humans run two near-complete reabsorption loops,** renal and biliary.
- **The human Vd that nine regulatory clearance factors inherit is a
  calculation,** reproducible to three figures from its source paper's own
  supplementary table given an assumed half-life.
- **18 corrections to the published record,** each traced to the table that
  contradicts it.
- **Half-life is the wrong structure-activity endpoint.** §11 proposes the
  dimensionless renal handling ratio `R = CL_renal/(fu·GFR)` instead, which
  divides out body size and binding and leaves the transport step.

## Conventions this repository holds itself to

Every number traces to a named table in a named document. Every database row
carries a `provenance` field and a PMID or DOI. Where a value was computed
here rather than read from a paper, the row says so. Where this project's own
method is known to be biased — the terminal-slope fitter on biphasic curves —
the affected rows carry a flag and the README says to read §8 item 9 before
using them. Corrections to earlier conclusions in this repository are kept in
the text rather than quietly edited out, because the failure mode they
illustrate (believing a sentence over a table) is the subject.

## Running things

Each strand installs separately; there is no repository-wide environment.

```bash
# the Chiu re-fit
cd chiu_replication && pip install -r requirements.txt && python model.py PFNA

# the lessons
cd tk_learning && pip install -r requirements.txt
python 02_explore.py && python test_lessons.py

# the review's analyses (numpy / pandas / matplotlib / scipy)
cd tk_review && python scripts/qsar_endpoint_table.py
```

## Licensing and attribution

`chiu_replication/data/` holds Chiu et al.'s MCSim input files unmodified,
under GPLv3 — see `chiu_replication/data/LICENSE-chiu-GPLv3`. The original
work is at [wachiuphd/2022-Bayes-PFAS-PK](https://github.com/wachiuphd/2022-Bayes-PFAS-PK).
Animal PK curves come from the EPA's
[CPHEA-Animal-PFAS-PK](https://github.com/USEPA/CPHEA-Animal-PFAS-PK)
database; `tk_learning/data/SOURCES.md` traces each one to its original
study. PDFs under `tk_review/papers/` are retrieved copies of
third-party publications, kept for extraction provenance and subject to
their publishers' terms.
