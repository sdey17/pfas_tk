# Transporters: what exists, what is wrong with it, and what to do next

Working document, written 2026-11-18 from a direct survey of the repository
rather than from the earlier narrative. Figures in §1 were counted from the
files; where the existing retrieval log disagrees, the files win and the
disagreement is named.

**Four priority full texts arrived on 2026-11-19** — Argoul 2026, Ebert 2020,
Nakagawa 2008 and Niu 2026, requests 1, 2, 3 and 6 — and are ingested. They
settle §4.3, overturn the framing of request 2, and give aim 1 its first
convergent answer (§8). They also correct a count this document reported: see
§3.6.

**Stage 0 is complete** (§7), and doing it corrected three claims in the first
draft of this document. Those corrections are kept in place rather than edited
away, marked ~~struck~~, because two of them are the kind of error that would
otherwise be repeated: a second-hand table counted as new data, and a count I
got wrong in my own survey.

---

## 1. What already exists

| Artefact | Size | State |
|---|---|---|
| `db/transporter_kinetics.csv` | **301 rows, 28 sources**, 29 DOIs | the record — everything read from a primary, including rows that are not transporter kinetics |
| `db/transporter_kinetics_working.csv` | **363 rows, 16 sources, 36 transporters, 22 compounds** | gated, normalised, split; built by `scripts/build_transporter_working_table.py`, never edited by hand |
| `db/membrane_partitioning.csv` | 9 rows | Ebert 2020; a different quantity, kept out of the transporter tables (§6 note) |
| `db/primary_2026/argoul2026_table1_human_reabsorption.csv` | 9 rows | the in vivo counterpart: net tubular reabsorption as % of filtration clearance |
| `db/transporter_mechanism.csv` | 83 rows | direction, sex and hormone regulation; Km on 24 |
| `db/primary_2026/*.csv` | 24 extracts | per-paper tables read from primaries |
| `papers/*.txt` | 29 transporter papers | **text only — no PDFs**, see §5 |
| `papers/SOURCES_transporters.md` | retrieval log | detailed, and now partly stale |
| `report/REPORT.md` §7 | narrative | direction-first, with three open conflicts |
| `scripts/saturation_margin_km_vs_kt.py` | analysis | the Km-versus-KT conflict, §4.1 below |

**Fill rates**, record and working table, after Stage 0:

| Field | Working table | Comment |
|---|---|---|
| `km_uM` | **58 / 363** | 10 studies, 8 compounds, 17 transporters |
| `vmax` | 38 / 363 | one basis after normalisation (§3.2) |
| **`clint_uL_min_mg`** | **51 / 363** | V_max/K_m, 7 studies; see below |

So the table is still mostly qualitative — and after the four new papers it is
*more* so, because Argoul 2026 contributes 120 substrate calls and no
constants at all. The quantitative core is about 58 rows.

The **CL_int column now holds 51 rows, 2 to 279 µL/min/mg protein**, and the
first pattern in it is worth stating carefully:

| Direction | n | median CL_int | range |
|---|---|---|---|
| secretion (basolateral uptake, apical efflux) | 23 | 24.2 | 4.6–279 |
| reabsorption (apical uptake) | 14 | **74.3** | 2–176 |
| hepatic uptake and other | 14 | 30.4 | 3.1–236 |

Reabsorptive transporters carry a ~3× higher median intrinsic clearance than
secretory ones. That is the direction that lengthens half-life, and it is the
first quantitative support in this compilation for the reabsorption-dominated
picture that Argoul's Table 1 shows in vivo (§8). **It is not yet an argument**:
the rows come from different expression systems with no REF/RAF scaling
(§2.3), so this is a hypothesis for Stage 3 to test, not a result.

---

## 2. Why this is not the binding problem again

The albumin and L-FABP work had one quantity, one identity, one sheet. That
does not transfer, for three reasons.

**2.1 The useful quantity is a clearance, not a constant.** For toxicokinetics
what matters is intrinsic clearance, CL_int = V_max / K_m, and neither factor
alone is interpretable. A low K_m with a negligible V_max moves nothing.

**2.2 Direction sets the sign.** The same CL_int lengthens or shortens
half-life depending on which membrane the transporter sits in:

- **Reabsorptive** (apical, lumen → cell; lengthens half-life): rat Oatp1a1,
  human OAT4, URAT1, ASBT, and NTCP/OSTα-β on the enterohepatic loop.
- **Secretory** (basolateral uptake then apical efflux; shortens half-life):
  OAT1, OAT3, and the apical pumps P-gp and BCRP.

A compilation that pools these without the sign is worse than no compilation.
`transporter_mechanism.csv` holds the direction calls; they have to be joined,
not assumed.

**2.3 The in vitro number needs a scaling factor to mean anything in vivo.**
A V_max measured in transfected HEK293 is per mg of that cell line's protein
at that expression level. Reaching a tissue clearance needs a relative
expression or activity factor, and that factor is where IVIVE either works or
quietly fails. The repository already contains the two papers that attempt it
(Lin 2024, Sakolish 2026) and 29 of the 230 rows are their composite-RPTEC
IVIVE outputs rather than measurements — those must not be pooled with measured
kinetics, and the working table now gates them out.

---

## 3. Defects found in the survey

**3.1 One real sync gap, and it was smaller than I first wrote. RESOLVED.**
Louisse 2024 was the only case. Its Table 1 holds fifteen rows, ten with a
K_m — but ~~ten~~ **six** were new:

| | `transporter_kinetics.csv` (before) | Table 1 | folded in |
|---|---|---|---|
| rows | 3 qualitative | 15 | 16 per-compound |
| with a K_m | 0 | 10 | **6** |

The other four are Louisse **2023**'s OAT4 values re-tabulated, digit-identical
to rows already in the record from the 2023 primary (60/4.5, 47/4.5, 58/8.5,
92/7.3), and the re-tabulation loses precision: it prints the standard error as
15 for all four where the primary gives 33 for PFHxS. They are stored with the
reason and **no constant**, which is the treatment Han 2021 and Yang 2026 get on
the binding side. Done by `scripts/add_louisse2024_oat_km.py`.

The extraction verifies against itself: efficiency is V_max/K_m, so the printed
third column is predicted by the first two, and all six reproduce to the
rounding (3.4/159 → 21.4 against a printed 21).

~~Yang 2010 and Han 2012 need the same check.~~ Both were checked and neither
needed a fold-in:

- **Yang 2010 was already synced** — all four rows, the three K_m values and
  the OATP1A2 negative result.
- **Han 2012 Table 6 is a compilation, not a measurement.** All thirteen of its
  values are already in the record attributed to the paper that measured them
  (Weaver 2010, Nakagawa 2008, Yang 2009). Table 6 also silently picks between
  conflicting primaries where the record keeps both: rat Oat1/PFOA 51
  (Nakagawa) against Weaver's 43.2, Oat3/PFOA 80.2 against 65.7, Oatp1a1/PFOA
  162.2 (Yang 2009) against Weaver's 126.4. Its extract additionally
  mis-assigned URAT1's 64.1 µM as a third OAT4 value, because Han's Table 6
  formats OAT4 and URAT1 as one row; the record had it right and the extract is
  corrected.

**3.1b A second-hand V_max must not meet a first-hand K_m. NEW.** Yang 2010's
OAT4 (37) and URAT1 (0.32 nmol/min/mg) V_max values are in this database only
as Louisse 2023's introduction quotes them — the retrieval log said so and the
table did not enforce it. Computing V_max/K_m from them builds a clearance out
of two papers, so the working table drops the V_max, keeps the K_m, and flags
the row. This costs two CL_int values, and they are the two that matter most
for the species question in Stage 4: human OAT4 and URAT1 reabsorption. **Yang
2010's own V_max is therefore a live request, not a nicety.**

**3.2 V_max units are not comparable, and one row is not a V_max. RESOLVED in
the working table.** The first two spellings are now one; the rest are carried
with a flag and excluded from the clearance column.

| Unit | Rows | Comment |
|---|---|---|
| `nmol/mg protein/min` | 26 | same unit as the next row, written differently |
| `nmol/min/mg protein` | 8 | → 34 rows share a basis once normalised |
| `nmol/min/10^6 cells` | 2 | needs mg protein per 10⁶ cells to join |
| `mg/h/kg BW^0.75` | 2 | a whole-body PBPK T_m, not an in vitro V_max |
| `uL/mg protein/min` | 1 | labelled V_max; it is a clearance |

**3.3 The retrieval log was stale in one place. RESOLVED.** It called the
Louisse 2024 K_m table "the main remaining numeric gap"; that line is replaced
and a dated reconciliation section appended to
`papers/SOURCES_transporters.md`. The rest of the log held up, including its
own record of the second-hand V_max in §3.1b.

**3.3b ~~29 sources~~ 27.** The §1 figure in the first draft of this document
was an error in my own survey, not a change from Stage 0: the pre-Stage-0 file
has 27 distinct study strings and 28 DOIs, and still does.

**3.4 Enzymes are in the transporter table.** Three carboxylesterase CES1 and
three CES2 rows, flagged in the `transporter` field as "(enzyme, NOT a
transporter)". They were included because they were on a requested list. They
need a scope decision, not silent carriage.

**3.6 I reported a K_m count that was 20% phantom. CORRECTED.** The Stage 0
summary gave 59 K_m values in the working table. Ten of those were numbers
parsed out of placeholder prose: `"within 1.41-2.78 (individual value only in
Fig 2b)"` yielded a K_m of 1.41, and `"not extractable from text (reported
only in Fig 1f)"` yielded a K_m of **1**, from the figure number — then
multiplied by four when its multi-compound cell was split. The true Stage 0
figure was **49**. The *gain* reported for Stage 0 (+6, from Louisse 2024) was
right; the absolute count was not. `num()` in the builder now requires a cell
to begin with a digit, so prose cannot masquerade as a measurement. All ten
belonged to Niu 2026 and have since been replaced by the real constants.

**3.7 Two defects inside the new papers.** Both were caught by the same method
this project uses throughout — checking a paper against its own internal
redundancy.

- **Niu 2026 transposes P-gp and BCRP** in its V_max/K_m paragraph and in its
  abstract. *Documented, not escalated — reviewer's decision, 2026-11-19; no
  author contact.* The figure values are what this database stores. Dividing the figures' own fitted constants gives, for Fig 3c
  (P-gp), HFPO-TA 227 / PFOA 212 / PFOS 139, and for Fig 3f (BCRP), PFOA 279 /
  HFPO-TA 233 / PFOS 24 — each set matching what the text attributes to the
  *other* transporter, on all six values. Three checks say the figures are
  right: the abstract's own K_m ranges (BCRP 1.29–7.92, P-gp 1.41–2.78) track
  the figures correctly; the identical arithmetic on OAT1 reproduces the text
  exactly; and the text cites panels 3b,d for constants that are in 3c,f. **The
  abstract's headline — highest V_max/K_m is PFOA by P-gp — should read BCRP.**
- **Nakagawa 2008's Table 1 mislabels its third column.** It is headed
  "K_m/V_max (µl/mg protein/min)" and the numbers are indeed K_m/V_max, but
  those are the units of V_max/K_m, a clearance. A reader taking the printed
  0.10 as a clearance is wrong by ~100×, and the ranking inverts: by the
  printed column rOAT3 is highest at 0.22; by actual clearance it is the lowest
  of the four.

**3.5 Composite and non-specific entities. RESOLVED in the working table.**
29 rows are "composite RPTEC proximal-tubule transport", which has no single
transporter identity to join a direction to, and are gated out. Multi-compound
cells like "PFBA, PFPeA (C4, C5)" are split one row per compound, each carrying
a `shared_across` count so one observation spread over seven compounds is never
read as seven observations; 21 cells split into 84 rows. One cell names no
specific compound at all ("6 PFCAs, 3 PFSAs and 4 fluorinated alternatives")
and is gated as UNSPECIFIED.

---

## 4. Conflicts already identified and still open

**4.1 Three orders of magnitude on the saturation question.**
`scripts/saturation_margin_km_vs_kt.py` is the most developed piece of work
here and it is worth reading before anything else. In vitro K_m values put
half-saturation 19–4,865× **above** human serum concentrations; the K_T values
that agency PBPK models run on put contaminated communities and occupational
cohorts **at or above** half-saturation. The script argues for the in vitro
side on three grounds, the strongest being that Han 2012's own Table 7
contains two rows marked "derived from an in vitro measurement" that are ~500×
higher than the fitted values *inside the same table*. This conclusion should
be carried into the transporter write-up, not re-derived.

**4.2 URAT1.** Yang 2010 reports PFOA K_m = 64.1 µM, but only with zero
extracellular chloride — a non-physiological condition, already noted in the
extract. Louisse 2023 found URAT1 transported nothing. Unresolved.

**4.3 ~~Argoul 2026 against everyone.~~ RESOLVED — it agrees with most of
them.** The full text shows the "against everyone" framing was wrong, built
from an abstract. Argoul agrees with Ryu 2024 on the chain-length threshold,
with Louisse 2023 that URAT1 transports no PFAS, and with Louisse 2023/2024 on
PFOA, PFNA and PFHxS. The genuine disagreements are narrow and specific:

| Pair | Argoul 2026 | The other study |
|---|---|---|
| OAT4 / PFBS | substrate | **not** a substrate (Louisse 2023) |
| OAT4 / PFDA | not a substrate | substrate, K_m 39 µM (Louisse 2023) |
| OAT1 / PFOA | **not** a substrate | substrate, K_m 185 µM (Louisse 2024) |
| OAT3 / PFBS | substrate | not detected (Louisse 2024) |
| rat Oatp1a1 / PFOS | no uptake | substantial transport (Zhao 2017) |
| OSTα/β / PFBS, PFHxS, PFOS | no uptake | all three transported (Zhao 2015) |
| rat Asbt / PFBS, PFHxS | substrates | no transport (Zhao 2015) |

Argoul says of the OAT1/PFOA clash only that "the reasons for these
discrepancies remain unclear". Note what the pattern is **not**: it is not a
chain-length or head-group effect, and the two labs disagree in both
directions. The most likely common factor is the mixture — Argoul dosed eleven
PFAS together at 0.2 or 2 µM each, Louisse dosed singly — and Argoul's own
observation that OAT3 and OAT4 ratios fall between 0.2 and 2 µM means the
transporters were near saturation, where competition between eleven
co-administered substrates would bite hardest. That is a testable explanation
and it belongs in the write-up.

**4.3b A real conflict survives, and it is about PFOS.** Argoul finds PFOS a
substrate of none of twelve transporters; Louisse 2023 reports OAT4/PFOS with
K_m 48 µM and V_max 2.2, and Zhao 2017 reports rat Oatp1a1/PFOS. This is the
one disagreement that changes a conclusion, because PFOS is the compound whose
reabsorption the PBPK models most depend on.

**DECISION (reviewer, 2026-11-19): hold both readings and let Stage 1's
inclusion rule decide.** Neither is demoted now. This is the right call for a
reason worth writing down: Louisse's row carries a K_m and a V_max and
Argoul's carries an uptake ratio, so the inclusion rule separates them on
evidence type before anyone has to adjudicate the biology. If the rule admits
K_m-with-V_max to the kinetics table and sends substrate calls to the
direction table, the conflict resolves itself structurally — PFOS appears as a
measured OAT4 substrate in one and a negative call in the other, both visible,
neither overwriting the other.

**4.4 Cheng 2006 contradicts itself.** One sentence of its abstract calls renal
Oatp1a1 "female-predominant", which conflicts with Cheng 2005 and with its own
androgen-stimulation conclusion. Needs the full text to resolve.

**4.5 Two numbers that must stay out of any kinetics table.** Worley & Fisher's
*fitted* male/female apical activity ratio of ~25,800×, and the Tm/K_T sex
ratios of exactly 91× in Loccisano and Kim. These are scalars absorbing an
observed clearance difference. REPORT.md §7 already says so.

---

## 5. No PDFs exist for this strand

Every transporter paper in `papers/` is a `.txt` extraction. The log records
that PMC returns an anti-bot interstitial for `/pdf/` and that Springer and
Elsevier return HTML landing pages, for **all** transporter papers tried. The
regulatory PDFs (ATSDR, EPA, OEHHA) came from other strands.

Consequence: anything in a figure rather than a table is currently unreachable.
Niu 2026's per-congener K_m and V_max exist "only inside Figures 1f and 2b,d".

---

## 6. Papers to request, in priority order

Abstract-only unless noted. Items 1–4 change conclusions; 5–9 fill the table;
10–13 settle the sex-difference argument.

| # | Paper | Why it matters |
|---|---|---|
| ~~1~~ | **Argoul 2026** — ✅ **OBTAINED 2026-11-19**, 66 rows ingested | Settled §4.3 and supplied half the answer in §8. Reports uptake ratios only, no constants |
| ~~2~~ | **Ebert 2020** — ✅ **OBTAINED**, 9 rows to `db/membrane_partitioning.csv` | **The framing here was wrong.** See the note below |
| ~~3~~ | **Niu 2026** — ✅ **OBTAINED**, 9 K_m/V_max read off Figures 2f and 3c,f | Also exposed a P-gp/BCRP transposition in the paper's own text (§3.7) |
| 4 | **Yang 2010**, *Toxicol Sci* 117:294, 10.1093/toxsci/kfq219 | **Now the top request.** K_m is in hand; the V_max values are second-hand, so §3.1b costs us the two CL_int values the species question needs most |
| 5 | **Yang 2009**, *Toxicol Lett* 190:163, 10.1016/j.toxlet.2009.07.011 | Oatp1a1 K_m and the C6–C10 K_i,app series, currently quoted via Weaver 2010's Discussion |
| ~~6~~ | **Nakagawa 2008** — ✅ **OBTAINED**, Table 1 supplies the four missing V_max values | Also a units defect worth knowing about (§3.7) |

**Request 2 was mis-specified, by me.** Ebert 2020 is a membrane partitioning
and permeability paper. Its entire transporter content is four passages: one
sentence of the abstract, one paragraph of the introduction, one paragraph of
the results, and one clause of the conclusions. "Reframes the whole strand"
was built from that single abstract sentence, read when only the abstract was
available, and it is far more weight than the paper asks to carry. It should
not have been request 2.

The paper's own conclusion is agnostic, and points where §8 points: "either of
the transport processes might dominate permeation, depending on compound
concentration or **amount and types of transporter proteins in the cell**."
That is expression-dependence, not a claim that transporters are unnecessary.

What it does bear on is one background assumption — that the uptake seen in
the **control** cells of transporter assays is passive diffusion. Ebert's
bilayer permeabilities reproduce effective permeabilities back-calculated from
the control cells of its refs 53, 55 and 59–61, which are Yang 2009, Weaver
2010 and three others already in this table, within a factor of 3 for PFSAs
and less clearly for PFCAs. That corroborates those papers' own subtraction
step. It says nothing about the saturable component they report, and
contradicts no K_m here.

Its value to this project is in the other strand: +0.63 log K_mem/w per
perfluorinated carbon, and PFSAs sorbing ~1.2 log units more strongly than
PFCAs of equal chain length — directly comparable to the per-carbon slope and
head-group offset already computed for albumin and L-FABP, which makes
membrane lipid a fourth binding phase for aim 3. Per-compound tables are in SI
S3.1.1 and S3.2.1, **not supplied**; that is a live request, and not one to be
met by digitising Figures 1 and 2.
| 7 | **Han 2008**, *Toxicol Lett* 181:81, 10.1016/j.toxlet.2008.07.002 | Rat hepatocyte uptake, male vs female K_m and V_max — the hepatic side |
| 8 | **Kimura 2017**, *Toxicol Lett* 277:18, 10.1016/j.toxlet.2017.05.012 | Caco-2 K_m 8.3 µM — an order of magnitude below every other K_m in the table, so either important or wrong |
| 9 | **Kummu 2015**, 10.1016/j.placenta.2015.07.119 | OAT4 and placental transfer; n = 4 |
| 10 | **Cheng 2006**, 10.1124/mol.106.025122 | Resolves the self-contradiction in §4.4 |
| 11 | **Cheng 2005**, 10.1124/dmd.105.003640 | Mouse renal Oatp sex differences |
| 12 | **Buist 2004**, 10.1124/dmd.32.6.620 | The direct rat-vs-mouse Oat comparison |
| 13 | **Kudo 2002**, 10.1016/s0009-2797(02)00006-6 | Sex-hormone regulation; the 70× male/female half-life |

Not obtainable: **Katakura 2007**, *J Health Sci* 53:77 — not indexed in
PubMed or Europe PMC, no DOI found. Its results survive only second-hand in
OEHHA PHG Table A6.3 and Weaver 2010's Discussion.

---

## 7. Proposed sequence

**Stage 0 — reconcile what is already here. COMPLETE, 2026-11-18.** Two
scripts, neither of which edits the record by hand:

- `scripts/add_louisse2024_oat_km.py` — folds in Louisse 2024 Table 1, six new
  K_m+V_max pairs, four carry-overs stored without constants, six negatives.
  Idempotent.
- `scripts/build_transporter_working_table.py` — gates, normalises, splits and
  computes CL_int into `db/transporter_kinetics_working.csv`. Six gates
  (NOT_MEASURED, NOT_A_TRANSPORTER, COMPOSITE, CARRIED_OVER, SECOND_HAND,
  NOT_VMAX) drop 52 of 230 rows; splitting adds 84 from 21 cells; CL_int comes
  out on 38.

It did not double the K_m count (52 → 58, not ~100), because the Louisse 2024
gap was six values rather than ten and the other two suspected gaps were not
gaps. It did deliver Stage 3's clearance column early, and it found the
second-hand-V_max defect in §3.1b, which no amount of new paper-fetching would
have surfaced.

**Stage 1 — decide the inclusion rule**, as `screen_fits_database.py` does for
binding. Candidate: a row earns a place if it reports **K_m with V_max**, or a
K_i with the probe substrate and its K_m. Uptake ratios and substrate/
non-substrate calls go to a separate direction table, not the kinetics table.

**Stage 2 — build the direction join.** Every transporter gets apical or
basolateral, reabsorptive or secretory, per species, with the source for the
localisation. Without this the kinetics table cannot be read.

**Stage 3 — compute CL_int = V_max/K_m** on a common basis, carrying the
expression system, and keep the IVIVE scaling explicitly separate so a reader
can see what was measured and what was scaled.

**Stage 4 — the species question**, which is the project's aim 1. The rat
Oatp1a1 reabsorption story has no clean human counterpart: Yang 2010 found
OATP1A2, the closest human orthologue, does **not** mediate saturable PFOA
uptake. Whether human OAT4 and URAT1 fill that role is exactly what §4.2 leaves
open, and item 1 of the request list is the best available evidence.

**Stage 5 — write it up** as a further section of `report/ALBUMIN_BINDING.md`,
or as its own document if it outgrows that.

---

## 8. Aim 1 has a convergent answer, from two independent routes

The four new papers do something the roadmap did not anticipate: two of them
reach the same conclusion about the species question by methods that share no
assumptions.

**Nakagawa 2008, from matched kinetics.** It measured human and rat OAT1 and
OAT3 *in the same expression system*, which is the only way to compare species
without also changing laboratory. Recomputing its clearances correctly (§3.7):

| | human | rat | human/rat |
|---|---|---|---|
| OAT1 | 9.99 | 8.26 | 1.21× |
| OAT3 | 8.13 | 4.56 | 1.78× |

against a species difference in PFOA half-life of roughly 50×. Its own
conclusion: species differences in renal elimination "are not attributable to
affinity differences in these OATs".

**Argoul 2026, from substrate specificity.** Ten PFAS against eight human and
four rat transporters. The same molecular-weight threshold — about 500 g/mol —
governs transporter interaction in both species: PFCAs up to nine carbons and
PFSAs up to six are substrates of human OAT3/OAT4 and rat Oat1/Oatp1a1, and
PFDA, PFOS and PFDS are substrates of none of the twelve, in either species.
Its conclusion: "transporter selectivity may be largely conserved, and
additional mechanisms may contribute to interspecies differences".

**So the transporters' handling of PFAS is not where the species difference
lives.** Affinity is within 2×; substrate specificity is shared. Both papers
name the same missing quantity — **transporter expression level** — and Argoul
says so explicitly: "data on the expression levels of these transporters in
human and rodent kidneys are needed".

**Niu 2026 supplies the human half of that number**, measured rather than
assumed: in human kidney cortex, OAT1 is 5.3 ± 1.9 pmol/mg protein, P-gp
2.1 ± 0.8, and BCRP below detection (n = 20). The rat counterpart is the
single most valuable thing still missing from this strand, and it reframes
Stage 4: the species question is now an **expression** question, not an
affinity one.

**What the in vivo side says.** Argoul's Table 1 assembles human unbound
fractions with the clearances from Abraham 2024:

| PFAS | net reabsorption, % of filtration clearance |
|---|---|
| PFOA | 99.0 |
| PFNA, PFDA | 98.5 |
| PFOS | 97.9 |
| PFHxS | 97.5 |
| PFBS | 97.3 |
| GenX | 77.1 |
| **PFHxA** | **9.6** |

Reabsorption is essentially complete for everything except PFHxA — and PFHxA
is the one compound Argoul finds to be an OAT1 substrate, OAT1 being
secretory. That is a clean, independent corroboration of the direction calls,
and it is the strongest single argument in this document that direction must
be joined before kinetics (Stage 2) rather than after.
