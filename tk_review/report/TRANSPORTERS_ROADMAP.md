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
| `db/transporter_kinetics_working.csv` | **58 rows** (tier A 52, B 2, C 4) | the constants. Built by `scripts/build_transporter_working_table.py`, never edited by hand |
| `db/transporter_direction.csv` | **305 rows** (tier D) | substrate calls, uptake ratios, inhibition — a different kind of evidence, in a different file |
| `scripts/screen_transporter_fits.py` | the inclusion screen | run before reading a paper; the counterpart of `screen_fits_database.py` |
| `db/transporter_localisation.csv` | **36 rows** | Stage 2's authority: membrane, direction and net process per (transporter, species, tissue) |
| `scripts/join_transporter_direction.py` | the join and its audit | read-only; reports coverage, conflicts and CL_int regrouped on direction |
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

Those three groupings were made from the free text each paper happened to
write. **§10 redoes them against the authority table and the gap narrows to
~1.9×** — 36.1 µL/min/mg for processes that lengthen half-life against 18.9
for those that shorten it, on comparable whole-cell assays. **It is still not
an argument**: the rows come from different expression systems with no
REF/RAF scaling (§2.3), so this is a hypothesis for Stage 3 to test.

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
PubMed or Europe PMC, no DOI found. Its results survive second-hand in OEHHA
PHG Table A6.3 and Weaver 2010's Discussion — and since the reviewer has
ruled the PHG directly reliable, **that is now a sufficient source** for its
two negatives (no Mrp2-mediated PFOA transport, no rat Npt2). It is removed
from the wanted list.

---

## 6b. Rat renal transporter expression — the papers exist

Found 2026-11-19 in answer to the question "is there a paper for rat renal
transporter expression levels?". According to PubMed:

| Paper | What it gives | Access |
|---|---|---|
| **Thakur A, Yue G, Ahire D, … Prasad B** (2024). "Sex and the Kidney Drug-Metabolizing Enzymes and Transporters: Are Preclinical Drug Disposition Data Translatable to Humans?" *Clin Pharmacol Ther* 116(1):235–246. [10.1002/cpt.3277](https://doi.org/10.1002/cpt.3277), PMID 38711199 | **78 DMETs in human, rat AND mouse kidney membrane fractions, both sexes.** 30 of 33 membrane proteins quantified in all three species differ significantly. No significant sex differences in humans; 28 in rats, 46 in mice | **PMC11218045, open access** |
| **Basit A, Radi Z, Vaidya VS, Karasu M, Prasad B** (2019). "Kidney Cortical Transporter Expression across Species Using Quantitative Proteomics." *Drug Metab Dispos* 47(8):802–808. [10.1124/dmd.119.086579](https://doi.org/10.1124/dmd.119.086579), PMID 31123036 | 19 transporters, kidney cortex, 5 species. All 19 higher in preclinical species than human except MDR1, OCT3, OCTN1. OAT1 cross-species variability 4.1× | paywalled |
| **Uchida Y, Toyohara T, Ohtsuki S, Moriyama Y, Abe T, Terasaki T** (2016). "Quantitative Targeted Absolute Proteomics for 28 Transporters in Brush-Border and Basolateral Membrane Fractions of Rat Kidney." *J Pharm Sci* 105(2):1011–1016. [10.1002/jps.24645](https://doi.org/10.1002/jps.24645), PMID 26367854 | **Rat, split by membrane — which maps straight onto §10's authority table.** BBM: bcrp 40.3, urat1 22.2, mate1 8.90, mrp4 4.69, mdr1a 3.22 fmol/µg. BLM: oat1 10.6, oat3 10.2, oct1 4.59 | paywalled |

**Thakur 2024 is the one to read first** — open access, all three species, both
sexes, one assay. Its sex result is strikingly parallel to the PFAS pattern:
no sex dependence in human kidney DMETs, 28 proteins sex-dependent in rat.

~~**The gap that remains**: none quantifies rat Oatp1a1.~~ **Wrong — Thakur
2024 quantifies it** (rat 1.9 ± 2.1 pmol/mg, not detected in human), and it
turned out to carry the sharpest result in this document. That claim was made
from abstracts before the full texts arrived; §11 has what they actually say.
Uchida's 28-protein panel genuinely does not include Oatp1a1, so that half of
the statement stands.

Two incidental finds while looking, both still useful: Kuo 2012
([10.1371/journal.pone.0039641](https://doi.org/10.1371/journal.pone.0039641))
independently confirms Oatp1a1 sits apically in the S3 segment, corroborating
the §10 authority row; and Cerrutti 2002
([10.1023/a:1015563021602](https://doi.org/10.1023/a:1015563021602)) reports
rat OAT1 protein in female basolateral membrane vesicles at 40% of male.
**Basit 2019 could not be obtained** and is no longer needed: Thakur 2024
covers the same ground with both sexes and three species.

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

**Stage 1 — the inclusion rule. COMPLETE, 2026-11-19.** Decided as proposed,
with one addition and one rejection. §9 states it and the evidence behind it.

**Stage 2 — the direction join. COMPLETE, 2026-11-19.** §10 states it.

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
2.1 ± 0.8, and BCRP below detection (n = 20).

**And the rat half now exists.** The reviewer supplied Thakur 2024 and
Uchida 2016 on 2026-11-19; both are ingested to
`db/transporter_expression.csv`. ~~None of them quantifies Oatp1a1~~ —
**Thakur 2024 does**, and that correction is the difference between a gap
and an answer. §11 states the result.

**OEHHA's Table A6.4 is the species comparison in vivo**, and the reviewer
has ruled that the PHG may be relied on directly — it carried dedicated
external peer review, so unlike the compilations in Table A6.3 it is stored
as OEHHA's own analysis
(`db/primary_2026/oehha2024_a64_species_reabsorption.csv`). Net tubular
reabsorption of PFOA, by species and sex:

| | CL_R (mL/kg-day) | % reabsorbed | % excreted |
|---|---|---|---|
| **human** | **0.06** | **99.8** | **0.2** |
| rat male | 18.2 | 93.7 | 6.3 |
| rat female | 666–1009 | **net secretion** | — |
| mouse male / female | 10 / 16 | 97.0 / 95.2 | 3.0 / 4.8 |
| macaque male / female | 15 / 32 | 91.2 / 81.2 | 8.8 / 18.8 |
| dog male / female | 43 / 50.8 | 59 / 52 | 41 / 48 |
| rabbit | 640–670 | **net secretion** | — |

**Two things to read off it, neither of which is the percentage.**

First, **the percentages compress the difference and the escaping fraction
does not.** Human 99.8% against male rat 93.7% looks like six points. The
fraction that *escapes* is 0.2% against 6.3% — **31-fold** — and the renal
clearances differ **303-fold**. The right denominator is what gets out.

Second, **the rat sex difference is qualitative.** Male rats reabsorb ~93%;
female rats are in net tubular *secretion*, in both Kemper 2003 and Ohmori
2003. Whatever drives the ~70× male/female half-life gap changes the **sign**
of net tubular handling, not its magnitude — and humans show no comparable
sex difference. Any mechanism proposed for the rat sex difference has to
account for a sign change, which rules out a simple scaling of one
transporter's abundance.

The one caveat worth attaching: OEHHA assumes f_u = 0.02 for *every* species.
Argoul 2026 measures human PFOA f_u at 0.0017, twelve-fold lower. Both still
land at ~99–100% reabsorption for human, because CL_R is tiny either way, so
the conclusion is robust to that disagreement — but the absolute net
reabsorption in mL/kg-day is not, and should not be quoted without the
assumption attached.

**Argoul's Table 1 is the independent check**, assembling human unbound
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

---

## 9. The inclusion rule (Stage 1)

One identity, as the binding side has one:

> **CL_int [µL/min/mg protein] = 1000 × V_max [nmol/min/mg protein] / K_m [µM]**

A source reaches the kinetics table only if it supplies what that identity
needs. Nothing else decides it — not how careful the paper is, not how many
compounds it covers, not whether it agrees with anyone.

| Tier | Meaning | Rows | Sources |
|---|---|---|---|
| **A** | reaches CL_int per mg protein | **52** | 8 |
| **B** | a clearance on a basis that will not join | 2 | Han 2008 only |
| **C** | affinity only — K_m, no clearance | 4 | Yang 2009/2010 only |
| **D** | no constant: substrate calls, uptake ratios, inhibition | 305 | 7 |

A, B and C live in `db/transporter_kinetics_working.csv`; **D is a separate
file**, `db/transporter_direction.csv`. That split is the rule made physical.
B and C stay beside A rather than in files of their own because they are six
rows between them, and because both are one measurement away from A (below).

**The addition: a directly reported clearance counts.** Kimura 2017 reports no
V_max, but prints "uptake clearance V_max/K_m = 55.0 µL/mg protein/min"
alongside K_m 8.3 µM. That is the quantity the identity produces, stated
outright, so the row is tier A. Its implied V_max, 55.0 × 8.3 = 456 pmol/mg
protein/min, sits inside the range Nakagawa (365–479) and Niu (192–674)
measured — a consistency check that supports admitting it. Note this also
settles the §6 item-8 question of whether Kimura's K_m of 8.3 µM, an order of
magnitude below every other K_m here, is "important or wrong": it is admitted,
flagged, and left for the write-up, because the clearance it implies is
unremarkable even though the K_m is not.

**The rejection: K_i and IC50 do not get in**, and the roadmap's original
candidate ("or a K_i with the probe substrate and its K_m") turns out to admit
nothing. An IC50 is as much a property of the probe substrate and its
concentration as of the inhibitor, so two IC50 values measured against
different probes are not comparable quantities. Converting one into a
probe-independent K_i needs the probe's own K_m, through Cheng–Prusoff.
**Of the 20 numeric IC50 values in this record, none carries it.** So the
clause admits zero rows and the quantity stays in tier D — which is the same
discipline the albumin and L-FABP sheets apply to probe-displacement
constants, and for the same reason.

While checking that, the `ic50` column turned out to have the defect `km` had:
28 populated cells, **20 numeric, 8 prose** ("significant inhibition only by
C6, C7, C8 at 10 µM"). The builder's `num()` guard already catches these.

**A substrate call is not a weak K_m.** Argoul 2026 contributes 120 rows and
no constants. That is a real measurement answering a real question — which
transporters touch which compounds, in which species — but no arithmetic turns
it into an affinity or a clearance, and keeping it in the same sheet as the
constants is how, eventually, somebody averages the two. Ryu 2024 adds 109
more of the same kind. **Together they are two thirds of the compilation**,
which is why the split matters more here than it did for binding.

### What the screen says is worth asking for

The request list now derives from the data rather than from memory, and it is
short:

| Source | Needs | Rows gained |
|---|---|---|
| **Yang 2010** | its own V_max (the record holds OAT4 37 and URAT1 0.32 only as Louisse 2023's introduction quotes them) | C → A, **+3** |
| **Yang 2009** | a V_max for Oatp1a1/PFOA | C → A, **+1** |
| **Han 2008** | mg protein per 10⁶ rat hepatocytes, to convert 64.8 and 47.6 µL/min/10⁶ cells | B → A, **+2** |

**Tier A is 52 rows; all three would make it 58. Nothing else in the record is
one measurement away.** Two things follow. First, the entire tier-C stratum is
Yang 2009/2010 — the only affinity-only evidence in the compilation comes from
the two papers still outstanding, which is a sharper argument for getting them
than "they are cited a lot". Second, and less comfortably: **the kinetics
table is nearly at its ceiling on the current literature.** Growth past 58
rows needs new measurements, not better retrieval.

---

## 10. The direction join (Stage 2)

`db/transporter_localisation.csv` is the authority: **36 rows, one per
(transporter, species, tissue)**, because direction is not a property of the
transporter alone.

- **OATP2B1 is apical in human kidney**, where it reabsorbs, **and sinusoidal
  in human liver**, where it pulls PFAS out of blood. Same gene, opposite jobs.
- **Oatp1a1** is apical in rat kidney and sinusoidal in rat liver.
- **OAT2 is basolateral in human kidney and apical in rat kidney** — the same
  gene on opposite membranes, so identical transport would secrete in one
  species and reabsorb in the other.

Nor is it a property of the compound. Membrane and polarity are fixed; what
varies per compound is only whether it is a substrate.

The renal rows come from the state-of-the-science review **Niu 2023**
(*Environ Health Perspect* 131:076002, [10.1289/EHP11885](https://doi.org/10.1289/EHP11885),
PMID 37418334) — the review Argoul cites — which states the polarity outright:
"OAT1 and OAT3 reside in the basolateral membrane … and therefore participate
in PFAS renal secretion. In contrast, OAT4, URAT1, and OATP2B1 are involved in
PFAS renal reabsorption because they are expressed in the apical membrane."
Tissue roles beyond the kidney follow **OEHHA's PHG** (March 2024, §4.3 and
Appendix 6), which names three: renal secretion/reabsorption; hepatic re-uptake
of biliary PFAS via bile-acid transporters, giving enterohepatic circulation;
and placental OAT4.

### The result that matters most is an absence

**OAT4 is human-specific.** Niu 2023: *"OAT4 is a human-specific transporter
that was previously demonstrated to transport PFOA; thus, it may be one of the
transporters responsible for species-specific differences in PFAS
toxicokinetics."* EFSA states it the same way: *"OAT4 (not expressed in rats)
and urate transporter 1 (URAT1), but not Oatp1a2 (the closest orthologous of
rat Oatp1a1), were shown to mediate PFOA."*

So the rat has **no orthologue of the apical reabsorptive transporter that
carries six of the measured human clearances in this compilation** (Louisse
2023, 39–92 µM). And OATP1A2, the obvious human counterpart to rat Oatp1a1,
does not transport PFOA (Yang 2010; Argoul 2026 confirms across ten PFAS).

**This is a different kind of answer from §8, and both hold.** §8 showed
affinity (Nakagawa, within 2×) and substrate specificity (Argoul, same MW
threshold) are conserved — *where a transporter exists in both species*. §10
adds that the reabsorptive machinery itself is **not the same set of genes**.
A presence/absence difference is not an affinity difference, and it is a far
better candidate for a ~50× half-life gap than anything in §8.

### Coverage and conflicts

| | |
|---|---|
| rows joined to a membrane and net process | **343 / 363 (94.5%)** |
| direction disagreements between a paper and the authority | **0** |
| unresolved | 20, all composite |

**Zero conflicts is itself a finding.** The papers that disagree sharply about
*which compounds* a transporter carries (§4.3) do not disagree at all about
*where it sits or which way it moves them*. Substrate identity is contested;
localisation is not.

The 20 unresolved rows all name several isoforms at once or an unidentified
one — "OATPs (Caco-2)", "OATP1B1, OATP1B3, OATP2B1", "hepatic uptake
transporter(s), OATP-like", "OAT4 and URAT1". There is no single membrane to
assign, so they are reported rather than forced.

**The join found two holes in the authority table**, which is what it is for:
OSTα/β needed a liver row (Zhao's rows describe efflux "out of
enterocyte/hepatocyte"), and rat Ostα/β needed a kidney row (Argoul tested it
renally). Both added and marked as such.

### CL_int regrouped on direction

| Net process | n | median | range |
|---|---|---|---|
| renal reabsorption (step 1) | 14 | **74.3** | 2.0 – 175.6 |
| hepatic uptake | 14 | 30.4 | 3.1 – 236.2 |
| renal secretion (step 1) | 17 | 18.9 | 4.6 – 60.2 |
| renal secretion (step 2, efflux) | 6 | 219.6 | 24.2 – 279.1 |

**The efflux row must not be read as a rate comparison.** Those six are P-gp
and BCRP from Niu 2026 in inside-out membrane vesicles, and that paper says
outright that "the specific rates of PFAS uptake cannot be directly compared
between these transporters" — different concentrations, and a vesicle
background of 1000–2000 pmol/mg against 100–400 in cells. They are not
evidence that apical efflux is an order of magnitude faster than basolateral
uptake.

On comparable whole-cell assays only:

| Effect on half-life | n | median CL_int |
|---|---|---|
| **lengthens** (reabsorption + hepatic uptake) | 28 | **36.1** |
| **shortens** (renal secretion) | 17 | **18.9** |

**~1.9×, not the ~3× §1 reported** from the free-text grouping. The direction
of the effect survives; its size roughly halves. Still a hypothesis for
Stage 3, for the REF/RAF reason in §2.3.

---

## 11. Measured expression, and the species answer

`db/transporter_expression.csv`, 40 rows from Thakur 2024 (*Clin Pharmacol
Ther* 116:235, [10.1002/cpt.3277](https://doi.org/10.1002/cpt.3277)), Uchida
2016 ([10.1002/jps.24645](https://doi.org/10.1002/jps.24645)) and Niu 2026.
Thakur measured 78 drug-metabolising enzymes and transporters in kidney
membrane fractions of human (n = 15), rat (n = 10) and mouse (n = 14), **both
sexes, one assay** — which is what makes the comparison legitimate.

### The secretory side is conserved; the reabsorptive side is not the same machinery

| Transporter | human | rat | mouse | |
|---|---|---|---|---|
| **Oat1** (basolateral, secretory) | 4.8 ± 1.2 | **32.8 ± 13.7** | 19.1 ± 7.1 | rat **6.8×** human |
| **Oat3** (basolateral, secretory) | 4.3 ± 1.1 | 6.3 ± 1.1 | 6.5 ± 3.6 | rat 1.5× human |
| **OAT4** (apical, reabsorptive) | **0.5 ± 0.2** | **not detected** | not detected | |
| **Oatp1a1** (apical, reabsorptive) | **not detected** | **1.9 ± 2.1** | 3.1 ± 2.8 | |
| URAT1 (apical) | 3.1 ± 1.2 | 2.2 ± 0.7 | 10.8 ± 1.1 | transports no PFAS in either |

pmol/mg membrane protein. **Each species reabsorbs with a transporter the
other does not have.** §10 reached that from Niu 2023 and EFSA; here it is
measured, in one assay, by proteomics — including Thakur's own listing of
Oat4 among proteins "absent in rats".

### The rat sex difference: a sign change, with a mechanism

> **"The abundances of Oatp1a1 and Oatp1a3 were 114- and 7-fold higher in
> male rats than in female rats, respectively (P < 0.05)."**

That single sentence explains the shape of the rat sex difference, not merely
its size. OEHHA PHG Table A6.4 has male rats reabsorbing 93.7% of filtered
PFOA while **female rats are in net tubular secretion**. A difference in
abundance alone cannot flip net handling from reabsorption to secretion; a
transporter that is effectively absent in one sex can. It also fits Cheng
2006 — "in kidney, gender-divergent Oatp expression is exclusively caused by
stimulation by androgens" — and Kudo 2002, where castration raised renal
clearance 14-fold.

And the human control: **no human renal DMET protein showed any sex
difference** (male n = 11, female n = 4), against 28 in rat and 40 in mouse.
Humans show no PFAS half-life sex difference either.

### Two §4.3 conflicts dissolve

**Asbt and OSTα/β were not detected in human *or* rat kidney** (mouse has
both). Zhao 2015 found transport by them in transfected cells; Argoul 2026
found none. **Both are right.** A transfected-cell result states what a
protein *can* carry, not whether it is present in the kidney to do it. That
retires two of the seven disagreements in §4.3, and it retires them as a
category error rather than as a winner-takes-all.

### A warning for Stage 3, from the same table

**Human renal BCRP is not detected** — by Thakur and independently by Niu
2026 ("below the detection limit"). Yet Niu 2026's BCRP/PFOA carries the
**highest CL_int in this entire compilation, 279 µL/min/mg**. Scaled by renal
expression it contributes essentially nothing in human kidney. This is the
clearest possible demonstration that Stage 3 must not rank transporters by
CL_int alone — and the vesicle caveat in §10 was already pointing at the same
row.

### The sources agree where they overlap

| | Thakur 2024 | Niu 2026 |
|---|---|---|
| human OAT1 | 4.8 | 5.3 pmol/mg |
| human P-gp | 1.2 | 2.1 pmol/mg |
| human BCRP | not detected | below detection |

Rat Oat1 is 32.8 pmol/mg in Thakur's total membrane fraction against Uchida's
10.6 fmol/µg in a purified basolateral fraction — different preparations,
same order.

### A tooling trap worth recording

`pdftotext` silently drops the micro sign from the Uchida PDF, turning every
"fmol/µg protein" into "fmol/mg" — in the table header, the caption *and* the
abstract alike, so the extracted text is internally consistent and gives no
hint of the loss. Reading the page as an image settles it. An audit of all 68
PDFs in this session found Uchida is **the only one affected**. The numbers
published earlier in this document came from the PubMed abstract and were
right; the PDF nearly overwrote them with values 1000× too small. **A quantity
extracted by `pdftotext` alone is not safe until the micro signs are
accounted for.**

---

## 12. Scaling factors, and the albumin problem (pre-Stage 3)

`db/ivive_scaling_factors.csv`. Two ways to turn a CL_int measured in
transfected cells into a clearance in a kidney, and **they disagree**:

| | REF (HEK293 ÷ kidney) | RAF (probe substrate) | apart |
|---|---|---|---|
| OAT1 | 0.16 | 0.64 (tenofovir) | 4× |
| OAT2 | 0.15 | 7.3 (acyclovir, ganciclovir) | 49× |
| **OAT3** | **0.37** | **4.1** (oseltamivir acid, benzylpenicillin) | **11×** |

Kumar 2021 ([10.1124/dmd.121.000367](https://doi.org/10.1124/dmd.121.000367),
PMC8232064) compared both over 26 OAT substrates and found **RAF the better
predictor** — 65% and 81% of plasma secretory clearances within 2- and 3-fold,
against 50% and 69% for REF. So the method with no proteomics in it wins, and
the choice changes a scaled OAT3 clearance by an order of magnitude.

**Two traps recorded in the table.** First, **REF's direction is not a safe
convention**: Kumar defines it as cells ÷ kidney (so REF < 1, and you divide),
while Tan 2024 reports REF 8.0 for OAT1 while stating its system expressed
*less* than in vivo — only consistent if theirs is kidney ÷ cells. Getting it
backwards is a factor of 40 on OAT1. Second, **REF is cell-line specific**:
Antonescu 2020 measured OAT1 at 1.32 fmol/µg in HEK293-Flp-In against 5.33 in
renal cortex, a ratio of 0.25 where Kumar's HEK293 gives 0.16.

### The albumin problem is bigger than the REF/RAF choice

Tan 2024 ([10.1021/acs.molpharmaceut.4c00504](https://doi.org/10.1021/acs.molpharmaceut.4c00504),
PMC11372837): with 4% human serum albumin present **and the unbound substrate
concentration held constant**, unbound intrinsic active uptake clearance rose
**3.5- to 122-fold for OAT1** and **up to 28-fold for OAT3**. The effect
appeared only for highly protein-bound substrates — f_u between **0.002 and
0.18** — and was absent for weakly bound ones. Mechanism: facilitated
dissociation, where albumin releases its ligand at the cell surface so the
transporter sees more than the bulk free fraction. The authors recommend
including HSA in routine OAT1/3 assays.

**PFOA's plasma f_u is 0.0017** (Argoul Table 1) — at or below the bottom of
that range, and PFOS, PFNA and PFHxS are lower still. PFAS are the extreme
case of the class this effect was demonstrated on.

**And no PFAS assay in this database is recorded as containing protein.** Of
the 52 tier-A rows across six studies, not one mentions albumin, serum or BSA.
Where the buffer is stated outright it is protein-free (Argoul: HBSS, and
Krebs-Henseleit for rat Oat1); for the rest this is an absence of mention
rather than a verified composition.

So the tier-A clearances are plausibly **lower bounds**, and understated
*most* for the longest-chain, most tightly bound compounds — the same
direction as the chain-length trend this compilation is trying to explain.
That is not a reason to discard the table; it is a caveat to state, and a
testable prediction: run one PFAS uptake assay with and without HSA.

**It also joins the strands.** The albumin association constants this project
has spent its effort on are exactly the input a facilitated-dissociation
correction needs. The binding sheet stops being only a distribution story.
