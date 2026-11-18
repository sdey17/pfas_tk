# Transporters: what exists, what is wrong with it, and what to do next

Working document, written 2026-11-18 from a direct survey of the repository
rather than from the earlier narrative. Figures in §1 were counted from the
files; where the existing retrieval log disagrees, the files win and the
disagreement is named.

**Stage 0 is complete** (§7), and doing it corrected three claims in the first
draft of this document. Those corrections are kept in place rather than edited
away, marked ~~struck~~, because two of them are the kind of error that would
otherwise be repeated: a second-hand table counted as new data, and a count I
got wrong in my own survey.

---

## 1. What already exists

| Artefact | Size | State |
|---|---|---|
| `db/transporter_kinetics.csv` | **230 rows, 27 sources**, 48 transporter entities, 61 chemical cells | the record — everything read from a primary, including rows that are not transporter kinetics |
| `db/transporter_kinetics_working.csv` | **241 rows, 15 sources, 34 transporters, 19 compounds** | gated, normalised, split; built by `scripts/build_transporter_working_table.py`, never edited by hand |
| `db/transporter_mechanism.csv` | 83 rows | direction, sex and hormone regulation; Km on 24 |
| `db/primary_2026/*.csv` | 24 extracts | per-paper tables read from primaries |
| `papers/*.txt` | 29 transporter papers | **text only — no PDFs**, see §5 |
| `papers/SOURCES_transporters.md` | retrieval log | detailed, and now partly stale |
| `report/REPORT.md` §7 | narrative | direction-first, with three open conflicts |
| `scripts/saturation_margin_km_vs_kt.py` | analysis | the Km-versus-KT conflict, §4.1 below |

**Fill rates**, record and working table, after Stage 0:

| Field | Record | Working table | Comment |
|---|---|---|---|
| `uptake_ratio` | 138 / 230 | — | substrate / not-substrate calls |
| `km` | **58 / 230** | **59 / 241** | all in µM or µmol/L; essentially already harmonised |
| `vmax` | 45 / 230 | **38 / 241** | five unit spellings in the record, one basis in the working table (§3.2) |
| `ic50` | 28 / 230 | — | inhibition, a different quantity again |
| **`clint_uL_min_mg`** | — | **38 / 241** | V_max/K_m, newly computed; see below |

So two thirds is still qualitative, and Stage 0 did **not** double the K_m
count as this document predicted: it moved it 52 → 58, because the ten Louisse
2024 values turned out to be six (§3.1). The quantitative core is about 60
rows, which is the size of a single good paper's output, not a compilation.

What Stage 0 did deliver is the **38-row CL_int column**, which was meant to be
Stage 3. It spans 2 to 236 µL/min/mg protein, and the highest value in it is
not a renal transporter at all but human NTCP with PFOS (236), with hepatic
OATP2B1/PFOS next (99) — a first indication that on a per-mg-protein basis the
hepatic uptake side is not the minor partner the renal framing assumes.

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

**4.3 Argoul 2026 against everyone.** Reports PFDA, PFOS and PFDS as substrates
of none of twelve transporters, against both Louisse and Weaver. Abstract only.

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
| 1 | **Argoul 2026**, *Arch Toxicol*, 10.1007/s00204-026-04508-7 | The only direct human-vs-rat comparison, 8 human + 4 rat transporters, 10 PFAS. Settles §4.3 and is the species-difference backbone |
| 2 | **Ebert 2020**, *ES&T* 54:5051, 10.1021/acs.est.0c00175 | Argues measured membrane permeabilities explain cellular uptake "with no need to postulate active uptake". If right, it reframes the whole strand |
| 3 | **Niu 2026** figures, 10.1021/envhealth.6c00088 | Per-congener K_m/V_max for BCRP, P-gp, OAT1 — the efflux side, currently ranges only |
| 4 | **Yang 2010**, *Toxicol Sci* 117:294, 10.1093/toxsci/kfq219 | K_m is in hand; the V_max values are second-hand from Louisse 2023's introduction. Also the URAT1 chloride condition (§4.2) |
| 5 | **Yang 2009**, *Toxicol Lett* 190:163, 10.1016/j.toxlet.2009.07.011 | Oatp1a1 K_m and the C6–C10 K_i,app series, currently quoted via Weaver 2010's Discussion |
| 6 | **Nakagawa 2008**, 10.1111/j.1742-7843.2007.00155.x | The four human/rat OAT1 and OAT3 K_m values everything else cites |
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
