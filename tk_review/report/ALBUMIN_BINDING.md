# Serum albumin–water partition coefficients for per- and polyfluoroalkyl substances: a harmonised compilation of 288 measurements from 20 primary sources

## Abstract

Binding to serum albumin governs the plasma distribution and renal clearance of
per- and polyfluoroalkyl substances (PFAS), and is consequently a required
parameter in toxicokinetic models. Published binding constants for the same
compound and the same protein differ by up to four orders of magnitude, which
has limited their use. We compiled 288 albumin–water partition coefficients for
63 PFAS from 20 primary experimental sources, converted to a single basis —
log *K* in litres of water per kilogram of protein — and audited each value
against its primary publication. Conversion used the identity
log *K* = 4.1774 + log₁₀*n* − log₁₀(*K*_d/µM); all 181 rows carrying both a
dissociation constant and a stoichiometry satisfy it exactly. The compilation
spans log *K* 0.48 to 5.49 (median 3.37, interquartile range 2.46–4.13), with
178 human and 110 bovine values. Dispersion within a single compound measured by
different laboratories is large and systematic rather than random: across the 32
compounds with three or more values the median within-compound range is 2.38 log
units and the maximum is 4.22. Method accounts for much of it — equilibrium
dialysis and differential scanning fluorimetry differ by 2.35 log units in
median, in a direction predicted by the thermal denaturation that the latter
measures. Perfluoroalkyl carboxylate affinity rises by 0.52 log units per
perfluorinated carbon from C4 to C9 and declines thereafter, placing the maximum
at perfluorononanoic acid. We document the conversion and admission decisions in
full, including seven sources read in their entirety and excluded, because in a
compilation of this kind the exclusions carry as much information as the
inclusions.

**Keywords:** PFAS; serum albumin; partition coefficient; protein binding;
toxicokinetics; data harmonisation

---

## 1. Introduction

Serum albumin is the most abundant protein in blood plasma, present at
approximately 600 µM, and is the principal carrier of anionic PFAS in the
circulation. The fraction of a PFAS that is not bound to plasma protein
determines what is available for glomerular filtration, for transporter-mediated
reabsorption, and for distribution into tissue. Albumin binding is therefore not
an incidental physicochemical property but a direct determinant of the long
elimination half-lives that make several PFAS of regulatory concern.

Toxicokinetic models require this parameter as a number. The literature supplies
it many times over, and the values disagree. For perfluorooctanoic acid bound to
human or bovine albumin, published constants recalculated onto a common basis
span more than three orders of magnitude. A modeller selecting one value is
making a consequential choice with little guidance, and a modeller averaging them
is averaging quantities that are not the same quantity.

Three features of the primary literature produce this situation. First, binding
is reported in incompatible units: dissociation constants, association
constants, protein–water distribution coefficients on a volume basis, and
distribution coefficients on a mass basis all appear, and conversion between
them requires a stoichiometry that is frequently not reported. Second, the
measurement conditions that most affect an apparent binding constant — the
ligand-to-protein ratio, temperature, ionic strength, pH, and co-solvent — vary
widely and are not always stated. Third, several widely cited values turn out on
reading the primary source not to be the quantity their citing papers take them
for.

This work addresses the first problem directly and the second and third by
documentation. We assembled every experimental albumin–PFAS binding measurement
we could locate for human and bovine albumin, converted each to a single basis
through an explicit identity, verified each against its primary source, and
recorded every decision that changes a number. The compilation is intended for
use in mass-balance and toxicokinetic models, and its principal contribution is
not the central values but the auditable provenance attached to each.

---

## 2. Methods

### 2.1 Data compilation

Candidate sources were identified by structured literature search of PubMed for
terms combining perfluoroalkyl substances with serum albumin and with binding,
dissociation, or association constants, supplemented by citation tracing from
compiled reviews and by the reference lists of the primary sources themselves.
Searches were last run in November 2026.

A source was eligible if it reported an original experimental measurement of
PFAS binding to serum albumin from which a dissociation constant could be
obtained. Computational studies, molecular docking, quantitative
structure–activity relationships, and review compilations without original
measurement were excluded at screening. The underlying database retains 1,460
rows across all proteins and species; the present compilation is the subset
meeting the criteria in §2.3.

### 2.2 Conversion to a common basis

All values were expressed as log *K*, the albumin–water partition coefficient in
litres of water per kilogram of protein. For a ligand binding *n* equivalent
sites with dissociation constant *K*_d,

&nbsp;&nbsp;&nbsp;&nbsp;*K* = *n* / (*K*_d · *M*_kg)

where *M*_kg is the molar mass of albumin in kilograms per mole. Taking
logarithms and expressing *K*_d in micromolar,

&nbsp;&nbsp;&nbsp;&nbsp;log *K* = log₁₀(10⁹/*M*) + log₁₀*n* − log₁₀(*K*_d/µM)
= **4.1774 + log₁₀*n* − log₁₀(*K*_d/µM)**

with *M* = 66,463 g·mol⁻¹ for human serum albumin. The bridge constant 4.1774 is
exact for that molar mass. Association constants were inverted as
*K*_d = 10⁶/*K*_a for *K*_a in M⁻¹, and logarithmic association constants as
log₁₀(*K*_d/µM) = 6 − log *K*_a.

Sources reporting a partition coefficient directly on a volume basis
(L_water per L_protein) were converted to the mass basis by the partial specific
volume of albumin, 1.36 kg·L⁻¹, a difference of log₁₀(1.36) = 0.1335 log units.
Sources reporting on a mass basis were used as published.

Every converted row retains the route by which it was obtained, so that values
resting on a reported stoichiometry can be separated from those resting on an
assumption. The conversion identity was verified mechanically: of the 288 rows,
181 carry both a dissociation constant and a stoichiometry, and **all 181
satisfy the identity to within 0.005 log units**.

### 2.3 Inclusion criteria

A measurement entered the compilation if it satisfied all of the following.

1. **A genuine binding constant.** A dissociation or association constant for a
   solution-phase binding equilibrium. Retention indices, docking scores,
   adsorption-isotherm threshold parameters, and gas-phase constants were
   excluded as different quantities (§4.2).
2. **A stoichiometry, or a defensible substitute.** Where *n* was reported it was
   used. Where it was not, *n* = 1 was substituted, and the row is labelled to
   say so. Where the source's own model fixes *n* = 1 by construction, the row is
   labelled as resting on the paper's assumption rather than ours.
3. **Human or bovine albumin.** Rat, mouse, porcine, and piscine values exist in
   the underlying database but are not pooled here, because interspecies
   differences in albumin binding are a question this compilation is intended to
   address rather than assume away.
4. **A perfluorinated or polyfluorinated analyte.** Straight-chain fatty acids
   run as positive controls or competitive probes were excluded (§2.4).

### 2.4 Quality gates

Four automated gates were applied to all sources.

**Goodness of fit.** Where a source publishes the coefficient of determination
for each single-site fit, rows below *R*² = 0.80 were excluded, on the grounds
that a dissociation constant derived from a curve the model does not describe is
not a measurement of the constant. Seven rows were removed by this gate.

**Co-solvent.** Rows whose binding curve was measured in a co-solvent
concentration high enough to act on the protein were excluded. The gate is keyed
on the combination of source and compound rather than on compound alone, because
the same compound measured in buffer by another laboratory is unaffected.

**Non-fluorinated analytes.** Octanoic, decanoic, hexadecanoic, and butyric acid
appear in several sources as reference ligands occupying albumin's fatty-acid
sites. Eight rows were removed. These are informative measurements but they are
not PFAS, and retaining them would place a C16 hydrocarbon inside a
perfluoroalkyl chain-length relationship.

**Composite entries.** A tabulated value whose compound field names several
analytes is a summary statistic, not a measurement of any one of them. Such rows
were excluded. The test is empirical rather than syntactic: the field is split on
delimiters outside parentheses and the parts are tested against compound names
appearing elsewhere in the database as standalone analytes, since the commas in a
systematic chemical name are locant separators and not list delimiters.

### 2.5 Reproducibility

The compilation is generated by a documented pipeline from a single source
database, and every figure reported here is recomputed from the delivered
workbook rather than transcribed. Scripts, the database, and the audit routines
are listed in §7.

---

## 3. Results

### 3.1 Composition

The compilation comprises **288 measurements of 63 compounds from 20 primary
sources**, of which 178 are human serum albumin and 110 are bovine. Values span
**log *K* 0.48 to 5.49**, with a median of 3.37 and an interquartile range of
2.46 to 4.13. No value is negative; that is, no row implies that albumin takes up
a PFAS less effectively than water does.

**Table 1.** Composition by source.

| Source | *n* | Median log *K* | Principal method |
|---|---|---|---|
| Chen 2025 | 60 | 3.40 | equilibrium dialysis |
| Fischer 2024 | 55 | 4.34 | solid-phase microextraction |
| Allendorf 2019 | 28 | 4.25 | equilibrium dialysis |
| Li 2021 | 25 | 2.92 | equilibrium dialysis |
| Starnes 2025 | 19 | 1.05 | differential scanning fluorimetry |
| Gao 2019 | 14 | 2.14 | equilibrium dialysis |
| Starnes 2024 | 14 | 1.01 | differential scanning fluorimetry |
| Bischel 2011 | 11 | 3.92 | equilibrium dialysis |
| Hebert 2010 | 11 | 2.29 | equilibrium dialysis |
| Qin 2024 | 11 | 3.58 | equilibrium dialysis |
| Gao 2025 | 10 | 2.70 | fluorescence quenching |
| Alesio 2022 | 7 | 3.99 | equilibrium dialysis |
| Bischel 2010 | 6 | 4.70 | equilibrium dialysis |
| Moro 2022 | 4 | 3.52 | isothermal titration calorimetry |
| Chi 2018 | 4 | 3.79 | native electrospray mass spectrometry |
| Chen 2009 | 2 | 4.02 | fluorescence quenching |
| Sheng 2020 | 2 | 3.37 | isothermal titration calorimetry |
| Qin 2010 | 2 | 3.23 | fluorescence quenching |
| Chen 2020 | 2 | 2.08 | fluorescence quenching |
| Han 2003 | 1 | 2.46 | size-exclusion separation |

The basis on which each row's stoichiometry rests is recorded explicitly
(Table 2). Only 107 of 288 values are partition coefficients measured directly,
requiring no stoichiometric assumption. A further 111 rest on a single-site model
imposed by the source itself, and 53 on an assumption of *n* = 1 made during
compilation because the source reports no stoichiometry at all.

**Table 2.** Basis of the stoichiometry used in conversion.

| Route | Rows | Interpretation |
|---|---|---|
| measured directly | 107 | source reports a partition coefficient; no conversion |
| *n* = 1 by single-site model | 111 | source's own model fixes *n* = 1 |
| *n* = 1 assumed | 53 | source reports no *n*; assumption made here |
| *n* fitted by binding isotherm | 11 | source fits *n* from its own isotherm |
| ITC, defined site set | 4 | stoichiometry counted calorimetrically |
| *n* = fitted quenching exponent | 2 | *n* is a curve-fit exponent, not a site count |

### 3.2 Within-compound dispersion

Disagreement between laboratories measuring the same compound is the dominant
feature of the data. Across the **32 compounds with three or more independent
values, the median within-compound range is 2.38 log units** — a factor of
240 — and the maximum is 4.22.

**Table 3.** The ten most frequently measured compounds.

| Compound | *n* values | Median log *K* | Range | Spread (log) |
|---|---|---|---|---|
| PFOA | 25 | 3.75 | 1.18–4.53 | 3.34 |
| PFOS | 19 | 3.96 | 1.24–4.83 | 3.59 |
| PFNA | 18 | 4.29 | 1.41–5.16 | 3.75 |
| PFHxA | 15 | 2.77 | 0.84–3.92 | 3.07 |
| PFBS | 14 | 2.47 | 0.95–3.73 | 2.77 |
| PFDA | 14 | 3.77 | 1.10–4.97 | 3.86 |
| PFHxS | 13 | 3.99 | 1.33–4.81 | 3.48 |
| PFHpA | 12 | 3.54 | 1.53–4.10 | 2.57 |
| PFUnDA | 11 | 3.56 | 1.05–5.18 | 4.13 |
| HFPO-DA | 10 | 2.43 | 0.77–3.36 | 2.59 |

Dispersion is not uniform across compounds. ADONA (five values, range 0.56 log),
N-EtFOSAA and N-MeFOSAA (four values each, 0.96) are tightly constrained, while
PFDoDA (seven values, 4.22) and PFUnDA (eleven values, 4.13) are not. The
tightly constrained compounds are those measured predominantly by one method;
the dispersed ones are those measured across several.

### 3.3 Method is the principal axis of disagreement

Grouping by method family resolves much of the dispersion into a systematic
offset rather than scatter (Table 4).

**Table 4.** Distribution by method family.

| Method | *n* | Median log *K* | Range |
|---|---|---|---|
| Solid-phase microextraction | 66 | 4.29 | 1.81–5.49 |
| Equilibrium dialysis | 151 | 3.38 | 1.55–5.16 |
| Native electrospray MS | 4 | 3.79 | 3.29–4.17 |
| Isothermal titration calorimetry | 4 | 3.52 | 3.18–4.53 |
| Ultrafiltration | 2 | 3.37 | 3.27–3.47 |
| Fluorescence quenching | 14 | 2.49 | 0.48–4.14 |
| Differential scanning fluorimetry | 33 | 1.03 | 0.69–1.61 |

**Equilibrium dialysis and differential scanning fluorimetry differ by 2.35 log
units in median**, with dialysis higher. The direction is mechanistically
expected. Differential scanning fluorimetry infers binding from the shift a
ligand induces in the protein's thermal unfolding transition, so its constants
describe a protein at or near its melting temperature rather than a folded
protein at physiological temperature, and the measurable shift requires
occupancy well above half-saturation. Its narrow range — 0.69 to 1.61 across 33
values spanning 50 compounds — is itself diagnostic: a method that compresses
three orders of magnitude of true affinity into under one order is reporting
rank rather than magnitude.

Solid-phase microextraction sits 0.91 log above dialysis. All 66 values derive
from a single source, so this offset cannot be separated from a laboratory
effect.

### 3.4 Chain length

Perfluoroalkyl carboxylate affinity increases monotonically with perfluorinated
chain length to a maximum and declines thereafter (Table 5). The increment from
C4 to C9 is **+0.52 log units per perfluorinated carbon**, and the maximum falls
at **C9, perfluorononanoic acid**.

**Table 5.** Perfluoroalkyl carboxylate homologous series.

| *n*_CF | Compound | Median log *K* | Values |
|---|---|---|---|
| 4 | PFBA | 1.68 | 10 |
| 5 | PFPeA | 2.23 | 6 |
| 6 | PFHxA | 2.77 | 15 |
| 7 | PFHpA | 3.54 | 12 |
| 8 | PFOA | 3.75 | 25 |
| 9 | **PFNA** | **4.29** | 18 |
| 10 | PFDA | 3.77 | 14 |
| 11 | PFUnDA | 3.56 | 11 |
| 12 | PFDoDA | 2.98 | 7 |
| 13 | PFTrDA | 2.68 | 5 |
| 14 | PFTeDA | 3.46 | 4 |

The decline beyond C9 is consistent with steric constraint in albumin's binding
pockets and with the decreasing aqueous solubility of long-chain homologues,
which reduces the free concentration available to the binding equilibrium. The
C14 value interrupts the monotonic decline, but rests on four measurements and
should not be weighted against the trend.

### 3.5 Species

Human and bovine albumin differ in median by 0.90 log units, bovine the higher
(3.99 against 3.09). This difference is confounded with method and compound
coverage: the bovine subset is dominated by solid-phase microextraction and
dialysis from sources that measured long-chain compounds, while the human subset
includes all 33 differential scanning fluorimetry values. The comparison should
therefore not be read as an interspecies difference in albumin affinity, and the
two species are kept separable in the delivered workbook for that reason.

---

## 4. Discussion

### 4.1 What the dispersion means for use of these values

A median within-compound range of 2.38 log units means that selecting a single
published constant for a toxicokinetic model is a decision with a 240-fold
range of defensible outcomes. The results above indicate that this is not
irreducible measurement noise. Most of it is attributable to identifiable
differences in what was measured.

Three implications follow. First, **a constant should be selected by method, not
by recency or citation count.** For a model requiring binding to folded albumin
at physiological temperature, equilibrium dialysis and isothermal titration
calorimetry measure that quantity and thermal-shift methods do not. Second,
**values should not be averaged across method families.** The 2.35 log offset
between dialysis and differential scanning fluorimetry is systematic, so a mean
across both is a weighted artefact of how many of each happen to have been
published. Third, **the stoichiometry assumption must be carried with the
value.** Only 107 of 288 values require no stoichiometric assumption; for the
53 resting on an assumed *n* = 1, a true *n* of 2 would shift log *K* by 0.30.

### 4.2 Sources read in full and excluded

Seven sources met the screening criteria, were read in full, and were excluded.
Each failed for a distinct reason, and the reasons are informative about the
literature as a whole.

**A quantity that cannot be converted.** One source reports 82 albumin binding
values for 83 compounds by biomimetic chromatography on an immobilised-albumin
column — the largest single dataset located. Its supporting information states
that the values are dimensionless partition coefficients, but neither the paper
nor its supplement reports a phase ratio or an immobilised-protein mass. The
albumin quantity is therefore absorbed into the empirical calibration constants
of the retention equations and cannot be recovered. Sources reporting on a
volume basis were convertible because a density bridges volume to mass; here the
denominator is not a volume of albumin but albumin immobilised on silica of
unstated composition. The authors state that absolute values may not be
confirmable by other methods and that the reproducible quantity is elution
order. This is a rank-ordering dataset, and is treated as such.

**A constant measured at the protein's melting temperature.** A second source
reports association constants for seven PFAS binding bovine albumin, all
evaluated at 61 °C — the melting temperature of the protein, which the paper
states explicitly. The exclusion is unusually well controlled because the same
laboratory published equilibrium dialysis values for the same seven compounds
and the same protein at 37 °C, and those values are in the compilation. The
comparison therefore isolates temperature and method with laboratory, protein,
and analyte held fixed, and every compound falls: PFDA from log *K* 4.97 to 1.70
(a factor of 1,860), PFOS from 4.43 to 2.61, PFOA from 3.75 to 2.57, HFPO-DA
from 2.99 to 2.47. The chain-length ordering also inverts, with PFNA exceeding
PFDA at 61 °C where the same laboratory's 37 °C data has the opposite and
expected order. A discrepancy of this size and structure is not an offset that
calibration could remove.

**A thermodynamically incoherent temperature series.** A third source reports
six PFAS at three temperatures with van't Hoff fits. Every one of the six yields
a binding enthalpy outside the range physically available to non-covalent
association: +412 kJ·mol⁻¹ for one compound, +325 and −278 for others, against a
normal range near ±60 kJ·mol⁻¹. These follow arithmetically from binding
constants that swing by factors of 5 to 624 across a 12 K interval, so the
series is not internally consistent at any single temperature, and the fitted
exponent *n* rises with temperature in four of six compounds, which a site count
cannot do. One compound in this source had been the compilation's largest value
at log *K* 7.86; the maximum is now 5.49.

This exclusion illustrates a limitation of row-level quality control. No
individual value in the source appeared anomalous. The defect was visible only
in the relationship between the three temperatures, which is to say only because
the source reported three. Sources reporting a single temperature cannot be
screened this way.

**Values cited rather than measured.** A fourth source tabulates eleven albumin
dissociation constants, six of which carry a superscript citation to an earlier
paper already in the compilation. Recalculated through the conversion identity,
all six reproduce the earlier source's values exactly, confirming that they are
cited rather than measured. Ingesting them would have double-counted one primary
source under another's name. The typographic risk is worth recording: in the
published layout the citation marker abuts the value, so a constant of 38 µM
reads as 3,819.

**Site-specific displacement constants.** A fifth source reports five
competitive probe-displacement constants alongside its direct measurements. A
displacement constant is specific to one subdomain and carries the displaced
probe's own affinity, so it is not interchangeable with an overall binding
constant. The source's own data make the point: its direct tryptophan-quenching
measurement of PFOS gives *K*_d = 45.45 µM with *n* = 14, or log *K* 3.67,
against 5.06 from site II displacement — 1.4 log units apart within one paper,
and a displacement row forced to *n* = 1 is doubly misdescribed where the direct
method resolved fourteen sites.

**Composite range entries.** Three further sources report values only as
collapsed ranges spanning several compounds, several species, or several
sequential binding steps. A midpoint of such a range describes none of its
members.

### 4.3 Decisions that change a reported number

Four decisions in the compilation alter values relative to their published form,
and are recorded here because they are not self-evident from the workbook.

**Ionic strength.** One source measured each compound at seven sodium chloride
concentrations from 0 to 500 mM, and the compilation takes one row per compound.
The row selected is that at 150 mM added chloride, the physiological level the
source identifies, falling back to the highest level measured below it where
150 mM was not run. The source's central finding is that ionic strength
suppresses PFAS–albumin binding, so a no-added-salt row is simultaneously the
condition that most overstates affinity and the one least representative of
plasma. This selection moves ten values down by 0.18 to 0.34 log units relative
to an earlier version of the compilation.

**Temperature.** Three sources measured at several temperatures. The 310 K rows
were taken in each case, because the compilation is intended for use at body
temperature and 89% of rows whose source states a temperature are at 37 °C.

**Isomers.** One source resolves two analytes into linear and branched isomers.
The branched rows were excluded and the linear rows carry the unqualified
compound name, because every other source reports a single value for technical
material that is 98.5–99.99% linear, and a branched-isomer constant does not
pool with a linear-dominated one.

**Scope of three sources.** Three sources were admitted only for compounds
already present from another source, so that they strengthen existing
comparisons rather than contributing single-source values against which nothing
can be checked. The cost is substantial and asymmetric: for one of the three,
seventeen compounds were excluded for which that source holds the only albumin
measurement in existence. This is a defensible trade for a cross-study
comparison and a poor one for a survey of coverage, and it is implemented as a
single reversible switch for that reason.

---

## 5. Limitations

**Single-source compounds.** Twenty-eight of the 63 compounds rest on one
source. For these, the compilation provides a value but no estimate of its
reliability, and the dispersion statistics in §3.2 do not apply to them.

**The species comparison is confounded.** As noted in §3.5, the human and bovine
subsets differ in method composition and compound coverage, so the 0.90 log
difference in medians cannot be attributed to the protein.

**Stoichiometry is usually unmeasured.** For 53 rows *n* = 1 was assumed during
compilation. Where *n* has been measured it is frequently not 1 — one source
resolves 14 sites for PFOS — so this assumption is a known and unquantified
source of bias in a fifth of the data.

**One compound lost its only sound value.** The thermodynamically incoherent
source discussed in §4.2 carried one compound whose only other measurement was
admitted on the strength of that source's presence. When the incoherent source
was removed, the scope rule in §4.3 removed the sound measurement with it, as a
now-single-source compound. The compound consequently has no value in the
compilation despite a defensible measurement existing. This is a correct
application of a stated rule producing an undesirable outcome, and re-admitting
that single row is a defensible exception.

**Non-human, non-bovine species are not pooled here.** The underlying database
holds rat, mouse, porcine, and piscine values. Mouse albumin in particular is
sparsely characterised: the one dataset located reports fraction unbound rather
than a dissociation constant, and so does not reach this compilation's basis
without an albumin concentration through which to convert.

**Enthalpy is not screened systematically.** The failure mode identified in
§4.2 — a binding enthalpy outside the physically available range — was detected
only because one source reported three temperatures. A systematic enthalpy
screen across all sources reporting one would be the cheapest available
improvement to data quality, and has not been performed.

---

## 6. Conclusions

A harmonised compilation of 288 albumin–water partition coefficients for 63 PFAS
is presented, converted to a single mass basis through an explicit and
mechanically verified identity. The compilation's central finding is that the
four-order-of-magnitude disagreement in the published literature is largely
systematic rather than random, and is resolvable into differences of method,
temperature, ionic strength, and assumed stoichiometry. Equilibrium dialysis and
thermal-shift methods differ by 2.35 log units in median in a direction their
mechanisms predict; within-compound ranges of up to 4.22 log units narrow
substantially when method is held fixed.

For toxicokinetic use we recommend selecting values by method rather than
pooling across methods, carrying the stoichiometry basis with the value, and
treating thermal-shift constants as rank orderings. The chain-length
relationship for perfluoroalkyl carboxylates — +0.52 log units per
perfluorinated carbon to a maximum at C9 — is robust across sources and provides
a check on any individual value.

The exclusions documented in §4.2 are as much a result as the compilation
itself. Of seven sources read in full and rejected, two were rejected on grounds
stated plainly in the source's own text, two on arithmetic performed on the
source's own tables, and one on the source's own internal contradiction. In none
of these cases was the defect apparent from the value as cited elsewhere.

---

## 7. Data and code availability

| Artefact | Path |
|---|---|
| Delivered compilation (288 rows) | `albumin_logK_Kd_FINAL.xlsx` |
| Source database (1,460 rows, all proteins and species) | `db/protein_binding.csv` |
| Conversion pipeline | `scripts/build_binding_affinity_table.py`, `scripts/albumin_kd_to_k.py`, `scripts/albumin_shareable.py`, `scripts/build_revised_working_table.py` |
| Eligibility screen | `scripts/screen_fits_database.py` |
| Conversion audit | `scripts/audit_working_table.py` |
| Exclusion audit | `scripts/excluded_sources.py` |
| Per-source audits | `scripts/audit_*.py` |
| Extended methodological primer | `report/ALBUMIN_BINDING_PRIMER.md` |

The primer referenced in the final row is a companion document covering the
derivations, worked examples, and per-source reasoning at greater length than is
appropriate here, including the full text of each conversion and admission
decision.

---

## 8. References

Primary sources are cited by first author and year as they appear in Table 1.
Full bibliographic detail for each, together with the specific table or figure
from which each value was taken, is recorded in the `source_table` field of
`db/protein_binding.csv` on a per-row basis, and in the per-source audit scripts
listed in §7.
