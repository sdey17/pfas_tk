# Protein–water partition coefficients for per- and polyfluoroalkyl substances: harmonised compilations for serum albumin (288 measurements, 20 sources) and liver fatty acid-binding protein (87 measurements, 9 sources)

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

We then applied the same
treatment to liver fatty acid-binding protein (L-FABP), the standing candidate
explanation for hepatic PFAS accumulation, compiling 87 coefficients for 32
PFAS from nine sources on a per-species bridge constant of 4.8475 for human
L-FABP against albumin's 4.1774. Across the 21 compounds measured against both
proteins the median difference is −0.12 log units, so L-FABP has no intrinsic
affinity advantage over the plasma carrier per unit protein mass — consistent
with the one in vivo test available, in which deleting L-FABP changes the mouse
liver partition coefficient by under 20% for 15 of 19 PFAS. Part of the
remaining L-FABP dispersion is structural: the displacement assays use probes
that occupy different sites of the two-site cavity, and sources using a
secondary-site probe and a primary-site probe are not estimating the same
quantity.

**Keywords:** PFAS; serum albumin; liver fatty acid-binding protein;
partition coefficient; protein binding; toxicokinetics; data harmonisation

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

Plasma protein is only half of the problem. Several PFAS accumulate
preferentially in liver, and the standing explanation for that is binding to
liver fatty acid-binding protein (L-FABP), the most abundant cytosolic protein
of the hepatocyte. The L-FABP literature has the same disorder as the albumin
literature and a quarter as many measurements, and it carries one additional
difficulty: the protein has two binding sites, and the displacement assays that
supply most of its constants use probes that occupy them differently. Section 5
applies the treatment of sections 2 to 4 to that protein, on the same basis, so
that the two can be compared directly — which turns out to matter, because the
comparison is the test of whether L-FABP binding can explain hepatic
accumulation at all.

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

## 5. Liver fatty acid-binding protein

### 5.1 A second protein, and a different bridge

Liver fatty acid-binding protein (L-FABP, *FABP1*) is the most abundant
cytosolic protein of the hepatocyte, 2–5% of cytosolic protein and reaching
0.1–0.4 mM. It binds long-chain fatty acids and their oxidation products, and
because several PFAS are structural analogues of those ligands it has been the
standing candidate explanation for the preferential accumulation of long-chain
PFAS in liver. The same harmonisation problem applies to it as to albumin, and
with fewer measurements to work from.

The conversion identity of §2.2 is unchanged, but its leading constant is not.
That constant is log₁₀(10⁹/*M*) for a protein of molar mass *M*, and L-FABP is
about one fifth the mass of albumin, so it must be recomputed per species:

| Protein | *M* (Da) | UniProt | Constant |
|---|---|---|---|
| human L-FABP | 14,208 | P07148 | 4.8475 |
| rat L-FABP | 14,273 | P02692 | 4.8455 |
| mouse L-FABP | 14,246 | P12710 | 4.8463 |
| human serum albumin, for contrast | 66,463 | P02768 | 4.1774 |

Passing an L-FABP dissociation constant through albumin's 4.1774 understates
log *K* by 0.67 log units, a factor of 4.7. The three L-FABP constants differ
from one another by 0.002 log and are interchangeable in practice; the
distinction that matters is the one from albumin.

### 5.2 Composition

The compilation holds **87 partition coefficients for 32 PFAS from nine primary
experimental sources**, spanning log *K* 1.35 to 5.85 with a median of 3.53 and
an interquartile range of 2.70–4.03. Sixty-eight values are human, from seven
sources; nineteen are rat, from two.

**Table 6.** Composition by source.

| Source | *n* | Species | Compounds | Median log *K* | Range | Method |
|---|---|---|---|---|---|---|
| Yang 2020 | 17 | human | 16 | 2.75 | 1.89–4.21 | probe displacement |
| Gao 2019 | 14 | human | 14 | 2.54 | 1.90–3.12 | equilibrium dialysis |
| Fan 2025 | 13 | rat | 13 | 3.72 | 2.72–4.12 | probe displacement |
| Zhang 2013 | 12 | human | 12 | 3.11 | 1.83–3.82 | probe displacement |
| Jia 2022 | 11 | human | 11 | 3.63 | 1.35–4.76 | microscale thermophoresis |
| Woodcroft 2010 | 6 | rat | 5 | 3.39 | 2.06–4.35 | probe displacement / ITC |
| Sheng 2018 | 6 | human | 6 | 4.18 | 3.66–4.32 | probe displacement |
| Sheng 2016 | 5 | human | 3 | 4.43 | 3.32–4.73 | probe displacement / ITC |
| Khazaee 2021 | 3 | human | 3 | 5.59 | 4.62–5.85 | equilibrium dialysis |

**Table 7.** Basis of the stoichiometry used in conversion.

| Route | *n* |
|---|---|
| *n* = 1 assumed (paper reports no stoichiometry) | 42 |
| *n* = 1 by the paper's own single-site model | 39 |
| ITC, one set of sites | 3 |
| Displacement, high-affinity site of a two-site fit | 2 |
| ITC, sequential, high-affinity site of a three-site fit | 1 |

Only four values rest on a stoichiometry the authors counted; the remaining 83
rest on *n* = 1, either assumed by us or fixed by the paper's own model. The
L-FABP literature is thus more dependent on the single-site assumption than the
albumin literature is, where 11 rows carried a fitted *n*.

Assay conditions are better reported here than for albumin: a protein
concentration is recoverable for all 87 rows, and a ligand concentration and
ligand-to-protein ratio for 64. Those ratios span 2 to 4,000 with a median of
300 — that is, most of this literature is measured at a large molar excess of
ligand over protein, the opposite of the physiological regime.

### 5.3 Per unit protein mass, L-FABP and albumin are not distinguishable

Twenty-one compounds have been measured against both proteins. The median
difference between them, taken per compound on the common log *K* basis, is
**−0.12 log units** — L-FABP marginally the weaker, by an amount far inside the
dispersion of either compilation.

This is the most consequential number in this section and it is easy to
misread, so it is worth stating what it does and does not say. It does not say
the two proteins are interchangeable in the body: albumin is present in plasma
at roughly 600 µM and L-FABP in hepatocyte cytosol at 100–400 µM, in different
compartments, so equal affinity per kilogram of protein still produces very
different capacities in context. What it does say is that there is no intrinsic
affinity advantage of L-FABP for PFAS over albumin. Any explanation of
preferential hepatic accumulation that rests on L-FABP binding PFAS
*more tightly* than the plasma carrier does not survive this comparison, and
§5.8 shows the in vivo evidence pointing the same way.

### 5.4 Dispersion, and a structural cause for part of it

Thirteen compounds carry three or more values. The median within-compound range
is 1.39 log units and the maximum is 3.30, for perfluorooctanoic acid across
nine sources (log *K* 2.55 to 5.85). The albumin compilation's equivalent
figures were 2.38 and 4.22 across 32 compounds, so L-FABP is better behaved,
but not by much and on a quarter as many compounds.

**Table 8.** Distribution by method family.

| Method | *n* | Median log *K* | Range |
|---|---|---|---|
| Fluorescent-probe displacement | 55 | 3.58 | 1.83–4.73 |
| Equilibrium dialysis | 17 | 2.84 | 1.90–5.85 |
| Microscale thermophoresis | 11 | 3.63 | 1.35–4.76 |
| Isothermal titration calorimetry | 4 | 4.19 | 3.32–4.43 |

Method matters, as it did for albumin, but here a second and more specific
axis is visible. The displacement assays do not all use the same probe, and
the probes do not all occupy the same site.

**Table 9.** Displacement rows by fluorescent probe.

| Probe | *n* | Median log *K* | Sources |
|---|---|---|---|
| 1,8-ANS | 38 | 3.69 | Sheng 2016, Sheng 2018, Yang 2020, Fan 2025 |
| NBD-lauric acid | 12 | 3.11 | Zhang 2013 |
| ANS | 5 | 3.06 | Woodcroft 2010 |

L-FABP has two binding sites, a completely internalised primary site and a
secondary site at the entrance to the cavity, and the probes report on them
differently. Luebker et al. (2002) state that their probe "is believed to bind
exclusively to the secondary binding site on L-FABP", and conclude that to
reach the primary site "a probe known to bind both sites, i.e. radiolabeled
oleic acid, would need to be used". Woodcroft et al. (2010) cite that
limitation explicitly as their reason for changing probe, and fit ANS to two
sites with 2:1 stoichiometry by Scatchard analysis. Zhang et al. (2013) use
NBD-lauric acid and state that it binds the primary site, which oleic and
palmitic acid displace.

Part of the spread in Table 9 is therefore structural rather than random: a
constant from a secondary-site probe and a constant from a primary-site probe
are not estimates of the same quantity. This is a different explanation from
the one that accounted for most of the albumin dispersion in §3.3, where the
separation was between an equilibrium measurement and a thermal-denaturation
one.

A floor on how much of the remaining spread is irreducible comes from one
laboratory measuring the same thing twice. Sheng 2016 and Sheng 2018 both
determined PFOA against human L-FABP by 1,8-ANS displacement, two years apart,
and obtained *K*_d of 2.36 and 8.03 µM — a factor of 3.4, tracking a factor of
2.1 in the probe's own constant (3.56 against 7.47 µM). Both are carried.

### 5.5 Chain length

**Table 10.** Perfluoroalkyl carboxylate homologous series, human values.

| Perfluorinated carbons | Compound | *n* | Median log *K* |
|---|---|---|---|
| 3 | PFBA | 1 | 1.90 |
| 4 | PFPeA | 1 | 2.01 |
| 5 | PFHxA | 3 | 2.20 |
| 6 | PFHpA | 2 | 2.31 |
| 7 | PFOA | 8 | 3.77 |
| 8 | PFNA | 6 | 4.12 |
| 9 | PFDA | 2 | 3.35 |
| 10 | PFUnDA | 2 | 3.47 |
| 11 | PFDoDA | 3 | 3.76 |
| 12 | PFTrDA | 1 | 2.35 |
| 13 | PFTeDA | 1 | 3.07 |
| 15 | PFHxDA | 1 | 2.79 |
| 17 | PFOcDA | 1 | 3.05 |

Affinity rises by **0.36 log units per perfluorinated carbon from C3 to C9**
and falls thereafter, placing the maximum at perfluorononanoic acid. This is
the same shape, and the same maximum, as the albumin series of §3.4, which rose
at 0.52 log units per carbon. Zhang et al. attribute the decline beyond C11 to
the binding cavity: the fluorocarbon chain is rigid, and beyond a length of
about 13 Å it must bend to be accommodated, which costs more energy than a
flexible hydrocarbon chain would. That explanation is specific to a small,
enclosed cavity and does not transfer to albumin, where the same inflection
occurs.

### 5.6 Species

The rat arm comes from two laboratories (Fan 2025, Woodcroft 2010) and the
human arm from five that measured by probe displacement, so the comparison can
be made with method held fixed. Restricting both sides to displacement:

| Comparison | Shared compounds | Median rat − human |
|---|---|---|
| rat vs Zhang 2013 | 10 | +0.50 |
| rat vs Sheng 2016 | 2 | −0.71 |
| rat vs Sheng 2018 | 3 | −0.62 |
| rat vs Yang 2020 | 6 | −0.10 |

The four human sources disagree among themselves by 1.21 log units — PFOA alone
is 2.36 µM in Sheng 2016 and 50.4 µM in Zhang 2013, a factor of 21 between two
human laboratories — and that spread is wider than any rat–human gap and
contains both signs. **The sign of the apparent species difference is set by the
choice of human comparator.** Pooling across all sources gives +0.40 log over
13 shared compounds, falling to −0.01 on the five compounds carrying four or
more independent human values, which is the signature of thin coverage rather
than biology.

Two published claims sit in opposition on this point and neither survives
cleanly. Modaresi et al. (2025) state that PFAS bind more strongly to human than
to rodent FABP; that holds against exactly one human source, on two compounds,
and reverses against the other three and against the pooled comparison. Fan et
al. (2025) state that their rat affinities are comparable to human values by the
same method, citing Zhang 2013; the agreement is near-exact against Yang 2020,
which they do not cite for it, while against Zhang 2013 the difference is the
largest of the four.

As with albumin, the practical consequence is that a rat-to-human extrapolation
factor for L-FABP binding cannot be estimated from this literature, because the
quantity it would correct for is smaller than the disagreement between
laboratories measuring one species.

### 5.7 Sources read in full and excluded

**Table 11.** Sources and rows examined and excluded, with reasons.

| Source or class | Rows | Reason |
|---|---|---|
| Sheng 2016 engineered variants | 22 | S39G, M74G, N111D and R122G are site-directed point mutants, three of them built to break binding. Pooling them placed PFOA in the sheet eight times across a 120-fold range. Wild type is retained. |
| Han 2021 | 13 | Its L-FABP constants are quoted from other papers, not measured; its own measurement is a thermal-shift Δ*T*_m. |
| Probes and natural ligands | 8 | DAUDA, ANS, 1,8-ANS, NBD-lauric acid, oleic and palmitic acid. Palmitic acid would otherwise have entered at log *K* 5.02, inside a PFAS chain-length relationship. |
| Luebker 2002 | 5 | Table 1 reports IC₅₀ only and fits no dissociation constant; two of four fluorochemicals have no finite value. Converting by Cheng–Prusoff would be the reader's arithmetic for a quantity the authors declined to publish. |
| Lower-affinity sites of multi-site fits | 4 | Retained in the database; only the high-affinity site enters the compilation, as for Moro 2022 on the albumin side. |
| Zhao 2023 | 4 | Docking scores and a QSAR over a 1,599-compound library; no measured constant. |
| Depth-only unique compounds | 2 | 6:2 FTCA and 6:2 FTSA appear only in Sheng 2018 and are held, so that adding a source does not widen the compound list. |
| Salt duplicates | 2 | Identical constants entered under a free acid and its salt. |
| Modaresi 2025 | — | Measures *K*_liver, a tissue partition coefficient by C18 solid-phase microextraction on liver homogenates, not a protein dissociation constant. Discussed in §5.8. |
| Jia 2025 | — | Its 44 constants are plotted and tabulated nowhere, supporting information included; and all 44 PFAS were incubated together as a mixture against 2 µM L-FABP, so each value is a competitive apparent constant under 43 competitors. |
| Zhao 2025 | — | Docking of 430 PFAS against albumin, PPARγ, transthyretin and GPR40, with QSAR and Kronecker-regularised least squares trained on the docking scores. No L-FABP, no measurement. |
| Cao & Ng 2025 | — | Molecular dynamics and docking for 96 photolithography PFAS; predicted affinities only. |
| Mann 2022; Mann & Berger 2023 | — | Engineered hL-FABP biosensors (L50W; acrylodan-conjugated C69S/F50C; split-protein cpGFP fusion). They report limits of detection in ppb, not dissociation constants, and the protein is a designed variant. |

One further source could not be assessed. Yang et al. (2026, *PNAS*) screened
50 human proteins against 481 chemicals by affinity-selection mass spectrometry
and report FABP1 as a selective target of legacy PFAS, with crystallography;
the full text is under a PubMed Central embargo until December 2026 and whether
it reports numerical constants is not determinable from the abstract.

### 5.8 The one in vivo test, and it is negative

Every value above was measured on isolated protein. Modaresi et al. (2025)
instead deleted the protein. Using C18 fibre solid-phase microextraction on
liver homogenates from wild-type and global L-FABP⁻/⁻ mice, with a 24-PFAS
mixture at 400 ng/L, they obtained liver partition coefficients for both
genotypes. **For 15 of the 19 PFAS tested, *K*_liver differed by less than 20%
between wild-type and knockout liver**, with 14 of 19 marginally above the 1:1
line. The four exceeding 20% all carry a sulfonate or sulfonamide head group:
PFDS +32%, PFNS +55%, 10:2 FTS +26%, FOSA −21%. No other FABP was upregulated
in the knockouts, so compensation does not account for the result.

This is not a contradiction of the binding data; it is what the binding data
predict once §5.3 is taken seriously. If L-FABP has no affinity advantage over
albumin per unit protein mass, then removing it from a tissue that contains
many other proteins should move the tissue partition coefficient very little,
which is what is observed. The tension that remains for toxicokinetic modelling
is not between the in vitro and in vivo results but in how the in vitro results
have been used: a strong binding constant for an abundant protein has been read
as evidence that the protein determines the accumulation, and that inference
does not follow from affinity alone.

---

## 6. Limitations

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

**The L-FABP compilation is thinner in every dimension.** Eighty-seven rows
against 288, 32 compounds against 63, and nine sources against 20. Nineteen of
its compounds rest on a single source, and only four of its 87 values carry a
stoichiometry the authors counted. Where the albumin dispersion statistics rest
on 32 compounds with three or more values, the L-FABP equivalents rest on 13.

**The L-FABP rat arm is two laboratories and the human arm is unbalanced.**
Nineteen rat values come from two sources; six of the seven human sources
contribute fewer than 18 rows each. The method-matched species comparison of
§5.6 rests on as few as two shared compounds for one of its four contrasts.

**Most L-FABP measurements are made far from physiological loading.** Where a
ligand-to-protein ratio is recoverable it has a median of 300 and reaches
4,000. Albumin's compilation contains measurements at ratios below 10⁻³; the
L-FABP literature contains none, so any non-linearity between the two regimes
is unobserved here rather than excluded.

**One L-FABP source could not be assessed.** Yang et al. (2026) is under a
PubMed Central embargo until December 2026 (§5.7).

**Enthalpy is not screened systematically.** The failure mode identified in
§4.2 — a binding enthalpy outside the physically available range — was detected
only because one source reported three temperatures. A systematic enthalpy
screen across all sources reporting one would be the cheapest available
improvement to data quality, and has not been performed.

---

## 7. Conclusions

Two harmonised compilations are presented on one basis: 288 albumin–water
partition coefficients for 63 PFAS, and 87 L-FABP–water coefficients for 32,
each converted to a mass basis through an explicit and mechanically verified
identity. The compilation's central finding is that the
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

For L-FABP the same treatment yields 87 coefficients for 32 PFAS from nine
sources, and one finding that bears directly on how the protein is used in
toxicokinetic models: across the 21 compounds measured against both proteins,
L-FABP and albumin are indistinguishable per unit protein mass, differing by a
median of −0.12 log units. The chain-length relationship has the same shape and
the same C9 maximum as albumin's, at +0.36 log units per perfluorinated carbon.
Part of the L-FABP dispersion has a structural cause absent from the albumin
literature — the displacement probes occupy different sites of a two-site
cavity, and the authors of the earliest of those studies say so explicitly. The
case for L-FABP as the determinant of hepatic PFAS accumulation does not rest on
affinity, since it has no affinity advantage; and the only study to delete the
protein in vivo finds liver partitioning essentially unchanged for 15 of 19
PFAS.

The exclusions documented in §4.2 and §5.7 are as much a result as the
compilations themselves. Of seven sources read in full and rejected, two were rejected on grounds
stated plainly in the source's own text, two on arithmetic performed on the
source's own tables, and one on the source's own internal contradiction. In none
of these cases was the defect apparent from the value as cited elsewhere.

---

## 8. Data and code availability

| Artefact | Path |
|---|---|
| Delivered albumin compilation (288 rows) | `albumin_logK_Kd_FINAL.xlsx` |
| Delivered L-FABP compilation (87 rows) | `lfabp_logK_Kd.xlsx`, `lfabp_rows_for_excel.tsv` |
| Source database (1,460 rows, all proteins and species) | `db/protein_binding.csv` |
| Conversion pipeline | `scripts/build_binding_affinity_table.py`, `scripts/albumin_kd_to_k.py`, `scripts/albumin_shareable.py`, `scripts/build_revised_working_table.py` |
| Eligibility screen | `scripts/screen_fits_database.py` |
| Conversion audit | `scripts/audit_working_table.py` |
| Exclusion audit | `scripts/excluded_sources.py` |
| L-FABP pipeline | `scripts/build_lfabp_working_table.py`, `scripts/lfabp_rows_to_tsv.py` |
| L-FABP audits | `scripts/audit_lfabp_species.py`, `scripts/audit_lfabp_anchors.py` |
| L-FABP report statistics | `scripts/lfabp_report_stats.py` |
| Per-source audits | `scripts/audit_*.py` |
| Extended methodological primer | `report/ALBUMIN_BINDING_PRIMER.md` |

The primer referenced in the final row is a companion document covering the
derivations, worked examples, and per-source reasoning at greater length than is
appropriate here, including the full text of each conversion and admission
decision.

---

## 9. References

Primary sources are cited by first author and year as they appear in Table 1.
Full bibliographic detail for each, together with the specific table or figure
from which each value was taken, is recorded in the `source_table` field of
`db/protein_binding.csv` on a per-row basis, and in the per-source audit scripts
listed in §7.
