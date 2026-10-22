# Mass-balance dosimetry for in vitro PFAS testing: a source-verified reassessment of two models against three uptake datasets

**Working manuscript**, 22 October 2026. Analysis, figures and tables are
reproducible from `tk_review/scripts/`; see *Data and code availability* for
the script behind each section and figure. Four figures and nine tables.
All model parameters and all exposure conditions are verified against primary
sources; the provenance of every reference is stated after the reference
list.

---

## Abstract

**Background.** Quantitative in vitro to in vivo extrapolation (QIVIVE) for
per- and polyfluoroalkyl substances (PFAS) requires knowing how much of an
administered dose actually reaches cells in a bioassay. Nominal medium
concentrations overstate it, sometimes by three orders of magnitude, because
PFAS bind extensively to serum albumin and sorb to vessel walls. Equilibrium
mass-balance models (IV-MBM) were developed to recover the freely dissolved and
cellular concentrations from first principles.

**Objective.** To reassess two published IV-MBM formulations — the
single-descriptor IV-MBM EQP v2.0 of Armitage et al. and the
concentration-dependent balance of Qin et al. — against the three published
datasets reporting PFAS cellular uptake in human liver cell lines, with all
model parameters and all exposure conditions verified against primary sources.

**Approach.** Both models were re-implemented from their published equations and
supplementary parameterisations. The three uptake datasets (Rosenmai et al.;
Amstutz et al.; Fragki et al.) were harmonised onto the fraction-of-dose axis.
Exposure conditions were read from the primary methods sections rather than
inherited from derived descriptions. Model residuals were decomposed by study,
compound, dose and chain length, and each candidate parameter was replaced with
the best available independent measurement.

**Principal findings.** (i) Reading the primary methods revealed that one of the
three assays was dosed in 2% rather than 10% serum, a five-fold error in the
dominant sorptive phase that is not apparent from secondary descriptions;
correcting it degraded the apparent fit of the dataset that had previously
looked best. (ii) On corrected media both models over-predict cellular uptake —
geometric-mean bias 3.9× (Armitage) and 7.5× (Qin) over 29 observations. (iii)
Every substitution of a best-available measured parameter made the fit worse,
in each case because the measurement moved the balance in the direction the
residual did not need. (iv) The residual separates into a level error and a
chain-length error with different causes, and the separation is clean: the
chain-length behaviour is a property of the partition coefficients and the
level is a property of the volumes and compositions, each invariant to changes
in the other. Reconstructing a third published formulation — the
measured-parameter linear balance of Fischer et al., using measured PFAS
albumin and membrane coefficients — removes the chain-length drift entirely
over nCF 5–7 while leaving a 5.6× over-prediction, so no parameterisation
tested here corrects the level. The over-steepening of the lipid term that
causes the drift is corroborated by two independent laboratories, two lipid
preparations and two analytical routes (+0.36 to +0.48 decades per CF₂
measured, against +0.78 predicted). (v) The same formalism applied to human serum reproduces measured
whole-serum partitioning for 26 PFAS with no fitted parameter, establishing
that the accounting is sound and localising the failure to the in vitro system.
(vi) All three assays operate 5,500–27,700× above the physiological
ligand:protein ratio, confirming and quantifying the regime concern raised by
Qin et al.

**Implication.** Improved bookkeeping within the well yields a better
description of the experiment but does not by itself yield a better oral
equivalent dose, because the experiments occupy a binding regime that does not
occur in vivo. We give specific reporting requirements that would make such
assays usable for QIVIVE.

**Keywords:** PFAS; in vitro dosimetry; mass balance; QIVIVE; serum albumin;
ligand:protein ratio; freely dissolved concentration

---

## 1. Introduction

### 1.1 The nominal concentration problem

In vitro bioactivity data are reported almost universally on a nominal basis:
the concentration weighed into the medium. For chemicals that distribute
appreciably out of the aqueous phase, that number is not the exposure the cells
experienced. Chemical partitions into serum lipids and proteins, sorbs to the
plastic of the well, is taken up by the cells themselves, and — for volatiles —
escapes to the headspace. The freely dissolved concentration that drives uptake
can be orders of magnitude below nominal [1,12].

This matters because the purpose of the measurement is extrapolation. If an
effect concentration measured in a dish is to be compared with a blood
concentration in a person, the two must be expressed on comparable bases. A
bioassay medium may contain anywhere from 0 to 20% serum; blood is blood. Using
a nominal effect concentration as a direct surrogate for a blood concentration
embeds a bias whose size depends on hydrophobicity and on assay composition
[1].

Equilibrium mass-balance models were developed to close this gap
[1,2,12,13,16]. They are now sufficiently mature that several independent
implementations exist, and a comparison of their features has been published
[17, cited via 1].

### 1.2 Why PFAS are a difficult case

The spLFER-based mass-balance approach was calibrated on neutral, hydrophobic
organics, for which a single descriptor — the octanol–water partition ratio —
predicts partitioning to most biological phases well. PFAS violate the premises
of that approach in three distinct ways.

First, they are surfactants: a perfluorinated tail with an ionisable head
group. Their behaviour at interfaces and in mixed phases is not well captured
by a bulk hydrophobicity descriptor, and the descriptor itself is difficult to
measure for them.

Second, they are essentially fully anionic at physiological pH. Ionisation is
handled in these models by scaling the neutral-form partition ratio by a
phase-specific penalty. For serum albumin, Armitage et al. set that penalty to
zero for acids — the anion is assumed to bind albumin exactly as well as the
neutral form would (ref. [1], SI Table S1). This is a modelling decision rather
than a measurement, and it is consequential: it is a major reason albumin ends
up carrying the great majority of the chemical in these calculations.

Third, PFAS bind albumin at discrete sites rather than dissolving into it
diffusely. Site binding saturates. Qin et al. [3] measured PFAS–albumin
isotherms over five orders of magnitude and found the binding strongly
concentration-dependent, which is the basis of their model and the reason the
linearity assumption common to spLFER approaches cannot simply be carried over.

### 1.3 Aims

This reassessment addresses four questions:

1. How well do the two models predict PFAS cellular uptake when both are
   implemented from their own published parameterisations?
2. Are the available uptake data mutually consistent enough to adjudicate
   between models?
3. Where does the residual lie, and can it be attributed to a parameter?
4. What do the answers imply for the use of these assays in QIVIVE?

We stress at the outset that question 4 was raised, and answered in principle,
by Qin et al. [3], who hypothesised that concentration-dependent protein
binding causes the binding regime in bioassays to differ substantially from
that in plasma, with consequences for QIVIVE. What follows quantifies that gap
for the specific uptake datasets used to parameterise and test these models,
and examines what it leaves achievable.

---

## 2. The mass-balance formalism

### 2.1 Governing relation

The well is treated as a set of phases at equilibrium. Conservation of the
administered amount, with each phase concentration expressed relative to the
freely dissolved concentration by a distribution ratio *D*, gives

> *n*<sub>total</sub> = *C*<sub>free</sub> · ( *V*<sub>w</sub> + Σ<sub>i</sub> *D*<sub>i</sub>*V*<sub>i</sub> )

so that *C*<sub>free</sub> = *n*<sub>total</sub> / *S*, where *S* ≡
*V*<sub>w</sub> + Σ *D*<sub>i</sub>*V*<sub>i</sub> is the total sorptive capacity
of the well in units of water-equivalent volume. The fraction of the
administered dose residing in cells follows as *x*<sub>cell</sub> =
*D*<sub>cell</sub>*V*<sub>cell</sub> / *S*.

Models differ only in how the *D*<sub>i</sub> are obtained. The accounting is
not in dispute.

### 2.2 Assumptions

The derivation requires: equilibrium; no degradation, metabolism or active
transport; no losses; additive, non-interacting phases; and — the load-bearing
one — linear sorption, i.e. constant activity coefficients, which holds only
while the chemical is dilute with respect to every sorptive phase.

For PFAS in these assays, equilibrium is defensible: Fragki et al. [7] report
cellular steady state within approximately 1 h against a 24 h exposure.
Metabolic stability is defensible for perfluoroalkyl acids, though not for
fluorotelomer alcohols. Active transport is assumed absent and cannot be
checked. The linearity assumption is the one that fails outright, and §7
quantifies by how much.

### 2.3 Choice of output metric

Two metrics are available: the cellular concentration and the fraction of the
dose in cells. The published uptake data are reported as amount per mass of
cell protein, so converting them to a concentration requires a cell volume,
while converting them to a fraction requires only protein per cell.

It should be noted — because it is easy to get wrong — that this does *not*
make the fraction axis immune to cell-volume uncertainty. The predicted
fraction is proportional to cell volume even though the measured fraction is
not, with the result that the ratio of predicted to measured scales with
assumed cell volume on *either* axis. What the fraction axis does buy is that
the three datasets' *measurements* become mutually comparable without any cell
volume being assumed. We use it for that reason.

---

## 3. The two models

### 3.1 IV-MBM EQP v2.0 (Armitage et al.)

The model predicts all required distribution ratios from the neutral-form
octanol–water partition ratio via sorbent-phase-specific linear free energy
relationships (spLFERs), with ionisation handled by phase-specific scaling
factors [1]. The cell is represented as a proximate composition — storage
lipid, membrane lipid, bulk protein and water — each contributing in proportion
to its volume fraction. Plastic sorption is included via one of two published
spLFERs.

We verified every parameter used here against the published supplementary
material: both spLFERs (SI §S1.3), all four ionisation scaling factors (SI
Table S1), and the cell-line compositions (SI Table S5). All agree exactly.

One parameter could not be verified. How bulk cell protein is converted into
sorptive capacity is stated nowhere in the main text or the supplementary
material; §S1 covers anionic phospholipids, lysosomal sequestration and the
spLFERs, then defers to ref. [2] for temperature and ionic-strength adjustments
only. The earlier paper does not contain it either, and the reason is stronger
than the absence of a documented value: in v1.0 there is no protein term at
all. Its mass balance represents the cell as a single lumped phase,
K<sub>CW</sub> = K<sub>OW</sub>·*f*<sub>L</sub>, where *f*<sub>L</sub> is "the
volume fraction of lipid equivalents (i.e., pseudo-octanol content)" after the
standard multimedia-model convention, assumed to be 0.05 for cells and varied
over 0.01–0.05 [2]. There is no membrane-lipid phase and no cell-water phase
either. The proximate composition of v2.0 is a new construction, so the
coefficient attached to its protein phase cannot be recovered from the earlier
record; and the 3.5% convention does not appear there. The implementation used
here therefore adopts the convention in EPA's `httk` package, under which 3.5%
of bulk protein is treated as octanol-equivalent.

Two further points follow from the earlier paper and bear on how much this
matters.

First, **v1.0 excludes PFAS by its own statement**: it "is intended for neutral
organic chemicals and should only be applied to ionogenic organic chemicals
(IOCs) if the compound is predominantly neutral at the test pH", and its
future-work list names the gap — "to expand the approach for neutral organic
chemicals to ionogenic organic chemicals that are predominantly charged in the
test system" [2]. The IOC machinery of v2.0 is that extension. Fully anionic
PFAS are among the first applications of it, not a case on which the formalism
was calibrated.

Second, the undocumented coefficient has **little leverage at its own value**.
With protein contributing 0.05–0.25% of the cell's sorptive capacity, the
modelled cell is 79% water for the shortest chain tested and 99.7% membrane
lipid for the longest: the cell is, in effect, a lipid droplet. The model's
output is therefore insensitive to the value of 0.035 over any plausible range.
What has leverage is the *choice* it encodes — octanol-equivalence for bulk
protein — and replacing that choice with a measured structural-protein affinity
enlarges the cell sink by 1.4–3.4× (§6.4). The distinction matters for reading
what follows: results that turn on the modelling choice are robust; results
that would turn on the particular constant do not arise, because none is
sensitive to it.

### 3.2 The concentration-dependent balance (Qin et al.)

Qin et al. [3] measured PFAS binding to bovine serum albumin, lipid, human
plasma, assay medium and cells by C18 solid-phase microextraction over five
orders of magnitude, and represent albumin binding as a saturable term plus a
linear non-specific term. Their biomaterial–water ratio is

> D<sub>i/w</sub> = D<sub>BSA/w</sub>·(V<sub>prot</sub>/V<sub>prot+lip</sub>)
> + D<sub>lip/w</sub>·(V<sub>lip</sub>/V<sub>prot+lip</sub>)

and the bioassay balance is

> C<sub>free</sub> = n<sub>tot</sub> /
> [(D<sub>cell/w</sub>·V<sub>cell,prot+lip</sub> + V<sub>cell,w</sub>)
> + (D<sub>medium/w</sub>·V<sub>medium,prot+lip</sub> + V<sub>medium,w</sub>)]

Two features of this form are easy to misread and are worth stating explicitly,
because both bear on the comparison in §3.3. First, the normalising volume is
the *sorptive* volume V<sub>prot+lip</sub>, not the volume of the biomaterial;
the authors give the units of their tabulated cell–water ratio as
L<sub>w</sub>/L<sub>cell,prot+lip</sub>, confirming it. Second, cell water does
appear in the balance, but as a non-sorbing phase at D = 1 pooled with the
medium water rather than inside D<sub>cell/w</sub>. The distinction the authors
draw against the earlier Henneberger formulation is about the *numerator*: they
"considered only proteins and lipids of the cell homogenate as the main
sorptive phase", so cell water lowers the free concentration but does not count
as PFAS in the cell. For the cell term itself the composition is protein and
lipid only, 95% protein by sorptive volume for their HEK293H line (8.14%
protein and 0.432% lipid of cell volume; water 91.6%).

The protein term is the measured non-specific protein–water ratio rather than
an octanol surrogate: the authors take D<sub>non-specific</sub> "as a surrogate
for partitioning to structural proteins", on the argument that cytoskeletal
protein binds organic acids more weakly than albumin does.

The saturable albumin term can be derived from mass action: for *n* independent
sites, the site concentration within pure albumin is *n*·ρ/*M*. With the
authors' own stated values (ρ = 1.36 kg/L, M = 66,463 Da) this evaluates to
*n* × 20,462 µM, which reproduces their published site-capacity column to
within 0.2% across all four measured compounds and implies *n* = 0.45 for PFBA
and 1.24 for PFOA. In the trace limit the term collapses to a constant,
recovering the standard bridge between a binding constant and a partition
coefficient.

Parameters were measured for four compounds (PFBA, PFOA, PFHxS, PFOS). Values
for PFHxA and PFHpA are interpolated between the two carboxylate anchors and
PFNA is extrapolated beyond both. Only the lipid term has a published
interpolation rule — a regression of log D<sub>lip/w</sub> on the number of
fluorinated carbons with slope +0.476 and intercept +0.172, which the authors
used to supply the PFBA value and which reproduces their tabulated carboxylate
entries exactly.

### 3.3 Structural comparison

**Table 1.** Structural differences between the two models.

| | IV-MBM EQP v2.0 | Qin et al. |
|---|---|---|
| source of *D* | predicted from log K<sub>OW</sub> | measured, then interpolated |
| albumin | linear, piecewise spLFER | saturable (site binding + non-specific) |
| albumin anion penalty | 0.0 (acids) | not applicable; measured on the anion |
| cell phases | storage lipid, membrane lipid, 3.5%·protein, water | protein, lipid; water non-sorbing |
| cell composition at short chain | ~79% water | ~95% protein (sorptive volume) |
| cell composition at long chain | ~99.7% membrane lipid | ~95% protein (sorptive volume) |
| dose dependence | linear by construction | non-linear |
| chemical coverage | any, from one descriptor | four measured compounds |

The two models hold opposite views of what a cell is. Their predicted uptake
nonetheless agrees to within a factor of ~2.5 in many cases, because albumin
dominates the denominator in both and masks the difference. Agreement between
them should not be read as corroboration of either cell representation.

### 3.4 The wider landscape, and why two models

Proença et al. [17] catalogue fourteen published in vitro distribution models
and classify them by the strategy used to obtain partitioning data:
experimental values for specific chemicals, experimental values for a specific
system, or QSPRs. Ten of the fourteen are specific to particular chemicals or
particular setups; four are generic. The models relevant here are the generic
ones, and they differ along axes that matter for PFAS:

**Table 1b.** Design axes across published in vitro distribution approaches.

| | source of partitioning | albumin treatment | cell representation | normalising volume |
|---|---|---|---|---|
| Armitage et al. v1.0 [2] | K<sub>OW</sub> × lipid-equivalent fraction | linear | one lumped pseudo-octanol phase | V<sub>cell</sub> |
| Armitage et al. v2.0 [1] | spLFERs from log K<sub>OW</sub> | linear, piecewise | lipid + 3.5%·protein + water | V<sub>cell</sub> |
| Fischer et al. [12] | **measured**, else PP-LFER / 3D-QSAR / COSMOmic | linear | protein + lipid + water, measured; D<sub>BSA/w</sub> for all protein | V<sub>cell</sub> |
| Qin et al. [3] | measured, PFAS-specific | **saturable** | protein + lipid; water non-sorbing | V<sub>prot+lip</sub> |
| Stadnicka-Michalak et al. [18] | regressions on log K<sub>OW</sub> | linear | kinetic, not equilibrium | — |
| VCBA / PBK-coupled [19] | varies | linear | kinetic, coupled to PBK | — |

Two of these — the kinetic approaches — answer a different question, since the
PFAS uptake data here reach steady state within ~1 h against a 24 h exposure
and an equilibrium treatment is adequate. Of the equilibrium models, the two
tested in full occupy opposite corners of the first two axes, which is why they
were selected. The third corner — measured biopartitioning with a *linear*
albumin term — is the Fischer et al. design, and §6.7 reconstructs it from
published PFAS measurements without reimplementing anyone's code.

Proença et al. also ran five of the generic models side by side on six test
chemicals in HepaRG cells at 1 µM, with and without 10% FBS. Two of their
observations bear directly on the results below. The Fischer et al. model
predicts the highest cellular concentration of the five, and the Armitage model
predicts "substantially lower" values, which they attribute to "both the
cellular compartments added (i.e. just lipid or protein and lipid) and QSPRs
used (log K<sub>OW</sub>-based QSPRs or PP-LFERs based on other descriptors)".
That ordering and that attribution are reproduced independently in §6.7.2
below. They also state the limit of their own exercise: "A comparison of model
predictive performance is only possible with adequate experimental data, and is
beyond the scope of this review" — and, in their discussion, that "especially
for the generic models more extensive experimental validation is long overdue".
The present work is an attempt at that validation for one chemical class.

Three further statements in that review anticipate findings reported below, and
are credited here rather than claimed. On albumin: "for anionic chemicals
sorption to albumin cannot be predicted with a simple partition since
saturation of albumin binding sites for these chemicals occurs at
concentrations relevant for toxicity testing", citing Henneberger et al. 2020
[21] — the argument of §7, which predates Qin et al. and which they sharpen by
noting that organic acids "appear to only have one binding site, which results
in the saturation at relatively low concentrations". On lipid: "a bulk solvent
such as octanol is not a good surrogate for anisotropic lipid bilayer
membranes, especially for partitioning of ionised chemicals", and membrane
partitioning "is more similar between neutral and charged species than
K<sub>OW</sub> based predictions will show" — the mechanism behind §6.7's
central result. On the applicability domain: "the main limitation remains
accurate QSPRs to predict partitioning of heavily charged chemicals". PFAS at
pH 7.4 are heavily charged, so the compounds modelled here sit on the boundary
of the domain the review identifies, which is the reason for testing rather
than assuming.

Two details in the last column of Table 1b are easy to lose and both caused
errors in earlier drafts of this work. The three equilibrium balances do not report
cell–water distribution on a common basis: Armitage and Fischer normalise by
total cell volume, Qin by the volume of cell protein plus lipid, and the two
differ by a factor of ~17 for a hepatic composition. And the treatment of cell
protein is not a detail of parameterisation but a structural choice: Armitage
gives it an octanol equivalence, which makes it negligible; Fischer gives it
the albumin surrogate, which makes it comparable to the lipid term; Qin gives
it a measured non-specific protein coefficient, which makes it dominant. A
model's cell is defined as much by that choice as by its composition.

---

## 4. The empirical basis

### 4.1 Three datasets

Three studies report PFAS cellular uptake in human liver cell lines on a
fraction-of-dose basis: Rosenmai et al. [5] (HepG2, 14 PFAS, three doses),
Amstutz et al. [6] (HepG2, 13 PFAS, one dose) and Fragki et al. [7]
(differentiated HepaRG, four PFAS, five doses). Twenty-nine observations
overlap the five compounds for which both models can be parameterised.

One provenance note on the third dataset. Its biokinetic measurements are in
Fragki et al. [7]; the companion paper from the same group on the same cell
system, Louisse et al. [8], reports hepatotoxic potencies and gene-expression
changes but no biokinetics, and is sometimes cited for the uptake data. All
uptake values used here are from ref. [7], and the exposure conditions in
Table 2 are from its methods section, not from ref. [8].

### 4.2 Exposure conditions, read from source

**Table 2.** Exposure conditions as stated in the primary methods sections.
Entries in bold differ from what a reader would infer from the cell-culture
section alone.

| | Rosenmai [5] | Amstutz [6] | Fragki [7] |
|---|---|---|---|
| cell line | HepG2 | HepG2 | differentiated HepaRG |
| growth medium serum | 10% FBS | 10% FBS | 10% FBS |
| **exposure medium serum** | 10% FBS | 10% FBS | **2% FBS** |
| co-solvent | **3% ethanol** | 0.1% DMSO | 0.5% DMSO |
| plate format | 12-well, 3 wells pooled | 96-well | 24-well |
| volume dosed → after t=0 sample | not stated | 250 → 200 µL | 550 → 500 µL |
| wash | 5× ice-cold PBS | 5× DPBS | 5× DPBS |
| uptake normalised to | nominal | **recovered** | nominal |

Two entries require comment. Fragki et al. state that exposure was conducted in
"assay medium (growth medium containing 2% FBS)"; the 10% figure belongs to the
growth and differentiation medium. Because serum fraction multiplies the volume
of the dominant sorptive phase, this is a five-fold difference in the parameter
to which the balance is most sensitive, and it is not recoverable from
descriptions that quote the culture conditions. Second, Amstutz et al.
normalise their intracellular fraction to *recovered* rather than
administered mass, so their denominator excludes material lost to the vessel;
their values are not strictly on the same axis as the other two.

Rosenmai et al. do not state a working volume; 1 mL per well is assumed, as in
all prior analyses of this dataset known to us.

### 4.3 Inter-study consistency

The three datasets disagree with one another before any model is involved. At
100 µM, cellular uptake of PFNA is reported as 0.31% [5], 0.30% [6] and 4.39%
[7] — a 14-fold spread. For PFHxA at the same dose the two HepG2 studies differ
by 3.6-fold (0.14% [5] vs 0.039% [6]) despite using the same cell line and
nominally similar protocols; part of this is attributable to the different
normalisation noted above.

Any model-performance statement must be read against this. A model disagreeing
with these data by a factor of four is disagreeing with a target that is itself
uncertain by a comparable factor.

---

## 5. Methods of this reassessment

Both models were implemented in Python from their published equations and
parameterisations (`mbm_published_equations.py`). Armitage parameters were
taken from ref. [1] and its supplementary material; Qin parameters from ref.
[3]. Where the Armitage implementation required a parameter not present in the
published record (cell protein sorption, §3.1), the `httk` convention was used
and is flagged throughout.

Assay geometry, cell number, serum fraction and dose series were taken from the
primary publications (Table 2). Measured uptake was placed on the
fraction-of-nominal-dose axis. Model performance is reported as geometric-mean
bias (predicted ÷ measured) and geometric-mean absolute error.

Residuals were decomposed by study, by compound, by dose and against
perfluorocarbon chain length. Four parameters were then replaced in turn with
the best available independent measurement: albumin partitioning from Fischer
et al. [9]; cell partitioning from the measured values in ref. [3]; cell protein
sorption from ref. [4]; and cell volume and composition from Fischer et al.'s
measurements on HepG2 [12]. A third formulation, the measured-parameter linear
balance of Fischer et al. [12], was then reconstructed from its published
equations and run on the same data with measured PFAS coefficients
throughout.

---

## 6. Results

### 6.1 Effect of the medium correction

Implementing Fragki et al. at the serum fraction stated in their methods rather
than the growth-medium value changes that dataset's geometric-mean bias from
1.54× to 14.92× (Qin) and from 2.48× to 10.80× (Armitage). Pooled across all 29
observations, bias moves from 3.43× to 7.52× (Qin) and 2.33× to 3.87×
(Armitage).

The direction is worth stating plainly: the dataset that had appeared to be the
models' best case was the one being supplied with five times more albumin than
it contained. Correcting an error degraded the apparent fit. This is the
expected behaviour when a model is accidentally compensated, and it is a
cautionary result for any benchmarking exercise that inherits assay conditions
from secondary sources.

The Armitage response (4.4×) matches the denominator arithmetic, since its
albumin term is linear. The Qin response (9.7×) exceeds it, because a saturable
sink amplifies an error in its own size: reducing albumin raises the free
concentration, which raises site occupancy, which lowers the apparent
distribution ratio, which reduces capacity further. The amplification is
largest at intermediate doses, where the 10% and 2% wells fall on opposite
sides of saturation, and vanishes at the top dose where both are saturated.

### 6.2 Overall performance

**Table 3.** Model performance on corrected media, 29 observations.

| | geometric-mean bias | geometric-mean \|error\| |
|---|---|---|
| IV-MBM EQP v2.0 | **3.87×** | **4.02×** |
| Qin et al. | 7.52× | 7.52× |

| by study | Armitage | Qin |
|---|---|---|
| Rosenmai [5] (n=15) | 1.97× | 3.72× |
| Amstutz [6] (n=4) | 3.77× | 18.95× |
| Fragki [7] (n=10) | 10.80× | 14.92× |

Both models over-predict throughout. The error is not uniform across
compounds: for PFBA–PFHpA, which only the two 10% FBS studies measure, Armitage
bias is 1.2–1.6×, while PFOA and PFNA — the compounds Fragki contributes —
carry almost the entire residual.

![**A**, the three uptake datasets on a common axis: measured fraction of dose
in cells against nominal concentration, for the five compounds both models can
be parameterised for. The datasets disagree with each other by up to 14× at
PFNA, so no single model can match all three. **B**, predicted against measured
for the 29 observations; open symbols are Fragki et al. as coded at 10% FBS and
filled symbols after the correction to the 2% at which the assay was actually
dosed (§6.1), with arrows joining the pairs. The grey diagonals are the 1:1
line and 10× over. **C**, where the balance puts the dose in each well, by
sorptive phase. Albumin carries 94–97% of the capacity at PFOA and PFNA in the
10% FBS assays; the cells carry a visible share only in Fragki's 2% FBS medium
at the long-chain end. **D**, the geometric-mean bias of each parameter
substitution, with 1× a perfect model (§6.1, §6.4). Every best-available
measurement moves at least one model further from the data. The last two bars
are the same substitution taken two ways: Fischer et al.'s HepG2 cell volume
alone, and that volume together with the cell composition they measured in the
same table.](../figures/fig20_mbm_three_datasets.png)

### 6.3 Dose dependence

Armitage is linear by construction and its fitted dose exponent is pinned at
1.00. Measured exponents run 0.85–1.42 and rise with chain length. Qin
reproduces that ordering, overshooting by 0.30–0.71 on the Rosenmai series and
by 0.07–0.20 on the Fragki series. This is the one respect in which the
additional machinery of the saturable model buys something the linear model
structurally cannot deliver, and it is also the one respect in which the medium
correction *improves* Qin: the Fragki exponents fall from 1.71 and 1.66 to 1.49
and 1.62 against a measured 1.42, indicating that part of the curvature
previously attributed to chemistry was an artefact of an oversized albumin
sink.

### 6.4 Parameter substitutions

**Table 4.** Effect of replacing each parameter with the best available
independent measurement. Bias is the geometric mean over 29 observations.

| substitution | direction of the measurement | Armitage | Qin |
|---|---|---|---|
| published models | — | 3.87× | 7.52× |
| albumin partitioning [9] | **below** both models | 8.08× | 7.75× |
| cell partitioning (measured, [3]) | mixed | 2.15× | 4.61× |
| cell protein sorption [4] | **above** both models | 9.81× | 12.24× |
| cell volume alone (measured HepG2 [12]) | — | 1.84× | 3.71× |
| cell volume **and** composition [12] | — | 9.10× | — |

Each measurement moves the balance in the direction the residual does not
require. A weaker albumin reduces competition and sends more chemical into the
cell; a stronger cell-protein affinity enlarges the cell sink and does the
same. Both models already over-predict.

The cell-partitioning substitution improves the aggregate but at the cost of
destroying the chain-length trend, over-predicting PFBA by 19.7× (Qin) and
23.6× (Armitage) while landing within a factor of 1.4 at PFNA. Aggregate
improvement achieved by averaging a large miss at one end of a series against a
good fit at the other is not evidence of a better parameter.

The cell-volume row requires care, and an earlier version of this work read it
wrongly. Taken alone it looks like the one substitution that helps: 3.87× to
1.84×, and it is also the only parameter that rescales every compound by the
same factor, which is what a chain-length-independent offset requires. But the
measured cell volume comes from Fischer et al. [12], who measured the
*composition* of the same HepG2 cells in the same table. Using the volume while
retaining Armitage's composition takes half of one measurement. Taking it whole
moves the bias to 9.10× — the measured hepatic lipid fraction is 6.9× larger
than the published parameterisation's while the cell volume is only 3.3×
smaller, so the cell sink roughly doubles.

**Every best-available measurement therefore worsens the fit, with no
exception. This is the central negative result of the reassessment: a parameter
can be the most influential term in a model and not be the term that is
wrong.** Albumin carries 94–97% of the sorptive capacity for the longer chains,
which makes it the obvious suspect and the obvious thing to tune; the best
available measurement of it moves the prediction further from the data. The
corollary is methodological: a substitution that improves a fit should be
checked for whether it imported only the convenient half of its source.

### 6.5 Chain-length structure of the residual

Armitage's over-prediction rises monotonically with perfluorocarbon chain
length; Qin's does not. Within the Rosenmai dataset alone — one laboratory, one
cell line, one medium, one vehicle — the Armitage residual climbs 0.86× to
9.33× between PFHxA and PFNA. Qin shows no trend within any single study.

The cause is not either distribution ratio alone. Where albumin dominates the
denominator, the predicted fraction in cells varies as the *ratio*
*D*<sub>cell</sub>/*D*<sub>albumin</sub>, so the residual's chain-length slope
is the error in the slope of that ratio.

**Table 5.** Chain-length dependence over nCF 5–8, decades per CF₂.

| | albumin | cell | **ratio** | residual, observed |
|---|---|---|---|---|
| IV-MBM EQP v2.0 | +0.386 | +0.769 | **+0.383** | **+0.345** |
| Qin et al. | +0.547 | +0.515 | **−0.032** | **−0.074** |
| measured albumin [9] | **+0.540** | — | — | — |

The ratio slope predicts the residual slope to within the uncertainty of a
four-point fit. Armitage's albumin spLFER rises more slowly with chain length
than the measurement, because the piecewise form flattens above log K<sub>OW</sub>
4.5 — a feature calibrated on neutral organics and pharmaceuticals, where a
binding cavity of finite depth does plateau. Its cell term rises faster than a
constant ratio would require, because a cell represented as lipid gains
affinity with chain length faster than albumin does. The two errors compound,
contributing roughly 40% and 60% of the gap respectively.

Qin's albumin slope sits within 0.007 of the measured value and his ratio is
flat. **The two models therefore fail in different ways: Qin has the
chain-length behaviour approximately right and the level wrong; Armitage has
the level closer and the chain-length behaviour wrong.**

Independent estimates of the cell distribution ratio's chain-length slope span
+0.177 (cell homogenate, [3]) to +0.782 (the Armitage parameterisation), with
the uptake-implied value (+0.407) and whole-tissue unbound fractions [11]
(+0.46 to +0.53) in between. The quantity is bounded, not known.

### 6.6 What does not explain the residual

Two candidate artefacts were tested and rejected.

*Wash losses.* All three studies wash five times with cold buffer, so it cannot
be a between-study variable; and the direction is wrong, since washing should
preferentially strip the shortest, most weakly bound chains, whereas
over-prediction rises with chain length.

*Cytotoxicity at the long-chain end.* Reduced protein content at high dose is
documented in two of the three studies, and Fragki et al. excluded PFNA at
100 µM from their own analysis on those grounds. If cell loss were suppressing
measured uptake, the chain-length trend would be strongest at the highest dose.
It is strongest at the lowest (+0.194 at 10 µM, +0.099 at 100 µM).

*Albumin secreted by the cells.* Two of the three datasets use hepatic lines,
which synthesise albumin, so the albumin sink could in principle exceed what
the serum supplement accounts for. Proença et al. [17] quantify this: albumin
excreted by cells in such assays is "equivalent to a supplementation of
0.0003–0.0009%" of FBS, four to five orders of magnitude below the 2–10% used
here. It cannot contribute.

*Co-solvent.* The three assays differ 30-fold in vehicle strength (Table 2).
A co-solvent lowers every water-to-phase distribution ratio by approximately
the same factor, and the predicted fraction is a ratio of sorptive terms, so
the effect cancels except against the aqueous term. It is therefore material
only where water carries a substantial share of the capacity — PFBA, where it
changes the residual from 1.42× to 1.11×, and negligible elsewhere (PFNA,
9.33× to 9.31×). The 30-fold difference in vehicle is, counter-intuitively, of
almost no consequence for this comparison.

### 6.7 A third parameterisation: measured biopartitioning

The design decision separating the Escher-school balances [12,14,15] from
IV-MBM EQP v2.0 is the use of measured — or multi-descriptor predicted —
membrane–water and protein–water partitioning in place of a single-descriptor
regression, retaining a linear albumin term. Fischer et al. state the premise
directly: distribution ratios derived "by log K<sub>OW</sub> correlation can
lead to a high uncertainty and inaccuracy", whereas methods using several
chemical properties "can substantially improve the quality of predicted
partition constants" [12]. Armitage et al. report the same for their own model
on ionisable organics [1]. Fischer et al. use polyparameter LFERs for neutral
species, a 3D-QSAR for anionic albumin–water and COSMOmic for anionic
lipid–water, preferring measured values where they exist. For PFAS measured
values now exist, so the reconstruction below uses them: albumin from Fischer
et al. [9], phospholipid membrane–water from Chen et al. [10]. PFBA is
excluded (no measured albumin value), so comparisons are over nCF 5–8 and are
directly comparable with the corresponding columns in §6.5.

Both measured series sit below the spLFER predictions — albumin by 2–7× and
membrane lipid by 8–90× over nCF 5–8, with the membrane discrepancy growing
systematically with chain length.

#### 6.7.1 Parameter sensitivity within the Armitage structure

Substituting each measured series into IV-MBM EQP v2.0 while leaving its
structure intact isolates the leverage of each coefficient.

**Table 6.** Coefficient substitutions within the Armitage structure, nCF 5–8.

| parameterisation | bias | slope (decades/CF₂) |
|---|---|---|
| Armitage, as published | 4.35× | +0.316 |
| + measured albumin only | 10.12× | +0.194 |
| + measured membrane only | 0.16× | −0.049 |
| + both measured | 0.40× | −0.156 |
| Qin, as published (saturable albumin) | 8.59× | −0.043 |

Read as a sensitivity analysis, this table says that the measured membrane
coefficient is where the chain-length drift lives: substituting it alone takes
the slope from +0.316 to −0.049 and no other single substitution comes close.

It does **not**, however, describe the Escher-school design, and an earlier
version of this section wrongly presented the last row as such. Within the
Armitage structure, cell protein sorbs by octanol-equivalence, so it
contributes almost nothing (§3.1); substituting a membrane coefficient 8–90×
smaller than the spLFER's therefore collapses the cell sink almost to its water
term, which is what produces the 0.40×. That is an artefact of mixing one
model's structure with another's coefficients, not a result about either.

#### 6.7.2 The Escher-school balance as published

Fischer et al.'s own form differs from the Armitage structure in four ways
[12]. Their cell–water ratio is

> D<sub>cell/w</sub> = VF<sub>prot,cell</sub>·D<sub>BSA/w</sub>
> + VF<sub>lip,cell</sub>·D<sub>lip/w</sub> + VF<sub>w,cell</sub>

so (i) **one protein surrogate is used everywhere** — cell protein receives
D<sub>BSA/w</sub>, not an octanol equivalence, which is the difference that
invalidates the reconstruction above; (ii) cell water is a phase of
D<sub>cell/w</sub>, unlike Qin's formulation; (iii) the normalising volume is
V<sub>cell</sub>, so the basis matches Armitage's and not Qin's; and (iv) there
is no plastic term. The medium protein volume is also different in kind: it is
total measured protein, not albumin. Fischer et al. measured 0.520% at 10%
dialysed FBS, 3.0× the albumin volume the serum convention gives.

Fischer et al. also supply the cell composition these datasets most need —
measured HepG2 (ARE-bla): 87.35% water, 9.46% protein, 3.19% lipid, 3.04 pg
per cell, with all cellular lipid treated as membrane lipid. This matters
because **Armitage's cell-line-specific entry is not a hepatic composition**.
His AREc2 is AREc32, which is MCF7-derived (Cellosaurus CVCL_1D32), and
carries 0.46% total lipid against the 3.19% measured for HepG2. Two of the
three datasets modelled here use HepG2 and the third HepaRG, so the published
parameterisation applies a breast-carcinoma composition throughout. (Armitage's
other entry, for the PPARγ line, reproduces Fischer et al.'s measured HEK293
numbers to the digits shown, which is evidently its source.)

**Table 7.** The Fischer et al. form, parameterised on measured PFAS
coefficients, nCF 5–8.

| cells | medium protein | bias | slope (decades/CF₂) |
|---|---|---|---|
| **measured HepG2 [12]** | **measured total [12]** | **5.56×** | **−0.100** |
| measured HepG2 [12] | albumin volume only | 17.09× | −0.033 |
| Armitage AREc2 (MCF7) | measured total [12] | 9.47× | −0.100 |
| Armitage AREc2 (MCF7) | albumin volume only | 27.64× | −0.044 |

Three results follow, and the first is a correction.

**The level is not fixed by measured parameters.** The faithful Escher-school
balance over-predicts by 5.56×, against 4.35× for Armitage as published. The
0.40× of Table 6 was a structural artefact, and the claim that measured
membrane partitioning moves the level toward unity does not survive. As it
stands, **no parameterisation tested in this work fixes the level**, which
strengthens rather than weakens the central negative result of §6.4.

**The chain-length drift is removed.** The slope falls from +0.316 to −0.100,
and over nCF 5–7 — where three of the four compounds sit — it is −0.007, i.e.
flat. The reason is visible in the coefficients: Fischer's measured albumin
series and Chen's measured membrane series rise with chain length at nearly
the same rate, so their ratio, which is what the predicted fraction depends
on (§6.5), is almost constant. Measured biopartitioning does solve the
chain-length problem; it simply does not solve the level problem.

**The two failures are cleanly separable.** The slope in Table 7 is unchanged
across a 6.9× change in cell lipid, a 3.3× change in cell volume and a 3×
change in medium protein — it is −0.100 in both rows with measured medium
protein and −0.03 to −0.04 in both without. Conversely the level moves by
5× across those same changes while the coefficients stay fixed. So:

> the chain-length behaviour is a property of the partition coefficients;
> the level is a property of the volumes and compositions.

Each error must be diagnosed with the corresponding class of parameter, and
neither can be diagnosed from the other. Both improvements available from
measurement work in the expected direction on the level — the measured
hepatic composition gains 1.7× over the MCF7 one, and the measured medium
protein gains 3.1× — and together they leave 5.56× standing.

Figure 2 places every parameterisation tested in this work on those two axes
alongside the lipid series that drives the horizontal one. The discriminator
for the slope is specifically the lipid coefficient: every variant whose
lipid–water ratio is measured falls at or below zero, Qin et al.'s included
since they take theirs from liposome measurements, and the only two that sit
on the positive side are the two that predict the lipid term from
log K<sub>OW</sub>.

![**A**, the lipid–water partition ratio against chain length. The orange line
is log D<sub>MW</sub> for the anion at pH 7.4 as IV-MBM EQP v2.0 evaluates it
(1.01·log K<sub>OW</sub> + 0.12, with the anion penalty d<sub>MW</sub> = 1.0);
the blue points and line are Chen et al.'s measurements on solid-supported
bilayers at pH 7.0 (SI Table S5, 12 perfluoroalkyl carboxylates); the dashed
green line is the regression Qin et al. fit to Ebert et al.'s liposome data
(Figure A1). The prediction rises more than twice as fast with chain length as
either measurement. **B**, every parameterisation tested in this work placed on
(chain-length slope, geometric-mean bias), recomputed over nCF 5–8 so both axes
are on one footing; PFBA is excluded throughout because it has no measured
albumin value. The star at (0, 1×) is a model with neither error. Substituting
measured coefficients moves points horizontally; changing cell and medium
volumes and compositions moves them vertically. Every variant whose lipid
coefficient is measured falls in the shaded band — Qin et al.'s included, since
they take theirs from liposome measurements — and the only two outside it are
the two that predict the lipid term from log K<sub>OW</sub>. Nothing tested
reaches the star.](../figures/fig22_mbm_separation.png)

Chen's coefficients carry the slope result in both Tables 6 and 7, so the
basis on which they are reported has to be established. It can be. Chen et al.
define K<sub>MW</sub> in their Supplementary equation S1 as

> K<sub>MW</sub> = (C<sub>spike</sub> − C<sub>PFAS,W</sub>)·V<sub>W</sub> /
> (C<sub>PFAS,W</sub>·M<sub>PL</sub>)

where C<sub>spike</sub> is the aqueous concentration in a lipid-free control at
the end of the incubation, C<sub>PFAS,W</sub> the aqueous concentration in the
sample, V<sub>W</sub> the water volume and M<sub>PL</sub> the quantity of
phospholipid membrane in the tube; the authors note that a unit conversion is
applied so that the reported coefficient is dimensionless. The measurement uses
a commercial solid-supported lipid membrane kit in which the membrane is
supplied as a stated *volume* — 0.067 to 2.166 µL of bilayer per tube in 500 µL
of PBS — so the lipid quantity is natively volumetric and the mass and volume
bases coincide at a phospholipid density of 1 g·cm<sup>−3</sup>. The reported
log K<sub>MW</sub> is therefore a dimensionless concentration ratio on a
volume-of-membrane to volume-of-water basis: the same basis as Armitage's
D<sub>MW</sub>, and the same density convention the cell-volume conversion used
here already adopts (10 ng per cell at 1 g·cm<sup>−3</sup>). Real phospholipid
densities of 1.01–1.06 g·cm<sup>−3</sup> bound the residual basis error at
≤0.03 log units, which is two orders of magnitude smaller than the
discrepancies in Tables 6 and 7. **The slope results therefore do not rest on
an unverified basis.**

Two independent checks support this reading. The tabulated values reproduce the
authors' own unweighted regression of log K<sub>MW</sub> on chain length — we
obtain +0.361 (R² = 0.983) for the perfluoroalkyl carboxylates and +0.374
(R² = 0.985) for the sulfonates against their stated +0.361 ± 0.014
(R² = 0.982) and +0.374 ± 0.019 (R² = 0.983) — confirming that the series was
transcribed correctly. And the absolute values sit in the published range for
PFAS membrane–water partitioning measured by the same technique.

The slope comparison is anchored at source: the measured series rises
+0.361 decades per CF<sub>2</sub> across the authors' full C4–C16 carboxylate
range (+0.44 over the nCF 5–8 window used in Table 6), against +0.78 for the
spLFER, which inherits +0.77 per CF<sub>2</sub> from the log K<sub>OW</sub>
series and multiplies it by 1.01. The spLFER over-steepens membrane
partitioning by a factor of about 2.6 per CF<sub>2</sub> added; compounded over
Chen's measured range this reaches three orders of magnitude, and the
predicted-to-measured ratio runs from 0.7 at PFBA to 0.001 at PFDoA.

A third, fully independent series agrees. Qin et al. [3] take their
lipid–water ratios from liposome measurements [20] and report the regression
they used to fill in PFBA: slope **+0.476** per fluorinated carbon for the
carboxylates. Their two sulfonate entries imply **+0.38**, against Chen's
**+0.374** for the same head group. So two laboratories, two lipid
preparations — phospholipid vesicles and solid-supported bilayers — and two
analytical routes place the lipid–water chain-length slope between +0.36 and
+0.48, while the spLFER inherits +0.77 from the log K<sub>OW</sub> series and
multiplies it by 1.01. **The over-steepening of the lipid term is not an
artefact of any one measurement.** It is the single most consistently
corroborated discrepancy in this reassessment, and it is the term that
represents almost the whole of the Armitage cell at long chain length
(Table 1).

One genuine caveat remains, and it is small. Chen et al. report K<sub>MW</sub>
at pH 7.0 while the balance is evaluated at pH 7.4, so a 0.4-unit
extrapolation is implicit. The authors measured at pH 6 and pH 7 and present
both, which bounds the pH sensitivity empirically, but they do not extend the
series to 7.4.

**The wider implication for the landscape** is that the two failure modes
documented here map onto two independent design axes. Measured biopartitioning
addresses the level; a saturable albumin term addresses the chain-length
behaviour. No published model currently occupies the corner that combines
both with a protein-realistic cell.

### 6.8 How much agreement did the source models claim?

A model's own validation sets a floor on what an external benchmark can
resolve, so it is worth asking what each model was shown to achieve before
being applied here.

Qin et al. [3] measured cell–water distribution ratios for their own HEK293H
line and compared them with the prediction from their Eq. 17. Reconstructing
that comparison from their published composition and parameters reproduces
their stated figures exactly:

**Table 8.** Qin's Eq. 17 against Qin's own measured D<sub>cell/w</sub>.

| | predicted | measured | Δ log | predicted/measured | % of log value |
|---|---|---|---|---|---|
| PFBA | 1.63 | 2.56 | −0.93 | 0.12× | −36% |
| PFOA | 3.69 | 3.27 | +0.42 | 2.65× | +13% |
| PFHxS | 3.54 | 3.23 | +0.31 | 2.05× | +10% |
| PFOS | 4.84 | 4.19 | +0.65 | 4.49× | +16% |

The paper describes this as agreement "within 15%" for PFOA, PFHxS and PFOS
and "around 35% lower" for PFBA. The final column shows what those percentages
are percentages of: the *logarithm* of the distribution ratio. On the quantity
that enters the mass balance the same residuals are a 2.0–4.5×
over-prediction for the three hydrophobic compounds and an 8.5×
under-prediction for PFBA. This is not a criticism of the measurement or of
the model's construction, and it does not affect the authors' conclusions,
which concern the free concentration in a medium where albumin carries more
than 98% of the capacity and the cell term is therefore nearly irrelevant. But
it has two consequences for the present work.

First, the cell-term residual is **ordered by chain length**: within the
carboxylates it gains +1.35 log units between nCF 3 and nCF 7, or +0.34 per
CF<sub>2</sub>. That is the same sign and nearly the same magnitude as the
drift found above in the Armitage cell term (+0.38 per CF<sub>2</sub>,
Table 5) — in a model whose cell parameters are measured rather than
predicted. The chain-length problem in the cell term is therefore not specific
to the spLFER approach.

Second, a 2–4.5× discrepancy in D<sub>cell/w</sub> is **inside the envelope
the model was published with**. An external benchmark cannot attribute a level
error smaller than that to the model rather than to its inputs, which bounds
how firmly the biases in Table 3 can be interpreted.

The same floor cannot be established for IV-MBM EQP v2.0 from the published
record, because the cell-protein sorption parameter is not documented there at
all (§3.1, §9.3).

---

## 7. The ligand:protein regime

The argument of this section is not new. Henneberger et al. [21] established
that for anionic chemicals "sorption to albumin cannot be predicted with a
simple partition since saturation of albumin binding sites for these chemicals
occurs at concentrations relevant for toxicity testing", as Proença et al. [17]
summarise it, and the same review notes that organic acids appear to have only
one binding site, which makes them saturate at relatively low concentrations.
Qin et al. [3] carried the argument into QIVIVE: concentration-dependent
protein binding makes the binding regime in bioassays differ substantially from
that in plasma. Their full text states the mechanism directly: the gap
between their nominal and free QIVIVE ratios arises "because protein binding is
highly concentration-dependent with a high-affinity specific binding at low
concentrations and lower non-specific binding at higher concentrations. Thus,
the low concentrations in human plasma are strongly bound to proteins, while
the high concentration needed to trigger effects in in vitro bioassays have a
higher fraction unbound." What is added here is the magnitude of the gap for
three datasets neither group analysed, and the demonstration that it is large
enough to invalidate the linearity assumption within a single dose series.
Figure 3 quantifies it.

Albumin offers approximately one high-affinity site per molecule for these
compounds (§3.2). The two 10% FBS assays contain 36.1 µM albumin and the 2% FBS
assay 7.2 µM, against nominal doses of 6.25–100 µM. The resulting nominal
ligand:albumin ratios are 0.28–2.77 and 0.87–13.85 mol/mol respectively.
Fischer et al. [9] worked at ≤ 0.004 mol/mol and estimate the ratio in human
serum at ≈ 0.0005. **The assays therefore operate 5,500–27,700× above the
physiological regime.**

The consequence for the models is direct. Solved self-consistently within the
wells, albumin site occupancy runs from 20% to 99% across the dose series, and
the apparent albumin distribution ratio falls to between 0.11 and 0.07 of its
trace-limit value — a sevenfold change in a quantity that a linear model treats
as a constant, within a single experiment (Fig. 3B).

This is why the Qin formulation requires its saturable term and why a linear
spLFER treatment cannot be correct here, whatever its merits for the neutral
organics it was calibrated on.

It also bounds what the exercise can deliver. Even an exactly correct model of
these wells would describe a binding regime that does not occur in vivo.
Improving the in vitro mass balance improves the description of the experiment;
it does not, on its own, improve the extrapolation, because the experiment is
in the wrong regime.

![**A**, the nominal ligand:albumin molar ratio against dose for the three
assays, with the physiological band marked. Albumin offers approximately one
high-affinity site per molecule for these compounds (§3.2), so a ratio of 1
means the sites are nominally filled. All three assays run three to four orders
of magnitude above the band. **B**, the consequence for the models: the
apparent albumin distribution ratio, relative to its trace limit, as a function
of freely dissolved concentration, with the positions the three assays occupy
marked. Within a single dose series the quantity a linear model treats as a
constant falls to between 0.11 and 0.07 of its trace value.](../figures/fig21_mbm_regime_gap.png)

---

## 8. A positive control: the same formalism in serum

The results above could be read as indicting the mass-balance approach itself.
A control argues otherwise.

Fischer et al. [9] measured, by one method in one laboratory, PFAS
distribution ratios to purified human serum albumin, to γ-globulin, and to
whole human serum, at a ligand:protein ratio of ≤ 0.004. The components and the
whole matrix are separate experiments, so predicting the latter from the former
is an out-of-sample test of the identical accounting used above.

For the 26 PFAS with all three measurements, the balance evaluated at
population-average serum composition and with no fitted parameter gives a
geometric-mean bias of 1.00× and a geometric-mean absolute error of 1.55×
(RMSE 0.262 log units; R² = 0.68 about the 1:1 line) for a quantity spanning
nearly two orders of magnitude. Allowing the two protein concentrations to
float returns 39.5 and 38.6 g/L against population values of 41 and 31 g/L and
does not reduce the RMSE at all — the model is not being supported by its
parameters. Globulin earns its term only above the ratio threshold the authors
themselves propose. Adding a phospholipid sink improves nothing.

**The accounting is therefore sound, and the failure is specific to the in
vitro system** — which has more phases, a vessel-wall term that is poorly
constrained, cell compositions the studies themselves did not measure, and a
ligand:protein regime in which the central linearity assumption does not
hold.

![The same mass-balance formalism applied to human serum, as a positive
control. Predicted against measured whole-serum–water distribution for the 26
PFAS in Fischer et al.'s Table S2, using their measured albumin coefficients,
physiological albumin and globulin concentrations, and no fitted parameter.
Geometric-mean bias 1.00×, geometric-mean |error| 1.55×, R² = 0.68. The
accounting is sound; the failure documented in §6 is specific to the in vitro
system.](../figures/fig19_fischer_serum_mbm.png)

---

## 9. Discussion

### 9.1 Implications for QIVIVE

Three findings bear on the use of these assays for reverse dosimetry.

First, the models over-predict cellular uptake by factors of 2–10 on data whose
own inter-study spread is a factor of 14. Neither model is currently accurate
enough to carry a QIVIVE calculation without an uncertainty statement of
comparable size.

Second, the residual is not attributable to any parameter for which a better
measurement exists. Substituting the best available value for albumin, for cell
partitioning and for cell protein each worsens the fit. The only parameter that
improves it is the one nobody measured.

Third, and most consequentially, the regime gap of §7 means that improving the
in vitro balance does not translate directly into a better oral equivalent
dose. Qin et al. [3] reached the same conclusion from the opposite direction,
showing that free-concentration-based QIVIVE ratios exceed nominal-based ones
precisely because of protein-content differences between blood and bioassays.
Our contribution is to show that the gap is large enough — three to four orders
of magnitude in ligand:protein ratio — that the saturation state of albumin
differs qualitatively, not merely quantitatively, between the two settings.

### 9.2 Reporting recommendations

Several of the difficulties encountered here are avoidable at no experimental
cost.

1. **Report the exposure medium composition, not the culture medium
   composition.** The single most consequential error in this work arose from
   a five-fold difference between the two in one study, stated in one sentence
   of its methods and absent from derived descriptions.
2. **Report the working volume.** One of the three studies does not, leaving a
   parameter that enters the balance linearly to be assumed.
3. **Report the cell volume, or cells per well together with protein per
   cell.** Cell volume is the only parameter that improved the fit here and
   the only one none of the studies measured.
4. **Report the ligand:protein molar ratio.** It determines whether a measured
   binding constant is in the trace regime and therefore whether it can be
   transferred to another setting at all.
5. **State the normalisation basis** — administered or recovered — since the
   two differ by whatever was lost to the vessel.
6. **Report the mass balance where it was measured.** Two of the three studies
   took a t = 0 medium sample explicitly for this purpose; neither reports the
   result. Doing so would constrain the vessel-wall term directly rather than
   by difference.

### 9.3 Limitations of this reassessment

The Armitage implementation depends on one parameter — the sorptive capacity
assigned to bulk cell protein — that is not documented in the published record
(§3.1). We used the `httk` convention. Having read both papers, we can now
bound what this costs. At the `httk` value the protein phase carries 0.05–0.25%
of the cell term, so no result here is sensitive to the number itself; what the
papers do not settle is whether the authors' own tool makes the same
octanol-equivalence *choice*, and that choice does have leverage (§6.4).
Deciding it would require the Excel/VBA implementation, which is available from
the authors on request. The structural argument is unaffected either way: a
cell whose protein term is small is represented as lipid, and a lipid
representation carries too steep a chain-length dependence — which three
independent measured lipid series now corroborate (§6.7).

Qin's parameters exist for four compounds; two of the five modelled here are
interpolated and one extrapolated. The extrapolated compound, PFNA, is among
the worst-fitting. The authors scope their own claim narrowly — they
demonstrated the model's validity "only for four anionic PFAS, one type of
bioassay medium and one cell line" — while arguing that prior work on other
organic anions supports transfer to other setups "if the lipid and protein
content has been quantified". Applying the model to three other datasets with
quantified composition is therefore the use the authors anticipate; extending
it along the chain-length axis to compounds they did not measure is not, and
the interpolation used here for PFHxA and PFHpA and the extrapolation for PFNA
have no published warrant beyond the lipid term.

The measured-membrane result of §6.7 rests on a coefficient reported at pH 7.0
and applied here at pH 7.4; the reporting basis itself is resolved (§6.7) but
the 0.4-unit pH extrapolation is not.

The uptake datasets are small, mutually inconsistent, and partly not on a
common normalisation basis. Rosenmai's working volume is assumed. Conclusions
about absolute bias are correspondingly weaker than conclusions about
*direction* and about *structure*, which is why we have emphasised the latter.

Finally, of the fourteen models Proença et al. [17] catalogue, three are
evaluated here — the two under assessment and a reconstruction of Fischer et
al. [12]. The remaining eleven are not, and the comparison here should not be
read as a ranking of the field. The Fischer reconstruction in particular is
built from that paper's published equations and its Table 2 system parameters,
not from the authors' own implementation, and substitutes measured PFAS
coefficients for the PP-LFER, 3D-QSAR and COSMOmic predictions they would use;
it is the design tested on the best available PFAS data, which is not the same
thing as the model as its authors would run it.

---

## 10. Conclusions

1. Both mass-balance models over-predict PFAS cellular uptake in human liver
   cell lines, by 3.9× (IV-MBM EQP v2.0) and 7.5× (Qin et al.) in geometric
   mean, once exposure conditions are read from primary sources.
2. Exposure media differ between the three available datasets in ways that
   dominate the model comparison and are not recoverable from secondary
   descriptions. One assay was dosed at one-fifth the serum concentration
   generally assumed for it.
3. Every substitution of a best-available measured parameter worsens the fit.
   Albumin is the most influential term in the balance and is not the source of
   the error.
4. The residual separates into a level component and a chain-length
   component, and the two are governed by different classes of parameter: the
   chain-length behaviour by the partition coefficients, the level by the
   volumes and compositions. Each is invariant to large changes in the other,
   so neither error can be diagnosed with the wrong class of parameter.
5. Measured biopartitioning removes the chain-length error and not the level
   error. A faithful reconstruction of the measured-parameter linear design is
   flat over nCF 5–7 and still over-predicts by 5.6×. The cause of the
   chain-length error is specifically that single-descriptor prediction
   over-steepens the lipid term: three independent measured series place the
   lipid–water chain-length slope at +0.36 to +0.48 decades per CF₂ against
   +0.78 predicted.
6. All three assays operate three to four orders of magnitude above the
   physiological ligand:protein ratio, with albumin between 20% and 99%
   saturated. The linearity assumption underlying spLFER-based balances does
   not hold in this regime.
7. The same formalism applied to human serum reproduces measured whole-serum
   partitioning for 26 PFAS with no fitted parameter. The accounting is sound;
   the in vitro system is the problem.
8. Consequently, better in vitro dosimetry yields a better description of the
   experiment, but not by itself a better oral equivalent dose. Closing that
   gap requires assays run nearer the physiological binding regime, not better
   models of assays that are not.

---

## Data and code availability

| | |
|---|---|
| model implementations (§3.1, §3.2) | `tk_review/scripts/mbm_published_equations.py` |
| exposure conditions, source-verified (§4.2) | `tk_review/scripts/mbm_medium_composition.py` |
| chain-length decomposition (§6.5) | `tk_review/scripts/mbm_chain_length_residual.py` |
| parameter substitutions (§6.4) | `mbm_fischer_albumin.py`, `mbm_cell_term.py`, `mbm_cell_protein.py` |
| coefficient substitutions, Armitage structure (§6.7.1) | `tk_review/scripts/mbm_measured_parameter_variant.py` |
| the Fischer et al. form (§6.7.2) | `tk_review/scripts/mbm_fischer_escher_balance.py` |
| Qin's model against Qin's measurement (§6.8) | `tk_review/scripts/mbm_qin_internal_check.py` |
| IV-MBM v1.0 against v2.0 (§3.1) | `tk_review/scripts/mbm_armitage_v1_vs_v2.py` |
| serum control (§8) | `tk_review/scripts/mbm_fischer_serum.py` |
| Figure 1 | `tk_review/scripts/make_mbm_correction_figure.py` |
| Figure 2 | `tk_review/scripts/make_mbm_separation_figure.py` |
| Figure 3 | `tk_review/scripts/make_mbm_regime_figure.py` |
| Figure 4 | `tk_review/scripts/make_mbm_fischer_serum_figure.py` |
| primary sources held locally | `literature/fulltext/` |

Extended derivations and the full approximation inventory are in
`report/MBM_PRIMER.md`; a non-technical exposition is in
`report/MBM_UNDERSTANDING.md`; the detailed working memo, including material
not carried into this manuscript, is `report/MBM_Qin_vs_Armitage.md`.

---

## References

1. Armitage JM, Sangion A, Parmar R, Looky AB, Arnot JA. Update and evaluation
   of a high-throughput in vitro mass balance distribution model: IV-MBM EQP
   v2.0. *Toxics* 2021;9:315. doi:10.3390/toxics9110315
2. Armitage JM, Wania F, Arnot JA. Application of mass balance models and the
   chemical activity concept to facilitate the use of in vitro toxicity data
   for risk assessment. *Environ Sci Technol* 2014;48(16):9770–9779.
   doi:10.1021/es501955g
3. Qin W, Henneberger L, Huchthausen J, König M, Escher BI. Role of
   bioavailability and protein binding of four anionic perfluoroalkyl
   substances in cell-based bioassays for quantitative in vitro to in vivo
   extrapolations. *Environ Int* 2023;173:107857.
   doi:10.1016/j.envint.2023.107857
4. Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L. Species difference?
   Bovine, trout, and human plasma protein binding of per- and polyfluoroalkyl
   substances. *Environ Sci Technol* 2024;58(23):9954–9966.
   doi:10.1021/acs.est.3c10824
5. Rosenmai AK, Ahrens L, le Godec T, Lundqvist J, Oskarsson A. Relationship
   between peroxisome proliferator-activated receptor alpha activity and
   cellular concentration of 14 perfluoroalkyl substances in HepG2 cells.
   *J Appl Toxicol* 2018;38:219–226. doi:10.1002/jat.3515
6. Amstutz VH, Cengo A, Gehres F, Sijm DTHM, Vrolijk MF. Investigating the
   cytotoxicity of per- and polyfluoroalkyl substances in HepG2 cells: a
   structure-activity relationship approach. *Toxicology* 2022;480:153312.
   doi:10.1016/j.tox.2022.153312
7. Fragki S, Louisse J, Bokkers B, Luijten M, Peijnenburg A, Rijkers D,
   Piersma AH, Zeilmaker MJ. New approach methodologies: a quantitative in
   vitro to in vivo extrapolation case study with PFASs. *Food Chem Toxicol*
   2022;172:113559. doi:10.1016/j.fct.2022.113559
8. Louisse J, Fragki S, Rijkers D, et al. Determination of in vitro hepatotoxic
   potencies of a series of perfluoroalkyl substances (PFASs) based on gene
   expression changes in HepaRG liver cells. *Arch Toxicol* 2023;97:1113–1131.
   doi:10.1007/s00204-023-03450-2
9. Fischer FC, et al. *Environ Sci Technol* 2024;58(2):1055–1063.
   doi:10.1021/acs.est.3c07415
10. Chen R, Muensterman D, Field J, Ng C. Deriving membrane-water and
    protein-water partition coefficients from in vitro experiments for per- and
    polyfluoroalkyl substances (PFAS). *Environ Sci Technol* 2025;59(1):82–91.
    doi:10.1021/acs.est.4c06734
11. Ryu S, Burchett W, Zhang S, Jia X, Modaresi SMS, Agudelo Areiza J,
    Rodrigues D, Zhu H, Sunderland EM, Fischer FC, Slitt AL. Unbound fractions
    of PFAS in human and rodent tissues: rat liver a suitable proxy for
    evaluating emerging PFAS? *Environ Sci Technol* 2024;58(33):14641–14650.
    doi:10.1021/acs.est.4c04050
12. Fischer FC, Henneberger L, König M, Bittermann K, Linden L, Goss K-U,
    Escher BI. Modeling exposure in the Tox21 in vitro bioassays. *Chem Res
    Toxicol* 2017;30(5):1197–1208. doi:10.1021/acs.chemrestox.7b00023
13. Fischer FC, Cirpka OA, Goss K-U, Henneberger L, Escher BI. Application of
    experimental polystyrene partition constants and diffusion coefficients to
    predict the sorption of neutral organic chemicals to multiwell plates.
    *Environ Sci Technol* 2018;52:13511–13522.
14. Henneberger L, Mühlenbrink M, Fischer FC, Escher BI. C18-coated solid-phase
    microextraction fibers for the quantification of partitioning of organic
    acids to proteins, lipids, and cells. *Chem Res Toxicol* 2019;32:168–178.
15. Huchthausen J, Mühlenbrink M, König M, Escher BI, Henneberger L.
    Experimental exposure assessment of ionizable organic chemicals in in
    vitro cell-based bioassays. *Chem Res Toxicol* 2020;33:1845–1854.
16. Endo S, Brown TN, Goss K-U. General model for estimating partition
    coefficients to organisms and their tissues using the biological
    compositions and polyparameter linear free energy relationships. *Environ
    Sci Technol* 2013;47:6630–6639.
17. Proença S, Escher BI, Fischer FC, Fisher C, Grégoire S, Hewitt NJ,
    Nicol B, Paini A, Kramer NI. Effective exposure of chemicals in in vitro
    cell systems: a review of chemical distribution models. *Toxicol In Vitro*
    2021;73:105133. doi:10.1016/j.tiv.2021.105133
18. Stadnicka-Michalak J, Tanneberger K, Schirmer K, Ashauer R. Measured and
    modeled toxicokinetics in cultured fish cells and application to in
    vitro–in vivo toxicity extrapolation. *PLoS One* 2014;9(3):e92303.
    doi:10.1371/journal.pone.0092303
19. Paini A, Sala Benito JV, Bessems J, Worth AP. From in vitro to in vivo:
    integration of the virtual cell based assay with physiologically based
    kinetic modelling. *Toxicol In Vitro* 2017;45(Pt 2):241–248.
    doi:10.1016/j.tiv.2017.06.015
20. Ebert A, Allendorf F, Berger U, Goss K-U, Ulrich N. Membrane/water
    partitioning and permeabilities of perfluoroalkyl acids and four of their
    alternatives and the effects on toxicokinetic behavior. *Environ Sci
    Technol* 2020;54(8):5051–5061. doi:10.1021/acs.est.0c00157
    — Source of the liposome–water ratios used by ref. [3]; cited here via
    [3], whose Table 2 and Figure A1 regression are the values used.

21. Henneberger L, Mühlenbrink M, Heinrich DJ, Teixeira A, Nicol B,
    Escher BI. Experimental validation of mass balance models for in vitro
    cell-based bioassays. *Environ Sci Technol* 2020;54(2):1120–1127.
    doi:10.1021/acs.est.9b06144
    — Cited via ref. [17] for the saturation of albumin binding by anionic
    chemicals at toxicologically relevant concentrations; not independently
    retrieved.

**Provenance of references.** Read at source, including supplementary
material where cited: 1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 17. Bibliographic record
verified: 4, 11, 18, 19. Cited via ref. [3], which is the source of the
values used: 20. Cited via ref. [17] for a statement attributed to it, not
independently retrieved: 21. Cited as methodological context from the
reference list of ref. 1, not independently retrieved: 13, 14, 15, 16.

Both models under assessment have now been read at source. For ref. [3] the
full text supplied the balance equations, the cell composition and the authors'
own cell-term validation, all used in §3.2 and §6.8. For ref. [1] the
supplementary information supplied every parameter except the cell-protein
sorption term; ref. [2] establishes that no such term existed in v1.0 and that
the 3.5% convention is not in the published record (§3.1, §9.3). Only the
main text of ref. [2] was read; its supplementary sections S1–S9 were not
retrieved, and the statements made here rest on its Eq. 1 and the surrounding
model description, which are in the main text. Ref. [12] was read in full; only
its main text was retrieved, and the cell and medium compositions used here are
from its Table 2, not its Tables S1–S2, which carry the full set of cell lines
and assays. The AREc32 provenance is from Cellosaurus (CVCL_1D32) and Wang,
Hayes & Wolf 2006, not from ref. [1], which does not state the line's
derivation.
