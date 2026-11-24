# Species differences in perfluoroalkyl acid elimination are carried by the renal tubule, not by plasma protein binding: an integrated analysis of albumin binding, liver fatty acid-binding protein and renal transporter expression

**[AUTHORS — to be supplied]**
**[AFFILIATIONS — to be supplied]**
**[CORRESPONDING AUTHOR — to be supplied]**

*Draft 1. Target journal not yet selected; see §Notes for the decisions that
depend on it.*

---

## Abstract

Perfluoroalkyl acid (PFAA) elimination half-lives differ between humans and
rodents by up to three orders of magnitude, and the difference is commonly
attributed to species differences in plasma protein binding. We test that
attribution quantitatively. Because the half-life of a one-compartment system
is ln2·V_d/CL and renal clearance is the product f_u × GFR × T, where T is the
fraction of the filtered load surviving the renal tubule, the human/rat
half-life ratio factorises exactly into four multiplicative terms —
distribution, binding, glomerular filtration and tubular handling — whose
logarithms partition the observed difference. Applying this to six PFAAs for
which volume of distribution and clearance are both measured in both species,
half-life ratios span 6- to 1792-fold, of which binding accounts for a mean of
3% of the log ratio and the tubular term for 80%. The binding term is not
merely small but bounded: across fourteen PFAAs measured in human and rat
plasma in a single laboratory, the rat/human unbound-fraction ratio spans only
0.76 to 2.03, so binding is arithmetically incapable of producing a
thousand-fold difference in half-life. The glomerular term is a species
constant (5.6-fold) and therefore cannot explain why the half-life ratio varies
between compounds, which it does by 300-fold. We then show why the tubular term
is the one that is not conserved. In kidney membrane fractions assayed together,
the secretory transporters OAT1 and OAT3 are present in human, rat and mouse and
differ only in abundance, whereas the reabsorptive transporters are not shared
at all: human OAT4 is undetectable in rat and mouse kidney, and rat Oatp1a1 is
undetectable in human. Since reabsorption determines the long half-lives, rodent
reabsorption data carries magnitude but not mechanism. We further find that
plasma protein binding cannot be used to derive an interspecies scaling factor
for an unrelated reason: across 288 albumin–water partition coefficients for 63
PFAAs from 20 sources, converted to a common basis and audited against the
primary publications, the median within-compound range is 2.38 log units, so the
disagreement between laboratories measuring one species exceeds the difference
between species. Liver fatty acid-binding protein, the standing explanation for
hepatic accumulation, is indistinguishable from albumin per unit protein mass
(median difference −0.12 log units across 21 compounds measured against both),
a conclusion supported independently by a knockout study and by tissue
homogenate binding data. We document eight errors in the primary and regulatory
literature, two of which propagate into published reabsorption estimates by
8- and 40-fold.

**Keywords:** perfluoroalkyl substances; toxicokinetics; interspecies
extrapolation; renal tubular reabsorption; serum albumin; liver fatty
acid-binding protein; organic anion transporters; elimination half-life

---

## 1. Introduction

Several perfluoroalkyl acids have human elimination half-lives measured in
years, and this persistence is the principal driver of their regulatory
significance. The same compounds are eliminated by rodents in days. Any
health assessment that uses rodent toxicology to set a human reference value
must cross that gap, and the size of the correction depends on what is
believed to cause it.

The most frequently invoked cause is plasma protein binding. PFAAs are
extensively bound to serum albumin; only the unbound fraction is available for
glomerular filtration; and binding is known to vary with chain length. The
inference that species differ in half-life because they differ in binding is
therefore mechanistically plausible, and it has the practical attraction that
binding is measurable in vitro, so an interspecies scaling factor could in
principle be obtained without animal work.

This paper tests that inference and finds it cannot be sustained. The argument
has three steps.

First, the attribution is **quantitatively testable**, because the half-life
decomposes exactly into separately measurable terms (§2.1). When the
decomposition is carried out, binding accounts for a few percent of the
observed species difference and the renal tubule accounts for most of it
(§3.5).

Second, the binding term is **bounded by direct measurement**. In the one
study that measured human and rodent plasma in a single laboratory, the
rat/human unbound-fraction ratio never exceeds about twofold in either
direction. A term bounded at twofold cannot explain a thousandfold difference,
and no improvement in the binding data can change this, because the bound is
on the quantity rather than on our knowledge of it (§3.5.3).

Third, the term that does carry the difference is the one that is
**mechanistically not conserved**. Measured in a single assay across three
species, the secretory transporters are shared and differ only in abundance,
while the reabsorptive transporters are different proteins in human and rat
(§3.3).

We also report two findings that bear on how this literature should be used.
The binding literature disagrees with itself, within one species, by more than
species differ from one another, which forecloses the estimation of a binding
scaling factor independently of the argument above (§3.1.2). And liver fatty
acid-binding protein, the standing explanation for the hepatic accumulation of
PFAAs, has no affinity advantage over the plasma carrier it is supposed to
outcompete (§3.2).

---

## 2. Methods

### 2.1 The governing identity and its decomposition

For a compound eliminated with first-order kinetics,

```
                        t½ = ln2 · V_d / CL                              (1)
```

Renal clearance is the product of the filterable fraction, the filtration
rate, and whatever the tubule does to the filtrate:

```
                        CL_R = f_u · GFR · T                             (2)
```

where `f_u` is the unbound fraction in plasma, `GFR` the glomerular filtration
rate, and `T` the dimensionless **tubular factor** — the fraction of the
filtered load that survives the tubule and is excreted. Writing `FR` for
fractional reabsorption and `S` for secretory clearance,

```
                   T = (1 − FR) + S / (f_u · GFR)                        (3)
```

so that T < 1 denotes net reabsorption, T = 1 filtration alone, and T > 1 net
secretion. Equation (3) contains two unknowns and is not separable from in
vivo data; T is therefore carried as a lumped quantity throughout, and
statements of the form "x% reabsorbed" should be read as "x% *net*
reabsorbed" (§5.1).

Combining (1) and (2) for one compound in two species gives an exact
factorisation of the half-life ratio:

```
    t½,human     V_d,human     f_u,rat      GFR_rat      T_rat
    ────────  =  ─────────  ×  ───────   ×  ───────   ×  ─────           (4)
    t½,rat       V_d,rat       f_u,human    GFR_human    T_human
                 \_________/   \_________/  \_________/  \_____/
                 DISTRIBUTION    BINDING    GLOMERULAR   TUBULE
```

On a logarithmic scale these four terms add, so the share of the observed
difference carried by each is well defined.

Two properties of (4) should be stated before the results. **The glomerular
term is a species constant**: it does not depend on the compound, so it
contributes the same amount to every half-life ratio and cannot explain
variation between compounds. And **the tubular term is a residual**: since
both clearances in this analysis are total rather than renal, the fourth
factor absorbs any species difference in non-renal elimination as well as the
tubule. The identity is therefore exact by construction, and the analysis
partitions explanatory burden rather than testing the identity itself. We
return to this in §5.2.

### 2.2 Data compilation

**Protein binding.** 288 albumin–water partition coefficients for 63 PFAAs
from 20 primary experimental sources, and 91 coefficients for 32 PFAAs against
liver fatty acid-binding protein (L-FABP) from ten sources, were compiled and
converted to a common basis — log *K*, litres of water per kilogram of protein
— through the identity

```
            log K = 4.1774 + log₁₀(n) − log₁₀(K_d / µM)                  (5)
```

with a per-species bridge constant of 4.8475 for human L-FABP. All 181 rows
carrying both a dissociation constant and a stoichiometry satisfy (5) exactly.
Each value was audited against its primary publication.

**Plasma unbound fractions.** Four independent human datasets were compiled
[Smeltz 2023; Fischer 2024; Argoul 2025; Ryu 2024a]. For interspecies ratios
only Ryu 2024a is used, being the only study measuring human, rat and two
mouse strains in one laboratory; in a ratio, a systematic method bias divides
out.

**Transporter expression and kinetics.** 305 rows of transporter localisation
and direction, 40 rows of quantitative abundance, and 52 intrinsic clearances
were compiled. Expression values are taken from a single assay across three
species [Thakur 2024] wherever possible.

**Clearance and volume of distribution.** For the decomposition, human values
are from Abraham 2024 — one male volunteer dosed with isotopically labelled
compounds, with V_d and CL both measured and the identity ln2·V_d/CL closing
at a consistent implied body weight of 82.9 kg across all fifteen compounds.
Rat values are male, measured, one study per compound [Kim 2016; Huang 2019;
Chengelis 2009; Ohmori 2003]. GFR is taken as 2570 mL/kg/day (human) and
14,400 mL/kg/day (rat).

### 2.3 Admission rules

Two rules are enforced programmatically rather than by care.

**No computed quantity mixes species.** The unbound fraction, filtration rate
and clearance entering any T must originate in one species; a cross-species
join raises an error rather than a warning, and every delivered row records
the binding matrix used. Cross-*strain* joins within a species are permitted
and recorded.

**A primary publication supersedes any compilation of it**, with the
compilation's value retained and marked superseded rather than deleted. This
rule produced the error findings of §3.6.

Where only total clearance is available, T is admitted as an upper bound,
since renal clearance is a component of total clearance. A bound below unity
establishes net reabsorption regardless of how the clearance divides between
routes.

### 2.4 Reproducibility

All tables are generated from the primary data by idempotent scripts, verified
byte-identical across consecutive runs. Data files, scripts and a notebook
reproducing every figure are listed in §Data and code availability.

---

## 3. Results

### 3.1 Serum albumin

#### 3.1.1 Chain length

Albumin affinity rises with chain length to a maximum at perfluorononanoic
acid (9 total carbons) at +0.52 log units per carbon from PFBA, and declines
thereafter (Table 1). The decline is consistent with steric constraint in
albumin's binding pockets and with the falling aqueous solubility of the
longest homologues.

**Table 1.** Perfluoroalkyl carboxylate homologous series, albumin.

| total carbons | compound | median log *K* | *n* |
|---|---|---|---|
| 4 | PFBA | 1.68 | 10 |
| 5 | PFPeA | 2.23 | 6 |
| 6 | PFHxA | 2.77 | 15 |
| 7 | PFHpA | 3.54 | 12 |
| 8 | PFOA | 3.75 | 25 |
| **9** | **PFNA** | **4.29** | **18** |
| 10 | PFDA | 3.77 | 14 |
| 11 | PFUnDA | 3.56 | 11 |
| 12 | PFDoDA | 2.98 | 7 |
| 13 | PFTriDA | 2.68 | 5 |
| 14 | PFTeDA | 3.46 | 4 |

#### 3.1.2 The literature disagrees with itself more than species differ

Across the 32 compounds with three or more independent values, the median
within-compound range is **2.38 log units** and the maximum is 4.22.
Selecting a single published constant for a toxicokinetic model is therefore a
decision with a 240-fold range of defensible outcomes.

Method is the principal axis of this disagreement: equilibrium dialysis and
differential scanning fluorimetry differ by 2.35 log units in median, in the
direction predicted by the thermal denaturation that the latter measures.
These are two quantities rather than two estimates of one, and averaging
across them is not meaningful.

Human and bovine albumin differ in median by 0.90 log units, but the
comparison is confounded with method and compound coverage and should not be
read as an interspecies difference in affinity.

> **Consequence.** An interspecies scaling factor for PFAA–albumin binding
> cannot be estimated from this literature, because the quantity it would
> correct for is smaller than the disagreement about it. This conclusion is
> independent of, and prior to, the quantitative argument of §3.5.

### 3.2 Liver fatty acid-binding protein

L-FABP does not enter equation (2); it is intracellular and hepatic. It is
examined here because it is the standing explanation for where PFAAs
accumulate instead of being excreted, and because that explanation is testable
against the albumin compilation on a common basis.

**Three independent lines agree that the explanation fails.**

First, **per unit protein mass the two proteins are not distinguishable.**
Across the 21 compounds measured against both, the median difference is
**−0.12 log units**, L-FABP marginally the weaker — far inside the dispersion
of either compilation. This does not imply the two are interchangeable in vivo,
since albumin is present in plasma at ~600 µM and L-FABP in hepatocyte cytosol
at 100–400 µM in a different compartment; but it does remove any *intrinsic
affinity advantage*, and an explanation requiring L-FABP to outcompete the
plasma carrier cannot proceed without one.

Second, **the one in vivo test is negative**: deleting L-FABP changes the
mouse liver partition coefficient by under 20% for 15 of 19 PFAAs.

Third, and independently, **tissue binding does not track L-FABP abundance.**
In whole-tissue homogenates, PFAA binding in liver was within twofold of
binding in heart and kidney for 14 of 16 compounds, despite the liver's
enrichment in L-FABP [Ryu 2024b], whose authors draw the same inference.

Part of the residual dispersion in the L-FABP literature is structural rather
than experimental: the protein has two binding sites, and displacement assays
using a secondary-site probe and a primary-site probe are not estimating the
same quantity. Stratifying on probe site is an available improvement.

As with albumin, the apparent rat–human difference in L-FABP affinity is set
by the choice of human comparator — five human sources disagree among
themselves by 1.21 log units, with PFOA alone differing 21-fold between two
human laboratories — so no scaling factor can be derived here either.

### 3.3 Renal transporters

#### 3.3.1 Localisation is agreed; substrate identity is not

Of 305 directional rows, 299 resolve to a reabsorption-or-secretion call
across 21 distinct transporters, with **zero genuine direction conflicts.**
Papers that disagree sharply about which compounds a transporter carries do
not disagree about where it sits. Localisation should be treated as settled
and substrate lists as contested.

#### 3.3.2 The reabsorptive machinery is not shared

**Table 2.** Renal transporter abundance in kidney membrane fractions,
pmol/mg membrane protein, single assay [Thakur 2024].

| transporter | membrane | net process | human | rat | mouse |
|---|---|---|---|---|---|
| Oat1 | basolateral | secretory | 4.8 ± 1.2 | **32.8 ± 13.7** | 19.1 ± 7.1 |
| Oat3 | basolateral | secretory | 4.3 ± 1.1 | 6.3 ± 1.1 | 6.5 ± 3.6 |
| **OAT4** | apical | **reabsorptive** | **0.5 ± 0.2** | **not detected** | **not detected** |
| **Oatp1a1** | apical | **reabsorptive** | **not detected** | **1.9 ± 2.1** | 3.1 ± 2.8 |
| URAT1 | apical | — | 3.1 ± 1.2 | 2.2 ± 0.7 | 10.8 ± 1.1 |

The secretory side is conserved in kind and differs in degree: both species
possess OAT1 and OAT3, the rat having 6.8-fold more OAT1. The reabsorptive
side is not the same machinery at all.

> Because reabsorption determines the long half-lives (§3.4), **the half of
> the renal system that governs persistence is the half that is not
> conserved.** A rat-to-human extrapolation of reabsorption is an
> extrapolation between two different proteins performing similar functions:
> defensible as an empirical scaling, not as a mechanistic one.

#### 3.3.3 Abundance must gate relevance

Applying an expression gate removes the three highest intrinsic clearances in
the compilation, all for BCRP, including 279 µL/min/mg for PFOA — the highest
value anywhere in it. Human renal BCRP is not detected by either laboratory
that has looked.

The same distinction resolves two standing disagreements: ASBT and OSTα/β were
found to transport PFAAs in transfected cells by one group and not by another,
and both results are correct, because a transfected-cell result reports what a
protein *can* carry rather than whether it is present to do so.

Two available scaling methods for OAT1 and OAT3 disagree 11-fold and in
opposite directions, and OAT4 has no scaling factor at all, so in vitro
clearance ratios cannot presently be taken in vivo (§5.3).

### 3.4 The tubular factor across species

Resolving the clearance, binding and filtration data onto equation (2) yields
T for 52 species–compound pairs. In the human, PFOA through PFUnDA exceed
**98% net reabsorbed** on every available binding source, and the C6 and C8
sulfonates reach 97.0%. The rat/human tubular factor is of order 10²
(22- to 161-fold across the available clearance sources), against a glomerular
factor of 6.8-fold.

Net secretion is not a corner case: mouse PFHxA secretes at T = 2.23, and the
female rat secretes PFOA at T = 24.6 against a male value of 0.42–1.07 — a
change of sign within one species, attributable to a 114-fold sex difference
in the abundance of the reabsorptive transporter Oatp1a1.

### 3.5 The half-life decomposition

#### 3.5.1 The partition

Applying equation (4) to the six PFAAs for which V_d and CL are measured in
both species gives Table 3.

**Table 3.** Human/rat half-life ratio and its four factors. Human: Abraham
2024, male, measured. Rat: male, measured, source per compound.

| compound | total C | t½ human (d) | t½ rat (d) | ratio | distribution | binding | glomerular | **tubule** |
|---|---|---|---|---|---|---|---|---|
| PFBS | 4 | 50.6 | 0.16 | 322 | 0.73 | 1.48 | 5.60 | **53.0** |
| PFHxA | 6 | 1.4 | 0.23 | 6 | 0.13 | 0.94 | 5.60 | **9.3** |
| PFHxS | 6 | 1635 | 32.8 | 50 | 0.56 | 1.09 | 5.60 | **14.6** |
| PFHpA | 7 | 152 | 0.08 | 1792 | 0.65 | 1.02 | 5.60 | **480.5** |
| PFOA | 8 | 2011 | 1.64 | 1228 | 1.08 | 1.22 | 5.60 | **166.8** |
| PFOS | 8 | 1211 | 28.7 | 42 | 0.40 | 1.45 | 5.60 | **13.1** |

**Table 4.** Share of the log₁₀ half-life ratio carried by each factor (%).

| compound | distribution | binding | glomerular | **tubule** |
|---|---|---|---|---|
| PFBS | −5 | 7 | 30 | **69** |
| PFHxA | −112 | −3 | 94 | **121** |
| PFHxS | −15 | 2 | 44 | **68** |
| PFHpA | −6 | 0 | 23 | **82** |
| PFOA | 1 | 3 | 24 | **72** |
| PFOS | −25 | 10 | 46 | **69** |
| **mean** | **−27** | **3** | **44** | **80** |

**Binding accounts for a mean of 3% of the log half-life ratio; the tubular
term accounts for 80%.** Distribution works *against* the observed difference
for five of six compounds: on volume of distribution alone the human would
eliminate these compounds faster than the rat, and the other terms must
overcome that.

#### 3.5.2 The glomerular term cannot explain the variation

The glomerular factor is 5.60 for every compound, being the ratio of two
species constants. The half-life ratios it must help explain vary from 6 to
1792 — a 300-fold spread between compounds. **A constant cannot account for a
varying quantity**, so whatever distinguishes PFHpA (1792-fold) from PFHxA
(6-fold) resides in the tubular term, which indeed varies from 9.3 to 480.

#### 3.5.3 The binding term is bounded, not merely small

Across all fourteen PFAAs measured in human and rat plasma in a single
laboratory, the **rat/human unbound-fraction ratio spans 0.76 to 2.03**
(median 1.16).

> This is the central quantitative result. The binding term is bounded at
> about twofold **by direct measurement**, while the half-life ratios to be
> explained span 6- to 1792-fold. Binding is arithmetically incapable of
> producing the observed species difference in either direction, and **no
> refinement of binding data can alter this, because the bound is on the
> quantity itself rather than on our knowledge of it.**

#### 3.5.4 Sensitivity

The conclusion does not depend on the single-subject human dataset. Repeating
the decomposition with population values from a Bayesian hierarchical analysis
of serum decay [Chiu 2022] shifts the partition but not the conclusion: the
binding share remains under 10% on every combination tested, and the tubular
residual ranges from 57% to 80% (Table 5).

**Table 5.** Sensitivity to the human data source.

| compound | human source | t½ (y) | ratio | binding % | tubule % |
|---|---|---|---|---|---|
| PFOA | Abraham 2024 (n = 1, measured) | 5.5 | 1228 | 3 | 72 |
| PFOA | Chiu 2022 (population) | 8.6 | 1915 | 3 | 57 |
| PFOS | Abraham 2024 | 3.3 | 42 | 10 | 69 |
| PFOS | Chiu 2022 | 9.2 | 117 | 8 | 60 |
| PFHxS | Abraham 2024 | 4.5 | 50 | 2 | 69 |
| PFHxS | Chiu 2022 | 22.0 | 245 | 2 | 62 |

### 3.6 Errors identified in the primary and regulatory literature

Applying the rule that a primary publication supersedes any compilation of it
(§2.3) identified eight discrepancies, of which two propagate into published
reabsorption estimates.

**Table 6.** Errors found by reading primary sources.

| # | error | magnitude |
|---|---|---|
| 1 | A regulatory toxicological profile misreports rat total clearances transcribed from the primary study: PFDA printed as 207/140 against the paper's 5.2/5.3 mL/kg/day, male PFOA as 135 against 50.5 | **40×**, 26×, 2.7× |
| 2 | A published health goal assumes an unbound fraction of 0.02 for every species and compound, uncited; measured values are 0.0025 (human) and 0.0030 (rat), so every derived reabsorption percentage is wrong in the denominator | ~**8×** |
| 3 | A units defect in a published K_m/V_max column | — |
| 4 | Transposition of two transporters between a paper's text and abstract | — |
| 5 | A citation for reported unbound fractions that does not resolve | — |
| 6 | A paper misquoting its own cited source in its Discussion | — |
| 7 | One compound reported twice, 2.8-fold apart, within one study | 2.8× |
| 8 | Silent loss of the micro sign during PDF text extraction, detected twice | **1000×** |

Recomputed with measured unbound fractions, the published figure of 93.7%
reabsorption of filtered PFOA in the male rat becomes **57.6%**; across three
rat studies the admissible range is T = 0.42 to 1.07, i.e. **at most half the
filtered load is returned, possibly none.**

Error 8 is a methods warning rather than a literature error and is the most
likely to recur elsewhere: the loss is undetectable by inspection and
detectable by arithmetic, since a value 1000-fold out of line with its
homologues is a dropped prefix until proven otherwise.

---

## 4. Discussion

### 4.1 Why the binding explanation has persisted

The inference from binding to half-life is mechanistically sound in direction
— bound compound is not filtered — and it fails only on magnitude. Three
features of the literature have obscured that failure.

The binding literature is large, and its internal dispersion (2.38 log units
within compounds) is of the same order as the species differences people have
sought in it; an analyst comparing single values from two species will
frequently find a difference of the size required, and it will be a method
artefact. The measurements that would settle it — one laboratory, several
species, one protocol — are rare, and the only such study available constrains
the rat/human ratio to within twofold. And the term that does carry the
difference is not directly measurable, so it has been inferred by subtraction
and absorbed into model-fitted parameters rather than confronted.

### 4.2 What the decomposition does and does not establish

It establishes that **binding cannot be the explanation**, because the binding
term is bounded at twofold and the effect to be explained is up to
1792-fold. That conclusion is robust to the choice of human dataset (§3.5.4)
and does not depend on the tubular attribution being correct.

It does **not** establish that the tubule is the sole alternative. The fourth
factor is a residual containing the tubule *and* any species difference in
non-renal elimination, and in the human a majority of PFOA elimination is not
renal. What the transporter data (§3.3) contributes is mechanistic support for
reading the residual as predominantly tubular: the reabsorptive machinery
differs qualitatively between the species, in the direction and roughly the
magnitude required.

A sharper test is available and not yet done (§6).

### 4.3 Implications for interspecies extrapolation

**Rodent data should be used for the magnitude of tubular handling, not for
its mechanism.** An assessment whose argument depends on which transporter
reabsorbs a compound cannot be supported by rat data, because the rat uses a
transporter the human does not possess. An assessment that depends only on how
much is reabsorbed can use rat data, with an uncertainty factor reflecting the
22- to 161-fold spread rather than a point estimate.

**"Rodent" is not a useful category here.** Against the human, the rat's
difference is carried by the tubule and the mouse's by the glomerulus. The two
rodents are fast for different reasons, and only one of them is a transporter
story.

**Sex may change the sign, not merely the magnitude.** Female rats net-secrete
PFOA while males net-reabsorb it, driven by a 114-fold difference in one
apical transporter. No comparable measurement exists in humans, and the
possibility that human subpopulations differ in the same qualitative way has
not been excluded.

### 4.4 Implications for regulatory practice

**A single unbound fraction for all compounds is not defensible.** The
measured range across the series is a factor of 780 (0.0004 for PFHpA to 0.31
for PFPA). Replacing a constant with a measured table requires no new data.

**Reabsorption percentages derived from assumed unbound fractions should be
recomputed.** The 8-fold error in the filtered-load denominator propagates
directly into the printed percentage.

**Short-chain compounds require a different model structure, not merely
different parameters.** For C4–C7 the tubular term can exceed unity, so a
model hard-coded for reabsorption-only kinetics will be wrong in sign for
precisely the compounds now entering regulatory attention.

---

## 5. Limitations

**5.1 The tubular factor is not separable.** Equation (3) has two unknowns and
one equation; every "x% reabsorbed" is "x% net reabsorbed". The one in vitro
route that could separate them is closed because OAT4 has no scaling factor.

**5.2 The fourth factor is a residual.** Both clearances in the decomposition
are total rather than renal, so the tubular term absorbs species differences
in biliary and faecal elimination. The identity is exact by construction; the
analysis partitions explanatory burden rather than testing the identity.

**5.3 In vitro transport was measured without plasma protein.** Transporters
respond to the free concentration, and a protein-free assay presents a free
fraction of 1 against an in vivo value near 0.0025 — a 400-fold discrepancy in
the quantity the protein sees. Rankings may survive; magnitudes are in doubt.

**5.4 The decomposition rests on six compounds and two species.** The human
side is a single volunteer, mitigated by a population sensitivity (§3.5.4);
the rat side is multi-source. Mouse and cynomolgus monkey cannot be included
for want of paired data (§6).

**5.5 No uncertainty propagation.** Table 3 reports point estimates; V_d and
CL carry experimental uncertainty that is not propagated into the shares.

**5.6 The binding compilations are dominated by method.** With a median
within-compound range of 2.38 log units, any single value drawn from them
carries an uncertainty larger than most effects reported here.

**5.7 The renal framework addresses at most 42% of human PFOA elimination.**

---

## 6. Conclusions and the experiments that would sharpen them

Plasma protein binding is not the reason humans eliminate perfluoroalkyl acids
more slowly than rats. The binding term is bounded at about twofold by direct
measurement and accounts for a few percent of a difference that reaches three
orders of magnitude; the glomerular term is a species constant and cannot
explain why the difference varies between compounds; and the term that does
carry it corresponds to the one part of the renal system whose machinery is
not shared between the species.

Four measurements would convert the central inference from strong to
decisive, and we state them in advance.

1. **An unbound fraction in cynomolgus plasma.** Fourteen clearances exist in
   this species and no binding measurement of any kind; one dialysis run
   across six compounds would add the bridge species to the decomposition.
2. **Direct renal clearances in place of total clearances** for the compounds
   in Table 3, which would convert the residual into a genuinely tubular term
   and test §4.2 directly.
3. **A transporter assay at physiological albumin**, to establish whether the
   in vitro clearance compilation transfers at all.
4. **An independent human unbound fraction for PFPeA and PFHxA**, the only two
   compounds whose tubular sign is presently in dispute.

A negative result on (1) — an unbound fraction showing the macaque to be a
poor human analogue — would be as informative as a positive one.

---

## Data and code availability

[TO BE COMPLETED — repository DOI or supplementary archive to be minted]

| artefact | path |
|---|---|
| Half-life decomposition (Tables 3–5) | `db/halflife_decomposition.csv`, `scripts/stage13_halflife_decomposition.py` |
| Tubular factors, 52 pairs | `db/tubular_factor_consolidated.csv` |
| Underlying clearance evidence, superseded values retained | `db/tubular_factor_evidence.csv` |
| Binding database, all proteins and species | `db/protein_binding.csv` |
| Transporter direction, abundance, clearance | `db/transporter_direction.csv`, `db/transporter_expression.csv`, `db/transporter_clearance_stage3.csv` |
| Human unbound-fraction recommendation | `db/human_fu_recommended.csv` |
| All figures, one runnable notebook | `notebook/pfas_tk_plots.ipynb` |

---

## References

[TO BE FORMATTED to the target journal's style. Full bibliographic detail for
all 72 sources is held in `db/citations.csv` with DOI and PMID.]

Two sources require disambiguation, both being "Ryu 2024":

**Ryu 2024a** — Ryu S, Burchett W, Zhang S, Modaresi SMS, Agudelo Areiza J,
Kaye E, Fischer FC, Slitt AL. Species-Specific Unbound Fraction Differences in
Highly Bound PFAS: A Comparative Study across Human, Rat, and Mouse Plasma and
Albumin. *Toxics* 2024;12(4):253. doi:10.3390/toxics12040253.

**Ryu 2024b** — Ryu S, Burchett W, Zhang S, Jia X, Modaresi SMS, Agudelo
Areiza J, Rodrigues D, Zhu H, Sunderland EM, Fischer FC, Slitt AL. Unbound
Fractions of PFAS in Human and Rodent Tissues: Rat Liver a Suitable Proxy for
Evaluating Emerging PFAS? *Environ Sci Technol* 2024;58(33):14641–14650.
doi:10.1021/acs.est.4c04050.
