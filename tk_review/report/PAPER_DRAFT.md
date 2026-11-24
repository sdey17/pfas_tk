# Different physiological terms carry the interspecies difference in perfluoroalkyl acid elimination in rat and mouse: a quantitative decomposition of half-life into protein binding, glomerular filtration and tubular handling

**Souvik Dey**

*[Affiliation to be supplied]*

*Draft 2 — internal. Formatted to Toxicological Sciences conventions.*

---

## Abstract

Perfluoroalkyl acid (PFAA) elimination half-lives differ between humans and
rodents by up to three orders of magnitude, and this difference is routinely
attributed to species differences in plasma protein binding. We tested that
attribution quantitatively. Because the half-life of a first-order system is
ln2·V_d/CL and renal clearance is f_u·GFR·T, where T is the fraction of the
filtered load surviving the renal tubule, the interspecies half-life ratio
factorizes exactly into four multiplicative terms whose logarithms partition
the observed difference. Applying this to eight PFAAs with measured volume of
distribution and clearance in each species pair, we find that **different
terms carry the difference in the two rodents.** For the rat, binding
contributes a mean of 8% of the log half-life ratio (bootstrap 95% CI 3–12%)
and the tubular term 65% (53–74%). For the mouse the pattern inverts: binding
contributes 46% (35–58%) and the tubular term is indistinguishable from zero
(−1%; −29 to +23%). The asymmetry has a direct cause: across 14 PFAAs measured
in one laboratory, rat/human unbound-fraction ratios span only 0.76–2.03,
whereas mouse/human ratios reach 18.6. Binding is therefore arithmetically
incapable of explaining the rat, and substantially explains the mouse. For
cynomolgus monkey, where no plasma unbound fraction has ever been measured, we
show that binding alone would require unbound fractions exceeding unity for
three of five compounds, which is impossible. We then identify why the rat's
tubular term is the one that differs: in kidney membrane fractions assayed
together, the secretory transporters OAT1 and OAT3 are present in all three
species, whereas the reabsorptive transporters are not shared at the level of
the genome rather than of detection: OAT4 is absent from the rat and mouse,
and the rodent Oatp1a1 family has no close human ortholog.
Separately, we show that liver fatty acid-binding protein has no affinity
advantage over albumin. That comparison proves to be method-confounded — 63
assays measured albumin only, 25 L-FABP only, and one measured both — so we
repeat it within the single fully paired design, where L-FABP binds 1.34-fold
more weakly (one-sided 95% limit, paired t(13) = -4.42), strengthening the
refutation. The within-method data also reveal a chain-length dependence the
pooled analysis concealed, +0.055 log units per carbon, so the disadvantage
vanishes above C11 — exactly where the hypothesis is invoked. Reading the primary sources
behind a widely used regulatory tabulation yields direct human renal
clearances for nine compounds, giving net tubular reabsorption of 98–99.5%
from C8 to C12, and reveals that the tabulation transcribed one study's total
clearance column as renal clearance, publishing a renal clearance for one
compound whose urinary clearance was below detection. For the three compounds
with direct renal clearances in both species, the renal clearance ratio
decomposes with no residual at all, and binding's share falls further, to 3%
against 68% for a tubular term that is now measured rather than inferred.
Reading one primary source at first hand adds a measured rat renal clearance
for PFOS and a fourth for PFOA, the latter widening the rat disagreement to
10.8-fold and the tubular factor with it to 15.9-fold, while its share of the
log ratio moves only from 58% to 74%. Two further results make the tubular
term partly mechanistic rather than merely residual. Within one study, one
strain and one assay, and with no sex difference in volume of distribution,
the tubular factor changes sign between the sexes for PFOA and PFHxS but not
for PFOS. And because reabsorption and secretion saturate in opposite
directions, the sign of the dose response discriminates between them: every
rat series rises with dose except PFOS, whose four series across two studies
all fall — in the one case where saturable plasma binding, which can only
raise clearance, is furthest from saturating. The practical consequence of
the main result is that no single interspecies scaling approach serves both
rodents; of the second, that none serves every compound.

**Key words:** perfluoroalkyl substances; toxicokinetics; interspecies
extrapolation; renal tubular reabsorption; serum albumin; organic anion
transporters; elimination half-life.

---

## Introduction

Several perfluoroalkyl acids (PFAAs) have human elimination half-lives
measured in years, while the same compounds are cleared by rodents within days
(Chang et al., 2012; Abraham et al., 2024). This gap is the central
quantitative problem in using rodent toxicology to derive human reference
values, and the size of the correction applied depends on what is believed to
cause it.

The most frequently invoked cause is plasma protein binding. PFAAs are
extensively bound to serum albumin, only the unbound fraction is available for
glomerular filtration, and binding varies systematically with perfluorinated
chain length (Forsthuber et al., 2020; Smeltz et al., 2023). The inference
that species differ in half-life because they differ in binding is
mechanistically plausible and practically attractive, since binding is
measurable in vitro and would permit an interspecies scaling factor to be
derived without animal experiments.

This inference is testable, because the half-life decomposes exactly into
separately measurable physiological terms. We carry out that decomposition and
reach three conclusions that differ from the prevailing account.

First, **the attribution is correct for the mouse and wrong for the rat.**
Binding accounts for roughly half of the human/mouse difference but under a
tenth of the human/rat difference, which is also the larger difference. A
single binding-based scaling approach cannot serve both.

Second, **the binding term for the rat is bounded, not merely small.** In the
one study measuring human and rodent plasma in a single laboratory, rat/human
unbound-fraction ratios never exceed about twofold in either direction, while
the half-life ratios to be explained reach 1800-fold.

Third, **the term that carries the rat difference corresponds to the one part
of the renal system whose molecular machinery is not conserved.** Measured in
a single assay across three species, secretory transporters are shared and
differ in abundance, while the apical reabsorptive transporters are different
proteins in human and rat.

We additionally report two findings bearing on how this literature should be
used: the dispersion within the albumin binding literature exceeds the
interspecies differences sought within it, and liver fatty acid-binding
protein (L-FABP), the standing explanation for hepatic PFAA accumulation, has
no measurable affinity advantage over the plasma carrier it is supposed to
outcompete.

---

## Materials and Methods

### Theoretical framework

For first-order elimination,

  t½ = ln2 · V_d / CL   (1)

Renal clearance is the product of the filterable fraction, the filtration
rate, and the net effect of the tubule:

  CL_R = f_u · GFR · T   (2)

where T is a dimensionless **tubular factor**, the fraction of the filtered
load that is excreted. Writing FR for fractional reabsorption and S for
secretory clearance,

  T = (1 − FR) + S/(f_u · GFR)   (3)

so T < 1 denotes net reabsorption, T = 1 filtration alone, and T > 1 net
secretion. Equation (3) contains two unknowns and is not separable from in
vivo data; T is therefore carried as a lumped quantity, and all statements of
fractional reabsorption are net.

Combining (1) and (2) for one compound in two species gives an exact
factorization:

  t½,human / t½,animal =
      (V_d,human/V_d,animal) × (f_u,animal/f_u,human)
      × (GFR_animal/GFR_human) × (T_animal/T_human)   (4)

denoted DISTRIBUTION × BINDING × GLOMERULAR × TUBULE. On a logarithmic scale
these terms add, so the share of the observed difference carried by each is
well defined.

Two properties of (4) constrain its interpretation. The glomerular term is a
ratio of species constants and is therefore identical for every compound,
so it cannot account for variation between compounds. And because the
clearances used are total rather than renal, the fourth term is a **residual**
that absorbs any species difference in non-renal elimination as well as
tubular handling; equation (4) is exact by construction, so the analysis
partitions explanatory burden rather than testing the identity.

### Data sources and selection

**Human kinetics.** Volume of distribution and total clearance are from
Abraham et al. (2024): a single male volunteer administered isotopically
labelled PFAAs as single oral doses, with V_d and CL both measured. The
identity ln2·V_d/CL closes at a consistent implied body weight of 82.9 kg
across all 15 compounds, an internal consistency check no other human dataset
permits. As a further check, the half-lives we compute from their V_d and CL
reproduce the half-lives that study reports directly: 50.6 against 51 days
(PFBS), 151.9 against 152 (PFHpA), 1.4 against 1.5 (PFHxA) and 294.9 against
295 (PFDoDA). No alternative human source reports measured V_d and CL for
more than three compounds.

**Animal kinetics.** Rows were restricted to male animals (matching the human
subject), measured rather than assumed values, and the lowest reported dose,
since several compounds show dose-dependent clearance. One study was used per
species–compound and is named with each value (Supplementary Table S1).

**Unbound fractions.** Plasma f_u is from Ryu et al. (2024a) throughout, the
only study measuring human, rat and two mouse strains in a single laboratory;
in a ratio, a systematic method bias divides out. Mouse values are the
geometric mean of CD-1 and C57BL/6.

**Glomerular filtration rates** were taken as 2570, 14 400 and 16 700
mL/kg/day for human, rat and mouse respectively.

**Protein binding compilations.** 288 albumin–water partition coefficients for
63 PFAAs from 20 primary sources and 91 L-FABP coefficients for 32 PFAAs from
10 sources were converted to a common basis, log K in L water per kg protein,
through log K = 4.1774 + log₁₀(n) − log₁₀(K_d/µM), with a bridge constant of
4.8475 for human L-FABP. All 181 rows carrying both a dissociation constant
and a stoichiometry satisfy this identity exactly. Each value was audited
against its primary publication.

**Transporter data.** 305 rows of localization and direction, 40 of
quantitative abundance and 52 intrinsic clearances were compiled; abundance is
from a single three-species assay (Thakur et al., 2024).

### Values excluded, and why

Two kinetic values were excluded, and the reasons form part of the result.

Mouse PFUnDA V_d is tabulated as 33 ± 60 mL/kg against 280 ± 80 for the female
of the same study (EFSA CONTAM Panel, 2020, Appendix C). A standard deviation
twice the mean, and a tenfold sex gap seen for no other compound in that
series, indicate a transcription error; 330 is the obvious reading but was not
substituted, because inferring a digit is not evidence.

Rat PFHxA V_d is 175 mL/kg in that same appendix and 930 mL/kg in our own
extraction from the primary study (Chengelis et al., 2009). The two are
probably V_ss and V_d,β, but neither source states so, and the compound was
dropped from the rat column rather than resolved by assumption.

### Statistical analysis

Binding comparisons report mean, SD, n and range per compound. The paired
albumin-minus-L-FABP difference was tested by paired t-test and Wilcoxon
signed-rank test; because a null result under high dispersion may reflect only
low power, two one-sided tests (TOST) against a ±0.5 log unit equivalence
margin are also reported, together with a one-sided 95% upper confidence limit,
which is the quantity relevant to a directional hypothesis. Confidence
intervals on decomposition shares are nonparametric bootstrap percentile
intervals over compounds (10 000 resamples). Analyses used Python 3.11 with
SciPy; all code and data are available (see Supplementary Material).

---

## Results

### Chain-length dependence and dispersion of albumin binding

Albumin affinity rises with chain length to a maximum at perfluorononanoic
acid (9 total carbons) and declines thereafter. Across the 32 compounds with
three or more independent values the median within-compound range is 2.38 log
units (maximum 4.22), and the pooled within-compound SD across compounds
measured against both proteins is **1.02 log units** (Table 1).

Equilibrium dialysis and differential scanning fluorimetry differ by 2.35 log
units in median, in the direction predicted by the thermal denaturation the
latter measures; these are not two estimates of one quantity.

> **Consequence.** Selecting a single published binding constant is a decision
> with a range of defensible outcomes exceeding two orders of magnitude, and
> the disagreement between laboratories measuring one species exceeds any
> difference reported between species. An interspecies scaling factor for
> PFAA–albumin binding cannot be estimated from this literature. This
> conclusion is independent of, and prior to, the quantitative analysis below.

### Liver fatty acid-binding protein has no affinity advantage over albumin

**Table 1.** Per-compound log K for albumin and L-FABP, with dispersion.
Difference is L-FABP minus albumin; negative values indicate weaker binding by
L-FABP. tC, total carbons.

| compound | tC | albumin n | albumin mean | SD | range | L-FABP n | L-FABP mean | SD | range | difference |
|---|---|---|---|---|---|---|---|---|---|---|
| PFEESA | 4 | 1 | 3.47 | — | 3.5–3.5 | 1 | 1.35 | — | 1.4–1.4 | −2.12 |
| PFBA | 4 | 9 | 1.91 | 0.96 | 0.5–3.6 | 1 | 1.90 | — | 1.9–1.9 | −0.01 |
| PFBS | 4 | 13 | 2.64 | 0.91 | 1.1–4.3 | 4 | 2.40 | 0.47 | 1.8–2.8 | −0.24 |
| PFPeA | 5 | 6 | 2.66 | 0.74 | 1.7–3.6 | 2 | 2.71 | 0.99 | 2.0–3.4 | +0.05 |
| PFHxS | 6 | 16 | 3.63 | 0.82 | 2.3–4.8 | 9 | 3.34 | 0.59 | 2.5–4.6 | −0.28 |
| PFHxA | 6 | 17 | 2.58 | 0.86 | 1.3–3.9 | 4 | 2.44 | 0.32 | 2.1–2.7 | −0.15 |
| PFHpA | 7 | 21 | 3.02 | 0.86 | 1.8–4.2 | 3 | 2.60 | 0.52 | 2.3–3.2 | −0.42 |
| PFOS | 8 | 20 | 3.98 | 0.75 | 2.1–4.8 | 6 | 3.87 | 0.90 | 2.9–5.6 | −0.11 |
| 6:2 Cl-PFESA | 8 | 4 | 4.09 | 1.26 | 2.4–5.0 | 2 | 3.29 | 0.64 | 2.8–3.8 | −0.79 |
| PFOA | 8 | 33 | 3.44 | 0.93 | 0.7–4.6 | 16 | 3.88 | 0.81 | 2.4–5.8 | +0.45 |
| HFPO-TA | 9 | 3 | 6.38 | 1.47 | 4.9–7.9 | 1 | 3.55 | — | 3.5–3.5 | −2.83 |
| PFNA | 9 | 26 | 4.03 | 0.83 | 2.3–5.2 | 13 | 4.00 | 0.60 | 2.8–4.7 | −0.03 |
| 8:2 Cl-PFESA | 10 | 1 | 3.63 | — | 3.6–3.6 | 1 | 4.12 | — | 4.1–4.1 | +0.49 |
| PFDA | 10 | 15 | 3.95 | 0.87 | 2.3–5.0 | 3 | 3.61 | 0.59 | 3.0–4.1 | −0.34 |
| PFUnDA | 11 | 10 | 3.91 | 1.09 | 2.3–5.2 | 3 | 3.65 | 0.46 | 3.1–4.0 | −0.26 |
| PFDoDA | 12 | 6 | 3.54 | 1.25 | 2.2–5.1 | 4 | 3.67 | 0.59 | 2.9–4.3 | +0.13 |
| PFTrDA | 13 | 4 | 3.36 | 1.37 | 1.8–4.5 | 1 | 2.35 | — | 2.4–2.4 | −1.01 |
| PFTeDA | 14 | 4 | 3.17 | 0.93 | 1.9–3.9 | 2 | 3.55 | 0.69 | 3.1–4.0 | +0.38 |
| PFHxDA | 16 | 1 | 3.16 | — | 3.2–3.2 | 1 | 2.79 | — | 2.8–2.8 | −0.37 |

Across all 19 compounds measured against both proteins the mean difference is
−0.39 log units (95% CI −0.79 to +0.01; paired t = −2.06, P = 0.055; Wilcoxon
P = 0.066), and equivalence within ±0.5 log units is **not** established
(TOST P = 0.29). The full set is dominated by compounds resting on a single
measurement on one side (PFEESA, HFPO-TA, PFTrDA), and restricting to the 10
compounds with n ≥ 3 against both proteins gives a mean difference of −0.126
log units (95% CI −0.310 to +0.058; t = −1.55, P = 0.156), for which
equivalence within ±0.5 log units **is** established (TOST P = 0.0006).

> **The test that bears on the hypothesis is one-sided.** The L-FABP
> explanation requires L-FABP to bind PFAAs *more tightly* than the plasma
> carrier. The one-sided 95% upper confidence limit on the difference is
> **+0.023 log units, or 1.06-fold**. The data therefore exclude any
> meaningful affinity advantage, which is a stronger statement than a failure
> to reject.

Two independent observations agree. Deletion of L-FABP changes the mouse liver
partition coefficient by under 20% for 15 of 19 PFAAs. And in whole-tissue
homogenates, binding in liver is within twofold of binding in heart and kidney
for 14 of 16 compounds despite the liver's enrichment in L-FABP (Ryu et al.,
2024b), whose authors draw the same inference.

The pooled within-compound SD is 1.02 log units for albumin and 0.65 for
L-FABP — in both cases larger than the difference between the proteins. The
dispersion, not the central difference, is the operative result.

**That comparison is method-confounded, and the confound had to be tested
rather than assumed.** The analysis above pairs the two proteins *by
compound*: for each compound it averages every albumin measurement and every
L-FABP measurement and differences the means. But those means are drawn from
almost disjoint sets of assays. Of the methods in the compilation, 63
measured albumin only, 25 measured L-FABP only, and exactly one measured
both. Each per-compound difference is therefore a difference between two
different method mixtures, and pairing on compound cannot remove a method
offset, because the pairing is on the wrong factor. Since method is the
principal axis of disagreement in this compilation (Figure 3), that is
precisely the condition under which such a mixture difference biases the
estimate.

Gao et al. provide the test. They measured 14 compounds against human serum
albumin, pooled human serum and human L-FABP in one laboratory by one method,
equilibrium dialysis — a fully paired within-method design in which the
method mixture is identical on both arms by construction.

**Table 2.** Within-method comparison, Gao et al., equilibrium dialysis. K_d
in µM, so a larger value is weaker binding; the difference is in log K
(association), negative where L-FABP binds more weakly. tC, total carbons.

| compound | tC | K_d albumin | K_d serum | K_d L-FABP | K_d ratio | Δlog K |
|---|---|---|---|---|---|---|
| PFBA | 4 | 419 | 482 | 879 | 2.10 | −0.322 |
| PFBS | 4 | 169 | 239 | 436 | 2.58 | −0.412 |
| PFPeA | 5 | 286 | 387 | 685 | 2.40 | −0.379 |
| PFHxA | 6 | 195 | 243 | 443 | 2.27 | −0.356 |
| PFHxS | 6 | 79 | 158 | 210 | 2.66 | −0.425 |
| PFHpA | 7 | 172 | 209 | 358 | 2.08 | −0.318 |
| 6:2 Cl-PFESA | 8 | 67 | 82 | 102 | 1.52 | −0.183 |
| PFOA | 8 | 115 | 166 | 197 | 1.71 | −0.234 |
| PFOS | 8 | 38 | 49 | 81 | 2.13 | −0.329 |
| PFNA | 9 | 73 | 84 | 90 | 1.23 | −0.091 |
| PFDA | 10 | 65 | 70 | 77 | 1.18 | −0.074 |
| PFUnDA | 11 | 79 | 75 | 54 | 0.68 | +0.165 |
| PFDoDA | 12 | 103 | 100 | 95 | 0.92 | +0.035 |
| PFTrDA | 13 | 265 | 284 | 317 | 1.20 | −0.078 |

![Within-method difference in log K between L-FABP and albumin against chain length, from the single fully paired design of Table 2. Negative values indicate weaker binding by L-FABP. The disadvantage is largest for short chains and closes near C12, the region in which the L-FABP hypothesis is invoked.](../figures/paper_fig2.png)

**The claim survives, and strengthens.** The mean within-method difference is
−0.214 log units (95% CI −0.319 to −0.110; paired t(13) = −4.42), sitting
between the two pooled estimates (−0.390 across all 19 compounds, −0.126 for
the n ≥ 3 subset) with the same sign and a far more decisive t. The one-sided
95% upper limit is **−0.128 log units**, which lies *below* zero. The
within-method data therefore do not merely exclude an affinity advantage for
L-FABP; they establish a disadvantage of at least **1.34-fold**. The pooled
headline limit of +0.023 is superseded by a positive finding.

**But the pooled analysis concealed a chain-length dependence, and it
qualifies the refutation.** The within-method difference scales strongly with
chain length: +0.0548 ± 0.0095 log units per total carbon, t(12) = 5.79,
r² = 0.736, surviving adjustment for head group at +0.0510 ± 0.0099. L-FABP
binds 2.13-fold more weakly at ≤ C8 but is indistinguishable from parity at
≥ C11. This bears directly on the hypothesis under test, because the L-FABP
explanation is invoked for *long-chain* hepatic accumulation: the refutation
is strongest exactly where the hypothesis is not applied and weakest where it
is. A claim that L-FABP has no affinity *advantage* is safe across the whole
range; a claim that it binds measurably more *weakly* is not safe above C11.

Two things are not claimed. The apparent sulfonate offset (−0.070 log units,
SE 0.060, t = −1.17) is not significant and adding head group raises r² only
from 0.736 to 0.765. And the top of the range is not monotonic — the
difference peaks at C11 (+0.165), falls at C12 (+0.035) and C13 (−0.078) —
with one compound per carbon there, so the linear fit must not be read as
continuing past C11.

### What binding could explain at most: an arithmetic bound

Before any decomposition, one claim can be made that depends on no residual,
no GFR value and no volume of distribution. Binding's contribution to a
half-life ratio is at most the unbound-fraction ratio itself. Measured ratios
and half-life ratios are therefore directly comparable.

**Table 3.** Binding's maximum possible share of the log half-life
ratio, per compound.

| comparison | n | f_u ratio | half-life ratio | binding ceiling, share of log ratio |
|---|---|---|---|---|
| human/rat | 8 | 1.02–1.62 | 12.8–1794 | **0.3–18.8%** (median 6.8%) |
| human/mouse | 8 | 2.03–18.56 | 3.3–1087 | **24.8–73.6%** (median 53.4%) |

At the largest rat half-life ratio, binding leaves at least 1109-fold
unexplained even if the largest observed binding ratio applied to it. This is
a ceiling rather than an estimate: it assumes the whole unbound-fraction
difference acts on the half-life and that nothing opposes it. **The rat
ceiling is already small, and the rat and mouse ceilings do not overlap.**

We lead with this because it is the most secure form of the central result.
The decomposition that follows adds the question of *which other term* takes
up the difference, and that question, as the Discussion sets out, is not
answered independently of the binding term itself.

### Decomposition of the human/rat half-life ratio

**Table 4.** Human/rat half-life ratio and its four factors, eight PFAAs with
measured V_d and CL in both species.

| compound | tC | t½ human (d) | t½ rat (d) | ratio | distribution | binding | glomerular | tubule |
|---|---|---|---|---|---|---|---|---|
| PFBA | 4 | 4.2 | 0.33 | 13 | 0.63 | 1.62 | 5.60 | 2.2 |
| PFBS | 4 | 50.6 | 0.16 | 322 | 0.73 | 1.48 | 5.60 | 53.0 |
| PFHxS | 6 | 1635 | 28.2 | 58 | 0.45 | 1.09 | 5.60 | 20.8 |
| PFHpA | 7 | 152 | 0.08 | 1794 | 0.65 | 1.02 | 5.60 | 480.5 |
| PFOA | 8 | 2011 | 1.64 | 1228 | 1.08 | 1.22 | 5.60 | 166.8 |
| PFOS | 8 | 1211 | 28.7 | 42 | 0.40 | 1.45 | 5.60 | 13.1 |
| PFNA | 9 | 1304 | 28.8 | 45 | 0.43 | 1.12 | 5.60 | 16.6 |
| PFDA | 10 | 862 | 46.4 | 19 | 0.53 | 1.62 | 5.60 | 3.9 |

![Share of the log half-life ratio carried by the binding and tubular terms, by compound, for the rat (left) and the mouse (right). Negative values indicate a factor acting against the observed difference. The tubular term dominates in the rat; binding dominates in the mouse. Values are those of Table 5.](../figures/paper_fig1.png)

**Table 5.** Share of the log₁₀ half-life ratio carried by each factor (%),
with bootstrap 95% confidence intervals of the mean.

| | distribution | binding | glomerular | **tubule** |
|---|---|---|---|---|
| **human/rat** mean (n = 8) | −14 | **8** | 42 | **65** |
| 95% CI | −20 to −8 | **3 to 12** | 32 to 53 | **53 to 74** |
| **human/mouse** mean (n = 8) | −12 | **46** | 67 | **−1** |
| 95% CI | −22 to −3 | **35 to 58** | 44 to 98 | **−29 to +23** |

For the rat, the tubular term carries 65% of the log half-life ratio and
binding 8%. Distribution works *against* the observed difference for seven of
eight compounds: on volume of distribution alone the human would eliminate
these compounds faster than the rat.

### The mouse inverts the pattern

For the mouse, binding carries 46% of the log ratio and the tubular term is
indistinguishable from zero (Table 5). The cause is direct: mouse/human
unbound-fraction ratios range from 2.0 (PFDoDA) to 18.6 (PFHxS), whereas
rat/human ratios across the same 14-compound panel span only 0.76–2.03
(Table 6). Both mouse strains agree closely except for PFOS, where CD-1 and
C57BL/6 differ fivefold.

**Table 6.** Plasma unbound fractions, single laboratory (Ryu et al., 2024a).

| compound | human | rat | mouse CD-1 | mouse C57BL/6 | mouse/human | rat/human |
|---|---|---|---|---|---|---|
| PFBA | 0.22900 | 0.37000 | 0.75700 | 0.76700 | 3.3 | 1.62 |
| PFHpA | 0.00445 | 0.00454 | 0.03310 | 0.03230 | 7.3 | 1.02 |
| PFOA | 0.00245 | 0.00298 | 0.01520 | 0.01420 | 6.0 | 1.22 |
| PFNA | 0.00187 | 0.00210 | 0.00455 | 0.00438 | 2.4 | 1.12 |
| PFDA | 0.00062 | 0.00100 | 0.00150 | 0.00127 | 2.2 | 1.62 |
| PFDoDA | 0.00011 | 0.00008 | 0.00022 | 0.00020 | 2.0 | 0.76 |
| PFHxS | 0.00069 | 0.00076 | 0.01280 | 0.01300 | 18.6 | 1.09 |
| PFOS | 0.00075 | 0.00109 | 0.01280 | 0.00256 | 7.6 | 1.45 |

> **The rat binding term is bounded, not merely small.** Across all 14 PFAAs
> in this panel the rat/human ratio spans 0.76–2.03. A term bounded at twofold
> cannot produce an 1800-fold difference in half-life, and no refinement of
> the binding data can alter that, because the bound is on the quantity rather
> than on knowledge of it. The same argument does not apply to the mouse,
> where the term reaches 18.6-fold and is correspondingly consequential.

### Cynomolgus monkey: the required binding is physically impossible

No plasma unbound fraction has ever been measured in cynomolgus plasma for any
PFAA, and our compilation of 1497 binding measurements contains no macaque
value of any kind. Rather than omit the species, we ask what unbound fraction
would be required for binding alone to account for the observed human/monkey
half-life difference (Table 7).

**Table 7.** Unbound fraction the cynomolgus monkey would require for binding
alone to explain its half-life difference from the human.

| compound | t½ human (d) | t½ monkey (d) | ratio | required f_u | multiple of human f_u |
|---|---|---|---|---|---|
| PFBA | 4.2 | 0.22 | 19 | **17.7** | 77× |
| PFHxA | 1.4 | 0.23 | 6 | **1.96** | 51× |
| PFHxS | 1635 | 153 | 11 | 0.0170 | 25× |
| PFOA | 2011 | 10.1 | 199 | 0.729 | 297× |
| PFOS | 1211 | 127 | 10 | 0.0095 | 13× |

Three of five required values exceed unity and are therefore impossible, and
all five lie far outside the range within which interspecies f_u ratios are
observed to fall. The conclusion that binding cannot explain the monkey is
reachable without the missing measurement.

### Direct human renal clearances, and a transcription error in a regulatory tabulation

The analysis above uses total clearances, because direct human renal
clearances were thought unavailable. Reading the three primary sources behind
EFSA's tabulation (EFSA CONTAM Panel, 2020) shows that they are available, and
that the tabulation misrepresents one of them.

**Table 8.** Human renal clearance from three primary studies (mL/kg/day).

| compound | Harada 2005 | Zhang 2013 mean | Zhang 2013 median | Fujii 2015 |
|---|---|---|---|---|
| PFHpA | — | 0.61 (n=31) | 0.41 | 0.674 ± 0.374 (n=10) |
| PFOA | **0.030 ± 0.013** (n=20) | 0.77 (n=66) | 0.19 | **0.044 ± 0.01** (n=10) |
| PFNA | — | 0.15 (n=50) | 0.094 | 0.038 ± 0.01 |
| PFDA | — | 0.096 (n=60) | 0.035 | 0.015 ± 0.01 |
| PFUnDA | — | 0.065 (n=63) | 0.042 | 0.005 ± 0.00 |
| PFDoDA | — | — | — | 0.005 ± 0.00 |
| PFTriDA | — | — | — | 0.006 ± 0.01 |
| PFTeDA | — | — | — | below detection |
| PFHxS | — | 0.027 (n=64) | 0.015 | — |
| PFOS | **0.015 ± 0.010** (n=20) | 0.037 (n=66) | 0.024 | — |

For PFOA, Harada (0.030) and Fujii (0.044) agree within 1.5-fold while Zhang's
mean is 18- to 26-fold higher. Zhang's own distribution is strongly
right-skewed (mean 0.77, median 0.19) and that study computed clearance from
whole blood converted to serum by a haematocrit factor, where the other two
measured serum directly. We treat Zhang as the outlier for PFOA rather than
averaging it in.

**A ninth literature error.** EFSA's per-compound renal clearance table
reproduces Zhang exactly — every value matches Zhang's Table 2, the apparent
"range" being the young-female against male/older-female groups — but for
Fujii it transcribes the **total** clearance column into a renal clearance
table:

| compound | EFSA "renal" | Fujii urinary | Fujii total |
|---|---|---|---|
| PFUnDA | 0.065 | 0.005 | 0.065 |
| PFDoDA | 0.070 | 0.005 | 0.070 |
| PFTriDA | 0.077 | 0.006 | 0.077 |
| PFTeDA | 0.224 | **below detection** | 0.224 |

For PFTeDA a renal clearance of 0.224 mL/kg/day is published for a urinary
clearance the authors could not detect. The figure of 0.005 that appears
elsewhere in the same EFSA document, and which made it appear
self-contradictory, is Fujii's correct urinary value: the summary sentence was
right and the table was wrong.

**Direct tubular factors.** Substituting measured renal for total clearance
converts the human upper bounds into point estimates (Table 9). Every direct
value falls below its former bound, as it must.

**Table 9.** Human tubular factors from measured renal clearance.

| compound | source | CL_R | f_u | filtered load | **T** | net reabsorbed | former bound |
|---|---|---|---|---|---|---|---|
| PFHpA | Fujii 2015 | 0.674 | 0.00445 | 11.44 | 0.0589 | 94.1% | T ≤ 0.051 |
| PFOA | Harada 2005 | 0.030 | 0.00245 | 6.30 | 0.0048 | 99.5% | T ≤ 0.019 |
| PFOA | Fujii 2015 | 0.044 | 0.00245 | 6.30 | 0.0070 | 99.3% | T ≤ 0.019 |
| PFNA | Fujii 2015 | 0.038 | 0.00187 | 4.81 | 0.0079 | 99.2% | T ≤ 0.032 |
| PFDA | Fujii 2015 | 0.015 | 0.00062 | 1.59 | 0.0094 | 99.1% | T ≤ 0.094 |
| PFUnDA | Fujii 2015 | 0.005 | 0.00042 | 1.08 | 0.0046 | 99.5% | T ≤ 0.266 |
| PFDoDA | Fujii 2015 | 0.005 | 0.00011 | 0.27 | 0.0185 | 98.2% | T ≤ 3.08 |
| PFOS | Harada 2005 | 0.015 | 0.00075 | 1.94 | 0.0078 | 99.2% | T ≤ 0.046 |
| PFHxS | Zhang 2013 | 0.027 | 0.00069 | 1.79 | 0.0151 | 98.5% | T ≤ 0.036 |

Human net tubular reabsorption is therefore 98–99.5% from C8 to C12 and 94%
for C7, tighter than the bounds permitted and systematically closer to
complete.

**Renal elimination is a minority route in humans and a majority route in
rodents.** Harada's Table 5 gives CL_R/CL_tot as 25.1% for human male PFOA
against 91.4% for male Wistar rat and 40.4% for male Japanese macaque. Fujii
gives human urinary 0.044 against a total of 0.096 (46%), with a biliary
clearance of 2.62 mL/kg/day — sixty times the urinary — of which an estimated
98% is reabsorbed. In the mouse, Fujii's urinary PFOA is 13.1 of a total 14.2
(92%). This is the measured basis for treating the human tubular term in
equation (4) as a residual, and it quantifies how large the non-renal
component of that residual is.

### A renal-only decomposition, with no residual term

The decomposition of Tables 2–4 uses total clearance, so its fourth factor is
a residual absorbing both tubular handling and any species difference in
non-renal elimination. With direct renal clearances now available in both
species, the renal clearance ratio decomposes exactly into three measured
terms and no residual:

  CL_R,rat / CL_R,human =
      (f_u,rat/f_u,human) × (GFR_rat/GFR_human) × (T_rat/T_human)   (5)

Here T is genuinely the tubular factor. Three compounds have a direct renal
clearance in male human and male rat (Table 10). PFNA has one in the human and
only in the female rat, and pairing across sexes is inadmissible in a species
where sex changes the sign of the tubular term.

**Table 10.** Renal-only decomposition, male human against male rat. Multiple
rows where sources disagree.

| compound | CL_R human (source) | CL_R rat (source) | ratio | T human | T rat | binding | glomerular | **tubule** |
|---|---|---|---|---|---|---|---|---|
| PFHpA | 0.674 (Fujii) | 983.1 (Ohmori) | 1459 | 0.0589 | 15.04 | 1.02 | 5.60 | **255.2** |
| **PFOA** | 0.030 (Harada) | **4.26 (Kim)** | 142 | 0.0048 | 0.0993 | 1.22 | 5.60 | **20.8** |
| **PFOA** | 0.044 (Fujii) | **4.26 (Kim)** | 97 | 0.0070 | 0.0993 | 1.22 | 5.60 | **14.2** |
| PFOA | 0.030 (Harada) | 18.2 (Kemper) | 607 | 0.0048 | 0.424 | 1.22 | 5.60 | 89.0 |
| PFOA | 0.030 (Harada) | 21.2 (Ohmori) | 707 | 0.0048 | 0.494 | 1.22 | 5.60 | **103.7** |
| PFOA | 0.030 (Harada) | 46.08 (Kudo) | 1536 | 0.0048 | 1.074 | 1.22 | 5.60 | **225.4** |
| PFOA | 0.044 (Fujii) | 18.2 (Kemper) | 414 | 0.0070 | 0.424 | 1.22 | 5.60 | **60.7** |
| PFOA | 0.044 (Fujii) | 21.2 (Ohmori) | 482 | 0.0070 | 0.494 | 1.22 | 5.60 | **70.7** |
| PFOA | 0.044 (Fujii) | 46.08 (Kudo) | 1047 | 0.0070 | 1.074 | 1.22 | 5.60 | **153.7** |
| **PFHxS** | 0.027 (Zhang) | **0.76 (Kim)** | 28 | 0.0151 | 0.0694 | 1.09 | 5.60 | **4.6** |
| **PFOS** | 0.012 (Harada) | **1.92 (Kim)** | 160 | 0.0062 | 0.1223 | 1.45 | 5.60 | **19.7** |
| **PFOS** | 0.037 (Zhang mean) | **1.92 (Kim)** | 52 | 0.0191 | 0.1223 | 1.45 | 5.60 | **6.4** |

Rows in bold are the primary analysis; the remainder are retained as a
sensitivity range. **PFOS enters as a measured value rather than a bound**,
and **a fourth rat PFOA renal clearance has been added**, both from Kim et
al. (2016) read at source.

Kim's values are weighted as primary because they are the only rat renal
clearances in this set taken from a primary table by the present authors,
with the clearance definition stated in that paper's Methods and four
internal consistency checks passing. Each alternative is weaker in a specific
way: Kemper (2003) is an unpublished study reached only through a regulatory
tabulation; the Ohmori value is a bar digitised from a figure whose two
independent digitisations disagree 2.1-fold; and the Kudo value reaches us
through Harada's Table 5 rather than its own.

For PFOA the four rat sources now span 10.8-fold, and the tubular factor with
them, from 14.2 to 225.4 — a 15.9-fold range. **Its share of the log ratio
moves only from 58% to 74%.** The magnitude of the tubular factor is badly
determined; the conclusion that it dominates is not. On the primary values
alone the tubular share runs 46% (PFHxS) to 76% (PFHpA).

**Table 11.** Share of the log renal clearance ratio, mean over the
stated source combinations.

| | n | binding | glomerular | **tubule** |
|---|---|---|---|---|
| renal-only, primary values (Kim-weighted) | 8 | **5.5%** | 37.5% | **57.0%** |
| renal-only, all sources | 14 | **4.4%** | 32.7% | **62.9%** |
| total clearance (equation 4, residual) | 8 | 8% | 42% | 65% |

**Removing the non-renal routes does not rescue binding.** Its share falls
from 8% to 4–6%, and the tubular term is confirmed at 57–63% as a measured
quantity rather than a residual. The glomerular and tubular entries should be
read with the GFR ambiguity below in mind; the binding entry is unaffected by
it. This is the strongest available form of the
central result, and it retires the principal limitation of the
total-clearance analysis for these three compounds.

PFHxS is the informative exception: its tubular factor differs only 4.6-fold
between the species and the glomerular term (52% of the log ratio) exceeds
the tubular one (46%). Not every compound is a tubular story.

**How much of the half-life this accounts for.** The renal route is 100% of
human PFHpA elimination, 25–46% of human PFOA and about 51% of human PFHxS,
against 61%, 42–91% and 8% respectively in the rat. Equation (5) therefore
accounts for all of the PFHpA difference, roughly half of the PFOA
difference, and only part of the PFHxS difference; the remainder is a
difference in route partition rather than in renal handling.

**PFDA's renal clearance is not recoverable from Ohmori's figure.** Figure 2
was re-examined at 300 dpi: the PFDA bars for both sexes sit flat on the
baseline carrying only their significance letters, yielding only
CL_R < 15 mL/kg/day and hence T < 1.04, a bound that straddles unity and
settles nothing. That paper's own Figure 3 caption states that PFNA-male,
PFDA-male and PFDA-female "are not distinguished on the graph", and its
published regression (CL_R = 0.48·CL_tot + 15.63, r² = 0.981) cannot recover
them: at PFDA's total clearance of 5.2 mL/kg/day it predicts a renal
clearance of 18.1, exceeding the total by 3.5-fold. The intercept is fitted
across totals reaching 3000 and swamps small values.

### Four further compounds, and why only two can be bounded

A targeted search was made for directly measured rat renal clearances of
PFBA, PFBS, PFDA and PFOS, the four compounds with the best remaining
prospect of entering equation (5). **None exists.** The reasons differ, and
the pattern of reasons is itself informative. PFBA (Chang 2008) and PFBS
(Olsen 2009; Chengelis 2009) were each studied with urine collected — the
methods of all three papers list urine among the matrices analysed, and
Olsen reports that "urine appeared to be a major route of elimination" — yet
none publishes a renal clearance. For PFOS no rat urine study was located at
all: the National Toxicology Program study that supplies the most recent rat
PFOS kinetics (Huang 2019) collected plasma, liver, kidney and brain only,
so no renal clearance can be derived from it at any dose. Harada (2005)
reports in its abstract that renal clearance sex differences "have been
reported in rats and Japanese macaques", but its Table 5 carries rat rows for
PFOA only; PFOS there is given for cynomolgus monkey and human. That lead
resolves to PFOA.

**The human side is the tighter constraint.** No human renal clearance exists
for PFBA or PFBS, so neither compound could enter equation (5) however
complete the rat data became. PFDA and PFOS do have human values (Table 8),
so the two compounds whose rodent measurements are closest to retrievable are
precisely the two that cannot be used, and the two that can be used are
limited by the rat. The practical consequence is that equation (5) is
extended not by new measurements but by one-sided bounds, since a bound on
CL_R propagates to a bound on every term.

Two bounds are available, of very different quality. For PFDA, Kudo (2001)
reports 0.2% of an intraperitoneal dose recovered in urine over 120 h in male
Wistar rats. Since A_urine = CL_R·AUC(0,τ) and Dose = CL_tot·AUC(0,∞),

  CL_R = (A_urine/Dose) · CL_tot / [AUC(0,τ)/AUC(0,∞)]   (6)

and the bracket evaluated on a one-compartment curve is 1 − exp(−ln2·τ/t½) =
0.083, so the 120 h window captured only 8.3% of the total exposure and the
raw 0.2% badly understates the renal share. For any multi-exponential curve
with the same *terminal* half-life the bracket is larger — it is a weighted
average of (1 − exp(−k·τ)) over rate constants all exceeding the terminal one
— so dividing by the one-compartment value overstates CL_R and equation (6)
returns an upper bound. For PFOS the bound is the identity CL_R ≤ CL_tot,
which is loose because the rat excretes PFOS substantially in bile.

**The PFDA route is cross-validated rather than assumed.** Four compounds
carry both a Kudo urinary fraction and an Ohmori renal clearance, so
equation (6) can be run where the answer is already known (Table 12). It
reads high in all four cases, by 1.05- to 1.80-fold. A one-directional
overstatement is what an upper bound tolerates, and it simultaneously denies
the route any claim to be a point estimate. For PFHpA-male and PFOA-female
the half-lives are 0.10 and 0.08 d, so τ = 5 d captures the entire exposure
and the correction is exactly unity; the 1.50- and 1.80-fold gaps there are
straight inter-study disagreements between Kudo's urinary recovery and
Ohmori's digitised bars — different doses and routes (Kudo intraperitoneal,
Ohmori intravenous) in the same strain and laboratory. That sets the floor on
how far any derived renal clearance can be trusted in this literature.

**Table 12.** Cross-validation of equation (6) against Ohmori's measured
renal clearances. The two rows with AUC fraction 1.000 isolate inter-study
disagreement from model error.

| sex | compound | urinary fraction | CL_tot | t½ (d) | AUC fraction | predicted CL_R | Ohmori CL_R | ratio |
|---|---|---|---|---|---|---|---|---|
| female | PFNA | 0.510 | 105.7 | 2.44 | 0.758 | 71.1 | 67.5 | 1.05 |
| male | PFOA | 0.550 | 50.5 | 5.63 | 0.460 | 60.4 | 45.0 | 1.34 |
| male | PFHpA | 0.920 | 1603.9 | 0.10 | 1.000 | 1475.6 | 983.1 | 1.50 |
| female | PFOA | 0.800 | 2233.2 | 0.08 | 1.000 | 1786.6 | 994.4 | 1.80 |
| male | PFNA | 0.020 | 6.9 | 29.57 | 0.111 | 1.25 | <15 | — |
| male | PFDA | 0.002 | 5.2 | 39.92 | 0.083 | 0.125 | <15 | — |

For male rat PFDA, equation (6) gives CL_R ≤ 0.125 mL/kg/day and T ≤ 0.0087,
a bound 120-fold tighter than the figure's T < 1.04 and, unlike it,
decisively below unity.

**Table 13.** Equation (5) extended as one-sided bounds. All rat inputs are
upper bounds, so the ratio and the tubular factor are upper bounds; the
binding and glomerular terms are exact. PFOS rows are summarised as ranges
over four human and three rat sources (12 combinations). Tubular shares of the log ratio are
−4% for the first PFDA row, 62–70% across the PFOS combinations, and
suppressed for the remaining two PFDA rows for the reason given below.

| compound | CL_R human | CL_R rat | ratio | T human | T rat | binding | glom. | **tubule** |
|---|---|---|---|---|---|---|---|---|
| PFDA | 0.015 (Fujii) | ≤ 0.125 | ≤ 8.3 | 0.0094 | ≤ 0.0087 | 1.62 | 5.60 | **≤ 0.92** |
| PFDA | 0.035 (Zhang med.) | ≤ 0.125 | ≤ 3.6 | 0.0220 | ≤ 0.0087 | 1.62 | 5.60 | **≤ 0.39** |
| PFDA | 0.096 (Zhang mean) | ≤ 0.125 | ≤ 1.3 | 0.0604 | ≤ 0.0087 | 1.62 | 5.60 | **≤ 0.14** |
| PFOS | 0.012–0.037 | ≤ 9.2–13.1 | ≤ 250–1092 | 0.0062–0.0191 | ≤ 0.59–0.83 | 1.45 | 5.60 | **≤ 31–135** |

Shares are suppressed where any single share exceeds 100% of the whole, the
symptom of the three logarithmic terms nearly cancelling: the percentages
diverge while the underlying factors remain exactly determined. This is the
same pathology as the difference form discussed in Methods, and in those rows
the factors, not the shares, are the interpretable quantities.

**PFDA breaks the pattern of the three compounds in equation (5).** Its
tubular factor is bounded below unity against all three human sources, so the
rat tubule is *not* more secretory than the human's for PFDA, and the renal
clearance difference is carried by glomerular filtration (5.60-fold, a species
constant) and binding (1.62-fold) alone. PFOS by contrast keeps the familiar
pattern, its tubular factor bounded above by 135 and exceeding 31 on every
one of the twelve source combinations while carrying 62–70% of the log ratio.
The tubular explanation is therefore compound-specific even within the
carboxylates: it is strong for PFHpA and PFOA, weak for PFHxS, and for PFDA
the available bound excludes it.

### What the glomerular term depends on, and what it does not

The glomerular factor is a ratio of species constants, so it is identical for
every compound and cannot explain between-compound variation. It can,
however, be wrong, and the project holds **two rat GFR values differing
1.82-fold** — 14400 and 7920 mL/kg/day — with a third near 12200 implied by
Harada et al. (2005). The tables above, like earlier work in this area, use
the larger value.

GFR enters equation (5) twice and in opposite directions: it multiplies the
glomerular term and divides the tubular one. The choice therefore moves
log10(1.82) = 0.26 log units from one to the other and changes nothing else.

**Table 14.** Mean share of the log renal clearance ratio under each
available rat GFR, over all compound and source combinations.

| rat GFR (mL/kg/day) | source | binding | glomerular | **tubule** |
|---|---|---|---|---|
| 7920 | Ohmori/Shim | **4.4%** | 21.4% | 74.2% |
| 12200 | implied by Harada 2005 | **4.4%** | 29.6% | 66.0% |
| 14400 | OEHHA/Han | **4.4%** | 32.7% | 62.9% |

**The binding share is invariant.** The central claim is therefore insensitive
to a GFR ambiguity that moves the glomerular share by half and the tubular
share by a sixth. The corollary is that the glomerular and tubular shares
must be reported as ranges rather than point values: glomerular 21–33%,
tubular 63–74%.

### Reabsorption or secretion? The sign of the dose response

Equation (3) has two unknowns and one equation, so a measured T cannot by
itself say whether net reabsorption or net secretion dominates. This has been
the principal limitation of the approach. It can be partly removed without
separating the two terms numerically, because both processes are
carrier-mediated and therefore saturable, and they saturate in opposite
directions:

- saturable **reabsorption** → at high dose less is reclaimed → T rises → clearance rises
- saturable **secretion** → at high dose less is pumped out → T falls → clearance falls

The *sign* of the dose response therefore discriminates. Dose series exist for
28 species–sex–compound combinations across five studies.

**One confounder matters, and it is one-directional.** Saturable plasma
protein binding raises f_u at high dose, and CL_R = f_u·GFR·T rises even with
T unchanged. But binding saturation can only push clearance *up*; it cannot
produce a fall. A rise is therefore consistent with either saturable
reabsorption or saturable binding, while a fall excludes reabsorption
dominance outright. The fall is the diagnostic observation; a rise is
diagnostic only where binding saturation can be excluded independently.

**Table 15.** Dose response of clearance. Direction assigned at ratio
thresholds of 1.15 and 0.87. The molar ratio is PFAS C_max at the high dose
against 0.60 mM serum albumin.

| compound | species | n series | rises | falls | PFAS:albumin at high dose |
|---|---|---|---|---|---|
| PFBS | rat | 2 | 2 | 0 | 0.70–1.25 |
| PFHxS | rat | 2 | 2 | 0 | 0.50–0.68 |
| PFOA | rat | 2 | 1 | 0 | not reported |
| **PFOS** | **rat** | **4** | **0** | **4** | **0.35–0.45** |
| PFHxS | mouse | 2 | 2 | 0 | not reported |
| PFOS | mouse | 2 | 1 | 0 | not reported |
| C7–C14 carboxylates | mouse | 14 | 8 | 3 | not reported |

![Ratio of clearance at the high dose to clearance at the low dose, for every dose series, on a logarithmic scale. Points right of unity indicate clearance rising with dose, the signature of saturating reabsorption; points left of unity indicate a fall, which saturable plasma binding cannot produce. Every rat series lies right of unity except the four PFOS series.](../figures/paper_fig3.png)

**All four rat PFOS series fall**, across two independent studies (ratios
0.24, 0.43, 0.66, 0.82). No other rat series falls. The confounder is
excluded exactly where it matters: PFOS's PFAS-to-albumin molar ratio at the
high dose is 0.35–0.45, the **lowest** of the six series for which serum
concentrations are reported, while the two compounds whose rises are
confounded sit at 0.50–1.25. Binding was furthest from saturation precisely
where the fall was observed.

**In the mouse the dose effect scales with chain length**: log10 of the
clearance ratio rises by 0.133 ± 0.027 per total carbon, t(12) = 4.89,
r² = 0.666, crossing unity at C9.2. Saturable reabsorption predicts this if
the reabsorptive carrier binds longer chains more tightly — the longer the
chain, the more reabsorption there is to saturate.

**T and the dose sign answer different questions, and the pair is more
informative than either.** T says which process *dominates* at the measured
dose; the sign says which *saturates first* as dose rises. PFOS is
net-reabsorbing (T = 0.122) yet its clearance falls, placing a smaller
secretory or non-renal component as the first to saturate. That is a coherent
picture, not a contradiction.

**PFOS separates from every other compound on three independent axes**: the
falling dose response, the absence of a sex difference in renal clearance,
and a tubular factor below unity in both sexes where PFOA's and PFHxS's cross
unity. Three unrelated kinds of evidence agree that PFOS is not handled like
the rest — a result in its own right, and a warning against treating the
chemical class uniformly.

The alternative is not dismissed. PFOS undergoes substantial enterohepatic
circulation, so a fall in *total* clearance may reflect saturable biliary
transport rather than saturable tubular secretion. Separating those requires
a dose series in renal clearance specifically, which does not exist for PFOS.
We also note that the rat tubular factor for PFHpA is 15.0, far above unity
and the only net-secreting value in the set, while mouse PFHpA clearance
falls with dose; those agree, but they are different species, and we do not
treat a cross-species agreement as corroboration.

### Renal transporter expression: the reabsorptive machinery is not shared

Of 305 directional rows, 299 resolve to a reabsorption-or-secretion
assignment across 21 transporters with zero genuine conflicts: papers that
disagree about which compounds a transporter carries do not disagree about
where it sits.

**Table 16.** Renal transporter abundance, kidney membrane fractions, single
assay (Thakur et al., 2024). Values are pmol/mg membrane protein, mean ± SD.

| transporter | membrane | net process | human (n=15) | rat (n=10) | mouse (n=14) |
|---|---|---|---|---|---|
| OAT1 | basolateral | secretory | 4.8 ± 1.2 | 32.8 ± 13.7 | 19.1 ± 7.1 |
| OAT3 | basolateral | secretory | 4.3 ± 1.1 | 6.3 ± 1.1 | 6.5 ± 3.6 |
| OAT2 | basolateral | secretory | 0.3 ± 0.1 | 3.1 ± 3.0 | 2.3 ± 0.6 |
| **OAT4** | apical | **reabsorptive** | **0.5 ± 0.2** | **absent** | **absent** |
| **Oatp1a1** | apical | **reabsorptive** | **no ortholog** | **1.9 ± 2.1** | 3.1 ± 2.8 |
| Oatp1a3 | apical | reabsorptive | no ortholog | 0.7 ± 0.6 | absent |
| URAT1 | apical | urate exchange | 3.1 ± 1.2 | 2.2 ± 0.7 | 10.8 ± 1.1 |

![Renal transporter abundance in kidney membrane fractions of human, rat and mouse (Thakur et al., 2024). The secretory transporters left of the divider are present in all three species and differ in degree. Of the reabsorptive transporters right of it, OAT4 is absent from the rodent genomes and the rodent Oatp1a family has no close human ortholog.](../figures/paper_fig4.png)

**The reabsorptive entries are genomic, not analytical.** OAT4 is absent from
the rat and mouse genomes, and the rodent Oatp1a1/Oatp1a3 family has no close
human ortholog, being only ~40% similar to human OATP1B1 (Thakur et al.,
2024, Figure 1). These are not proteins that fell below a detection
threshold; the genes are not present. The claim is therefore one of orthology
rather than of abundance, and is correspondingly firmer.

The secretory side is conserved in kind and differs in degree; the
reabsorptive side is not the same machinery at all. Because reabsorption
determines the long half-lives, the half of the renal system governing
persistence is the half that is not shared.

**Two transporters independently predict the sex reversal reported above.**
In the rat, the apical reabsorptive Oatp1a1 is 114-fold more abundant in
males, while the basolateral Oat2, which feeds secretion, is 20-fold more
abundant in females; Thakur et al. state explicitly that higher female Oat2
drives the higher PFOA excretion of female rats. Both differences act in the
same direction — males reabsorb more, females secrete more — and that is the
molecular counterpart of the sign change in the tubular factor between sexes
obtained from Kim's renal clearances (Table 17). Two laboratories, two
methods, one conclusion.

Two caveats attach to this table. Abundances come from label-free global
proteomics (data-independent acquisition with the total protein approach)
rather than targeted assays with stable-isotope standards, so the absolute
values carry more uncertainty than the quoted standard deviations suggest,
and Thakur et al. flag proteins below 0.5 pmol/mg as needing confirmation in
enriched fractions — human OAT4 sits at exactly 0.5. The absence of rodent
OAT4 does not depend on this, being a genomic fact rather than a measurement.

Expression also gates relevance. Applying an expression filter removes the
three highest intrinsic clearances in the compilation, all BCRP, including the
single highest value anywhere in it (279 µL/min/mg for PFOA); human renal BCRP
is not detected by either laboratory that has examined it. The same
distinction reconciles conflicting reports for ASBT and OSTα/β, transport by
which was observed in transfected cells but which are absent from human and
rat kidney: a transfected-cell result reports what a protein can carry, not
whether it is present to do so.

### Tubular factors and sex

Resolving the clearance, binding and filtration data onto equation (2) gives T
for 52 species–compound pairs. In the human, PFOA through PFUnDA exceed 98%
net reabsorbed on every available binding source. Net secretion is not a
corner case: mouse PFHxA secretes at T = 2.23, and the female rat secretes
PFOA at T = 24.6 against a male value of 0.42–1.07 — a change of sign within
one species, concurrent with a 114-fold sex difference in abundance of the
reabsorptive transporter Oatp1a1.

**One study settles this without any cross-study assembly.** Kim et al.
(2016) measured renal clearance in both sexes of one rat strain by one assay
and report no significant sex difference in volume of distribution for any of
the three compounds studied. The sex difference therefore cannot be a
distribution effect; it sits in the tubular term, and that term changes sign
(Table 17).

**Table 17.** Tubular factors by sex from a single study (Kim et al., 2016),
male rat GFR 14400 mL/kg/day.

| compound | CL_R male | CL_R female | T male | T female | female/male |
|---|---|---|---|---|---|
| PFOA | 4.26 | 323.60 | 0.0993 | **7.541** | 76 |
| PFHxS | 0.76 | 63.87 | 0.0694 | **5.836** | 84 |
| PFOS | 1.92 | 1.72 | 0.1223 | 0.1096 | 0.9 |

For PFOA and PFHxS the same compound is net-reabsorbed in the male and
net-secreted in the female within one experiment. PFOS, which shows no sex
difference in renal clearance, is the internal control establishing that the
contrast is not an artefact of the assay.

### Errors identified in primary and secondary sources

Applying the rule that a primary publication supersedes any compilation of it
identified nine discrepancies, three of which propagate into published
estimates.

ATSDR (2021), Table 3-6, misreports rat total clearances transcribed from
Ohmori et al. (2003): PFDA as 207/140 against the paper's 5.2/5.3 mL/kg/day
(40- and 26-fold), and male PFOA as 135 against 50.5 (2.7-fold). OEHHA (2024),
Table A6.4, assumes an unbound fraction of 0.02 for every species and
compound, uncited; measured values are 0.0025 (human) and 0.0030 (rat), so
every reabsorption percentage derived from it is wrong in the denominator by
roughly eightfold. Recomputed with measured values, the published figure of
93.7% reabsorption of filtered PFOA in the male rat becomes 57.6%; across
three rat studies the admissible range is T = 0.42–1.07, i.e. at most half the
filtered load is returned and possibly none.

A ninth, identified above, is EFSA's transcription of Fujii's total
clearance column into a renal clearance table, which publishes a renal
clearance of 0.224 mL/kg/day for a compound whose urinary clearance was below
detection.

The remaining six comprise a units defect in a published K_m/V_max column, a
transposition of two transporters between a paper's text and abstract, an
unresolvable citation for reported unbound fractions, a paper misquoting its
own cited source, one compound reported twice 2.8-fold apart within a single
study, and silent loss of the micro sign during PDF text extraction, detected
twice and worth 1000-fold. The last is a methods hazard rather than a
literature error: it is undetectable by inspection and detectable by
arithmetic, since a value 1000-fold out of line with its homologues is a
dropped prefix until shown otherwise.

---

## Discussion

### Why the binding explanation persisted, and where it holds

The inference from binding to half-life is sound in direction and fails only
on magnitude, and then only for some species. Three features of the literature
obscured this. The albumin literature is large and its internal dispersion
(pooled SD 1.02 log units) is of the same order as the species differences
sought within it, so an analyst comparing single values from two species will
frequently find a difference of the required size that is a method artifact.
Studies measuring several species in one laboratory are rare. And the term
that carries the rat difference is not directly measurable, so it has been
inferred by subtraction or absorbed into fitted model parameters.

Our results do not reject the binding explanation outright; they localize it.
For the mouse it is largely correct. For the rat — the species with the larger
difference, and the species more often used in chronic toxicology — it is
quantitatively excluded.

### The binding and tubular shares are one finding, not two

In equation (4) the tubular term is a **residual**: it is what the half-life
ratio leaves after the distribution, binding and glomerular terms are
removed. An error in the unbound-fraction ratio therefore enters the binding
term with one sign and the tubular term with exactly the opposite sign, and
the two log shares move with correlation −1 by construction.

This is demonstrable. Multiplying every unbound-fraction ratio by a factor k
and recomputing gives, for the human/rat set, mean binding shares of −9.4%,
7.5% and 24.5% at k = 0.5, 1.0 and 2.0, with mean tubular shares of 81.8%,
64.8% and 47.8% — and a **sum that is invariant at 72.3% throughout**. Every
log unit added to binding is subtracted from the residual, exactly.

It follows that "binding is small" and "the tubular term is large" are the
same statement, and reporting them as two findings overstates the evidence.
The result is not weakened by this; it is relocated. The secure claim is the
measured unbound-fraction ratio and the arithmetic ceiling it implies
(*What binding could explain at most*). Equation (5), where T is measured rather than inferred, is the
only place in this work where the two terms are genuinely independent — which
is the reason extending it matters, and the reason the renal-only compounds
carry more weight than their number suggests.

### What the decomposition establishes, and what it does not

It establishes that binding cannot explain the human/rat difference, because
the term is bounded at twofold by direct measurement while the effect reaches
1800-fold. That conclusion does not depend on the tubular attribution being
correct.

It does not establish that the tubule is the sole alternative. The fourth term
is a residual containing tubular handling and any species difference in
non-renal elimination, and in humans a majority of PFOA elimination is not
renal. The transporter data provide mechanistic support for reading the
residual as predominantly tubular — the reabsorptive machinery differs
qualitatively between the species, in the required direction — but a direct
test requires renal rather than total clearances.

### Implications for interspecies extrapolation

**No single scaling approach serves both rodents.** A binding-based correction
is approximately right for the mouse and wrong by orders of magnitude for the
rat. Assessments should state which rodent their correction was derived from.

**Rodent data support the magnitude of tubular handling, not its mechanism.**
An argument depending on which transporter reabsorbs a compound cannot be
supported by rat data, because the rat uses a transporter humans lack.

**Sex may change the sign rather than the magnitude.** Female rats net-secrete
PFOA while males net-reabsorb it. No comparable human measurement exists, and
the possibility that human subpopulations differ qualitatively has not been
excluded.

**A single unbound fraction for all compounds is not defensible.** The measured
range across the series is a factor of 780. Replacing a constant with a
measured table requires no new data.

---

## Limitations

The tubular factor is not separable into reabsorption and secretion at a
given dose; equation (3) has two unknowns and one equation, and the one in
vitro route that could separate them is closed because OAT4 has no scaling
factor to the kidney, while the two available methods for OAT1 and OAT3
disagree 11-fold and in opposite directions. The dose-response sign test
(*Reabsorption or secretion?*, below) partly circumvents this by asking which process saturates
first rather than which dominates, but it rests on *total* clearance for most
compounds, so a falling response cannot be attributed to tubular secretion
rather than to saturable biliary elimination without a renal-clearance dose
series. None exists for PFOS, the compound where the question is sharpest.

The fourth factor is a residual rather than a measured tubular term for the
five compounds lacking paired renal clearances; for PFHpA, PFOA and PFHxS it
has been replaced by the residual-free decomposition of equation (5), which
gives the same answer. Its non-renal content is quantified rather than
assumed: renal elimination is
25–46% of total human PFOA clearance (Harada 2005; Fujii 2015) against 91% in
the male rat and 92% in the mouse, so the residual is substantially non-renal
in the human and substantially tubular in the rodents. Direct human renal
clearances are available for nine compounds (Table 7) but matching male rat
values for only three, which is what limits equation (5) to that set. A
targeted search established that no directly measured rat renal clearance
exists for PFBA, PFBS, PFDA or PFOS, and that for PFBA and PFBS the human
value is missing as well, so those two cannot be added by any rodent
measurement. PFDA and PFOS enter only as the one-sided bounds of Table 12.
The PFDA bound inherits the 1.05- to 1.80-fold spread of Table 11, and the
PFOS bound is loose by construction because it rests on CL_R <= CL_tot in a
species with substantial biliary excretion of PFOS; neither should be read as
an estimate.

**The unbound fractions rest unequally across species, and the rat worst of
all.** Four plasma unbound-fraction sources are available: Ryu et al.
(2024a), two Argoul et al. (2025) experiments, and Fischer et al. (2024).
They cover the human with four, the mouse with three, and **the rat with one**.
The rat binding factor therefore cannot be varied on its numerator at all,
only on the human denominator, so the sensitivity range below understates rat
uncertainty by an unknown amount. The mouse, testable properly, retains high
ratios under all twelve of its source pairings, so the rat–mouse asymmetry is
not an artefact of a single laboratory.

**The binding range quoted in the abstract and in *What binding could explain at most* is
source-specific.** The rat/human range of 0.76–2.03 is the Ryu/Ryu pairing.
Across all four available human sources the rat binding factor reaches 4.89
(PFOA against Fischer et al., 2024), which lifts the worst-case binding
ceiling from 18.8% to 25.5% of the log half-life ratio. The conclusion
survives — binding remains a minority of the difference at every pairing, and
the Fischer comparison rests on only three shared compounds — but a claim
that binding is bounded at twofold is not defensible. Bounded at roughly
fivefold across every available source is.

All in vitro transport measurements were made without plasma protein, against
an in vivo unbound fraction near 0.0025 — a 400-fold discrepancy in the
quantity transporters respond to. Rankings may survive; magnitudes are in
doubt.

The human kinetic data derive from a single volunteer. A sensitivity analysis
against population values from a Bayesian hierarchical analysis of serum decay
shifts the partition but not the conclusion (binding below 10% of the log
ratio in every case, tubular residual 57–80%). The animal side draws on
several studies per species, and strains are mixed within the mouse column.

Decomposition shares are reported with bootstrap intervals over compounds but
do not propagate the experimental uncertainty in individual V_d and CL values.

The cynomolgus calculation assumes a GFR twice the human value per kilogram;
the conclusion would require that factor to be wrong by orders of magnitude to
change.

---

## Conclusions

The interspecies difference in PFAA elimination is carried by different
physiological terms in different species. For the rat, plasma protein binding
is bounded at twofold and contributes under a tenth of a difference reaching
three orders of magnitude, while the tubular term contributes about two
thirds. For the mouse the pattern inverts and binding contributes about half.
For the cynomolgus monkey, binding alone would require unbound fractions
exceeding unity.

The term carrying the rat difference corresponds to the one component of the
renal system whose molecular machinery is not conserved between the species:
the apical reabsorptive transporters, which are different proteins in human
and rat. Rodent reabsorption data therefore carry magnitude but not mechanism.

For two compounds the tubular term is no longer only a magnitude. Within a
single study, in one rat strain and one assay and with no sex difference in
volume of distribution, the tubular factor for PFOA and PFHxS changes sign
between the sexes — net reabsorption in the male, net secretion in the female
— while PFOS, which has no sex difference in renal clearance, does not. And
the sign of the dose response separates saturating reabsorption from
saturating secretion without requiring the two to be resolved numerically:
every rat series rises with dose except PFOS, whose four series across two
studies all fall, in the one case where saturable plasma binding can be
excluded as the cause. The tubular term is therefore partly mechanistic here,
not merely residual.

The tubular explanation is nonetheless compound-specific rather than general.
Within the compounds for which the renal clearance ratio can be decomposed
without a residual, the tubular term carries about two thirds of the
difference for PFHpA and PFOA, is exceeded by the glomerular term for PFHxS,
and for PFDA is bounded below unity, meaning the rat tubule is not the more
secretory of the two. Any extrapolation that treats rodent tubular handling
as a uniform property of the chemical class will therefore misstate at least
some compounds, in a direction that depends on chain length and head group
rather than on a single scaling factor.

Five measurements would sharpen these conclusions, and we state them in
advance: an unbound fraction in cynomolgus plasma, which would convert the
reductio of Table 5 into a decomposition; direct renal clearances in place of
total clearances, which would convert the residual into a tubular term; a
transporter assay at physiological albumin; an independent human unbound
fraction for PFPeA and PFHxA, the only compounds whose tubular sign is
presently in dispute; -- newly identified here as the binding constraint
rather than the rodent data -- human renal clearances for PFBA and PFBS,
without which neither compound can enter the residual-free decomposition
however completely the rat is characterised; and a second laboratory's
unbound fractions in **rat** plasma, the only species for which this work
rests on a single source with no possibility of substitution on the animal
side. A renal-clearance dose series for PFOS would additionally separate
saturable tubular secretion from saturable biliary elimination, the one
alternative the mechanism test cannot exclude. A negative result on the first — an unbound fraction
showing the macaque to be a poor human analogue — would be as informative as a
positive one.

---

## Supplementary Material

Supplementary Table S1 (per-compound kinetic values and sources), the complete
binding compilations, transporter tables, and all analysis code are available
in the project repository. Analyses are reproducible by idempotent scripts
verified byte-identical across consecutive runs.

| content | file |
|---|---|
| Species decomposition (Tables 2–5) | `db/species_decomposition.csv`, `scripts/stage15_species_decomposition.py` |
| Binding statistics (Table 1) | `db/binding_statistics.csv`, `scripts/stage14_binding_statistics.py` |
| Tubular factors, 52 pairs | `db/tubular_factor_consolidated.csv` |
| Binding compilation, all proteins and species | `db/protein_binding.csv` |
| Transporter abundance and kinetics | `db/transporter_expression.csv`, `db/transporter_clearance_stage3.csv` |
| Renal-only decomposition (Tables 9–10) | `db/renal_only_decomposition.csv`, `scripts/stage17_renal_only_decomposition.py` |
| Rat renal clearance search, bounds and cross-validation (Tables 11–12) | `db/rat_renal_clearance_search.csv`, `db/kudo_ohmori_crossvalidation.csv`, `db/renal_only_bounds.csv`, `scripts/stage18_rat_renal_clearance_search.py` |
| Open provenance flags | `db/renal_clearance_provenance_flags.csv` |
| Within-method L-FABP comparison | `db/lfabp_within_method.csv`, `scripts/stage19_lfabp_within_method.py` |
| Kim 2016 read at source; GFR sensitivity | `db/primary_2026/kim2016_table4_rat.csv`, `db/renal_only_v2.csv`, `db/gfr_sensitivity.csv`, `scripts/stage20_kim2016_primary.py` |
| Dose-response mechanism test | `db/dose_sign_test.csv`, `db/tubular_mechanism_verdict.csv`, `scripts/stage21_dose_sign_test.py` |
| Arithmetic bound and f_u source sensitivity | `db/binding_arithmetic_bound.csv`, `db/fu_source_sensitivity.csv`, `scripts/stage22_fu_sensitivity_and_bound.py` |
| Figures | `notebook/pfas_tk_plots.ipynb` |

---

## References

Abraham, K., Mertens, H., Richter, L., Mielke, H., Schwerdtle, T., and
Monien, B. H. (2024). Kinetics of 15 per- and polyfluoroalkyl substances
(PFAS) after single oral application as a mixture — a pilot investigation in
a male volunteer. *Environ. Int.* **193**, 109047.

ATSDR (2021). *Toxicological Profile for Perfluoroalkyls*. Agency for Toxic
Substances and Disease Registry, Atlanta, GA.

Chang, S.-C., Das, K., Ehresman, D. J., Ellefson, M. E., Gorman, G. S.,
Hart, J. A., Noker, P. E., Tan, Y.-M., Lieder, P. H., Lau, C., Olsen, G. W.,
and Butenhoff, J. L. (2008). Comparative pharmacokinetics of perfluorobutyrate
in rats, mice, monkeys, and humans and relevance to human exposure via
drinking water. *Toxicol. Sci.* **104**, 40–53.

Chang, S., Noker, P. E., Gorman, G. S., Gibson, S. J., Hart, J. A., Ehresman,
D. J., and Butenhoff, J. L. (2012). Comparative pharmacokinetics of
perfluorooctanesulfonate (PFOS) in rats, mice, and monkeys. *Reprod. Toxicol.*
**33**, 428–440.

Chengelis, C. P., Kirkpatrick, J. B., Myers, N. R., Shinohara, M., Stetson,
P. L., and Sved, D. W. (2009). Comparison of the toxicokinetic behavior of
perfluorohexanoic acid (PFHxA) and nonafluorobutane-1-sulfonic acid (PFBS) in
cynomolgus monkeys and rats. *Reprod. Toxicol.* **27**, 400–406.

EFSA CONTAM Panel (2020). Risk to human health related to the presence of
perfluoroalkyl substances in food. *EFSA J.* **18**, 6223.

Fujii, Y., Niisoe, T., Harada, K. H., Uemoto, S., Ogura, Y., Takenaka, K.,
and Koizumi, A. (2015). Toxicokinetics of perfluoroalkyl carboxylic acids with
different carbon chain lengths in mice and humans. *J. Occup. Health* **57**,
1–12.

Forsthuber, M., Kaiser, A. M., Granitzer, S., Hassl, I., Hengstschläger, M.,
Stangl, H., and Gundacker, C. (2020). Albumin is the major carrier protein for
PFOS, PFOA, PFHxS, PFNA and PFDA in human plasma. *Environ. Int.* **137**,
105324.

Harada, K., Inoue, K., Morikawa, A., Yoshinaga, T., Saito, N., and Koizumi,
A. (2005). Renal clearance of perfluorooctane sulfonate and perfluorooctanoate
in humans and their species-specific excretion. *Environ. Res.* **99**,
253–261.

Huang, M. C., Dzierlenga, A. L., Robinson, V. G., Waidyanatha, S., DeVito,
M. J., Eifrid, M. A., Granville, C. A., Gibbs, S. T., and Blystone, C. R.
(2019). Toxicokinetics of perfluorobutane sulfonate (PFBS),
perfluorohexane-1-sulphonic acid (PFHxS), and perfluorooctane sulfonic acid
(PFOS) in male and female Hsd:Sprague Dawley SD rats after intravenous and
gavage administration. *Toxicol. Rep.* **6**, 645–655.

Kim, S.-J., Heo, S.-H., Lee, D.-S., Hwang, I. G., Lee, Y.-B., and Cho, H.-Y.
(2016). Gender differences in pharmacokinetics and tissue distribution of 3
perfluoroalkyl and polyfluoroalkyl substances in rats. *Food Chem. Toxicol.*
**97**, 243–255.

Kudo, N., Suzuki, E., Katakura, M., Ohmori, K., Noshiro, R., and
Kawashima, Y. (2001). Comparison of the elimination between perfluorinated
fatty acids with different carbon chain length in rats. *Chem. Biol. Interact.*
**134**, 203–216.

Kudo, N., Katakura, M., Sato, Y., and Kawashima, Y. (2002). Sex
hormone-regulated renal transport of perfluorooctanoic acid. *Chem. Biol.
Interact.* **139**, 301–316.

OEHHA (2024). *Public Health Goals for Perfluorooctanoic Acid and
Perfluorooctane Sulfonic Acid in Drinking Water*. Office of Environmental
Health Hazard Assessment, California EPA, Sacramento, CA.

Ohmori, K., Kudo, N., Katayama, K., and Kawashima, Y. (2003). Comparison of
the toxicokinetics between perfluorocarboxylic acids with different carbon
chain length. *Toxicology* **184**, 135–140.

Olsen, G. W., Chang, S.-C., Noker, P. E., Gorman, G. S., Ehresman, D. J.,
Lieder, P. H., and Butenhoff, J. L. (2009). A comparison of the
pharmacokinetics of perfluorobutanesulfonate (PFBS) in rats, monkeys, and
humans. *Toxicology* **256**, 65–74.

Ryu, S., Burchett, W., Zhang, S., Modaresi, S. M. S., Agudelo Areiza, J.,
Kaye, E., Fischer, F. C., and Slitt, A. L. (2024a). Species-specific unbound
fraction differences in highly bound PFAS: A comparative study across human,
rat, and mouse plasma and albumin. *Toxics* **12**, 253.

Ryu, S., Burchett, W., Zhang, S., Jia, X., Modaresi, S. M. S., Agudelo
Areiza, J., Rodrigues, D., Zhu, H., Sunderland, E. M., Fischer, F. C., and
Slitt, A. L. (2024b). Unbound fractions of PFAS in human and rodent tissues:
Rat liver a suitable proxy for evaluating emerging PFAS? *Environ. Sci.
Technol.* **58**, 14641–14650.

Smeltz, M., Wambaugh, J. F., and Wetmore, B. A. (2023). Plasma protein binding
evaluations of per- and polyfluoroalkyl substances for category-based
toxicokinetic assessment. *Chem. Res. Toxicol.* **36**, 1133–1144.

Zhang, Y., Beesoon, S., Zhu, L., and Martin, J. W. (2013). Biomonitoring of
perfluoroalkyl acids in human urine and estimates of biological half-life.
*Environ. Sci. Technol.* **47**, 10619–10627.

Sundström, M., Chang, S.-C., Noker, P. E., Gorman, G. S., Hart, J. A.,
Ehresman, D. J., Bergman, Å., and Butenhoff, J. L. (2012). Comparative
pharmacokinetics of perfluorohexanesulfonate (PFHxS) in rats, mice, and
monkeys. *Reprod. Toxicol.* **33**, 441–451.

Thakur, A., Yue, G., Ahire, D., Mettu, V. S., Al Maghribi, A., Ford, K.,
Peixoto, L., Leeder, J. S., and Prasad, B. (2024). Sex and the kidney
drug-metabolizing enzymes and transporters: are preclinical drug disposition
data translatable to humans? *Clin. Pharmacol. Ther.* **116**, 235–246.

*[Remaining references to be completed from `db/citations.csv`, which holds
87 entries with DOI or PMID for all but one.]*
