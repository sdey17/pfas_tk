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
species, whereas the reabsorptive transporters are not shared — human OAT4 is
undetectable in rat and mouse, and rodent Oatp1a1 is undetectable in human.
Separately, we show that liver fatty acid-binding protein has no affinity
advantage over albumin (one-sided 95% upper limit 1.06-fold), removing the
standard explanation for hepatic accumulation. Reading the primary sources
behind a widely used regulatory tabulation yields direct human renal
clearances for nine compounds, giving net tubular reabsorption of 98–99.5%
from C8 to C12, and reveals that the tabulation transcribed one study's total
clearance column as renal clearance, publishing a renal clearance for one
compound whose urinary clearance was below detection. For the three compounds
with direct renal clearances in both species, the renal clearance ratio
decomposes with no residual at all, and binding's share falls further, to 3%
against 68% for a tubular term that is now measured rather than inferred. The practical
consequence of the main result is that no single interspecies scaling
approach serves both rodents.

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

### Decomposition of the human/rat half-life ratio

**Table 2.** Human/rat half-life ratio and its four factors, eight PFAAs with
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

**Table 3.** Share of the log₁₀ half-life ratio carried by each factor (%),
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
indistinguishable from zero (Table 3). The cause is direct: mouse/human
unbound-fraction ratios range from 2.0 (PFDoDA) to 18.6 (PFHxS), whereas
rat/human ratios across the same 14-compound panel span only 0.76–2.03
(Table 4). Both mouse strains agree closely except for PFOS, where CD-1 and
C57BL/6 differ fivefold.

**Table 4.** Plasma unbound fractions, single laboratory (Ryu et al., 2024a).

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
half-life difference (Table 5).

**Table 5.** Unbound fraction the cynomolgus monkey would require for binding
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

**Table 7.** Human renal clearance from three primary studies (mL/kg/day).

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
converts the human upper bounds into point estimates (Table 8). Every direct
value falls below its former bound, as it must.

**Table 8.** Human tubular factors from measured renal clearance.

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
clearance in male human and male rat (Table 9). PFNA has one in the human and
only in the female rat, and pairing across sexes is inadmissible in a species
where sex changes the sign of the tubular term.

**Table 9.** Renal-only decomposition, male human against male rat. Multiple
rows where sources disagree.

| compound | CL_R human (source) | CL_R rat (source) | ratio | T human | T rat | binding | glomerular | **tubule** |
|---|---|---|---|---|---|---|---|---|
| PFHpA | 0.674 (Fujii) | 983.1 (Ohmori) | 1459 | 0.0589 | 15.04 | 1.02 | 5.60 | **255.2** |
| PFOA | 0.030 (Harada) | 18.2 (Kemper) | 607 | 0.0048 | 0.424 | 1.22 | 5.60 | **89.0** |
| PFOA | 0.030 (Harada) | 21.2 (Ohmori) | 707 | 0.0048 | 0.494 | 1.22 | 5.60 | **103.7** |
| PFOA | 0.030 (Harada) | 46.08 (Kudo) | 1536 | 0.0048 | 1.074 | 1.22 | 5.60 | **225.4** |
| PFOA | 0.044 (Fujii) | 18.2 (Kemper) | 414 | 0.0070 | 0.424 | 1.22 | 5.60 | **60.7** |
| PFOA | 0.044 (Fujii) | 21.2 (Ohmori) | 482 | 0.0070 | 0.494 | 1.22 | 5.60 | **70.7** |
| PFOA | 0.044 (Fujii) | 46.08 (Kudo) | 1047 | 0.0070 | 1.074 | 1.22 | 5.60 | **153.7** |
| PFHxS | 0.027 (Zhang) | 0.76 (Kim) | 28 | 0.0151 | 0.069 | 1.09 | 5.60 | **4.6** |

All eight combinations of the disagreeing human and rat sources are shown
rather than a selected pair; the partition in Table 10 is the mean over all
eight. For PFOA the rat sources span 2.5-fold and the human two 1.5-fold, so its
tubular factor ranges from 61 to 225 across the six combinations — but its
*share* of the log ratio is insensitive to that, staying between 68% and 74%.
Across all three compounds the tubular share spans 46% (PFHxS) to 76%
(PFHpA).

**Table 10.** Share of the log renal clearance ratio, mean over all eight
source combinations.

| | binding | glomerular | **tubule** |
|---|---|---|---|
| renal-only (equation 5, no residual) | **3%** | 29% | **68%** |
| total clearance (equation 4, residual) | 8% | 42% | 65% |

**Removing the non-renal routes does not rescue binding.** Its share falls
from 8% to 3%, and the tubular term is confirmed at 68% as a measured
quantity rather than a residual. This is the strongest available form of the
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

**PFDA could not be added, and the attempt is recorded.** Ohmori's Figure 2
was re-examined at 300 dpi: the PFDA bars for both sexes sit flat on the
baseline carrying only their significance letters. That paper's own Figure 3
caption states that PFNA-male, PFDA-male and PFDA-female "are not
distinguished on the graph", and its published regression
(CL_R = 0.48·CL_tot + 15.63, r² = 0.981) cannot recover them: at PFDA's total
clearance of 5.2 mL/kg/day it predicts a renal clearance of 18.1, exceeding
the total by 3.5-fold. The intercept is fitted across totals reaching 3000
and swamps small values.

### Renal transporter expression: the reabsorptive machinery is not shared

Of 305 directional rows, 299 resolve to a reabsorption-or-secretion
assignment across 21 transporters with zero genuine conflicts: papers that
disagree about which compounds a transporter carries do not disagree about
where it sits.

**Table 6.** Renal transporter abundance, kidney membrane fractions, single
assay (Thakur et al., 2024). Values are pmol/mg membrane protein, mean ± SD.

| transporter | membrane | net process | human (n=15) | rat (n=10) | mouse (n=14) |
|---|---|---|---|---|---|
| OAT1 | basolateral | secretory | 4.8 ± 1.2 | 32.8 ± 13.7 | 19.1 ± 7.1 |
| OAT3 | basolateral | secretory | 4.3 ± 1.1 | 6.3 ± 1.1 | 6.5 ± 3.6 |
| **OAT4** | apical | **reabsorptive** | **0.5 ± 0.2** | **not detected** | **not detected** |
| **Oatp1a1** | apical | **reabsorptive** | **not detected** | **1.9 ± 2.1** | 3.1 ± 2.8 |
| URAT1 | apical | — | 3.1 ± 1.2 | 2.2 ± 0.7 | 10.8 ± 1.1 |

The secretory side is conserved in kind and differs in degree; the
reabsorptive side is not the same machinery. Because reabsorption determines
the long half-lives, the half of the renal system governing persistence is the
half that is not shared.

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

The tubular factor is not separable into reabsorption and secretion;
equation (3) has two unknowns and one equation, and the one in vitro route
that could separate them is closed because OAT4 has no scaling factor to the
kidney, while the two available methods for OAT1 and OAT3 disagree 11-fold and
in opposite directions.

The fourth factor is a residual rather than a measured tubular term for the
five compounds lacking paired renal clearances; for PFHpA, PFOA and PFHxS it
has been replaced by the residual-free decomposition of equation (5), which
gives the same answer. Its non-renal content is quantified rather than
assumed: renal elimination is
25–46% of total human PFOA clearance (Harada 2005; Fujii 2015) against 91% in
the male rat and 92% in the mouse, so the residual is substantially non-renal
in the human and substantially tubular in the rodents. Direct human renal
clearances are available for nine compounds (Table 7) but matching male rat
values for only three, which is what limits equation (5) to that set.

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

Four measurements would sharpen these conclusions, and we state them in
advance: an unbound fraction in cynomolgus plasma, which would convert the
reductio of Table 5 into a decomposition; direct renal clearances in place of
total clearances, which would convert the residual into a tubular term; a
transporter assay at physiological albumin; and an independent human unbound
fraction for PFPeA and PFHxA, the only compounds whose tubular sign is
presently in dispute. A negative result on the first — an unbound fraction
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
| Figures | `notebook/pfas_tk_plots.ipynb` |

---

## References

Abraham, K., Mertens, H., Richter, L., Mielke, H., Schwerdtle, T., and
Monien, B. H. (2024). Kinetics of 15 per- and polyfluoroalkyl substances
(PFAS) after single oral application as a mixture — a pilot investigation in
a male volunteer. *Environ. Int.* **193**, 109047.

ATSDR (2021). *Toxicological Profile for Perfluoroalkyls*. Agency for Toxic
Substances and Disease Registry, Atlanta, GA.

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

OEHHA (2024). *Public Health Goals for Perfluorooctanoic Acid and
Perfluorooctane Sulfonic Acid in Drinking Water*. Office of Environmental
Health Hazard Assessment, California EPA, Sacramento, CA.

Ohmori, K., Kudo, N., Katayama, K., and Kawashima, Y. (2003). Comparison of
the toxicokinetics between perfluorocarboxylic acids with different carbon
chain length. *Toxicology* **184**, 135–140.

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

Thakur, A., et al. (2024). Quantitative proteomic analysis of renal
transporters in human, rat and mouse kidney. *Clin. Pharmacol. Ther.* **116**,
235. [AUTHOR LIST TO BE COMPLETED]

*[Remaining references to be completed from `db/citations.csv`, which holds
72 entries with DOI and PMID.]*
