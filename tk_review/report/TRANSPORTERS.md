# PFAS renal transporters and the tubular axis

### What the measurements support, what they bound, and what they do not settle

Companion to `ALBUMIN_BINDING.md`. Working data are in `db/`; the
stage-by-stage record of how this was assembled, including every defect found
along the way and every claim later withdrawn, is `TRANSPORTERS_ROADMAP.md`.

Evidence base: **81 catalogued sources**, 58 transporter kinetic rows that
survive the inclusion rule (§1), a 36-row transporter localisation authority
table, measured plasma unbound fractions for 14 PFAS in four plasma matrices,
and in vivo renal or total clearances for three species consolidated into one
table of 44 verdicts over 74 evidence rows.

---

## Summary of findings

### Established — robust to every source choice tested

1. **The reabsorptive machinery is not the same genes in human and rat.**
   Human OAT4 is not detected in rat kidney; rat Oatp1a1 is not detected in
   human. **Each species reabsorbs PFAS with a transporter the other lacks** —
   a presence/absence difference, not a difference of degree.
2. **Transporter affinity is conserved** where a transporter exists in both:
   1.2 to 1.8-fold on OAT1 and OAT3, measured in one expression system.
3. **Substrate specificity is conserved** — the same ~500 g/mol threshold
   governs transporter interaction in both species.
4. **Abundance gates relevance.** The highest intrinsic clearance in the entire
   compilation belongs to a transporter that two laboratories cannot detect in
   human kidney.
5. **The rat sex difference is a change of sign, not of degree.** Female rats
   *net-secrete* PFOA and PFNA; males do not. Established twice over, by routes
   sharing no assumptions: Oatp1a1 is 114-fold more abundant in males, and the
   in vivo renal fractions put the female above the secretion threshold and the
   male below it.
6. **The two rodents are fast for different reasons.** Against the human, the
   rat's difference is carried by the tubule and the mouse's by the glomerulus.
   "Rodents clear PFOA faster than humans" is one sentence covering two
   mechanisms.
7. **Long-chain PFAS are almost completely reabsorbed in the human.** C8 through
   C11 exceed **98%** reabsorbed, and the C6 and C8 sulfonates **97%**, on every
   unbound-fraction source available. (The C4 sulfonate PFBS is not in this
   group: it sits at 87%, and short chains generally are unsettled — finding 11.)
8. **A QSAR for PFAS renal clearance cannot be a QSAR for albumin binding
   alone.** Binding explains 61% of the log-variance in the human clearance
   series and 83% in the mouse, but the tubular residual still spans 466- and
   823-fold. Aims 1 and 2 do not separate.

### Bounded — direction settled, magnitude not

9. **The rat/human tubular factor is of order 10²** — 22 to 161-fold across the
   four human and three rat clearance sources. Every combination leaves the
   tubule at least three-fold above the 6.8-fold glomerular factor, so the
   dominance is secure and the number is not. **It should be quoted as an order
   of magnitude.**
10. **The male rat reabsorbs far less filtered PFOA than the regulatory figure
    of 93.7%**, but how much less is unsettled: three rat studies give T = 0.42
    to 1.07 on OEHHA's glomerular filtration rate and 0.77 to 1.95 on the rate
    Ohmori 2003 cites. **At most half the filtered load is returned, possibly
    none.**

### Not established

11. **The sign of human tubular handling for the short chains (C5–C7).** Three
    laboratories' human unbound fractions disagree by up to 11-fold, and the
    filtered load scales with that term, so human PFHxA is either 52% reabsorbed
    or secreted at three times its filtered load depending on which is used.
    *This is now the binding constraint on the whole strand.*
12. **Whether human OAT4 reabsorbs more effectively in situ than rat Oatp1a1.**
    Per mg of transporter protein they are within 1.3-fold; per kidney they
    cannot be compared, because neither has a scaling factor.
13. **Whether any measured clearance survives contact with plasma protein.**
    Every assay in this compilation was run without it.
14. **Anything about the monkey.** No plasma unbound fraction has ever been
    measured in cynomolgus plasma for any PFAS, although 14 clearances exist.

### Errors found in the regulatory and secondary literature

Eight, listed with evidence in §9. The two with regulatory consequence:

- **ATSDR 2021 Table 3-6 misreports rat PFDA total clearance by forty-fold**
  (207 against the paper's 5.2 mL/kg·day) and male rat PFOA by 2.7-fold.
- **OEHHA PHG Table A6.4 assumes an unbound fraction of 0.02 for every
  species, uncited.** The measured values are 0.0025 (human) and 0.0030 (rat),
  so the filtered loads and therefore the printed reabsorption percentages are
  wrong by roughly eight-fold in the denominator.

---

## 1. The data, and the rule that decides what counts

The record holds 301 rows from 28 sources. One identity decides which can be
compared:

> **CL_int [µL/min/mg protein] = 1000 × V_max [nmol/min/mg] / K_m [µM]**

A source reaches the kinetics table only if it supplies what that identity
needs. Four tiers result:

| Tier | | rows |
|---|---|---|
| **A** | reaches CL_int per mg protein | **52** |
| B | a clearance on a basis that will not join (per 10⁶ cells) | 2 |
| C | affinity only — K_m, no clearance | 4 |
| D | substrate calls, uptake ratios, inhibition | 305 |

A, B and C are in `db/transporter_kinetics_working.csv`; **D is a separate
file**, because a substrate call is not a weak K_m. Argoul 2026 contributes 120
uptake-ratio rows and Ryu 2024 another 109 — together two thirds of the
compilation — and no arithmetic turns an uptake ratio into an affinity.

**K_i and IC50 do not get in.** An IC50 belongs as much to the probe substrate
and its concentration as to the inhibitor. Converting one to a probe-independent
K_i needs the probe's own K_m, and **none of the 20 numeric IC50 values in this
record carries it.** This is the same discipline the albumin and L-FABP sheets
apply to probe-displacement constants.

---

## 2. Direction, which has to come before kinetics

`db/transporter_localisation.csv`, 36 rows, one per **(transporter, species,
tissue)** — because direction is not a property of the transporter alone.
OATP2B1 is apical in human kidney, where it reabsorbs, and sinusoidal in human
liver, where it pulls PFAS out of blood. **OAT2 is basolateral in human kidney
and apical in rat**, so identical transport would secrete in one species and
reabsorb in the other.

Renal polarity follows Niu 2023 (*Environ Health Perspect* 131:076002), which
states it outright; tissue roles beyond the kidney follow OEHHA's PHG.

**343 of 363 rows join to a membrane and a net process, with zero direction
conflicts.** That is itself a result: the papers that disagree sharply about
*which compounds* a transporter carries do not disagree at all about *where it
sits*. Substrate identity is contested; localisation is not.

---

## 3. Expression, measured in one assay across three species

Thakur 2024 (*Clin Pharmacol Ther* 116:235) quantified 78 enzymes and
transporters in kidney membrane fractions of human (n=15), rat (n=10) and mouse
(n=14), both sexes. pmol/mg membrane protein:

| | human | rat | mouse | |
|---|---|---|---|---|
| **Oat1** (basolateral, secretory) | 4.8 ± 1.2 | **32.8 ± 13.7** | 19.1 ± 7.1 | rat 6.8× |
| **Oat3** (basolateral, secretory) | 4.3 ± 1.1 | 6.3 ± 1.1 | 6.5 ± 3.6 | rat 1.5× |
| **OAT4** (apical, reabsorptive) | **0.5 ± 0.2** | **not detected** | not detected | |
| **Oatp1a1** (apical, reabsorptive) | **not detected** | **1.9 ± 2.1** | 3.1 ± 2.8 | |
| URAT1 (apical) | 3.1 ± 1.2 | 2.2 ± 0.7 | 10.8 ± 1.1 | transports no PFAS |

Independent agreement where sources overlap: human OAT1 4.8 (Thakur) against
5.3 (Niu 2026); P-gp 1.2 against 2.1; BCRP not detected in both.

**Two consequences.**

**The secretory side is conserved in kind and differs in degree** — both species
have OAT1 and OAT3, the rat simply has more. **The reabsorptive side is not the
same machinery at all.** Niu 2023 calls OAT4 "a human-specific transporter …
it may be one of the transporters responsible for species-specific differences
in PFAS toxicokinetics"; EFSA says "OAT4 (not expressed in rats) and urate
transporter 1 (URAT1), but not Oatp1a2 (the closest orthologous of rat
Oatp1a1), were shown to mediate PFOA". Here it is measured rather than cited.

**Two long-standing disagreements dissolve.** ASBT and OSTα/β were not detected
in human *or* rat kidney. Zhao 2015 found transport by both in transfected
cells; Argoul 2026 found none. **Both are right** — a transfected-cell result
states what a protein *can* carry, not whether it is present to do it.

---

## 4. Clearance, and why abundance gates it

52 tier-A rows carry an intrinsic clearance. Eleven can be scaled to kidney —
human OAT1 and OAT3, by relative expression or relative activity factor — and
the two methods **disagree 11-fold on OAT3 and in opposite directions**, so
they are reported separately and never averaged.

**The expression gate removes the compilation's three highest values**, all
BCRP from Niu 2026, including **279 µL/min/mg for PFOA — the highest anywhere
in it.** Human renal BCRP is not detected by Thakur or by Niu. The fastest
transporter in the table is one two laboratories cannot find in human kidney.

**Reabsorption against secretion**, within one laboratory (Louisse 2023 for
OAT4, 2024 for OAT1/OAT3), measured and unscaled:

| | OAT4 (reabsorb) | OAT3 | OAT1 | OAT4/OAT3 |
|---|---|---|---|---|
| C6 PFHpA | 75.0 | 8.7 | 21.4 | 8.6 |
| C6 PFHxS | 79.3 | 7.5 | — | 10.6 |
| C7 PFOA | 95.7 | 16.7 | 30.8 | 5.7 |
| C8 PFNA | 146.6 | 23.8 | — | 6.2 |
| C8 PFOS | 45.8 | — | 53.5 | **0.9** |
| C9 PFDA | 153.8 | — | — | |

Both sides rise with chain length — OAT3 2.7-fold C6→C8, OAT4 2.1-fold C6→C9 —
so chain length does not by itself tip the balance. **What stands out is the
level, and it splits by head group**: for every carboxylate, OAT4 carries
3.1–8.6 times the clearance of a secretory transporter. PFHxS behaves the same.
**PFOS is the exception** — OAT1 matches OAT4, and OAT3 does not transport it.

**The caveat that bounds this.** OAT4 cannot be scaled and OAT1/OAT3 can.
Applying the expression factor lifts OAT1 to 118–334 µL/min/mg, which would
overturn the ranking. The comparison is suggestive, not decisive, and **an OAT4
abundance in the expressing cells would settle it.**

---

## 5. The species question

### What has to be explained

The male rat clears PFOA renally **303 times faster** than the human. That
number is observed and factors exactly in two:

> **303× = (filtered load) × (fraction escaping the tubule)**

### Three candidate explanations, two eliminated

**Affinity is not it.** Nakagawa 2008 measured human and rat OAT1 and OAT3 in
the *same* expression system — the only way to compare species without also
changing laboratory. Correcting a units error in its Table 1 (a K_m/V_max column
labelled with the units of V_max/K_m, which inverts the ranking), the intrinsic
clearances are:

| | human | rat | human/rat |
|---|---|---|---|
| OAT1 | 9.99 | 8.26 | 1.21× |
| OAT3 | 8.13 | 4.56 | 1.78× |

Within two-fold, against a ~50-fold difference in half-life. The paper's own
conclusion: species differences "are not attributable to affinity differences
in these OATs".

**Substrate specificity is not it.** Argoul 2026 tested ten PFAS against eight
human and four rat transporters. The same ~500 g/mol threshold holds in both;
PFDA, PFOS and PFDS are substrates of none of the twelve, in either species.

**What is left is the machinery itself** — and it differs in kind, not degree
(§3). That is a presence/absence difference, and a far better candidate for a
50-fold gap than anything in the first two.

### The sex difference says a reabsorptive transporter can do this

Within the rat: same GFR, same binding, same everything except sex. Escape goes
from partial reabsorption in males to **net secretion** in females — a 37 to
48-fold difference in renal clearance, **with a change of sign**. (The escape
fraction is given as a direction rather than a percentage because §6 shows the
male figure is not settled.) What differs:

| | male/female |
|---|---|
| **Oatp1a1** (apical, reabsorptive) | **114×** |
| Oatp1a3 (apical, reabsorptive) | 7× |
| Oat1 (basolateral, **secretory**) | 2× |

**Note the direction of the Oat1 difference.** Males have *more* secretory
transporter and still clear PFOA 37–48× more slowly. The reabsorptive side is
not merely contributing — it is overriding a secretory difference pointing the
other way. A difference in abundance alone cannot flip net handling from
reabsorption to secretion; a transporter effectively absent in one sex can.

It fits the endocrinology: Cheng 2006 finds renal gender-divergent Oatp
expression "exclusively caused by stimulation by androgens", and Kudo 2002
found castration raised renal clearance 14-fold. And the human control is
clean — **no human renal enzyme or transporter showed any sex difference**
(n = 11 male, 4 female), against 28 in rat. Humans show no PFAS half-life sex
difference either.

**This is the strongest mechanistic result in the strand**, because it is the
one case where the mechanism is unambiguous and the prediction is qualitative.

### A second, independent route to the same sex difference

Everything above comes from transporter expression. Ohmori 2003's in vivo
clearances reach the same conclusion with no transporter data in them at all.
Net secretion requires the renal fraction of elimination to exceed
(f_u × GFR)/CL_tot, and Kudo 2001 measured that fraction in the same strain
and laboratory:

| | threshold | measured renal fraction | |
|---|---|---|---|
| **PFOA female** | 1.9% | **80%** | **net secretion** |
| PFOA male | 85.0% | 55% | fails the test |
| **PFNA female** | 28.6% | **51%** | **net secretion** |
| PFNA male | impossible | 2% | reabsorption certain |

**The female rat kidney secretes PFOA and PFNA; the male rat kidney reabsorbs
them.** The 70-fold sex difference in rat PFOA half-life is the difference
between a secreting kidney and a reabsorbing one. The female PFOA verdict
survives a ten-fold error in the unbound fraction; the PFNA one survives only
1.8-fold.

---

## 6. The tubular axis, and what the regulatory figures get wrong

### The quantity, and why there is no third term

Renal clearance factors exactly:

> **CL_R = f_u × GFR × T**, so **T = CL_R / (f_u × GFR)**
>
> | | |
> |---|---|
> | T < 1 | net tubular **reabsorption**, of (1 − T) of the filtered load |
> | T = 1 | filtration alone |
> | T > 1 | net tubular **secretion** |

T is reported per species throughout, not only as a between-species ratio: a
ratio of 45× hides whether both species reabsorb or one secretes, and those are
different mechanisms.

**T is deliberately left lumped.** Writing FR for fractional reabsorption and S
for secretory clearance, T = 1 − FR + S/(f_u × GFR), and **only the combination
is identifiable from in vivo data.** Splitting it needs an independent
measurement of FR or S *in the same species*; scaled in vitro CL_int is not
that measurement, because the expression gate and the REF/RAF choice move it
about five-fold and the albumin effect (§8) is unquantified for PFAS. A
three-term model fitted on it would reintroduce exactly the failure mode of the
assumed f_u. **The in vitro work constrains T rather than decomposing it.**

### Three tiers of evidence, and two of them need no renal clearance

| tier | evidence | gives |
|---|---|---|
| **direct** | a study in that species reports renal clearance | T outright |
| **upper** | only total clearance; CL_R ≤ CL_tot always | T ≤ CL_tot/(f_u·GFR) — below 1, reabsorption holds whatever the renal share |
| **lower** | a measured renal *fraction* | secretion once f_ren > (f_u·GFR)/CL_tot |

The second and third tiers were what made the rat tractable: it had direct
renal clearance for one compound of fourteen, and now has a verdict for nine.

Two rules are enforced by the build rather than by care: **no row mixes
species** (the unbound fraction, the filtration rate and the clearance must all
come from one), and **a paper beats any compilation of it**, with the
compilation's number retained and marked superseded.

### The assumed unbound fraction, and what measuring it does

The figures in circulation — human 99.8% of filtered PFOA reabsorbed, male rat
93.7% — come from OEHHA PHG Table A6.4, computed as f_u × GFR − CL_R with
**f_u assumed to be 0.02 in every species**, uncited. Ryu 2024 measured it:
14 PFAS by presaturation equilibrium dialysis in human, rat and two mouse
strains.

| | f_u measured | f_u assumed | ratio |
|---|---|---|---|
| human | 0.00245 | 0.02 | 8.2× |
| rat | 0.00298 | 0.02 | 6.7× |
| mouse (C57BL/6) | 0.0142 | 0.02 | 1.4× |

**The assumption is wrong in value and approximately right in structure for the
rat–human pair** — those two differ only 1.22-fold — and **wrong in structure
for the mouse**, which sits 5.8-fold above the human. Across all fourteen
compounds the pattern holds: rat/human has a geometric mean of 1.20× and lies
within two-fold for **13 of 14**; mouse/human is 3.7–4.1× and within two-fold
for only **3 of 14**. So holding f_u equal is a sound approximation for the rat
and a bad one for the mouse, **across the whole class** — not a lucky break
on PFOA.

The mouse error is also chain-length-dependent rather than constant: 3.4× at C4,
a 4.8–7.3× plateau over C5–C8, falling to ~2× by C9–C12. **It is largest over
the compounds of greatest regulatory interest.**

### The two rodents differ from the human for different reasons

| PFOA | observed CL_R ratio | glomerular | tubular |
|---|---|---|---|
| human vs **rat** | 303× | 6.8× | **~10²** |
| human vs **mouse** | 85× | **37.7×** | 2.3× |

The mouse differs from the human in *both* glomerular factors — f_u is 5.8×
higher and GFR 6.5× higher, and 5.8 × 6.5 = 37.7 — so the glomerulus accounts
for nearly the whole difference. The rat differs in *one*: its f_u matches the
human, so its glomerular term is almost pure GFR and the tubule must carry the
rest.

**A PBPK model calibrated on the mouse tubule and applied to the rat, or the
reverse, imports the wrong term.** This is the first result in the strand that
separates the two rodents mechanistically rather than quantitatively.

Across the human–mouse pairs the glomerular factor dominates six of seven, and
for four of the six long chains the tubular factor is *below one* — the tubule
partly cancels the mouse's faster clearance rather than amplifying it. For PFDA
the cancellation is nearly exact (13.4× glomerular, 0.06× tubular, 0.9×
observed), which is why PFDA has almost the same renal clearance in mouse and
human despite a 13-fold higher filtered load.

### Why the rat percentage should not be quoted

Three successive corrections each loosened it, and the pattern is itself the
finding.

| source for rat male PFOA CL_R | T (OEHHA GFR) | T (Ohmori GFR) |
|---|---|---|
| Kemper 2003, unpublished — **OEHHA's choice** | 0.42 | 0.77 |
| Ohmori 2003, digitised by Han 2012 | 0.49 | 0.90 |
| **Kudo 2002** — the only one cross-checked against its own total clearance | **1.07** | **1.95** |

The filtration clearance sits *inside* that span, so **the sign of the male rat
tubular term depends on which primary source is taken.** And the rat GFR is
itself not agreed: OEHHA uses 10 mL/min·kg, Ohmori 5.5 citing Shim 1984 — a
1.82-fold fork in a quantity treated throughout as known.

**The human tubule reabsorbs over 98% of filtered PFOA under every combination
of these sources** — 98.1% on the loosest, 99.3% on the tightest. The direction of the species contrast is secure; the rat
number is not, and the PHG's 93.7% is the most reabsorptive of three available
answers rather than a consensus.

---

## 7. Chain length, and the bridge to the QSAR aim

### The question, posed exactly

Because CL_R = f_u × GFR × T and GFR is a species constant, **all variation
across a compound series within one species sits in binding and tubular
handling.** If binding carried almost all of it, a PFAS QSAR for renal
clearance would be a QSAR for albumin binding and the two aims would separate.

| | clearance spread | binding spread | **tubular spread** | binding r² |
|---|---|---|---|---|
| human (13 compounds) | 3516× | 2180× | **466×** | **0.61** |
| mouse (7 compounds) | 52154× | 604× | **823×** | **0.83** |

**Binding is the larger term in both species — and not the only one.** The
tubular residual spans 466- and 823-fold, far too wide to treat as a constant.
**So aim 2 needs the transporter layer that aim 1 has built, and the two aims
do not separate.**

### What is settled along the series, and what is not

![Each bar spans the range of T across that species' own sources. A bar wholly below the dashed line establishes net reabsorption; a bar crossing it means the sign is not determined by the data. The long chains are settled; the short chains are not, and the spread is binding uncertainty, not clearance uncertainty.](fig_T_vs_chain_length.png)

**The long chains are settled, and one measurement settled them.** Abraham 2024
followed a single subject for 450 days with plasma, urine **and faeces** all
measured, reporting per compound the urinary, total and faecal clearances, an
unbound fraction, and the filtration clearance computed from *that subject's
own* creatinine clearance of 113.7 mL/min. **The clearances and the filtration
rate are certainly that subject's**, which is unique in this strand.

**But the binding term is not.** ✅ *Resolved from the papers.* Abraham
Table 4 carries the footnote *"as published by Smeltz et al. (2023)"*, and its
Methods state filtration clearance *"was calculated as product of the fraction
unbound in plasma (f_up, available from Smeltz et al. 2023) and the GFR
measured in"* the volunteer. **Abraham is therefore not an independent binding
series** — its clearances and GFR are the subject's own, its binding is EPA's.
This review's earlier claim that Abraham's T "needs no cross-study join" was
wrong: it is a join, within species, between Abraham's clearances and Smeltz's
binding. There are **four** independent human binding datasets, not five:

| | method | compounds |
|---|---|---|
| **Smeltz 2023** (US EPA) | ultracentrifugation | **67** |
| Fischer 2024 | C18 solid-phase microextraction | 29 |
| Argoul 2025 | Dianorm equilibrium dialysis | 16 |
| Ryu 2024 | presaturation equilibrium dialysis | 14 |

**The short chains have no determined sign.** Three laboratories measured human
plasma f_u and disagree by up to 11-fold, worst at the short chains:

On Abraham's own urinary clearances and the subject's own GFR of 113.7
mL/min, so that only the binding term changes between rows:

| human PFHxA | f_u | T | |
|---|---|---|---|
| **Smeltz 2023** (EPA; Abraham's source) | 0.0068 | **3.47** | secretion |
| **Fischer 2024** (derived) | 0.0095 | **2.48** | secretion |
| **Argoul 2025** | 0.0200 | **1.18** | secretion |
| Ryu 2024 | 0.0380 | 0.62 | reabsorption |

**Three of the four now give net secretion, and Ryu is the outlier.** An
earlier draft of this report presented 0.48 (Ryu) first and 3.48 as the
alternative; that ordering was wrong. PFPeA behaves the same way.

**The disagreement is confined to exactly two compounds.** PFPeA and PFHxA.
Every PFAS of seven carbons and above gives net reabsorption on all four
datasets, so **no conclusion about PFOA, PFNA, PFDA, PFUnDA, PFHxS or PFOS
depends on which binding source is chosen.**

### The uncertainty is binding, not clearance

Abraham's urinary clearance of PFHxA is **47.1 mL/kg·day**; Argoul 2026 reports
**47.0** from an unrelated population study — agreement to better than 1%.
**The renal clearances agree across laboratories where they overlap; the
unbound fractions do not.** That is the single most useful thing this
comparison says about where effort should go.

~~Binding affinity rises with chain length, so within the carboxylate series
f_u must fall monotonically — a constraint on the data independent of any
laboratory's protocol. Ryu's series has zero inversions in 8 pairs, Argoul's
1–2 and Abraham's 4, so the test prefers Ryu.~~

**That test was invalid, and is withdrawn.** Fischer 2024 measured
albumin–water distribution across ηpfc 4–13 and its own series **rises to
PFUnDA and then falls** — log D_HSA/w of 5.18 at ηpfc 10, 5.05 at 11, 4.54 at
12, 3.87 at 13 — with the mechanism stated: *"for ηpfc > 8 PFAS, reduced
molecular flexibility limits their accessibility to the HSA binding pockets."*
Smeltz 2023 reports the same independently — binding lower for PFCAs of chain
length ≥ 11. **So "f_u must fall monotonically" is false above about eight
perfluorinated carbons, which is exactly where this review applied it.**

Argoul's PFUdA value remains an outlier — 40-fold above Ryu's and 130-fold
above Fischer's — but on the evidence of two other laboratories, not on a
premise that turns out to be wrong.

**The position adopted here is therefore:** Abraham's *faecal* term to tighten
the long chains (a mass balance, not a binding assay, and nothing here
questions it); Ryu's f_u for the filtered load in the tables, as the series
with the widest compound coverage, **with the explicit note that it sits at the
high end of the short-chain range**; and the short chains reported as
sign-undetermined.

### Two inequalities that survive every choice

1. **The C6→C8 carboxylate drop in T is real** — at least 26-fold on Abraham's
   own numbers and 73-fold on Ryu's. The human kidney returns about half of
   filtered PFHxA and essentially all filtered PFOA.
2. **Mouse PFHxA is net secretion** (T = 2.23 on Ryu's binding, 2.71 on
   Argoul's own) — its renal clearance exceeds its filtered load whichever
   source is used. This is the only direct secretory verdict in the mouse
   column.

*An earlier draft read those two together as a clean human-reabsorbs /
mouse-secretes species contrast tied to OAT4. That is withdrawn: on the human
subject's own binding the human secretes PFHxA too. More generally, **every
sign verdict near T = 1 is a verdict about a binding measurement**, because
the filtered load scales with f_u and nothing else in the expression is in
doubt.*

---

## 7b. Which human unbound fraction to use

`scripts/stage9_human_fu_resolution.py`, `scripts/stage9b_build_recommended_fu.py`
→ **`db/human_fu_recommended.csv`**

### The recommendation is by purpose, because the four datasets are not interchangeable

| use | source | why |
|---|---|---|
| **human absolute values — the default** | **Smeltz 2023 (UC)** | 67 PFAS so the series has no gaps; US EPA's own; the method choice is *argued* — the paper runs UC and RED side by side on 22 compounds and shows RED failing to equilibrate for high-mass, high-logK_ow PFAS; Bayesian credible intervals exist per compound; and Abraham's in vivo mass balance already runs on it |
| **individual variability** | **Fischer 2024** | measures albumin and globulin binding *separately* and returns f_u from a person's own serum protein concentrations, R² = 0.99 |
| **species ratios** | **Ryu 2024, and only Ryu** | the only dataset measuring human, rat and two mouse strains in one laboratory; in a *ratio* a systematic method bias divides out |

**Do not use Ryu for human absolute values.** It runs at the high end — 4× above
Smeltz for PFHxA, 15× above Fischer for PFPeA — and those are precisely the two
compounds in dispute.

Two caveats travel with Smeltz: **PFBS is 0.0046 in its Table 2 and 0.0128 in
its Table 3**, a 2.8× internal discrepancy in the same assay; and its PFOS
value is an average over three salt forms.

### On a variable f_u — the instinct is right, and there are three different things it can mean

1. **Variable across compounds.** Settled, and the single largest correction to
   the regulatory position: OEHHA's 0.02 against a measured range of 0.0004
   (PFHpA) to 0.31 (PFPA) — **a factor of 780 across the series.** Already
   implemented throughout this review.
2. **Variable across individuals.** This is the one worth adding. Fischer's
   equation returns f_u from measured albumin and globulin, and f_u varies up
   to **2.5× across 1,919 NHANES participants**, widest below 500 g/mol —
   which is the short-chain range. `human_fu_recommended.csv` carries the band
   per compound (2.5× below 500 g/mol, 1.6× from 500–600, 1.3× above).
   **Note the direction before using it:** a higher f_u means a higher filtered
   load and faster elimination, so the individuals with the *lowest* f_u carry
   the longest half-lives. For a health-protective value the low tail is the
   sensitive one, not the high.
3. **Variable with concentration.** Real — Qin 2024 shows binding becomes
   non-linear above a PFAS:protein molar ratio near 0.1 — but at environmental
   exposures the ratio is ~0.0005, so this bites for *comparing in vitro
   datasets*, not for human dosimetry. Lowest priority of the three.

### What the choice actually decides

**Two compounds.** PFPeA and PFHxA. Ten of the twelve PFAS with enough data
give the same tubular sign on every binding dataset, and all ten are
reabsorption. The decision is far narrower than the 780-fold spread in f_u
makes it look.

---

## 8. A caveat that may matter more than the scaling

Every transport assay in this compilation was run **without plasma protein**.
Of the 52 tier-A rows across six studies, not one mentions albumin, serum or
BSA; where the buffer is stated outright it is protein-free.

Tan 2024 showed that with 4% human serum albumin present *and the unbound
concentration held constant*, unbound intrinsic active uptake clearance rose
**3.5 to 122-fold for OAT1** and up to 28-fold for OAT3 — but only for highly
protein-bound substrates, f_u between 0.002 and 0.18. The mechanism proposed is
facilitated dissociation at the cell surface, so the transporter sees more than
the bulk free fraction.

**PFOA's plasma f_u is 0.0015 — at or below the bottom of that range**, and
PFOS, PFNA and PFHxS are lower still. PFAS are the extreme case of exactly the
class this was demonstrated on.

If it applies, the clearances in §4 are **lower bounds**, understated *most* for
the longest-chain, most tightly bound compounds — the same direction as the
chain-length trend. That is uncomfortable, because part of what we would
attribute to transport could be an assay artefact.

**It also joins this document to its companion.** The albumin association
constants in `ALBUMIN_BINDING.md` are exactly the input a facilitated-
dissociation correction requires. The binding work stops being only a
distribution story. The test is simple: one PFAS uptake assay, with and without
albumin.

---

## 9. Errors found in the source literature

Each was caught by an internal inconsistency in the source itself or between a
source and its own primary, not by preference for one number over another.

**1. ATSDR 2021 Table 3-6 misreports Ohmori 2003 by up to forty-fold.**

| | Ohmori Table 1 | ATSDR | |
|---|---|---|---|
| PFHpA male / female | 1603.9 / 3070.5 | 1604 / 3071 | ✓ |
| **PFOA male** | **50.5** | **135** | 2.7× |
| PFOA female | 2233.2 | 2233 | ✓ |
| PFNA male / female | 6.9 / 105.7 | 6.9 / 106 | ✓ |
| **PFDA male / female** | **5.2 / 5.3** | **207 / 140** | **40× / 26×** |

Five of eight are right. **The PFDA error was load-bearing**: 207 was the sole
reason rat PFDA registered as uninformative, and at the true value it closes as
certain reabsorption.

**2. OEHHA PHG Table A6.4 assumes f_u = 0.02 for every species, uncited.** The
measured values differ from it by 6–8-fold for human and rat (§6), so every
filtered load in that table — and therefore every printed reabsorption
percentage — is wrong by roughly that factor in the denominator. The *ordering*
of species survives; the numbers do not.

**3. Nakagawa 2008 Table 1 carries a units defect.** A column of K_m/V_max
values is labelled with the units of V_max/K_m, which inverts the ranking.
Corrected, the human/rat intrinsic clearances are 1.21× (OAT1) and 1.78×
(OAT3) — the basis for finding affinity conserved.

**4. Niu 2026 transposes P-gp and BCRP** between its text and abstract on one
hand and its figures on the other. Verified across all six values; the
abstract's own K_m ranges and OAT1 arithmetic confirm the figures are right.

**5. Argoul 2026's citation for its unbound fractions does not resolve.** Its
reference list carries the LC-MS/MS method paper twice, as 2024 and 2025, with
identical titles and DOI. The actual f_u source is *Environ Toxicol Pharmacol*
**120**:104879.

**6. Abraham 2024 misquotes Smeltz in its own Discussion.** Its §4.5 says the
PFCA unbound fractions "ranged from 0.3108 for PFBA". **0.3108 is Smeltz's
value for PFPA** — perfluoropropanoate, C3. PFBA is 0.0930, which is what
Abraham's own Table 4 correctly uses. A text error only; no arithmetic depends
on it.

**7. Smeltz 2023 reports PFBS twice, 2.8-fold apart.** Its Table 2 gives
f_up = 0.0046 and its Table 3 gives 0.0128 for the same ultracentrifugation
assay. Abraham used the Table 2 value. Both are recorded in
`db/smeltz2023_fu.csv`; neither is preferred here, because nothing in the paper
distinguishes them.

**8. `pdftotext` silently drops the micro sign**, and the loss is internally
consistent so it cannot be detected from the text alone. It converted Uchida
2016's `fmol/µg` to `fmol/mg` — a 1000-fold error caught only because it would
have placed rat OAT1 500-fold below human — and rendered Ohmori 2003's dose of
48.64 µmol/kg as "48.64 mmol/kg", which for PFOA would be 20 g/kg. **Any
quantitative extraction from PDF in this project is checked against the
arithmetic, not against the character stream.**

Two of these have regulatory consequence (1 and 2); the rest would have
propagated into this review had they not been caught.

---

## 10. What is open, and what single measurement would close it

| Measurement | What it decides | Status |
|---|---|---|
| ~~A human f_u series that settles the short-chain disagreement~~ | **Decided (§7b): Smeltz 2023 as the default, Fischer for individual variability, Ryu for species ratios only.** The residual dispute is PFPeA and PFHxA alone | **closed** |
| **One equilibrium-dialysis experiment on cynomolgus plasma** | opens a third species for the entire framework; the monkey's PFOA half-life is closest to the human, and 14 clearances already exist | **highest-value single experiment** |
| **One PFAS uptake assay ± albumin** | whether every clearance in this compilation is a lower bound (§8) | not done for PFAS |
| **OAT4 abundance in the expressing cells** | whether OAT4 really dominates human renal handling, or scaled OAT1 overtakes it (§4) | does not exist |
| Rat f_u by sex for PFOA | the one input the sex-difference result assumes rather than measures; sex-invariant for PFHxS (1.10×), unmeasured for PFOA | — |
| The eight renal clearances in Ohmori 2003 Fig. 2 | replaces four rat bounds with direct values and adds the first female-rat rows | plotted as a bar chart, no numbers in the text |
| Rat Oatp1a1 REF | whether the two reabsorptive sides can be compared per kidney at all (§5) | does not exist |

**Four of the first five are each a single experiment**, and each decides
something this report has to leave open.

### A literature search against these four items

`scripts/stage7_open_items_literature.py` → `db/open_items_literature.csv`

**Item 11 is not exhausted — there are five human f_u series, not three.**

| | human PFOA f_u | method |
|---|---|---|
| **Fischer 2024** (Harvard) | **0.00061** | C18 SPME, HSA and globulin separately |
| Abraham 2024 | 0.00100 | in vivo, one subject |
| Argoul 2025 | 0.00170 | equilibrium dialysis |
| Ryu 2024 | 0.00245 | presaturation equilibrium dialysis |

plus **Smeltz 2023** (US EPA), **67 PFAS by ultracentrifugation with Bayesian
uncertainty bounds** — for a regulatory conversation, the obvious reference
dataset, and this review did not have it.

**And the question changes shape.** Fischer predicts f_unbound from an
individual's *own* albumin and globulin concentrations (R² = 0.99 against
measured serum distribution) and finds it varies **up to 2.5-fold across 1919
NHANES participants** — most for compounds below 500 g/mol, which is exactly
the short-chain range in dispute. Compounds under 500 g/mol are >98%
albumin-bound; globulins matter only above 600.

**So OEHHA's assumed f_u is wrong in kind as well as magnitude.** It is not
that 0.02 is the wrong constant — it is that there is no constant to be right
about, and the variability is largest precisely where §7 found the sign
undetermined.

Qin 2024 adds a quantified mechanism: PFAS binding to human plasma is
**concentration-dependent**, specific binding dominating below a PFAS:protein
molar ratio of 0.1 and nonspecific binding above 1. Fischer independently holds
its own ratios below 0.004–0.009, citing a physiological ratio of 0.0005 and
BSA non-linearity from 0.39 upward. **Checked against the disagreement: Ryu
dosed at a ratio near 0.008 and Abraham's subject carried near 2 × 10⁻⁶, so
both are inside the linear range and concentration does *not* explain the
5.6-fold PFHxA gap.** The mechanism is real and now quantified; it is not this
gap's cause.

**Item 12 is half-closed, and I had it wrong.** I previously reported that no
absolute quantification of OAT4 exists. The *kidney* side does, and is
replicated: **Prasad 2016 gives 0.5 ± 0.2 pmol/mg in human cortex (n = 41) and
Thakur 2024 gives 0.5 ± 0.2 (n = 15)** — two laboratories, nine years apart,
identical to the reported precision. What is missing is only the **in vitro**
side: OAT4 abundance in the cell line that supplied the kinetics. That is one
proteomics run on an existing line, not a new study.

**Item 13 is confirmed open.** Nothing PFAS-specific found; Tan 2024 remains
the only albumin-mediated-uptake demonstration and it used drugs. Fischer 2024
proposes the same experiment in its closing paragraph and names renal tubular
reabsorption as the target.

**Item 14 is confirmed absent — and the prize is larger than stated.** No
plasma protein binding measurement exists for any macaque, for any PFAS. But
**Butenhoff 2004** measured cynomolgus PFOA in **serum, urine and faeces** at
steady state across three oral doses plus an IV arm — the same mass-balance
design as Abraham 2024 in the human. So one equilibrium-dialysis run on macaque
plasma would not merely add a species with a bounded T; it would add one with
the full three-tier treatment **including the faecal term**. This search raises
rather than lowers that experiment's value.

**The concrete next step is a reading task, not a search:** ingesting Fischer
2024 and Smeltz 2023 would take the human series from three to five and let
§7's monotonicity test run on all five. It does not pre-empt the OEHHA
question — it sharpens it, from *which of three series* to *which of five, and
whether a single value should be used at all.*

---

### What the strand does not attempt

T is a **renal** accounting and does not describe the longest chains. Abraham's
mass balance puts 58–61% of the elimination of PFNA, PFDA, PFUnDA and PFDoDA in
the faeces; human PFDoDA has a total clearance 3.1 times its filtered load, so
its elimination cannot be mostly glomerular. **Anything above C11 belongs in a
mass-balance model with a faecal term**, which is the third aim, not this axis.

---

## 11. First QSAR fit on the tubular endpoint

`scripts/qsar0_build_modelling_table.py`, `scripts/qsar1_censored_fit.py`
→ `db/qsar/qsar_modelling_table.csv`, `db/qsar/qsar_censored_fit.csv`

**log₁₀ T ~ chain length + head group**, per species, by maximum likelihood
**with censoring** — 30 of the 39 usable rows are one-sided bounds, and fitting
them as point values would bias every coefficient toward the censoring
direction.

| species | n | censored | nC effect | sulfonate vs carboxylate | σ |
|---|---|---|---|---|---|
| human | 19 | 15 | ×0.39 per carbon (z = −1.65) | ×0.05 (z = −1.67) | 0.80 |
| rat | 11 | 7 | ×0.72 per carbon (z = −0.56) | ×0.03 (z = −1.25) | 1.17 |
| mouse | 10 | 3 | ×0.60 per carbon (z = −1.60) | ×0.12 (z = −1.74) | 0.70 |

*(The rat improved from 10 rows and σ = 1.75 when Ohmori's Fig. 2 was
digitised — §10a — taking it from one direct renal clearance to four.)*

**The sign is the same in all three species and the point estimates are large
— but no linear coefficient reaches significance** (largest |z| = 1.74). With
10–19 rows, most of them bounds, that is what should be expected. The honest
statement is that the fit is *consistent with* a strong chain-length effect
and does not *establish* one; three species agreeing in sign happens by chance
one time in four.

### The shape is not linear, and the curvature is decided by the disputed rows

Adding nC² improves the human fit by ΔAIC = −7.08 and the rat by −7.43, and
in the quadratic model the coefficients do become significant (human nC
z = 2.71, nC² z = −2.60, sulfonate z = −4.18). The fitted vertices are
**maxima** at nC 5.3 (human), 6.9 (rat), 5.7 (mouse).

That is mechanically what should happen. T divides CL_R by f_u, and f_u has a
*minimum* near C7 (§7b), so 1/f_u has a maximum there and **T inherits an
inverted U** — the mid-short chains, C5–C7, are the compounds the tubule
returns least of.

**In the human it does not survive a sensitivity check.** The vertex sits
exactly on PFPeA and PFHxA — the two compounds whose tubular sign depends on
which binding dataset is used:

| fit | ΔAIC, quadratic vs linear |
|---|---|
| human, all 19 rows | **−7.08** quadratic justified |
| human, without PFPeA and PFHxA | **+7.78** quadratic rejected |
| **rat, all 11 rows** | **−12.04** quadratic justified |
| **rat, without the Fig. 2 additions** | **−6.10** **still justified** |

*(Compare ΔAIC within each row only; absolute AICs are not comparable across
different datasets.)*

**So the choice of human f_u source decides whether the human chain-length
relationship is monotone or humped** — the functional form of the model. The
binding question of §7b is upstream of the model structure, not a detail
beside it.

**But the rat says the shape is real.** Its curvature survives removing the
newly digitised rows (they strengthen it from −6.10 to −12.04 rather than
creating it), it rests on four *direct* renal clearances rather than contested
binding, and the rat has no four-way f_u dispute. **Same shape, very different
evidentiary standing — and the species without the binding problem is the one
that supports it.**

---

## 12. Questions for discussion

Each question carries the finding that provokes it, so it can be taken up
without reading the sections above.

### 12.1 On the species difference

**1.1 — If the two rodents are fast for different reasons, which one is the
model?** Against the human, the rat's difference in PFOA renal clearance is
carried by the tubule (~10²-fold) and the mouse's by the glomerulus (37.7× of
an 85× total). A PBPK model calibrated on the mouse tubule and applied to the
rat imports the wrong term. **Does our current practice of treating "rodent"
as one extrapolation source survive this?**

**1.2 — What does it mean that each species reabsorbs with a transporter the
other lacks?** Human OAT4 is not detected in rat or mouse kidney; rat Oatp1a1
is not detected in human. This is presence/absence, not a difference of
degree. **Is a scaling factor even the right instrument for a
presence/absence difference, or does it require a structural change to the
model?**

**1.3 — The rat sex difference is a change of sign, not of degree.** Female
rats *net-secrete* PFOA and PFNA; males reabsorb them. Established twice
independently — Oatp1a1 is 114× more abundant in males, and the measured
renal fractions put females above the secretion threshold and males below.
**If sex can flip the sign within one species, what is the defensible way to
use a single-sex rodent study for human extrapolation?**

**1.4 — Which species should we be arguing about?** No plasma unbound fraction
has ever been measured in any macaque, for any PFAS — yet cynomolgus PFOA has
serum, urine *and* faecal clearances at steady state (Butenhoff 2004), and
the monkey's PFOA half-life is the closest to the human. **One equilibrium
dialysis experiment would open the species we most need. Is that worth
commissioning or collaborating on?**

---

### 12.2 On the unbound fraction

**2.1 — Is a single f_u defensible at all any more?** OEHHA Table A6.4 assumes
f_u = 0.02 for every species and every compound, uncited. Measured human
values run from **0.0004 (PFHpA) to 0.3108 (PFPA) — a factor of 780.**
Measured human and rat PFOA differ only 1.22×, but mouse and human differ
5.8×. **Which of our standing documents still carry a constant f_u, and what
is the triage order for fixing them?**

**2.2 — Variable in which sense?** Three distinct meanings get conflated:

| sense | status | consequence |
|---|---|---|
| across **compounds** | settled, 780× range | largest single correction available |
| across **individuals** | measurable; Fischer 2024 predicts f_u from a person's own albumin and globulin, 2.5× across NHANES | affects susceptibility, not the central estimate |
| with **concentration** | real above a PFAS:protein ratio ≈0.1; environmental exposure sits near 0.0005 | matters for comparing in vitro datasets, not for human dosimetry |

**Which of these does a PHG actually need, and at which point in the
calculation does each enter?**

**2.3 — Which tail is protective?** A higher f_u means a larger filtered load
and *faster* elimination, so the individuals with the **lowest** f_u carry the
longest half-lives and the highest body burdens. **When we build an
f_u distribution rather than a point value, are we confident the health-
protective percentile is the low one — and is that how our current defaults
are oriented?**

**2.4 — Low albumin as a susceptibility factor.** Fischer finds f_u varies most
for compounds below 500 g/mol, driven by serum albumin, and notes that the
lowest albumin levels track diabetes and fatty liver disease. **Is
albumin-deficient status a PFAS susceptibility factor we should be naming
explicitly — and does it cut the way intuition suggests?**

**2.5 — Four datasets, four methods, and the sign of the answer depends on the
choice.** For human PFHxA, T is 3.47 on Smeltz, 2.48 on Fischer, 1.18 on
Argoul and 0.62 on Ryu — secretion on three, reabsorption on one. **Is there a
principled basis for choosing, or should short-chain conclusions simply be
reported as sign-undetermined?** (This review recommends Smeltz as the default
and Ryu for species ratios only; §7b of the report sets out why.)

---

### 12.3 On the transporters themselves

**3.1 — Does any measured transport survive contact with plasma protein?**
Every transport assay in this compilation was run **without albumin**. Tan
2024 showed that with 4% human serum albumin present and unbound concentration
held constant, unbound uptake clearance rose **3.5 to 122-fold at OAT1** — but
for drugs, not PFAS, and PFAS sit at or below the bottom of the f_u range
where the effect was largest. **If this holds for PFAS, every clearance in
the literature is a lower bound. How much of our mechanistic picture is
conditional on an experiment nobody has run?**

**3.2 — Abundance versus affinity.** The highest intrinsic clearance in the
entire compilation belongs to a transporter that two laboratories cannot
detect in human kidney. Affinity is conserved across species to within
1.2–1.8× where a transporter exists in both. **Are we over-weighting affinity
data because it is abundant, and under-weighting expression because it is
not?**

**3.3 — Secretion is not the minor term we treat it as.** Mouse PFHxA is net
*secretion* (T = 2.2–2.7), female rat PFOA and PFNA are net secretion, and on
three of four binding datasets human PFHxA is too. **Our framing is almost
entirely about reabsorption. Is that a property of PFAS or an artefact of
having studied PFOA first?**

**3.4 — What is the minimum transporter dataset for a defensible IVIVE?** The
kidney side of the OAT4 scaling factor exists and is replicated to the
reported precision by two laboratories (0.5 ± 0.2 pmol/mg). The in vitro side
does not exist, and it is one proteomics run on a cell line that already
exists. **Is this the kind of gap we can close ourselves, or fund?**

---

### 12.4 On how this meets the regulatory numbers

**4.1 — Two secondary-source errors with regulatory reach.** ATSDR Table 3-6
misreports rat PFDA total clearance by **forty-fold** (207 against the paper's
5.2 mL/kg·day) and male rat PFOA by 2.7×. OEHHA's own A6.4 percentages rest on
an uncited f_u. **What is our process for catching this class of error, given
both were found by arithmetic rather than by reading?**

**4.2 — A number that should stop being quoted.** The male rat's PFOA
reabsorption has been reported here as 93.7%, then 58%, then "at most half and
possibly none", then "at most a quarter" — each change driven by a defensible
source choice, not by new data. Three rat renal clearances straddle the
filtration boundary, and the rat GFR itself is disputed 1.82×. **Is the
honest move to publish the range, or to stop reporting a rat percentage at
all?** The *direction* relative to the human has never wavered.

**4.3 — Where the renal framework stops.** 58–61% of human elimination of
PFNA, PFDA, PFUnDA and PFDoDA is faecal, measured. Human PFDoDA clears at 3.1×
its filtered load. **Above C11 the renal accounting does not describe the
compound. Do our long-chain assessments acknowledge that, and what is the
minimum mass-balance treatment that would?**

**4.4 — What would change our minds?** Stated in advance, the cheapest
experiments that would overturn conclusions in this review: a cynomolgus
plasma f_u (opens a third species); a PFAS uptake assay ± albumin (could make
every clearance a lower bound); OAT4 abundance in the expressing cell line
(decides whether OAT4 dominates human handling); and rat f_u by sex for PFOA
(the one input the sex-difference result assumes rather than measures).
**Are any of these within reach of a collaboration?**

---

### 12.5 For the QSAR work specifically

**5.1 — A QSAR for PFAS renal clearance cannot be a QSAR for albumin binding
alone.** Binding explains 61% of the log-variance in the human clearance
series and 83% in the mouse — but the tubular residual still spans 466- and
823-fold. **Does that settle whether the transporter layer belongs inside the
QSAR, or alongside it?**

**5.2 — Most of the endpoint data is censored.** 28 of 44 rows are one-sided
bounds. Fitting them as point values biases every coefficient; dropping them
loses two-thirds of the data. **Is a censored likelihood something our
modelling workflow can accommodate, or do we need to build that?**

**5.3 — The descriptors are collinear by construction.** In a homologous
series, molecular weight is linear in carbon number (r = 0.91 here), so a
model using both fits one variable twice. Independent information lives in the
head group, ether oxygens, branching and pKa — and only head group has more
than a couple of levels in our set. **Does this argue for a deliberately small
model, or for widening the compound set before modelling at all?**

**5.4 — Three laboratories now report a U-shaped f_u against molecular
weight**, with the minimum near 400–500 g/mol (measured vertex: PFHpA, C7,
f_u = 0.0004), and a mechanism — above ~8 perfluorinated carbons the molecule
is too rigid to enter the albumin pocket. **A monotonic chain-length
descriptor cannot represent that. What is the right functional form?**

---

## Sources for the central claims

Bibliographic details verified against PubMed. The full catalogue is
`db/source_catalogue.csv` (81 sources).

### In vivo clearance and binding — the tubular axis

| Source | Supplies |
|---|---|
| **Ryu S, Burchett W, Zhang S, Modaresi SMS, Agudelo Areiza J, Kaye E, Fischer FC, Slitt AL** (2024) *Toxics* **12**:253, [10.3390/toxics12040253](https://doi.org/10.3390/toxics12040253) | measured plasma f_u, 14 PFAS × human, rat, 2 mouse strains |
| **Abraham K** et al. (2024) *Environ Int* **193**:109047, [10.1016/j.envint.2024.109047](https://doi.org/10.1016/j.envint.2024.109047) | human urinary, total **and faecal** clearance + f_u + creatinine clearance, one subject, 15 PFAS |
| **Ohmori K, Kudo N, Katayama K, Kawashima Y** (2003) *Toxicology* **184**:135–140, [10.1016/s0300-483x(02)00573-5](https://doi.org/10.1016/s0300-483x(02)00573-5) | rat total and renal clearance, 4 PFCAs × both sexes; the T > 1 secretion test, in 2003 |
| **Kudo N** et al. (2001) *Chem Biol Interact* **134**:203–216 | rat urinary recovery by compound and sex — the renal-fraction tier |
| **Kudo N, Suzuki E, Katakura M, Ohmori K, Noshiro R, Kawashima Y** (2002) *Chem Biol Interact* **139**:301–316 | rat renal clearance, probenecid and castration controls |
| **Kim S-J** et al. (2016) *Drug Metab Dispos* **44**:1341–1348 | rat renal clearance and f_u by sex for PFHxS |
| **Argoul CML** et al. (2026) *Arch Toxicol* **100**:5151–5158, [10.1007/s00204-026-04508-7](https://doi.org/10.1007/s00204-026-04508-7) | human reabsorption table; mouse renal, total and urinary partition; substrate specificity human vs rat |
| **Argoul CML** et al. (2025) *Environ Toxicol Pharmacol* **120**:104879, [10.1016/j.etap.2025.104879](https://doi.org/10.1016/j.etap.2025.104879) | measured plasma f_u, human and mouse, two experiments |
| **OEHHA** (2024) *PFOA and PFOS Public Health Goals*, Appendix 6 Table A6.4 | the species reabsorption axis this report recomputes |
| **ATSDR** (2021) *Toxicological Profile for Perfluoroalkyls*, Table 3-6 | secondary compilation; three values corrected in §9 |

### Transporters — localisation, expression and kinetics

| Source | Supplies |
|---|---|
| **Niu S** et al. (2023) *Environ Health Perspect* **131**:076002, [10.1289/EHP11885](https://doi.org/10.1289/EHP11885) | renal localisation, both species |
| **Thakur A** et al. (2024) *Clin Pharmacol Ther* **116**:235–246, [10.1002/cpt.3277](https://doi.org/10.1002/cpt.3277) | expression, 3 species, both sexes; the 114× Oatp1a1 sex difference |
| **Nakagawa H** et al. (2008) *Basic Clin Pharmacol Toxicol* **103**:1–8, [10.1111/j.1742-7843.2007.00155.x](https://doi.org/10.1111/j.1742-7843.2007.00155.x) | matched human/rat OAT kinetics in one expression system |
| **Louisse J** et al. (2023) *Arch Toxicol* **97**:685–696, [10.1007/s00204-022-03428-6](https://doi.org/10.1007/s00204-022-03428-6) | OAT4 kinetics |
| **Louisse J** et al. (2024) *Toxicology* **509**:153961, [10.1016/j.tox.2024.153961](https://doi.org/10.1016/j.tox.2024.153961) | OAT1/OAT3 kinetics |
| **Niu S** et al. (2026) *Environ Health (Wash)* **4**:2150–2157, [10.1021/envhealth.6c00088](https://doi.org/10.1021/envhealth.6c00088) | efflux kinetics, human cortex expression |
| **Weaver YM** et al. (2010) *Toxicol Sci* **113**:305–314, [10.1093/toxsci/kfp275](https://doi.org/10.1093/toxsci/kfp275) | rat renal kinetics |
| **Uchida Y** et al. (2016) *J Pharm Sci* **105**:1011–1016, [10.1002/jps.24645](https://doi.org/10.1002/jps.24645) | rat kidney, membrane-resolved |
| **Kumar AR** et al. (2021) *Drug Metab Dispos* **49**:470–478, [10.1124/dmd.121.000367](https://doi.org/10.1124/dmd.121.000367) | IVIVE scaling factors; recommends REF over RAF |
| **Tan SPF** et al. (2024) *Mol Pharm* **21**:4603–4617, [10.1021/acs.molpharmaceut.4c00504](https://doi.org/10.1021/acs.molpharmaceut.4c00504) | albumin-mediated uptake — the §8 caveat |

---

## Reproducing this

Every table in this report is rebuilt from the primary data by an idempotent
script; none of the record is edited by hand.

```bash
cd tk_review
python3 -I scripts/build_transporter_working_table.py     # the 58-row kinetic record
python3 -I scripts/stage6_consolidate_tubular.py          # the consolidated T table
python3 -I scripts/stage6_abraham_partition.py            # the human mass balance
python3 -I scripts/stage6_chain_length_bridge.py          # the binding/tubule decomposition
python3 -I scripts/stage6_figure_T_vs_chain.py            # the figure
```

Outputs land in `db/`; `db/tubular_factor_consolidated.csv` is the single
authoritative verdict table and `db/tubular_factor_evidence.csv` carries every
underlying source row, including the three superseded values from §9.
