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

Six, listed with evidence in §9. The two with regulatory consequence:

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
measured, reporting per compound the urinary, total and faecal clearances, the
unbound fraction, and the filtration clearance computed from *that subject's
own* creatinine clearance. T therefore needs no cross-study join and no assumed
GFR — unique in this strand. Subtracting the measured faecal term instead of
charging it to the kidney:

| | faecal share of total | bound before | bound after |
|---|---|---|---|
| PFNA | 58% | ≤0.014 | **≤0.0089** |
| PFDA | 61% | ≤0.093 | **≤0.0106** |
| PFUnDA | 60% | ≤0.266 | **≤0.0039** |

PFUnDA moves from the loosest long-chain bound to the tightest. **There is no
long-chain loss of reabsorption in the human:** the apparent rise in T toward
C12 was the faecal term, as predicted before the data were available.

**The short chains have no determined sign.** Three laboratories measured human
plasma f_u and disagree by up to 11-fold, worst at the short chains:

| | Abraham f_up | Ryu f_u | ratio | T (Abraham) | T (Ryu) | |
|---|---|---|---|---|---|---|
| PFPeA | 0.0440 | 0.0888 | 2.0 | **1.78** | **0.68** | **sign differs** |
| PFHxA | 0.0068 | 0.0380 | 5.6 | **3.48** | **0.48** | **sign differs** |
| PFBS | 0.0046 | 0.0365 | 7.9 | 0.129 | 0.013 | agree |

Human PFHxA either reabsorbs half its filtered load or secretes at three times
it, **decided entirely by a binding measurement.**

### The uncertainty is binding, not clearance

Abraham's urinary clearance of PFHxA is **47.1 mL/kg·day**; Argoul 2026 reports
**47.0** from an unrelated population study — agreement to better than 1%.
**The renal clearances agree across laboratories where they overlap; the
unbound fractions do not.** That is the single most useful thing this
comparison says about where effort should go.

Binding affinity rises with chain length, so within the carboxylate series f_u
must fall monotonically — a constraint on the data independent of any
laboratory's protocol:

| series | pairs | inversions |
|---|---|---|
| **Ryu f_u** | 8 | **0** |
| Argoul f_u | 3–6 | 1–2, at C9→C10 and C10→C11 |
| Abraham f_up | 7 | **4** |

Abraham's f_up falls to 0.0004 at C7 then rises **37-fold** to 0.0148 at C11,
which no binding mechanism supports. **The position adopted here is therefore
split, and deliberately so:** Abraham's measured *faecal* term (a mass balance,
not a binding assay) to tighten the long chains; Ryu's *f_u* for the filtered
load, because its series alone passes the monotonicity test; and the short
chains reported as sign-undetermined, because no test available here separates
0.48 from 3.48.

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

**6. `pdftotext` silently drops the micro sign**, and the loss is internally
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
| **A human f_u series that settles the short-chain disagreement** | the sign of T for C5–C7 in humans — currently 0.48 or 3.48 for PFHxA (§7) | *referred to OEHHA; and see below — two further datasets exist* |
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
