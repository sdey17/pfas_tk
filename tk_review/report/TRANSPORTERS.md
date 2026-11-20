# PFAS renal transporters: what the measurements support

Companion to `ALBUMIN_BINDING.md`. Working files are in `db/`; the stage-by-stage
record of how this was assembled, including the defects found along the way, is
`TRANSPORTERS_ROADMAP.md`.

---

## What this establishes, and what it does not

**Five things the measurements support.**

1. **The reabsorptive machinery is not the same genes in human and rat.** Human
   OAT4 is not detected in rat kidney; rat Oatp1a1 is not detected in human.
   Each species reabsorbs PFAS with a transporter the other does not have.
2. **Transporter affinity is conserved** where a transporter exists in both —
   1.2 to 1.8-fold on OAT1 and OAT3, measured in one expression system.
3. **Substrate specificity is conserved**: the same molecular-weight threshold,
   about 500 g/mol, governs transporter interaction in both species.
4. **The rat sex difference has a mechanism, and it is a sign change.** Rat
   Oatp1a1 is 114-fold more abundant in males; male rats reabsorb filtered PFOA
   while females are in net tubular secretion.
5. **Abundance gates relevance.** The highest intrinsic clearance in this whole
   compilation belongs to a transporter that two laboratories cannot detect in
   human kidney.
6. **The human–rat contrast in renal reabsorption is real and larger than the
   regulatory figure.** Measured unbound fractions put human reabsorption of
   filtered PFOA at **99.0%** and the male rat at **58%**, not the 99.8% and
   93.7% computed from an assumed binding value.

**Three things it does not.**

1. Whether human OAT4 reabsorbs more effectively *in situ* than rat Oatp1a1.
   Per mg of transporter protein they are within 1.3-fold; per kidney they
   cannot be compared, because neither has a scaling factor.
2. ~~The species contrast in renal reabsorption.~~ **Resolved 2026-11-19**,
   and it moved to the supported list as item 6 — but the widely quoted
   figures were wrong. See §6.
3. Whether any of the measured clearances survive contact with plasma protein.
   Every assay in this compilation was run without it.

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
from 6.3% in males to **net secretion** in females — 37 to 48-fold, with a sign
change. What differs:

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

---

## 6. The species contrast in reabsorption is not established

The figures everyone quotes — human 99.8% of filtered PFOA reabsorbed, male rat
93.7% — come from OEHHA's PHG Table A6.4, computed as f_u × GFR − CL_R with
**f_u assumed to be 0.02 in every species**. The PHG states it in one clause,
with no citation and no species justification, and draws a headline conclusion
from it: *"this mechanism would explain the exceptionally long human T½
relative to animals."*

**The assumption is not even constant across compounds.** Measured human plasma
f_u spans 52-fold across nine PFAS, from 0.08% (PFHxS) to 4.2% (GenX). The
assumed 2% is about right for the short chains (PFHxA 2.0%, PFBS 1.7%) and **12
to 25-fold too high for the long chains the PHG is about.**

**All three species have now been measured.** Ryu 2024 (*Toxics* 12:253,
presaturation equilibrium dialysis, 14 PFAS) reports plasma f_u in human, rat
and two mouse strains; Argoul 2025 reports human and mouse independently.
Recomputing OEHHA's table with measured values rather than 0.02:

| | f_u | filtered | CL_R | **reabsorbed** | OEHHA printed |
|---|---|---|---|---|---|
| human | 0.00245 | 6.30 | 0.06 | **99.0%** | 99.88% |
| **male rat** | 0.00298 | 42.9 | 18.2 | **57.6%** | **93.68%** |
| male mouse | 0.0142 | 237 | 10 | **95.8%** | 97.01% |

**The male rat does not reabsorb 94% of filtered PFOA. It reabsorbs about
58%.** The regulatory figure is wrong by a wide margin, because the assumed
f_u is eight-fold too high for the rat.

**But the assumption is approximately right in *structure* for this pair**, and
that is why the species conclusion survives. Rat and human plasma f_u for PFOA
differ only **1.22-fold** (0.00298 against 0.00245); Ryu's own equivalence
testing puts human and rat within two-fold, mean absolute fold difference 1.32.
So the decomposition of the observed 303-fold clearance difference holds:

| | glomerular | tubular |
|---|---|---|
| f_u assumed equal (OEHHA) | 5.6× | 54× |
| **f_u measured** | **6.8×** | **45×** |

**The mouse is where the assumption fails structurally, not the rat.** Mouse
f_u is 5.8-fold the human (Ryu) and 6.7-fold (Argoul, independently) — so for
mouse-versus-human, holding f_u equal shifts the difference from glomerular to
tubular and reverses the reading.

**The corrected position is narrower and firmer than either earlier one.** The
human–rat contrast in renal reabsorption **is** established, and it is larger
in character than the PHG states: the human reabsorbs essentially all filtered
PFOA, the male rat a little over half. **The percentages OEHHA prints are not
reliable. The attribution of the difference to the tubule is.**

### All fourteen compounds, and the split inverts between the rodents

Ryu measures fourteen PFAS, so the binding layer can be run whole rather than
on PFOA alone. Two results follow.

**First, the rat approximation generalises.** Rat/human f_u has a geometric
mean of **1.20×** across all fourteen and lies within two-fold for **13 of
14**. Mouse/human is 3.7–4.1× and within two-fold for only 3 of 14. OEHHA's
constant-f_u assumption is therefore a sound approximation for the rat and a
bad one for the mouse, **across the whole class** — not a lucky break on PFOA.
The mouse/human error is also chain-length-dependent, rising from 3.4× at C4
to a 4.8–7.3× plateau over C5–C8 and falling to ~2× by C9–C12: largest exactly
over the compounds of greatest regulatory interest.

**Second, and this is the substantive finding, the two rodents are fast for
different reasons.** For PFOA:

| | observed CL_R ratio | glomerular | tubular |
|---|---|---|---|
| human vs **rat** | 303× | 6.8× | **~10²** |
| human vs **mouse** | 85× | **37.7×** | 2.3× |

The rat tubular factor is an order of magnitude rather than a figure. Three rat
studies report different PFOA renal clearances and four human sources report
clearances spanning 2.9-fold; across every combination the ratio runs **22× to
161×**. **The dominance holds under all of them — the smallest value is still
three-fold above the 6.8× glomerular factor — but the number is not
determined**, and quoting a single figure for it would misrepresent the
evidence.

The mouse differs from the human in *both* glomerular factors — f_u 5.8× higher
and GFR 6.5× higher, and 5.8 × 6.5 = 37.7 — so the glomerulus accounts for
nearly the whole difference. The rat differs in *one*, its f_u matching the
human, so the tubule must carry 45× of the 303×. **"Rodents clear PFOA faster
than humans" is one sentence covering two mechanisms**, and a PBPK model
calibrated on the mouse tubule and applied to the rat imports the wrong term.

Across the human–mouse pairs the glomerular factor dominates six of seven, and
for four of the six long chains the tubular factor is *below one* — the tubule
partly cancels the mouse's faster clearance rather than amplifying it.

**Coverage, stated plainly.** The f_u layer is complete (14 × 4 matrices, one
laboratory). The reabsorbed-fraction layer needs renal clearance, which exists
for 7 of 14 in the human and 7 of 14 in the mouse. The rat began at 1 of 14 and
is now at **6 of 14**, by a route that needs no new experiment: since
CL_R ≤ CL_total always, a rat total clearance gives an upper bound on T, and
**where that bound falls below 1 net reabsorption holds whatever the renal
share.** It does for PFBA (≥84% reabsorbed), PFNA (≥77%), PFHxS (≥18%) and PFOS
(≥17%) on the loosest bound available for each — and PFHxS has a direct renal
clearance too (Kim 2016), giving 93%. Four compounds stay uninformative because
their total clearance exceeds the filtered load, which permits either sign.

**One correction this forces.** The male rat's PFOA reabsorption is not 58%.
Three rat studies give renal clearances of 18.2, 21.2 and 46.1 mL/kg·day
against a filtration clearance of 42.9 — **so the sign of the rat tubular term
depends on which primary source is taken.** Kudo 2002's value, the highest and
the only one internally cross-checked against its own total clearance, puts the
male rat marginally into net *secretion*. OEHHA's choice (Kemper 2003,
unpublished) is the most reabsorptive of the three. The defensible statement is
that **the male rat tubule returns at most half the filtered load and possibly
none** — which still contrasts sharply with the human (T ≈ 0.007–0.010) and the
mouse (T = 0.015), and so leaves the species conclusion intact while removing
the printed percentage.

### The same rule applied to every species

One two-tier rule — direct renal clearance where a study in that species
reports it, total clearance as an upper bound otherwise — across all three
species, with the no-contamination constraint enforced by an assertion in the
build rather than by care. Every row's f_u, GFR and clearance come from one
species. Strain is tracked separately and reported, not blocked: Wistar against
Wistar-Han is flagged as a substrain difference, anything else as a strain
mismatch.

| species | with a clearance | conclusive | of which direct |
|---|---|---|---|
| human | **13 / 14** | **12 / 14** | 2 |
| rat | 9 / 14 | 4 / 14 | 1 |
| mouse | 7 / 14 | **7 / 14** | 7 |

The human jump is almost entirely one source — **Abraham 2024**, total
clearance for fifteen PFAS in a single 67-year-old male volunteer, from a
two-compartment fit to 450 days of plasma decay. **n = 1**, and eleven of the
twelve human conclusions rest on it. Worth stating plainly rather than
presenting 13/14 as settled coverage.

It does carry one independent confirmation, though: Abraham's human PFOA
clearance gives T = 0.00662 against Argoul's 0.00667 — two unconnected studies,
one a single-volunteer decay curve and one a population estimate, agreeing to
three significant figures on a quantity spanning three orders of magnitude
across the table.

**Two results from the wider sweep.** Mouse **PFHxA is net secretion**
(T = 2.23), the only direct secretory verdict anywhere in the table — the
mouse's renal clearance of PFHxA exceeds its filtered load. And human
**PFDoDA** has a total clearance 3.1× its filtered load, which means its
elimination cannot be mostly glomerular: for a C12 carboxylate the likely answer
is that most of it is not renal at all. **T is a renal accounting and it does
not describe the longest chains** — a limit of the framework, not a gap in the
data. Anything above C11 belongs in a mass-balance model with a faecal term.

**The monkey is excluded, and the reason points at the next experiment.** The
project holds 22 cynomolgus rows with 14 clearances, but **no plasma unbound
fraction has ever been measured in cynomolgus plasma for any PFAS.** The only
macaque GFR available is for the *Japanese* macaque, a different species, so
borrowing it is refused. One equilibrium-dialysis experiment on cynomolgus
plasma would open a third species for the entire framework — and the monkey is
the species whose PFOA half-life is closest to the human. **That is the
highest-value single experiment this strand has identified.**

### Ohmori 2003, read from the paper

That paper turned out to settle more than expected, and to correct three values
this project had been carrying.

**Two ATSDR values are wrong, one by forty-fold.** ATSDR Table 3-6 attributes a
male rat PFOA total clearance of 135 mL/kg·day to Ohmori; the paper says
**50.5**. It attributes PFDA clearances of 207 and 140; the paper says **5.2
and 5.3**. Five of the eight values agree. The PFDA error was the sole reason
rat PFDA sat as uninformative, and at the true value the compound closes as
certain reabsorption.

**A third evidence tier, and it needs no renal clearance measurement.** With
f_ren the fraction of elimination that is renal, T = f_ren · CL_tot/(f_u · GFR),
so **net secretion is established as soon as f_ren exceeds (f_u · GFR)/CL_tot**.
Where a study in the same species measures urinary recovery above that
threshold, the verdict follows. Kudo 2001 measured it in Wistar rats from
Ohmori's own laboratory:

| | threshold | measured | verdict |
|---|---|---|---|
| PFHpA male | 4.1% | **92%** | **net secretion**, T ≈ 23 |
| PFOA male | 85.0% | 55% | undetermined |
| **PFOA female** | 1.9% | **80%** | **net secretion**, T ≈ 42 |
| **PFNA female** | 28.6% | **51%** | **net secretion**, T ≈ 1.8 |
| PFNA male | impossible | 2% | reabsorption certain |
| PFDA both | impossible | 0.2% | reabsorption certain |

**So the rat sex difference is a change of sign, not of degree.** The female rat
*net-secretes* PFOA and PFNA; the male rat does not — male PFOA needs an 85%
renal fraction against a measured 55%, and male PFNA cannot secrete at all.
**The 70-fold sex difference in rat PFOA half-life is the difference between a
secreting kidney and a reabsorbing one.** §14 of the working record reached the
same conclusion from Oatp1a1 expression (114×, with a sign flip); this reaches
it from in vivo clearance alone, with no transporter data in it. It survives a
ten-fold error in the pooled f_u for PFOA, though only a 1.8-fold error for
PFNA.

**Ohmori reached this framework in 2003.** The paper concludes tubular secretion
for some PFCAs "because CL_R of free PFCA was higher than the glomerular
filtration rate" — which is T > 1 — with probenecid inhibition as support. Its
own binding data cannot be used (">98% bound, data not shown", ultrafiltration
at 300 µM, sixty-fold Ryu's concentration and near albumin saturation), so the
test above uses Ryu's measured f_u.

**And the rat GFR is not agreed.** Ohmori uses **5.5 mL/min/kg** citing Shim
1984; OEHHA uses **10**. A 1.82× fork, and on Ohmori's figure the male rat
reabsorbs **at most a quarter** of filtered PFOA rather than at most half. That
is the third successive correction to loosen this number, which is itself the
result: **the rat reabsorption figure has no settled value, only a settled
direction relative to the human.** The human tubule reabsorbs over 99% on every
combination of sources.

The remaining gaps are lab measurements, not search problems: cynomolgus plasma
f_u, rat f_u by sex for PFOA, and the eight renal clearances in Ohmori's Fig. 2,
which are plotted as a bar chart with no numbers in the text.

Ryu's plasma is equal-sex-pooled Wistar-Han; the rat clearance it is paired
with is male Sprague-Dawley. One rat study does measure rat f_u by sex — Kim
2016, for PFHxS — and the answer is favourable: **0.076% male against 0.069%
female, a 1.10-fold difference, against a 7-fold sex difference in that
compound's clearance.** So for PFHxS the sex difference is not a binding
difference and pooling is defensible. PFOA, where the rat sex difference is
largest, is still unmeasured by sex. The mouse pairing is strain-matched (C57BL/6) but not
sex-matched. Mouse f_u for PFOS differs **5-fold between CD-1 and C57BL/6**, so
"mouse f_u" is not one number for that compound.

**On the two laboratories:** the apparent consistent offset was an artefact of
looking at PFOA. Across fifteen overlapping pairs there is no net offset
(geometric mean 0.90×) and the spread is thirteen-fold, concentrated at the
long chains. A test that trusts neither laboratory settles it: within the
carboxylate series f_u must fall monotonically with chain length, and **Ryu's
series is monotone in all three species while Argoul's inverts at C9→C10 and
again at C10→C11.** Argoul's PFUdA value is forty times Ryu's and nearly as
high as Argoul's own C6. The figures above use Ryu's f_u, and that choice is
argued rather than assumed.

---

## 7. A caveat that may matter more than the scaling

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

## 8. What would change these conclusions

| Measurement | What it decides |
|---|---|
| **OAT4 abundance in the expressing cells** | whether OAT4 really dominates human renal handling, or whether scaled OAT1 overtakes it (§4) |
| **One PFAS uptake assay ± albumin** | whether every clearance here is a lower bound (§7) |
| Rat Oatp1a1 REF | whether the reabsorptive sides can be compared per kidney at all (§5) |
| Yang 2010's own V_max | adds 3 rows to tier A; changes no conclusion |

The first three are each a single experiment, and each decides something this
document currently has to leave open.

---

## Sources for the central claims

Bibliographic details verified against PubMed.

| Source | Supplies |
|---|---|
| **Niu S** et al. (2023) *Environ Health Perspect* **131**:076002, [10.1289/EHP11885](https://doi.org/10.1289/EHP11885) | renal localisation, both species |
| **Thakur A** et al. (2024) *Clin Pharmacol Ther* **116**:235–246, [10.1002/cpt.3277](https://doi.org/10.1002/cpt.3277) | expression, 3 species, both sexes |
| **Argoul CML** et al. (2026) *Arch Toxicol* **100**:5151–5158, [10.1007/s00204-026-04508-7](https://doi.org/10.1007/s00204-026-04508-7) | substrate specificity, human vs rat |
| **Argoul CML** et al. (2025) *Environ Toxicol Pharmacol* **120**:104879, [10.1016/j.etap.2025.104879](https://doi.org/10.1016/j.etap.2025.104879) | measured plasma f_u, human and mouse |
| **Nakagawa H** et al. (2008) *Basic Clin Pharmacol Toxicol* **103**:1–8, [10.1111/j.1742-7843.2007.00155.x](https://doi.org/10.1111/j.1742-7843.2007.00155.x) | matched human/rat OAT kinetics |
| **Louisse J** et al. (2023) *Arch Toxicol* **97**:685–696, [10.1007/s00204-022-03428-6](https://doi.org/10.1007/s00204-022-03428-6) | OAT4 kinetics |
| **Louisse J** et al. (2024) *Toxicology* **509**:153961, [10.1016/j.tox.2024.153961](https://doi.org/10.1016/j.tox.2024.153961) | OAT1/OAT3 kinetics |
| **Niu S** et al. (2026) *Environ Health (Wash)* **4**:2150–2157, [10.1021/envhealth.6c00088](https://doi.org/10.1021/envhealth.6c00088) | efflux kinetics, human cortex expression |
| **Weaver YM** et al. (2010) *Toxicol Sci* **113**:305–314, [10.1093/toxsci/kfp275](https://doi.org/10.1093/toxsci/kfp275) | rat renal kinetics |
| **Tan SPF** et al. (2024) *Mol Pharm* **21**:4603–4617, [10.1021/acs.molpharmaceut.4c00504](https://doi.org/10.1021/acs.molpharmaceut.4c00504) | albumin-mediated uptake |
| **Kumar AR** et al. (2021) *Drug Metab Dispos* **49**:470–478, [10.1124/dmd.121.000367](https://doi.org/10.1124/dmd.121.000367) | IVIVE scaling factors |
| **Uchida Y** et al. (2016) *J Pharm Sci* **105**:1011–1016, [10.1002/jps.24645](https://doi.org/10.1002/jps.24645) | rat kidney, membrane-resolved |
| **OEHHA** (2024) *PFOA and PFOS Public Health Goals*, Appendix 6 | species reabsorption, in vivo clearances |

One citation correction worth recording: Argoul 2026 cites "Argoul et al. 2025"
for its unbound fractions, but its reference list contains that year only as a
duplicate of a different paper — the LC-MS/MS method paper, same title, same
DOI as its 2024 entry. The f_u source is the *Environ Toxicol Pharmacol* paper
listed above.
