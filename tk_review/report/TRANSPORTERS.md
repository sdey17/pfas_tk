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

Two laboratories agree with a consistent offset — Ryu runs 1.3–1.9× higher than
Argoul throughout, within the spread of equilibrium-dialysis methods, with the
species ordering and ratios preserved.

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
