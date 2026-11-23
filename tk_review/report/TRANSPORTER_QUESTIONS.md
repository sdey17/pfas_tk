# Transporters, unbound fraction and PFAS toxicokinetics

### Questions for discussion

Each question carries the finding that provokes it, so it can be taken up
without reading the full review. Evidence and sources are in
`TRANSPORTERS.md`; the working record, including everything withdrawn along
the way, is `TRANSPORTERS_ROADMAP.md`.

A note on the framing used throughout. Renal clearance factors exactly:

> **CL_R = f_u × GFR × T**, where **T = CL_R / (f_u × GFR)**
>
> T < 1 net reabsorption · T = 1 filtration alone · T > 1 net secretion

GFR is a species constant, so for any compound the whole argument reduces to
two terms: **how much is unbound**, and **what the tubule does with what is
filtered**.

---

## 1. On the species difference

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

## 2. On the unbound fraction

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

## 3. On the transporters themselves

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

## 4. On how this meets the regulatory numbers

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

## 5. For the QSAR work specifically

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
