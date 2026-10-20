# PFAS toxicokinetics — executive summary

**Three questions.** Why do PFAS serum half-lives differ so much between humans,
rats and mice? What did each reported half-life assume? Is exposure level related
to half-life or to volume of distribution?

**One equation underneath all of it.**

```
t½ = ln2 · Vd / CL
```

Half-life is not a property. It is a ratio of how much of the body a chemical
occupies to how fast the body removes it. Any two terms fix the third. Nearly
every disagreement in this literature turns out to be about **which two were
measured and which one was assumed** — and in the most consequential case, the
assumed one is the one everybody quotes.

---

## The seven findings

### 1. The species difference is a clearance difference, and it is a female one

Six datasets measure *both* terms in *both* sexes inside single experiments —
three species, two compounds, four laboratories, no cross-study comparison
needed:

| | Vd male/female | clearance female/male |
|---|---|---|
| **span across all six** | **1.33 – 2.18×** | **0.78 – 44.3×** |
| | a 1.6× spread | a 56× spread |

The volume of distribution sex ratio is male-higher every single time and never
leaves a narrow band. The clearance ratio ranges over nearly two orders of
magnitude **and changes sign between species**. Argoul 2026 reaches the same
conclusion inside one experiment: across 11 PFAS dosed as a single cocktail,
clearance spans **5,254×** while Vss spans **7.7×**. (That study's *absolute*
clearances run ~2.5× below the one-compound-at-a-time literature, but a shared
scaling cancels from a ratio, so this contrast is unaffected — see REPORT.md
§11.3a.)

Taking both limbs from primary sources, the gap is concentrated in females —
female rat ÷ female mouse PFOA clearance is **373–496×**, against **7.0×** in
males. Tatum-Gibbs 2011 replicates the whole structure in a second compound with
strains matched (rat sex ratio 21.9×, mouse ~1×).

### 2. One mechanistic axis explains the ordering: fractional renal reabsorption

```
CL_renal = fu · GFR · (1 − FR)
```

Human 99.94% → mouse 97/95.2% → male rat 93.7% → macaque 91.2/81.2% → dog
59/52% → **female rat and rabbit in net secretion**. It reproduces rodent
half-lives within 1.7× over a 37× span with no free parameters, and an
unrelated dataset recomputes mouse PFOA at **95.9%**, inside the 95.2–97.0%
already in use.

**But the axis has two factors, and this review previously treated it as one.**
Humans have the longest half-life yet reabsorb *less* in absolute terms
(51 mL/d/kg) than male rats (270) or mice (318–324). Of the 333× male-mouse-to-
human renal clearance gap, the escape fraction carries **50×** and the six-fold
lower human GFR carries **6.5×** — 68% and 32% on a log scale. Reabsorption is
the larger term, which is why the single-axis framing works; but a third of the
difference is filtration rate, which no transporter story explains.

**The transporter is identified by elimination, and it is not the same protein
in humans.** A candidate must be both sex-divergent and able to carry PFOA. In
the rat only **Oatp1a1** is both (23× male-predominant, androgen-induced,
transports C8–C10); Oat2 is strongly sex-divergent but carries no PFOA, and
Oat1/Oat3 carry it but are not sex-divergent. In humans, **OATP1A2 — the closest
orthologue of rat Oatp1a1 — does not transport PFOA at all**; human apical
reabsorption runs through OAT4 and URAT1, neither androgen-regulated. The
reabsorbed *fraction* transfers across species. The *mechanism* does not.

### 3. Exposure is related to neither half-life nor volume of distribution

Vd slopes against dose run −0.21 to +0.25 with inconsistent sign. Four
independent dose slopes for half-life cluster at **+0.08 to +0.12** against a
ceiling of 1. The between-person association does not survive age adjustment.

One real exception: the **female rat**, slope +0.21 over a 3,200× dose range
(saturable secretion); the male rat is −0.02.

### 4. The human volume of distribution that regulation rests on is an assumption

Nine of 42 adopted regulatory clearance factors are computed from Thompson
2010's 170 mL/kg, and EPA's own appendix shows ten human studies *assigning* it.
Its supplementary table is headed "Input data for the **calibration** of the Vd
parameter", with a final column headed **"calculated Vd"** — and the published
values reproduce to three figures from intake, serum and an assumed elimination
rate taken from a study in *the same two communities*.

The consequence is sharper than "circular". Forming `CL = ln2·Vd/t½` makes the
half-life **cancel exactly**, leaving `CL = Dose/Serum = 0.132–0.138 mL/kg-day`.
But Vd does not cancel: it is proportional to whatever half-life is assumed, from
**168 mL/kg at 2.3 y to 277 at 3.8 y**. Adopting "170" alongside a different
half-life silently contradicts the data it came from. The paper's own corrigendum
is that error, corrected by a factor of 0.6 — exactly the ratio of the two rate
constants.

Meanwhile every direct measurement (74, 121, 113–199 mL/kg) falls *below* the
assigned range, and the one population fit (430) sits far above it.

### 5. Five sixths of the field is empty

39 chemicals × 11 species × 4 parameters = 1,716 cells. **1,466 have no data at
all (85%)**; a further 94 rest on a single study. Even PFOS, the best-covered
compound in the world, fills 29 of 44. By species: human 71, rat 61, mouse 45,
monkey 28 — then **dog 2**. Mining four more compilations moved coverage from
143 to 156 cells, so the hole is real, not a search artefact.

### 6. Humans run two near-complete reabsorption loops, not one

The renal axis in finding 2 leaves a residual: predicted human half-life comes
out **4.3× short** of observation. The missing term is enterohepatic. Harada 2007
sampled serum and bile from four gallstone-surgery patients and measured a
**biliary resorption rate of 0.97**. Bile PFOS (27.9 ng/mL) actually *exceeds*
serum (23.2), so bile is a real excretion route — but 97% of what is secreted is
reabsorbed from the gut.

| loop | fraction reabsorbed | source |
|---|---|---|
| renal (PFOA) | **99.94%** | Han 2012 Table 4 |
| biliary (PFOS) | **97%** | Harada 2007, n = 4 |

Both are near-unity, both lengthen the human half-life, and **both are
interruptible** — which is why bile-acid sequestrants work. Delaere 2025:
treated PFOS half-life **1.2 y** (n = 19) against **7.3 y** under observation
(n = 9), a 6.1× difference; PFHxS 2.5 y against 9.4 y. Genuis 2010 and Møller
2024 agree in direction.

It also reconciles the faecal-route dispute without either side being wrong. At
97% resorption, **gross** biliary flux is large while **net** faecal elimination
is small — so Andersson, measuring faeces under ongoing intake, saw a large
signal, and Abraham, following a labelled bolus, saw almost nothing leave that
way. Different quantities, one reconciling constant.

The caveat is that n = 4, and the number is now load-bearing.

### 7. Half-life is the wrong endpoint for a structure-activity model

Because `t½ = ln2·Vd/CL`, half-life carries a glomerular filtration term that
differs **6.5×** between mouse and human with no change in chemistry. A
structural model fitted to it is being asked to absorb body size. The
replacement is the dimensionless **renal handling ratio**

```
R = CL_renal / (fu · GFR)
```

which divides out GFR and the unbound fraction and leaves the transport step.
`R < 1` is net reabsorption, `R > 1` net secretion, and `R = 1 − FR`, so it is
finding 2's axis on a scale that neither saturates near 1 nor runs unboundedly
negative.

On the nine compounds where both terms were measured in one experiment (Argoul
2026, male mouse), R spans **875×**, and three results constrain any QSAR:
chain length alone does **not** order it (within the carboxylates the series is
non-monotonic, and PFHxA crosses into net secretion); the head group carries
~3× at matched chain length (PFNA vs PFOS, both 8 fluorinated carbons); and
ether oxygens move it **1.4 log units** at constant chain length, non-monotonically
(PFHxA secreted, GenX reabsorbed, PFO2OA secreted).

The same construction sharpens finding 4's successor problem. R is linear in
`1/fu`, so the three human free fractions in the literature give three different
answers for human PFOA — log10 R of −3.93 (fu 0.10), −3.23 (fu 0.02, assumed),
−1.72 (fu 0.00061, measured) — and the mouse-to-human gap in the *transport*
step is correspondingly **354×, 71× or 2.2×**. Under the measured value the
species difference is almost entirely **binding**, which would mean the
transporter literature is explaining a quantity that barely differs between the
species. Stated as a sensitivity, not a result, because the mouse and human
values come from different methods at different ligand:protein ratios — the
exact artefact diagnosed above. One single-method `fu` panel would settle it.

---

## What this corrects in the published record

| | correction |
|---|---|
| **EPA** | the rat Oatp1a1 male/female ratio quoted as "2.5-fold" is a *different transporter's* number (OAT-K); the primary value is **23×** |
| **OEHHA** | its adaptation of the source reabsorption table altered two values (human 99.94→99.8, male rat 93.7→93.2) |
| **Provenance** | the `fu = 0.02` assumption behind the axis is Han 2012's, not OEHHA's, and Han calls its own values "rough estimates" |
| **Andersson 2025** | transposes Thompson's PFOA and PFOS volumes (it is PFOA 170, PFOS 230) — and an earlier claim in this repository that Zhang 2013 did so was wrong |
| **Cheng 2006** | its abstract calls renal Oatp1a1 "female-predominant", then reports androgens *increase* it and concludes androgens are the exclusive cause — self-contradictory, and against Cheng 2005 |
| **This review** | its own CPHEA refits compress sex ratios ~13× (see below) |

**The self-correction matters most.** The terminal-slope fits in
`db/cphea_fitted_halflives.csv` select the window with the best adjusted R²,
which on a biphasic curve is the shallow tail. Against a published value for the
same experiment the fit returns **1.44 d vs 0.08 d** for the female rat and
compresses a **71× sex ratio to 5.6×**. The bias runs one way, so every
conclusion above survives and several strengthen — but that column must not be
read as comparable to published half-lives. It is flagged, not deleted, because
it remains useful for within-curve comparison.

---

## The one thing still genuinely unresolved

Whether reabsorptive transport saturates at real human exposures depends on which
parameter family you believe. There are now **six independent in vitro Km values
for PFOA against human transporters** — 47–310 µM, from two laboratories,
agreeing within a factor of 7 — against a PBPK transport-affinity constant of
**0.133 µM**, which is **354–2,336× lower than anything ever measured in a
cell**. The in vitro values put human serum far below half-saturation; the fitted
KT puts three of four human populations at or above it, which would make
clearance dose-dependent in contaminated communities and mean a single clearance
factor cannot transfer between exposure settings.

The weight of evidence has moved decisively onto the in vitro side: the KT values
are fitted to plasma curves rather than measured, the mouse row has standard
errors exceeding its estimates, and the dose-response evidence agrees with the in
vitro answer. But **nobody has measured a human KT**, and that single number
would close it outright.

**One formerly-open conflict is now resolved.** The ~100× disagreement over
PFOA's plasma free fraction was not a contradiction between two measurements.
Han 2003's ">90% bound" is a *calculation* from Kd and albumin concentration, and
it is a floor that Fischer 2024's measured 0.00061 satisfies. The underlying Kd
gap is a ligand:protein ratio artefact — Han titrated at 1.7:1 to 60:1, Fischer at
≤0.004:1, and human serum sits at 10⁻⁵–10⁻³:1. The defect was in the
inheritance: PBPK models read ">90%" as "≈90%" and used a free fraction about
**164× too high**.

---

## What exists now

| | |
|---|---|
| `db/combined/` | 2,351 rows in four schemas; 20 chemicals, 8 species groups; every row carries its provenance and primary source. Also one Excel workbook. |
| `db/primary_2026/` | 24 per-paper extractions from the 25 full texts obtained during this work |
| `db/qsar/` | per-compound structure descriptors paired with the renal handling ratio (finding 7) |
| `scripts/` | 34 runnable scripts — every figure and table regenerates |
| `figures/` | 16 figures; `fig00_master.png` is the whole argument in one panel pair |
| `report/REPORT.md` | the full review, ~2,000 lines |
| `SUMMARY.md` | illustrated summary, 29 citations |
| `WANTED.md` | the papers still unobtainable, and what each would change |

Nothing is quoted from memory. Where a value came from another paper's citation
rather than the original, the row says so.
