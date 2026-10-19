# Qin 2023 vs Armitage 2021 (IV-MBM EQP v2.0) on PFBA, PFHxA, PFHpA, PFOA, PFNA

**This version supersedes the earlier httk-based analysis in this file's
history.** Armitage's Supplementary Section S1, the Fragki main paper and the
Rosenmai supporting information are now in hand, and two of the earlier
conclusions were artefacts of substituting EPA's `httk` for the published
parameterisation. What changed, and why, is in §2.

Reproduce: `python3 tk_review/scripts/mbm_published_equations.py` and
`python3 tk_review/scripts/make_mbm_figure.py`.
Outputs: `db/mbm_published_equations.csv`, `figures/fig10_mbm_qin_vs_armitage.png`.

## 1. The two equations, as published

**Qin** (confirmed by his own text):

    C_free = C_nom x V_tot / ( V_w
                             + D_BSA/w  x V_protein,medium
                             + D_lip/w  x V_lipid,medium
                             + D_SP/w   x V_protein,cell
                             + D_lip/w  x V_lipid,cell )

`D_SP/w = D_non-specific`: *"D_cell/w was predicted by the MBM (Eq. 17) with
D_non-specific (Table 1) as a surrogate for partitioning to structural proteins
and D_lip/w (Table 2) for lipid."* There is **no cell-water term**: *"we
considered only proteins and lipids of the cell homogenate as the main sorptive
phase for D_cell/w, which is different from Henneberger et al. where the cells
were recognized as an entity containing water, protein and lipid and all of
them contributed."*

**Armitage** Eq. 2, with the SI S1 spLFERs:

    C_W = M_T / ( V_W + D_AW V_A + D_SaW V_Sa + D_SlW V_Sl
                      + D_DsW V_Ds + D_CW V_C + D_PlW A_Pl )

    log K_MW = 1.01 log Kow,N + 0.12
    log K_SaW = 1.08 log Kow,N - 0.70   (log Kow,N < 4.5)
              = 0.37 log Kow,N + 2.56   (log Kow,N >= 4.5)
    d_OW = d_PlW = 3.1,  d_MW = 1.0,  d_SaW = 0.0 for acids   (Table S1)

## 2. What the real SI changed, against the httk stand-in

| | httk `armitage.R` | Armitage 2021 SI | effect |
|---|---|---|---|
| cell composition | one generic cell, membrane lipid **0.04** | **cell-line specific**; the HepG2-derived AREc2 line is 0.00046 / **0.00414** / 0.0514 / 0.944 (Table S4) | **9.1×** |
| albumin spLFER | `0.71 logKow + 0.42` | piecewise `1.08 logKow − 0.70` / `0.37 logKow + 2.56` | 1.6× |
| ionisation offsets | d_OW = d_PlW = 3.5 | **3.1** | small |
| anionic phospholipids, lysosomes | present | absent from the published model | small |

The cell-composition term dominates, because Armitage's PFAS cell is almost
entirely membrane lipid and the hepatic line has ten times less of it than the
generic default. On Rosenmai PFOA at 100 µM: published 0.660% of dose in cells,
with httk's generic cell 6.007%, measured 0.330%.

Two independent checks that the composition table is the right one to use:
Armitage's PPARγ line (HEK293H — Qin's own cell line) is listed at **0.08**
protein by volume against Qin's measured **0.0814**; and Armitage's AREc2 at
10 ng/cell × 5.14% protein × 1.35 g/mL gives **694 pg protein per cell**
against the project spreadsheet's **675 pg/cell**.

## 3. Why the two models behave differently

**They put the PFAS in different places, and nearly agree by coincidence.**
At a common C_w = 1 µM:

| | nCF | Qin D_cell/w | Armitage D_cell/w | ratio | Qin D_BSA/w | Armitage D_SaW | ratio |
|---|---|---|---|---|---|---|---|
| PFBA | 3 | 2.4 | 1.2 | 1.99 | 520 | 138 | 3.76 |
| PFHxA | 5 | 25.4 | 10.1 | 2.52 | 4 952 | 6 503 | 0.76 |
| PFHpA | 6 | 83.2 | 56.4 | 1.47 | 12 044 | 26 501 | 0.45 |
| PFOA | 7 | 272 | 337 | 0.81 | 24 298 | 51 615 | 0.47 |
| PFNA | 8 | 892 | 2 036 | 0.44 | 45 849 | 100 531 | 0.46 |

The two cell ratios land within 0.44–2.5× of each other — but for unrelated
reasons:

- **Qin's cell is 92–95% protein** (the measured D_non-specific) and 5–8%
  lipid, with no water at all.
- **Armitage's cell is 79% water at PFBA**, flipping to **99.7% membrane
  lipid** by PFNA. Its bulk-protein term contributes 0.1–0.2%, because protein
  is treated as 3.5% of octanol and the anion then takes the −3.1 penalty.

So the one phase Qin says carries nearly all the PFAS is the one phase
Armitage's cell effectively ignores. The numerical agreement across these five
compounds is a coincidence of two wrong-looking routes crossing.

**The real divergence is now in albumin, not the cell.** Armitage's piecewise
spLFER rises 1.08 per log Kow below 4.5 — 0.83 per CF2 — then flattens to 0.285
per CF2 above it. Qin's measured affinity rises smoothly at ~0.49 per CF2. The
kink falls between PFBA and PFHxA, so Armitage reads 3.8× *low* at PFBA and
2.1–2.2× *high* from PFHpA up. Since albumin holds 19–98% of the dose in these
assays, that is what moves the predictions.

## 4. Dose dependence: the one test Armitage cannot pass

Fitting C_cell ∝ C_nom^b:

| series | measured | Qin | Armitage |
|---|---|---|---|
| Rosenmai PFBA | 0.85 | 1.15 | 1.00 |
| Rosenmai PFHxA | 0.87 | 1.58 | 1.00 |
| Rosenmai PFOA | 1.14 | 1.82 | 1.00 |
| Rosenmai PFNA | 1.31 | 1.81 | 1.00 |
| Fragki PFOA | 1.42 | 1.71 | 1.00 |
| Fragki PFNA | 1.42 | 1.66 | 1.00 |

Armitage is linear by construction and is pinned at exactly 1.00. The measured
exponent rises monotonically with chain length, 0.85 → 1.42, and Qin reproduces
that ordering while overshooting the magnitude by 0.3–0.7. This is the only
place where Qin's extra machinery buys something Armitage structurally cannot
deliver, and it is why the two models' error swaps sign with dose:

| | ≤25 µM | 100 µM |
|---|---|---|
| Armitage | 2.73× | 2.29× |
| Qin | **1.58×** | **8.95×** |

## 5. Performance

Fraction of dose in cells, n = 29:

| | Armitage / measured | Qin / measured |
|---|---|---|
| PFBA | 1.2–1.7× | 1.8–3.4× |
| PFHxA | 0.8–3.3× | 2.5–40.5× |
| PFHpA | 1.0–4.3× | 2.8–50.5× |
| PFOA | 1.9–5.8× | 1.7–9.5× |
| PFNA | 1.0–13.5× | 0.5–7.3× |
| **overall bias** | **2.33×** | 3.43× |
| **overall \|error\|** | **2.42×** | 3.95× |

Implemented from its own SI, Armitage is the better predictor overall and is
remarkably good on Rosenmai's short chains (0.8–1.7×). Qin wins at low dose and
loses badly at 100 µM, where its albumin term has saturated too far.

Qin's parameters are measured only for PFBA and PFOA; PFHxA and PFHpA are
interpolated between those two carboxylate anchors and PFNA is extrapolated
beyond both (open symbols in the figure). Its two worst compounds, PFHxA and
PFHpA, are interpolated ones.

## 6. Why Rosenmai stands out: the axis label, not the measurement

Rosenmai's Abstract says **"4–2500 ng mg⁻¹ protein"**; Results §3.3 and the
Fig. 2A axis say **nmol mg⁻¹ protein**. The full Results sentence settles that
these are the same numbers:

> "The cellular PFAS concentrations, expressed as nmol mg⁻¹ protein, for the 14
> PFASs at treatment concentrations of 10, 30 and 100 μM are shown in Figure
> 2(A). The concentrations ranged 18–2490 nmol mg⁻¹ protein for the PFCAs,
> 4–494 nmol mg⁻¹ protein for the PFSAs, and were 81 and 435 nmol mg⁻¹ protein
> for FOSA."

The union of those three sub-ranges is **4 to 2490** — minimum 4 from the
PFSAs, maximum 2490 from the PFCAs. The Abstract's range is **4 to 2500**.
Same minimum, same maximum to rounding. One set of numbers is being reported
under two different units, so **exactly one label is a typo and no reading
makes both correct**. (The Abstract is therefore not an independent vote for
ng/mg — it is the same data. Which label is wrong has to be settled on
physical grounds.)

Four checks say the numbers are ng/mg.

**(a) Rosenmai against Amstutz.** Same cell line, same dose, same stated unit:

| at 100 µM | Rosenmai | Amstutz | ratio | Rosenmai ÷ MW | ratio |
|---|---|---|---|---|---|
| PFHxA | 164 | 0.289 | 568× | 0.522 | 1.81× |
| PFHpA | 297 | 0.326 | 911× | 0.816 | 2.50× |
| PFOA | 525 | 1.815 | 289× | 1.268 | 0.70× |
| PFNA | 458 | 2.237 | 205× | 0.987 | 0.44× |

**(b) Rosenmai against itself.** Its uptake-% table and its Fig. 2A values are
the same measurement normalised differently, so the medium volume reconciling
them must be one number for all 15 rows. Read as ng/mg it is 1.03–2.20 mL for
the three pooled wells; read as nmol/mg it is 423–645 mL. A Corning 12-well is
22.1 mm across and ~17.5 mm deep, so one well brim-full is 6.7 mL and the whole
plate is 81 mL: the ng/mg reading needs 5–11% of three wells, the nmol/mg
reading needs **5–8× the entire plate**. It is the absolute magnitude that
discriminates, and it does so without anyone knowing the true working volume.

*A check that does NOT work, recorded so it is not repeated.* The implied
volume is constant within each compound across doses (PFOA 1.39 / 1.41 /
1.44 mL over a 10× range). An earlier version of this section called that
tightness decisive. It is not evidence at all: the two readings differ by a
constant molecular weight per compound, so they scale identically with dose and
the dose consistency is the same under either. Nor can reading one value off
the figure settle it — PFOA at 100 µM reads 525, which under ng/mg implies
0.47 mL/well and under nmol/mg implies 195 mL/well; without an independently
known working volume that observation only re-derives the volume, it does not
test the unit.

*A residual that survives the fix.* Between compounds the implied volume runs
2.14× (PFBA 2.11 mL down to PFNA 1.15 mL), decreasing monotonically with chain
length. Rosenmai's Table S2 internal-standard recoveries of 78–98% span only
1.26×, so something else is also going on in that paper.

**(c) Scale, by sub-class.** Back-calculating every compound in Rosenmai's
Table 1 at the PFOA-implied volume gives PFCAs 12–1510 ng/mg and PFSAs
8–96 ng/mg, against the stated 18–2490 and 4–494. The same scale. Read as
nmol/mg the stated minima are 330× (PFCA) and 520× (PFSA) above what the
uptake table supports. The PFSA minimum of 4 is the sharpest of these, because
PFOS has the lowest uptake of any compound in Table 1 (0.04–0.05%) and so pins
the bottom of the range.

The FOSA entry corroborates that these are the Fig. 2A numbers: it has only
*two* values, 81 and 435, because 100 µM was dropped for FOSA on cytotoxicity
grounds — so Fig. 2A carries two FOSA bars, not three.

**(d) Physical ceiling.** As nmol/mg, Fig. 2A puts PFOA at 253 mM and PFNA at
266 mM inside the cell — a quarter of the cell by mass.

**(e) A second label error in the same paper.** The supporting information's
Figure S2 gives total protein per well in **ng/cm²**, with bars at a few
hundred. A 12-well is 3.5 cm², so that reads as a few micrograms of protein per
well, about 40× too little for 1.8 × 10⁵ HepG2 cells; as µg/cm² it is ~20× too
much. Whichever way it is read it cannot be right as printed, which establishes
that unit labels in this paper are not reliable.

Reading the axis as nmol/mg multiplies C_cell by exactly the molecular weight,
**214× (PFBA) to 464× (PFNA)** — 2.33–2.67 decades — and because the factor
tracks MW it also tilts the chain-length trend by 0.34 decades.

On a fraction-of-dose basis Rosenmai does not stand out at all. The measurement
is sound; the conversion to concentration is where the outlier is created.

## 7. Which axis to plot

Every source agrees on **protein per cell** — Armitage's AREc2 gives 694
pg/cell, the project spreadsheet 675 pg/cell — and disagrees **3.3× on cell
volume** (10 pL vs 3.04 pL). Since all three datasets report per mg protein,
the *fraction of dose in cells* is insensitive to that disagreement while
*C_cell* carries all of it. Model discrimination should be done on the fraction.

## 8. Geometry, now sourced rather than assumed

| | plate | medium | cells | FBS |
|---|---|---|---|---|
| Rosenmai | 12-well × 3 pooled | 3 × 1 mL (assumed; §6b implies ~0.7 mL/well) | 5.43 × 10⁵ | 10% |
| Amstutz | 96-well | 250 µL | 4 × 10⁴ | 10% |
| Fragki | 24-well | **550 µL** (stated) | seeded 55 000/well, differentiated | 10% |

Fragki's main paper confirms the plate and volume I had assumed. The cell count
at exposure is still not stated for differentiated HepaRG; 2.2 × 10⁵ is used.
Because the measured fraction and the predicted fraction both scale with cell
number, the ratio is insensitive to it.

## 9. Remaining caveats

- Fragki's 6.25 µM values in the project sheet (PFOA 0.0562, PFNA 0.281) still
  disagree with supplementary Table 2 (0.09226, 0.22203); the 25 and 100 µM
  values match to 1–3%. On the supplementary values Fragki PFOA gives b = 1.23
  rather than 1.42, so the direction is robust.
- Armitage 2021 did run PFOS through IV-MBM EQP v2.0 (workbook sheets AEID767,
  AEID1325), which would make a direct validation target, but the sheets' two-
  row headers defeat clean extraction and the PFAS rows could not be read
  reliably. Not pursued.
- AREc2's water content of 0.944 is high for a cell and its protein of 0.0514
  low; it is used because it is the composition Armitage's own model was run
  with for a HepG2-derived line, not because it is independently plausible.

## 10. The fraction-of-dose axis, and the project spreadsheet on it

Figure 10 is now built on fraction of dose in cells (§7). Converting the
spreadsheet's nmol/mg columns with the agreed 675 pg protein per cell:

| over 29 points | bias | \|error\| |
|---|---|---|
| project sheet Armitage | **0.90×** | 2.86× |
| published-equation Armitage | 2.33× | **2.42×** |
| published-equation Qin | 3.43× | 3.95× |
| project sheet Qin | **57×** | 57× |

The sheet's Armitage is essentially unbiased and runs 0.39× of the
published-equation version on average — the two bracket the data.

The sheet's Qin is a different matter. Its bias and its absolute error are
equal to four figures, so every point is high by nearly the same factor: a
scale error, not scatter. On the fraction axis the failure is obvious —
it places **6–82% of the dose inside cells that occupy 0.18% of the well
volume**, where Rosenmai and Amstutz both measure 0.04–0.33%. Two candidates,
in order of size:

1. **D_cell/w used as a lumped ratio** rather than composed as
   `D_SP/w × V_protein,cell + D_lip/w × V_lipid,cell`. For PFBA the lumped
   route gives 31 per litre of whole cell against the composed route's 2.4 —
   13×.
2. **Dividing by C_nominal rather than C_free.** Qin's ratios are defined
   against the freely dissolved concentration. Small for PFBA (~1.25×) but
   large by PFOA, where albumin holds 89% of the dose.

Neither is visible in the sheet's output columns, so these are checks to run
rather than a decomposition.

### Rosenmai on the fraction axis

This is the cleanest form of the unit problem, because the axis is unit-free
and the comparison is internal to one paper. At 100 µM:

| | Table 1 uptake | Fig 2A as ng/mg | Fig 2A as nmol/mg |
|---|---|---|---|
| PFBA | 0.100% | 0.072% | 15.4% |
| PFHxA | 0.140% | 0.064% | 20.1% |
| PFHpA | 0.210% | 0.100% | 36.3% |
| PFOA | 0.330% | 0.155% | 64.2% |
| PFNA | 0.310% | 0.121% | 56.0% |

Across all 15 rows the nmol/mg reading contradicts Rosenmai's own uptake table
by **141–215×**, and would put 15–64% of the dose inside 0.18% of the well
volume. The ng/mg reading sits within 1.4–2.1× of the uptake table, the
residual being the unstated working volume (§6b implies ~0.7 mL/well against
the 1 mL assumed).

## 11. The wash artefact hypothesis: refuted on three counts

Earlier sections floated five ice-cold PBS/DPBS washes as a chain-length-
dependent loss that would bias measured C_cell low. It does not survive.

**(a) It is not a between-study variable.** All three studies wash five times:
Rosenmai "washed 5 times with ice-cold phosphate-buffered saline", Amstutz
"washed... with DPBS" following Rosenmai's protocol, and Fragki "Cells were
washed 5 times with Dulbecco's Phosphate Buffered Saline". A shared protocol
cannot explain differences between them.

**(b) The direction is backwards.** Washing should strip the shortest, least
bound, fastest-desorbing chains hardest, so the models should over-predict
*short* chains most. They do the opposite — Armitage's over-prediction rises
monotonically with chain length:

| | nCF | Armitage/measured | Qin/measured |
|---|---|---|---|
| PFBA | 3 | 1.42× | 2.36× |
| PFHxA | 5 | 1.20× | 8.78× |
| PFHpA | 6 | 1.55× | 9.03× |
| PFOA | 7 | 2.68× | 3.80× |
| PFNA | 8 | 3.86× | 1.51× |

Armitage: **+0.091 decades per CF2, r = +0.83**. Qin shows no trend
(−0.035, r = −0.20).

**(c) The magnitude is two orders out.** Fragki's time course puts steady state
at ~1 h, so the cell↔medium exchange half-time is ~10–30 min. Five quick rinses
give 1–2 min of contact:

| exchange t½ | 1 min contact | 2 min | 5 min |
|---|---|---|---|
| 10 min | 6.7% lost | 12.9% | 29.3% |
| 20 min | 3.4% | 6.7% | 15.9% |
| 30 min | 2.3% | 4.5% | 10.9% |

A 3–5× loss needs 1.5–2.3 exchange half-lives of contact — 15–45 minutes of
washing. The hypothesis is dead at the order-of-magnitude level.

The one thing it was invoked to explain, Amstutz reading 3.6–4.8× below
Rosenmai for PFHxA and PFHpA while agreeing to 1.0–1.4× for PFOA and PFNA, is
better explained by short-chain analytical recovery: Rosenmai's own Table S2
reports internal-standard recoveries of 78–98% that vary by compound.

---

## 12. Testing the albumin parameter against the best measurement

Added 2026-10-16. Reproduce with
`python3 tk_review/scripts/mbm_fischer_albumin.py`; full detail in
`report/REPORT.md` §6.1b–6.1c.

Section 4 of this memo found the divergence between the two models to be in
albumin, not in the cell. That raised an obvious question: can either model be
improved by giving it a *measured* albumin partition coefficient instead of
its own? Fischer 2024 is the only albumin dataset taken at a ligand:protein
ratio below 0.004, so it is the natural reference. The assays run on fetal
bovine serum, so its BSA column applies.

Converting a binding constant into a partition coefficient is exact in the
trace limit — `D = 1000·n·Ka / M_alb`, so `log D = log Ka + log n − 1.823`
(mass basis) or `− 1.688` (the volume basis both models use). Qin's own
albumin term is this expression: identifying `Cmax = n × 20,527 µM` reproduces
Qin's two Cmax columns to 0.3% on all four of their measured compounds.

**Result: the substitution makes both models worse.**

| variant | bias | \|err\| | ≤25 µM | 100 µM | dose swing |
|---|---|---|---|---|---|
| Qin, published | 3.43× | 3.95× | 1.58× | 5.52× | 3.49× |
| Qin + Fischer, constant D | 4.74× | 4.77× | 4.78× | 4.72× | **0.99×** |
| Qin + Fischer, rescaled Langmuir | 9.68× | 9.68× | 7.66× | 11.18× | 1.46× |
| Armitage, published | **2.33×** | **2.42×** | 2.73× | 2.12× | 0.78× |
| Armitage + Fischer | 5.09× | 5.10× | 5.79× | 4.71× | 0.81× |

Fischer's value is below both models' for every compound (0.14–0.61× of Qin,
0.17–0.61× of Armitage), so substituting it shrinks the albumin sink and
pushes more dose into the cell — and both models already over-predict.

**Two findings survive the negative result.**

*Qin's dose dependence is an albumin-saturation artefact.* Its published error
swings 3.49× across the dose range; with a constant albumin coefficient the
swing is 0.99×. The dose exponent Qin produces (1.15–1.82 against a measured
0.85–1.42, section 4 of this memo) therefore says nothing about PFAS uptake —
it is the Langmuir term coming off saturation.

*Albumin cannot carry the level error.* Inverting the balance for the D_alb
that would make each model unbiased gives values **above** each model's
current one, while the measurement sits **below** it. For PFHxA the needed
value exceeds anything ever measured for that compound by 8.6–13.2×, even
allowing six binding sites. So albumin is where the sensitivity is (it holds
19–94% of sorptive capacity) but not where the error is, and the remaining
candidates are the cell term, the plastic term, and the measured C_cell values.

---

## 13. The cell term: slope exonerated, level convicted

Added 2026-10-19. Reproduce with `python3 tk_review/scripts/mbm_cell_term.py`;
full detail in `report/REPORT.md` §6.1d.

Section 12 ruled albumin out as the source of the over-prediction. The cell
term is the other large compartment, and it is where the two models disagree
most about mechanism while agreeing most about value.

**Composition.** Armitage's cell is 79% water at PFBA and 99.7% membrane lipid
at PFNA, with protein contributing 0.1–0.2% throughout (only 3.5% of protein
is sorptive, and that fraction is treated as octanol). Qin's cell is 92–95%
protein at every chain length, with no water phase. They differ on the cell
partition by 0.44–9.5×, yet their predicted fraction of dose in the cell
agrees to 0.6–2.5×, because the albumin term damps the difference. Agreement
on the output is not corroboration of the mechanism.

**The slope.** Against Qin's own Table 2 both models look 2.9–4.4× too steep in
chain length. That reading is wrong. Inverting the balance for the cell
partition the uptake data require gives +0.549 decades per CF₂; the models
give +0.515 and +0.782; Ryu 2024's whole-tissue f_unbound gives 0.46–0.53.
Qin's Table 2 alone gives +0.177, and at PFBA it is 20× above what the uptake
requires. Three independent lines against one: the models' slope is fine.

Substituting Qin's measured D_cell/w confirms it. Mean bias improves (Qin
3.43 → 1.99, Armitage 2.33 → 1.25) while absolute error worsens (3.95 → 4.80,
2.42 → 3.09), with PFBA over-predicted 20–24× and PFNA under-predicted to
0.27–0.31×. Centring the fit while destroying the trend is what substituting
an outlier looks like.

**The level.** Qin's cell partition runs 1.5–9.2× above the uptake requirement
and Armitage's 0.25–6.2×, both worst mid-series at PFHxA–PFHpA. Those factors
are conditional: redoing the inversion under Qin's, Armitage's and Fischer's
albumin moves the implied level 2.6–8.3×, though the slope survives all three
(+0.72 to +0.82). Read the level as "several-fold too high, worst mid-chain".

**What to fix.** Not the lipid or octanol descriptor — the chain-length
dependence is already right. The two remaining candidates are the sorptive
capacity assigned to cell protein (Qin's measured D_non-specific versus
Armitage's "3.5% of protein behaves as octanol" — wildly different guesses at
the same quantity, and the substance of each model's cell term) and the cell
volume the result is normalised by, which is uncertain to 3.1× across the
three compositions in play.

---

## 14. Candidate 1 eliminated; the cell volume is what is left

Added 2026-10-19. Reproduce with `python3 tk_review/scripts/mbm_cell_protein.py`;
detail in `report/REPORT.md` §6.1e. This closes the mass-balance line.

Two papers made the test possible. **Qin et al. 2024** (ES&T 58:9954, the same
group as the 2023 model) measured specific and nonspecific protein–water
distribution ratios for 16 PFAS by SPME over four orders of magnitude.
**Chen et al. 2025** (ES&T 59:82) independently writes the Kd → partition
coefficient bridge of §6.1b as its Equation 4.

**Candidate 1 was cell protein.** Qin's 2023 interpolation turns out to sit
within 1.1–2.5× of his own later measurement with the same slope, so the model
was not damaged by interpolating three of the five compounds. Armitage's
octanol surrogate is another matter: 0.035·D_ow runs **457–7,820× below** the
measured protein affinity and 1.5× too steep. That is why Armitage's cell is
~100% membrane lipid — the convention has defined its protein term out of
existence.

But substituting the measured value makes both models **worse** (Qin 3.43 →
5.92, Armitage 2.33 → 6.10), because the measurement is *higher* than either
model assumes and both already over-predict.

**The pattern across all three substitutions:**

| parameter | best measurement | direction | effect |
|---|---|---|---|
| albumin (Fischer 2024) | below both models | smaller competing sink | worse |
| cell D_cell/w (Qin 2023 T2) | mixed | — | centres, disperses |
| cell protein (Qin 2024) | above both models | bigger cell sink | worse |

All three move dose into the cell. No partition coefficient fixes an
over-prediction here, because the measured values do not point that way.

**The cell volume does.** Rerunning the published models at 3.04 pL/cell
(project sheet, HepG2) instead of 10 ng/cell (AREc2) gives Qin bias 1.60×
(|err| 3.95 → 2.26) and **Armitage 1.08×**. It is the one knob that rescales
every compound by the same factor, which is what a level offset independent of
chain length requires — and the assays modelled are HepG2 and HepaRG, not
AREc2. None of the three uptake papers reports its cell volume, so that is
where this stops.

**One result worth keeping regardless.** Qin 2024 reports both binding regimes
for the same protein and compound: specific exceeds nonspecific by 3.0–5.0×
across PFBA–PFOA, with specific dominating below ν ≈ 0.01 and nonspecific
above ν ≈ 1. That is the direct measurement of the ligand:protein regime
effect §6.1a inferred, and it quantifies the penalty for carrying an in vitro
constant into a physiological calculation without asking which regime produced
it.
