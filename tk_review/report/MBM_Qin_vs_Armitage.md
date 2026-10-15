# Qin 2023 vs Armitage 2021 (IV-MBM EQP v2.0) on PFBA, PFHxA, PFHpA, PFOA, PFNA

Reproduce with `python3 tk_review/scripts/mbm_qin_vs_armitage.py` (19/19
internal checks against the papers' own published numbers) and
`python3 tk_review/scripts/make_mbm_figure.py`.
Outputs: `db/mbm_qin_vs_armitage.csv`, `figures/fig10_mbm_qin_vs_armitage.png`.

## 1. Where the Armitage parameters came from

Armitage 2021 puts the biopartitioning spLFERs in Supplementary Section S1,
which is not in the copy on hand (MDPI returns 403 for `/s1`). They were taken
instead from EPA's own implementation of the same model, `httk` `R/armitage.R`
(`USEPA/CompTox-ExpoCast-httk`, function `armitage_eval`). That file reproduces
Armitage 2021 Eq. 5 (`0.97 log Kow - 6.94`, Kramer), Eq. 6
(`0.56 log Kow - 4.635`, Fischer) and Eq. 3 (acid ion trapping,
`10^(cell pH - bulk pH)`) exactly, and its mass-balance denominator is Eq. 2 —
so it is the same parameterisation, not a reconstruction.

| | Qin 2023 | Armitage 2021 / httk |
|---|---|---|
| albumin | **measured** D_BSA/w, saturable (Eq. 16) | spLFER `0.71 log Kow_neutral + 0.42`, linear |
| anion penalty on albumin | measured | **none** (`SFbsa_acidic = 0`) |
| cell | **measured** lumped D_cell/w | composed: storage lipid + membrane lipid + structural protein + water |
| anion penalty on cell | measured anionic log Kow | generic: `SFkow 3.5`, `SFmw 1` |
| dose dependence | yes | none — linear everywhere |
| plastic | not in the SPME balance | Kramer or Fischer spLFER |

**Volume basis.** Qin's Table 2 header reads `log Dcell/w [Lw/Lcell,prot+lip]`
and `log Dplasma/w [Lw/Lprot+lip]` — his distribution ratios are per litre of
*sorptive phase*, not per litre of cell or plasma. Armitage's `DR_kcw` is per
litre of whole cell (its composition sums to 1.0 and includes 0.90 water). Qin's
cell term is therefore scaled by the cell's protein+lipid volume fraction,
`F_PROTLIP = 0.1265` (675 pg protein+lipid in a 3.04 pL cell), before the two
share an axis — a −0.90 decade shift on Qin's cell only. The same reading is
what reproduces Qin's own Table 2: Eq. 17 with D_BSA/w = 10^4.60, D_lip/w =
10^3.52 and plasma at 93.3% / 6.7% protein:lipid by volume gives log D_plasma/w
= 4.573 against his published MBM value of 4.56, with 99.4% of the binding on
protein against his 99.1%.

Both sets of albumin constants sit on the same volume basis. Qin reports Cmax
both per litre of BSA and per mole; the ratio is 2.049e4 µmol/L for all four
compounds, i.e. 1364 g/L, which is 1/0.733 mL g⁻¹ — the partial specific volume
httk uses for `Valb`. The two models' albumin terms are therefore directly
comparable, and the script asserts this.

## 2. Parameter provenance per compound

Qin measured PFBA, PFOA, PFHxS and PFOS. Of the five compounds here only
**PFBA and PFOA are measured**; PFHxA and PFHpA are interpolated between those
two carboxylate anchors and **PFNA (8 fluorinated carbons) is extrapolated
beyond both**. Only log D_lip/w has a published rule (Fig. A1,
`0.476·nCF + 0.172`, verified to ±0.02 on both carboxylates and shown to fail
by 1.1 decades on the sulfonate PFHxS, which is why sulfonates are excluded as
anchors). Kd, Cmax, D_non-specific, the Freundlich constants and D_cell/w have
no published rule and are two-point interpolations. Open symbols in the figure
mark every extrapolated point.

log Kow is likewise reported by Qin only for PFBA and PFOA; the gaps are filled
at +0.770 per CF2 (OPERA) and +0.695 per CF2 (COSMOtherm), both carried in the
CSV.

## 3. Why the two models diverge: it is the cell, and it is a slope

Albumin is **not** the driver — the two treatments agree to within 0.4–2.6×
across all five compounds. The cell does not:

| | nCF | Qin log D_cell/w | Armitage log D_cell/w | ratio |
|---|---|---|---|---|
| PFBA | 3 | 1.66 | 0.46 | **15.8×** |
| PFHxA | 5 | 2.02 | 1.88 | 1.4× |
| PFHpA | 6 | 2.19 | 2.65 | 0.35× |
| PFOA | 7 | 2.37 | 3.43 | 0.09× |
| PFNA | 8 | 2.55 | 4.22 | **0.02×** |

(Qin's column is his measured value plus the −0.90 decade prot+lip conversion,
so both are per litre of whole cell.)

Qin's *measured* cell/water ratio rises **+0.18 per CF2**. Armitage's
*constructed* one rises **+0.75 per CF2** — four times steeper. They cross at
nCF 5.2, just above PFHxA, and diverge in opposite directions on either side:
a **700-fold swing** across five compounds. The prot+lip conversion moves the
crossover (it was at nCF 6.9 when Qin's ratio was wrongly applied to whole-cell
volume) but cannot touch the slope, which is what the disagreement is about.

The reason Armitage is that steep is that for an anionic PFAS its cell is
effectively a single compartment. The storage-lipid and structural-protein
terms both run on `Kow_anion = 10^(log Kow - 3.5)` and contribute 0.0–0.1% of
D_cell/w; the membrane-lipid term runs on `K_mw = 10^(1.01 log Kow + 0.12)`
with only a 1-decade anion penalty and carries **98.7–99.9%** of it for
anything from PFHxA up. So the whole Armitage PFAS cell limb reduces to

    D_cell/w  ~=  0.033 x 10^(1.01 log Kow - 0.88)

and inherits the chain-length slope of log Kow itself. At the short end the
anion penalty leaves nothing but the cell's water content (D_cell/w = 2.9 for
PFBA), which is why Armitage predicts ~0.1% of dose in cells for PFBA
regardless of anything else.

## 4. Second difference: dose dependence

Qin's albumin sink saturates, so his predicted cell fraction climbs with dose;
Armitage's is a constant. Over 10 → 100 µM Qin's lift grows with chain length
(PFBA 1.4×, PFOA 7.1×, PFNA 6.9×) while Armitage's is
1.00× for every compound.

The measured Rosenmai dose series flips sign at the same place: PFBA
(0.14→0.13→0.10%) and PFHxA (0.19→0.19→0.14%) **fall** with dose, PFHpA, PFOA
(0.24→0.27→0.33%) and PFNA (0.15→0.22→0.31%) **rise**. So Qin has the direction
right for the long chains and badly overshoots the magnitude (6.6× predicted
against 1.4× observed for PFOA); Armitage cannot produce the trend at all.

## 5. Performance against the data

Both models over-predict nearly everywhere, and each is best where the other is
worst — Armitage at the short end, Qin at the long end:

| | Armitage / measured | Qin / measured |
|---|---|---|
| PFBA | **0.8–1.2×** | 10–20× |
| PFHxA | 4.0–17× | 3.1–52× |
| PFHpA | 6.6–29× | 1.6–31× |
| PFOA | 7.2–14× | **0.4–2.7×** |
| PFNA | 2.1–29× | **0.1–1.0×** |

Armitage still holds the short end and Qin the long end, but the prot+lip
correction flips the *sign* of Qin's error above PFHpA: he now under-predicts
PFOA and PFNA rather than over-predicting them. Qin's best compounds remain the
extrapolated ones, and his worst (PFBA) is measured.

A candidate artefact runs the same way: Rosenmai and Amstutz both washed cells
five times with ice-cold PBS/DPBS before lysis, and the shortest, least-bound,
fastest-desorbing chains are exactly the ones those washes should strip. Amstutz
could not detect PFBA at all, which does rule out Qin's ~10%-of-dose prediction
for it, but the short-chain measurements are a lower bound in a chain-length-
dependent way.

## 6. Why Rosenmai stands out: the axis label, not the measurement

Rosenmai's own units are inconsistent. The Abstract says **"4–2500 ng mg⁻¹
protein"**; Results §3.3 and the Fig. 2A axis say **"18–2490 nmol mg⁻¹
protein"**. Three independent routes say ng/mg is the correct one.

**(a) Rosenmai's own uptake table.** Its Table 1 percentages are unit-free, so
they fix the absolute amount in the cells independently of any axis label.
Rosenmai pooled *three* wells per treatment, giving 543,400 cells and 0.3668 mg
protein per measured sample. Taking the plotted PFBA values of 18 / 45 / 126 at
10 / 30 / 100 µM and asking what medium volume reconciles them with the uptake
table:

| reading | 10 µM | 30 µM | 100 µM | per well |
|---|---|---|---|---|
| as **ng/mg** | 2.20 mL | 1.98 mL | 2.16 mL | 0.66–0.73 mL |
| as nmol/mg | 472 mL | 423 mL | 462 mL | 141–157 mL |

The ng/mg reading lands on a single working volume consistent across three
doses to **±6%**, and 0.7 mL/well is an ordinary 12-well volume. The nmol/mg
reading needs 423–472 mL in a 12-well plate.

**(b) Cross-study.** Fragki's HepaRG PFOA at 100 µM is 2.79 nmol/mg. Rosenmai's
uptake-derived PFOA at 100 µM is 3.40 nmol/mg — agreement to 1.2×.

**(c) Physical ceiling.** Read as nmol/mg, Fig. 2A puts PFOA at 253 mM and PFNA
at 266 mM inside the cell — a quarter of the cell by mass.

Reading the axis as nmol/mg multiplies C_cell by exactly the molecular weight:
**214× for PFBA rising to 464× for PFNA**, i.e. 2.33–2.67 decades, and because
the factor scales with MW it also tilts the chain-length trend by 0.34 decades.

On a fraction-of-dose basis Rosenmai does not stand out at all. At 100 µM it
agrees with Amstutz to **1.35× for PFOA and 1.03× for PFNA** (the two agree less
well for the short chains, 3.6× and 4.8×, consistent with the wash artefact
above). The measurement is sound; the conversion to concentration is where the
outlier is created.

## 7. Caveats carried in the code

- Rosenmai's working volume is not stated. 1.0 mL/well × 3 pooled wells is used
  in the model runs; the ng/mg reconciliation in §6 independently implies
  ~0.70 mL/well. The difference is 1.4× on the fractions and does not affect
  any conclusion here.
- Rosenmai pooled three wells per treatment, so the measured unit is 543,400
  cells and 0.3668 mg protein, not one well. Uptake % is invariant to the
  pooling; the absolute nmol and the protein mass are not.
- Fragki's plate format and working volume are not in the supplementary file on
  hand; 24-well / 500 µL is assumed. This affects Fragki's *fraction* but not
  its C_cell, which comes straight from nmol/mg.
- Fragki's µg/L column is a whole-liver scaling (126 mg protein/g liver), not a
  cell concentration. The nmol/mg column is the one used here.
- The protein → cell-volume factor (222 mg protein per mL of cell, from
  675 pg protein in a 3.04 pL cell) is common-mode: it moves Rosenmai and Fragki together and does not touch
  Amstutz, so it cannot generate between-study scatter, only a uniform offset.
- Amstutz's intracellular fraction is normalised to **recovered** mass
  (cells + medium at 24 h), not nominal; its mass balance was 93–106%, so the
  two bases are within a few percent here. The model is scored on the matching
  basis for each study.
- `F_PROTLIP = 0.1265` is a ≤1.6× lever on Qin's cell term alone (httk's cell is
  0.10 non-water; a 1.35 g/mL protein density on 675 pg/cell gives 0.166). It
  shifts Qin's curve vertically and does not change the chain-length slope that
  the Qin/Armitage disagreement consists of.

## 8. Reconciliation with the project spreadsheet

`python3 tk_review/scripts/mbm_vs_project_sheet.py`; the sheet's rows are in
`db/mbm_project_sheet.csv`, the diff in `db/mbm_sheet_reconciliation.csv`.

**The sheet's Qin column is linear.** C_cell / C_test is exactly constant
across every dose series in both the Qin and the Armitage column. Armitage
should be — it is a linear model. Qin should not: the saturable albumin term
(Eqs. 12–16) is what that paper is about, and without it the Qin column is
Armitage's structure carrying Qin's constants. The measured data are not linear
either; C_cell / C_test rises 3.1–3.5× over Fragki's dose range.

**Rosenmai vs Amstutz settles the unit question with no modelling at all.**
Same cell line, same dose, same stated unit:

| at 100 µM | Rosenmai | Amstutz | ratio | Rosenmai ÷ MW | ratio |
|---|---|---|---|---|---|
| PFHxA | 164 | 0.289 | 568× | 0.522 | 1.81× |
| PFHpA | 297 | 0.326 | 911× | 0.816 | 2.50× |
| PFOA | 525 | 1.815 | 289× | 1.268 | 0.70× |
| PFNA | 458 | 2.237 | 205× | 0.987 | 0.44× |

**And Rosenmai is internally inconsistent unless the axis is ng/mg.** Its
uptake-% table and its Fig. 2A values are the *same* measurement differently
normalised, so the medium volume implied by reconciling them must be one number
for all 15 rows. Read as ng/mg it is 1.03–2.20 mL for the three pooled wells —
and within each compound it is constant to 4–25% (PFOA 1.39 / 1.41 / 1.44 mL
across a 10× dose range; PFBA 1.98 / 2.16 / 2.20). Read as nmol/mg it is
423–645 mL in a 12-well plate.

**The dose exponent is the sharpest test of the two models.** Fitting
C_cell ∝ C_nom^b:

| series | n | b |
|---|---|---|
| Fragki PFNA | 5 | 1.42 |
| Fragki PFOA | 5 | 1.42 |
| Rosenmai PFNA | 3 | 1.28 |
| Rosenmai PFOA | 3 | 1.13 |
| Rosenmai PFHpA | 3 | 1.02 |
| Rosenmai PFBA | 3 | 0.85 |
| Rosenmai PFHxA | 3 | 0.80 |

Armitage requires b = 1.00 exactly and cannot represent any of this. Qin's
measured Freundlich slopes imply b = 1/n_Fr ≈ 1.62. The data sit between, and
b rises monotonically with chain length — a trend neither model reproduces,
since Qin's Freundlich slopes are nearly chain-length independent (−0.38 to
−0.39) and so predict ~1.62 for every compound. (Fitting Fragki PFOA on the
three supplementary-table values alone gives b = 1.23, so the direction is
robust to the low-dose discrepancy below.)

**One data discrepancy.** Fragki's supplementary Table 2 lists only 6.25 / 25 /
100 µM. The sheet's 25 and 100 µM values match it to 1–3%; its 6.25 µM values
do not (PFOA 0.0562 against 0.09226; PFNA 0.281 against 0.22203), and 12.5 /
50 µM are not in the supplementary at all. Worth checking against the main
paper, which is not on hand here.

**Caveat on the sheet's D_cell/w column.** It divides by C_nominal, not by the
freely dissolved C_w that Qin's ratios are defined against. For PFBA albumin
sequesters only ~30% of the dose so the error is ~1.4×, but by PFOA the free
fraction in 10% FBS is a few percent, so the same formula understates D_cell/w
by 20–50× at the long-chain end.
