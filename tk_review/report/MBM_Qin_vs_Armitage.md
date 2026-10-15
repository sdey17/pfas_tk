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
Fig. 2A axis say **"18–2490 nmol mg⁻¹ protein"**. Four independent checks say
ng/mg.

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
the three pooled wells, and *within each compound* constant to 4–25% (PFOA
1.39 / 1.41 / 1.44 mL across a 10× dose range). Read as nmol/mg it is
423–645 mL in a 12-well plate.

**(c) Physical ceiling.** As nmol/mg, Fig. 2A puts PFOA at 253 mM and PFNA at
266 mM inside the cell — a quarter of the cell by mass.

**(d) A second label error in the same paper.** The supporting information's
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
