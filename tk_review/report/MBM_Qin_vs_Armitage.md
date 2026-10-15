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
| PFBA | 3 | 2.56 | 0.46 | **125×** |
| PFHxA | 5 | 2.92 | 1.88 | 11× |
| PFHpA | 6 | 3.09 | 2.65 | 2.8× |
| PFOA | 7 | 3.27 | 3.43 | 0.7× |
| PFNA | 8 | 3.45 | 4.22 | **0.2×** |

Qin's *measured* cell/water ratio rises **+0.18 per CF2**. Armitage's
*constructed* one rises **+0.94 per CF2** — five times steeper. They cross
between PFHpA and PFOA and diverge in opposite directions on either side: a
600-fold swing across five compounds.

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
(PFBA 1.3×, PFHxA 3.4×, PFHpA 5.3×, PFOA 6.6×, PFNA 6.7×) while Armitage's is
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
| PFBA | **0.8–1.2×** | 71–135× |
| PFHxA | 3.9–17× | 23–343× |
| PFHpA | 6.5–29× | 12–212× |
| PFOA | 7.6–17× | 3.4–20× |
| PFNA | 2.6–28× | **0.9–7×** |

Uncomfortably, Qin's single best compound (PFNA) is the one whose parameters are
extrapolated, and his worst (PFBA) is measured.

A candidate artefact runs the same way: Rosenmai and Amstutz both washed cells
five times with ice-cold PBS/DPBS before lysis, and the shortest, least-bound,
fastest-desorbing chains are exactly the ones those washes should strip. Amstutz
could not detect PFBA at all, which does rule out Qin's ~10%-of-dose prediction
for it, but the short-chain measurements are a lower bound in a chain-length-
dependent way.

## 6. Why Rosenmai stands out: the axis label, not the measurement

Rosenmai's own units are inconsistent. The Abstract says **"4–2500 ng mg⁻¹
protein"**; Results §3.3 and the Fig. 2A axis say **"18–2490 nmol mg⁻¹
protein"**. Two independent routes say ng/mg is the correct one.

1. Rosenmai's Table 1 uptake percentages are unit-free. Back-calculating from
   them, with the paper's 18×10⁴ cells/well and a 12-well working volume, gives
   **0.54 mg protein per 10⁶ cells** — a standard HepG2 value, confirming the
   geometry — and cell concentrations of **31–1480 ng/mg**, inside the
   Abstract's range and nowhere near 18–2490 nmol/mg.
2. Fragki's HepaRG PFOA at 100 µM is 2.79 nmol/mg. Rosenmai's uptake-derived
   PFOA at 100 µM is 3.40 nmol/mg — agreement to 1.2×.

Reading the axis as nmol/mg multiplies C_cell by exactly the molecular weight:
**214× for PFBA rising to 464× for PFNA**, i.e. 2.33–2.67 decades, and because
the factor scales with MW it also tilts the chain-length trend by 0.34 decades.
At 100 µM it puts PFOA at 253 mM and PFNA at 266 mM inside the cell — a quarter
of the cell by mass, which is not physical.

On a fraction-of-dose basis Rosenmai does not stand out at all. At 100 µM it
agrees with Amstutz to **1.35× for PFOA and 1.03× for PFNA** (the two agree less
well for the short chains, 3.6× and 4.8×, consistent with the wash artefact
above). The measurement is sound; the conversion to concentration is where the
outlier is created.

## 7. Caveats carried in the code

- Rosenmai's working volume is not stated; 1.0 mL (standard 12-well) is assumed
  and flagged. The derived 0.54 mg protein/10⁶ cells is the check on it.
- Fragki's plate format and working volume are not in the supplementary file on
  hand; 24-well / 500 µL is assumed. This affects Fragki's *fraction* but not
  its C_cell, which comes straight from nmol/mg.
- Fragki's µg/L column is a whole-liver scaling (126 mg protein/g liver), not a
  cell concentration. The nmol/mg column is the one used here.
- The protein → cell-volume factor (180 mg protein per mL of cell) is
  common-mode: it moves Rosenmai and Fragki together and does not touch
  Amstutz, so it cannot generate between-study scatter, only a uniform offset.
- Amstutz's intracellular fraction is normalised to **recovered** mass
  (cells + medium at 24 h), not nominal; its mass balance was 93–106%, so the
  two bases are within a few percent here. The model is scored on the matching
  basis for each study.
