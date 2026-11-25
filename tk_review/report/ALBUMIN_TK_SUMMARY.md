# Can albumin binding carry a toxicokinetic model?

**A short summary of three attempts and what each established**

Souvik Dey · Office of Environmental Health Hazard Assessment, California
Environmental Protection Agency, Sacramento, CA

---

## The question

Plasma protein binding is the standard explanation for why perfluoroalkyl
acids are eliminated so slowly, and albumin is the protein almost always
named. If that explanation carries the kinetics, then an albumin binding
constant should be enough to predict a half-life. We tested that directly,
in three increasingly permissive forms, on the nine human compounds with a
measured volume of distribution and total clearance.

The governing identity throughout is

> t½ = ln2 · V_d / (f_u · GFR · T)

where f_u is the unbound fraction in plasma and T is the tubular factor, the
share of the filtered load surviving the renal tubule. Albumin binding enters
only through f_u, by the single-site relation f_u = 1/(1 + n[Alb]K_a) with
[Alb] = 0.6 mM.

---

## Attempt 1 — binding alone, nothing fitted

Take K_a, convert to f_u, assume T = 1, predict. No parameters.

**It fails, and the arithmetic says why before the model is run.** A binding
constant has units of M⁻¹ and a half-life has units of time; with no volume
and no clearance scale supplied, the relation can be *calibrated* but not
*derived*. Supplying V_d and GFR fixes the units but not the problem:

| quantity | median | range |
|---|---|---|
| f_u from K_a vs measured f_u | **10.6-fold** error | 0.2 to 1395 |

The binding-derived unbound fraction disagrees with the measured one by a
median factor of ten and, at the extreme, by three orders of magnitude. The
half-life predictions are correspondingly short, because T = 1 credits the
kidney with clearing everything it filters.

**What this establishes:** binding data alone cannot supply a human unbound
fraction to better than about tenfold, so it cannot supply a half-life.

---

## Attempt 2 — a fitted regression on the binding constant

Abandon the mechanism and calibrate: fit log t½ = a + b·log K_a within a
species, on one binding source, scored leave-one-out so each compound is
predicted from a fit that excluded it.

| species | binding source | n | LOO median | LOO worst | fitted slope |
|---|---|---|---|---|---|
| human | Gao (dialysis) | 9 | **2.4-fold** | 9.7 | +2.56 |
| human | Chen (dialysis) | 9 | 3.4 | 87.1 | +0.65 |
| human | Li (fluorescence) | 9 | 5.7 | 173.6 | +0.49 |
| rat | Fan (calorimetry) | 8 | 13.8 | 95.3 | **+1.45** |
| rat | Fan (ultrafiltration) | 8 | 14.7 | 39.2 | **−1.30** |

This buys accuracy and costs the mechanism. The calibration absorbs V_d, GFR
and T together, so the systematic shortfall of Attempt 1 disappears by
construction — but a and b have no physiological reading, they differ between
species, and because they are estimated *from* half-lives the model cannot be
carried to a species whose half-lives are the unknown, which is the usual
reason for wanting it.

**The slope is a diagnostic, and it fails.** If V_d, GFR and T do not track
binding, then t½ ∝ 1/f_u ∝ K_a and the slope must be **+1**. The human arms
give +0.49, +0.65 and +2.56 — a fivefold span on the same nine compounds,
differing only in which laboratory supplied K_a. Worse, the two rat arms come
from **the same paper by two methods** and differ in *sign*: −1.30 by
ultrafiltration against +1.45 by calorimetry. A relationship whose direction
depends on the assay cannot be reporting a property of the binding.

**Two coverage facts belong with the table.** No usable mouse albumin binding
constant exists in this compilation, so the mouse arm cannot be built under
this constraint. One mouse albumin row is present — Fischer et al. (2025) —
but it is a PBTK parameterisation carrying no numeric K_a or K_d and naming
its chemicals in prose ("9 PFAA with 4–10 perfluorinated carbons") rather than
per compound, so it cannot enter a per-compound regression. And Fan's two methods disagree by up to 63-fold on a single
compound within that one paper (PFBS: 2.78 µM by ultrafiltration against
176.30 µM by calorimetry).

**What this establishes:** a fitted binding–half-life relationship is as much
a property of the assay as of the chemistry, and is not transferable between
species.

---

## Attempt 3 — binding-derived f_u inside the mechanistic model

Keep the mechanism, keep V_d and GFR, fit *one* parameter — a constant
tubular factor — and let albumin binding supply f_u in place of a measured
value. This is the only form in which binding data earned their place.

Scored on the four long-chain compounds (PFOA, PFOS, PFNA, PFDA) whose
tubular factor is determined, leave-one-out:

| model | f_u from | k | median | worst |
|---|---|---|---|---|
| M1d | **albumin K_a** | 1 | **2.1-fold** | 3.6 |
| M1b | measured (Smeltz 2023) | 1 | 2.2-fold | 2.4 |
| M1a | measured, T = 1 | 0 | 62.4 | 133 |

**Binding-derived and measured unbound fractions perform the same here**, at
this sample size. That is the one positive result in this line of work:
albumin binding is a usable stand-in for a measured f_u *once a tubular
factor is fitted*. It is not an improvement on a measured f_u, and it does
nothing about T — the step from 62.4-fold to about 2-fold is bought entirely
by the fitted tubular constant, not by the binding.

**The number of binding sites is not identifiable.** Varying n over 1, 2 and
3 moves the median by 2.3%. Under tight binding f_u ≈ 1/(n[Alb]K_a) is
proportional to 1/n while the fitted T is proportional to n, so n cancels in
the product f_u·T. It survives only through the weakest binder, PFBA, which
moves 1.29-fold; no other compound moves more than 4.9%. These data therefore
cannot be used to infer stoichiometry.

---

## Bottom line

1. **Albumin binding cannot set the scale.** Alone it predicts no half-life;
   with V_d and GFR supplied it still misses the measured unbound fraction by
   a median tenfold.
2. **A fitted binding–half-life relation is assay-dependent to the point of
   sign reversal** within a single paper, and cannot cross species.
3. **Binding-derived f_u does work inside the mechanistic model** — matching a
   measured f_u on the long chains with one fitted parameter — but the
   parameter doing the work is the tubular factor, not the binding.
4. **Stoichiometry is unrecoverable** once T is fitted, by algebra rather than
   by noise.

This is consistent with the main finding of the larger analysis: binding is
arithmetically incapable of carrying the interspecies difference in the rat,
and the tubular term is where the difference lives. The modelling reaches the
same place from the forward direction.

---

## Caveats that bound all three attempts

- **Every human half-life and volume of distribution is from one volunteer**
  (Abraham et al. 2024). No fold error above is an estimate of population
  accuracy.
- **Four compounds per stratum** in Attempt 3, so at most two residual
  degrees of freedom per fitted model. The comparison is real; its size is
  not estimable at that n.
- Attempt 3's tubular factors are **renal-only upper bounds**, corrected for
  a measured faecal route; urinary clearance was below detection for every
  long-chain compound, so they are not measurements.
- f_u is Smeltz 2023 except PFDoDA, which Smeltz does not report and which
  takes Ryu's value.

---

## Scripts

| script | what it does |
|---|---|
| `scripts/stage25_binding_only_tk.py` | Attempt 1: binding-only mechanistic, two-link decomposition |
| `scripts/stage26_albumin_only_model.py` | Attempt 2: fitted log t½ ~ log K_a, slope diagnostic |
| `scripts/stage30_rerun_corrected.py` | Attempt 3: M1d within the stratified model comparison |
| `scripts/stage29_pfdoda_and_fu_consensus.py` | the f_u consensus and faecal correction both rely on |
