# In vitro mass balance for PFAS: a primer from first principles

**What this is.** Everything in `MBM_Qin_vs_Armitage.md` rests on one equation
and about ten approximations. This document derives the equation, names every
approximation, shows where each one breaks, and then reads the five papers
against it — two models and three datasets.

Read it with `scripts/mbm_published_equations.py` open. Every number quoted
here comes from that file, and every code block below runs as written.

---

## Part I — The one equation

### 1. What a partition coefficient actually is

Put a chemical in contact with two phases, water and something else. At
equilibrium the chemical potential is equal in both:

```
μ_A = μ_B           μ_i = μ_i° + RT ln(γ_i x_i)
```

Rearranging, the ratio of mole fractions is

```
x_A / x_B = (γ_B / γ_A) · exp( −(μ_A° − μ_B°) / RT )
```

The exponential is a constant at fixed temperature. So the ratio is constant
**if and only if the activity coefficients γ are constant**. Convert mole
fractions to concentrations and you have the partition coefficient

```
K_AB = C_A / C_B                    (constant)
```

> **Approximation 1 (dilute / linear regime).** γ is constant only when the
> chemical is a trace. Add enough of it and the sorbent changes — sites fill,
> the phase swells, the chemical starts interacting with itself. Then K is no
> longer a constant and becomes an *apparent* ratio that falls as you add more.
>
> This is the single most important approximation in the whole subject. Qin's
> model is built around relaxing it; Armitage's is not.

**Notation.** `K` is used for a true constant; `D` for a *distribution ratio*,
which is what you measure when the chemical exists in more than one form (here:
neutral and anionic) or when the sorbent is a mixture (a cell, serum). `D` is
what appears in every balance below.

### 2. The well as a set of phases

A well in a 96- or 24-well plate contains, from the model's point of view:

| phase | symbol | what it is |
|---|---|---|
| medium water | V_w | the aqueous phase, where "free" concentration lives |
| serum albumin | V_alb | from the FBS in the medium |
| serum lipid | V_slip | from the FBS |
| cells | V_cell | protein + lipid + cell water |
| plastic | A_plastic | the well wall and base, a surface not a volume |

Each has a distribution ratio to water, `D_i = C_i / C_w`.

### 3. Deriving the balance

Total moles are conserved and distributed over the phases:

```
n_total = C_w·V_w + C_alb·V_alb + C_slip·V_slip + C_cell·V_cell + C_pl·A_pl
```

Substitute `C_i = D_i · C_w` in every term — this is the step that uses
Approximation 1, five times:

```
n_total = C_w · ( V_w + D_alb·V_alb + D_slip·V_slip + D_cell·V_cell + D_pl·A_pl )
```

Define the **total sorptive capacity**

```
S ≡ V_w + Σ_i D_i · V_i                              [litres of water-equivalent]
```

and the balance is a one-liner:

```
C_w = n_total / S                                    (★)
```

Everything else follows. The concentration in the cells is `C_cell = D_cell·C_w`,
and the **fraction of the nominal dose that ends up inside cells** is

```
x_cell = C_cell·V_cell / n_total = D_cell·V_cell / S     (★★)
```

That is the whole model. Both papers are (★★) with different recipes for the
`D_i`.

### 4. The approximation ledger

Everything that had to be true to get from the well to (★):

| # | approximation | breaks when | does it bite here? |
|---|---|---|---|
| A1 | linear sorption, γ constant | sorbent saturates | **yes** — Qin measures K_d = 0.317 µM for PFOA against doses of 6.25–100 µM |
| A2 | equilibrium reached | slow uptake | no — Fragki's time course reaches steady state in ~1 h, exposure is 24 h |
| A3 | no metabolism or degradation | metabolised compounds | no for PFCA/PFSA; **yes** for fluorotelomer alcohols |
| A4 | no active transport | transporter-expressing cells | unknown — HepG2 express few PFAS transporters, but this is assumed, not shown |
| A5 | phases are independent and additive | one sorbent alters another | assumed |
| A6 | volumes additive and constant | — | fine |
| A7 | no losses (evaporation, tubing, pipettes) | volatiles | fine for ionised PFAS |
| A8 | the aqueous phase is pure water | co-solvents present | **3% ethanol / 0.5% / 0.1% DMSO** — see §6 |
| A9 | the cell is a well-mixed bag of phases | organelles, lysosomal trapping, pH gradients | assumed; Armitage's published model has no lysosome term |
| A10 | one `D_cell` describes all cells equally | different cell lines | **yes** — three assays, two cell lines, one composition used for all |

A1, A8 and A10 are the ones that have actually cost us something in this
project. A4 is the one nobody can check.

### 5. Two axes, and a trap worth naming

There are two ways to report the answer:

```
C_cell  = D_cell · n_total / S                      [concentration]
x_cell  = D_cell · V_cell / S                       [fraction of dose]
```

Read those as written and `x_cell` is the one carrying `V_cell`. That is
correct — and it is the opposite of what a careless reading of "which axis
needs the cell volume" suggests, because the question has two sides.

**On the model side.** `x_cell ∝ V_cell` directly. `C_cell` contains `V_cell`
only inside `S`, where the cell term is 0.1–18% of the total, so predicted
`C_cell` is nearly volume-free.

**On the measurement side**, converting a published "nmol per mg protein":

```
measured x_cell = (nmol/mg) · (mg protein in the well) / (nmol added)       no volume
measured C_cell = (nmol/mg) · (mg protein per cell) / (litres per cell)     needs volume
```

**The volume dependence does not disappear — it changes sides.** And the
quantity actually tested, predicted ÷ measured, is *identical* on the two
axes. Rosenmai PFOA at 100 µM, holding protein per cell fixed at 694 pg and
changing only the volume:

| ng/cell | meas x_cell | pred x_cell | **ratio** | meas C | pred C | **ratio** |
|---|---|---|---|---|---|---|
| 10.0 | 0.330% | 0.660% | **2.00×** | 182 µM | 364 µM | **2.00×** |
| 3.04 | 0.330% | 0.201% | **0.61×** | 599 µM | 366 µM | **0.61×** |

> **Choosing the fraction axis buys no protection against the cell-volume
> uncertainty.** Both ratios scale as (volume per cell)/(protein per cell).
> That is exactly why §14 of the memo can move Armitage's bias from 3.87× to
> 1.84× by swapping 10 ng/cell for 3.04 pL/cell — it is a real lever, not an
> artefact of the axis.

**What the fraction axis does buy** is that the *measurements* become mutually
comparable with no volume assumption at all. Rosenmai, Amstutz and Fragki can
be put on one plot and read against each other directly (fig. 20A). On the
concentration axis they could not be, because each would need a cell volume
nobody measured.

**A cancellation that does work — sometimes.** Cell *number* per well is also
uncertain (HepG2 divide between seeding and exposure; Fragki never state the
count for differentiated HepaRG). Whether it cancels depends on what the paper
reported:

| | reports | measured x_cell | ratio over a 4× change in N |
|---|---|---|---|
| Fragki | nmol/mg protein | built from N, so ∝ N | **8.69× → 7.89×**, cancels |
| Rosenmai | uptake % directly | a fixed published number | **1.00× → 3.97×**, does not |

So the memo's §8 remark — that the ratio is insensitive to cell number — holds
for Fragki, which is where it is made, and **not** for Rosenmai, whose Table 1
is already a percentage. For Rosenmai the assumed 5.43 × 10⁵ (three wells at
the *seeded* 1.8 × 10⁵) propagates straight into the comparison.

---

## Part II — Three results you can derive in two lines each

These are the levers that explain most of what we found. Each is a direct
consequence of (★★); none needs a model.

### 6. Why a co-solvent almost cancels

A co-solvent makes water a better solvent, so **every** water→phase partition
coefficient falls by roughly the same factor `k`. Put `D_i → D_i/k` into (★★)
and work it through. Split the capacity into its water part and the rest,
`S = V_w + Σ` with `Σ ≡ Σ_i D_i V_i`, and write `w ≡ V_w/S`, so `Σ = S(1−w)`:

```
x'_cell = (D_cell/k)·V_cell / ( V_w + Σ/k )

x'_cell     (D_cell/k)·V_cell        S              1          S
-------  =  ----------------- × ------------  =   ---  ×  ----------
x_cell       V_w + Σ/k           D_cell·V_cell      k      V_w + Σ/k
```

Substitute `V_w = wS` and `Σ = S(1−w)`; every `S` cancels:

```
x'_cell        1              1                  1                   1
-------  =    ---  ×  ----------------  =  -------------  =  ---------------
x_cell         k       w + (1−w)/k           k·w + 1 − w       1 + (k−1)·w
```

which is

```
x'_cell / x_cell = 1 / ( 1 + (k−1)·w )               (★★★)
```

Read it off: the numerator falls by `k`, and so does the denominator — but only
on the `(1−w)` of it that carries a `D`. The mismatch is exactly water's share.

Three sanity checks, all immediate: `k = 1` gives 1 (no co-solvent, no effect);
`w = 0` gives 1 (nothing dissolved in water, total cancellation); `w = 1` gives
`1/k` (no sorbents at all, the full effect). And it is **exact**, not a
linearisation — checked against the model to machine precision for all five
compounds.

- `w → 0` (long chains, albumin dominates): no effect at all.
- `w → 1` (short chains, water dominates): the full `1/k`.

Check it against the real run. Rosenmai, Armitage, 3% ethanol ⇒ k ≈ 1.35:

| | w | (★★★) predicts | model gives |
|---|---|---|---|
| PFBA | 0.80 | 1/(1+0.35·0.80) = **0.781** | 1.11/1.42 = **0.782** |
| PFNA | 0.005 | 1/(1+0.35·0.005) = **0.998** | 9.31/9.33 = **0.998** |

Exact to five decimals.

> **The approximation inside the result.** (★★★) assumes one `k` for every
> phase. Real cosolvency powers σ differ by sorbent, so strictly
> `D_cell → D_cell/k₁` and `D_alb → D_alb/k₂`. In the albumin-dominated limit
> the ratio becomes `≈ k₂/k₁`, so what survives is the *difference* in σ, not
> σ itself. A full unit of σ difference at 3% ethanol
> (σ_cell = 5 vs σ_alb = 4) is worth 10^(0.03×1) = 1.07, i.e. 7%. Still small
> against 5×, but it is not identically zero the way the single-`k` version is.

**An error that multiplies every sorptive term cancels on this axis. Only an
error in the *relative* size of two sinks survives.**

### 7. Why albumin's share sets the sensitivity

Let `a ≡ D_alb·V_alb / S` be albumin's share of the capacity. Multiply the
albumin coefficient by `f` (an error):

```
S' = S·(1 − a + f·a)        so      x'_cell / x_cell = 1 / (1 − a + f·a)   (★4)
```

Two things follow, and it is worth being careful about which is which.

**(i) How much a wrong albumin *coefficient* costs you.** For a 10×-too-low
one (`f = 0.1`), (★4) is `1/(1 − 0.9a)`. In Fragki's well:

| | a at 0% FBS | a at 2% | a at 10% | cost at 2% | cost at 10% |
|---|---|---|---|---|---|
| PFBA | 0.000 | 0.046 | 0.195 | 1.04× | 1.21× |
| PFHxA | 0.000 | 0.687 | 0.917 | 2.62× | 5.71× |
| PFOA | 0.000 | 0.884 | **0.974** | 4.89× | **8.09×** |
| PFNA | 0.000 | 0.791 | 0.948 | 3.47× | 6.80× |

At 0% serum the albumin parameter is irrelevant and the number is exactly 1.00
by construction. **A model validated on serum-free data says nothing about its
behaviour at 10% FBS** — and the same model is 4–8× more sensitive to its
albumin parameter in a 10% FBS well than in a 2% one.

**(ii) How much a wrong albumin *volume* costs you**, which is the Fragki case.
Coding 10% FBS where the paper says 2% leaves the coefficient alone and
multiplies `V_alb` by 5 — equivalent to `f = 0.2` going from the true well to
the coded one. (★4) with the true `a` then predicts the whole error:

| | a (true, 2% FBS) | (★4) with f = 0.2 | actual Armitage |
|---|---|---|---|
| PFBA | 0.046 | 1.19× | 1.18× |
| PFHxA | 0.687 | 3.75× | 3.75× |
| PFOA | 0.884 | **4.53×** | **4.54×** |
| PFNA | 0.791 | 4.14× | 4.17× |

Three lines of algebra reproduce a 4.5× modelling error to two decimals. That
is what it means to say the balance is an albumin calculation with a cell
attached.

> Note `a` is evaluated in the *reference* system — the one you are coming
> from. Use the coded system's `a` instead and you get the inverse factor, not
> the same number, because (★4) is not symmetric in `f`.

### 8. Why the residual's chain-length slope is a *ratio*

Once albumin dominates (`a` = 0.94–0.97 from PFOA up), `S ≈ D_alb·V_alb`, so

```
x_cell ≈ D_cell·V_cell / (D_alb·V_alb)
log x_cell ≈ log(D_cell/D_alb) + log(V_cell/V_alb)
```

The second term has no chain-length dependence. Differentiate with respect to
chain length `n` and subtract the same for the measurement:

```
d(log x_pred/x_meas)/dn  =  d log(D_cell/D_alb)_model/dn − d log(D_cell/D_alb)_true/dn   (★5)
```

**The residual's slope is the error in the slope of the ratio** — not in either
coefficient alone. A model can have both coefficients wrong and no residual
trend (if the errors track), or both nearly right and a strong trend (if they
don't). This is exactly what separates the two models; §13 below.

---

## Part III — The two models

Both are (★★). They differ entirely in how they produce the `D_i`.

### 9. Armitage 2021 (IV-MBM EQP v2.0)

**Strategy: predict every `D` from one descriptor, log K_ow.**

#### 9a. Ionisation

PFAS at pH 7.4 are essentially fully anionic. Armitage splits each `D` into a
neutral and an ionic part:

```
D_ow = [ F_n·10^logKow + F_i·10^(logKow − d_OW) ] / 10^(−k_salt·c_salt)
```

with `F_n` the neutral fraction, `F_i = 1 − F_n`, and `d_OW` the penalty the
anion pays. The denominator is a salting-out correction,
`k_salt = 0.04·logKow + 0.114` at `c_salt = 0.15 M`.

For PFAS `F_n = 0`, so every term reduces to the anionic branch. The penalties
in SI Table S1:

| phase | penalty `d` |
|---|---|
| octanol | 3.1 |
| membrane lipid | 1.0 |
| **serum albumin** | **0.0** |
| plastic | 3.1 |

> **Approximation A11.** `d_SaW_acid = 0` says the anion binds albumin exactly
> as well as the neutral form would. That is a modelling choice, not a
> measurement, and it is a large part of why albumin ends up holding 94–97% of
> the well in this model.

#### 9b. The spLFERs

"Sorbent-phase-specific linear free energy relationships" — one regression per
phase, from SI S1:

```
log K_mw = 1.01·logKow + 0.12                      membrane lipid–water
log K_pl = 0.97·logKow − 6.94                      plastic–water
log K_sa = 1.08·logKow − 0.70     (logKow < 4.5)   albumin, lower branch
log K_sa = 0.37·logKow + 2.56     (logKow ≥ 4.5)   albumin, upper branch
```

> **Approximation A12.** log K_ow is a sufficient descriptor. For a
> perfluorinated surfactant with an anionic head this is the weakest link in
> the model: PFAS are neither simply hydrophobic nor simply lipophilic, and
> their K_ow is itself hard to measure.
>
> **Approximation A13.** The albumin spLFER is *piecewise*, flattening to a
> slope of 0.37 above log K_ow 4.5 — calibrated on neutral organics and
> pharmaceuticals, where a plateau is physically sensible because the binding
> cavity has finite depth. Whether PFAS–albumin binding plateaus there is an
> open question, and §13 shows the consequences of assuming it does.

#### 9c. The cell as a sum of phases

```
D_cell/w = f_sl·D_ow + f_ml·D_mw + 0.035·f_prot·D_ow + f_water·1
```

Four terms: storage lipid behaves as octanol; membrane lipid as a phospholipid
membrane; **3.5% of the protein behaves as octanol and the rest is inert**; the
cell's own water behaves as water.

The AREc2 composition (SI Tables S3–S6), used for all three assays:

```
f_sl = 0.00046    f_ml = 0.00414    f_prot = 0.0514    f_water = 0.944
```

> **Approximation A14, the big one.** "3.5% of bulk protein behaves as octanol"
> is a convention imported from the neutral-organics literature. Measured
> against Qin 2024's protein–water ratios it is **457–7,820× too low**.
>
> **Provenance caveat, added on audit.** The 0.035 in this implementation comes
> from EPA's `httk`, and the code flags it as *not found in Armitage's SI S1*.
> So "Armitage's convention" — as §13–14 and §16 of the memo call it — is an
> attribution this project has not verified at source. The *consequence* is a
> property of the implementation either way: with a protein term this small the
> cell is lipid, and a lipid cell has too steep a chain-length slope. But whose
> convention it is remains open, and it matters, because if Armitage's own
> model treats cell protein differently then the §16 verdict lands on `httk`
> rather than on the paper. **Resolving this needs SI S1 §2, which is not in
> hand.** The
> consequence: the protein term contributes 0.1–0.2% of Armitage's cell, so
> the cell is 79% water at PFBA and **99.7% membrane lipid at PFNA**. The
> convention has defined protein out of the cell.

### 10. Qin 2023

**Strategy: measure the `D`s, and let albumin saturate.**

#### 10a. The Langmuir, derived

Take a protein with `n` identical, independent binding sites. For one site:

```
S + L ⇌ SL          K_a = [SL] / ([S]·[L])          K_d ≡ 1/K_a
```

Sites are conserved, `[S]_tot = [S] + [SL]`, so the occupancy is

```
θ = [SL]/[S]_tot = K_a[L] / (1 + K_a[L]) = [L] / (K_d + [L])
```

Now put it on a per-volume-of-albumin basis. The concentration of **sites
inside pure albumin** is just the molar concentration of albumin times `n`:

```
C_max = n · ρ_alb / M_alb = n · 1364 g/L ÷ 66,463 g/mol = n × 20,523 µM
```

So bound ligand per litre of albumin is `C_b = C_max·θ`, and the distribution
ratio is

```
D_specific = C_b / C_w = C_max / (K_d + C_w)
```

Add a linear, non-saturating term for everything that is not a defined site:

```
D_BSA/w = C_max / (K_d + C_w) + D_non-specific              Qin Eq. 16   (★6)
```

**Check the derivation against Qin's own table.** He reports C_max both per
litre of BSA and per mole of BSA; the ratio must be 20,523 µM:

| | C_max [µmol/L_BSA] | C_max [mol/mol] | ratio |
|---|---|---|---|
| PFBA | 9,220 | 0.450 | 20,489 |
| PFOA | 25,400 | 1.24 | 20,484 |

0.2% from the derived value. The derivation is right, and `n` = 0.45 for PFBA
and 1.24 for PFOA — **less than one site for PFBA**.

#### 10a-bis. How badly Approximation 1 breaks here

A1 does not break *in* Qin's model. **Qin's model is what you build because A1
broke in the experiment** — he measured the isotherm, found it curved, and
wrote (★6) to carry the curvature. Armitage keeps A1; Qin abandons it. Here is
how much was at stake.

**Albumin has about one site.** From the C_max identity above, `n` = 0.45
(PFBA), 1.24 (PFOA), 0.92 (PFHxS), 0.39 (PFOS). So albumin's high-affinity
capacity for PFAS is bounded at roughly **one mole per mole of albumin**.

**Now count the albumin actually in the well.** Armitage's serum albumin is
24 g/L, so the medium holds `24 × f_FBS` g/L:

| | FBS | albumin | doses | **nominal PFAS : albumin** |
|---|---|---|---|---|
| Rosenmai | 10% | 2.40 g/L = **36.1 µM** | 10–100 µM | 0.28 – 2.77 mol/mol |
| Amstutz | 10% | 2.40 g/L = **36.1 µM** | 100 µM | 2.77 mol/mol |
| Fragki | 2% | 0.48 g/L = **7.2 µM** | 6.25–100 µM | 0.87 – **13.85 mol/mol** |

Against `n ≈ 1`, a ratio of 2.8 means roughly three times more PFAS than there
are sites. A ratio of 13.9 means fourteen times. For scale, Fischer 2024 worked
at ≤ 0.004 mol/mol and estimated the ratio in human serum at ≤ 0.0005 — so
**these assays run 5,500–27,700× above the physiological ligand:protein
ratio.** A1 is not marginally violated here; it is nowhere near true.

**What that does to the "constant".** Rearranging (★6), the specific term is
just its own trace limit times a saturation penalty:

```
D_specific = C_max/(K_d + C_w) = D_trace,spec · 1 / (1 + C_w/K_d)
```

That factor `1/(1 + C_w/K_d)` **is** the activity-coefficient correction of
Approximation 1, written for a site-binding sorbent. γ constant means "the
sorbent the Nth molecule meets is the same one the first molecule met". Fill
the sites and it is not.

Between the trace limit and full saturation, `D_BSA/w` falls to the
non-specific floor:

| | K_d (µM) | D at trace | D at the floor | **collapse** |
|---|---|---|---|---|
| PFBA | 18.3 | 546 | 43 | 12.8× |
| PFOA | 0.317 | 85,100 | 5,010 | 17.0× |
| PFHxS | 0.0744 | 256,000 | 2,950 | **86.6×** |
| PFOS | 0.00624 | 1,350,000 | 69,200 | 19.6× |

**And inside the actual wells**, solving the balance self-consistently:

| study | cmpd | dose | C_free | C_free/K_d | occupancy | D / D_trace |
|---|---|---|---|---|---|---|
| Rosenmai | PFOA | 10 µM | 0.081 µM | 0.26 | 20% | 0.81× |
| Rosenmai | PFOA | 30 µM | 0.419 µM | 1.3 | 57% | 0.46× |
| Rosenmai | PFOA | 100 µM | 5.30 µM | 16.7 | **94%** | **0.11×** |
| Fragki | PFOA | 6.25 µM | 0.367 µM | 1.2 | 54% | 0.50× |
| Fragki | PFOA | 100 µM | 22.9 µM | 72.4 | **99%** | **0.07×** |

**A parameter that is supposed to be a constant moves 7-fold inside a single
dose series.** That is the whole justification for Qin's extra machinery.

**Four consequences that run through the entire comparison:**

1. **The dose exponent (§4 of the memo).** Armitage is linear by construction,
   pinned at b = 1.00. Measured b runs 0.85–1.42 and rises with chain length.
   Only a saturating model can produce that, and Qin does.
2. **Qin's error swings with dose** where Armitage's does not — his saturation
   is fitted to BSA isotherms, and the well's saturation need not match it.
3. **A saturating sink amplifies an error in its own size.** Cut the albumin
   5× and Armitage's prediction moves by exactly the denominator arithmetic,
   4.54× at every dose. Qin's moves by up to 14×, because less albumin → more
   free → higher occupancy → lower D_BSA → less capacity → more free:

   | Fragki PFOA, 10% → 2% FBS | occupancy | D_BSA ratio | **Qin x_cell** | **Armitage x_cell** |
   |---|---|---|---|---|
   | 6.25 µM | 13% → 54% | 0.56× | **7.89×** | 4.54× |
   | 25 µM | 48% → 93% | 0.23× | **14.17×** | 4.54× |
   | 100 µM | 94% → 99% | 0.63× | 4.56× | 4.54× |

   The feedback is largest **mid-range**, where the 10% well is still linear
   and the 2% well has already tipped into saturation. At the top dose both
   are saturated and it vanishes.
4. **A literature binding constant is meaningless without its ν.** An apparent
   K_a measured at ν = 3 is not the number measured at ν = 0.004 — Qin 2024
   later measured both regimes and found specific binding exceeds non-specific
   by 3.0–5.0×, with the changeover around ν ≈ 0.01 to ν ≈ 1. This is why the
   binding database in `ALBUMIN_BINDING.md` records the ligand:protein ratio
   for every row it can.

```python
import sys; sys.path.insert(0, "tk_review/scripts")
import mbm_published_equations as P, mbm_qin_vs_armitage as Q

p = Q.qin_params("PFOA", 7)
trace = p["Cmax"] / p["Kd"] + 10 ** p["logDnonspec"]
a = P.ASSAYS["Fragki"]; V = P.volumes(a)
alb = P.SI["conc_ser_alb"] * a["FBS"] / Q.BSA_MW * 1e6      # uM albumin
for d in a["doses"]:
    r = P.qin_run("PFOA", d, V)
    theta = r["Cw"] / (p["Kd"] + r["Cw"])
    print(f"{d:6g} uM   nominal PFAS:alb {d/alb:5.2f}   C_free {r['Cw']:7.3f} uM"
          f"   occupancy {theta:6.1%}   D/D_trace {r['Dbsa']/trace:5.2f}x")
```

#### 10b. The trace limit, and the bridge to the binding literature

When `C_w ≪ K_d`, (★6) collapses to a constant:

```
D_trace = C_max/K_d + D_ns = n·(ρ_alb/M_alb)·K_a + D_ns
```

Dropping the non-specific term and converting to a **mass** basis (divide by
ρ_alb):

```
log D [L_w/kg_protein] = log K_a + log n − log₁₀(M_protein/1000)          (★7)
```

For albumin (M = 66,500) the constant is **1.823**. This is the bridge that
makes the entire binding literature — reported as K_a in M⁻¹ — comparable with
a partition coefficient in L/kg. Without it, the two cannot even be compared.

> **Approximation A15.** (★7) needs `n`, the stoichiometry, which is almost
> never reported with a K_a. It moves log D by up to 0.78 across the plausible
> range.

#### 10c. Qin's cell

```
D_cell/w = D_SP/w · f_protein + D_lip/w · f_lipid       — and NO water term
```

`D_SP/w` is a *measured* non-specific protein–water ratio. The cell is 92–95%
protein at every chain length.

> **Approximation A16.** Omitting cell water means Qin's cell cannot fall below
> its protein term, which is why his PFBA cell partition (log 0.38) sits above
> Armitage's (log 0.08) — Armitage's is almost pure cell water there.

### 11. The two models, side by side

| | Armitage 2021 | Qin 2023 |
|---|---|---|
| source of `D` | predicted from log K_ow | measured, then interpolated |
| albumin | linear spLFER, piecewise | **saturating Langmuir** |
| albumin anion penalty | 0 | n/a (measured on the anion) |
| cell | lipid + 3.5% protein-as-octanol + water | protein + lipid, no water |
| cell at PFBA | 79% water | 92% protein |
| cell at PFNA | **99.7% membrane lipid** | 95% protein |
| dose dependence | linear by construction | non-linear |
| compounds parameterised | any, from log K_ow | PFBA, PFOA, PFHxS, PFOS only |

> **Approximation A17.** Qin measured four compounds. PFHxA and PFHpA are
> interpolated between the two carboxylate anchors (PFBA, PFOA) and **PFNA is
> extrapolated beyond both**. Only `log D_lip/w` has a published interpolation
> rule (Fig. A1, `0.476·nCF + 0.172`); everything else is a two-point
> interpolation with no basis in the paper.

### 12. Where the dose actually goes — run this

```python
import sys; sys.path.insert(0, "tk_review/scripts")
import mbm_published_equations as P, mbm_qin_vs_armitage as Q

a = P.ASSAYS["Rosenmai"]; V = P.volumes(a)
print(f"{'cmpd':7s}{'albumin':>10s}{'ser.lipid':>11s}{'water':>9s}"
      f"{'cells':>9s}{'plastic':>9s}{'  a (alb share)':>16s}")
for nm, _nC, nCF, _mw, ko, _kc, _s in Q.COMPOUNDS:
    lk = ko if ko is not None else Q.KOW_OPERA(nCF)
    D = P.armitage_D(lk, V["comp"])
    terms = {"albumin":   D["Dsa"] * V["Valb"],
             "ser.lipid": D["Dow"] * V["Vslip"],
             "water":     V["Vw"],
             "cells":     D["Dcw"] * V["Vcells"],
             "plastic":   D["Dpl"] * V["sarea"] * 1000.0}
    S = sum(terms.values())
    print(f"{nm:7s}" + "".join(f"{terms[k]/S:>{w}.1%}" for k, w in
          (("albumin",10),("ser.lipid",11),("water",9),("cells",9),("plastic",9)))
          + f"{terms['albumin']/S:>16.3f}")
```

Running it gives `a` = 0.196, 0.918, 0.977, 0.982, 0.971 for PFBA→PFNA, with
water at 80.3% for PFBA and 0.5% for PFNA. **The balance is a water calculation
for the shortest chain and an albumin calculation for everything else** — and
the switch happens between PFBA and PFHxA, in one step of two CF₂. That abrupt
switch is why PFBA behaves like a different compound in every analysis in this
project, and why it has so much leverage on any regression across the series.

### 13. Why the two models fail differently

Apply (★5). Over nCF 5–8, in decades per CF₂:

| | albumin | cell | **ratio** | residual, observed |
|---|---|---|---|---|
| Armitage | +0.386 | +0.769 | **+0.383** | **+0.345** |
| Qin | +0.547 | +0.515 | **−0.032** | **−0.074** |
| Fischer 2024, *measured* BSA | **+0.540** | — | — | — |

(★5) predicts the residual slope from the ratio slope, and it does — to within
the noise on a four-point fit.

**Armitage is wrong from both ends and the errors stack.** Its albumin rises at
+0.386 where Fischer measures +0.540 (A13, the piecewise flattening); its cell
rises at +0.769 where a flat ratio needs ≈ +0.540 (A14, the cell is lipid).
40% of the gap is the spLFER, 60% the cell.

**Qin gets both right** — albumin within 0.007 of Fischer's measurement — and
has essentially no chain-length residual, despite being *worse* on the level.

> **Qin has the right chain-length physics and the wrong level. Armitage has
> the better level and the wrong chain-length physics.**

---

## Part IV — The three data papers

All three measure the same quantity on the same axis: **moles found in the
cells ÷ moles added to the well**. They differ in almost everything else.

### 14. Rosenmai et al. 2018 — *J Appl Toxicol* 38:219, doi 10.1002/jat.3515

**What.** 14 PFAS in HepG2; cellular concentration by LC-MS/MS alongside a
PPARα reporter assay; the point of the paper is the *relationship* between
cellular concentration and receptor activation.

**For us.** §2.5: 12-well plates, 18 × 10⁴ cells/well, 24 h at 10/30/100 µM,
5× ice-cold PBS wash, lysates from **three wells pooled**. Table 1 is the
uptake percentage — the 15 points this project uses. Recovery 78–98%.

**Gotchas.**
- Stocks are in **ethanol**, and §2.4 holds it at **3% v/v** — the largest
  co-solvent of the three. (§6 shows why it nearly cancels.)
- The **working volume is never stated.** 1 mL/well is an assumption.
- The 3% ethanol is stated for the PPARα assay; §2.5 does not repeat it.
- Uptake is < 1% for everything except PFDoDA, PFTeDA and FOSA.

**Read:** §2.1, §2.2, §2.4 (the ethanol sentence), §2.5, §3.3 + Table 1.

### 15. Amstutz et al. 2022 — *Toxicology* 480:153312, doi 10.1016/j.tox.2022.153312

**What.** 13 PFAS in HepG2; MTT viability and DCFH reactive-oxygen assays at
0–800 µM, plus intracellular PFAS at 100 µM. The point is a structure–activity
ranking for cytotoxicity.

**For us.** §2.5: 96-well, 4 × 10⁴ cells, 250 µL dosed, 50 µL taken at t = 0
for a mass balance, 24 h on the remaining 200 µL, 5× DPBS. Exposure is in
**the culture medium — DMEM + 10% FBS**, 0.1% DMSO carrier. Explicitly "based
on the study of Rosenmai et al. (2018), with small modifications".

**Gotcha.** Amstutz reports an *intracellular fraction* normalised to
**recovered** mass, not to nominal. Its denominator excludes whatever went to
the plastic, so its numbers are not quite on the same axis as the other two.

**Read:** §2.2, §2.5. Skip MTT/DCFH entirely.

### 16. Fragki et al. 2022 — *Food Chem Toxicol* 172:113559, doi 10.1016/j.fct.2022.113559

**What.** A QIVIVE case study: PFOA, PFNA, PFHxS, PFOS in differentiated
HepaRG, with cell-associated concentrations fed through a PBK model to oral
equivalent doses. The biokinetics are a means to that end.

**For us.** The dose series in nmol/mg protein, 6.25–100 µM, Supplementary
Table 2. §2.3: 24-well, 550 µL dosed, 50 µL at t = 0, 5× DPBS. The time course
reaches steady state at **~1 h** — which is what licenses A2.

**Gotchas.**
- **§2.2: exposure is in "assay medium (growth medium containing 2% FBS)"**,
  not the 10% of the growth medium. This single sentence was worth ~4× on
  every Fragki prediction.
- Cells are differentiated in **1.7% DMSO for twelve days** and held there —
  hepatocyte-like and triglyceride-loaded, nothing like the AREc2 composition
  they are modelled with (A10, A14).
- The authors themselves **excluded PFNA at 100 µM** from their QIVIVE for low
  protein content suggesting cytotoxicity. It is the largest value in the set.

**Read:** §2.2 (last sentence), §2.3, §3.1, Supplementary Table 2.

### 17. What makes the three mutually inconsistent

| | Rosenmai | Amstutz | Fragki |
|---|---|---|---|
| cells | HepG2 | HepG2 | **differentiated HepaRG** |
| exposure FBS | 10% | 10% | **2%** |
| co-solvent | **3% ethanol** | 0.1% DMSO | 0.5% DMSO |
| plate | 12-well ×3 pooled | 96-well | 24-well |
| volume | **not stated** | 250 → 200 µL | 550 → 500 µL |
| normalised to | nominal | **recovered** | nominal |

And the measurements themselves disagree: at 100 µM PFHxA, Rosenmai reads
0.14% and Amstutz 0.039% — **3.6× apart on the same cell line and the same
axis**. No model error explains that one.

---

## Part V — The same equation, in serum

Fischer et al. 2024 (*ES&T* 58:1055, doi 10.1021/acs.est.3c07415) applies (★)
to serum instead of a well. Serum = water + albumin + globulin:

```
D_serum/w = VF_w + D_HSA/w·c_alb + D_glob/w·c_glob            [L_w / L_serum]
f_unbound = VF_w / D_serum/w
```

**Watch the units.** `D_HSA/w` is in L_water per **kg** of albumin, `c_alb` in
kg per litre of serum, so the product is L_water per litre of serum. The
protein **density cancels**: writing it as Fischer does, `D_vol·VF` with
`VF = c/ρ` and `D_vol = D_mass·ρ`, the ρ drops out. It survives only in `VF_w`.

This is the one place in the whole subject where the balance is *validated*:
26 PFAS have D_HSA/w, D_glob/w **and** D_serum/w measured separately, and
predicting the third from the first two with no fitted parameter gives a
geometric-mean bias of 1.00× and a typical error of 1.55×
(`ALBUMIN_BINDING.md` §13a). The structure works. It is the in vitro *well*
that does not.

---

## Check your understanding

**1.** In (★★), `x_cell = D_cell·V_cell/S`. If you doubled the number of cells
in the well, would the measured fraction of dose in cells double? Under what
condition?

**2.** The co-solvent result (★★★) says the effect scales with `w = V_w/S`.
Rosenmai's PFBA has `w = 0.80`. What does that tell you about which compound
you should use to *detect* a co-solvent problem experimentally?

**3.** Qin reports `n` = 0.45 binding sites for PFBA. A site count below one is
physically odd. Give two readings of it — one that is a real property of the
system, one that is an artefact of the fit.

**4.** Armitage's albumin anion penalty is `d_SaW_acid = 0` while its octanol
penalty is 3.1. What would happen to albumin's share `a` if the albumin
penalty were also 3.1, and what would that do to the model's predicted uptake?

**5.** Both models were given the AREc2 cell composition for Fragki's
differentiated HepaRG. HepaRG carry visible triglyceride droplets. Which of
Armitage's four cell terms would that change most, and would it make the
over-prediction better or worse?

**6.** (★5) says the residual's chain-length slope is the error in the slope of
`D_cell/D_alb`. Suppose someone measured a new albumin series and found a slope
of +0.70 per CF₂ instead of Fischer's +0.540. Would that make Armitage's
chain-length problem better or worse — and what would it do to Qin's?

<details>
<summary>Answers</summary>

**1.** Only if the cells are a negligible share of `S`. `x_cell = D_cell·V_cell/
(V_w + … + D_cell·V_cell)`. When the cell term is small the denominator barely
moves and `x_cell ∝ V_cell` — doubling cells doubles uptake. When cells hold a
large share, the denominator grows too and the fraction saturates toward 1. In
these assays the cell share is 0.1–18%, so the proportionality roughly holds —
which is also why the *ratio* predicted/measured is insensitive to the cell
count, and why §8 of the memo says the comparison survives not knowing it.

**2.** The opposite of what you would guess. The compound where water dominates
is where a co-solvent shows up — so you test it with **PFBA**, not PFNA. For a
long chain the effect cancels to 0.2%. Generally: to detect an error that
multiplies all the `D`s, look where the one term without a `D` is large.

**3.** (a) Real: `n` is an *apparent* mean site count over a heterogeneous
population of sites with different affinities; a weak, partly-occupied class of
sites fitted as one homogeneous class returns `n < 1`. (b) Artefact: `C_max`
and `K_d` are strongly correlated in a Langmuir fit, so if the measured
concentration range does not reach saturation, the pair is only determined
along a ridge and `n` is not separately identifiable. For PFBA `K_d` = 18.3 µM
against doses up to 100 µM, so saturation *is* sampled — which favours reading
(a).

**4.** `D_sa` would fall by 10^3.1 ≈ 1,260×, so `a` would collapse from ~0.95
to near zero and the cell and plastic would carry the well. Free concentration
would rise by up to the same factor and the model would over-predict uptake
catastrophically. This is the sense in which `d_SaW_acid = 0` is doing enormous
work: it is the assumption that makes albumin the dominant sink.

**5.** Storage lipid, `f_sl`, which is 0.00046 in AREc2 — essentially zero.
Triglyceride droplets are storage lipid, so the real `f_sl` could be orders
above it. That raises `D_cell`, raises `x_cell`, and makes the over-prediction
**worse**. Like every other correction in this project, the direction is
unhelpful — which is itself the finding.

**6.** Armitage **better**, Qin **worse**. Armitage's ratio slope is
`cell − albumin` = 0.769 − 0.386 = +0.383; a steeper true albumin slope means
Armitage's +0.386 is further below truth, so more of the residual is explained
and a corrected albumin would flatten more of it. Qin's +0.547 would now be
0.15 *below* a true +0.70, opening a ratio-slope error it does not currently
have. The point of the exercise: the ratio is what matters, so "better albumin
data" helps or hurts depending entirely on which model you are holding.

</details>

---

## Where to go next

| file | what |
|---|---|
| `report/MBM_Qin_vs_Armitage.md` | the full comparison, §15 the medium correction, §16 the chain-length residual |
| `report/ALBUMIN_BINDING.md` §12–13a | the K_d ↔ partition bridge (★7) and the serum balance that works |
| `scripts/mbm_published_equations.py` | both models, from the published equations |
| `scripts/mbm_medium_composition.py` | §6 and §7 of this primer, run against the data |
| `scripts/mbm_chain_length_residual.py` | §8 and §13 of this primer |
