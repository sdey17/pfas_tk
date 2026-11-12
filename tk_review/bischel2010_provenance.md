# Where to verify each Bischel 2010 row

Bischel BJ, MacManus-Spencer LA, Luthy RG. "Noncovalent interactions of
long-chain perfluoroalkyl acids with serum albumin." *Environ Sci Technol*
2010;**44**(13):5263–5269. doi [10.1021/es101334s](https://doi.org/10.1021/es101334s)

Companion: `bischel2010_for_excel.tsv` (8 rows).

---

## The four measured rows — Table 1, page 5265

> **TABLE 1. Percent Bound and log K_PW for PFNA Binding to 500 µM Albumin
> Determined by Equilibrium Dialysis and LC-MS/MS**
>
> | albumin | conditions | percent bound | log K_PW | [ligand]:[protein] mole ratio |
> |---|---|---|---|---|
> | HSA | 50 mM sodium phosphate, pH 7.4 | 99.95 ± 0.01% | 4.93 ± 0.05 | 5 ± 1 × 10⁻⁴ |
> | BSA | 50 mM sodium phosphate, pH 7.4 | 99.92 ± 0.01% | 4.74 ± 0.05 | 6 ± 1 × 10⁻⁴ |
> | BSA | 50 mM sodium phosphate, pH 7.4; 9 g/L NaCl | 99.89 ± 0.01% | 4.56 ± 0.05 | 10 ± 1 × 10⁻⁴ |

The fourth is **not in any table**. It is one sentence of Results prose, p5265:

> "Confirming results obtained for PFNA with essentially fatty acid free BSA
> (99.93 ± 0.01% bound for a 6 ± 1 × 10⁻⁴ [L]:[P] mole ratio, corresponding to
> a log K_PW of 4.80 ± 0.08) suggest that trace fatty acids in the essentially
> γ-globulin free BSA used in this and subsequent tests had little influence on
> the binding affinities."

`Protein conc = 500` is in the Table 1 title and again in Methods p5264
("prepared with either 500 µM albumin (human or bovine)").

**`PFAS conc` — two reported, two derived.** The tables give loading only as
a mole ratio, but the Results text on p5265 states the molarities outright:

> "equilibrium dialysis was conducted with 500 µM HSA or BSA and **0.24 ±
> 0.08 µM or 0.32 ± 0.07 µM PFNA**, respectively, in 50 mM sodium phosphate
> buffer"

So HSA takes 0.24 and BSA 0.32 as reported. The fatty-acid-free control runs
at the same 6×10⁻⁴ ratio as that BSA test and takes the same 0.32. The NaCl
row is not covered by the sentence, so its 0.50 µM is ratio × 500 µM.

**`ligand:protein ratio` is then the quotient, not Table 1's rounded figure.**
Table 1 prints 5, 6 and 10 ± 1 × 10⁻⁴; the sheet shows 0.00048, 0.00064 and
0.001, which is the two concentrations divided. The two agree inside the
paper's own ±1×10⁻⁴, and the pipeline prefers the quotient wherever it does
not contradict a stated value.

### Why the sheet's log K is 0.135 below the paper's

Eq 2, p5265, defines K_PW on a **volume** basis:

> K_PW = C_P/C_W = f_bound / (F_albumin · [P] · (1 − f_bound))

with "A partial specific volume of 0.733 mL/g was used in calculations for BSA
and HSA." `F_albumin · [P]` is a volume fraction, so K_PW is L water per **litre
of protein**. This project's column is per **kg**, so log10(1/0.733) = 0.135 is
subtracted.

Check it yourself from the percent-bound column alone, at 500 µM albumin
(= 0.03323 g/mL):

| condition | recompute, volume basis | recompute, mass basis | paper prints |
|---|---|---|---|
| HSA | **4.914** | 4.779 | 4.93 |
| BSA | **4.710** | 4.575 | 4.74 |
| BSA + NaCl | **4.571** | 4.437 | 4.56 |
| BSA, fatty-acid free | **4.768** | 4.633 | 4.80 |

Volume basis to within 0.02 log; mass basis off by 0.13–0.15. So:

| sheet log K | = paper − 0.135 |
|---|---|
| 4.795 | 4.93 (HSA) |
| 4.605 | 4.74 (BSA) |
| 4.665 | 4.80 (BSA, fatty-acid free) |
| 4.425 | 4.56 (BSA + 9 g/L NaCl) |

---

## The four isotherm rows — Table 2, page 5266

> **TABLE 2. Association Constants (Ka) and Binding Stoichiometries (n) for
> PFOA and PFNA Binding to BSA Determined by Equilibrium Dialysis**
>
> | compound | Ka,1 (M⁻¹) | n1 | R² | [ligand]:[protein] range | model |
> |---|---|---|---|---|---|
> | PFOA | 0.20 (±0.14) × 10⁶ | 4.3 (±2.0) | 0.772 | 0.04–5 | one-class |
> | PFOA | 1.4 (±1.9) × 10⁶ | 1.4 (±0.5) | 0.913 | 0.04–70 | two-class |
> | PFNA | 1.1 (±0.2) × 10⁶ | 4.6 (±0.3) | 0.943 | 0.02–5 | one-class |
> | PFNA | 3.3 (±3.2) × 10⁶ | 2.9 (±0.7) | 0.953 | 0.02–120 | two-class |

`PFAS conc` again comes from the ratio, and here it is free: [P] is 1 µM,
so the mole ratio *is* the micromolar concentration. `Protein conc = 1` is
in the Results text immediately above, p5265:
"equilibrium dialysis was conducted … over a wide range of PFNA and PFOA
concentrations and 1 µM BSA; [L]:[P] mole ratios ranged from 0.02 to 120".

### How Ka becomes the sheet's two columns

`Kd (uM) = 1e6 / Ka`, then `log K = 4.1774 + log10(n) − log10(Kd)`,
where 4.1774 = log10(1e9/66463) carries albumin's molar mass.

| Table 2 Ka | n | → Kd (µM) | → log K |
|---|---|---|---|
| 0.20 × 10⁶ | 4.3 | 5 | 4.112 |
| 1.4 × 10⁶ | 1.4 | 0.7143 | 4.470 |
| 1.1 × 10⁶ | 4.6 | 0.9091 | 4.882 |
| 3.3 × 10⁶ | 2.9 | 0.303 | 5.158 |

---

## What is deliberately left out

**The nanoESI-MS results** (Figure 2 and Results, p5266–5267; abstract: "BSA
complexes with up to eight occupied binding sites"). In the database but not in
the sheet: they carry no Kd or Ka, so nothing converts. Table 3 p5266 lists a
nanoESI-MS Ka of 1.3 × 10⁵ M⁻¹ for PFOA — but that is Bischel's own summary
table comparing methods, and the gas-phase number is not an aqueous binding
constant.

**SI Tables S2, S3 and S4** — seventeen refits of the same two datasets Table 2
already contributes. Reported by `scripts/bischel2010_refit_spread.py`, never
ingested: counting one experiment seventeen times would inflate the spread this
project is trying to explain.
