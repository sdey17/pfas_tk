# Wu et al. 2024: why the HSA constants are not usable

**Wu Y, Bao J, Liu Y, Wang X, Lu X, Wang K.** In Vitro and In Silico
Analysis of the Bindings between Legacy and Novel Per- and
Polyfluoroalkyl Substances and Human Serum Albumin. *Toxics*
2024;12(1):46. doi [10.3390/toxics12010046](https://doi.org/10.3390/toxics12010046),
PMID 38251003. Open access.

Reproduce with `python3 tk_review/scripts/audit_wu2024_hsa.py`. Every
number below is arithmetic on the paper's own Tables 1–6; nothing is
refitted.

## What prompted this

On the per-compound albumin figures Wu sat about four logs below every
other human-HSA source, on all ten of its compounds at once:

| compound | Wu 2024 | other human HSA (median) | gap | n others |
|---|---|---|---|---|
| PFBA | −0.84 | 1.90 | −2.74 | 2 |
| PFBS | 0.09 | 2.98 | −2.90 | 3 |
| 6:2 Cl-PFESA | 0.11 | 3.17 | −3.06 | 2 |
| PFOS | 0.71 | 3.93 | −3.22 | 4 |
| 8:2 Cl-PFESA | 0.14 | 3.63 | −3.50 | 1 |
| HFPO-DA | −0.26 | 3.36 | −3.62 | 1 |
| PFOA | 0.03 | 3.65 | −3.62 | 8 |
| ADONA | −0.03 | 3.70 | −3.73 | 2 |
| PFHxS | 0.24 | 4.06 | −3.82 | 3 |

A *uniform* offset across ten chemically distinct compounds is the
signature of a method or an analysis, not of chemistry. Four
independent checks, all of which fail.

## Check 1 — the accessible-fluorophore fraction exceeds 1

The modified Stern–Volmer form the paper uses is

    F0 / (F0 − F) = 1 / (K'_SV · f_a · C) + 1 / f_a

so the intercept of their fitted line is `1/f_a`. HSA has **one**
tryptophan — the paper says so itself in the introduction — so `f_a`,
the fraction of fluorophore accessible to the quencher, cannot exceed
1.0. From the printed intercepts:

| | PFOA | PFOS | PFHxS | PFPrS | PFBS | PFBA | 9Cl | 11Cl | HFPO-DA | NaDONA |
|---|---|---|---|---|---|---|---|---|---|---|
| f_a | 1.57 | **4.07** | 2.11 | 1.42 | 1.66 | 1.00 | 1.41 | 1.65 | 1.24 | 1.48 |

Ten out of ten above 1. Re-deriving `Kq = K'_SV/τ0` from those slopes
reproduces the paper's printed `Kq` exactly, so this is not a
transcription error on our side — the model simply does not fit, and
the static-quenching conclusion rests on it.

## Check 2 — K_A is in M^−n, and inverting it needs n = 1

The double-logarithmic (Lineweaver–Burk) equation is

    lg[(F0 − F)/F] = lg K_A + n·lg C

Dimensionally `K_A` carries units of **M^−n**. Our Kd→K bridge
(`log K = 4.1774 + log n − log Kd(µM)`) needs a dissociation constant
in molar, i.e. it is only valid at `n = 1`. Wu's fitted n runs
**0.47 to 0.87**.

Re-express each fit as the 1:1 constant it implies *inside* the
measured window — `K_eff = K_A · C*^(n−1)` at the geometric-mean
concentration C* = 315 µM — and compare against the paper's own,
independently fitted Stern–Volmer constant:

| | K_A (M^−n) | K_eff (M^−1) | K'_SV (M^−1) | ratio |
|---|---|---|---|---|
| PFOA | 103.3 | 1352 | 1261 | 1.07 |
| PFOS | 394.6 | 1131 | 308 | 3.67 |
| PFHxS | 157.9 | 1362 | 840 | 1.62 |
| PFPrS | 79.5 | 1350 | 1466 | 0.92 |
| PFBS | 116.5 | 1323 | 1129 | 1.17 |
| PFBA | 20.4 | 1471 | 3022 | 0.49 |
| 9Cl-PF3ONS | 117.1 | 1084 | 1083 | 1.00 |
| 11Cl-PF3OUdS | 127.6 | 1320 | 1148 | 1.15 |
| HFPO-DA | 59.0 | 1249 | 1580 | 0.79 |
| NaDONA | 93.1 | 1313 | 1318 | 1.00 |

Within a factor of 2 on eight of ten, against a constant fitted from a
different equation. That agreement is what says the re-expression is
the right one.

## Check 3 — re-expressed, the between-compound result disappears

| | as ingested | re-expressed | shift |
|---|---|---|---|
| PFOA | 0.02 | 1.31 | +1.28 |
| PFOS | 0.71 | 1.23 | +0.52 |
| PFHxS | 0.24 | 1.31 | +1.07 |
| PFPrS | −0.11 | 1.31 | +1.42 |
| PFBS | 0.09 | 1.30 | +1.21 |
| PFBA | −0.84 | 1.35 | +2.19 |
| 9Cl-PF3ONS | 0.11 | 1.21 | +1.11 |
| 11Cl-PF3OUdS | 0.13 | 1.30 | +1.16 |
| HFPO-DA | −0.26 | 1.27 | +1.53 |
| NaDONA | −0.03 | 1.30 | +1.32 |

**2.2 logs of reported spread become 0.13.** The paper's headline
findings — that PFSAs bind more strongly than PFCAs, that affinity
rises with chain length, that the novel alternatives bind more weakly
than the legacy compounds — are all differences in the fitted exponent
`n`, not in affinity.

Note the correction is *in the right direction but far too small*: it
moves PFOA from 0.02 to 1.31, still about 2.3 logs below the human-HSA
pack. The dimensional error accounts for roughly a third of the gap.
The assay accounts for the rest.

## Check 4 — lg K_A is a restatement of n

`K_A` is the value of `(F0−F)/F` extrapolated to C = 1 mol/L, which is
**3.5 decades beyond** the highest concentration measured. If the
quenching actually observed at C* is similar for every compound, then
`lg K_A = y* − n·lg C*`, so `lg K_A` must be linear in `n` with slope
`−lg C* = 3.50`.

Observed: **slope 3.23, R² = 0.9921.**

And the quenching actually measured at C* is nearly identical across
the set: `(F0−F)/F` = 0.34 to 0.46, a **1.36-fold** spread — while the
reported K_A spans **19.4-fold**. Nineteen-fold of reported affinity
range out of 1.36-fold of measured signal, manufactured by
extrapolation.

## The assay

From §2.1 and §2.2: each PFAS was prepared as a 1 mM stock in
**anhydrous ethanol**, and 5 mL of 1 mM HSA was mixed with 1, 3 and
5 mL of it.

| mix | PFAS | HSA | L:P | ethanol |
|---|---|---|---|---|
| 5:1 | 167 µM | 833 µM | 0.20 | **16.7 % v/v** |
| 5:3 | 375 µM | 625 µM | 0.60 | **37.5 % v/v** |
| 5:5 | 500 µM | 500 µM | 1.00 | **50.0 % v/v** |

Against the <1 % that albumin binding work normally tolerates. Two
further consequences:

- **The protein is not held constant.** HSA falls 1.67-fold across the
  series while PFAS rises. This is a dilution series of a fixed total,
  not a titration, and `F0` is never defined against a
  concentration-matched blank. *Tested:* a pure dilution artefact would
  give n = 1.47; the fitted n is 0.47–0.87, so dilution alone does not
  account for it.
- **Three concentrations, one replicate, no error bars.** A
  two-parameter line through three points reaches R² = 0.99 for free.

## Disposition

The ten rows **stay** in `db/protein_binding.csv` — a real paper,
correctly transcribed; the database records what the literature says.
They are **held out of the albumin figures** and of any pooled
statistic, by name, via `EXCLUDED` in
`scripts/albumin_chainlength_data.py`, and are written to
`figures/albumin_chainlength_excluded.csv` so the holdout is visible
rather than silent.

## The independent check: does it clear the physiological bound?

`binding_kd_physiological_bound.py` asks a question that owes nothing to
any of the four checks above: at 600 µM albumin, what Kd is needed to
reproduce the >99% bound fraction PFAS actually shows in plasma? At
n = 6 sites that is **Kd ≲ 36 µM**.

The script did **not** catch Wu, and structurally could not: its `A62`
table is a hard-coded transcription of OEHHA 2024 Table A6.2 plus two
Starnes DSF rows, and it never opens `db/protein_binding.csv`. (It does
carry a `Wu et al. 2009` — a *different* paper, PFOA Kd 26 µM, which
passes at 99.3%. Easy to mistake for coverage.) Wu 2024 is now scored
in a separate `FLAGGED`/`WU2024` block, kept apart so the OEHHA
transcription stays a transcription:

| compound | Kd (µM) | bound, n=1 | bound, n=6 |
|---|---|---|---|
| PFOS | 2 534 | 19.1 % | 58.7 % |
| PFOA | 9 676 | 5.8 % | 27.1 % |
| PFBA | 49 102 | 1.2 % | 6.8 % |

All ten fail. The *tightest*, PFOS, is weaker than every row in Table
A6.2 except Klevens & Ellenbogen 1954; PFOA is 3.1× weaker than that
floor and PFBA 15.8×.

And the part that matters for the diagnosis: **re-expressed, it still
fails.** K_eff ≈ 1084–1471 M⁻¹ → Kd 680–922 µM → 80–84 % bound at n = 6
— the same class as the DSF thermal-transition rows, which also fail.
So the dimensional error is not the whole gap, confirming from an
independent direction what Check 3 implied: the analysis accounts for
about a third of it, the assay for the rest.

## Still to do

Wu 2024 also feeds `logk_converted_vs_measured.py`,
`make_logk_conversion_figure.py`, `kd_k_interconversion.py` and
`binding_all_data_protein_specific.py`. Those have not been
re-generated.
