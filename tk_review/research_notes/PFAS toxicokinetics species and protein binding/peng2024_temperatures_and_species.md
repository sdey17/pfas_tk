# Two findings: Peng's "outlier" was ours, and the species gap is a lab gap

## 1. Peng 2024 — one experiment at two temperatures

**Peng Y, et al.** *Toxics* 2024;12(1):43.
doi [10.3390/toxics12010043](https://doi.org/10.3390/toxics12010043), PMID 38250999.
Fix: `scripts/fix_peng2024_temperatures.py`.

Peng contributed two human-albumin values per compound — for PFNA, log K
5.07 and 3.21 — and the selection rule was treating them as two
independent candidates. They are the same binding constant at two
temperatures.

**First hypothesis, tested and rejected.** It looked like Wu 2024: a
double-log K_A in M^−n alongside a Stern–Volmer constant in M^−1.
Re-expressing the double-log value inside the measured window
(`K_eff = K_A·C*^(n−1)`, C* = 6 µM) made the gap **worse** on four of
five compounds — PFNA 1.90 → 2.42 log, PFOA 0.94 → 1.51. Not the Wu
case; the dimensional story does not apply.

**What it actually is.** Peng's Table 2 reports Kb and n at 298, 304 and
310 K. Table 3's "no probe" row is the *same 298 K number* re-listed as
the baseline for the site-competition experiment — identical to four
significant figures for PFNA, PFOA, PFO3DA, PFHpA and DFSA. Not a
second measurement.

**Three faults, all in our ingestion, none in the paper:**

1. The temperature column carried all three values as one lump,
   `"25, 31, 37 (298, 304, 310 K)"`, on every row — so neither row said
   which temperature it was.
2. The 298 K row **lost its n**. Table 2 gives n at every temperature;
   without it the conversion assumed n = 1, wrong by log(1.51) = 0.18
   for PFNA.
3. Nothing recorded that the two rows are one experiment, so the
   selection could choose between them as independent evidence.

**Effect on PFNA's admissible spread:** 509-fold → 768-fold (once the
missing n was restored) → **134-fold** (once the 25 °C row was ruled
inadmissible where the same source has a 37 °C one).

**A caution that survives the fix.** The temperature dependence is
extreme and not consistent in direction: PFNA falls 79-fold from 298 to
310 K while HFPO-TA rises 624-fold and PFHpA 165-fold. The van't Hoff
enthalpies Peng derives reach **+412 kJ/mol** (HFPO-TA) and +325
(PFHpA). Non-covalent protein–ligand binding is normally within about
±50 kJ/mol; 412 kJ/mol is covalent-bond territory. The numbers are
transcribed as reported and the temperatures are now right, but this
series should not be read as a reliable thermodynamic result.

## 2. The human/bovine difference is a laboratory difference

The selection gave human PFOA fu = 0.00083 and bovine 0.00157, about
2-fold apart — but human came from Fischer and bovine from Allendorf,
so species was confounded with laboratory.

**Fischer 2024 measured 27 compounds on both human and bovine albumin,
same laboratory, same method.** That removes the confound:

| | human − bovine, log |
|---|---|
| median | **−0.02** |
| mean | +0.08 |
| sd | 0.27 |
| median absolute | **0.09** → fu differs **1.23×** |

The five compounds under comparison:

| PFAS | human | bovine | diff | fu human | fu bovine |
|---|---|---|---|---|---|
| PFHxA | 3.20 | 2.83 | +0.37 | 0.01553 | 0.03566 |
| PFHpA | 4.10 | 3.93 | +0.17 | 0.00198 | 0.00293 |
| PFOA | 4.48 | 4.20 | +0.28 | 0.00083 | 0.00157 |
| PFNA | 4.49 | 4.54 | −0.05 | 0.00081 | 0.00072 |

**Within one laboratory, human and bovine albumin are nearly
indistinguishable for PFAS binding — 1.2× in fu.** Against a
between-laboratory spread of 1.24 log on the same compound by the same
method, which is **17× in fu**.

So the between-laboratory term is about **fourteen times larger than
the between-species term**. The 2-fold human/bovine gap in the selected
values is Fischer-versus-Allendorf, not human-versus-bovine.

**This does not license pooling.** The species separation in the figures
and the selection rule stays: the finding is that for *albumin binding*
the two proteins behave alike in the one study that measured both, not
that species can be mixed in general. Nothing here speaks to rat
albumin, to L-FABP, or to the renal transporters where the real
species differences in PFAS toxicokinetics live.
