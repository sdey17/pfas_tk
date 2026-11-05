# The L-FABP outlier is not the one that looked like it

Reproduce with `python3 tk_review/scripts/audit_khazaee2021_lfabp.py`.

**Khazaee M, Christie E, Cheng W, Michalsen M, Field J, Ng C.** *Toxics*
2021;9(3):45. doi [10.3390/toxics9030045](https://doi.org/10.3390/toxics9030045),
PMID 33652875. Full text via PubMed Central (PMC7996760).

**Gao K, Zhuang T, Liu X, et al.** *Environ Sci Technol* 2019;53(11):6529-38.
doi [10.1021/acs.est.9b00715](https://doi.org/10.1021/acs.est.9b00715), PMID 31099564.

## What prompted it

On the L-FABP figure Khazaee's PFOA sits at log K 5.85, more than two
logs above the next human value, from a single Kd of 0.099 µM. It
looked like the L-FABP set's Wu 2024.

## The paper's own defence is real, but does not cover this case

Khazaee §3.3 argues that equilibrium dialysis systematically returns
lower Kd than fluorescence displacement, backed by a search over 37
studies: *"the low Kd values measured here are not an artifact of this
study but rather a consistent outcome of the EqD approach."*

True, and useful. But **Gao 2019 is also equilibrium dialysis**, on the
same protein, for the same compound, and the two differ 1,990-fold:

| | method | [protein] | ligand | reported Kd | log K |
|---|---|---|---|---|---|
| Khazaee 2021 | equilibrium dialysis | 1 µM | 0.05–5 µM | **0.099 µM** | 5.85 |
| Gao 2019 | equilibrium dialysis | 6 µM | 2–10 µM | **197 µM** | 2.55 |

An EqD-versus-FD argument cannot explain a gap between two EqD studies.

## The real difference: where each titration sat

A Kd is identifiable only if the ligand concentrations bracket it.

- **Khazaee**: ligand spans **0.5× to 50×** its reported Kd.
- **Gao**: ligand spans **0.010× to 0.051×** — never within twenty-fold.

Equilibrium dialysis measures the *free* concentration, so the test is
whether changing the Kd changes that measurement by more than the
assay's precision (Khazaee states 4–18 % RSD whole-method):

| | bound at any point | free conc. shift if Kd were 10× tighter | verdict |
|---|---|---|---|
| Khazaee, PFOA | 20–91 % | up to **89 %** | identifiable |
| Khazaee, PPAR-δ PFHxS | 20–91 % | up to **89 %** | identifiable |
| Gao, PFOA | **2.9 %** | **19.7 %** | marginal — only just clears |
| Gao, PFBA | **0.7 %** | **5.6 %** | **not identifiable** |

## Conclusion, which reverses the starting assumption

**Khazaee's constants stay, unflagged.** The better-designed titration
of the two is the one that disagrees with the pack.

What the audit actually found is that **Gao 2019's L-FABP constants are
weakly identified, and progressively worse the shorter the chain** —
the titration is fixed at 2–10 µM while the reported Kd climbs from
77 µM (PFDA) to 879 µM (PFBA). At PFBA only 0.7 % of the ligand is
bound and a ten-fold change in Kd moves the measurement 5.6 %, inside
the method's own noise.

**Gao supplies the only published PFBA value for L-FABP in this
database**, and it is the least identifiable number in the set. The
short-chain end of the L-FABP figure rests on it.

## Still worth knowing about Khazaee, none of it fatal

- **Duplicates only** — five ratios, n = 2, no replication beyond that.
- **Negative binding coefficients.** Some dialysate concentrations came
  out *above* the non-binding control, which is impossible. The paper
  attributes this to 15- to 3000-fold dilutions and says a decision
  tree handled them; the tree is not in the main text.
- **The sorption check does not establish what it is cited for.**
  Recovery of **75–235 %** is quoted as showing "no level of detectable
  sorption". A window that wide cannot exclude a two-fold loss.
- **The MD screen that chose the pairs did not predict the outcomes:**
  PFBA predicted to bind PPAR-α and did not; PFBS predicted for PPAR-δ
  and did not; PFBA bound PPAR-δ unpredicted; and for I-FABP *both*
  predicted binders showed no binding. Recorded in
  `db/insilico_models.csv`.

## Disposition

No change to the figures. The flag that was going to go on Khazaee
goes on Gao 2019 instead, and the next thing worth doing is deciding
whether the short-chain L-FABP points should carry a marker for
"reported Kd lies far outside the concentration range tested".
