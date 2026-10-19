# PFAS binding to serum albumin: a standalone primer and review

**What this is.** A self-contained introduction to PFAS–albumin binding, written
for someone starting from zero, followed by a review of what the literature
actually establishes and what our own analysis added. It assumes no prior
knowledge of binding chemistry.

**Why it exists separately.** The albumin material in `REPORT.md` (§6.1–6.3) is
written as evidence inside a larger argument about species differences. This
document builds the subject from the ground up instead.

Every code block runs as-is on plain Python 3 — no imports from this
repository, no packages to install. Copy and run them.

**Contents**

- Part I — [Foundations](#part-i--foundations): §1–4
- Part II — [Why the literature disagrees](#part-ii--why-the-literature-disagrees): §5–8
- Part III — [What is actually known](#part-iii--what-is-actually-known): §9–11a
- Part IV — [Using a binding constant](#part-iv--using-a-binding-constant): §12–14
- [Check your understanding](#check-your-understanding)
- [Reading list](#reading-list)

---

## Part I — Foundations

### 1. Why anyone cares

The chain of reasoning that makes albumin binding a toxicokinetics question,
in four steps:

1. PFAS circulate in blood. Most of what is in blood is **not floating free** —
   it is stuck to proteins.
2. Only the **free** (unbound) fraction can be filtered by the kidney, cross a
   membrane, or reach a receptor. Bound PFAS is, for that moment, inert and
   un-excretable.
3. So the unbound fraction, written **f_u**, sets how fast the body can clear
   the compound. A smaller f_u means slower clearance and a longer half-life.
4. Albumin is the most abundant protein in plasma and the one most PFAS stick
   to. So "how tightly does PFAS X bind albumin?" is upstream of "why does PFOA
   have a 3.8-year half-life in humans?"

That is the motivation. Whether the chain actually holds up quantitatively is
§13, and the answer is "partly".

### 2. What albumin is

Human serum albumin (HSA) is a single protein, 66,463 g/mol, about 585 amino
acids, shaped like a heart in three domains. Facts worth carrying:

- **It is abundant.** Plasma holds 35–50 g/L of it. At 66,463 g/mol that is
  **526–752 µmol/L**; 600 µM is the usual working figure. Nothing else in
  plasma comes close.
- **Its day job is carrying fatty acids.** Fatty acids are oily and do not
  dissolve in blood; albumin ferries them. It has **seven** recognised
  fatty-acid sites, plus two classic drug sites (Sudlow I and II), plus a
  single tryptophan (Trp214) whose fluorescence makes a convenient probe.
- **PFAS look like fatty acids.** A perfluorinated chain with an acid head is
  structurally a fatty acid with fluorines where hydrogens would be. It binds
  in the same pockets. That is the whole mechanism in one sentence.
- **Different species have their own version** — bovine (BSA), rat (RSA),
  porcine (PSA). They are similar but not identical, which matters when a
  study uses BSA as a cheap stand-in for HSA (§10).

The key structural consequence: **albumin does not have "a" binding site.** It
has several, of different affinity. Most reported numbers lump them into one
effective constant, and a lot of the apparent disagreement in the literature
is really disagreement about which site got measured (§8).

### 3. The vocabulary, built from one reaction

Everything starts with one reversible reaction — ligand L (the PFAS) meets
protein P (albumin) and sticks:

```
L + P  ⇌  LP
```

At equilibrium the rates balance, and the ratio of concentrations is fixed:

**K_d — the dissociation constant**

```
K_d = [L][P] / [LP]          units: concentration, usually µM or M
```

K_d is the **ligand concentration at which half the sites are occupied**. That
is the intuition to keep. It has units of concentration because it *is* a
concentration.

- **Small K_d = tight binding** (it takes very little ligand to fill the sites)
- **Large K_d = weak binding**

**K_a — the association constant** is simply its reciprocal:

```
K_a = 1 / K_d                units: M⁻¹  (per molar)
```

**Large K_a = tight binding.** Papers use one or the other more or less
arbitrarily, which is the first thing that makes the literature hard to read.
A paper reporting K_a = 10⁶ M⁻¹ and one reporting K_d = 1 µM are saying the
same thing.

```python
# Converting between them. Run this.
def ka_from_kd(kd_uM):
    """Kd in micromolar -> Ka in M^-1."""
    return 1.0 / (kd_uM * 1e-6)

def kd_from_ka(ka_M):
    """Ka in M^-1 -> Kd in micromolar."""
    return 1.0 / ka_M * 1e6

print(f"Kd =  1 uM   -> Ka = {ka_from_kd(1.0):.3g} M^-1")
print(f"Kd = 100 uM  -> Ka = {ka_from_kd(100.0):.3g} M^-1")
print(f"Ka = 2.7e5   -> Kd = {kd_from_ka(2.7e5):.3g} uM")
```

**n — the number of binding sites**, sites per protein molecule. Usually fitted
alongside K_d, often badly. A value below 1 means the fit found less than one
site's worth of capacity, which tells you it is an *effective lumped
parameter*, not a count of pockets.

**f_u — the unbound fraction**, the number TK actually wants: what share of the
total PFAS in plasma is free.

### 4. Three views of the same number

This is the part worth slowing down on, because the same physical fact gets
written three ways and people talk past each other.

#### View 1 — occupancy (the chemist's view)

What fraction of the *protein's sites* are filled?

```
occupancy = [L_free] / (K_d + [L_free])
```

At physiological PFAS levels this number is tiny — albumin is almost entirely
empty.

#### View 2 — the unbound fraction (the toxicokineticist's view)

What fraction of the *ligand* is free? Starting from K_d = [L][P]/[LP] and
substituting, the free-ligand term cancels and you get something clean:

```
f_bound = [P] / ([P] + K_d)        and     f_u = 1 - f_bound
```

where [P] is the concentration of *binding sites* = n × [albumin].

**The cancellation only works when there is far more protein than ligand**
([L] ≪ [P]). In plasma that holds easily. In most laboratory assays it does
not, and that is §6.

```python
# View 2. How tight must binding be to explain what we see in blood?
ALB_uM = 600.0          # albumin in human plasma, 40 g/L / 66463 g/mol

def f_bound(kd_uM, n=1, alb_uM=ALB_uM):
    return (n * alb_uM) / (n * alb_uM + kd_uM)

def kd_needed(target_bound, n=1, alb_uM=ALB_uM):
    """Kd that would produce a given bound fraction."""
    return n * alb_uM * (1.0 / target_bound - 1.0)

print("Kd required to reach a given bound fraction:")
for f in (0.95, 0.99, 0.999):
    print(f"  {100*f:5.1f}% bound ->  n=1: {kd_needed(f,1):7.1f} uM"
          f"   n=6: {kd_needed(f,6):7.1f} uM")

print("\nMeasured PFOA (Chen 2025, Kd = 2.57 uM):")
for n in (1, 3, 6):
    print(f"  n={n}: {100*f_bound(2.57, n):.3f}% bound")
```

PFOA is measured as >99% bound in human plasma. Run the block: that requires
K_d below about 6 µM at one site. **Any reported constant much weaker than
that cannot reproduce what is seen in blood** — which is a free, powerful
filter on a messy literature.

#### View 3 — the partition coefficient (the modeller's view)

A mass-balance or PBPK model does not want an affinity, it wants a **capacity**:
how many litres of water's worth of PFAS one kilogram of albumin can hold. That
is a partition coefficient, written **D** (or K_protein/water):

```
D = (amount per kg protein) / (amount per L water)
```

Views 1–3 are the same measurement in different clothes. §12 gives the exact
bridge between them.

> **Pause and check.** Before reading on: if K_d for compound A is 1 µM and for
> compound B is 100 µM, which binds more tightly, and by how much? And which
> would you expect to be cleared from blood faster, all else equal?
>
> <sub>A binds 100× more tightly. B should clear faster, because more of it is
> free at any moment. "All else equal" is doing a lot of work — §13.</sub>

---

## Part II — Why the literature disagrees

### 5. How people actually measure this

Six methods dominate. They do not measure quite the same thing, and knowing
which was used tells you most of what you need about a number's reliability.

| method | what it does | what it really gives | watch out for |
|---|---|---|---|
| **Equilibrium dialysis** | separates free from bound across a membrane | the reference standard; a true K_d | slow; PFAS stick to the apparatus |
| **Solid-phase microextraction (SPME)** | a coated fibre samples only free ligand | free concentration over a wide range | needs careful calibration |
| **Fluorescence quenching** | watches Trp214's glow dim as ligand binds | binding *near that one tryptophan* | reports only the site it can see |
| **Probe displacement** | measures how much ligand it takes to evict a dye from a known pocket | site-resolved constants | depends on the probe's own affinity |
| **Differential scanning fluorimetry (DSF)** | measures how much the ligand stabilises the protein against melting | an apparent affinity **at the melting temperature**, ~60–70 °C | not a 37 °C equilibrium constant |
| **Ultrafiltration / nanoESI-MS** | physically separates, or counts complexes | varies | non-specific losses |

### 6. The spread, and what causes it

For a single PFAS–protein pair the reported association constants span
**four to five orders of magnitude**:

| chemical | protein | n studies | K_a range (M⁻¹) | spread |
|---|---|---|---|---|
| PFBS | HSA | 5 | 116 – 6.5 × 10⁶ | 55,800× |
| PFBA | HSA | 4 | 20.4 – 1.1 × 10⁶ | 54,000× |
| PFOS | HSA | 6 | 395 – 7.6 × 10⁶ | 19,300× |
| PFOA | HSA | 9 | 103 – 1.98 × 10⁶ | 19,200× |

The natural guess is "different methods". **That is not the main cause** —
equilibrium dialysis alone spans 3.2 × 10² to ~10⁶, and fluorescence quenching
spans 10² to 2 × 10⁶ by itself. The real drivers are four things, and all four
are knowable from a paper's methods section.

### 7. The four axes that move a binding constant

#### Axis 1 — the ligand:protein ratio (the big one)

Write **ν** for the PFAS:protein molar ratio.

| situation | [PFAS] | ν |
|---|---|---|
| human background, 4 ng/mL | 0.0097 µM | 1.6 × 10⁻⁵ |
| high-exposure community | 0.48 µM | 8 × 10⁻⁴ |
| fluorochemical worker | 12 µM | 0.02 |
| **typical DSF assay** | 1000 µM | **1.7** |
| **typical fluorescence assay** | 80 µM vs 2 µM HSA | **40** |

Laboratory assays run **three to six orders of magnitude** above the body.

Why it matters: at low ν, ligand finds the few tight ("specific") sites. At
high ν those fill up and the rest smears onto weak, non-saturable ("non-specific")
surface. The same protein and the same compound give **two different constants**
depending on which regime you are in.

This used to be an inference. **Qin et al. 2024 measured both, separately:**

| compound | specific log D | non-specific log D | ratio |
|---|---|---|---|
| PFBA | 2.44 | 1.97 | 3.0× |
| PFHxA | 3.26 | 2.69 | 3.7× |
| PFHpA | 4.12 | 3.47 | 4.5× |
| PFOA | 4.58 | 3.88 | 5.0× |

Specific binding dominates below ν ≈ 0.01 (the body); non-specific takes over
above ν ≈ 1 (the assay). **A constant carried from assay to body without asking
which regime produced it is wrong by 3–5×**, and that is the gap within one
method in one laboratory — a floor, not a ceiling.

```python
# Axis 1, made concrete. Saturation drives the apparent constant down.
def apparent_D(D_trace, Cw_uM, kd_uM):
    """A partition coefficient measured at free concentration Cw."""
    return D_trace / (1.0 + Cw_uM / kd_uM)

D0, kd = 1e5, 2.57          # trace-limit D, and PFOA's Kd in uM
print("free PFAS (uM)   apparent D   fraction of the true value")
for cw in (0.01, 0.1, 1, 10, 100, 1000):
    d = apparent_D(D0, cw, kd)
    print(f"{cw:>12}   {d:10.0f}   {d/D0:.3f}")
```

#### Axis 2 — temperature

DSF reads affinity off a **melting transition**, around 60–70 °C. Its own
authors say so. DSF constants rank congeners usefully and should not be carried
into a kinetic model as 37 °C affinities.

#### Axis 3 — ionic strength

Gao et al. 2025 is the only chain-length series measured across added salt.
Affinity falls **monotonically** with ionic strength, by roughly **1.5–2×**
between a 20 mM buffer and the ~150 mM of plasma, and the shift is independent
of chain length and head group. A constant measured in dilute buffer overstates
binding before any other correction.

#### Axis 4 — which site

See §8. This one is routinely invisible because most papers report one number
without saying which pocket produced it.

### 8. The site problem

Chen & Guo 2009 measured five PFAS at **three** sites with three probes. The
full grid, in K_a (M⁻¹):

| compound | Trp214 | Sudlow I | Sudlow II |
|---|---|---|---|
| PFBA | no binding | 1.1 × 10⁶ | no binding |
| PFBS | no binding | 2.2 × 10⁶ | 6.5 × 10⁶ |
| PFOA | **2.7 × 10⁵** | binds, not quantified | binds, not quantified |
| PFOS | **2.2 × 10⁴** | binds, not quantified | 7.6 × 10⁶ |
| PFDoA | binds, not quantified | no binding | 1.2 × 10⁶ |

Read across a row: **PFOS binds Sudlow II 345× more tightly than Trp214.** Same
paper, same protein, same compound.

The bold cells are the ones regulatory summaries carry — and the authors
themselves identify the tryptophan site as **the weakest of the three**.
Quoting it makes PFOS look like a weak binder (K_d 45 µM) when its best site
gives 0.13 µM.

**A K_d quoted without its site is as incomplete as one quoted without its n.**

The same paper also supplies the control the field usually lacks. For the same
three probes the authors measured plain, non-fluorinated fatty acids: octanoic
acid binds Sudlow II at 1.1 × 10⁷ M⁻¹ — **tighter than every PFAS in the grid**.
PFAS occupy the fatty-acid pockets without outcompeting the fatty acids, which
is why the "PFAS displace endogenous ligands" story needs care.

> **Pause and check.** A paper reports "K_a = 2.2 × 10⁴ M⁻¹ for PFOS–HSA by
> fluorescence quenching". Name three things you would want to know before
> using it, and say roughly what each could do to the number.
>
> <sub>Which site (up to 345×), the ligand:protein ratio (3–5×, or more), the
> buffer ionic strength (1.5–2×). Temperature too, if the method were DSF.</sub>

---

## Part III — What is actually known

### 9. Structure: chain length and head group

**Chain length is not monotonic.** This is the finding that most surprises
people, and three independent datasets now agree on it.

Chen et al. 2025, equilibrium dialysis, HSA at 1 µM — the cleanest series
available, K_d in µM (**lower = tighter**):

| perfluorinated carbons | PFCA | K_d (µM) | | PFSA | K_d (µM) |
|---|---|---|---|---|---|
| 3 | PFBA | 84.9 | | | |
| 4 | PFPeA | 16.3 | | PFBS | 15.6 |
| 5 | PFHxA | 14.0 | | PFPeS | 17.2 |
| 6 | PFHpA | 6.2 | | PFHxS | **1.31** |
| 7 | **PFOA** | **2.57** | | PFHpS | 1.68 |
| 8 | PFNA | 3.75 | | **PFOS** | **0.94** |
| 9 | PFDA | 3.64 | | PFNS | 6.56 |
| 10 | PFUnDA | 4.73 | | PFDS | 1.88 |
| 11 | PFDoDA | 15.8 | | | |
| 12 | PFTrDA | 31.4 | | PFDoS | 1.29 |
| 13 | PFTeDA | 13.1 | | | |

Binding tightens from PFBA to about PFOA, then **weakens again**. The usual
explanation: beyond a certain length the chain no longer fits the pocket and
has to bend.

Two consequences:

- **A structure–activity model linear in carbon number is wrong at both ends.**
- The compounds that bind albumin most tightly (C6–C10) are not the longest
  ones. Long-chain persistence is not an albumin story.

**The head group is worth about a log unit.** At matched chain length the
sulfonate binds harder than the carboxylate — PFHxS 1.31 vs PFHpA 6.2 µM at
six fluorinated carbons; PFOS 0.94 vs PFNA 3.75 at eight.

**Albumin stops being the main carrier at long chains.** Fischer et al. 2024
measured γ-globulin alongside albumin: globulin binds 2–3 log units weaker at
short chains, catches up steadily, and for PFTeDA binds *more* tightly than
albumin. Treating "protein binding" as "albumin binding" degrades exactly where
the persistent compounds live.

### 10. Species: albumin is conserved, and that is the problem

Starnes et al. 2024 measured four albumins by one method, K_d in mM
(lower = tighter):

| compound | human | rat | bovine | porcine | rat/human |
|---|---|---|---|---|---|
| PFOA | 0.79 | 0.40 | 0.99 | 0.37 | **0.51×** |
| PFOS | 0.69 | 0.42 | 0.86 | 0.38 | 0.61× |
| PFBA | 2.64 | 1.90 | 3.06 | 1.28 | 0.72× |

The differences are statistically significant and **biologically trivial**:
highest-to-lowest never exceeds ~4×, usually under 2×.

Now the punchline. **Rat albumin binds PFOA about twice as tightly as human
albumin — and the rat clears PFOA roughly a thousand times faster.** The
direction is wrong and the magnitude is off by two and a half orders of
magnitude. Albumin affinity cannot be the cause of the species difference in
elimination.

Three independent matrices now agree:

- **Albumin** (Starnes 2024, Fischer 2024): ≤4× across four species.
- **Whole tissue** (Ryu 2024, 16 PFAS × 5 tissues × 4 organisms): log f_unbound
  correlates across every species and tissue pair at R² ≥ 0.97.
- **Whole plasma** (Qin 2024): the human-versus-trout difference is explained
  by **composition** — human plasma is 42.25 mL/L protein and 4.46 mL/L lipid,
  trout 15.46 and 7.08 — not by protein identity. BSA is a good surrogate for
  human plasma.

**Species differ in what their plasma and tissues are made of, not in how their
proteins bind.**

> **A gap worth knowing about.** There is **no in vitro measurement of PFAS
> binding to mouse serum albumin** — not in our 865-row database, not in
> PubMed. The mouse is the animal most of the toxicology runs on.

### 11. The best current numbers

If you need one dataset, use **Chen et al. 2025 Table S5**: 60 PFAS, equilibrium
dialysis, HSA fixed at 1.0 µM with ligand at 0.0625–2 µM, so ν = 0.06–2. That
is one to two orders of magnitude closer to physiological than the DSF and
fluorescence datasets, and it is the widest coverage at that ratio.

The sanity check that matters: **56 of its 60 compounds reproduce >99% plasma
binding** at six sites. The whole dataset lands in the physiologically possible
band, which is not true of the DSF or 1954 surfactant literature.

For comparison, scoring the older literature against the same requirement:

| source | K_d (µM) | predicted bound (n=6) | |
|---|---|---|---|
| Chen 2025 (PFOS) | 0.94 | 100% | ✓ |
| Chen 2025 (PFOA) | 2.57 | 99.9% | ✓ |
| Chen & Guo 2009 (PFOS, Sudlow II) | 0.13 | 100% | ✓ |
| Bischel 2010 | 0.7–5 | ~100% | ✓ |
| Qin 2010 | 21 | 99.4% | ✓ |
| Beesoon & Martin 2015 (n-PFOA) | 100 | 97.3% | borderline |
| Han 2003 | 400 | 90.0% | ✗ |
| Starnes 2024 (DSF) | 790 | 82.0% | ✗ |
| Klevens & Ellenbogen 1954 | 3100 | 53.7% | ✗ |

### 11a. The two best datasets appear to disagree. They do not.

Chen 2025 and Fischer 2024 are the two cleanest albumin datasets here, both at
defensible ligand:protein ratios. Converting Chen's PFOA (K_d = 2.57 µM) at one
site gives log D = 3.77. Fischer measured **4.48**. About 5× apart — which
looked, at first, like a real problem.

It is not a disagreement. **They measure the two different binding regimes**,
and each lands on the matching term of the one study that measured both.

Putting all four on the same basis (mass, L_water/kg_protein; Qin shifted
−0.135 from his volume basis):

| compound | Fischer | Qin **specific** | gap | | Chen | Qin **non-specific** | gap |
|---|---|---|---|---|---|---|---|
| PFHxA | 3.20 | 3.13 | +0.07 | | 3.03 | 2.56 | +0.48 |
| PFHpA | 4.10 | 3.99 | +0.11 | | 3.39 | 3.34 | **+0.05** |
| PFOA | 4.48 | 4.45 | **+0.03** | | 3.77 | 3.75 | **+0.02** |
| PFHxS | 4.77 | 4.89 | −0.12 | | 4.06 | 3.45 | +0.62 |
| PFOS | 4.62 | 5.14 | −0.52 | | 4.20 | 4.04 | +0.17 |

For PFOA the two matches are **0.03 and 0.02 log units** — 7% and 5%. Three
laboratories, three methods, two regimes, and each pair lands on the other.

**Why it is not a stoichiometry artefact.** The competing explanation was that
Chen's conversion assumed n = 1 when it is really ~5, which would close the gap
arithmetically (§12: log D = log K_a + log n − 1.823). But Chen *currently* sits
on Qin's non-specific term to 0.02 log units. Raising his n to 5 would lift him
onto Fischer and in doing so **break an agreement that presently holds to 5%**.
You cannot have both. And n ≈ 1 is what §12's backwards conversion found
independently (n_implied = 1.01 for PFOA, 1.12 for PFHpA).

**So which do you use?** It depends entirely on the question:

- **A PBTK or plasma-binding calculation → the specific constant** (Fischer, or
  Qin's specific term). The body runs at ν < 0.01.
- **An in vitro mass-balance model of a dosed assay → the non-specific
  constant** (Chen, or Qin's non-specific term). The assay runs at ν > 1.

Getting it backwards costs about 5× for carboxylates and more for sulfonates.

**Two honest caveats.** The sulfonates reconcile less cleanly — Qin's own
specific-to-non-specific separation is much larger for them (PFHxS 27×, PFOS
13×) than the Fischer–Chen gap, so something else is also moving, most likely
that a single-site fit over ν = 0.06–2 cannot resolve a very tight sulfonate
site. And **Chen's PFNA breaks his own series**: his K_d puts PFNA (3.75 µM)
*weaker* than PFOA (2.57 µM), where Fischer and Qin both put it at or above.
That is the single largest residual in the table (−0.87 log units).

Reproduce: `python3 scripts/chen_fischer_reconciliation.py`

---

## Part IV — Using a binding constant

### 12. The bridge from K_d to a partition coefficient

A model wants capacity, not affinity. Here is the exact conversion.

Albumin at total molar concentration P, with n equivalent sites of dissociation
constant K_d. Langmuir occupancy gives bound ligand per litre of solution:

```
C_bound = n · P · C_w / (K_d + C_w)
```

A partition coefficient is bound ligand per unit **sorbent** over free ligand
per unit water. On a mass basis the sorbent is P·M_alb/1000 kg protein per
litre, so:

```
D = [C_bound / (P·M_alb/1000)] / C_w  =  1000·n / [M_alb·(K_d + C_w)]
```

**P cancels.** The partition coefficient does not depend on how much albumin is
in the tube — only on how many sites each albumin carries and how tightly they
hold. In the trace limit (C_w ≪ K_d):

> **log D = log K_a + log n − 1.823**  (mass basis, L_water/kg_protein)
>
> **log D = log K_a + log n − 1.688**  (volume basis, L_water/L_protein)

```python
import math
M_ALB = 66463.0          # g/mol

def logD_from_ka(ka_M, n=1, basis="mass"):
    shift = math.log10(M_ALB / 1000.0)
    if basis == "volume":
        shift -= math.log10(1.0 / 0.733)      # albumin density 1.364 kg/L
    return math.log10(ka_M) + math.log10(n) - shift

ka = 1.0 / (2.57e-6)      # Chen 2025 PFOA
for n in (1, 3, 6):
    print(f"n={n}: log D = {logD_from_ka(ka, n):.2f} (mass) "
          f"{logD_from_ka(ka, n, 'volume'):.2f} (volume)")
```

**Three groups write this same relation.** Our derivation; Qin's mass-balance
model, whose albumin term `C_max/(K_d + C_w) + D_ns` is this expression (its two
C_max columns recover the pure-albumin molar concentration to 0.3%); and Chen
2025's Equation 4, published independently. It is not a modelling choice — it is
the definition of a partition coefficient applied to a saturable site.

**What the conversion needs that affinity papers rarely report.** K_a and D are
strictly proportional, log-log slope exactly 1. Everything a model needs is in
the **intercept**, which is set by the molar-mass basis and by **n**. Moving n
from 1 to 6 shifts log D by 0.78 with no new measurement.

**Run it backwards and it becomes a diagnostic.** Ask what stoichiometry would
reconcile the tightest reported K_a with Fischer's directly measured log D:
**n = 1.01 for PFOA, 1.12 for PFHpA.** Two completely different kinds of
measurement agree at one binding site, with nothing fitted.

### 13. What happened when we used it

Short version: **albumin binding turned out not to explain what we hoped.**

Two published mass-balance models (Qin 2023, Armitage 2021) predict how much
PFAS gets into cells in an in vitro assay. Both over-predict against 29
measured points. Albumin carries 19–94% of the sorptive capacity in those
assays, so it looked like the obvious culprit.

We substituted Fischer's measured albumin partition coefficient — the only one
measured at a physiological ligand:protein ratio — into both models. **Both got
worse** (Qin 3.43× → 4.74×, Armitage 2.33× → 5.09×). Fischer's value is *lower*
than the models assume, which shrinks the albumin sink and pushes more dose into
the cell, and the models already over-predicted.

Inverting the balance for the albumin coefficient that *would* fix them gives
values **above** anything measured. The data and the residual point in opposite
directions. **Albumin is where the sensitivity is and not where the error is.**

The same thing happened with cell protein, from the other direction. Tracking
it down (`MBM_Qin_vs_Armitage.md` §12–14) ended at the **cell volume**, which is
uncertain to 3.1× and is the only parameter that rescales every compound by the
same factor.

One result from that work is worth keeping regardless of the models:
**Armitage's convention that 3.5% of cell protein behaves as octanol gives a
protein affinity 457–7,820× below the measured value.** A coefficient calibrated
on neutral organics does not describe an anionic surfactant on protein.

### 14. What is still open

1. **Mouse albumin has never been measured.** The largest single gap.
2. **n is almost never reported** with a K_d, and it moves log D by up to 0.78.
   Where it has been checked, n ≈ 1 is what makes the datasets consistent
   (§11a, §12) — but that is three compounds, not a rule.
3. **Nobody has measured a PFAS–albumin constant at a genuinely physiological
   ν** (< 0.01) across a wide compound set. Chen's 0.06–2 is the closest.
4. **The sulfonates do not reconcile as cleanly as the carboxylates** (§11a),
   and nobody has explained why their specific-to-non-specific separation is
   so much larger.
5. **Chen's PFNA inversion** (§11a) is unexplained and is his largest residual
   against the other two datasets.

*Resolved since the first draft:* the apparent Chen-versus-Fischer conflict,
which turned out to be the specific/non-specific regime difference (§11a).

---

## Check your understanding

Work these before looking at the answers. They are ordered by difficulty.

**1.** A study reports K_a = 1 × 10⁵ M⁻¹ for PFHxA–HSA. What is K_d in µM?

**2.** Using the 600 µM albumin figure and one binding site, what bound fraction
does that K_a predict? Is it consistent with PFHxA being highly bound in plasma?

**3.** Two papers report PFOS–HSA constants: 2.2 × 10⁴ M⁻¹ and 7.6 × 10⁶ M⁻¹.
They are the same paper. What explains the 345× difference?

**4.** A DSF study reports K_d = 790 µM for PFOA. Why should this not go into a
PBPK model, and what is it still good for?

**5.** Convert Chen 2025's PFOS (K_d = 0.94 µM) to log D on the mass basis at
n = 1 and n = 6. Which would you use for a mass-balance model, and what would
you need to know to decide?

**6.** Rat albumin binds PFOA ~2× more tightly than human albumin. The rat
clears PFOA ~1000× faster. Reconcile these.

**7.** Why does substituting a *better* albumin measurement make a mass-balance
model *worse*?

<details>
<summary>Answers</summary>

**1.** K_d = 1/10⁵ M = 10 µM. (`kd_from_ka(1e5)`)

**2.** f_bound = 600/(600+10) = 98.4%. Consistent with "highly bound" but it
would not reach >99% at one site — you would need n ≈ 2 or more.

**3.** Different **sites** on the same protein. 2.2 × 10⁴ is the tryptophan
site, which the authors call the weakest of three; 7.6 × 10⁶ is Sudlow site II.
Neither is wrong; quoting one without naming the site is.

**4.** DSF measures affinity at the protein's melting transition (~60–70 °C),
not at 37 °C, and 790 µM predicts only 82% bound at n=6 where PFOA is observed
>99% bound. It is still good for **ranking** congeners and comparing species by
one consistent method.

**5.** log D = log(1/0.94e-6) + log n − 1.823 → 4.20 at n=1, 4.98 at n=6. To
choose you need the stoichiometry, which Chen fits as B_max per compound — and
you need to know whether your model's albumin term is on a mass or volume basis
(a further 0.135 log units).

**6.** Binding affinity is not the rate-limiting step. Clearance is dominated by
renal **transporters** and tubular reabsorption, which differ enormously between
species. Albumin binding modulates how much is available to be filtered, but a
2× difference cannot produce a 1000× one — and it points the wrong way.

**7.** Because the measured value is *lower* than the model assumed, which
shrinks the albumin sink and sends more dose to the cell — and the model already
over-predicted cell uptake. A parameter can be the most *sensitive* one in a
model and still not be the one that is *wrong*.

</details>

---

## Reading list

**Start here**

- **Chen et al. 2025**, *Environ Sci Technol* 59:82–91,
  [10.1021/acs.est.4c06734](https://doi.org/10.1021/acs.est.4c06734) — 60 PFAS,
  equilibrium dialysis at ν = 0.06–2, plus membrane partitioning. The best
  single dataset. Its Eq. 4 is the K_d → D bridge.
- **Qin et al. 2024**, *Environ Sci Technol* 58:9954–9966,
  [10.1021/acs.est.3c10824](https://doi.org/10.1021/acs.est.3c10824) — separates
  specific from non-specific binding, and settles the species question on whole
  plasma. Open access.

**The site problem**

- **Chen & Guo 2009**, *Arch Toxicol* 83:255–261,
  [10.1007/s00204-008-0359-x](https://doi.org/10.1007/s00204-008-0359-x) — the
  5 × 3 compound-by-site grid, and the fatty-acid controls.

**Species**

- **Starnes et al. 2024**, *Toxicol Sci* 199:132–149,
  [10.1093/toxsci/kfae028](https://doi.org/10.1093/toxsci/kfae028) — four
  albumins, one method.
- **Ryu et al. 2024**, *Environ Sci Technol* 58:14641–14650,
  [10.1021/acs.est.4c04050](https://doi.org/10.1021/acs.est.4c04050) — whole
  tissue, four organisms, two mouse strains.

**Structure and the physiological ratio**

- **Fischer et al. 2024**, *Environ Sci Technol*,
  [10.1021/acs.est.3c07415](https://doi.org/10.1021/acs.est.3c07415) — log D at
  ν ≤ 0.004, 27 compounds, three protein targets.

**Context on how much the conditions matter**

- **Gao et al. 2025**, *Environ Sci Technol* 59:26716–26728,
  [10.1021/acs.est.5c08779](https://doi.org/10.1021/acs.est.5c08779) — the ionic
  strength series.
- **Starnes et al. 2025**, *Toxicol Sci* 203:67–78,
  [10.1093/toxsci/kfae124](https://doi.org/10.1093/toxsci/kfae124) — 65 PFAS by
  DSF; read as a ranking, not as affinities.

---

## Where the data and code live

| file | what |
|---|---|
| `db/protein_binding.csv` | 865 rows, 66 study-years. **Group by study before averaging** — one study supplies 320 rows |
| `scripts/check_protein_binding.py` | validation; reproduces published statistics from three papers |
| `scripts/binding_kd_physiological_bound.py` | §4 View 2 and the §11 scoring table |
| `scripts/kd_to_partition.py` | §12, the full derivation and diagnostics |
| `scripts/binding_species_and_structure.py` | §9–10 |
| `scripts/chen_fischer_reconciliation.py` | §11a |
| `scripts/add_chen2025_tableS5.py` | the §9/§11 dataset |
| `report/REPORT.md` §6.1–6.3 | the same material as evidence in the species argument |
| `report/MBM_Qin_vs_Armitage.md` §12–14 | §13 in full |
