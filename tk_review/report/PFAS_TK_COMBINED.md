# PFAS toxicokinetics: albumin, L-FABP and renal transporters

### A combined review, written to be read start to finish

**What this is.** Three strands of the project — binding to serum albumin,
binding to liver fatty acid-binding protein, and renal transporters —
assembled into one document and connected through the quantity they jointly
determine. Part I builds the governing equation from scratch, assuming no
prior renal physiology. Parts II–V are the evidence.

**What is deliberately not here.** The QSAR fit is set aside for this round,
for reasons given in §14.3. Nothing else has been dropped: every finding from
the component reports is carried over, including the ones that cut against the
project's own conclusions.

**Figures.** All twelve live in one notebook, `notebook/pfas_tk_plots.ipynb`,
each reading one CSV and editable in place. §18 explains it. Every
chain-length axis in this document and that notebook uses **total carbons**.

---

# Part I — Foundations

## 1. What a clearance actually is

Start with the only definition you need:

> **Clearance is the volume of plasma that would have to be completely
> stripped of a substance, per unit time, to account for the rate at which
> the body is removing it.**

It is a *virtual volume per time*, not a real one. The body does not clear
some volumes completely and leave others untouched; it removes a little from
all the blood passing through. Clearance is the bookkeeping device that turns
"a little from a lot" into one number.

```
                         rate of removal  (mass / time)
          clearance  =  ────────────────────────────────
                         plasma concentration  (mass / volume)
```

The units work out as volume per time — mL/day, or, scaled to body weight,
**mL/kg/day**, which is what this project uses throughout so that a 25 g mouse
and a 70 kg human can appear in the same table.

**Why this definition and not "fraction removed per day"?** Because clearance
is additive across organs. The kidney clears, the liver clears, the gut
clears; total clearance is their sum. Fractions do not add like that. This
additivity is what lets §12 say that 58–61% of human PFOA elimination is *not*
renal — a statement that only makes sense because the routes are separable.

---

## 2. What the kidney does, in three steps

A nephron handles a substance in up to three ways. All three must be
understood before the equation makes sense.

```
   blood in
      │
      │   ① GLOMERULAR FILTRATION
      │      Blood is pushed against a sieve. Water and small solutes pass
      ▼      into the tubule. Cells and PROTEINS do not.
   ┌──────────────┐
   │  glomerulus  │──────────► filtrate enters the tubule
   └──────────────┘                      │
      │                                  │
      │  blood continues                 │   ② TUBULAR REABSORPTION
      │  past the tubule                 │      Transporters in the tubule
      │                                  │      wall pull substance BACK
      │◄─────────────────────────────────┤      into the blood.
      │                                  │
      │   ③ TUBULAR SECRETION            │
      ├─────────────────────────────────►│      Transporters push additional
      │                                  │      substance INTO the tubule,
      ▼                                  ▼      from blood that was never
   blood out                           urine    filtered.
```

**① Filtration.** The glomerulus is a sieve with a molecular-weight cut-off
far below albumin's 66,500 Da. So **a PFAS molecule bound to albumin cannot be
filtered** — the albumin will not fit through. Only the *free* molecules pass.
This single fact is why the binding half of this document exists.

**② Reabsorption.** Downstream of the sieve, transporter proteins on the
tubule's apical (urine-facing) membrane can pull material back out of the
filtrate and return it to the blood. For PFAS this is the dominant process,
and it is why half-lives are years rather than hours.

**③ Secretion.** Other transporters, on the basolateral (blood-facing)
membrane, take up substance from blood that was never filtered and hand it
across the cell into the urine. This *adds* to excretion.

Reabsorption and secretion run in opposite directions, in the same tubule, at
the same time. **Nothing in an in vivo measurement separates them** — a point
we return to in §4.3, because it constrains everything that follows.

---

## 3. Building the equation

### 3.1 Step one — how much gets filtered

Let

* `C` = total plasma concentration of the PFAS (bound + free), mass/volume
* `f_u` = the **unbound fraction**, the dimensionless share of `C` that is not
  stuck to protein. `f_u = 0.0025` means 0.25% is free and 99.75% is on
  albumin.
* `GFR` = **glomerular filtration rate**, the volume of plasma water the
  glomeruli sieve per unit time.

The free concentration is `f_u × C`. The filtrate is produced at rate `GFR`.
So the amount entering the tubule per unit time — the **filtered load** — is

```
          filtered load  =  GFR × f_u × C          (mass / time)
```

Now apply the definition of clearance. If every filtered molecule left in the
urine and nothing else happened, the removal rate would be `GFR × f_u × C`,
and dividing by `C`:

```
          CL_filtration-only  =  f_u × GFR         (volume / time)
```

**This is the single most useful quantity in the document.** It is the
clearance a compound would have if the kidney were nothing but a sieve. Call
it the filtered load clearance. It is the yardstick everything is measured
against.

### 3.2 Step two — what the tubule does to it

Let `FR` be the fraction of the filtered load that the tubule reabsorbs
(`0 ≤ FR ≤ 1`), and let `S` be the clearance contributed by secretion
(volume/time, defined the same way as any other clearance).

Material reaching the urine per unit time:

```
   urinary rate  =  (filtered load) × (1 − FR)   +   (secreted)
                 =  GFR × f_u × C × (1 − FR)     +   S × C
```

Divide by `C` to get renal clearance:

```
          CL_R  =  f_u × GFR × (1 − FR)  +  S
```

### 3.3 Step three — collecting the tubule into one number

Factor out the filtered load clearance:

```
          CL_R  =  f_u × GFR × [ (1 − FR)  +  S / (f_u × GFR) ]
                   └────┬────┘   └──────────────┬─────────────┘
                     the sieve            everything the
                                           tubule does
```

Name the bracket `T`, the **tubular factor**:

> ## CL_R = f_u × GFR × T
>
> ## T = CL_R / (f_u × GFR) = (1 − FR) + S / (f_u × GFR)

That is the equation the whole project turns on. Reading it:

| `T` | what it means | mechanism |
|---|---|---|
| **T < 1** | the kidney excretes **less** than it filtered | net **reabsorption** |
| **T = 1** | the kidney excretes exactly what it filtered | filtration alone |
| **T > 1** | the kidney excretes **more** than it filtered | net **secretion** |

`T` is dimensionless. It is the fraction of the filtered load that survives
the tubule and leaves the body.

### 3.4 Why this particular grouping, and not some other

Three reasons, and they are worth being explicit about:

1. **It isolates what we do not know.** `f_u` is measurable in a test tube.
   `GFR` is a textbook physiological constant per species. `CL_R` is
   measurable in vivo. **`T` is the only term that is not directly
   measurable** — so putting all the ignorance in one place, and solving for
   it, is the honest move.
2. **It is dimensionless, so it compares across species.** A rat clears PFOA
   at thousands of mL/kg/day and a human at hundredths. Those numbers cannot
   be compared. Their `T` values can.
3. **It separates the two reasons a species can be fast.** A species can clear
   quickly because it filters more (`f_u × GFR` large) or because its tubule
   gives back less (`T` large). **These are different mechanisms with
   different implications**, and §11.3 shows the rat and the mouse differ from
   the human by one each.

---

## 4. A worked example, with real numbers

### 4.1 Where the GFR number comes from

Human GFR is conventionally about **125 mL/min**. Converting to the project's
units for a 70 kg adult:

```
    125 mL/min × 1440 min/day  =  180,000 mL/day
    180,000 mL/day ÷ 70 kg     =   2,571 mL/kg/day
```

which is the 2570 used throughout. The number is not arbitrary and it is worth
being able to reconstruct it in your head.

### 4.2 Human PFOA, step by step

```
    f_u   = 0.00245           Ryu 2024, human plasma
    GFR   = 2570 mL/kg/day
    CL    = 0.120 mL/kg/day   US EPA 2024, TOTAL clearance
```

**Filtered load clearance:**

```
    f_u × GFR = 0.00245 × 2570 = 6.30 mL/kg/day
```

So if the human kidney were only a sieve, PFOA would be cleared at
6.3 mL/kg/day.

**The tubular factor:**

```
    T = 0.120 / 6.30 = 0.0191
```

**What that means:** the kidney excretes about **1.9% of what it filtered**.
It gives back **98.1%**.

### 4.3 Two things to notice, and they matter

**First — this is an upper bound, not a measurement.** The 0.120 is *total*
clearance, which includes biliary and faecal routes. Renal clearance can only
be smaller, so `T` can only be smaller. The honest statement is **T ≤ 0.0191**,
and therefore **at least 98.1% reabsorbed**. This is the "UPPER tier" of §10.1,
and it is logically clean: if the *whole* clearance is below the filtered load,
the tubule must be returning material, however the clearance divides between
routes.

**Second — `T` cannot be split.** Our one equation,

```
          T = (1 − FR) + S / (f_u × GFR)
```

has **two unknowns**, `FR` and `S`. One equation, two unknowns: only the
combination is identifiable. A `T` of 0.0191 is equally consistent with

* 98.1% reabsorption and no secretion, or
* 99.5% reabsorption plus a little secretion, or
* any other pair on that line.

To split them you need an independent measurement of `FR` or `S` **in the same
species**. Scaled in vitro intrinsic clearance is *not* that measurement — the
expression gate and the choice of scaling factor move it about five-fold
(§9.3), and the effect of albumin in the assay is unquantified (§9.6).

**So wherever this document says "net reabsorption", it means `T < 1`. It does
not mean secretion is absent.**

---

## 5. Run it yourself

Self-contained; no dependencies beyond the standard library.

```python
# ---------------------------------------------------------------
# The tubular factor, from first principles.
# ---------------------------------------------------------------

def filtered_load_clearance(fu, gfr):
    """The clearance a compound would have if the kidney were only a sieve."""
    return fu * gfr


def tubular_factor(cl_renal, fu, gfr):
    """T = CL_R / (fu * GFR).  <1 reabsorption, 1 filtration, >1 secretion."""
    return cl_renal / filtered_load_clearance(fu, gfr)


def percent_reabsorbed(T):
    """Only meaningful when T < 1."""
    return (1 - T) * 100


GFR = {"human": 2570.0, "rat": 14400.0, "mouse": 16700.0}   # mL/kg/day

# --- human PFOA -------------------------------------------------
fu, cl = 0.00245, 0.120
T = tubular_factor(cl, fu, GFR["human"])
print(f"filtered load clearance : {filtered_load_clearance(fu, GFR['human']):.2f} mL/kg/day")
print(f"tubular factor T        : {T:.4f}")
print(f"at least                : {percent_reabsorbed(T):.1f}% reabsorbed")

# --- now the same compound in the rat ----------------------------
# Ohmori 2003, male Wistar, renal clearance measured DIRECTLY
fu_rat, cl_rat = 0.00298, 21.2
T_rat = tubular_factor(cl_rat, fu_rat, GFR["rat"])
print(f"\nrat male PFOA  T        : {T_rat:.3f}  "
      f"({'reabsorption' if T_rat < 1 else 'secretion'})")

# --- the comparison that matters ---------------------------------
print(f"\nrat T is {T_rat / T:.0f}x the human T")
print(f"but the filtered loads differ by only "
      f"{filtered_load_clearance(fu_rat, GFR['rat']) / filtered_load_clearance(fu, GFR['human']):.1f}x")
```

Expected output:

```
filtered load clearance : 6.30 mL/kg/day
tubular factor T        : 0.0191
at least                : 98.1% reabsorbed

rat male PFOA  T        : 0.494  (reabsorption)

rat T is 26x the human T
but the filtered loads differ by only 6.8x
```

**That last pair of lines is the whole species argument in miniature.** The rat
is faster than the human at PFOA, and only a small part of it is because the
rat filters more. Most of it is that the rat's tubule gives back far less.

### 5.1 Try changing one thing

```python
# The FOUR real human measurements of PFOA binding, same clearance.
for source, fu_val in [("Smeltz 2023 (UC)",   0.001),
                       ("Fischer 2024",       0.001311),
                       ("Argoul 2025",        0.0017),
                       ("Ryu 2024",           0.00245)]:
    T = tubular_factor(0.120, fu_val, GFR["human"])
    print(f"{source:20} fu={fu_val:<9} T={T:.4f}  "
          f"{percent_reabsorbed(T):.2f}% reabsorbed")
```

```
Smeltz 2023 (UC)     fu=0.001     T=0.0467  95.33% reabsorbed
Fischer 2024         fu=0.001311  T=0.0356  96.44% reabsorbed
Argoul 2025          fu=0.0017    T=0.0275  97.25% reabsorbed
Ryu 2024             fu=0.00245   T=0.0191  98.09% reabsorbed
```

`T` moves **in exact inverse proportion** to `f_u`: halve the unbound fraction
and you double the tubular factor.

**For PFOA this does not matter and it is important to see why.** The four
laboratories span 2.5-fold in `f_u`, and `T` duly spans 2.5-fold — but all
four land far below 1, so all four say *net reabsorption*, and all four say
*more than 95% of the filtered load is returned*. The conclusion is robust to
the disagreement.

**For PFHxA the same 2.5-fold spread flips the sign** (§13.3), because there
the values straddle `T = 1` instead of sitting two orders of magnitude below
it. **Whether a binding disagreement matters depends entirely on where the
compound sits relative to the filtration line** — which is why §8's 240-fold
dispersion is reported as a property of the literature rather than as a
uniform uncertainty on every answer.

---

## 6. Check your understanding

Work these before going on. Answers follow, but the value is in attempting
them.

**Q1.** A compound has `f_u = 0.5` and a measured renal clearance exactly
equal to `0.5 × GFR`. What is `T`, and what is the tubule doing?

**Q2.** Two laboratories measure the same compound in the same species. One
reports `f_u = 0.01`, the other `f_u = 0.001`. The renal clearance is agreed.
By what factor do their `T` values differ, and which one makes the compound
look like it is being secreted?

**Q3.** A compound's *total* clearance gives `T ≤ 0.4`. A colleague says
"we can't conclude anything until we know the renal share." Are they right?

**Q4.** Rat PFOA: `f_u = 0.0030`, GFR = 14,400 mL/kg/day, renal clearance
1054 mL/kg/day (female, Kudo 2002). Compute `T`. What is the tubule doing, and
what does that imply about the sex difference given the male value of 0.49?

**Q5.** Why can't you use a scaled in vitro OAT4 clearance to work out `FR`?

<details>
<summary><b>Answers</b></summary>

**A1.** `T = (0.5 × GFR) / (0.5 × GFR) = 1`. The tubule is doing nothing on
net — excretion exactly matches filtration. Note this does *not* mean no
transport is happening; it means reabsorption and secretion are balanced.

**A2.** Tenfold, because `T` is inversely proportional to `f_u`. The
laboratory reporting the **lower** `f_u` (0.001) computes the smaller filtered
load and therefore the **larger** `T` — so that one makes the compound look
secreted. This is exactly the PFHxA situation in §13.3.

**A3.** No, they are wrong, and the reason is worth internalising. Renal
clearance is a *part* of total clearance, so `CL_R ≤ CL_total`, so the true
`T` is **at most** 0.4. Since 0.4 < 1, net reabsorption is established no
matter how the clearance divides between routes. Knowing the renal share would
sharpen the number, not change the conclusion.

**A4.** Filtered load `= 0.00298 × 14,400 = 42.9 mL/kg/day`.
`T = 1054 / 42.9 = 24.6`. That is far above 1, so the female rat is **net
secreting** PFOA — excreting about 24 times what it filters. Against the male's
0.49 (net reabsorption), **sex changes the sign, not just the magnitude**. See
§9.7 for the mechanism: a 114-fold difference in one apical transporter.

**A5.** Because `T` lumps `FR` and `S`, and an in vitro clearance tells you
about one transporter rather than about `FR`. To convert it you need a scaling
factor from the expression system to the whole kidney — OAT4 has none (§9.5) —
and the two available methods for OAT1/OAT3 disagree 11-fold and in opposite
directions (§9.3). You would also be using a protein-free assay to predict
behaviour in plasma that is 99.75% bound (§9.6).

</details>

---

## 7. A map of the rest

Now that the equation is in place, here is where each strand attaches to it.

```
        CL_R   =    f_u      ×     GFR      ×      T
                     │                             │
                     │                             │
         ┌───────────┘                             └──────────┐
         │                                                    │
   PART II: BINDING                                 PART III: TRANSPORT
   §8  serum albumin  ──── sets f_u                 §9  renal transporters
   §8b L-FABP ──── does not enter the equation          are the MECHANISM of T
       at all; it is here because it is the
       standing explanation for where PFAS
       goes INSTEAD of being excreted
                             │                              │
                             └──────────┬───────────────────┘
                                        │
                              PART IV: INTEGRATION
                              §10–13  what T actually is,
                                      per compound and species
```

---

# Part II — Binding

## 8. Serum albumin

`report/ALBUMIN_BINDING.md` · `db/protein_binding.csv` ·
`db/albumin_kd_to_k.csv` · notebook figures 1–4

Albumin is the most abundant protein in plasma, present at roughly **600 µM**,
and it is the principal carrier of anionic PFAS in the circulation. From §3.1,
the fraction *not* bound to it is what gets filtered — so albumin binding is
not an incidental physicochemical property but a direct determinant of the
long half-lives that make PFAS a regulatory problem.

**288 albumin–water partition coefficients for 63 PFAS from 20 primary
sources**, converted to one basis and audited against the primary publication
in each case.

### 8.1 The common basis, and why the arithmetic is not the problem

Published values arrive as dissociation constants, association constants,
volume-basis partition coefficients and mass-basis partition coefficients.
Converting between them requires a stoichiometry that is frequently not
reported. Everything here is put on **log *K*, litres of water per kilogram
of protein**, through an explicit identity:

```
        log K  =  4.1774  +  log₁₀(n)  −  log₁₀(K_d / µM)
```

All **181 rows carrying both a K_d and a stoichiometry satisfy it exactly**.
State that plainly in any discussion: the conversion arithmetic is auditable
and correct. The uncertainty is in the measurements, not the algebra.

| | |
|---|---|
| range | log *K* 0.48 to 5.49 |
| median | 3.37 (IQR 2.46–4.13) |
| species split | 178 human, 110 bovine |

### 8.2 The dispersion is systematic, not random

Across the 32 compounds with three or more values the **median within-compound
range is 2.38 log units**; the maximum is 4.22.

> **A modeller selecting one published constant is making a choice with a
> 240-fold range of defensible outcomes. A modeller averaging them is
> averaging quantities that are not the same quantity.**

Three features of the primary literature produce this, and they are worth
separating because only one of them is fixable by better measurement:

1. **Incompatible units**, requiring a stoichiometry that is often unreported.
   *Fixable by documentation* — and done, above.
2. **Measurement conditions that move the answer** — ligand-to-protein ratio,
   temperature, ionic strength, pH, co-solvent — varying widely and not always
   stated. *Partly fixable.*
3. **Several widely cited values are not the quantity their citing papers take
   them for.** *Only fixable by reading the primary source.* §14 lists eight
   instances found by doing exactly that.

**Method is the principal axis.** Equilibrium dialysis and differential
scanning fluorimetry differ by **2.35 log units in median**, in the direction
the physics predicts, because DSF measures thermal denaturation rather than an
equilibrium. **These are not competing estimates of one quantity; they are two
quantities.**

*(Notebook figure 3 shows the family-by-family spread among the K_d-based
assays and prints its own numbers. It cannot show the DSF comparison, because
DSF reports a melting temperature and so contributes no K_d — the cell says so
rather than implying otherwise.)*

### 8.3 Chain length

Affinity rises with chain length to a maximum at **PFNA (9 total carbons)**,
then declines. The rise is **+0.52 log units per carbon from PFBA to PFNA**.

| total carbons | compound | median log *K* | n |
|---|---|---|---|
| 4 | PFBA | 1.68 | 10 |
| 5 | PFPeA | 2.23 | 6 |
| 6 | PFHxA | 2.77 | 15 |
| 7 | PFHpA | 3.54 | 12 |
| 8 | PFOA | 3.75 | 25 |
| **9** | **PFNA** | **4.29** | **18** |
| 10 | PFDA | 3.77 | 14 |
| 11 | PFUnDA | 3.56 | 11 |
| 12 | PFDoDA | 2.98 | 7 |
| 13 | PFTriDA | 2.68 | 5 |
| 14 | PFTeDA | 3.46 | 4 |

The decline past PFNA is attributed to steric constraint in albumin's pockets
and to falling aqueous solubility, which reduces the free concentration
available to the binding equilibrium. Fischer's statement is the crispest:
*"for ηpfc > 8 PFAS, reduced molecular flexibility limits their accessibility
to the HSA binding pockets."* The C14 value interrupts the decline but rests
on four measurements and should not be weighted against the trend.

> **A documentation defect, recorded rather than propagated.** Table 5 of
> `ALBUMIN_BINDING.md` labels this column *n*_CF — perfluorinated carbons —
> but the values in it are **total** carbons (it lists PFBA as 4; PFBA has 3
> perfluorinated carbons). Table 10, for L-FABP, genuinely uses perfluorinated
> carbons. The two tables therefore *appear* to put the maximum at different
> chain lengths when both put it at PFNA. The table above is relabelled
> correctly.
>
> **The conversion is not a constant**, which is what makes this easy to get
> wrong: a carboxylate `CF₃(CF₂)ₙCOOH` has one carbon that is not
> perfluorinated, so total = perfluorinated + 1; a sulfonate `CF₃(CF₂)ₙSO₃H`
> has none, so total = perfluorinated. Add one to everything and every
> sulfonate shifts a place. The notebook encodes this in exactly one function,
> `total_carbons()`.

### 8.4 Species: the comparison cannot be made

Human and bovine albumin differ in median by **0.90 log units**, bovine the
higher (3.99 against 3.09).

**This should not be read as a species difference.** It is confounded with
method and compound coverage: the bovine subset is dominated by
microextraction and dialysis from sources that measured long-chain compounds,
while the human subset contains all 33 differential scanning fluorimetry
values. The two are kept separable in the delivered data for that reason, and
pooling them would be a mistake.

> **The general form of this problem, which recurs in §8b.4:** the
> disagreement *within* a species is larger than the difference *between*
> species. **A species extrapolation factor for binding cannot be estimated
> from this literature**, because the quantity it would correct for is smaller
> than the disagreement about it.

---

## 8b. Liver fatty acid-binding protein

`report/ALBUMIN_BINDING.md` §5 · notebook figures 5–6

**L-FABP does not appear in the equation of §3.** It is intracellular and
hepatic; it does not set `f_u` and it does not act on the tubule. It is in
this document because it is the standing explanation for where PFAS goes
*instead* of being excreted — several PFAS accumulate preferentially in liver,
and L-FABP binding is the usual account of why.

**91 coefficients for 32 PFAS from ten sources**, on the same basis, using a
per-species bridge constant of **4.8475** for human L-FABP against albumin's
4.1774.

### 8b.1 The comparison that decides the hypothesis

Twenty-one compounds have been measured against both proteins. The median
difference is **−0.12 log units**, L-FABP marginally the **weaker** — far
inside the dispersion of either compilation.

**This is the most consequential number in the section and it is easy to
misread, so state both halves:**

- **It does not say the two proteins are interchangeable in the body.**
  Albumin is at ~600 µM in plasma, L-FABP at 100–400 µM in hepatocyte
  cytosol, in different compartments. Equal affinity per kilogram of protein
  still produces very different capacity in context.
- **It does say there is no intrinsic affinity advantage.** Any explanation of
  preferential hepatic accumulation that rests on L-FABP binding PFAS *more
  tightly* than the plasma carrier does not survive this comparison.

![albumin against L-FABP](fig_binding_chain_profile.png)

*Rendered by `scripts/fig_binding_chain_profile.py` from
`db/binding_chain_profile.csv`; the same comparison is notebook figure 5,
which rebuilds it from the per-measurement files instead of the summary.*

### 8b.2 The in vivo test points the same way

The one available in vivo test is a knockout: deleting L-FABP changes the
mouse liver partition coefficient by **under 20% for 15 of 19 PFAS**.

An explanation carrying the accumulation should not survive its own deletion
that quietly. Two independent lines — a binding comparison and a gene
deletion — agree, which is more than most claims in this literature have.

### 8b.3 A structural cause for part of the scatter

L-FABP has **two binding sites**, and most of its constants come from probe
displacement assays. **Sources using a secondary-site probe and sources using
a primary-site probe are not estimating the same quantity.** Luebker's assay,
for example, displaces DAUDA, which binds the secondary site — and the paper
says so.

So part of this dispersion is *structural* and fixable by stratifying on probe
site, rather than being irreducible experimental noise. That is a concrete,
cheap improvement available to anyone revisiting this dataset.

Within-compound dispersion: median range **1.40 log units**, maximum **3.30**
for PFOA across nine sources (log *K* 2.55 to 5.85). Better behaved than
albumin, but on half as many compounds.

### 8b.4 Chain length, same shape, same maximum

Affinity rises **0.36 log units per carbon from PFBA to PFNA** and falls
thereafter — **the same shape and the same maximum as albumin**, which rose at
0.52.

Zhang attributes the L-FABP decline to the binding cavity: the fluorocarbon
chain is rigid, and beyond about 13 Å it must bend to be accommodated, which
costs more energy than a flexible hydrocarbon chain would. **That explanation
is specific to a small enclosed cavity and does not transfer to albumin** —
yet the same inflection occurs there. Either two different mechanisms produce
the same inflection at the same chain length, or the real cause is a property
of the ligand rather than of either protein. Worth raising; not settled.

### 8b.5 Species, and two published claims that do not survive

Restricting both sides to displacement assays, so method is held fixed:

| comparison | shared compounds | median rat − human |
|---|---|---|
| rat vs Zhang 2013 | 10 | **+0.50** |
| rat vs Sheng 2016 | 2 | **−0.71** |
| rat vs Sheng 2018 | 3 | −0.62 |
| rat vs Yang 2020 | 6 | −0.10 |
| rat vs Yang 2026 | 1 | −0.06 |

The five human sources disagree among themselves by **1.21 log units** — PFOA
alone is 2.36 µM in Sheng 2016 and 50.4 µM in Zhang 2013, **a factor of 21
between two human laboratories.** That spread is wider than any rat–human gap
and contains both signs.

> **The sign of the apparent species difference is set by the choice of human
> comparator.** Pooled across all sources: +0.38 log over 13 shared compounds,
> falling to **−0.01** on the five compounds carrying four or more independent
> human values — the signature of thin coverage rather than biology.

Two published claims sit in opposition here and **neither survives cleanly**:

- **Modaresi 2025** states PFAS bind human FABP more strongly than rodent
  FABP. That holds against exactly one human source, on two compounds, and
  reverses against the other three and against the pooled comparison.
- **Fan 2025** states their rat affinities are comparable to human values by
  the same method, citing Zhang 2013. The agreement is near-exact against
  **Yang 2020, which they do not cite for it**; against Zhang 2013, which they
  do cite, the difference is the largest of the four.

As with albumin, **a rat-to-human extrapolation factor for L-FABP binding
cannot be estimated from this literature.**

---

# Part III — Transport

## 9. Renal transporters

`report/TRANSPORTERS.md` · `db/transporter_direction.csv` ·
`db/transporter_expression.csv` · `db/transporter_clearance_stage3.csv` ·
notebook figures 7–8

Transporters are the **mechanism of `T`**. Reabsorptive transporters on the
apical membrane push `T` below 1; secretory transporters on the basolateral
membrane push it above.

### 9.1 Localisation is agreed; substrate identity is not

**299 of the 305 rows in `db/transporter_direction.csv` resolve to a
reabsorption-or-secretion call, across 21 distinct transporters, and the
number of genuine direction conflicts is zero.** Five rows carry no parseable
direction and one is a combined two-transporter row. The apparent conflicts
are free-text variants of the same call — "secretion (basolateral)" against
"secretion (basolateral uptake into the proximal tubule cell)" — plus
"not operative here" qualifiers.

That is a result, not a formality: **papers that disagree sharply about which
compounds a transporter carries do not disagree at all about where it sits.**
When reading this literature, treat localisation as settled and substrate
lists as contested.

> *Note for the record:* `report/TRANSPORTERS.md` states this as "343 of 363
> rows". No current data file has 363 rows and no script reproduces that
> figure, so the counts above were recomputed from the delivered file. The
> substantive claim is confirmed; the older row counts are stale.

### 9.2 Expression, measured in one assay across three species

Thakur 2024 quantified 78 enzymes and transporters in kidney membrane
fractions of human (n=15), rat (n=10) and mouse (n=14), both sexes. **This is
the single most useful table in the strand, because it is one assay** — the
species comparison is not assembled across laboratories, so it is not exposed
to the between-laboratory problem that §8.4 and §8b.5 describe.

pmol/mg membrane protein:

| transporter | membrane | net process | human | rat | mouse |
|---|---|---|---|---|---|
| **Oat1** | basolateral | secretory | 4.8 ± 1.2 | **32.8 ± 13.7** | 19.1 ± 7.1 |
| **Oat3** | basolateral | secretory | 4.3 ± 1.1 | 6.3 ± 1.1 | 6.5 ± 3.6 |
| **OAT4** | apical | **reabsorptive** | **0.5 ± 0.2** | **not detected** | **not detected** |
| **Oatp1a1** | apical | **reabsorptive** | **not detected** | **1.9 ± 2.1** | 3.1 ± 2.8 |
| URAT1 | apical | — | 3.1 ± 1.2 | 2.2 ± 0.7 | 10.8 ± 1.1 |

*(URAT1 transports no PFAS in this compilation.)* Independent agreement where
sources overlap: human OAT1 4.8 (Thakur) against 5.3 (Niu 2026); P-gp 1.2
against 2.1; BCRP not detected in both.

**The two consequences:**

- **The secretory side is conserved in kind and differs in degree.** Both
  species have OAT1 and OAT3; the rat simply has 6.8× more OAT1.
- **The reabsorptive side is not the same machinery at all.** Human has OAT4
  and no Oatp1a1. Rat has Oatp1a1 and no OAT4.

> **This is the most important fact in the transporter strand**, because
> reabsorption is what sets the long half-lives, and it is the half of the
> system that is **not shared**. A rat-to-human extrapolation of reabsorption
> is an extrapolation between two different proteins that happen to do
> similar jobs — defensible as an empirical scaling, indefensible as a
> mechanistic one.

Niu 2023 calls OAT4 *"a human-specific transporter … it may be one of the
transporters responsible for species-specific differences in PFAS
toxicokinetics"*; EFSA says *"OAT4 (not expressed in rats) and urate
transporter 1 (URAT1), but not Oatp1a2 (the closest orthologous of rat
Oatp1a1), were shown to mediate PFOA"*. Here it is **measured** rather than
cited.

**Two long-standing disagreements dissolve.** ASBT and OSTα/β were not
detected in human *or* rat kidney, yet Zhao 2015 found transport by both in
transfected cells while Argoul 2026 found none. **Both are right.** A
transfected-cell result states what a protein *can* carry, not whether it is
present to do it.

> **Capability versus presence resolves more apparent conflicts in this
> literature than any other single idea.** Carry it into any paper you read.

### 9.3 Abundance gates relevance

52 rows carry an intrinsic clearance. Eleven can be scaled to the kidney —
human OAT1 and OAT3, by relative expression (REF) or relative activity factor
(RAF). **The two methods disagree 11-fold on OAT3 and in opposite
directions**, so they are reported separately and never averaged.

**The expression gate removes the compilation's three highest values**, all
BCRP from Niu 2026, including **279 µL/min/mg for PFOA — the highest intrinsic
clearance anywhere in the dataset.** Human renal BCRP is not detected by
Thakur or by Niu.

> **The fastest transporter in the dataset is one that two laboratories cannot
> find in human kidney.** This is the cleanest illustration of why an in vitro
> clearance ranking is not a mechanistic ranking.

### 9.4 Reabsorption against secretion, within one laboratory

Louisse (2023 for OAT4, 2024 for OAT1/OAT3), measured and unscaled, µL/min/mg:

| compound | total C | OAT4 (reabsorptive) | OAT3 | OAT1 | OAT4 / OAT3 |
|---|---|---|---|---|---|
| PFHpA | 7 | 75.0 | 8.7 | 21.4 | 8.6 |
| PFHxS | 6 | 79.3 | 7.5 | — | 10.6 |
| PFOA | 8 | 95.7 | 16.7 | 18.9 | 5.7 |
| PFNA | 9 | 146.6 | 23.8 | — | 6.2 |
| PFOS | 8 | 45.8 | — | 53.5 | **—** |
| PFDA | 10 | 153.8 | — | — | |

Both sides rise with chain length — OAT3 2.7-fold from PFHpA to PFNA, OAT4
2.1-fold from PFHpA to PFDA — **so chain length does not by itself tip the
balance.**

What stands out is the level, and **it splits by head group**: for every
carboxylate, OAT4 carries **3.1–8.6 times** the clearance of a secretory
transporter. PFHxS behaves like the carboxylates. **PFOS is the exception** —
OAT1 matches OAT4, and OAT3 does not transport it at all.

### 9.5 The caveat that bounds §9.4, and it is large

**OAT4 cannot be scaled to the kidney and OAT1/OAT3 can.** Applying the
expression factor lifts OAT1 to **118–334 µL/min/mg**, which would **reverse
the ordering** in the table above.

> **So the in vitro ratio is not evidence about the in vivo balance. It is
> evidence about the proteins.** Any statement of the form "reabsorption
> dominates secretion in the human tubule" that cites §9.4 is over-reading it.

This is also the concrete reason `T` cannot be split (§4.3): the one in vitro
quantity that might have given us `FR` independently is the one that has no
scaling factor.

### 9.6 Does any of it survive contact with plasma protein?

**Every in vitro assay in this compilation was run without albumin.**

The concern is quantitative, not theoretical. Transporters see the *free*
concentration. An assay run protein-free presents the transporter with a free
fraction of 1, against an in vivo value near 0.0025 for PFOA — **a 400-fold
difference in the quantity the protein actually responds to.**

This does not necessarily invalidate the *rankings*; if binding affects all
compounds similarly, much of it divides out. It does call the *magnitudes*
into question, and **any IVIVE that scales these clearances to a kidney
without a binding correction is overstating them.**

### 9.7 The rat sex difference is a change of sign

Female rats **net-secrete** PFOA and PFNA; males do not. Established twice
over, by routes sharing no assumptions:

- **Oatp1a1 is 114-fold more abundant in male rats** than in female — the
  reabsorptive transporter, present in the sex that retains.
- **In vivo renal fractions** put the female above the secretion threshold and
  the male below it.

Worked through with the numbers from §6 Q4: female `T = 24.6`, male
`T = 0.42–1.07`. **The same compound, the same species, the same tubule — and
the sign of the tubular term inverts.**

> This is the clearest mechanism-to-outcome link anywhere in the project, and
> the template for what every other species lacks: an abundance measurement
> and an in vivo measurement that agree, in the same animal, about the same
> transporter.

---

# Part IV — Integration

## 10. How a tubular factor gets established

`db/tubular_factor_consolidated.csv` (52 rows) ·
`db/tubular_factor_evidence.csv` (87 rows) · notebook figure 9

### 10.1 Three tiers of evidence

| tier | what is available | what it yields |
|---|---|---|
| **DIRECT** | renal clearance measured | `T` itself |
| **UPPER** | only total clearance | `T ≤ CL_tot/(f_u·GFR)` |
| **LOWER** | measured renal fraction | secretion established once `f_ren > (f_u·GFR)/CL_tot` |

**The UPPER tier is the workhorse, and it is logically clean.** Renal
clearance is a part of total clearance, so a bound computed from the total is
a true upper bound on `T`. **If the whole clearance is less than the filtered
load, the tubule must be returning material, no matter how the clearance
divides between routes.** (This is §6 Q3.)

**The rule the project enforces:** a value may be combined across studies only
**within one species**. `pick_fu()` raises `SystemExit` on a cross-species
join rather than warning, and every delivered row records which binding matrix
it used in a `fu_matrix` column. Cross-*strain* joins within a species are
permitted and recorded, not hidden.

![T against chain length](fig_T_vs_chain_length.png)

*Rendered by `scripts/stage6_figure_T_vs_chain.py` from
`db/tubular_factor_consolidated.csv` and
`db/human_tubular_sign_by_fu_source.csv`; notebook figure 9 is the same plot.
**Two different spreads are drawn and must not be conflated:** a coloured bar
spans that species' own* clearance *sources under a single unbound fraction;
a grey whisker spans the four human* f_u *datasets with the clearance held
fixed. The grey whiskers are the wider ones — which is §8.2 and §13 in one
picture.*

### 10.2 What is established

**Long-chain PFAS are almost completely reabsorbed in the human.** PFOA
through PFUnDA exceed **98% reabsorbed** on every unbound-fraction source
available; PFHxS and PFOS reach **97.0%**. PFBS is the exception at **87%**,
and short chains generally are unsettled (§13.3).

**The rat/human tubular factor is of order 10²** — **22 to 161-fold** across
the four human and three rat clearance sources. Every combination leaves the
tubule at least three-fold above the **6.8-fold glomerular factor**, so the
dominance is secure and the number is not.

> **Quote it as an order of magnitude.** A point estimate implies a precision
> the evidence does not have.

### 10.3 The two rodents are fast for different reasons

Against the human, **the rat's difference is carried by the tubule and the
mouse's by the glomerulus.**

"Rodents clear PFOA faster than humans" is one sentence covering two
mechanisms — and the distinction is load-bearing, because **only one of them
is a transporter story.** A mouse-based argument about reabsorptive transport
is not supported by the mouse's own kinetics.

### 10.4 What is not established

- **The sign of human tubular handling for C5–C7.** §13.3.
- **Whether human OAT4 reabsorbs more effectively in situ than rat Oatp1a1.**
  Per mg of transporter protein they are within 1.3-fold; **per kidney they
  cannot be compared, because neither has a scaling factor** (§9.5).
- **Whether any measured in vitro clearance survives contact with plasma
  protein** (§9.6). Given §8.2's 240-fold dispersion in binding, this is not a
  minor omission.
- **Anything about the monkey.** Fourteen clearances exist across six
  compounds; **no unbound fraction has ever been measured in cynomolgus plasma
  for any PFAS**, and the 1,497-row binding collection contains **no macaque
  measurement of any kind** (§16).

---

## 11. Where the renal framework stops

**58–61% of human PFOA elimination is not renal.**

Say this early in any renal discussion, so the tubular work is not asked to
carry more than it can. The renal axis explains the **species differences** in
elimination well, and the **absolute** human half-life only partially. Biliary
and faecal routes, and menstrual loss in women of reproductive age, carry the
rest.

This is a scope statement, not a weakness — but **a mass-balance model that
closes only on renal clearance will be wrong by more than a factor of two in
the human**, and that should be stated before someone discovers it downstream.

---

## 12. Is the gap clearance, or binding?

`db/clearance_fu_coverage.csv` · `scripts/stage11b_clearance_fu_coverage.py` ·
notebook figure 11

The project holds **1,497 protein-binding rows** (767 human, 190 bovine, 162
mouse, 153 rat, 31 rainbow trout) across albumin, FABP, γ-globulin and tissue
homogenates. The obvious hope is that they extend the renal endpoint.

**They do not unlock a single additional species-compound**, and the reason is
structural rather than a matter of digging harder. **Binding supplies `f_u`,
and `f_u` is not the scarce term.**

| species | compounds with a clearance | already have f_u | gap |
|---|---|---|---|
| human | 16 | 15 | 1 (DONA) |
| rat | 12 | 9 | 3 (PFOS isomers, clearances read "NR") |
| mouse | 12 | 10 | 2 (PFTriDA, PFTeDA) |
| **monkey** | **6** | **0** | **6** |

The only places binding could help are mouse PFTriDA/PFTeDA and the monkey,
and the collection has neither: **no mouse plasma, serum or albumin binding
for C13–C14** — PFTriDA appears only as human and bovine albumin, and
PFTeDA's mouse rows are *tissue homogenates*, which do not give a plasma
unbound fraction — and **no macaque binding for any PFAS at all.**

Over all 46 species-compounds that have a clearance, the number the binding
collection could unlock is **zero**.

> **The bottleneck is renal clearance, and no amount of binding work moves
> it.** A negative result, but a load-bearing one: it redirects effort away
> from the largest dataset in the project and toward the smallest.

---

## 13. The human unbound fraction

`db/human_fu_recommended.csv` (28 compounds) · notebook figure 10

Four **independent** human datasets exist. A fifth that looks independent is
not: **Abraham 2024's binding is Smeltz's**, by its own Table 4 footnote
("as published by Smeltz et al. (2023)"). Treating it as a fifth source would
double-count the EPA measurement.

| dataset | method | compounds |
|---|---|---|
| **Smeltz 2023** | ultracentrifugation | 67 (US EPA) |
| **Fischer 2024** | solid-phase microextraction | 29 |
| **Argoul 2025** | Dianorm dialysis | 16 |
| **Ryu 2024** | presaturation dialysis | 14 |

### 13.1 The recommendation is by purpose, not a single winner

| use | source | why |
|---|---|---|
| **human absolute values (the default)** | **Smeltz 2023 (UC)** | 67 PFAS so the series has no gaps; US EPA's own; the method choice is *argued* — UC and RED are run side by side on 22 compounds and RED is shown failing to equilibrate for high-mass, high-logK_ow PFAS; per-compound Bayesian credible intervals; and Abraham's in vivo mass balance already runs on it |
| **individual variability** | **Fischer 2024** | measures albumin and globulin binding *separately* and returns f_u from a person's own serum protein concentrations, R² = 0.99 |
| **species ratios** | **Ryu 2024, and only Ryu** | the only dataset measuring human, rat and two mouse strains in one laboratory; in a *ratio*, a systematic method bias divides out |

**Do not use Ryu for human absolute values.** It runs at the high end — 4×
above Smeltz for PFHxA, 15× above Fischer for PFPeA — and those are precisely
the two compounds in dispute.

Two caveats travel with Smeltz: **PFBS is 0.0046 in its Table 2 and 0.0128 in
its Table 3**, a 2.8× internal discrepancy within the same assay; and its PFOS
value is an average over three salt forms.

### 13.2 "Variable f_u" means three different things

1. **Variable across compounds.** Settled, and the single largest correction
   to the regulatory position: OEHHA's flat 0.02 against a measured range of
   **0.0004 (PFHpA) to 0.31 (PFPA) — a factor of 780** across the series.
   Already implemented throughout this review.
2. **Variable across individuals.** The one worth adding. Fischer's equation
   returns f_u from measured albumin and globulin; f_u varies up to **2.5×
   across 1,919 NHANES participants**, widest below 500 g/mol — the
   short-chain range. The band is carried per compound in
   `human_fu_recommended.csv` (2.5× below 500 g/mol, 1.6× for 500–600, 1.3×
   above).
   > **Check the direction before using it.** A higher f_u means a larger
   > filtered load and therefore **faster** elimination, so the individuals
   > with the **lowest** f_u carry the **longest** half-lives. For a
   > health-protective value the low tail is the sensitive one, not the high.
3. **Variable with concentration.** Real — Qin 2024 shows binding goes
   non-linear above a PFAS:protein molar ratio near 0.1 — but environmental
   exposures sit near **0.0005**. This matters for *comparing in vitro
   datasets*, not for human dosimetry. Lowest priority of the three.

### 13.3 What the choice actually decides: two compounds

Human PFHxA, the same clearance, four binding datasets:

| f_u source | T | verdict |
|---|---|---|
| Smeltz 2023 | 3.47 | net secretion |
| Fischer 2024 | 2.48 | net secretion |
| Argoul 2025 | 1.18 | net secretion |
| **Ryu 2024** | **0.62** | **net reabsorption** |

**Ten of the twelve PFAS with enough data give the same tubular sign on every
binding dataset, and all ten are reabsorption.** The disagreement is confined
to **PFPeA and PFHxA**.

So the decision is far narrower than the 780-fold spread in f_u makes it look
— compare the PFOA case in §5.1, where the same kind of disagreement changes
nothing because the compound sits two orders of magnitude below the
filtration line. **Whether a binding disagreement matters depends on where the
compound sits relative to `T = 1`.**

It is not cosmetic, though: those two compounds also decide whether the human
chain-length relationship is monotone or humped, which is a question about
model *structure*, not parameters.

---

# Part V — Quality, gaps and discussion

## 14. Errors found in the source and regulatory literature

Eight, documented with evidence in `report/TRANSPORTERS.md` §9.

| # | error | size |
|---|---|---|
| 1 | **ATSDR 2021 Table 3-6 misreports Ohmori 2003.** Rat PFDA total clearance printed as 207/140 against the paper's 5.2/5.3; male PFOA 135 against 50.5 | **40×**, 26×, 2.7× |
| 2 | **OEHHA PHG Table A6.4 assumes f_u = 0.02 for every species, uncited.** Measured values are 0.0025 (human) and 0.0030 (rat), so the filtered loads — and every reabsorption percentage derived from them — are wrong in the denominator | ~**8×** |
| 3 | Nakagawa 2008 Table 1 carries a units defect in a K_m/V_max column | — |
| 4 | Niu 2026 transposes P-gp and BCRP between its text and abstract | — |
| 5 | Argoul 2026's citation for its unbound fractions does not resolve | — |
| 6 | Abraham 2024 misquotes Smeltz in its own Discussion (§4.5) | — |
| 7 | Smeltz 2023 reports PFBS twice, 2.8-fold apart (Table 2 vs Table 3) | 2.8× |
| 8 | **`pdftotext` silently drops the micro sign** — caught twice here (fmol/µg→fmol/mg; 48.64 µmol/kg→mmol) | **1000×** |

### 14.1 Error 8 is the one most likely to bite someone else

It is not a literature error but a methods warning. The micro sign vanishes on
extraction, and the loss is **undetectable by eye and detectable by
arithmetic**: a value 1000-fold out of line with its neighbours is a dropped
prefix until proven otherwise. Check magnitudes against the series, not
character counts against the PDF.

### 14.2 A number that should stop being quoted

**The male rat's PFOA reabsorption of 93.7%.** Three rat studies give
`T = 0.42 to 1.07` on OEHHA's GFR and 0.77 to 1.95 on the rate Ohmori cites.
**At most half the filtered load is returned, possibly none.**

### 14.3 Why the QSAR is set aside this round

Reported here rather than buried, because the diagnosis is itself a finding.

A per-species censored regression of `log₁₀ T` on chain length and head group
was fitted (46 usable rows, 31 of them one-sided bounds, Tobit likelihood).
Results:

| species | n | censored | nC effect | σ |
|---|---|---|---|---|
| human | 19 | 15 | ×0.39 per carbon (z = −1.65) | 0.80 |
| rat | 11 | 7 | ×0.72 per carbon (z = −0.56) | 1.17 |
| mouse | 16 | 9 | ×0.57 per carbon (z = −1.93) | 0.67 |

**The sign agrees across all three species and no linear coefficient reaches
significance** (largest |z| = 1.93). Adding six mouse rows moved its z from
−1.60 to −1.93 — **a 60% increase in n bought about 0.3 of a z.**

A quadratic improves the human (ΔAIC −7.08) and rat (−12.04) but not the mouse
(+0.86), putting a **maximum** near 5–7 carbons — mechanically expected, since
`T` divides by `f_u` and `f_u` has a minimum there, so `T` inherits an
inverted U. **But in the human the curvature does not survive removing PFPeA
and PFHxA** (ΔAIC +7.78 without them) — the two compounds of §13.3. In the
rat it does survive (−12.04 → −6.10 without the digitised rows), and the rat
has no four-way binding dispute.

> **So the choice of human f_u source decides whether the human chain-length
> relationship is monotone or humped.** The binding question is upstream of
> the model structure, not a detail beside it.

**It is not wrong, it is underpowered — and what it is short of is measured
renal clearances, not cleverer descriptors.** That is why §15 matters more
than another modelling round. The fit and all its sensitivity analyses remain
in `db/qsar/` and `report/TRANSPORTERS.md` §11–11a; nothing has been deleted.

---

## 15. Is there more direct renal clearance in the literature?

`db/renal_clearance_leads.csv` · `scripts/stage12_renal_clearance_leads.py` ·
notebook figure 12

**Short answer: yes, about ten compounds' worth, and it has never been
ingested.** The human DIRECT tier currently holds **two** compounds (PFBS and
PFHxA, both from Argoul). Everything else human is an upper bound from total
clearance.

### 15.1 What was found

EFSA's 2020 Opinion tabulates human renal clearances for ten compounds, from
four primary studies:

- **Zhang 2013** — 86 paired blood/urine samples, adults in China,
  isomer-specific. The widest-coverage human source.
  [10.1021/es401905e](https://doi.org/10.1021/es401905e)
- **Fujii 2015** — renal, biliary and faecal clearances in humans
- **Harada 2005** — PFOA and PFOS, young/old × male/female, n = 5 per cell.
  [10.1016/j.envres.2004.12.003](https://doi.org/10.1016/j.envres.2004.12.003)
- **Harada 2007** — follow-up

*(Bibliographic details for Zhang 2013 and Harada 2005 retrieved from
PubMed.)*

### 15.2 Why they are not yet in the pipeline

Two reasons, each disqualifying on its own:

1. **The secondary source contradicts itself.** Two passages of the *same*
   EFSA document give different values for the same compounds — PFHpA 0.61
   against 0.67, and **PFUnDA 0.064 against 0.005, a 13-fold gap.** That is
   not a rounding difference.
2. **Some values exceed independently reported total clearances.** Renal
   clearance cannot exceed total clearance in the same subject. EFSA's PFOA
   renal range (0.13–0.77 mL/kg·day) sits **3- to 18-fold above** Abraham
   2024's measured *total* clearance for PFOA (0.0417) and Argoul's (0.042).
   Different cohorts, so not a strict contradiction — but a gap that large
   between one group's renal part and another group's whole must be
   understood before either is used.

### 15.3 Would they change anything?

Provisional `T`, computed on the project's own human f_u and GFR, **clearly
marked as not-yet-admitted**:

| compound | total C | CL_R (mL/kg·d) | provisional T | currently held as |
|---|---|---|---|---|
| PFHxS | 6 | 0.027–0.039 | 0.015–0.022 | bound |
| PFHpA | 7 | 0.61 | 0.053 | bound |
| PFOA | 8 | 0.13–0.77 | 0.021–0.122 | bound |
| PFOS | 8 | 0.037–0.106 | 0.019–0.055 | bound |
| PFNA | 9 | 0.062–0.25 | 0.013–0.052 | bound |
| PFDA | 10 | 0.066–0.096 | 0.042–0.060 | bound |
| PFUnDA | 11 | 0.064–0.065 | 0.059–0.060 | bound |
| PFDoDA | 12 | 0.070 | 0.259 | bound (currently uninformative) |
| PFTriDA | 13 | 0.077 | no human f_u | absent |
| PFTeDA | 14 | 0.224 | no human f_u | absent |

**Eight would upgrade from a bound to a measurement** — which is the scarce
thing. **All eight point the same way the project already concludes.** So the
*direction* is not at stake; what these buy is evidence tier and magnitude,
which is exactly what §14.3 is short of.

**And the honest limitation:** none of these covers **PFPeA or PFHxA**, the
two compounds whose sign is actually disputed. **This literature strengthens
the part that was never in doubt and leaves the open question open.**

### 15.4 Recommended next step

Three PDFs would land most of it: **Zhang 2013, Fujii 2015, Harada 2005.**
The pattern that has worked in this project is for primary PDFs to be read
directly rather than taken from a secondary tabulation — it is how the Ohmori,
Smeltz and Abraham numbers were landed, and **in each case reading the primary
changed something.**

---

## 16. What is open, and what would close it

Stated in advance, so that the answers count whichever way they come out.

| gap | what would close it | cost | value |
|---|---|---|---|
| **The monkey has no binding data at all** — 14 clearances across 6 compounds, 0 unbound fractions, and no macaque measurement anywhere in the 1,497-row collection | one dialysis or ultracentrifugation run on cynomolgus plasma across six compounds | **one afternoon of bench work** | **highest.** Converts 14 orphan clearances into 14 tubular factors in the species physiologically closest to the human |
| **Eight human compounds sit at bounds rather than measurements** | read Zhang 2013, Fujii 2015, Harada 2005 primary (§15) | three PDFs | high |
| **The sign for PFPeA and PFHxA is undecided** | an independent human f_u measurement for those two compounds | one assay, two compounds | high — it is the only thing in dispute, and it decides the model's functional form (§14.3) |
| **No in vitro transport assay has been run at physiological albumin** | one OAT4 or OAT1 assay with albumin present | one experiment | high — it tells us whether the whole in vitro compilation transfers (§9.6) |
| **OAT4 has no scaling factor** | a relative-expression or relative-activity measurement for the OAT4 expression system | one calibration | medium — it is the single missing piece preventing the in vitro ratio in §9.4 from being taken in vivo |
| **The two scaling methods disagree 11-fold** | an agreed convention, or an experiment that discriminates them | discussion, then possibly one experiment | medium — may be worth more than new measurements |

**What is explicitly *not* on this list: more binding measurements on
compounds that already have them.** The project has 1,497 and they unlock
nothing further (§12).

---

## 17. Discussion questions, with suggested directions

Each carries a **suggested direction** — a position to argue from, not a
settled answer. Where the evidence does not support a position, that is said.

### 17.1 On the species difference

**Q1. If the two rodents are fast for different reasons — the rat at the
tubule, the mouse at the glomerulus — which is the right model for human
reabsorption, and does "rodent" mean anything as a category here?**

> **Suggested direction.** Neither is a model of human *mechanism*, and the
> category "rodent" should be retired from this argument. The human reabsorbs
> with OAT4, which neither rodent has; the rat reabsorbs with Oatp1a1, which
> the human lacks (§9.2).
>
> The practical rule: **use the rat for the magnitude of the tubular effect
> and never for its mechanism.** If a PHG argument depends on *which*
> transporter does the reabsorbing, the rat cannot support it. If it depends
> only on *how much* is reabsorbed, the rat is usable — with an uncertainty
> factor reflecting the 22–161-fold spread (§10.2), not a point estimate.

**Q2. The rat sex difference is a change of sign, not of degree. If sex can
flip the sign within one species, what is the defensible way to pick a single
human value?**

> **Suggested direction.** This is the strongest argument in the project
> *against* point estimates, and it should be used as such. Female rats
> net-secrete PFOA at `T = 24.6`; males net-reabsorb at `T = 0.42–1.07`; the
> mechanism is a 114-fold difference in one apical transporter (§9.7).
> Nothing about human physiology guarantees we sit at a stable point on that
> axis.
>
> Two defensible positions, and they are **not the same**:
> - **Health-protective:** take the slowest-eliminating configuration — lowest
>   f_u, strongest reabsorption. Defensible for a PHG.
> - **Central-estimate:** report the range and refuse the point. Defensible
>   for a dosimetry model feeding an epidemiological analysis.
>
> What is *not* defensible is a single value presented without saying which of
> these it is. **Ask the room which one the current number is meant to be** —
> often the answer is that nobody has decided.

**Q3. Which species should we actually be arguing about?**

> **Suggested direction.** The monkey is the highest-value and cheapest gap in
> the entire project (§16). Six compounds have clearance, none has binding,
> and the binding collection contains no macaque measurement of any kind.
>
> **If one experiment comes out of this meeting, argue for this one.**

### 17.2 On the unbound fraction

**Q4. Is a single f_u defensible at all any more?**

> **Suggested direction.** No, and this is the clearest correction available.
> The measured range across the series is **0.0004 to 0.31 — a factor of
> 780** — and OEHHA's 0.02 is not close to the middle of it for the compounds
> that matter most. For PFOA the measured human value is 0.0025, eight-fold
> below the assumption, so the filtered load used in the PHG is eight-fold too
> high and every reabsorption percentage derived from it is wrong in the
> denominator (§14, error 2).
>
> The correction requires **no new data** — only replacing one constant with a
> table that already exists (`db/human_fu_recommended.csv`).

**Q5. Variable f_u — in which of the three senses, and which does a PHG
need?**

> **Suggested direction.** Separate them explicitly; they have different costs
> and different readiness (§13.2):
> - **Across compounds** — do it now. Settled, and large (780×).
> - **Across individuals** — do it next. 2.5×, widest exactly where the short
>   chains sit. **But get the direction right:** higher f_u means faster
>   elimination, so the sensitive individuals are the **low**-f_u tail. A
>   model taking the upper percentile of f_u to be protective is protective of
>   the wrong person.
> - **With concentration** — not needed for human dosimetry at environmental
>   exposures (molar ratio ~0.0005 against a threshold near 0.1). It matters
>   only when reconciling in vitro datasets run at high loading.

**Q6. Low albumin as a susceptibility factor — is it real, and is it
actionable?**

> **Suggested direction.** Mechanistically it must be real: f_u is set by
> albumin and globulin concentrations and Fischer's equation reproduces it at
> R² = 0.99. **But the direction cuts against intuition.** Low albumin raises
> f_u, which raises the filtered load and *speeds* elimination. So
> hypoalbuminaemia predicts a **lower** PFAS body burden at steady state —
> while simultaneously raising the free concentration available to tissues at
> any given total concentration.
>
> **These two effects oppose each other, and which dominates depends on
> whether the toxicological driver is total body burden or free
> concentration.** That is genuinely open, and a good question to put to the
> room, because the answer determines whether low albumin is a susceptibility
> factor or a protective one.

### 17.3 On the transporters

**Q7. Does any measured transport survive contact with plasma protein?**

> **Suggested direction.** Unknown, and it is the largest unquantified caveat
> in the strand (§9.6). Transporters see the free concentration; an assay run
> protein-free presents a free fraction of 1 against an in vivo 0.0025 — a
> **400-fold** difference in the quantity the protein responds to.
>
> It may not invalidate the *rankings* — if binding affects all compounds
> similarly, much divides out — but it does call the *magnitudes* into
> question. **Ask whether anyone has run even one transporter assay at
> physiological albumin.** If not, that is the second experiment to argue for.

**Q8. Abundance versus affinity — which should gate inclusion in a model?**

> **Suggested direction.** Abundance gates, and it should gate *first*. The
> highest intrinsic clearance in the compilation (BCRP, 279 µL/min/mg for
> PFOA) belongs to a transporter two laboratories cannot detect in human
> kidney (§9.3).
>
> **A transfected-cell result states what a protein *can* carry; it says
> nothing about whether the protein is *present* to do it.** Applying that one
> rule dissolves two long-standing disagreements (ASBT, OSTα/β — Zhao found
> transport, Argoul found none, both right).
>
> Recommend an explicit two-gate rule for any model built from this
> compilation: **(1) is it expressed in the tissue, at what abundance;
> (2) does it transport the compound.** Papers routinely report (2) and omit
> (1).

**Q9. Is secretion the minor term we treat it as?**

> **Suggested direction.** No. Mouse PFHxA net-secretes at `T = 2.2`; rat
> PFHpA and female rat PFOA run well above 1. The framing of PFAS renal
> handling as "reabsorption, with some filtration" is a **long-chain
> generalisation that does not survive contact with the short chains.**
>
> Practical consequence: a model hard-coding reabsorption-only kinetics will
> be **wrong in sign** for the short chains — exactly the compounds now
> entering regulatory attention. **The structure of the model, not just its
> parameters, is chain-length dependent.**

**Q10. What is the minimum transporter dataset for a defensible IVIVE?**

> **Suggested direction.** Four things, and the strand currently has at most
> three for any transporter:
> 1. abundance in the target tissue, in the target species, by a quantitative
>    method (Thakur-style) — **not mRNA**;
> 2. intrinsic clearance for the compound, in a system where the expressed
>    amount is known;
> 3. a scaling factor from the expression system to the kidney — **OAT4 has
>    none** (§9.5);
> 4. the assay run at physiological binding, or a stated correction (§9.6).
>
> The two available scaling methods disagree 11-fold and in opposite
> directions, so item 3 is unsolved even where it exists. **Arguing for an
> agreed scaling convention may be worth more than more measurements.**

### 17.4 On meeting the regulatory numbers

**Q11. Two secondary-source errors have regulatory reach. What is the process
response?**

> **Suggested direction.** The specific corrections matter, but the process
> point is larger: **both errors are of a kind that reading the primary source
> catches immediately and that no amount of cross-checking between secondary
> sources would ever catch**, because the secondary sources copy each other.
>
> Suggest a standing rule: **any number entering a PHG calculation is read
> from the primary publication, and the primary is recorded alongside it.**
> This project found eight errors by applying exactly that rule, and in three
> cases the error was in a document that had already been through formal
> review.

**Q12. Where does the renal framework stop?**

> **Suggested direction.** Say it early (§11), so the tubular work is not
> asked to carry more than it can. The renal axis explains **species
> differences** well and the **absolute** human half-life only partially —
> 58–61% of human PFOA elimination is not renal. A mass-balance model closing
> only on renal clearance will be wrong by more than a factor of two.

**Q13. What would change our minds?**

> **Suggested direction.** Commit to §16's table before the data arrives. Note
> in particular that a *negative* result on the cynomolgus binding run — say,
> an f_u so high that the monkey turns out to be a poor human analogue — would
> be just as valuable as a positive one, and the project should say so now
> rather than after.

---

## 18. The notebook, the data and the code

### 18.1 One notebook, twelve figures

`notebook/pfas_tk_plots.ipynb`, generated by
`scripts/build_plots_notebook.py`.

```bash
cd tk_review
python3 -I scripts/build_plots_notebook.py        # writes the .ipynb
python3 -I scripts/build_plots_notebook.py --py   # same code as a flat script
cd notebook && jupyter lab pfas_tk_plots.ipynb
```

**Why a generator rather than a hand-written notebook.** A `.ipynb` is JSON
with the source chopped into cells: horrible to edit and impossible to review
in a diff. The cells live in the generator as ordinary Python strings, so the
diff is readable, and `--py` dumps the same code as a flat script that is
executed to prove the notebook runs before it ships.

| § | fig | what it shows |
|---|---|---|
| 8.2 | 1 | albumin dispersion by compound — the two-orders-of-magnitude choice |
| 8.3 | 2 | albumin affinity against total carbons — the PFNA maximum |
| 8.2 | 3 | albumin binding by assay family — method as a systematic axis |
| 8.2 | 4 | log *K* against assay loading — a hypothesis allowed to fail |
| 8b.1 | 5 | albumin against L-FABP — no affinity advantage |
| 8b.3 | 6 | L-FABP by compound — the two-site scatter |
| 9.2 | 7 | transporter abundance, three species — the reabsorptive side is not shared |
| 9.4 | 8 | OAT4 against the secretory transporters |
| 10 | 9 | `T` against total carbons, three species, with the f_u whiskers |
| 13.3 | 10 | where the choice of binding dataset flips the sign |
| 12 | 11 | is the gap clearance, or binding? |
| 15.3 | 12 | the candidates not yet admitted |

Each figure is one cell, reads one CSV, uses matplotlib and the standard
library only — **no pandas, no project imports** — so a cell can be lifted out
whole. Every function takes `save=`:
`fig2_chain_length(save="chain.png")`.

### 18.2 Editing

* **Carbon axis** — `total_carbons()` in the setup cell is the **only** place
  the perfluorinated→total rule is written. To switch back, make it
  `return int(n_cf)`.
* **Palette** — `C1`, `C2`, `C3` in the setup cell. They are validated for
  colour-vision deficiency (worst all-pairs ΔE 9.2 CVD, 24.0 normal); if you
  substitute your own, check them rather than eyeballing.
* **Figure 2's summary data** — `db/binding_chain_profile.csv` is a 24-row
  hand-editable file with commented provenance.
* **Thresholds** — e.g. `fig1_albumin_dispersion(min_n=2)`,
  `fig4_loading(["PFBA", "PFHxA", "PFDA"])`.

### 18.3 Nothing has been deleted

The project contains ~30 older figure scripts writing into `tk_review/figures/`,
covering strands this document does not use — mass balance, Qin 2024
comparisons, van 't Hoff thermodynamics, coverage matrices. **None has been
removed.** Section 6 of the notebook indexes them by theme with their file
names, so any can be pulled in. Run one with
`python3 -I scripts/<name>.py`.

> **Do not run `for f in scripts/*.py`.** That re-executes one-shot ingestion
> scripts and duplicates rows in `db/protein_binding.csv`.

### 18.4 The data files

| file | rows | what it is |
|---|---|---|
| `db/protein_binding.csv` | 1,497 | every binding measurement, all proteins and species |
| `db/albumin_kd_to_k.csv` | 119 convertible | albumin K_d → log *K*, with saturation and ratio diagnostics |
| `db/binding_chain_profile.csv` | 24 | figure 2's summary source — hand-editable |
| `figures/albumin_chainlength_data.csv` | 309 | per-measurement albumin values behind figures 1, 2, 5 |
| `figures/lfabp_by_compound_data.csv` | 87 | per-measurement L-FABP values behind figures 5, 6 |
| `db/transporter_direction.csv` | 305 | transporter × membrane × net process |
| `db/transporter_expression.csv` | 40 | Thakur and Uchida abundance, pmol/mg |
| `db/transporter_clearance_stage3.csv` | 52 | intrinsic clearances, tier-graded, with scaling columns |
| `db/tubular_factor_consolidated.csv` | 52 | **the authoritative `T` table** |
| `db/tubular_factor_evidence.csv` | 87 | every underlying clearance row, superseded values included |
| `db/human_fu_recommended.csv` | 28 | the f_u recommendation, with per-compound variability bands |
| `db/clearance_fu_coverage.csv` | 46 | which species-compounds have clearance, f_u, binding |
| `db/renal_clearance_leads.csv` | 10 | §15's candidates, **not admitted** |

### 18.5 Rebuilding the analysis

```bash
cd tk_review
python3 -I scripts/stage11_fujii_mouse_longchain.py    # must precede stage6
python3 -I scripts/stage11b_clearance_fu_coverage.py
python3 -I scripts/stage6_consolidate_tubular.py
python3 -I scripts/stage9b_build_recommended_fu.py
python3 -I scripts/stage12_renal_clearance_leads.py
python3 -I scripts/build_plots_notebook.py
```

Run in that order; each reads the outputs above it. All are idempotent and
verified byte-identical across two consecutive full runs.

---

## 19. One paragraph, if you read nothing else

PFAS renal elimination is the product of three terms — the free fraction, the
filtration rate, and what the tubule does — and the identity
`CL_R = f_u × GFR × T` is the whole framework. The binding literature that
sets the first term disagrees with itself by **240-fold**, more than species
differ, so **no species extrapolation factor for binding can be estimated from
it**. The transporter literature that explains the third term shows the
reabsorptive machinery is **not the same protein in human and rat**, so rodent
reabsorption data carries magnitude but not mechanism. L-FABP, the standing
explanation for hepatic accumulation, binds PFAS **no more tightly than
albumin does**, so that explanation needs replacing. The project's 1,497
binding measurements **cannot extend the renal endpoint any further**, because
the binding term is not the one that is missing — renal clearance is, and
roughly ten human compounds' worth of it sits in the literature un-ingested.
Three PDFs would land most of it, and one afternoon of bench work on
cynomolgus plasma would open the species that matters most.
