# Why PFAS half-lives differ between humans, rats and mice

A literature review built around one question: **is the species difference in
PFAS elimination a difference in how much is given, how much is held, or how
much is let go?**

The answer is the third, and it is not close. Everything else in this report is
either the evidence for that, the assumptions that obscured it, or the places
where the evidence runs out.

---

## 0. How to read this

Section 1 builds the foundation — one equation that governs the whole subject.
Everything after is an application of it. Sections 2–4 answer the three
questions you asked, in order. Sections 5–7 are the supporting evidence
(assumptions, binding, transporters). Section 8 lists what this review corrects
in the existing repository, section 9 what nobody knows yet. Section 10 is the
practical summary; section 11 takes the same machinery and asks which *endpoint*
a structure-activity model would have to predict, which turns out to change the
species answer too.

Code blocks are runnable from `tk_review/`. Every number traces to
`db/*.csv`, and every row there carries a PMID or DOI and the table it came
from.

There are questions at the end of several sections. They are there because the
reasoning is easy to follow and easy to half-follow, and the difference matters.

---

## 1. The foundation: one equation, three quantities, two of them free

Elimination kinetics rests on an identity:

```
t½ = ln2 · Vd / CL
```

- **Vd**, volume of distribution (L/kg) — the apparent volume the compound
  occupies. It is not an anatomical volume. It is the number that relates total
  body burden to plasma concentration.
- **CL**, clearance (L/kg/day) — the volume of plasma cleared of compound per
  unit time.
- **t½**, half-life — *not an independent quantity*. It is what the other two
  produce.

This matters more than it looks. **Any two of the three fix the third.** So a
paper reporting a half-life has either measured two things or measured one and
assumed another. Which it did decides how much to trust the number — and most
of the disagreement in this literature is an argument about an assumption, not
about data.

```bash
cd tk_review
python3 scripts/human_parameter_reconciliation.py
```

That script lays each major human source side by side and completes whichever
term it left implicit. The result is section 5.

### Two designs, two failure modes

| design | what it measures | what it must assume | fails when |
|---|---|---|---|
| **serum decay** — follow concentration after exposure stops | the elimination rate directly, from the slope of log C against t | that exposure actually stopped; that background is negligible or subtracted | exposure continues; levels approach background |
| **mass balance** — measure excretion and divide | clearance | a volume of distribution; *and* that you counted every excretion route | Vd is wrong; a route is missed |

**Serum-decay designs never need a Vd.** Mass-balance designs need both a Vd and
a complete excretion accounting, each uncertain by severalfold, and both
multiply straight into the answer. That single structural distinction sorts the
literature better than any judgement about authors or journals.

> **Check your understanding.** A paper reports a human PFOA half-life of 1.3 y
> from urinary clearance and an assumed Vd of 170 mL/kg. Someone objects that
> the true Vd is 430 mL/kg. Does the half-life go up or down, and by how much?
> *(Answer: up, by the Vd ratio — 2.5×, to 3.3 y. Half-life scales linearly with
> Vd. This is exactly what happened with Zhang 2013; see §8.)*

---

## 2. Your citation list: what it can and cannot support

The attached `PFAS_TK.xlsx` holds 81 citations. Each was resolved to a PMID
(title-matched against PubMed; the match check caught three wrong
identifications the automated search had made) and classified by evidence type.

```bash
python3 scripts/parse_xlsx_list.py && python3 scripts/resolve_pmids.py
python3 scripts/verify_pmids.py && python3 scripts/finalize_source_list.py
# -> db/source_catalogue.csv, 81 rows, 81 verified PMIDs
```

| category | n | can yield a TK parameter? |
|---|---|---|
| mechanistic / hazard only | 47 | no |
| in vivo toxicokinetics | 11 | yes |
| in vitro protein binding | 9 | yes |
| transporter or enzyme | 8 | yes |
| receptor binding | 3 | no |
| membrane biophysics | 3 | no |

**42 of 81 cannot contribute a half-life, Vd, clearance or binding constant.**
They are immunotoxicity, hepatic steatosis, oxidative-stress and
gene-expression studies — a hazard-characterisation corpus, which is a
different thing from a toxicokinetic one.

A further **21 are repeat-dose toxicity studies that measured serum or tissue
PFAS** without reporting kinetics. Those are flagged `supplies_internal_dose`,
because a measured internal concentration at a known administered dose is
exactly the dose↔serum pairing the exposure question needs, even with no
half-life attached.

The genuinely load-bearing entries are a small set: Kudo 2002 (sex-hormone
regulation of renal PFOA transport), Yang 2009 (Oatp1a1 reabsorption),
Nakagawa 2008 (OAT secretion), Han 2004 (α2u-globulin), Chang 2009 and
Chengelis 2009 and Lau 2020 and Das 2008 (rodent TK), and the albumin-binding
cluster (Luo 2012, Chen 2015, Wang 2016, Zhang 2009, Weiss-Errico 2018).

**So the review had to be built mostly from outside your list.** That is not a
criticism of the list; it is a statement about what it is for.

---

## 3. The species question: one axis, and the female rat falls off it

### 3.1 It is clearance, not distribution

Split the identity for every chemical × sex pair where both rodents were fitted:

```
ln(t_mouse/t_rat) = ln(Vd_mouse/Vd_rat) − ln(CL_mouse/CL_rat)
```

```bash
python3 scripts/rat_mouse_decomposition.py   # -> db/mouse_rat_decomposition.csv
```

| term | correlation with the half-life gap |
|---|---|
| volume of distribution | **r = +0.00** |
| clearance | **r = +0.87** |

Vd ratios stay inside ~4× for every chemical and sex. Clearance ratios reach
50×. Across all 46 species × sex × chemical groups in the database, **Vd spans
12.6× while clearance spans 7,384×**. Distribution is not where the action is.

### 3.2 The gap is specific to female rats — for long chains

| | female ÷ male clearance, **rat** | female ÷ male clearance, **mouse** |
|---|---|---|
| PFOA | **44×** | 1.8× |
| PFNA | **21×** | 1.1× |
| PFHxS | **20×** | 0.8× |
| PFOS | 0.8× | 1.0× |
| PFBA | 3.6× | **2.9×** |

Male mouse and male rat agree within ~4× for every chemical. The headline
15–31× "species" differences are **all in females, and all in long chains**.

The last row matters and I nearly missed it. **The mouse does show a PFBA sex
difference** — Chang 2008 measured 4.6–5.3×, as large as the rat's 5.2–6.2×. So
the mouse is perfectly capable of sex-dependent renal PFAS handling; what it
lacks is the sex dependence *for long-chain compounds*. Any explanation that
makes the mouse constitutionally sex-indifferent is wrong. And PFBA's sex
difference cannot be Oatp1a1-mediated, because C4 compounds do not interact with
Oatp1a1 at all (Yang 2009: no inhibition at 1 mM).

**How reliable are these ratios?** Less than their precision suggests, in three
specific ways:

- They are composites, `ln2·Vd,ss/CL`, not measured half-lives. For rat PFOA the
  measured male/female ratio is **8.5–71×** across five primary studies
  (Kudo 2002 **71×**, now read from the primary Table 2 — 5.68 d male against
  0.08 d female, and the paper's own text says "1/70"; Ohmori 2003 70×;
  Kemper 2003 39–63×; Kim 2016 8.5–12×) against the composite's 20×. The EPA
  fit's female-rat Vd,ss of 0.649 L/kg is 2–3× every other rat value in the same
  file, and measured female half-lives cluster at 1.9–4.6 h rather than the
  implied 16.6 h.
- **The mouse side rests on one study for PFHxS and PFNA** — Sundström 2012 and
  Tatum-Gibbs 2011 — while the rat side is replicated four to six times. The
  mouse PFNA male value of 227 d sits 3.3× above the measured ceiling of
  34.3–68.9 d, with a tenfold-wide credible interval. **For PFOA this caveat is
  now retired:** Argoul 2026 independently re-measured the female mouse and
  lands within 1.33× of Lou 2009 (§3.2a).
- A **route confound attacks the female rat specifically**: IV clearance
  significantly exceeded gavage clearance in female rats for six compounds
  (apparent bioavailability > 100%). EPA's PFHxS assessment responded by
  discarding all IV rat data; the fits producing these ratios pool IV and oral.

The direction and rough magnitude survive all three. The second decimal does not.

### 3.2a The same answer from primary measurements on both limbs

```bash
python3 scripts/primary_sex_species_decomposition.py
```

Three full texts now on disk measure both terms of the identity in both sexes,
inside single experiments, so strain, dose, route, assay and model are held
fixed: **Kudo 2002** Table 2 (Wistar rat, one IV dose, both sexes), **Lou 2009**
Table 2 (CD-1 mouse, both sexes, one model) and **Argoul 2026** Table 1 (female
CD-1 mouse, 11 PFAS in one cocktail).

| PFOA | rat (Kudo 2002) | mouse (Lou 2009) |
|---|---|---|
| clearance, female ÷ male | **44.3×** | **0.83×** |
| Vd, male ÷ female | 1.64× | 1.67× |
| half-life, male ÷ female | 71× | 1.39× |
| share of the log half-life ratio carried by clearance | 89% | 26% |

**The Vd sex ratio is the same in both species — 1.64× and 1.67× — while the
clearance sex ratio differs by 53×.** Sex-dependent distribution is conserved,
small, and in the same direction everywhere. Sex-dependent clearance is
species-specific and enormous. Note the mouse sign: the female mouse clears PFOA
slightly *slower* than the male, the opposite direction to the rat.

The mouse row is the same point seen from the other side. Where the clearance
difference nearly vanishes (0.83×), what remains of the half-life difference is
**74% Vd**. Vd matters only where clearance does not.

Taking the species ratio from primary sources on both sides:

| PFOA clearance | rat | mouse | ratio |
|---|---|---|---|
| **female** | 2,233 mL/kg-day | 4.5–6.0 | **373–496×** |
| **male** | 50.4 mL/kg-day | 7.2 | **7.0×** |

The species gap is **53× larger in females than in males**. The dose mismatch
runs the safe way: the rat was dosed at 20.14 mg/kg and the mice at 0.08–10, and
female rat PFOA clearance *falls* as dose rises (§4.2), so a dose-matched
comparison would widen the female gap, not narrow it.

**The conserved Vd ratio now rests on six datasets, not two.** Sundström 2012
measured PFHxS in rat, mouse *and* monkey, both sexes in all three, in one
laboratory — so the same split can be run four more times:

| dataset | compound | Vd M/F | clearance F/M |
|---|---|---|---|
| rat (Kudo 2002) | PFOA | 1.64× | **44.3×** |
| rat (Sundström 2012, Table 2) | PFHxS | 2.18× | **7.95×** |
| monkey (Sundström 2012, Table 5) | PFHxS | 1.35× | 1.45× |
| mouse (Lou 2009) | PFOA | 1.67× | **0.83×** |
| mouse 1 mg/kg (Sundström 2012, Table 3) | PFHxS | 1.34× | **0.91×** |
| mouse 20 mg/kg (Sundström 2012, Table 3) | PFHxS | 1.33× | **0.78×** |

**Vd male/female spans 1.33–2.18× — a 1.6× spread. Clearance female/male spans
0.78–44.3× — a 56× spread.** Three species, two compounds, four laboratories.
The Vd ratio is male-higher every single time and never leaves a narrow band;
the clearance ratio ranges over nearly two orders of magnitude and changes sign
between species. Sundström's authors reach the distribution half of this
independently, noting that in all three species "the mean Vdss suggested
predominantly extracellular distribution".

**The mouse's inverted sex difference is real, not noise.** The mouse clearance
ratio is below 1 in all three mouse rows — PFOA 0.83×, PFHxS 0.91× and 0.78×.
Three independent estimates, two compounds, two laboratories, all on the same
side of unity. The female mouse genuinely clears these compounds slightly *more
slowly* than the male, which is the opposite sign to the rat. It is a small
effect and nothing in this review turns on it, but it is not a sampling artefact
and it is one more thing any mouse-mechanism account has to accommodate.

**A data-quality note on the rat PFHxS row.** Sundström's Table 1 (24-hour
follow-up) carries male IV parameters estimated from **a single rat** and female
from two, with its own footnotes saying so; Table 2 (10 weeks, N = 4/sex) is the
usable one and is what the table above uses. Even there the female β-phase could
not be estimated and a one-compartment model was substituted. Regulatory
tabulations of rat PFHxS inherit these numbers without the footnotes.

**The rat PFOA sex ratio now has two primary sources that agree.** Ohmori 2003
measured all four of PFHpA, PFOA, PFNA and PFDA in both sexes:

| | C7 PFHpA | C8 PFOA | C9 PFNA | C10 PFDA |
|---|---|---|---|---|
| male t½ (d) | 0.10 | **5.63** | 29.5 | 39.9 |
| female t½ (d) | 0.05 | **0.08** | 2.44 | 58.6 |
| ratio M/F | 2.0× | **70.4×** | 12.1× | **0.68×** |

Its PFOA ratio of **70.4×** sits against Kudo 2002's 71× — two independent
laboratories, same strain-free conclusion. Note the last column: **PFDA is the
one compound where the female half-life exceeds the male**, so the "female rats
eliminate PFAS faster" generalisation fails at C10 even within the rat.

Ohmori also states this review's §3.1 conclusion outright, in 2003: *"Distribution
volumes in steady state (Vss) were not much different between PFCAs and between
sexes"*, and in the Discussion, *"The difference in Vss between PFCAs was not so
significant as the difference in t½."* It further reports a correlation between
total and renal clearance of **r² = 0.981** — total clearance in the rat
essentially *is* renal clearance. (Its Table 1 is an image in the PDF, so the
numeric Vss values could not be extracted; the half-lives above are from the
abstract.)

**The whole pattern replicates in a second compound.** Tatum-Gibbs 2011 is the
only strain-matched rat-versus-mouse PFNA experiment, run in one laboratory with
both sexes of both species:

| PFNA half-life | rat (Sprague-Dawley) | mouse (CD-1) |
|---|---|---|
| male | 30.6 d | 34.3–68.9 d |
| female | **1.4 d** | 25.8–68.4 d |
| sex ratio M/F | **21.9×** | ~1.0–1.3× |
| species ratio, same sex | — | male 1.1–2.3×, **female 18–49×** |

Rat sex difference large, mouse sex difference absent, species gap concentrated
in females — the identical structure to PFOA, in a different compound, with the
strain confound removed. Its abstract also records that mouse PFNA elimination
is "non-linear with exposure dose" while the rat's is "by and large linear",
which is the §4.2 dose-dependence pattern with the species roles reversed from
PFOA and is not explained here.

Two incidental confirmations from this paper: its measured mouse ceiling of
**68.9 d** is the figure against which §3.2's composite of 227 d sits 3.3× too
high, and it independently cites Lou 2009 as "21.7 days for males and 15.6 days
for females", confirming that those two numbers are the sexes and not a dose
range.

**And the Vd-versus-clearance question settles inside one experiment.** Argoul
2026 dosed 11 PFAS as a single cocktail to one sex of one strain:

| | span across 11 PFAS |
|---|---|
| plasma clearance | 1.3 → 6,830 mL/kg-day = **5,254×** |
| Vss | 0.062 → 0.48 L/kg = **7.7×** |

PFHxA is the one exception at Vss 4.0 L/kg, and the authors attribute 3.8 L/kg
of that to a deep peripheral compartment, putting it at 0.12 L/kg without it.
Even counting PFHxA at face value, clearance varies over a range 81× wider than
Vss. §3.1 reached this conclusion across studies; it holds with every
cross-study confound removed.

This study's *absolute* clearances run about 2.5× below the
one-compound-at-a-time literature (§11.3a), but the comparison above is a ratio
within one experiment, and a shared scaling cancels from a ratio. The 5,254×
against 7.7× stands whatever the common factor turns out to be.

### 3.3 The axis that unifies all of it

OEHHA's PHG appendix (Table A6.4, pp. 327–328, adapted from Han et al. 2012)
tabulates GFR, measured renal clearance and net tubular reabsorption of PFOA
across seven species. Every group sits on one scale:

| | % of filtered PFOA reabsorbed |
|---|---|
| human | 99.8% |
| mouse male / female | 97.0% / 95.2% |
| **rat male** | **93.2%** |
| macaque male / female | 91.2% / 81.2% |
| dog male / female | 59% / 52% |
| **rat female** | **net secretion** |
| rabbit male / female | net secretion |

**The female rat does not reabsorb less — it secretes.** Its renal clearance
(666–1,009 mL/kg-day) exceeds its entire filtered load (288 mL/kg-day). The
mouse sits at 95–97% in *both* sexes. The 31× PFOA "species difference" is one
group crossing from net reabsorption into net secretion.

Because `t½ ∝ 1/(1−FR)`, small changes near FR = 1 produce enormous changes in
half-life. Going from 97% to 99.8% reabsorbed cuts the escaping fraction
fifteenfold.

```bash
python3 scripts/reabsorption_axis.py   # -> db/reabsorption_axis.csv
```

Using **only** measured renal clearance and fitted Vd — no free parameters —
the identity reproduces every rodent group to within **1.7×** across a 37×
half-life span:

| | Vd (L/kg) | CL_renal (L/kg/d) | predicted t½ | observed t½ | ratio |
|---|---|---|---|---|---|
| rat female | 0.649 | 0.838 | 0.5 d | 0.7 d | 0.78 |
| rat male | 0.289 | 0.0197 | 10.2 d | 14.2 d | 0.72 |
| mouse female | 0.299 | 0.0160 | 12.9 d | 21.5 d | 0.60 |
| mouse male | 0.252 | 0.0100 | 17.5 d | 25.6 d | 0.68 |

**The axis, read at its origin.** OEHHA's Table A6.4 says it is adapted from
Han et al. 2012, and that review is now on disk, so the chain can be checked:

```bash
python3 scripts/han2012_axis_at_source.py
```

The adaptation preserved nine of eleven rows and altered two: **human 99.94% →
99.8%** and **male rat 93.7% → 93.2%**. Neither changes a conclusion; the
primary figures are used from here on. Han also calls its own values "rough
estimates" in the body text, a caveat the citation chain dropped.

More importantly, Han states a point this review had missed. **Humans have the
longest PFOA half-life but not the largest reabsorption.** In absolute terms
humans reabsorb 51 mL/d/kg against the male rat's 270 and the mouse's 318–324.
What is extreme about humans is the *fraction*. That is because the axis has two
factors, not one:

```
CL_renal = [fu · GFR] × [1 − FR]
            filtration   escape fraction
```

and both vary across species. Decomposing the 333× male-mouse-to-human gap in
renal clearance:

| factor | contribution | basis |
|---|---|---|
| escape fraction | **50×** | 3.0% escapes in mouse vs 0.06% in human |
| filtration term | **6.5×** | GFR 16.7 vs 2.57 L/d/kg |

On a log scale reabsorption carries **68%** and GFR **32%**. So the single-axis
framing is justified — reabsorption is the bigger term — but a third of the
human/rodent difference is simply that humans filter blood about six times more
slowly per kg, which no transporter story explains and which the allometric
scaling in §3.6 would partly absorb.

**An independent dataset reproduces the axis.** Everything above traces to one
compilation (Han 2012, via OEHHA Table A6.4). Argoul 2026 Table 4 reports the
three quantities the axis needs — unbound fraction, free filtration clearance,
measured renal clearance — for ten PFAS in female mice from its own laboratory,
so `FR = 1 − CL_renal/(fu·GFR)` can be recomputed from scratch:

```bash
python3 scripts/argoul_reabsorption_check.py
```

| chemical | fu | FR recomputed |
|---|---|---|
| PFDA | 0.42% | +0.997 |
| PFHxS | 1.3% | +0.996 |
| PFNA | 0.35% | +0.983 |
| **PFOA** | 0.87% | **+0.959** |
| PFOS | 0.25% | +0.950 |
| PFBA | 77% | +0.898 |
| GenX | 26% | +0.897 |
| PFO2OA | 10% | −0.533 |
| PFHxA | 25% | −1.712 |

Mouse PFOA comes out at **95.9%**, inside the 95.2–97.0% the table above already
uses, from a completely separate dataset. Two compounds fall past zero into net
tubular secretion — the same end of the axis OEHHA assigns to the female rat and
the rabbit. So the secretion end is not a rat-female peculiarity; the mouse
reaches it too, and which compounds get there is chain- and head-group-specific
rather than species-specific.

The relationship to binding is strong but **one-directional** (Spearman
ρ = −0.60, n = 9). Every compound bound above 98.5% sits above 95% reabsorption,
without exception. The loosely bound compounds spread from +0.90 to −1.71, so a
high unbound fraction is compatible with either. Tight binding is close to
sufficient for near-complete reabsorption; loose binding predicts nothing on its
own. That asymmetry is convenient, because it is the tightly bound end — where
humans sit — that behaves predictably.

### 3.4 Two honest caveats

**The "99.94%" is softer than it looks.** It depends on an unbound fraction
assumed to be 0.02 for every species — **Han 2012's assumption**, stated in its
own Table 4 footnotes, which OEHHA inherited rather than introduced (this report
previously attributed it to OEHHA). Fischer et al. 2024 measured human
serum PFOA f_unbound at **0.00061** by solid-phase microextraction — 33× lower.
Recomputed, the figure is 96.2%. The qualitative claim survives; the headline
number does not. The *prediction* above is unaffected, because `fu · GFR · (1−FR)`
is just `CL_renal` rearranged and f_unbound cancels.

**In the rat, the mechanism is established.** The sex difference is renal and
androgen-driven at the apical reabsorptive step, and the evidence is
interventional rather than correlational: castration raises male PFOA renal
clearance **14-fold** to the female level and testosterone reverses it
(Kudo 2002); oestradiol raises male 24-h urinary excretion from 9 ± 4% to
61 ± 19% of dose (Ylinen 1989); probenecid works only in castrated males
(Vanden Heuvel 1992). Rat renal Oatp1a1 is apical, male-dominant,
testosterone-induced (Lu 1996), and transports PFOA (Km 126–162 µM;
Weaver 2010, Yang 2009). Worley 2015's PBPK puts the whole difference in one
parameter, and Gotoh 2002 showed the mechanism can generate > 250× using a
non-PFAS Oatp1a1 substrate with GFR and protein binding held equal.

Kudo 2002's full text adds the cleanest single piece of that evidence.
Probenecid collapses renal PFOA clearance in male, castrated-male and female
rats onto a common value with **no significant difference between the three
groups** (Table 3: 0.032 / 0.447 / 0.732 → 0.010 / 0.013 / 0.016 mL/min/kg).
The entire sex difference is carried by organic anion transport; what remains
when transport is blocked is sex-indifferent.

**Which transporter, though?** The full text makes this an elimination argument
rather than an assertion. A candidate must be *both* sex-divergent *and* able to
carry PFOA:

| transporter | sex difference (Kudo 2002) | transports PFOA? | verdict |
|---|---|---|---|
| Oat1 | none among males, slightly lower in females | yes, Km 43.2 µM (Weaver 2010) | transports, not sex-divergent |
| Oat3 | flat except ~2× after ovariectomy | yes, C8–C9 (Weaver 2010) | transports, not sex-divergent |
| Oat2 | **7.5× female-predominant** | **no** (Weaver 2010 Fig. 9A) | sex-divergent, cannot carry it |
| OAT-K | 2.5× male-predominant | untested | untested |
| oatp2 | none | no | out |
| **Oatp1a1 (oatp1)** | **23× male-predominant**, androgen-induced | **yes**, C8–C10 | **the only candidate meeting both** |

Worth recording that **Kudo 2002's own conclusion was different**: its multiple
regression put OAT2 and OAT3 as positive correlates of renal clearance, with
OAT3's coefficient the larger, and the paper concludes "both OAT2 and OAT3 are
responsible for the urinary elimination" of PFOA. On oatp1 it is explicitly
agnostic — "it remains unclear whether oatp1 is responsible for the
re-absorption of PFOA" — because the literature was then split on whether oatp1
moved anions in or out. Direct transport assays a decade later resolved both
halves against it: Oat2 carries no PFOA, so its correlation across hormonal
manipulations is co-regulation rather than mechanism, and Oat3, which does carry
PFOA, is not sex-divergent. The paper everyone cites for the female-rat
mechanism proposed a mechanism its own successors refuted, and the surviving
explanation is the one it declined to endorse.

Two further confirmations now sit at source rather than second-hand:

- **Oat2 carries no PFAS in humans either.** Louisse 2024 tested PFHpA, PFOA,
  PFNA, PFBS, PFHxS and PFOS against OAT1-, OAT2- and OAT3-transduced human
  cells and found **no transport by OAT2 for any of them**. With Weaver 2010
  (rat Oat2) and Nakagawa 2008, that is three independent negatives across two
  species. The transporter with the largest sex difference is the one that
  cannot carry the compound, and that is no longer a rat-specific observation.
- **The rat-versus-mouse difference really is in Oat2.** Buist & Klaassen 2004,
  now read directly: "the most notable species differences in Oat mRNA
  expression were **a lack of Oat2 female predominance in mouse kidney** and a
  less dramatic Oat3 male predominance in mouse liver." The one clear species
  difference is in the transporter that does not transport PFOA — which is
  precisely why it cannot explain the species gap.

And **Yang 2009**, the source of the rat Oatp1a1 kinetics this report had been
quoting through Weaver 2010's Discussion, confirms them directly: Km =
**162.2 ± 20.2 µM**, no inhibition by C4 or C5 at 1 mM, and an in vitro–to–in
vivo relationship of `log(CL_total) = 2.233·log(Ki,app) − 2.627` with
**R² = 0.982** across C6–C10 against male-rat total clearance. It adds one
methodological caveat: total PFOA uptake into Oatp1a1-expressing cells does
*not* saturate even at 1 mM, because a large linear passive-diffusion component
runs alongside — the saturable Km emerges only after subtracting the
vector-control rate.

**Three competing hypotheses are dead** and should be stated as such:

- **α2u-globulin** — Han 2004 measured Kd ~10⁻³ M and concluded it "cannot
  adequately explain the sex-dependent elimination of PFOA in rats". The mouse's
  lack of α2u-globulin is a red herring.
- **L-FABP / I-FABP** — Modaresi 2025's double-knockout mouse showed no change
  in PFOS half-life, clearance or Vd.
- **Oat2** — the most sex-divergent rat renal transporter (male = 13% of
  female, i.e. 7.5× female-predominant; Kudo 2002 Fig. 3B) does not transport
  PFOA at all. Nakagawa 2008 concluded this, and Weaver 2010 Fig. 9A — full text
  on disk — shows no significant net Oat2-mediated uptake of C7, C8, C9 or C10
  in Oat2-expressing CHO cells.

**The rat Oatp1a1 fold-value contradiction is resolved, and EPA has it wrong.**
EPA's PFHxA review states the rat male/female renal Oatp1a1 mRNA ratio is
**2.5-fold**; EPA's PFOA assessment says **5–20-fold**. Neither traces to
anything in Kudo 2002's abstract, and this report previously refused to quote any
single value. The primary text gives all three numbers, and they belong to
different transporters (p. 310, Figs. 3–4):

| | male ÷ female renal mRNA |
|---|---|
| **oatp1 (Oatp1a1)** | **23×** |
| OAT-K | 2.5× |
| OAT2 | 0.13× (7.5× the other way) |

So **23× is the Oatp1a1 figure**. EPA's 2.5-fold is this paper's **OAT-K**
ratio, a different transporter, apparently misattributed; EPA's 5–20-fold
understates the primary value. In oestradiol-treated males and in ovariectomised
females, oatp1 mRNA was undetectable at 21 PCR cycles — the regulation is close
to on/off, not graded, which is why a 23× mRNA swing can support a 71×
half-life difference.

**In the human, the reabsorptive transporter is a different protein.** This is
the most consequential thing the new full texts changed. Yang 2010 tested the
three candidate human apical transporters and found that **OATP1A2 — the closest
human orthologue of rat Oatp1a1 — does *not* mediate saturable PFOA uptake at
all.** What does: **OAT4** (Km 172.3 µM at pH 6, 310.3 µM at pH 7.4) and
**URAT1** (Km 64.1 µM), with chain-length-dependent inhibition across C4–C12 for
both.

So the axis is real at the level of fractional reabsorption, and §3.3 shows it
predicts half-lives across seven species — but it is **implemented by different
proteins in rat and human**. The rat's reabsorptive step runs through an
androgen-induced transporter; the human's runs through two that are not
androgen-regulated. That is a mechanistic reason to expect what §5.3 reports
epidemiologically: no human sex difference of the rat's kind, and whatever human
sex difference exists arising from something else. It also means rat-to-human
extrapolation of the *reabsorption mechanism* — as opposed to the reabsorbed
fraction — has no molecular warrant.

One caveat on URAT1 that resolves a contradiction §4.3 flagged. Yang's 64.1 µM
was measured **in the absence of extracellular chloride**; the paper says URAT1
uptake "was greatly enhanced by an outward Cl⁻ gradient". URAT1 is an exchanger,
so under physiological chloride that transport is much weaker. Louisse 2023
found URAT1 transported nothing, and the two results are probably not in
conflict at all — they are the same transporter under different counter-ion
conditions. Han 2012's own Table 6 lists a human URAT1 row with **no Km value
entered**, which suggests its authors reached the same judgement about the
number's physiological relevance. Yang's conclusion that URAT1 "may contribute
significantly to the long half-life of PFO in humans" is weaker than it reads.

**The Cheng 2005 / Cheng 2006 contradiction is resolved, against Cheng 2006.**
This report previously flagged the direction of mouse renal Oatp1a1 as
unresolved, because Cheng 2006's abstract calls it "female-predominant". With
Cheng 2005's full text in hand and Cheng 2006's abstract recovered in full, the
2006 sentence is simply wrong, and it contradicts its own paper:

- **Cheng 2005** (same laboratory, primary data, C57BL/6, n = 10/sex) states it
  twice — abstract: "In kidney, expression of Oatp1a1, 3a1, and 4c1 was higher in
  males than in females"; Discussion: "Gender differences in Oatp1a1 expression
  were observed in mouse liver and kidney, with higher expression in males than
  in females." It also confirms the transporter's role: Oatp1a1 "is localized to
  the apical membrane domain of proximal tubules in kidney, where it reabsorbs
  organic anions from the lumen", and that male predominance "is
  androgen-dependent (Lu et al., 1996; Isern et al., 2001)".
- **Cheng 2006** says in its background sentence that renal Oatp1a1 is
  female-predominant, then reports that "**androgens increased Oatp1a1 mRNA in
  liver and kidney**" and concludes that "in kidney, gender-divergent Oatp
  expression is exclusively caused by stimulation by androgens". If androgens
  raise it and androgens are the exclusive cause, it cannot be female-predominant.

- **Cheng & Klaassen 2009** — the same laboratory, five years later, *citing
  Cheng 2006 itself* — writes: "male-predominant Oatps (**Oatp1a1** and 3a1) and
  female-predominant Mrp3 in mouse kidneys are due to stimulatory effect of
  androgens and estrogens, respectively (Cheng et al., 2006; Maher et al.,
  2006)." The authors themselves read their 2006 paper as male-predominant.

So mouse renal Oatp1a1 is **male-predominant and androgen-driven**, same
direction as the rat. Three lines of evidence, one of them the authors' own
later reading of the contradictory paper. Cheng 2006 adds one clean mechanistic detail worth keeping:
male-pattern growth hormone raises hepatic Oatp1a1 but **not** renal, so the
renal sex difference is purely androgenic.

**That deepens the mouse puzzle rather than solving it.** Mouse renal Oatp1a1 is
*also* androgen-dependent, so its regulation is not rat-specific and cannot by
itself generate the species difference. The only clean measured rat-vs-mouse renal comparison (Buist &
Klaassen 2004) found the species difference in **Oat2 — the wrong transporter**.
**No mouse orthologue transport kinetics exist for any PFAS.** This remains the
single largest hole in the mechanism. And one quantitative comparison is still
impossible: Kudo 2002 puts the rat renal Oatp1a1 male/female ratio at 23×, but
Cheng 2005 reports the mouse renal direction without a fold value — the only
renal number it quantifies is Oatp3a1 at **3.3×**. So whether the mouse's
androgen-driven Oatp1a1 difference is 23× like the rat's or much smaller is
unknown, and that single number might settle the whole species question. A
strain caveat sharpens it: Cheng used C57BL/6, while the PFAS kinetics come from
CD-1.

One thing this review previously listed as missing is not: **a sex-resolved
mouse PFOA half-life pair has been measured.** Lou 2009 Table 2 reports female
15.6 d (95% CI 14.7–16.5) and male 21.7 d (19.5–24.1), with Vd 0.135 against
0.226 L/kg and ke 0.00185 against 0.00133 /h. The "15.6–21.7 d range" EPA's PFOA
assessment quotes from this paper is not a dose range — it is the two sexes. The
pair exists, it is tight, and it says the mouse sex difference in PFOA half-life
is **1.39×** against the rat's 71×.

One plausible reconciliation, flagged as inference rather than finding: mouse
elimination is multi-route and compensating. Furukawa 2024's Abcb4-null mice
lose biliary PFOA but gain renal clearance via upregulated Abcb1, and
Abcb1a/b-null mice lose renal clearance. A redundant system is structurally
less sensitive to any single androgen-regulated step.

**So: the phenomenology is one axis and is quantitative. The mechanism is
established for the male rat and absent for the mouse.** That gap is itself a
result, and it is the highest-value experiment in §9.

### 3.5 PFOS is not exempt for the reason it looks

The obvious reading of the PFOS row — sulfonate head group escapes a pathway
that acts on carboxylates — **is wrong, and measured data refute it.** Rat
Oatp1a1 transports PFOS with ~7× *higher* affinity than PFHxS (Km 37 vs 256 µM),
and Zhao 2017 concluded outright that "rOATP1A1 is not the major transporter for
the renal elimination of PFOS".

What PFOS and PFDA actually share is that **renal excretion has become
negligible** — Kudo 2001 recovered 0.2% of a PFDA dose in urine over 120 h — so
there is nothing left for a reabsorptive transporter to modulate. The sex
difference needs a *window*: substantial renal elimination **and** adequate
Oatp1a1 affinity. C8–C9 carboxylates and the C6 sulfonate sit inside it; C4–C5
are not Oatp1a1 ligands at all; C10 carboxylate and C8 sulfonate are renally
silent.

> **Check your understanding.** Male mouse PFOS half-life is 35.5 d, male rat
> 57.1 d — the mouse is *faster*. But for PFOA the mouse is *slower* (25.6 vs
> 14.2 d). Why does the direction flip between two compounds in the same pair of
> species? *(Because the rat's androgen-regulated renal pathway acts on PFOA and
> is irrelevant for PFOS — not because PFOS fails to bind the transporter, but
> because PFOS barely leaves by the kidney at all. Remove that pathway and the
> two species sit within ~2× of each other.)*

### 3.6 Does body size explain any of it? Partly, and not for PFOA

If the species differences were a size effect, allometric scaling would absorb
them. Argoul 2026 tests this directly, scaling its female-mouse clearances to
human and comparing against measured human clearances (Abraham 2024):

| chemical | mouse CL | human CL measured | measured ÷ predicted |
|---|---|---|---|
| PFOS | 2.68 | 0.088 | **1.01** |
| PFNA | 2.1 | 0.067 | **0.98** |
| PFHxS | 1.34 | 0.053 | **1.22** |
| PFBA | 783 | 21.8 | **0.86** |
| PFDA | 1.78 | 0.148 | 2.56 |
| **PFOA** | 4.54 | **0.042** | **0.29** |
| PFHxA | 6830 | 56.9 | 0.26 |
| GenX | 274 | 43.0 | 4.86 |

All in mL/kg-day; predictions use the exponent −0.4138 fitted to these data,
steeper than either textbook value (−0.25, −0.33).

For four of eight compounds allometry lands within 1.25× — which is a real
result, and more than this review expected. **PFOA is one of the failures**, and
it fails by 3.4× in the direction of humans clearing it more slowly than size
alone predicts. So for the one compound the species question is usually asked
about, body size explains none of the gap: the residual is exactly the kind of
term §3.3's reabsorption axis supplies.

Two cautions. The comparison is female mouse against mixed-sex human, and §3.2a
shows sex costs a factor of 0.83 in the mouse, so that is not the explanation for
a 3.4× miss. And the human column is Abraham's — the low side of the 6.7×
human-clearance disagreement in §5.1. Argoul adopting **0.042 mL/kg-day** for
human PFOA is itself a data point in that dispute: it is 6.7× below OEHHA's
0.280, and the ratio is the disagreement, not a coincidence. Use the OEHHA value
instead and PFOA's measured-over-predicted becomes ≈1.9, overshooting rather
than undershooting. The allometric verdict for PFOA therefore depends on which
side of §5.1 you stand, and is not independent of it.

---

## 4. Does exposure relate to half-life or to volume of distribution?

**To neither, materially.** Six independent lines agree, and they fail in
different ways, which is what makes the convergence persuasive.

### 4.1 Vd does not track dose or concentration

```bash
python3 scripts/vd_vs_dose.py
```

Slope of ln(Vd) on ln(serum), within species × sex across chemicals:

| group | slope | r | concentration range |
|---|---|---|---|
| primate male / female | −0.21 / −0.15 | −0.76 / −0.89 | 280–380× |
| rat male / female | −0.13 / −0.02 | −0.45 / −0.06 | 21–46× |
| mouse male / female | −0.14 / +0.25 | −0.40 / +0.30 | 14–74× |

Inconsistent in sign, weak, and the same on dose rather than serum. The
albumin-saturation prediction (saturated plasma binding → more compound in
tissue → larger apparent Vd) **fails from both ends**: Kudo 2007 measured 52% of
a 0.041 mg/kg IV PFOA dose in rat liver versus 27% at 16.56 mg/kg, so hepatic
uptake saturates first and apparent Vd *decreases* with dose. And human serum
albumin is ~410 µM while human PFOS at 5 ng/mL occupies **0.0024% of binding
sites**.

### 4.2 Half-life is nearly dose-independent — except in one group

Over 8–25× dose in the NTP rat studies (Huang 2019), half-life moves 12–20% at
most — while *apparent* clearance moves 1.5–2.9×. In four of six series the
clearance change is matched within a few percent by the opposite change in
bioavailability (PFHxS male: CL +1.87×, F 98%→52% = −1.88×). **The apparent
dose-dependence of clearance in those studies is a bioavailability artefact.**
EPA's PFHxA IRIS review tested dose-dependence formally and rejected it; Lau
2020 states PFBS clearance "was linear with exposure doses".

**The exception is the female rat, and it is large.** Kemper 2003 dosed PFOA at
four levels over a 250× range in both sexes:

| dose (mg/kg) | male t½ (h) | female t½ (h) | ratio |
|---|---|---|---|
| 0.1 | 202 | 3.2 | 63× |
| 1 | 138 | 3.5 | 39× |
| 5 | 174 | 4.6 | 38× |
| 25 | 157 | **16.2** | **10×** |

Over that 250× dose rise the male half-life moves 0.78× (slope −0.05) while the
female moves **5.1× (slope +0.29)**. The sex ratio collapses from 63× to 10×.

This is the clearest dose-dependence anywhere in the dataset, and it is
sex-specific — which makes mechanistic sense once the directions are kept
straight:

| group | regime | what saturates | effect of raising dose |
|---|---|---|---|
| male rat | reabsorption-dominated | reabsorption | clearance **↑** (slope +0.11) |
| female rat | secretion-dominated | secretion | clearance **↓** (slope −0.29) |

Saturating *reabsorption* lets more escape, so clearance rises. Saturating
*secretion* lets less be pumped out, so clearance falls. Both are observed, in
the predicted directions, in the same species and compound.

**Two consequences.** First, pooling doses blurs the sex difference, which is
part of why the EPA composite gives 20× where matched-dose measurements give
~70×. Second — and this is the one that matters for extrapolation — dose
dependence is strong only where a high-capacity secretory pathway exists to
saturate. Male rats, mice, monkeys and humans are all
reabsorption-dominated, and there the effect is the +0.1 of §4.4.

### 4.3 Saturation is unreachable at human exposure

Measured transporter Km values sit **2,200× to 82,000×** above current human
serum concentrations (NHANES geometric means: PFOA 1.56 ng/mL = 3.8 nM; PFOS
4.72 ng/mL = 9.4 nM). The tightest margin anywhere — C8 Health Project Little
Hocking mean 227.6 ng/mL against human URAT1 Km 64.1 µM = 26,542 ng/mL — is
**~117×**, at which the Michaelis–Menten term is 99.1% linear. Tubular water
reabsorption concentrates luminal fluid ~100-fold, closing two of the four
orders; a 20–800× margin remains.

**Two caveats, pulling opposite ways — and the second is serious.**

```bash
python3 scripts/saturation_margin_km_vs_kt.py
```

*The tightest margin above uses a number that may not apply in vivo.* Yang
2010's human URAT1 Km of 64.1 µM was measured **in the absence of extracellular
chloride** (§3.4). Under physiological chloride that transport is much weaker,
and Louisse 2023 found none at all. Drop URAT1 and the binding margin is set by
OAT4 at 172.3–310.3 µM, which *widens* the tightest margin from ~117× to
310–840×. This caveat strengthens the section's conclusion.

*The PBPK models agencies rely on imply the opposite conclusion.* Han 2012
Table 7 compiles the transport-affinity constants that published PFOA PBPK
models actually run on. The human value is **KT = 0.055 mg/L = 0.133 µM**.
There are now **six independent in vitro Km values for PFOA against human
transporters**, from two laboratories:

| transporter | Km (µM) | source |
|---|---|---|
| OAT4 | 47 | Louisse 2023, via Louisse 2024 Table 1 |
| URAT1 | 64.1 | Yang 2010 (zero extracellular chloride) |
| OAT3 | 90 | Louisse 2024 Table 1 |
| OAT4 | 172.3 | Yang 2010 (pH 6.0) |
| OAT1 | 185 | Louisse 2024 Table 1 |
| OAT4 | 310.3 | Yang 2010 (pH 7.4) |

They span 47–310 µM — a factor of 6.6 — while the fitted KT sits
**354–2,336× below all of them.** At that value the margins invert:

| human serum PFOA | ÷ PBPK KT (0.133 µM) | ÷ lowest in vitro Km (64.1 µM) |
|---|---|---|
| general population, ~4 ng/mL | 0.1× | 0.0002× |
| Lubeck WV, 68 ng/mL | 1.2× | 0.003× |
| Little Hocking OH, 448 ng/mL | **8.1×** | 0.017× |
| occupational, ~1000 ng/mL | **18.2×** | 0.038× |

So on the PBPK parameters, reabsorption is already **saturated** in the
contaminated communities — which would make clearance dose-dependent there and
mean a single clearance factor cannot transfer between exposure settings. That
is the opposite of this section's conclusion, from the same literature.

**The in vitro side should be believed, and the case is now stronger than when
this section was first written.** Six measurements from two laboratories agree
within a factor of 7; the KT values are fitted to plasma curves, not measured.
Three further reasons, worth stating rather than assuming. (a) The KT values are not measurements: Han's
footnote says "Tmc and KT are obtained by fitting PFOA plasma elimination
curves", and the two rows in the same table marked "derived from an in vitro
measurement" carry KT = 67 mg/L — about 1,200× higher, and squarely inside the
in vitro Km range. (b) The fitted values are poorly identified: the mouse row is
Lou 2009 Table 4, where Tm = 860.9 ± 1298.3 and KT = 0.0015 ± 0.0022 both have
standard errors exceeding the estimate, and its Tmc/KT of 1.4 × 10⁷ L/d/kg is
four orders off every other species in the table. (c) The dose-response evidence
agrees with the in vitro side: §4.2 finds human half-life essentially
dose-independent and §4.6 finds the between-person association does not survive
age adjustment, whereas saturation at 68–448 ng/mL would make half-life rise
with exposure across exactly that range.

**So the answer stands, but its basis is narrower than the section originally
claimed.** Saturation is unreachable at human exposures according to every
direct measurement of the transporters; the contrary implication of the PBPK
constants is an artefact of fitting a saturable model to data that do not
constrain its parameters. What would settle it is a measured human KT, and
nobody has one.

### 4.4 The four slopes that do exist all land in the same place

| contrast | slope d ln(CL)/d ln(C) |
|---|---|
| male rat PFOA, 250× dose range (this repo) | +0.11 |
| human cross-cohort PFOA (this repo) | +0.09 |
| Abraham tracer dose vs cohorts | +0.08 |
| three PFBS studies over 1,200× | +0.12 |

Against a mechanistic ceiling of 1. Applied to the 187× human/rat serum gap,
+0.1 predicts a 1.6× half-life difference; the observed gap is 81×. **Dose
explains ~13% on a log scale, and that is generous.**

### 4.5 The confounder that manufactures the opposite conclusion

This is the most useful methodological finding in the review. Unsubtracted
background inflates a fitted half-life, and the inflation scales with how close
the cohort sits to background:

| serum ÷ background | inflation | source |
|---|---|---|
| 1.7 | **+150%** (5.0 → 2.0 y) | Nilsson 2022a, PFOA |
| 5.5 | +57% | Lyu 2026 |
| 13 | +20% | Xu 2020 |
| 49 | +15% | Batzella 2024 |
| 260 | +0.7% | Li 2022, PFHxS |

**Any cross-cohort regression of half-life on exposure that pools corrected and
uncorrected studies produces a spurious negative slope that is
indistinguishable from saturable reabsorption.** Low-exposure cohorts get
inflated half-lives; high-exposure cohorts do not. That alone can generate the
"higher exposure → shorter half-life" pattern with no mechanism behind it.

The limiting case is Gasiorowski 2022's randomised *control* arm: at serum PFOS
~11 ng/mL, 52 weeks produced −0.01 ng/mL (p = 0.96) — an effectively infinite
apparent half-life — with PFOA actually *rising* (p = 0.02).

And every between-person direct test gives the wrong sign for saturation:
Batzella 2024 (n = 5,860) finds PFOA half-life rising from 1.92 y in the lowest
baseline quartile to 2.85 y in the highest. The authors name the artefact
themselves — people with faster excretion have lower baseline concentrations
*because* they excreted more between exposure ending and first sampling, and
baseline was 41–75 months after cessation.

### 4.6 The between-person association does not survive adjustment

This project's longest-standing open question was whether Li 2022's
between-person result — higher initial level going with *longer* half-life,
the opposite sign to saturation — holds once age is controlled. Age is the
obvious confounder, because it drives both half-life and accumulated burden.

Individual data are not published, so no re-fit is possible. But the answer is
already in the supplement, in a table nobody quotes. **Table S14 reports partial
R² for each determinant of individual half-life — mutually adjusted by
construction**, since partial R² is the variance a term explains given the
others in the model.

```bash
python3 scripts/li2022_age_adjustment.py   # -> db/li2022_age_adjustment.csv
```

| determinant | PFOA | PFOS | PFHxS | range across 8 compounds |
|---|---|---|---|---|
| **Age** | 9.40% | 28.39% | 16.31% | **3.62 – 28.39%** |
| eGFR | 7.05% | 4.05% | 2.95% | 1.25 – 8.16% |
| Calprotectin | 6.78% | 5.43% | 7.71% | 0.24 – 7.71% |
| Gender | 0.87% | 3.02% | 5.39% | 0.87 – 6.16% |
| **initial PFAS tertile** | **0.70%** | **3.18%** | **2.32%** | **0.70 – 5.32%** |

**Initial PFAS level is the weakest determinant for 8 of 8 compounds**, beaten
by age by a median factor of 8.9× (range 2.0–37.4×). Adjusted, it explains under
1% of the variance in individual PFOA half-life.

A second, independent check points the same way. Comparing the two unadjusted
contrasts on the same cohort (n = 114), the tertile effect is only **6–37% the
magnitude of the age effect** for every compound — so a modest age imbalance
across tertiles accounts for all of it, and accumulated burden rising with age
makes that imbalance expected rather than hypothetical.

**Verdict: the between-person contradiction dissolves.** The within-person and
within-study evidence for mild concentration dependence stands alone and
unopposed, and the apparent between-person reversal was an age confound, exactly
as predicted. Note the sample for the adjusted table is n = 54, those with
complete data on all six determinants, against n = 114 for the unadjusted
tertile table.

> **Check your understanding.** Section 4.4 says higher dose → faster clearance
> (positive slope). Section 4.5 says the literature often shows higher exposure →
> *shorter* half-life, and calls that an artefact. Are those in conflict?
> *(No — they agree in sign. Faster clearance means shorter half-life. The point
> of 4.5 is that the observed between-cohort effect is far larger than +0.1 and
> is inflated by the background artefact, so it cannot be read as the mechanism.)*

---

## 5. Why published half-lives disagree: the assumption audit

`db/human_halflife_extended.csv` holds 142 rows from 39 studies, each carrying
the design fields that determine trustworthiness. `db/exposure_vs_halflife.csv`
(184 rows) pairs chemical, serum concentration, exposure amount and half-life
on one line, human and animal together.

### 5.1 The disagreement is about clearance, not Vd

```bash
python3 scripts/human_parameter_reconciliation.py
```

PFOA, each source's own triple with the missing term completed:

| source | Vd (mL/kg) | t½ (y) | CL (mL/kg/d) | derived |
|---|---|---|---|---|
| Abraham 2024 (measured, n=1) | **121** | 5.51 | 0.0417 | CL |
| Chiu 2022 (fitted) | 430 | 3.14 | 0.260 | CL |
| EPA 2024 (Vd assumed) | 170 | 2.70 | 0.120 | — |
| OEHHA 2024 (CL measured) | 398 | 2.70 | **0.280** | Vd |

Spread: Vd 3.6×, half-life 2.0×, **clearance 6.7×**.

Chiu's fitted parameters imply CL = 0.260 mL/kg-day against OEHHA's measured
0.280 — agreement to 7%. **So Chiu and OEHHA agree on clearance**, and the
~400 mL/kg Vd is what that shared clearance implies at a 2.7–3.1 y half-life.

Abraham et al. 2024 is the only **direct measurement** of human Vd: one
volunteer swallowed ~4 µg each of fifteen ¹³C-labelled PFAS and was followed
450 days in plasma, urine and faeces, giving `Vd = D_abs/C₀` with nothing
assumed. It gets 121 mL/kg and a 5.5 y PFOA half-life, hence a clearance **6.7×
below OEHHA's**. It needs no assumption about exposure cessation and no
background subtraction, because the label separates dose from background — it
is methodologically the strongest design in the human literature. It is also
n = 1, at tracer dose, with a half-life above every cohort estimate.

### 5.1a The donation trial adjudicates it — and sides with the measurement

Gasiorowski 2022 randomised 285 firefighters to donate plasma (up to 800 mL
every 6 weeks, mean 6.4 donations), donate whole blood (~470 mL every 12 weeks,
mean 4.3), or be observed. It removed a **known volume** of plasma from people
at **known serum concentrations** and measured the resulting fall. That is a
mass-balance experiment that fixes Vd without assuming a half-life:

```
Vd(L) = V_plasma removed × C_mean / ΔC_serum
```

The trial's own observation arm supplies the correction that makes it work —
subtracting it nets out continuing background intake and natural elimination.
**The authors never did this calculation.**

```bash
python3 scripts/gasiorowski_vd.py   # -> db/gasiorowski_vd.csv
```

| chemical | blood arm | plasma arm | Abraham 2024 | Chiu 2022 |
|---|---|---|---|---|
| PFOS | 113–124 | 182–199 | **152** | 320 |
| PFHxS | 71–78 | 171–187 | **125** | 290 |

Two arms of one trial, analysed independently, bracket Abraham's directly
measured values. **Chiu's fitted Vd sits above the entire range for both
compounds.**

Known directions of error: using the 800 mL ceiling as the plasma mean
*overstates* volume removed and so *overstates* the plasma-arm figure, which is
the higher of the two; counting only the plasma fraction of whole blood
*understates* the blood-arm figure. Body weight is reconstructed from the
reported BMI of 27.9 (weight is not reported) over a 84–92 kg range. The PFOA
row is not reliable — baseline was 1.2 ng/mL against a 1 ng/mL reporting limit,
and the authors' own sensitivity analysis more than halves the effect.

**So three independent human approaches — a labelled-dose study, a faecal/urinary
mass balance, and a randomised removal trial — all land at 74–199 mL/kg, while
the population-model fit gives 290–430.** The disagreement is real and it is now
three-to-one against the fitted value.

What it does *not* settle is clearance. For PFOS, Abraham's half-life (3.32 y)
and Chiu's (3.36 y) agree almost exactly while their Vd differ 2.1×, so the
clearance disagreement is entirely a Vd disagreement. For PFOA they differ on
both. **A population model that reproduces observed serum trajectories with a Vd
two to four times the measured one is fitting something, and what it is fitting
is the open question.**

### 5.1b EPA's own appendix says the human Vd literature never measured anything

The decisive confirmation is in a table nobody cites. EPA's 2024 appendix volume
tabulates the human PFOA Vd literature, and the table's **title** is the
finding:

> *Table B-26. Summary of PFOA Volume of Distribution Values **Assigned** in
> Human Studies*

```bash
python3 scripts/extract_epa_table_b26.py   # -> db/epa_table_b26_human_vd.csv
```

Ten entries. **All assigned, none measured.** Six assign exactly **170 mL/kg** —
Thompson 2010's calibrated value. The full range is 170–200 mL/kg, a spread of
1.18×.

| | Vd (mL/kg) | how obtained |
|---|---|---|
| Andersson 2025 | 74 | **measured**, urinary/faecal mass balance |
| Gasiorowski 2022 (derived, §5.1a) | 113–199 | **measured**, randomised removal |
| Abraham 2024 | 121 | **measured**, labelled oral dose |
| *EPA Table B-26, ten human studies* | *170–200* | ***assigned*** |
| Thompson 2010 | 170 | calibrated against an assumed half-life |
| Chiu 2022 | 430 | fitted by population model |

**The apparent tight agreement of the human Vd literature is an artefact of
everyone adopting the same assumed constant.** Every direct measurement falls
*below* the assigned range; the one population fit sits more than twice *above*
it. The assigned consensus sits in neither place, and nothing in that table
tests it.

Two details worth keeping: Shin 2011 is the only entry with a sex-specific Vd
(181 male, 198 female), and Gomis 2017's 200 mL/kg is explicitly an *average of
human and animal* values — a different circularity from Thompson's.

**One caveat on comparability, which cuts in the strengthening direction.**
These are not all the same kind of volume. Abraham reports a *terminal* volume
(`Vd = D_abs/C₀`, back-extrapolated from the terminal phase, so Vz or V_area),
while the EPA animal fits report **Vdss** and a **terminal beta-phase**
half-life — established by recovering Zurlinden 2025's supplementary tables,
where the variables are named `Vdss [pop]` and `halft_beta`. For a
two-compartment drug Vz ≥ Vss, so the steady-state volume implied by Abraham's
data is **at or below** 121 mL/kg. Treating his number as directly comparable to
the assigned 170–200 therefore *understates* the gap rather than inflating it.
Chiu's one-compartment fit has no such distinction, since Vz = Vss = V there.

### 5.1c The number ten studies assign, read at its source

```bash
python3 scripts/thompson_vd_circularity.py
```

Thompson et al. 2010 is where the 170 mL/kg comes from. Nine of the 42 adopted
regulatory clearance factors in §5.4 are computed from it, and Table B-26 above
shows ten human studies assigning it. The full text and its five supplementary
files are now on disk, and three things follow.

**It is not a measurement.** The supplementary table is headed, in the authors'
own words, *"Input data for the calibration of the Vd parameter from two US
communities"*, and its last column is headed **"calculated Vd"**. The paper's
Eq. 2c is `Vd = DP/(CP·kP)`. With the assigned `kP = 0.0008 /day`:

| community | dose | serum | Vd computed here | published |
|---|---|---|---|---|
| Little Hocking | 62 ng/kg-day | 448 ng/mL | 173.0 mL/kg | 173 |
| Lubeck | 9 ng/kg-day | 68 ng/mL | 165.4 mL/kg | 165 |

Both reproduce to three figures. The 170 mL/kg is an arithmetic consequence of
an assumed elimination rate, not an observation about distribution.

**The circularity is real, and worse than circular.** Inverting Eq. 2c, both
communities imply the same **2.37-year** half-life — and the paper says plainly
that it took that half-life from Bartell 2010 *because Bartell measured it in
the same two communities used for the calibration*. Now watch what happens when
an agency combines the adopted Vd with a half-life to get clearance:

```
CL = ln2·Vd/t½ = k·Vd = k · DP/(CP·k) = DP/CP
```

The half-life **cancels exactly**. Thompson's data support one clearance — the
intake-to-serum ratio, 0.132–0.138 mL/kg-day — and no other. But Vd does *not*
cancel; it is directly proportional to whatever half-life is assumed:

| assumed t½ | self-consistent Vd |
|---|---|
| 2.3 y (main text) | 168 mL/kg |
| 2.37 y (kP as used) | 173 mL/kg |
| 2.5 y (corrigendum) | 182 mL/kg |
| 3.3 y (Olsen 2007) | 241 mL/kg |
| 3.8 y (kP = 0.0005 as printed) | 277 mL/kg |

So adopting "170 mL/kg" while assuming any other half-life is not conservative
or neutral — it silently contradicts the exposure data the number was built
from. The paper's own **appended corrigendum** is this exact error: the authors
switched kP from 0.0005 to 0.0008 /day, applied the new Vd to the ng/day intakes
but not the ng/kg-day intakes, and had to correct the latter by a factor of
**0.6** — which is 0.0005/0.0008.

**The transposition is settled, and only one source has it backwards.** The
abstract and Table 1 give **PFOA = 170** and **PFOS = 230** mL/kg. This
repository already had them that way round. **Andersson 2025** is the source
that transposes them, stating Thompson yielded "230 mL/kg body weight for PFOA
and 170" for PFOS; Zhang 2013 assigns 170 to PFOA correctly, per Table B-26.

**And the PFOS value was never calibrated against anything.** The paper says
PFOS's 230 is "based on adjustment of the PFOA value" — scaled up by 20–50% on
the strength of a cited monkey model, landing on 170 × 1.35. No PFOS serum or
intake data entered it. §5.5 already noted the arithmetic; the full text
confirms there is nothing underneath it.

### 5.2 The assumptions that move the number, with directions

1. **Ongoing unmonitored exposure** → inflates half-life. The dominant bias.
2. **Background not subtracted** → inflates, scaling with proximity to
   background (quantified in §4.5).
3. **Isomers pooled** → not one-directional. Xu 2020 has every branched PFOS
   isomer faster than linear; Nilsson 2022b has branched slower; Rosato's pooled
   1m-PFOS is 5.86 y against L-PFOS 3.13 y.
4. **Which average** → AM vs GM is worth 10–15%; and Lyu 2026 shows
   `ln2/mean(k)` and `mean(ln2/k)` differ by 15–23% in the same data.
5. **Follow-up short relative to half-life** → unstable either way.
6. **Age and sex** → real modifiers. Lewis-Michl 2025 gives a clean 1.8-fold
   gradient within one cohort (1.96 y at 0–17 → 3.55 y at 60+). Paediatric
   estimates exist only from Batzella 2024 (1.64 y) and Lewis-Michl (1.96 y),
   both 30–45% below their adult strata, both cautioning that growth dilution is
   not separated from intrinsic elimination.
7. **Excretion routes counted** → renal-only underestimates clearance, hence
   overestimates half-life. Shi 2016 is the cleanest demonstration: F-53B total
   half-life 15.3 y against a renal-clearance-only 280 y **in the same paper**,
   an 18-fold gap.
8. **The assumed Vd** → scales mass-balance estimates linearly.

### 5.3 The faecal route is contested, not settled

| | Andersson 2025 (matched serum/urine/faeces, n=147) | Abraham 2024 (labelled bolus, n=1) |
|---|---|---|
| PFOS | faeces dominant, **4:1** | **faecal not detected** |
| PFOA | urine ~1.7:1 | faecal 17% of total |
| PFNA–PFDoDA | — | faecal 58–61% |

A candidate resolution neither paper tests: Andersson measured faecal PFOS in
people with ongoing dietary intake, so unabsorbed dietary PFOS passing through
would be scored as elimination — which a fully absorbed labelled bolus cannot
produce. If that is right, corrections to PFOS half-lives based on Andersson's
4:1 are not justified. Cholestyramine trials push the other way (Møller 2024,
n = 45: 63% PFOS lowering in 12 weeks vs 3% control).

**A controlled animal measurement now quantifies the split.** Argoul 2026
modelled urine and faeces simultaneously with plasma, over 119 days, in female
mice (Table 2):

| | renal | faecal |
|---|---|---|
| PFBA, PFHxA, GenX, PFO2OA | 97–100% | 0–2.4% |
| **PFOA** | **79.1%** | 20.9% |
| PFDS | 70% | 30% |
| **PFOS** | **46.3%** | 53.7% |
| **PFHxS** | **44%** | 56% |
| PFNA | 28.6% | 71.4% |
| PFDA | 7.3% | 92.7% |

**And the rat agrees, from 25 years earlier.** Kudo 2001 followed urine and
faeces for four PFCAs in Wistar rats over 120 h: urinary recovery falls from
**92% (PFHpA, C7) → 55% (PFOA, C8) → 2.0% (PFNA, C9) → 0.2% (PFDA, C10)** in
males, while faecal recovery rises with chain length. The paper states outright
that "feces must become a major route of elimination when PFCA was hardly
eliminated in urine as was observed with PFNA in male rats and PFDA in both
sexes." Two laboratories, two species, two decades apart, same monotonic
hand-off from renal to faecal as the chain lengthens. This is now one of the
better-replicated findings in the review.

In the rat the hand-off is also **sex-dependent**: female rats put 51% of a PFNA
dose in urine against the male's 2.0%, so C9 is a renal compound in females and
a faecal one in males. That is the same androgen-regulated step as §3.4, seen in
the route split rather than the half-life.

This does not adjudicate the human dispute — it is a different species, and the
mouse is the species §3.4 flags as having multi-route redundancy. But it shows
the route split is **strongly chain-length dependent within one species and one
experiment**, moving monotonically from 100% renal at C4 to 93% faecal at C10.
Any human accounting that assigns a single route fraction across compounds is
wrong in a predictable direction. It also means the Andersson-versus-Abraham
disagreement for PFOS, where our two human sources sit at 4:1 faecal and
faecal-not-detected, straddles a compound that in mice is genuinely close to
50:50 — so the dispute is about a compound where neither extreme is plausible.

**The dispute largely dissolves once the biliary loop is quantified.** EPA's 2016
PFOS Health Effects Support Document carries a number this review had never
seen: Harada 2007 sampled serum and bile from four gallstone-surgery patients and
measured a **biliary resorption rate of 0.97**. Bile PFOS (27.9 ng/mL) actually
*exceeds* serum (23.2), so bile is a genuine excretion route — but **97% of what
is secreted is reabsorbed from the gut.**

That is the enterohepatic twin of §3.3's renal axis, and it reframes the whole
section. Humans run **two near-complete reabsorption loops**, not one:

| loop | fraction reabsorbed | source |
|---|---|---|
| renal (PFOA) | **99.94%** | Han 2012 Table 4 |
| biliary (PFOS) | **97%** | Harada 2007 |

Both are near-unity, both are the reason the human half-life is long, and
crucially **both are interruptible** — which is why bile-acid sequestrants work:

| | PFOS half-life | PFHxS |
|---|---|---|
| **Delaere 2025**, treated (cholestyramine and/or plasma donation, n=19) | **1.2 y** | **2.5 y** |
| **Delaere 2025**, observation (n=9) | **7.3 y** | **9.4 y** |
| ratio | **6.1×** | **3.8×** |

plus Genuis 2010 (4 g/day cholestyramine, 20 weeks: PFOS 23 → 14.4 ng/g) and
Møller 2024 (63% lowering in 12 weeks against 3% control). A sequestrant that
blocks a 97% resorption loop should accelerate elimination several-fold, and it
does.

This also reconciles Andersson and Abraham without either being wrong. With 97%
resorption, **gross** biliary flux is large while **net** faecal elimination is
small. Andersson measured faeces in people with ongoing intake and saw a large
signal; Abraham followed a labelled bolus and saw almost nothing leave by that
route. Those are measurements of different quantities, and the 0.97 reconciles
them.

*Two caveats on Delaere.* Its own abstract states that "the study did not conduct
statistical comparisons" and — more seriously — that "the calculations only
included data from participants whose serum PFOS and PFHxS concentrations
decreased." Conditioning on a decrease biases apparent half-lives downward in
both arms. The direction of the treatment effect is credible; the absolute
half-lives are not.

**And two further routes appear in no agency's clearance accounting at all.**

*Breastfeeding.* Mondal 2014 (C8 Science Panel, n = 633) finds each month of
breastfeeding lowers maternal serum by **3% for PFOA, 3% PFOS, 2% PFNA, 1%
PFHxS** — and raises the infant's by 6% and 4%. For a woman breastfeeding a
year, that is a third of her PFOA burden leaving by a route no clearance factor
counts. The mother's excretion is the infant's dose.

*Menstrual blood loss.* This now has a number attached. Wong 2014 fitted a
population PK model to six NHANES cycles and concluded menstruation accounts for
**about 30%** of the male/female PFOS half-life difference; Verner & Longnecker
2015 revised the annual serum loss from 432 to 868 mL/year, implying **>30%**
(both via EPA's 2016 PFOS document). Upson 2022 reviews the epidemiology: postmenopausal
women carry higher PFAS than premenopausal women, concentrations rise with years
since menopause and with hysterectomy, and a life-stage PBPK bias analysis
reproduces the apparent PFAS–menopause association purely from the loss of this
excretion route.

This matters for §3's framing. The rat sex difference is transporter-mediated
and androgen-driven; the **human** sex difference may be mostly mechanical —
blood leaving the body. Wallis 2023 states the position plainly: sex differences
from hormone-mediated transporters "have not been directly observed in humans."
Two species, two different reasons for the same-looking pattern.

### 5.4 Agencies disagree about clearance more than about half-life

```bash
python3 scripts/agency_clearance_comparison.py
```

| | adopted human clearance factor | spread |
|---|---|---|
| PFOA | 0.059 → 0.280 mL/kg-day | **4.7×** |
| PFOS | 0.016 → 0.390 mL/kg-day | **24×** |

EPA and OEHHA cite the **same** Li et al. half-life studies and land 2.3× apart,
because EPA computes `CL = Vd·ln2/t½` from Thompson's assumed 170 mL/kg while
OEHHA regresses intake against serum directly and states this "obviat[es] the
need for a T½ estimate".

Thompson et al. 2010 is worth knowing about: it is an Australian
exposure-reconstruction paper, and its Vd is *calibrated*, not measured —
`Vd = Dose·t½/(Css·ln2)` using Emmett 2006 serum, an **assumed** 2.3-y
half-life and fa = 0.91. And PFOS 230 mL/kg is literally 170 × 1.35, the 1.35
taken from Andersen 2006's monkey *central-compartment* volumes. EPA then uses
that Vd to compute clearance — circular. The measured value is 121.

Neither **OEHHA 2024 nor EPA 2024 cites Chiu 2022 at all** (verified by grep
over both complete documents, 38,696 and 13,747 lines).

### 5.5 The states agree with each other by repeating one derivation

Adding thirteen US state and further international bodies (261 rows,
`db/state_international_regulatory.csv`) widens the picture rather than
settling it.

```bash
python3 scripts/all_agency_clearance.py   # -> db/all_agency_clearance.csv
```

For **PFOA**, five jurisdictions — New Jersey DWQI, Michigan MPART,
Pennsylvania DEP, the Drexel advisory group that supports it, and New Hampshire
DES — all adopt **0.140 mL/kg-day**. That looks like consensus. It is one
derivation performed five times: every one computes `CL = Vd·ln2/t½` from
**Thompson 2010's Vd of 170 mL/kg** and **Bartell 2010's 2.3 y half-life**. The
same pattern repeats for PFOS at 0.128–0.130, where the Vd is Thompson's 230.

Across all 42 adopted values:

| how the number was produced | n |
|---|---|
| renal clearance only | 11 |
| **computed from Thompson 2010's assumed Vd** | **9** |
| animal-derived | 7 |
| adopted wholesale from another agency | 7 |
| measured: intake vs serum at steady state | 4 |
| other or unstated | 4 |

Recall from §5.4 what Thompson's Vd actually is: *calibrated* from two US water
communities using an **assumed** 2.3-year half-life, not measured. So the
apparent agreement among states is not independent confirmation — it is one
assumption propagating, and §5.1a shows the measured Vd is 121 mL/kg, not 170.

The jurisdictions that do **not** use it diverge sharply, and in both
directions: OEHHA's measured intake-vs-serum clearance is ~2× the
Thompson-derived cluster, while its renal-clearance-only figure is 2–8× below
it. **Minnesota MDH is the one state that broke ranks**, adopting OEHHA's
0.280/0.390 wholesale and landing twice as far from its neighbours as any
methodological disagreement between them.

Two further observations worth recording:

- **New York derived its PFOA clearance (0.092) from Macon et al. 2011, a mouse
  study.** Using animal data to set a human clearance factor is a different kind
  of assumption from the others here.
- **Four states set their human PFHxS clearance (0.086–0.090) from Sundström
  2012's cynomolgus monkey Vd.** The human PFHxS factor in Michigan, Minnesota,
  New Hampshire and New York rests on a monkey volume of distribution.
- The extraction also caught an apparent error: **ITRC's PFOS entry uses
  Bartell 2010's 840-day half-life, which is the PFOA value.** Recorded as
  found, not corrected.

---

## 6. In vitro binding to albumin and related proteins

`db/protein_binding.csv` — 273 rows, 54 studies, 30 protein targets, 272 with a
DOI, computational values flagged.

```bash
python3 scripts/build_binding_table.py
# -> db/binding_summary_ka.csv (association constants)
# -> db/binding_summary_logD.csv (protein/water distribution)
```

### 6.1 A single affinity number for a PFAS–protein pair does not exist

| chemical | protein | n | Ka range (M⁻¹) | spread |
|---|---|---|---|---|
| PFBS | HSA | 4 | 116 – 6.5 × 10⁶ | 55,800× |
| PFBA | HSA | 3 | 20.4 – 1 × 10⁶ | 49,100× |
| PFOA | HSA | 6 | 103 – 1.98 × 10⁶ | 19,200× |
| PFOA | BSA | 6 | 322 – 1 × 10⁶ | 3,110× |
| PFOS | BSA | 4 | 1.16 × 10³ – 2.44 × 10⁵ | 210× |

**The spread is not fluorescence-versus-dialysis.** Equilibrium dialysis alone
spans 3.2 × 10² to ~10⁶, and fluorescence quenching spans 10² to 2 × 10⁶ by
itself. The drivers are the **ligand:protein ratio regime** and the **fitting
model**: Alesio & Bothun 2022 obtained Ka from 0.07 to 6.16 × 10⁶ M⁻¹ from a
*single dataset* using three models. Differential scanning fluorimetry is
systematically low and its own authors say so, because Kd comes off a thermal
ramp rather than 37 °C.

Reported stoichiometry is a titration curve read at five points, not a
disagreement: n ≈ 1 (quenching), 1–3 (low-ratio dialysis), 6–9 (microdesalting),
8 (nanoESI-MS at 4:1), 45 (dialysis at 1.2 mM). **Physiological PFAS:albumin is
≤ 0.0005**; most in vitro work runs three to six orders of magnitude above it.

*Any table printing one number per pair is misleading. Use the range.*

### 6.2 The one dataset measured at a realistic ratio

Fischer et al. 2024 (solid-phase microextraction, PFAS:protein ≤ 0.004), log
D_protein/water in L_water/kg:

| chemical | HSA | γ-globulin |
|---|---|---|
| PFHxA | 3.20 | 1.73 |
| PFHpA | 4.10 | 1.81 |
| PFOA | 4.48 | 2.16 |
| PFNA | 4.49 | 2.81 |
| PFDA | 4.73 | 3.67 |
| PFHxS | 4.77 | 1.95 |
| PFOS | 4.62 | 3.33 |

Albumin dominates below ~7 perfluorinated carbons (D_glob/D_HSA < 0.1 for
MW < 500) and globulins catch up as the chain lengthens, reaching parity around
ηpfc 12–13.

Measured serum unbound fractions: **PFOA 0.061%, PFOS 0.042%, PFHxS 0.041%,
PFHpS 0.035%** — with interindividual variation of only 1.6–2.5×.

**The apparent conflict with Han 2003 is resolved, and it was never a
contradiction between two measurements.**

```bash
python3 scripts/binding_conflict_han_vs_fischer.py
```

Three things fall out of reading Han 2003 directly:

- **">90% bound" is a calculation, not a measurement.** Han's own abstract says
  so: "*On the basis of these binding parameters and the estimated plasma
  concentration of serum albumin, greater than 90% of PFOA would be bound*". The
  measured quantities are Kd = **0.3–0.4 mM** and **n = 6–9 sites** per albumin.
- **It is a floor, and Fischer's value satisfies it.** f_unbound = 0.00061 is
  99.94% bound, which is indeed >90%. The error entered downstream, when PBPK
  models read ">90%" as "≈90%" and used f_unbound ≈ 0.1 — about **164× too
  high**. That is the real defect, and it is in the inheritance, not in Han.
- **The underlying disagreement is in Kd, and it is a ratio artefact.** Han
  titrated **50–60 µM albumin with 0.1–3 mM PFOA** — a ligand:protein ratio of
  **1.7:1 to 60:1**. Fischer worked at **≤0.004:1**. Human serum sits at
  1.6 × 10⁻⁵ to 4 × 10⁻³ : 1, so Han's lowest ratio is **414× above even an
  occupational serum**. At 1.7–60:1 the high-affinity site is saturated
  instantly and the fitted constant is dominated by the 6–9 low-affinity sites;
  reproducing Fischer's f_unbound would need Kd ≈ 3.3 µM, about **91–121×**
  below what Han measured. The two studies measured different things, and only
  one measured the regime humans are in.

A third value sits between them and agrees with neither extreme: **Ohmori 2003**
reports rat plasma protein binding "over 98% for all PFCAs tested"
(f_unbound < 0.02).

Han 2003 also records that **ultrafiltration failed outright** for this assay —
"PFOA completely nonspecifically bound to the membrane" — which is why
microdesalting columns were used, and is a standing warning about method choice
in PFAS binding work.

**A species-matched set now exists for the mouse.** Argoul 2026 Table 4 reports
equilibrium-dialysis unbound fractions in CD-1 mouse plasma for ten PFAS — the
first set in this review measured in the same animals whose clearance is
reported alongside:

| chemical | mouse fu | chemical | mouse fu |
|---|---|---|---|
| PFOS | 0.25% | PFO2OA | 10% |
| PFNA | 0.35% | PFHxA | 25% |
| PFDA | 0.42% | GenX | 26% |
| PFOA | **0.87%** | PFBS | 29% |
| PFHxS | 1.3% | PFBA | 77% |

Mouse PFOA at 0.87% sits **14× above** Fischer's human 0.061% and **11× below**
Han's human ceiling of 10%. Two things follow. First, the Han-versus-Fischer gap
is not a species artefact — the mouse value lands between them, so it cannot be
invoked to rescue either. Second, a PBPK model that borrows a human f_unbound
for a mouse compartment, or the reverse, is off by about an order of magnitude
before anything else goes wrong; the ordering across compounds is nonetheless
the same in both species, with the long-chain sulfonates most tightly bound.

### 6.3 Binding does not explain the species difference

Two results settle this, and both are negative:

- **Starnes 2024** (Table 5): rat albumin binds PFOA ~2× *more* tightly than
  human albumin, while the rat half-life is ~1000× *shorter*. Wrong direction,
  wrong magnitude.
- **Han et al. 2004** explicitly reject α2u-globulin (Kd ~10⁻³ M) as the cause
  of male-rat PFOA retention — the authors of the α2u-globulin paper in your
  own citation list.

Also worth knowing: PFAS bind albumin at only **2–10% the affinity of
similar-chain fatty acids** (Hebert & MacManus-Spencer 2010), so nearly all in
vitro work uses fatty-acid-free albumin while in vivo binding competes against a
higher-affinity endogenous excess.

PFOA's albumin site is genuinely contested — crystallography puts the
high-affinity site in subdomain IIIA (Maso 2021) while warfarin displacement
puts it at Sudlow site I (Peng 2024, 96.1% reduction; Chen 2015). PFOS is
consistent across methods (FA3/4 + FA5, PDB 4E99).

---

### 6.3a What the in vitro literature says about species and structure

`db/protein_binding.csv` holds 273 rows from 54 studies, but most of it cannot
answer either question, because §6.1 showed a single affinity number per pair
does not exist. What can answer them is the handful of studies that hold method
constant and vary one thing. Run
`python3 scripts/binding_species_and_structure.py`.

**Species: albumin affinity is close to species-invariant.** Starnes 2024
measured four albumins by one method (differential scanning fluorimetry), Kd in
mM, lower = tighter:

| compound | nC | human | rat | bovine | porcine | rat/human | max/min |
|---|---|---|---|---|---|---|---|
| PFBA | 3 | 2.64 | 1.90 | 3.06 | 1.28 | 0.72× | 2.4× |
| PFBS | 4 | 1.68 | 1.41 | 0.80 | 1.25 | 0.84× | 2.1× |
| GenX | 4 | 1.57 | 1.64 | 2.55 | 2.12 | 1.04× | 1.6× |
| PFHxA | 5 | 1.71 | 1.33 | 2.17 | 0.95 | 0.78× | 2.3× |
| 6:2 FTSA | 6 | 0.41 | 0.33 | 1.39 | 0.52 | 0.80× | 4.2× |
| PFOA | 7 | 0.79 | 0.40 | 0.99 | 0.37 | **0.51×** | 2.7× |
| PFOS | 8 | 0.69 | 0.42 | 0.86 | 0.38 | 0.61× | 2.3× |

The spread across four species never exceeds ~4× and is usually under 2×.
Albumin is a conserved protein and behaves like one. Rat binds *tighter* than
human throughout (geometric mean 0.74×, PFOA 0.51×) while eliminating PFOA
~1000× faster — the wrong direction, as §6.3 says. Fischer 2024 agrees
independently at a physiological ligand:protein ratio: bovine-minus-human log D
is within 0.3 log units for most of 27 compounds, median −0.14. **BSA is a fair
stand-in for HSA**, which matters because much of the literature uses BSA.

**Structure: three findings, and only the first is the textbook one.** Fischer
2024, log D_protein/water, one method, 10 chain lengths:

1. *Chain length is not monotonic.* Carboxylates rise from 2.16 (nC4) to a
   maximum 5.18 around nC10–11, then **fall** to 3.87 at nC13. Slope +0.45 per
   CF₂ up to nC10, −0.44 after. A model linear in carbon number is wrong at both
   ends.
2. *The head group is worth about a log unit.* At matched chain length the
   sulfonate binds harder — PFHxS 4.77 vs PFHpA 4.10 at nC6; PFHpS 5.49 vs PFOA
   4.48 at nC7. Slopes are similar (+0.75 vs +0.79 per CF₂); the series are
   offset.
3. *Albumin stops being the main carrier.* γ-globulin binds 2–3 log units weaker
   at short chains and catches up steadily; for PFTeDA (nC13) globulin binds
   **more** tightly than albumin. Treating "protein binding" as "albumin binding"
   gets worse as the chain lengthens — and the long chains are the persistent
   ones.

**The gap.** There is **no in vitro PFAS–mouse serum albumin binding study**, in
this database or in PubMed: a search returns one record, which measures mouse
albumin *concentration*, not binding. So the species question can be answered
for human, rat, bovine and porcine, and not for the mouse — the animal §10's
endpoint is built from. See `WANTED.md`.

### 6.4 How much of the species difference is binding? The arithmetic

§6.2 settles what the unbound fraction is and §11.4 shows what the choice does
to the transport endpoint. The arithmetic joining them is worth stating plainly,
because it is what "does binding affect elimination?" actually resolves to.

Only unbound PFAS is filtered, so the clearance free filtration alone would
produce is fu·GFR, and the measured renal gap factorises:

```
CL_renal(mouse)     fu(mouse)     GFR(mouse)      R(mouse)
---------------  =  ---------  ×  ----------  ×  ---------
CL_renal(human)     fu(human)     GFR(human)      R(human)
   120×              BINDING       BODY SIZE      TRANSPORT
```

PFOA, female mouse (CL_renal 3.59 mL/d/kg, fu 0.87%) against human (CL_renal
0.03 mL/d/kg, Han 2012 Table 4). The 120× gap is fixed; what the fu assumption
changes is how it is **apportioned**:

| human fu | binding | body size | transport | dominant term |
|---|---|---|---|---|
| 0.00061 (measured) | **14.3×** | 3.9× | 2.2× | **binding** |
| 0.02 (assumed) | 0.4× | 3.9× | 70.7× | transport |
| 0.10 (inherited) | 0.1× | 3.9× | 353.5× | transport |

Under the measured value **binding is the largest single term**: humans hold
PFOA ~14× more tightly than mice, which alone keeps 14× less of it in front of
the glomerulus. Under the assumed values the apportionment inverts — binding
runs the *wrong way* (humans would bind less tightly than mice) and transport
must carry 71–354× to balance. That decides whether the species difference is a
binding problem or a transporter problem.

Two cautions. The mouse and human fu come from different methods (equilibrium
dialysis vs solid-phase microextraction) at different ligand:protein ratios —
precisely the artefact §6.1 diagnoses — so the 14× is the least secure number
here. And the mouse limb inherits the cocktail study's ~2.5× low clearance bias
(§11.3a), which lands on the transport term: correcting it raises 2.2× to ~5×,
still far below 71×.

**Binding alone does not explain long human half-lives.** At fu = 0.00061, free
filtration would clear 1.57 mL/d/kg, which with Vd = 170 mL/kg is a half-life of
**75 days** — weeks, not years. Measured human renal clearance is 0.03 mL/d/kg,
**52× lower still**, because tubular reabsorption recovers **98.1%** of what is
filtered. Binding sets how much is presented to the glomerulus; reabsorption
disposes of almost all of that. Binding is the species-variable factor,
reabsorption the near-complete one in every species measured — a model that gets
binding right and reabsorption wrong still misses the half-life by two orders of
magnitude.

Reproduce with `python3 scripts/binding_elimination_decomposition.py`.

---

## 7. Transporters: direction first, then kinetics

`db/transporter_kinetics.csv` — 217 rows, 27 studies, 48 transporter entities.

**Direction decides the sign of every inference:**

- **Reabsorptive (lengthen half-life):** rat Oatp1a1 (apical), human OAT4,
  human URAT1, ASBT; NTCP/OSTα-β on the enterohepatic loop. Weaver 2010:
  *"Oatp1a1 is the major player in the reabsorption of PFCAs."*
- **Secretory (shorten half-life):** OAT1, OAT3 (basolateral uptake into the
  tubule cell), and the apical pumps P-gp and BCRP.

Representative measured Km, from primary tables:

| transporter | species | PFAS | Km (µM) | source |
|---|---|---|---|---|
| Oatp1a1 | rat | PFOA / PFNA / PFDA | 126.4 / 20.5 / 28.5 | Weaver 2010 Tables 1–2 |
| Oatp1a1 | rat | PFOA | 162.2 | Yang 2009 |
| OAT1 / OAT3 | human | PFOA | 48.0 / 49.1 | Nakagawa 2008 |
| OAT1 / OAT3 | rat | PFOA | 51.0 / 80.2 | Nakagawa 2008 |
| OAT4 | human | PFOA | 172.3 (pH 6) / 310.3 (pH 7.4) | Yang 2010 |
| URAT1 | human | PFOA | 64.1 | Yang 2010 |
| OAT4 | human | 6 PFAS | 39–92 | Louisse 2023 Table 1 |

Three conflicts left open rather than silently resolved: Louisse 2023 found
**URAT1 transported nothing**, contradicting Yang 2010; Argoul 2026 reports
PFDA/PFOS as substrates of nothing, against both Louisse and Weaver; and
Cheng 2006's abstract calls renal Oatp1a1 "female-predominant", contradicting
Cheng 2005 and its own androgen conclusion.

**Two numbers that must not be cited as transporter measurements.** Worley &
Fisher's rat PBPK required a *fitted* male/female apical activity ratio of
~25,800×, and Loccisano (PFOA) and Kim (PFHxS) both report Tm/KT sex ratios of
exactly 91×. These are scalars absorbing an observed clearance difference, not
measurements of transporter expression. Relatedly, EPA's own two documents give
the rat male/female Oatp1a1 mRNA ratio as 2.5-fold and as 5–20-fold, and neither
traces to Kudo 2002's abstract — **no single fold-value should be quoted without
the primary full text.**

---

## 8. What this review corrects in the existing repository

1. **`SYNTHESIS.md`'s reconciliation of Zhang 2013 does not survive.** It
   rescales Zhang's clearance with Chiu's *fitted* Vd of 430 mL/kg, gets 3.29 y,
   and concludes Zhang and Chiu "agree to within 5%". Abraham's **measured**
   121 mL/kg gives ≈0.93 y instead. The reconciliation was an artefact of
   choosing the largest published Vd.
2. **Zhang 2013's headline PFOA numbers are 2.1 y (young females) and 2.6 y
   (males/older females), arithmetic means** — not the 1.3 y GM currently in
   `studies.csv`. Zhang's own Discussion says the value "is similar to the actual
   PFOA serum elimination half-lives reported for highly exposed populations
   (i.e., 2.3 or 2.9 years)". The short-half-life argument built on 1.3 y is
   weaker than the existing files imply.
3. **~~Zhang 2013 and~~ Andersson 2025 transposes Thompson 2010's Vd.**
   *Resolved* (§5.1c). Thompson gives PFOA = 170 and PFOS = 230 mL/kg. This
   repository already had them the right way round. **Andersson 2025** is the
   single source with them reversed ("230 mL/kg body weight for PFOA and 170");
   Zhang 2013 assigns 170 to PFOA correctly, so naming it here was wrong.
4. **The `APPRAISAL.md` framing of Vd** treats Chiu's 0.43 L/kg as the
   comparator against "commonly assumed" 0.17–0.20. Abraham's measurement
   supports the assumed values, not the fitted one.
5. **An earlier claim in this session was overstated.** I reported that
   OEHHA's implied Vd (398 mL/kg) and Chiu's fitted Vd (430) were a convergence
   of two independent routes. They are not independent: Chiu's clearance and
   OEHHA's clearance agree, so at a common half-life they necessarily imply a
   common Vd. The real finding is agreement on *clearance*, and Abraham
   disagrees with both by 6.7×.
6. **`species_exposure.csv` mixes time units silently** — Chiu's human
   half-lives in years, EPA animal fits in days, no units column.
   `db/master_exposure_halflife.csv` carries an explicit `native_time_unit`.
7. **EPA's PFHxA human Vd is wrong by 6.1×** — it assumed human Vd = monkey Vd
   = 730 mL/kg; the measured human value is 119. EPA's *reasoning* (Vd roughly
   conserved across mammals) is vindicated; its number is not.
8. **EPA's rat Oatp1a1 fold-value is a different transporter's number.** The
   PFHxA review's "2.5-fold" male/female renal Oatp1a1 mRNA ratio is Kudo 2002's
   **OAT-K** ratio. The Oatp1a1 figure in that paper is **23×** (§3.4). EPA's
   PFOA assessment's 5–20-fold understates it.
9. **This review's own CPHEA refits are biased, and the direction is now
   known.** `db/cphea_fitted_halflives.csv` holds 186 terminal log-linear fits
   computed here, selecting the window with the best adjusted R². On a biphasic
   curve that rule picks the shallow terminal tail, not the dominant elimination
   phase a published half-life describes. Kudo 2002 Table 2 makes this
   measurable, because CPHEA study 2990271 *is* that experiment, dose and all:

   | | fitted here | Kudo 2002 published | ratio |
   |---|---|---|---|
   | female rat | 1.44 d | 0.08 d | **17.9×** |
   | male rat | 8.08 d | 5.68 d | 1.4× |
   | M/F ratio | 5.6× | **71×** | understated 13× |

   Both sexes inflate, the fast-eliminating one far more, so the bias does not
   cancel in a ratio — it **compresses ratios toward 1**. Every conclusion in §3
   argues these ratios are large, so all of them survive and several
   strengthen. But the fitted column must not be read as comparable to published
   half-lives. It now carries a `tail_selection_flag`, and
   `scripts/validate_cphea_fits.py` reproduces the check. Argoul 2026 makes the
   same methodological point from the other direction, preferring mean residence
   time precisely because "the terminal half-life does not reliably reflect the
   overall persistence" when more than one compartment is present.
10. **One gap this review listed as open was not open.** §3.4 stated that no
   sex-resolved mouse PFOA half-life pair had been measured anywhere. Lou 2009
   Table 2 reports it — female 15.6 d, male 21.7 d — and the "15.6–21.7 d range"
   EPA quotes from that paper is the two sexes, not a dose range. Tatum-Gibbs
   2011 cites it the same way, independently.
11. **The fu = 0.02 assumption behind the reabsorption axis is Han 2012's, not
   OEHHA's.** §3.4 attributed it to OEHHA. It is stated in the footnotes to Han's
   own Table 4; OEHHA inherited it. Han also calls those values "rough
   estimates", a caveat the citation chain dropped.
12. **OEHHA's adaptation of Han's Table 4 altered two values** — human 99.94% →
   99.8% and male rat 93.7% → 93.2%. Neither changes a conclusion. §3.3 now
   quotes the primary figures.
13. **§3.3's "single axis" is two factors, and the report said one.** Han 2012
   points out that humans reabsorb *less* PFOA in absolute terms (51 mL/d/kg)
   than male rats (270) or mice (318–324); what is extreme is the fraction. Of
   the 333× male-mouse-to-human renal clearance gap, the escape fraction carries
   50× and the ~6× lower human GFR carries 6.5× — 68% and 32% on a log scale.
   The axis framing survives because reabsorption is the larger term, but a third
   of the species difference is filtration rate, not transport.
14. **§4.3's conclusion was stated more strongly than its basis supports.**
   Saturation is unreachable on every in vitro transporter measurement, but the
   PBPK constants agencies' models run on (Han 2012 Table 7, human KT = 0.133 µM)
   are 483–2,336× lower and imply saturation *is* reached in contaminated
   communities. §4.3 now gives both and argues for the in vitro side rather than
   assuming it.
15. **The want list had the wrong identifier for the dog study, for the whole
   project.** `WANTED.md` paired "Hanhijärvi et al. 1988, beagle dog" with
   doi:10.1007/978-3-642-71248-7_96. That DOI is **Kojo, Hanhijärvi, Ylinen &
   Kosma 1986, "Toxicity and Kinetics of Perfluoro-octanoic Acid in the Wistar
   Rat"** — a 28-day rat study with no dog data in it at all (Hanhijärvi is a
   co-author of both, which is how the two were conflated). Repeated retrieval
   attempts were therefore chasing the wrong paper. The dog study is a separate
   publication, indexed as EPA HERO 5412773.
16. **A PFOS calculation inherited the PFOA volume.** EPA's 2016 PFOS document
   records that Zhang's Tianjin intake estimates used "a volume of distribution
   of **170 mL/kg** (Thompson et al. 2010; Egeghy and Lorber 2011)" — for
   **PFOS**, whose value in Thompson is 230. Whether the 170 came from
   misreading Thompson or from Egeghy & Lorber, a PFOS calculation ran on the
   PFOA number, a 1.35× error in the conservative direction.
17. **"PFAS are metabolically inert" is not universal.** §1 and most of the
   literature treat PFAS as non-metabolised. Yi 2022 shows 6:2 Cl-PFESA
   undergoing reductive dechlorination to 6:2 H-PFESA (13.6% in rat liver,
   reductive conditions only) and notes it is the **second** perfluoroalkyl acid
   reported to biotransform in mammals. The premise holds for the legacy
   compounds this review centres on; it does not hold for the chlorinated ethers.
18. **The human reabsorptive transporter is not the rat's.** §3.4 implied the
   Oatp1a1 mechanism extends to humans. Yang 2010 shows **OATP1A2, the closest
   human orthologue, does not transport PFOA at all**; human apical reabsorption
   runs through OAT4 and URAT1, neither androgen-regulated. The reabsorbed
   *fraction* transfers across species; the *mechanism* does not.

---

## 9. What nobody knows

### 9.0 The shape of the hole, measured

Before the specific gaps, the general one. Merging every extraction in `db/`
onto a single grid of chemical × species × parameter, and counting *distinct
studies* per cell:

```bash
python3 scripts/build_coverage_matrix.py   # -> db/coverage_matrix.csv
```

| | cells | share |
|---|---|---|
| two or more studies | 156 | 9% |
| **a single study only** | **94** | **5%** |
| **no data at all** | **1,466** | **85%** |

39 chemicals × 11 species × 4 parameters = 1,716 cells. **Five sixths are
empty.** Even PFOS, the best-covered compound in the world, fills 29 of its 44
cells; PFOA 27.

By species, the collapse is steep: human 71 cells with any data, rat 61, mouse
45, monkey 28 — then pig 13, cattle 9, bird 8, fish 6, sheep/goat 4, rabbit 3,
**dog 2**.

### 9.0a The gaps are real, not an artefact of searching

That emptiness could be an indictment of this review's retrieval rather than of
the literature. It was worth testing, and the test is the strongest single
result about the shape of the field.

Four regulatory compilations were mined in full — EPA's 2024 Appendix B, NJ
DWQI's three MCL support documents, the remaining state and national documents,
and a final targeted sweep — yielding **806 new extracted rows citing 332
distinct primary studies**, of which roughly **125 had a first author not
present anywhere else in the collection**.

Coverage moved from 143 cells with two or more studies to **156**. Twelve cells.

**The compilations deepen provenance; they do not widen coverage.** Hundreds of
studies, concentrated on the same handful of chemical × species combinations
that were already well covered. Some of what they add exists *only* through
them — unpublished Clewell 2006 analyses, Clark data used by EPA 2021, NHANES
survey rounds, and Hanhijärvi 1988's dog work, a book chapter not indexed in
PubMed.

So the empty five sixths of the grid is a property of the field, not of the
search. **Nobody has measured these things.**

Coverage collapses away from four compounds and two species. PFOA, PFOS, PFHxS,
PFBS, PFBA and PFHxA account for most of what exists. The replacement
chemistries that are displacing them — HFPO-DA aside — are nearly empty:
6:2 Cl-PFESA has 2 of 44 cells, cC6O4 2, ADONA 3, and the diPAPs 1 each.

Two late additions moved cells without changing that picture. **EFSA 2020
Appendix C** (`scripts/parse_efsa_appendix_c.py`, 109 records) carries mouse
clearance and Vd for PFUnDA, PFDoDA, PFTrDA and PFTeDA from Fujii 2015 that
appear nowhere else, and deepened twelve cells from single-study to corroborated
— but opened none that were empty. Three of its rows were **refused** rather
than parsed: the web rendering lost `<br>` separators, so the PFOS row offers 18
half-lives against 20 doses, with one cell merging three values and another
repeating one. Aligning those would have meant guessing which half-life belongs
to which dose.

And **Wallis 2023** supplies what are, in its own words, "the first and possibly
only estimates of human elimination half-lives" for three fluoroethers, from the
GenX Exposure Study after Cape Fear discharges were controlled:

| compound | human half-life | 95% CI |
|---|---|---|
| PFO4DA | 127 d | 86–243 |
| Nafion byproduct 2 | 296 d | 176–924 |
| PFO5DoA | 379 d | 199–3,870 |

These sit between the short-chain carboxylates (PFHxA ~32 d, PFHpA ~62 d) and
the long-chain PFAAs (800–1,200 d) — consistent with ether oxygens shortening
half-life, which is the design rationale for the replacement chemistries. They
are the only human numbers testing it, from 44 people and two blood draws, and
the widest interval spans twentyfold. (Note a transcription error in the paper
itself: its Results text repeats Nafion byproduct 2's confidence interval for
PFO4DA; the abstract's 86–243 is the consistent one.)

This is the context for every ratio in this report. **When a cell holds one
study, a headline number is that study**, which is how the 31× mouse/rat PFOA
ratio came to rest on Lou 2009 alone (§3.2).

### 9.1 Specific gaps

Ranked by how much a single experiment would settle.

1. **No sex-resolved mouse PFOA half-life pair has ever been measured.** The
   31× mouse/rat ratio that frames this whole question depends on it. This is
   the single cheapest decisive experiment in the review.
2. **No study has measured renal Oatp1a1 protein in rat and mouse side by
   side**, and no Slco1a1-null animal has been dosed with PFAS. Without these
   the rat-vs-mouse mechanism cannot be closed. **No mouse transporter kinetics
   exist for any PFAS**, with any transporter.
3. **Abraham 2024 at a second dose level in the same volunteer** would settle
   human PFOA concentration-dependence outright, and would resolve whether the
   6.7× clearance disagreement is a dose effect or an individual.
4. ~~No controlled-removal study reports a Vd~~ — **done, §5.1a.** The
   Gasiorowski 2022 reanalysis gives PFOS 113–199 and PFHxS 71–187 mL/kg,
   bracketing Abraham's measured values and excluding Chiu's fitted ones. What
   remains open is *why* a population model needs a Vd two to four times the
   measured one to reproduce observed serum trajectories.
5. ~~Li 2022's tertile analysis re-run with age adjustment~~ — **done, §4.6.**
   The adjusted answer was already published as partial R² in Table S14:
   initial level is the weakest of six determinants for all eight compounds.
6. **Rat-vs-mouse tissue partition coefficients for PFOA/PFOS** are missing, as
   is any rat-vs-mouse biliary comparison under a common protocol. MRP4 appears
   never to have been tested for PFAS.
7. **Huang 2021's corrigendum content is unobtainable**, so which NTP numbers
   were corrected is unknown; affected rows are flagged.

### Three standing caveats on the animal data

- **"Mouse" means CD-1.** Every mouse PK dataset retrieved is CD-1, and Buist
  2004 found the *direction* of renal Oat3 sex bias differs between 129J and
  C57BL/6. Mouse renal transporter sex bias is not a species constant.
- **Transporter abundance is not a fixed parameter in these studies.** PFAS
  downregulate their own transporters: PFOA cuts hepatic Oatp1a1/1a4/1b2 and
  PFDA cuts all four while raising serum bile acids 300%, PPARα-dependently
  (Cheng & Klaassen 2008).
- **Non-linearity is not confined to the female rat.** Mouse PFOA is non-linear
  above ~40 mg/kg single or 5 mg/kg/d; mouse PFNA is non-linear while rat PFNA
  is linear. Fits dominated by high-dose data "would overestimate PFAS clearance
  at environmentally relevant doses" — the EPA animal PK authors' own warning
  about the dataset this review leans on.

---

## 10. The practical summary

**For a central human estimate**, Rosato 2024's pooled general-population PFOA
figure of 2.35 y (2.20–2.51) remains the best anchor — explicit cessation
criterion, risk-of-bias scoring, and the subgroup where studies actually agree
(I² = 43%). Batzella 2024 (n = 5,860, 2.36 y, 2.33–2.40) now corroborates it
from a cohort 29× larger than anything prior.

**For dosimetry**, the clearance factor matters more than the half-life, and it
is the least settled quantity in the field: 4.7× spread for PFOA, 24× for PFOS
across agencies. State which you used and why.

**For the species question**, the answer is renal handling, and dose is not a
competing explanation at any exposure humans experience.

**For in vitro binding**, report ranges and the ligand:protein ratio, never a
single Ka.

---

## 11. Toward a QSAR: choosing the endpoint before fitting anything

Everything above treats the species question as a physiology problem. The same
machinery answers a different question — which structures behave which way —
but only if the endpoint is chosen correctly first. This section argues that
the usual endpoint is the wrong one, proposes a replacement, and reports what
the data already on disk say about it.

### 11.1 Half-life is not a QSAR endpoint

From the identity that opens this report,

    t½ = ln2 · Vd / CL

half-life is a composite of three terms that structure affects separately:

| term | what sets it | is it chemistry? |
|---|---|---|
| plasma binding (fu) | albumin affinity, chain length, head group | yes |
| glomerular filtration (GFR) | body size and renal physiology | **no** |
| tubular transport | transporter affinity and direction | yes |

GFR alone differs ~6.5× between mouse and human (§3.4) with no change in
chemistry whatsoever. Fitting descriptors to half-life therefore asks a
structural model to absorb a body-size term, and it cannot. This is the same
failure as fitting descriptors to clearance in absolute units.

### 11.2 The proposed endpoint: the renal handling ratio

    R = CL_renal / (fu · GFR)

R is the measured renal clearance divided by the clearance that free
filtration alone would produce. It is dimensionless, and it has a mechanistic
reading at every value:

- R < 1 — net tubular **reabsorption**; less leaves than was filtered
- R = 1 — pure filtration; no net transport
- R > 1 — net tubular **secretion**; more leaves than was filtered

R divides GFR out, so a mouse and a human are directly comparable, and it
divides fu out, so the binding step is not counted twice. What remains is the
transport step — the part a structural model could plausibly learn. `log10 R`
is the working scale: signed, continuous, and spanning ~3 log units across the
compounds below.

R is a re-expression of the fractional reabsorption axis used throughout §3:
**FR = 1 − R**, so FR = 99.94% is log10 R = −3.2. The reason to prefer R is
arithmetic rather than conceptual. FR crowds against a ceiling at 1 (the
interesting human-to-mouse range compresses into 95–99.99%) and runs
unboundedly negative under secretion, where PFHxA reaches FR = −1.71. In log R
the same range is an ordinary 3-log-unit axis with no special points.

`scripts/qsar_endpoint_table.py` builds the table and runs every number below;
output in `db/qsar/qsar_endpoint_table.csv`.

### 11.3 What the one complete dataset says

Argoul 2026 is the only source that measures CL_renal **and** fu for many
compounds in one experiment, one species, one laboratory — which is the
condition a structure-activity comparison requires and almost nothing else in
this literature satisfies. Male mouse, nine compounds with both terms:

| compound | fluorinated C | head group | ether O | fu % | R | log10 R | handling |
|---|---|---|---|---|---|---|---|
| PFDA | 9 | carboxylate | 0 | 0.42 | 0.0031 | −2.51 | reabsorbed |
| PFHxS | 6 | sulfonate | 0 | 1.30 | 0.0045 | −2.34 | reabsorbed |
| PFNA | 8 | carboxylate | 0 | 0.35 | 0.0171 | −1.77 | reabsorbed |
| PFOA | 7 | carboxylate | 0 | 0.87 | 0.0413 | −1.38 | reabsorbed |
| PFOS | 8 | sulfonate | 0 | 0.25 | 0.0496 | −1.30 | reabsorbed |
| PFBA | 3 | carboxylate | 0 | 77.0 | 0.102 | −0.99 | reabsorbed |
| GenX | 5 | ether-carboxylate | 1 | 26.0 | 0.104 | −0.98 | reabsorbed |
| PFO2OA | 5 | ether-carboxylate | 2 | 10.0 | 1.533 | +0.19 | **secreted** |
| PFHxA | 5 | carboxylate | 0 | 25.0 | 2.712 | +0.43 | **secreted** |

875-fold span. "Fluorinated C" counts carbons bearing fluorine, so it excludes
a PFCA's carboxyl carbon and includes every carbon of a PFSA. That definition
was chosen because it makes PFOS and PFNA the same size on the chain axis — a
prediction the table can falsify.

Three findings, all of which constrain what a QSAR may assume:

**Chain length alone does not order the endpoint.** Within the carboxylates the
series is non-monotonic: C3 −0.99, C5 **+0.43**, C7 −1.38, C8 −1.77, C9 −2.51.
PFHxA breaks it, and not marginally — it crosses from reabsorption into
secretion. Any model using carbon number as its sole descriptor is already
falsified on this dataset.

**The head group carries about 3-fold at matched chain length.** PFNA and PFOS
both have 8 fluorinated carbons; log10 R = −1.77 vs −1.30, a 2.9-fold
difference. So the descriptor definition survives as an approximation — the two
are closer to each other than either is to its own chain neighbours — but the
head group is not negligible.

**Ether oxygens are a strong, non-monotonic descriptor.** Three compounds with
five fluorinated carbons each, differing only in ether substitution:

| ether O | compound | log10 R | handling |
|---|---|---|---|
| 0 | PFHxA | +0.43 | secreted |
| 1 | GenX | −0.98 | reabsorbed |
| 2 | PFO2OA | +0.19 | secreted |

Chain length is held constant and the endpoint moves 1.4 log units, with the
single-ether compound retained and both of its neighbours secreted. Whatever
the mechanism, it is not a monotonic function of ether count, and the
replacement chemicals sit on both sides of the divide.

### 11.3a How far the one dataset can be trusted

Everything in §11.3 rests on a single experiment, and that experiment is
unusual in three ways at once: eleven PFAS dosed as one cocktail, at
0.019–1.55 mg/kg, fitted by nonlinear mixed effects. Each of those could bias
the absolute numbers, so the comparison against the one-compound-at-a-time
literature is worth making explicitly. All comparators below are female mice,
matching Argoul. Built by `scripts/argoul_vs_single_compound.py` into
`db/argoul_vs_single_compound.csv`; drawn as Figure 4 by
`scripts/make_review_figures.py`.

**Clearance and volume are both low, by similar factors.**

| compound | comparator | dose mg/kg | CL ratio | Vd ratio |
|---|---|---|---|---|
| PFOA | Fujii 2015 | 0.128 | 0.38× | 0.59× |
| PFHxS | Sundström 2012 | 1 | 0.49× | 0.65× |
| PFHxS | Sundström 2012 | 20 | 0.34× | 0.42× |
| **geometric mean** | | | **0.40×** | **0.54×** |

Only comparisons with **one attributable dose** are included. Three otherwise
usable comparators pooled two dose levels into a single reported parameter and
are left out: Lou fitted one model across 1 and 10 mg/kg; Tatum-Gibbs reported
a half-life *range* across both doses rather than a value per arm; and the EPA
PFHxA volume is a Vd,β computed from Gannon, whose mouse arm ran 2 *and*
100 mg/kg. Keeping them would mean a dose column that cannot be read as a dose.
The cost is PFOA's second comparator, PFNA entirely, and PFHxA entirely — the
last of which was the one outlier in the other direction (5.13× **high**, and
the only compound Argoul places in net secretion, so its kinetics differ in
kind regardless).

**How much of a check is that, really?** In most of these studies Vd is not
measured apart from clearance — non-compartmental analysis gives
Vz = CL/λz and Vss = CL·MRT, and Lou's clearance is ke·Vd derived in this
review. So for Lou, Sundström and the EPA PFHxA value, "clearance low *and*
volume low" is largely one observation seen twice. Fujii 2015 is the
exception: Vd from Dose/C(0) by back-extrapolation, CL from Dose/AUC, which
are independent, and there CL is 0.38× and Vd 0.59× — both low, independently.
One row is thin evidence, so the claim defended here is the narrower one: the
**rate constant** agrees across the half-life comparisons and the **scale** does not.

`db/argoul_vs_single_compound.csv` holds the comparison itself — values,
ratios, and the doses: Argoul's IV and oral dose for each compound, the
comparator's dose and route, and the gap, which runs 1.6× (Fujii on PFOA,
essentially dose-matched) to 526× (Chang 2012 PFOS at 20 mg/kg). Full
references, DOIs, PMIDs and the route each value took into this review — some
arrive through an agency compilation rather than the primary paper — are in
`db/argoul_vs_single_compound_sources.csv`, joined on the `source` column.

**The half-life, which is their ratio, is not low at all.** Argoul reports MRT
rather than a terminal half-life; MRT = Vss/CL holds in its own table to within
6% for nine of ten compounds, so ln2·MRT is the half-life a one-compartment
system with that Vss and CL would show. Against a measured *terminal* half-life
that is an approximation, and for a compound with a deep second compartment a
lower bound, so read the ratios as indicative.

| compound | comparators | dose mg/kg | ratio range |
|---|---|---|---|
| PFBA | Chang 2008 × 3 | 10 / 30 / 100 | 0.97–1.07× |
| PFHxS | Sundström 2012 × 2 | 1 / 20 | 1.20–1.30× |
| PFOS | Chang 2012 × 2 | 1 / 20 | 0.79–0.98× |
| **all seven** | | | **geometric mean 1.04×** |

Same dose rule as above: Lou (PFOA) and Tatum-Gibbs (PFNA) are excluded because
neither can be tied to one dose. Including their three rows would give ten
comparisons at 0.99× instead of seven at 1.04×, so nothing here turns on the
exclusion — it buys a dose column that means what it says.

Every comparison above is against a study that dosed that compound **alone**.
The EPA pooled fit (Zurlinden 2025) is deliberately not among them and is not
in `db/argoul_vs_single_compound.csv`: it is a Bayesian hierarchical
re-analysis of these same studies, so placing it beside them counts the same
animals twice, and it is a model output rather than a measurement — for PFBA it
returns 0.186 d where every Chang arm sits at 0.116–0.128. For the record its
six ratios would move the geometric mean to 0.95× (span 0.54–1.74×), which does
not change the conclusion.

That exclusion has one cost: **PFHxA has no single-compound female-mouse
comparator**, so it drops out of the half-life comparison entirely and survives
only in the volume comparison. PFBS and PFDS have no comparator either, and
Argoul reports no clearance for PFBS at all.

**What the shape of that result means.** A factor that multiplies Vd and CL
together cancels from ln2·Vd/CL, so "both terms twofold low, half-life
unchanged" is the signature of a shared scaling rather than a disagreement
about elimination. Three candidates fit, and the data in hand cannot separate
them:

- **Dose — the obvious candidate, and it does not survive testing.** §4.4's
  saturable-reabsorption account predicts clearance rising with dose, and the
  cocktail was dosed below its comparators, so this looked like the answer.
  Three things undercut it, all in `scripts/argoul_dose_matching.py`:

  1. *The dose is tuned per compound, not swept.* Argoul ran 0.019–1.55 mg/kg,
     but that range tracks persistence (short-chain near 1.3–1.55, long-chain
     at 0.019–0.080) rather than being a dose-response design. So "10–100×
     below its comparators" is only true compound by compound, with a
     different factor each time.
  2. *Matched that way, dose does not order the discrepancies.* Fujii dosed
     PFOA at 0.128 mg/kg against the cocktail's 0.080 — a gap of 1.6×,
     essentially dose-matched — and has the **worst** agreement of the four
     (0.38×). Lou at 1 mg/kg and above, a gap of ≥12×, has the **best**
     (0.75×). If dose drove the gap, the ordering would run the other way.
  3. *The slope quoted earlier was circular.* It was fitted through all three
     PFHxS points including Argoul's, which carries **79% leverage** on that
     fit — so "the cocktail sits at the bottom of the trend rather than off it"
     was close to guaranteed by construction. Fitting Sundström's two doses
     alone and extrapolating down to 0.094 mg/kg predicts 2.04 mL/kg/d against
     1.30 observed: dose closes part of the gap (0.49× → 0.64×) and about 1.6×
     survives it.

  Dose is a contributor, not the explanation.
- **Cocktail competition.** Eleven PFAS competing for the same reabsorptive
  transporters would *raise* clearance, not lower it, so it argues against the
  observed direction — unless the competition is for plasma binding sites
  instead, which would raise fu and lower the apparent volume.
- **Mixed-effects shrinkage.** Shrinkage toward a population mean compresses
  spread rather than shifting the centre, so it is the weakest of the three.

**What this costs §11.3, and what it does not.** R is built from an absolute
renal clearance, so a 2.5-fold low bias in that clearance is a 2.5-fold low bias
in R — a flat −0.40 log unit offset on every compound. Because it is flat it
changes nothing §11.3 actually claims: the ranking, the 2.9-log-unit span, the
non-monotonicity in chain length, and the head-group and ether contrasts are
all *differences between compounds*, and a shared factor cancels from every one
of them. It does move the compounds relative to the fixed divide at R = 1, but
not far enough to matter: the nearest reabsorbed compound (GenX, −0.98) sits
0.98 log units below the divide, so the bias would have to be 9.6× rather than
2.5× to reclassify anything. What it does cost is §11.4, where an absolute
mouse value is set against an absolute human one — read those gaps as carrying
a factor of about 2.5 on top of the fu uncertainty, in the direction of
understating the mouse.

### 11.4 The sensitivity that decides the species question

R is linear in 1/fu, so the endpoint inherits the binding uncertainty of §6.2
in full. Human PFOA, with CL_renal = 0.03 mL/d/kg and GFR = 2,570 mL/d/kg both
fixed, under the three unbound fractions the literature carries:

| fu | source | log10 R | gap to female mouse PFOA |
|---|---|---|---|
| 0.10 | PBPK models reading ">90% bound" as "≈90% bound" | −3.93 | 354× |
| 0.02 | Han 2012's stated assumption | −3.23 | 71× |
| 0.00061 | Fischer, measured at physiological ligand:protein | −1.72 | **2.2×** |

This is the sharpest consequence of the binding adjudication in §6.2. Under the
assumed fu, the mouse-to-human difference in the transport step is ~71× and the
species gap is a transport-biology problem. Under the measured fu it is ~2× and
the species gap is almost entirely a **binding** problem — which would mean the
transporter work in §7, including the Oatp1a1 question, is explaining a
quantity that barely differs between the species.

The caveat is load-bearing and is the reason this is stated as a sensitivity
rather than a result: Argoul's mouse fu and Fischer's human fu come from
different methods at different ligand:protein ratios, which is precisely the
artefact §6.2 diagnoses. Cross-method fu comparison is not yet legitimate.

That makes one inexpensive experiment decisive for both the species question
and the QSAR: **measure fu for these compounds in mouse, rat and human plasma
by a single method at physiological ligand:protein ratio.** No animals are
required. It is listed in §9 as a gap; §11 raises it to the top of the list,
because every value of R in this section is proportional to it.

### 11.5 What the QSAR cannot yet be fitted on

Nine compounds in one species from one laboratory is a hypothesis generator,
not a training set. The table above is the complete set of PFAS for which
CL_renal and fu were measured together anywhere in this review. The per-compound
material that exists alongside it, and what is missing:

| asset | compounds | what it gives | what it lacks |
|---|---|---|---|
| Argoul 2026 (§3.2a) | 11 | CL, Vss, MRT, F, fu, R | one species, one sex, n small |
| Louisse 2024 + Yang 2010 (§7) | 5 | human OAT1/3/4 Km | no matched fu or CL_renal |
| Han 2012 Table 6 (§7) | 4 | rat Oat1/Oat3/Oatp1a1 Km | rat only; mixed methods |
| Kudo 2001, Ohmori 2003 (§3) | 4 | chain-length elimination, both sexes | no fu, no CL_renal |
| Fischer (§6) | 11 | log D at realistic ratio | no in vivo pairing |

The binding and transporter columns are per-compound and could feed a QSAR
directly; the in vivo column is the bottleneck. Expanding that column — a
second species with CL_renal and fu measured together for the same compound
set — is what turns this from an argument into a model.

## Appendix: the databases

| file | rows | what it answers |
|---|---|---|
| `source_catalogue.csv` | 81 | your citation list, PMID-verified and classified |
| `exposure_vs_halflife.csv` | 184 | chemical, serum, exposure, half-life on one line |
| `master_exposure_halflife.csv` | 46 | dose/serum/half-life/Vd/CL per species × sex |
| `human_halflife_extended.csv` | 142 | human estimates with design and modelling assumptions |
| `regulatory_values.csv` | 173 | what each agency adopted, from which study, under which assumptions |
| `agency_clearance_comparison.csv` | 16 | clearance factors normalised to mL/kg-day |
| `protein_binding.csv` | 273 | in vitro binding, with method and ratio regime |
| `binding_summary_ka.csv` | 71 | association constants, as ranges |
| `binding_summary_logD.csv` | 111 | protein/water distribution at realistic ratio |
| `transporter_kinetics.csv` | 217 | Km/Vmax/IC50 with reabsorption-or-secretion direction |
| `transporter_mechanism.csv` | 83 | expression, hormonal regulation, species comparison |
| `animal_halflife_measured.csv` | 217 | measured animal half-lives by strain, sex, dose, route (187 measured, 30 flagged modelled) |
| `vd_clearance.csv` | 143 | Vd by species, flagged measured/fitted/assumed, with provenance |
| `dose_dependence.csv` | 32 | multi-dose studies, half-life at each level |
| `mouse_rat_decomposition.csv` | 14 | the species gap split into Vd and CL terms |
| `sex_decomposition.csv` | 21 | the same split for the sex difference |
| `reabsorption_axis.csv` | 5 | renal reabsorption and the half-life it predicts |
| `qsar/qsar_endpoint_table.csv` | 11 | per-compound structure descriptors paired with the renal handling ratio R (§11) |

Figures in `figures/`; retrieval logs in `papers/SOURCES_*.md`; per-subtopic
research notes in `research_notes/`.
