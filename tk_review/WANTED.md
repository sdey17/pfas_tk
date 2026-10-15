# Papers still wanted

**Nearly everything is in.** Twenty-five full texts were supplied between
2026-10-01 and 2026-10-03 and are now in `papers/`, extracted into
`db/primary_2026/`, and folded into the report. Everything received has been
removed from this list.

**Received and removed (25):** Argoul/Gayrard 2026 · Thompson 2010 (+5
supplementary + corrigendum) · Kudo 2002 · Kudo 2001 · Lou 2009 · Han 2012 ·
Han 2003 · Yang 2010 · Yang 2009 · Tatum-Gibbs 2011 · Sundström 2012 ·
Cheng 2005 · Cheng 2009 · Buist & Klaassen 2004 · Ohmori 2003 · Louisse 2024 ·
Shi 2016 · Jia 2022 · Maso 2021 · Zurlinden 2025 · Huang 2021 corrigendum ·
Yi 2022 · Delaere 2025 · Kojo 1986 · EPA 2016 PFOS HESD.
Also removed: **Weaver 2010** (listed in error; on disk all along) and
**Cheng 2006** (needed no full text — its own abstract contradicts itself, and
Cheng 2005 and 2009 settle the direction).

## Binding: the one gap that is in the literature, not in this file

A PubMed search for PFAS binding to **mouse serum albumin** returns exactly one
record, and it measures mouse albumin *concentration* under microplastic
co-exposure, not PFAS binding to it (Zhou 2025, doi:10.1016/j.scitotenv.2025.179381).
There is no in vitro PFAS–mouse-albumin binding study to obtain. This is not a
retrieval failure and should not be chased; it is a hole in the field, directly
under this review's own mouse limb. See `scripts/binding_species_and_structure.py`
part 3.

Worth obtaining if cheap, found 2026-10-15 and not yet in `db/protein_binding.csv`:

| study | why | id |
|---|---|---|
| Gao 2025 | BSA binding for 10 PFAS resolved by chain length AND head group, plus ionic-strength and pH dependence; experimental, not docking | PMID 41344664, doi:10.1021/acs.est.5c08779 |
| Gallagher 2023 | QSAR over 24 PFAS HSA affinities; useful as a structure-activity comparator, not as data | PMID 37513249, doi:10.3390/molecules28145375 |
| Zhao 2023 | review of PFAS binding models across albumin, FABP, OATs and membranes | PMID 37126916, doi:10.1016/j.envint.2023.107951 |

The rest of that search returned molecular-docking papers, which
`db/protein_binding.csv` already flags and excludes from the experimental
analysis.

---

Identifiers below are copied from the retrieval logs in `papers/SOURCES_*.md`,
not written from memory. Every item here was attempted and failed.

---

## Still blocks something

### 1. A measured human KT — *no paper known to supply it*

Not a retrieval request, but the single most valuable missing number. §4.3 rests
on a judgement between two parameter families. There are now **six independent
in vitro Km values for PFOA against human transporters** (47–310 µM, Yang 2010
and Louisse 2024, two laboratories, agreeing within a factor of 7), against a
PBPK transport-affinity constant of **0.133 µM** — **354–2,336× lower**. Every
KT in Han 2012 Table 7 is *fitted to plasma curves*, not measured. The in vitro
side is now strongly favoured, but a directly measured human renal reabsorptive
KT would close the question outright.

### 2. Kemper 2003 — *PFOA: toxicokinetics in the rat*, DuPont Haskell
Unpublished; EPA docket.

The rat PFOA dose series behind §4.2's female dose-dependence — the one place in
the whole dataset where dose genuinely moves half-life (female t½ 3.2 → 16.2 h
across 0.1 → 25 mg/kg). Grey literature, read only through ATSDR's tabulation,
yet it underpins Wambaugh 2013 and Worley & Fisher 2015.

---

## Would fill a hole in the coverage grid

85% of the 1,716 chemical × species × parameter cells have no data.

| paper | identifier | fills |
|---|---|---|
| **Hanhijärvi et al. 1988**, beagle dog — *"A proposed species difference in the renal excretion of perfluoro-octanoic acid in the beagle dog and rat"* | **EPA HERO 5412773**. ⚠️ **NOT** doi:10.1007/978-3-642-71248-7_96 — that DOI, which this list carried until 2026-10-03, is **Kojo, Hanhijärvi, Ylinen & Kosma 1986**, a Wistar *rat* study with no dog data. Hanhijärvi co-authors both, which is how they were conflated; every previous retrieval attempt was chasing the wrong paper | **the entire dog column** — dog has 2 of 44 cells. Han 2012 Table 4's compiled values (CLR 50.8 F / 43 M mL/d/kg, 52%/59% reabsorbed) are the current substitute |
| **Kerstner-Wood 2003**, SRI contract report | no DOI/PMID | the main source for neat human plasma f_unbound of PFOS/PFOA/PFHxS. Less critical now that Han 2003 is in hand and the binding conflict is resolved (§6.2), but still a primary source nobody has read |

---

## Not worth your time

- **Griffith & Long 1980** — the obvious place to expect dog data. Checked: rat
  and mouse acute toxicity only, no TK.
- **ATSDR 2021 for dog values** — contains "dog"/"beagle"/"canine" only inside
  its literature-search Boolean strings.
- **cC6O4** — already captured via Fustinoni 2023.
- **Katakura 2007** — the only study testing Npt2 and Mrp2 for PFOA. Not indexed
  in PubMed or Europe PMC; *J Health Sci* has no located DOI. Known second-hand
  through OEHHA Table A6.3 and Weaver 2010's Discussion.
- Anything already in `papers/` — 186 full-text extractions.

---

## If you can only chase one

**Hanhijärvi 1988.** Everything else on this list changes a number or
corroborates a claim. That one is a whole species: the dog is the only mammal on
the reabsorption axis with just two data cells, it sits in the interesting middle
of the axis (52–59% reabsorbed, between the macaque and the female rat), and two
agencies tabulate it inconsistently with no way to adjudicate. It is also the
only item here that probably cannot be solved online — an interlibrary loan is
the realistic route.
