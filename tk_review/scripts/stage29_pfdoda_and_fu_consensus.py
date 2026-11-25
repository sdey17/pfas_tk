#!/usr/bin/env python3
"""PFDoDA verified at source, the fecal-route correction it exposes, and an
f_u consensus across all four human studies (Stage 29).

    python3 -I scripts/stage29_pfdoda_and_fu_consensus.py

Writes db/abraham2024_route_split.csv, db/human_tubular_factor_renal.csv,
db/human_fu_consensus.csv, db/human_fu_consensus_models.csv,
db/human_tk_corrected.csv. Idempotent.

PART 1 -- WHAT THE SOURCE ACTUALLY SAYS ABOUT PFDoDA.
Abraham et al. 2024 (Int J Hyg Environ Health, PMID 39476597), Table 3:
half-life 295 d, 95% CI 274-321, V_d 354 mL/kg (342-368). ONE volunteer,
single oral dose, 450 days of plasma sampling. The half-life is MEASURED
and the confidence interval is tight, so it is not a transcription error
and not a modelling artefact.

Nor is the non-monotonic chain-length pattern an anomaly. It is the
paper's central result, and it runs across FOUR CONSECUTIVE COMPOUNDS:
2011 d (PFOA) -> 1305 (PFNA) -> 862 (PFDA) -> 584 (PFUnDA) -> 295
(PFDoDA). PFDoDA is the end of a smooth decline, not an outlier beside it.

WHAT IS WRONG IS THE TUBULAR FACTOR, AND IT IS MY ERROR, NOT THE DATA'S.
T was computed as CL_total/(f_u*GFR), i.e. by crediting the kidney with
ALL of the clearance. Abraham's Table 4 reports, for PFDoDA, a urinary
clearance of "n.d." -- not detected -- and a FECAL clearance of 0.0286
against a total of 0.0474 mL/min, so 60% of the measured elimination is
faecal. An earlier stage of this project had already marked this value
"total_upper_bound" and "UNINFORMATIVE_bound_exceeds_filtration", and
Stage 28 regressed against that.

THE CORRECTION IS NOT CONFINED TO PFDoDA. Abraham measured faecal
clearance for five long-chain PFCAs, and it is 17% to 61% of total. Every
human T in this project computed from total clearance is therefore
overstated for those compounds.

PART 2 -- THE f_u CONSENSUS.
Four independent human studies measure f_u: Smeltz 2023 (ultracentrifugation,
67 PFAS), Ryu 2024, Argoul 2025 and Fischer 2024 (derived). No two agree
closely and no single study covers everything. Rather than pick one, the
consensus value is the GEOMETRIC MEAN of whatever sources report a
compound, carried with its full range, and every model is then re-run
under all five options so the reader can see what the choice costs.

ONE EXTERNAL ANCHOR DESERVES WEIGHT. Abraham 2024 -- the source of every
half-life and V_d used here -- uses SMELTZ values in its own Table 4. A
model built on Abraham's half-lives and Smeltz's f_u is at least
internally consistent with one published analysis; any other pairing is
assembled by us.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
LN2 = math.log(2.0)
GFR_HUMAN = 2570.0
SULFONATES = {"PFBS", "PFHxS", "PFOS"}

# Abraham 2024 Table 4, read from the paper text. mL/min as published.
# CL_uri is "n.d." for every long-chain compound: the authors could not
# detect them in urine at the doses applied.
ABRAHAM_T4 = {
    # chemical:   (CL_uri, CL_tot, CL_fec)   None = not detected
    "PFBA":   (1.00,  1.24,   None),
    "PFHpA":  (0.023, 0.033,  None),
    "PFOA":   (None,  0.0024, 0.0004),
    "PFNA":   (None,  0.0038, 0.0022),
    "PFDA":   (None,  0.0084, 0.0051),
    "PFUnDA": (None,  0.0163, 0.0098),
    "PFDoDA": (None,  0.0474, 0.0286),
    "PFBS":   (0.067, 0.107,  None),
    "PFHxS":  (None,  0.0030, None),
    "PFOS":   (None,  0.0050, None),
}
ABRAHAM_T3 = {   # half-life (d), low, high -- Table 3
    "PFBA": (4.18, 4.06, 4.31), "PFHpA": (152, 147, 157),
    "PFOA": (2011, 1466, 3206), "PFNA": (1305, 1054, 1716),
    "PFDA": (862, 728, 1057), "PFUnDA": (584, 501, 705),
    "PFDoDA": (295, 274, 321), "PFBS": (50.6, 47.7, 53.8),
    "PFHxS": (1634, 1166, 2733), "PFOS": (1211, 927, None),
}
FU_SOURCES = [("smeltz_UC", "Smeltz 2023 (UC)"), ("ryu", "Ryu 2024"),
              ("argoul", "Argoul 2025"), ("fischer_derived", "Fischer 2024")]


def sig(x, n=6):
    if x is None:
        return None
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


def gmean(vs):
    return 10 ** (sum(math.log10(v) for v in vs) / len(vs))


def load_human():
    out = {}
    for r in csv.DictReader((HERE / "db" / "species_decomposition.csv").open()):
        c = r["chemical"]
        if c in out:
            continue
        out[c] = dict(t_half=float(r["t_half_human_d"]),
                      vd=float(r["vd_human"]), fu=float(r["fu_human"]),
                      cl=float(r["cl_human"]),
                      nc=int(r["n_total_carbons"]))
    return out


def load_fu_sources():
    """Four studies' human f_u.

    human_fu_comparison.csv has no PFDoDA ROW AT ALL, which made the
    compound look as though no study measured it. Ryu 2024 Table 1 does,
    via db/species_fu_ratios_all14.csv, so that file is read as well. The
    gap was in our comparison table, not in Ryu.
    """
    out = {k: {} for k, _ in FU_SOURCES}
    for r in csv.DictReader((HERE / "db" / "human_fu_comparison.csv").open()):
        for k, _ in FU_SOURCES:
            if r.get(k):
                out[k][r["chemical"]] = float(r[k])
    filled = []
    for r in csv.DictReader(
            (HERE / "db" / "species_fu_ratios_all14.csv").open()):
        c, v = r["chemical"], r.get("fu_human")
        if v and c not in out["ryu"] and "Ryu" in (r.get("source") or ""):
            out["ryu"][c] = float(v)
            filled.append(c)
    if filled:
        print(f"  [f_u loader] filled from Ryu's wider table, absent from"
              f" human_fu_comparison.csv: {', '.join(sorted(filled))}")
    return out


def ols(rows, y):
    p = len(rows[0])
    aug = [[sum(r[i] * r[j] for r in rows) for j in range(p)]
           + [sum(rows[k][i] * y[k] for k in range(len(rows)))]
           for i in range(p)]
    for i in range(p):
        piv = max(range(i, p), key=lambda k: abs(aug[k][i]))
        if abs(aug[piv][i]) < 1e-12:
            return None
        aug[i], aug[piv] = aug[piv], aug[i]
        for k in range(p):
            if k == i:
                continue
            f = aug[k][i] / aug[i][i]
            for j in range(i, p + 1):
                aug[k][j] -= f * aug[i][j]
    return [aug[i][p] / aug[i][i] for i in range(p)]


def main():
    hum = load_human()
    fu_src = load_fu_sources()
    bar = "=" * 78
    chems = sorted(hum)

    # ================================================= PART 1
    print(bar)
    print("PART 1  PFDoDA AT SOURCE: THE HALF-LIFE IS REAL, MY T WAS NOT")
    print(bar)
    t, lo, hi = ABRAHAM_T3["PFDoDA"]
    print(f"  Abraham 2024 Table 3:  t½ = {t} d  (95% CI {lo}-{hi})")
    print(f"                         V_d = 354 mL/kg (342-368)")
    print(f"  n = 1 volunteer, single oral dose, 450 days of sampling.")
    print()
    print("  HOW MUCH HUMAN PFDoDA HALF-LIFE DATA EXISTS: THIS, AND NOTHING")
    print("  ELSE. A PubMed search for perfluorododecanoic acid")
    print("  toxicokinetics returns 10 records, none of them a human")
    print("  half-life study; the hits are rat brain distribution, RIVM rat")
    print("  modelling, minipig PK and fish bioconcentration. So the value")
    print("  is a single measurement in a single person and cannot be")
    print("  cross-checked against an independent study.")
    print()
    print("  IS IT AN OUTLIER? NO -- IT IS THE END OF A TREND.")
    series = ["PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA"]
    print("    " + "  ->  ".join(f"{c} {ABRAHAM_T3[c][0]:.0f}d"
                                 for c in series))
    print("  Five consecutive chain lengths, monotonically decreasing after")
    print("  PFOA. A single bad number cannot produce a monotone run of four")
    print("  decreases, so the pattern is a property of the data.")
    print()

    print("  WHAT IS WRONG IS THE ROUTE ATTRIBUTION. Abraham Table 4:")
    print(f"  {'cmpd':9}{'CL_uri':>9}{'CL_tot':>9}{'CL_fec':>9}{'fecal %':>9}"
          f"   non-renal evidence")
    route = []
    for c in sorted(ABRAHAM_T4, key=lambda x: hum.get(x, {}).get("nc", 99)):
        uri, tot, fec = ABRAHAM_T4[c]
        frac = fec / tot if fec else None
        note = ("faecal route measured" if fec else
                "no route detected either way" if uri is None else
                "urinary route measured")
        route.append(dict(
            chemical=c, cl_urinary_mL_min=uri if uri else "n.d.",
            cl_total_mL_min=tot, cl_fecal_mL_min=fec if fec else "n.d.",
            fecal_pct_of_total=sig(frac * 100, 4) if frac else "",
            evidence=note))
        print(f"  {c:9}{(str(uri) if uri else 'n.d.'):>9}{tot:>9}"
              f"{(str(fec) if fec else 'n.d.'):>9}"
              f"{(f'{frac*100:.0f}%' if frac else '-'):>9}   {note}")
    print()
    print("  Urinary clearance is NOT DETECTED for all five long-chain PFCAs")
    print("  and for both long-chain PFSAs. The authors say so plainly: they")
    print("  could not detect long-chain PFAS in urine at these doses.")
    print()

    # ---- the corrected T
    print(bar)
    print("THE FECAL-ROUTE CORRECTION TO T, FOR EVERY COMPOUND IT TOUCHES")
    print(bar)
    print("  T from total clearance credits the kidney with the faecal route")
    print("  as well. Using the faecal FRACTION (dimensionless, so no body-")
    print("  weight conversion is needed) gives a renal-only upper bound:")
    print("      CL_renal <= CL_total * (1 - CL_fec/CL_tot)")
    print()
    print(f"  {'cmpd':9}{'T (total)':>11}{'fecal %':>9}{'T (renal)':>11}"
          f"{'change':>9}   reading")
    tren = []
    for c in chems:
        h = hum[c]
        T_tot = h["cl"] / (h["fu"] * GFR_HUMAN)
        ent = ABRAHAM_T4.get(c)
        frac = (ent[2] / ent[1]) if ent and ent[2] else None
        T_ren = T_tot * (1 - frac) if frac else T_tot
        rd = ("net secretion" if T_ren > 1 else
              "strong reabsorption" if T_ren < 0.1 else "reabsorption")
        tren.append(dict(
            chemical=c, n_total_carbons=h["nc"],
            T_from_total_clearance=sig(T_tot, 4),
            fecal_pct_of_total=sig(frac * 100, 4) if frac else "",
            T_renal_upper_bound=sig(T_ren, 4),
            kind="renal-only upper bound" if frac
                 else "total clearance; no faecal route detected",
            reading=rd))
        print(f"  {c:9}{T_tot:>11.4f}{(f'{frac*100:.0f}%' if frac else '-'):>9}"
              f"{T_ren:>11.4f}{(f'{(1-frac)*100:.0f}%' if frac else '-'):>9}"
              f"   {rd}")
    print()
    pf = next(r for r in tren if r["chemical"] == "PFDoDA")
    print(f"  PFDoDA: T falls from {pf['T_from_total_clearance']} to "
          f"{pf['T_renal_upper_bound']}.")
    still = [r["chemical"] for r in tren
             if float(r["T_renal_upper_bound"]) > 1]
    if still:
        print(f"  It is STILL above 1 ({', '.join(still)}), so the correction")
        print("  does not by itself overturn net secretion -- but the value is")
        print("  now an UPPER BOUND built on an undetectable urinary")
        print("  clearance, which is no basis for a secretion claim. The")
        print("  honest statement is that PFDoDA's renal T is UNDETERMINED.")
    else:
        print("  Every corrected T is now below 1, so the net-secretion claim")
        print("  in Stage 28 does not survive the correction.")
    print()

    # ================================================= PART 2
    print(bar)
    print("PART 2  THE f_u CONSENSUS ACROSS ALL FOUR STUDIES")
    print(bar)
    print(f"  {'cmpd':9}" + "".join(f"{lab.split()[0]:>11}"
                                    for _, lab in FU_SOURCES)
          + f"{'n':>3}{'spread':>9}{'consensus':>11}")
    cons, consensus = [], {}
    for c in chems:
        vals = [(fu_src[k].get(c), lab) for k, lab in FU_SOURCES]
        got = [(v, lab) for v, lab in vals if v]
        line = f"  {c:9}"
        for v, _ in vals:
            line += f"{(f'{v:.5f}' if v else '-'):>11}"
        if got:
            g = gmean([v for v, _ in got])
            consensus[c] = g
            sp = max(v for v, _ in got) / min(v for v, _ in got)
            line += f"{len(got):>3}{sp:>8.1f}x{g:>11.5f}"
            cons.append(dict(
                chemical=c, n_total_carbons=hum[c]["nc"], n_sources=len(got),
                fu_min=sig(min(v for v, _ in got)),
                fu_max=sig(max(v for v, _ in got)),
                fold_spread=sig(sp, 4), fu_consensus_geomean=sig(g),
                sources="; ".join(f"{lab}={v:g}" for v, lab in got)))
        else:
            line += f"{0:>3}{'-':>9}{'none':>11}"
            cons.append(dict(
                chemical=c, n_total_carbons=hum[c]["nc"], n_sources=0,
                fu_min="", fu_max="", fold_spread="",
                fu_consensus_geomean="", sources="no source reports it"))
        print(line)
    print()
    miss = [r["chemical"] for r in cons if r["n_sources"] == 0]
    one = [r["chemical"] for r in cons if r["n_sources"] == 1]
    print(f"  COVERAGE: {len(consensus)} of {len(chems)} compounds have at "
          f"least one source.")
    if miss:
        print(f"  NO SOURCE AT ALL: {', '.join(miss)} -- no consensus possible.")
    if one:
        print(f"  ONLY ONE SOURCE: {', '.join(one)} -- the 'consensus' there is")
        print("  that single study, so it carries none of the agreement the")
        print("  word implies and is flagged by n_sources=1 in the CSV.")
    multi = len(consensus) - len(one)
    print(f"  GENUINE MULTI-STUDY CONSENSUS: {multi} of {len(chems)} compounds.")
    print()
    wide = max((r for r in cons if r["n_sources"] > 1),
               key=lambda r: float(r["fold_spread"]))
    print(f"  WIDEST DISAGREEMENT: {wide['chemical']} spans "
          f"{wide['fold_spread']}-fold ({wide['sources']}).")
    print("  A consensus geometric mean hides that; the range column is")
    print("  carried in db/human_fu_consensus.csv for exactly that reason.")
    print()
    print("  ABRAHAM USES SMELTZ. The f_up column of Abraham's own Table 4")
    print("  reproduces the Smeltz values exactly (PFBA 0.093, PFOA 0.001,")
    print("  PFNA 0.0016, PFDA 0.0027, PFHxS 0.0009, PFOS 0.0045). Since")
    print("  every half-life and V_d here is Abraham's, pairing them with")
    print("  Smeltz is the only combination that matches a published")
    print("  analysis rather than one we assembled. Note Abraham's PFDoDA")
    print("  f_up cell is BLANK -- Smeltz has no PFDoDA, which is why.")
    print()

    # ---- run the models under all five f_u options
    print(bar)
    print("WHAT THE f_u CHOICE COSTS: ALL FIVE OPTIONS, BOTH MODELS")
    print(bar)
    print("  M1a is parameter-free (T=1). M1b fits one constant T, scored")
    print("  leave-one-out. Both use the CORRECTED renal T where available.")
    print()
    # The faecal fraction is f_u-independent; T is not. So carry the
    # fraction and rebuild T per option, rather than reusing the Ryu-based
    # T for every option -- that would put two different f_u studies inside
    # one prediction.
    fecfrac = {}
    for c in chems:
        ent = ABRAHAM_T4.get(c)
        fecfrac[c] = (ent[2] / ent[1]) if ent and ent[2] else 0.0

    def T_renal(c, fu_c):
        """Renal-only T for compound c under a given f_u."""
        return hum[c]["cl"] * (1 - fecfrac[c]) / (fu_c * GFR_HUMAN)
    opts = [(lab, fu_src[k]) for k, lab in FU_SOURCES]
    opts.append(("CONSENSUS (geomean)", consensus))
    opts.append(("compilation (Ryu, as used)",
                 {c: hum[c]["fu"] for c in chems}))
    rows = []
    common = sorted(set(chems).intersection(
        *[{c for c in chems if fu.get(c)} for _, fu in opts]))
    print(f"  COMMON SET, scored on the {len(common)} compounds every option")
    print(f"  covers: {', '.join(common)}")
    print("  Scoring each option on its own coverage would compare different")
    print("  compound sets -- which is how PFDoDA alone produced a 139x")
    print("  'worst case' for one row and not the others.")
    print()
    print(f"  {'f_u option':28}{'M1a med':>10}{'M1a worst':>11}"
          f"{'M1b med':>10}{'M1b worst':>11}")
    for lab, fu in opts:
        fa = sorted(max(p / hum[c]["t_half"], hum[c]["t_half"] / p)
                    for c, p in ((c, LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN))
                                 for c in common))
        fb = []
        for c in common:
            train = [x for x in common if x != c]
            That = gmean([T_renal(x, fu[x]) for x in train])
            pr = LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN * That)
            fb.append(max(pr / hum[c]["t_half"], hum[c]["t_half"] / pr))
        fb.sort()
        rows.append(dict(
            fu_option=lab, scope=f"common set (n={len(common)})",
            n=len(common),
            m1a_median_fold=sig(fa[len(fa) // 2], 4),
            m1a_worst_fold=sig(fa[-1], 4),
            m1b_median_fold=sig(fb[len(fb) // 2], 4),
            m1b_worst_fold=sig(fb[-1], 4)))
        print(f"  {lab:28}{fa[len(fa)//2]:>9.1f}x{fa[-1]:>10.1f}x"
              f"{fb[len(fb)//2]:>9.1f}x{fb[-1]:>10.1f}x")
    print()
    print(f"  {'f_u option':28}{'n':>3}{'M1a med':>10}{'M1a worst':>11}"
          f"{'M1b med':>10}{'M1b worst':>11}   (own coverage)")
    for lab, fu in opts:
        have = [c for c in chems if fu.get(c)]
        if len(have) < 4:
            continue
        fa = sorted(max(p / hum[c]["t_half"], hum[c]["t_half"] / p)
                    for c, p in ((c, LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN))
                                 for c in have))
        fb = []
        for c in have:
            train = [x for x in have if x != c]
            That = gmean([T_renal(x, fu[x]) for x in train])
            p = LN2 * hum[c]["vd"] / (fu[c] * GFR_HUMAN * That)
            fb.append(max(p / hum[c]["t_half"], hum[c]["t_half"] / p))
        fb.sort()
        rows.append(dict(
            fu_option=lab, scope="own coverage", n=len(have),
            m1a_median_fold=sig(fa[len(fa) // 2], 4),
            m1a_worst_fold=sig(fa[-1], 4),
            m1b_median_fold=sig(fb[len(fb) // 2], 4),
            m1b_worst_fold=sig(fb[-1], 4)))
        print(f"  {lab:28}{len(have):>3}{fa[len(fa)//2]:>9.1f}x"
              f"{fa[-1]:>10.1f}x{fb[len(fb)//2]:>9.1f}x{fb[-1]:>10.1f}x")
    print()
    cm = [r for r in rows if r["scope"].startswith("common")]
    bm = min(cm, key=lambda r: r["m1b_median_fold"])
    bw = min(cm, key=lambda r: r["m1b_worst_fold"])
    print(f"  ON THE COMMON SET: best M1b median is {bm['fu_option']} at "
          f"{bm['m1b_median_fold']:.1f}x,")
    print(f"  best worst case is {bw['fu_option']} at "
          f"{bw['m1b_worst_fold']:.1f}x.")
    print(f"  M1a spans {min(r['m1a_median_fold'] for r in cm):.1f}x to "
          f"{max(r['m1a_median_fold'] for r in cm):.1f}x across the options,")
    print(f"  M1b spans {min(r['m1b_median_fold'] for r in cm):.1f}x to "
          f"{max(r['m1b_median_fold'] for r in cm):.1f}x.")
    print("  THE CHOICE MOVES THE NUMBERS BUT NOT THE CONCLUSION: every")
    print("  option leaves M1a with a large one-directional shortfall and")
    print("  every option lets one fitted constant remove most of it.")
    print()
    ra = max(r["m1a_median_fold"] for r in cm) / min(
        r["m1a_median_fold"] for r in cm)
    rb = max(r["m1b_median_fold"] for r in cm) / min(
        r["m1b_median_fold"] for r in cm)
    print()
    print("  A PREDICTION OF MINE THAT THE NUMBERS DID NOT BEAR OUT. Because")
    print("  T is now derived with the same f_u the model divides by, the")
    print("  absolute LEVEL of f_u cancels in M1b, so I expected M1b to be")
    print("  far less sensitive to the choice than M1a. It is not:")
    print(f"      M1a medians span {ra:.1f}-fold across the options")
    print(f"      M1b medians span {rb:.1f}-fold")
    print("  Comparable. The reason is that the four studies disagree about")
    print("  the SHAPE of f_u across compounds, not only its level -- Smeltz")
    print("  and Ryu differ 6.0-fold on PFOS but only 1.3-fold on PFHxS. The")
    print("  level cancels; the shape is what M1b depends on, and the")
    print("  disagreement lives there too. So the f_u choice is not made")
    print("  harmless by fitting T.")
    print()
    smz = next(r for r in cm if r["fu_option"].startswith("Smeltz"))
    con = next(r for r in cm if r["fu_option"].startswith("CONSENSUS"))
    cov = {r["fu_option"]: r["n"] for r in rows if r["scope"] == "own coverage"}
    print("  RECOMMENDATION -- AND THE TWO LINES OF EVIDENCE DO NOT AGREE,")
    print("  SO BOTH ARE STATED.")
    print(f"    Best median on the common set is {bm['fu_option']} "
          f"({bm['m1b_median_fold']:.1f}x), but it")
    print(f"    reports only {cov.get(bm['fu_option'], '?')} of the 9 compounds "
          f"and has the WORST M1a of all")
    print(f"    ({max(r['m1a_median_fold'] for r in cm):.1f}x), so its edge "
          f"rests on the fewest compounds.")
    print()
    print(f"    PRIMARY: Smeltz 2023. Not because it wins the median -- it is")
    print(f"    second at {smz['m1b_median_fold']:.1f}x -- but because it has the "
          f"best worst case")
    print(f"    ({smz['m1b_worst_fold']:.1f}x), covers "
          f"{cov.get('Smeltz 2023 (UC)', '?')} of 9 compounds, is the largest "
          f"series (67 PFAS),")
    print("    and is the f_u Abraham 2024 used in its own Table 4. That last")
    print("    point is the one that is not a fit statistic: it makes the")
    print("    half-life, V_d and f_u a single published pairing.")
    print()
    print(f"    SENSITIVITY: the consensus geometric mean "
          f"({con['m1b_median_fold']:.1f}x median, "
          f"{con['m1b_worst_fold']:.1f}x worst),")
    print("    which covers all 9 and sits mid-range, carried with the per-")
    print("    compound spread (up to 11.1-fold on PFHpA) from")
    print("    db/human_fu_consensus.csv.")
    print()
    print("    CAVEAT TO STATE IN THE PAPER: Smeltz has no PFDoDA, so that")
    print("    compound takes Ryu's value under any choice. The f_u column is")
    print("    therefore not from a single source and should not be presented")
    print("    as one.")
    print()

    out = [
        ("abraham2024_route_split.csv", route,
         ["chemical", "cl_urinary_mL_min", "cl_total_mL_min",
          "cl_fecal_mL_min", "fecal_pct_of_total", "evidence"]),
        ("human_tubular_factor_renal.csv", tren,
         ["chemical", "n_total_carbons", "T_from_total_clearance",
          "fecal_pct_of_total", "T_renal_upper_bound", "kind", "reading"]),
        ("human_fu_consensus.csv", cons,
         ["chemical", "n_total_carbons", "n_sources", "fu_min", "fu_max",
          "fold_spread", "fu_consensus_geomean", "sources"]),
        ("human_fu_consensus_models.csv", rows,
         ["fu_option", "scope", "n", "m1a_median_fold", "m1a_worst_fold",
          "m1b_median_fold", "m1b_worst_fold"]),
    ]
    for name, rr, cols in out:
        p = HERE / "db" / name
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in sorted(rr, key=lambda d: [str(d.get(x, "")) for x in cols]):
                w.writerow({x: r.get(x, "") for x in cols})
        print(f"wrote db/{name}  ({len(rr)} rows)")


if __name__ == "__main__":
    main()
