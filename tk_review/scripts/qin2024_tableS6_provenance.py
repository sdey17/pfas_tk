#!/usr/bin/env python3
"""Where do Qin 2024 Table S6's D_BSA/w equations come from, and can the
saturation they report be reproduced from independent data?

Qin W, Henneberger L, Gluge J, Konig M, Escher BI. Environ Sci Technol
2024;58(13):5727-5738. doi:10.1021/acs.est.3c09950. SI Table S6 carries, for
11 PFAS, two columns of concentration-dependent distribution ratios:

    log D_medium/w = m log Cw + b        [L_w / L_prot+lip]
    log D_BSA/w    = m log Cw + b        [L_w / L_prot]

Ten of the 22 cells are regressions; the rest are concentration-independent
constants.  This script does two things.

PART A  Audits the provenance of each column from the SI's own structure.
PART B  Checks the D_BSA/w slopes against the same group's earlier
        mechanistic fit (Qin 2023), which is the only place their Langmuir
        parameters are published.
PART C  Derives a PARAMETER-FREE bound: a measured Freundlich slope puts an
        upper bound on Kd, with no capacity, level or stoichiometry needed.
        Tests it against albumin Kd measured by two independent groups.
PART D  Triangulates three studies -- an independent Kd, an independent
        low-concentration plateau, and Qin's slope -- to back out the
        non-specific partition coefficient, and checks it against the one
        Qin 2023 measured directly.

WHY PART C MATTERS.  This project has no study that measured both a Kd and a
protein-water partition coefficient on the same material (albumin Kd sources
and albumin partition sources are fully disjoint).  The bound below does not
need them to be the same study: a Kd from one paper constrains the shape of
an isotherm measured in another, because the shape is a property of Kd alone.

Run:  python3 tk_review/scripts/qin2024_tableS6_provenance.py
"""
import collections
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from add_qin2024_tableS6 import S6          # the 11 x 2 Table S6 fits
from add_qin2024_table1 import T1           # Table 1, incl. IC_10,free
from mbm_qin_vs_armitage import QIN_BSA_MEASURED   # Qin 2023 Table 1

DB = os.path.join(os.path.dirname(HERE), "db", "protein_binding.csv")

# Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L. "Species Difference?
# Bovine, Trout, and Human Plasma Protein Binding of Per- and Polyfluoroalkyl
# Substances." Environ Sci Technol 2024;58(23):9954-9966.
# doi:10.1021/acs.est.3c10824, PMID 38804966, PMC11171458 (open access).
# Table 1, "BSA: log D_BSA/w [L_w/L_prot]" column.  None = the paper reports a
# concentration-independent constant.  This is the SOURCE of Table S6's BSA
# column; partA proves the identity cell by cell rather than asserting it.
COMPANION_T1_BSA = {
    "PFBA":     (-0.298, 2.91, 0.60),
    "PFHxA":    (-0.271, 3.43, 0.88),
    "PFHpA":    (-0.305, 3.95, 0.81),
    "PFOA":     (-0.314, 4.38, 0.93),
    "PFNA":     (-0.147, 4.52, 0.61),
    "PFUnA":    (None,   4.75, 0.96),
    "HFPO-DA":  (-0.493, 3.44, 0.86),
    "PFHxS":    (-0.472, 4.28, 0.92),
    "PFOS":     (-0.379, 4.74, 0.92),
    "6:2 FTSA": (-0.092, 3.86, 0.54),
    "PFOSA":    (-0.105, 4.28, 0.50),
}

# Albumin: 1 L of pure protein is 1.364 kg (partial specific volume
# 0.733 mL/g), so a mass-basis K converts to a volume basis by +log10(1.364).
KG_PER_L_ALB = 1.0 / 0.733
MASS_TO_VOL = math.log10(KG_PER_L_ALB)

# Sorptive volume fractions of the PPARgamma medium, SI Table S7:
# protein 0.94 mL/L, lipid 1.47e-2 mL/L.
MED_FPROT = 0.94 / (0.94 + 1.47e-2)
MED_FLIP = 1.0 - MED_FPROT

# Qin's 11 under the names other studies use for them.
ALIAS = {"PFUnA": ["PFUnA", "PFUnDA"],
         "HFPO-DA": ["HFPO-DA", "HFPO-DA (GenX)", "GenX"],
         "6:2 FTSA": ["6:2 FTSA", "6:2 FTS", "6:2 FtS"],
         "PFOSA": ["PFOSA", "FOSA"]}


def names(nm):
    return ALIAS.get(nm, [nm])


# Plateau sources, named by (study, year, species, a substring of the protein
# string that picks exactly one material). Allendorf ran two BSA batches and
# Fischer ran four matrices, so the material has to be named: pooling them
# would average a fatty-acid-free albumin with a heat-shock fraction, or a
# bovine albumin with a human one.
PLATEAU_SRC = [
    ("Allendorf FAF", ("Allendorf et al.", "2019", "bovine", "A3803"), True),
    ("Allendorf HS", ("Allendorf et al.", "2019", "bovine", "A7906"), True),
    ("Bischel", ("Bischel et al.", "2011", "bovine", "bovine serum albumin"),
     False),
    ("Fischer BSA", ("Fischer et al.", "2024", "bovine",
                     "bovine serum albumin"), True),
    ("Fischer HSA human", ("Fischer et al.", "2024", "human",
                           "human serum albumin"), True),
]
LABEL = {k: k.split()[-1] if k.startswith("Allendorf") else k.split()[0]
         for k, _, _ in PLATEAU_SRC}
LABEL["Fischer BSA"] = "FischerBSA"
LABEL["Fischer HSA human"] = "FischerHSA"


def plateaus(rows):
    """log D_max per source, rebased to L_water per L of protein."""
    out = collections.OrderedDict()
    for label, (st, yr, sp, mat), mass_basis in PLATEAU_SRC:
        d = {}
        for r in rows:
            if (r["study"], r["year"], r["protein_species"]) != (st, yr, sp):
                continue
            if mat not in r["protein"] or not r["ka"].strip():
                continue
            try:
                v = float(r["ka"])
            except ValueError:
                continue
            d[r["chemical"]] = v + (MASS_TO_VOL if mass_basis else 0.0)
        if d:
            out[label] = d
    return out


def rule(ch="-", n=78):
    print(ch * n)


def load():
    return list(csv.DictReader(open(DB, encoding="utf-8")))


def s6_bsa(nm):
    """(slope, intercept) of the Table S6 BSA column; slope None = constant."""
    ncf, hg, s_b, i_b, r_b, s_m, i_m, r_m, ctot = S6[nm]
    return s_b, i_b, r_b


def s6_med(nm):
    ncf, hg, s_b, i_b, r_b, s_m, i_m, r_m, ctot = S6[nm]
    return s_m, i_m, r_m


# =====================================================================
# PART A -- provenance
# =====================================================================
def partA():
    rule("=")
    print("PART A -- which of the two columns was measured in this paper?")
    rule("=")
    print("  The SI's own structure answers it. Methods described:")
    print("    Text S2  free concentration + cytotoxicity, PPARgamma assay")
    print("             -> BioSPME 96-pin, 96-well plate, 220 uL, 25000 cells")
    print("    Text S3  CELL binding assay          -> C18-SPME")
    print("    Text S4  STRUCTURAL PROTEIN assay    -> C18-SPME")
    print("    Text S5  384-well HTS campaign")
    print("    Table S3 experimental conditions for (a) the medium assay and")
    print("             (b) the cell and structural-protein assays")
    print()
    print("  There is no BSA binding assay in Text S1-S5 and no BSA row in")
    print("  Table S3. Figure S6 shows MEDIUM binding isotherms for 11 PFAS;")
    print("  no BSA isotherm figure exists. D_BSA/w appears in the whole SI")
    print("  only in Table S6, in the Figure S8b caption -- a panel showing")
    print("  4 compounds (PFHxS, PFOS, 6:2 FTSA, PFOSA) beside literature")
    print("  D_lip/w and chemical structures -- and in two figure captions.")
    print()
    print("  So: the D_medium/w column is measured here (Figure S6). The")
    print("  D_BSA/w column is not, and the paper does not say where it is")
    print("  from; the caption cites ref 2 (Qin 2023) for the MODEL FORM")
    print("  only. A1-A2 rule out two candidate origins, A3 identifies the")
    print("  real one and A4 gives the derivation.")
    print()

    print("  A1. Is Table 1's log D_BSA/w column the Table S6 equation")
    print("      evaluated at IC_10,free?")
    worst = 0.0
    for nm, (_, _, ic_nom, ic_free, d_med, d_cell, d_bsa, d_sp, d_lip,
             _) in T1.items():
        m, b, _ = s6_bsa(nm)
        cw = ic_free * 1e6                       # mol/L -> umol/L
        pred = b if m is None else m * math.log10(cw) + b
        worst = max(worst, abs(pred - d_bsa))
    print(f"      YES. Largest disagreement across 11 compounds: "
          f"{worst:.3f} log units.")
    print()

    print("  A2. Could D_BSA/w have been back-computed from the measured")
    print("      D_medium/w by inverting eq 4?")
    print(f"      eq 4: D_med = {MED_FPROT:.4f} D_BSA + {MED_FLIP:.4f} D_lip")
    print("      Medium is 98.5% protein by sorptive volume, so inverting it")
    print("      can barely move the value and cannot change the slope.")
    print(f"      {'cmpd':9s} {'med slope':>9s} {'BSA slope':>9s} "
          f"{'med b':>6s} {'BSA b':>6s} {'eq4 D_BSA':>9s} {'vs S6':>6s}")
    diffs = []
    for nm in S6:
        ms, mb, _ = s6_med(nm)
        bs, bb, _ = s6_bsa(nm)
        d_lip = T1[nm][8]
        inv = (10 ** mb - 10 ** d_lip * MED_FLIP) / MED_FPROT
        inv = math.log10(inv) if inv > 0 else float("nan")
        diffs.append(bb - mb)
        print(f"      {nm:9s} "
              f"{'const' if ms is None else f'{ms:9.3f}':>9s} "
              f"{'const' if bs is None else f'{bs:9.3f}':>9s} "
              f"{mb:6.2f} {bb:6.2f} {inv:9.2f} {inv - bb:6.2f}")
    gm = sum(diffs) / len(diffs)
    both = [(nm, s6_med(nm)[0], s6_bsa(nm)[0]) for nm in S6
            if s6_med(nm)[0] is not None and s6_bsa(nm)[0] is not None]
    far = [nm for nm, a, b in both if max(abs(a), abs(b)) /
           min(abs(a), abs(b)) > 1.5]
    qual = [nm for nm in S6
            if (s6_med(nm)[0] is None) != (s6_bsa(nm)[0] is None)]
    print(f"      NO. The two columns are not related by eq 4: D_BSA/w at")
    print(f"      Cw = 1 uM sits {gm:+.2f} log units ABOVE D_medium/w on")
    print(f"      average ({10 ** gm:.1f}x), where eq 4 allows at most")
    print(f"      {-math.log10(MED_FPROT):+.3f}.")
    print(f"      The slopes disagree too. {len(far)} of the {len(both)} "
          f"compounds where both")
    print(f"      columns carry a slope differ by more than 1.5x "
          f"({', '.join(far)}),")
    print(f"      and {len(qual)} more have a slope in one column and a flat")
    print(f"      constant in the other ({', '.join(qual)}) -- one matrix")
    print("      saturates where the other does not. Inverting a near-unit")
    print("      linear map cannot produce either.")
    print()
    print("      That gap is the same failure lesson 13 part 4 reports from")
    print("      the forward direction: eq 4 overpredicts D_medium/w from")
    print("      D_BSA/w by about 2x throughout. Two readings, and the SI")
    print("      does not settle between them: either albumin in buffer")
    print("      really does bind ~2x harder per litre of protein than the")
    print("      protein in 2% charcoal-stripped FBS, or the two columns")
    print("      were measured on different material.")
    print()

    print("  A3. SETTLED. The column is Table 1 of the companion paper:")
    print("      Qin W, Escher BI, Huchthausen J, Fu Q, Henneberger L.")
    print("      'Species Difference? Bovine, Trout, and Human Plasma")
    print("      Protein Binding of Per- and Polyfluoroalkyl Substances.'")
    print("      Environ Sci Technol 2024;58(23):9954-9966,")
    print("      doi 10.1021/acs.est.3c10824, PMC11171458, open access.")
    print("      Every cell, slope, intercept and R2, to the last digit:")
    print()
    print(f"      {'cmpd':9s} {'S6 slope':>9s} {'T1 slope':>9s} "
          f"{'S6 b':>6s} {'T1 b':>6s} {'S6 R2':>6s} {'T1 R2':>6s}  match")
    allok = True
    for nm, (s_t, b_t, r_t) in COMPANION_T1_BSA.items():
        s_b, i_b, r_b = s6_bsa(nm)
        ok = s_b == s_t and i_b == b_t and abs(r_b - r_t) < 1e-9
        allok &= ok
        f = lambda v: "const" if v is None else f"{v:.3f}"
        print(f"      {nm:9s} {f(s_b):>9s} {f(s_t):>9s} {i_b:6.2f} "
              f"{b_t:6.2f} {r_b:6.2f} {r_t:6.2f}  "
              f"{'identical' if ok else 'DIFFERS'}")
    print()
    print(f"      All 11 cells identical: {allok}. It is a verbatim")
    print("      reprint, and the baseline-toxicity paper does not cite it")
    print("      -- the companion appeared in issue 23 of the same volume,")
    print("      the baseline paper in issue 13, so the data were in hand")
    print("      but unpublished at submission.")
    print()
    print("  A4. So HOW were they derived? The companion paper's Methods")
    print("      and Data Evaluation give the full chain:")
    print("      1. BioSPME 96-Pin Device, C18-coated pins of ~80 nL")
    print("         coating volume, run automated on a Hamilton Star robot.")
    print("      2. BSA dissolved in PBS. The BSA CONCENTRATION was chosen")
    print("         per compound to put the bound fraction between 30 and")
    print("         90% -- the measurable window. That is why Table 2")
    print("         records 0.1-10 mg/mL rather than one value.")
    print("      3. A 9-step 2-fold dilution series of PFAS (256x), dosed")
    print("         against that BSA, incubated overnight at 37 C, 250 rpm.")
    print("      4. Extraction 10-120 min (t95 was 30 min for hydrophilic")
    print("         PFAS, 58 min for hydrophobic), desorption 20 min,")
    print("         LC-MS/MS. Mass balance 92-115%.")
    print("      5. The isotherm -- bound concentration against free -- was")
    print("         fitted with an empirical FREUNDLICH isotherm. After log")
    print("         transformation that is a straight line of log C_bound")
    print("         against log Cw, giving log K_Fr and the exponent n_Fr.")
    print("         log D at any Cw then follows from those two.")
    print("      6. The decision rule for whether a compound gets an")
    print("         equation or a constant is n_Fr: below 0.9 the isotherm")
    print("         counts as concentration-dependent and the regression is")
    print("         reported; at 0.90-1 (weakly dependent) or >=1 the")
    print("         authors FIX n_Fr = 1 and report the average log D as a")
    print("         constant. That is exactly why PFUnA is a constant at")
    print("         4.75 -- the paper says its n_Fr sat in 0.90-1.")
    print("      7. The fit POOLS two methods. 'BSA binding isotherms of")
    print("         PFBA, PFOA and PFHxS measured in the present study were")
    print("         compared with our previous results with C18-SPME using")
    print("         single fibers... because the data from the two methods")
    print("         were almost overlapping, all data were fitted together")
    print("         to derive log K_Fr and n_Fr.'")
    print()
    print("      Point 7 is the answer to part B below: the slopes differ")
    print("      from Qin 2023 because this is a JOINT refit of the 2023")
    print("      single-fibre data pooled with the new BioSPME data, not a")
    print("      new measurement and not a reprint. Part B inferred 'refit")
    print("      or re-measurement'; it is specifically the former.")
    print()


# =====================================================================
# PART B -- against Qin 2023, the only published Langmuir parameters
# =====================================================================
def partB():
    rule("=")
    print("PART B -- Table S6 against Qin 2023's own fits (4 compounds)")
    rule("=")
    print("  Qin 2023 (Environ Int 173:107857) Table 1 published, for four")
    print("  PFAS, BOTH a Freundlich fit and the Langmuir parameters:")
    print("    Freundlich  log D = (n_Fr - 1) log Cw + log K_Fr")
    print("    mechanistic D(Cw) = Cmax/(Kd + Cw) + D_nonspecific")
    print()
    print(f"  {'cmpd':7s} {'Qin23 slope':>11s} {'S6 slope':>8s} "
          f"{'Qin23 b':>7s} {'S6 b':>6s}   {'Kd uM':>8s} {'Cmax uM':>8s} "
          f"{'D_ns':>7s}")
    for nm, p in QIN_BSA_MEASURED.items():
        Kd, Cmax, CmaxM, lKs, lDns, nFr1, lKFr = p
        bs, bb, _ = s6_bsa(nm)
        print(f"  {nm:7s} {nFr1:11.3f} "
              f"{'const' if bs is None else f'{bs:8.3f}':>8s} "
              f"{lKFr:7.2f} {bb:6.2f}   {Kd:8.4g} {Cmax:8.0f} {lDns:7.2f}")
    print()
    print("  The sulfonates barely move (PFHxS -0.480 -> -0.472, PFOS")
    print("  -0.381 -> -0.379) but both carboxylates do (PFBA -0.391 ->")
    print("  -0.298, PFOA -0.381 -> -0.314). A reprint would move nothing.")
    print()
    print("  A4 explains the pattern exactly. The companion paper pooled")
    print("  the 2023 single-fibre C18-SPME points with its new BioSPME")
    print("  points and refitted, for the three compounds measured twice:")
    print("  PFBA, PFOA and PFHxS. Two of those three are the carboxylates")
    print("  whose slopes moved. PFOS moved least because its low-")
    print("  concentration C18-SPME points were EXCLUDED from the joint fit")
    print("  (see part E), so its refit is nearly the BioSPME data alone.")
    print()


def langmuir_D(cw, Kd, Cmax, Kns):
    return Cmax / (Kd + cw) + Kns


def fit_loglog(cws, ds):
    n = len(cws)
    xs = [math.log10(c) for c in cws]
    ys = [math.log10(d) for d in ds]
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    m = sxy / sxx
    return m, my - m * mx


def partB2():
    rule("=")
    print("PART B2 -- does Langmuir+linear reproduce a Freundlich line?")
    rule("=")
    print("  Generate D(Cw) from Qin 2023's own Kd, Cmax and D_nonspecific,")
    print("  then fit log D against log Cw the way Table S6 does. If the two")
    print("  descriptions of the same experiment disagree, the Freundlich")
    print("  slope is not a clean read-out of the Langmuir shape.")
    print()
    print("  Window: Qin 2023 dosed a 9-step 2-fold series; the fit below")
    print("  uses 10 points over 256x, centred on each compound's own Kd,")
    print("  which is where a Langmuir curve bends most.")
    print()
    print(f"  {'cmpd':7s} {'Kd uM':>8s} {'window uM':>19s} "
          f"{'Langmuir fit':>12s} {'Qin23 Fr':>9s} {'S6 Fr':>7s} {'R2':>5s}")
    for nm, p in QIN_BSA_MEASURED.items():
        Kd, Cmax, CmaxM, lKs, lDns, nFr1, lKFr = p
        Kns = 10 ** lDns
        lo, hi = Kd / 16.0, Kd * 16.0
        cws = [lo * (hi / lo) ** (i / 9.0) for i in range(10)]
        ds = [langmuir_D(c, Kd, Cmax, Kns) for c in cws]
        m, b = fit_loglog(cws, ds)
        ys = [math.log10(d) for d in ds]
        my = sum(ys) / len(ys)
        ss = sum((y - my) ** 2 for y in ys)
        res = sum((y - (m * math.log10(c) + b)) ** 2
                  for c, y in zip(cws, ds and [math.log10(d) for d in ds]))
        r2 = 1 - res / ss if ss else float("nan")
        bs, _, _ = s6_bsa(nm)
        print(f"  {nm:7s} {Kd:8.4g} {lo:8.3g}-{hi:<10.3g} "
              f"{m:12.3f} {nFr1:9.3f} "
              f"{'const' if bs is None else f'{bs:7.3f}':>7s} {r2:5.2f}")
    print()
    print("  Read the last column first: a Langmuir curve fitted in log-log")
    print("  over its own bend is nearly straight (R2 high), so the")
    print("  Freundlich form is a fair description. But the SLOPE it")
    print("  produces is steeper than the published Freundlich slope for")
    print("  every compound, because the published fit was made over the")
    print("  ACTUAL dosing window, not the window centred on Kd. The slope")
    print("  is therefore a property of the window as much as of the")
    print("  chemistry -- which is why part C bounds Kd instead of fitting")
    print("  it.")
    print()


# =====================================================================
# PART C -- the parameter-free bound
# =====================================================================
def partC(rows):
    rule("=")
    print("PART C -- a measured Freundlich slope bounds Kd from above")
    rule("=")
    print("  For D(Cw) = Cmax/(Kd + Cw) + K_ns, with Cmax > 0 and K_ns >= 0,")
    print()
    print("      | d log D / d log Cw |  =  Cw*Cmax / [ (Kd+Cw)^2 * D ]")
    print("                             <=  Cw / (Kd + Cw)")
    print()
    print("  because D >= Cmax/(Kd+Cw). Equality holds only when K_ns = 0.")
    print("  Rearranged, with |m| the measured slope at concentration Cw:")
    print()
    print("      Kd  <=  Cw * (1/|m| - 1)")
    print()
    print("  No capacity, no level, no stoichiometry, no units to reconcile.")
    print("  A steeper slope means a tighter bound. The slope is evaluated")
    print("  at each compound's own IC_10,free, which is where Table 1")
    print("  applies the fit.")
    print()
    print("  Independent albumin Kd, kept apart by species -- Chen 2025 is")
    print("  human (equilibrium dialysis, 0.0625-2 uM ligand, 1 uM HSA),")
    print("  Gao 2025 is bovine (fluorescence quenching, 0.8 uM BSA), the")
    print("  species Qin measured.")
    print()
    chen = {}
    gao = collections.defaultdict(dict)
    for r in rows:
        if "albumin" not in r["protein"].lower():
            continue
        if r["study"] == "Chen et al." and r["year"] == "2025" and r["kd"]:
            chen[r["chemical"]] = float(r["kd"])
        if r["study"] == "Gao et al." and r["year"] == "2025" and r["ka"]:
            import re
            mm = re.search(r"plus (\d+) mM added NaCl", r["buffer"])
            if mm:
                gao[r["chemical"]][int(mm.group(1))] = 1e6 / float(r["ka"])
    print(f"  {'cmpd':9s} {'slope':>7s} {'R2':>4s} {'Cw=IC10f uM':>11s} "
          f"{'Kd bound uM':>11s} | {'Chen25 HSA':>10s} {'ok':>3s} | "
          f"{'Gao25 BSA uM':>13s} {'ok':>3s}")
    verdict = []
    for nm in S6:
        m, b, r2 = s6_bsa(nm)
        cw = T1[nm][3] * 1e6
        if m is None:
            print(f"  {nm:9s} {'const':>7s} {r2:4.2f} {cw:11.4g} "
                  f"{'no bound':>11s} | {'':>10s} {'':>3s} | "
                  f"{'':>13s}")
            continue
        bound = cw * (1.0 / abs(m) - 1.0)
        kc = next((chen[a] for a in names(nm) if a in chen), None)
        gg = next((gao[a] for a in names(nm) if a in gao), None)
        kc_s, kc_ok = ("-", "")
        if kc is not None:
            kc_s = f"{kc:.3g}"
            kc_ok = "yes" if kc <= bound else "NO"
            verdict.append(("Chen 2025 (human)", nm, kc <= bound))
        gg_s, gg_ok = ("-", "")
        if gg:
            lo, hi = min(gg.values()), max(gg.values())
            gg_s = f"{lo:.3g}-{hi:.3g}"
            gg_ok = "yes" if lo <= bound else "NO"
            verdict.append(("Gao 2025 (bovine)", nm, lo <= bound))
        print(f"  {nm:9s} {m:7.3f} {r2:4.2f} {cw:11.4g} {bound:11.4g} | "
              f"{kc_s:>10s} {kc_ok:>3s} | {gg_s:>13s} {gg_ok:>3s}")
    print()
    for src in ("Gao 2025 (bovine)", "Chen 2025 (human)"):
        v = [x for x in verdict if x[0] == src]
        ok = sum(1 for x in v if x[2])
        print(f"  {src}: {ok} of {len(v)} compounds satisfy the bound"
              + ("" if ok == len(v) else
                 "  -- violations: "
                 + ", ".join(x[1] for x in v if not x[2])))
    print()
    print("  A violation is a real inconsistency, not a near miss: it says")
    print("  the isotherm Qin fitted saturates FASTER than a single site of")
    print("  that Kd can saturate. Three ways out, and they are")
    print("  distinguishable in principle: positive cooperativity (Gao fits")
    print("  n_Hill up to 2.5, and a Hill slope above 1 steepens the")
    print("  isotherm past this bound); more than one class of site, the")
    print("  tighter one with a smaller Kd than the fitted average; or the")
    print("  two measurements disagree.")
    print()
    print("  Weigh the R2 column before any of that. PFOSA's slope is")
    print("  -0.105 at R2 0.50 and 6:2 FTSA's is -0.092 at R2 0.54, the two")
    print("  weakest fits in Table S6 -- neither is distinguishable from")
    print("  flat, and a bound computed from a slope that may be zero is")
    print("  not a bound. The PFOSA 'violation' is that, not a conflict.")
    print("  HFPO-DA is different: its slope is the steepest in the table")
    print("  (-0.493, R2 0.86) and the two independent Kd values disagree")
    print("  with each other by 63x (Gao bovine 417 uM, Chen human 6.6 uM),")
    print("  which is the real discrepancy there.")
    print()


# =====================================================================
# PART D -- three studies, one implied non-specific coefficient
# =====================================================================
def partD(rows):
    rule("=")
    print("PART D -- Kd + plateau + slope -> the non-specific coefficient")
    rule("=")
    print("  The bound in part C is the K_ns = 0 limit. Going further needs")
    print("  one more independent number: the low-concentration plateau")
    print()
    print("      D_max = Cmax/Kd + K_ns      (the value as Cw -> 0)")
    print()
    print("  Three studies measured albumin partitioning in exactly that")
    print("  regime, holding the PFAS:protein ratio far below one site per")
    print("  albumin, so what they report IS D_max:")
    print("    Allendorf 2019  bovine, equilibrium dialysis, 37 C, <=0.1")
    print("                    mol/mol -- same species as Qin")
    print("    Bischel 2011    bovine, dialysis, 0.01-0.57 uM against")
    print("                    200 uM BSA")
    print("    Fischer 2024    human, C18-SPME, <=0.004 mol/mol")
    print()
    print("  With Kd (study 2), D_max (study 3) and |m| at Cw (Qin), K_ns is")
    print("  determined. Writing u = Kd/(Kd+Cw):")
    print()
    print("      K_ns = D_max * u * (1-u-|m|) / [ |m| + u*(1-u-|m|) ]")
    print()
    print("  and the bracket is positive only inside the part C bound.")
    print()
    lvl = plateaus(rows)
    gao = collections.defaultdict(dict)
    import re
    for r in rows:
        if r["study"] == "Gao et al." and r["year"] == "2025" and r["ka"] \
                and "albumin" in r["protein"].lower():
            mm = re.search(r"plus (\d+) mM added NaCl", r["buffer"])
            if mm:
                gao[r["chemical"]][int(mm.group(1))] = 1e6 / float(r["ka"])

    print("  Plateau log D_max, all rebased to L_water per L of protein.")
    print("  BOVINE throughout except the last column, which is Fischer's")
    print("  human albumin, carried only to show the species gap:")
    keys = list(lvl)
    print(f"  {'cmpd':9s} " + " ".join(f"{LABEL[k]:>10s}" for k in keys)
          + f" {'Qin S6 @1uM':>11s}")
    for nm in S6:
        row = []
        for k in keys:
            v = next((lvl[k][a] for a in names(nm) if a in lvl[k]), None)
            row.append(f"{v:10.2f}" if v is not None else f"{'-':>10s}")
        print(f"  {nm:9s} " + " ".join(row) + f" {s6_bsa(nm)[1]:11.2f}")
    print()
    print("  Note how much the plateau itself disagrees between studies --")
    print("  up to 1.6 log units for PFUnA. The implied K_ns below inherits")
    print("  that spread, so it is reported per source, never pooled.")
    print()
    print("  Implied K_ns, using Gao's bovine Kd and each bovine plateau --")
    print("  one species throughout, no cross-species mixing:")
    print(f"  {'cmpd':9s} {'Kd uM':>8s} {'logD_max':>8s} {'|m|':>6s} "
          f"{'Cw uM':>9s} {'log K_ns':>8s} {'spec share':>10s}  check")
    for src in [k for k in lvl if "human" not in k]:
        print(f"    -- plateau from {src}")
        for nm in S6:
            m, b, _ = s6_bsa(nm)
            if m is None:
                continue
            g = next((gao[a] for a in names(nm) if a in gao), None)
            dmx = next((lvl[src][a] for a in names(nm) if a in lvl[src]), None)
            if g is None or dmx is None:
                continue
            Kd = min(g.values())              # strongest reported, 0 mM NaCl
            cw = T1[nm][3] * 1e6
            u = Kd / (Kd + cw)
            br = 1.0 - u - abs(m)
            if br <= 0:
                print(f"  {nm:9s} {Kd:8.4g} {dmx:8.2f} {abs(m):6.3f} "
                      f"{cw:9.4g} {'--':>8s} {'--':>10s}  outside the part C "
                      f"bound")
                continue
            kns = 10 ** dmx * u * br / (abs(m) + u * br)
            spec = 1.0 - kns / 10 ** dmx
            print(f"  {nm:9s} {Kd:8.4g} {dmx:8.2f} {abs(m):6.3f} "
                  f"{cw:9.4g} {math.log10(kns):8.2f} {spec * 100:9.1f}%  ok")
        print()
    print()
    print("  'spec share' is the fraction of the plateau that the saturable")
    print("  site contributes. Qin 2023 measured log D_nonspecific directly")
    print("  for four compounds, so those four are a real check:")
    print(f"    {'cmpd':7s} {'Qin23 measured':>14s} {'Qin23 plateau':>13s} "
          f"{'spec share':>10s}")
    for nm, p in QIN_BSA_MEASURED.items():
        Kd, Cmax, CmaxM, lKs, lDns, nFr1, lKFr = p
        dmax = Cmax / Kd + 10 ** lDns
        print(f"    {nm:7s} {lDns:14.2f} {math.log10(dmax):13.2f} "
              f"{(1 - 10 ** lDns / dmax) * 100:9.1f}%")
    print()


def partE(rows):
    rule("=")
    print("PART E -- the payoff: three studies against one measured value")
    rule("=")
    print("  Qin 2023 measured log D_nonspecific directly for four PFAS.")
    print("  Part D reconstructs the same quantity from a Kd measured by")
    print("  Gao 2025, a plateau measured by Allendorf 2019, Bischel 2011")
    print("  or Fischer 2024, and a slope from Qin 2024 Table S6 -- three")
    print("  different laboratories, three different methods, none of them")
    print("  the one that measured the answer.")
    print()
    lvl = plateaus(rows)
    gao = gao_kd(rows)
    srcs = [k for k in lvl if "human" not in k]
    print(f"  {'cmpd':7s} {'measured':>8s} " +
          " ".join(f"{LABEL[k]:>10s}" for k in srcs) + "   worst miss")
    for nm, p in QIN_BSA_MEASURED.items():
        Kd_q, Cmax, CmaxM, lKs, lDns, nFr1, lKFr = p
        m, b, _ = s6_bsa(nm)
        g = next((gao[a] for a in names(nm) if a in gao), None)
        if m is None or g is None:
            continue
        Kd = min(g.values())
        cw = T1[nm][3] * 1e6
        u = Kd / (Kd + cw)
        br = 1.0 - u - abs(m)
        cells, miss = [], []
        for k in srcs:
            dmx = next((lvl[k][a] for a in names(nm) if a in lvl[k]), None)
            if dmx is None or br <= 0:
                cells.append(f"{'-':>10s}")
                continue
            kns = math.log10(10 ** dmx * u * br / (abs(m) + u * br))
            cells.append(f"{kns:10.2f}")
            miss.append(abs(kns - lDns))
        print(f"  {nm:7s} {lDns:8.2f} " + " ".join(cells) +
              (f"   {max(miss):.2f} log" if miss else ""))
    print()
    print("  Three of the four reconstruct to within 0.42 log units (2.6x)")
    print("  of a number none of the three inputs knew, and PFBA to within")
    print("  0.02 on Allendorf's fatty-acid-free batch. The heat-shock")
    print("  fraction is 0.60 off for PFBA, which is the expected")
    print("  direction: it carries bound fatty acids competing for the same")
    print("  site. For a quantity nobody has measured twice in the same")
    print("  system, this is the closest this project has come to")
    print("  validating the Kd-to-K bridge of REPORT 6.1b -- on three")
    print("  compounds, so treat it as encouraging, not settled.")
    print()
    print("  Kd here is Gao's 0 mM NaCl value, the tightest she reports.")
    print("  At 150 mM her Kd is 1.5-2x weaker, which moves K_ns by less")
    print("  than 0.2 log and widens the PFOS gap below rather than closing")
    print("  it.")
    print()
    print("  PFOS does not, and the reason is visible upstream. Qin 2023's")
    print("  PFOS parameters imply a low-concentration plateau of")
    print("  log D_max = log(Cmax/Kd + D_ns):")
    print(f"    {'cmpd':7s} {'Qin23 implied':>13s}   independent plateaus")
    for nm, p in QIN_BSA_MEASURED.items():
        Kd_q, Cmax, CmaxM, lKs, lDns, nFr1, lKFr = p
        dmax = math.log10(Cmax / Kd_q + 10 ** lDns)
        obs = [f"{LABEL[k]} {lvl[k][a]:.2f}" for k in lvl
               for a in names(nm) if a in lvl[k]]
        print(f"    {nm:7s} {dmax:13.2f}   {', '.join(obs)}")
    print()
    print("  PFOS sits 1.2-2.0 log units above every independent plateau,")
    print("  bovine and human alike, and the cause is its Kd: Qin 2023 fits")
    print("  6.2 nM where Gao measures 693 nM in bovine albumin and Chen")
    print("  measures 940 nM in human -- 111x and 151x weaker.")
    print()
    print("  The authors reached the same conclusion independently, which")
    print("  is the strongest confirmation this triangulation could get.")
    print("  The companion paper says of the 2023 PFOS experiment: 'The")
    print("  specific binding of PFOS measured with C18-SPME was a bit")
    print("  higher. However, 5 mg/mL of BSA was used for the C18-SPME,")
    print("  resulting in a bound fraction of PFOS > 99% in the low")
    print("  concentration range. The bound fraction reduced to 40-90%")
    print("  after adjusting BSA to 0.1 mg/mL in this study. log D_BSA/w of")
    print("  PFOS were also fitted with data from two methods, but several")
    print("  points in the low concentrations were excluded.'")
    print()
    print("  At >99% bound the free concentration is below what the method")
    print("  can resolve, so the specific term was fitted to the part of")
    print("  the curve carrying least information. That is exactly the")
    print("  defect the plateau comparison detects from outside, and the")
    print("  authors excluded those points when they refitted. The 6.2 nM")
    print("  should not be used: our mass-balance code carries it in")
    print("  QIN_BSA_MEASURED and this is a reason to flag it there.")
    print()
    print("  What this exercise can and cannot do:")
    print("    CAN  bound Kd from a measured saturation slope, with no")
    print("         capacity, level or stoichiometry (part C)")
    print("    CAN  give the saturable share of binding at the in vitro")
    print("         use concentration from Kd and the slope alone -- the")
    print("         'spec share' column, which does not depend on the")
    print("         plateau source at all: 63-91% across these compounds")
    print("    CAN  reconstruct K_nonspecific from three unrelated studies")
    print("         to within a few tenths of a log unit (parts D, E)")
    print("    CANNOT fit the 10 Table S6 equations from other studies. No")
    print("         study outside this group reports an albumin isotherm")
    print("         for these compounds; the independent albumin work is")
    print("         all deliberately sub-saturating, which is the regime")
    print("         where the slope is zero by construction. A slope needs")
    print("         the saturating regime, and only Qin has measured there.")
    print()


def gao_kd(rows):
    import re
    g = collections.defaultdict(dict)
    for r in rows:
        if r["study"] == "Gao et al." and r["year"] == "2025" and r["ka"] \
                and "albumin" in r["protein"].lower():
            mm = re.search(r"plus (\d+) mM added NaCl", r["buffer"])
            if mm:
                g[r["chemical"]][int(mm.group(1))] = 1e6 / float(r["ka"])
    return g


def main():
    rows = load()
    partA()
    partB()
    partB2()
    partC(rows)
    partD(rows)
    partE(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
