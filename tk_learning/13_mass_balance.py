"""
LESSON 13 -- Mass-balance models: what they are, and the one job they do.

A mass-balance model (MBM) is bookkeeping for a single well on a plate.
It is not a toxicokinetic model, it does not know what an oral dose is,
and it cannot convert one. Its whole job is:

    "Of the amount I pipetted into this well, how much is actually
     dissolved free in water, and how much is sitting in the cells?"

That is it. Everything else -- plasma concentration, oral dose, half-life
-- belongs to models that run AFTER the MBM has handed over an answer.

WHY IT HAS TO EXIST
-------------------
The number a bioassay paper reports is the NOMINAL concentration:

    C_nom = (moles dosed) / (total volume of the well)

Read that definition again. It is a concentration of nothing in
particular. It assumes the chemical is uniformly dissolved in the whole
well, which is false the moment the medium contains serum albumin or the
plate contains cells. The cells do not respond to C_nom. They respond to
what reaches them.

For PFAS this is not a rounding error. In the assay below, 80% of the
dosed PFOA is stuck to albumin in the medium and 2% is in the cells, so
the free concentration is ~5.6x lower than the nominal one -- and the
discrepancy grows with chain length, because the longer compounds bind
harder. A chain-length trend in a bioassay can therefore be partly an
artefact of the plate, which is exactly what an MBM is for: it separates
"this compound is more potent" from "more of this compound stayed in
solution".

THREE WAYS TO SEE THE SAME EQUATION
-----------------------------------
(1) ACCOUNTING. Mass is conserved. The dose is distributed over phases:

        n_tot = C_w*V_w + C_prot*V_prot + C_lip*V_lip + ...

    and each phase is in equilibrium with the water, C_i = D_i * C_w.
    Substitute, and C_w factors out.

(2) EFFECTIVE VOLUME (the "sink"). Factoring C_w out gives

        n_tot = C_w * [ V_w + D_prot*V_prot + D_lip*V_lip + ... ]
              = C_w * S                                    <- call it S

    S has units of litres. It is the volume of PURE WATER that would
    hold the whole dose at the same free concentration. A well with
    albumin in it behaves like a much bigger glass of water. This is the
    most useful view: every MBM in the literature is a recipe for S.

(3) DILUTION FACTOR. Dividing,

        C_free / C_nom = V_tot / S

    which contains no toxicology at all. It is a property of the plate,
    the medium and the chemical's partitioning -- not of the cells'
    biology. Two labs running the same chemical at the same C_nom in
    different media are dosing different experiments.

WHERE IT SITS IN IVIVE  (the answer to "does it convert to oral dose?")
-----------------------------------------------------------------------
No. It is step 1 of 3, and the three steps use completely different
models and different data:

   step 1   IN VITRO MASS BALANCE          C_nom  ->  C_free (or C_cell)
            Qin 2024 eq 3/6, Armitage 2014, Fischer 2017.
            Inputs: well volumes, medium and cell composition,
            partition coefficients. No physiology.

   step 2   EQUIVALENCE ASSUMPTION         C_free,in vitro == C_free,in vivo
            Not a model -- an assumption, and the weakest link in the
            chain. Usually "free in the well == free in plasma", or
            "cell concentration == liver concentration".

   step 3   REVERSE DOSIMETRY              C_plasma  ->  oral dose
            A TK/PBK model. Inputs: clearance, volume of distribution,
            half-life, absorption. THIS is the step that produces a
            mg/kg/day. Fragki 2023 does this for PFAS.

The MBM converts a CONCENTRATION into another CONCENTRATION inside the
same well. The oral dose appears only in step 3. If you skip step 1 you
hand step 3 a nominal concentration -- an accounting fiction -- and ask a
PBK model to reproduce it as if it were a plasma concentration. The error
is silent and compound-specific.

WHICH EQUATION IS "THE MBM" IN QIN 2024?
----------------------------------------
Qin et al. 2024 (Environ Sci Technol 58:5727) numbers nine equations. The
MBM is eq 3, and eq 6 is its fully predictive version:

    eq 1   C_nom = n_tot / V_tot                     definition, not a model
    eq 2   the mass balance in general form
    eq 3   THE MBM. C_free,medium from C_nom, with D_medium/w, D_cell/w
           and the protein+lipid volumes of medium and cells
    eq 4   D_medium/w  =  D_BSA/w * f_prot  +  D_lip/w * f_lip    (medium)
    eq 5   D_cell/w    =  D_SP/w  * f_prot  +  D_lip/w * f_lip    (cells)
    eq 6   eqs 4 and 5 inserted into eq 3 -> the MBM with no measured
           D_medium/w or D_cell/w required, only surrogates
    eq 7   D_lip/w definition
    eq 8   BASELINE TOXICITY MODEL, a DIFFERENT model: the critical
           membrane burden, 69 mmol/L_lip, converted to a free water
           concentration
    eq 9   eqs 6 + 8 chained -> the nominal concentration a pure baseline
           toxicant would need

Note what eqs 4 and 5 are: little mass-balance models nested inside the
big one. Each is a volume-weighted mean of two surrogate measurements. So
an MBM has TWO failure modes, and they are independent:

    the PARTITIONING can be wrong   (D_BSA/w, D_lip/w, D_SP/w)
    the VOLUMES can be wrong        (f_prot, f_lip, V_cell, V_medium)

Section 6.7 of ../tk_review/report/MBM_REVIEW.md is about exactly this
split. Part 5 below is a closure test that isolates the volume term.

And note what Qin 2024 does NOT contain: any oral dose at all. It stops
at a toxic ratio. The MBM is complete without ever leaving the plate.

Run:  python3 tk_learning/13_mass_balance.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tk_review", "scripts"))
# Table 1 of Qin 2024, as ingested into the binding database. Importing it
# rather than retyping it keeps the lesson and the database from drifting.
from add_qin2024_table1 import T1

# --- the actual assay, from Qin 2024 SI ---------------------------------
# Text S2: the free-concentration and cytotoxicity experiment behind Table 1
# was run in a 96-well plate, 25000 HEK293H cells per well, 220 uL total.
# (The 384-well numbers in Tables S7/S8 are the separate HTS campaign.)
V_WELL = 220e-6          # L of medium
N_CELLS = 25_000

# Table S7, PPARgamma row: 98% OptiMEM + 2% csFBS
MED_F_PROT = 0.94e-3     # L protein per L medium
MED_F_LIP = 1.47e-5      # L lipid   per L medium
MED_F_W = 0.999          # L water   per L medium

# Table S8, HEK293H row
CELL_V_PER_1E6 = 2.87e-3 * 1e-3   # L per 1e6 cells
CELL_F_PROT = 33.9e-3    # L protein per L cell
CELL_F_LIP = 5.3e-3      # L lipid   per L cell
CELL_F_W = 0.961         # L water   per L cell


def rule(ch="-", n=72):
    print(ch * n)


# =======================================================================
# PART 1. The simplest possible well: water only.
# =======================================================================
def part1():
    rule("=")
    print("PART 1 -- a well with nothing in it but water")
    rule("=")
    n_tot = 1e-9                      # mol dosed
    V = V_WELL
    C_nom = n_tot / V
    # The only phase is water, so D = 1 and S = V.
    S = V
    C_free = n_tot / S
    print(f"  dosed            {n_tot:.3g} mol into {V * 1e6:.0f} uL")
    print(f"  C_nom            {C_nom * 1e6:.3f} umol/L")
    print(f"  sink S           {S:.3e} L   (= the water volume)")
    print(f"  C_free           {C_free * 1e6:.3f} umol/L")
    print(f"  C_free / C_nom   {C_free / C_nom:.3f}")
    print()
    print("  With one non-sorbing phase, C_nom IS the free concentration.")
    print("  This is the only case where the nominal number means anything")
    print("  physical, and it is not the case any bioassay runs in.")
    print()


# =======================================================================
# PART 2. Add one sorbing phase. The sink grows; nothing else changes.
# =======================================================================
def part2():
    rule("=")
    print("PART 2 -- add albumin: one phase, one partition coefficient")
    rule("=")
    n_tot = 1e-9
    V = V_WELL
    V_prot = V * MED_F_PROT
    V_w = V - V_prot
    print(f"  medium           {V * 1e6:.0f} uL, of which protein "
          f"{V_prot * 1e9:.1f} nL ({MED_F_PROT * 100:.3f}% by volume)")
    print(f"  {'log D_BSA/w':>12s}  {'S (L)':>10s}  {'S/V_w':>7s}  "
          f"{'C_free/C_nom':>12s}  {'% of dose on protein':>21s}")
    for logD in (0, 1, 2, 3, 3.98, 5):
        D = 10 ** logD
        S = V_w + D * V_prot
        frac_prot = D * V_prot / S
        print(f"  {logD:12.2f}  {S:10.3e}  {S / V_w:7.2f}  "
              f"{V / S:12.3f}  {frac_prot * 100:20.1f}%")
    print()
    print("  Read the third column: that is how many times bigger the well")
    print("  looks than its own water volume. 0.094% of the well is protein,")
    print("  but at log D 3.98 (PFOA on BSA) that sliver holds 90% of the")
    print("  dose, because 10^3.98 = 9550 beats 1/0.00094 = 1064.")
    print()
    print("  RULE OF THUMB: a phase matters when D_i > 1 / f_i, i.e. when")
    print("  its partition coefficient beats the reciprocal of its volume")
    print("  fraction. For 2% serum medium that threshold is log D ~ 3.0.")
    print()


# =======================================================================
# PART 3. The real thing: Qin 2024 eq 3, medium and cells.
# =======================================================================
def sink(d_med_w, d_cell_w, n_cells=N_CELLS, v_well=V_WELL):
    """Qin 2024 eq 3 denominator, broken into its four terms.

    S = D_medium/w * V_medium,prot+lip + V_medium,w
      + D_cell/w   * V_cell,prot+lip   + V_cell,w

    Note which volume each D is multiplied by: the SORPTIVE volume only.
    Water is in the sink at D = 1, but it is not inside either D. That is
    Qin's convention and it differs from Fischer 2017, which folds water
    into D_cell/w. The two are not interchangeable, and mixing them up
    double-counts or drops the water.
    """
    v_cell = n_cells / 1e6 * CELL_V_PER_1E6
    terms = {
        "medium protein+lipid": d_med_w * v_well * (MED_F_PROT + MED_F_LIP),
        "medium water":         v_well * MED_F_W,
        "cell protein+lipid":   d_cell_w * v_cell * (CELL_F_PROT + CELL_F_LIP),
        "cell water":           v_cell * CELL_F_W,
    }
    return sum(terms.values()), terms, v_well + v_cell


def part3():
    rule("=")
    print("PART 3 -- Qin 2024 eq 3, for PFOA, with the measured D values")
    rule("=")
    _, _, ic_nom, ic_free, d_med, d_cell, d_bsa, d_sp, d_lip, _ = T1["PFOA"]
    S, terms, V_tot = sink(10 ** d_med, 10 ** d_cell)
    print(f"  measured log D_medium/w {d_med:.2f}   "
          f"log D_cell/w {d_cell:.2f}")
    print(f"  V_tot {V_tot * 1e6:.2f} uL     sink S {S:.4e} L"
          f"     S / V_tot = {S / V_tot:.2f}")
    print()
    print("  where the dose actually is:")
    for k, v in terms.items():
        print(f"    {k:22s} {v:.4e} L   {v / S * 100:6.2f}%")
    print()
    print(f"  C_free / C_nom  = V_tot / S = {V_tot / S:.4f}")
    print(f"  C_nom  / C_free =            {S / V_tot:.2f}x")
    print()
    print(f"  measured IC_10,nom  {ic_nom:.3g} mol/L")
    print(f"  measured IC_10,free {ic_free:.3g} mol/L")
    print(f"  measured ratio      {ic_nom / ic_free:.2f}x")
    print()
    print("  The cells hold 2% of the dose. The medium protein holds 80%.")
    print("  This is the single most important fact about PFAS bioassays:")
    print("  the serum in the medium, not the cells, is the main sink, so")
    print("  the serum percentage is an exposure variable.")
    print()


# =======================================================================
# PART 4. Where the two D values come from: eqs 4 and 5.
# =======================================================================
def part4():
    rule("=")
    print("PART 4 -- eqs 4 and 5: each D is itself a little mass balance")
    rule("=")
    print("  D_i/w = D_protein/w * f_prot + D_lip/w * f_lip,  f summing to 1")
    print("  over the SORPTIVE volume only (protein + lipid, no water).")
    print()
    fm = MED_F_PROT / (MED_F_PROT + MED_F_LIP)
    fc = CELL_F_PROT / (CELL_F_PROT + CELL_F_LIP)
    print(f"  medium: protein is {fm * 100:.1f}% of the sorptive volume, "
          f"lipid {100 - fm * 100:.1f}%")
    print(f"  cell:   protein is {fc * 100:.1f}% of the sorptive volume, "
          f"lipid {100 - fc * 100:.1f}%")
    print()
    print(f"  {'cmpd':9s} {'logD_BSA':>8s} {'logD_lip':>8s} "
          f"{'pred med':>8s} {'meas med':>8s} {'x':>5s}   "
          f"{'logD_SP':>7s} {'pred cell':>9s} {'meas cell':>9s} {'x':>5s}")
    for nm, (_, _, _, _, d_med, d_cell, d_bsa, d_sp, d_lip, _) in T1.items():
        p_med = math.log10(10 ** d_bsa * fm + 10 ** d_lip * (1 - fm))
        p_cell = math.log10(10 ** d_sp * fc + 10 ** d_lip * (1 - fc))
        print(f"  {nm:9s} {d_bsa:8.2f} {d_lip:8.2f} "
              f"{p_med:8.2f} {d_med:8.2f} {10 ** abs(p_med - d_med):5.2f}   "
              f"{d_sp:7.2f} {p_cell:9.2f} {d_cell:9.2f} "
              f"{10 ** abs(p_cell - d_cell):5.2f}")
    print()
    print("  A weighted mean of two numbers must lie BETWEEN them. So if a")
    print("  measured D sits outside [D_protein, D_lip], no choice of")
    print("  composition can reproduce it -- the surrogates are wrong, and")
    print("  no amount of retuning volume fractions will help. Count those:")
    bad = []
    for nm, (_, _, _, _, d_med, d_cell, d_bsa, d_sp, d_lip, _) in T1.items():
        if not min(d_sp, d_lip) <= d_cell <= max(d_sp, d_lip):
            bad.append(f"{nm} (cell {d_cell:.2f} outside "
                       f"[{min(d_sp, d_lip):.2f}, {max(d_sp, d_lip):.2f}])")
    for b in bad:
        print(f"    {b}")
    print(f"  -> {len(bad)} of {len(T1)} compounds, for the CELL term.")
    print("  This is a composition-independent failure. It is the finding")
    print("  section 6.7 of the review turns on.")
    print()


# =======================================================================
# PART 5. Closure test: do the volumes work, given measured D values?
# =======================================================================
def part5():
    rule("=")
    print("PART 5 -- closure test on the VOLUMES, all 11 compounds")
    rule("=")
    print("  Feed eq 3 the MEASURED D_medium/w and D_cell/w and compare the")
    print("  predicted C_nom/C_free with the measured IC_10,nom/IC_10,free.")
    print("  This is NOT an independent validation -- D_medium/w was itself")
    print("  measured as bound/free in this same system, so the comparison")
    print("  is close to a consistency check. That is precisely what makes")
    print("  it useful: the partitioning cancels, so what is left under")
    print("  test is the VOLUME and COMPOSITION term.")
    print()
    print(f"  {'cmpd':9s} {'nCF':>3s} {'meas ratio':>10s} {'pred ratio':>10s} "
          f"{'pred/meas':>9s}  {'cells hold':>10s}")
    res = []
    for nm, (ncf, _, ic_nom, ic_free, d_med, d_cell, _, _, _, _) in T1.items():
        S, terms, V_tot = sink(10 ** d_med, 10 ** d_cell)
        pred = S / V_tot
        meas = ic_nom / ic_free
        res.append(pred / meas)
        print(f"  {nm:9s} {ncf:3d} {meas:10.2f} {pred:10.2f} "
              f"{pred / meas:9.2f}  "
              f"{terms['cell protein+lipid'] / S * 100:9.2f}%")
    gm = math.exp(sum(math.log(x) for x in res) / len(res))
    print()
    print(f"  geometric mean pred/meas {gm:.2f}x, "
          f"range {min(res):.2f}-{max(res):.2f}x")
    print()
    print("  Within a factor of ~2 across three orders of magnitude of")
    print("  partitioning. So the volumes and compositions are roughly")
    print("  right -- and the trouble in PART 4 is the surrogate D values,")
    print("  not the bookkeeping.")
    print()
    print("  But look at the direction. The predicted sink is too SMALL in")
    print("  9 of 11 cases, by about 1.3x on average. A one-sided miss is")
    print("  not noise; it says eq 3 is missing a sink, or that the D fed")
    print("  in is not the right average over the dose range. Two")
    print("  candidates, not separable from Table 1 alone:")
    print("    - the plate. 96-well polystyrene sorbs PFAS, Armitage 2014")
    print("      carries an explicit plastic term, Qin eq 3 has none.")
    print("    - D_medium/w here is the Freundlich fit evaluated at one")
    print("      point, IC_10,free; the dose-response spans 10 dilutions")
    print("      and D rises as concentration falls.")
    print()


# =======================================================================
# PART 6. What the MBM does not do.
# =======================================================================
def part6():
    rule("=")
    print("PART 6 -- the thing an MBM cannot tell you")
    rule("=")
    _, _, ic_nom, ic_free, d_med, d_cell, _, _, d_lip, _ = T1["PFOA"]
    S, _, V_tot = sink(10 ** d_med, 10 ** d_cell)
    print(f"  PFOA IC_10,nom   {ic_nom * 1e6:9.2f} umol/L   "
          f"<- what the bioassay paper reports")
    print(f"  PFOA IC_10,free  {ic_free * 1e6:9.2f} umol/L   "
          f"<- what the MBM delivers ({S / V_tot:.1f}x lower)")
    print(f"       C_cell      {ic_free * 10 ** d_cell * 1e6:9.2f} umol/L   "
          f"per L of cell protein+lipid")
    print(f"       C_membrane  {ic_free * 10 ** d_lip * 1e6:9.2f} umol/L   "
          f"per L of membrane lipid (eq 7)")
    print()
    print("  Four concentrations, all of the same experiment, spanning "
          f"{10 ** d_cell:.0f}x.")
    print("  The MBM's output is whichever of these the next step needs.")
    print("  It has said nothing about:")
    print("    - how much of an oral dose reaches plasma (absorption)")
    print("    - how fast it leaves (clearance, half-life)")
    print("    - how it is distributed in a body (Vd, PBK compartments)")
    print("    - whether a free concentration in a well is the right")
    print("      counterpart of a free concentration in plasma")
    print("  Those are steps 2 and 3. The MBM is finished at the plate edge.")
    print()
    rule("=")
    print("EXERCISES")
    rule("=")
    print("  1. In PART 2, find the log D at which albumin and water hold")
    print("     the dose equally. Do it algebraically first, then check it")
    print("     by editing the loop. What does the answer depend on?")
    print()
    print("  2. Qin's assay uses 2% serum. Re-run PART 3 with the 10% FBS")
    print("     generic medium (f_prot 3.0e-3, f_lip 7e-5). How much does")
    print("     the PFOA IC_10,free move? Which published IC values are")
    print("     therefore not comparable across labs?")
    print()
    print("  3. PART 5 sends the measured D values in and gets the measured")
    print("     ratio back out. Write down what would have to be true for")
    print("     that test to FAIL. (Hint: only one group of inputs is")
    print("     actually being tested.)")
    print()
    print("  4. Fischer 2017 defines D_cell/w to INCLUDE cell water, Qin")
    print("     excludes it. Modify sink() to Fischer's convention without")
    print("     changing the physical answer. What has to change in the")
    print("     tabulated D_cell/w for the two to agree? (Qin's HEK293H")
    print(f"     cell is {CELL_F_W * 100:.1f}% water by volume.)")
    print()


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
