"""
LESSON 14 -- Where log K = 4.1774 + log n - log Kd comes from, and whether
K and D are the same thing.

Two questions, and the second is not the easy one.

PART A derives the bridge from nothing but a definition and a box of
liquid. No equation is quoted; each one is built and then checked
numerically against a simulated experiment.

PART B asks whether log K and log D are the same quantity in different
units. The short answer is "almost, and the three ways they are not will
each cost you a factor".

Run:  python3 tk_learning/14_kd_and_partition.py
"""
import math

M_ALB = 66463.0          # g/mol, serum albumin
RHO_ALB = 1.0 / 0.733    # kg/L of the pure albumin phase


def rule(ch="=", n=72):
    print(ch * n)


# =====================================================================
# PART A
# =====================================================================
def a1_two_definitions():
    rule()
    print("A1. TWO QUANTITIES THAT DESCRIBE THE SAME TUBE")
    rule()
    print("  Put a chemical in a tube with water and albumin. It ends up")
    print("  partly free in the water and partly stuck to the protein.")
    print("  Two traditions describe that, and they are asking different")
    print("  questions.")
    print()
    print("  THE CHEMIST'S QUESTION: how tightly does one site hold one")
    print("  molecule? Write the binding as a reaction,")
    print()
    print("      P  +  L  <->  PL          Kd = [P][L] / [PL]")
    print()
    print("  Kd has units of concentration, and it has a physical reading:")
    print("  the free concentration at which HALF the sites are occupied.")
    print("  Small Kd = tight binding. It says nothing about how much")
    print("  protein is in the tube.")
    print()
    print("  THE ENGINEER'S QUESTION: if I know the free concentration, how")
    print("  much chemical is sitting in the protein phase? Write a ratio,")
    print()
    print("      K = (amount of chemical per unit of protein) / Cw")
    print()
    print("  K has units of L_water per kg_protein. It is the quantity a")
    print("  mass-balance model needs, because a model multiplies it by how")
    print("  much protein is present to get an amount.")
    print()
    print("  Neither is more fundamental. The bridge between them exists")
    print("  because both are statements about the same tube.")
    print()


def a2_the_bridge():
    rule()
    print("A2. BUILDING THE BRIDGE, ONE LINE AT A TIME")
    rule()
    print("  Let the albumin be at molar concentration P (mol per litre of")
    print("  SOLUTION), with n equivalent sites each of dissociation")
    print("  constant Kd.")
    print()
    print("  (1) What fraction of sites are occupied at free concentration")
    print("      Cw? Rearrange Kd = [P][L]/[PL] for the occupied fraction:")
    print()
    print("          theta = Cw / (Kd + Cw)")
    print()
    print("      Sanity: at Cw = Kd, theta = 1/2. That is the definition.")
    print()
    print("  (2) Each albumin carries n sites, so bound ligand per litre of")
    print("      solution is")
    print()
    print("          C_bound = n * P * Cw / (Kd + Cw)")
    print()
    print("  (3) Now change the denominator of the bookkeeping. K is per")
    print("      unit PROTEIN, not per litre of solution. One litre of")
    print("      solution contains P moles of albumin, which weigh")
    print("      P * M_alb grams, i.e. P * M_alb / 1000 kg. So")
    print()
    print("          amount per kg protein = C_bound / (P * M_alb / 1000)")
    print()
    print("  (4) Divide by Cw to make it a ratio, and watch what happens:")
    print()
    print("          K = [n P Cw/(Kd+Cw)] / [P M_alb/1000] / Cw")
    print("            = 1000 n / [M_alb (Kd + Cw)]")
    print()
    print("      *** P HAS CANCELLED. ***")
    print()
    print("      That is the whole reason this conversion is possible. The")
    print("      partition coefficient does not care how much albumin you")
    print("      put in the tube -- only how many sites each albumin has")
    print("      and how tightly they hold. A3 checks that numerically,")
    print("      because it is the step people disbelieve.")
    print()
    print("  (5) The trace limit. If Cw << Kd the sites are nearly empty,")
    print("      (Kd + Cw) -> Kd, and K stops depending on dose:")
    print()
    print("          K0 = 1000 n / (M_alb Kd)")
    print()
    print("      A partition coefficient is only a CONSTANT in this limit.")
    print("      Outside it, K is a function of concentration.")
    print()


def a3_p_cancels():
    rule()
    print("A3. CHECKING THAT P CANCELS, BY SIMULATING THE TUBE")
    rule()
    print("  Forget the algebra. Build the experiment: dose a fixed total")
    print("  amount into tubes holding different amounts of albumin, solve")
    print("  for equilibrium numerically, and compute K straight from its")
    print("  definition. If step (4) is right, K comes out the same.")
    print()
    Kd, n, L_tot = 2.57, 1.0, 0.10      # uM, sites, uM total ligand
    print(f"  Kd = {Kd} uM, n = {n:.0f}, total ligand {L_tot} uM "
          f"in every tube.")
    print()
    K0 = 1e9 * n / (M_ALB * Kd)
    print(f"  {'albumin uM':>11s}{'Cw uM':>10s}{'% bound':>9s}"
          f"{'K measured':>13s}{'K0/(1+Cw/Kd)':>14s}")
    for P in (0.1, 1.0, 10.0, 100.0):
        cw = free_concentration(L_tot, P, n, Kd)
        bound = n * P * cw / (Kd + cw)              # uM per L solution
        kg_prot = P * 1e-6 * M_ALB / 1000.0         # kg per L solution
        K = (bound / kg_prot) / cw                  # L_water per kg protein
        print(f"  {P:>11.1f}{cw:>10.4f}{100 * bound / L_tot:>8.1f}%"
              f"{K:>13.1f}{K0 / (1 + cw / Kd):>14.1f}")
    print()
    print(f"  Read the last two columns together. K moves a little --")
    print(f"  5643 to 5849, about 3.6% -- and the final column accounts for")
    print(f"  ALL of it: it is the exact expression from step (4), with no")
    print(f"  P in it anywhere. The drift is the Cw term, not an albumin")
    print(f"  term, because adding albumin pulls Cw down and a lower Cw")
    print(f"  means emptier sites.")
    print()
    print(f"  So P really has cancelled. Over a 1000-fold range of albumin,")
    print(f"  with the bound fraction going from 3.6% to 97.5%, the only")
    print(f"  thing that moved K is where the free concentration landed.")
    print(f"  Push the dose up (exercise 1) and that term stops being a")
    print(f"  rounding error.")
    print()
    print(f"  The Cw -> 0 limit of this column is K0 = {K0:,.0f}, which is")
    print(f"  what the formula in A4 returns.")
    print()


def free_concentration(L_tot, P, n, Kd, lo=1e-15, hi=None):
    """Solve L_tot = Cw + n*P*Cw/(Kd+Cw) for Cw, by bisection."""
    hi = L_tot if hi is None else hi
    f = lambda cw: cw + n * P * cw / (Kd + cw) - L_tot
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def a4_the_constant():
    rule()
    print("A4. WHERE 4.1774 COMES FROM -- it is not a fitted number")
    rule()
    print("  Take the trace result and put Kd in micromolar, which is how")
    print("  every paper reports it:")
    print()
    print("      Kd(M) = Kd(uM) x 1e-6")
    print("      K0 = 1000 n / (M_alb x Kd(uM) x 1e-6) = 1e9 n / "
          "(M_alb Kd(uM))")
    print()
    per_kg = 1e9 / M_ALB
    per_L = 1e9 * RHO_ALB / M_ALB
    print(f"      1e9 / M_alb = 1e9 / {M_ALB:.0f} = {per_kg:,.0f}")
    print()
    print(f"      K0 = {per_kg:,.0f} x n / Kd(uM)")
    print(f"      log K0 = {math.log10(per_kg):.4f} + log n - log Kd(uM)")
    print()
    print("  So the constant is not a calibration. READ IT:")
    print()
    print(f"      {per_kg:,.0f} = the number of MICROMOLES OF ALBUMIN in one")
    print("                 KILOGRAM of albumin.")
    print()
    print("  which makes the whole formula say something you can state in")
    print("  words:")
    print()
    print("      K0  =  (how many micromoles of SITES one kg of protein")
    print("              offers)  /  (the Kd those sites have)")
    print()
    print("      capacity divided by affinity. That is all it is.")
    print()
    print("  And on a volume basis, one LITRE of albumin weighs 1.364 kg")
    print("  (from the 0.733 mL/g partial specific volume), so it holds")
    print(f"  {per_L:,.0f} umol, and")
    print()
    print(f"      log D0 = {math.log10(per_L):.4f} + log n - log Kd(uM)")
    print()
    print(f"  The two constants differ by log10(1.364) = "
          f"{math.log10(RHO_ALB):.4f}. That is a DENSITY,")
    print("  not a unit typo, and it is the first of the three traps in")
    print("  part B.")
    print()


def a5_trace_limit():
    rule()
    print("A5. THE ASSUMPTION YOU WILL BREAK FIRST")
    rule()
    print("  K0 is the Cw -> 0 limit. The exact expression from step (4)")
    print("  keeps the dose in it:")
    print()
    print("      K(Cw) = K0 / (1 + Cw/Kd)")
    print()
    Kd, n = 2.57, 1.0
    K0 = 1e9 * n / (M_ALB * Kd)
    print(f"  For PFOA on human albumin, Kd = {Kd} uM, n = {n:.0f}:")
    print(f"  K0 = {K0:,.0f} L_w/kg, i.e. log K0 = {math.log10(K0):.2f}")
    print()
    print(f"  {'Cw uM':>8s}{'Cw/Kd':>8s}{'K(Cw)':>12s}{'log K':>8s}"
          f"{'fraction of K0':>16s}")
    for cw in (0.01, 0.1, 1.0, 2.57, 10.0, 100.0):
        K = K0 / (1 + cw / Kd)
        print(f"  {cw:>8g}{cw / Kd:>8.2f}{K:>12,.0f}"
              f"{math.log10(K):>8.2f}{1 / (1 + cw / Kd):>15.3f}x")
    print()
    print("  A 'partition coefficient' quoted without the concentration it")
    print("  was measured at is therefore incomplete. Two labs reporting")
    print("  3.77 and 3.07 for the same compound may both be right and")
    print("  simply have dosed a decade apart.")
    print()


# =====================================================================
# PART B
# =====================================================================
def b1_same_or_not():
    rule()
    print("B1. ARE log K AND log D THE SAME THING?")
    rule()
    print("  Mostly yes -- and three times no. Both are")
    print()
    print("      (chemical in the sorbent phase) / (chemical free in water)")
    print()
    print("  so they are the same KIND of quantity, and in this project's")
    print("  tables they are interchangeable in that sense. But a number")
    print("  carried from one paper to another has to survive three things")
    print("  that have nothing to do with 'units' in the trivial sense.")
    print()
    print("  TRAP 1 -- THE BASIS. Per kilogram or per litre of protein?")
    print(f"    They differ by log10(1.364) = {math.log10(RHO_ALB):.3f}, "
          f"a factor of 1.36.")
    print("    Fischer 2024 and Allendorf 2019 report per kg. Qin 2024")
    print("    reports per L. Small enough to go unnoticed, big enough to")
    print("    matter when you are arguing about a factor of 2.")
    print()
    print("  TRAP 2 -- WHAT COUNTS AS THE SORBENT. 'The protein' can mean")
    print("    purified albumin, total plasma protein, or protein + lipid")
    print("    together. Qin's D_medium/w and D_cell/w are per litre of")
    print("    PROTEIN PLUS LIPID; his D_BSA/w is per litre of protein")
    print("    alone. Same symbol D, different denominators:")
    print()
    print("      D_BSA/w      [L_w / L_protein]      purified albumin")
    print("      D_plasma/w   [L_w / L_prot+lip]     whole plasma")
    print("      D_cell/w     [L_w / L_prot+lip]     whole cell")
    print("      K_HSA        [L_w / kg_protein]     purified albumin, mass")
    print()
    print("    These are not convertible into one another without the")
    print("    composition of the matrix. That is not a unit conversion,")
    print("    it is a different measurement.")
    print()
    print("  TRAP 3 -- THE ONE THAT IS ACTUALLY A DIFFERENT QUANTITY.")
    print("    In solution chemistry the two letters are not synonyms:")
    print()
    print("      K  partition coefficient of ONE chemical species")
    print("      D  distribution ratio, summing ALL species at a given pH")
    print()
    print("    For a neutral compound they coincide. For an acid they do")
    print("    not, and the gap is set by the pKa:")
    print()
    print("      D = K_neutral x f_neutral  +  K_anion x f_anion")
    print()
    b1_speciation()


def b1_speciation():
    print(f"    {'compound':14s}{'pKa':>6s}{'% anionic at pH 7.4':>21s}"
          f"   does the K/D distinction bite?")
    for name, pka in (("PFOA", -0.2), ("PFOS", -3.3), ("PFOSA", 8.77),
                      ("hexaflumuron", 9.11), ("flubendiamide", 9.03)):
        frac = 1.0 / (1.0 + 10 ** (pka - 7.4))
        verdict = ("no -- fully anionic, one species" if frac > 0.99 else
                   "YES -- two species in play")
        print(f"    {name:14s}{pka:>6.2f}{100 * frac:>20.1f}%   {verdict}")
    print()
    print("    This is why the distinction usually sleeps for PFAS: at pH")
    print("    7.4 the carboxylates and sulfonates are >99.9% anionic, so")
    print("    there is only one species and D = K. It wakes up for the")
    print("    neutral-ish ones. Qin measured PFOSA's pKa at 8.77 for")
    print("    exactly this reason -- 4.1% of it is anionic at pH 7.4, and")
    print("    the anion binds differently from the neutral form, which is")
    print("    what makes its isotherm curve at all.")
    print()
    print("  SO: in this project's albumin tables, log K and log D are the")
    print("  same quantity and you may compare them AFTER fixing the basis")
    print("  (trap 1) and checking the sorbent (trap 2). Across the wider")
    print("  literature, treat the letters as a hint about the author's")
    print("  tradition, never as a definition -- always read what the")
    print("  denominator is.")
    print()


def exercises():
    rule()
    print("EXERCISES")
    rule()
    print("  1. In A3, K came out the same at every albumin concentration.")
    print("     Raise L_tot from 0.10 to 100 uM and re-run. K is no longer")
    print("     constant across the tubes. Which step of A2 broke, and why")
    print("     does the answer now depend on P after all?")
    print()
    print("  2. Qin fits n = 0.39 sites per albumin for PFOS. A site count")
    print("     below one is impossible to take literally. Write down two")
    print("     readings of what a fitted n < 1 means, and say what")
    print("     experiment would distinguish them.")
    print()
    print("  3. You find log K = 4.48 (per kg, human albumin) and Kd = 2.57")
    print("     uM for PFOA in two different papers. Compute the n needed")
    print("     to reconcile them. Then explain why you should NOT publish")
    print("     that n. (tk_review/scripts/kd_k_interconversion.py part 5.)")
    print()
    print("  4. A paper reports D_plasma/w = 4.02 per litre of protein +")
    print("     lipid. Human plasma is 42.25 mL/L protein and 4.46 mL/L")
    print("     lipid. What is the largest value the albumin-only")
    print("     coefficient could have, and what assumption did you make?")
    print()


if __name__ == "__main__":
    print(__doc__.split("Run:")[0].rstrip())
    print()
    a1_two_definitions()
    a2_the_bridge()
    a3_p_cancels()
    a4_the_constant()
    a5_trace_limit()
    b1_same_or_not()
    exercises()
