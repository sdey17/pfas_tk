#!/usr/bin/env python3
"""Converting a binding constant into a mass-balance partition coefficient.

A mass-balance model does not want an affinity, it wants a capacity: how many
litres of water's worth of PFAS one litre (or kilogram) of albumin can hold.
The two are the same measurement in different clothes, but the conversion
carries a factor that depends on the molar mass of the protein and on the
number of sites, so "proportional" is not "equal" and a model cannot be tuned
on proportionality alone.  This script does four things:

  1  derives the conversion and states the assumptions it needs
  2  checks it against Qin 2023, whose own D_BSA/w equation IS this conversion
  3  converts every albumin binding constant in db/protein_binding.csv and
     puts them on one axis with the two mass-balance models' own values
  4  shows what an error in D_albumin/water does to predicted C_free and
     C_cell, which is the reason absolute values are needed

Run:  python3 tk_review/scripts/kd_to_partition.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")

M_ALB = 66463.0            # g/mol, Qin's value for BSA (HSA is 66437)
PARTIAL_SPEC_VOL = 0.733   # mL/g -> density of the pure albumin phase
RHO_ALB = 1000.0 / PARTIAL_SPEC_VOL / 1000.0      # kg/L = 1.364
P_ALB_PHASE = RHO_ALB * 1000.0 / M_ALB * 1e6      # uM of albumin inside pure albumin

LINE = "=" * 78


# -------------------------------------------------------------------- 1
def derivation():
    print(LINE)
    print("1.  THE CONVERSION")
    print(LINE)
    print(f"""
  Albumin at total molar concentration P (mol per litre of SOLUTION), with n
  equivalent independent sites of dissociation constant Kd.  Site occupancy is
  Langmuir, so the ligand bound per litre of solution is

      C_bound = n P Cw / (Kd + Cw)                                       (1)

  A partition coefficient is bound ligand per unit of SORBENT over free ligand
  per unit of water.  On a mass basis the sorbent is P*M_alb/1000 kg of protein
  per litre, so

      D_pw = [C_bound / (P M_alb / 1000)] / Cw
           = 1000 n / [M_alb (Kd + Cw)]                                  (2)

  P cancels.  That is the whole trick: the partition coefficient does not
  depend on how much albumin is present, only on how many sites each albumin
  carries and how tightly they hold.

  In the trace limit Cw << Kd, which is where a partition coefficient is a
  constant at all,

      D_pw0 = 1000 n / (M_alb Kd)  =  1000 n Ka / M_alb                  (3)

  With M_alb = {M_ALB:.0f} g/mol:

      D_pw0 [L_water/kg_protein] = {1000.0/M_ALB:.6f} x n x Ka[M^-1]
      log D_pw0                  = log Ka + log n - {math.log10(M_ALB/1000.0):.4f}

  Qin's equation is on a VOLUME basis, L_water per L_albumin.  Multiply by the
  density of the pure albumin phase, {RHO_ALB:.4f} kg/L from the {PARTIAL_SPEC_VOL} mL/g partial
  specific volume:

      D_vw0 [L_water/L_albumin]  = {1000.0*RHO_ALB/M_ALB:.6f} x n x Ka[M^-1]
      log D_vw0                  = log Ka + log n - {math.log10(M_ALB/1000.0/RHO_ALB):.4f}
                                 = log D_pw0 + {math.log10(RHO_ALB):.4f}

  So Ka and D are strictly proportional and the log-log slope is exactly 1.
  The INTERCEPT is what a mass-balance model needs, and it is set by two things
  the affinity measurement usually does not report: the molar mass basis, and n.
  Going from n = 1 to n = 6 moves log D by {math.log10(6):.2f} without changing any
  measurement.  That is the first reason proportionality is not enough.

  THREE ASSUMPTIONS, each of which fails somewhere in this literature:
    (a) trace ligand, Cw << Kd.  Equation (2) shows the exact correction:
        D_pw(Cw) = D_pw0 / (1 + Cw/Kd).  A constant fitted at Cw = 10 Kd is
        {1/(1+10):.3f} of the trace value, i.e. 11x low.
    (b) one class of equivalent sites.  Albumin has several of different
        affinity; a single (n, Kd) pair is a lumped fit whose n is an
        effective capacity, not a count.
    (c) the sorbed ligand is all site-bound.  Non-specific surface
        partitioning adds a term that does NOT saturate, which is exactly why
        Qin carries a separate D_nonspecific.
""")


# -------------------------------------------------------------------- 2
def check_against_qin():
    """Qin 2023's D_BSA/w = Cmax/(Kd + Cw) + D_ns is equation (2) plus (c)."""
    import mbm_qin_vs_armitage as Q
    print(LINE)
    print("2.  VALIDATION: Qin 2023's own equation is this conversion")
    print(LINE)
    print(f"""
  Qin writes  D_BSA/w = Cmax / (Kd + Cw) + D_nonspecific   (their Eq. 16),
  on a volume basis.  Comparing with equation (2) on a volume basis,
  D_vw = 1000 rho n / [M_alb (Kd + Cw)], the identification is

      Cmax [umol per L_albumin] = n x (molar concentration of albumin
                                       inside the pure albumin phase)

  Qin reports Cmax both ways, so their ratio must recover that concentration,
  {P_ALB_PHASE:.0f} uM at {RHO_ALB:.3f} kg/L and {M_ALB:.0f} g/mol.  Checking all four measured
  compounds:
""")
    print(f"  {'cmpd':8s}{'Cmax umol/L':>13s}{'Cmax mol/mol':>14s}{'ratio uM':>11s}"
          f"{'vs expected':>13s}")
    ok = True
    for name, (kd, cmax_v, cmax_m, lks, ldns, nfr1, lkfr) in Q.QIN_BSA_MEASURED.items():
        ratio = cmax_v / cmax_m
        dev = ratio / P_ALB_PHASE
        ok &= abs(dev - 1) < 0.01
        print(f"  {name:8s}{cmax_v:>13.3g}{cmax_m:>14.3f}{ratio:>11.0f}{dev:>12.3f}x")
    print(f"\n  expected {P_ALB_PHASE:.0f} uM -- agreement within 1%: {ok}")
    print("""
  So Qin's "Cmax in mol per mol" column IS n, read straight off the fit:
""")
    print(f"  {'cmpd':8s}{'n (sites/albumin)':>20s}{'Kd uM':>10s}"
          f"{'log D_vw0 site':>16s}{'log D_ns':>10s}{'log D_vw0 total':>17s}")
    for name, (kd, cmax_v, cmax_m, lks, ldns, nfr1, lkfr) in Q.QIN_BSA_MEASURED.items():
        site = cmax_v / kd
        tot = site + 10 ** ldns
        print(f"  {name:8s}{cmax_m:>20.3f}{kd:>10.4g}{math.log10(site):>16.2f}"
              f"{ldns:>10.2f}{math.log10(tot):>17.2f}")
    print("""
  Two things worth noticing. Qin's fitted n is BELOW one site per albumin for
  PFBA (0.45) and PFOS (0.39) and only 1.24 for PFOA, so the Langmuir capacity
  is sub-stoichiometric -- these are lumped effective parameters, not counts.
  And the non-specific term is not a rounding correction: for PFOS it is
  10^4.84, within a factor of 3 of the saturable term.
""")
    return Q


# -------------------------------------------------------------------- 3
def literature_to_partition(Q):
    print(LINE)
    print("3.  EVERY ALBUMIN CONSTANT IN THE FILE, ON ONE AXIS")
    print(LINE)
    print("""
  log D_albumin/water on a MASS basis (L_water per kg_protein), the basis
  Fischer 2024 reports and the one a tissue-composition model wants.  Each
  literature Ka is converted by equation (3) at n = 1 and n = 6 to bracket the
  stoichiometry, which is almost never reported alongside the constant.
""")
    rows = list(csv.DictReader(open(DB)))
    alb = [r for r in rows
           if "serum albumin" in r["protein"]
           and not r["is_computational"].startswith("yes")]

    def logd_from_ka(ka, n):
        return math.log10(1000.0 * n * ka / M_ALB)

    # measured log D, Fischer 2024
    fischer = {}
    for r in alb:
        if (r["study"], r["year"]) == ("Fischer et al.", "2024") \
                and r["ka_units"].startswith("log D"):
            fischer.setdefault(r["chemical"], {})[r["protein_species"]] = float(r["ka"])

    # Qin's own trace-limit D, converted to the mass basis for comparison
    qin_logd = {}
    for name, (kd, cmax_v, cmax_m, lks, ldns, nfr1, lkfr) in Q.QIN_BSA_MEASURED.items():
        tot_vol = cmax_v / kd + 10 ** ldns
        qin_logd[name] = math.log10(tot_vol / RHO_ALB)

    wanted = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
              "PFBS", "PFHxS", "PFOS"]
    print(f"  {'cmpd':8s}{'Fischer HSA':>12s}{'Qin BSA':>9s}{'Armitage':>9s}"
          f"{'  lit. Ka -> log D (n=1)':>26s}{'  n=6':>8s}  sources")
    out = []
    for c in wanted:
        kas = []
        for r in alb:
            if r["chemical"] != c:
                continue
            if "M^-1" not in r["ka_units"] or not r["ka"]:
                continue
            try:
                kas.append((float(r["ka"]), f"{r['study']} {r['year']}"))
            except ValueError:
                pass
        fi = fischer.get(c, {}).get("human")
        qi = qin_logd.get(c)
        arm = armitage_logd(c, Q)
        if kas:
            lo, hi = min(kas)[0], max(kas)[0]
            rng = f"{logd_from_ka(lo,1):.2f} - {logd_from_ka(hi,1):.2f}"
            n6 = f"{logd_from_ka(hi,6):.2f}"
        else:
            rng, n6 = "-", "-"
        print(f"  {c:8s}{(f'{fi:.2f}' if fi else '-'):>12s}"
              f"{(f'{qi:.2f}' if qi else '-'):>9s}"
              f"{(f'{arm:.2f}' if arm is not None else '-'):>9s}"
              f"{rng:>26s}{n6:>8s}  {len(kas)} value(s)")
        out.append((c, fi, qi, arm, kas))
    print("""
  Read the columns against each other, not down them.  Where all four exist the
  spread between the two mass-balance models and the one physiological-ratio
  measurement is the quantity that matters for a mass balance, and it is not
  small: for PFOS, Qin's albumin holds 24x more than Fischer measured.

  A more useful way to use equation (3) is backwards.  Fischer 2024 measured
  log D directly at a physiological ligand:protein ratio, so for each compound
  we can ask what stoichiometry would be needed to reconcile the tightest
  reported binding constant with that measurement:

      n_implied = D_Fischer x M_alb / (1000 x Ka_tightest)

  A value near 1 means the constant and the partition coefficient are telling
  the same story at one site; far from 1 means they are not reconcilable by
  stoichiometry and something else differs.
""")
    print(f"  {'cmpd':8s}{'Fischer log D':>14s}{'tightest Ka':>13s}"
          f"{'n implied':>11s}   source of the tightest Ka")
    for c, fi, qi, arm, kas in out:
        if fi is None or not kas:
            continue
        ka, src = max(kas)
        ni = 10 ** fi * M_ALB / (1000.0 * ka)
        print(f"  {c:8s}{fi:>14.2f}{ka:>13.3g}{ni:>11.2f}   {src}")
    print("""
  PFOA lands at n = 1.0 and PFHpA at 1.1: the tightest measured constant and
  the physiological-ratio partition coefficient agree at one binding site, with
  no fitting.  That is the strongest evidence in this file that equation (3) is
  the right bridge.  Where n_implied runs to tens, the two measurements
  disagree about the chemistry, not about the stoichiometry -- and in those
  cases a mass-balance model built on the binding constant and one built on the
  partition coefficient will not agree however carefully either is converted.
""")
    return out, fischer, qin_logd


def logkow(chem, Q):
    """Armitage drives everything off log Kow; same source as the MBM scripts."""
    for (name, nC, nCF, mw, ko, kc, ksrc) in Q.COMPOUNDS:
        if name == chem:
            return ko if ko is not None else Q.KOW_OPERA(nCF)
    return None


def armitage_logd(chem, Q):
    """Armitage's albumin partition, already L_water per kg_albumin."""
    lk = logkow(chem, Q)
    if lk is None:
        return None
    gKsa = (1.08 * lk - 0.70) if lk < 4.5 else (0.37 * lk + 2.56)
    salt = 10 ** ((0.04 * lk + 0.114) * 0.15)
    return math.log10(10 ** gKsa / salt)     # d_SaW = 0 for acids, so Fi term is 1


# -------------------------------------------------------------------- 4
def sensitivity(Q):
    print(LINE)
    print("4.  WHY ABSOLUTE VALUES ARE NEEDED: error propagation in the balance")
    print(LINE)
    import mbm_published_equations as P
    print("""
  Qin's balance is   Cfree = Cnom Vtot / S,   S = Vw + D_alb V_alb + R
  where R collects every non-albumin sorptive term.  Write a = D_alb V_alb / S
  for the share of total sorptive capacity held by albumin.  If D_alb is wrong
  by a factor f (f < 1 = underestimate), the denominator becomes S[(1-a) + f a]
  and every downstream concentration is wrong by

      Cfree_wrong / Cfree_true = 1 / [(1 - a) + f a]                     (4)

  and C_cell, being proportional to Cfree, is wrong by the same factor.  Two
  consequences.  The error is BOUNDED by 1/(1-a) however badly D_alb is
  underestimated -- a model cannot put more than all the dose in the cells.
  And it depends on a, which varies with the compound and with the serum
  content of the medium.  So a model that is off by a CONSTANT factor in D_alb
  is off by a VARIABLE factor in C_cell, which is why a proportional
  correction fitted on one set of points does not transfer, and why absolute
  values are required.

  Albumin's share of capacity, a, under Qin's own parameters:
""")
    print(f"  {'assay':10s}{'cmpd':8s}{'dose uM':>9s}{'a = alb share':>15s}"
          f"{'10x low':>9s}{'3x low':>8s}{'3x high':>9s}{'10x high':>10s}"
          f"{'  ceiling':>10s}")
    rowsout = []
    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        for nm in ("PFBA", "PFOA", "PFNA"):
            dose = a["doses"][-1]
            q = P.qin_run(nm, dose, V)
            alpha = q["x_alb"]
            def err(f):
                return 1.0 / ((1 - alpha) + f * alpha)
            print(f"  {study:10s}{nm:8s}{dose:>9.1f}{alpha:>15.3f}"
                  f"{err(0.1):>9.2f}{err(1/3):>8.2f}{err(3.0):>9.2f}"
                  f"{err(10.0):>10.2f}{1/(1-alpha):>10.1f}x")
            rowsout.append((study, nm, dose, alpha, err(0.1), 1/(1-alpha)))
    amin = min(r[3] for r in rowsout)
    amax = max(r[3] for r in rowsout)
    bycmp = {}
    for study, nm, dose, alpha, e10, ceil in rowsout:
        bycmp.setdefault(nm, []).append(alpha)
    print(f"""
  a runs {amin:.2f} to {amax:.2f}, and the axis it runs along is the COMPOUND, not the
  assay:""")
    for nm, al in bycmp.items():
        print(f"    {nm:6s} a = {min(al):.3f} - {max(al):.3f} across the three assays"
              f"  (spread {max(al)/min(al):.3f}x)")
    print(f"""
  All three assays use 10% FBS and the same cell line, so they carry almost the
  same albumin capacity and a barely moves between them.  The chain length
  moves it from 0.19 to 0.94.  So a single multiplicative error in D_alb
  produces a C_cell error of {min(r[4] for r in rowsout):.1f}x for PFBA and {max(r[4] for r in rowsout):.1f}x for PFNA in the SAME
  experiment: one fitted correction factor cannot absorb a discrepancy that
  varies with chain length.

  Where the assay axis does bite is serum content, which the wider literature
  varies from 0 to 10%.  Holding the compound fixed at PFOA and sweeping it:
""")
    print(f"  {'FBS %':>7s}{'a':>10s}{'10x low':>10s}{'ceiling':>10s}")
    base = P.ASSAYS["Amstutz"]
    for fbs in (0.0, 0.005, 0.01, 0.02, 0.05, 0.10):
        a2 = dict(base, FBS=fbs)
        V2 = P.volumes(a2)
        q2 = P.qin_run("PFOA", 100.0, V2)
        al = q2["x_alb"]
        e10 = 1.0 / ((1 - al) + 0.1 * al)
        print(f"  {100*fbs:>7.1f}{al:>10.3f}{e10:>10.2f}{1/(1-al) if al < 1 else float('inf'):>10.1f}x")
    print("""
  At 0% serum the albumin parameter is irrelevant and C_cell is set entirely
  by the cell and the plastic.  At 10% it carries almost nine tenths of the
  capacity.  A model validated against serum-free data says nothing about how
  it will behave at 10% FBS, and vice versa.

  The ceiling column is the number to check before spending effort on the
  albumin parameter at all: it is the largest error any albumin
  mis-specification can produce in that cell.  Near 1, the balance is
  insensitive to D_alb and the error is somewhere else; large, and the balance
  is essentially an albumin calculation with a cell attached.

  PRACTICAL CONSEQUENCE.  The three assays in this review are all 10% FBS and
  the compounds span a = 0.19 to 0.94, so fitting a single correction to the
  pooled 29 points will absorb the albumin error for the long chains and
  mis-correct the short ones.  If the albumin parameter is to be improved it
  has to be improved as an absolute value, per compound, on the mass or volume
  basis the model actually uses -- which is what equation (3) and section 3
  provide.
""")
    return rowsout


def main():
    derivation()
    Q = check_against_qin()
    literature_to_partition(Q)
    sensitivity(Q)


if __name__ == "__main__":
    main()
