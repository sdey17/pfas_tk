#!/usr/bin/env python3
"""Fit a serum mass-balance model to Fischer 2024 -- one species, one lab.

WHY THIS DATASET.  A mass-balance model needs an ABSOLUTE partition
coefficient for every sink in the matrix and an independent measurement of
what the matrix as a whole does.  Section 6.1b showed why most of our 955-row
binding table cannot supply that: the f_unbound values are 94% Ryu tissue
homogenates, the purified-albumin f_unbound are circular (the same experiment
that gave the K), and mixing a human K with a mouse f_unbound is not a test
of the model, it is a test of the species difference.

Fischer et al. 2024 (Environ Sci Technol 58:1055-1063, doi 10.1021/acs.est.3c07415)
is the one dataset that closes the loop inside a single species:

    components   log D_HSA/w   (Lw per kg HSA)      28 PFAS
                 log D_glob/w  (Lw per kg globulin) 27 PFAS
    whole matrix log D_serum/w (Lw per L serum)     27 PFAS

measured by the same C18-fibre depletion method, in the same lab, at a
PFAS:protein molar ratio <= 0.004 -- the trace regime the body is in.  The
components and the whole matrix are separate experiments, so predicting one
from the other is a real out-of-sample test, not a rearrangement.

THE MODEL (Fischer eqs 2-3, and the standard serum balance).  Serum is water
plus albumin plus globulin.  Each sink multiplies the free concentration:

    D_serum/w = VF_w + D_HSA/w . c_alb + D_glob/w . c_glob          (Lw/Lserum)
    f_unbound = VF_w / D_serum/w

with c in kg protein per L serum and VF_w the water volume fraction.  Note
the protein DENSITY cancels: Fischer writes D_vol . VF, and VF = c/rho, so
D_mass . rho . c/rho = D_mass . c.  Density survives only in VF_w.

WHAT IS FITTED.  Fischer fixes the composition at the NHANES 2017-18 average
and checks the prediction (their Figure 3B).  That is test 1 below.  Tests 2-4
go further: let the composition float and ask what serum the data actually
behave like, whether globulin earns its term, and whether the residual has
chain-length structure a two-protein model cannot carry.

Run:  python3 tk_review/scripts/mbm_fischer_serum.py
"""
import csv
import math
import os

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(os.path.dirname(HERE), "db")
FIG = os.path.join(os.path.dirname(HERE), "figures")

# ---------------------------------------------------------------- composition
# Fischer section 3.4: "mean concentrations were 41 g/L albumin and 31 g/L
# globulins" in NHANES 2017-18 (n = 1919).  Ranges 21-52 and 20-55 g/L.
C_ALB_NHANES = 0.041          # kg/L serum
C_GLOB_NHANES = 0.031         # kg/L serum
RHO_ALB = 1.36                # kg/L, Fischer section 2.7
RHO_GLOB = 1.20               # kg/L, Fischer section 2.7
NHANES_ALB_RANGE = (0.021, 0.052)
NHANES_GLOB_RANGE = (0.020, 0.055)

# Fischer section 3.3, read off the text for four compounds.
FU_REPORTED = {"PFOA": 6.1e-4, "PFHxS": 4.1e-4, "PFHpS": 3.5e-4, "PFOS": 4.2e-4}

PROT = {"human serum albumin (HSA)": "HSA",
        "gamma-globulin": "glob",
        "whole human serum (surrogate)": "serum",
        "bovine serum albumin (BSA)": "BSA"}

CLASS = {"PFPeA": "carboxylate", "PFHxA": "carboxylate", "PFHpA": "carboxylate",
         "PFOA": "carboxylate", "PFNA": "carboxylate", "PFDA": "carboxylate",
         "PFUnDA": "carboxylate", "PFDoDA": "carboxylate",
         "PFTriDA": "carboxylate", "PFTeDA": "carboxylate",
         "PFBS": "sulfonate", "PFPeS": "sulfonate", "PFHxS": "sulfonate",
         "PFHpS": "sulfonate", "PFOS": "sulfonate", "PFNS": "sulfonate",
         "PFDS": "sulfonate",
         "FHxSA": "sulfonamide", "FOSA": "sulfonamide", "FOSAA": "sulfonamide",
         "N-MeFOSAA": "sulfonamide", "N-EtFOSAA": "sulfonamide",
         "FDSA": "sulfonamide",
         "5:3 FTCA": "fluorotelomer", "6:2 FtS": "fluorotelomer",
         "7:3 FTCA": "fluorotelomer", "8:2 FtS": "fluorotelomer",
         "10:2 FtS": "fluorotelomer", "ADONA": "ether"}

# Chen 2025 names -> Fischer names, for the optional membrane sensitivity.
CHEN_ALIAS = {"PFTrDA": "PFTriDA", "PFOSA": "FOSA", "6:2 FTS": "6:2 FtS",
              "8:2 FTS": "8:2 FtS", "10:2 FTS": "10:2 FtS"}


# ---------------------------------------------------------------------- data
def load():
    """Fischer 2024 rows, pivoted to one record per compound."""
    rows = [r for r in csv.DictReader(open(os.path.join(DB, "protein_binding.csv")))
            if r["study"].startswith("Fischer") and r["year"] == "2024"]
    d = {}
    for r in rows:
        key = PROT.get(r["protein"])
        if key is None:
            continue
        rec = d.setdefault(r["chemical"], {"cmpd": r["chemical"],
                                           "nCF": int(r["n_perfluorinated_carbons"]),
                                           "cls": CLASS.get(r["chemical"], "?")})
        rec[key] = float(r["ka"])
        if r["unbound_fraction"]:
            rec["fu_rep"] = float(r["unbound_fraction"])
    return d


def load_membrane():
    """Chen 2025 log K_membrane/water, for the optional third-sink test."""
    out = {}
    for r in csv.DictReader(open(os.path.join(DB, "protein_binding.csv"))):
        if "Chen" not in r["study"] or "log K_MW" not in r["ka_units"]:
            continue
        name = CHEN_ALIAS.get(r["chemical"], r["chemical"])
        out[name] = float(r["ka"])
    return out


# --------------------------------------------------------------------- model
def vf_water(c_alb, c_glob):
    return 1.0 - c_alb / RHO_ALB - c_glob / RHO_GLOB


def predict(recs, c_alb, c_glob, c_lip=0.0, logKmw=None):
    """log10 D_serum/w predicted from the component coefficients."""
    vw = vf_water(c_alb, c_glob)
    out = []
    for r in recs:
        d = vw + 10 ** r["HSA"] * c_alb + 10 ** r["glob"] * c_glob
        if c_lip and logKmw is not None and r["cmpd"] in logKmw:
            d += 10 ** logKmw[r["cmpd"]] * c_lip
        out.append(math.log10(d))
    return np.array(out)


def fit(recs, free_glob=True, c_lip_free=False, logKmw=None):
    """Least squares on the log10 residual, composition free and positive."""
    meas = np.array([r["serum"] for r in recs])

    def resid(p):
        c_alb = 10 ** p[0]
        c_glob = 10 ** p[1] if free_glob else 0.0
        c_lip = 10 ** p[2] if c_lip_free else 0.0
        return predict(recs, c_alb, c_glob, c_lip, logKmw) - meas

    p0 = [math.log10(C_ALB_NHANES), math.log10(C_GLOB_NHANES), math.log10(0.0023)]
    n = 3 if c_lip_free else 2
    lo = [-4.0, -6.0, -6.0][:n]
    hi = [math.log10(0.30)] * n
    r = least_squares(resid, p0[:n], bounds=(lo, hi))
    c_alb = 10 ** r.x[0]
    c_glob = 10 ** r.x[1] if free_glob else 0.0
    c_lip = 10 ** r.x[2] if c_lip_free else 0.0
    return c_alb, c_glob, c_lip, r


def report(pred, meas, label, npar):
    res = pred - meas
    rmse = float(np.sqrt(np.mean(res ** 2)))
    bias = float(np.mean(res))
    mae = float(np.mean(np.abs(res)))
    ss_res = float(np.sum(res ** 2))
    ss_tot = float(np.sum((meas - meas.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot
    n = len(meas)
    aic = n * math.log(ss_res / n) + 2 * npar
    print(f"  {label:42s}{10**bias:>8.2f}{10**mae:>8.2f}{rmse:>8.3f}"
          f"{r2:>8.3f}{aic:>9.1f}{n:>5d}")
    return dict(bias=bias, mae=mae, rmse=rmse, r2=r2, aic=aic, n=n, res=res)


# ---------------------------------------------------------------------- main
def main():
    data = load()
    recs = sorted((r for r in data.values()
                   if "HSA" in r and "glob" in r and "serum" in r),
                  key=lambda r: (r["cls"], r["nCF"]))
    meas = np.array([r["serum"] for r in recs])
    print(__doc__.split("Run:")[0].rstrip())
    print("=" * 84)
    print(f"DATA.  {len(recs)} PFAS have all three of D_HSA/w, D_glob/w and D_serum/w.")
    print("=" * 84)
    dropped = [c for c, r in data.items()
               if not ("HSA" in r and "glob" in r and "serum" in r)]
    print(f"\n  dropped for an N.D. in one column: {', '.join(sorted(dropped))}")
    print(f"\n  {'cmpd':11s}{'nCF':>4s}{'class':>14s}{'D_HSA':>8s}{'D_glob':>8s}"
          f"{'D_serum':>9s}{'glob/HSA':>10s}")
    for r in recs:
        print(f"  {r['cmpd']:11s}{r['nCF']:>4d}{r['cls']:>14s}{r['HSA']:>8.2f}"
              f"{r['glob']:>8.2f}{r['serum']:>9.2f}"
              f"{10**(r['glob']-r['HSA']):>10.3f}")

    # ------------------------------------------------------------- test 1
    print("\n" + "=" * 84)
    print("TEST 1.  Fischer's own prediction: composition fixed at the NHANES average")
    print("=" * 84)
    print(f"""
  c_alb = {C_ALB_NHANES*1000:.0f} g/L, c_glob = {C_GLOB_NHANES*1000:.0f} g/L, VF_w = {vf_water(C_ALB_NHANES, C_GLOB_NHANES):.4f}
  Nothing is fitted.  This is their Figure 3B, reproduced from Table S2.
""")
    print(f"  {'model':42s}{'bias':>8s}{'|err|':>8s}{'RMSE':>8s}{'R2':>8s}"
          f"{'AIC':>9s}{'n':>5s}")
    p_fixed = predict(recs, C_ALB_NHANES, C_GLOB_NHANES)
    s_fixed = report(p_fixed, meas, "albumin + globulin, NHANES average", 0)
    p_albonly = predict(recs, C_ALB_NHANES, 0.0)
    s_albonly = report(p_albonly, meas, "albumin only, NHANES average", 0)

    print(f"\n  {'cmpd':11s}{'meas':>8s}{'pred':>8s}{'pred/meas':>11s}"
          f"{'alb share':>11s}{'glob share':>12s}{'water share':>13s}")
    vw = vf_water(C_ALB_NHANES, C_GLOB_NHANES)
    for r, p in zip(recs, p_fixed):
        a = 10 ** r["HSA"] * C_ALB_NHANES
        g = 10 ** r["glob"] * C_GLOB_NHANES
        tot = vw + a + g
        print(f"  {r['cmpd']:11s}{r['serum']:>8.2f}{p:>8.2f}"
              f"{10**(p-r['serum']):>11.2f}{a/tot:>11.1%}{g/tot:>12.1%}"
              f"{vw/tot:>13.2%}")

    # ------------------------------------------------------------ test 1b
    print("\n" + "=" * 84)
    print("TEST 1b.  Reconciling with Fischer's own reported fit, and the basis question")
    print("=" * 84)
    print("""
  Fischer report "slope = 0.97, 95%-CI = 0.94-1.01, R2 = 0.99" for Figure 3B.
  That is a NO-INTERCEPT regression of measured on predicted in log space.  On
  data where log D runs 2.0-3.9 and no point is anywhere near zero, forcing the
  line through the origin makes R2 almost a foregone conclusion: it is measured
  against sum(y^2), not against the variance of y.  Both parameterisations below
  reproduce it.  The honest predictive numbers are the ones in test 1.

  THE BASIS QUESTION.  Fischer's section 2.7 says D was "transformed to
  volumetric units assuming protein densities of 1.36 and 1.20 kg/L", and
  section 3.4 then states VF_HSA = 4.1%, VF_glob = 3.1%, VF_w = 92.8% -- which
  are the g/L concentrations read straight off as volume percentages.  Taken
  literally that multiplies albumin's term by rho twice over.  Reading it as a
  mass balance instead (D in Lw/kg times c in kg/L) makes the density cancel and
  gives VF_w = 94.4%.  The two differ by 1.36x on the albumin term.
""")
    print(f"  {'parameterisation':34s}{'c_alb':>8s}{'c_glob':>8s}"
          f"{'slope0':>9s}{'R2(no-ic)':>11s}{'bias':>8s}{'RMSE':>8s}{'R2 vs 1:1':>11s}")
    for lab, calb, cglob in (
            ("mass basis, density cancels", C_ALB_NHANES, C_GLOB_NHANES),
            ("VF read literally, x density", C_ALB_NHANES * RHO_ALB,
             C_GLOB_NHANES * RHO_GLOB)):
        pp = predict(recs, calb, cglob)
        b0 = float((pp * meas).sum() / (pp * pp).sum())
        sse = float(((meas - b0 * pp) ** 2).sum())
        r2_0 = 1 - sse / float((meas ** 2).sum())
        res0 = pp - meas
        r2_11 = 1 - float((res0 ** 2).sum()) / float(((meas - meas.mean()) ** 2).sum())
        print(f"  {lab:34s}{calb*1000:>8.1f}{cglob*1000:>8.1f}{b0:>9.3f}"
              f"{r2_0:>11.4f}{10**res0.mean():>8.2f}"
              f"{float(np.sqrt((res0**2).mean())):>8.3f}{r2_11:>11.3f}")
    print("""
  Fischer's published 0.97 (CI 0.94-1.01) sits between the two, so the
  comparison cannot settle which they ran; the 1.36x is inside the model's own
  1.55x error.  What does lean one way is test 2: letting the albumin term float
  returns ~40 g/L, the mass concentration, not the 56 g/L the literal reading
  implies.""")

    # ------------------------------------------------------------- test 2
    print("\n" + "=" * 84)
    print("TEST 2.  Let the composition float: what serum do the data behave like?")
    print("=" * 84)
    ca, cg, _cl, _r = fit(recs, free_glob=True)
    p_fit = predict(recs, ca, cg)
    ca1, _cg1, _cl1, _r1 = fit(recs, free_glob=False)
    p_fit1 = predict(recs, ca1, 0.0)
    print(f"""
  Fitted (2 parameters): c_alb = {ca*1000:5.1f} g/L, c_glob = {cg*1000:5.1f} g/L
  Fitted (1 parameter):  c_alb = {ca1*1000:5.1f} g/L, globulin term removed
  NHANES average:        c_alb = {C_ALB_NHANES*1000:5.1f} g/L, c_glob = {C_GLOB_NHANES*1000:5.1f} g/L
  Serum product used:    total protein 40-90 g/L (Sigma H4522, Fischer section 2.1)
""")
    print(f"  {'model':42s}{'bias':>8s}{'|err|':>8s}{'RMSE':>8s}{'R2':>8s}"
          f"{'AIC':>9s}{'n':>5s}")
    report(p_fixed, meas, "fixed NHANES (0 fitted)", 0)
    s_fit1 = report(p_fit1, meas, "albumin only, c_alb fitted (1)", 1)
    s_fit = report(p_fit, meas, "albumin + globulin, both fitted (2)", 2)

    # ------------------------------------------------------------- test 3
    print("\n" + "=" * 84)
    print("TEST 3.  Does globulin earn its term?")
    print("=" * 84)
    hi_glob = [r for r in recs if r["glob"] - r["HSA"] > -1.0]
    lo_glob = [r for r in recs if r["glob"] - r["HSA"] <= -1.0]
    print(f"""
  Fischer's own rule: globulin matters once D_glob/D_HSA > 0.1.  That splits
  the set into {len(hi_glob)} compounds above the line and {len(lo_glob)} below it.
  If the rule is right, dropping globulin should hurt the first group and
  leave the second alone.
""")
    print(f"  {'subset':42s}{'bias':>8s}{'|err|':>8s}{'RMSE':>8s}{'R2':>8s}"
          f"{'AIC':>9s}{'n':>5s}")
    for nm, sub in (("D_glob/D_HSA > 0.1", hi_glob), ("D_glob/D_HSA <= 0.1", lo_glob)):
        if len(sub) < 3:
            continue
        m = np.array([r["serum"] for r in sub])
        report(predict(sub, ca, cg), m, f"{nm}: with globulin", 2)
        report(predict(sub, ca, 0.0), m, f"{nm}: without globulin", 2)

    # ------------------------------------------------------------- test 4
    print("\n" + "=" * 84)
    print("TEST 4.  Residual structure: what the two-protein model cannot carry")
    print("=" * 84)
    res = p_fit - meas
    print(f"\n  {'class':16s}{'n':>4s}{'bias':>9s}{'|err|':>9s}  compounds")
    for cls in ("carboxylate", "sulfonate", "sulfonamide", "fluorotelomer", "ether"):
        idx = [i for i, r in enumerate(recs) if r["cls"] == cls]
        if not idx:
            continue
        sub = res[idx]
        names = ", ".join(recs[i]["cmpd"] for i in idx)
        print(f"  {cls:16s}{len(idx):>4d}{10**sub.mean():>9.2f}"
              f"{10**np.abs(sub).mean():>9.2f}  {names}")
    ncf = np.array([r["nCF"] for r in recs], float)
    sl, ic = np.polyfit(ncf, res, 1)
    rho = float(np.corrcoef(ncf, res)[0, 1])
    print(f"""
  Residual vs chain length: slope {sl:+.3f} log units per CF2, Pearson r = {rho:+.2f}.
  Residual vs log D_glob/D_HSA: r = """
          f"{float(np.corrcoef([r['glob']-r['HSA'] for r in recs], res)[0,1]):+.2f}")
    worst = sorted(zip(recs, res), key=lambda t: -abs(t[1]))[:6]
    print(f"\n  {'worst misses':16s}{'nCF':>5s}{'pred/meas':>11s}")
    for r, e in worst:
        print(f"  {r['cmpd']:16s}{r['nCF']:>5d}{10**e:>11.2f}")
    print("""
  INDEPENDENT CHECK.  Fischer's own section 3.3 names exactly four compounds
  that the binding equation misses: "Three fluorotelomers (5:3 FTCA, 7:3 FTCA,
  6:2 FtS) were underpredicted ... In contrast, PFHpS was overpredicted."  The
  four largest residuals above are PFHpS over-predicted and those same three
  fluorotelomers under-predicted, and nothing else comes close.  The
  implementation here therefore reproduces theirs, outlier for outlier.

  Their reading of the fluorotelomers is a missing sink (alpha-globulins,
  phospholipids); test 7 puts the phospholipid half of that to the test.  Their
  reading of PFHpS is fatty-acid competition.  A second possibility sits in
  their own data: PFHpS has the largest HSA-BSA gap in the set (5.49 vs 4.88),
  and Fischer flag PFHpS as one of four compounds where D_BSA/w is a poor proxy
  for D_HSA/w.  Whichever albumin number is the outlier, the 5.1x miss is driven
  by a single measurement, not by the structure of the balance.""")

    # ------------------------------------------------------------- test 5
    print("\n" + "=" * 84)
    print("TEST 5.  Unbound fraction")
    print("=" * 84)
    print(f"""
  f_unbound = VF_w / D_serum/w.  Three routes to it:
    (a) from the MEASURED D_serum/w                -- Fischer's Figure 3A
    (b) from the two-protein model at NHANES average -- Fischer's eq 3
    (c) the four values stated in Fischer section 3.3
""")
    print(f"  {'cmpd':11s}{'(a) measured':>14s}{'(b) model':>12s}"
          f"{'(c) stated':>12s}{'b/a':>7s}{'c/a':>7s}")
    for r in recs:
        if r["cmpd"] not in FU_REPORTED:
            continue
        fa = vw / 10 ** r["serum"]
        fb = vf_water(ca, cg) / 10 ** predict([r], ca, cg)[0]
        fc = FU_REPORTED[r["cmpd"]]
        print(f"  {r['cmpd']:11s}{fa:>14.2e}{fb:>12.2e}{fc:>12.2e}"
              f"{fb/fa:>7.2f}{fc/fa:>7.2f}")
    print("""
  Column (c) does not reconcile with column (a).  f_unbound = VF_w/D_serum/w
  applied to Table S2 gives 0.14% for PFOA where section 3.3 states 0.061%, and
  0.08% for PFOS where it states 0.042%; PFHxS and PFHpS agree to within 20%.
  The stated RANGE (0.01-1.2%) does reproduce exactly from Table S2 -- PFDS at
  the bottom, PFHxA at the top -- so the formula is right and the four
  individual numbers in the running text are the part that does not tie out.
  Everything below uses column (a), computed from Table S2, not the text.""")

    # ------------------------------------------------------------- test 6
    print("\n" + "=" * 84)
    print("TEST 6.  Interindividual variability from NHANES serum composition")
    print("=" * 84)
    print(f"""
  Hold each compound's D_HSA/w and D_glob/w fixed and move only the serum.
  NHANES 2017-18: albumin {NHANES_ALB_RANGE[0]*1000:.0f}-{NHANES_ALB_RANGE[1]*1000:.0f} g/L (CV 8.5%), globulin {NHANES_GLOB_RANGE[0]*1000:.0f}-{NHANES_GLOB_RANGE[1]*1000:.0f} g/L (CV 14%).
  How much f_unbound moves depends entirely on whether the two travel together.

    (A) INDEPENDENT.  Sweep the four corners of the box.  Albumin and globulin
        are free to be high at the same time, so every compound sees the full
        swing of its dominant sink.
    (B) ANTI-CORRELATED.  Fischer's Figure S12B/C: individuals with low globulin
        have high albumin.  Impose it at its limit -- total protein fixed at
        {(C_ALB_NHANES+C_GLOB_NHANES)*1000:.0f} g/L, albumin swept {NHANES_ALB_RANGE[0]*1000:.0f}-{NHANES_ALB_RANGE[1]*1000:.0f} g/L, globulin taking the rest.

  (A) is the upper bound, (B) the lower; the real NHANES joint distribution is
  between them.  The interesting part is that they differ BY COMPOUND.
""")
    print(f"  {'cmpd':11s}{'nCF':>4s}{'glob/HSA':>10s}{'f_u avg':>11s}"
          f"{'(A) fold':>10s}{'(B) fold':>10s}{'(B) f_u min':>13s}{'(B) f_u max':>13s}")
    varrows = []
    tot_prot = C_ALB_NHANES + C_GLOB_NHANES
    for r in recs:
        indep = []
        for calb in NHANES_ALB_RANGE:
            for cglob in NHANES_GLOB_RANGE:
                indep.append(vf_water(calb, cglob) / 10 ** predict([r], calb, cglob)[0])
        anti = []
        for calb in NHANES_ALB_RANGE:
            cglob = max(tot_prot - calb, 0.0)
            anti.append(vf_water(calb, cglob) / 10 ** predict([r], calb, cglob)[0])
        favg = vf_water(C_ALB_NHANES, C_GLOB_NHANES) / 10 ** predict(
            [r], C_ALB_NHANES, C_GLOB_NHANES)[0]
        print(f"  {r['cmpd']:11s}{r['nCF']:>4d}{10**(r['glob']-r['HSA']):>10.3f}"
              f"{favg:>11.2e}{max(indep)/min(indep):>10.2f}"
              f"{max(anti)/min(anti):>10.2f}{min(anti):>13.2e}{max(anti):>13.2e}")
        varrows.append((r, min(indep), favg, max(indep),
                        min(anti), max(anti)))
    byanti = sorted(varrows, key=lambda t: t[5] / t[4])
    print(f"""
  Under (A) every compound lands at 2.6-2.8x, because the corner sweep lets both
  proteins rise together and the balance is nearly linear in total protein.  That
  is an artefact of the sweep, not a result.

  Under (B) the fold range runs {byanti[0][5]/byanti[0][4]:.2f}x ({byanti[0][0]['cmpd']}) to {byanti[-1][5]/byanti[-1][4]:.2f}x ({byanti[-1][0]['cmpd']}), and the ordering is
  mechanistic: a compound that binds ONLY albumin sees the full albumin swing,
  while a compound that binds albumin and globulin about equally is buffered,
  because what it loses on one it gains on the other.

  Fischer report a maximum of 2.5x, for PFHxS, and 1.6x for 10:2 FtS (Figure 4);
  p99/p1 of 1.6x for PFHxS and 1.3x for PFTriDA.  Model (B) gives PFHxS
  {dict((t[0]['cmpd'], t[5]/t[4]) for t in varrows).get('PFHxS', float('nan')):.2f}x -- the albumin-only limit {NHANES_ALB_RANGE[1]/NHANES_ALB_RANGE[0]:.2f}x almost exactly -- and PFTriDA
  {dict((t[0]['cmpd'], t[5]/t[4]) for t in varrows).get('PFTriDA', float('nan')):.2f}x.  PFHxS lands on their number exactly, because an albumin-only
  compound is already at the albumin-only limit and no joint distribution can
  push it past 2.5x.  10:2 FtS does not: model (B) gives 1.12x against their
  1.6x, and (A) gives 2.74x, so the real NHANES correlation is partial, nearer
  the anti-correlated end.  The anti-correlation is therefore the MECHANISM
  that buffers the long-chain PFAS, but its perfect-correlation limit
  over-buffers them; the true spread sits between the two columns.

  FOR TK.  This is a floor on interindividual variability in the free
  concentration that is NOT pharmacokinetic at all -- same dose, same clearance,
  same everything, 2.5x apart in free PFHxS because one person's albumin is
  21 g/L and another's is 52.  Clearance scales with the free fraction, so this
  propagates straight into half-life.
""")

    # ------------------------------------------------------- test 7 (optional)
    kmw = load_membrane()
    shared = [r for r in recs if r["cmpd"] in kmw]
    print("\n" + "=" * 84)
    print("TEST 7.  Sensitivity: does a phospholipid sink change anything?")
    print("=" * 84)
    print(f"""
  CROSSES STUDIES, not species.  Chen et al. 2025 measured log K_membrane/water
  on solid-supported phosphatidylcholine bilayers for {len(shared)} of these {len(recs)} compounds.
  Human serum carries roughly 2-3 g/L phospholipid, inside lipoproteins.
  Chen's published basis is not fully pinned down in our table (the SI equation
  defining it is not in hand), so this is a sensitivity, not a result.
""")
    if len(shared) >= 10:
        m3 = np.array([r["serum"] for r in shared])
        print(f"  {'model':42s}{'bias':>8s}{'|err|':>8s}{'RMSE':>8s}{'R2':>8s}"
              f"{'AIC':>9s}{'n':>5s}")
        report(predict(shared, ca, cg), m3, "2 proteins, composition from test 2", 2)
        report(predict(shared, ca, cg, 0.0023, kmw), m3,
               "+ phospholipid at 2.3 g/L (not fitted)", 2)
        ca3, cg3, cl3, _ = fit(shared, free_glob=True, c_lip_free=True, logKmw=kmw)
        report(predict(shared, ca3, cg3, cl3, kmw), m3,
               "+ phospholipid, all three fitted", 3)
        print(f"\n  fitted: c_alb = {ca3*1000:.1f} g/L, c_glob = {cg3*1000:.1f} g/L, "
              f"c_lipid = {cl3*1000:.2f} g/L")

    # ------------------------------------------------------------------ write
    dst = os.path.join(DB, "mbm_fischer_serum.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["compound", "nCF", "class", "logD_HSA_w_Lw_per_kg",
                    "logD_glob_w_Lw_per_kg", "logD_serum_w_meas_Lw_per_L",
                    "logD_serum_w_pred_NHANES", "logD_serum_w_pred_fitted",
                    "ratio_pred_meas_NHANES", "ratio_pred_meas_fitted",
                    "albumin_share_of_D", "globulin_share_of_D",
                    "f_unbound_from_measured_D", "f_unbound_model_avg",
                    "f_unbound_indep_min", "f_unbound_indep_max",
                    "fold_range_independent", "f_unbound_anticorr_min",
                    "f_unbound_anticorr_max", "fold_range_anticorrelated"])
        for r, pf, pt, (rr, fmin, favg, fmax, amin, amax) in zip(
                recs, p_fixed, p_fit, varrows):
            a = 10 ** r["HSA"] * C_ALB_NHANES
            g = 10 ** r["glob"] * C_GLOB_NHANES
            tot = vw + a + g
            w.writerow([r["cmpd"], r["nCF"], r["cls"], r["HSA"], r["glob"],
                        r["serum"], round(pf, 3), round(pt, 3),
                        round(10 ** (pf - r["serum"]), 3),
                        round(10 ** (pt - r["serum"]), 3),
                        round(a / tot, 4), round(g / tot, 4),
                        f"{vw / 10 ** r['serum']:.3e}", f"{favg:.3e}",
                        f"{fmin:.3e}", f"{fmax:.3e}", round(fmax / fmin, 2),
                        f"{amin:.3e}", f"{amax:.3e}", round(amax / amin, 2)])
    print(f"\n  wrote {dst} ({len(recs)} rows)")
    return dict(recs=recs, meas=meas, p_fixed=p_fixed, p_fit=p_fit,
                ca=ca, cg=cg, res=res, varrows=varrows,
                s_fixed=s_fixed, s_fit=s_fit, s_fit1=s_fit1, s_albonly=s_albonly)


if __name__ == "__main__":
    main()
