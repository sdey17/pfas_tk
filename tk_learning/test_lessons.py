"""
Regression tests: do the lessons still produce the numbers they explain?

Run this after changing a lesson, after a library upgrade, and before
walking anyone through the material. It takes a few seconds -- it
deliberately avoids the MCMC in lessons 06-07 and re-derives everything
else from the data.

The tolerances are loose on purpose. These are not unit tests of an
implementation; they are a check that the STORY each lesson tells still
matches the arithmetic. If a number moves by 1% nothing is wrong. If it
moves by 20%, some prose in a lesson has quietly become false.

Usage:  python test_lessons.py
"""
import sys

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

from tk import (load, available, iv_1comp, oral_1comp, iv_2comp,
                half_life, clearance, steady_state)

FAILS = []


def check(label, got, want, tol=0.02, unit=""):
    """Relative comparison, except when `want` is 0."""
    ok = abs(got - want) <= tol * abs(want) if want else abs(got) <= tol
    mark = "ok  " if ok else "FAIL"
    print(f"   [{mark}] {label:<46s} {got:>12.5g} {unit:<9s} "
          f"(expected {want:g}, +/-{100*tol:g}%)")
    if not ok:
        FAILS.append(label)
    return ok


print("Regression tests for the tk_learning lessons\n")

# ----------------------------------------------------------------------
print("tk.py -- the model equations")
check("half_life(0.06)", half_life(0.06), 11.5525, 1e-4, "d")
check("clearance(0.06, 0.2)", clearance(0.06, 0.2), 0.012, 1e-9, "L/kg/d")
check("steady_state(0.1, 0.0108)", steady_state(0.1, 0.0108), 9.259, 1e-3, "mg/L")
check("iv_1comp at one half-life",
      iv_1comp(half_life(0.06), 10, 0.18, 0.06), 10 / 0.18 / 2, 1e-9, "mg/L")
# an oral curve must have the same AUC as the IV curve when F = 1
g = np.linspace(0, 600, 60001)
auc_oral = np.trapezoid(oral_1comp(g, 1.0, 0.25, 0.05, 1.0), g)
check("oral AUC == IV AUC (F=1)", auc_oral, 1.0 / clearance(0.05, 0.25),
      1e-3, "mg/L*d")

# ----------------------------------------------------------------------
print("\ndata files -- schema and self-consistency")
import os, glob
SCHEMA = ["study", "author", "sex", "route", "dose_mgkg", "dose_mg", "bw_kg",
          "time_d", "conc_mgL", "conc_sd", "n_animals", "animal_id", "dataset"]
for f in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       "data", "*.csv"))):
    stem = os.path.basename(f)[:-4]
    if stem.count("_") != 2:
        continue
    chem, sex, species = stem.split("_")
    d = pd.read_csv(f)
    check(f"{stem}: schema", float(list(d.columns) == SCHEMA), 1.0, 1e-9)
    # the sex column must agree with the filename -- the whole point of
    # having it is that a concatenated frame stays unambiguous
    check(f"{stem}: sex column matches name",
          float(sorted(set(d.sex.astype(str))) == [sex]), 1.0, 1e-9)

print("\nlesson 03/04 -- the monkey, one compartment")
mk = load("PFOA_Male_primate")
one = mk[(mk.animal_id == 2054) & (mk.conc_mgL > 0)].sort_values("time_d")
t, C, DOSE = one.time_d.values, one.conc_mgL.values, 10.0
check("n observations, monkey 2054", len(t), 14, 1e-9)

slope, intercept = np.polyfit(t, np.log(C), 1)
check("log-linear on everything: half-life", half_life(-slope), 8.59, 0.01, "d")
check("log-linear on everything: Vd", DOSE / np.exp(intercept), 0.187, 0.01, "L/kg")

m = t >= 7
check("terminal-only half-life (t>=7d)",
      half_life(-np.polyfit(t[m], np.log(C[m]), 1)[0]), 10.16, 0.01, "d")

popt, pcov = curve_fit(
    lambda tt, a, b: np.log(iv_1comp(tt, DOSE, np.exp(a), np.exp(b))),
    t, np.log(C), p0=[np.log(0.2), np.log(0.07)])
Vd, k = np.exp(popt)
check("NLS == log-linear (same objective)", half_life(k), 8.59, 0.01, "d")
se = np.sqrt(np.diag(pcov))
se_CL = np.sqrt(pcov[0, 0] + pcov[1, 1] + 2 * pcov[0, 1])
check("95% interval on k  (x/)", np.exp(1.96 * se[1]), 1.16, 0.02)
check("95% interval on Vd (x/)", np.exp(1.96 * se[0]), 1.42, 0.02)
check("95% interval on CL (x/)", np.exp(1.96 * se_CL), 1.35, 0.02)
check("CL is wider than k, tighter than Vd",
      float(se[1] < se_CL < se[0]), 1.0, 1e-9)

# ----------------------------------------------------------------------
print("\nlesson 05 -- oral dosing")
# A lesson that hardcodes a dataset label breaks SILENTLY when the data
# are regenerated and a label changes (e.g. "6.0 mg/kg" -> "6 mg/kg").
# Guard every literal label in every lesson, not just the ones below.
import re
HERE = os.path.dirname(os.path.abspath(__file__))
labels = set()
for src in sorted(glob.glob(os.path.join(HERE, "[0-9][0-9]_*.py"))):
    labels |= set(re.findall(r'"(\d{6,7}-[^"]*mg/kg-[a-z]+)"', open(src).read()))
known = set()
for name in available():
    known |= set(load(name).dataset.unique())
check("every hardcoded dataset label exists",
      float(bool(labels) and labels <= known), 1.0, 1e-9)
for miss in sorted(labels - known):
    FAILS.append(f"unknown dataset label {miss!r}")

# lesson 05 section B: PFOA, male rats, 6 mg/kg gavage, study 5916078
_g = (load("PFOA_Male_rat")
      .query('dataset == "5916078-6 mg/kg-gavage" and conc_mgL > 0')
      .groupby("time_d", as_index=False).conc_mgL.mean())
check("lesson 05 gavage arm: n time points", len(_g), 10, 1e-9)
_p, _ = curve_fit(
    lambda tt, lv, lk, lka: np.log(oral_1comp(tt, 6.0, np.exp(lv),
                                              np.exp(lk), np.exp(lka))),
    _g.time_d.values, np.log(_g.conc_mgL.values),
    p0=[np.log(0.25), np.log(0.05), np.log(2.0)], maxfev=20000)
_Vd, _k, _ka = np.exp(_p)
check("lesson 05 Vd/F", _Vd, 0.1593, 1e-3, "L/kg")
check("lesson 05 half-life", half_life(_k), 11.57, 1e-3, "d")
check("lesson 05 absorption is fast (flip-flop safe)",
      float(_ka > 10 * _k), 1.0, 1e-9)

# ----------------------------------------------------------------------
print("\nlesson 08/09 -- ODEs and volumes")


def rhs2(t_, A, k10, k12, k21):
    return [-(k10 + k12) * A[0] + k21 * A[1], k12 * A[0] - k21 * A[1]]


tt = np.linspace(0, 60, 13)
num = solve_ivp(rhs2, (0, 60), [10.0, 0.0], t_eval=tt, args=(0.10, 0.50, 0.30),
                rtol=1e-10, atol=1e-13).y[0] / 0.14
ana = iv_2comp(tt, 10.0, 0.14, 0.10, 0.50, 0.30)
check("solve_ivp == analytical 2-compartment",
      float(np.max(np.abs(num / ana - 1))), 0.0, 1e-8)


def predict(t_, a, b, c, d):
    V1_, k10_, k12_, k21_ = np.exp([a, b, c, d])
    s = solve_ivp(rhs2, (0.0, float(t_.max())), [DOSE, 0.0], t_eval=t_,
                  args=(k10_, k12_, k21_), rtol=1e-9, atol=1e-11)
    return np.log(np.maximum(s.y[0] / V1_, 1e-12))


p, _ = curve_fit(predict, t, np.log(C), p0=np.log([0.14, 0.08, 0.4, 0.2]),
                 maxfev=20000)
V1, k10, k12, k21 = np.exp(p)
sm = k10 + k12 + k21
disc = np.sqrt(sm ** 2 - 4 * k10 * k21)
alpha, beta = (sm + disc) / 2, (sm - disc) / 2
CL = V1 * k10
Vss, Vz = V1 * (1 + k12 / k21), CL / beta
A_ = DOSE / V1 * (alpha - k21) / (alpha - beta)
B_ = DOSE / V1 * (k21 - beta) / (alpha - beta)
AUC, AUMC = A_ / alpha + B_ / beta, A_ / alpha ** 2 + B_ / beta ** 2

check("2-compartment terminal half-life", half_life(beta), 15.94, 0.02, "d")
check("V1", V1, 0.1261, 0.02, "L/kg")
check("Vss", Vss, 0.1651, 0.02, "L/kg")
check("Vz", Vz, 0.4052, 0.02, "L/kg")
check("Vextrap", DOSE / B_, 2.124, 0.03, "L/kg")
check("ordering V1 < Vss < Vz < Vextrap",
      float(V1 < Vss < Vz < DOSE / B_), 1.0, 1e-9)
check("NCA Vss == compartmental Vss", DOSE * AUMC / AUC ** 2, Vss, 1e-6, "L/kg")
check("NCA CL == compartmental CL", DOSE / AUC, CL, 1e-6, "L/kg/d")
check("Vz == V1*(1+k12/(k21-beta))", V1 * (1 + k12 / (k21 - beta)), Vz, 1e-9)
check("ln2*Vz/CL == terminal half-life", np.log(2) * Vz / CL, half_life(beta),
      1e-9, "d")

# ----------------------------------------------------------------------
print("\nlesson 08 -- mass balance in the 2-compartment ODE")
def _rhs2(t, state, k10, k12, k21):
    A1, A2 = state
    return [-(k10 + k12) * A1 + k21 * A2,  k12 * A1 - k21 * A2]

_k10, _k12, _k21 = 0.10, 0.50, 0.30
# every INTERNAL transfer term must appear twice with opposite signs, so
# summing the right-hand side leaves only elimination. One call proves it
# for all t; no integration needed.
_d1, _d2 = _rhs2(0.0, [7.0, 3.0], _k10, _k12, _k21)
check("transfer terms cancel: sum(rhs) == -k10*A1",
      _d1 + _d2, -_k10 * 7.0, 1e-12, "mg/kg/d")
_s = solve_ivp(_rhs2, (0, 60), [10.0, 0.0], t_eval=np.linspace(0, 60, 20001),
               args=(_k10, _k12, _k21), rtol=1e-10, atol=1e-13)
_tot = _s.y[0] + _s.y[1]
check("total at t=0 equals the dose", _tot[0], 10.0, 1e-12, "mg/kg")
check("body burden is monotone (bolus, no input)",
      float(np.all(np.diff(_tot) <= 1e-12)), 1.0, 1e-9)
# the general ledger: with an INPUT the total rises, so carry a state that
# accumulates what left, and assert body + eliminated == delivered.
def _rhs_inf(t, state, k, dose, dur):
    A, A_elim = state
    return [(dose / dur if t < dur else 0.0) - k * A, k * A]
_g = np.linspace(0, 60, 2401)
_si = solve_ivp(_rhs_inf, (0, 60), [0.0, 0.0], t_eval=_g, args=(0.06, 10.0, 2.0),
                rtol=1e-11, atol=1e-14, max_step=0.25)
_delivered = np.minimum(_g, 2.0) * (10.0 / 2.0)
check("ledger closes: body + eliminated == delivered",
      float(np.max(np.abs(_si.y[0] + _si.y[1] - _delivered))), 0.0, 1e-6, "mg/kg")

# ----------------------------------------------------------------------
print("\nlesson 08 D2 -- stiffness")
# A stiff system must give the SAME answer on an explicit and an implicit
# solver; only the cost differs. If they disagree, one of them is wrong.
def _rhs_two_speed(t, A, Vb, Vf, Vs, Qf, Qs, CL):
    A_bl, A_f, A_s = A
    C_bl = A_bl / Vb
    out_f, out_s = A_f / Vf, A_s / Vs
    return [Qf * out_f + Qs * out_s - (Qf + Qs) * C_bl - CL * C_bl,
            Qf * (C_bl - out_f), Qs * (C_bl - out_s)]
_p = (0.08, 0.05, 0.60, 8.0, 0.5, 0.02)
_te = np.linspace(0, 120, 50)
_sol = {m: solve_ivp(_rhs_two_speed, (0, 120), [10.0, 0.0, 0.0], t_eval=_te,
                     args=_p, rtol=1e-10, atol=1e-13, method=m)
        for m in ("LSODA", "RK45")}
_tot = {m: v.y.sum(0) for m, v in _sol.items()}
check("stiff system: LSODA and RK45 agree",
      float(np.max(np.abs(_tot["LSODA"] / _tot["RK45"] - 1))), 0.0, 1e-6)
check("implicit solver needs far fewer rhs calls",
      float(_sol["RK45"].nfev > 10 * _sol["LSODA"].nfev), 1.0, 1e-9)
_Vb, _Vf, _Vs, _Qf, _Qs, _CL = _p
_J = np.array([[-(_Qf + _Qs + _CL) / _Vb, _Qf / _Vf, _Qs / _Vs],
               [_Qf / _Vb, -_Qf / _Vf, 0.0],
               [_Qs / _Vb, 0.0, -_Qs / _Vs]])
_ev = np.sort(np.abs(np.linalg.eigvals(_J)))
check("stiffness ratio is large", float(_ev[-1] / _ev[0] > 1e3), 1.0, 1e-9)
# and the lesson 08 section D fit is NOT stiff -- the prose says ratio 3.7
_s = 0.1398 + 0.0156 + 0.0506
_d = np.sqrt(_s ** 2 - 4 * 0.1398 * 0.0506)
check("section D alpha/beta ratio", (_s + _d) / (_s - _d), 3.73, 0.01)

# ----------------------------------------------------------------------
print("\nlesson 10 -- sex and species")
def cl_nca(df, ds):
    g = (df[(df.dataset == ds) & (df.conc_mgL > 0)]
         .groupby("time_d", as_index=False).conc_mgL.mean())
    if len(g) < 4:
        return None
    tt_, cc = g.time_d.values, g.conc_mgL.values
    term = tt_ >= tt_.max() / 3
    if term.sum() < 3:
        return None
    kk = -np.polyfit(tt_[term], np.log(cc[term]), 1)[0]
    if kk <= 0:
        return None
    r = df[df.dataset == ds].iloc[0]
    return dict(study=r.study, dose=r.dose_mgkg, route=r.route,
                CL=r.dose_mgkg / (np.trapezoid(cc, tt_) + cc[-1] / kk))


M = pd.DataFrame([x for x in (cl_nca(load("PFOA_Male_rat"), d)
                              for d in load("PFOA_Male_rat").dataset.unique()) if x])
F = pd.DataFrame([x for x in (cl_nca(load("PFOA_Female_rat"), d)
                              for d in load("PFOA_Female_rat").dataset.unique()) if x])
pairs = M.merge(F, on=["study", "dose", "route"], suffixes=("_M", "_F"))
ratio = np.exp(np.log(pairs.CL_F / pairs.CL_M).mean())
check("matched male/female pairs found", len(pairs), 7, 1e-9)
check("females faster in every pair",
      float((pairs.CL_F > pairs.CL_M).all()), 1.0, 1e-9)
check("geometric mean CL ratio (F/M)", ratio, 28.3, 0.05, "x")

# ----------------------------------------------------------------------
print("\nlesson 11 -- continuous dosing")
DWI, K, VD, DAYS = np.exp(-4.3955), np.log(2) / 3.14, 0.43, 365.25
BAF = DWI * DAYS / (K * VD)
check("serum:water ratio at steady state", BAF, 47.3, 0.01, "x")
css = BAF * 0.1
up = solve_ivp(lambda t_, y: [DWI * DAYS * 0.1 / VD - K * y[0]], (0, 100),
               [0.0], t_eval=[3.14, 15.7], rtol=1e-10).y[0]
check("50% of Css after one half-life", up[0] / css, 0.5, 1e-3)
check("~97% of Css after five half-lives", up[1] / css, 0.969, 0.01)

# ----------------------------------------------------------------------
print("\nlesson 12 -- PBPK")
BW = 0.25
V = dict(liver=.034 * BW, kidney=.007 * BW, blood=.074 * BW, rest=.76 * BW)
QC = 14.0 * BW ** 0.75
Q = dict(liver=.174 * QC, kidney=.141 * QC, rest=.685 * QC)
CLint = 0.02


def rhs_pbpk(t_, A, P):
    Cb = A[0] / V["blood"]
    o = [A[1] / V["liver"] / P, A[2] / V["kidney"] / P, A[3] / V["rest"] / P]
    el = CLint * o[1]
    return [Q["liver"] * o[0] + Q["kidney"] * o[1] + Q["rest"] * o[2]
            - sum(Q.values()) * Cb,
            Q["liver"] * (Cb - o[0]), Q["kidney"] * (Cb - o[1]) - el,
            Q["rest"] * (Cb - o[2])]


tp = np.linspace(0, 200, 8001)
sol = solve_ivp(rhs_pbpk, (0, 200), [10 * BW, 0, 0, 0], t_eval=tp, args=(1.0,),
                rtol=1e-11, atol=1e-14, method="LSODA")
tot = sol.y.sum(axis=0)
check("mass balance holds at t=0", tot[0], 10 * BW, 1e-9, "mg")
mask = (tp > 60) & (tot > tot[0] * 1e-12)
k_num = -np.polyfit(tp[mask], np.log(tot[mask]), 1)[0]
Vtot = sum(V.values())
CL_ws = Q["kidney"] * CLint / (Q["kidney"] + CLint)
check("finite flows: well-stirred formula", k_num, CL_ws / Vtot, 0.005, "/d")
check("well-stirred CL is below CLint", float(CL_ws < CLint), 1.0, 1e-9)

# flows to infinity must recover the one-compartment answer
Qbig = {k: v * 10000 for k, v in Q.items()}
Q_save = dict(Q)
Q.update(Qbig)
sol2 = solve_ivp(rhs_pbpk, (0, 200), [10 * BW, 0, 0, 0], t_eval=tp, args=(1.0,),
                 rtol=1e-11, atol=1e-14, method="LSODA")
Q.update(Q_save)
tot2 = sol2.y.sum(axis=0)
m2 = (tp > 60) & (tot2 > tot2[0] * 1e-12)
check("infinite flows -> k = CLint/Vtot",
      -np.polyfit(tp[m2], np.log(tot2[m2]), 1)[0], CLint / Vtot, 1e-4, "/d")

# ----------------------------------------------------------------------
print()
if FAILS:
    print(f"{len(FAILS)} CHECK(S) FAILED:")
    for f in FAILS:
        print(f"   - {f}")
    print("\nA failure means a number quoted in a lesson no longer matches\n"
          "what the code computes. Fix the prose or fix the code -- but do\n"
          "not ship a lesson that explains a number it does not produce.")
    sys.exit(1)
print("All checks passed. The lessons' numbers match their prose.")
