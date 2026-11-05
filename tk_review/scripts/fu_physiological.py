#!/usr/bin/env python3
"""Pick one albumin coefficient per compound, then compute fu.

WHY A RULE RATHER THAN A MEDIAN.  fu = 1/(1 + K c_p), and at
physiological albumin every PFAS here is over 99 % bound, so fu is a
small difference between large numbers: an error in log K lands on it
almost undiluted, and the 1.2-log between-laboratory spread in this
literature is a 16-fold spread in fu.  Taking a median across the
literature would average a measured coefficient together with one
extrapolated from a titration that never reached its own Kd, and throw
away the three columns this project built precisely to tell them
apart.  So the rule RANKS candidates and reports which one won and
why, alongside the full admissible spread.

ADMISSIBILITY -- hard filters, applied first.

  species        must match the target. Human, bovine and rat albumin
                 are different proteins and are never pooled.
  not held out   Wu 2024's HSA constants and Han 2021's quoted values
                 are excluded for documented defects; see
                 scripts/excluded_sources.py.
  has a value    a log K must exist.

PREFERENCE -- the surviving candidates are ranked, best first.

  1  a coefficient the PAPER measured beats one converted here from a
     Kd. A conversion needs a stoichiometry the paper usually does not
     report and a molar mass, and inherits the Kd's error one for one.
  2  a Freundlich-type coefficient evaluated at an operating point is
     demoted below a trace-limit one: it is concentration-dependent by
     construction, so it is not the same quantity.
  3  a STATED basis beats an unstated one. Mass and volume differ by
     0.135 log and an unstated basis is that much ambiguity.
  4  among converted values, one whose Kd sits INSIDE the tested
     concentration range beats one outside it, and identifiable beats
     marginal beats not identifiable.
  5  among converted values, a stoichiometry the paper counted or
     fixed beats one the reader assumed, which beats a fluorescence
     quenching exponent.

WHAT IS REPORTED.  The winner, the reason it won, and the range of fu
implied by EVERY admissible candidate. The second is the honest
answer; the first is what a model would have to use.

Run:  python3 tk_review/scripts/fu_physiological.py
Writes db/fu_physiological.csv
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "albumin_lfabp_reported_conc.xlsx")
OUT = os.path.join(HERE, "..", "db", "fu_physiological.csv")

M_ALB = 66463.0
VBAR = 0.733e-3                       # L/g, serum albumin
ALBUMIN_G_L = [("35 g/L", 35.0), ("40 g/L", 40.0), ("50 g/L", 50.0)]
COMPOUNDS = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]
HELD_OUT = {"Wu 2024", "Han 2021"}


def fu(logK, g_per_L, basis):
    c_p = g_per_L * VBAR if basis == "volume" else g_per_L / 1000.0
    return 1.0 / (1.0 + (10.0 ** logK) * c_p)


def basis_of(text):
    t = (text or "").lower()
    if t.startswith("volume"):
        return "volume"
    if t.startswith("mass"):
        return "mass"
    return "unstated"


def rank(row):
    """Lower sorts better. Each element is one preference, in order."""
    route = str(row["route"] or "")
    reported = route.startswith("measured")
    freundlich = "equation at Cw" in route
    basis = basis_of(row["basis_text"])
    rng = str(row["range_text"] or "")
    ident = ("identifiable" if "- identifiable" in rng else
             "marginal" if "marginal" in rng else
             "not" if "not identifiable" in rng else "none")
    geo = ("inside" if rng.startswith("inside") else
           "edge" if rng.startswith("edge") else
           "far" if rng.startswith("far") else
           "outside" if rng.startswith("outside") else "none")
    n_how = (0 if "n counted" in route else
             1 if "fixed by the paper" in route else
             2 if "assumed by the reader" in route else
             3 if "quenching-plot exponent" in route else 2)
    return (
        0 if (reported and not freundlich) else 1 if reported else 2,
        0 if basis != "unstated" else 1,
        {"inside": 0, "edge": 1, "none": 2, "outside": 3, "far": 4}[geo],
        {"identifiable": 0, "none": 1, "marginal": 2, "not": 3}[ident],
        n_how,
    )


def why(row):
    route = str(row["route"] or "")
    bits = []
    if route.startswith("measured") and "equation at Cw" not in route:
        bits.append("the paper measured the coefficient")
    elif route.startswith("measured"):
        bits.append("measured, but a Freundlich fit at an operating "
                    "point")
    else:
        bits.append("converted here from a Kd")
    bits.append(f"{basis_of(row['basis_text'])} basis")
    rng = str(row["range_text"] or "")
    if rng and rng != "None":
        bits.append(rng)
    return "; ".join(bits)


def load():
    from openpyxl import load_workbook
    ws = load_workbook(SRC)["Coefficients"]
    ix = {h: i for i, h in enumerate(c.value for c in ws[1])}
    rows = []
    for r in ws.iter_rows(min_row=2, values_only=True):
        if r[ix["PFAS"]] not in COMPOUNDS:
            continue
        if "albumin" not in str(r[ix["Protein"]]).lower():
            continue
        if r[ix["log K (L water / kg protein)"]] is None:
            continue
        src = str(r[ix["Source"]]).split(" Cw")[0].strip()
        if src in HELD_OUT:
            continue
        rows.append(dict(
            pfas=r[ix["PFAS"]], species=r[ix["Species"]], source=src,
            logK=float(r[ix["log K (L water / kg protein)"]]),
            route=r[ix["Route"]],
            basis_text=r[ix["log K basis"]],
            range_text=r[ix["Kd within tested range?"]]))
    return rows


def main():
    rows = load()
    out = []
    for species in ("human", "bovine"):
        here = [r for r in rows if r["species"] == species]
        if not here:
            continue
        print("=" * 78)
        print(f"{species.upper()} SERUM ALBUMIN")
        print("=" * 78)
        for c in COMPOUNDS:
            cand = sorted((r for r in here if r["pfas"] == c), key=rank)
            if not cand:
                print(f"\n  {c}: NO admissible candidate")
                continue
            best = cand[0]
            b = basis_of(best["basis_text"])
            # Candidates the rule cannot separate. Picking one of
            # these is a coin flip, and a coin flip between
            # Allendorf's two BSA batches is 0.62 log, so it has to
            # be said rather than hidden by sort order.
            tied = [r for r in cand if rank(r) == rank(best)]
            print(f"\n  {c}   {len(cand)} admissible")
            print(f"    chosen  {best['source']:16} log K "
                  f"{best['logK']:.2f}   {why(best)}")
            if len(tied) > 1:
                lk = [r["logK"] for r in tied]
                tfu = [fu(r["logK"], 40.0, basis_of(r["basis_text"]))
                       for r in tied]
                print(f"    TIED    {len(tied)} candidates the rule "
                      f"cannot separate: log K "
                      f"{min(lk):.2f} to {max(lk):.2f}")
                print(f"            fu {min(tfu):.5f} to {max(tfu):.5f}"
                      f" ({max(tfu)/min(tfu):.1f}-fold) on the tie "
                      f"alone -- the choice between them is arbitrary")
            lo = min(fu(r["logK"], 40.0, basis_of(r["basis_text"]))
                     for r in cand)
            hi = max(fu(r["logK"], 40.0, basis_of(r["basis_text"]))
                     for r in cand)
            line = "    fu     "
            for lab, g in ALBUMIN_G_L:
                line += f" {lab} {fu(best['logK'], g, b):.5f} "
            print(line)
            print(f"    the whole admissible set at 40 g/L spans fu "
                  f"{lo:.5f} to {hi:.5f}  ({hi/lo:.0f}-fold)")
            for r in cand[1:]:
                print(f"      also  {r['source']:16} log K "
                      f"{r['logK']:6.2f}  fu "
                      f"{fu(r['logK'], 40.0, basis_of(r['basis_text'])):.5f}"
                      f"   {why(r)}")
            tfu = [fu(r["logK"], 40.0, basis_of(r["basis_text"]))
                   for r in tied]
            rec = dict(species=species, pfas=c, n_admissible=len(cand),
                       n_tied_at_top=len(tied),
                       tie_fu_lo_40=round(min(tfu), 6),
                       tie_fu_hi_40=round(max(tfu), 6),
                       chosen_source=best["source"],
                       chosen_logK=round(best["logK"], 3),
                       chosen_basis=b, chosen_because=why(best),
                       fu_lo_40=round(lo, 6), fu_hi_40=round(hi, 6),
                       spread_fold=round(hi / lo, 1))
            for lab, g in ALBUMIN_G_L:
                rec[f"fu_{lab.split()[0]}"] = round(
                    fu(best["logK"], g, b), 6)
            out.append(rec)
        print()

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")


if __name__ == "__main__":
    main()
