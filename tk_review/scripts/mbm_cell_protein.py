#!/usr/bin/env python3
"""Candidate 1: the sorptive capacity assigned to CELL PROTEIN.

Section 6.1d left two candidates for the models' over-prediction. This is the
first: Qin's measured non-specific protein affinity against Armitage's
convention that 3.5% of protein behaves as octanol. They are guesses at the
same quantity and they are the substance of each model's cell term -- Qin's
cell is 92-95% protein, and Armitage's protein term is the one part of its
cell that is NOT lipid.

What makes the test possible now is Qin et al. 2024 (ES&T 58:9954), the same
group, which measured specific and nonspecific protein-water distribution
ratios for 16 PFAS by SPME over four orders of magnitude. The 2023
mass-balance model had measured D_nonspecific for only PFBA and PFOA and
interpolated the rest; 2024 measures all five of ours.

Run:  python3 tk_review/scripts/mbm_cell_protein.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_published_equations as P
import mbm_qin_vs_armitage as Q
import mbm_fischer_albumin as F
import mbm_cell_term as C

DB = os.path.join(os.path.dirname(HERE), "db")
LINE = "=" * 78
ORDER = ["PFBA", "PFHxA", "PFHpA", "PFOA", "PFNA"]

# Qin 2024 Table 2, BSA column, log D_protein/water at pH 7.4 (volume basis).
QIN24 = {           # chem: (specific, nonspecific)
    "PFBA":  (2.44, 1.97),
    "PFHxA": (3.26, 2.69),
    "PFHpA": (4.12, 3.47),
    "PFOA":  (4.58, 3.88),
    "PFNA":  (None, 4.61),      # specific not resolvable
}


def main():
    print(__doc__.split("Run:")[0].rstrip())

    # ------------------------------------------------------------------ 1
    print(LINE)
    print("1.  THREE ESTIMATES OF THE SAME QUANTITY")
    print(LINE)
    print("""
  log D_protein/water on the volume basis both models use.  "Qin 2023 model"
  is what the published mass balance actually runs on; "Qin 2024 measured" is
  the same laboratory's later direct measurement; "Armitage" is 0.035 x D_ow,
  the octanol convention, expressed as an effective protein partition
  coefficient so the three are comparable.
""")
    print(f"  {'cmpd':8s}{'nCF':>4s}{'Qin23 model':>13s}{'Qin24 meas':>12s}"
          f"{'Armitage':>10s}{'  23 vs 24':>11s}{'Arm vs 24':>11s}  Qin23 provenance")
    rows = []
    for nm in ORDER:
        nCF, lk = C._ncf_lk(nm)
        p = Q.qin_params(nm, nCF)
        q23 = p["logDnonspec"]
        q24 = QIN24[nm][1]
        arm = math.log10(P.P_NLOM * P.armitage_D(lk, P.CELL_COMP["AREc2"])["Dow"])
        print(f"  {nm:8s}{nCF:>4d}{q23:>13.2f}{q24:>12.2f}{arm:>10.2f}"
              f"{10**(q23-q24):>11.2f}{10**(arm-q24):>11.1e}  {p['source']}")
        rows.append((nm, nCF, q23, q24, arm))
    s23 = C._slope([(r[1], r[2]) for r in rows])
    s24 = C._slope([(r[1], r[3]) for r in rows])
    sar = C._slope([(r[1], r[4]) for r in rows])
    print(f"""
  decades per CF2:  Qin 2023 model {s23:+.3f}   Qin 2024 measured {s24:+.3f}   Armitage {sar:+.3f}

  Two things stand out.  Qin's 2023 interpolation and his 2024 measurement
  agree closely in slope and to {min(10**abs(r[2]-r[3]) for r in rows):.1f}-{max(10**abs(r[2]-r[3]) for r in rows):.1f}x in level, so the model was not
  mis-parameterised by the interpolation -- the three compounds it had to
  interpolate came out near enough.

  Armitage's octanol surrogate is a different story.  It runs
  {min(10**(r[4]-r[3]) for r in rows):.1e} to {max(10**(r[4]-r[3]) for r in rows):.1e} times the measured protein affinity -- that is
  {1/max(10**(r[4]-r[3]) for r in rows):.0f}x to {1/min(10**(r[4]-r[3]) for r in rows):.0f}x TOO LOW, three to four orders of magnitude -- and its
  slope is {sar/s24:.1f}x too steep because D_ow inherits the full octanol-water slope
  while real protein binding is shallower.  It never catches up: the gap is
  worst at PFBA and still a factor of {1/max(10**(r[4]-r[3]) for r in rows):.0f} at PFNA.

  That explains the composition finding of section 6.1d.  Armitage's cell is
  ~100% membrane lipid not because protein does not bind PFAS, but because
  the octanol convention has defined its protein term out of existence.  A
  partition coefficient calibrated on neutral organics does not describe an
  anionic surfactant binding to protein.
""")

    # ------------------------------------------------------------------ 2
    print(LINE)
    print("2.  DOES THE MEASURED VALUE BEAT THE CONVENTION IN THE BALANCE?")
    print(LINE)
    variants = [
        ("Qin, published", "qin", None),
        ("Qin + Qin 2024 cell protein", "qin", "q24"),
        ("Armitage, published", "arm", None),
        ("Armitage + Qin 2024 cell protein", "arm", "q24"),
    ]
    series = {v[0]: [] for v in variants}
    bycmp = {v[0]: {} for v in variants}
    bydose = {v[0]: {"low": [], "high": []} for v in variants}
    out = []
    for study, a in P.ASSAYS.items():
        V = P.volumes(a)
        for nm in ORDER:
            nCF, lk = C._ncf_lk(nm)
            for d in a["doses"]:
                m = P.measured_x_cell(study, nm, d, V, a)
                if m is None:
                    continue
                rec = dict(study=study, compound=nm, dose_uM=d,
                           qin2024_logD_nonspecific=QIN24[nm][1],
                           meas_x_cell_pct=round(m * 100, 4))
                for label, kind, mode in variants:
                    x = _run(kind, nm, lk, d, V,
                             10 ** QIN24[nm][1] if mode else None)
                    rec[label.lower().replace(" ", "_").replace(",", "")
                        .replace("+", "plus") + "_x_cell_pct"] = round(x * 100, 4)
                    series[label].append((x, m))
                    bycmp[label].setdefault(nm, []).append((x, m))
                    bydose[label]["low" if d <= 25 else "high"].append((x, m))
                out.append(rec)
    print(f"\n  {'variant':36s}{'bias':>8s}{'|err|':>8s}{'<=25 uM':>10s}{'100 uM':>9s}")
    for label, _k, _m in variants:
        b, e, n = F.stats(series[label])
        bl, _, _ = F.stats(bydose[label]["low"])
        bh, _, _ = F.stats(bydose[label]["high"])
        print(f"  {label:36s}{b:>8.2f}{e:>8.2f}{bl:>10.2f}{bh:>9.2f}")
    print(f"\n  bias by compound:\n")
    print(f"  {'variant':36s}" + "".join(f"{c:>9s}" for c in ORDER))
    for label, _k, _m in variants:
        line = f"  {label:36s}"
        for c in ORDER:
            b, _, _ = F.stats(bycmp[label].get(c, []))
            line += f"{b:>9.2f}" if b else f"{'-':>9s}"
        print(line)

    # ------------------------------------------------------------------ 3
    print("\n" + LINE)
    print("3.  WHAT QIN 2024 SETTLES THAT THE MODELS CANNOT")
    print(LINE)
    print(f"\n  {'cmpd':8s}{'specific':>10s}{'nonspecific':>13s}{'ratio':>8s}"
          f"   which one applies where")
    for nm in ORDER:
        sp, ns = QIN24[nm]
        if sp is None:
            print(f"  {nm:8s}{'n/a':>10s}{ns:>13.2f}{'-':>8s}"
                  f"   specific term not resolvable")
            continue
        print(f"  {nm:8s}{sp:>10.2f}{ns:>13.2f}{10**(sp-ns):>8.1f}"
              f"   physiological / in vitro")
    print("""
  This is the measurement sections 6.1a and 6.1b could only infer.  The SAME
  protein and the SAME compound give two constants differing 3-7x, and which
  one you get depends entirely on the PFAS:protein molar ratio: specific
  binding dominates below nu ~ 0.01, which is the physiological regime, and
  nonspecific binding takes over above nu ~ 1, which is where in vitro assays
  are dosed.  Qin's own mass balance correctly uses the nonspecific value,
  because the assays it models are dosed in excess.  A PBTK model of a human
  must not.
""")

    # ------------------------------------------------------------------ 4
    print(LINE)
    print("4.  THE PATTERN ACROSS ALL THREE SUBSTITUTIONS, AND WHAT IS LEFT")
    print(LINE)
    print("""
  Three parameters have now been replaced by the best available measurement:

    albumin      Fischer 2024      measured BELOW both models  -> smaller sink
    cell D_cell/w Qin 2023 Table 2  mixed                      -> centres, disperses
    cell protein Qin 2024          measured ABOVE both models  -> bigger sink

  The first two push dose INTO the cell by shrinking the competing sink; the
  third pushes it in by enlarging the cell's own.  All three therefore move
  the prediction the same way, and all three make models that already
  over-predict worse.  No partition coefficient in this system is going to
  fix an over-prediction, because the measured values do not point that way.

  That leaves the one parameter that scales x_cell directly, with no
  partition coefficient involved: the cell volume.  Section 6.1d put it at
  3.1x uncertain.  Armitage's AREc2 is 10 ng/cell, i.e. ~10 pL at 1 mg/mL;
  the project spreadsheet has 3.04 pL/cell for HepG2 with the same protein
  per cell.  Rerunning the published models on the smaller cell:
""")
    print(f"  {'variant':46s}{'bias':>8s}{'|err|':>8s}{'<=25 uM':>10s}{'100 uM':>9s}")
    for label, kind in (("Qin, published", "qin"), ("Armitage, published", "arm")):
        for ngc, tag in ((10.0, "10 ng/cell (AREc2, as published)"),
                         (3.04, "3.04 ng/cell (project sheet, HepG2)")):
            pairs, lo, hi = [], [], []
            for study, a in P.ASSAYS.items():
                a2 = dict(a)
                comp = dict(P.CELL_COMP[a["line"]], ng_cell=ngc)
                V = P.volumes(a2)
                V = dict(V, Vcells=a["cells"] * ngc * 1e-12,
                         Vprot_cell=a["cells"] * ngc * 1e-12 * comp["prot"],
                         Vlip_cell=a["cells"] * ngc * 1e-12 * (comp["sl"] + comp["ml"]))
                for nm in ORDER:
                    nCF, lk = C._ncf_lk(nm)
                    for d in a["doses"]:
                        m = P.measured_x_cell(study, nm, d, V, a2)
                        if m is None:
                            continue
                        x = _run(kind, nm, lk, d, V)
                        pairs.append((x, m))
                        (lo if d <= 25 else hi).append((x, m))
            b, e, n = F.stats(pairs)
            bl, _, _ = F.stats(lo)
            bh, _, _ = F.stats(hi)
            print(f"  {label.split(',')[0] + ' @ ' + tag:46s}{b:>8.2f}{e:>8.2f}"
                  f"{bl:>10.2f}{bh:>9.2f}")
    print("""
  The cell volume moves the bias almost exactly as its ratio predicts,
  because x_cell is very nearly proportional to V_cells while the cell holds
  a small share of the dose.  That is the point: it is the one knob that
  rescales every compound by the SAME factor, which is what an unexplained
  level offset that does not vary with chain length actually needs.

  It is also not a free parameter.  Both cell volumes come from real
  measurements of different cell lines, and the assays being modelled are
  HepG2 and HepaRG, not AREc2.  Resolving it needs the cell volume of the
  lines actually used, which none of the three uptake papers reports.  That
  is where this line of work stops until such a measurement exists.
""")

    dst = os.path.join(DB, "mbm_cell_protein.csv")
    keys = []
    for r in out:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    print(f"  wrote {dst} ({len(out)} rows)")


def _run(kind, nm, lk, dose, V, Dprot=None):
    """x_cell with the CELL PROTEIN partition optionally overridden."""
    ntot = dose * V["Vbm"]
    comp = V["comp"]
    if kind == "arm":
        D = P.armitage_D(lk, comp)
        if Dprot is None:
            Dcw = D["Dcw"]
        else:
            Dcw = (comp["sl"] * D["Dow"] + comp["ml"] * D["Dmw"]
                   + comp["prot"] * Dprot + comp["water"])
        den = (V["Vw"] + D["Dsa"] * V["Valb"] + D["Dow"] * V["Vslip"]
               + Dcw * V["Vcells"] + D["Dpl"] * V["sarea"] * 1000.0)
        Cw = ntot / den
        return Dcw * Cw * V["Vcells"] / ntot
    Cw = dose * 1e-3
    for _ in range(500):
        Dbsa, Dsp, Dlip, _p = P.qin_D(nm, Cw)
        sp = Dsp if Dprot is None else Dprot
        den = (V["Vw"] + Dbsa * V["Valb"] + Dlip * V["Vslip"]
               + sp * V["Vprot_cell"] + Dlip * V["Vlip_cell"])
        new = ntot / den
        if abs(math.log10(max(new, 1e-30)) - math.log10(max(Cw, 1e-30))) < 1e-12:
            Cw = new
            break
        Cw = math.sqrt(Cw * new)
    Dbsa, Dsp, Dlip, _p = P.qin_D(nm, Cw)
    sp = Dsp if Dprot is None else Dprot
    n_cell = sp * V["Vprot_cell"] * Cw + Dlip * V["Vlip_cell"] * Cw
    return n_cell / ntot


if __name__ == "__main__":
    main()
