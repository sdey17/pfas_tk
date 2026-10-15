"""
Diff the project spreadsheet's mass-balance predictions against this repo's
implementation, in the spreadsheet's own units (nmol/mg protein).

Four tests the spreadsheet's numbers answer on their own, before any model:
  1  are the sheet's predictions linear in dose?  (Qin's should not be)
  2  Rosenmai vs Amstutz at 100 uM -- same cell line, dose and stated unit
  3  three-way agreement once Rosenmai is divided by MW
  4  the measured dose exponent b in C_cell ~ C_nom^b, against what each
     model requires

Run:  python3 tk_review/scripts/mbm_vs_project_sheet.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mbm_qin_vs_armitage as M

DB = os.path.join(os.path.dirname(HERE), "db")
SHEET = os.path.join(DB, "mbm_project_sheet.csv")
MW = {c[0]: c[3] for c in M.COMPOUNDS}
NCF = {c[0]: c[2] for c in M.COMPOUNDS}
LKOW = {c[0]: (c[4] if c[4] is not None else None) for c in M.COMPOUNDS}


def logkow(name):
    v = LKOW[name]
    return v if v is not None else M.KOW_OPERA(NCF[name])


def to_nmol_per_mg(C_cell_uM):
    """whole-cell uM -> nmol per mg protein (uM = nmol/mL; / mg protein per mL)"""
    return C_cell_uM / M.MG_PROTEIN_PER_ML_CELL


def loglog_slope(pairs):
    n = len(pairs)
    lx = [math.log10(x) for x, _ in pairs]
    ly = [math.log10(y) for _, y in pairs]
    sx, sy = sum(lx), sum(ly)
    sxx = sum(v * v for v in lx)
    sxy = sum(a * b for a, b in zip(lx, ly))
    return (n * sxy - sx * sy) / (n * sxx - sx ** 2)


def load():
    rows = []
    for r in csv.DictReader(open(SHEET)):
        r["Ctest_uM"] = float(r["Ctest_uM"])
        for k in ("Ccell_meas_nmolmg", "Ccell_qin_nmolmg_sheet",
                  "Ccell_arm_nmolmg_sheet"):
            r[k] = float(r[k])
        rows.append(r)
    return rows


# Fragki's supplementary Table 2 gives only 6.25 / 25 / 100 uM. The sheet's
# 25 and 100 uM values match it to 1-3%; its 6.25 uM values do not, and
# 12.5 / 50 uM are not in the supplementary at all (they must come from the
# main paper, which is not on hand).
FRAGKI_SUPPL_T2 = {"PFOA": {6.25: 0.09226, 25.0: 0.52306, 100.0: 2.79309},
                   "PFNA": {6.25: 0.22203, 100.0: 15.69467}}


def main():
    rows = load()

    print("=" * 78)
    print("TEST 1  Are the sheet's predictions linear in dose?")
    print("        C_cell / C_test is constant iff the model is linear.")
    print("=" * 78)
    series = [("Fragki", "PFOA"), ("Fragki", "PFNA"), ("Rosenmai", "PFOA"),
              ("Rosenmai", "PFNA"), ("Rosenmai", "PFBA")]
    hdr = "%-26s" % "" + "".join("%12s" % (s[:4] + "/" + c) for s, c in series)
    print(hdr)
    for lab, key in (("sheet Qin", "Ccell_qin_nmolmg_sheet"),
                     ("sheet Armitage", "Ccell_arm_nmolmg_sheet"),
                     ("MEASURED", "Ccell_meas_nmolmg")):
        line = "%-26s" % lab
        for st, cp in series:
            v = sorted((r["Ctest_uM"], r[key] / r["Ctest_uM"])
                       for r in rows if r["study"] == st and r["compound"] == cp)
            line += "%11.2fx" % (v[-1][1] / v[0][1])
        print(line)
    print()
    print("  Both sheet columns are exactly linear. Qin's model is not: the")
    print("  saturable albumin term (Eqs. 12-16) is the point of that paper.")
    print("  The measured data are not linear either.")

    print()
    print("=" * 78)
    print("TEST 2  Rosenmai vs Amstutz at 100 uM: same cell line (HepG2),")
    print("        same dose, same stated unit (nmol/mg protein).")
    print("=" * 78)
    g = lambda s, c, d: [r["Ccell_meas_nmolmg"] for r in rows
                         if r["study"] == s and r["compound"] == c
                         and r["Ctest_uM"] == d][0]
    print("%-7s %10s %10s %9s | %12s %8s" %
          ("cmpd", "Rosenmai", "Amstutz", "ratio", "Rosenmai/MW", "ratio"))
    for c in ("PFHxA", "PFHpA", "PFOA", "PFNA"):
        R, A = g("Rosenmai", c, 100), g("Amstutz", c, 100)
        print("%-7s %10.1f %10.3f %8.0fx | %12.3f %7.2fx"
              % (c, R, A, R / A, R / MW[c], R / MW[c] / A))
    print()
    print("  As plotted, two HepG2 studies disagree by 205-911x. Dividing")
    print("  Rosenmai by molecular weight brings them to within 0.44-2.5x.")

    print()
    print("=" * 78)
    print("TEST 3  Three-way agreement at 100 uM once Rosenmai is corrected")
    print("=" * 78)
    for c in ("PFOA", "PFNA"):
        vals = [g("Fragki", c, 100), g("Amstutz", c, 100), g("Rosenmai", c, 100) / MW[c]]
        print("  %-5s Fragki(HepaRG) %7.2f  Amstutz %7.3f  Rosenmai/MW %7.3f"
              "   spread %.1fx" % (c, vals[0], vals[1], vals[2],
                                   max(vals) / min(vals)))
    print()
    print("  PFOA closes to 2.2x. PFNA does not: Fragki's HepaRG value is 7-16x")
    print("  the two HepG2 studies, which is a cell-line question, not a unit one.")

    print()
    print("=" * 78)
    print("TEST 4  Dose exponent b in C_cell ~ C_nom^b.")
    print("        Armitage requires b = 1.00 exactly. Qin's measured")
    print("        Freundlich slope n_Fr-1 implies b = 1/n_Fr.")
    print("=" * 78)
    print("%-20s %6s %6s %8s %10s" % ("series", "n", "b", "Qin b", "Armitage b"))
    for st, cp in (("Fragki", "PFOA"), ("Fragki", "PFNA"),
                   ("Rosenmai", "PFNA"), ("Rosenmai", "PFOA"),
                   ("Rosenmai", "PFHpA"), ("Rosenmai", "PFBA"),
                   ("Rosenmai", "PFHxA")):
        v = sorted((r["Ctest_uM"], r["Ccell_meas_nmolmg"])
                   for r in rows if r["study"] == st and r["compound"] == cp)
        b = loglog_slope(v)
        nfr = M.qin_params(cp, NCF[cp])["nFr_minus1"]
        print("%-20s %6d %6.2f %8.2f %10.2f"
              % (st + " " + cp, len(v), b, 1 / (1 + nfr), 1.00))
    # robustness: Fragki on the three supplementary-table values only
    v = sorted(FRAGKI_SUPPL_T2["PFOA"].items())
    print("%-20s %6d %6.2f %8s %10s"
          % ("Fragki PFOA (supplT2)", len(v), loglog_slope(v), "", ""))
    print()
    print("  b rises monotonically with chain length, 0.80 (PFHxA) to 1.42")
    print("  (PFNA). Armitage cannot produce b != 1 at all. Qin can, but its")
    print("  Freundlich slopes are nearly chain-length independent (-0.38 to")
    print("  -0.39), so it predicts b ~ 1.62 for every compound and misses the")
    print("  trend in the other direction.")

    print()
    print("=" * 78)
    print("This repo's predictions, in the sheet's units (nmol/mg protein)")
    print("=" * 78)
    print("%-9s %-6s %6s %9s %9s | %9s %9s %9s | %7s"
          % ("study", "cmpd", "dose", "meas", "meas/MW", "qin_nl", "qin_lin",
             "armitage", "src"))
    out = []
    for r in rows:
        cp, st, d = r["compound"], r["study"], r["Ctest_uM"]
        key = {"Rosenmai": "Rosenmai2017", "Amstutz": "Amstutz2022",
               "Fragki": "Fragki2023"}[st]
        V = M.volumes(M.ASSAYS[key])
        qn = to_nmol_per_mg(M.qin_predict(cp, NCF[cp], d, V, "mechanistic")["C_cell"])
        ql = to_nmol_per_mg(M.qin_predict(cp, NCF[cp], d, V, "linear")["C_cell"])
        am = to_nmol_per_mg(M.armitage_predict(logkow(cp), d, V)["C_cell"])
        corr = r["Ccell_meas_nmolmg"] / (MW[cp] if st == "Rosenmai" else 1.0)
        src = M.qin_params(cp, NCF[cp])["source"]
        print("%-9s %-6s %6g %9.3f %9.3f | %9.3f %9.3f %9.3f | %7s"
              % (st, cp, d, r["Ccell_meas_nmolmg"], corr, qn, ql, am, src[:7]))
        out.append(dict(study=st, compound=cp, dose_uM=d,
                        meas_sheet_nmolmg=r["Ccell_meas_nmolmg"],
                        meas_MW_corrected_nmolmg=round(corr, 4),
                        qin_nonlinear_nmolmg=round(qn, 4),
                        qin_linear_nmolmg=round(ql, 4),
                        armitage_nmolmg=round(am, 4),
                        qin_sheet_nmolmg=r["Ccell_qin_nmolmg_sheet"],
                        armitage_sheet_nmolmg=r["Ccell_arm_nmolmg_sheet"],
                        qin_param_source=src,
                        qin_over_meas=round(qn / corr, 2),
                        armitage_over_meas=round(am / corr, 2)))
    dst = os.path.join(DB, "mbm_sheet_reconciliation.csv")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    print()
    print("wrote %s  (%d rows)" % (dst, len(out)))


if __name__ == "__main__":
    main()
