#!/usr/bin/env python3
"""Fraction unbound from a binding constant: the derivation, and a test.

THE DERIVATION, in four steps.

1.  One ligand, one class of site.

        L + P  <->  PL          Kd = [L][P] / [PL]

    The fraction unbound is free ligand over total ligand:

        fu = [L] / ([L] + [PL])

    Substitute [PL] = [L][P]/Kd and the free-ligand term cancels:

        fu = [L] / ([L] + [L][P]/Kd) = 1 / (1 + [P]/Kd)

        fu = Kd / (Kd + [P])                                   (1)

    [P] here is FREE protein. In the trace limit, [L] << [P], free
    protein is approximately total protein, and for n equivalent
    sites per molecule the site concentration is n[P]_total:

        fu = Kd / (Kd + n[P]_total)                            (2)

    The trace limit is not a detail: it is the same assumption that
    makes a partition coefficient a constant at all. Outside it, fu
    depends on how much ligand is present and is not a property of
    the pair.

2.  The same thing in partition-coefficient form. A protein-water
    partition coefficient is

        K = C_protein / C_water

    with C_protein per kg (or per L) of protein and C_water molar.
    Bound ligand per litre of solution is K * C_water * c_p, where
    c_p is the protein content of the solution in the SAME basis as
    K -- kg/L if K is per kg, L/L if K is per L. Then

        fu = C_water / (C_water + K C_water c_p)

        fu = 1 / (1 + K c_p)                                   (3)

3.  The two forms are the same statement. In the trace limit the
    moles of sites per kg of protein are 1000n/M, each occupied with
    probability [L]/Kd, so

        K = 1000 n / (M Kd[M])  =  1e9 n / (M Kd[uM])          (4)

    which is this project's bridge, log K = log10(1e9/M) + log n
    - log Kd(uM). Substituting (4) into (3) returns (2) exactly.

4.  WHICH BASIS, which is where this goes wrong in practice. A
    coefficient per KILOGRAM of protein needs c_p in kg/L; one per
    LITRE of protein needs the protein's volume fraction, which is
    its mass concentration times its partial specific volume
    (0.733 mL/g for serum albumin). Getting that wrong costs
    log10(1/0.733) = 0.135 log -- a 1.4-fold error in fu, in the
    direction of over-predicting binding. Bischel 2011 reports on
    the VOLUME basis and says so; this script follows it.

THE TEST.  Bischel et al. 2011 is the only study in this database
that reports, for the same compounds in the same table, BOTH a
measured partition coefficient AND a measured fraction bound at a
stated albumin concentration. So equation (3) can be checked against
measurement with no cross-study assumption at all: same paper, same
experiment, same 200 uM BSA.

Bischel SM, Macmanus-Spencer LA, Zhang C, Luthy RG.  Strong
associations of short-chain perfluoroalkyl acids with serum albumin
and investigation of binding mechanisms.  Environ Toxicol Chem
2011;30(11):2423-30.  doi 10.1002/etc.647, PMID 21842491.

Ryu et al. 2024 cannot serve as this test, which is worth saying
plainly because it looks like it should. Its 320 rows are the largest
block of measured unbound fractions here, but they are TISSUE
HOMOGENATE -- "total protein plus phospholipid, not an isolated
protein" -- at a 1:5 dilution, with no albumin concentration and no
Kd. Predicting them needs the protein and lipid composition of each
tissue, which is a separate model and not a check on this equation.

Run:  python3 tk_review/scripts/funbound_from_binding.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "db", "funbound_predicted_vs_measured.csv")

M_ALB = 66463.0        # g/mol, serum albumin
V_BAR = 0.733e-3       # L/g, partial specific volume of serum albumin
BSA_UM = 200.0         # Bischel's albumin concentration

PHYSIOLOGICAL = [("plasma albumin, 35 g/L", 35.0),
                 ("plasma albumin, 40 g/L", 40.0),
                 ("plasma albumin, 50 g/L", 50.0)]


def fu_from_logK(logK, g_per_L, basis="volume"):
    """Equation (3). basis says what K is per: 'mass' -> kg of
    protein, 'volume' -> litres of protein."""
    c_p = (g_per_L * V_BAR) if basis == "volume" else (g_per_L / 1000.0)
    return 1.0 / (1.0 + (10.0 ** logK) * c_p)


def fu_from_Kd(kd_uM, n, prot_uM):
    """Equation (2)."""
    return kd_uM / (kd_uM + n * prot_uM)


def bischel():
    rows = []
    for r in csv.DictReader(open(DB, encoding="utf-8")):
        if not (r["study"].startswith("Bischel") and r["year"] == "2011"):
            continue
        if not (r["unbound_fraction"].strip() and r["ka"].strip()):
            continue
        rows.append((r["chemical"], float(r["ka"]),
                     float(r["unbound_fraction"]), r["ka_units"]))
    return rows


def main():
    print("=" * 76)
    print("CHECK: equation (3) against Bischel 2011's own measurements")
    print("=" * 76)
    # 200e-6 mol/L * 66463 g/mol = 13.3 g/L. No further factor:
    # mol/L times g/mol is already g/L.
    g = BSA_UM * 1e-6 * M_ALB                   # g/L of BSA
    print(f"  {BSA_UM:.0f} uM BSA = {g:.2f} g/L"
          f" = {g*V_BAR*1000:.2f} mL albumin per L"
          f"   (vbar {V_BAR*1e3:.3f} mL/g)")
    print(f"  Bischel reports log K on the VOLUME basis, so c_p = "
          f"{g*V_BAR:.5f} L/L\n")
    print(f"  {'PFAS':9}{'log K':>7}{'fu measured':>13}"
          f"{'fu predicted':>14}{'pred/meas':>11}{'  on mass basis':>16}")
    out, ratios = [], []
    for chem, logK, fu_meas, units in bischel():
        fu_vol = fu_from_logK(logK, g, "volume")
        fu_mass = fu_from_logK(logK, g, "mass")
        ratios.append(fu_vol / fu_meas)
        out.append(dict(pfas=chem, log_K=logK, fu_measured=fu_meas,
                        fu_predicted=round(fu_vol, 5),
                        ratio=round(fu_vol / fu_meas, 3),
                        fu_predicted_mass_basis=round(fu_mass, 5)))
        print(f"  {chem:9}{logK:7.2f}{fu_meas:13.3f}{fu_vol:14.4f}"
              f"{fu_vol/fu_meas:11.2f}{fu_mass:16.4f}")
    ratios.sort()
    med = ratios[len(ratios) // 2]
    print(f"\n  median predicted/measured {med:.2f}, "
          f"range {min(ratios):.2f} to {max(ratios):.2f}")

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"  wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")

    print("\n" + "=" * 76)
    print("CONSISTENCY: equation (2) and equation (3) are one statement")
    print("=" * 76)
    print("""
  Step 3 of the derivation says the Kd form and the partition form
  must agree exactly, through K = 1e9 n /(M Kd). If they ever differ
  here, the bridge constant is wrong somewhere.""")
    print(f"\n  {'Kd (uM)':>9}{'n':>4}{'log K':>9}{'fu via (2)':>13}"
          f"{'fu via (3)':>13}{'difference':>13}")
    for kd, n in ((0.1, 1), (1.0, 1), (10.0, 1), (100.0, 1),
                  (10.0, 3), (10.0, 6)):
        logK = math.log10(1e9 * n / (M_ALB * kd))
        a = fu_from_Kd(kd, n, 600.0)            # 600 uM albumin
        gl = 600.0 * 1e-6 * M_ALB               # the same, as g/L
        b = fu_from_logK(logK, gl, "mass")
        print(f"  {kd:9.4g}{n:4d}{logK:9.3f}{a:13.6f}{b:13.6f}"
              f"{abs(a-b):13.2e}")

    print("\n" + "=" * 76)
    print("WHAT RYU 2024 WOULD NEED, and why it is a different job")
    print("=" * 76)
    print("""
  Ryu's 320 rows are measured unbound fractions, but in TISSUE
  HOMOGENATE diluted 1:5 -- total protein plus phospholipid, no
  albumin concentration, no Kd. Two things stand between them and
  equation (3).

  First, the value reported is the DILUTED one and has to be undone
  before it means anything about tissue:

      fu_tissue = (1/D) / ( (1/fu_diluted - 1) + 1/D )

  at dilution factor D = 5. Second, predicting it needs the protein
  AND lipid content of each tissue, and a coefficient for each, which
  is a composition model rather than a check on this equation.""")
    print(f"\n  the undilution alone, at D = 5:")
    print(f"  {'fu measured (diluted)':>24}{'fu in tissue':>16}")
    for fud in (0.3, 0.1, 0.01, 0.001, 1e-4):
        D = 5.0
        fut = (1 / D) / ((1 / fud - 1) + 1 / D)
        print(f"  {fud:>24.5g}{fut:16.5g}")

    print("\n" + "=" * 76)
    print("WHY fu IS A BRUTAL TEST OF A BINDING CONSTANT")
    print("=" * 76)
    print("""
  At physiological albumin almost every PFAS is >99 % bound, so fu is
  a small difference between two large numbers and an error in log K
  lands on it almost undiluted. One log of disagreement between
  laboratories is one log in fu.""")
    print(f"\n  {'log K':>7}" + "".join(f"{lab.split(',')[1]:>14}"
                                        for lab, _ in PHYSIOLOGICAL))
    for lk in (3.0, 3.5, 4.0, 4.5, 5.0, 5.5):
        line = f"  {lk:7.1f}"
        for _lab, gl in PHYSIOLOGICAL:
            line += f"{fu_from_logK(lk, gl):14.5f}"
        print(line)
    print("\n  The albumin literature for a single compound spans about")
    print("  1.2 log (median, same method, two laboratories), which is")
    print("  a 16-fold spread in fu at any albumin concentration.")


if __name__ == "__main__":
    main()
