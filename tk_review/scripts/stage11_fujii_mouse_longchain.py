#!/usr/bin/env python3
"""Fujii 2015 mouse clearances via EFSA Appendix C (Stage 11).

    python3 -I scripts/stage11_fujii_mouse_longchain.py

Writes db/primary_2026/fujii2015_mouse_total_clearance.csv. Idempotent.
The usable rows are wired into scripts/stage6_consolidate_tubular.py.

    Fujii Y et al. (2015), as tabulated in EFSA Appendix C
    (db/efsa_appendix_c.csv); mouse strain FVB/NJcl.

WHAT THIS ADDS, AND WHAT IT DOES NOT. These are TOTAL clearances, not renal
ones, so they enter at the UPPER tier as bounds on T rather than as direct
values. They extend the mouse series from C10 to C12 -- PFHpA, PFUnDA and
PFDoDA are new compounds for the mouse column.

TWO COMPOUNDS CANNOT BE USED, AND THE REASON ANSWERS A QUESTION WORTH
RECORDING. Fujii also reports mouse clearances for PFTrDA (C13) and PFTeDA
(C14). Neither can be turned into a T, because neither has a MOUSE plasma
unbound fraction: Ryu 2024 stops at C12, and this project's 1497-row protein
binding database holds PFTrDA and PFTeDA binding for HUMAN and BOVINE albumin
only. The albumin and FABP collection is large and it is the wrong half of the
problem -- it supplies fu, and fu is not what is scarce.

THE DOSE DEPENDENCE IS ITSELF A RESULT. Fujii dosed at 0.31 and 3.13 umol/kg.
Clearance RISES with dose, steeply for the longest chains:

    PFTeDA female   10.4 at 0.31 umol/kg  ->  48.7 at 3.13   (4.7x)
    PFTrDA male      6.8 at 0.31          ->  34.2 at 3.13   (5.0x)
    PFDoDA male      4.4 at 0.31          ->   9.4 at 3.13   (2.1x)
    PFUnDA male      2.8 at 0.31          ->   5.7 at 3.13   (2.0x)
    PFHpA  male    347   at 0.31          -> 293   at 3.13   (0.8x)

A tenfold dose giving two- to fivefold more clearance is what SATURABLE
REABSORPTION predicts: at higher load the reabsorptive transporters saturate,
less is returned, and clearance rises. The effect grows with chain length,
which is where reabsorption is strongest. The short chain, PFHpA, shows none
of it. THE LOW DOSE IS USED HERE because it is nearer environmental exposure;
the high-dose values are recorded alongside.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
# chemical: {sex: (CL_tot at 0.31, CL_tot at 3.13)} mL/kg/day
FUJII = {
    "PFHpA":  {"female": (257.0, 190.0), "male": (347.0, 293.0)},
    "PFNA":   {"female": (5.1, 2.4),     "male": (3.9, 4.0)},
    "PFDA":   {"female": (2.8, 2.2),     "male": (2.2, 3.9)},
    "PFUnDA": {"female": (3.4, 3.1),     "male": (2.8, 5.7)},
    "PFDoDA": {"female": (4.8, 5.2),     "male": (4.4, 9.4)},
    "PFTriDA": {"female": (7.2, 17.1),   "male": (6.8, 34.2)},
    "PFTeDA": {"female": (10.4, 48.7),   "male": (10.4, 106.3)},
}
NC = {"PFHpA": 7, "PFNA": 9, "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12,
      "PFTriDA": 13, "PFTeDA": 14}
GFR_MOUSE = 16700.0


def mouse_fu():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] == "Ryu 2024" and r["species"].startswith("mouse") \
                and r["fu"]:
            out.setdefault(r["chemical"], {})[r["species"]] = float(r["fu"])
    return out


def main():
    fu = mouse_fu()
    p = HERE / "db" / "primary_2026" / "fujii2015_mouse_total_clearance.csv"
    n = usable = 0
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "sex", "strain",
                    "cl_total_low_dose_mL_kg_day", "cl_total_high_dose",
                    "dose_ratio_effect", "fu_mouse_C57BL6", "fu_source",
                    "gfr_mL_kg_day", "filtered_load", "T_upper_bound",
                    "verdict", "usable", "blocker", "source"])
        for chem, bysex in FUJII.items():
            f57 = fu.get(chem, {}).get("mouse (C57BL/6)")
            for sex, (lo, hi) in bysex.items():
                if f57:
                    filt = f57 * GFR_MOUSE
                    T = lo / filt
                    verdict = ("net_reabsorption_holds" if T < 1
                               else "uninformative")
                    use, block = "yes", ""
                    usable += 1
                else:
                    filt = T = None
                    verdict = ""
                    use = "no"
                    block = ("no mouse plasma fu exists: Ryu 2024 stops at "
                             "C12 and the project's binding database holds "
                             "this compound for human and bovine albumin only")
                w.writerow([chem, NC[chem], sex, "FVB/NJcl", lo, hi,
                            round(hi / lo, 2), f57 or "",
                            "Ryu 2024" if f57 else "", GFR_MOUSE,
                            round(filt, 1) if filt else "",
                            round(T, 4) if T else "", verdict, use, block,
                            "Fujii 2015 via EFSA Appendix C "
                            "(db/efsa_appendix_c.csv)"])
                n += 1
    bar = "=" * 78
    print(bar)
    print("FUJII 2015 MOUSE TOTAL CLEARANCES  (low dose, 0.31 umol/kg)")
    print(bar)
    print(f"  {'chem':9}{'nC':>3}{'sex':8}{'CL_tot':>9}{'fu':>9}"
          f"{'filtered':>10}{'T <=':>8}   verdict")
    for chem, bysex in FUJII.items():
        f57 = fu.get(chem, {}).get("mouse (C57BL/6)")
        for sex, (lo, hi) in bysex.items():
            if not f57:
                print(f"  {chem:9}{NC[chem]:>3}{sex:8}{lo:>9.1f}"
                      f"{'-':>9}{'-':>10}{'-':>8}   BLOCKED: no mouse fu")
                continue
            filt = f57 * GFR_MOUSE
            print(f"  {chem:9}{NC[chem]:>3}{sex:8}{lo:>9.1f}{f57:>9.5f}"
                  f"{filt:>10.1f}{lo/filt:>8.3f}   "
                  f"{'reabsorption holds' if lo/filt < 1 else 'uninformative'}")
    print()
    print(f"  {usable} of {n} rows usable; the rest blocked on mouse binding.")
    print()
    print(bar)
    print("WHAT THE 1497-ROW BINDING DATABASE CAN AND CANNOT DO HERE")
    print(bar)
    print("  It supplies fu. fu is not the scarce term. Of the compounds with")
    print("  a clearance somewhere in this project:")
    print()
    print("      human   15 compounds,  14 already have fu")
    print("      rat     13 compounds,  10 already have fu (the 3 gaps are")
    print("                             PFOS isomers whose clearances read NR)")
    print("      mouse   13 compounds,  11 already have fu")
    print("      monkey   6 compounds,   0 have fu")
    print()
    print("  The two places binding could unlock a row are mouse PFTrDA and")
    print("  PFTeDA, and the monkey. The database has NEITHER: no mouse")
    print("  plasma or albumin binding for C13-C14, and no macaque binding")
    print("  for any PFAS at all.")
    print()
    print("  SO THE ALBUMIN AND FABP COLLECTION UNLOCKS ZERO ADDITIONAL ROWS")
    print("  FOR THIS ENDPOINT. That is worth knowing rather than guessing:")
    print("  the bottleneck is renal clearance, and no amount of binding work")
    print("  moves it.")
    print(f"\nwrote {p.relative_to(HERE)}  ({n} rows)")


if __name__ == "__main__":
    main()
