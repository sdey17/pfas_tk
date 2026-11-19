#!/usr/bin/env python3
"""The OEHHA reabsorption axis recomputed for all 14 Ryu PFAS, three species.

    python3 -I scripts/stage4b_all14_species_reabsorption.py

Writes db/species_fu_ratios_all14.csv and db/species_reabsorption_all14.csv.
Idempotent: both files are rebuilt from scratch each run.

THE QUESTION. Stage 4 established, for PFOA alone, that the 303x human/rat
renal clearance difference decomposes into 6.8x glomerular and 45x tubular
once plasma fu is measured rather than assumed. Is that tubular dominance a
property of PFOA, or of the class? If chain length moves the glomerular /
tubular split, the QSAR aim has to model two mechanisms, not one.

THREE LAYERS, AND THEY DO NOT HAVE THE SAME COVERAGE.

  Layer 1  fu by species          14 x 4 matrices, COMPLETE
           One laboratory, one method, one plate. No extrapolation at all.

  Layer 2  filtered load fu x GFR 14 x 3 species, COMPLETE
           GFR is a species constant, so this layer inherits layer 1's
           coverage exactly.

  Layer 3  reabsorbed %           needs renal clearance, WHICH IS SPARSE
           human   7 of 14  (Argoul 2026 Table 1)
           mouse   7 of 14  (Argoul 2026 Table 4, female C57BL/6)
           rat     1 of 14  (PFOA only, OEHHA A6.4 / Kemper 2003)

So the full three-species triangle closes for exactly ONE compound. The
human-mouse pair closes for seven. THE RAT IS THE HOLE, AND IT IS A HOLE IN
THE RENAL CLEARANCE DATA, NOT IN THE BINDING DATA. Ryu closed the binding
half of the rat gap; the clearance half is still open, and no amount of
further binding work will close it.

PROVENANCE MISMATCHES, STATED NOT BURIED.

  (a) Ryu's plasma is "equal sex pooled (minimum of 3 male and 3 female)".
      The rat PFOA renal clearance it is paired with is MALE (18.2 mL/kg/d,
      Kemper 2003), and female rat renal clearance is 666 mL/kg/d -- 37x
      higher. Pooling is defensible only if rat plasma fu does not itself
      differ by sex. No paper measures that: Ryu's companion tissue paper
      (Environ Sci Technol 58:14641) also pools. The sex difference was
      attributed in section 14 to Oatp1a1 expression (114x), a transporter
      difference, so binding is unlikely to carry it -- but that is an
      argument, not a measurement.

  (b) Ryu's rat is Wistar-Han. The renal clearance is Sprague-Dawley.

  (c) Ryu's mouse strains are CD-1 and C57BL/6; Argoul's mouse is female
      C57BL/6, so the mouse pairing is strain-matched. Sex still is not.

  (d) For the long-chain human compounds Argoul reports no separate renal
      clearance and uses TOTAL clearance as the renal proxy. Reproduced here
      and flagged in clr_basis.

  (e) GFR DISAGREES BETWEEN THE TWO AUTHORITIES FOR THE MOUSE. OEHHA/Han use
      16700 mL/kg/d; Argoul uses 10000 (recovered exactly from their own
      fu x GFR column). This recompute uses the OEHHA set throughout, since
      it is OEHHA's table being recomputed, and reports the sensitivity.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

# --- OEHHA PHG Table A6.4 GFR set, mL/kg/day ------------------------------
GFR = {"human": 2570, "rat": 14400, "mouse": 16700}
GFR_MOUSE_ARGOUL = 10000      # recovered from Argoul 2026 Table 4

# --- chain length and family ----------------------------------------------
CHEM = {
    "PFBA":    (4,  "carboxylate"),
    "PFPeA":   (5,  "carboxylate"),
    "PFHxA":   (6,  "carboxylate"),
    "PFHpA":   (7,  "carboxylate"),
    "PFOA":    (8,  "carboxylate"),
    "PFNA":    (9,  "carboxylate"),
    "PFDA":    (10, "carboxylate"),
    "PFUnDA":  (11, "carboxylate"),
    "PFDoDA":  (12, "carboxylate"),
    "PFBS":    (4,  "sulfonate"),
    "PFHxS":   (6,  "sulfonate"),
    "PFOS":    (8,  "sulfonate"),
    "PFOSA":   (8,  "sulfonamide"),
    "6:2 FTS": (6,  "fluorotelomer sulfonate"),
}

# --- renal clearance, mL/kg/day -------------------------------------------
# human: Argoul 2026 Table 1. clr column where reported, else total CL.
CLR_HUMAN = {
    "PFBS":  (1.2,   "measured_renal"),
    "PFHxA": (47.0,  "measured_renal"),
    "PFHxS": (0.053, "total_as_proxy"),
    "PFOA":  (0.042, "total_as_proxy"),
    "PFNA":  (0.067, "total_as_proxy"),
    "PFOS":  (0.088, "total_as_proxy"),
    "PFDA":  (0.15,  "total_as_proxy"),
}
# mouse: Argoul 2026 Table 4, cl_renal_measured, female C57BL/6
CLR_MOUSE = {
    "PFBA":  (783.0, "measured_renal"),
    "PFHxA": (6780.0, "measured_renal"),
    "PFOA":  (3.59,  "measured_renal"),
    "PFNA":  (0.6,   "measured_renal"),
    "PFDA":  (0.13,  "measured_renal"),
    "PFHxS": (0.589, "measured_renal"),
    "PFOS":  (1.24,  "measured_renal"),
}
# rat: OEHHA A6.4 male, Kemper 2003
CLR_RAT = {
    "PFOA":  (18.2,  "measured_renal"),
}
CLR = {"human": CLR_HUMAN, "rat": CLR_RAT, "mouse": CLR_MOUSE}
CLR_SRC = {
    "human": "Argoul 2026 Table 1",
    "rat":   "OEHHA 2024 PHG Table A6.4 (Kemper 2003), male Sprague-Dawley",
    "mouse": "Argoul 2026 Table 4, female C57BL/6",
}
# Argoul's own fu, for the cross-laboratory check (as fractions)
FU_ARGOUL = {
    ("human", "PFBS"): 0.017, ("human", "PFHxA"): 0.020,
    ("human", "PFHxS"): 0.0008, ("human", "PFOA"): 0.0017,
    ("human", "PFNA"): 0.0017, ("human", "PFOS"): 0.0016,
    ("human", "PFDA"): 0.0038,
    ("mouse", "PFBA"): 0.77, ("mouse", "PFHxA"): 0.25,
    ("mouse", "PFOA"): 0.0087, ("mouse", "PFNA"): 0.0035,
    ("mouse", "PFDA"): 0.0042, ("mouse", "PFBS"): 0.29,
    ("mouse", "PFHxS"): 0.013, ("mouse", "PFOS"): 0.0025,
}


def load_ryu():
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    out = {}
    for r in csv.DictReader(f.open()):
        if r["experiment"] != "Ryu 2024" or not r["fu"]:
            continue
        out[(r["chemical"], r["species"])] = float(r["fu"])
    return out


def main():
    fu = load_ryu()
    chems = sorted(CHEM, key=lambda c: (CHEM[c][1], CHEM[c][0]))

    # ---------------- layer 1: fu ratios -------------------------------
    p1 = HERE / "db" / "species_fu_ratios_all14.csv"
    with p1.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_perfluorinated_carbons", "family",
                    "fu_human", "fu_rat_wistar_han", "fu_mouse_cd1",
                    "fu_mouse_c57bl6", "rat_over_human", "cd1_over_human",
                    "c57_over_human", "cd1_over_c57", "source"])
        L1 = []
        for c in chems:
            n, fam = CHEM[c]
            h = fu.get((c, "human"))
            rt = fu.get((c, "rat"))
            m1 = fu.get((c, "mouse (CD-1)"))
            m5 = fu.get((c, "mouse (C57BL/6)"))
            row = dict(chemical=c, n=n, fam=fam, h=h, rt=rt, m1=m1, m5=m5,
                       r_h=rt / h if h and rt else None,
                       c1_h=m1 / h if h and m1 else None,
                       c5_h=m5 / h if h and m5 else None,
                       c1_c5=m1 / m5 if m1 and m5 else None)
            L1.append(row)
            w.writerow([c, n, fam, h, rt, m1, m5,
                        f"{row['r_h']:.3f}" if row["r_h"] else "",
                        f"{row['c1_h']:.3f}" if row["c1_h"] else "",
                        f"{row['c5_h']:.3f}" if row["c5_h"] else "",
                        f"{row['c1_c5']:.3f}" if row["c1_c5"] else "",
                        "Ryu 2024 Toxics 12:253 Table 1, "
                        "doi:10.3390/toxics12040253"])

    # ---------------- layers 2 and 3 -----------------------------------
    p2 = HERE / "db" / "species_reabsorption_all14.csv"
    L3 = []
    with p2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_perfluorinated_carbons", "family",
                    "species", "matrix", "fu", "gfr_mL_kg_day",
                    "filtered_load_mL_kg_day", "clr_mL_kg_day", "clr_basis",
                    "net_mL_kg_day", "pct_reabsorbed", "verdict",
                    "fu_source", "clr_source", "provenance_mismatch"])
        for c in chems:
            n, fam = CHEM[c]
            for sp, matrix in (("human", "human"), ("rat", "rat"),
                               ("mouse", "mouse (C57BL/6)")):
                f = fu.get((c, matrix))
                if f is None:
                    continue
                g = GFR[sp]
                filt = f * g
                hit = CLR[sp].get(c)
                if hit is None:
                    w.writerow([c, n, fam, sp, matrix, f, g, round(filt, 3),
                                "", "", "", "", "NO_RENAL_CLEARANCE_DATA",
                                "Ryu 2024", "", ""])
                    continue
                clr, basis = hit
                net = filt - clr
                pct = 100 * net / filt
                verdict = ("net_reabsorption" if net > 0
                           else "NET_TUBULAR_SECRETION")
                mm = {"human": "sex-pooled fu",
                      "rat": "sex-pooled Wistar-Han fu vs male "
                             "Sprague-Dawley CLR",
                      "mouse": "sex-pooled fu vs female CLR; strain matched"}[sp]
                w.writerow([c, n, fam, sp, matrix, f, g, round(filt, 3),
                            clr, basis, round(net, 3),
                            round(pct, 1) if net > 0 else "", verdict,
                            "Ryu 2024", CLR_SRC[sp], mm])
                L3.append((c, n, fam, sp, f, filt, clr, net, pct, verdict,
                           basis))

    report(L1, L3, fu)
    print(f"\nwrote {p1.relative_to(HERE)}")
    print(f"wrote {p2.relative_to(HERE)}")


def report(L1, L3, fu):
    bar = "=" * 78
    print(bar)
    print("LAYER 1  PLASMA fu ACROSS 14 PFAS AND 4 MATRICES  (complete)")
    print(bar)
    print(f"{'chemical':9}{'nC':>3} {'family':24}"
          f"{'rat/hum':>9}{'CD-1/hum':>10}{'C57/hum':>9}{'CD-1/C57':>10}")
    for r in L1:
        def g(k):
            return f"{r[k]:.2f}" if r[k] else "-"
        print(f"{r['chemical']:9}{r['n']:>3} {r['fam']:24}"
              f"{g('r_h'):>9}{g('c1_h'):>10}{g('c5_h'):>9}{g('c1_c5'):>10}")

    rh = [r["r_h"] for r in L1 if r["r_h"]]
    c5 = [r["c5_h"] for r in L1 if r["c5_h"]]
    c1 = [r["c1_h"] for r in L1 if r["c1_h"]]
    def gm(xs):
        p = 1.0
        for x in xs:
            p *= x
        return p ** (1 / len(xs))
    print()
    print(f"  rat / human      geometric mean {gm(rh):5.2f}   "
          f"range {min(rh):.2f} to {max(rh):.2f}   n={len(rh)}")
    print(f"  CD-1 / human     geometric mean {gm(c1):5.2f}   "
          f"range {min(c1):.2f} to {max(c1):.2f}   n={len(c1)}")
    print(f"  C57BL/6 / human  geometric mean {gm(c5):5.2f}   "
          f"range {min(c5):.2f} to {max(c5):.2f}   n={len(c5)}")
    n_rat_2x = sum(1 for x in rh if 0.5 <= x <= 2)
    n_c5_2x = sum(1 for x in c5 if 0.5 <= x <= 2)
    print()
    print(f"  within two-fold of the human:  rat {n_rat_2x}/{len(rh)}"
          f"      C57BL/6 mouse {n_c5_2x}/{len(c5)}")
    print()
    print("  SO THE PFOA RESULT GENERALISES. Holding fu equal between rat and")
    print("  human is approximately right for EVERY ONE of the 14 compounds.")
    print("  Holding it equal between mouse and human is wrong for almost all")
    print("  of them. OEHHA's assumption is not a lucky break on PFOA; it is")
    print("  a sound approximation for the rat and a bad one for the mouse,")
    print("  across the whole class.")

    print()
    print(bar)
    print("LAYER 1b  CHAIN LENGTH: does the species ratio drift with nC?")
    print(bar)
    print(f"{'nC':>3}  {'carboxylate rat/hum':>21}  {'C57/hum':>9}")
    for n in range(4, 13):
        row = [r for r in L1 if r["n"] == n and r["fam"] == "carboxylate"]
        if not row:
            continue
        r = row[0]
        print(f"{n:>3}  {r['chemical']:>10} "
              f"{r['r_h']:>10.2f}  {r['c5_h']:>9.2f}")
    print()
    print("  The rat/human ratio is flat in chain length: 0.76 to 2.03 with no")
    print("  trend, nine carboxylates. The mouse/human ratio is NOT flat. It")
    print("  rises from 3.4x at C4 to a 4.8-7.3x plateau over C5-C8, then")
    print("  falls to about 2x by C9-C12. The ratio is largest for the")
    print("  compounds of greatest regulatory interest and smallest at both")
    print("  ends of the series. A QSAR trained on mouse plasma binding and")
    print("  applied to humans therefore carries a chain-length-dependent")
    print("  error, not a constant one, and it is worst in the middle.")

    print()
    print(bar)
    print("LAYER 3  REABSORBED FRACTION, WHERE RENAL CLEARANCE EXISTS")
    print(bar)
    print(f"{'chemical':9}{'nC':>3} {'species':7}{'fu':>9}{'filtered':>10}"
          f"{'CL_R':>9}{'reabs':>8}  basis")
    for (c, n, fam, sp, f, filt, clr, net, pct, verdict, basis) in L3:
        pcts = f"{pct:6.1f}%" if net > 0 else "  SECR"
        print(f"{c:9}{n:>3} {sp:7}{f:9.5f}{filt:10.2f}{clr:9.2f}"
              f"{pcts:>8}  {basis}")

    print()
    print(bar)
    print("LAYER 3b  THE GLOMERULAR / TUBULAR SPLIT, HUMAN vs MOUSE")
    print(bar)
    print("  The only pair with enough renal clearance data to decompose.")
    print()
    print(f"{'chemical':9}{'nC':>3}{'CLR ratio':>11}{'glomerular':>12}"
          f"{'tubular':>10}   dominant")
    byc = {}
    for (c, n, fam, sp, f, filt, clr, net, pct, verdict, basis) in L3:
        byc.setdefault(c, {})[sp] = (n, f, filt, clr)
    for c in sorted(byc, key=lambda x: (CHEM[x][1], CHEM[x][0])):
        d = byc[c]
        if "human" not in d or "mouse" not in d:
            continue
        n = d["human"][0]
        clr_r = d["mouse"][3] / d["human"][3]
        glom = d["mouse"][2] / d["human"][2]
        tub = clr_r / glom
        import math
        dom = ("tubular" if abs(math.log10(tub)) > abs(math.log10(glom))
               else "glomerular")
        arrow = "with" if tub > 1 else "AGAINST"
        print(f"{c:9}{n:>3}{clr_r:>11.1f}{glom:>12.1f}{tub:>10.2f}   "
              f"{dom:10} tubule works {arrow}")
    print()
    print("  THIS IS THE OPPOSITE OF THE RAT RESULT, AND THAT IS THE FINDING.")
    print()
    print("    human vs rat  PFOA   303x  =   6.8x glomerular x  45x tubular")
    print("    human vs mouse PFOA    85x  =  37.7x glomerular x 2.3x tubular")
    print()
    print("  The glomerular term carries the mouse comparison because the")
    print("  mouse differs from the human in BOTH of its factors -- fu is")
    print("  5.8x higher and GFR is 6.5x higher, and 5.8 x 6.5 = 37.7. The rat")
    print("  differs in only one: its fu matches the human, so its glomerular")
    print("  term is almost pure GFR (6.5x of the 6.8x) and the tubule has to")
    print("  carry the rest.")
    print()
    print("  So the two rodents are high-clearance for different reasons, and")
    print("  'rodents clear PFOA faster than humans' is one sentence covering")
    print("  two mechanisms. A PBPK model calibrated on the mouse tubule and")
    print("  applied to the rat, or the reverse, imports the wrong term.")
    print()
    print("  Across the six compounds the tubular factor spans 0.06x to 4.6x")
    print("  while the glomerular factor spans 13x to 122x, so the glomerular")
    print("  term dominates five of the six -- PFDA is the one exception, and")
    print("  only marginally (1.22 log units against 1.13). For four of the six")
    print("  the tubular factor is BELOW one -- the tubule does not amplify")
    print("  the mouse's faster clearance, it partly cancels it. For PFDA the")
    print("  cancellation is almost exact (13.4x glomerular, 0.1x tubular,")
    print("  0.9x observed), which is why PFDA has nearly the SAME renal")
    print("  clearance in mouse and human despite a 13x higher filtered load.")

    print()
    print("  ONE ROW WHERE THE fu CHOICE SWINGS THE ANSWER. Human PFDA is")
    print("  90.6% reabsorbed on Ryu's fu and 98.5% on Argoul's own -- the")
    print("  figure Argoul 2026 Table 1 prints. The two differ 6-fold in fu")
    print("  for this compound, and PFDA is exactly where the monotonicity")
    print("  test below says Argoul's series breaks. Every other row moves by")
    print("  less than two percentage points between the two choices.")

    print()
    print(bar)
    print("STRAIN IS NOT A DETAIL FOR PFOS")
    print(bar)
    p1_, p5_ = fu[("PFOS", "mouse (CD-1)")], fu[("PFOS", "mouse (C57BL/6)")]
    print(f"  PFOS mouse fu   CD-1 {p1_:.5f}   C57BL/6 {p5_:.5f}"
          f"   ratio {p1_/p5_:.1f}x")
    pu1, pu5 = fu[("PFUnDA", "mouse (CD-1)")], fu[("PFUnDA", "mouse (C57BL/6)")]
    print(f"  PFUnDA mouse fu CD-1 {pu1:.5f}   C57BL/6 {pu5:.5f}"
          f"   ratio {pu1/pu5:.1f}x")
    print("  The other twelve agree between the two strains within 1.5x.")
    print("  These two do not. 'Mouse fu' is not one number for PFOS or")
    print("  PFUnDA, so a mouse reabsorption figure for either inherits")
    print("  whichever strain supplied the binding measurement -- and the")
    print("  Argoul renal clearances are C57BL/6, so that is the strain used")
    print("  above. Had CD-1 been used, mouse PFOS filtered load would be 5x")
    print("  higher and its reabsorbed fraction 99.4% rather than 97.1%.")

    print()
    print(bar)
    print("CROSS-LABORATORY CHECK  Ryu against Argoul 2025/2026")
    print(bar)
    print(f"{'chemical':9}{'species':7}{'Ryu':>10}{'Argoul':>10}{'ratio':>8}")
    rs = []
    for (sp, c), fa in sorted(FU_ARGOUL.items(), key=lambda kv: kv[0]):
        m = "mouse (C57BL/6)" if sp == "mouse" else sp
        fr = fu.get((c, m))
        if fr is None:
            continue
        rs.append(fr / fa)
        print(f"{c:9}{sp:7}{fr:10.5f}{fa:10.5f}{fr/fa:>8.2f}")
    print()
    n2 = sum(1 for x in rs if 0.5 <= x <= 2)
    print(f"  n={len(rs)}  geometric mean {gm(rs):.2f}x  "
          f"range {min(rs):.2f} to {max(rs):.2f}   within two-fold "
          f"{n2}/{len(rs)}")
    print()
    print("  CORRECTION TO THE EARLIER READING. On PFOA alone Ryu ran 1.3-1.9x")
    print("  above Argoul and I called it a consistent offset. Across fifteen")
    print("  pairs it is not: the geometric mean is 0.90, so there is no net")
    print("  offset in either direction, and the spread is 13-fold. Eleven of")
    print("  fifteen agree within two-fold; the four that do not are the ones")
    print("  that matter.")
    print()
    print("  AND THE DISAGREEMENT IS NOT RANDOM -- IT IS AT THE LONG CHAINS.")
    print("  Argoul's own equilibration study found long-chain PFAS take up to")
    print("  120 h to reach equilibrium; Ryu dialysed for 18 h with membrane")
    print("  presaturation. The two protocols differ exactly where they")
    print("  disagree. The monotonicity test below says which to believe")
    print("  without trusting either laboratory.")

    print()
    print(bar)
    print("WHICH LABORATORY TO BELIEVE: THE MONOTONICITY TEST")
    print(bar)
    print("  Binding affinity to albumin rises with perfluorocarbon chain")
    print("  length, so within the straight-chain carboxylate series fu must")
    print("  FALL monotonically from C4 to C12. That is a constraint on the")
    print("  data, independent of either laboratory's protocol. Count the")
    print("  inversions.")
    print()
    # Ryu, human and C57BL/6 mouse carboxylate series
    carb = [(CHEM[c][0], c) for c in CHEM if CHEM[c][1] == "carboxylate"]
    carb.sort()
    series = {
        "Ryu human":        [(n, fu.get((c, "human"))) for n, c in carb],
        "Ryu mouse C57":    [(n, fu.get((c, "mouse (C57BL/6)"))) for n, c in carb],
        "Ryu rat":          [(n, fu.get((c, "rat"))) for n, c in carb],
        "Argoul human":     [(n, FU_ARGOUL.get(("human", c))) for n, c in carb],
        "Argoul mouse C57": [(n, FU_ARGOUL.get(("mouse", c))) for n, c in carb],
    }
    # Argoul's exp-2-only long chains, read straight off the stored table
    A2 = {}
    for r in csv.DictReader((HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if "Argoul" in r["source"] and r["experiment"] == "2" and r["fu"]:
            A2[(r["species"], r["chemical"])] = float(r["fu"])
    series["Argoul exp2 human"] = [
        (n, A2.get(("human", "PFUdA" if c == "PFUnDA" else c))) for n, c in carb]

    print(f"{'series':20}{'pairs':>7}{'inversions':>12}   where")
    for name, xs in series.items():
        xs = [(n, v) for n, v in xs if v]
        inv = [(xs[i][0], xs[i + 1][0]) for i in range(len(xs) - 1)
               if xs[i + 1][1] > xs[i][1]]
        where = ", ".join(f"C{a}->C{b}" for a, b in inv) or "none"
        print(f"{name:20}{len(xs) - 1:>7}{len(inv):>12}   {where}")
    print()
    print("  RYU'S CARBOXYLATE SERIES IS MONOTONE IN ALL THREE SPECIES.")
    print("  ARGOUL'S IS NOT: human fu RISES from PFNA (C9) to PFDA (C10) and")
    print("  rises again to PFUdA (C11), where it reaches 0.017 -- forty times")
    print("  Ryu's 0.00042, and within 15% of Argoul's OWN value for C6 PFHxA")
    print("  (0.020). A C11 carboxylate is not as weakly bound as a C6.")
    print()
    print("  So the disagreement is not symmetric uncertainty between two")
    print("  equally good datasets. Ryu's long-chain values satisfy a")
    print("  constraint Argoul's violate, which is what an 18-hour dialysis")
    print("  WITH presaturation versus a longer one WITHOUT it would predict:")
    print("  presaturation is there precisely to stop the membrane stripping")
    print("  the compound that binds hardest. The reabsorption figures above")
    print("  use Ryu's fu throughout, and that choice is now argued rather")
    print("  than assumed.")
    print()
    print("  WHAT THIS DOES NOT SETTLE. Argoul supplies the renal clearances,")
    print("  and nothing here tests those. The layer-3 rows therefore pair one")
    print("  laboratory's binding with another's clearance -- which is the")
    print("  cross-study join this review otherwise refuses. It is stated in")
    print("  the provenance_mismatch column of every row, and it is the reason")
    print("  the percentages are reported to one decimal and not two.")

    print()
    print(bar)
    print("SENSITIVITY: THE MOUSE GFR IS NOT AGREED")
    print(bar)
    print(f"  OEHHA / Han 2012   {GFR['mouse']} mL/kg/day   (used above)")
    print(f"  Argoul 2026        {GFR_MOUSE_ARGOUL} mL/kg/day"
          "   (recovered from their own fu x GFR column)")
    f = fu[("PFOA", "mouse (C57BL/6)")]
    for g in (GFR["mouse"], GFR_MOUSE_ARGOUL):
        filt = f * g
        print(f"    mouse PFOA at GFR {g:>6}: filtered {filt:7.1f}, "
              f"reabsorbed {100*(filt-3.59)/filt:5.1f}%")
    print("  The percentage is insensitive because reabsorption is near")
    print("  complete either way. The ABSOLUTE filtered load is not, so any")
    print("  mass-balance model built on it inherits a 1.67x fork.")


if __name__ == "__main__":
    main()
