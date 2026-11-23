#!/usr/bin/env python3
"""The human unbound fraction, resolved (Stage 9).

    python3 -I scripts/stage9_human_fu_resolution.py

Writes db/smeltz2023_fu.csv, db/human_fu_comparison.csv and
db/human_tubular_sign_by_fu_source.csv. Idempotent.

    Smeltz M, Wambaugh JF, Wetmore BA (2023) Chem Res Toxicol 36:870-881,
    doi:10.1021/acs.chemrestox.3c00003
    Abraham K et al. (2024) Environ Int 193:109047,
    doi:10.1016/j.envint.2024.109047

1. THE ABRAHAM/SMELTZ IDENTITY IS RESOLVED, AND IT WAS NOT A FABRICATION.
   Abraham 2024 Table 4 carries the footnote "* as published by Smeltz et al.
   (2023)", and its Methods state that filtration clearance "was calculated as
   product of the fraction unbound in plasma (fup, available from Smeltz et
   al. 2023) and the GFR measured in" the volunteer.

   SO ABRAHAM IS NOT AN INDEPENDENT BINDING SERIES. Its clearances and its
   GFR are the subject's own; its binding is EPA's. This review's claim that
   Abraham's T "needs no cross-study join" was WRONG and is corrected here.
   There are FOUR independent human binding datasets, not five:

       Smeltz 2023    ultracentrifugation,            67 PFAS   (US EPA)
       Ryu 2024       presaturation equilibrium dialysis, 14 PFAS
       Argoul 2025    Dianorm equilibrium dialysis,   16 PFAS
       Fischer 2024   C18 solid-phase microextraction, 29 PFAS

2. AND ABRAHAM MISQUOTES SMELTZ IN ITS OWN DISCUSSION. The text says the PFCA
   values "ranged from 0.3108 for PFBA". 0.3108 is Smeltz's value for PFPA,
   perfluoropropanoate, C3. PFBA is 0.0930, which is what Abraham's Table 4
   correctly uses. A text error only; the arithmetic is unaffected.

3. SMELTZ HAS AN INTERNAL DISCREPANCY OF ITS OWN. PFBS is 0.0046 in Table 2
   and 0.0128 in Table 3's UC column, a 2.8-fold difference in the same assay
   in the same paper. PFOS is 0.0045 against 0.0049, which is rounding (Table
   2 states it is an average of three salt forms). Abraham used the Table 2
   values. Both are recorded.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
SRC = ("Smeltz 2023, Chem Res Toxicol 36:870-881, "
       "doi:10.1021/acs.chemrestox.3c00003")

# ---- Smeltz Table 2, ultracentrifugation (the method the paper prefers) --
# chemical, n_carbons, family, fup
T2 = [
    ("PFPA", 3, "carboxylate", 0.3108), ("PFBA", 4, "carboxylate", 0.0930),
    ("PFPeA", 5, "carboxylate", 0.0440), ("PFHxA", 6, "carboxylate", 0.0068),
    ("PFHpA", 7, "carboxylate", 0.0004), ("PFOA", 8, "carboxylate", 0.0010),
    ("PFNA", 9, "carboxylate", 0.0016), ("PFDA", 10, "carboxylate", 0.0027),
    ("PFUnDA", 11, "carboxylate", 0.0148),
    ("PFTrDA", 13, "carboxylate", 0.0433),
    ("PFTeDA", 14, "carboxylate", 0.0271),
    ("PFBS", 4, "sulfonate", 0.0046), ("PFHxS", 6, "sulfonate", 0.0009),
    ("PFHpS", 7, "sulfonate", 0.0006), ("PFOS", 8, "sulfonate", 0.0045),
    ("PFDS", 10, "sulfonate", 0.0528),
    ("4:2FTS", 6, "fluorotelomer sulfonate", 0.0050),
    ("6:2FTS", 8, "fluorotelomer sulfonate", 0.0143),
    ("8:2FTS", 10, "fluorotelomer sulfonate", 0.0396),
    ("3:3FTCA", 6, "fluorotelomer carboxylate", 0.0518),
    ("5:3FTCA", 8, "fluorotelomer carboxylate", 0.0059),
    ("6:3FTCA", 9, "fluorotelomer carboxylate", 0.0023),
    ("7:3FTCA", 10, "fluorotelomer carboxylate", 0.0051),
    ("4H-PFBA", 4, "polyfluoro carboxylate", 0.1714),
    ("7H-PFHpA", 7, "polyfluoro carboxylate", 0.0026),
    ("8H-PFOA", 8, "polyfluoro carboxylate", 0.0021),
    ("9H-PFNA", 9, "polyfluoro carboxylate", 0.0009),
    ("11H-PFUnDA", 11, "polyfluoro carboxylate", 0.0015),
]
# Table 3: UC vs RED, and whether RED reached equilibrium
T3 = {
    "PFBA": (0.0930, 0.1453, True), "PFPeA": (0.0440, 0.0455, True),
    "PFHxA": (0.0068, 0.0106, True), "PFHpA": (0.0004, 0.0013, True),
    "PFOA": (0.0010, 0.0011, True), "PFNA": (0.0016, 0.0004, True),
    "PFDA": (0.0027, 0.0001, True), "PFUnDA": (0.0148, 0.0001, False),
    "PFTrDA": (0.0433, 0.000048, False), "PFTeDA": (0.0271, 0.0016, False),
    "PFBS": (0.0128, 0.0121, True), "PFHxS": (0.0009, 0.0004, False),
    "PFHpS": (0.0006, None, True), "PFOS": (0.0049, 0.0003, True),
    "PFOSA": (0.0039, 0.0008, True), "PFDS": (0.0528, None, False),
    "4:2FTS": (0.0050, 0.0079, True), "8:2FTS": (0.0396, 0.0012, True),
    "3:3FTCA": (0.0518, 0.0761, True),
}

# ---- Abraham 2024 Table 4: the subject's own clearances, mL/min ----------
GFR_SUBJ = 113.7           # mL/min, the volunteer's measured creatinine cl.
ABRAHAM = {   # chemical: (CL_uri, CL_tot, CL_fec or None)
    "PFBA": (1.00, 1.24, None), "PFPeA": (8.89, 8.45, None),
    "PFHxA": (2.68, 3.24, None), "PFHpA": (0.023, 0.033, None),
    "PFOA": (None, 0.0024, 0.0004), "PFNA": (None, 0.0038, 0.0022),
    "PFDA": (None, 0.0084, 0.0051), "PFUnDA": (None, 0.0163, 0.0098),
    "PFBS": (0.067, 0.107, None), "PFHxS": (None, 0.0030, None),
    "PFOS": (None, 0.0050, None), "6:2FTS": (0.97, 1.13, None),
}

# ---- the four independent human binding datasets ------------------------
RYU, ARGOUL, FISCHER = {}, {}, {}


def load_others():
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["species"] != "human" or not r["fu"]:
            continue
        c, v = r["chemical"], float(r["fu"])
        if r["experiment"] == "Ryu 2024":
            RYU[c] = v
        elif "Argoul" in r["source"] and r["experiment"] == "1":
            ARGOUL[c] = v
    # Fischer: the derived whole-serum series
    p = HERE / "db" / "fischer2024_fu_derived.csv"
    if p.exists():
        for r in csv.DictReader(p.open()):
            FISCHER[r["chemical"]] = float(r["fu_derived"])


def main():
    load_others()
    smeltz = {c: v for c, _n, _f, v in T2}

    # ---------- 1. Smeltz as a record ----------
    p1 = HERE / "db" / "smeltz2023_fu.csv"
    with p1.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "family", "fup_ultracentrifugation",
                    "fup_RED", "RED_reached_equilibrium",
                    "table3_UC_if_different", "note", "source"])
        for c, n, fam, v in T2:
            t3 = T3.get(c)
            diff = ""
            note = ""
            if t3 and abs(t3[0] - v) > 1e-9:
                diff = t3[0]
                note = (f"Table 2 gives {v}, Table 3 gives {t3[0]} for the "
                        f"same UC assay, a {max(t3[0],v)/min(t3[0],v):.1f}x "
                        f"internal discrepancy")
            if c == "PFOS":
                note += ("; Table 2 value is an average of three salt forms"
                         if not note else
                         "; Table 2 value is an average of three salt forms")
            w.writerow([c, n, fam, v,
                        t3[1] if t3 and t3[1] is not None else "",
                        ("yes" if t3[2] else "no") if t3 else "",
                        diff, note, SRC + ", Tables 2 and 3"])

    # ---------- 2. the four-source comparison ----------
    chems = [c for c, _n, _f, _v in T2
             if c in RYU or c in ARGOUL or c in FISCHER]
    p2 = HERE / "db" / "human_fu_comparison.csv"
    with p2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "family", "smeltz_UC", "ryu",
                    "argoul", "fischer_derived", "n_sources", "min", "max",
                    "fold_spread"])
        comp = []
        for c, n, fam, v in T2:
            vals = {"smeltz": smeltz.get(c), "ryu": RYU.get(c),
                    "argoul": ARGOUL.get(c), "fischer": FISCHER.get(c)}
            present = [x for x in vals.values() if x]
            if len(present) < 2:
                continue
            w.writerow([c, n, fam, vals["smeltz"], vals["ryu"],
                        vals["argoul"], vals["fischer"], len(present),
                        min(present), max(present),
                        round(max(present) / min(present), 1)])
            comp.append((c, n, fam, vals, max(present) / min(present)))

    # ---------- 3. the sign of T, by binding source ----------
    p3 = HERE / "db" / "human_tubular_sign_by_fu_source.csv"
    signs = []
    with p3.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "cl_renal_mL_min",
                    "cl_renal_basis", "fu_source", "fu", "cl_filt_mL_min",
                    "T", "verdict"])
        for c, (uri, tot, fec) in ABRAHAM.items():
            if uri is not None:
                clr, basis = uri, "measured urinary clearance"
            elif fec is not None:
                clr, basis = tot - fec, "total minus measured faecal (bound)"
            else:
                clr, basis = tot, "total clearance (bound)"
            for name, d in (("Smeltz 2023 (used by Abraham)", smeltz),
                            ("Ryu 2024", RYU), ("Argoul 2025", ARGOUL),
                            ("Fischer 2024 (derived)", FISCHER)):
                fu = d.get(c)
                if not fu:
                    continue
                filt = fu * GFR_SUBJ
                T = clr / filt
                w.writerow([c, dict((x, n) for x, n, _f, _v in T2).get(c, ""),
                            round(clr, 5), basis, name, fu, round(filt, 4),
                            round(T, 4),
                            "net secretion" if T > 1 else "net reabsorption"])
                signs.append((c, name, T, basis))

    report(comp, signs, smeltz)
    for p in (p1, p2, p3):
        print(f"wrote {p.relative_to(HERE)}")


def report(comp, signs, smeltz):
    bar = "=" * 78
    print(bar)
    print("A.  THE CONTRADICTION, AND EXACTLY WHERE ABRAHAM REPORTS IT")
    print(bar)
    print("  Abraham 2024 section 4.5 'Renal elimination', page 10, and")
    print("  Fig. 7. The paper states it in its own words:")
    print()
    print('    "Remarkably, for MPFPeA and MPFHxA (shortest half-lives of')
    print('     0.5 and 1.4 days, respectively), the total clearance was')
    print('     higher than the filtration clearance (MPFPeA: 8.5 vs. 5.0')
    print('     mL/min; MPFHxA: 3.2 vs. 0.77 mL/min), INDICATING AN')
    print('     ADDITIONAL (NETTO) TUBULAR SECRETION."')
    print()
    print('    "In contrast, for the small molecules with longer half-lives,')
    print('     MPFBA and PFBS, the total clearance was much lower than the')
    print('     filtration clearance (MPFBA 1.2 vs. 10.6 mL/min; PFBS: 0.11')
    print('     vs. 0.52), INDICATING A CONSIDERABLE (NETTO) TUBULAR')
    print('     REABSORPTION."')
    print()
    print("  Fig. 7 plots CLFilt against CLtot for those four compounds with")
    print("  net secretion in green and net reabsorption in red. So the")
    print("  contradiction is not between the papers -- IT IS WITHIN ONE")
    print("  PAPER, ACROSS COMPOUNDS, AND ABRAHAM PUTS IT IN A FIGURE.")
    print()
    print("  Abraham also flags its own fragility, in the same paragraph:")
    print('    "these considerations are based on the data of Smeltz et al.')
    print('     (2023; NOT YET CONFIRMED BY OTHER STUDIES) and the')
    print('     investigation of a single person, and therefore have the')
    print('     CHARACTER OF A HYPOTHESIS."')
    print()
    print("  That sentence is the honest frame for the whole question, and it")
    print("  is the paper's own.")

    print()
    print(bar)
    print("B.  THE SIGN OF T, BY WHICH BINDING SOURCE IS USED")
    print(bar)
    print("  On Abraham's own clearances and the subject's own GFR of 113.7")
    print("  mL/min, so only the binding term changes between rows.")
    print()
    byc = {}
    for c, name, T, basis in signs:
        byc.setdefault(c, []).append((name, T))
    order = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
             "PFUnDA", "PFBS", "PFHxS", "PFOS", "6:2FTS"]
    print(f"  {'chem':8}{'Smeltz':>9}{'Fischer':>9}{'Argoul':>9}{'Ryu':>9}"
          f"   sign agreed?")
    for c in order:
        if c not in byc:
            continue
        d = dict((n.split(" ")[0], T) for n, T in byc[c])
        cells = [d.get(k) for k in ("Smeltz", "Fischer", "Argoul", "Ryu")]
        got = [x for x in cells if x is not None]
        agreed = all((x > 1) == (got[0] > 1) for x in got)
        row = "".join(f"{(f'{x:.2f}' if x is not None else '-'):>9}"
                      for x in cells)
        verdict = ("yes, " + ("SECRETION" if got[0] > 1 else "reabsorption")
                   if agreed else "NO -- sign depends on the source")
        print(f"  {c:8}{row}   {verdict}")
    print()
    print("  THE DISAGREEMENT IS CONFINED TO TWO COMPOUNDS, BOTH SHORT:")
    print("  PFPeA (C5) and PFHxA (C6). Nothing else in the series.")
    print("  Every compound of C7 and above agrees on the sign whichever")
    print("  binding dataset is used, and every one of them is reabsorption.")
    print("  So the choice of fu source does NOT affect any conclusion about")
    print("  PFOA, PFNA, PFDA, PFUnDA, PFHxS or PFOS.")
    print()
    print("  AND FOR PFHxA THREE OF THE FOUR NOW SAY SECRETION. Section 7")
    print("  reported this as 0.48 against 3.48, a two-way tie with Ryu's")
    print("  reading given first. On the subject's own GFR it is 3.47")
    print("  (Smeltz), 2.48 (Fischer), 1.18 (Argoul), 0.62 (Ryu). Ryu is the")
    print("  outlier, not the default.")

    print()
    print(bar)
    print("C.  WHERE THE FOUR DATASETS AGREE AND DISAGREE")
    print(bar)
    print(f"  {'chem':9}{'nC':>3}{'Smeltz':>9}{'Ryu':>9}{'Argoul':>9}"
          f"{'Fischer':>9}{'spread':>8}")
    for c, n, fam, vals, sp in comp:
        if fam not in ("carboxylate", "sulfonate"):
            continue
        row = "".join(
            f"{(f'{vals[k]:.4f}' if vals[k] else '-'):>9}"
            for k in ("smeltz", "ryu", "argoul", "fischer"))
        print(f"  {c:9}{n:>3}{row}{sp:>7.1f}x")
    print()
    print("  Agreement is good for the long chains and poor for the short")
    print("  ones -- the opposite of what methodological intuition suggests,")
    print("  and the reason the disputed compounds are the short ones.")

    print()
    print(bar)
    print("D.  THE U-SHAPE IS NOW REPORTED BY THREE LABORATORIES")
    print(bar)
    print("  Smeltz Table 2, perfluorinated carboxylates, ultracentrifugation:")
    print()
    print(f"    {'nC':>4}{'compound':>10}{'fup':>10}")
    for c, n, fam, v in T2:
        if fam != "carboxylate" or c.startswith(("4H", "7H", "8H", "9H", "11H")):
            continue
        mark = "  <-- minimum" if c == "PFHpA" else ""
        print(f"    {n:>4}{c:>10}{v:>10.4f}{mark}")
    print()
    print("  Smeltz calls it: 'The U-shaped trend noted between molecular")
    print("  weight and PPB (Table 2)', and italicises the PFCAs of at least")
    print("  11 carbons 'that exhibit reduced binding'. Fischer 2024 reports")
    print("  the same turn with a mechanism -- reduced molecular flexibility")
    print("  above 8 perfluorinated carbons blocks access to the albumin")
    print("  pocket -- and Argoul 2025 reports 'a U-shaped relationship")
    print("  between fu and molecular weight ... with the lowest values near")
    print("  400-500 g/mol'.")
    print()
    print("  THREE INDEPENDENT LABORATORIES, THREE METHODS, SAME U. The")
    print("  monotonic-binding assumption this review used in section 7 is")
    print("  refuted by all three, not just by Fischer.")


if __name__ == "__main__":
    main()
