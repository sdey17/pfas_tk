#!/usr/bin/env python3
"""The tubular factor T for every species, under one rule (Stage 4e).

    python3 -I scripts/stage4e_all_species_tubular.py

Writes db/tubular_factor_all_species.csv. Idempotent.

THE RULE, ENFORCED IN CODE AND NOT ONLY IN PROSE. For every row, the unbound
fraction, the glomerular filtration rate and the clearance must all come from
the SAME species. assert_same_species() refuses to build a row otherwise, so
species contamination cannot enter by oversight. Cross-STUDY joins within a
species are allowed; cross-species joins are not, in either direction.

    T  =  CL_R / (fu x GFR)        T<1 reabsorption of (1-T) of the filtered
                                   load; T=1 filtration alone; T>1 secretion

TWO TIERS.

  DIRECT   a study in that species reports renal clearance. T outright.
  BOUNDED  only total clearance is reported. CL_R <= CL_total always, so
           T <= CL_total/(fu x GFR). Where that bound is below one, net
           reabsorption holds whatever the renal share -- no assumption. Where
           it exceeds one it permits either sign and is reported as
           uninformative rather than dropped. Where several studies in the
           species disagree, the LOOSEST bound is used: the highest total
           clearance, so the conclusion holds whichever source is right.

WHAT THIS ADDS OVER STAGE 4D. Stage 4d did the rat. Abraham 2024 measured
human total clearance for fifteen PFAS in one volunteer, which takes the human
from 7 of 14 to 13 of 14. The mouse gains cross-checks but no new compounds.
The monkey is excluded, for two independent reasons given at the end.

STRAIN IS TRACKED SEPARATELY FROM SPECIES. Stage 4c found mouse fu for PFOS
differs five-fold between CD-1 and C57BL/6, so each mouse row records which
strain supplied the fu and which supplied the clearance, and flags a mismatch.
Strain mismatch is not species contamination and does not block a row; it is
reported.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

# ---- glomerular filtration rate, mL/kg/day, per species -----------------
GFR = {
    "human": (2570.0, "OEHHA PHG Table A6.4"),
    "rat":   (14400.0, "OEHHA PHG Table A6.4"),
    "mouse": (16700.0, "OEHHA PHG Table A6.4; Argoul 2026 uses 10000"),
}

# Ryu's matrices, by species
FU_MATRIX = {
    "human": ["human"],
    "rat":   ["rat"],
    "mouse": ["mouse (C57BL/6)", "mouse (CD-1)"],
}
STRAIN_OF_MATRIX = {"human": "", "rat": "Wistar-Han",
                    "mouse (C57BL/6)": "C57BL/6", "mouse (CD-1)": "CD-1"}

# Ryu's names -> this project's names
ALIAS = {"PFUdA": "PFUnDA", "PFDoA": "PFDoDA", "PFPA": "PFPeA",
         "PFHA": "PFHpA", "PFUDA": "PFUnDA", "6:2 FtS": "6:2 FTS"}

# ---- DIRECT renal clearances, mL/kg/day ---------------------------------
# (species, strain, sex, chemical, CL_R, source, note)
DIRECT = [
    ("human", "", "both", "PFBS", 1.2, "Argoul 2026 Table 1", ""),
    ("human", "", "both", "PFHxA", 47.0, "Argoul 2026 Table 1", ""),
    ("rat", "Wistar", "male", "PFOA", 21.2, "Ohmori 2003",
     "digitized from the paper by Han 2012; the paper reports CL_R for "
     "PFHpA, PFOA, PFNA and PFDA in both sexes but only this value has "
     "been extracted"),
    ("rat", "Sprague-Dawley", "male", "PFOA", 18.2,
     "Kemper 2003 (unpublished), as used by OEHHA PHG Table A6.4", ""),
    ("rat", "Wistar", "male", "PFOA", 0.032 * 1440, "Kudo 2002 Table 3",
     "0.032 mL/min/kg at 5 h after IV; 91% of its own total clearance"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 0.76, "Kim 2016",
     "8% of its own total clearance of 9.01"),
    ("mouse", "C57BL/6", "female", "PFBA", 783.0,
     "Argoul 2026 Table 4", ""),
    ("mouse", "C57BL/6", "female", "PFHxA", 6780.0,
     "Argoul 2026 Table 4", ""),
    ("mouse", "C57BL/6", "female", "PFOA", 3.59,
     "Argoul 2026 Table 4", ""),
    ("mouse", "C57BL/6", "female", "PFNA", 0.6,
     "Argoul 2026 Table 4", ""),
    ("mouse", "C57BL/6", "female", "PFDA", 0.13,
     "Argoul 2026 Table 4", ""),
    ("mouse", "C57BL/6", "female", "PFHxS", 0.589,
     "Argoul 2026 Table 4", ""),
    ("mouse", "C57BL/6", "female", "PFOS", 1.24,
     "Argoul 2026 Table 4", ""),
]

# ---- TOTAL clearances, mL/kg/day ----------------------------------------
TOTAL = [
    # human: Abraham 2024, one 67-year-old male volunteer, 2-compartment fit
    # to 450 days of plasma decay after a known high exposure.
    ("human", "", "male", "PFBA", 21.89, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFPeA", 146.6, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFHxA", 56.89, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFHpA", 0.584, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFOA", 0.0417, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFNA", 0.0659, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFDA", 0.148, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFUnDA", 0.286, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFDoDA", 0.832, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFBS", 1.877, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFHxS", 0.0530, "Abraham 2024", "n=1"),
    ("human", "", "male", "PFOS", 0.0870, "Abraham 2024", "n=1"),
    ("human", "", "male", "6:2 FTS", 19.92, "Abraham 2024", "n=1"),
    ("human", "", "both", "PFOA", 0.120, "US EPA 2024", ""),
    ("human", "", "both", "PFOA", 0.096, "Fujii 2015", ""),
    ("human", "", "both", "PFOS", 0.069, "ATSDR 2021", "derived"),
    ("human", "", "both", "PFNA", 0.154, "ATSDR 2021", "derived"),
    ("human", "", "male", "PFHxS", 0.064, "ATSDR 2021", "derived"),
    # rat
    ("rat", "Sprague-Dawley", "male", "PFOA", 47.39, "Kim 2016", "IV"),
    ("rat", "Wistar", "male", "PFOA", 50.4, "Kudo 2002", "IV"),
    ("rat", "Wistar", "male", "PFNA", 6.9,
     "Ohmori 2003 via ATSDR 2021 Table 3-6", ""),
    ("rat", "Wistar", "male", "PFDA", 207.0,
     "Ohmori 2003 via ATSDR 2021 Table 3-6",
     "ATSDR prints SD 0.054, implausible"),
    ("rat", "Wistar", "male", "PFHpA", 1604.0,
     "Ohmori 2003 via ATSDR 2021 Table 3-6", ""),
    ("rat", "Sprague-Dawley", "male", "PFOS", 9.24, "Kim 2016", "IV"),
    ("rat", "Sprague-Dawley", "male", "PFOS", 13.1,
     "Huang 2019 (NTP)", "beta phase"),
    ("rat", "Sprague-Dawley", "male", "PFOS", 11.3, "Chang 2012",
     "2 mg/kg"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 6.7, "Sundstrom 2012", ""),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 9.01, "Kim 2016", "IV"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 4.73,
     "Huang 2019 (NTP)", "beta phase"),
    ("rat", "Sprague-Dawley", "male", "PFBS", 828.0,
     "Huang 2019 (NTP)", ""),
    ("rat", "Sprague-Dawley", "male", "PFBS", 946.0, "Chengelis 2009", ""),
    ("rat", "Sprague-Dawley", "male", "PFBS", 2856.0, "Olsen 2009", ""),
    ("rat", "Sprague-Dawley", "male", "PFBA", 851.0, "Chang 2008", "IV"),
    ("rat", "Sprague-Dawley", "male", "PFHxA", 2784.0,
     "Chengelis 2009", ""),
    # mouse
    ("mouse", "CD-1", "male", "PFOS", 4.7, "Chang 2012",
     "1 and 20 mg/kg alike"),
    ("mouse", "CD-1", "male", "PFHxS", 4.8, "Sundstrom 2012",
     "20 mg/kg; 2.9 at 1 mg/kg"),
    ("mouse", "FVB/NJcl", "male", "PFOA", 14.2, "Fujii 2015", ""),
]

# ---- species that cannot be built, and why ------------------------------
EXCLUDED = {
    "monkey (cynomolgus)": [
        "no plasma unbound fraction has been measured in cynomolgus macaque "
        "plasma for any PFAS -- Ryu 2024 covers human, rat and two mouse "
        "strains only, and this project holds no monkey value from any source",
        "the only macaque GFR in this project is 8500 mL/kg/day for the "
        "JAPANESE macaque (Macaca fuscata, OEHHA Table A6.4), and the "
        "clearance data are cynomolgus (Macaca fascicularis). Borrowing the "
        "one for the other is species contamination and is refused.",
    ],
    "monkey (rhesus)": ["no clearance and no unbound fraction"],
}


def strain_note(fu_strain, cl_strain):
    """Wistar-Han is a Wistar-derived stock, so that pairing is a substrain
    difference and not a strain difference. Everything else that differs is
    reported as a strain mismatch. Neither blocks a row: strain is not
    species, and both are tracked separately."""
    if not fu_strain or not cl_strain or fu_strain == cl_strain:
        return ""
    fam = {"Wistar", "Wistar-Han"}
    if {fu_strain, cl_strain} <= fam:
        return f"substrain: fu {fu_strain} vs clearance {cl_strain}"
    return f"STRAIN: fu {fu_strain} vs clearance {cl_strain}"


def assert_same_species(label, *species_tags):
    """Refuse to build a row from more than one species. The hard rule."""
    uniq = {s for s in species_tags if s}
    if len(uniq) != 1:
        raise SystemExit(f"SPECIES CONTAMINATION in {label}: {sorted(uniq)}")
    return uniq.pop()


def ryu_fu():
    """fu keyed by (matrix, chemical). Matrix carries the species."""
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] != "Ryu 2024" or not r["fu"]:
            continue
        chem = ALIAS.get(r["chemical"], r["chemical"])
        out[(r["species"], chem)] = float(r["fu"])
    return out


def species_of_matrix(m):
    return m.split(" (")[0]


def build():
    fu = ryu_fu()
    rows = []
    for tier, table in (("DIRECT", DIRECT), ("BOUNDED", TOTAL)):
        for sp, strain, sex, chem, cl, src, note in table:
            chem = ALIAS.get(chem, chem)
            g, gsrc = GFR[sp]
            # pick the fu matrix for this species, preferring a strain match
            cands = [m for m in FU_MATRIX[sp] if (m, chem) in fu]
            if not cands:
                continue
            match = [m for m in cands if STRAIN_OF_MATRIX[m] == strain]
            m = (match or cands)[0]
            # THE HARD RULE
            assert_same_species(f"{sp}/{chem}/{src}",
                                sp, species_of_matrix(m))
            f = fu[(m, chem)]
            filt = f * g
            T = cl / filt
            rows.append(dict(
                species=sp, strain=strain, sex=sex, chemical=chem, tier=tier,
                fu=f, fu_matrix=m, fu_strain=STRAIN_OF_MATRIX[m],
                gfr=g, gfr_src=gsrc, filt=filt, cl=cl, src=src, note=note,
                T=T,
                strain_mismatch=strain_note(STRAIN_OF_MATRIX[m], strain)))
    return rows, fu


def verdict_of(r):
    if r["tier"] == "DIRECT":
        if r["T"] < 1:
            return "net_reabsorption", round(100 * (1 - r["T"]), 1)
        return "NET_TUBULAR_SECRETION", ""
    if r["T"] < 1:
        return "net_reabsorption_holds_on_upper_bound", round(
            100 * (1 - r["T"]), 1)
    return "UNINFORMATIVE_bound_exceeds_filtration", ""


def main():
    rows, fu = build()
    p = HERE / "db" / "tubular_factor_all_species.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "strain", "sex", "chemical", "tier",
                    "fu", "fu_matrix", "fu_source", "gfr_mL_kg_day",
                    "gfr_source", "filtered_load_mL_kg_day",
                    "clearance_mL_kg_day", "clearance_kind", "T", "verdict",
                    "min_pct_reabsorbed", "clearance_source",
                    "strain_mismatch", "note"])
        for r in sorted(rows, key=lambda x: (x["species"], x["chemical"],
                                             x["tier"])):
            v, mn = verdict_of(r)
            w.writerow([r["species"], r["strain"], r["sex"], r["chemical"],
                        r["tier"], r["fu"], r["fu_matrix"],
                        "Ryu 2024 Toxics 12:253", r["gfr"], r["gfr_src"],
                        round(r["filt"], 3), r["cl"],
                        "renal" if r["tier"] == "DIRECT" else
                        "total_upper_bound", round(r["T"], 5), v, mn,
                        r["src"], r["strain_mismatch"], r["note"]])
    report(rows, fu)
    print(f"\nwrote {p.relative_to(HERE)}  ({len(rows)} rows)")


def report(rows, fu):
    bar = "=" * 78
    CHEMS = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
             "PFUnDA", "PFDoDA", "PFBS", "PFHxS", "PFOS", "PFOSA", "6:2 FTS"]

    # ---- one summary per (species, chemical) ---------------------------
    # DIRECT rows: report the full range across sources. A verdict is given
    # only when every source agrees on the sign; otherwise DISPUTED, because
    # picking the most cautious source would manufacture a verdict the
    # literature does not support.
    # BOUNDED rows: use the HIGHEST clearance, the loosest bound, so a
    # reabsorption claim holds whichever source is right.
    cell = {}
    for r in rows:
        cell.setdefault((r["species"], r["chemical"]), []).append(r)

    def summarise(rs):
        d = [x for x in rs if x["tier"] == "DIRECT"]
        if d:
            Ts = sorted(x["T"] for x in d)
            lo, hi = Ts[0], Ts[-1]
            if hi < 1:
                return ("DIRECT", lo, hi, "net_reabsorption",
                        100 * (1 - hi))
            if lo > 1:
                return ("DIRECT", lo, hi, "NET_SECRETION", None)
            return ("DIRECT", lo, hi, "DISPUTED", None)
        w = max(rs, key=lambda x: x["cl"])
        if w["T"] < 1:
            return ("BOUNDED", w["T"], w["T"], "net_reabsorption_holds",
                    100 * (1 - w["T"]))
        return ("BOUNDED", w["T"], w["T"], "uninformative", None)

    summ = {k: summarise(v) for k, v in cell.items()}

    print(bar)
    print("T BY SPECIES AND COMPOUND")
    print(bar)
    print("  D = direct renal clearance; a range means several sources, and a")
    print("      verdict is given only where they agree on the sign.")
    print("  B = total clearance as an upper bound, on the HIGHEST source, so")
    print("      a reabsorption claim holds whichever source is right.")
    print("  ?   the bound permits either sign; nothing is concluded.")
    print("  !   sources disagree across the filtration boundary.")
    print()
    hdr = ("human", "rat", "mouse")
    print(f"{'chemical':9}" + "".join(f"{h:>24}" for h in hdr))
    for c in CHEMS:
        line = f"{c:9}"
        for sp in hdr:
            t = summ.get((sp, c))
            if t is None:
                line += f"{'-':>24}"
                continue
            tag, lo, hi, v, mn = t
            rng = f"{lo:.3f}" if abs(hi - lo) < 1e-9 else f"{lo:.2f}-{hi:.2f}"
            if v == "DISPUTED":
                body = f"{tag[0]} T={rng} !"
            elif v == "NET_SECRETION":
                body = f"{tag[0]} T={rng} SECR"
            elif v == "uninformative":
                body = f"{tag[0]} T={rng} ?"
            else:
                body = f"{tag[0]} T={rng} >={mn:.1f}%"
            line += f"{body:>24}"
        print(line)

    print()
    print(bar)
    print("COVERAGE")
    print(bar)
    print(f"  {'species':8}{'clearance':>11}{'conclusive':>12}"
          f"{'of which direct':>17}{'disputed':>10}{'uninformative':>15}")
    for sp in hdr:
        got = [c for c in CHEMS if (sp, c) in summ]
        con = [c for c in got if summ[(sp, c)][3] in
               ("net_reabsorption", "NET_SECRETION",
                "net_reabsorption_holds")]
        dr = [c for c in got if summ[(sp, c)][0] == "DIRECT"
              and summ[(sp, c)][3] != "DISPUTED"]
        dis = [c for c in got if summ[(sp, c)][3] == "DISPUTED"]
        un = [c for c in got if summ[(sp, c)][3] == "uninformative"]
        print(f"  {sp:8}{len(got):>8}/14{len(con):>9}/14{len(dr):>17}"
              f"{len(dis):>10}{len(un):>15}")
    print()
    print("  The human moves from 7 of 14 to 13 of 14 on Abraham 2024 alone,")
    print("  twelve of them conclusive. The rat moves from 1 to 9 with a")
    print("  clearance but only 4 conclusive, because its short-chain total")
    print("  clearances exceed the filtered load. The mouse stays at 7, all")
    print("  direct: the only mouse total clearances in the project are for")
    print("  compounds Argoul already measured renally, so they serve as")
    print("  cross-checks instead.")

    print()
    print(bar)
    print("THE HUMAN PFOA T IS ITSELF A RANGE, AND IT WIDENS THE RAT RATIO")
    print(bar)
    hp = sorted(x["cl"] for x in cell[("human", "PFOA")])
    filt_h = cell[("human", "PFOA")][0]["filt"]
    print(f"  human PFOA filtered load  {filt_h:.2f} mL/kg/day")
    for x in sorted(cell[("human", "PFOA")], key=lambda y: y["cl"]):
        print(f"    CL {x['cl']:>7.4f}  T {x['cl']/filt_h:.5f}"
              f"   {x['src']} ({x['tier'].lower()})")
    print(f"  plus Argoul 2026's 0.042 -> T 0.00667, used in sections 16-17")
    print()
    Th_lo, Th_hi = 0.042 / filt_h, max(hp) / filt_h
    print(f"  so T_human(PFOA) spans {Th_lo:.5f} to {Th_hi:.5f}, a "
          f"{Th_hi/Th_lo:.1f}-fold range across human sources alone")
    print()
    rr = [x["T"] for x in cell[("rat", "PFOA")] if x["tier"] == "DIRECT"]
    print(f"  T_rat(PFOA, male) spans {min(rr):.3f} to {max(rr):.3f} "
          f"(section 17)")
    print(f"  T_rat / T_human therefore spans "
          f"{min(rr)/Th_hi:.0f}x to {max(rr)/Th_lo:.0f}x")
    print()
    print("  SECTIONS 16 AND 17 SAID 45x AND THEN 45-160x. WITH EVERY HUMAN")
    print("  SOURCE IN THE PROJECT IT IS WIDER STILL. The low end moves")
    print("  because the US EPA and Fujii human clearances are two to three")
    print("  times Argoul's. Every combination still leaves the tubule")
    print("  dominant over the 6.8x glomerular factor, by a factor of at")
    print("  least three, so the conclusion is unchanged and the number")
    print("  should be quoted as an order of magnitude, not a figure.")

    print()
    print(bar)
    print("THE MOUSE CROSS-CHECKS PASS")
    print(bar)
    print("  A total clearance must not fall BELOW the renal clearance of the")
    print("  same compound and species. Three mouse pairs can be tested:")
    print()
    for c in ("PFOS", "PFHxS", "PFOA"):
        d = [r for r in rows if r["species"] == "mouse" and r["chemical"] == c
             and r["tier"] == "DIRECT"]
        b = [r for r in rows if r["species"] == "mouse" and r["chemical"] == c
             and r["tier"] == "BOUNDED"]
        if not d or not b:
            continue
        ok = "PASS" if b[0]["cl"] >= d[0]["cl"] else "FAIL"
        print(f"  {c:7} renal {d[0]['cl']:>7.2f} ({d[0]['strain']}, "
              f"Argoul)   total {b[0]['cl']:>6.1f} ({b[0]['strain']}, "
              f"{b[0]['src']})   {ok}")
    print()
    print("  All three pass, and the strains differ across each pair, so this")
    print("  is also a weak check that mouse strain does not dominate the")
    print("  renal numbers for these compounds. For PFOS it is the compound")
    print("  where Stage 4c found a 5-fold STRAIN difference in fu, so the")
    print("  check is worth having: CD-1 total 4.7 against C57BL/6 renal")
    print("  1.24 is consistent, which it would not have to be.")

    print()
    print(bar)
    print("ONE HUMAN COMPOUND THE RENAL ACCOUNTING CANNOT HOLD")
    print(bar)
    r = max(cell[("human", "PFDoDA")], key=lambda x: x["cl"])
    print(f"  PFDoDA  fu {r['fu']:.6f}  filtered {r['filt']:.3f} mL/kg/day")
    print(f"          total clearance {r['cl']} mL/kg/day, which is "
          f"{r['T']:.1f}x the filtered load")
    print()
    print("  So human PFDoDA elimination CANNOT be mostly glomerular. Either")
    print("  the tubule secretes it, or -- far likelier for a C12 carboxylate")
    print("  -- most of its elimination is not renal at all. Kudo 2001 found")
    print("  faecal elimination dominant for PFDA and above in the rat, and")
    print("  section 10's biliary row says the same for the human. THE T")
    print("  FRAMEWORK IS A RENAL ACCOUNTING AND IT DOES NOT DESCRIBE THE")
    print("  LONGEST CHAINS. That is a limit of the framework, not a gap in")
    print("  the data, and PFUnDA at T <= 0.27 is already close to the edge.")

    print()
    print(bar)
    print("SPECIES EXCLUDED, AND WHY")
    print(bar)
    for sp, reasons in EXCLUDED.items():
        print(f"  {sp}")
        for why in reasons:
            print(f"    - {why}")
        print()
    print("  This project holds 22 cynomolgus rows with 14 clearances, so the")
    print("  monkey is excluded by the UNBOUND FRACTION, not by the")
    print("  clearance. One equilibrium-dialysis experiment on cynomolgus")
    print("  plasma would open a third species for the whole framework, and")
    print("  the monkey is the species whose PFOA half-life is closest to the")
    print("  human. THAT IS THE HIGHEST-VALUE SINGLE EXPERIMENT THIS STRAND")
    print("  HAS IDENTIFIED.")

    print()
    print(bar)
    print("WHAT OHMORI 2003 WOULD CLOSE")
    print(bar)
    print("  Ohmori K, Kudo N, Katayama K, Kawashima Y (2003). Toxicology")
    print("  184(2-3):135-140. PMID 12499116.")
    print()
    print("  The abstract states that renal clearance WAS DETERMINED for")
    print("  PFHpA, PFOA, PFNA and PFDA in male and female Wistar rats, and")
    print("  reports r^2 = 0.981 between total and renal clearance. Only the")
    print("  PFOA value has reached this project, through Han 2012's")
    print("  digitization. The paper's own numbers would:")
    print("    - close rat PFHpA and PFDA, both currently uninformative")
    print("    - replace the rat PFNA bound with a direct value")
    print("    - give the first female-rat rows in the table")
    print("    - come from WISTAR rats, matching Ryu's Wistar-Han fu better")
    print("      than the Sprague-Dawley sources do")
    print()
    print("  AND IT SETTLES AN ERROR. ATSDR 2021 Table 3-6 attributes a male")
    print("  rat PFOA total clearance of 135 mL/kg/day to Ohmori, but")
    print("  Ohmori's own Vss of 196 mL/kg and half-life of 5.63 d give")
    print("  ln2 x Vd / t = 24.1. Han's digitized CL_R of 21.2 is 88% of")
    print("  that, exactly as the paper's r^2 = 0.981 implies. ATSDR's 135 is")
    print("  5.6x too high and is almost certainly a transcription or unit")
    print("  error. The female value, 2233 against 1741 derived, is fine.")
    print()
    print("  Nothing above depends on the 135: it was used only as a loosest")
    print("  bound, where it is uninformative either way.")


if __name__ == "__main__":
    main()
