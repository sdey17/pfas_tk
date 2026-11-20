#!/usr/bin/env python3
"""One authoritative tubular-factor table, three tiers (Stage 6, part 1).

    python3 -I scripts/stage6_consolidate_tubular.py

Writes db/tubular_factor_evidence.csv (long form, every source) and
db/tubular_factor_consolidated.csv (one verdict per species, sex, compound).
Idempotent. SUPERSEDES db/species_reabsorption_all14.csv,
db/rat_tubular_factor.csv, db/rat_tubular_factor_ohmori.csv and
db/tubular_factor_all_species.csv, which were built stage by stage and now
disagree with each other.

WHY THIS EXISTS. Section 19 read Ohmori 2003 from the paper and found ATSDR's
rat PFDA total clearance wrong by forty-fold. The Stage 4e table was built on
ATSDR's value and still says rat PFDA is uninformative, where the true value
makes it certain reabsorption. Two of this project's own files disagree. A
stage-by-stage record is right for the roadmap and wrong for the database, so
the database gets one table.

    T  =  CL_R / (fu x GFR)      T<1 net reabsorption of (1-T) of the
                                 filtered load; T=1 filtration alone;
                                 T>1 net secretion

THREE TIERS, EACH A DIFFERENT KIND OF EVIDENCE.

  DIRECT   a study in that species reports renal clearance. T outright. With
           several sources, the range; a verdict only where they agree on the
           sign, else DISPUTED.
  UPPER    only total clearance. CL_R <= CL_total, so T <= CL_tot/(fu x GFR).
           Below one, net reabsorption holds whatever the renal share. Built
           on the HIGHEST total clearance so it holds whichever source is
           right.
  LOWER    a measured renal FRACTION. T = f_ren x CL_tot/(fu x GFR), so
           secretion is established once f_ren > (fu x GFR)/CL_tot.

TWO RULES THE BUILD ENFORCES RATHER THAN ASSERTS.

  SPECIES PURITY. assert_same_species() refuses any row whose fu, GFR and
  clearance are not all from one species. Tested by deliberate violation.

  PRIMARY OVER SECONDARY. Where a value exists both in a paper and in a
  compilation of that paper, the paper wins and the compilation's value is
  recorded as superseded with the discrepancy. This is what Stage 4e got
  wrong.
"""
import csv, math, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

GFR = {"human": 2570.0, "rat": 14400.0, "mouse": 16700.0}
GFR_SRC = "OEHHA PHG Table A6.4 / Han 2012"
GFR_ALT = {"rat": (7920.0, "Ohmori 2003 citing Shim 1984, 5.5 mL/min/kg"),
           "mouse": (10000.0, "Argoul 2026, recovered from its fu x GFR column")}

FU_MATRIX = {"human": ["human"], "rat": ["rat"],
             "mouse": ["mouse (C57BL/6)", "mouse (CD-1)"]}
STRAIN_OF = {"human": "", "rat": "Wistar-Han",
             "mouse (C57BL/6)": "C57BL/6", "mouse (CD-1)": "CD-1"}
ALIAS = {"PFUdA": "PFUnDA", "PFDoA": "PFDoDA", "PFHA": "PFHpA"}

CHEMS = [("PFBA", 4, "carboxylate"), ("PFPeA", 5, "carboxylate"),
         ("PFHxA", 6, "carboxylate"), ("PFHpA", 7, "carboxylate"),
         ("PFOA", 8, "carboxylate"), ("PFNA", 9, "carboxylate"),
         ("PFDA", 10, "carboxylate"), ("PFUnDA", 11, "carboxylate"),
         ("PFDoDA", 12, "carboxylate"), ("PFBS", 4, "sulfonate"),
         ("PFHxS", 6, "sulfonate"), ("PFOS", 8, "sulfonate"),
         ("PFOSA", 8, "sulfonamide"), ("6:2 FTS", 6, "fluorotelomer")]
NC = {c: n for c, n, _f in CHEMS}
FAM = {c: f for c, _n, f in CHEMS}

# ---- renal clearances, mL/kg/day: (species, strain, sex, chem, val, src) --
RENAL = [
    ("human", "", "both", "PFBS", 1.2, "Argoul 2026 Table 1"),
    ("human", "", "both", "PFHxA", 47.0, "Argoul 2026 Table 1"),
    ("rat", "Sprague-Dawley", "male", "PFOA", 18.2,
     "Kemper 2003 (unpublished) via OEHHA PHG Table A6.4"),
    ("rat", "Wistar", "male", "PFOA", 21.2,
     "Ohmori 2003 Fig 2, digitized by Han 2012"),
    ("rat", "Wistar", "male", "PFOA", 46.08, "Kudo 2002 Table 3"),
    ("rat", "Sprague-Dawley", "female", "PFOA", 666.0,
     "Kemper 2003 (unpublished) via OEHHA PHG Table A6.4"),
    ("rat", "Wistar", "female", "PFOA", 1009.0,
     "Ohmori 2003 Fig 2, digitized by Han 2012"),
    ("rat", "Wistar", "female", "PFOA", 1054.1, "Kudo 2002 Table 3"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 0.76, "Kim 2016"),
    ("mouse", "C57BL/6", "female", "PFBA", 783.0, "Argoul 2026 Table 4"),
    ("mouse", "C57BL/6", "female", "PFHxA", 6780.0, "Argoul 2026 Table 4"),
    ("mouse", "C57BL/6", "female", "PFOA", 3.59, "Argoul 2026 Table 4"),
    ("mouse", "C57BL/6", "female", "PFNA", 0.6, "Argoul 2026 Table 4"),
    ("mouse", "C57BL/6", "female", "PFDA", 0.13, "Argoul 2026 Table 4"),
    ("mouse", "C57BL/6", "female", "PFHxS", 0.589, "Argoul 2026 Table 4"),
    ("mouse", "C57BL/6", "female", "PFOS", 1.24, "Argoul 2026 Table 4"),
]

# ---- total clearances, mL/kg/day ----------------------------------------
TOTAL = [
    ("human", "", "male", "PFBA", 21.89, "Abraham 2024"),
    ("human", "", "male", "PFPeA", 146.6, "Abraham 2024"),
    ("human", "", "male", "PFHxA", 56.89, "Abraham 2024"),
    ("human", "", "male", "PFHpA", 0.584, "Abraham 2024"),
    ("human", "", "male", "PFOA", 0.0417, "Abraham 2024"),
    ("human", "", "male", "PFNA", 0.0659, "Abraham 2024"),
    ("human", "", "male", "PFDA", 0.148, "Abraham 2024"),
    ("human", "", "male", "PFUnDA", 0.286, "Abraham 2024"),
    ("human", "", "male", "PFDoDA", 0.832, "Abraham 2024"),
    ("human", "", "male", "PFBS", 1.877, "Abraham 2024"),
    ("human", "", "male", "PFHxS", 0.0530, "Abraham 2024"),
    ("human", "", "male", "PFOS", 0.0870, "Abraham 2024"),
    ("human", "", "male", "6:2 FTS", 19.92, "Abraham 2024"),
    ("human", "", "both", "PFOA", 0.120, "US EPA 2024"),
    ("human", "", "both", "PFOA", 0.096, "Fujii 2015"),
    ("human", "", "both", "PFOA", 0.042, "Argoul 2026 Table 1"),
    ("human", "", "both", "PFOS", 0.069, "ATSDR 2021 (derived)"),
    ("human", "", "both", "PFOS", 0.088, "Argoul 2026 Table 1"),
    ("human", "", "both", "PFNA", 0.154, "ATSDR 2021 (derived)"),
    ("human", "", "both", "PFNA", 0.067, "Argoul 2026 Table 1"),
    ("human", "", "male", "PFHxS", 0.064, "ATSDR 2021 (derived)"),
    ("human", "", "both", "PFHxS", 0.053, "Argoul 2026 Table 1"),
    ("human", "", "both", "PFDA", 0.15, "Argoul 2026 Table 1"),
    # rat: Ohmori's own table, PRIMARY
    ("rat", "Wistar", "male", "PFHpA", 1603.9, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "female", "PFHpA", 3070.5, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "male", "PFOA", 50.5, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "female", "PFOA", 2233.2, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "male", "PFNA", 6.9, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "female", "PFNA", 105.7, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "male", "PFDA", 5.2, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "female", "PFDA", 5.3, "Ohmori 2003 Table 1"),
    ("rat", "Wistar", "male", "PFOA", 50.4, "Kudo 2002"),
    ("rat", "Sprague-Dawley", "male", "PFOA", 47.39, "Kim 2016"),
    ("rat", "Sprague-Dawley", "male", "PFOS", 9.24, "Kim 2016"),
    ("rat", "Sprague-Dawley", "male", "PFOS", 13.1, "Huang 2019 (NTP)"),
    ("rat", "Sprague-Dawley", "male", "PFOS", 11.3, "Chang 2012"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 6.7, "Sundstrom 2012"),
    ("rat", "Sprague-Dawley", "female", "PFHxS", 53.35, "Sundstrom 2012"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 9.01, "Kim 2016"),
    ("rat", "Sprague-Dawley", "male", "PFHxS", 4.73, "Huang 2019 (NTP)"),
    ("rat", "Sprague-Dawley", "male", "PFBS", 828.0, "Huang 2019 (NTP)"),
    ("rat", "Sprague-Dawley", "male", "PFBS", 946.0, "Chengelis 2009"),
    ("rat", "Sprague-Dawley", "male", "PFBS", 2856.0, "Olsen 2009"),
    ("rat", "Sprague-Dawley", "male", "PFBA", 851.0, "Chang 2008"),
    ("rat", "Sprague-Dawley", "male", "PFHxA", 2784.0, "Chengelis 2009"),
    # mouse
    ("mouse", "CD-1", "male", "PFOS", 4.7, "Chang 2012"),
    ("mouse", "CD-1", "male", "PFHxS", 4.8, "Sundstrom 2012"),
    ("mouse", "FVB/NJcl", "male", "PFOA", 14.2, "Fujii 2015"),
    ("mouse", "C57BL/6", "female", "PFBA", 783.0, "Argoul 2026 Table 2"),
    ("mouse", "C57BL/6", "female", "PFHxA", 6830.0, "Argoul 2026 Table 2"),
    ("mouse", "C57BL/6", "female", "PFOA", 4.5, "Argoul 2026 Table 2"),
    ("mouse", "C57BL/6", "female", "PFNA", 2.1, "Argoul 2026 Table 2"),
    ("mouse", "C57BL/6", "female", "PFDA", 1.8, "Argoul 2026 Table 2"),
    ("mouse", "C57BL/6", "female", "PFHxS", 1.3, "Argoul 2026 Table 2"),
    ("mouse", "C57BL/6", "female", "PFOS", 2.7, "Argoul 2026 Table 2"),
]

# values a compilation attributes to a paper, superseded by the paper itself
SUPERSEDED = [
    ("rat", "male", "PFOA", 135.0, "ATSDR 2021 Table 3-6 citing Ohmori 2003",
     50.5, "Ohmori 2003 Table 1"),
    ("rat", "male", "PFDA", 207.0, "ATSDR 2021 Table 3-6 citing Ohmori 2003",
     5.2, "Ohmori 2003 Table 1"),
    ("rat", "female", "PFDA", 140.0, "ATSDR 2021 Table 3-6 citing Ohmori 2003",
     5.3, "Ohmori 2003 Table 1"),
]

# ---- measured renal fraction of elimination -----------------------------
F_REN = [
    ("rat", "male", "PFHpA", 0.92, "Kudo 2001, % of dose in urine, Wistar"),
    ("rat", "male", "PFOA", 0.55, "Kudo 2001, % of dose in urine, Wistar"),
    ("rat", "male", "PFNA", 0.020, "Kudo 2001, % of dose in urine, Wistar"),
    ("rat", "male", "PFDA", 0.002, "Kudo 2001, % of dose in urine, Wistar"),
    ("rat", "female", "PFOA", 0.80, "Kudo 2001, % of dose in urine, Wistar"),
    ("rat", "female", "PFNA", 0.51, "Kudo 2001, % of dose in urine, Wistar"),
    ("mouse", "female", "PFBA", 1.00, "Argoul 2026 Table 2, C57BL/6"),
    ("mouse", "female", "PFHxA", 0.993, "Argoul 2026 Table 2, C57BL/6"),
    ("mouse", "female", "PFOA", 0.791, "Argoul 2026 Table 2, C57BL/6"),
    ("mouse", "female", "PFNA", 0.286, "Argoul 2026 Table 2, C57BL/6"),
    ("mouse", "female", "PFDA", 0.073, "Argoul 2026 Table 2, C57BL/6"),
    ("mouse", "female", "PFHxS", 0.44, "Argoul 2026 Table 2, C57BL/6"),
    ("mouse", "female", "PFOS", 0.463, "Argoul 2026 Table 2, C57BL/6"),
]


def assert_same_species(label, *tags):
    u = {t for t in tags if t}
    if len(u) != 1:
        raise SystemExit(f"SPECIES CONTAMINATION in {label}: {sorted(u)}")
    return u.pop()


def load_fu():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] != "Ryu 2024" or not r["fu"]:
            continue
        out[(r["species"], ALIAS.get(r["chemical"], r["chemical"]))] = \
            float(r["fu"])
    return out


def pick_fu(fu, sp, chem, strain):
    cands = [m for m in FU_MATRIX[sp] if (m, chem) in fu]
    if not cands:
        return None, None
    exact = [m for m in cands if STRAIN_OF[m] == strain]
    m = (exact or cands)[0]
    assert_same_species(f"{sp}/{chem}", sp, m.split(" (")[0])
    return fu[(m, chem)], m


def build(fu, gfr_map):
    ev = []
    for tier, table in (("DIRECT", RENAL), ("UPPER", TOTAL)):
        for sp, strain, sex, chem, val, src in table:
            chem = ALIAS.get(chem, chem)
            f, mat = pick_fu(fu, sp, chem, strain)
            if f is None:
                continue
            g = gfr_map[sp]
            filt = f * g
            ev.append(dict(species=sp, strain=strain, sex=sex, chemical=chem,
                           tier=tier, fu=f, fu_matrix=mat, gfr=g, filt=filt,
                           value=val, src=src, T=val / filt))
    return ev


def consolidate(ev, fu, gfr_map):
    # group the evidence
    byk = {}
    for e in ev:
        byk.setdefault((e["species"], e["sex"], e["chemical"]), []).append(e)
    # merge 'both' into each sex so a population estimate informs both
    out = {}
    for (sp, sex, chem), es in byk.items():
        out[(sp, sex, chem)] = list(es)
    for (sp, sex, chem), es in list(byk.items()):
        if sex != "both":
            continue
        for other in ("male", "female"):
            out.setdefault((sp, other, chem), []).extend(es)

    fren = {(s, x, c): (v, src) for s, x, c, v, src in F_REN}
    res = []
    for (sp, sex, chem), es in sorted(out.items()):
        if sex == "both" and (sp, "male", chem) in out:
            continue
        direct = [e for e in es if e["tier"] == "DIRECT"]
        upper = [e for e in es if e["tier"] == "UPPER"]
        f, mat = pick_fu(fu, sp, chem, "")
        filt = f * gfr_map[sp]
        row = dict(species=sp, sex=sex, chemical=chem, n_carbons=NC[chem],
                   family=FAM[chem], fu=f, fu_matrix=mat, gfr=gfr_map[sp],
                   filt=filt)
        if direct:
            Ts = sorted(e["T"] for e in direct)
            row.update(tier="DIRECT", T_lo=Ts[0], T_hi=Ts[-1],
                       n_sources=len(direct),
                       sources="; ".join(sorted({e["src"] for e in direct})))
            if Ts[-1] < 1:
                row.update(verdict="net_reabsorption",
                           min_pct=100 * (1 - Ts[-1]))
            elif Ts[0] > 1:
                row.update(verdict="net_secretion", min_pct=None)
            else:
                row.update(verdict="DISPUTED", min_pct=None)
            res.append(row)
            continue
        # no direct: try the lower bound from a renal fraction
        hit = fren.get((sp, sex, chem))
        if upper:
            loosest = max(upper, key=lambda e: e["value"])
            tight = min(upper, key=lambda e: e["value"])
            thresh = filt / loosest["value"]
            row.update(tier="UPPER", T_lo=None, T_hi=loosest["T"],
                       n_sources=len(upper),
                       sources="; ".join(sorted({e["src"] for e in upper})),
                       cltot_loosest=loosest["value"],
                       cltot_tightest=tight["value"], threshold=thresh)
            if loosest["T"] < 1:
                row.update(verdict="net_reabsorption_holds",
                           min_pct=100 * (1 - loosest["T"]))
            elif hit and hit[0] > thresh:
                row.update(tier="LOWER", verdict="net_secretion_established",
                           min_pct=None, f_ren=hit[0], f_ren_src=hit[1],
                           T_lo=hit[0] * loosest["value"] / filt)
            else:
                row.update(verdict="undetermined", min_pct=None,
                           f_ren=hit[0] if hit else None,
                           f_ren_src=hit[1] if hit else None)
            res.append(row)
    return res


def main():
    fu = load_fu()
    ev = build(fu, GFR)
    cons = consolidate(ev, fu, GFR)

    p1 = HERE / "db" / "tubular_factor_evidence.csv"
    with p1.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "strain", "sex", "chemical", "n_carbons",
                    "family", "tier", "fu", "fu_matrix", "gfr_mL_kg_day",
                    "filtered_load_mL_kg_day", "clearance_mL_kg_day",
                    "clearance_kind", "T", "clearance_source"])
        for e in sorted(ev, key=lambda x: (x["species"], x["chemical"],
                                           x["sex"], x["tier"])):
            w.writerow([e["species"], e["strain"], e["sex"], e["chemical"],
                        NC[e["chemical"]], FAM[e["chemical"]], e["tier"],
                        e["fu"], e["fu_matrix"], e["gfr"],
                        round(e["filt"], 3), e["value"],
                        "renal" if e["tier"] == "DIRECT" else "total",
                        round(e["T"], 5), e["src"]])
        for sp, sex, chem, bad, badsrc, good, goodsrc in SUPERSEDED:
            w.writerow([sp, "", sex, chem, NC[chem], FAM[chem], "SUPERSEDED",
                        "", "", "", "", bad, "total",
                        "", f"{badsrc} -- SUPERSEDED by {goodsrc} = {good}, "
                        f"a {bad/good:.1f}x discrepancy"])

    p2 = HERE / "db" / "tubular_factor_consolidated.csv"
    with p2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "sex", "chemical", "n_carbons", "family",
                    "fu", "fu_matrix", "gfr_mL_kg_day",
                    "filtered_load_mL_kg_day", "tier", "T_low", "T_high",
                    "verdict", "min_pct_reabsorbed",
                    "secretion_threshold_renal_fraction",
                    "measured_renal_fraction", "n_clearance_sources",
                    "clearance_sources"])
        for r in cons:
            w.writerow([r["species"], r["sex"], r["chemical"], r["n_carbons"],
                        r["family"], r["fu"], r["fu_matrix"], r["gfr"],
                        round(r["filt"], 3), r["tier"],
                        "" if r.get("T_lo") is None else round(r["T_lo"], 5),
                        "" if r.get("T_hi") is None else round(r["T_hi"], 5),
                        r["verdict"],
                        "" if r.get("min_pct") is None
                        else round(r["min_pct"], 1),
                        "" if r.get("threshold") is None
                        else round(r["threshold"], 4),
                        "" if r.get("f_ren") is None else r["f_ren"],
                        r["n_sources"], r["sources"]])

    report(ev, cons, fu)
    print(f"\nwrote {p1.relative_to(HERE)}  ({len(ev)+len(SUPERSEDED)} rows)")
    print(f"wrote {p2.relative_to(HERE)}  ({len(cons)} rows)")


def report(ev, cons, fu):
    bar = "=" * 78
    print(bar)
    print("CONSOLIDATED VERDICTS")
    print(bar)
    sym = {"net_reabsorption": "reabs", "net_reabsorption_holds": "reabs>=",
           "net_secretion": "SECR", "net_secretion_established": "SECR*",
           "DISPUTED": "DISPUTED", "undetermined": "-"}
    print(f"  {'species':7}{'sex':7}{'chem':9}{'nC':>3}{'T':>16}"
          f"{'tier':8}  verdict")
    for r in sorted(cons, key=lambda x: (x["species"], x["sex"],
                                         x["n_carbons"], x["chemical"])):
        lo, hi = r.get("T_lo"), r.get("T_hi")
        if r["tier"] == "DIRECT":
            t = f"{lo:.3f}" if abs(hi - lo) < 1e-9 else f"{lo:.2f}-{hi:.2f}"
        elif r["tier"] == "LOWER":
            t = f">={lo:.1f}"
        else:
            t = f"<={hi:.3f}"
        extra = (f" {r['min_pct']:.0f}%" if r.get("min_pct") is not None
                 else "")
        print(f"  {r['species']:7}{r['sex']:7}{r['chemical']:9}"
              f"{r['n_carbons']:>3}{t:>16}{r['tier']:8}  "
              f"{sym[r['verdict']]}{extra}")

    print()
    print(bar)
    print("WHAT CHANGED FROM THE STAGE 4E TABLE")
    print(bar)
    prev = {}
    f = HERE / "db" / "tubular_factor_all_species.csv"
    for x in csv.DictReader(f.open()):
        prev[(x["species"], x["chemical"])] = x["verdict"]
    moved = 0
    for r in sorted(cons, key=lambda x: (x["species"], x["chemical"])):
        if r["sex"] == "female":
            continue
        old = prev.get((r["species"], r["chemical"]))
        if old is None:
            print(f"  NEW   {r['species']:6}{r['chemical']:8}"
                  f"-> {r['verdict']}")
            moved += 1
            continue
        oldv = ("net_reabsorption" if "reabsorption" in old
                else "net_secretion" if "SECRETION" in old
                else "undetermined")
        newv = ("net_reabsorption" if "reabsorption" in r["verdict"]
                else "net_secretion" if "secretion" in r["verdict"]
                else r["verdict"])
        if oldv != newv:
            print(f"  MOVED {r['species']:6}{r['chemical']:8}"
                  f"{oldv:18} -> {newv}")
            moved += 1
    print(f"\n  {moved} cells changed. The PFDA and PFHpA moves are the")
    print("  Ohmori correction; the rest are the sex split and the new tier.")

    print()
    print(bar)
    print("THREE SUPERSEDED VALUES, KEPT IN THE EVIDENCE FILE")
    print(bar)
    for sp, sex, chem, bad, badsrc, good, goodsrc in SUPERSEDED:
        print(f"  {sp} {sex} {chem}: {badsrc.split(' citing')[0]} says {bad}, "
              f"the paper says {good}  ({bad/good:.1f}x)")

    print()
    print(bar)
    print("A VALIDATION THAT FALLS OUT: ARGOUL IS INTERNALLY CONSISTENT")
    print(bar)
    print("  Argoul 2026 reports mouse renal clearance (Table 4), total")
    print("  clearance (Table 2) and the urinary % of dose (Table 2)")
    print("  separately. CL_R / CL_tot must equal the urinary fraction.")
    print()
    print(f"  {'chem':8}{'CL_R':>9}{'CL_tot':>9}{'ratio':>8}"
          f"{'stated':>9}  agree")
    fren = {(s, x, c): v for s, x, c, v, _ in F_REN}
    for chem in ("PFBA", "PFHxA", "PFOA", "PFNA", "PFDA", "PFHxS", "PFOS"):
        r = [e for e in ev if e["species"] == "mouse" and e["tier"] == "DIRECT"
             and e["chemical"] == chem]
        t = [e for e in ev if e["species"] == "mouse" and e["tier"] == "UPPER"
             and e["chemical"] == chem and "Argoul" in e["src"]]
        if not r or not t:
            continue
        ratio = r[0]["value"] / t[0]["value"]
        stated = fren[("mouse", "female", chem)]
        ok = "yes" if abs(ratio - stated) < 0.02 else "NO"
        print(f"  {chem:8}{r[0]['value']:>9.2f}{t[0]['value']:>9.1f}"
              f"{ratio:>8.3f}{stated:>9.3f}  {ok}")
    print()
    print("  All seven agree to within two percentage points, across four")
    print("  orders of magnitude of clearance. Argoul's three tables are")
    print("  mutually consistent, which is the strongest internal check any")
    print("  source in this strand has passed.")


if __name__ == "__main__":
    main()
