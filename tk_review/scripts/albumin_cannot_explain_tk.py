#!/usr/bin/env python3
"""Can albumin binding explain the sex and species differences in PFAS TK?

No, and this script is the arithmetic.  Three tests, each designed so that
albumin binding is held fixed or measured, and each with enough of a gap
that no plausible model closes it.

  1  SAME PROTEIN.  Male and female rats of one strain carry the same
     serum albumin.  Binding is identical BY CONSTRUCTION -- not measured
     to be equal, incapable of differing.  Any half-life difference
     between them is therefore something other than albumin.  This test
     needs no binding measurement and no model.

  2  MATCHED BINDING, DIFFERENT HEAD GROUP.  PFNA and PFOS both have
     eight perfluorinated carbons and sit within noise of each other on
     albumin.  Albumin affinity is dominated by chain length, so if
     albumin set the sex difference these two should behave alike.

  3  THE SPECIES BUDGET.  How large is the measured rat-human albumin
     difference, against the observed half-life difference?  The binding
     term is taken from species_albumin_comparison.py, which uses only
     WITHIN-LABORATORY pairs, because this project has measured the
     laboratory term at 45% of the albumin variance and the species term
     at 2% (what_drives_albumin_spread.py).

WHY THE GAP CANNOT BE MODELLED AWAY.  Binding reaches the half-life
through distribution: t1/2 = ln2 * Vd / CL, and albumin acts on Vd.  Grant
the most generous coupling -- t1/2 strictly proportional to binding
strength -- and a 1.3x binding difference buys 1.3x of half-life.  The
observations are 47x and 7000x.  Argoul 2026 makes the same point from
inside one experiment: across 11 PFAS dosed as one cocktail to one sex of
one strain, clearance spans 5254x while Vss spans 7.7x.  Vd is the small
term, and binding only moves Vd.

WHAT DOES EXPLAIN IT is in db/transporter_kinetics.csv: Oatp1a1, apical
and male-predominant in rat kidney, driving REABSORPTION rather than
secretion, and not transporting the short chains.  Part 4 prints the rows.

Run:  python3 -I tk_review/scripts/albumin_cannot_explain_tk.py
"""
import collections
import csv
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db")
ALB = os.path.join(HERE, "..", "notebook", "albumin_logK_Kd_revised.csv")

HOURS = {"h": 1.0, "hr": 1.0, "hours": 1.0,
         "d": 24.0, "day": 24.0, "days": 24.0,
         "y": 8766.0, "year": 8766.0, "years": 8766.0}

# Perfluorinated carbons and head group, for the matched-binding test.
CHEM = {"PFBA": (3, "carboxylate"), "PFPeA": (4, "carboxylate"),
        "PFHxA": (5, "carboxylate"), "PFHpA": (6, "carboxylate"),
        "PFOA": (7, "carboxylate"), "PFNA": (8, "carboxylate"),
        "PFDA": (9, "carboxylate"), "PFBS": (4, "sulfonate"),
        "PFHxS": (6, "sulfonate"), "PFOS": (8, "sulfonate")}

# species_albumin_comparison.py, part 1: Starnes 2024, four albumins on
# one instrument, 7 compounds.  Median dlogK against human.
RAT_VS_HUMAN_LOGK = 0.109          # rat binds this much tighter
BOVINE_VS_HUMAN_LOGK = -0.098      # for scale: bovine is the same size


def hours(value, units):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v * HOURS[units.lower().strip()] \
        if (units or "").lower().strip() in HOURS else None


def animal():
    with open(os.path.join(DB, "animal_halflife_measured.csv"),
              newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def human_pfoa():
    out = []
    with open(os.path.join(DB, "human_halflife_extended.csv"),
              newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["chemical"] != "PFOA":
                continue
            try:
                out.append(float(r["halflife_y"]) * 8766.0)
            except (TypeError, ValueError):
                pass
    return out


def albumin_logk():
    by = collections.defaultdict(list)
    with open(ALB, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                by[r["PFAS"]].append(
                    float(r["log K (L water / kg protein)"]))
            except (TypeError, ValueError):
                pass
    return by


def rule(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main():
    rows = animal()
    alb = albumin_logk()

    rule("1.  SAME PROTEIN: male vs female rat")
    print("  Albumin is identical between the sexes, so this difference")
    print("  cannot be a binding difference. No model is needed.\n")
    print("  %-7s %3s %-12s %10s %10s %9s  %s"
          % ("PFAS", "nC", "head group", "M med h", "F med h", "M/F",
             "n M,F"))
    sex = []
    for c in sorted(CHEM, key=lambda c: (CHEM[c][1], CHEM[c][0])):
        m = [h for h in (hours(r["halflife"], r["halflife_units"])
                         for r in rows
                         if r["chemical"] == c and r["species"] == "rat"
                         and r["sex"] == "M") if h]
        f = [h for h in (hours(r["halflife"], r["halflife_units"])
                         for r in rows
                         if r["chemical"] == c and r["species"] == "rat"
                         and r["sex"] == "F") if h]
        if not m or not f:
            continue
        mm, ff = st.median(m), st.median(f)
        sex.append((c, CHEM[c][0], CHEM[c][1], mm / ff, len(m), len(f)))
        print("  %-7s %3d %-12s %10.1f %10.1f %8.1fx  %d,%d"
              % (c, CHEM[c][0], CHEM[c][1], mm, ff, mm / ff, len(m),
                 len(f)))
    worst = max(sex, key=lambda s: s[3])
    print(f"\n  >>> Up to {worst[3]:.0f}x ({worst[0]}) on one protein.")
    print("      Whatever sets this, it is not albumin affinity.")
    print("      Thin cells: PFHpA and PFDA rest on 1-2 animals a sex;")
    print("      PFOA (15,12), PFHxS (15,14) and PFOS (11,9) do not.")

    rule("2.  MATCHED BINDING: PFNA vs PFOS, both 8 fluorinated carbons")
    a, b = st.median(alb["PFNA"]), st.median(alb["PFOS"])
    ra = dict((c, r) for c, _, _, r, _, _ in sex)
    print("  Albumin affinity tracks chain length. These two share a chain")
    print("  length and sit within noise of each other on albumin:\n")
    print("    %-6s  albumin log K %5.2f  (n=%d)" % ("PFNA", a,
                                                     len(alb["PFNA"])))
    print("    %-6s  albumin log K %5.2f  (n=%d)" % ("PFOS", b,
                                                     len(alb["PFOS"])))
    print("    apart by %.2f log = %.2fx -- inside the ~1 log spread of"
          % (abs(a - b), 10 ** abs(a - b)))
    print("    either compound, so they are not resolvably different.\n")
    print("    rat M/F half-life ratio   PFNA %.1fx   PFOS %.1fx"
          % (ra["PFNA"], ra["PFOS"]))
    print("    apart by %.0fx" % (ra["PFNA"] / ra["PFOS"]))
    print("\n  >>> Indistinguishable on albumin, %.0fx apart in vivo."
          % (ra["PFNA"] / ra["PFOS"]))
    print("      The difference between them is the HEAD GROUP, which is")
    print("      what a transporter sees and albumin largely does not.")

    rule("3.  THE SPECIES BUDGET: rat vs human PFOA")
    fem = [h for h in (hours(r["halflife"], r["halflife_units"])
                       for r in rows if r["chemical"] == "PFOA"
                       and r["species"] == "rat" and r["sex"] == "F") if h]
    mal = [h for h in (hours(r["halflife"], r["halflife_units"])
                       for r in rows if r["chemical"] == "PFOA"
                       and r["species"] == "rat" and r["sex"] == "M") if h]
    hu = human_pfoa()
    f, m, h = st.median(fem), st.median(mal), st.median(hu)
    print("  observed half-life")
    print("    female rat  %10.1f h   (n=%d)" % (f, len(fem)))
    print("    male rat    %10.1f h   (n=%d)" % (m, len(mal)))
    print("    human       %10.0f h   = %.2f y  (n=%d)"
          % (h, h / 8766, len(hu)))
    print("    human / female rat  %8.0fx" % (h / f))
    print("    human / male rat    %8.0fx" % (h / m))
    print("\n  measured albumin difference, rat vs human")
    print("    %+.3f log = %.2fx  (Starnes 2024, one instrument, 7 "
          "compounds)" % (RAT_VS_HUMAN_LOGK, 10 ** RAT_VS_HUMAN_LOGK))
    print("    for scale, bovine vs human is %.2fx -- the same size, and"
          % (10 ** abs(BOVINE_VS_HUMAN_LOGK)))
    print("    nobody proposes bovine albumin to explain a TK difference.")
    print("\n  the budget, granting binding the most generous coupling")
    print("    needed   %8.0fx      observed human / female rat" % (h / f))
    print("    supplied %8.2fx      if t1/2 were proportional to binding"
          % (10 ** RAT_VS_HUMAN_LOGK))
    print("    short by %8.0fx      = %.3f%% of the gap explained"
          % ((h / f) / 10 ** RAT_VS_HUMAN_LOGK,
             100 * 10 ** RAT_VS_HUMAN_LOGK / (h / f)))

    rule("4.  WHAT IS LEFT: the renal transporter rows")
    seen, n = set(), 0
    with open(os.path.join(DB, "transporter_kinetics.csv"),
              newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            e = (r["sex_or_hormone_effect"] or "").strip()
            key = (r["transporter"], e[:60])
            if not e or key in seen or "not assessed" in e.lower():
                continue
            if not any(w in e.lower() for w in
                       ("predominant", "sex-different", "sex difference",
                        "cannot explain")):
                continue
            seen.add(key)
            n += 1
            print("  %-26s %-12s %s"
                  % (r["transporter"][:26],
                     r["species_orthologue"][:12],
                     r["direction_reabsorption_or_secretion"][:28]))
            for i in range(0, min(len(e), 300), 66):
                print("     %s" % e[i:i + 66])
    print(f"\n  {n} distinct sex/hormone findings.  The pattern they make:")
    print("    Oatp1a1  male-predominant, APICAL, drives REABSORPTION,")
    print("             does not transport the short chains -- which is")
    print("             exactly where part 1 shows no sex difference.")
    print("    Oat1/Oat3 are basolateral SECRETION, and where Oat1 is")
    print("             male-predominant that would speed males up -- the")
    print("             WRONG direction for a longer male half-life, which")
    print("             is why the field settled on Oatp1a1 instead.")
    print("    Oat2 IS strongly female-predominant, the right direction,")
    print("             but later direct transport work shows it does not")
    print("             carry PFOA: excluded on substrate, not on sex.")
    print("    The mouse replicates it: renal Oatp1a1 is male-predominant")
    print("             there too, so the argument is not rat-specific.")
    print("\n  A reabsorption transporter explains all three tests at once:")
    print("    it differs by sex within a species (test 1), it reads the")
    print("    head group as well as the chain (test 2), and its renal")
    print("    expression differs between species (test 3).  Albumin")
    print("    explains none of them.")

    rule("WHAT THIS DOES NOT SHOW")
    print("  Albumin still sets the free fraction, and a mass-balance")
    print("  model needs it.  The claim here is narrow and specific:")
    print("  albumin affinity is not the knob that sets the DIFFERENCES")
    print("  in half-life between sexes or between species.")
    print("  The rat-human binding term is 7 compounds from one")
    print("  held-out DSF study plus a Han 2003 pair that is really one")
    print("  shared number; species_albumin_comparison.py calls that")
    print("  equivalence provisional, and it stays provisional here.")
    print("  Nothing above rests on its exact value: at 1.3x, or 3x, or")
    print("  10x, the gap in part 3 is still three orders wide.")


if __name__ == "__main__":
    main()
