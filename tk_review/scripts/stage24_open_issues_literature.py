#!/usr/bin/env python3
"""What the literature can and cannot settle for the six open issues (Stage 24).

    python3 -I scripts/stage24_open_issues_literature.py

Writes db/open_issues_literature.csv and db/starnes_vs_ryu_rat.csv. Idempotent.

TWO OF THE SIX MOVE; FOUR DO NOT. The searches behind this stage are recorded
in OPEN_ISSUES below, one row per issue, so a later reader can see what was
looked for as well as what was found.

ISSUE 1, THE RAT UNBOUND FRACTION. The claim that no independent rat binding
data exist was too strong. The project already holds rat albumin binding from
six studies, and ONE of them, Starnes 2024, measured rat AND human albumin
for the same eight compounds in one laboratory by one assay. That is a
genuine independent cross-species binding comparison.

IT DISAGREES WITH RYU, AND THE DISAGREEMENT IS THE FINDING. Starnes gives a
geometric mean K_d ratio (rat/human) of 0.68 over the five compounds
both sources cover -- rat albumin binds PFAS MORE
tightly, which should leave LESS free -- while Ryu's plasma unbound fractions
give a ratio of 1.32, rat plasma leaving MORE free. They point in opposite
directions, and agree on only one of the five compounds both cover.

THAT DOES NOT OVERTURN RYU, AND THE REASON IS QUANTITATIVE. The two measure
different things: K_d to isolated albumin by differential scanning
fluorimetry against f_u in whole plasma by equilibrium dialysis. Starnes's
absolute K_d values are in the millimolar range, and a one-site model
f_u = 1/(1 + [albumin]/K_d) at a physiological 0.6 mM albumin turns them into
unbound fractions near 0.4, where the measured values are near 0.003 -- two
orders of magnitude out. A method whose absolute scale cannot reproduce
plasma behaviour cannot adjudicate a plasma measurement. What it can do is
remove the claim that the direction is unchallenged, which is how it is now
reported.

A UNIT TRAP WAS AVOIDED HERE. db/protein_binding.csv carries Starnes rows
from TWO papers: the cross-species study in mM and an earlier human-only
study in uM. Comparing across them silently inflates the ratio a
thousand-fold. Only rows whose kd_units is mM are used below.

ISSUE 4, PFOS SECRETION OR BILIARY. Harada 2007 measured biliary clearance
in humans directly (serum and bile from patients with biliary drainage) and
concluded that biliary excretion is "comparable in both rats and humans".
That bounds the size of the alternative: human biliary clearance for PFOS is
2.98 mL/kg/day gross, and after 97% biliary reabsorption about 0.089 net,
against a renal clearance of 0.012 to 0.037. The non-renal route is
therefore the same size as or larger than the renal one for PFOS, so a fall
in TOTAL clearance with dose can plausibly be biliary. What the literature
does NOT supply is a biliary or renal dose series in the rat, which is the
only thing that would separate the two. The alternative is bounded, not
eliminated.

THE OTHER FOUR. Searches found nothing that moves them, and the reasons
differ: for two of them no study of the right design appears to exist, and
for two the limitation is structural rather than evidential.
"""
import csv
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent

ALBUMIN_mM = 0.60     # physiological serum albumin, for the scale check

OPEN_ISSUES = [
    dict(issue="1. rat unbound fraction rests on one source",
         searched="PubMed for rat plasma protein binding / unbound fraction "
                  "of PFAS from any laboratory; the project's own "
                  "protein_binding.csv for rat albumin",
         found="PARTLY MOVES. No second rat PLASMA f_u exists. But rat "
               "albumin binding is held from six studies, and Starnes 2024 "
               "measured rat and human albumin for the same 8 compounds in "
               "one assay -- an independent cross-species comparison that "
               "the paper did not use.",
         effect="The direction DISAGREES with Ryu (Kd ratio 0.68 against f_u "
                "ratio 1.32, agreeing on 1 of 5 compounds). Starnes cannot "
                "adjudicate, because its absolute Kd implies f_u near 0.4 "
                "where the measured value is near 0.003. Report as an "
                "unresolved cross-method disagreement, not as support.",
         status="partly moves"),
    dict(issue="2. human kinetics from one volunteer",
         searched="PubMed for controlled human PFAS kinetic studies beyond "
                  "Abraham 2024",
         found="NO MOVE. The project already carries the population "
               "alternatives it can (Zhang 2013 n=66, Fujii 2015 n=10, "
               "Harada 2005 n=20). No further controlled-dosing human study "
               "was located; deliberate human dosing is not repeatable.",
         effect="none; the existing population sensitivity analysis is the "
                "available mitigation",
         status="no move"),
    dict(issue="3. no within-compound uncertainty propagation",
         searched="not a literature question",
         found="STRUCTURAL, NOT EVIDENTIAL. Kim 2016 reports SEMs and is now "
               "ingested with them (db/primary_2026/kim2016_table4_rat.csv), "
               "so three compounds could carry propagated error. The rest "
               "would need SEMs gathered from each primary source.",
         effect="doable with work already in hand for three compounds; no "
                "new literature required",
         status="no move (structural)"),
    dict(issue="4. PFOS: tubular secretion or biliary elimination",
         searched="PubMed for PFOS biliary excretion, enterohepatic "
                  "circulation and rat bile studies",
         found="PARTLY MOVES. Harada 2007 measured human biliary clearance "
               "directly and states biliary excretion is comparable in rats "
               "and humans. The project holds the extract but has never used "
               "it.",
         effect="Bounds the alternative: human PFOS biliary clearance 2.98 "
                "gross, about 0.089 net of 97% reabsorption, against a renal "
                "clearance of 0.012-0.037 mL/kg/day. The non-renal route is "
                "the same size or larger, so a falling TOTAL clearance can "
                "plausibly be biliary. Still no rat dose series in bile or "
                "urine, so the two are bounded, not separated.",
         status="partly moves"),
    dict(issue="5. PFBA and PFBS lack a human renal clearance",
         searched="PubMed for human urinary or renal clearance of PFBA and "
                  "PFBS",
         found="NO MOVE. Nothing located. Olsen 2009 observed urine as an "
               "elimination pathway in its human subjects but published no "
               "clearance, and Chang 2008's human arm reports serum "
               "half-lives only.",
         effect="none; these two remain blocked on the human side, as Stage "
                "18 concluded",
         status="no move"),
    dict(issue="6. mouse column mixes strains",
         searched="PubMed for single-strain mouse PFAS kinetic series across "
                  "several compounds",
         found="NO MOVE WITHIN THE PRESENT SET. Fujii 2015 is one strain "
               "(FVB/NJcl) across 7 carboxylates and already supplies the "
               "mouse dose series; the mixing arises where Chang 2012 (CD-1) "
               "and Sundstrom 2012 are brought in for the sulfonates.",
         effect="could be narrowed by restricting the mouse column to "
                "Fujii's strain at the cost of losing the sulfonates; a "
                "presentational choice, not new evidence",
         status="no move"),
]

# Harada 2007, from the project's extract; human, n = 4 serum/bile pairs.
HARADA2007 = {
    "PFOA": dict(biliary_cl=1.06, biliary_reabsorption=0.89,
                 renal_cl_low=0.030, renal_cl_high=0.044),
    "PFOS": dict(biliary_cl=2.98, biliary_reabsorption=0.97,
                 renal_cl_low=0.012, renal_cl_high=0.037),
}


def starnes_vs_ryu():
    """Matched rat/human albumin Kd against matched rat/human plasma f_u."""
    d = {}
    for r in csv.DictReader((HERE / "db" / "protein_binding.csv").open()):
        # ONLY the cross-species paper's rows: its Kd is in mM. The earlier
        # human-only Starnes study reports uM, and mixing them is a
        # thousand-fold error.
        if (r["study"] == "Starnes et al."
                and "ALBUMIN" in r["protein"].upper()
                and r["kd_units"] == "mM" and r["kd"]
                and r["kd"] != "nonbinding"):
            d[(r["chemical"], r["protein_species"])] = float(r["kd"])
    fu = {}
    for r in csv.DictReader(
            (HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r["experiment"] == "Ryu 2024" and r["fu"]:
            fu[(r["chemical"], r["species"].split(" (")[0])] = float(r["fu"])
    out = []
    for chem in sorted({c for c, _ in d}):
        kr, kh = d.get((chem, "rat")), d.get((chem, "human"))
        fr, fh = fu.get((chem, "rat")), fu.get((chem, "human"))
        if not (kr and kh):
            continue
        row = dict(chemical=chem, kd_rat_mM=kr, kd_human_mM=kh,
                   kd_ratio_rat_over_human=round(kr / kh, 3))
        if fr and fh:
            row.update(fu_rat=fr, fu_human=fh,
                       fu_ratio_rat_over_human=round(fr / fh, 3),
                       directions_agree=("yes" if (kr / kh > 1) == (fr / fh > 1)
                                         else "no"))
            # one-site scale check at physiological albumin
            row["fu_implied_by_kd_rat"] = round(1 / (1 + ALBUMIN_mM / kr), 4)
        else:
            row.update(fu_rat="", fu_human="", fu_ratio_rat_over_human="",
                       directions_agree="", fu_implied_by_kd_rat="")
        out.append(row)
    return out


def gm(v):
    return 10 ** (sum(math.log10(x) for x in v) / len(v))


def main():
    bar = "=" * 78
    rows = starnes_vs_ryu()
    paired = [r for r in rows if r["directions_agree"]]

    print(bar)
    print("ISSUE 1: AN INDEPENDENT RAT COMPARISON EXISTS, AND IT DISAGREES")
    print(bar)
    print(f"  {'cpd':10}{'Kd rat':>9}{'Kd human':>10}{'Kd r/h':>9}"
          f"{'fu r/h':>9}{'agree':>8}{'fu from Kd':>12}")
    for r in rows:
        fr = (f"{r['fu_ratio_rat_over_human']:.2f}"
              if r["fu_ratio_rat_over_human"] != "" else "--")
        ag = r["directions_agree"] or "--"
        fk = (f"{r['fu_implied_by_kd_rat']:.3f}"
              if r["fu_implied_by_kd_rat"] != "" else "--")
        print(f"  {r['chemical'][:9]:10}{r['kd_rat_mM']:>9.2f}"
              f"{r['kd_human_mM']:>10.2f}{r['kd_ratio_rat_over_human']:>9.2f}"
              f"{fr:>9}{ag:>8}{fk:>12}")
    print()
    n_ag = sum(1 for r in paired if r["directions_agree"] == "yes")
    print(f"  directional agreement: {n_ag} of {len(paired)}")
    # both means over the SAME five compounds, or they are not comparable
    print(f"  geometric mean Kd  rat/human = "
          f"{gm([r['kd_ratio_rat_over_human'] for r in paired]):.2f}"
          f"   (<1: rat albumin binds TIGHTER, so less should be free)")
    print(f"  geometric mean f_u rat/human = "
          f"{gm([r['fu_ratio_rat_over_human'] for r in paired]):.2f}"
          f"   (>1: rat plasma leaves MORE free)")
    print()
    print("  THE TWO POINT OPPOSITE WAYS. But the scale check in the last")
    print("  column shows why neither overturns the other: a one-site model")
    print(f"  at {ALBUMIN_mM} mM albumin turns Starnes's Kd into unbound")
    print("  fractions near 0.4, where the measured plasma values are near")
    print("  0.003 -- two orders of magnitude out. Differential scanning")
    print("  fluorimetry ranks affinities; it does not reproduce plasma.")
    print("  The honest report is an unresolved cross-method disagreement,")
    print("  which is still more than the paper says now.")
    print()

    print(bar)
    print("ISSUE 4: THE BILIARY ALTERNATIVE IS NOW BOUNDED")
    print(bar)
    print(f"  {'cpd':7}{'biliary CL':>12}{'reabsorbed':>12}{'net biliary':>13}"
          f"{'renal CL':>16}   non-renal is")
    for chem, h in sorted(HARADA2007.items()):
        net = h["biliary_cl"] * (1 - h["biliary_reabsorption"])
        lo, hi = h["renal_cl_low"], h["renal_cl_high"]
        print(f"  {chem:7}{h['biliary_cl']:>12.2f}"
              f"{h['biliary_reabsorption']*100:>11.0f}%{net:>13.3f}"
              f"{f'{lo:.3f}-{hi:.3f}':>16}   "
              f"{net/hi:.1f}-{net/lo:.1f}x the renal route")
    print()
    print("  Harada 2007, human, n = 4 serum/bile pairs, and it states that")
    print("  biliary excretion is comparable in rats and humans. For PFOS the")
    print("  net biliary route is 2 to 7 times the renal one, so a fall in")
    print("  TOTAL clearance with dose can plausibly be saturating biliary")
    print("  transport rather than saturating tubular secretion.")
    print()
    print("  WHAT IS STILL MISSING is a dose series in bile or in urine")
    print("  specifically, in the rat. Without it the two remain bounded")
    print("  rather than separated, exactly as Stage 21 said -- but the size")
    print("  of the alternative is now a number rather than a worry.")
    print()

    print(bar)
    print("ALL SIX ISSUES")
    print(bar)
    for d in OPEN_ISSUES:
        print(f"  [{d['status'].upper()}] {d['issue']}")
        print(f"     found : {d['found'][:200]}")
        print()

    p = HERE / "db" / "starnes_vs_ryu_rat.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["chemical"]):
            w.writerow(r)
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} rows)")

    p = HERE / "db" / "open_issues_literature.csv"
    cols = ["issue", "searched", "found", "effect", "status"]
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for d in sorted(OPEN_ISSUES, key=lambda x: x["issue"]):
            w.writerow({k: d[k] for k in cols})
    print(f"wrote {p.relative_to(HERE)}  ({len(OPEN_ISSUES)} rows)")


if __name__ == "__main__":
    main()
