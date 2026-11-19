#!/usr/bin/env python3
"""Screen a transporter source for whether it can EVER reach CL_int.

Run this BEFORE reading a paper, and after ingesting one. It is the
transporter counterpart of screen_fits_database.py, and it works the same
way: there is one identity the strand is built on, and a source earns a place
in the kinetics table only if it supplies what that identity needs.

    CL_int [uL/min/mg protein] = 1000 * Vmax [nmol/min/mg protein] / Km [uM]

So the screening question for a new paper is exactly two questions. Does it
report a Km, in uM, for a NAMED transporter and a NAMED compound? And does it
report a Vmax on a basis that converts to nmol/min/mg protein -- or, failing
that, does it report Vmax/Km itself in uL/min/mg protein? Nothing else about
the paper decides this. Not how careful it is, not whether it agrees with
anyone, not how many compounds it covers.

THE FOUR TIERS, which build_transporter_working_table.py assigns:

    A   reaches CL_int per mg protein. Km + Vmax, or a directly reported
        Vmax/Km in uL/min/mg. These are the only rows that may be compared
        with each other as clearances.
    B   reports a clearance on a basis that will not join. Han 2008's
        uL/min per 10^6 hepatocytes is the only instance in this record;
        joining it needs mg protein per 10^6 cells, which the paper does not
        give. Kept, visible, and out of the pooled comparison.
    C   affinity only: a Km and no clearance anywhere. Still fully usable for
        the saturation question -- whether half-saturation sits above or
        below human serum concentrations is a statement about Km alone and
        needs no Vmax. Not usable for clearance.
    D   no constant: substrate calls, uptake ratios, inhibition. Written to
        db/transporter_direction.csv. This is the direction evidence, and it
        is evidence of a different kind, not weaker evidence of the same
        kind.

WHY IC50 AND Ki SIT IN D RATHER THAN C. An IC50 is as much a property of the
probe substrate and the concentration it was used at as of the inhibitor
being measured, so two IC50 values from different probes are not comparable
quantities. Turning one into a probe-independent Ki needs the probe's own Km,
through Cheng-Prusoff. Of the 20 numeric IC50 values in this record, NONE
carries the probe's Km. So none can be converted, and pooling them would be
the same error as comparing binding constants measured against different
fluorescent probes -- which is exactly what the albumin and L-FABP work spent
its effort avoiding.

WHY A SUBSTRATE CALL IS NOT A WEAK Km. Argoul 2026 contributes 120 rows and
no constants: it reports whether uptake exceeded twice control at 0.2 and
2 uM. That is a real measurement and it answers a real question -- which
transporters touch which compounds, in which species -- but no arithmetic
turns it into an affinity or a clearance. Keeping it in the same sheet as the
constants is how, eventually, somebody averages the two.

Usage:
    python3 -I scripts/screen_transporter_fits.py           # every source
    python3 -I scripts/screen_transporter_fits.py Yang Niu  # just these
"""
import collections, csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, "..", "db", "transporter_kinetics_working.csv")
DIR = os.path.join(HERE, "..", "db", "transporter_direction.csv")

# Sources already ruled on. They pass or fail on their merits, but the
# decision is made and the screen should not re-propose them every run.
DECIDED = {
    "Sakolish": "excluded: IVIVE-predicted relative clearances, not "
                "measurements, on a composite RPTEC entity with no single "
                "transporter identity to join a direction to",
    "Lin": "excluded: same composite RPTEC IVIVE output as Sakolish",
    "Liu": "out of scope: carboxylesterase CES1/CES2, an enzyme, not a "
           "transporter. In the record because it was on a requested list",
    "Worley": "excluded: Vmax_basoC and Vmax_apicalC are whole-body PBPK "
              "fitted parameters in mg/h/kg BW^0.75, not in vitro rates",
    "Loccisano": "excluded: PBPK Tm/KT, a fitted scalar, not a measurement",
    "Kim": "excluded: PBPK Tm/KT, a fitted scalar, not a measurement",
    "Ebert": "different quantity: membrane/water partitioning and passive "
             "permeability, held in db/membrane_partitioning.csv. Its "
             "transporter content is one intro paragraph and one comparison "
             "against the CONTROL cells of other people's assays",
    "Argoul": "tier D by design: 120 uptake-ratio rows, no constants. The "
              "species-difference backbone, and none of it is a clearance",
    "Ryu": "tier D: substrate classification across 14 PFAS",
    "Kudo": "tier D: in vivo clearance and mRNA, no in vitro constants",
    "Cheng": "tier D: mRNA expression and sex predominance",
    "Buist": "tier D: mRNA expression and sex predominance",
    "Kummu": "tier D: placental transfer correlation, n = 4",
}

# What each tier-C or tier-B source would need to reach tier A. This is the
# request list, derived from the data rather than from memory.
#
# (rows that would move, what is needed). The row count is STORED, not parsed
# back out of the sentence: an earlier version of this script recovered it
# with a digit search over the prose and picked up the "37" in "OAT4 37",
# reporting that three obtained papers would take tier A to 92 rather than
# 58. That is the same mistake the working-table builder made when it read a
# Km of 1 out of "reported only in Fig 1f", committed here within the hour.
BLOCKER = {
    ("Yang", "2010"): (3, "its own Vmax. The record holds OAT4 37 and URAT1 0.32 "
                      "ONLY as Louisse 2023's introduction quotes them, and "
                      "pairing a quoted Vmax with a Km read from the primary "
                      "builds a clearance out of two papers"),
    ("Yang", "2009"): (1, "a Vmax for Oatp1a1/PFOA. The paper reports Km 162.2 "
                      "and no Vmax in what has been obtained"),
    ("Han", "2008"): (2, "mg protein per 10^6 rat hepatocytes, which converts its "
                     "clearances (64.8 and 47.6 uL/min/10^6 cells) onto the "
                     "pooled basis. A standard value would do it, but it is "
                     "not in this paper"),
}


def load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def short(study):
    return (study or "?").split()[0].rstrip(",")


def main():
    want = [a.lower() for a in sys.argv[1:]]
    kin, dirn = load(KIN), load(DIR)
    rows = kin + dirn

    by = collections.defaultdict(list)
    for r in rows:
        by[(short(r["study"]), r["year"])].append(r)

    print("=" * 78)
    print("CL_int [uL/min/mg protein] = 1000 * Vmax [nmol/min/mg] / Km [uM]")
    print("A source reaches the kinetics table only if it supplies both.")
    print("=" * 78)

    order = {"A": 0, "B": 1, "C": 2, "D": 3}
    best = {}
    for key, rs in by.items():
        best[key] = min((r["tier"] for r in rs), key=lambda t: order[t])

    for tier, label in (("A", "TIER A - reaches CL_int per mg protein"),
                        ("B", "TIER B - clearance on a basis that will not join"),
                        ("C", "TIER C - affinity only (Km, no clearance)"),
                        ("D", "TIER D - direction evidence, no constant")):
        group = sorted(k for k, v in best.items() if v == tier)
        if not group:
            continue
        print(f"\n{label}   [{len(group)} sources]")
        for key in group:
            name, year = key
            if want and not any(w in name.lower() for w in want):
                continue
            rs = by[key]
            n = collections.Counter(r["tier"] for r in rs)
            counts = " ".join(f"{t}:{n[t]}" for t in "ABCD" if n[t])
            print(f"   {name + ' ' + year:22s} {len(rs):3d} rows  ({counts})")
            if key in BLOCKER:
                n_gain, why = BLOCKER[key]
                print(f"       NEEDS: {why}")
                print(f"       -> {n_gain} row(s) would move to tier A")
            if name in DECIDED and tier != "A":
                print(f"       {DECIDED[name]}")

    print("\n" + "=" * 78)
    print("WORKLIST - what would move rows up a tier")
    print("=" * 78)
    gain = 0
    for key, (n_gain, why) in sorted(BLOCKER.items()):
        if key in by:
            print(f"   {key[0]} {key[1]} (+{n_gain}): {why}")
            gain += n_gain
    a = sum(1 for r in kin if r["tier"] == "A")
    print(f"\n   tier A now {a} rows; all three obtained would make it "
          f"{a + gain}.")
    print(f"   Nothing else in the record is one measurement away from "
          f"tier A.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
