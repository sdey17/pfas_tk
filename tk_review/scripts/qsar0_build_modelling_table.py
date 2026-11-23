#!/usr/bin/env python3
"""QSAR Stage 0: the modelling table, from the transporter strand's endpoint.

    python3 -I scripts/qsar0_build_modelling_table.py

Writes db/qsar/qsar_modelling_table.csv and db/qsar/qsar_descriptors.csv.
Idempotent.

WHAT THIS REPLACES. db/qsar/qsar_endpoint_table.csv was a prototype: 11 rows,
mouse only, from Argoul 2026, with the endpoint already correctly chosen as
R = CL_R / (fu x GFR) -- the same quantity the transporter strand calls T. The
strand has since built that endpoint for three species over 44 verdicts with
graded evidence, so the prototype is superseded rather than extended, and kept
on disk for provenance.

THE ENDPOINT, AND WHY IT IS THE RIGHT ONE FOR A QSAR

Renal clearance factors exactly as CL_R = fu x GFR x T. GFR is a species
constant and fu has its own well-populated QSAR literature, so T is the part
no existing model predicts. The transporter strand showed why that matters:
binding explains 61% of the log-variance in the human clearance series and 83%
in the mouse, but the tubular residual still spans 466-fold and 823-fold.
A QSAR FOR RENAL CLEARANCE THAT MODELS ONLY BINDING LEAVES THAT ON THE TABLE.

THREE THINGS THIS TABLE DOES THAT A NAIVE ENDPOINT TABLE WOULD NOT

1. It carries the EVIDENCE TIER, so censored rows are usable. Most rows are
   one-sided bounds, not point values. A model fitted on bounds as though they
   were points is biased; a model that drops them loses most of the data. Both
   are avoidable -- censored regression uses them correctly -- but only if the
   censoring direction travels with the row, which it does here in
   `censoring`.

2. It carries the fu SOURCE and a flag for whether the sign survives changing
   it. Two compounds (PFPeA, PFHxA) flip between reabsorption and secretion
   depending on which human binding dataset is used. Those rows should not
   silently anchor a model.

3. It keeps fu and CL_R as separate columns beside T, so the same table
   supports three different modelling targets without rebuilding: T (the
   tubular question), fu (the binding question, already well served), and
   CL_R (the composite, which is what regulators actually want).
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
OUT = HERE / "db" / "qsar"

# ---- descriptors, assigned here from structure ---------------------------
# chemical: (n_C_total, n_perfluorinated_C, head_group, n_ether_O, n_branch, MW)
D = {
    "PFBA":   (4, 3, "carboxylate", 0, 0, 214.04),
    "PFPeA":  (5, 4, "carboxylate", 0, 0, 264.05),
    "PFHxA":  (6, 5, "carboxylate", 0, 0, 314.05),
    "PFHpA":  (7, 6, "carboxylate", 0, 0, 364.06),
    "PFOA":   (8, 7, "carboxylate", 0, 0, 414.07),
    "PFNA":   (9, 8, "carboxylate", 0, 0, 464.08),
    "PFDA":   (10, 9, "carboxylate", 0, 0, 514.08),
    "PFUnDA": (11, 10, "carboxylate", 0, 0, 564.09),
    "PFDoDA": (12, 11, "carboxylate", 0, 0, 614.10),
    "PFBS":   (4, 4, "sulfonate", 0, 0, 300.10),
    "PFHxS":  (6, 6, "sulfonate", 0, 0, 400.12),
    "PFOS":   (8, 8, "sulfonate", 0, 0, 500.13),
    "PFOSA":  (8, 8, "sulfonamide", 0, 0, 499.14),
    "6:2 FTS": (8, 6, "fluorotelomer sulfonate", 0, 0, 428.17),
}
# h_logD and pKa where Ryu 2024 supplies them
EXTRA = {}


def load_extra():
    f = HERE / "db" / "ryu2024_funbound.csv"
    if not f.exists():
        return
    for r in csv.DictReader(f.open()):
        c = r["chemical"]
        if c in EXTRA:
            continue
        EXTRA[c] = (r.get("h_logD", ""), r.get("pKa", ""))


def main():
    load_extra()
    OUT.mkdir(exist_ok=True)

    # ---- descriptor table -------------------------------------------
    pd = OUT / "qsar_descriptors.csv"
    with pd.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_C_total", "n_perfluorinated_C",
                    "head_group", "n_ether_O", "n_branch", "mw_g_mol",
                    "h_logD", "pKa", "descriptor_source"])
        for c, (nc, npf, hg, no, nb, mw) in sorted(D.items(),
                                                   key=lambda kv: kv[1][0]):
            ld, pk = EXTRA.get(c, ("", ""))
            w.writerow([c, nc, npf, hg, no, nb, mw, ld, pk,
                        "structure, assigned here; h_logD and pKa from "
                        "Ryu 2024 where available"])

    # ---- modelling table --------------------------------------------
    cons = list(csv.DictReader(
        (HERE / "db" / "tubular_factor_consolidated.csv").open()))
    robust = {}
    for r in csv.DictReader(
            (HERE / "db" / "human_fu_recommended.csv").open()):
        robust[r["chemical"]] = r["tubular_sign_robust_to_fu_choice"]

    pm = OUT / "qsar_modelling_table.csv"
    n = usable = 0
    with pm.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "species", "sex", "chemical", "n_C_total", "n_perfluorinated_C",
            "head_group", "n_ether_O", "n_branch", "mw_g_mol", "h_logD",
            "pKa", "fu", "fu_matrix", "gfr_mL_kg_day",
            "filtered_load_mL_kg_day", "T_low", "T_high", "T_point",
            "censoring", "evidence_tier", "verdict",
            "log10_T_point", "usable_for_regression",
            "sign_robust_to_fu_source", "clearance_sources"])
        for r in cons:
            c = r["chemical"]
            if c not in D:
                continue
            nc, npf, hg, no, nb, mw = D[c]
            ld, pk = EXTRA.get(c, ("", ""))
            lo = r["T_low"] or ""
            hi = r["T_high"] or ""
            tier = r["tier"]
            # censoring: DIRECT with one source -> a point; UPPER -> right-
            # censored from above (T <= hi); LOWER -> left-censored
            if tier == "DIRECT" and lo and hi and abs(
                    float(lo) - float(hi)) < 1e-9:
                point, cens = float(lo), "none (point estimate)"
            elif tier == "DIRECT":
                point, cens = (float(lo) + float(hi)) / 2, \
                    "interval (several direct sources)"
            elif tier == "UPPER":
                point, cens = float(hi), "right-censored: T <= value"
            elif tier == "LOWER":
                point, cens = float(lo), "left-censored: T >= value"
            else:
                point, cens = None, "unknown"
            import math
            use = "yes" if point and point > 0 else "no"
            if r["verdict"] == "undetermined":
                use = "no (bound permits either sign)"
            if use == "yes":
                usable += 1
            w.writerow([
                r["species"], r["sex"], c, nc, npf, hg, no, nb, mw, ld, pk,
                r["fu"], r["fu_matrix"], r["gfr_mL_kg_day"],
                r["filtered_load_mL_kg_day"], lo, hi,
                round(point, 5) if point else "", cens, tier, r["verdict"],
                round(math.log10(point), 4) if point and point > 0 else "",
                use, robust.get(c, ""), r["clearance_sources"][:90]])
            n += 1

    print(f"wrote {pd.relative_to(HERE)}  ({len(D)} compounds)")
    print(f"wrote {pm.relative_to(HERE)}  ({n} rows, {usable} usable)")
    report(pm)


def report(pm):
    import collections
    rows = list(csv.DictReader(pm.open()))
    bar = "=" * 78
    print()
    print(bar)
    print("WHAT A MODEL CAN ACTUALLY BE FITTED ON")
    print(bar)
    cens = collections.Counter(r["censoring"] for r in rows)
    print("  censoring structure:")
    for k, v in cens.most_common():
        print(f"    {v:>3}  {k}")
    print()
    sp = collections.Counter(r["species"] for r in rows
                             if r["usable_for_regression"] == "yes")
    print("  usable rows by species:", dict(sp))
    print()
    print("  MOST OF THE DATA IS CENSORED, AND THAT IS THE DEFINING FEATURE")
    print("  OF THIS DATASET. Fitting ordinary least squares to the bounds as")
    print("  if they were measurements would bias every coefficient toward")
    print("  the censoring direction. A Tobit or interval-censored likelihood")
    print("  uses them correctly and is the right default here.")
    print()
    print(bar)
    print("THE HONEST LIMIT: n IS SMALL AND THE DESCRIPTORS ARE COLLINEAR")
    print(bar)
    import math
    use = [r for r in rows if r["usable_for_regression"] == "yes"]
    nc = [float(r["n_C_total"]) for r in use]
    mw = [float(r["mw_g_mol"]) for r in use]
    m1, m2 = sum(nc)/len(nc), sum(mw)/len(mw)
    num = sum((a-m1)*(b-m2) for a, b in zip(nc, mw))
    den = math.sqrt(sum((a-m1)**2 for a in nc)*sum((b-m2)**2 for b in mw))
    print(f"  usable rows: {len(use)}")
    print(f"  carbon number vs molecular weight:  r = {num/den:.4f}")
    print()
    print("  CHAIN LENGTH AND MOLECULAR WEIGHT ARE THE SAME DESCRIPTOR HERE.")
    print("  For a homologous series they must be -- MW is linear in nC. Any")
    print("  model using both is fitting one variable twice. The descriptors")
    print("  that carry INDEPENDENT information are the head group, the ether")
    print("  oxygen count, branching, and pKa; and of those only head group")
    print("  has more than a couple of levels in this set.")
    print()
    print("  So the realistic first model is small: log10 T against chain")
    print("  length and head group, fitted with a censored likelihood,")
    print("  per species. Anything larger than that is fitting noise, and")
    print("  the table is built to make that limit visible rather than to")
    print("  hide it.")


if __name__ == "__main__":
    main()
