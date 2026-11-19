#!/usr/bin/env python3
"""Stage 3: intrinsic clearance, scaled where it can be and labelled where not.

    python3 -I scripts/build_stage3_clearance.py

Writes db/transporter_clearance_stage3.csv and prints the analysis.

FIVE LAYERS, KEPT SEPARATE ON PURPOSE. The point of this stage is not to
produce one number per row; it is to make visible which parts of a renal
clearance estimate are measured and which are assumed.

  1. CL_int per mg protein, unscaled -- the 52 tier-A rows, untouched.
  2. Scaled to kidney by REF and by RAF, SEPARATELY and never merged.
  3. Direction from db/transporter_localisation.csv, so a clearance carries
     its sign.
  4. An expression gate from db/transporter_expression.csv, so a transporter
     that is not detected in the tissue cannot contribute to that species'
     picture however fast it looks in vitro.
  5. The albumin caveat from db/ivive_scaling_factors.csv, as a stated range
     rather than an applied multiplier.

SCALING IS AVAILABLE FOR THREE TRANSPORTERS AND ONE SPECIES. Kumar 2021's
factors are human OAT1, OAT2 and OAT3 in OAT-overexpressing HEK293 against
human kidney. Both REF and RAF are defined in vitro over in vivo, so an in
vitro clearance is DIVIDED by the factor:

    REF   OAT1 0.16 -> x6.25     RAF   OAT1 0.64 -> x1.56
          OAT2 0.15 -> x6.67           OAT2 7.3  -> x0.137
          OAT3 0.37 -> x2.70           OAT3 4.1  -> x0.244

For OAT3 the two methods disagree about the SIGN of the correction: REF
raises the clearance 2.7-fold, RAF lowers it 4.1-fold. Both are reported.
Anything that is not human OAT1/OAT3 in a HEK293-type system stays unscaled,
with the reason recorded -- which is most of the reabsorptive side, and that
is the honest limit of this stage rather than an oversight.

WHAT THIS STAGE DELIBERATELY DOES NOT DO: convert to a whole-organ clearance.
That needs membrane protein per gram of kidney and kidney mass per kg body
weight. The first is not in this database; Scotcher 2017 gives microsomal and
cytosolic protein per gram of human kidney cortex, which is the wrong scalar
for transporters quantified in a total membrane fraction. Inventing it would
put a fabricated number under every row downstream, so the stage stops at
uL/min/mg protein and names the missing factor.
"""
import collections, csv, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, "..", "db", "transporter_kinetics_working.csv")
LOC = os.path.join(HERE, "..", "db", "transporter_localisation.csv")
EXP = os.path.join(HERE, "..", "db", "transporter_expression.csv")
OUT = os.path.join(HERE, "..", "db", "transporter_clearance_stage3.csv")

sys.path.insert(0, HERE)
from join_transporter_direction import canon, tissue_of   # noqa: E402

# transporter -> (REF, RAF), both in vitro / in vivo, human, HEK293.
FACTORS = {"OAT1": (0.16, 0.64), "OAT2": (0.15, 7.3), "OAT3": (0.37, 4.1)}
HEK = ("HEK293", "HEK 293")
ALBUMIN = ("OAT1/OAT3 uptake of highly protein-bound substrates rose 3.5-122x "
           "(OAT1) and up to 28x (OAT3) with 4% HSA present at constant "
           "unbound concentration (Tan 2024). PFOA fu = 0.0017, at or below "
           "the range tested. This assay carried no protein, so the value is "
           "plausibly a lower bound")
FIELDS = ["study", "year", "transporter", "gene", "species", "tissue",
          "chemical", "n_perfluorinated_carbons", "head_group",
          "km_uM", "vmax_nmol_min_mg", "clint_uL_min_mg",
          "membrane", "net_process", "effect_on_half_life",
          "tissue_abundance_pmol_mg", "expression_gate",
          "ref_factor", "clint_ref_scaled", "raf_factor", "clint_raf_scaled",
          "scaling_note", "albumin_caveat", "expression_system"]


def load(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main():
    kin = [r for r in load(KIN) if r["tier"] == "A"]
    loc, exp = load(LOC), load(EXP)

    auth = {}
    for r in loc:
        auth[(canon(r["transporter"]) or r["transporter"], r["species"],
              r["tissue"])] = r
    # measured tissue abundance, per (canonical transporter, species)
    abund = {}
    for r in exp:
        if r["units"] != "pmol/mg membrane protein" or r["sex"] != "both":
            continue
        c = canon(r["transporter"]) or r["transporter"]
        abund[(c, r["species"])] = (r["value"], r["note"])

    out = []
    for r in kin:
        c = canon(r["transporter"]) or r["transporter"]
        sp = (r["species_orthologue"] or "").strip()
        t = tissue_of(r)
        a = auth.get((c, sp, t))
        if a is None:
            cand = [v for k, v in auth.items() if k[0] == c and k[1] == sp]
            a = cand[0] if len(cand) == 1 else None
        cl = float(r["clint_uL_min_mg"])
        sysm = r["expression_system"] or ""

        # --- expression gate ---
        gate, ab_val, ab_note = "not assessed", "", ""
        if t == "kidney proximal tubule":
            if (c, sp) in abund:
                ab_val, ab_note = abund[(c, sp)]
                gate = "PASS" if ab_val else "FAIL - not detected in this tissue"
            else:
                gate = "unknown - not in the Thakur 2024 panel"
        else:
            gate = ("not assessed - Thakur 2024 measured kidney only, and "
                    "this row is not renal")

        # --- scaling ---
        ref = raf = ref_v = raf_v = ""
        if sp == "human" and c in FACTORS and any(h in sysm for h in HEK):
            ref_f, raf_f = FACTORS[c]
            ref, raf = ref_f, raf_f
            ref_v, raf_v = round(cl / ref_f, 1), round(cl / raf_f, 1)
            note = ("scaled: human {} in a HEK293-type system. Kumar 2021's "
                    "factors come from a different HEK293 line, so these are "
                    "estimates, and REF and RAF disagree {:.0f}-fold here"
                    .format(c, max(ref_f, raf_f) / min(ref_f, raf_f)))
        elif sp != "human":
            note = ("NOT SCALED: the available factors are human kidney "
                    "references; no rodent equivalent exists")
        elif c not in FACTORS:
            note = ("NOT SCALED: no REF or RAF published for {}. RAF would "
                    "need a transporter-selective probe, which does not "
                    "exist for it".format(c))
        else:
            note = ("NOT SCALED: factors are for OAT-overexpressing HEK293 "
                    "and this measurement is in {}".format(sysm[:50]))

        d = {f: "" for f in FIELDS}
        d.update({
            "study": r["study"].split()[0] + " " + r["year"],
            "year": r["year"], "transporter": r["transporter"],
            "gene": r["gene"], "species": sp, "tissue": t,
            "chemical": r["chemical"],
            "n_perfluorinated_carbons": r["n_perfluorinated_carbons"],
            "head_group": r["head_group"], "km_uM": r["km_uM"],
            "vmax_nmol_min_mg": r["vmax_nmol_min_mg"],
            "clint_uL_min_mg": cl,
            "membrane": a["membrane"] if a else "",
            "net_process": a["net_process"] if a else "unresolved",
            "effect_on_half_life": a["effect_on_half_life"] if a else "",
            "tissue_abundance_pmol_mg": ab_val,
            "expression_gate": gate,
            "ref_factor": ref, "clint_ref_scaled": ref_v,
            "raf_factor": raf, "clint_raf_scaled": raf_v,
            "scaling_note": note,
            "albumin_caveat": (ALBUMIN if c in ("OAT1", "OAT3")
                               else "not assessed for this transporter"),
            "expression_system": sysm,
        })
        out.append(d)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader(); w.writerows(out)

    # ----------------------------- report -----------------------------
    sc = [d for d in out if d["clint_ref_scaled"] != ""]
    print("=" * 76)
    print(f"STAGE 3   {len(out)} tier-A rows")
    print("=" * 76)
    print(f"  scaled to kidney   {len(sc):3d}  (human OAT1/OAT3 in HEK293)")
    print(f"  unscaled           {len(out) - len(sc):3d}")
    for why, n in collections.Counter(
            d["scaling_note"].split(": ", 1)[-1].split(".")[0]
            for d in out if d["clint_ref_scaled"] == "").most_common():
        print(f"       {n:3d}  {why}")

    print("\n" + "=" * 76)
    print("THE SCALED ROWS, both methods")
    print("=" * 76)
    print(f"{'transporter':8s} {'compound':8s} {'study':16s} "
          f"{'measured':>9s} {'REF':>8s} {'RAF':>8s}  process")
    for d in sorted(sc, key=lambda d: (d["transporter"], d["chemical"])):
        print(f"{d['transporter']:8s} {d['chemical']:8s} {d['study']:16s} "
              f"{d['clint_uL_min_mg']:9.1f} {d['clint_ref_scaled']:8.1f} "
              f"{d['clint_raf_scaled']:8.1f}  {d['net_process'][:26]}")

    print("\n" + "=" * 76)
    print("EXPRESSION GATE")
    print("=" * 76)
    for g, n in collections.Counter(
            d["expression_gate"].split(" - ")[0] for d in out).most_common():
        print(f"   {n:3d}  {g}")
    failed = [d for d in out if d["expression_gate"].startswith("FAIL")]
    if failed:
        print("\n   rows removed from their species' renal picture:")
        for d in sorted(failed, key=lambda d: -d["clint_uL_min_mg"]):
            print(f"      {d['clint_uL_min_mg']:7.1f} uL/min/mg  "
                  f"{d['transporter']:6s} {d['species']:5s} "
                  f"{d['chemical']:8s} {d['study']}")

    print("\n" + "=" * 76)
    print("DIRECTION, measured CL_int, whole-cell assays, expression-gated")
    print("=" * 76)
    keep = [d for d in out
            if not d["expression_gate"].startswith("FAIL")
            and "step 2" not in d["net_process"]]
    g = collections.defaultdict(list)
    for d in keep:
        g[d["effect_on_half_life"].split(" (")[0]].append(
            d["clint_uL_min_mg"])
    for k, v in sorted(g.items(), key=lambda x: -len(x[1])):
        print(f"   {k:12s} n={len(v):3d}  median {statistics.median(v):7.1f}"
              f"  range {min(v):g}-{max(v):g}")
    print(f"\n   {len(out) - len(keep)} rows excluded: failed the expression "
          f"gate or are apical efflux measured in vesicles")

    print("\n" + "=" * 76)
    print("REABSORPTION AGAINST SECRETION, PER COMPOUND")
    print("=" * 76)
    print("  Human, Louisse 2023 (OAT4) and 2024 (OAT1/OAT3): same group,")
    print("  same system, so the comparison is within-laboratory. Measured")
    print("  CL_int, unscaled, uL/min/mg protein.")
    nc = {"PFHpA": 6, "PFOA": 7, "PFNA": 8, "PFDA": 9, "PFHxS": 6,
          "PFOS": 8, "PFBS": 4, "HFPO-TA": 8}
    by = collections.defaultdict(dict)
    for d in out:
        if d["species"] == "human":
            by[d["transporter"]][d["chemical"]] = d["clint_uL_min_mg"]
    print(f"\n  {'compound':9s} {'OAT4':>7s} {'OAT3':>7s} {'OAT1':>7s}   "
          f"{'OAT4/OAT3':>9s} {'OAT4/OAT1':>9s}")
    ratios = []
    for c in sorted(nc, key=lambda c: nc[c]):
        o4 = by.get("OAT4", {}).get(c)
        o3 = by.get("OAT3", {}).get(c)
        o1 = by.get("OAT1", {}).get(c)
        if o4 is None and o3 is None and o1 is None:
            continue
        r43 = f"{o4 / o3:9.1f}" if (o4 and o3) else " " * 9
        r41 = f"{o4 / o1:9.1f}" if (o4 and o1) else " " * 9
        if o4 and o3:
            ratios.append(o4 / o3)
        fmt = lambda v: f"{v:7.1f}" if v else "      -"
        print(f"  C{nc[c]} {c:6s} {fmt(o4)} {fmt(o3)} {fmt(o1)}   "
              f"{r43} {r41}")
    if ratios:
        print(f"\n  BOTH SIDES RISE WITH CHAIN LENGTH -- OAT3 2.7-fold from C6")
        print("  to C8, OAT4 2.1-fold from C6 to C9 -- so chain length does")
        print("  not by itself tip the balance: the reabsorption-to-secretion")
        print(f"  ratio stays between {min(ratios):.1f} and {max(ratios):.1f} "
              f"across C6-C8.")
        HEAD = {"PFHpA": "carboxylate", "PFOA": "carboxylate",
                "PFNA": "carboxylate", "PFDA": "carboxylate",
                "PFHxS": "sulfonate", "PFOS": "sulfonate"}
        cx, su = [], []
        for c in by.get("OAT4", {}):
            for sec in ("OAT3", "OAT1"):
                v = by.get(sec, {}).get(c)
                if v:
                    (cx if HEAD.get(c) == "carboxylate" else su).append(
                        (by["OAT4"][c] / v, c, sec))
        if cx:
            print(f"\n  WHAT STANDS OUT IS THE LEVEL, AND IT SPLITS BY HEAD")
            print(f"  GROUP. For every CARBOXYLATE measured, OAT4 carries")
            print(f"  {min(r for r, _, _ in cx):.1f}-{max(r for r, _, _ in cx):.1f} "
                  f"times the intrinsic clearance of a secretory")
            print("  transporter -- the in vitro shape of the 97-99%")
            print("  reabsorption OEHHA Table A6.4 reports in vivo.")
        if su:
            print("\n  The SULFONATES do not behave alike:")
            for r, c, sec in sorted(su, key=lambda x: -x[0]):
                print(f"     {c:6s} OAT4/{sec} = {r:.1f}")
            print("  PFHxS is reabsorption-dominated like the carboxylates;")
            print("  PFOS is the one compound whose secretion keeps pace with")
            print("  its reabsorption, and OAT3 does not transport it at all.")
        print("\n  AND OAT4 IS THE HUMAN-SPECIFIC ONE (section 11): the")
        print("  transporter that dominates human renal PFAS handling by")
        print("  intrinsic clearance is the one the rat does not have.")
        print("\n  THE CAVEAT THAT DECIDES HOW FAR THIS GOES. OAT4 cannot be")
        print("  scaled and OAT1/OAT3 can. Applying REF lifts OAT1 to")
        print("  118-334 uL/min/mg, which would overturn the ranking. So the")
        print("  unscaled comparison is suggestive, not decisive, and the one")
        print("  measurement that would settle it is an OAT4 abundance in the")
        print("  expressing cells -- a REF for OAT4.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
