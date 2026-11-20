#!/usr/bin/env python3
"""Is the renal axis QSAR-able? (Stage 6, part 2)

    python3 -I scripts/stage6_chain_length_bridge.py

Reads db/tubular_factor_evidence.csv. Writes db/renal_axis_decomposition.csv.
Idempotent. Prints the analysis.

THE QUESTION AIM 2 ACTUALLY ASKS. Renal clearance factors exactly:

    CL_R  =  fu x GFR x T

GFR is a species constant, so across a compound series within one species all
the variation sits in two terms: PLASMA BINDING (fu) and TUBULAR HANDLING (T).
If fu carries almost all of it, then a PFAS QSAR for renal clearance is really
a QSAR for albumin binding and the tubule can be treated as near-constant. If
T carries a comparable share, the QSAR needs transporter descriptors too, and
the two aims cannot be separated.

This is the first question in the strand that can be asked of the whole
series at once, because Stage 6 part 1 put every compound on the same basis.

WHAT THE BOUNDS DO AND DO NOT ALLOW. Most human T values are upper bounds from
total clearance, and the bound LOOSENS WITH CHAIN LENGTH because faecal
elimination grows with chain length while the bound charges all of it to the
kidney. So any apparent rise in T at the long end is an artefact and is
reported as such. No human renal-fraction measurement exists in this project
to tighten it. Two human values are direct (PFHxA and PFBS, Argoul 2026) and
all seven mouse values are direct, so the comparisons that matter are built on
those.
"""
import csv, math, pathlib, statistics as st

HERE = pathlib.Path(__file__).resolve().parent.parent
NC = {"PFBA": 4, "PFPeA": 5, "PFHxA": 6, "PFHpA": 7, "PFOA": 8, "PFNA": 9,
      "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12, "PFBS": 4, "PFHxS": 6,
      "PFOS": 8, "PFOSA": 8, "6:2 FTS": 6}
FAM = {c: ("sulfonate" if c in ("PFBS", "PFHxS", "PFOS") else
           "sulfonamide" if c == "PFOSA" else
           "fluorotelomer" if c == "6:2 FTS" else "carboxylate")
       for c in NC}


def load():
    f = HERE / "db" / "tubular_factor_evidence.csv"
    out = []
    for r in csv.DictReader(f.open()):
        if not r["T"]:
            continue
        out.append(dict(species=r["species"], sex=r["sex"],
                        chemical=r["chemical"], kind=r["clearance_kind"],
                        fu=float(r["fu"]), filt=float(r["filtered_load_mL_kg_day"]),
                        cl=float(r["clearance_mL_kg_day"]), T=float(r["T"]),
                        src=r["clearance_source"]))
    return out


def pearson_log(xs, ys):
    n = len(xs)
    lx = [math.log10(v) for v in xs]
    ly = [math.log10(v) for v in ys]
    mx, my = sum(lx) / n, sum(ly) / n
    num = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    den = math.sqrt(sum((a - mx) ** 2 for a in lx) *
                    sum((b - my) ** 2 for b in ly))
    return num / den


def decomp(label, d, kind_note):
    cl = [v["cl"] for v in d.values()]
    fg = [v["filt"] for v in d.values()]
    T = [v["T"] for v in d.values()]
    r = pearson_log(fg, cl)
    print(f"  {label}   n={len(d)}   ({kind_note})")
    print(f"    clearance spread       {max(cl)/min(cl):>9.0f}x"
          f"   sd log10 {st.pstdev([math.log10(x) for x in cl]):.2f}")
    print(f"    filtered-load spread   {max(fg)/min(fg):>9.0f}x"
          f"   sd log10 {st.pstdev([math.log10(x) for x in fg]):.2f}")
    print(f"    T spread (residual)    {max(T)/min(T):>9.0f}x"
          f"   sd log10 {st.pstdev([math.log10(x) for x in T]):.2f}")
    print(f"    log-log r = {r:+.3f},  r2 = {r*r:.3f}"
          f"  -> binding explains {100*r*r:.0f}% of the variance in clearance")
    print()
    return dict(n=len(d), r2=r * r,
                sd_cl=st.pstdev([math.log10(x) for x in cl]),
                sd_fg=st.pstdev([math.log10(x) for x in fg]),
                sd_T=st.pstdev([math.log10(x) for x in T]))


def main():
    ev = load()
    bar = "=" * 78

    def series(sp, srcfrag, kind):
        out = {}
        for e in ev:
            if (e["species"] == sp and srcfrag in e["src"]
                    and e["kind"] == kind):
                out[e["chemical"]] = e
        return out

    H = series("human", "Abraham 2024", "total")
    M = series("mouse", "Argoul 2026 Table 4", "renal")

    print(bar)
    print("A.  HOW MUCH OF THE CLEARANCE SERIES IS BINDING, AND HOW MUCH IS")
    print("    THE TUBULE?")
    print(bar)
    hs = decomp("HUMAN", H, "Abraham 2024 total clearance, upper bounds on T")
    ms = decomp("MOUSE", M, "Argoul 2026 renal clearance, direct")
    print("  SO BINDING IS THE LARGER TERM IN BOTH SPECIES, BUT NOT THE ONLY")
    print("  ONE. It accounts for 61% of the log-variance in the human series")
    print("  and 83% in the mouse. The tubular residual still spans 466-fold")
    print("  and 823-fold across the series, which is far too wide to treat")
    print("  as a constant.")
    print()
    print("  A QSAR FOR PFAS RENAL CLEARANCE THEREFORE CANNOT BE A QSAR FOR")
    print("  ALBUMIN BINDING ALONE. Aim 2 needs the transporter layer, which")
    print("  is what aim 1 has just built. The two aims are not separable and")
    print("  this is the quantity that says so.")

    print()
    print(bar)
    print("B.  T AGAINST CHAIN LENGTH -- AND WHY THE RIGHT ARM IS AN ARTEFACT")
    print(bar)
    print(f"  {'nC':>3}{'compound':>9}{'human T':>14}{'mouse T':>14}")
    for c in ("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
              "PFUnDA", "PFDoDA"):
        h = H.get(c)
        m = M.get(c)
        hs_ = ("-" if h is None else
               (f"{h['T']:.4f}" if h["kind"] == "renal"
                else f"<={h['T']:.4f}"))
        ms_ = ("-" if m is None else
               (f"{m['T']:.4f}" if m["kind"] == "renal"
                else f"<={m['T']:.4f}"))
        print(f"  {NC[c]:>3}{c:>9}{hs_:>14}{ms_:>14}")
    print()
    print("  The human column looks like a V with its floor at C8 and a sharp")
    print("  rise to 3.08 at C12. THE RISE IS NOT A TUBULAR EFFECT. Every")
    print("  human value except PFHxA and PFBS is an upper bound computed from")
    print("  TOTAL clearance, and faecal elimination grows with chain length,")
    print("  so the bound charges an increasingly large non-renal term to the")
    print("  kidney. Kudo 2001 measured faecal dominance from PFDA upward in")
    print("  the rat. No human renal-fraction measurement exists in this")
    print("  project, so the right arm cannot be tightened and no tubular")
    print("  claim is made for C9 and above in the human.")

    print()
    print(bar)
    print("C.  WHAT THE DIRECT VALUES DO ESTABLISH")
    print(bar)
    print("  A direct value above a bound, or a bound below a direct value,")
    print("  gives a real inequality. Three such pairs exist.")
    print()
    hx, pb = H["PFHxA"], None
    for e in ev:
        if e["species"] == "human" and e["chemical"] == "PFBS" \
                and e["kind"] == "renal":
            pb = e
    hxd = [e for e in ev if e["species"] == "human"
           and e["chemical"] == "PFHxA" and e["kind"] == "renal"][0]
    po = H["PFOA"]
    pba = H["PFBA"]
    print(f"  1. C6 -> C8 CARBOXYLATE, HUMAN.  T(PFHxA) = {hxd['T']:.3f} "
          f"measured directly,")
    print(f"     T(PFOA) <= {po['T']:.4f}.  So T falls by AT LEAST "
          f"{hxd['T']/po['T']:.0f}-fold from C6 to C8.")
    print("     The human kidney returns about half of filtered PFHxA and")
    print("     essentially all filtered PFOA. This is the chain-length")
    print("     effect, and it is real rather than bounded away.")
    print()
    print(f"  2. HEAD GROUP BEATS CHAIN LENGTH AT THE SHORT END. Both of the")
    print(f"     human direct values sit at the short end and they disagree")
    print(f"     by {hxd['T']/pb['T']:.0f}-fold IN THE DIRECTION OPPOSITE TO "
          f"CHAIN LENGTH:")
    print(f"       PFBS  C4 sulfonate    T = {pb['T']:.4f}  "
          f"({100*(1-pb['T']):.1f}% reabsorbed)")
    print(f"       PFHxA C6 carboxylate  T = {hxd['T']:.4f}  "
          f"({100*(1-hxd['T']):.1f}% reabsorbed)")
    print("     The SHORTER compound is the better reabsorbed, by a wide")
    print("     margin, because it is a sulfonate. Any QSAR on chain length")
    print("     alone gets this pair backwards.")
    print()
    print(f"  3. AND THE HALF-LIFE DIFFERENCE IS NOT IN T. PFBS clears from")
    print(f"     human plasma in about a month and PFOA in years, yet their")
    print(f"     T values are {pb['T']:.4f} and <={po['T']:.4f} -- the same to")
    print(f"     within a factor of 2. The difference is the FILTERED LOAD:")
    print(f"     fu x GFR is {pb['filt']:.1f} for PFBS and {po['filt']:.1f} "
          f"for PFOA, a {pb['filt']/po['filt']:.0f}-fold gap.")
    print("     For this pair the tubule is not the discriminator; binding")
    print("     is. For the C6/C8 carboxylate pair it is the other way. THAT")
    print("     IS WHY NEITHER TERM CAN BE DROPPED.")

    print()
    print(bar)
    print("D.  RETRACTED: THERE IS NO ESTABLISHED SPECIES CONTRAST AT PFHxA")
    print(bar)
    print("  An earlier version of this section read human PFHxA (T = 0.48,")
    print("  reabsorption) against mouse PFHxA (T = 2.23, secretion) as the")
    print("  one clean sign flip between species in the direct data, and tied")
    print("  it to OAT4.")
    print()
    print("  scripts/stage6_abraham_partition.py withdraws it. The human")
    print("  value depends on whose binding measurement is used:")
    print()
    print("    fu 0.038 (Ryu 2024)        -> T = 0.48   reabsorption")
    print("    fu 0.020 (Argoul 2026)     -> T = 0.92   borderline")
    print("    fu 0.0068 (Abraham 2024,   -> T = 3.48   SECRETION")
    print("               the same subject as the clearance)")
    print()
    print("  On the subject's own binding the human secretes PFHxA as the")
    print("  mouse does, and the contrast disappears. The mouse verdict is")
    print("  unaffected -- it is secretion on both mouse binding sources, 2.23")
    print("  on Ryu's and 2.71 on Argoul's own -- but a one-sided result is")
    print("  not a species contrast.")
    print()
    print("  THE LESSON GENERALISES BEYOND THIS COMPOUND. Every sign verdict")
    print("  near T = 1 is a verdict about a binding measurement, because the")
    print("  filtered load scales with fu and nothing else in the expression")
    print("  is in doubt. The compounds whose verdicts are safe are the ones")
    print("  far from the line, which is every chain of C8 and above.")

    # ---- write the decomposition ------------------------------------
    p = HERE / "db" / "renal_axis_decomposition.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "chemical", "n_carbons", "family",
                    "fu", "filtered_load_mL_kg_day", "clearance_mL_kg_day",
                    "clearance_kind", "T", "T_is_bound", "source"])
        for sp, d in (("human", H), ("mouse", M)):
            for c, e in sorted(d.items(), key=lambda kv: NC[kv[0]]):
                w.writerow([sp, c, NC[c], FAM[c], e["fu"],
                            round(e["filt"], 3), e["cl"], e["kind"],
                            round(e["T"], 5),
                            "no" if e["kind"] == "renal" else "yes",
                            e["src"]])
    print(f"\nwrote {p.relative_to(HERE)}")
    print(f"  human r2 {hs['r2']:.3f}, mouse r2 {ms['r2']:.3f}")


if __name__ == "__main__":
    main()
