#!/usr/bin/env python3
"""Reabsorption or secretion? The sign of the dose response (Stage 21).

    python3 -I scripts/stage21_dose_sign_test.py

Writes db/dose_sign_test.csv and db/tubular_mechanism_verdict.csv. Idempotent.

THE PROBLEM THIS ADDRESSES. The tubular factor is not separable:
T = (1 - FR) + S/(f_u*GFR) has two unknowns and one equation, so a measured T
cannot say whether net reabsorption or net secretion dominates. The paper
records this as its principal limitation. It need not be.

THE TEST. Both processes are carrier-mediated and therefore saturable, but
they saturate in OPPOSITE directions:

    saturable REABSORPTION -> at high dose less is reclaimed -> T rises -> CL rises
    saturable SECRETION    -> at high dose less is pumped out -> T falls -> CL falls

So the SIGN of d(CL)/d(dose) discriminates, without needing to separate the
two terms numerically. Dose series exist for 24 species-sex-compound
combinations across five studies.

THE ASYMMETRY THAT MAKES THIS WORK. One confounder matters: saturable plasma
protein binding raises f_u at high dose, and CL_R = f_u*GFR*T rises even with
T unchanged. But binding saturation can only push clearance UP. It cannot
produce a fall. Therefore:

    CL rises  -> saturable reabsorption OR saturable binding  (ambiguous)
    CL falls  -> saturable secretion OR saturable non-renal   (reabsorption
                 dominance is EXCLUDED; binding cannot do it)

A fall is the diagnostic observation. A rise is diagnostic only where binding
saturation can be ruled out, which the serum concentrations allow for the NTP
series: a PFAS:albumin molar ratio well below 1 leaves binding unsaturated.

WHAT IS NOT CLAIMED. A fall in TOTAL clearance can also come from saturable
biliary or hepatic transport rather than saturable tubular secretion. For
PFOS, which undergoes substantial enterohepatic circulation, that alternative
is live and is recorded rather than dismissed.
"""
import csv
import math
import pathlib
import statistics as st

HERE = pathlib.Path(__file__).resolve().parent.parent

# Serum albumin molarity, mM. Rodent serum albumin ~30-40 g/L, MW 65.9 kDa.
ALBUMIN_mM = 0.60
ALBUMIN_RANGE_mM = (0.45, 0.70)

# Cmax at the low and high dose, mM, from db/dose_dependence.csv
# serum_conc_range, parsed here so the molar ratio is reproducible.
CMAX_mM = {
    ("Huang et al. (NTP)", "PFBS", "male"): (0.053, 0.750),
    ("Huang et al. (NTP)", "PFBS", "female"): (0.028, 0.422),
    ("Huang et al. (NTP)", "PFHxS", "male"): (0.08, 0.41),
    ("Huang et al. (NTP)", "PFHxS", "female"): (0.06, 0.30),
    ("Huang et al. (NTP)", "PFOS", "male"): (0.01, 0.21),
    ("Huang et al. (NTP)", "PFOS", "female"): (0.02, 0.27),
}

# Rows in dose_dependence.csv that are NOT dose escalations.
EXCLUDE_REASON = {
    "10 (as 2x5)": "repeat-dose (2 mg/kg x 5) against a single 20 mg/kg "
                   "dose, not a dose escalation; total dose differs by 2-fold "
                   "while the dosing SCHEDULE also differs, so the sign "
                   "carries no information about saturation",
}

TOTAL_CARBONS = {"PFHpA": 7, "PFOA": 8, "PFNA": 9, "PFDA": 10, "PFUnDA": 11,
                 "PFDoDA": 12, "PFTriDA": 13, "PFTeDA": 14,
                 "PFBS": 4, "PFHxS": 6, "PFOS": 8}

# From Stage 20: Kim 2016, one study, one strain, both sexes.
KIM_SEX = {"PFOA": (0.0993, 7.5410), "PFHxS": (0.0694, 5.8361),
           "PFOS": (0.1223, 0.1096)}

# Male rat tubular factors, OEHHA/Han GFR, from db/renal_only_v2.csv.
# PFOA, PFHxS and PFOS are Kim 2016 read at source (the weighted primary);
# PFHpA is Ohmori 2003 Fig. 2 digitised, the only source for it.
RAT_T_MALE = {"PFHpA": 15.0376, "PFOA": 0.0993,
              "PFHxS": 0.0694, "PFOS": 0.1223}


def collect():
    rows = []
    for r in csv.DictReader((HERE / "db" / "dose_dependence.csv").open()):
        if not (r["clearance_low"] and r["clearance_high"]):
            continue
        skip = next((v for k, v in EXCLUDE_REASON.items()
                     if k in r["dose_low"]), None)
        if skip:
            continue
        lo, hi = float(r["clearance_low"]), float(r["clearance_high"])
        rows.append(dict(
            study=f"{r['study']} {r['year']}".replace(" et al.", ""),
            chemical=r["chemical"], species=r["species"], sex=r["sex"],
            dose_low=r["dose_low"], dose_high=r["dose_high"],
            dose_units=r["dose_units"], cl_low=lo, cl_high=hi,
            cl_ratio=round(hi / lo, 3),
            cmax_key=(r["study"], r["chemical"], r["sex"])))
    for r in csv.DictReader(
            (HERE / "db" / "primary_2026" /
             "fujii2015_mouse_total_clearance.csv").open()):
        lo = float(r["cl_total_low_dose_mL_kg_day"])
        hi = float(r["cl_total_high_dose"])
        rows.append(dict(
            study="Fujii 2015", chemical=r["chemical"], species="mouse",
            sex=r["sex"], dose_low="0.31", dose_high="3.13",
            dose_units="umol/kg", cl_low=lo, cl_high=hi,
            cl_ratio=round(hi / lo, 3), cmax_key=None))
    return rows


def main():
    rows = collect()
    bar = "=" * 78

    print(bar)
    print("THE SIGN OF THE DOSE RESPONSE, 5 STUDIES")
    print(bar)
    print(f"  {'study':16}{'cmpd':8}{'sp':7}{'sex':8}{'CL lo':>9}{'CL hi':>9}"
          f"{'ratio':>8}  direction")
    for r in sorted(rows, key=lambda x: (x["chemical"], x["species"],
                                         x["sex"], x["study"])):
        d = ("RISES" if r["cl_ratio"] > 1.15 else
             "falls" if r["cl_ratio"] < 0.87 else "flat")
        r["direction"] = d
        print(f"  {r['study'][:15]:16}{r['chemical']:8}{r['species'][:6]:7}"
              f"{r['sex'][:7]:8}{r['cl_low']:>9.3g}{r['cl_high']:>9.3g}"
              f"{r['cl_ratio']:>8.2f}  {d}")
    print()
    nr = sum(1 for r in rows if r["direction"] == "RISES")
    nf = sum(1 for r in rows if r["direction"] == "falls")
    print(f"  {nr} rise, {nf} fall, {len(rows)-nr-nf} flat, of {len(rows)}.")
    print()

    print(bar)
    print("RULING THE BINDING CONFOUNDER IN OR OUT, WHERE Cmax IS REPORTED")
    print(bar)
    print(f"  Serum albumin taken as {ALBUMIN_mM} mM "
          f"(range {ALBUMIN_RANGE_mM[0]}-{ALBUMIN_RANGE_mM[1]} considered).")
    print(f"  A PFAS:albumin molar ratio near or above 1 at the high dose")
    print(f"  means binding saturation could itself raise clearance.")
    print()
    print(f"  {'cmpd':8}{'sex':8}{'Cmax lo':>9}{'Cmax hi':>9}"
          f"{'ratio:alb':>11}{'CL ratio':>10}  reading")
    conf = []
    for r in rows:
        if r["cmax_key"] not in CMAX_mM:
            continue
        clo, chi = CMAX_mM[r["cmax_key"]]
        mr = chi / ALBUMIN_mM
        if r["direction"] == "falls":
            read = "FALL: binding cannot cause this"
        elif mr >= 0.5:
            read = "rise CONFOUNDED (binding may saturate)"
        else:
            read = "rise, binding unsaturated"
        print(f"  {r['chemical']:8}{r['sex'][:7]:8}{clo:>9.3f}{chi:>9.3f}"
              f"{mr:>11.2f}{r['cl_ratio']:>10.2f}  {read}")
        r["cmax_low_mM"], r["cmax_high_mM"] = clo, chi
        r["pfas_to_albumin_at_high_dose"] = round(mr, 3)
        r["confounder_reading"] = read
        conf.append(r)
    print()
    print("  PFOS is the decisive row. Its molar ratio at the high dose is")
    print(f"  the LOWEST of the six "
          f"({min(r['pfas_to_albumin_at_high_dose'] for r in conf if r['chemical']=='PFOS'):.2f}"
          f"-{max(r['pfas_to_albumin_at_high_dose'] for r in conf if r['chemical']=='PFOS'):.2f}),"
          f" so binding is furthest from saturation,")
    print("  and its clearance nevertheless FELL. The one confounder that")
    print("  could manufacture a rise was weakest exactly where the fall was")
    print("  observed.")
    print()

    print(bar)
    print("THE MOUSE SERIES: THE DOSE EFFECT TURNS ON WITH CHAIN LENGTH")
    print(bar)
    fj = [r for r in rows if r["study"] == "Fujii 2015"]
    xs = [TOTAL_CARBONS[r["chemical"]] for r in fj]
    ys = [math.log10(r["cl_ratio"]) for r in fj]
    mx, my = st.fmean(xs), st.fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    icept = my - slope * mx
    res = [y - (icept + slope * x) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(e * e for e in res) / (len(xs) - 2) / sxx)
    r2 = 1 - sum(e * e for e in res) / sum((y - my) ** 2 for y in ys)
    print(f"  log10(CL ratio) = {icept:+.3f} {slope:+.4f} x total carbons")
    print(f"  slope {slope:+.4f} +/- {se:.4f}  t({len(xs)-2}) = {slope/se:.2f}"
          f"   r2 = {r2:.3f}")
    print(f"  crosses 1 at {-icept/slope:.1f} total carbons")
    print()
    print("  Saturable reabsorption predicts exactly this if the reabsorptive")
    print("  carrier binds longer chains more tightly: the longer the chain,")
    print("  the more reabsorption there is to saturate, so the larger the")
    print("  dose effect. Below ~C9 there is little to saturate and the")
    print("  ratio sits at or under 1.")
    print()
    print("  NOTE ON COVERAGE. The four PFTriDA and PFTeDA rows are marked")
    print("  unusable in the source file because no mouse plasma f_u exists")
    print("  for them, which blocks a tubular FACTOR. It does not block a")
    print("  SIGN test, which needs only the rodent clearances. They are")
    print("  included here and carry the top of the trend.")
    print()

    print(bar)
    print("T AND THE DOSE SIGN ARE JOINTLY INFORMATIVE")
    print(bar)
    print("  T says which process DOMINATES at the measured dose. The dose")
    print("  sign says which process SATURATES FIRST as dose rises. They are")
    print("  different questions, and the pair is more informative than")
    print("  either alone. A compound can be net-reabsorbing (T < 1) and")
    print("  still show a falling clearance if its smaller secretory")
    print("  component is the one that saturates.")
    print()
    print(f"  {'cmpd':8}{'sp':7}{'T male':>9}{'dominates':>14}"
          f"{'dose sign':>11}  joint reading")
    verdicts = []
    for chem in sorted({r["chemical"] for r in rows},
                       key=lambda c: (TOTAL_CARBONS[c], c)):
        for sp in ("rat", "mouse"):
            sub = [x for x in rows
                   if x["chemical"] == chem and x["species"] == sp]
            if not sub:
                continue
            up = sum(1 for r in sub if r["direction"] == "RISES")
            dn = sum(1 for r in sub if r["direction"] == "falls")
            conf_up = sum(1 for r in sub if r["direction"] == "RISES"
                          and "CONFOUNDED" in r.get("confounder_reading", ""))
            T = RAT_T_MALE.get(chem) if sp == "rat" else None
            dom = ("secretion" if T and T > 1 else
                   "reabsorption" if T else "unmeasured")
            sign = ("rises" if up > dn else "falls" if dn > up else "mixed/flat")
            if T is None:
                joint = f"no T for {sp}; sign alone: {sign}"
            elif T < 1 and sign == "rises":
                joint = "net reabsorption, and reabsorption saturates first"
            elif T < 1 and sign == "falls":
                joint = ("net reabsorption, but a SECRETORY or non-renal "
                         "component saturates first")
            elif T > 1 and sign == "falls":
                joint = "net secretion, and secretion saturates first"
            elif T > 1 and sign == "rises":
                joint = "net secretion, but reabsorption saturates first"
            else:
                joint = "indeterminate"
            tstr = f"{T:.3f}" if T else "--"
            print(f"  {chem:8}{sp:7}{tstr:>9}{dom:>14}{sign:>11}  {joint}")
            if conf_up:
                print(f"  {'':8}{'':7}{'':9}{'':14}{'':11}  "
                      f"({conf_up} of {up} rises confounded by possible "
                      f"binding saturation)")
            sexd = KIM_SEX.get(chem)
            sexnote = ((f"T {sexd[0]:.3f}(M) vs {sexd[1]:.3f}(F), "
                        + ("crosses 1" if sexd[0] < 1 < sexd[1]
                           else "no sex difference"))
                       if sexd and sp == "rat" else "")
            verdicts.append(dict(
                chemical=chem, total_carbons=TOTAL_CARBONS[chem], species=sp,
                n_dose_series=len(sub), n_rise=up, n_fall=dn,
                n_rise_confounded=conf_up,
                T_male_rat=(f"{T:.4f}" if T else ""),
                dominant_process=dom, dose_sign=sign,
                kim_T_male=(sexd[0] if sexd and sp == "rat" else ""),
                kim_T_female=(sexd[1] if sexd and sp == "rat" else ""),
                sex_evidence=sexnote, joint_reading=joint))
    print()
    print("  PFHpA POINTS THE SAME WAY TWICE, BUT ACROSS SPECIES. Its RAT")
    print("  tubular factor is 15.0, far above 1, so in the rat it is")
    print("  net-SECRETED -- the only compound here that is. Its MOUSE")
    print("  clearance falls with dose, the signature of saturating")
    print("  secretion. The two agree, but they are different species, and")
    print("  this project does not treat a cross-species agreement as a")
    print("  check. It is suggestive and nothing more; the rat dose series")
    print("  for PFHpA that would make it a check does not exist.")
    print()
    print("  PFOS SEPARATES FROM EVERYTHING ELSE ON THREE INDEPENDENT AXES.")
    print("  Its rat clearance FALLS with dose in four observations across")
    print("  two studies, where every other rat compound rises. It is the")
    print("  only compound with no sex difference in renal clearance. And its")
    print("  tubular factor stays below 1 in both sexes while PFOA's and")
    print("  PFHxS's cross 1. Three different kinds of evidence agree that")
    print("  PFOS is not handled like the others, which is a result in its")
    print("  own right and a warning against treating the class uniformly.")
    print()
    print("  THE ALTERNATIVE IS NOT DISMISSED. PFOS undergoes substantial")
    print("  enterohepatic circulation, so a fall in TOTAL clearance may")
    print("  reflect saturable biliary transport rather than saturable")
    print("  tubular secretion. Separating those needs a dose series in")
    print("  renal clearance specifically, which does not exist for PFOS.")
    print()

    out = []
    p = HERE / "db" / "dose_sign_test.csv"
    cols = ["study", "chemical", "species", "sex", "dose_low", "dose_high",
            "dose_units", "cl_low", "cl_high", "cl_ratio", "direction",
            "cmax_low_mM", "cmax_high_mM", "pfas_to_albumin_at_high_dose",
            "confounder_reading"]
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["chemical"], x["species"],
                                             x["sex"], x["study"])):
            w.writerow(r)
    out.append((p, len(rows)))

    p = HERE / "db" / "tubular_mechanism_verdict.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(verdicts[0].keys()))
        w.writeheader()
        for r in sorted(verdicts, key=lambda x: (x["total_carbons"],
                                                 x["chemical"],
                                                 x["species"])):
            w.writerow(r)
    out.append((p, len(verdicts)))

    for p, n in out:
        print(f"wrote {p.relative_to(HERE)}  ({n} rows)")


if __name__ == "__main__":
    main()
