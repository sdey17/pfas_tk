#!/usr/bin/env python3
"""Abraham 2024's own clearance partition: one subject, no cross-study join.

    python3 -I scripts/stage6_abraham_partition.py

Writes db/abraham2024_clearance_partition.csv. Idempotent.

    Abraham K et al. (2024). Environ Int 193:109047. PMID 39476597,
    doi:10.1016/j.envint.2024.109047

WHY THIS SUPERSEDES THE HUMAN COLUMN OF STAGE 6 PART 1. A single 67-year-old
male volunteer, 82 kg, took a known oral dose of fifteen mostly 13C-labelled
PFAS and was followed for 450 days with plasma, URINE AND FAECES all measured.
The study reports, per compound and in the same subject:

    CLuri    urinary (renal) clearance
    CLtot    total clearance
    CLfec    faecal clearance
    fup      plasma unbound fraction
    CLfilt   glomerular filtration clearance = fup x the subject's OWN
             measured creatinine clearance of 113.7 mL/min

So T = CLuri / CLfilt needs NO cross-study join and no assumed GFR -- every
term is from one person in one experiment. That is unique in this strand, and
it does two things Stage 6 part 1 could not.

1. IT TIGHTENS THE LONG CHAINS INSTEAD OF LOOSENING THEM. Part 1 bounded T by
   CLtot/(fu x GFR), which charges faecal elimination to the kidney, so the
   bound loosened as chain length grew and I reported the apparent rise in T
   at C11-C12 as an artefact that could not be removed. Abraham MEASURES the
   faecal term -- 58 to 61% of total clearance for PFNA through PFDoDA -- so
   the bound becomes (CLtot - CLfec)/CLfilt and the artefact disappears. The
   prediction was right and is now demonstrated rather than asserted.

2. IT PUTS THE SHORT CHAINS IN DOUBT, WHICH PART 1 DID NOT SEE. Abraham's own
   fup for PFHxA is 0.0068 against Ryu's 0.038, a 5.6-fold gap, and the
   filtered load scales with it. On Ryu's fu the human REABSORBS half of
   filtered PFHxA; on Abraham's it SECRETES it at three times the filtered
   load. The sign of T for the short chains is not determined by the data.

Units: the paper reports clearances in mL/min for an 82 kg subject. Converted
here by x 1440 / 82 = x 17.561 to the mL/kg/day used throughout the strand.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
KG, MINS = 82.0, 1440.0
CONV = MINS / KG                      # mL/min -> mL/kg/day
CREAT_CL = 113.7                      # mL/min, the subject's own

# chemical: (CLuri, CLtot, CLfec, fup, CLfilt, uri_flag)  all mL/min
# uri_flag: "measured", "not_quantifiable", "not_detected"
P = {
    "PFBA":    (1.00,  1.24,  None,   0.093,  10.57, "measured"),
    "PFPeA":   (8.89,  8.45,  None,   0.044,   5.00, "measured"),
    "PFHxA":   (2.68,  3.24,  None,   0.0068,  0.77, "measured"),
    "PFHpA":   (0.023, 0.033, None,   0.0004,  0.045, "measured"),
    "PFOA":    (None,  0.0024, 0.0004, 0.001,  0.11, "not_quantifiable"),
    "PFNA":    (None,  0.0038, 0.0022, 0.0016, 0.18, "not_quantifiable"),
    "PFDA":    (None,  0.0084, 0.0051, 0.0027, 0.31, "not_detected"),
    "PFUnDA":  (None,  0.0163, 0.0098, 0.0148, 1.68, "not_detected"),
    "PFDoDA":  (None,  0.0474, 0.0286, None,   None, "not_detected"),
    "PFBS":    (0.067, 0.107, None,   0.0046,  0.52, "measured"),
    "PFHxS":   (None,  0.0030, None,  0.0009,  0.10, "not_detected"),
    "PFOS":    (None,  0.0050, None,  0.0045,  0.51, "not_detected"),
    "6:2 FTS": (0.97,  1.13,  None,   0.0143,  1.63, "measured"),
    "HFPO-DA (GenX)": (1.43, 2.45, None, None, None, "measured"),
    "DONA":    (0.178, 0.191, None,   None,    None, "measured"),
}
NC = {"PFBA": 4, "PFPeA": 5, "PFHxA": 6, "PFHpA": 7, "PFOA": 8, "PFNA": 9,
      "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12, "PFBS": 4, "PFHxS": 6,
      "PFOS": 8, "6:2 FTS": 6, "HFPO-DA (GenX)": 6, "DONA": 7}
ORDER = ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
         "PFUnDA", "PFDoDA", "PFBS", "PFHxS", "PFOS", "6:2 FTS",
         "HFPO-DA (GenX)", "DONA"]


def ryu_human():
    out = {}
    f = HERE / "db" / "plasma_unbound_fraction.csv"
    for r in csv.DictReader(f.open()):
        if r["experiment"] == "Ryu 2024" and r["species"] == "human" \
                and r["fu"]:
            out[r["chemical"]] = float(r["fu"])
    return out


def main():
    ryu = ryu_human()
    bar = "=" * 78
    rows = []
    for c in ORDER:
        uri, tot, fec, fup, filt, flag = P[c]
        # renal clearance, or the tightest bound on it
        if uri is not None:
            cl_r, kind = uri, "measured"
        elif fec is not None:
            cl_r, kind = tot - fec, "bound_total_minus_faecal"
        else:
            cl_r, kind = tot, "bound_total"
        T = cl_r / filt if filt else None
        rows.append(dict(chem=c, nc=NC[c], uri=uri, tot=tot, fec=fec,
                         fup=fup, filt=filt, cl_r=cl_r, kind=kind, T=T,
                         flag=flag))

    print(bar)
    print("A.  ABRAHAM'S OWN T: ONE SUBJECT, ONE EXPERIMENT, NO JOIN")
    print(bar)
    print("  Clearances in mL/min as the paper reports them. CLfilt is the")
    print("  subject's own fup x creatinine clearance of 113.7 mL/min.")
    print()
    print(f"  {'chem':16}{'nC':>3}{'CLuri':>9}{'CLtot':>9}{'CLfec':>9}"
          f"{'fup':>9}{'CLfilt':>8}{'T':>9}  verdict")
    for r in rows:
        if r["T"] is None:
            print(f"  {r['chem']:16}{r['nc']:>3}"
                  f"{(r['uri'] if r['uri'] is not None else '-'):>9}"
                  f"{r['tot']:>9}"
                  f"{(r['fec'] if r['fec'] is not None else '-'):>9}"
                  f"{'-':>9}{'-':>8}{'-':>9}  no fup reported")
            continue
        pre = "" if r["kind"] == "measured" else "<="
        v = ("SECRETION" if r["T"] > 1 and r["kind"] == "measured"
             else "secretion?" if r["T"] > 1
             else f"reabs {100*(1-r['T']):.1f}%")
        print(f"  {r['chem']:16}{r['nc']:>3}"
              f"{(r['uri'] if r['uri'] is not None else 'nd'):>9}"
              f"{r['tot']:>9}"
              f"{(r['fec'] if r['fec'] is not None else '-'):>9}"
              f"{r['fup']:>9}{r['filt']:>8}{pre+format(r['T'],'.4f'):>9}  {v}")

    print()
    print(bar)
    print("B.  THE LONG-CHAIN ARTEFACT IS GONE, AND THE PREDICTION HELD")
    print(bar)
    print("  Part 1 charged all non-renal elimination to the kidney, so its")
    print("  bound loosened with chain length. Abraham measures the faecal")
    print("  term, so it can be subtracted:")
    print()
    print(f"  {'chem':9}{'nC':>3}{'CLfec/CLtot':>13}"
          f"{'part 1 T bound':>16}{'Abraham T bound':>17}")
    import csv as _c
    prev = {}
    f = HERE / "db" / "tubular_factor_evidence.csv"
    for x in _c.DictReader(f.open()):
        if x["species"] == "human" and "Abraham" in x["clearance_source"]:
            prev[x["chemical"]] = float(x["T"])
    for r in rows:
        if r["fec"] is None or r["T"] is None:
            continue
        old = prev.get(r["chem"])
        print(f"  {r['chem']:9}{r['nc']:>3}{100*r['fec']/r['tot']:>12.0f}%"
              f"{('<='+format(old,'.4f')) if old else '-':>16}"
              f"{'<='+format(r['T'],'.4f'):>17}")
    print()
    print("  PFUnDA moves from a bound of 0.266 to one of 0.0039 -- from the")
    print("  loosest long-chain bound in the human series to the TIGHTEST,")
    print("  and PFDA from 0.093 to 0.0106. The apparent rise in T toward C12")
    print("  was the faecal term, as part 1 predicted but could not show.")
    print("  There is NO LONG-CHAIN LOSS OF REABSORPTION in the human: C8")
    print("  through C11 all exceed 98% reabsorbed on this subject's own")
    print("  numbers, and the trend with chain length is flat or slightly")
    print("  tightening rather than rising.")
    print()
    print("  PFOA is the exception that proves the mechanism: its bound moves")
    print("  the OTHER way, from 0.0066 to 0.0182, because only 17% of its")
    print("  clearance is faecal so there is little to subtract -- while its")
    print("  fup differs 2.4-fold between the two sources, which dominates.")
    print("  The verdict is unchanged either way, both being far below 1.")

    print()
    print(bar)
    print("C.  BUT THE SHORT CHAINS NOW HAVE NO DETERMINED SIGN")
    print(bar)
    print("  Three laboratories measured human plasma fu. They disagree most")
    print("  at the short chains, and the filtered load scales with fu, so the")
    print("  sign of T goes with the choice:")
    print()
    print(f"  {'chem':9}{'Abraham fup':>13}{'Ryu fu':>10}{'ratio':>8}"
          f"{'T(Abraham)':>12}{'T(Ryu)':>9}   sign")
    GFR_H = 2570.0
    for c in ("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFBS", "6:2 FTS"):
        r = [x for x in rows if x["chem"] == c][0]
        if r["fup"] is None or c not in ryu:
            continue
        cl_kg = r["cl_r"] * CONV
        T_ab = r["T"]
        T_ryu = cl_kg / (ryu[c] * GFR_H)
        sign = ("agree" if (T_ab > 1) == (T_ryu > 1) else "DISAGREE")
        print(f"  {c:9}{r['fup']:>13.4f}{ryu[c]:>10.4f}"
              f"{ryu[c]/r['fup']:>8.1f}{T_ab:>12.3f}{T_ryu:>9.3f}   {sign}")
    print()
    print("  PFHxA: 0.48 on Ryu's fu, 3.48 on Abraham's own. Reabsorption or")
    print("  secretion, decided by a 5.6-fold disagreement in a binding")
    print("  measurement. PFPeA is the same story.")
    print()
    print("  AND THE CLEARANCES DO NOT DISAGREE AT ALL. Abraham's urinary")
    print(f"  clearance of PFHxA is {2.68*CONV:.1f} mL/kg/day; Argoul 2026")
    print("  reports 47.0 from an unrelated population study. They agree to")
    print("  better than 1%. SO THE UNCERTAINTY IS ENTIRELY IN THE BINDING,")
    print("  NOT IN THE CLEARANCE -- which is the single most useful thing")
    print("  this comparison says about where effort should go.")

    print()
    print(bar)
    print("D.  WHICH fu SERIES TO BELIEVE: THE SECTION 16 TEST, RERUN")
    print(bar)
    print("  Binding rises with chain length, so fu must FALL monotonically")
    print("  across the carboxylate series. Section 16 used this to prefer Ryu")
    print("  over Argoul at the long chains. Abraham is a third series:")
    print()
    for name, d in (("Abraham fup", {r["chem"]: r["fup"] for r in rows
                                     if r["fup"]}),
                    ("Ryu fu", ryu)):
        ser = [(NC[c], d[c]) for c in ("PFBA", "PFPeA", "PFHxA", "PFHpA",
                                       "PFOA", "PFNA", "PFDA", "PFUnDA",
                                       "PFDoDA") if c in d]
        inv = [(ser[i][0], ser[i+1][0]) for i in range(len(ser)-1)
               if ser[i+1][1] > ser[i][1]]
        where = ", ".join(f"C{a}->C{b}" for a, b in inv) or "none"
        print(f"  {name:14}{len(ser)-1} pairs, {len(inv)} inversions: {where}")
    print()
    print("  ABRAHAM'S SERIES INVERTS FOUR TIMES, RYU'S NOT AT ALL. Abraham's")
    print("  fup falls to 0.0004 at C7 and then rises 37-fold to 0.0148 at")
    print("  C11, which no binding mechanism supports. So the monotonicity")
    print("  test, applied impartially, prefers Ryu's BINDING -- while")
    print("  Abraham remains the only source whose clearance partition needs")
    print("  no cross-study join.")
    print()
    print("  THE HONEST POSITION IS THEREFORE SPLIT, AND IT IS NOT A FUDGE:")
    print("    use Abraham's measured FAECAL term, which no one else has, to")
    print("      tighten the long chains -- it is a mass balance, not a")
    print("      binding assay, and nothing above questions it;")
    print("    use Ryu's fu for the filtered load, because its series passes")
    print("      the monotonicity test and Abraham's does not;")
    print("    and report the short chains as sign-undetermined, because on")
    print("      that combination PFHxA sits at T = 0.48 while Abraham's own")
    print("      numbers give 3.48, and no test available here separates")
    print("      them.")

    p = HERE / "db" / "abraham2024_clearance_partition.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "cluri_mL_min", "cluri_status",
                    "cltot_mL_min", "clfec_mL_min", "faecal_pct_of_total",
                    "fup", "clfilt_mL_min", "creatinine_clearance_mL_min",
                    "cl_renal_used_mL_min", "cl_renal_basis",
                    "cl_renal_mL_kg_day", "T_abraham_own", "source"])
        for r in rows:
            w.writerow([r["chem"], r["nc"],
                        r["uri"] if r["uri"] is not None else "", r["flag"],
                        r["tot"], r["fec"] if r["fec"] is not None else "",
                        round(100*r["fec"]/r["tot"], 1) if r["fec"] else "",
                        r["fup"] if r["fup"] is not None else "",
                        r["filt"] if r["filt"] is not None else "",
                        CREAT_CL, round(r["cl_r"], 6), r["kind"],
                        round(r["cl_r"] * CONV, 4),
                        round(r["T"], 5) if r["T"] is not None else "",
                        "Abraham 2024 Tables 1/3/4, Environ Int 193:109047, "
                        "doi:10.1016/j.envint.2024.109047; n=1, 67 y male, "
                        "82 kg"])
    print(f"\nwrote {p.relative_to(HERE)}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
