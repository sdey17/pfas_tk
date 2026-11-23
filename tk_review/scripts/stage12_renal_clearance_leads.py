#!/usr/bin/env python3
"""Literature leads for DIRECT human renal clearance (Stage 12).

    python3 -I scripts/stage12_renal_clearance_leads.py

Writes db/renal_clearance_leads.csv. Idempotent. NOTHING HERE ENTERS THE
T PIPELINE -- these are candidates awaiting primary verification, and the
file records exactly what would be needed to admit each one.

WHY THIS FILE EXISTS. The human DIRECT tier of db/tubular_factor_
consolidated.csv currently holds TWO compounds (PFBS and PFHxA, both from
Argoul). Everything else human is an UPPER bound derived from total
clearance. A search of the literature shows that is not because the
measurements do not exist -- it is because they were never ingested.

WHAT WAS FOUND. EFSA's 2020 Opinion (doi 10.2903/j.efsa.2020.6223) tabulates
human renal clearances for TEN compounds, citing four primary studies:

    Zhang et al. 2013   PMID 23980546, doi 10.1021/es401905e
                        86 paired blood/urine samples, adults in China,
                        isomer-specific. The widest-coverage human source.
    Fujii et al. 2015   renal, biliary and faecal clearances in humans
    Harada et al. 2005  PMID 16194675, doi 10.1016/j.envres.2004.12.003
                        PFOA and PFOS, young/old x M/F, N=5 per cell
    Harada et al. 2007  follow-up

WHY THEY ARE NOT INGESTED HERE. Two reasons, both disqualifying on their own.

  (1) THE SECONDARY SOURCE CONTRADICTS ITSELF. Two passages of the SAME EFSA
      document give different values for the same compounds:
          PFHpA    0.61  in the per-compound tabulation
                   0.67  in the chain-length summary sentence
          PFUnDA   0.064 in the per-compound tabulation
                   0.005 in the chain-length summary sentence   <- 13x
      A 13-fold disagreement inside one document is not a rounding
      difference. Until the primary tables are read, neither value is
      usable, and a project rule that prefers primary sources exists
      precisely for this case.

  (2) SOME VALUES EXCEED INDEPENDENTLY REPORTED TOTAL CLEARANCES. Renal
      clearance cannot be larger than total clearance in the same subject.
      EFSA's PFOA renal range (0.13-0.77 mL/kg/day) is 3- to 18-fold ABOVE
      Abraham 2024's measured TOTAL clearance for PFOA (0.0417) and Argoul's
      (0.042). Different cohorts, so not a contradiction in the strict
      sense -- but a difference that large between the renal part and other
      people's whole must be understood before either is used.

THE PROVISIONAL T BELOW IS A DECISION AID, NOT A RESULT. It asks: IF the
EFSA values survive primary checking, would they change any conclusion? That
is what decides whether chasing three PDFs is worth anyone's afternoon.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
GFR_HUMAN = 2570.0

# compound: (CL_renal low, high, primary study as EFSA cites it)
EFSA = {
    "PFHpA":  (0.61, 0.61,   "Zhang 2013; Fujii 2015"),
    "PFOA":   (0.13, 0.77,   "Harada 2005; Harada 2007; Zhang 2013; Fujii 2015"),
    "PFNA":   (0.062, 0.25,  "Zhang 2013; Fujii 2015"),
    "PFDA":   (0.066, 0.096, "Zhang 2013; Fujii 2015"),
    "PFUnDA": (0.064, 0.065, "Zhang 2013; Fujii 2015"),
    "PFDoDA": (0.070, 0.070, "Fujii 2015"),
    "PFTriDA": (0.077, 0.077, "Fujii 2015"),
    "PFTeDA": (0.224, 0.224, "Fujii 2015"),
    "PFHxS":  (0.027, 0.039, "Zhang 2013"),
    "PFOS":   (0.037, 0.106, "Harada 2005; Harada 2007; Zhang 2013"),
}
CONFLICT = {"PFHpA": "summary sentence of the same document says 0.67",
            "PFUnDA": "summary sentence of the same document says 0.005 (13x)"}
NC = {"PFHpA": 7, "PFOA": 8, "PFNA": 9, "PFDA": 10, "PFUnDA": 11,
      "PFDoDA": 12, "PFTriDA": 13, "PFTeDA": 14, "PFHxS": 6, "PFOS": 8}


def human_fu():
    out = {}
    for r in csv.DictReader(
            (HERE / "db" / "plasma_unbound_fraction.csv").open()):
        if r["experiment"] == "Ryu 2024" and r["species"] == "human" and r["fu"]:
            out[r["chemical"]] = float(r["fu"])
    return out


def existing_human():
    """What the pipeline already has for the human, by compound and tier."""
    out = {}
    f = HERE / "db" / "tubular_factor_consolidated.csv"
    for r in csv.DictReader(f.open()):
        if r["species"] == "human":
            out.setdefault(r["chemical"], set()).add(r["tier"])
    return out


def main():
    fu, have = human_fu(), existing_human()
    bar = "=" * 78
    print(bar)
    print("CANDIDATE HUMAN RENAL CLEARANCES FOUND IN THE LITERATURE")
    print(bar)
    print("  Source: EFSA 2020 Opinion, which tabulates them from four")
    print("  primary studies. NOT ingested -- see the header for why.")
    print()
    print(f"  {'chem':9}{'nC':>3}{'CL_R low-high':>16}{'fu':>10}"
          f"{'filtered':>10}{'T low-high':>18}  now in pipeline as")
    rows = []
    for chem, (lo, hi, src) in sorted(EFSA.items(), key=lambda kv: NC[kv[0]]):
        f = fu.get(chem)
        tiers = have.get(chem, set())
        cur = "+".join(sorted(tiers)) if tiers else "absent"
        if f:
            filt = f * GFR_HUMAN
            tlo, thi = lo / filt, hi / filt
            tt = f"{tlo:.4f} - {thi:.4f}"
            verdict = ("net_reabsorption" if thi < 1 else
                       "net_secretion" if tlo > 1 else "straddles_1")
        else:
            filt = tlo = thi = None
            tt, verdict = "no human fu", "undetermined"
        print(f"  {chem:9}{NC[chem]:>3}{f'{lo:g} - {hi:g}':>16}"
              f"{(f'{f:g}' if f else '-'):>10}"
              f"{(f'{filt:.2f}' if filt else '-'):>10}{tt:>18}  {cur}")
        rows.append([chem, NC[chem], lo, hi, "mL/kg/day", f or "",
                     round(filt, 3) if filt else "",
                     round(tlo, 5) if tlo else "",
                     round(thi, 5) if thi else "", verdict, cur,
                     "EFSA 2020 Opinion doi 10.2903/j.efsa.2020.6223 "
                     "(secondary tabulation)", src,
                     CONFLICT.get(chem, ""), "primary PDF not yet read"])
    print()
    print(bar)
    print("WOULD THIS CHANGE ANY CONCLUSION?")
    print(bar)
    reab = [r for r in rows if r[9] == "net_reabsorption"]
    new = [r for r in rows if r[10] == "absent"]
    upg = [r for r in rows if "UPPER" in r[10] and "DIRECT" not in r[10]]
    print(f"  {len(rows)} candidate compounds, {len(reab)} of which give")
    print(f"  net reabsorption on the project's own human fu.")
    print(f"  {len(new)} are compounds the human column does not have at all.")
    print(f"  {len(upg)} would be UPGRADED from an upper bound to a direct")
    print("  measurement -- which is the scarce thing, since the human DIRECT")
    print("  tier currently holds exactly two compounds.")
    print()
    if len(reab) == len(rows):
        print("  EVERY candidate points the same way the project already")
        print("  concludes: net reabsorption across the human series. So the")
        print("  DIRECTION is not at stake. What these would buy is the")
        print("  EVIDENCE TIER and the MAGNITUDE -- turning bounds into")
        print("  measurements, which is exactly what the chain-length")
        print("  analysis is short of.")
    print()
    print("  WORTH NOTING FOR THE SHORT CHAINS: none of these candidates")
    print("  covers PFPeA or PFHxA, the two compounds whose SIGN is disputed")
    print("  (section 7b). So this literature does NOT settle the open")
    print("  question; it strengthens the part that was never in doubt.")
    print()
    p = HERE / "db" / "renal_clearance_leads.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_carbons", "cl_renal_low", "cl_renal_high",
                    "units", "fu_human_Ryu2024", "filtered_load_mL_kg_day",
                    "T_provisional_low", "T_provisional_high",
                    "provisional_verdict", "current_pipeline_status",
                    "secondary_source", "primary_studies_cited",
                    "internal_conflict", "admission_blocker"])
        w.writerows(rows)
    print(f"wrote {p.relative_to(HERE)}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
