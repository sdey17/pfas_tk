#!/usr/bin/env python3
"""Stage 4: the species question, as an accounting rather than a model.

    python3 -I scripts/stage4_species_accounting.py

Read-only. Aim 1 asks why PFOA persists ~50-fold longer in humans than in
rats. Sections 8, 10, 11 and 13 each supply a piece. This stage puts them
against the observed difference and asks how much of it they account for.

IT IS AN ACCOUNTING, NOT A MODEL, AND THE DISTINCTION IS THE POINT. Turning
the measured ratios below into a predicted renal clearance needs membrane
protein per gram of kidney, which Stage 3 declined to invent. So this stage
decomposes the OBSERVED difference exactly, where the arithmetic is OEHHA's
and solid, then lists each candidate mechanism with the size of the species
difference it carries, and says which are measured and which are assumed. It
does not multiply them together to produce a prediction, because two of the
factors are not independent and one of them cannot be put on a common scale.

THE SEX DIFFERENCE IS THE VALIDATION CASE. Within the rat, one transporter
differs 114-fold between sexes and the direction of net tubular handling
flips. If the framework explains that, it has earned the right to be applied
across species. If it did not, nothing else here would be worth saying.
"""
import sys

# ---- OEHHA PHG Table A6.4, fu assumed 0.02 for every species -------------
FU = 0.02
CASES = {   # label: (GFR mL/kg-day, CLR mL/kg-day)
    "human": (2570, 0.06),
    "rat male (Kemper)": (14400, 18.2),
    "rat male (Ohmori)": (14400, 21.2),
    "rat female (Kemper)": (14400, 666),
    "rat female (Ohmori)": (14400, 1009),
    "mouse male": (16700, 10),
    "mouse female": (16700, 16),
}

# ---- measured species ratios, rat / human --------------------------------
#  factor, rat, human, units, what it is, how certain
FACTORS = [
    ("GFR", 14400, 2570, "mL/kg-day",
     "glomerular filtration rate", "MEASURED (OEHHA Table A6.4)"),
    ("OAT1 abundance", 32.8, 4.8, "pmol/mg membrane protein",
     "basolateral uptake, step 1 of secretion",
     "MEASURED in one assay (Thakur 2024)"),
    ("OAT3 abundance", 6.3, 4.3, "pmol/mg membrane protein",
     "basolateral uptake, step 1 of secretion",
     "MEASURED in one assay (Thakur 2024)"),
    ("OAT1 CL_int for PFOA", 8.26, 9.99, "uL/min/mg protein",
     "per-molecule activity of the same transporter",
     "MEASURED in one system (Nakagawa 2008)"),
    ("OAT3 CL_int for PFOA", 4.56, 8.13, "uL/min/mg protein",
     "per-molecule activity of the same transporter",
     "MEASURED in one system (Nakagawa 2008)"),
]


def main():
    print("=" * 76)
    print("1. WHAT HAS TO BE EXPLAINED")
    print("=" * 76)
    print(f"{'case':22s} {'filtered':>9s} {'CL_R':>8s} {'reabsorbed':>11s} "
          f"{'escaping':>9s}")
    esc = {}
    for k, (gfr, clr) in CASES.items():
        f = FU * gfr
        e = clr / f
        esc[k] = e
        tag = f"{100 * (1 - e):10.2f}%" if e < 1 else "   SECRETION"
        print(f"{k:22s} {f:9.1f} {clr:8.2f} {tag} {100 * e:8.2f}%")

    h, r = esc["human"], esc["rat male (Kemper)"]
    gfr_ratio = CASES["rat male (Kemper)"][0] / CASES["human"][0]
    clr_ratio = CASES["rat male (Kemper)"][1] / CASES["human"][1]
    print()
    print(f"  Male rat clears PFOA renally {clr_ratio:.0f}x faster than human,")
    print(f"  and that number factors EXACTLY in two:")
    print(f"     {clr_ratio:.0f}x  =  {gfr_ratio:.1f}x (filtered load, i.e. GFR "
          f"at equal fu)")
    print(f"           x  {r / h:.0f}x (the fraction that escapes the tubule)")
    print()
    print("  So roughly a sixth of the species difference is glomerular and")
    print("  needs no transporter at all. The rest is tubular, and that is")
    print("  what sections 10, 11 and 13 are about.")

    print("\n" + "=" * 76)
    print("2. THE VALIDATION CASE: the rat sex difference")
    print("=" * 76)
    for a, b in (("rat male (Kemper)", "rat female (Kemper)"),
                 ("rat male (Ohmori)", "rat female (Ohmori)")):
        print(f"  {a.split(' (')[1][:-1]:8s}  escape {100*esc[a]:6.2f}% (male) "
              f"-> {100*esc[b]:7.2f}% (female),  "
              f"{esc[b] / esc[a]:5.1f}x, and the sign flips")
    print()
    print("  Same species, same GFR, same assumed fu. What differs:")
    print("     Oatp1a1 (apical, REABSORPTIVE)   114x higher in males")
    print("     Oatp1a3 (apical, reabsorptive)     7x higher in males")
    print("     Oat1    (basolateral, SECRETORY)   2x higher in males")
    print()
    print("  Note the direction of the Oat1 difference: males have MORE")
    print("  secretory transporter and still clear PFOA 37-48x more slowly.")
    print("  The reabsorptive side is not merely contributing, it is")
    print("  overriding a secretory difference pointing the other way.")
    print()
    print("  A 114-fold change in one apical transporter turns net")
    print("  reabsorption of 94% into net secretion. That is the framework")
    print("  working on a case where the mechanism is unambiguous, and it is")
    print("  the reason to take the species comparison below seriously.")

    print("\n" + "=" * 76)
    print("3. THE SPECIES FACTORS, MEASURED")
    print("=" * 76)
    print(f"{'factor':24s} {'rat':>9s} {'human':>9s} {'rat/human':>10s}  basis")
    for name, rv, hv, units, what, cert in FACTORS:
        print(f"{name:24s} {rv:9.1f} {hv:9.1f} {rv / hv:10.2f}  {cert}")
    print()
    sec = (32.8 / 4.8) * (8.26 / 9.99)
    print(f"  Combining abundance with per-molecule activity, rat OAT1")
    print(f"  secretory capacity is about {sec:.1f}x the human. OAT3 is")
    print(f"  {(6.3 / 4.3) * (4.56 / 8.13):.2f}x -- rat has more protein but "
          f"lower activity, and they")
    print("  nearly cancel.")

    print("\n" + "=" * 76)
    print("4. THE REABSORPTIVE SIDE CANNOT BE PUT ON A RATIO AT ALL")
    print("=" * 76)
    print("  human OAT4     0.5 pmol/mg      rat: NOT DETECTED")
    print("  rat Oatp1a1    1.9 pmol/mg      human: NOT DETECTED")
    print()
    print("  These are different genes. There is no rat/human ratio to form,")
    print("  and no REF for either, so neither can be scaled. What can be")
    print("  compared is intrinsic clearance per mg in each species' own")
    print("  system, and they are strikingly similar for PFOA:")
    print("     human OAT4     95.7 uL/min/mg  (Louisse 2023)")
    print("     rat Oatp1a1    73.6 uL/min/mg  (Weaver 2010)")
    print("  a 1.3-fold difference, in the direction of the human reabsorbing")
    print("  slightly faster per mg of transporter protein.")

    print("\n" + "=" * 76)
    print("5. WHAT THIS DOES AND DOES NOT SETTLE")
    print("=" * 76)
    print("  SETTLED, and quantitative:")
    print(f"     - {gfr_ratio:.1f}x of the {clr_ratio:.0f}x is GFR, measured, "
          f"and needs no transporter.")
    print("     - Transporter AFFINITY is not the explanation: Nakagawa's")
    print("       matched human/rat OAT1 and OAT3 differ 1.2-1.8x (section 8).")
    print("     - Substrate SPECIFICITY is not the explanation: the same")
    print("       ~500 g/mol threshold holds in both species (Argoul 2026).")
    print("     - The reabsorptive machinery is NOT THE SAME GENES, measured")
    print("       by proteomics (section 11).")
    print("     - A reabsorptive transporter CAN dominate net handling, shown")
    print("       by the rat sex difference above.")
    print()
    print("  NOT SETTLED, and honestly so:")
    print("     - Whether human OAT4 reabsorbs MORE EFFECTIVELY in situ than")
    print("       rat Oatp1a1. Per mg they are within 1.3-fold, but neither")
    print("       has a REF, so per-kidney they cannot be compared.")
    print("     - The fu assumption. OEHHA assumes 0.02 for every species.")
    print("       Argoul 2026 measures human PFOA fu at 0.0017, twelvefold")
    print("       lower; Starnes 2024 (DSF) finds rat albumin binds PFOA")
    print("       about 2x TIGHTER than human, which would put rat fu lower")
    print("       still. If rat fu is lower than human fu, the filtered-load")
    print("       term above is overstated for the rat and MORE of the")
    print("       species difference has to come from the tubule, not less.")
    print("     - No measured rat plasma fu for PFOA is in this database.")

    print("\n" + "=" * 76)
    print("6. CAN THE PROJECT'S OWN BINDING DATA FIX THE fu ASSUMPTION? NO.")
    print("=" * 76)
    print("  The albumin sheet is human and bovine by the reviewer's own")
    print("  decision, and Fan 2025 and Han 2003 were set aside as rat-only.")
    print("  For the albumin sheet that was right. For this stage it removes")
    print("  the one species comparison that would matter, so the rat rows")
    print("  were re-examined here. They do not help, and the reason is")
    print("  worth recording.")
    print()
    print("  Fan 2025 measured the SAME rat serum albumin by TWO methods:")
    import csv, collections, os, statistics
    here = os.path.dirname(os.path.abspath(__file__))
    db = os.path.join(here, "..", "db", "protein_binding.csv")
    rsa = collections.defaultdict(dict)
    with open(db, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if (r["study"].startswith("Fan")
                    and "albumin" in (r["protein"] or "").lower()
                    and r["kd"].strip()):
                m = ("UF" if "ultrafiltration" in r["method"]
                     else "ITC" if "calorimetry" in r["method"] else "?")
                rsa[r["chemical"]][m] = float(r["kd"])
    folds = []
    for c, v in sorted(rsa.items()):
        if "UF" in v and "ITC" in v:
            f = v["ITC"] / v["UF"]
            folds.append(f)
            print(f"     {c:9s} ultrafiltration {v['UF']:7.2f} uM   "
                  f"ITC {v['ITC']:7.2f} uM   {f:5.1f}x apart")
    if folds:
        print(f"\n  Median disagreement {statistics.median(folds):.0f}-fold, "
              f"range {min(folds):.1f}-{max(folds):.1f}-fold, within one")
        print("  paper and one protein preparation. The rat-versus-human")
        print("  difference this stage would need to detect is about 2-fold.")
        print("  The binding data is four times noisier than the signal, so")
        print("  it cannot settle the question either way.")
        print()
        print("  Starnes 2024 points the same direction but by a method this")
        print("  project treats cautiously (DSF): rat albumin about 2x")
        print("  tighter than human for PFOA. Suggestive, not usable.")
        print()
        print("  WHAT WOULD SETTLE IT is not another albumin Kd. It is a")
        print("  measured unbound fraction in rat PLASMA, by the method")
        print("  Argoul 2026 used for human and mouse. That single number")
        print("  would turn the filtered-load term in section 1 from an")
        print("  assumption into a measurement, for both species.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
