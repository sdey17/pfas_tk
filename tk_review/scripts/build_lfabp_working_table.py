#!/usr/bin/env python3
"""Build the L-FABP partition-coefficient sheet, on the same footing as albumin.

Produces lfabp_logK_Kd.xlsx with the same eleven columns as
albumin_logK_Kd_FINAL.xlsx, so the two can sit side by side.

THE BRIDGE IS NOT ALBUMIN'S. L-FABP is about a fifth of albumin's mass, so

    log K = log10(1e9 / M) + log10(n) - log10(Kd / uM)

has a different leading constant for every species:

    human L-FABP  M = 14208 Da (UniProt P07148)  ->  4.8475
    rat   L-FABP  M = 14273 Da (UniProt P02692)  ->  4.8455
    mouse L-FABP  M = 14246 Da (UniProt P12710)  ->  4.8463
    (albumin, for contrast)     66463 Da         ->  4.1774

Running an L-FABP Kd through albumin's 4.1774 understates log K by 0.67, a
factor of 4.7. The per-species masses come from add_lfabp_reported_conc.py,
which established them.

BOTH SPECIES ARE CARRIED, which is the opposite of the albumin sheet's choice.
There, human and bovine were kept separable because the species question was
confounded with method and coverage. Here rat is the only non-human arm and it
comes from a single laboratory, so the sheet carries a Species column and the
audit reports the comparison, but no row is pooled across species.

GATES, mirroring the albumin sheet:

  holdout        Han 2021. Its hL-FABP constants are quoted from other papers,
                 not measured -- the row notes say so, and its PFBS value is
                 Zhang 2013's. Its own measurements are thermal-shift dTm,
                 which is not a dissociation constant.
  computational  Zhao 2023's four rows are docking scores and QSAR summaries
                 over a 1599-compound library, with no measured constant.
  not a PFAS     the fluorescent probes these assays displace -- DAUDA,
                 1,8-ANS -- and oleic acid, L-FABP's natural ligand. These are
                 the reference ligands, not analytes.
  several        a cell naming several compounds is a summary. Tested the
  compounds      data-driven way, as in the albumin builder: split on commas
                 outside brackets and count how many parts appear elsewhere as
                 standalone compounds, because an IUPAC name's commas are
                 locants. "2H,2H,3H,3H-perfluorooctanoic acid" is one compound;
                 "PFOS, N-EtFOSA, Wyeth-14643, N-EtFOSE, PFOA" is five.
  composite Kd   Woodcroft 2010 reports two and three site-specific Kd values
                 in one cell ("3.1e-6, 26.2e-6, 52.6e-6"). Those are real
                 counted stoichiometries -- n = 2 by displacement, n = 3 by
                 ITC, among the few measured n in this literature -- but a
                 multi-site series cannot enter a one-row-per-measurement
                 sheet without a decision about which site. Held out pending
                 a read of the primary.
  salt duplicate Yang 2020 lists PFOA and PFOA ammonium salt, PFOS and PFOS
                 potassium salt, at IDENTICAL Kd: one measurement under two
                 names. The free-acid name is kept and the salt row dropped.
                 PFBS and PFBS potassium salt differ (116.4 vs 125.25) and are
                 kept as two.

ROUTE LABELS follow the albumin sheet's vocabulary, assigned from what each
paper actually reports:

  n = 1 by single-site model   the paper's own model fixes n = 1
                               (Fan 2025, Zhang 2013, Khazaee 2021)
  n fitted by binding isotherm the paper fits n from its own displacement
                               isotherm (Sheng 2016). NOTE 7 of its 22 fitted
                               values are BELOW 1, which a site count cannot
                               be, so these are fitted stoichiometries and not
                               counted sites; the label says "fitted" for that
                               reason.
  n = 1 assumed                the paper reports no n (Yang 2020, Gao 2019,
                               Jia 2022)

Run:  python3 -I scripts/build_lfabp_working_table.py
"""
import collections, csv, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "lfabp_logK_Kd.xlsx")

MW = {"human": 14208.0, "rat": 14273.0, "mouse": 14246.0}
HOLDOUT = {"Han 2021": "constants quoted from other papers, not measured here"}
COMPUTATIONAL = {"Zhao 2023"}
# Reference ligands and natural ligands, not analytes. These assays work by
# displacing a fluorescent probe (DAUDA, 1,8-ANS, NBD-lauric acid) and are
# calibrated against L-FABP's own fatty-acid ligands (oleic, palmitic), so all
# of these appear in the tables as controls. Palmitic acid in particular reads
# as a strong binder (log K 5.02 from Zhang 2013) and would be the third
# highest value in the sheet, inside a PFAS chain-length relationship.
NOT_ANALYTE = re.compile(
    r"^(DAUDA|1,8-ANS|ANS\b|NBD|oleic acid|palmitic acid|stearic acid|"
    r"lauric acid|myristic acid|arachidonic acid)", re.I)
PAPER_SINGLE_SITE = {"Fan 2025", "Zhang 2013", "Khazaee 2021"}

# Rows whose Kd lies outside the concentration series the paper actually
# tested, keyed on (source, compound). A fitted constant beyond the top of the
# ligand range is an extrapolation of the model, not a measurement.
#
#   Khazaee 2021 PFHxS, Kd 1.7 uM. Its protein was ~77 nM (derived in
#   fix_khazaee2021_protein_conc.py) and its ligand series ran 3.8-384 nM, so
#   1.7 uM is 22x the protein concentration and 4.4x above the highest ligand
#   concentration used. Its PFOA (0.099 uM) and PFOS (0.18 uM) sit INSIDE that
#   series and are kept.
OUTSIDE_RANGE = {("Khazaee 2021", "pfhxs")}
FITTED = {"Sheng 2016"}

COLS = ["PFAS", "Species", "Kd (uM)", "n used",
        "log K (L water / kg protein)", "Route", "Source",
        "PFAS conc (uM)", "Protein conc (uM)", "ligand:protein ratio",
        "f_bound (%)"]
TO_UM = {"uM": 1.0, "nM": 1e-3, "mM": 1e3, "M": 1e6}

def num(s):
    m = re.match(r"\s*([0-9.]+(?:[eE][-+]?\d+)?)", str(s or ""))
    return float(m.group(1)) if m else None

def source_of(r):
    st = (r.get("study") or "").replace("et al.", "").strip().rstrip(",")
    return f"{st} {r.get('year') or ''}".strip()

def names_several(name, known):
    txt = str(name or "")
    if "," not in re.sub(r"\([^)]*\)", "", txt):
        return False
    parts, depth, cur = [], 0, ""
    for ch in txt:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            parts.append(cur); cur = ""
        else:
            cur += ch
    parts.append(cur)
    return sum(1 for q in parts
               if re.sub(r"\s*\([^)]*\)", "", q).strip().lower() in known) >= 2

def salt_base(name):
    return re.sub(r",?\s*(potassium|ammonium|sodium)\s*salt", "",
                  str(name or ""), flags=re.I).strip()

def main():
    raw = list(csv.DictReader(open(DB, encoding="utf-8")))
    lf = [r for r in raw
          if "liver fatty acid" in (r.get("protein") or "").lower()
          and "albumin" not in (r.get("protein") or "").lower()]
    known = {str(r["chemical"]).strip().lower() for r in lf
             if "," not in str(r["chemical"])}

    kept, dropped = [], collections.Counter()
    seen_salt = {}
    for r in lf:
        src, sp, ch = source_of(r), r.get("protein_species"), r["chemical"]
        if src in HOLDOUT:
            dropped["holdout: " + src] += 1; continue
        if src in COMPUTATIONAL:
            dropped["computational: " + src] += 1; continue
        if sp not in MW:
            dropped[f"species not carried ({sp})"] += 1; continue
        if NOT_ANALYTE.match(ch.strip()):
            dropped["not an analyte (probe or natural ligand)"] += 1; continue
        if names_several(ch, known):
            dropped["several compounds in one cell"] += 1; continue
        if "," in str(r.get("kd") or ""):
            dropped["composite Kd (multi-site in one cell)"] += 1; continue
        if (src, ch.strip().lower()) in OUTSIDE_RANGE:
            dropped["Kd outside the concentration series tested"] += 1; continue

        kd, ku = r.get("kd"), (r.get("kd_units") or "").strip()
        ka = r.get("ka")
        if kd and num(kd) is not None:
            f = TO_UM.get(ku.split()[0] if ku else "")
            if f is None:
                dropped[f"Kd units not understood ({ku[:16]!r})"] += 1; continue
            kd_uM = num(kd) * f
        elif ka and num(ka) is not None:
            kd_uM = 1e6 / num(ka)
        else:
            dropped["no numeric constant"] += 1; continue

        # one measurement under two names: free acid and its salt
        base = salt_base(ch)
        key = (src, sp, base, round(kd_uM, 6))
        if key in seen_salt:
            dropped["salt form of a row already kept"] += 1; continue
        seen_salt[key] = True

        n_raw = num(r.get("n_binding_sites"))
        if src in FITTED and n_raw:
            n, route = n_raw, "n fitted by binding isotherm"
        elif src in PAPER_SINGLE_SITE:
            n, route = 1.0, "n = 1 by single-site model"
        elif n_raw:
            n, route = n_raw, "n fitted by binding isotherm"
        else:
            n, route = 1.0, "n = 1 assumed"

        bridge = math.log10(1e9 / MW[sp])
        kept.append({
            "PFAS": base, "Species": sp, "Kd (uM)": round(kd_uM, 4),
            "n used": n,
            "log K (L water / kg protein)":
                round(bridge + math.log10(n) - math.log10(kd_uM), 3),
            "Route": route, "Source": src,
            "PFAS conc (uM)": num(r.get("ligand_conc_range")) or "",
            "Protein conc (uM)": num(r.get("protein_conc")) or "",
            "ligand:protein ratio": num(r.get("ligand_over_protein_max")) or "",
            "f_bound (%)": num(r.get("unbound_fraction")) or ""})

    kept.sort(key=lambda d: (d["Species"], d["Source"], d["PFAS"]))
    print(f"L-FABP sheet: {len(kept)} rows from "
          f"{len({d['Source'] for d in kept})} sources, "
          f"{len({d['PFAS'] for d in kept})} compounds")
    for sp in ("human", "rat", "mouse"):
        v = [d["log K (L water / kg protein)"] for d in kept if d["Species"] == sp]
        if v:
            v.sort()
            print(f"   {sp:6s} n={len(v):3d}  log K {v[0]:.2f} to {v[-1]:.2f}  "
                  f"median {v[len(v)//2]:.2f}")
    print("\nroutes:")
    for rt, n in collections.Counter(d["Route"] for d in kept).most_common():
        print(f"   {n:4d}  {rt}")
    print("\ndropped:")
    for k, n in sorted(dropped.items(), key=lambda x: -x[1]):
        print(f"   {n:4d}  {k}")

    import openpyxl
    from openpyxl.styles import Font
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Coefficients"
    ws.append(COLS)
    for c in ws[1]:
        c.font = Font(bold=True)
    for d in kept:
        ws.append([d[c] for c in COLS])
    ws.freeze_panes = "A2"
    for i, c in enumerate(COLS, 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = \
            max(11, min(34, len(c) + 3))
    nb = wb.create_sheet("Notes")
    nb.append(["L-FABP partition coefficients"])
    nb.append([])
    nb.append(["bridge", "log K = log10(1e9/M) + log10(n) - log10(Kd/uM)"])
    for sp, m in sorted(MW.items()):
        nb.append([f"  {sp} L-FABP", f"M = {m:.0f} Da -> constant {math.log10(1e9/m):.4f}"])
    nb.append(["  albumin, for contrast", "M = 66463 Da -> constant 4.1774"])
    nb.append([])
    nb.append(["EXCLUDED"])
    for k, n in sorted(dropped.items(), key=lambda x: -x[1]):
        nb.append([f"  {k}", n])
    nb.append([])
    nb.append(["Species are NOT pooled. Rat comes from a single laboratory "
               "(Fan 2025), so a rat-human difference cannot be separated "
               "from a laboratory effect."])
    wb.save(OUT)
    print(f"\nwrote {OUT}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
