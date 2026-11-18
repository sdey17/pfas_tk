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
The sheet carries a Species column and the audit reports the comparison, but
no row is pooled across species.

GATES, in the order the builder applies them:

  holdout        Han 2021. Its hL-FABP constants are quoted from other papers,
                 not measured -- the row notes say so, and its PFBS value is
                 Zhang 2013's. Its own measurements are thermal-shift dTm,
                 which is not a dissociation constant.
  computational  Zhao 2023's rows are docking scores and QSAR summaries over a
                 1599-compound library, with no measured constant.
  ENGINEERED     Sheng 2016 measured WILD TYPE hL-FABP and FOUR SITE-DIRECTED
  PROTEIN        POINT MUTANTS (S39G, M74G, N111D, R122G), and the mutants are
                 not the protein this sheet is about. Three of them were built
                 precisely to break binding: N111D loses it outright (Kd 284 uM
                 for PFOA against 6.49 uM for wild type, a 44x loss) and R122G
                 halves the site count (ITC N = 0.554 and 0.562, which is what
                 a sub-unity "stoichiometry" in this source means). Pooling
                 them put PFOA into the sheet EIGHT times spanning 2.36 to
                 284 uM -- a 120-fold spread that is five different proteins,
                 not replication. 22 rows go; the 7 wild-type rows stay.
  not a PFAS     the fluorescent probes these assays displace -- DAUDA,
                 1,8-ANS, ANS, NBD-lauric acid -- and L-FABP's natural fatty
                 acid ligands. These are the reference ligands, not analytes.
                 Palmitic acid in particular reads as a strong binder (log K
                 5.02 from Zhang 2013) and would be the third highest value in
                 the sheet, sitting inside a PFAS chain-length relationship.
  several        a cell naming several compounds is a summary. Tested the
  compounds      data-driven way, as in the albumin builder: split on commas
                 outside brackets and count how many parts appear elsewhere as
                 standalone compounds, because an IUPAC name's commas are
                 locants.
  lower-affinity A multi-site fit gives one constant per site. Only the
  site of a      HIGH-AFFINITY site reaches the sheet, which is the rule the
  multi-site fit albumin sheet already follows -- Moro 2022 sits in it as
                 "ITC, sequential site 1 of 2, 37 C". The remaining sites stay
                 in the database as the record. This UNDERSTATES total
                 capacity for the few compounds with resolved extra sites:
                 Woodcroft 2010's PFOA by ITC carries 1/3.13 + 1/26.25 +
                 1/52.56 = 0.377 uM^-1 across three sites against 0.319 for
                 site 1 alone, so the sheet's value is 0.07 log low there.
                 Named rather than corrected, because a one-row-per-constant
                 sheet cannot carry a site series.
  salt duplicate Yang 2020 lists PFOA and PFOA ammonium salt, PFOS and PFOS
                 potassium salt, at IDENTICAL Kd: one measurement under two
                 names. The free-acid name is kept and the salt row dropped.
                 PFBS and PFBS potassium salt differ (116.4 vs 125.25) and are
                 kept as two.
  outside the    Khazaee 2021's PFHxS, Kd 1.7 uM against a ligand series
  range tested   topping out at 0.384 uM. A fit beyond the top of the series
                 is an extrapolation of the model, not a measurement. Its PFOA
                 (0.099 uM) and PFOS (0.18 uM) sit inside the series and stay.

ROUTE LABELS follow the albumin sheet's vocabulary, assigned from what each
paper actually reports. The one correction worth naming:

  SHENG 2016's FLUORESCENCE n IS THE PROBE'S, NOT THE PFAA'S. Table 2's n
  column (1.01, 1.19, 1.07) sits under the "1,8-ANS" header and counts ANS
  molecules per protein. The PFAA numbers in that table are IC50 and a Kd
  derived from it by a competition relation that assumes 1:1 at the probe's
  site. So those rows are "n = 1 by single-site model", not a fitted PFAA
  stoichiometry. Numerically this moves log K by 0.004; the label was simply
  describing the wrong molecule.

  Sheng's ITC n, by contrast, IS the PFAA's and IS counted: Table 3's N (sites)
  column gives 0.999 for PFOA, 1.20 for PFNA, 0.987 for PFHxS on wild type.
  Together with Woodcroft 2010 that makes two sources, not one, with counted
  L-FABP stoichiometries. Sheng fits Origin's One Set of Sites model, so those
  rows take the albumin sheet's existing label for exactly that case --
  "ITC, one set of sites, <T> C", as Moro 2022 carries it there.

ROUTE VOCABULARY IS THE ALBUMIN SHEET'S, not a parallel one. Every label used
here either already appears in albumin_logK_Kd_FINAL.xlsx ("n = 1 assumed",
"n = 1 by single-site model", "n fitted by binding isotherm", "ITC, one set of
sites, <T> C", "ITC, sequential site 1 of N, <T> C") or extends that pattern
to a case albumin has no instance of ("displacement, site 1 of 2, <T> C").

CONCENTRATION COLUMNS. These are NOT parsed out of the database's free-text
condition fields, because the leading number in those fields is not reliably
the quantity wanted and the units column does not always match: Zhang's
"400 nM hL-FABP with 50 nM NBD-lauric acid" read naively gives 400 uM, a
1000-fold error, and Sheng's ligand field begins with the PROTEIN
concentration. Every value below was read from the primary and carries its
provenance in CONDITIONS. A blank means the paper does not state it; the
Notes sheet lists which papers those are. ligand:protein is always computed as
PFAS conc / protein conc, never stored, which is the albumin sheet's rule.

Run:  python3 -I scripts/build_lfabp_working_table.py
"""
import collections, csv, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "lfabp_logK_Kd.xlsx")

MW = {"human": 14208.0, "rat": 14273.0, "mouse": 14246.0}
HOLDOUT = {"Han 2021": "constants quoted from other papers, not measured here"}
COMPUTATIONAL = {"Zhao 2023"}
ENGINEERED = re.compile(r"ENGINEERED POINT VARIANT|\b(S39G|M74G|N111D|R122G)\b")
NOT_ANALYTE = re.compile(
    r"^(DAUDA|1,8-ANS|ANS\b|NBD|oleic acid|palmitic acid|stearic acid|"
    r"lauric acid|myristic acid|arachidonic acid)", re.I)
LOWER_SITE = re.compile(r"site\s+([2-9])\s+of", re.I)
HIGH_SITE = re.compile(r"site\s+1\s+of\s+(\d)", re.I)

# Papers whose own model fixes n = 1. Woodcroft 2010 belongs here for the
# three PFCAs it fits to a ONE-site competitive model (PFHpA, PFHxA, PFPeA);
# its two-site and sequential-ITC rows never reach this branch because the
# resolved-site test above them catches those first.
PAPER_SINGLE_SITE = {"Fan 2025", "Zhang 2013", "Khazaee 2021", "Woodcroft 2010"}
ITC_COUNTED = {"Sheng 2016"}        # Table 3 N (sites) is the PFAA's own
OUTSIDE_RANGE = {("Khazaee 2021", "pfhxs")}

# --- assay concentrations, read from each primary --------------------------
# (source, method class) -> (protein_uM, ligand_uM, provenance)
# None means the paper does not state it.
KHAZAEE_MW = {"pfoa": 414.07, "pfos": 500.13, "pfhxs": 400.11}
KHAZAEE_NG_ML = 153.5               # top of the stated ligand series, ng/mL
KHAZAEE_RATIO = 5.0                 # the paper's own top ligand:protein ratio

CONDITIONS = {
    ("Fan 2025", "fluorescence"): (
        1.0, None,
        "1 uM rL-FABP in 100 uL per well, 1,8-ANS at saturation; the PFAS "
        "series is 'different concentrations' in SM Text S1 and is not given "
        "in the paper"),
    ("Gao 2019", "dialysis"): (
        6.0, 12.0,
        "0.006 mM (6 uM) hL-FABP; six-point ligand series 0.002-0.012 mM "
        "(2-12 uM), top of series quoted"),
    ("Jia 2022", "MST"): (
        5.0, None,
        "10 uM RED-NHS-labelled hL-FABP mixed 1:1 with the ligand series, so "
        "5 uM in the capillary; SI Text S7 describes the serial dilution but "
        "never states its concentrations"),
    ("Sheng 2016", "fluorescence"): (
        0.25, None,
        "0.25 uM (250 nM) hL-FABP with 1,8-ANS at 12.25 uM; per-compound "
        "ligand series taken from the figures where plotted"),
    ("Woodcroft 2010", "fluorescence"): (
        1.0, None,
        "1 uM rat L-FABP pre-bound with 70 uM ANS in 2.5 ml; Fig. 3 plots "
        "log[competitor] from -6 to -2 (1 uM to 10 mM) but the per-compound "
        "endpoints are not tabulated"),
    ("Woodcroft 2010", "ITC"): (
        162.0, 4050.0,
        "0.162 mM L-FABP in the sample cell, 4.05 mM PFOA in the syringe"),
    ("Yang 2020", "fluorescence"): (
        1.0, 300.0,
        "1 uM hL-FABP with 1,8-ANS at 40 uM; 300 uM screening concentration "
        "from 30 mM DMSO stocks, dose-response series not tabulated"),
    ("Zhang 2013", "fluorescence"): (
        0.4, None,
        "400 nM hL-FABP with 50 nM NBD-lauric acid, a 1:8 probe:protein "
        "ratio; the per-compound PFC series is SI Figure S3, not tabulated"),
}

# Ligand series a figure in the primary does give, per (source, compound).
LIGAND_FROM_FIGURE = {
    ("Sheng 2016", "pfoa"): (1000.0,
                             "Fig. 2c, log[PFOA] plotted at -1, 0, 1, 2, 3 "
                             "(0.1-1000 uM); top of series"),
}

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


def method_class(r):
    m = (r.get("method") or "").lower()
    if "isothermal titration" in m or "itc" in m:
        return "ITC"
    if "equilibrium dialysis" in m or "dialysis" in m:
        return "dialysis"
    if "thermophoresis" in m or "mst" in m:
        return "MST"
    return "fluorescence"


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


def concentrations(src, mcls, compound, row):
    """(protein_uM, ligand_uM) for one row, from the primary. None = unstated."""
    key = compound.strip().lower()
    if src == "Khazaee 2021":                       # per compound, from ng/mL
        mw = KHAZAEE_MW.get(key)
        if mw:                      # unrounded, so the ratio comes out at 5.0
            lig = KHAZAEE_NG_ML / mw
            return lig / KHAZAEE_RATIO, lig
        return None, None
    if (src, mcls) == ("Sheng 2016", "ITC"):        # varies by row, SI Table S2
        p = re.search(r"([\d.]+)\s*uM protein in the cell",
                      row.get("protein_conc") or "")
        l = re.search(r"([\d.]+)\s*mM in the syringe",
                      row.get("ligand_conc_range") or "")
        return (float(p.group(1)) if p else None,
                float(l.group(1)) * 1e3 if l else None)
    prot, lig, _ = CONDITIONS.get((src, mcls), (None, None, ""))
    fig = LIGAND_FROM_FIGURE.get((src, key))
    if fig:
        lig = fig[0]
    return prot, lig


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
        prot_desc = r.get("protein") or ""
        if src in HOLDOUT:
            dropped["holdout: " + src] += 1; continue
        if src in COMPUTATIONAL:
            dropped["computational: " + src] += 1; continue
        if ENGINEERED.search(prot_desc):
            dropped["engineered point mutant, not the native protein"] += 1
            continue
        if sp not in MW:
            dropped[f"species not carried ({sp})"] += 1; continue
        if NOT_ANALYTE.match(ch.strip()):
            dropped["not an analyte (probe or natural ligand)"] += 1; continue
        if names_several(ch, known):
            dropped["several compounds in one cell"] += 1; continue
        site = r.get("binding_site_identified") or ""
        if LOWER_SITE.search(site):
            dropped["lower-affinity site of a multi-site fit"] += 1; continue
        if "," in str(r.get("kd") or ""):
            dropped["composite Kd (several sites in one cell)"] += 1; continue
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

        base = salt_base(ch)
        key = (src, sp, base, round(kd_uM, 6))
        if key in seen_salt:
            dropped["salt form of a row already kept"] += 1; continue
        seen_salt[key] = True

        mcls = method_class(r)
        n_raw = num(r.get("n_binding_sites"))
        hs = HIGH_SITE.search(site)
        if hs:                                  # site 1 of a resolved series
            n = 1.0
            tag = "ITC" if mcls == "ITC" else "displacement"
            if "sequential" in (r.get("n_binding_sites") or "").lower():
                tag += ", sequential"     # "ITC, sequential site 1 of 3, 30 C"
            T = (r.get("temperature_c") or "").strip()
            route = (f"{tag} site 1 of {hs.group(1)}" if tag.endswith("sequential")
                     else f"{tag}, site 1 of {hs.group(1)}")
            route += f", {T} C" if T else ""
        elif src in ITC_COUNTED and mcls == "ITC" and n_raw:
            # Sheng fits Origin's One Set of Sites model and reports its N, so
            # this is the albumin sheet's "ITC, one set of sites, 37 C" case
            # (Moro 2022), word for word, with its own temperature.
            T = (r.get("temperature_c") or "").strip()
            n = n_raw
            route = "ITC, one set of sites" + (f", {T} C" if T else "")
        elif src in PAPER_SINGLE_SITE or (src in ITC_COUNTED
                                          and mcls == "fluorescence"):
            n, route = 1.0, "n = 1 by single-site model"
        elif n_raw:
            n, route = n_raw, "n fitted by binding isotherm"
        else:
            n, route = 1.0, "n = 1 assumed"

        p_uM, l_uM = concentrations(src, mcls, ch, r)
        # The albumin sheet PREFERS the ratio the paper itself reports and
        # computes one only as a fallback. Allendorf 2019 is why: its 28 rows
        # carry the measured analyte:BSA molar ratio at equilibrium from SI
        # Table SI-2.3, deliberately in place of the 0.1 design ceiling and in
        # place of any quotient of the two concentration columns, which hold
        # c_ref and a g/L protein loading. Same rule here. On the present
        # L-FABP rows the two agree wherever both exist -- Gao's 2, Sheng's
        # 13/15 from SI Table S2, Khazaee's 5 -- so this changes no value; it
        # makes the RULE the same one, which is what matters if a source that
        # reports a ratio its concentrations do not reproduce is added later.
        reported = num(r.get("ligand_over_protein_max"))
        computed = l_uM / p_uM if (p_uM and l_uM) else None
        ratio = round(reported if reported is not None
                      else computed, 4) if (reported is not None
                                            or computed is not None) else ""
        p_uM = round(p_uM, 5) if p_uM is not None else None
        l_uM = round(l_uM, 5) if l_uM is not None else None
        bridge = math.log10(1e9 / MW[sp])
        kept.append({
            "PFAS": base, "Species": sp, "Kd (uM)": round(kd_uM, 4),
            "n used": n,
            "log K (L water / kg protein)":
                round(bridge + math.log10(n) - math.log10(kd_uM), 3),
            "Route": route, "Source": src,
            "PFAS conc (uM)": l_uM if l_uM is not None else "",
            "Protein conc (uM)": p_uM if p_uM is not None else "",
            "ligand:protein ratio": ratio,
            "f_bound (%)": num(r.get("unbound_fraction")) or ""})

    kept.sort(key=lambda d: (d["Species"], d["Source"], d["PFAS"]))
    print(f"L-FABP sheet: {len(kept)} rows from "
          f"{len({d['Source'] for d in kept})} sources, "
          f"{len({d['PFAS'] for d in kept})} compounds")
    for sp in ("human", "rat", "mouse"):
        v = sorted(d["log K (L water / kg protein)"]
                   for d in kept if d["Species"] == sp)
        if v:
            labs = sorted({d["Source"] for d in kept if d["Species"] == sp})
            print(f"   {sp:6s} n={len(v):3d}  log K {v[0]:.2f} to {v[-1]:.2f}  "
                  f"median {v[len(v)//2]:.2f}  ({len(labs)} sources)")
    print("\nroutes:")
    for rt, n in collections.Counter(d["Route"] for d in kept).most_common():
        print(f"   {n:4d}  {rt}")
    print("\nconcentration columns filled:")
    for c in ("PFAS conc (uM)", "Protein conc (uM)", "ligand:protein ratio"):
        print(f"   {sum(1 for d in kept if d[c] != ''):4d}/{len(kept)}  {c}")
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
    nb.append(["L-FABP partition coefficients"]); nb.append([])
    nb.append(["bridge", "log K = log10(1e9/M) + log10(n) - log10(Kd/uM)"])
    for sp, m in sorted(MW.items()):
        nb.append([f"  {sp} L-FABP",
                   f"M = {m:.0f} Da -> constant {math.log10(1e9/m):.4f}"])
    nb.append(["  albumin, for contrast", "M = 66463 Da -> constant 4.1774"])
    nb.append([])
    nb.append(["ASSAY CONCENTRATIONS", "read from each primary; "
               "ligand:protein is computed as PFAS conc / protein conc"])
    for (s, m), (p, l, why) in sorted(CONDITIONS.items()):
        nb.append([f"  {s} ({m})",
                   f"protein {p if p is not None else 'not stated'} uM, "
                   f"ligand {l if l is not None else 'NOT STATED'} uM", why])
    nb.append(["  Khazaee 2021 (dialysis)",
               "per compound: ligand = 153.5 ng/mL / MW, protein = ligand / 5",
               "the paper prints the series as 0.33-153.5 ng/mL across "
               "ligand:protein mole ratios 0.05-5 and prints no molar "
               "protein concentration"])
    for (s, c), (v, why) in sorted(LIGAND_FROM_FIGURE.items()):
        nb.append([f"  {s} / {c.upper()}", f"ligand {v} uM", why])
    nb.append([])
    nb.append(["EXCLUDED"])
    for k, n in sorted(dropped.items(), key=lambda x: -x[1]):
        nb.append([f"  {k}", n])
    nb.append([])
    nb.append(["Multi-site fits contribute their HIGH-AFFINITY site only, as "
               "in the albumin sheet. Total capacity is therefore understated "
               "where a paper resolves further sites: Woodcroft 2010's PFOA "
               "by ITC is 0.07 log low on that account."])
    nb.append(["Species are NOT pooled. Both the rat rows (Fan 2025, "
               "Woodcroft 2010) and the human rows come from more than one "
               "laboratory, so a species difference is no longer confounded "
               "with a single laboratory on either side."])
    nb.append(["Sheng 2016 contributes WILD-TYPE hL-FABP only. Its four "
               "engineered point mutants (S39G, M74G, N111D, R122G) are "
               "excluded: three were built to break binding, and pooling them "
               "put PFOA in the sheet eight times across a 120-fold range."])
    wb.save(OUT)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
