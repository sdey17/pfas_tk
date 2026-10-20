#!/usr/bin/env python3
"""Consolidate every extraction in this review into four combined dataset files.

The review accumulated ~30 CSVs in db/ plus 13 per-paper extractions in
db/primary_2026/. They use different column names, different units and
different provenance conventions because they were built by different passes.
This script merges them into four files with one schema each, normalised units
and a provenance column on every row:

  combined/tk_parameters.csv   one row per chemical x species x sex x source:
                               half-life, clearance, Vd, with units normalised
  combined/binding.csv         protein binding and unbound fractions
  combined/transporters.csv    transporter kinetics, expression and direction
  combined/regulatory.csv      what each agency adopted and from where

Nothing is invented here. Where a source gave a range, both ends are kept in
`value_low`/`value_high`. Where a unit conversion was applied, the original
value and unit are preserved alongside.

Run:  python3 scripts/build_combined_datasets.py
"""
import csv
import os
import sys
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db")
PRIM = os.path.join(DB, "primary_2026")
OUT = os.path.join(DB, "combined")

TK_FIELDS = [
    "chemical", "species", "strain", "sex", "parameter", "value", "value_low",
    "value_high", "units", "dose", "dose_units", "route", "method",
    "n_animals", "study", "year", "pmid_or_doi", "source_table", "provenance",
    "notes",
]


def rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def num(x):
    """First number in a string, or None. Keeps ranges out of `value`."""
    if x is None:
        return None
    m = re.match(r"^\s*(-?\d+(?:\.\d+)?)(?:[eE][-+]?\d+)?\s*$", str(x))
    return float(m.group(1)) if m else None


def rng(x):
    """(low, high) if the string is a range like '34.3-68.9', else (None, None)."""
    if not x:
        return None, None
    m = re.match(r"^\s*(-?\d+(?:\.\d+)?)\s*[-–]\s*(-?\d+(?:\.\d+)?)\s*$", str(x))
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


def rec(**kw):
    r = {f: "" for f in TK_FIELDS}
    r.update({k: v for k, v in kw.items() if v is not None})
    return r


# ----------------------------------------------------------------- TK table
def build_tk():
    out = []

    # --- existing review tables -------------------------------------------
    for r in rows(os.path.join(DB, "animal_halflife_measured.csv")):
        out.append(rec(
            chemical=r.get("chemical"), species=r.get("species"),
            strain=r.get("strain"), sex=r.get("sex"), parameter="halflife",
            value=num(r.get("halflife")), units=r.get("halflife_units"),
            dose=r.get("dose"), dose_units=r.get("dose_units"),
            route=r.get("route"), method=r.get("model"),
            n_animals=r.get("n"), study=r.get("study"), year=r.get("year"),
            pmid_or_doi=r.get("pmid") or r.get("doi"),
            source_table=r.get("source_table"), provenance="animal_halflife_measured.csv",
            notes=r.get("notes")))

    for r in rows(os.path.join(DB, "vd_clearance.csv")):
        for par, val, un in (("vd", r.get("vd"), r.get("vd_units")),
                             ("clearance", r.get("clearance"), r.get("clearance_units")),
                             ("halflife", r.get("halflife"), r.get("halflife_units"))):
            if not val:
                continue
            lo, hi = rng(val)
            out.append(rec(
                chemical=r.get("chemical"), species=r.get("species"),
                strain=r.get("strain"), sex=r.get("sex"), parameter=par,
                value=num(val), value_low=lo, value_high=hi, units=un,
                method=r.get("vd_method_measured_fitted_assumed") if par == "vd" else "",
                n_animals=r.get("n"), study=r.get("study"), year=r.get("year"),
                pmid_or_doi=r.get("pmid") or r.get("doi"),
                source_table=r.get("source_table"), provenance="vd_clearance.csv",
                notes=r.get("vd_provenance_cited_from") if par == "vd" else ""))

    for r in rows(os.path.join(DB, "human_halflife_extended.csv")):
        out.append(rec(
            chemical=r.get("chemical"), species="human", sex=r.get("sex"),
            parameter="halflife", value=num(r.get("halflife")),
            units=r.get("halflife_units") or "y", method=r.get("design"),
            n_animals=r.get("n"), study=r.get("study"), year=r.get("year"),
            pmid_or_doi=r.get("pmid") or r.get("doi"),
            provenance="human_halflife_extended.csv", notes=r.get("assumptions")))

    for r in rows(os.path.join(DB, "cphea_fitted_halflives.csv")):
        out.append(rec(
            chemical=r.get("chemical"), species=r.get("species"),
            strain=r.get("strain"), sex=r.get("sex"), parameter="halflife",
            value=num(r.get("halflife_days")), units="d", dose=r.get("dose"),
            dose_units=r.get("dose_units"), route=r.get("route"),
            method="terminal log-linear fit computed in this review",
            study="EPA CPHEA " + (r.get("study") or ""),
            provenance="cphea_fitted_halflives.csv",
            notes=("TAIL-SELECTION FLAG: " + r["tail_selection_flag"]
                   if r.get("tail_selection_flag") else
                   "see REPORT.md section 8 item 9 before using")))

    # --- the nine primary full texts --------------------------------------
    for r in rows(os.path.join(PRIM, "argoul2026_mouse_tk.csv")):
        base = dict(chemical=r["chemical"], species="mouse", strain="CD-1",
                    sex="female", route="IV and oral", dose=r.get("iv_dose_mgkg"),
                    dose_units="mg/kg", method="NLME, 119-day sampling",
                    study="Argoul et al.", year="2026",
                    pmid_or_doi="42162713 / 10.1016/j.envres.2026.124802",
                    source_table="Table 1", provenance="primary_2026/argoul2026_mouse_tk.csv")
        for par, key, un in (("clearance", "cl_mL_kg_day", "mL/kg-day"),
                             ("vd_ss", "vss_L_kg", "L/kg"),
                             ("mrt", "mrt_d", "d"),
                             ("bioavailability", "bioavailability_pct", "%")):
            if r.get(key):
                lo, hi = rng(r[key])
                out.append(rec(parameter=par, value=num(r[key]), value_low=lo,
                               value_high=hi, units=un, **base))

    for r in rows(os.path.join(PRIM, "kudo2002_rat_tk.csv")):
        if r["parameter"] in ("total_clearance", "renal_clearance_untreated",
                              "renal_clearance_probenecid"):
            un = "mL/min/kg"
        elif r["parameter"] == "total_clearance_per_day":
            un = "mL/kg-day"
        elif r["parameter"] == "halflife":
            un = "d"
        else:
            un = "mL/kg"
        for sex in ("male", "female"):
            out.append(rec(
                chemical="PFOA", species="rat", strain="Wistar", sex=sex,
                parameter=r["parameter"], value=num(r[sex]), units=un,
                dose="20.14", dose_units="mg/kg (48.63 umol/kg)", route="iv",
                method="measured, both sexes in one experiment",
                study="Kudo et al.", year="2002",
                pmid_or_doi="11879818 / 10.1016/s0009-2797(02)00006-6",
                source_table=r["source_table"],
                provenance="primary_2026/kudo2002_rat_tk.csv", notes=r.get("note")))

    for r in rows(os.path.join(PRIM, "lou2009_mouse_tk.csv")):
        un = {"vd": "L/kg", "ke": "1/h", "halflife": "d",
              "clearance_derived": "mL/kg-day"}.get(r["parameter"], "")
        for sex in ("male", "female"):
            out.append(rec(
                chemical="PFOA", species="mouse", strain="CD-1", sex=sex,
                parameter=f"{r['parameter']} ({r['matrix']})",
                value=num(r[sex]), units=un, dose="1 and 10", dose_units="mg/kg",
                route="oral gavage", method="one-compartment model",
                study="Lou et al.", year="2009",
                pmid_or_doi="19005225 / 10.1093/toxsci/kfn234",
                source_table=r["source_table"],
                provenance="primary_2026/lou2009_mouse_tk.csv",
                notes=f"95% CI {r[sex + '_ci']}" if r.get(sex + "_ci") else ""))

    for r in rows(os.path.join(PRIM, "tatumgibbs2011_pfna_rat_mouse.csv")):
        lo, hi = rng(r["halflife_d"])
        out.append(rec(
            chemical="PFNA", species=r["species"], strain=r["strain"],
            sex=r["sex"], parameter="halflife", value=num(r["halflife_d"]),
            value_low=lo, value_high=hi, units="d", dose=r["dose_mgkg"],
            dose_units="mg/kg", route=r["route"],
            method=f"dose linearity: {r['dose_linearity']}",
            study="Tatum-Gibbs et al.", year="2011",
            pmid_or_doi="10.1016/j.tox.2011.01.003",
            source_table=r["source"],
            provenance="primary_2026/tatumgibbs2011_pfna_rat_mouse.csv",
            notes=r.get("note")))

    for r in rows(os.path.join(PRIM, "kudo2001_chain_length_elimination.csv")):
        if r["pct_dose_urine_120h"]:
            out.append(rec(
                chemical=r["chemical"], species="rat", strain="Wistar",
                sex=r["sex"], parameter="pct_dose_in_urine_120h",
                value=num(r["pct_dose_urine_120h"]), units="% of dose",
                route="intraperitoneal", method="cumulative recovery to 120 h",
                study="Kudo et al.", year="2001",
                pmid_or_doi="11311214 / 10.1016/s0009-2797(01)00155-7",
                source_table=r["source"],
                provenance="primary_2026/kudo2001_chain_length_elimination.csv",
                notes=r.get("route_verdict")))

    for r in rows(os.path.join(PRIM, "sundstrom2012_pfhxs_three_species.csv")):
        base = dict(chemical="PFHxS", species=r["species"], strain=r["strain"],
                    sex=r["sex"], dose=r["dose_mgkg"], dose_units="mg/kg",
                    route=r["route"], n_animals=r["n"],
                    study="Sundstrom et al.", year="2012",
                    pmid_or_doi="21856411 / 10.1016/j.reprotox.2011.07.004",
                    source_table=r["source_table"],
                    provenance="primary_2026/sundstrom2012_pfhxs_three_species.csv",
                    notes=r.get("note"))
        for par, key, un, meth in (
                ("halflife", "halflife_d", "d", r["halflife_type"]),
                ("clearance", "clearance_mL_d_kg", "mL/kg-day", ""),
                ("vd", "vd_mL_kg", "mL/kg", r["vd_type"]),
                ("pct_dose_in_urine_24h", "pct_dose_urine_24h", "% of dose", "")):
            if r.get(key):
                out.append(rec(parameter=par, value=num(r[key]), units=un,
                               method=meth, **base))

    for r in rows(os.path.join(PRIM, "ohmori2003_rat_chain_length.csv")):
        for sex, key in (("male", "halflife_male_d"), ("female", "halflife_female_d")):
            out.append(rec(
                chemical=r["chemical"], species="rat", sex=sex,
                parameter="halflife", value=num(r[key]), units="d", route="iv",
                method="terminal phase", study="Ohmori et al.", year="2003",
                pmid_or_doi="12499116", source_table="abstract",
                provenance="primary_2026/ohmori2003_rat_chain_length.csv",
                notes=r.get("note")))

    for r in rows(os.path.join(PRIM, "delaere2025_firefighter_removal.csv")):
        out.append(rec(
            chemical=r["chemical"], species="human", parameter="halflife",
            value=num(r["apparent_halflife_y"]), units="y",
            method=f"one-compartment first-order; {r['group']}",
            n_animals=r["n"], study="Delaere et al.", year="2025",
            pmid_or_doi="Environ Int 202:109609",
            source_table="abstract",
            provenance="primary_2026/delaere2025_firefighter_removal.csv",
            notes="CAVEAT: only participants whose concentrations decreased were "
                  "included, and no statistical comparison was made"))

    for r in rows(os.path.join(PRIM, "shi2016_clpfesa_human.csv")):
        out.append(rec(
            chemical=r["chemical"], species="human", parameter=r["parameter"],
            value=num(r["median"]), units=r["units"], method="one-compartment",
            study="Shi et al.", year="2016",
            pmid_or_doi="26866980 / 10.1021/acs.est.5b05849",
            source_table="abstract",
            provenance="primary_2026/shi2016_clpfesa_human.csv",
            notes=f"range {r['range']}; {r['note']}" if r["range"] else r["note"]))

    for r in rows(os.path.join(PRIM, "han2012_table4_reabsorption_axis.csv")):
        for par, key, un in (("gfr", "gfr_L_d_kg", "L/d/kg"),
                             ("renal_clearance", "clr_mL_d_kg", "mL/d/kg"),
                             ("pct_renal_reabsorption", "pct_reabsorption", "%")):
            if r.get(key):
                out.append(rec(
                    chemical="PFOA", species=r["species"], sex=r["sex"],
                    parameter=par, value=num(r[key]), units=un,
                    method="compilation; fu assumed 0.02 for all species",
                    study="Han et al.", year="2012",
                    pmid_or_doi="21985250 / 10.1021/tx200363w",
                    source_table="Table 4",
                    provenance="primary_2026/han2012_table4_reabsorption_axis.csv",
                    notes=r.get("note")))

    os.makedirs(OUT, exist_ok=True)
    out = [r for r in out if r["value"] != "" or r["value_low"] != ""]
    with open(os.path.join(OUT, "tk_parameters.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TK_FIELDS)
        w.writeheader()
        w.writerows(out)
    return out


def binding_rows():
    """The measurement rows that db/binding_affinity.csv is written from."""
    sys.path.insert(0, HERE)
    from build_binding_affinity_table import build as _build
    return _build()[1]


# ------------------------------------------------------------- binding table
BIND_FIELDS = ["chemical", "species", "protein", "measure", "value", "units",
               "method", "study", "year", "pmid_or_doi", "provenance", "notes"]


def build_binding():
    """One row per binding measurement, from db/protein_binding.csv.

    This used to read binding_summary_ka.csv and binding_summary_logD.csv,
    which are per-chemical AGGREGATES carrying ka_M-1_geomean / studies, not
    value / study. Every .get() missed and 311 of 327 rows came out with
    species, value, units, method, study, year and the identifier all blank.
    The per-measurement table is the right source and is shared with
    db/binding_affinity.csv so the two cannot drift apart.
    """
    out = []
    for r in binding_rows():
        out.append({
            "chemical": r["pfas"], "species": r["species"],
            "protein": r["protein"], "measure": r["measure"],
            "value": r["value"], "units": r["units"],
            "method": r["_method"], "study": r["_study"], "year": r["_year"],
            "pmid_or_doi": r["_id"], "provenance": "protein_binding.csv",
            "notes": r["location"]})
    for r in rows(os.path.join(PRIM, "han2003_pfoa_albumin_binding.csv")):
        out.append({
            "chemical": "PFOA", "species": "rat and human",
            "protein": r["protein"], "measure": r["parameter"],
            "value": r["value"], "units": r["units"], "method": r["method"],
            "study": "Han et al.", "year": "2003",
            "pmid_or_doi": "12807361 / 10.1021/tx034005w",
            "provenance": "primary_2026/han2003_pfoa_albumin_binding.csv",
            "notes": f"ligand:protein {r['ligand_protein_molar_ratio']}; {r['note']}"})

    for r in rows(os.path.join(PRIM, "argoul2026_mouse_tk.csv")):
        if r.get("fu_pct"):
            out.append({
                "chemical": r["chemical"], "species": "mouse (CD-1, female)",
                "protein": "plasma (whole)", "measure": "unbound fraction fu",
                "value": r["fu_pct"], "units": "%", "method": "equilibrium dialysis",
                "study": "Argoul et al.", "year": "2026",
                "pmid_or_doi": "42162713 / 10.1016/j.envres.2026.124802",
                "provenance": "primary_2026/argoul2026_mouse_tk.csv",
                "notes": "Table 4; PFBA value from Ryu 2024"})
    with open(os.path.join(OUT, "binding.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BIND_FIELDS)
        w.writeheader()
        w.writerows(out)
    return out


# --------------------------------------------------------- transporter table
TR_FIELDS = ["transporter", "species", "chemical", "measure", "value", "units",
             "direction", "condition", "study", "year", "pmid_or_doi",
             "provenance", "notes"]


def build_transporters():
    out = []
    for r in rows(os.path.join(DB, "transporter_kinetics.csv")):
        out.append({
            "transporter": r.get("transporter", ""), "species": r.get("species", ""),
            "chemical": r.get("chemical", ""), "measure": r.get("parameter", "Km"),
            "value": r.get("value", ""), "units": r.get("units", ""),
            "direction": r.get("direction", ""), "condition": "",
            "study": r.get("study", ""), "year": r.get("year", ""),
            "pmid_or_doi": r.get("pmid", "") or r.get("doi", ""),
            "provenance": "transporter_kinetics.csv", "notes": r.get("notes", "")})
    for r in rows(os.path.join(PRIM, "han2012_table6_transporter_km.csv")):
        if not r["km_uM"]:
            continue
        out.append({
            "transporter": r["transporter"], "species": r["species"],
            "chemical": r["chemical"], "measure": "Km", "value": r["km_uM"],
            "units": "uM", "direction": "", "condition": "",
            "study": "Han et al.", "year": "2012",
            "pmid_or_doi": "21985250 / 10.1021/tx200363w",
            "provenance": "primary_2026/han2012_table6_transporter_km.csv",
            "notes": r["note"]})
    for r in rows(os.path.join(PRIM, "yang2010_human_apical_transporters.csv")):
        out.append({
            "transporter": r["transporter"], "species": "human",
            "chemical": r["substrate"], "measure": "Km",
            "value": r["km_uM"] or "no saturable uptake", "units": "uM",
            "direction": "apical / reabsorptive", "condition": r["condition"],
            "study": "Yang et al.", "year": "2010",
            "pmid_or_doi": "20639259 / 10.1093/toxsci/kfq219",
            "provenance": "primary_2026/yang2010_human_apical_transporters.csv",
            "notes": r["note"]})
    for r in rows(os.path.join(PRIM, "louisse2024_human_oat_km.csv")):
        out.append({
            "transporter": r["transporter"], "species": "human",
            "chemical": r["chemical"], "measure": "Km",
            "value": r["km_uM"] or "no transport observed", "units": "uM",
            "direction": "basolateral / secretory" if r["transporter"] in
                         ("OAT1", "OAT2", "OAT3") else "apical / reabsorptive",
            "condition": f"Vmax {r['vmax_nmol_min_mg']} nmol/min/mg; "
                         f"efficiency {r['efficiency_uL_min_mg']} uL/min/mg"
                         if r["vmax_nmol_min_mg"] else "",
            "study": "Louisse et al.", "year": "2024",
            "pmid_or_doi": "10.1016/j.tox.2024.153961",
            "provenance": "primary_2026/louisse2024_human_oat_km.csv",
            "notes": f"SE {r['km_se']}" if r["km_se"] else r["source"]})

    for r in rows(os.path.join(PRIM, "cheng2005_mouse_oatp_sex.csv")):
        out.append({
            "transporter": r["transporter"], "species": "mouse (C57BL/6)",
            "chemical": "(endogenous organic anions)",
            "measure": f"mRNA sex predominance, {r['tissue']}",
            "value": r["fold"] or r["sex_predominance"],
            "units": "fold" if r["fold"] else "direction only",
            "direction": "apical / reabsorptive" if r["transporter"] == "Oatp1a1" else "",
            "condition": f"androgen-dependent: {r['androgen_dependent']}"
                         if r["androgen_dependent"] else "",
            "study": "Cheng et al.", "year": "2005",
            "pmid_or_doi": "15843488 / 10.1124/dmd.105.003640",
            "provenance": "primary_2026/cheng2005_mouse_oatp_sex.csv",
            "notes": r["note"]})

    for r in rows(os.path.join(PRIM, "kudo2002_transporter_mrna.csv")):
        for key, meas in (("male_over_female_fold", "renal mRNA, male/female fold"),
                          ("female_over_male_fold", "renal mRNA, female/male fold")):
            if r.get(key):
                out.append({
                    "transporter": r["transporter"], "species": "rat",
                    "chemical": "PFOA", "measure": meas, "value": r[key],
                    "units": "fold", "direction": r["direction"], "condition": "",
                    "study": "Kudo et al.", "year": "2002",
                    "pmid_or_doi": "11879818 / 10.1016/s0009-2797(02)00006-6",
                    "provenance": "primary_2026/kudo2002_transporter_mrna.csv",
                    "notes": f"transports PFOA: {r['transports_pfoa_direct_assay']}; "
                             f"{r['hormonal_response']}"})
    for r in rows(os.path.join(PRIM, "han2012_table7_pbpk_tm_kt.csv")):
        for key, meas, un in (("tmc_mg_h_kg", "PBPK transport maximum Tmc", "mg/h/kg"),
                              ("kt_mg_L", "PBPK transport affinity KT", "mg/L")):
            out.append({
                "transporter": "(lumped renal reabsorption)", "species": r["species"],
                "chemical": "PFOA", "measure": meas, "value": r[key], "units": un,
                "direction": "reabsorptive", "condition": "fitted to plasma curves",
                "study": "Han et al.", "year": "2012",
                "pmid_or_doi": "21985250 / 10.1021/tx200363w",
                "provenance": "primary_2026/han2012_table7_pbpk_tm_kt.csv",
                "notes": f"{r['tmc_unit_note']} {r['kt_note']}".strip()})
    with open(os.path.join(OUT, "transporters.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TR_FIELDS)
        w.writeheader()
        w.writerows(out)
    return out


# --------------------------------------------------------- regulatory table
def build_regulatory():
    src = rows(os.path.join(DB, "regulatory_values.csv"))
    agency = rows(os.path.join(DB, "all_agency_clearance.csv"))
    fields = sorted({k for r in src + agency for k in r} | {"provenance"})
    out = []
    for r in src:
        r = dict(r); r["provenance"] = "regulatory_values.csv"; out.append(r)
    for r in agency:
        r = dict(r); r["provenance"] = "all_agency_clearance.csv"; out.append(r)
    with open(os.path.join(OUT, "regulatory.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(out)
    return out


if __name__ == "__main__":
    tk = build_tk()
    bd = build_binding()
    tr = build_transporters()
    rg = build_regulatory()
    print(f"combined/tk_parameters.csv   {len(tk):5d} rows")
    print(f"combined/binding.csv         {len(bd):5d} rows")
    print(f"combined/transporters.csv    {len(tr):5d} rows")
    print(f"combined/regulatory.csv      {len(rg):5d} rows")
    print(f"{'':29}{len(tk) + len(bd) + len(tr) + len(rg):5d} rows total")
