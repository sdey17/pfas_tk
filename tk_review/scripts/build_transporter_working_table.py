#!/usr/bin/env python3
"""Gate, normalise and split db/transporter_kinetics.csv into a working table.

    python3 -I scripts/build_transporter_working_table.py

db/transporter_kinetics.csv is the record: everything read from a primary goes
in, including rows that are not measured transporter kinetics at all. This
script is the filter, and it is the transporter counterpart of
build_lfabp_working_table.py. It never writes to the record.

It does four things and reports each one.

1. GATES. A row is dropped unless it is a measurement of a named transporter
   handling a named compound. Five gates, each counted below:

     NOT_MEASURED      IVIVE predictions, mRNA levels, in vivo clearances,
                       regression correlations, and Km values that are PBPK
                       model parameters rather than observations
     NOT_A_TRANSPORTER the carboxylesterase CES1/CES2 rows, which are in the
                       record because they were on a requested list
     COMPOSITE         "composite RPTEC proximal-tubule transport", which has
                       no single transporter identity to join a direction to
     CARRIED_OVER      a constant the paper prints but another paper measured
     NOT_VMAX          a whole-body PBPK Tm (mg/h/kg BW^0.75) or a clearance
                       mislabelled Vmax; the Km beside it is still kept
     SECOND_HAND       a constant this database holds only as another paper
                       quoted it. Yang 2010's OAT4 and URAT1 Vmax values are
                       read from Louisse 2023's introduction, not from Yang
                       2010 itself, and pairing one with a Km read from the
                       primary would build a clearance out of two papers. The
                       Vmax is dropped, the Km kept, and the row flagged. This
                       is why Yang 2010's own Vmax is on the request list

2. NORMALISES Vmax UNITS. `nmol/mg protein/min` and `nmol/min/mg protein` are
   the same basis written two ways; they become one spelling.
   `nmol/min/10^6 cells` is a different basis and is kept, flagged, and left
   out of the clearance column -- joining it needs mg protein per 10^6 cells,
   which none of these papers report.

3. SPLITS MULTI-COMPOUND CELLS. `chemical` sometimes holds a list
   ("PFBS, PFHxS, PFOS") or a Louisse-style classification
   ("substrates: PFHpA, PFOA | non-substrates: PFNA, ..."). Each becomes one
   row per compound, carrying a `shared_across` count so that a single
   observation spread over seven compounds is never read as seven
   observations. Cells that name no specific compound ("6 PFCAs, 3 PFSAs and
   4 fluorinated alternatives") cannot be split and are gated as UNSPECIFIED.

4. ASSIGNS AN EVIDENCE TIER (Stage 1's inclusion rule) and computes CL_int.

     A  reaches intrinsic clearance per mg protein. Either Km and Vmax are
        both reported on a convertible basis, or the paper reports Vmax/Km
        itself in uL/min/mg protein (Kimura 2017 does, and only that).
     B  reports a clearance on a basis that cannot be joined to A. Han 2008's
        uL/min per 10^6 hepatocytes is the only case; joining it needs mg
        protein per 10^6 cells, which that paper does not give.
     C  affinity only: a Km and no clearance. Still answers the saturation
        question in the report, which does not need Vmax -- half-saturation
        against serum concentration is a statement about Km alone.
     D  no constant: a substrate call, an uptake ratio, or an inhibition
        measurement. These are the direction evidence, written to a separate
        file, and they never enter a clearance or affinity comparison.

   WHY IC50 AND Ki ARE TIER D AND NOT TIER C. An IC50 is a property of the
   probe substrate and its concentration as much as of the inhibitor, so
   pooling IC50 values measured against different probes compares nothing.
   Converting one to a probe-independent Ki needs the probe's own Km
   (Cheng-Prusoff), and NOT ONE of the 20 numeric IC50 values in this record
   carries it. They stay as evidence that a compound interacts with a
   transporter, which is a direction claim, not an affinity one.

       CL_int [uL/min/mg protein] = 1000 * Vmax [nmol/min/mg] / Km [uM]

   after converting Vmax to nmol/min/mg; pmol/mg protein/min, which Nakagawa
   2008 and Niu 2026 use, is the same basis a thousand times finer.

   because nmol/umol = 1e-3 and 1 L = 1e6 uL. Louisse 2024 prints its own
   efficiency column and this reproduces it (3.4/159 -> 21.4 against a printed
   21), which is the check that the conversion is right.

WHAT THIS DOES NOT DO. It does not decide the inclusion rule -- that is Stage
1, and a row with a Km and no Vmax is still here. It does not assert the
direction; `direction_reabsorption_or_secretion` is carried as the paper gave
it and still has to be reconciled against db/transporter_mechanism.csv, which
is Stage 2.
"""
import collections, csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "db", "transporter_kinetics.csv")
OUT = os.path.join(HERE, "..", "db", "transporter_kinetics_working.csv")
DIR_OUT = os.path.join(HERE, "..", "db", "transporter_direction.csv")

NOT_MEASURED = re.compile(
    r"IVIVE|predicted|mRNA|in vivo renal clearance|multiple regression|"
    r"correlation of|used in PBPK model|in vitro-in silico", re.I)
NOT_A_TRANSPORTER = re.compile(r"\bCES\d\b|enzyme, NOT a transporter", re.I)
COMPOSITE = re.compile(r"composite", re.I)
CARRIED_OVER = re.compile(r"NOT MEASURED HERE", re.I)
NOT_VMAX = re.compile(r"BW\^0\.75|clearance, not Vmax", re.I)
SECOND_HAND = re.compile(r"as quoted|as reported|cited (?:in|by)", re.I)
# A clearance the paper states outright instead of a Vmax. Kimura 2017 prints
# "uptake clearance Vmax/Km = 55.0 uL/mg protein/min"; Han 2008 prints
# "uptake clearance 64.8 +/- 15.7 uL/min/10^6 cells". Matched narrowly, on the
# vmax cell and the source note only, never across the whole row.
DIRECT_CL = re.compile(
    r"(?:uptake\s+)?clearance\s*(?:Vmax/Km\s*)?=?\s*([0-9]+(?:\.[0-9]+)?)"
    r"(?:\s*\+/-\s*[0-9]+(?:\.[0-9]+)?)?\s*"
    r"u?L/(?:mg\s+protein/min|min/mg\s+protein|min/10\^6\s+cells)", re.I)
PER_CELLS = re.compile(r"10\^6\s+cells", re.I)
VMAX_CANON = {"nmol/mg protein/min": "nmol/min/mg protein",
              "pmol/min/mg protein": "pmol/mg protein/min"}
# unit -> factor converting it to nmol/min/mg protein. pmol is a thousandth
# of a nmol; Nakagawa 2008 and Niu 2026 both report in pmol/mg protein/min,
# which is the same basis as Louisse's nmol/min/mg, a thousand times finer.
VMAX_TO_NMOL = {"nmol/min/mg protein": 1.0, "pmol/mg protein/min": 0.001}
VMAX_SCALABLE = set(VMAX_TO_NMOL)

# a compound token: PFOA, PFHxS, 6:2 FTS, HFPO-TA, F53B, N-EtFOSA, PFO3TDA ...
TOKEN = re.compile(r"^(?:[A-Za-z0-9]+:[0-9]+\s+)?[A-Za-z0-9][A-Za-z0-9:\-]*$")
UNSPECIFIED = re.compile(r"^\d+\s+PF|several|various|^any$", re.I)


def num(s):
    """The number a cell STARTS with: '185 +/- 97', '47 (15)', '1.29-7.92'.

    It must start with one, after an optional < > ~ or =. Cells that open
    with prose are placeholders, not measurements, and a loose search finds
    numbers inside them that were never constants: 'within 1.41-2.78
    (individual value only in Fig 2b)' yielded a Km of 1.41, and 'not
    extractable from text (reported only in Fig 1f)' yielded a Km of 1, from
    the figure number. Ten such phantom Km values were in the working table
    before this guard, four of them multiplied out by a multi-compound split.
    """
    t = (s or "").replace(",", "").strip().lstrip("<>~=").lstrip()
    m = re.match(r"-?\d+(?:\.\d+)?", t)
    return float(m.group()) if m else None


def compounds(cell):
    """Every compound named in a `chemical` cell, with how the cell was read."""
    c = (cell or "").strip()
    if not c:
        return [], "empty"
    if UNSPECIFIED.match(c):
        return [], "unspecified"
    # Louisse style: "substrates: A, B | non-substrates: C, D"
    if "substrates:" in c:
        out = []
        for part in c.split("|"):
            kind = ("not a substrate" if "non-substrate" in part.lower()
                    else "substrate")
            names = part.split(":", 1)[1] if ":" in part else ""
            for n in names.split(","):
                n = n.strip()
                if n and n.lower() != "none" and TOKEN.match(n):
                    out.append((n, kind))
        return out, "classification"
    # "PFOA (C8) and PFNA (C9) transported; PFHpA (C7) not"
    if ";" in c or " transported" in c:
        out = []
        for part in re.split(r";", c):
            kind = ("not a substrate" if re.search(r"\bnot\b", part)
                    else "substrate")
            for n in re.findall(r"\b(PF[A-Za-z0-9]+|\d:\d\s?FT[A-Za-z]*)\b",
                                part):
                out.append((n, kind))
        return out, "prose"
    # plain list, with parenthetical carbon ranges stripped
    bare = re.sub(r"\([^)]*\)", "", c).strip()
    names = [n.strip() for n in bare.split(",") if n.strip()]
    if len(names) > 1 and all(TOKEN.match(n) for n in names):
        return [(n, None) for n in names], "list"
    if TOKEN.match(bare):
        return [(bare, None)], "single"
    return [(c, None)], "single"


def main():
    with open(SRC, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    dropped = collections.Counter()
    flags = collections.Counter()
    out = []

    for r in rows:
        note = " ".join((r.get(k) or "") for k in
                        ("assay_volume_note", "source_table", "measurement_type"))
        mt = r.get("measurement_type") or ""
        if NOT_A_TRANSPORTER.search(r["transporter"] + note):
            dropped["NOT_A_TRANSPORTER"] += 1; continue
        if COMPOSITE.search(r["transporter"]):
            dropped["COMPOSITE"] += 1; continue
        if CARRIED_OVER.search(note):
            dropped["CARRIED_OVER"] += 1; continue
        if NOT_MEASURED.search(mt):
            dropped["NOT_MEASURED"] += 1; continue

        km, vmax = num(r["km"]), num(r["vmax"])
        vu = (r["vmax_units"] or "").strip()
        if vmax is not None and SECOND_HAND.search(r["vmax"]):
            if km is None:
                dropped["SECOND_HAND"] += 1; continue
            flags["Vmax discarded (second-hand, quoted by another paper), "
                  "Km kept"] += 1
            vmax, vu = None, ""
        if num(r["ic50"]) is not None and SECOND_HAND.search(r["source_table"]):
            flags["Ki/IC50 held only as another paper quotes it"] += 1
        if vmax is not None and NOT_VMAX.search(vu):
            if km is None:
                dropped["NOT_VMAX"] += 1; continue
            flags["Vmax discarded (not an in vitro rate), Km kept"] += 1
            vmax, vu = None, ""
        vu = VMAX_CANON.get(vu, vu)
        if vu and vu not in VMAX_SCALABLE and vmax is not None:
            flags[f"Vmax on a basis that cannot be joined: {vu}"] += 1

        names, how = compounds(r["chemical"])
        if not names:
            dropped["UNSPECIFIED" if how == "unspecified" else "NO_COMPOUND"] += 1
            continue
        if len(names) > 1:
            flags[f"multi-compound cell split ({how})"] += 1

        for name, kind in names:
            d = dict(r)
            d["chemical"] = name
            d["shared_across"] = len(names)
            d["chemical_cell_read_as"] = how
            d["km_uM"] = km if km is not None else ""
            vmax_n = (vmax * VMAX_TO_NMOL[vu]
                      if (vmax is not None and vu in VMAX_SCALABLE) else None)
            d["vmax_nmol_min_mg"] = (round(vmax_n, 6)
                                     if vmax_n is not None else "")
            d["vmax_units"] = vu
            d["substrate_call"] = kind or r["substrate_or_inhibitor"]
            clint = (round(1000.0 * vmax_n / km, 1)
                     if (km and vmax_n is not None and km > 0) else "")
            # a clearance the paper states outright, where we have no Vmax
            direct, per_cells = "", False
            if clint == "":
                for cell in (r["vmax"], r["source_table"]):
                    m = DIRECT_CL.search(cell or "")
                    if m:
                        direct = float(m.group(1))
                        per_cells = bool(PER_CELLS.search(m.group(0)))
                        break
            if direct != "" and not per_cells:
                clint = direct
            d["clint_uL_min_mg"] = clint
            d["clint_uL_min_10e6cells"] = direct if per_cells else ""
            d["tier"] = ("A" if clint != "" else
                         "B" if per_cells else
                         "C" if km else "D")
            if len(names) > 1:
                d["chemical_cell_original"] = r["chemical"]
            out.append(d)

    fields = list(rows[0].keys()) + [
        "chemical_cell_original", "chemical_cell_read_as", "shared_across",
        "km_uM", "vmax_nmol_min_mg", "clint_uL_min_mg",
        "clint_uL_min_10e6cells", "tier", "substrate_call"]
    fields = list(dict.fromkeys(fields))
    # Stage 1's inclusion rule, made physical: constants and direction
    # evidence are different quantities and live in different files. A row
    # that carries no constant cannot enter a clearance or affinity
    # comparison, and keeping it in the same sheet is how it eventually does.
    kin = [d for d in out if d["tier"] in ("A", "B", "C")]
    dirn = [d for d in out if d["tier"] == "D"]
    for path, rowset in ((OUT, kin), (DIR_OUT, dirn)):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for d in rowset:
                w.writerow({f: d.get(f, "") for f in fields})

    print(f"read    {len(rows):4d} rows from the record")
    print(f"dropped {sum(dropped.values()):4d}")
    for g, n in dropped.most_common():
        print(f"           {n:3d}  {g}")
    print(f"kept    {len(out):4d} rows, split by evidence tier")
    print()
    for f, n in flags.most_common():
        print(f"   flag {n:3d}  {f}")

    tiers = collections.Counter(d["tier"] for d in out)
    cl = [d for d in out if d["clint_uL_min_mg"] != ""]
    km = [d for d in out if d["km_uM"] != ""]
    rel = lambda p: os.path.relpath(p, os.path.join(HERE, ".."))
    print()
    print(f"  {rel(OUT)}   {len(kin)} rows")
    print(f"    tier A  {tiers['A']:3d}  CL_int per mg protein")
    print(f"    tier B  {tiers['B']:3d}  clearance on another basis "
          f"(per 10^6 cells)")
    print(f"    tier C  {tiers['C']:3d}  affinity only (K_m, no clearance)")
    print(f"  {rel(DIR_OUT)}        {len(dirn)} rows")
    print(f"    tier D  {tiers['D']:3d}  substrate calls, uptake ratios, "
          f"inhibition")
    print()
    print(f"K_m present      {len(km):3d}  "
          f"({len({d['study'][:18] for d in km})} studies, "
          f"{len({d['chemical'] for d in km})} compounds, "
          f"{len({d['transporter'] for d in km})} transporters)")
    if cl:
        v = sorted(d["clint_uL_min_mg"] for d in cl)
        print(f"CL_int           {len(cl):3d}  {v[0]:g} to {v[-1]:g} "
              f"uL/min/mg protein, "
              f"{len({d['study'][:18] for d in cl})} studies")
    return 0


if __name__ == "__main__":
    sys.exit(main())
