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

4. COMPUTES CL_int where Km and Vmax are both present on a common basis.

       CL_int [uL/min/mg protein] = 1000 * Vmax [nmol/min/mg] / Km [uM]

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

NOT_MEASURED = re.compile(
    r"IVIVE|predicted|mRNA|in vivo renal clearance|multiple regression|"
    r"correlation of|used in PBPK model|in vitro-in silico", re.I)
NOT_A_TRANSPORTER = re.compile(r"\bCES\d\b|enzyme, NOT a transporter", re.I)
COMPOSITE = re.compile(r"composite", re.I)
CARRIED_OVER = re.compile(r"NOT MEASURED HERE", re.I)
NOT_VMAX = re.compile(r"BW\^0\.75|clearance, not Vmax", re.I)
SECOND_HAND = re.compile(r"as quoted|as reported|cited (?:in|by)", re.I)
VMAX_CANON = {"nmol/mg protein/min": "nmol/min/mg protein"}
VMAX_SCALABLE = {"nmol/min/mg protein"}

# a compound token: PFOA, PFHxS, 6:2 FTS, HFPO-TA, F53B, N-EtFOSA, PFO3TDA ...
TOKEN = re.compile(r"^(?:[A-Za-z0-9]+:[0-9]+\s+)?[A-Za-z0-9][A-Za-z0-9:\-]*$")
UNSPECIFIED = re.compile(r"^\d+\s+PF|several|various|^any$", re.I)


def num(s):
    """First number in a cell like '185 +/- 97', '47 (15)', '1.29-7.92'."""
    m = re.search(r"-?\d+(?:\.\d+)?", (s or "").replace(",", ""))
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
        if vu and vu != "nmol/min/mg protein" and vmax is not None:
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
            d["vmax_nmol_min_mg"] = vmax if (
                vmax is not None and vu in VMAX_SCALABLE) else ""
            d["vmax_units"] = vu
            d["substrate_call"] = kind or r["substrate_or_inhibitor"]
            d["clint_uL_min_mg"] = (
                round(1000.0 * vmax / km, 1)
                if (km and vmax is not None and vu in VMAX_SCALABLE and km > 0)
                else "")
            if len(names) > 1:
                d["chemical_cell_original"] = r["chemical"]
            out.append(d)

    fields = list(rows[0].keys()) + [
        "chemical_cell_original", "chemical_cell_read_as", "shared_across",
        "km_uM", "vmax_nmol_min_mg", "clint_uL_min_mg", "substrate_call"]
    fields = list(dict.fromkeys(fields))
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for d in out:
            w.writerow({f: d.get(f, "") for f in fields})

    print(f"read    {len(rows):4d} rows from the record")
    print(f"dropped {sum(dropped.values()):4d}")
    for g, n in dropped.most_common():
        print(f"           {n:3d}  {g}")
    print(f"wrote   {len(out):4d} rows to {os.path.relpath(OUT, os.path.join(HERE,'..'))}")
    print()
    for f, n in flags.most_common():
        print(f"   flag {n:3d}  {f}")

    km = [d for d in out if d["km_uM"] != ""]
    cl = [d for d in out if d["clint_uL_min_mg"] != ""]
    print()
    print(f"K_m present      {len(km):3d} / {len(out)}  "
          f"({len({d['study'][:18] for d in km})} studies, "
          f"{len({d['chemical'] for d in km})} compounds, "
          f"{len({d['transporter'] for d in km})} transporters)")
    print(f"CL_int computed  {len(cl):3d} / {len(out)}")
    if cl:
        v = sorted(d["clint_uL_min_mg"] for d in cl)
        print(f"   CL_int {v[0]:g} to {v[-1]:g} uL/min/mg protein")
    return 0


if __name__ == "__main__":
    sys.exit(main())
