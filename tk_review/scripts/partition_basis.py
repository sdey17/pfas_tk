#!/usr/bin/env python3
"""Per kg of protein, or per litre of it? Classify every coefficient.

A protein-water partition coefficient can be reported on either of two
bases and they are not interchangeable:

    MASS    L of water per KG of protein     c_p = g/L / 1000
    VOLUME  L of water per L of protein      c_p = g/L * vbar

For serum albumin the partial specific volume vbar is 0.733 mL/g, so
the two differ by a constant

    log10(1 / 0.733) = 0.135 log

in every direction that matters: read a volume-basis coefficient as if
it were mass-basis and you over-predict binding by 1.36-fold.  That is
small next to the between-laboratory spread in this literature, which
is why it is easy to lose -- and it is systematic, which is why it
should not be.  It is the difference between 0.0054 and 0.0074 on
Bischel's PFOA fraction unbound, against a measured 0.009.

WHERE EACH ROW'S BASIS COMES FROM.

  the paper        a MEASURED coefficient carries whatever basis its
                   paper used.  Most say; the units strings in this
                   database were transcribed with that wording intact,
                   so the classifier reads them rather than guessing.

  this project     a CONVERTED coefficient was computed here from a Kd
                   through log K = log10(1e9/M) + log n - log Kd(uM).
                   M is a molar mass in g/mol, so the result is per KG
                   of protein: mass basis, by construction, whatever
                   the source paper used for anything else.

  nobody           some papers state neither, and two in this database
                   say so in as many words ("basis not stated in the
                   paper, so not converted").  Those rows carry a live
                   0.135-log ambiguity and the column says so instead
                   of picking one.

NOT EVERY COEFFICIENT HERE IS PROTEIN-WATER.  Phospholipid membrane
and whole-serum coefficients are reported per litre of the sorbing
phase and have no mass basis to convert to, because a membrane has no
molar mass.  They are classified as volume and flagged by phase, so
they are never handed to a per-kg calculation.

Run:  python3 tk_review/scripts/partition_basis.py
Writes db/partition_basis.csv
"""
import collections
import csv
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "db", "partition_basis.csv")

VBAR_ALBUMIN = 0.733e-3                      # L/g
BASIS_GAP = math.log10(1.0 / 0.733)          # 0.135 log

IS_PARTITION = re.compile(r"log ?[KD]_|partition|distribution coefficient",
                          re.I)
PER_KG = re.compile(r"per\s+kg|/\s*kg|kg[_ ]?(protein|albumin)", re.I)
PER_L = re.compile(r"per\s+L[_ ]?(water|protein|albumin|serum|lipid|"
                   r"prot\+lip|membrane)|volume[- ]of[- ]", re.I)
SAYS_MASS = re.compile(r"MASS basis", re.I)
SAYS_VOLUME = re.compile(r"VOLUME basis", re.I)
NOT_STATED = re.compile(r"basis not stated|BASIS NOT", re.I)
NON_PROTEIN = re.compile(r"membrane|lipid|IAM|serum|plasma|prot\+lip|"
                         r"medium|cell", re.I)

FIELDS = ["study", "year", "chemical", "protein", "species", "value",
          "units", "phase", "basis", "basis_evidence",
          "shift_if_misread"]


def classify(units):
    """(basis, evidence). Reads the paper's own wording first."""
    if SAYS_VOLUME.search(units):
        return "volume", "the paper's units say VOLUME basis"
    if SAYS_MASS.search(units):
        return "mass", "the paper's units say MASS basis"
    if NOT_STATED.search(units):
        return "not stated", "the paper does not state a basis"
    if PER_KG.search(units):
        return "mass", "units are per kg of protein"
    if PER_L.search(units):
        return "volume", "units are per litre of the sorbing phase"
    return "not stated", "no basis wording in the units"


def main():
    with open(DB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    out = []
    for r in rows:
        units = (r["ka_units"] or r["kd_units"]).strip()
        if not IS_PARTITION.search(units):
            continue
        val = (r["ka"] or r["kd"]).strip()
        basis, why = classify(units)
        phase = ("non-protein sorbing phase"
                 if NON_PROTEIN.search(units) else "protein")
        out.append(dict(
            study=r["study"], year=r["year"], chemical=r["chemical"],
            protein=r["protein"][:46], species=r["protein_species"],
            value=val, units=units[:110], phase=phase, basis=basis,
            basis_evidence=why,
            shift_if_misread=("" if basis == "not stated"
                              else f"{BASIS_GAP:.3f} log")))

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: (r["study"], r["year"],
                                               r["chemical"])))

    print(f"wrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")
    print(f"  {len(out)} reported partition / distribution coefficients\n")
    print(f"  the 0.135-log question: mass basis or volume basis?")
    c = collections.Counter(r["basis"] for r in out)
    for k in ("mass", "volume", "not stated"):
        if c[k]:
            print(f"    {c[k]:5d}  {k}")
    p = collections.Counter(r["phase"] for r in out)
    print(f"\n  what is doing the sorbing")
    for k, v in p.most_common():
        print(f"    {v:5d}  {k}")

    print(f"\n  rows whose basis nobody states -- a live {BASIS_GAP:.3f} "
          f"log either way")
    unk = [r for r in out if r["basis"] == "not stated"]
    for (s, y, u), n in collections.Counter(
            (r["study"], r["year"], r["units"][:62]) for r in unk
    ).most_common():
        print(f"    {n:5d}  {s} {y}")
        print(f"           {u}")

    print(f"\n  CONVERTED coefficients are not in this count. A value")
    print(f"  computed here from a Kd goes through log10(1e9/M) with M")
    print(f"  in g/mol, so it is per KG by construction -- mass basis,")
    print(f"  whatever the source paper used.")


if __name__ == "__main__":
    main()
