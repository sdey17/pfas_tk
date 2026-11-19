#!/usr/bin/env python3
"""Stage 2: join the kinetics and direction tables onto the localisation table.

    python3 -I scripts/join_transporter_direction.py

Read-only. Reports three things:

  1. COVERAGE  -- how many rows get a membrane and a net process, and which
     transporter strings do not resolve to an entry in the authority table.
  2. CONFLICTS -- rows whose paper-stated direction disagrees with the
     authority, and rows for a transporter the authority says does not
     transport PFAS at all.
  3. CL_int BY PROCESS -- the clearance column regrouped on the authority's
     direction rather than on the free text each paper happened to write.

TISSUE IS RESOLVED FROM THE PAPER, NOT GUESSED. Several transporters appear
in two tissues with opposite polarity -- OATP2B1 is apical in kidney and
sinusoidal in liver -- so the join needs to know which one a row is about.
It reads that from the paper's own words: a row whose stated direction or
expression system mentions hepatic, liver, sinusoidal or biliary is a liver
row, ileal or intestinal is an intestine row, placental is placenta, and
everything else is kidney, which is what the renal transporter papers are.
Where that rule picks a tissue the authority has no row for, the row is
reported as unresolved rather than silently assigned.
"""
import collections, csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOC = os.path.join(HERE, "..", "db", "transporter_localisation.csv")
KIN = os.path.join(HERE, "..", "db", "transporter_kinetics_working.csv")
DIR = os.path.join(HERE, "..", "db", "transporter_direction.csv")

# transporter string -> canonical name in the authority table. Case is not
# meaningful in these names (OAT1 and Oat1 are the same gene in two species,
# and the species column carries that), but the spelling variants are.
def canon(t):
    k = (t or "").strip().lower()
    k = k.replace("alpha", "a").replace("beta", "b").replace("/", "")
    k = re.sub(r"[\s\-_()]", "", k)
    table = {
        "oat1": "OAT1", "oat2": "OAT2", "oat3": "OAT3", "oat4": "OAT4",
        "urat1": "URAT1", "asbt": "ASBT", "ntcp": "NTCP", "bcrp": "BCRP",
        "pgp": "P-gp", "mrp2": "MRP2", "bsep": "BSEP", "oct2": "OCT2",
        "ostab": "OSTalpha/beta", "ostaostb": "OSTalpha/beta",
        "oatp1a1": "Oatp1a1", "oatp1a2": "OATP1A2", "oatp1a5": "Oatp1a5",
        "oatp1b1": "OATP1B1", "oatp1b3": "OATP1B3", "oatp2b1": "OATP2B1",
        "npt2napiiia": "Npt2 (NaPi-IIa)", "npt2napiia": "Npt2 (NaPi-IIa)",
        "npt2": "Npt2 (NaPi-IIa)",
    }
    return table.get(k)


LIVER = re.compile(r"hepatic|liver|sinusoidal|biliary|canalicular|hepatocyte",
                   re.I)
GUT = re.compile(r"ileum|ileal|intestin|enterocyte|enterohepatic", re.I)
PLACENTA = re.compile(r"placenta", re.I)


def tissue_of(row):
    blob = " ".join((row.get(f) or "") for f in
                    ("direction_reabsorption_or_secretion",
                     "expression_system", "measurement_type"))
    if PLACENTA.search(blob):
        return "placenta"
    if GUT.search(blob) and not LIVER.search(blob):
        return "ileum"
    if LIVER.search(blob):
        return "liver"
    return "kidney proximal tubule"


SAYS_REABS = re.compile(r"reabsorp", re.I)
SAYS_SECR = re.compile(r"secretion", re.I)


def main():
    with open(LOC, newline="", encoding="utf-8") as fh:
        loc = list(csv.DictReader(fh))
    # Both sides go through canon(). The authority deliberately spells the
    # rat entries Oat1/Oat3/Asbt and the human ones OAT1/OAT3/ASBT, which is
    # the convention the literature uses, so keying on the raw string would
    # miss every rodent row.
    auth = {}
    for r in loc:
        c = canon(r["transporter"]) or r["transporter"]
        auth[(c, r["species"], r["tissue"])] = r
    # a species-agnostic fallback key, for rows whose species cell is blank
    rows = []
    for path, kind in ((KIN, "kinetics"), (DIR, "direction")):
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                r["_kind"] = kind
                rows.append(r)

    resolved, unresolved = [], collections.Counter()
    for r in rows:
        c = canon(r["transporter"])
        sp = (r["species_orthologue"] or "").strip()
        t = tissue_of(r)
        hit = auth.get((c, sp, t)) if c else None
        if hit is None and c:
            # the same transporter in that species, any tissue, if unique
            cand = [v for k, v in auth.items() if k[0] == c and k[1] == sp]
            if len(cand) == 1:
                hit = cand[0]
        if hit is None:
            unresolved[f"{r['transporter'][:44]} [{sp}]"] += 1
        else:
            r["_auth"] = hit
            resolved.append(r)

    print("=" * 76)
    print(f"COVERAGE   {len(resolved)} / {len(rows)} rows carry a membrane "
          f"and a net process")
    print("=" * 76)
    if unresolved:
        print(f"\n{sum(unresolved.values())} rows did not resolve "
              f"({len(unresolved)} distinct strings):")
        for k, n in unresolved.most_common():
            print(f"   {n:4d}  {k}")
        print("   (composite or unidentified entities; no single membrane to "
              "assign)")

    print("\n" + "=" * 76)
    print("CONFLICTS")
    print("=" * 76)
    flipped, dead = [], collections.Counter()
    for r in resolved:
        a = r["_auth"]
        said = r["direction_reabsorption_or_secretion"] or ""
        proc = a["net_process"]
        if said:
            if SAYS_REABS.search(said) and "secretion" in proc:
                flipped.append((r, a, "paper says reabsorption, authority "
                                      "says secretion"))
            elif SAYS_SECR.search(said) and "reabsorption" in proc:
                flipped.append((r, a, "paper says secretion, authority says "
                                      "reabsorption"))
        if a["transports_pfas"] == "no" and r["substrate_call"] in (
                "substrate", "likely substrate"):
            dead[f"{a['transporter']} [{a['species']}]"] += 1
    if flipped:
        for r, a, why in flipped[:12]:
            print(f"   {a['transporter']:14s} [{a['species']}] "
                  f"{r['chemical'][:12]:12s} {r['study'].split()[0]:10s} {why}")
    else:
        print("   no direction disagreements between any paper and the "
              "authority table")
    if dead:
        print("\n   rows calling a compound a substrate of a transporter the "
              "authority marks as NOT transporting PFAS:")
        for k, n in dead.most_common():
            print(f"      {n:3d}  {k}")

    print("\n" + "=" * 76)
    print("CL_int REGROUPED ON THE AUTHORITY'S DIRECTION")
    print("=" * 76)
    import statistics
    g = collections.defaultdict(list)
    for r in resolved:
        v = r.get("clint_uL_min_mg") or ""
        if v:
            g[r["_auth"]["net_process"]].append((float(v), r))
    tot = sum(len(v) for v in g.values())
    print(f"{'net process':44s} {'n':>3s} {'median':>8s} {'range':>18s}")
    for proc, vals in sorted(g.items(), key=lambda x: -len(x[1])):
        v = sorted(x[0] for x in vals)
        print(f"{proc[:44]:44s} {len(v):3d} {statistics.median(v):8.1f} "
              f"{v[0]:7.1f} - {v[-1]:7.1f}")
    print(f"{'':44s} {tot:3d}  rows carry a clearance")
    print()
    print("  READ THE SECOND-STEP ROW WITH CARE. The six 'renal secretion")
    print("  (step 2 of 2)' rows are P-gp and BCRP from Niu 2026, measured in")
    print("  INSIDE-OUT MEMBRANE VESICLES, and that paper states outright that")
    print("  'the specific rates of PFAS uptake cannot be directly compared")
    print("  between these transporters' because the concentrations differed")
    print("  and the vesicle background (1000-2000 pmol/mg) far exceeds the")
    print("  cell background (100-400). They are NOT evidence that apical")
    print("  efflux is an order of magnitude faster than basolateral uptake.")
    cells = [x for proc, v in g.items() for x, r in v
             if "step 2" not in proc]
    if cells:
        print(f"  Excluding them, the clearance range across whole-cell")
        print(f"  assays is {min(cells):g} to {max(cells):g} uL/min/mg over "
              f"{len(cells)} rows.")

    print()
    e = collections.defaultdict(list)
    for r in resolved:
        v = r.get("clint_uL_min_mg") or ""
        if v:
            k = r["_auth"]["effect_on_half_life"].split(" (")[0]
            e[k].append(float(v))
    print(f"{'effect on half-life':24s} {'n':>3s} {'median':>8s}   "
          f"whole-cell assays only")
    e2 = collections.defaultdict(list)
    for r in resolved:
        v = r.get("clint_uL_min_mg") or ""
        if v and "step 2" not in r["_auth"]["net_process"]:
            e2[r["_auth"]["effect_on_half_life"].split(" (")[0]].append(
                float(v))
    for k, v in sorted(e.items(), key=lambda x: -len(x[1])):
        w = e2.get(k, [])
        tail = (f"{len(w):3d} {statistics.median(w):8.1f}" if w else "")
        print(f"{k:24s} {len(v):3d} {statistics.median(v):8.1f}   {tail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
