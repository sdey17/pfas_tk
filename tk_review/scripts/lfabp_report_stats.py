#!/usr/bin/env python3
"""Every number quoted in report/LFABP_BINDING.md, computed from the sheet.

Read-only. Run after build_lfabp_working_table.py so the report and the
workbook cannot drift apart.

    python3 -I scripts/lfabp_report_stats.py
"""
import collections, csv, os, re, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
XL = os.path.join(HERE, "..", "lfabp_logK_Kd.xlsx")
AL = os.path.join(HERE, "..", "albumin_logK_Kd_FINAL.xlsx")
DB = os.path.join(HERE, "..", "db", "protein_binding.csv")
L = "log K (L water / kg protein)"

# method family per source, from the primary
METHOD = {
    "Fan 2025": "fluorescent-probe displacement",
    "Sheng 2016": "fluorescent-probe displacement / ITC",
    "Sheng 2018": "fluorescent-probe displacement",
    "Woodcroft 2010": "fluorescent-probe displacement / ITC",
    "Yang 2020": "fluorescent-probe displacement",
    "Zhang 2013": "fluorescent-probe displacement",
    "Gao 2019": "equilibrium dialysis",
    "Khazaee 2021": "equilibrium dialysis",
    "Jia 2022": "microscale thermophoresis",
}
PROBE = {
    "Fan 2025": "1,8-ANS", "Sheng 2016": "1,8-ANS", "Sheng 2018": "1,8-ANS",
    "Woodcroft 2010": "ANS", "Yang 2020": "1,8-ANS",
    "Zhang 2013": "NBD-lauric acid",
}
NCF = {}   # compound -> n perfluorinated carbons, from the database
HEAD = set()  # compounds the database calls a plain carboxylate


def load(path, species):
    import openpyxl
    ws = openpyxl.load_workbook(path)["Coefficients"]
    hdr = [c.value for c in ws[1]][:11]
    out = []
    for r in ws.iter_rows(min_row=2, max_col=11):
        v = [c.value for c in r]
        if v[0] is None:
            continue
        d = dict(zip(hdr, v))
        if d["Species"] in species:
            out.append(d)
    return out


def q(vals):
    v = sorted(vals)
    n = len(v)
    return v[n // 4], statistics.median(v), v[(3 * n) // 4]


def main():
    R = load(XL, {"human", "rat", "mouse"})
    A = load(AL, {"human", "bovine"})
    for row in csv.DictReader(open(DB, encoding="utf-8")):
        if row["n_perfluorinated_carbons"].strip().isdigit():
            NCF.setdefault(row["chemical"].strip(),
                           int(row["n_perfluorinated_carbons"]))
            if row["head_group"].strip().lower().startswith("carboxylate"):
                HEAD.add(row["chemical"].strip())

    lk = [d[L] for d in R]
    lo, med, hi = q(lk)
    print("=" * 72)
    print(f"HEADLINE  {len(R)} rows | {len({d['Source'] for d in R})} sources | "
          f"{len({d['PFAS'] for d in R})} compounds")
    print(f"  log K {min(lk):.2f} to {max(lk):.2f}, median {med:.2f}, "
          f"IQR {lo:.2f}-{hi:.2f}")
    for sp in ("human", "rat"):
        v = [d[L] for d in R if d["Species"] == sp]
        print(f"  {sp:6s} n={len(v):3d} median {statistics.median(v):.2f} "
              f"({min(v):.2f} to {max(v):.2f}) from "
              f"{len({d['Source'] for d in R if d['Species']==sp})} sources")

    print("\nTABLE 1  composition by source")
    print(f"{'source':16s} {'n':>3s} {'species':7s} {'cpds':>4s} "
          f"{'median':>7s} {'range':>13s}  method")
    for s in sorted({d["Source"] for d in R}):
        sub = [d for d in R if d["Source"] == s]
        v = sorted(x[L] for x in sub)
        sp = "/".join(sorted({x["Species"] for x in sub}))
        print(f"{s:16s} {len(sub):3d} {sp:7s} {len({x['PFAS'] for x in sub}):4d} "
              f"{statistics.median(v):7.2f} {v[0]:6.2f}-{v[-1]:6.2f}  "
              f"{METHOD.get(s,'?')}")

    print("\nTABLE 2  route (basis of the stoichiometry used)")
    for rt, n in collections.Counter(d["Route"] for d in R).most_common():
        print(f"  {n:3d}  {rt}")

    print("\nTABLE 3  within-compound dispersion (compounds with >= 3 values)")
    by = collections.defaultdict(list)
    for d in R:
        by[d["PFAS"]].append(d)
    multi = {c: v for c, v in by.items() if len(v) >= 3}
    rng = []
    print(f"{'compound':10s} {'n':>2s} {'min':>6s} {'max':>6s} {'range':>6s}  sources")
    for c, v in sorted(multi.items(), key=lambda x: -len(x[1])):
        lv = sorted(x[L] for x in v)
        rng.append(lv[-1] - lv[0])
        print(f"{c:10s} {len(v):2d} {lv[0]:6.2f} {lv[-1]:6.2f} "
              f"{lv[-1]-lv[0]:6.2f}  {len({x['Source'] for x in v})}")
    if rng:
        print(f"  -> {len(multi)} compounds, median range "
              f"{statistics.median(rng):.2f} log, max {max(rng):.2f}")

    print("\nTABLE 4  by method family")
    fam = collections.defaultdict(list)
    for d in R:
        f = METHOD.get(d["Source"], "?")
        if "ITC" in str(d["Route"]):
            f = "isothermal titration calorimetry"
        elif "/" in f:
            f = f.split(" / ")[0]
        fam[f].append(d[L])
    for f, v in sorted(fam.items(), key=lambda x: -len(x[1])):
        print(f"  {len(v):3d}  {statistics.median(v):5.2f}  "
              f"({min(v):.2f} to {max(v):.2f})  {f}")

    print("\nTABLE 4b  by probe, displacement rows only")
    pr = collections.defaultdict(list)
    for d in R:
        p = PROBE.get(d["Source"])
        if p and "ITC" not in str(d["Route"]):
            pr[p].append(d[L])
    for p, v in sorted(pr.items(), key=lambda x: -len(x[1])):
        print(f"  {len(v):3d}  {statistics.median(v):5.2f}  "
              f"({min(v):.2f} to {max(v):.2f})  {p}")

    print("\nTABLE 5  PFCA homologous series (human, median per compound)")
    series = []
    for c, v in sorted(by.items()):
        n = NCF.get(c)
        hum = [x[L] for x in v if x["Species"] == "human"]
        if n and hum and c in HEAD:
            series.append((n, c, len(hum), statistics.median(hum)))
    series.sort()
    for n, c, k, m in series:
        print(f"  C{n:<3d} {c:8s} n={k:2d}  median log K {m:.2f}")
    rise = [(n, m) for n, c, k, m in series if n <= 9]
    if len(rise) > 2:
        xs = [x for x, _ in rise]; ys = [y for _, y in rise]
        mx, my = statistics.mean(xs), statistics.mean(ys)
        b = (sum((x-mx)*(y-my) for x, y in rise)
             / sum((x-mx)**2 for x in xs))
        print(f"  -> C3-C9 slope {b:+.2f} log per perfluorinated carbon "
              f"(n = {len(rise)} compounds)")

    print("\nL-FABP vs ALBUMIN, compounds measured against both")
    al = collections.defaultdict(list)
    for d in A:
        al[d["PFAS"]].append(d[L])
    shared = sorted(set(al) & set(by))
    diffs = []
    for c in shared:
        lf = statistics.median(x[L] for x in by[c])
        ab = statistics.median(al[c])
        diffs.append(lf - ab)
    if diffs:
        print(f"  {len(shared)} shared compounds; median L-FABP minus albumin "
              f"{statistics.median(diffs):+.2f} log")
        print(f"  albumin: {len(A)} rows, median {statistics.median([d[L] for d in A]):.2f}")
        print(f"  L-FABP : {len(R)} rows, median {med:.2f}")

    print("\nCONCENTRATION REPORTING")
    for c in ("PFAS conc (uM)", "Protein conc (uM)", "ligand:protein ratio"):
        f = sum(1 for d in R if d[c] not in (None, ""))
        print(f"  {f:3d}/{len(R)}  {c}")
    rat = sorted(d["ligand:protein ratio"] for d in R
                 if isinstance(d["ligand:protein ratio"], (int, float)))
    if rat:
        print(f"  ligand:protein ratio {min(rat):g} to {max(rat):g}, "
              f"median {statistics.median(rat):g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
