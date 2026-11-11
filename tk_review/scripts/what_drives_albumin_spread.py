#!/usr/bin/env python3
"""What actually drives the spread in albumin log K?  Rank the candidates.

THE HYPOTHESIS THIS KILLS.  Albumin coefficients for one compound scatter
by two to four log across the literature.  Binding is saturable, and Qin
measured the saturation directly at about -0.3 log per decade of free
concentration, so the obvious explanation is that the labs that disagree
were dosing at different loadings.  If true, the spread is an artefact of
the operating point and correcting for it would reconcile the field.

It is not true, and this script is why.  Four parts, each a test the
hypothesis has to survive:

  A  the cross-laboratory slope, per compound, against Qin's own
  B  how much scatter correcting for concentration actually removes
  C  whether METHOD and concentration are collinear -- fluorescence
     quenching runs at high micromolar, dialysis runs low, so a slope
     fitted across methods may be measuring the method and calling it
     saturation.  Tested by refitting inside one method family.
  D  a variance decomposition: compound, then species, method, and the
     laboratory itself, ranked by how much each removes.

WHAT IT FINDS.  The pooled slope matches Qin's almost exactly, which looks
like confirmation -- until C shows method medians separated by nearly five
decades of concentration, and the slope inside a method family falls to
about a third of Qin's.  D then puts the laboratory at roughly 45% of the
variance against 2% for species, with method adding nothing once the
laboratory is known.  Concentration is real and second order; the dominant
term is whose bench it was measured on.

The within-method cells are small (n = 4) over narrow concentration spans,
so C is underpowered and attenuates a true slope.  It is evidence that the
cross-laboratory slope is confounded, not proof that saturation is absent.
Qin's own within-laboratory -0.30, measured over four orders on one bench,
stands.

Run:  python3 -I tk_review/scripts/what_drives_albumin_spread.py
"""
import collections
import csv
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "notebook", "albumin_logK_Kd_revised.csv")
LOGK = "log K (L water / kg protein)"
CONC = "PFAS conc (uM)"

FAM = [(("quenching-plot exponent",), "fluorescence"),
       (("counted by ITC",), "ITC"),
       (("single-site model",), "single-site fit"),
       (("assumed by the reader",), "reader n=1"),
       (("measured",), "measured directly")]

# Qin's own slopes, his Table 1 Freundlich fits (slope = n_Fr - 1).
QIN = {"PFBA": -0.30, "PFHxA": -0.27, "PFHpA": -0.30, "PFOA": -0.31,
       "PFNA": -0.15, "HFPO-DA": -0.49, "PFHxS": -0.47, "PFOS": -0.38,
       "6:2 FTSA": -0.09, "PFOSA": -0.10}


def fam(route):
    for words, label in FAM:
        if any(w in route for w in words):
            return label
    return "other"


def load():
    out = []
    with open(CSV, newline="", encoding="utf-8") as fh:
        for d in csv.DictReader(fh):
            try:
                lk = float(d[LOGK])
            except (TypeError, ValueError):
                continue
            try:
                c = float(d[CONC])
            except (TypeError, ValueError):
                c = None
            out.append({"pfas": d["PFAS"], "logk": lk, "src": d["Source"],
                        "fam": fam(d["Route"]), "sp": d["Species"],
                        "x": math.log10(c) if c and c > 0 else None})
    return out


def fit(pairs):
    """slope, r2, n, residual SD -- None if it cannot be fitted."""
    if len(pairs) < 4:
        return None
    xs = [x for x, _ in pairs]
    ys = [y for _, y in pairs]
    mx, my = st.mean(xs), st.mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    res = [y - (a + b * x) for x, y in zip(xs, ys)]
    sdy, sdr = st.pstdev(ys), st.pstdev(res)
    return b, (1 - sdr ** 2 / sdy ** 2 if sdy else 0.0), len(pairs), sdy, sdr


def main():
    rows = load()
    withx = [r for r in rows if r["x"] is not None]
    print(f"{len(rows)} albumin coefficients, {len(withx)} with a "
          f"concentration\n")

    print("=" * 68)
    print("A + B.  Per compound: does concentration explain the spread?")
    print("=" * 68)
    by = collections.defaultdict(list)
    for r in withx:
        by[r["pfas"]].append((r["x"], r["logk"]))
    print("%-9s %4s   %7s %8s   %8s %9s %7s" % (
        "PFAS", "n", "slope", "Qin's", "SD raw", "SD resid", "r2"))
    slopes, raws, resids = [], [], []
    for p in sorted(by, key=lambda p: -len(by[p])):
        f = fit(by[p])
        if not f:
            continue
        b, r2, n, sdy, sdr = f
        q = f"{QIN[p]:+.2f}" if p in QIN else "   -"
        print("%-9s %4d   %+7.2f %8s   %8.2f %9.2f %7.2f"
              % (p, n, b, q, sdy, sdr, r2))
        slopes.append(b); raws.append(sdy); resids.append(sdr)
    print("\n  median slope %+.2f   (Qin's own median %+.2f) -- they agree"
          % (st.median(slopes), st.median(list(QIN.values()))))
    print("  median SD %.2f -> %.2f after correcting: only %.0f%% removed"
          % (st.median(raws), st.median(resids),
             100 * (1 - st.median(resids) / st.median(raws))))

    print()
    print("=" * 68)
    print("C.  Is that slope the METHOD in disguise?")
    print("=" * 68)
    byf = collections.defaultdict(list)
    for r in withx:
        byf[r["fam"]].append(r["x"])
    print("  concentration each method is run at:")
    for f in sorted(byf, key=lambda f: -st.median(byf[f])):
        v = byf[f]
        print("    %-18s n=%-4d median %9.3g uM" % (f, len(v),
                                                    10 ** st.median(v)))
    span = max(st.median(v) for v in byf.values()) - \
        min(st.median(v) for v in byf.values())
    print(f"  -> method medians span {span:.1f} decades of concentration.")
    print("     Method and concentration are badly collinear.\n")
    print("  slope refitted INSIDE one method family:")
    inside = []
    for f in sorted(byf):
        g = [(r["x"], r["logk"]) for r in withx if r["fam"] == f]
        res = fit(g)
        if not res:
            print("    %-18s n=%-4d cannot fit (no spread in x)"
                  % (f, len(g)))
            continue
        b, r2, n, _, _ = res
        inside.append(b)
        print("    %-18s n=%-4d slope %+.2f   r2 %.2f" % (f, n, b, r2))
    print("\n  within compound AND within method:")
    cell, out = collections.defaultdict(list), []
    for r in withx:
        cell[(r["pfas"], r["fam"])].append((r["x"], r["logk"]))
    for (p, f), g in sorted(cell.items()):
        res = fit(g)
        if res:
            out.append(res[0])
            print("    %-9s %-18s n=%-3d slope %+.2f   r2 %.2f"
                  % (p, f, res[2], res[0], res[1]))
    if out:
        print("\n  median of those %d slopes: %+.2f, against Qin's %+.2f"
              % (len(out), st.median(out), st.median(list(QIN.values()))))
        print("  Held at one method, most of the slope goes away.")
        print("  (n = 4 per cell over narrow spans: underpowered, so read")
        print("   this as 'confounded', not as 'saturation is absent'.)")

    print()
    print("=" * 68)
    print("D.  Variance decomposition: what removes the most scatter?")
    print("=" * 68)

    def within(key):
        g = collections.defaultdict(list)
        for r in rows:
            g[(r["pfas"], key(r))].append(r["logk"])
        v = n = 0
        for c in g.values():
            if len(c) >= 2:
                m = st.mean(c)
                v += sum((x - m) ** 2 for x in c)
                n += len(c) - 1
        return (math.sqrt(v / n) if n else float("nan")), n

    base, _ = within(lambda r: 0)
    print("  pooled SD inside a compound: %.2f log\n" % base)
    print("  %-32s %7s %9s %6s" % ("also holding fixed", "SD", "removed",
                                   "df"))
    for name, key in [("nothing", lambda r: 0),
                      ("species", lambda r: r["sp"]),
                      ("method family", lambda r: r["fam"]),
                      ("species + method", lambda r: (r["sp"], r["fam"])),
                      ("the laboratory", lambda r: r["src"]),
                      ("laboratory + method",
                       lambda r: (r["src"], r["fam"]))]:
        sd, n = within(key)
        print("  %-32s %7.2f %8.0f%% %6d"
              % (name, sd, 100 * (1 - sd / base), n))
    lab, _ = within(lambda r: r["src"])
    both, _ = within(lambda r: (r["src"], r["fam"]))
    print("\n  Species removes almost nothing.  The laboratory removes the")
    print("  most.  Method and laboratory overlap heavily -- knowing both")
    print("  (%.2f) beats knowing the laboratory alone (%.2f) by only %.2f"
          % (both, lab, lab - both))
    print("  log, so method is largely a label for which bench it was")
    print("  measured on rather than an independent term.")
    print("  Even knowing BOTH, %.2f log (%.0fx) remains inside one"
          % (both, 10 ** both))
    print("  compound, one laboratory and one method.")
    print()
    print("  bischel2010_refit_spread.py accounts for a good part of that")
    print("  last residue: refitting ONE of Bischel's datasets, with the")
    print("  bench and the operator held fixed, moves PFOA's log K by 0.55")
    print("  log on its own.")


if __name__ == "__main__":
    main()
