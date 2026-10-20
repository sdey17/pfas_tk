#!/usr/bin/env python3
"""Where is the albumin chain-length optimum for perfluorocarboxylates?

Fischer 2024 puts the maximum at 10 perfluorinated carbons (PFUnDA).
Chen 2025 puts it at 7 (PFOA).  Bischel 2011 also says 7.  Starnes 2025
says 6.  Gao 2025 and Qin 2024 never see a maximum at all.  This script
sorts out which of those are measurements of the optimum and which are
measurements of where a method stops working.

Series are read from db/protein_binding.csv except where a study only
supplies a summary range, which is hard-coded with its source noted.

Run:  python3 scripts/albumin_chain_length_optimum.py
      python3 scripts/albumin_chain_length_optimum.py --figure
"""
import argparse
import collections
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")
OUT = os.path.join(HERE, "..", "figures")

# Chen 2025 SI Table S5 also reports a membrane-water partition coefficient
# for every compound, measured in the same lab on the same chemicals.  It is
# the control for "does this compound behave normally at long chain length".
CHEN_KMW = {3: 1.63, 4: 2.02, 5: 2.42, 6: 2.85, 7: 3.28, 8: 3.75,
            9: 4.18, 10: 4.50, 11: 4.82, 12: 5.11, 13: 5.36, 15: 5.73}

# Bischel et al. 2011, Environ Toxicol Chem 30:2423-2430, abstract: BSA,
# C4-C12 carboxylates (nCF 3-11), equilibrium dialysis + nanoESI-MS +
# fluorescence at low PFAA:albumin mole ratios, log K_protein-water 3.3-4.3,
# with affinity "decreasing from the C8 to C12 PFCAs" -- so a maximum at
# C8 = nCF 7, and a total span of 1.0 log unit over 8 CF2 groups.
NAME = {3: "PFBA", 4: "PFPeA", 5: "PFHxA", 6: "PFHpA", 7: "PFOA", 8: "PFNA",
        9: "PFDA", 10: "PFUnDA", 11: "PFDoDA", 12: "PFTrDA", 13: "PFTeDA",
        15: "PFHxDA"}

BISCHEL11 = dict(span_log=1.0, n_cf2=8, peak_ncf=7,
                 note="abstract only; no per-compound table in hand")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load():
    with open(CSV) as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def series(rows, study, year, species=None, protein_sub=None, head="carboxylate",
           field="ka", agg="max", invert=False, islog=False):
    """Collapse to one value per nCF, on a log scale where bigger is tighter.

    invert=True for K_D-like fields where smaller is tighter; islog=True when
    the stored value is already a logarithm (Fischer and Qin report log D, not
    an association constant); agg='median' for studies with replicate rows."""
    acc = collections.defaultdict(list)
    for r in rows:
        if r["study"] != study or r["year"] != year:
            continue
        if r["head_group"] != head:
            continue
        if species is not None and r["protein_species"] != species:
            continue
        if protein_sub is not None and protein_sub not in r["protein"]:
            continue
        v, n = num(r[field]), num(r["n_perfluorinated_carbons"])
        if v is None or n is None or v <= 0:
            continue
        acc[int(n)].append(v)
    out = {}
    for n, vs in acc.items():
        vs.sort()
        v = (vs[len(vs) // 2] if agg == "median"
             else (min(vs) if invert else max(vs)))
        if islog:
            out[n] = v
        else:
            out[n] = -math.log10(v) if invert else math.log10(v)
    return out


# ---------------------------------------------------------------------------
# Shape statistics
# ---------------------------------------------------------------------------

def peak(s):
    return max(s, key=lambda n: s[n])


def interior(s):
    """True when the maximum has at least one measured point on each side."""
    p = peak(s)
    return min(s) < p < max(s)


def sign_changes(s):
    """Number of sign changes in the first difference. A single-peaked curve
    has exactly 1 (up then down); 0 means still rising or still falling."""
    ks = sorted(s)
    d = [s[b] - s[a] for a, b in zip(ks, ks[1:])]
    sg = [1 if x > 0 else -1 for x in d if abs(x) > 1e-12]
    return sum(1 for a, b in zip(sg, sg[1:]) if a != b)


def limb_slope(s, lo=4, hi=7):
    """log-units per CF2 on the common rising limb nCF lo..hi."""
    pts = [(n, s[n]) for n in sorted(s) if lo <= n <= hi]
    if len(pts) < 3:
        pts = [(n, s[n]) for n in sorted(s) if n <= hi]
    if len(pts) < 2:
        return None, len(pts)
    mx = sum(n for n, _ in pts) / len(pts)
    my = sum(v for _, v in pts) / len(pts)
    sxy = sum((n - mx) * (v - my) for n, v in pts)
    sxx = sum((n - mx) ** 2 for n, _ in pts)
    return (sxy / sxx if sxx else None), len(pts)


def own_limb(s):
    """log-units per CF2 from the shortest compound to one CF2 below the
    series' own peak -- each series' own rising limb, not a fixed window.
    A fixed window misreads the curved limbs (globulin accelerates toward its
    peak, Fischer's albumins decelerate), which is why this is the fairer
    statistic even though it uses a different window for each series."""
    p = peak(s)
    pts = [(n, s[n]) for n in sorted(s) if n <= p - 1]
    if len(pts) < 2:
        return None, 0
    (n0, v0), (n1, v1) = pts[0], pts[-1]
    return (v1 - v0) / (n1 - n0), len(pts)


def window_span(s, lo=4, hi=9):
    """Dynamic range over a window every series covers, so the comparison is
    like for like rather than a function of how far each series reaches."""
    vs = [s[n] for n in s if lo <= n <= hi]
    return (max(vs) - min(vs), len(vs)) if len(vs) >= 3 else (None, len(vs))


def spearman(pairs):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    xs = rank([a for a, _ in pairs])
    ys = rank([b for _, b in pairs])
    mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs); syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")



def vertex(s, half=3):
    """Quadratic vertex through the peak and +/- half CF2 either side."""
    p = peak(s)
    pts = [(n, s[n]) for n in sorted(s) if abs(n - p) <= half]
    if len(pts) < 3:
        return None
    n = len(pts)
    sx = sum(x for x, _ in pts); sx2 = sum(x * x for x, _ in pts)
    sx3 = sum(x ** 3 for x, _ in pts); sx4 = sum(x ** 4 for x, _ in pts)
    sy = sum(y for _, y in pts); sxy = sum(x * y for x, y in pts)
    sx2y = sum(x * x * y for x, y in pts)
    m = [[sx4, sx3, sx2, sx2y], [sx3, sx2, sx, sxy], [sx2, sx, n, sy]]
    for i in range(3):                               # Gaussian elimination
        pv = max(range(i, 3), key=lambda r: abs(m[r][i]))
        if abs(m[pv][i]) < 1e-12:
            return None
        m[i], m[pv] = m[pv], m[i]
        for r in range(3):
            if r == i:
                continue
            f = m[r][i] / m[i][i]
            for c in range(i, 4):
                m[r][c] -= f * m[i][c]
    a, b = m[0][3] / m[0][0], m[1][3] / m[1][1]
    return None if a >= 0 else -b / (2 * a)


def fmt(s):
    return "  ".join(f"{n}:{s[n]:.2f}" for n in sorted(s))


def banner(t):
    print()
    print(t)
    print("-" * len(t))


# ---------------------------------------------------------------------------

def build(rows):
    """Every albumin PFCA chain-length series in the file, on a log scale
    where bigger always means tighter."""
    S = {}
    S["Fischer 2024 HSA"] = dict(
        s=series(rows, "Fischer et al.", "2024", "human", "albumin", islog=True),
        method="SPME depletion", units="log D (L_w/kg_protein)",
        conc="0.2-42 ug/L, PFAS:protein <= 0.004")
    S["Fischer 2024 BSA"] = dict(
        s=series(rows, "Fischer et al.", "2024", "bovine", "albumin", islog=True),
        method="SPME depletion", units="log D",
        conc="0.2-42 ug/L, PFAS:protein <= 0.004")
    S["Fischer 2024 globulin"] = dict(
        s=series(rows, "Fischer et al.", "2024", protein_sub="globulin",
                 islog=True),
        method="SPME depletion", units="log D",
        conc="0.2-42 ug/L, PFAS:protein <= 0.009")
    S["Fischer 2024 whole serum"] = dict(
        s=series(rows, "Fischer et al.", "2024",
                 protein_sub="whole human serum", islog=True),
        method="SPME depletion", units="log D", conc="0.2-42 ug/L")
    S["Allendorf 2019 BSA"] = dict(
        s=series(rows, "Allendorf et al.", "2019", protein_sub="A3803",
                 islog=True),
        method="equilibrium dialysis, 37 C", units="log D (L_w/kg_protein)",
        conc="compound:BSA molar ratio below 0.1")
    S["Bischel 2011 BSA"] = dict(
        s=series(rows, "Bischel et al.", "2011", islog=True),
        method="equilibrium dialysis", units="log K_PW (volume basis)",
        conc="BSA 200 uM, PFAA 0.009-0.68 uM, ratio <= 0.0034")
    S["Chen 2025 HSA"] = dict(
        s=series(rows, "Chen et al.", "2025", field="kd", invert=True),
        method="equilibrium dialysis", units="-log K_D (M)",
        conc="0.0625-2 uM, HSA 1.0 uM")
    S["Gao 2025 BSA"] = dict(
        s=series(rows, "Gao et al.", "2025", agg="median"),
        method="Trp fluorescence quenching", units="log Ka_Hill",
        conc="uM-range titration")
    S["Qin 2024 BSA specific"] = dict(
        s={n: v for n, v in _qin(rows, "specific").items()},
        method="SPME, specific term", units="log D", conc="4 decades of isotherm")
    S["Qin 2024 BSA nonspecific"] = dict(
        s={n: v for n, v in _qin(rows, "nonspecific").items()},
        method="SPME, nonspecific term", units="log D", conc="4 decades of isotherm")
    S["Starnes 2025 HSA"] = dict(
        s=series(rows, "Starnes et al.", "2025"),
        method="DSF", units="log Ka", conc="uM-mM dye assay")
    S["Starnes 2024 HSA"] = dict(
        s=series(rows, "Starnes et al.", "2024", "human", field="kd", invert=True),
        method="DSF", units="-log K_D (mM)", conc="uM-mM dye assay")
    S["Alesio 2024 BSA"] = dict(
        s=series(rows, "Alesio & Bothun", "2024"),
        method="DSF", units="log Ka", conc="uM-mM dye assay")
    return {k: v for k, v in S.items() if len(v["s"]) >= 3}


def _qin(rows, which):
    out = {}
    for r in rows:
        if r["study"] != "Qin et al." or r["year"] != "2024":
            continue
        # "nonspecific" contains "specific", so match the full column label
        if r["head_group"] != "carboxylate":
            continue
        if f", {which} column" not in r["source_table"]:
            continue
        v, n = num(r["ka"]), num(r["n_perfluorinated_carbons"])
        if v is not None and n is not None:
            out[int(n)] = v              # already log D
    return out


# ---------------------------------------------------------------------------

def ryu_slope(lo=4, hi=13):
    """Mean log10 f_unbound across Ryu 2024's tissues and species, regressed
    on nCF. Negative because a falling unbound fraction means tighter binding.
    Computed rather than quoted so the anchor cannot drift."""
    acc = collections.defaultdict(list)
    for r in load():
        if not r["study"].startswith("Ryu") or r["head_group"] != "carboxylate":
            continue
        fu, n = num(r["unbound_fraction"]), num(r["n_perfluorinated_carbons"])
        if fu and n and fu > 0 and lo <= n <= hi:
            acc[int(n)].append(math.log10(fu))
    pts = [(n, sum(v) / len(v)) for n, v in sorted(acc.items())]
    mx = sum(a for a, _ in pts) / len(pts)
    my = sum(b for _, b in pts) / len(pts)
    sxy = sum((a - mx) * (b - my) for a, b in pts)
    sxx = sum((a - mx) ** 2 for a, _ in pts)
    return sxy / sxx


def section1(S):
    banner("1. Two different questions, often run together")
    print("'Does albumin binding rise with chain length?' and 'where does it")
    print("stop rising?' are not the same question, and most of the datasets")
    print("in this file answer only the first. A series whose largest value")
    print("sits at its longest compound has not located an optimum; it has")
    print("placed a lower bound on one.")
    print()
    print("So the first cut is purely structural: does the maximum have a")
    print("measured point on BOTH sides of it?")
    print()
    print(f"{'dataset':26s} {'nCF range':>11s} {'pts':>4s} {'peak':>5s} "
          f"{'interior?':>10s} {'sign chg':>9s}")
    for k, v in S.items():
        s = v["s"]
        print(f"{k:26s} {f'{min(s)}-{max(s)}':>11s} {len(s):4d} "
              f"{peak(s):5d} {('YES' if interior(s) else 'no'):>10s} "
              f"{sign_changes(s):9d}")
    print()
    print("'sign chg' counts reversals in the first difference. A clean")
    print("single-peaked curve gives exactly 1. A 0 means the series is still")
    print("climbing when it runs out. Anything above 1 means the series")
    print("wanders, and a binding curve does not wander -- the measurement does.")
    print()
    interiors = [k for k, v in S.items() if interior(v["s"])]
    print(f"{len(interiors)} of {len(S)} datasets bracket a maximum at all:")
    for k in interiors:
        s = S[k]["s"]
        print(f"  {k:26s} peak at nCF {peak(s):2d}, {sign_changes(s)} sign "
              f"change(s), {len(s)} points, span "
              f"{max(s.values()) - min(s.values()):.2f} log units")
    print()
    print("One of those gets struck out immediately. Alesio 2024 has three")
    print("points (nCF 7, 8, 9), so one sign change is automatic rather than")
    print("informative, and its PFDA sits 293x below its PFNA -- no binding")
    print("curve drops two and a half orders of magnitude for one CF2. It is")
    print("dropped from the adjudication and kept in the table only so the")
    print("exclusion is on the record.")
    print()
    print("(Bischel 2011 and Allendorf 2019 are now per-compound series in")
    print(" the table above; the abstract-only summary below is superseded.)")
    print("Bischel et al. 2011 (BSA, C4-C12 carboxylates, equilibrium dialysis")
    print("+ nanoESI-MS + fluorescence at low PFAA:albumin mole ratios) adds a")
    print(f"further bracketed maximum, from its abstract only: peak at nCF "
          f"{BISCHEL11['peak_ncf']},")
    print(f"total span {BISCHEL11['span_log']:.1f} log units over "
          f"{BISCHEL11['n_cf2']} CF2 groups.")
    print()
    print("So the field splits three ways:")
    print("  peak at nCF 10-11   Fischer (4 sorbents, one method)")
    print("  peak at nCF 6-8     Chen 2025, Bischel 2011, Starnes 2025")
    print("  no peak reached     Gao 2025, Qin 2024 (both terms), Starnes 2024")
    print()
    print("The third group is not neutral. Gao is still climbing at nCF 9 and")
    print("Qin's nonspecific term at nCF 10, which is incompatible with a")
    print("maximum at 6-8: both would have to be falling by then, and neither")
    print("is. Qin's specific term stops at nCF 7, so it only says the")
    print("optimum is at 7 or beyond.")



def section2(S):
    banner("2. The compression cannot be the explanation")
    print("Section 11b of the albumin primer established the compression:")
    print("the shortfall (Fischer - Chen) grows with Fischer at slope 0.82,")
    print("r2 = 0.90 on the carboxylates, which is the same as saying Chen")
    print("reads about 0.18 of any true affinity difference. Four fifths of")
    print("the dynamic range is lost. That is a real finding, and it is the")
    print("wrong tool for this job.")
    print()
    print("An optimum is an argmax, and argmax survives any strictly")
    print("increasing transform. If Chen were exactly a + b*Fischer with")
    print("b > 0, Chen's peak would sit at Fischer's nCF, whatever b is.")
    print()
    f = S["Fischer 2024 HSA"]["s"]
    c = S["Chen 2025 HSA"]["s"]
    for b in (1.0, 0.82, 0.5, 0.25):
        t = {n: 1.0 + b * v for n, v in f.items()}
        print(f"  Fischer rescaled by b = {b:4.2f} -> peak still at nCF {peak(t):2d}")
    print()
    print(f"  Fischer HSA actual peak  nCF {peak(f)}")
    print(f"  Chen    HSA actual peak  nCF {peak(c)}")
    print()
    print("So the two disagree in SHAPE, not in scale, and the residuals from")
    print("the Chen-on-Fischer line have to carry chain-length structure.")
    print("They do:")
    print()
    common = sorted(set(f) & set(c))
    mx = sum(f[n] for n in common) / len(common)
    my = sum(c[n] for n in common) / len(common)
    sxy = sum((f[n] - mx) * (c[n] - my) for n in common)
    sxx = sum((f[n] - mx) ** 2 for n in common)
    syy = sum((c[n] - my) ** 2 for n in common)
    sl = sxy / sxx
    r2 = sxy ** 2 / (sxx * syy)
    print(f"  Chen = {my - sl * mx:.2f} + {sl:.2f} * Fischer   "
          f"r2 = {r2:.2f}  on {len(common)} carboxylates (nCF {min(common)}-{max(common)})")
    print()
    print(f"  The {sl:.2f} is the same compression seen from the other side")
    print("  (0.18 = 1 - 0.82, fitted on nine compounds rather than ten), so")
    print(f"  the two agree. What is new is the r2 of {r2:.2f}: Chen's curve is")
    print("  not a shrunken copy of Fischer's, because almost none of the")
    print("  variance Chen does have is Fischer's variance. That leftover")
    print("  variance is the shape argument.")
    print(f"{'nCF':>4s} {'Fischer':>8s} {'Chen':>8s} {'fitted':>8s} {'resid':>8s}")
    res = []
    for n in common:
        fit = my + sl * (f[n] - mx)
        res.append((n, c[n] - fit))
        print(f"{n:4d} {f[n]:8.2f} {c[n]:8.2f} {fit:8.2f} {c[n] - fit:8.2f}")
    lo = [r for n, r in res if n <= 8]
    hi = [r for n, r in res if n >= 11]
    print()
    print(f"  mean residual, nCF <= 8  {sum(lo) / len(lo):+.2f}")
    print(f"  mean residual, nCF >= 11 {sum(hi) / len(hi):+.2f}")
    print()
    print("The residuals are not scatter. They are positive across nCF 6-10")
    print("and negative at 11-12, so they carry a chain-length trend with a")
    print("sign change in it. Whichever study is right, the disagreement is")
    print("about the long-chain end specifically, and no rescaling of either")
    print("curve will reconcile the two.")


def section3(S):
    banner("3. Three shape statistics, and what each one does and does not say")
    use = {k: v for k, v in S.items() if len(v["s"]) >= 5}
    print("Restricting to series with at least 5 points, which drops Alesio")
    print("2024 and Starnes 2024.")
    print()
    print("(a) SMOOTHNESS. A single binding cavity filling and then being")
    print("    outgrown gives one reversal. More than one means the series is")
    print("    not tracing a cavity.")
    print()
    print(f"{'dataset':26s} {'pts':>4s} {'sign chg':>9s} {'peak nCF':>9s}")
    for k, v in sorted(use.items(), key=lambda kv: sign_changes(kv[1]["s"])):
        s = v["s"]
        mark = f"{peak(s):9d}" if interior(s) else f"{'>=' + str(peak(s)):>9s}"
        print(f"{k:26s} {len(s):4d} {sign_changes(s):9d} {mark}")
    print()
    print("    This one discriminates cleanly. Every series that places the")
    print("    peak at nCF 10 or never reaches it is smooth (0 or 1 reversal).")
    print("    The two that place it early wander 3 and 4 times.")
    print()
    print("(b) DYNAMIC RANGE over nCF 4-9, a window every series covers, so")
    print("    the comparison does not just reward reaching further.")
    print()
    print(f"{'dataset':26s} {'span 4-9':>9s} {'per CF2':>8s} {'peak nCF':>9s}")
    for k, v in use.items():
        sp, n = window_span(v["s"])
        if sp is None:
            continue
        s = v["s"]
        mark = f"{peak(s):9d}" if interior(s) else f"{'>=' + str(peak(s)):>9s}"
        print(f"{k:26s} {sp:9.2f} {sp / 5:8.2f} {mark}")
    print()
    print("    Physical anchors for what the per-CF2 increment should be:")
    print(f"      Chen's own membrane-water K_MW, nCF 3-7   "
          f"{(CHEN_KMW[7] - CHEN_KMW[3]) / 4:+.2f}")
    print(f"      Ryu 2024 tissue f_unbound, nCF 4-13        "
          f"{ryu_slope():+.2f} (binding up)")
    print("      octanol-water for PFCAs                     ~+0.6")
    print()
    print("(c) RISING LIMB, measured on each series' own limb (shortest")
    print("    compound to one CF2 below its own peak) rather than a fixed")
    print("    window, because the limbs are curved in different directions.")
    print()
    print(f"{'dataset':26s} {'limb/CF2':>9s} {'pts':>4s} {'peak nCF':>9s}")
    for k, v in use.items():
        sl, n = own_limb(v["s"])
        if sl is None:
            continue
        s = v["s"]
        mark = f"{peak(s):9d}" if interior(s) else f"{'>=' + str(peak(s)):>9s}"
        print(f"{k:26s} {sl:9.2f} {n:4d} {mark}")
    print(f"{'Bischel 2011 BSA':26s} "
          f"{BISCHEL11['span_log'] / BISCHEL11['n_cf2']:9.2f} {'n/a':>4s} "
          f"{BISCHEL11['peak_ncf']:9d}")
    print()
    pairs = []
    for k, v in use.items():
        sl, _ = own_limb(v["s"])
        if sl is not None:
            pairs.append((sl, peak(v["s"])))
    pairs.append((BISCHEL11["span_log"] / BISCHEL11["n_cf2"],
                  BISCHEL11["peak_ncf"]))
    print(f"    Spearman rho between limb slope and peak position across the")
    print(f"    {len(pairs)} series: {spearman(pairs):+.2f}.")
    print()
    print("    This test is weaker than (a) and it is worth saying why rather")
    print("    than leaning on the rho. The purified-protein series behave as")
    print("    the hypothesis predicts -- Fischer's two albumins and its")
    print("    globulin climb at 0.50 to 0.58 per CF2 and peak at nCF 10-11,")
    print("    while Chen (0.38), Starnes (0.10) and Bischel (0.12) climb")
    print("    shallowly and peak at 6-7. But Fischer's WHOLE SERUM breaks it:")
    print("    a shallow limb and a late peak. That is explicable (the")
    print("    denominator there is total serum protein, most of which is not")
    print("    albumin, so the short-chain end is propped up by everything")
    print("    else and the limb flattens) but an explicable exception is")
    print("    still an exception. Treat (c) as corroboration, not evidence.")
    print()
    print("Taken together: smoothness separates the two camps cleanly, range")
    print("and limb slope point the same way with one named exception, and no")
    print("statistic in this section places the optimum below nCF 9.")



def section4(S):
    banner("4. Fischer's turnover survives the test a method artefact fails")
    print("Fischer measured four sorbents in one experiment with one method.")
    print("If the fall after nCF 10 were the SPME fiber or the chemistry")
    print("giving out, it would give out at a fixed log D -- a ceiling is a")
    print("level. Instead it happens at a fixed chain length across sorbents")
    print("whose levels differ by more than a log unit:")
    print()
    print(f"{'sorbent':26s} {'peak nCF':>9s} {'log D at peak':>14s} "
          f"{'fall to nCF 13':>15s}")
    for k in ("Fischer 2024 HSA", "Fischer 2024 BSA", "Fischer 2024 globulin",
              "Fischer 2024 whole serum"):
        s = S[k]["s"]
        p = peak(s)
        fall = s[p] - s.get(13, s[max(s)])
        print(f"{k:26s} {p:9d} {s[p]:14.2f} {fall:15.2f}")
    print()
    lv = [S[k]['s'][peak(S[k]['s'])] for k in
          ("Fischer 2024 HSA", "Fischer 2024 globulin", "Fischer 2024 whole serum")]
    print(f"The peak LEVEL spans {max(lv) - min(lv):.2f} log units across the three")
    print("sorbents; the peak POSITION spans one carbon. A ceiling in the")
    print("measurement cannot produce that pattern.")
    print()
    print("The one-carbon spread is itself informative rather than noise.")
    ga = S["Fischer 2024 globulin"]["s"]
    ha = S["Fischer 2024 HSA"]["s"]
    print(f"  albumin  peaks at nCF {peak(ha)} ({ha[peak(ha)]:.2f}), "
          f"nCF 11 is {ha.get(11, float('nan')):.2f}")
    print(f"  globulin peaks at nCF {peak(ga)} ({ga[peak(ga)]:.2f}), "
          f"nCF 10 is {ga.get(10, float('nan')):.2f}")
    print("Two proteins, same tube, same fiber, optima one CF2 apart. That is")
    print("what two binding cavities of slightly different depth look like,")
    print("and it is not what a shared artefact looks like.")
    print()
    print("The two albumins are the replicate check on the fall itself:")
    h, b = S["Fischer 2024 HSA"]["s"], S["Fischer 2024 BSA"]["s"]
    for n in sorted(set(h) & set(b)):
        print(f"   nCF {n:2d}   HSA {h[n]:.2f}   BSA {b[n]:.2f}   "
              f"diff {h[n] - b[n]:+.2f}")
    d13 = abs(h[13] - b[13])
    print(f"They agree to {d13:.2f} log units at nCF 13, the point furthest down")
    print("the falling limb, while that point sits")
    print(f"{h[peak(h)] - h[13]:.2f} log units below the HSA peak. A "
          f"{10 ** (h[peak(h)] - h[13]):.0f}-fold fall reproduced")
    print("to within a few percent on a second protein is not scatter.")


def section5(S):
    banner("5. Concentration headroom: which studies could see nCF 12-15 at all")
    print("The 2026 Starnes review states the problem outright for its own")
    print("method family: 'the poor aqueous solubility of many PFAS")
    print("congeners, particularly long-chain PFAAs (C12 and above), prevents")
    print("accurate measurement of their binding in the assay.' C12 is nCF 11.")
    print("That is precisely the region where the datasets disagree.")
    print()
    print("What each study actually put in the tube:")
    print()
    for k in ("Fischer 2024 HSA", "Chen 2025 HSA", "Gao 2025 BSA",
              "Starnes 2025 HSA", "Alesio 2024 BSA"):
        if k in S:
            print(f"  {k:24s} {S[k]['conc']}")
    print()
    top_f = 42e-6 / 714.0            # 42 ug/L of PFTeDA, g/L over g/mol
    top_c = 2e-6                     # 2 uM
    print(f"  For PFTeDA (nCF 13, MW 714) the top nominal aqueous level is")
    print(f"    Fischer  {top_f * 1e9:6.1f} nM   ({42} ug/L)")
    print(f"    Chen     {top_c * 1e9:6.1f} nM   (2 uM = {2e-6 * 714 * 1000:.2f} mg/L)")
    print(f"    ratio    {top_c / top_f:6.0f}x more concentrated in Chen")
    print()
    print("Measured aqueous solubilities for nCF >= 11 PFCAs are scarce and")
    print("the modelled values scatter over orders of magnitude, so no sharp")
    print("solubility line can be drawn here, and this section does not try.")
    print("What can be said is the ordering. Fischer works at tens of")
    print("nanomolar, 34x below Chen's top level and two to four orders below")
    print("the DSF assays, which need micromolar to millimolar to move a")
    print("melting temperature at all. The study with the most headroom is")
    print("the one that reports a smooth, replicated turnover; the studies")
    print("with the least are the ones whose long-chain tails misbehave.")
    print("That is an ordering, not a proof, and it is only worth anything")
    print("because the misbehaviour is visible directly:")
    print()
    c = S["Chen 2025 HSA"]["s"]
    print()
    print("Chen's own long-chain tail misbehaves in a way a real binding")
    print("curve cannot:")
    print(f"   nCF 10 {c[10]:.2f}   11 {c[11]:.2f}   12 {c[12]:.2f}   "
          f"13 {c[13]:.2f}   15 {c.get(15, float('nan')):.2f}")
    print("Affinity weakens to nCF 12 and then TIGHTENS again at 13 and 15.")
    print("The printed uncertainties make that a real reversal, not rounding:")
    print("  K_D(nCF 12) = 31.4 +/- 5.1 uM against K_D(nCF 13) = 13.1 +/- 3.0.")
    print("A binding cavity does not do that. A sorption or solubility")
    print("artefact switching on partway along the series does.")
    print()
    print("And it is specific to the protein measurement. Chen measured a")
    print("membrane-water partition coefficient on the same compounds in the")
    print("same study, and it is perfectly behaved to nCF 15:")
    ks = sorted(CHEN_KMW)
    print("   " + "  ".join(f"{n}:{CHEN_KMW[n]:.2f}" for n in ks))
    d = [(CHEN_KMW[b] - CHEN_KMW[a]) / (b - a) for a, b in zip(ks, ks[1:])]
    print(f"   per-CF2 increments {min(d):+.2f} to {max(d):+.2f}, no sign change, "
          f"no turnover.")
    print("So the long-chain compounds were present and quantifiable in that")
    print("lab. Whatever goes wrong goes wrong in the dialysis, not in the")
    print("chemicals.")


def section6(S):
    banner("6. What the optimum is NOT: it is invisible in bulk tissue")
    print("This is where the question most often gets muddled. Albumin has a")
    print("cavity, and a cavity has a size, so a cavity can have an optimum.")
    print("A lipid bilayer does not, and neither does 'tissue'.")
    print()
    print("Chen's membrane-water coefficient (section 5) rises monotonically")
    print("to nCF 15. Ryu 2024's measured unbound fractions, across five")
    print("tissues in four species, fall monotonically to nCF 13 -- binding")
    print("keeps getting stronger with every CF2, with no hint of a maximum:")
    print()
    rows = load()
    acc = collections.defaultdict(list)
    for r in rows:
        if not r["study"].startswith("Ryu") or r["head_group"] != "carboxylate":
            continue
        fu, n = num(r["unbound_fraction"]), num(r["n_perfluorinated_carbons"])
        if fu and n and fu > 0:
            acc[int(n)].append(math.log10(fu))
    ks = sorted(acc)
    print("   " + "  ".join(f"{n}:{sum(acc[n]) / len(acc[n]):.2f}" for n in ks))
    print(f"   mean log10 f_unbound; slope {ryu_slope():+.2f} per CF2 over "
          f"nCF 4-13, no sign change")
    print()
    print("So three statements that sound contradictory are all true:")
    print("  - tissue binding has no chain-length optimum (Ryu, to nCF 13)")
    print("  - membrane partitioning has no chain-length optimum (Chen, to 15)")
    print("  - albumin binding does have one (Fischer, at nCF 10)")
    print()
    print("They are consistent because past the optimum the long-chain PFCAs")
    print("do not stop binding, they stop binding to ALBUMIN and go on")
    print("partitioning into lipid instead. Anything measured on whole tissue")
    print("or whole serum sums the two and sees the sum keep rising. That is")
    print("why an optimum found in a purified-protein assay is not refuted by")
    print("an in vivo series that shows none, and why the reverse inference --")
    print("'half-lives keep rising, so albumin affinity must too' -- does not")
    print("follow either.")


def section7(S):
    banner("7. Where the optimum is")
    print("Fitting a quadratic through each bracketed peak and its three")
    print("neighbours either side puts the vertex at:")
    print()
    print(f"{'dataset':26s} {'argmax nCF':>11s} {'vertex nCF':>11s}")
    for k, v in S.items():
        if not interior(v["s"]):
            continue
        vt = vertex(v["s"])
        print(f"{k:26s} {peak(v['s']):11d} "
              f"{(f'{vt:.1f}' if vt else 'n/a'):>11s}")
    print()
    h = S["Fischer 2024 HSA"]["s"]
    b = S["Fischer 2024 BSA"]["s"]
    print("Taking the two Fischer albumins as the measurement that survives")
    print("sections 3 to 5, and reading the top of the curve rather than a")
    print("single point:")
    print()
    print(f"{'nCF':>4s} {'HSA':>7s} {'BSA':>7s} {'mean':>7s} "
          f"{'below peak':>11s}")
    pk = max((h[n] + b[n]) / 2 for n in set(h) & set(b))
    for n in sorted(set(h) & set(b)):
        m = (h[n] + b[n]) / 2
        print(f"{n:4d} {h[n]:7.2f} {b[n]:7.2f} {m:7.2f} {m - pk:11.2f}")
    print()
    for thr in (0.15, 0.25, 0.50):
        flat = [n for n in sorted(set(h) & set(b)) if (h[n] + b[n]) / 2 >= pk - thr]
        print(f"Within {thr:.2f} log units of the top: nCF {flat}  "
              f"({', '.join(NAME[n] for n in flat)})")
    flat = [n for n in sorted(set(h) & set(b)) if (h[n] + b[n]) / 2 >= pk - 0.25]
    wide = [n for n in sorted(set(h) & set(b)) if (h[n] + b[n]) / 2 >= pk - 0.50]
    print()
    print("So the honest answer is a plateau, not a point.")
    print()
    print("  THE OPTIMUM FOR PERFLUOROCARBOXYLATES ON SERUM ALBUMIN SITS AT")
    print(f"  nCF {min(flat)}-{max(flat)} ({NAME[min(flat)]} to {NAME[max(flat)]}, "
          f"the C{min(flat) + 1} and C{max(flat) + 1} acids),")
    print(f"  centred on nCF 10 ({NAME[10]}) at log D ~ {pk:.1f} L_water per kg")
    print(f"  protein. Widening the tolerance to half a log unit brings in")
    print(f"  nCF {min(wide)} ({NAME[min(wide)]}) and nothing else.")
    print()
    print("Confidence, stated plainly:")
    print("  - that the optimum exists and is not an artefact: strong. The")
    print(f"    two albumins agree to {abs(h[10] - b[10]):.2f} log at the peak and "
          f"{abs(h[13] - b[13]):.2f} at nCF 13,")
    print("    the far end of a 20-fold fall, and the position is invariant")
    print("    across three sorbents whose levels differ by 1.3 log units.")
    print("  - that it sits at nCF 10 rather than 9 or 11: moderate. The top")
    print("    is flat and no study publishes per-compound error bars on this")
    print("    stretch; nCF 9-11 cannot be separated.")
    print("  - that it is NOT at nCF 6-8: strong, and this is the part that")
    print("    matters. Fischer's two albumins peak at 10; Gao's fluorescence")
    print("    quenching is still climbing at nCF 9; Qin's nonspecific SPME")
    print("    term is still climbing at nCF 10. Those are three different")
    print("    methods in three different labs, and a maximum at nCF 7 needs")
    print("    all three to be wrong in the same direction. Chen's own")
    print("    quadratic vertex, fitted through its peak and neighbours, lands")
    print(f"    at {vertex(S['Chen 2025 HSA']['s']):.1f} rather than at its argmax of 7.")
    print()
    print("Open, and worth saying: the sulfonate optimum is NOT located.")
    print("Fischer's PFSA series ends at nCF 10 (PFDS) and is still at or")
    print("near its top there, so all that can be said is that it lies at")
    print("nCF >= 9. No dataset in this file carries a C12 sulfonate.")


def section8(S):
    banner("8. What it changes")
    h = S["Fischer 2024 HSA"]["s"]
    pk = h[peak(h)]
    print("1. The regulated compounds are on the rising limb, not at the top.")
    print()
    print(f"{'compound':10s} {'nCF':>4s} {'log D':>7s} {'below optimum':>14s}")
    for n, name in ((7, "PFOA"), (8, "PFNA"), (9, "PFDA"), (10, "PFUnDA")):
        if n in h:
            print(f"{name:10s} {n:4d} {h[n]:7.2f} {h[n] - pk:14.2f}")
    print()
    print(f"PFOA sits {pk - h[7]:.2f} log units -- a factor of {10 ** (pk - h[7]):.0f} -- below the")
    print("albumin optimum. Every PFAS with a human health advisory is on the")
    print("rising limb, which is the regime where affinity is most sensitive")
    print("to chain length and therefore where a QSAR has the most to say.")
    print()
    sl, _ = limb_slope(h, 4, 9)
    print(f"2. A QSAR fitted on nCF <= 9 will extrapolate straight through the")
    print(f"   optimum. Fischer's HSA limb over nCF 4-9 runs at {sl:+.2f} per CF2;")
    print("   extending that line past the peak:")
    for n in (11, 12, 13):
        if n in h:
            pred = h[4] + sl * (n - 4)
            print(f"     nCF {n}: linear extrapolation {pred:.2f} vs measured "
                  f"{h[n]:.2f}, over by {10 ** (pred - h[n]):.0f}x")
    print("   Any chain-length descriptor trained on the compounds that have")
    print("   toxicity values will do this unless the training set reaches")
    print("   nCF 12, which almost none do.")
    print()
    print("3. It does not rescue the albumin hypothesis for the mass balance.")
    print("   The earlier finding stands: substituting Fischer's measured")
    print("   D_alb makes both published mass-balance models worse, and the")
    print("   inverse problem demands values above anything measured. Moving")
    print(f"   to the optimum would only help compounds at nCF 10-11, and the")
    print("   compounds that failed -- PFHxA at nCF 5, PFOA at 7 -- are on the")
    print("   rising limb where the optimum is irrelevant. PFHxA needed log D")
    print(f"   3.8-4.0 against {h.get(5, float('nan')):.2f} measured; the optimum is "
          f"five carbons away")
    print("   and cannot be borrowed.")


# Hand-placed label offsets for the one panel where points collide.
NUDGE = {"Starnes 2025 HSA": (-8, 16), "Bischel 2011 BSA": (14, -17),
         "Fischer 2024 HSA": (-22, 15), "Fischer 2024 BSA": (24, -18),
         "Fischer 2024 whole serum": (0, 22), "Fischer 2024 globulin": (0, 22),
         "Gao 2025 BSA": (16, 14), "Qin 2024 BSA nonspecific": (-14, -20),
         "Chen 2025 HSA": (0, -20)}


def figure(S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
    SURFACE = "#fcfcfb"
    INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
    RULE = "#d8d7d3"
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
        "font.size": 9, "axes.edgecolor": RULE, "axes.linewidth": 0.8,
        "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
        "ytick.color": INK2, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.fontsize": 8,
    })

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.5, 4.3))

    # -- panel 1: every series, each centred on its own mean so shapes compare
    order = ["Fischer 2024 HSA", "Fischer 2024 BSA", "Qin 2024 BSA nonspecific",
             "Gao 2025 BSA", "Chen 2025 HSA", "Starnes 2025 HSA",
             "Qin 2024 BSA specific"]
    style = {"Fischer 2024 HSA": (BLUE, "-", 2.0, "o"),
             "Fischer 2024 BSA": (BLUE, "--", 1.4, "s"),
             "Qin 2024 BSA nonspecific": (AQUA, "-", 1.4, "^"),
             "Qin 2024 BSA specific": (AQUA, "--", 1.4, "v"),
             "Gao 2025 BSA": ("#7a5cc6", "-", 1.4, "D"),
             "Chen 2025 HSA": (ORANGE, "-", 2.0, "o"),
             "Starnes 2025 HSA": (MUTED, "-", 1.2, "x")}
    for k in order:
        if k not in S:
            continue
        s = S[k]["s"]
        ks = sorted(s)
        m = sum(s.values()) / len(s)
        col, ls, lw, mk = style[k]
        ax1.plot(ks, [s[n] - m for n in ks], ls, color=col, lw=lw, marker=mk,
                 ms=4, label=k)
        p = peak(s)
        if interior(s):
            ax1.plot([p], [s[p] - m], marker="*", ms=13, color=col,
                     mec=SURFACE, mew=0.8, zorder=5)
    ax1.axvspan(9.5, 11.5, color=BLUE, alpha=0.08, lw=0, zorder=0)
    ax1.set_xlabel("perfluorinated carbons (nCF)")
    ax1.set_ylabel("log affinity, centred on each series' own mean")
    ax1.set_title("Every albumin PFCA series, shapes overlaid\n"
                  "stars mark a bracketed maximum", fontsize=9.5, color=INK)
    ax1.legend(loc="lower left", fontsize=7.4)

    # -- panel 2: each series' own rising limb vs where it says the peak is
    pts = []
    for k, v in S.items():
        if len(v["s"]) < 5:
            continue
        sl, _ = own_limb(v["s"])
        if sl is None:
            continue
        pts.append((sl, peak(v["s"]), interior(v["s"]), k))
    pts.append((BISCHEL11["span_log"] / BISCHEL11["n_cf2"],
                BISCHEL11["peak_ncf"], True, "Bischel 2011 BSA"))
    pts.sort()
    ax2.axvspan(0.40, 0.65, color=AQUA, alpha=0.10, lw=0, zorder=0)
    for i, (sl, pk, inter, k) in enumerate(pts):
        col = BLUE if inter else MUTED
        ax2.scatter([sl], [pk], s=58, color=col, zorder=3,
                    marker="o" if inter else "^")
        lab = {"Fischer 2024 HSA": "Fischer HSA",
               "Fischer 2024 BSA": "Fischer BSA"}.get(
            k, k.replace("Fischer 2024 ", "Fischer ")
                .replace("Qin 2024 BSA ", "Qin ")
                .replace(" BSA", "").replace(" HSA", ""))
        dx, dy = NUDGE.get(k, (0, (11, -18, 24, -31)[i % 4]))
        ax2.annotate(lab, (sl, pk), textcoords="offset points",
                     xytext=(dx, dy),
                     ha="center" if abs(dx) < 12 else ("right" if dx < 0 else "left"),
                     fontsize=7.2, color=INK2,
                     arrowprops=dict(arrowstyle="-", lw=0.5, color="#c9c8c4",
                                     shrinkA=0, shrinkB=5))
    ax2.set_ylim(5.0, 12.4)
    ax2.set_xlim(0.0, 0.82)
    ax2.annotate("physically expected\n+0.4 to +0.6 per CF2", (0.525, 5.15),
                 ha="center", va="bottom", fontsize=7.6, color=INK2)
    ax2.set_xlabel("rise per CF2 on each series' own rising limb")
    ax2.set_ylabel("apparent peak nCF  (triangles = still rising)")
    ax2.set_title("Shallow rising limb, early apparent peak\n"
                  "with whole serum the named exception",
                  fontsize=9.5, color=INK)

    # -- panel 3: Fischer's four sorbents
    ax3.axvspan(9.5, 11.5, color=BLUE, alpha=0.08, lw=0, zorder=0)
    for k, col in (("Fischer 2024 HSA", BLUE), ("Fischer 2024 BSA", "#6fa8e8"),
                   ("Fischer 2024 globulin", ORANGE),
                   ("Fischer 2024 whole serum", AQUA)):
        s = S[k]["s"]
        ks = sorted(s)
        ax3.plot(ks, [s[n] for n in ks], "-o", color=col, lw=1.7, ms=4,
                 label=k.replace("Fischer 2024 ", ""))
        p = peak(s)
        ax3.plot([p], [s[p]], marker="*", ms=14, color=col, mec=SURFACE,
                 mew=0.8, zorder=5)
    ax3.set_xlabel("perfluorinated carbons (nCF)")
    ax3.set_ylabel(r"log D (L$_{water}$ per kg protein)")
    ax3.set_title("One method, four sorbents, 1.3 log units apart\n"
                  "all turn over at nCF 10-11", fontsize=9.5, color=INK)
    ax3.legend(loc="upper left", fontsize=8)
    ax3.set_ylim(1.4, 5.8)

    fig.suptitle("Where the albumin chain-length optimum is, and why the "
                 "datasets disagree about it", fontsize=11.5, color=INK, y=1.02)
    fig.tight_layout()
    path = os.path.join(OUT, "fig17_albumin_chain_length_optimum.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"wrote {os.path.abspath(path)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", action="store_true")
    args = ap.parse_args()
    rows = load()
    S = build(rows)
    if args.figure:
        figure(S)
        return
    print("=" * 78)
    print("THE ALBUMIN CHAIN-LENGTH OPTIMUM FOR PERFLUOROCARBOXYLATES")
    print("=" * 78)
    section1(S)
    section2(S)
    section3(S)
    section4(S)
    section5(S)
    section6(S)
    section7(S)
    section8(S)
    print()


if __name__ == "__main__":
    main()
