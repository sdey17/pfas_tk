#!/usr/bin/env python3
"""Fischer 2024 as a fifth human unbound-fraction series (Stage 8).

    python3 -I scripts/stage8_ingest_fischer2024_fu.py

Appends to db/plasma_unbound_fraction.csv and writes
db/fischer2024_fu_derived.csv. Idempotent.

    Fischer FC, Ludtke S, Thackray C, Pickard HM, Haque F, Dassuncao C,
    Endo S, Schaider L, Sunderland EM (2024). "Binding of Per- and
    Polyfluoroalkyl Substances (PFAS) to Serum Proteins: Implications for
    Toxicokinetics in Humans." Environ Sci Technol 58:1055-1063.
    PMID 38166384, PMC11149785. doi:10.1021/acs.est.3c07415

THE DATA WERE ALREADY IN THIS PROJECT, IN THE WRONG STRAND. The albumin
strand ingested Fischer's Table S2 in full -- 29 compounds x 4 matrices (HSA,
BSA, gamma-globulin, whole human serum) as log D values -- into
db/protein_binding.csv. The transporter strand never used it, which is why
sections 7 and 10 described three human fu series when the project held a
fourth.

TWO TIERS, BECAUSE THE DERIVATION DOES NOT REPRODUCE THE PAPER EXACTLY.

  STATED   four f_unbound values Fischer gives in its own Results text.
           Verified, ingested as measurements.

  DERIVED  the whole series, from the paper's own equation 3 rearranged,
           f_unbound = VF_w / D_serum/w with VF_w = 0.928.

The derivation reproduces the four stated values only to within 1.1 to 2.2
fold, and the ratio is not constant, so it is not a simple volume-fraction
mismatch -- Fischer's surrogate serum carried 40-90 g/L total protein against
the NHANES average used for VF. The derived series is therefore written to a
SEPARATE file and not into the measurement record, with the discrepancy
carried on every row.

                 derived   paper states   ratio
    PFOA         0.00131      0.00061      2.15
    PFHxS        0.00047      0.00041      1.13
    PFOS         0.00075      0.00042      1.80
    PFHpS        0.00039      0.00035      1.11

AND THE SERIES OVERTURNS AN ARGUMENT THIS REVIEW MADE. Section 7 used
"binding rises with chain length, so fu must fall monotonically" to prefer
Ryu's series over Argoul's at the long chains. FISCHER SHOWS THAT PREMISE IS
FALSE ABOVE ABOUT 8 PERFLUORINATED CARBONS: its own log D_HSA/w rises to 5.18
at PFUnDA and then FALLS to 4.54 and 3.87 for PFTriDA and PFTeDA, and the
paper gives the mechanism -- "for eta_pfc > 8 PFAS, reduced molecular
flexibility limits their accessibility to the HSA binding pockets". Smeltz
2023 reports the same independently: "binding was lower for the PFCAs with a
carbon chain length of >= 11".

So the monotonicity test was not valid at the chain lengths where it was used.
Argoul's PFUdA value is still an outlier -- 40-fold above Ryu's and 130-fold
above Fischer's derived value -- but it is an outlier on the evidence of two
other laboratories, not on a premise that turns out to be wrong.
"""
import csv, pathlib

HERE = pathlib.Path(__file__).resolve().parent.parent
VF_W = 0.928

# f_unbound stated in Fischer 2024 Results, "Distribution of PFAS in human serum"
STATED = {"PFOA": 0.00061, "PFHxS": 0.00041, "PFOS": 0.00042, "PFHpS": 0.00035}
SRC = ("Fischer 2024, Environ Sci Technol 58:1055-1063, "
       "doi:10.1021/acs.est.3c07415")
METHOD = ("C18 solid-phase microextraction fibre depletion + HPLC-MS/MS, "
          "surrogate human serum, 37 C, PBS pH 7.4, 48 h; PFAS:protein molar "
          "ratio <=0.004 (albumin) and <=0.009 (globulin)")
FIELDS = ["chemical", "mw_g_mol", "species", "experiment", "fu_percent", "fu",
          "range_low_percent", "range_high_percent", "n_replicates", "method",
          "source"]
# Fischer names compounds by perfluorinated carbons; map to this project's names
ALIAS = {"PFUdA": "PFUnDA", "PFDoA": "PFDoDA", "6:2 FtS": "6:2 FTS",
         "FOSA": "PFOSA"}


def fischer_logd():
    out = {}
    f = HERE / "db" / "protein_binding.csv"
    for r in csv.DictReader(f.open()):
        if r.get("study") != "Fischer et al." or r.get("year") != "2024":
            continue
        if not r.get("ka"):
            continue
        chem = ALIAS.get(r["chemical"], r["chemical"])
        out.setdefault(chem, {})[r["protein"].split(" (")[0]] = (
            float(r["ka"]), r.get("n_perfluorinated_carbons", ""),
            r.get("head_group", ""))
    return out


def main():
    d = fischer_logd()

    # ---- tier 1: the stated values, into the measurement record ------
    p = HERE / "db" / "plasma_unbound_fraction.csv"
    existing = list(csv.DictReader(p.open()))
    have = {(r["chemical"], r["species"], r["experiment"]) for r in existing}
    new = []
    for chem, fu in sorted(STATED.items()):
        key = (chem, "human", "Fischer 2024")
        if key in have:
            continue
        new.append({"chemical": chem, "mw_g_mol": "", "species": "human",
                    "experiment": "Fischer 2024", "fu_percent": fu * 100,
                    "fu": fu, "range_low_percent": "",
                    "range_high_percent": "", "n_replicates": "4",
                    "method": METHOD,
                    "source": SRC + ", Results text (f_unbound in surrogate "
                                    "human serum)"})
    if new:
        with p.open("a", newline="") as fh:
            csv.DictWriter(fh, fieldnames=FIELDS).writerows(new)
        print(f"appended {len(new)} stated Fischer rows to "
              f"{p.relative_to(HERE)}")
    else:
        print("stated Fischer rows already present (idempotent)")

    # ---- tier 2: the derived series, to its own file -----------------
    p2 = HERE / "db" / "fischer2024_fu_derived.csv"
    with p2.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chemical", "n_perfluorinated_carbons", "head_group",
                    "logD_HSA_w", "logD_globulin_w", "logD_serum_w",
                    "fu_derived", "fu_stated_in_paper",
                    "derived_over_stated", "derivation", "caveat", "source"])
        rows = []
        for chem in sorted(d, key=lambda c: (list(d[c].values())[0][2],
                                             int(list(d[c].values())[0][1] or 0))):
            m = d[chem]
            s = m.get("whole human serum")
            if not s:
                continue
            logd, npfc, head = s
            fu = VF_W / 10 ** logd
            st = STATED.get(chem)
            w.writerow([chem, npfc, head,
                        m.get("human serum albumin", ("",))[0],
                        m.get("gamma-globulin", ("",))[0], logd,
                        round(fu, 6), st if st else "",
                        round(fu / st, 2) if st else "",
                        "fu = VF_w / D_serum/w, VF_w = 0.928 "
                        "(Fischer eq. 2 and 3 rearranged)",
                        "DERIVED, not measured: reproduces the paper's four "
                        "stated f_unbound values only to within 1.1-2.2x, "
                        "and not by a constant factor. Do not use as a "
                        "measurement.", SRC + ", Table S2"])
            rows.append((chem, npfc, head, fu))

    report(rows)
    print(f"\nwrote {p2.relative_to(HERE)}  ({len(rows)} rows)")


def report(rows):
    bar = "=" * 78
    print()
    print(bar)
    print("THE SHORT-CHAIN PICTURE WITH FOUR HUMAN SOURCES")
    print(bar)
    fu_other = {
        "PFPeA": {"Ryu 2024": 0.0888, "Abraham 2024": 0.044},
        "PFHxA": {"Ryu 2024": 0.0380, "Abraham 2024": 0.0068,
                  "Argoul 2026": 0.0200},
        "PFHpA": {"Ryu 2024": 0.00445, "Abraham 2024": 0.0004},
    }
    fis = {c: fu for c, _n, _h, fu in rows}
    CLR = {"PFHxA": 47.0}          # mL/kg/day, Argoul 2026 and Abraham agree
    GFR = 2570.0
    for chem in ("PFPeA", "PFHxA", "PFHpA"):
        print(f"\n  {chem}")
        srcs = dict(fu_other[chem])
        srcs["Fischer 2024 (derived)"] = fis.get(chem)
        for n, v in sorted(srcs.items(), key=lambda kv: kv[1] or 0):
            if v is None:
                continue
            line = f"    {n:26}fu {v:.5f}"
            if chem in CLR:
                T = CLR[chem] / (v * GFR)
                line += (f"   T = {T:5.2f}   "
                         f"{'SECRETION' if T > 1 else 'reabsorption'}")
            print(line)
    print()
    print("  FOR PFHxA THE FOUR SOURCES SPLIT TWO AND TWO, AND THE SPLIT IS")
    print("  NOT RANDOM: the two lowest unbound fractions -- Abraham's, from")
    print("  the same subject as the clearance, and Fischer's -- BOTH give")
    print("  net secretion, and they agree with each other within 1.4-fold.")
    print("  Ryu sits 4-fold above Fischer and 15-fold above Fischer for")
    print("  PFPeA, which is the largest short-chain disagreement in the set.")
    print()
    print("  SO THE BALANCE OF EVIDENCE HAS MOVED. Section 7 reported human")
    print("  PFHxA at T = 0.48 on Ryu's binding, with Abraham's 3.48 as the")
    print("  alternative. With four sources the two that agree point to")
    print("  secretion. THIS DOES NOT SETTLE IT -- Fischer's value is derived")
    print("  rather than measured and carries a 1.1-2.2x reproduction error,")
    print("  which spans the boundary on its own -- but it is no longer a")
    print("  two-way tie, and the report should not present Ryu's reading as")
    print("  the default.")

    print()
    print(bar)
    print("AND THE LONG-CHAIN MONOTONICITY TEST WAS INVALID")
    print(bar)
    print("  Section 7 preferred Ryu over Argoul at the long chains because")
    print("  Ryu's carboxylate series falls monotonically and Argoul's")
    print("  inverts. Fischer's own albumin series does this:")
    print()
    carb = [(int(n), c, None) for c, n, h, _f in rows if h == "carboxylate"]
    import csv as _c
    hs = {}
    for r in _c.DictReader((HERE / "db" / "protein_binding.csv").open()):
        if (r.get("study") == "Fischer et al." and r.get("year") == "2024"
                and r["protein"].startswith("human serum albumin")
                and r["ka"]):
            hs[r["chemical"]] = float(r["ka"])
    ser = sorted(((int(n), c) for c, n, h, _f in rows if h == "carboxylate"))
    print(f"    {'eta_pfc':>8}{'compound':>10}{'log D_HSA/w':>14}")
    for n, c in ser:
        if c in hs:
            mark = ""
            print(f"    {n:>8}{c:>10}{hs[c]:>14.2f}{mark}")
    print()
    print("  IT RISES TO PFUnDA AND THEN FALLS. Fischer gives the mechanism:")
    print("  'for eta_pfc > 8 PFAS, reduced molecular flexibility limits")
    print("  their accessibility to the HSA binding pockets'. Smeltz 2023")
    print("  reports the same independently -- binding lower for PFCAs of")
    print("  chain length >= 11.")
    print()
    print("  SO 'fu MUST FALL MONOTONICALLY' IS FALSE ABOVE ABOUT 8")
    print("  PERFLUORINATED CARBONS, WHICH IS EXACTLY WHERE SECTION 7 APPLIED")
    print("  IT. The conclusion about Argoul's PFUdA survives on other")
    print("  grounds -- it is 40x Ryu's value and 130x Fischer's derived one,")
    print("  so two independent laboratories place it as an outlier -- but")
    print("  the argument given was wrong and is withdrawn.")


if __name__ == "__main__":
    main()
