#!/usr/bin/env python3
"""Per-congener L-FABP and HSA Kd for Gao 2019 and Zhang 2013.

Both papers' per-congener tables were recorded in this database as
unretrievable.  Fan et al. 2025 Table 1 reprints them, so they are
transcribed here FROM FAN, as a secondary source, and every row says so.

  Gao K, Zhuang T, Liu X, Fu J, Zhang J, Fu J, Wang L, Zhang A, Liang
  Y, Song M, Jiang G.  "Prenatal Exposure to Per- and Polyfluoroalkyl
  Substances (PFASs) and Association between the Placental Transfer
  Efficiencies and Dissociation Constant of Serum Proteins-PFAS
  Complexes."  Environ Sci Technol 2019;53(11):6529-6538.
  doi 10.1021/acs.est.9b00715, PMID 31099564.
  Equilibrium dialysis, 8 kDa MWCO, 4 h -- read from the primary, whose
  full text is in papers/.

  Zhang L, Ren X-M, Guo L-H.  "Structure-based investigation on the
  interaction of perfluorinated compounds with human liver fatty acid
  binding protein."  Environ Sci Technol 2013;47(19):11293-11301.
  doi 10.1021/es4026722, PMID 24006842.
  NBD-lauric acid fluorescence displacement (Fan's footnote c).

WHY A SECONDARY TRANSCRIPTION IS ACCEPTABLE HERE, exceptionally: both
columns can be checked against something independent.

  Gao PFBA.  Gao's own main text states "they were 482 +/- 45, 419 +/-
  37, and 879 +/- 51 uM, respectively" for Kd-SP, Kd-HP and Kd-LP.
  Fan's HSA column gives PFBA 419 +/- 37 and its hL-FABP column gives
  879 +/- 51.  Exact, on both proteins.

  Gao HSA, five more.  This database already held PFHxA 195, PFOA 115,
  PFNA 73, PFDoA 103 and PFOS 38 for Gao, transcribed independently
  from Jia et al. 2022.  Fan's HSA column reproduces all five exactly.
  Two independent secondary transcriptions agreeing is close to a
  primary read.

  Zhang.  Han et al. 2021 quoted six of Zhang's values (PFBA >1000,
  PFBS 1034, PFHxS 85.7, PFOA 50.4, PFOS 18.5, PFNA 16.2) and this
  database flagged them as "quoted from another study, not measured
  here" without being able to check them.  Fan's column c reproduces
  all six.  Those Han rows can now be read against Zhang's own.

WHAT IS STILL MISSING: neither paper's assay concentrations reach
Fan's table, and Gao's are in its SI rather than its main text.  The
concentration columns stay empty rather than carrying a guess.

Run:  python3 scripts/add_gao2019_zhang2013_via_fan.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

VIA = ("Fan et al. 2025 Table 1 (Environ Int 201:109591, "
       "doi 10.1016/j.envint.2025.109591), which reprints it")

GAO = dict(study="Gao et al.", year="2019", pmid="31099564",
           doi="10.1021/acs.est.9b00715")
ZHANG = dict(study="Zhang et al.", year="2013", pmid="24006842",
             doi="10.1021/es4026722")

GAO_METHOD = ("equilibrium dialysis, rapid equilibrium dialysis device, "
              "8 kDa molecular weight cutoff (Thermo Fisher), 4 h to "
              "equilibrium, protein precipitated with 90/10 "
              "acetonitrile/water + 0.1% formic acid, then LC-MS/MS "
              "-- method read from the primary, values from Fan")
ZHANG_METHOD = ("NBD-lauric acid fluorescence displacement against "
                "hL-FABP expressed in E. coli (per Fan's footnote c)")

# (compound, nCF, head, Zhang hL-FABP, Gao hL-FABP, Gao HSA)
# "" is an empty cell in Fan's table; "ND" is printed as ND there.
T1 = [
    ("PFBA",   3,  "carboxylate", "ND",            "879 +/- 51",  "419 +/- 37"),
    ("PFPeA",  4,  "carboxylate", "ND",            "685 +/- 23",  "286 +/- 17"),
    ("PFHxA",  5,  "carboxylate", "ND",            "443 +/- 25",  "195 +/- 14"),
    ("PFHpA",  6,  "carboxylate", "336.0 +/- 40.3", "358 +/- 16", "172 +/- 23"),
    ("PFOA",   7,  "carboxylate", "50.4 +/- 4.0",  "197 +/- 13",  "115 +/- 16"),
    ("PFNA",   8,  "carboxylate", "16.2 +/- 0.6",  "90 +/- 8",    "73 +/- 21"),
    ("PFDA",   9,  "carboxylate", "12.9 +/- 0.6",  "77 +/- 14",   "65 +/- 8"),
    ("PFUnDA", 10, "carboxylate", "10.6 +/- 1.2",  "54 +/- 5",    "79 +/- 13"),
    ("PFDoDA", 11, "carboxylate", "12.3 +/- 0.7",  "95 +/- 11",   "103 +/- 12"),
    ("PFTeDA", 13, "carboxylate", "60.5 +/- 6.7",  "317 +/- 20",  "265 +/- 21"),
    ("PFBS",   4,  "sulfonate",   "1034 +/- 110",  "436 +/- 33",  "169 +/- 15"),
    ("PFHxS",  6,  "sulfonate",   "85.7 +/- 11.8", "210 +/- 14",  "79 +/- 12"),
    ("PFOS",   8,  "sulfonate",   "18.5 +/- 0.6",  "81 +/- 7",    "38 +/- 5"),
]

# Rows this database already holds for Gao's HSA, transcribed from Jia
# et al. 2022.  They are replaced, and the agreement is asserted below.
JIA_CHECK = {"PFHxA": 195.0, "PFOA": 115.0, "PFNA": 73.0, "PFOS": 38.0}

TK = ("yes - Gao measured HSA and L-FABP side by side and related the "
      "ratio to placental transfer efficiency, which is the same "
      "liver-versus-blood competition Fan et al. 2025 later measured "
      "in rats")


def blank(fields):
    return dict.fromkeys(fields, "")


def num(v):
    return v.split("+/-")[0].strip() if v and v != "ND" else ""


def rows(fields):
    out = []
    for chem, ncf, head, zh, gl, gh in T1:
        for src, prot, meth, val, tbl in (
            (GAO, "human serum albumin (HSA)", GAO_METHOD, gh,
             "Kd-HP (HSA-PFAS complex)"),
            (GAO, "liver fatty acid-binding protein (hL-FABP)",
             GAO_METHOD, gl, "Kd-LP (L-FABP-PFAS complex)"),
            (ZHANG, "liver fatty acid binding protein (L-FABP), "
                    "recombinant, expressed in E. coli",
             ZHANG_METHOD, zh, "hL-FABP Kd"),
        ):
            if not val:
                continue
            r = blank(fields)
            r.update(src)
            r.update(
                chemical=chem, n_perfluorinated_carbons=str(ncf),
                head_group=head, protein=prot, protein_species="human",
                method=meth, kd=num(val),
                kd_units=("uM" if num(val) else
                          "ND, no binding detected (as printed)"),
                is_computational="no", tk_relevance_noted=TK,
                protein_conc="", protein_conc_units="",
                ligand_conc_range="", ligand_conc_units="",
                source_table=(
                    f"{tbl} = {val} uM | SECONDARY: this paper's "
                    f"per-congener table was not retrieved; the value is "
                    f"transcribed from {VIA}"))
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, existing = reader.fieldnames, list(reader)

    # Check the Jia-sourced Gao HSA values against Fan's column BEFORE
    # replacing them. A disagreement would mean one transcription is
    # wrong and must stop the run.
    fan_hsa = {c: float(gh.split("+/-")[0]) for c, _, _, _, _, gh in T1}
    for chem, old in JIA_CHECK.items():
        assert abs(fan_hsa[chem] - old) < 1e-9, (
            f"{chem}: Fan says {fan_hsa[chem]}, the existing Jia-sourced "
            f"row says {old} -- transcriptions disagree, not proceeding")
    print(f"cross-check passed: Fan's Gao-HSA column matches all "
          f"{len(JIA_CHECK)} values this database already held from Jia "
          f"2022 ({', '.join(sorted(JIA_CHECK))})")

    before = len(existing)
    drop = {("Gao et al.", "2019"), ("Zhang et al.", "2013")}
    existing = [r for r in existing
                if (r["study"], r["year"]) not in drop]
    print(f"removed {before - len(existing)} existing Gao 2019 / "
          f"Zhang 2013 rows")

    new = rows(fields)
    print(f"adding  {len(new)} rows -> {len(existing) + len(new)} total")
    for key, label in ((("Gao et al.", "human serum albumin"), "Gao HSA"),
                       (("Gao et al.", "liver fatty"), "Gao hL-FABP"),
                       (("Zhang et al.", "liver fatty"), "Zhang hL-FABP")):
        got = [r for r in new
               if r["study"] == key[0] and key[1] in r["protein"]]
        with_kd = sum(1 for r in got if r["kd"])
        print(f"   {label:<14} {len(got):>3} rows, {with_kd} with a Kd")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(existing + new)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    main()
