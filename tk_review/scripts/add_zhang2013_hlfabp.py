#!/usr/bin/env python3
"""Zhang et al. 2013, read at source: 17 PFCs and 2 fatty acids on hL-FABP.

Zhang L, Ren X-M, Guo L-H.  "Structure-based investigation on the
interaction of perfluorinated compounds with human liver fatty acid
binding protein."  Environ Sci Technol 2013;47(19):11293-11301.
doi 10.1021/es4026722, PMID 24006842.  Table 1, read from the PDF.

This replaces the 13 rows transcribed from Fan et al. 2025 Table 1 two
commits ago.  Fan's column reproduces Zhang's Table 1 EXACTLY on all
ten shared compounds plus the three NDs -- the secondary transcription
was sound -- but these are now primary, and the primary carries six
things Fan's column did not:

  * IC50 alongside every Kd;
  * two more compounds, PFHxDA (C15) and PFOcDA (C17), which matter
    because they are what turns the chain-length trend over;
  * the two fatty acids used to validate the assay, oleic (Kd 4.93 uM)
    and palmitic (0.67 uM), which set the scale PFAS are weak against;
  * both fluorotelomer alcohols, not detected -- the paper's cleanest
    structural result, since 8:2 FTOH has the same chain as PFOS;
  * the probe's own constants, NBD-lauric acid Kd 280 +/- 20 nM at
    n = 1.12 +/- 0.04, i.e. ONE site;
  * the assay concentrations, which were blank before: 400 nM hL-FABP
    with 50 nM NBD-lauric acid, a probe:protein ratio of 1:8.

THE CHAIN-LENGTH TURN.  Affinity rises from C4 to C11 and falls after:
PFUnA is tightest at 10.6 uM, then PFDoA 12.3, PFTeDA 60.5, PFHxDA
115.4.  The paper attributes the fall above C11 to the binding cavity
running out of room.  Allendorf saw the same turn on albumin at a
similar place, and Fan saw a "U" on rat L-FABP, so three independent
datasets agree the trend is not monotonic.

Run:  python3 scripts/add_zhang2013_hlfabp.py [--dry-run]
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

STUDY, YEAR = "Zhang et al.", "2013"
DOI, PMID = "10.1021/es4026722", "24006842"
PROTEIN = ("liver fatty acid binding protein (hL-FABP), recombinant, "
           "expressed in E. coli, 98% pure by SDS-PAGE")
METHOD = ("NBD-lauric acid fluorescence displacement; the probe's own "
          "Kd (280 +/- 20 nM, n = 1.12 +/- 0.04, one site) measured "
          "first by titration, then PFC added over a concentration "
          "range and Kd derived from the IC50 of the displacement "
          "curve; n = 3")
PROT_CONC = ("400 nM hL-FABP with 50 nM NBD-lauric acid, a 1:8 "
             "probe:protein molar ratio", "nM")
LIG_CONC = ("a range of concentrations per compound; the displacement "
            "curves are SI Figure S3, not transcribed", "")

# Table 1: (chemical, nCF or None, head group, IC50 uM, Kd uM, note)
# "" for both means the table prints ND.
T1 = [
    ("oleic acid", None, "fatty acid (endogenous ligand)",
     "0.88 +/- 0.03", "4.93 +/- 0.17"),
    ("palmitic acid", None, "fatty acid (endogenous ligand)",
     "0.12 +/- 0.02", "0.67 +/- 0.11"),
    ("PFBA", 3, "carboxylate", "", ""),
    ("PFPeA", 4, "carboxylate", "", ""),
    ("PFHxA", 5, "carboxylate", "", ""),
    ("PFHpA", 6, "carboxylate", "60.0 +/- 7.2", "336.0 +/- 40.3"),
    ("PFOA", 7, "carboxylate", "9.0 +/- 0.7", "50.4 +/- 4.0"),
    ("PFNA", 8, "carboxylate", "2.9 +/- 0.1", "16.2 +/- 0.6"),
    ("PFDA", 9, "carboxylate", "2.3 +/- 0.1", "12.9 +/- 0.6"),
    ("PFUnDA", 10, "carboxylate", "1.9 +/- 0.2", "10.6 +/- 1.2"),
    ("PFDoDA", 11, "carboxylate", "2.2 +/- 0.1", "12.3 +/- 0.7"),
    ("PFTeDA", 13, "carboxylate", "10.8 +/- 1.2", "60.5 +/- 6.7"),
    ("PFHxDA", 15, "carboxylate", "20.6 +/- 1.4", "115.4 +/- 7.8"),
    ("PFOcDA", 17, "carboxylate", "11.1 +/- 0.9", "62.2 +/- 5.2"),
    ("PFBS", 4, "sulfonate", "185 +/- 20", "1034 +/- 110"),
    ("PFHxS", 6, "sulfonate", "15.3 +/- 2.1", "85.7 +/- 11.8"),
    ("PFOS", 8, "sulfonate", "3.3 +/- 0.1", "18.5 +/- 0.6"),
    ("6:2 FTOH", None, "fluorotelomer alcohol", "", ""),
    ("8:2 FTOH", None, "fluorotelomer alcohol", "", ""),
]

# Fan et al. 2025 Table 1 column c, the secondary transcription this
# replaces. Asserted against the primary before anything is written:
# if they ever disagree, one of the two readings is wrong and the run
# must stop rather than quietly overwrite.
FAN_C = {"PFHpA": 336.0, "PFOA": 50.4, "PFNA": 16.2, "PFDA": 12.9,
         "PFUnDA": 10.6, "PFDoDA": 12.3, "PFTeDA": 60.5,
         "PFBS": 1034.0, "PFHxS": 85.7, "PFOS": 18.5}

TK = ("yes - hL-FABP is the liver counterpart to albumin in blood. "
      "Affinity rises C4 to C11 then falls, which the paper attributes "
      "to the cavity running out of room; the same turn appears in "
      "Allendorf on albumin and in Fan on rat L-FABP")


def blank(fields):
    return dict.fromkeys(fields, "")


def rows(fields):
    out = []
    for chem, ncf, head, ic50, kd in T1:
        r = blank(fields)
        nd = not kd
        r.update(
            study=STUDY, year=YEAR, pmid=PMID, doi=DOI, chemical=chem,
            n_perfluorinated_carbons="" if ncf is None else str(ncf),
            head_group=head, protein=PROTEIN, protein_species="human",
            method=METHOD, temperature_c="",
            kd="" if nd else kd.split("+/-")[0].strip(),
            kd_units=("ND, not detected (as printed in Table 1)"
                      if nd else "uM"),
            n_binding_sites=("1 (the probe binds one site, n = 1.12 "
                             "+/- 0.04; the paper notes the cavity "
                             "holds two fatty acids)"),
            binding_site_identified=("primary fatty-acid site; Arg122 "
                                     "and Asn111 implicated"),
            protein_conc=PROT_CONC[0], protein_conc_units=PROT_CONC[1],
            ligand_conc_range=LIG_CONC[0], ligand_conc_units=LIG_CONC[1],
            buffer="", is_computational="no", tk_relevance_noted=TK,
            source_table=(
                f"Table 1 | " + ("ND, not detected" if nd else
                                 f"Kd {kd} uM | IC50 {ic50} uM")
                + " | read from the PDF, primary"))
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    primary = {c: float(kd.split("+/-")[0]) for c, _, _, _, kd in T1 if kd}
    for chem, secondary in FAN_C.items():
        got = primary.get(chem)
        assert got is not None and abs(got - secondary) < 1e-9, (
            f"{chem}: the primary says {got}, Fan's Table 1 said "
            f"{secondary} -- the two readings disagree, not proceeding")
    print(f"cross-check passed: the primary reproduces all "
          f"{len(FAN_C)} values transcribed from Fan's Table 1")

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, existing = reader.fieldnames, list(reader)

    before = len(existing)
    existing = [r for r in existing
                if not (r["study"] == STUDY and r["year"] == YEAR
                        and r["doi"] == DOI)]
    print(f"removed {before - len(existing)} secondary {STUDY} {YEAR} rows")

    new = rows(fields)
    kd = sum(1 for r in new if r["kd"])
    print(f"adding  {len(new)} primary rows ({kd} with a Kd, "
          f"{len(new) - kd} ND) -> {len(existing) + len(new)} total\n")
    for r in new:
        print(f"   {r['chemical']:<14}"
              f"{(r['kd'] + ' uM') if r['kd'] else 'ND':>12}")

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
