#!/usr/bin/env python3
"""Muensterman et al. 2026 SI Table S13: biomimetic partition coefficients.

Environ Sci Technol 60(1):1207-1217, doi 10.1021/acs.est.5c11716.
Two rows per compound, from the same table:

  log K(HSA)  human serum albumin-water, the protein limb
  log K(IAM)  immobilised artificial membrane-water, the lipid limb

Both are CHROMATOGRAPHIC SURROGATES, not equilibrium measurements. Retention
on a protein- or phospholipid-coated column is calibrated against a set of
pharmaceuticals of known partition coefficient, and the PFAS value is read
off that calibration. That makes this the largest protein-water data set in
the file by a wide margin, and the least like the others in kind. The SI
says so itself: "Calibration is typically reliable within its defined range
... and caution must be exercised when extrapolating beyond its dynamic
range", and several PFAS here sit well outside it -- log K(IAM) reaches 9.1
for PFODA against a calibrant span topping out near 5.

log k(AGP), the third column, is a retention FACTOR on an alpha-1-acid
glycoprotein column, not a partition coefficient. It is carried in the
source note rather than as a value, so it cannot be mistaken for one.

Only the target-PFAS block is parsed. The 37 suspect PFAS lower in the same
table have names that wrap across two lines, which splits their numeric
columns; parsing those without risking a misalignment is not possible from
the text layer, and a misaligned row is worse than a missing one.

CAUTION, recorded on every row: a Correction to this paper was published
(Environ Sci Technol 60(23):17051, doi 10.1021/acs.est.6c07668). Its
content is not in hand, so these values are the originals as published and
have NOT been checked against it.

Run:  python3 scripts/add_muensterman2026_biomimetic.py TABLE_S13.txt [--dry-run]
"""
import argparse
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "db", "protein_binding.csv")

NUM = r"(?:-?\d+(?:\.\d+)?|\*\*|N/A)"
ROW = re.compile(r"^\s*(?:\S[\S ]*?)??([A-Za-z0-9][\w:\-+()' ]*?)\s+"
                 + r"\s+".join([f"({NUM})"] * 7) + r"\s*$")

# nCF and head group for the compounds we already carry, so the new rows join
# the existing chain-length series instead of sitting apart from them.
KNOWN = {
    "PFBA": (3, "carboxylate"), "PFPeA": (4, "carboxylate"),
    "PFHxA": (5, "carboxylate"), "PFHpA": (6, "carboxylate"),
    "PFOA": (7, "carboxylate"), "PFNA": (8, "carboxylate"),
    "PFDA": (9, "carboxylate"), "PFUdA": (10, "carboxylate"),
    "PFDoA": (11, "carboxylate"), "PFTrDA": (12, "carboxylate"),
    "PFTeDA": (13, "carboxylate"), "PFHxDA": (15, "carboxylate"),
    "PFODA": (17, "carboxylate"),
    "PFPrS": (3, "sulfonate"), "PFBS": (4, "sulfonate"),
    "PFPeS": (5, "sulfonate"), "PFHxS": (6, "sulfonate"),
    "PFHpS": (7, "sulfonate"), "PFOS": (8, "sulfonate"),
    "PFNS": (9, "sulfonate"), "PFDS": (10, "sulfonate"),
    "PFDoS": (12, "sulfonate"), "PFTrDS": (13, "sulfonate"),
    "HFPO-DA": ("", "ether carboxylate (polyfluoroether)"),
    "ADONA": ("", "ether carboxylate (polyfluoroether)"),
    "FOSA": (8, "sulfonamide"), "FHxSA": (6, "sulfonamide"),
    "6:2 FTS": (6, "fluorotelomer sulfonate"),
    "8:2 FTS": (8, "fluorotelomer sulfonate"),
    "4:2 FTS": (4, "fluorotelomer sulfonate"),
    "10:2 FTS": (10, "fluorotelomer sulfonate"),
}

# Rename to the acronyms the rest of the file uses, so a compound does not
# appear twice under two spellings.
ALIAS = {"PFUdA": "PFUnDA", "PFDoA": "PFDoDA", "HFPO-DA": "HFPO-DA (GenX)",
         "6:2 FTS": "6:2 FtS", "8:2 FTS": "8:2 FtS", "10:2 FTS": "10:2 FtS"}

METHOD = ("biomimetic chromatography, rapid-gradient; retention on a "
          "protein- or phospholipid-coated column calibrated against "
          "pharmaceuticals of known partition coefficient (SURROGATE, not "
          "an equilibrium measurement)")
CAUTION = ("| CAUTION: a Correction to this paper was published (Environ Sci "
           "Technol 60(23):17051, doi 10.1021/acs.est.6c07668); its content "
           "is not in hand and these are the ORIGINAL published values")

SORBENT = {
    "HSA": ("human serum albumin (HSA), immobilised stationary phase", "human",
            "log K_HSA (albumin-water, biomimetic); basis not stated in the "
            "paper, so not converted"),
    "IAM": ("immobilised artificial membrane (IAM, phospholipid stationary "
            "phase)", "not applicable (synthetic phospholipid)",
            "log K_IAM (phospholipid membrane-water, biomimetic); basis not "
            "stated in the paper, so not converted"),
}


def val(x):
    return "" if x in ("**", "N/A") else x


# The left-hand class column bleeds into the name on the line where a class
# label starts, so it is stripped explicitly rather than by guesswork.
CLASS_LABELS = ["Me/Et FASAA", "Me/Et FASA", "n:2 UFTCA", "n:2 FTCA",
                "n:3 FTCA", "n:2 FTOH", "n:2 FTS", "n:2 PAP", "Cl-PFESA",
                "Cl-PFECA", "Cl-PFSA", "Cy-PFSA", "Zwitterionic PFAS",
                "Fatty acids", "FASAA", "FASA", "FASE", "PFECA", "PFESA",
                "PFCA", "PFSA", "Me/Et", "Et"]
# Three names wrap across two lines, which leaves only the parenthetical on
# the line carrying the numbers.
WRAPPED = {"FHET)": "6:2 FTOH", "FOET)": "8:2 FTOH", "FDET)": "10:2 FTOH"}


def clean_name(name):
    name = name.strip()
    for lab in CLASS_LABELS:
        if name == lab:
            return ""
        if name.startswith(lab + " "):
            name = name[len(lab):].strip()
            break
    return WRAPPED.get(name, name)


def parse(path):
    out = []
    for line in open(path, errors="replace"):
        m = ROW.match(line.rstrip())
        if not m:
            continue
        name = clean_name(m.group(1))
        if not name or name.lower().startswith(("class", "acronym", "table",
                                                "pfas")):
            continue
        g = [val(m.group(i)) for i in range(2, 9)]
        if not (g[2] or g[4]):          # no partition coefficient at all
            continue
        out.append((name, dict(tr_iam=g[0], chi=g[1], k_iam=g[2],
                               tr_hsa=g[3], k_hsa=g[4],
                               tr_agp=g[5], k_agp=g[6])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table", help="text of SI Table S13 (target-PFAS block)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    parsed = parse(args.table)
    print(f"parsed {len(parsed)} compound rows")

    with open(CSV) as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, list(reader)

    new = []
    for name, d in parsed:
        chem = ALIAS.get(name, name)
        ncf, head = KNOWN.get(name, ("", ""))
        for tag, key in (("HSA", "k_hsa"), ("IAM", "k_iam")):
            if not d[key]:
                continue
            protein, species, units = SORBENT[tag]
            r = dict.fromkeys(fields, "")
            r.update(
                study="Muensterman et al.", year="2026", pmid="41473938",
                doi="10.1021/acs.est.5c11716", chemical=chem,
                n_perfluorinated_carbons=str(ncf), head_group=head,
                protein=protein, protein_species=species,
                method=METHOD, ka=d[key], ka_units=units,
                n_binding_sites="not applicable (partitioning, not site binding)",
                protein_conc="immobilised stationary phase",
                ligand_conc_range="10 uL injection of a standard",
                temperature_c="25",
                buffer="pH 7.4 mobile phase; isopropanol gradient 0-30%",
                is_computational="no",
                tk_relevance_noted=(
                    "the largest single set of protein-water and "
                    "membrane-water coefficients for PFAS, but a "
                    "chromatographic surrogate rather than an equilibrium "
                    "measurement"),
                source_table=(f"SI Table S13 | tR(IAM) {d['tr_iam']}, CHI_IAM "
                              f"{d['chi']}, tR(HSA) {d['tr_hsa']}, log k(AGP) "
                              f"{d['k_agp'] or 'n/a'} (an AGP retention "
                              f"FACTOR, not a partition coefficient) "
                              + CAUTION))
            new.append(r)

    print(f"{len(new)} rows to add")
    if args.dry_run:
        for r in new[:6]:
            print(f"   {r['chemical']:10s} {r['protein'][:34]:34s} {r['ka']}")
        return
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows + new)
    print(f"wrote {os.path.abspath(CSV)}  ({len(rows) + len(new)} rows)")


if __name__ == "__main__":
    main()
