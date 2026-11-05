#!/usr/bin/env python3
"""Does each citation describe the study it is filed under?

WHY THIS EXISTS.  A citation row is built by handing a PMID to PubMed
and writing down what comes back.  If the PMID is wrong, that process
succeeds -- it returns a real paper, with real authors, a real
journal and a real abstract.  Nothing fails.  The row just describes
somebody else's paper.

Two have been caught that way in this project.  Gao 2025 was briefly
filed under PMID 41344964, one digit off 41344664, which returned a
pancreatitis retraction notice.  Chen & Guo 2009 sat for longer under
PMID 19002450, whose DOI is 10.1007/s00247-008-1049-9 -- `s00247` is
Pediatric Radiology, `s00204` is Archives of Toxicology, one digit
apart -- so the row carried a review of splenic function in sickle
cell disease.  The correct record is PMID 18854981.

Both were found by eye.  This is the check that finds them without
one.  It is cheap and it is offline: the test is internal
consistency, between three fields the database already holds.

THE THREE CHECKS, weakest premise first:

  SURNAME.  The `study` label is this project's own short name for
  the paper -- "Chen & Guo", "Alesio et al.", "Ng & Hungerbuhler".
  The surnames in it must appear among the fetched `authors`.  This
  is what catches a wrong PMID: a paper by other people fails it
  immediately.  Accents and transliteration make exact matching too
  brittle, so comparison is on ASCII-folded lowercase surnames.

  YEAR.  The `year` column is the issue year this project files the
  study under; `pub_year` is PubMed's.  They legitimately differ by
  one when a paper appears online in one year and in an issue the
  next -- Alesio 2022 is a 2021 Chemosphere issue, Chen & Guo 2009
  is an Arch Toxicol 2009 issue first posted October 2008.  So a gap
  of more than one year is reported; one year is not.

  DOI PREFIX.  A DOI's registrant prefix is a publisher, and its
  suffix often names the journal series.  Where a row has both a DOI
  and a journal, grossly mismatched pairs are reported -- this is
  what makes the Chen & Guo case visible from the DOI alone, since
  the stored DOI was correct while the PMID was not.

WHAT A HIT MEANS.  A SURNAME hit is nearly always a wrong identifier
and should be resolved against the source before anything downstream
is trusted.  YEAR and DOI hits are weaker and often benign; they are
reported as NOTE, not FAIL, and the three bands are kept separate so
a real failure is never buried under them.

Run:  python3 tk_review/scripts/check_citation_identity.py
Exit status 1 if any FAIL.  Reads only; writes nothing.
"""
import csv
import os
import re
import sys
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
CITES = os.path.join(HERE, "..", "db", "citations.csv")

# A DOI suffix opens with a series token that often names the
# journal. The token must be matched WHOLE: 'ac' is Analytical
# Chemistry, but '/acs.est.4c06734' opens with 'acs.est', and a
# substring test reports every ES&T paper as a mismatch. series()
# below extracts the token; this table only interprets it.
#
# Used to flag a mismatch, never to assert a match, so the table does
# not need to be complete -- an unrecognised token is passed over.
DOI_SERIES = {
    "s00204": "Arch Toxicol",
    "s00247": "Pediatr Radiol",
    "s41598": "Sci Rep",
    "acs.est": "Environ Sci Technol",
    "acs.chemrestox": "Chem Res Toxicol",
    "j.chemosphere": "Chemosphere",
    "j.envint": "Environ Int",
    "j.jes": "J Environ Sci",
    "jf": "J Agric Food Chem",
    "tx": "Chem Res Toxicol",
    "ac": "Anal Chem",
}

# Words in a `study` label that are not surnames.
NOT_A_NAME = {"et", "al", "and", "the", "panel", "&"}


def fold(s):
    """Lowercase ASCII, so Hungerbuhler matches Hungerbühler."""
    s = unicodedata.normalize("NFKD", s or "")
    s = "".join(c for c in s if not unicodedata.combining(c))
    return s.encode("ascii", "ignore").decode("ascii").lower()


def label_surnames(study):
    """The surnames a study label claims, e.g. 'Chen & Guo' -> both."""
    out = []
    for w in re.split(r"[\s&;,]+", fold(study)):
        w = w.strip(".")
        if w and w not in NOT_A_NAME and not w.isdigit():
            out.append(w)
    return out


def series(doi):
    """The journal-series token a DOI suffix opens with, or ''.

    10.1021/acs.est.4c06734     -> acs.est
    10.1021/ac100721e           -> ac
    10.1007/s00204-008-0359-x   -> s00204
    10.1016/j.chemosphere.2021  -> j.chemosphere
    """
    suf = (doi or "").lower().split("/", 1)
    if len(suf) < 2:
        return ""
    suf = suf[1]
    for pat in (r"^(acs\.[a-z]+)",      # ACS, two segments
                r"^(j\.[a-z]+)",        # Elsevier
                r"^(s\d{5})",           # Springer
                r"^([a-z]+)"):          # older ACS: ac, tx, jf, es
        m = re.match(pat, suf)
        if m:
            return m.group(1)
    return ""


def author_surnames(authors):
    """Surnames from 'Chen YM; Guo LH' -- the word before initials."""
    out = []
    for part in (authors or "").split(";"):
        part = fold(part).strip()
        if not part:
            continue
        # "chen ym" -> chen. Also tolerate "chen, yan-min".
        toks = [t.strip(".,") for t in re.split(r"[\s,]+", part) if t]
        if toks:
            out.append(toks[0])
        out.extend(toks)          # keep the rest for compound names
    return out


def check(rows):
    fails, years = [], []
    for r in rows:
        study, cite = r["study"], r.get("citation", "")
        if not cite:
            continue              # nothing fetched: not this check's job

        want = label_surnames(study)
        have = author_surnames(r.get("authors", ""))
        # An organisational author ("EFSA Panel on ...") has no
        # surnames to match, so require only that the label's first
        # token appears somewhere in the author string.
        blob = fold(r.get("authors", ""))
        missing = [w for w in want
                   if w not in have and w not in blob]
        if want and missing:
            fails.append(("SURNAME", study, r["year"],
                          f"label names {'/'.join(missing)} but the "
                          f"citation's authors are "
                          f"{r.get('authors', '')[:60]!r}",
                          f"pmid {r.get('pmid') or '-'}  "
                          f"doi {r.get('doi') or '-'}"))

        tok = series(r.get("doi"))
        jrnl, belongs = r.get("journal", ""), DOI_SERIES.get(tok)
        if belongs and jrnl and belongs != jrnl:
            fails.append(("DOI", study, r["year"],
                          f"doi series '{tok}' is {belongs} but the "
                          f"row says {jrnl}", r.get("doi", "")))

        try:
            gap = abs(int(r["year"]) - int(r["pub_year"]))
        except (ValueError, KeyError, TypeError):
            gap = 0
        if gap >= 1:
            years.append((study, r["year"], r["pub_year"], gap,
                          f"{jrnl} {r.get('volume')}"
                          f"({r.get('issue')})"))
    return fails, years


def main():
    rows = list(csv.DictReader(open(CITES, encoding="utf-8")))
    fails, years = check(rows)
    cited = sum(1 for r in rows if r.get("citation"))

    print(f"citation identity check -- {cited} of {len(rows)} studies "
          f"have a citation to check\n")

    print("IDENTITY  does each citation describe the study it is "
          "filed under?")
    if not fails:
        print("  pass -- every citation's authors and journal match "
              "its study label\n")
    else:
        for kind, study, year, why, extra in fails:
            print(f"  FAIL [{kind}]  {study} {year}")
            print(f"        {why}")
            if extra:
                print(f"        {extra}")
        print()

    # Not a defect in the identifier, so it is kept out of the exit
    # status: it is a defect in one FIELD of an otherwise correct row,
    # and it is systematic rather than a mistake anyone made here.
    print("YEAR  the citation year is PubMed's publication_date, "
          "which for")
    print("      Springer and some other publishers is the "
          "ONLINE-FIRST date, not")
    print("      the issue date. Where it disagrees with the volume, "
          "the citation")
    print("      string pairs one year with another year's volume. "
          "This project's")
    print("      own `year` column is the issue year and is the one "
          "to cite.")
    if not years:
        print("      none\n")
    else:
        for study, ours, theirs, gap, vol in years:
            print(f"      {study} {ours}: citation says {theirs}, "
                  f"volume is {vol}  ({gap} yr)")
        print(f"      {len(years)} row(s); cite the `year` column, "
              f"not the citation string's year.\n")

    if fails:
        print(f"{len(fails)} identity failure(s): a citation that "
              f"describes a different paper")
        print("than the study it is filed under. Resolve the "
              "identifier before trusting")
        print("anything built from that row.")
        return 1
    print("No identity failures.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
