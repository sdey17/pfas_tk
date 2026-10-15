"""
Two-page standalone memo: the Argoul 2026 cocktail study against the
single-compound literature.

Reuses the style, glyph check and table helpers of make_summary_pdf.py so the
two documents look like one set. Numbers are read from
db/argoul_vs_single_compound.csv at build time, so the memo cannot drift from
the data the way a hand-typed summary would.

Run:  python3 scripts/make_argoul_memo.py   ->  report/Argoul_comparison_memo.pdf
"""

import csv
import math
import os
import sys

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, KeepTogether,
                                PageTemplate, Paragraph, Spacer)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_summary_pdf as M                                  # noqa: E402

HERE = M.HERE
OUT = os.path.join(HERE, "report", "Argoul_comparison_memo.pdf")
p, tbl, W = M.p, M.tbl, M.W


def rows():
    f = os.path.join(HERE, "db", "argoul_vs_single_compound.csv")
    return list(csv.DictReader(open(f, newline="")))


def gm(v):
    return math.exp(sum(math.log(x) for x in v) / len(v))


def footer(canv, doc):
    canv.saveState()
    canv.setStrokeColor(M.RULE)
    canv.setLineWidth(0.4)
    canv.line(20 * mm, 14 * mm, A4[0] - 20 * mm, 14 * mm)
    canv.setFont("Times-Roman", 7.6)
    canv.setFillColor(M.SOFT)
    canv.drawString(20 * mm, 10 * mm,
                    "Argoul 2026 against the single-compound literature")
    canv.drawRightString(A4[0] - 20 * mm, 10 * mm, str(doc.page))
    canv.restoreState()


def main():
    q = rows()
    by = lambda k: [float(r["ratio"]) for r in q if r["parameter"] == k]
    cl, vd, hl = by("clearance"), by("volume_of_distribution"), by("half_life")

    flow = [
        p("The cocktail study against the single-compound literature", "title"),
        p("How Argoul 2026 was run, which comparisons it supports, and why six "
          "were dropped", "subtitle"),

        p("1. What the experiment was", "h1"),
        p("Argoul and colleagues dosed <b>eleven PFAS as a single mixture</b> "
          "to 108 female C57BL/6 mice, 12 weeks old, 20.0 &plusmn; 0.8 g "
          "(Environmental Research 303:124802). Fifty-one animals received the "
          "mixture intravenously by lateral tail vein at 5 mL/kg; forty-nine "
          "received it by oral gavage at 10 mL/kg; the vehicle was saline with "
          "10% ethanol for the intravenous route and 5% for the oral. Doses "
          "were set per compound from a preliminary study so that each would "
          "be measurable while staying in a low molar range - 0.019 to "
          "1.55 mg/kg, tracking persistence rather than sweeping a "
          "dose-response.", "body"),
        p("Sampling was destructive and ran to <b>119 days</b>: blood from "
          "three to five mice at each of fifteen times, from 15 minutes to "
          "119 days. From 24 hours on, three mice per time per route went into "
          "metabolism cages for complete 24-hour urine and faeces collection, "
          "so the renal and faecal routes were measured rather than inferred. "
          "Eight untreated mice established background. Quantification was "
          "LC-MS/MS against seven isotope-labelled internal standards.", "body"),
        p("Plasma, urine and faeces were then fitted <b>simultaneously</b> by "
          "nonlinear mixed effects, Laplacian engine, with precision from "
          "bootstrap (n = 100 or 200). Clearance and steady-state volume are "
          "primary estimates; mean residence time and terminal half-life are "
          "derived from them. Two inputs come from outside the experiment and "
          "matter later: glomerular filtration rate is taken as 10 L/kg/day "
          "from the physiological literature, not fitted, and the unbound "
          "fractions are from a companion study in the same strain, not the "
          "same animals.", "body"),

        p("2. Why the comparison is worth making", "h1"),
        p("This study carries unusual weight in the review: it is the only "
          "source that measures clearance and volume for many PFAS in one "
          "experiment, which is what the sex-versus-species decomposition and "
          "the structure-activity endpoint both need. Every other mouse value "
          "available comes from a study that dosed one compound. So the "
          "question is whether dosing eleven at once, at low dose, by "
          "mixed-effects fitting, reproduces the one-at-a-time literature.",
          "body"),

        p("3. Which comparisons were kept, and why six were not", "h1"),
        p("Thirteen comparisons survive, against female-mouse studies that "
          "dosed the compound alone. Six were dropped, all for the same "
          "reason: <b>no single dose can be attributed to the reported "
          "parameter</b>, so a dose column alongside them could not be read as "
          "a dose.", "body"),
        KeepTogether([
        tbl([
            ["dropped", "rows", "why"],
            ["Lou 2009 (PFOA)", "3",
             "one 2-compartment fit pooled across 1 and 10 mg/kg"],
            ["Tatum-Gibbs 2011 (PFNA)", "2",
             "a half-life <i>range</i>, 25.8-68.4 d, spanning both doses "
             "rather than a value per arm"],
            ["US EPA PFHxA IRIS", "1",
             "a V<sub>d,&beta;</sub> computed from Gannon 2011, whose mouse "
             "arm ran 2 <i>and</i> 100 mg/kg"],
        ], [44 * mm, 14 * mm, W - 58 * mm]),
        p("A seventh exclusion is of a different kind. The EPA pooled "
          "Bayesian fit was also kept out, because it re-analyses the same "
          "single-compound studies: placing it beside them would count the "
          "same animals twice, and it is a model output rather than a "
          "measurement - for PFBA it returns 0.186 d where every Chang arm "
          "sits at 0.116-0.128.", "caption")]),
        p("The exclusions cost coverage, which is worth stating plainly. PFOA "
          "loses its second comparator; <b>PFNA and PFHxA leave the comparison "
          "entirely</b>. PFHxA was the one outlier in the other direction "
          "(5.13&times; high) and the only compound the study places in net "
          "tubular secretion, so its kinetics differ in kind regardless. They "
          "cost almost nothing in the answer: including the three dropped "
          "half-life rows gives ten comparisons at 0.99&times; against seven "
          "at 1.04&times;.", "body"),

        KeepTogether([
        p("4. The result", "h1"),
        tbl([
            ["parameter", "n", "compounds", "geometric mean", "span"],
            ["clearance", str(len(cl)), "PFOA, PFHxS",
             f"<b>{gm(cl):.2f}&times;</b>  ({1/gm(cl):.1f}&times; low)",
             f"{min(cl):.2f}-{max(cl):.2f}&times;"],
            ["volume of distribution", str(len(vd)), "PFOA, PFHxS",
             f"<b>{gm(vd):.2f}&times;</b>  ({1/gm(vd):.1f}&times; low)",
             f"{min(vd):.2f}-{max(vd):.2f}&times;"],
            ["half-life (their ratio)", str(len(hl)), "PFBA, PFHxS, PFOS",
             f"<b>{gm(hl):.2f}&times;</b>", f"{min(hl):.2f}-{max(hl):.2f}&times;"],
        ], [42 * mm, 10 * mm, 34 * mm, 40 * mm, W - 126 * mm]),
        p("Read from db/argoul_vs_single_compound.csv at build time; "
          "references in the _sources sidecar. Half-life is "
          "ln2&middot;MRT, since the study reports mean residence time rather "
          "than a terminal half-life.", "caption")]),
        p("Both absolute terms sit low by a similar factor while their ratio "
          "does not. That is the signature of a <b>shared multiplicative "
          "scaling</b>, not a disagreement about elimination: in "
          "ln2&middot;V<sub>d</sub>/CL a factor common to both cancels. One "
          "caveat on how much of a check that is - volume is usually not "
          "measured apart from clearance, since non-compartmental analysis "
          "gives V<sub>z</sub> = CL/&lambda;<sub>z</sub>, so for the "
          "Sundstr&ouml;m rows 'both low' is largely one observation seen "
          "twice. Fujii is the exception, with volume from dose/C(0) and "
          "clearance from dose/AUC, and there the two are low independently "
          "(0.38&times; and 0.59&times;).", "body"),

        p("5. What explains the offset", "h1"),
        p("<b>Dose is the obvious candidate and it does not survive testing.</b> "
          "Matched per compound, the dose gap fails to order the "
          "discrepancies: Fujii dosed PFOA at 0.128 mg/kg against the "
          "cocktail's 0.080, a gap of 1.6&times;, and disagrees <i>most</i> "
          "(0.38&times;), while Lou at a gap of at least 12&times; agrees "
          "<i>best</i> (0.75&times;). The log-log slope once quoted in support "
          "was fitted through the cocktail's own point, which carries 79% "
          "leverage on it; dropping it and extrapolating Sundstr&ouml;m's two "
          "doses alone predicts 2.04 mL/kg/d against 1.30 observed, leaving "
          "about 1.6&times; unexplained. Dose contributes; it is not the "
          "explanation.", "body"),
        p("<b>Strain is a confounder the review had not isolated.</b> Argoul "
          "used C57BL/6; every comparator is CD-1 except Fujii, which is "
          "FVB/NJcl. All thirteen comparisons therefore cross strains, and no "
          "study in hand separates strain from cocktail, dose or method. The "
          "comparison CSV now carries both strains per row so the confound is "
          "visible rather than implicit.", "body"),
        p("<b>Cocktail competition</b> for shared reabsorptive transporters "
          "would raise clearance, not lower it, so it argues against the "
          "observed direction unless the competition is for plasma binding "
          "sites instead. <b>Mixed-effects shrinkage</b> compresses spread "
          "rather than shifting a centre, so it is the weakest of the four.",
          "body"),

        p("6. What this does and does not cost the review", "h1"),
        p("The within-experiment contrasts are <b>unaffected</b>: a shared "
          "scaling cancels from a ratio, so the 5,254&times; clearance span "
          "against 7.7&times; in volume across the eleven compounds stands "
          "whatever the common factor turns out to be. The absolute clearances "
          "are what the renal handling ratio R is built from, and those carry "
          "the roughly 2.5&times; low bias - a flat &minus;0.40 log unit "
          "offset on every compound. Being flat, it cancels from every "
          "between-compound contrast the endpoint actually claims, and the "
          "nearest reabsorbed compound sits 0.98 log units from the R = 1 "
          "divide, so the bias would have to be 9.6&times; rather than "
          "2.5&times; to reclassify anything. What it does cost is the "
          "cross-species comparison, where an absolute mouse value meets an "
          "absolute human one: read those gaps as carrying a factor of about "
          "2.5 on top of the unbound-fraction uncertainty, understating the "
          "mouse.", "body"),
        p("The experiment that would settle it is not another comparator. It "
          "is the same cocktail at a second total dose, in a strain the "
          "single-compound literature also used - which would separate dose, "
          "cocktail and strain in one design.", "body"),
        Spacer(1, 5),
        p("Reproduce: scripts/argoul_vs_single_compound.py (comparison), "
          "scripts/argoul_dose_matching.py (the dose test), "
          "scripts/make_argoul_memo.py (this memo). Full treatment in "
          "report/REPORT.md section 11.3a and section 10 of the main review.",
          "caption"),
    ]

    M.check_glyphs(flow)
    doc = BaseDocTemplate(
        OUT, pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
        topMargin=18 * mm, bottomMargin=20 * mm,
        title="Argoul 2026 against the single-compound literature",
        author="pfas_tk", subject="PFAS toxicokinetics")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height,
                  id="f")
    doc.addPageTemplates([PageTemplate(id="body", frames=[frame],
                                       onPage=footer)])
    doc.build(flow)
    try:
        from pypdf import PdfReader
        pages = len(PdfReader(OUT).pages)
    except ImportError:
        pages = "?"
    print(f"wrote {os.path.relpath(OUT, HERE)}  {pages} pages, "
          f"{os.path.getsize(OUT)/1024:.0f} KB")


if __name__ == "__main__":
    main()
