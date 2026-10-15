"""
Build the printable review: report/PFAS_TK_review.pdf

A mini literature review, written to stand on its own and be printed: the
question, the evidence base, nine substantive sections, what the work corrects
in the published record, what is still open, and a numbered reference list.

Two rules this file enforces mechanically rather than by care.

  * Every number in a table is read from db/ at build time, or transcribed from
    report/REPORT.md with the section noted. The document cannot drift from
    the data it describes.
  * In-text citations are keys into REFS below, numbered by order of first
    appearance at render time, so a reference cannot be mis-numbered and an
    unused or undefined one fails the build.

reportlab's built-in fonts are WinAnsi only; anything outside that set renders
as a solid black box, so there are no Unicode sub/superscripts anywhere (the
<sub>/<super> markup is used instead) and check_glyphs() fails the build
rather than shipping boxes.

Needs: reportlab.  Run:  python3 scripts/make_summary_pdf.py
"""

import csv
import os
import re
import sys

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "report", "PFAS_TK_review.pdf")
FIG = os.path.join(HERE, "figures", "fig00_master.png")
FIG2 = os.path.join(HERE, "figures", "fig02_sex_and_species_decomposition.png")
FIG3 = os.path.join(HERE, "figures",
                    "fig03_structure_activity_and_sensitivity.png")
FIG4 = os.path.join(HERE, "figures",
                    "fig04_argoul_vs_single_compound.png")

INK = colors.HexColor("#111111")
SOFT = colors.HexColor("#55534f")
RULE = colors.HexColor("#d8d7d3")
ACCENT = colors.HexColor("#2a78d6")

SAFE_EXTRA = set("\u00b5\u00d7\u00b7\u00b0\u2013\u2014\u2018\u2019"
                 "\u201c\u201d\u2026\u00bd\u00e9\u00f6\u00e5\u00f8"
                 "\u00c5")

REFS = {
"argoul2026": "Argoul CML, Toutain P-L, Picard-Hagen N, Mselli-Lakhal L, Dauwe Y, Roques BB, Lacroix MZ, Gayrard V (2026). Nonlinear mixed-effects modeling of the intravenous and oral kinetics of eleven perfluoroalkyl substances in female mice. <i>Environmental Research</i> 303:124802. doi:10.1016/j.envres.2026.124802",
"han2012": "Han X, Nabb DL, Russell MH, Kennedy GL, Rickard RW (2012). Renal elimination of perfluorocarboxylates (PFCAs). <i>Chemical Research in Toxicology</i> 25(1):35-46. doi:10.1021/tx200363w",
"kudo2001": "Kudo N, Suzuki E, Katakura M, Ohmori K, Noshiro R, Kawashima Y (2001). Comparison of the elimination between perfluorinated fatty acids with different carbon chain length in rats. <i>Chemico-Biological Interactions</i> 134(2):203-216. doi:10.1016/s0009-2797(01)00155-7",
"kudo2002": "Kudo N, Katakura M, Sato Y, Kawashima Y (2002). Sex hormone-regulated renal transport of perfluorooctanoic acid. <i>Chemico-Biological Interactions</i> 139(3):301-316. doi:10.1016/s0009-2797(02)00006-6",
"lou2009": "Lou I, Wambaugh JF, Lau C, Hanson RG, Lindstrom AB, Strynar MJ, Zehr RD, Setzer RW, Barton HA (2009). Modeling single and repeated dose pharmacokinetics of PFOA in mice. <i>Toxicological Sciences</i> 107(2):331-341. doi:10.1093/toxsci/kfn234",
"tatum2011": "Tatum-Gibbs K, Wambaugh JF, Das KP, Zehr RD, Strynar MJ, Lindstrom AB, Delinsky A, Lau C (2011). Comparative pharmacokinetics of perfluorononanoic acid in rat and mouse. <i>Toxicology</i> 281(1-3):48-55. doi:10.1016/j.tox.2011.01.003",
"thompson2010": "Thompson J, Lorber M, Toms L-ML, Kato K, Calafat AM, Mueller JF (2010). Use of simple pharmacokinetic modeling to characterize exposure of Australians to perfluorooctanoic acid and perfluorooctane sulfonic acid. <i>Environment International</i> 36(4):390-397. doi:10.1016/j.envint.2010.02.008, with its corrigendum, <i>Environment International</i> 36(6):652. doi:10.1016/j.envint.2010.05.008",
"yang2010": "Yang C-H, Glover KP, Han X (2010). Characterization of cellular uptake of perfluorooctanoate via organic anion-transporting polypeptide 1A2, organic anion transporter 4, and urate transporter 1 for their potential roles in mediating human renal reabsorption of perfluorocarboxylates. <i>Toxicological Sciences</i> 117(2):294-302. doi:10.1093/toxsci/kfq219",
"sundstrom2012": "Sundstr&ouml;m M, Chang S-C, Noker PE, Gorman GS, Hart JA, Ehresman DJ, Bergman &Aring;, Butenhoff JL (2012). Comparative pharmacokinetics of perfluorohexanesulfonate (PFHxS) in rats, mice, and monkeys. <i>Reproductive Toxicology</i> 33(4):441-451. doi:10.1016/j.reprotox.2011.07.004",
"cheng2005": "Cheng X, Maher J, Chen C, Klaassen CD (2005). Tissue distribution and ontogeny of mouse organic anion transporting polypeptides (Oatps). <i>Drug Metabolism and Disposition</i> 33(7):1062-1073. doi:10.1124/dmd.105.003640",
"cheng2006": "Cheng X, Maher J, Lu H, Klaassen CD (2006). Endocrine regulation of gender-divergent mouse organic anion-transporting polypeptide (Oatp) expression. <i>Molecular Pharmacology</i> 70(4):1291-1297. doi:10.1124/mol.106.025122 (abstract only; full text paywalled)",
"huang2019": "Huang MC, Dzierlenga AL, Robinson VG, Waidyanatha S, DeVito MJ, Eifrid MA, Granville CA, Gibbs ST, Blystone CR (2019). Toxicokinetics of perfluorobutane sulfonate, perfluorohexane-1-sulphonic acid, and perfluorooctane sulfonic acid in male and female Hsd:Sprague Dawley SD rats after intravenous and gavage administration. <i>Toxicology Reports</i> 6:645-655. doi:10.1016/j.toxrep.2019.06.016, with its corrigendum, <i>Toxicology Reports</i> 8:365. doi:10.1016/j.toxrep.2021.02.001",
"han2003": "Han X, Snow TA, Kemper RA, Jepson GW (2003). Binding of perfluorooctanoic acid to rat and human plasma proteins. <i>Chemical Research in Toxicology</i> 16(6):775-781. doi:10.1021/tx034005w",
"yang2009": "Yang C-H, Glover KP, Han X (2009). Organic anion transporting polypeptide (Oatp) 1a1-mediated perfluorooctanoate transport and evidence for a renal reabsorption mechanism of Oatp1a1 in renal elimination of perfluorocarboxylates in rats. <i>Toxicology Letters</i> 190(2):163-171. doi:10.1016/j.toxlet.2009.07.011",
"louisse2024": "Louisse J, Pedroni L, van den Heuvel JJMW, Rijkers D, Leenders L, Noorlander A, Punt A, Russel FGM, Koenderink JB, et al. (2024). In vitro and in silico characterization of the transport of selected perfluoroalkyl carboxylic acids and perfluoroalkyl sulfonic acids by human organic anion transporter 1 (OAT1), OAT2 and OAT3. <i>Toxicology</i> 509:153961. doi:10.1016/j.tox.2024.153961",
"louisse2023": "Louisse J, et al. (2023). Perfluoroalkyl substances (PFASs) are substrates of the renal human organic anion transporter 4 (OAT4). <i>Archives of Toxicology</i>. doi:10.1007/s00204-022-03428-6",
"ohmori2003": "Ohmori K, Kudo N, Katayama K, Kawashima Y (2003). Comparison of the toxicokinetics between perfluorocarboxylic acids with different carbon chain length. <i>Toxicology</i> 184(2-3):135-140. doi:10.1016/s0300-483x(02)00573-5",
"cheng2009": "Cheng X, Klaassen CD (2009). Tissue distribution, ontogeny, and hormonal regulation of xenobiotic transporters in mouse kidneys. <i>Drug Metabolism and Disposition</i> 37(11):2178-2185. doi:10.1124/dmd.109.027177",
"buist2004": "Buist SCN, Klaassen CD (2004). Rat and mouse differences in gender-predominant expression of organic anion transporter (Oat1-3; Slc22a6-8) mRNA levels. <i>Drug Metabolism and Disposition</i> 32(6):620-625. doi:10.1124/dmd.32.6.620",
"shi2016": "Shi Y, Vestergren R, Xu L, Zhou Z, Li C, Liang Y, Cai Y (2016). Human exposure and elimination kinetics of chlorinated polyfluoroalkyl ether sulfonic acids (Cl-PFESAs). <i>Environmental Science &amp; Technology</i> 50(5):2396-2404. doi:10.1021/acs.est.5b05849",
"maso2021": "Maso L, Trande M, Liberi S, Moro G, Daems E, Linciano S, et al. (2021). Unveiling the binding mode of perfluorooctanoic acid to human serum albumin. <i>Protein Science</i> 30(4):830-841. doi:10.1002/pro.4036",
"zurlinden2025": "Zurlinden TJ, Dzierlenga MW, Kapraun DF, Ring C, Bernstein AS, Schlosser PM, Morozov V (2025). Estimation of species- and sex-specific PFAS pharmacokinetics in mice, rats, and non-human primates using a Bayesian hierarchical methodology. <i>Toxicology and Applied Pharmacology</i> 499:117336. doi:10.1016/j.taap.2025.117336. Data: github.com/USEPA/CPHEA-Animal-PFAS-PK",
"weaver2010": "Weaver YM, Ehresman DJ, Butenhoff JL, Hagenbuch B (2010). Roles of rat renal organic anion transporters in transporting perfluorinated carboxylates with different chain lengths. <i>Toxicological Sciences</i> 113(2):305-314. doi:10.1093/toxsci/kfp275",
"zhao2017": "Zhao W, Zitzow JD, Weaver Y, Ehresman DJ, Chang SC, Butenhoff JL, Hagenbuch B (2017). Organic anion transporting polypeptides contribute to the disposition of perfluoroalkyl acids in humans and rats. <i>Toxicological Sciences</i> 156(1):84-95. doi:10.1093/toxsci/kfw236",
"abraham2024": "Abraham K, Mertens H, Richter L, Mielke H, et al. (2024). Kinetics of 15 PFAS after a single oral dose in one adult volunteer: terminal half-lives, clearances and derived volumes of distribution. <i>Environment International</i>. doi:10.1016/j.envint.2024.109047",
"fischer2024": "Fischer FC, et al. (2024). Protein binding of PFAS measured by solid-phase microextraction at environmentally relevant PFAS:protein ratios. <i>Environmental Science &amp; Technology</i>. doi:10.1021/acs.est.3c07415, and Fischer FC, et al. (2025). doi:10.1021/acs.est.5c05473",
"andersson2025": "Andersson AG, et al. (2025). The relative importance of fecal and urinary excretion of perfluorooctane sulfonic acid and perfluorooctanoic acid after high exposure: an observational study in Ronneby, Sweden. <i>Environmental Research</i> 285:122487. doi:10.1016/j.envres.2025.122487",
"li2022": "Li Y, Andersson A, Xu Y, Pineda D, Nilsson CA, Lindh CH, Jakobsson K, Fletcher T (2022). Determinants of serum half-lives for perfluoroalkyl substances after end of exposure to contaminated drinking water, Ronneby cohort. (Held as a structured abstract; full text paywalled.)",
"chiu2022": "Chiu WA, et al. (2022). Bayesian estimation of human population toxicokinetics of PFOA, PFOS, PFHxS and PFNA from studies of contaminated drinking water. <i>Environmental Health Perspectives</i> 130(12):127001. doi:10.1289/EHP10103",
"butenhoff2004": "Butenhoff JL, Kennedy GL, Hinderliter PM, Lieder PH, Jung R, Hansen KJ, Gorman GS, Noker PE, Thomford PJ (2004). Pharmacokinetics of perfluorooctanoate in cynomolgus monkeys. <i>Toxicological Sciences</i> 82(2):394-406. doi:10.1093/toxsci/kfh302",
"fujii2015": "Fujii Y, Niisoe T, Harada KH, Uemoto S, Ogura Y, Takenaka K, Koizumi A (2015). Toxicokinetics of perfluoroalkyl carboxylic acids with different carbon chain lengths in mice and humans. <i>Journal of Occupational Health</i> 57(1):1-12. doi:10.1539/joh.14-0136-OA. (Values reached via OEHHA 2024 PHG Table 4.8.1 and ATSDR 2021 Table 3-6.)",
"starnes2024": "Starnes HM, Jackson TW, Rock KD, Belcher SM (2024). Quantitative cross-species comparison of serum albumin binding of per- and polyfluoroalkyl substances from five structural classes. <i>Toxicological Sciences</i> 199(1):132-149. doi:10.1093/toxsci/kfae028",
"gallagher2023": "Gallagher A, Kar S, Sep&uacute;lveda MS (2023). Computational modeling of human serum albumin binding of per- and polyfluoroalkyl substances employing QSAR, read-across, and docking. <i>Molecules</i> 28(14):5375. doi:10.3390/molecules28145375. (Computational; its affinity data are the [starnes2024] set.)",
"gao2025": "Gao Y, Kibbey TCG, Lee MJ, O'Carroll DM (2025). Predicting salt effects on the binding of PFAS to protein. <i>Environmental Science &amp; Technology</i> 59(49):26716-26728. doi:10.1021/acs.est.5c08779. (Abstract only; per-compound affinities not obtained.)",
"zhao2023": "Zhao L, Teng M, Zhao X, Li Y, Sun J, Zhao W, Ruan Y, Leung KMY, Wu F (2023). Insight into the binding model of per- and polyfluoroalkyl substances to proteins and membranes. <i>Environment International</i> 175:107951. doi:10.1016/j.envint.2023.107951. (Review.)",
"chang2008": "Chang S-C, Das K, Ehresman DJ, Ellefson ME, Gorman GS, Hart JA, Noker PE, Tan Y-M, Lieder PH, Lau C, Olsen GW, Butenhoff JL (2008). Comparative pharmacokinetics of perfluorobutyrate in rats, mice, monkeys, and humans and relevance to human exposure via drinking water. <i>Toxicological Sciences</i> 104(1):40-53. doi:10.1093/toxsci/kfn057",
"chang2012": "Chang S-C, Noker PE, Gorman GS, Gibson SJ, Hart JA, Ehresman DJ, Butenhoff JL (2012). Comparative pharmacokinetics of perfluorooctanesulfonate (PFOS) in rats, mice, and monkeys. <i>Reproductive Toxicology</i> 33(4):428-440. doi:10.1016/j.reprotox.2011.07.002",
"dzierlenga2020": "Dzierlenga AL, Robinson VG, Waidyanatha S, DeVito MJ, Eifrid MA, Gibbs ST, Granville CA, Blystone CR (2020). Toxicokinetics of PFHxA, PFOA and PFDA in male and female Hsd:Sprague Dawley SD rats following intravenous or gavage administration. <i>Xenobiotica</i> 50(6):722-732. doi:10.1080/00498254.2019.1683776",
"kim2016": "Kim S-J, Heo S-H, Lee D-S, Hwang IG, Lee Y-B, Cho H-Y (2016). Gender differences in pharmacokinetics and tissue distribution of 3 perfluoroalkyl and polyfluoroalkyl substances in rats. <i>Food and Chemical Toxicology</i> 97:243-255. doi:10.1016/j.fct.2016.09.017",
"iwabuchi2017": "Iwabuchi K, Senzaki N, Mazawa D, Sato I, Hara M, Ueda F, Liu W, Tsuda S (2017). Tissue toxicokinetics of perfluoro compounds with single and chronic low doses in male rats. <i>Journal of Toxicological Sciences</i> 42(3):301-317. doi:10.2131/jts.42.301",
"kemper2003": "Kemper RA (2003). Perfluorooctanoic acid: toxicokinetics in the rat. Unpublished report, DuPont-7473; US EPA public docket AR-226-1499. Haskell Laboratory, E.I. du Pont de Nemours. (Not obtained; its data are used as digitised by the EPA database [zurlinden2025].)",
"bartell2010": "Bartell SM, Calafat AM, Lyu C, Kato K, Ryan PB, Steenland K (2010). Rate of decline in serum PFOA concentrations after granular activated carbon filtration at two public water systems in Ohio and West Virginia. <i>Environmental Health Perspectives</i> 118(2):222-228. doi:10.1289/ehp.0901252",
"gasiorowski2022": "Gasiorowski R, Forbes MK, Silver G, Krastev Y, Hamdorf B, Lewis B, et al. (2022). Effect of plasma and blood donations on levels of perfluoroalkyl and polyfluoroalkyl substances in firefighters in Australia: a randomized clinical trial. <i>JAMA Network Open</i> 5(4):e226257. doi:10.1001/jamanetworkopen.2022.6257",
"harada2007": "Harada K, et al. (2007). Biliary excretion and enterohepatic recirculation of perfluorooctane sulfonate in humans. Values used here are as tabulated by US EPA (2016), <i>Health Effects Support Document for Perfluorooctane Sulfonate</i>; the original was not obtained.",
"delaere2025": "Delaere I, Harris K, Gaskin S, Tefera Y, Mitchell K, Springer D, Mills S (2025). Changes in serum perfluorooctane sulfonic acid and perfluorohexane sulfonic acid concentrations in firefighters accessing a voluntary PFAS reduction treatment program. <i>Environment International</i> 202:109609. doi:10.1016/j.envint.2025.109609",
"genuis2010": "Genuis SJ, et al. (2010). Human elimination of perfluorinated compounds under cholestyramine administration. (Cited via [andersson2025] and the EPA assessment; the original was not obtained.)",
"moller2024": "M&oslash;ller S, et al. (2024). Cholestyramine cross-over trial for PFAS lowering. <i>Environment International</i>. doi:10.1016/j.envint.2024.108471 (cited via [andersson2025])",
"zhang2013": "Zhang Y, Beesoon S, Zhu L, Martin JW (2013). Biomonitoring of perfluoroalkyl acids in human urine and estimates of biological half-life. <i>Environmental Science &amp; Technology</i> 47(18):10619-10627. doi:10.1021/es401905e",
"yi2022": "Yi S, et al. (2022). Biotransformation of 6:2 chlorinated polyfluoroalkyl ether sulfonate in the rat. (Held as a structured abstract; see db/primary_2026/yi2022_clpfesa_rat.csv.)",
"epa2024pfoa": "US EPA (2024). <i>Final Human Health Toxicity Assessment for Perfluorooctanoic Acid (PFOA)</i>, EPA-815R24006, and its Appendix (Table B-26).",
"epa2024pfos": "US EPA (2024). <i>Final Human Health Toxicity Assessment for Perfluorooctane Sulfonic Acid (PFOS)</i>.",
"epa2025pfhxs": "US EPA (2025). <i>IRIS Toxicological Review of Perfluorohexanesulfonic Acid (PFHxS)</i>, Table 3-3.",
"epa2023pfhxa": "US EPA (2023). <i>IRIS Toxicological Review of Perfluorohexanoic Acid (PFHxA)</i>.",
"oehha2024": "California OEHHA (2024). <i>Public Health Goals for PFOA and PFOS in Drinking Water</i>, Appendix Tables A6.3, A6.4 and 4.8.1.",
"atsdr2021": "ATSDR (2021). <i>Toxicological Profile for Perfluoroalkyls</i>, Tables 3-5 and 3-6.",
"efsa2020": "EFSA CONTAM Panel (2020). Risk to human health related to the presence of perfluoroalkyl substances in food. <i>EFSA Journal</i> 18(9):6223. doi:10.2903/j.efsa.2020.6223",
"njdwqi2018": "New Jersey Drinking Water Quality Institute (2017-2018). Health-based maximum contaminant level support documents for PFOA, PFOS and PFNA.",
}

# ---------------------------------------------------------------------------
# citation machinery: numbers assigned by order of first appearance
# ---------------------------------------------------------------------------
_ORDER = []


def c(*keys):
    """Render an in-text citation marker for one or more reference keys."""
    nums = []
    for k in keys:
        if k not in REFS:
            raise KeyError(f"undefined reference key: {k}")
        if k not in _ORDER:
            _ORDER.append(k)
        nums.append(str(_ORDER.index(k) + 1))
    return "[" + ", ".join(nums) + "]"


# ---------------------------------------------------------------------------
# styles
# ---------------------------------------------------------------------------
def _styles():
    ss = getSampleStyleSheet()
    base = dict(fontName="Times-Roman", textColor=INK)
    return {
        "title": ParagraphStyle("t", parent=ss["Title"], fontName="Times-Bold",
                                fontSize=22, leading=26, textColor=INK,
                                alignment=0, spaceAfter=4),
        "subtitle": ParagraphStyle("st", fontName="Times-Roman", fontSize=12.5,
                                   leading=16.5, textColor=SOFT, spaceAfter=16),
        "h1": ParagraphStyle("h1", fontName="Times-Bold", fontSize=14,
                             leading=16.5, textColor=INK, spaceBefore=15,
                             spaceAfter=6),
        "h2": ParagraphStyle("h2", fontName="Times-Bold", fontSize=10.8,
                             leading=13.5, textColor=INK, spaceBefore=11,
                             spaceAfter=4),
        "body": ParagraphStyle("b", fontSize=9.6, leading=12.8,
                               alignment=TA_JUSTIFY, spaceAfter=6, **base),
        "lead": ParagraphStyle("l", fontSize=10.6, leading=14.6,
                               alignment=TA_JUSTIFY, spaceAfter=7, **base),
        "abs": ParagraphStyle("ab", fontSize=9.8, leading=13.2,
                              alignment=TA_JUSTIFY, spaceAfter=6,
                              leftIndent=6 * mm, rightIndent=6 * mm, **base),
        "caption": ParagraphStyle("c", fontName="Times-Italic", fontSize=8.2,
                                  leading=10.8, textColor=SOFT, spaceBefore=3,
                                  spaceAfter=9),
        "eq": ParagraphStyle("e", fontName="Times-Italic", fontSize=11.5,
                             leading=15, textColor=INK, alignment=1,
                             spaceBefore=5, spaceAfter=7),
        # Reference entries are set ragged-right, not justified. A DOI is an
        # unbreakable token, so justifying pushes it to the next line and
        # stretches the one before it into wide gaps -- which is what
        # happened to the Kudo 2001 and Moller 2024 entries.
        "ref": ParagraphStyle("r", fontSize=8.5, leading=10.8,
                              spaceAfter=3.5, leftIndent=7 * mm,
                              firstLineIndent=-7 * mm, **base),
        "cell": ParagraphStyle("tc", fontName="Times-Roman", fontSize=8.3,
                               leading=10.2, textColor=INK),
        "cellh": ParagraphStyle("th", fontName="Times-Bold", fontSize=8.1,
                                leading=10, textColor=INK),
    }


S = _styles()
W = A4[0] - 40 * mm


def p(t, k="body"):
    return Paragraph(t, S[k])


def tbl(rows, widths, head=True):
    data = []
    for i, r in enumerate(rows):
        st = "cellh" if (head and i == 0) else "cell"
        data.append([x if isinstance(x, Paragraph) else Paragraph(str(x), S[st])
                     for x in r])
    t = Table(data, colWidths=widths, hAlign="LEFT",
              repeatRows=1 if head else 0)
    cmds = [("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("LINEBELOW", (0, 0), (-1, -2), 0.3, RULE),
            ("LINEBELOW", (0, -1), (-1, -1), 0.9, INK)]
    if head:
        cmds += [("LINEABOVE", (0, 0), (-1, 0), 0.9, INK),
                 ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
                 ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f4f3f1"))]
    t.setStyle(TableStyle(cmds))
    return t


# ---------------------------------------------------------------------------
# numbers read from the database at build time
# ---------------------------------------------------------------------------
def _rows(path):
    return list(csv.DictReader(open(os.path.join(HERE, path), newline="")))


def qsar_rows():
    r = [x for x in _rows("db/qsar/qsar_endpoint_table.csv")
         if x["log10R_gfr_argoul"]]
    r.sort(key=lambda x: float(x["log10R_gfr_argoul"]))
    return r


def axis_rows():
    return _rows("db/primary_2026/han2012_table4_reabsorption_axis.csv")


def inventory():
    import glob
    comb = os.path.join(HERE, "db", "combined")
    return {
        "texts": len(glob.glob(os.path.join(HERE, "papers", "*.txt"))),
        "extracts": len(glob.glob(os.path.join(HERE, "db", "primary_2026",
                                               "*.csv"))),
        "figures": len(glob.glob(os.path.join(HERE, "figures", "*.png"))),
        "scripts": len(glob.glob(os.path.join(HERE, "scripts", "*.py"))),
        "rows": sum(len(list(csv.DictReader(open(os.path.join(comb, f)))))
                    for f in sorted(os.listdir(comb)) if f.endswith(".csv")),
        "report_lines": sum(1 for _ in open(os.path.join(HERE, "report",
                                                         "REPORT.md"))),
    }


# ---------------------------------------------------------------------------
# page furniture
# ---------------------------------------------------------------------------
def _footer(canv, doc):
    canv.saveState()
    canv.setStrokeColor(RULE)
    canv.setLineWidth(0.4)
    canv.line(20 * mm, 14 * mm, A4[0] - 20 * mm, 14 * mm)
    canv.setFont("Times-Roman", 7.6)
    canv.setFillColor(SOFT)
    canv.drawString(20 * mm, 10 * mm,
                    "Why PFAS half-lives differ between humans, rats and mice")
    canv.drawRightString(A4[0] - 20 * mm, 10 * mm, str(doc.page))
    canv.restoreState()


def _cover(canv, doc):
    canv.saveState()
    canv.setFillColor(ACCENT)
    canv.rect(0, A4[1] - 11 * mm, A4[0], 11 * mm, stroke=0, fill=1)
    canv.restoreState()


# ---------------------------------------------------------------------------
# the document
# ---------------------------------------------------------------------------
def story():
    inv = inventory()
    out = []

    # ===================== cover and abstract =====================
    out += [
        Spacer(1, 22 * mm),
        p("Why PFAS half-lives differ between humans, rats and mice", "title"),
        p("A review of the evidence, the assumptions behind the published "
          "numbers, and what is still unresolved", "subtitle"),
        tbl([
            ["Evidence base", f"{inv['texts']} full texts read; "
                              f"{inv['extracts']} per-paper extractions; "
                              f"{inv['rows']:,} database rows; "
                              f"{len(REFS)} references"],
            ["Underlying review", f"report/REPORT.md, ~{inv['report_lines']:,} "
                                  f"lines, {inv['scripts']} runnable scripts, "
                                  f"{inv['figures']} figures"],
            ["Repository", "github.com/sdey17/pfas_tk"],
        ], [32 * mm, W - 32 * mm], head=False),
        Spacer(1, 5 * mm),
        KeepTogether([
        p("Contents", "h1"),
        tbl([
            ["1", "Introduction: one identity, three questions"],
            ["2", "Evidence base and method"],
            ["3", "Where the species difference lives: distribution or "
                  "clearance?"],
            ["4", "The mechanism: fractional renal reabsorption"],
            ["5", "Transporter identity, and why the mechanism does not "
                  "transfer"],
            ["6", "The second loop: enterohepatic recirculation"],
            ["7", "Plasma protein binding, and a conflict that was not one"],
            ["8", "Does exposure level change elimination?"],
            ["9", "The human volume of distribution"],
            ["10", "Toward structure-activity: the renal handling ratio"],
            ["11", "What this review corrects in the published record"],
            ["12", "Coverage: how much of the field is empty"],
            ["13", "Open questions, and what would close them"],
            ["14", "Conclusions"],
            ["", "References"],
        ], [8 * mm, W - 8 * mm], head=False)]),
        NextPageTemplate("body"),
        PageBreak(),
        p("Abstract", "h1"),
        p("Published serum half-lives for the same per- and polyfluoroalkyl "
          "substance differ by orders of magnitude between humans and "
          "laboratory animals, between rat and mouse, between male and female, "
          "and sometimes between papers describing the same experiment. This "
          "review asks where that variance actually sits. Because half-life is "
          "not a property but a ratio - the volume a chemical apparently "
          "occupies divided by the rate the body clears it - the question is "
          "answerable, and the answer is unambiguous: across six datasets that "
          "measure both terms in both sexes within single experiments, the "
          "volume of distribution varies 1.6-fold while clearance varies "
          "56-fold and changes sign between species. The variance is in "
          "clearance.", "abs"),
        p("One mechanistic axis then orders every species for which data "
          "exist. Fractional renal reabsorption runs from 99.94% in humans "
          "through 97% in the male mouse and 93.7% in the male rat to net "
          "tubular secretion in the female rat and the rabbit. Rodent "
          "half-lives follow from it within 1.7-fold over a 37-fold span, "
          "though that agreement establishes that renal clearance is nearly "
          "all of total clearance in rodents rather than validating an "
          "independent model. Two refinements matter. The axis is two "
          "factors rather than one: of the 333-fold male-mouse-to-human renal "
          "clearance gap, the escape fraction carries 50-fold and the "
          "six-fold lower human glomerular filtration rate carries 6.5-fold. "
          "And humans run a second, enterohepatic loop of comparable "
          "completeness, which closes the residual the renal axis leaves and "
          "explains why bile-acid sequestrants shorten human half-lives "
          "several-fold.", "abs"),
        p("The review also audits the provenance of the numbers regulation "
          "uses. The human volume of distribution behind nine adopted "
          "clearance factors is a calculation rather than a measurement, "
          "reproducible to three figures from its source paper's own "
          "supplementary table given an assumed elimination rate drawn from a "
          "study in the same two communities; and a widely-inherited plasma "
          "free fraction is roughly 164-fold too high because a calculated "
          "lower bound was read as a measurement. Eighteen such corrections "
          "are documented, one of them to this work's own method.", "abs"),
        p("Finally, half-life is argued to be the wrong endpoint for any "
          "structure-activity model, because it carries a body-size term. A "
          "dimensionless renal handling ratio is proposed in its place and "
          "evaluated on the only dataset that supports the comparison - a "
          "cocktail study whose absolute clearances and volumes both run "
          "two- to threefold below the one-compound-at-a-time literature while "
          "their ratio does not, a shared offset that leaves every "
          "between-compound contrast the endpoint is built on intact.", "abs"),
    ]

    # ===================== 1. introduction =====================
    out += [
        p("1. Introduction: one identity, three questions", "h1"),
        p("Per- and polyfluoroalkyl substances persist in the body for a long "
          "time, and how long is the single number that drives their "
          "regulation: it converts an external dose into an internal one, and "
          "it is what any cross-species extrapolation has to carry. The "
          "published values do not agree. Human PFOA half-life estimates span "
          "roughly seventeen-fold, from about half a year to eight and a half; "
          "the male rat and the female rat differ by about seventy-fold in the "
          "same experiment " + c("kudo2002", "ohmori2003") + "; and the rat "
          "and the mouse, two rodents of similar size, differ in a direction "
          "that reverses between sexes " + c("tatum2011") + ".", "lead"),
        p("Disagreements of that size usually mean a quantity is being "
          "compared across things that are not comparable. The way into the "
          "problem is that half-life is not an independent property. It is "
          "fixed by an identity:", "body"),
        p("t<sub>1/2</sub> = ln2 &middot; V<sub>d</sub> / CL", "eq"),
        p("where V<sub>d</sub> is the apparent volume of distribution - the "
          "number relating total body burden to plasma concentration, not an "
          "anatomical volume - and CL is clearance, the volume of plasma "
          "irreversibly cleared per unit time. Any two of the three terms fix "
          "the third. A study that measures serum decay obtains the rate "
          "constant directly and never needs a volume; a study that works from "
          "mass balance needs a volume and a complete excretion accounting, "
          "each uncertain severalfold, and both multiply into the answer. The "
          "practical consequence runs through this entire review: most "
          "published disagreements are not disagreements about measurements. "
          "They are differences in <b>which two terms were measured and which "
          "one was assumed</b> - and in the most consequential case, the "
          "assumed one is the number everybody quotes.", "body"),
        p("Three questions follow, and this review answers them in order. "
          "First, is the species difference a difference in how much is held "
          "or in how much is let go - distribution or clearance? Second, what "
          "did each reported half-life assume? Third, does exposure level "
          "itself change elimination, as saturable kinetics would predict?",
          "body"),
    ]

    # ===================== 2. evidence base =====================
    out += [
        p("2. Evidence base and method", "h1"),
        p(f"The review rests on {inv['texts']} full texts read directly, "
          f"{inv['extracts']} per-paper extractions made from those texts, and "
          f"a consolidated database of {inv['rows']:,} rows across four "
          "schemas - toxicokinetic parameters, protein binding, transporter "
          "kinetics, and adopted regulatory values. Animal serum time courses "
          "come from the US EPA's compiled per-record database "
          + c("zurlinden2025") + ", which digitised or transcribed them from "
          "the primary studies - principally the monkey intravenous series "
          + c("butenhoff2004", "chang2012") + ", the rodent dose series "
          + c("dzierlenga2020", "huang2019", "kemper2003") + " and a "
          "low-dose chronic study " + c("iwabuchi2017") + ". Human values "
          "come from the primary papers and from agency assessments "
          + c("epa2024pfoa", "epa2024pfos", "epa2023pfhxa", "oehha2024",
              "atsdr2021", "efsa2020") + ".", "body"),
        p("Three rules govern what is in the database, and they are the reason "
          "several of this review's conclusions differ from the secondary "
          "literature.", "body"),
        p("<b>Read the table, not the sentence.</b> Where a paper's abstract "
          "and its own tables disagree, the table wins and the discrepancy is "
          "recorded. Two of the largest corrections in section 11 are of "
          "exactly this kind, and so is the correction this review makes to "
          "its own earlier work.", "body"),
        p("<b>Mark what was computed.</b> Every row records whether a value "
          "was measured, fitted, assumed, or computed in this work, and "
          "carries a PMID or DOI and the table it came from. A clearance "
          "derived from an assumed volume is not the same kind of object as a "
          "clearance measured from urine, and conflating the two is how the "
          "problem in section 9 propagated.", "body"),
        p("<b>Flag your own biases rather than deleting them.</b> This review "
          "refitted the EPA raw curves and found its own fitting procedure "
          "biased on biphasic data. The affected column is kept, flagged, and "
          "documented, because it remains valid for the within-curve "
          "comparisons it is used for - and because deleting it would hide a "
          "failure mode that applies to much of the published literature too.",
          "body"),
    ]
    return out


def story_2():
    out = []

    # ===================== 3. where the difference lives =====================
    out += [
        KeepTogether([
        p("3. Where the species difference lives: distribution or clearance?",
          "h1"),
        p("The identity in section 1 makes this testable rather than "
          "rhetorical, provided both terms are measured in the same animals. "
          "Six datasets satisfy that condition - three species, two compounds, "
          "four laboratories, no cross-study comparison required "
          + c("kudo2002", "ohmori2003", "sundstrom2012", "kim2016",
              "dzierlenga2020", "huang2019") + ".", "body"),
        tbl([
            ["", "V<sub>d</sub> male/female", "clearance female/male"],
            ["span across all six datasets", "<b>1.33 - 2.18&times;</b>",
             "<b>0.78 - 44.3&times;</b>"],
            ["", "a 1.6&times; spread", "a 56&times; spread"],
        ], [62 * mm, (W - 62 * mm) / 2, (W - 62 * mm) / 2])]),
        Spacer(1, 3),
        p("The distribution ratio is male-higher in every dataset and never "
          "leaves a narrow band, which is what one would expect of a quantity "
          "set mostly by body composition and plasma protein concentration. "
          "The clearance ratio ranges over nearly two orders of magnitude "
          "<b>and changes sign between species</b>. That asymmetry is the "
          "central finding of this review, and everything after it is an "
          "attempt to explain the clearance term.", "body"),
        p("A single recent experiment reproduces the result internally. "
          "Argoul and colleagues dosed eleven PFAS as one cocktail in female "
          "mice and fitted all of them simultaneously " + c("argoul2026")
          + "; across those eleven compounds, clearance spans 5,254-fold while "
          "steady-state volume of distribution spans 7.7-fold. Same animals, "
          "same assay, same model. Note that this varies <i>chemistry</i> "
          "at fixed species and sex, where the table above varies sex at "
          "fixed chemistry, so it is corroboration by analogy rather than "
          "replication - though it is hard to explain why clearance should be "
          "the volatile term along both axes unless it simply is the volatile "
          "term. This study's <i>absolute</i> clearances run about 2.5&times; "
          "below the one-compound-at-a-time literature, for reasons set out "
          "in section 10; the comparison here is a ratio within one "
          "experiment, from which a shared scaling cancels, so the "
          "5,254-fold against 7.7-fold is unaffected.", "body"),
        p("Taken from primary sources rather than compilations, the species "
          "gap is also concentrated in one sex. Female rat divided by female "
          "mouse PFOA clearance is 373-496&times;, against 7.0&times; in males "
          + c("kudo2002", "lou2009") + " - this comparison <i>is</i> "
          "cross-study, rat and mouse from different laboratories, and "
          "carries the usual caveat. Tatum-Gibbs and colleagues "
          "replicate the whole structure in a second compound with strains "
          "matched: a rat sex ratio of 21.9&times; against a mouse ratio near "
          "unity " + c("tatum2011") + ". Any account of 'the species "
          "difference' that does not mention sex is describing an average over "
          "two quite different animals.", "body"),
    ]

    # ===================== 4. the axis =====================
    ax = axis_rows()
    rows = [["species", "sex", "GFR (L/d/kg)", "renal CL (mL/d/kg)",
             "reabsorbed"]]
    for sp in ("human", "mouse", "rat", "Japanese macaque", "dog", "rabbit"):
        for r in ax:
            if r["species"] != sp:
                continue
            pct = r["pct_reabsorption"]
            rows.append([r["species"], r["sex"], r["gfr_L_d_kg"],
                         r["clr_mL_d_kg"],
                         f"<b>{pct}%</b>" if pct else "<b>net secretion</b>"])
    out += [
        p("4. The mechanism: fractional renal reabsorption", "h1"),
        p("PFAS are filtered freely at the glomerulus to the extent they are "
          "unbound, and then substantially recovered from the tubular "
          "filtrate. Renal clearance is therefore what escapes:", "body"),
        p("CL<sub>renal</sub> = f<sub>u</sub> &middot; GFR &middot; (1 - FR)",
          "eq"),
        p("where f<sub>u</sub> is the unbound fraction in plasma and FR the "
          "fraction reabsorbed. Han and colleagues assembled this across "
          "species " + c("han2012") + "; the table below is read from that "
          "paper's Table 4 at source rather than from any adaptation of it.",
          "body"),
        KeepTogether([
            tbl(rows, [32 * mm, 17 * mm, 23 * mm, 32 * mm, W - 104 * mm]),
            p("Read at source into "
              "db/primary_2026/han2012_table4_reabsorption_axis.csv. Two "
              "values differ from a widely-cited adaptation of this table; see "
              "section 11.", "caption"),
        ]),
        p("The ordering is the half-life ordering. A caveat on how much that "
          "demonstrates, because it is easy to overclaim: the reabsorbed "
          "fraction is <i>derived</i> from the measured renal clearance by "
          "rearranging the equation above, so predicting a half-life from it "
          "restates the identity t<sub>1/2</sub> = ln2 &middot; V<sub>d</sub> "
          "/ CL<sub>renal</sub> rather than testing a reabsorption model "
          "independently. What the agreement to 1.7&times; over a 37&times; "
          "span does establish is that <b>renal clearance accounts for nearly "
          "all of total clearance in these rodents</b> - a real and "
          "non-trivial result, and the reason the human residual in section 6 "
          "stands out so sharply. The axis orders the species and names the "
          "quantity that varies; it does not independently predict them.",
          "body"),
        p("The reabsorbed fraction does get one independent check. A separate "
          "dataset, different laboratory and different animals, recomputes "
          "mouse PFOA reabsorption at 95.9% " + c("argoul2026") + ", inside "
          "the 95.2-97.0% range already in use.", "body"),
        p("Elimination hands off from kidney to gut as the chain "
          "lengthens", "h2"),
        p("The axis is renal, and for the shorter carboxylates that is nearly "
          "the whole story: PFHpA leaves the male rat 92% in urine over five "
          "days. By PFNA that falls to 2.0% and by PFDA to 0.2%, with faeces "
          "becoming the major route " + c("kudo2001") + ". So a renal "
          "clearance is close to a total clearance at C7 and is a small "
          "minority of it at C10, and any human accounting that assigns one "
          "route fraction across compounds is wrong in a predictable "
          "direction. The handoff is also sex-dependent: female rats excrete "
          "51% of a PFNA dose in urine where males excrete 2.0% "
          + c("kudo2001") + ", which is the same sex difference as section 3 "
          "seen through the route split rather than the rate.", "body"),
        p("The axis has two factors, not one", "h2"),
        p("It is tempting to read the table as 'humans reabsorb more', and "
          "this review initially did. The absolute numbers say otherwise: "
          "humans recover 51 mL/d/kg of filtrate, against 270 in the male rat "
          "and 318-324 in the mouse. Humans reabsorb <i>less</i> in absolute "
          "terms and still clear far more slowly, because they filter far less "
          "to begin with. Decomposing the 333&times; male-mouse-to-human renal "
          "clearance gap gives 50&times; from the escape fraction and "
          "6.5&times; from the six-fold lower human GFR - 68% and 32% on a log "
          "scale. Reabsorption is the larger term, which is why the "
          "single-axis framing works at all; but a third of the difference is "
          "filtration rate, and no transporter story explains that third.",
          "body"),
    ]

    # ===================== 5. transporters =====================
    out += [
        KeepTogether([
        p("Figure 2", "h2"),
        Image(FIG2, width=W, height=W / 2.565),
        p("<b>A:</b> 21 matched female/male pairs - eight compounds across "
          "three species - with <i>both</i> ratios written in the same "
          "direction, so the contrast cannot be an artefact of how each was "
          "expressed. The volume of distribution stays inside a 4.7&times; "
          "band straddling 1.0; clearance spans 103&times; and crosses it in "
          "both directions. This is a wider set than the six datasets in the "
          "table above, and it says the same thing more strongly. <b>B:</b> "
          "the male-mouse-to-human renal clearance gap split into its two "
          "multiplicative factors; because they multiply, their logarithms "
          "add, and the escape fraction takes 68% of the distance. "
          "Regenerate with scripts/make_review_figures.py.", "caption")]),
        KeepTogether([
        p("5. Transporter identity, and why the mechanism does not transfer",
          "h1"),
        p("If reabsorption is the controlling step, some apical transporter "
          "must perform it, and in the rat the sex difference constrains which "
          "one. A candidate has to satisfy two conditions simultaneously: it "
          "must be sex-divergent in the right direction, and it must actually "
          "carry PFOA. Those two conditions are jointly much more restrictive "
          "than either alone. The table covers the transporters for which "
          "both properties have been measured, so it is an elimination among "
          "tested candidates rather than a proof that no untested protein "
          "qualifies.", "body"),
        tbl([
            ["transporter", "sex-divergent?", "transports PFOA?", "verdict"],
            ["Oatp1a1", "<b>yes</b>, 23&times; male-predominant, "
                        "androgen-induced", "<b>yes</b>, C8-C10",
             "<b>the only candidate satisfying both</b>"],
            ["Oat2", "yes, strongly", "<b>no</b> - three negative reports, "
                     "two species", "excluded by transport"],
            ["Oat1 / Oat3", "no", "yes", "excluded by regulation"],
        ], [26 * mm, 50 * mm, 46 * mm, W - 122 * mm])]),
        Spacer(1, 3),
        p("The transport evidence is from direct uptake assays "
          + c("yang2009", "weaver2010", "louisse2024") + " and the regulation "
          "evidence from expression studies under hormonal manipulation "
          + c("kudo2002", "cheng2005", "cheng2009", "buist2004") + "; the "
          "organic anion transporting polypeptides as a family have been "
          "shown to contribute to PFAA disposition in both humans and rats "
          + c("zhao2017") + ". The "
          "23&times; figure is from the primary source; a frequently quoted "
          "'2.5-fold' belongs to a different transporter (see section 11).",
          "body"),
        p("The surprise is in the human", "h2"),
        p("The reabsorbed <i>fraction</i> transfers across species - that is "
          "what section 4 shows. The <i>mechanism</i> does not. OATP1A2, the "
          "closest human orthologue of rat Oatp1a1, does not mediate saturable "
          "PFOA uptake at all " + c("yang2010") + ". Human apical reabsorption "
          "instead runs through OAT4 and URAT1 " + c("yang2010", "louisse2023")
          + ", neither of which is androgen-regulated - which is consistent "
          "with the absence in humans of anything like the rat's seventyfold "
          "sex difference. Extrapolating a rat mechanism to humans and "
          "extrapolating a rat reabsorbed fraction to humans are therefore "
          "very different acts, and only the second is supported.", "body"),
    ]

    # ===================== 6. the second loop =====================
    out += [
        p("6. The second loop: enterohepatic recirculation", "h1"),
        p("The renal axis leaves a residual. Predicted human half-life comes "
          "out about 4.3&times; short of observation - visible as the single "
          "off-line point in Figure 1. Three candidate explanations exist: the "
          "human reabsorbed fraction is underestimated, the human clearance "
          "figure is wrong, or a second elimination route is being "
          "reabsorbed too.", "body"),
        p("The third is correct, and the number was already in an agency "
          "document. Harada and colleagues sampled serum and bile from four "
          "gallstone-surgery patients and measured a biliary resorption rate "
          "of 0.97 " + c("harada2007") + ". Bile PFOS (27.9 ng/mL) actually "
          "<i>exceeds</i> serum (23.2), so bile is a genuine excretion route - "
          "but 97% of what is secreted is recovered from the gut.", "body"),
        tbl([
            ["loop", "fraction reabsorbed", "source"],
            ["renal (PFOA)", "<b>99.94%</b>", "Han 2012 Table 4 "
                                              + c("han2012")],
            ["biliary (PFOS)", "<b>97%</b>", "Harada 2007, n = 4 "
                                             + c("harada2007")],
        ], [38 * mm, 38 * mm, W - 76 * mm]),
        Spacer(1, 3),
        p("Humans therefore run two near-complete reabsorption loops rather "
          "than one. Both are near unity, both lengthen the half-life, and - "
          "the testable part - <b>both are interruptible</b>. If a 97% "
          "resorption loop is blocked, elimination should accelerate several "
          "fold, and it does. In a treatment programme using cholestyramine "
          "and plasma donation, apparent PFOS half-life was 1.2 y in nineteen "
          "treated participants against 7.3 y in nine observed ones, and "
          "PFHxS 2.5 y against 9.4 y " + c("delaere2025") + ". Supporting "
          "evidence comes from a cholestyramine series " + c("genuis2010")
          + " and a cross-over trial reporting 63% lowering in twelve weeks "
          "against 3% in controls " + c("moller2024") + ". Plasma donation "
          "alone, which removes burden without touching either loop, produces "
          "a much smaller effect " + c("gasiorowski2022") + ".", "body"),
        p("The loop also settles an apparent contradiction in the human "
          "excretion literature. One study measuring faeces under ongoing "
          "intake finds faecal elimination dominant for PFOS "
          + c("andersson2025") + "; another, following a single labelled "
          "dose, barely detects it " + c("abraham2024") + ". At 97% "
          "resorption, <b>gross</b> biliary flux is large while <b>net</b> "
          "faecal elimination is small, so the two studies measured different "
          "quantities and both are right. The caveat is that the constant "
          "doing this work rests on four patients.", "body"),
    ]
    return out


def story_3():
    out = []

    # ===================== 7. binding =====================
    out += [
        p("7. Plasma protein binding, and a conflict that was not one", "h1"),
        p("The unbound fraction f<sub>u</sub> enters the reabsorption equation "
          "directly, so its value propagates into everything above. The "
          "literature appeared to contain a hundredfold disagreement: a "
          "long-standing result reports PFOA 'over 90% bound' to plasma "
          "protein " + c("han2003") + ", implying f<sub>u</sub> near 0.1, "
          "while recent measurements at environmentally realistic "
          "ligand:protein ratios give 0.00061 " + c("fischer2024") + ".",
          "body"),
        p("Reading both papers at source dissolves the conflict. The '>90%' "
          "is not a measurement but a <i>calculation</i> from a dissociation "
          "constant and an albumin concentration, and as a lower bound it is "
          "satisfied by 0.00061 as comfortably as by 0.1. The underlying "
          "difference in dissociation constant is a titration artefact: the "
          "older work titrated 50-60 &micro;M albumin with 0.1-3 mM PFOA, a "
          "ligand:protein ratio of 1.7:1 to 60:1, against 0.004:1 or below in "
          "the newer work - while human serum sits at 10<super>-5</super> to "
          "10<super>-3</super>:1. Near-saturating ratios understate affinity, "
          "and structural work on the albumin complex is consistent with "
          "multiple sites of differing affinity " + c("maso2021") + ".",
          "body"),
        p("The defect was in the inheritance rather than in either "
          "measurement. Physiologically based models read '>90% bound' as "
          "'about 90% bound' and adopted a free fraction roughly "
          "<b>164&times; too high</b>. Section 10 shows how far that single "
          "substitution propagates.", "body"),
        p("Species and structure, from the studies that hold method fixed",
          "h2"),
        p("Most of the 285-row binding file cannot answer a comparative "
          "question, for the reason above: a single affinity number per pair "
          "does not exist. What can are the few studies that fix the method "
          "and vary one thing.", "body"),
        KeepTogether([
        p("On <b>species</b>, " + c("starnes2024") + " measured four albumins "
          "by one assay. Dissociation constant in mM, lower is tighter:",
          "body"),
        tbl([
            ["compound", "F-carbons", "human", "rat", "bovine", "porcine",
             "rat / human"],
            ["PFBA", "3", "2.64", "1.90", "3.06", "1.28", "0.72&times;"],
            ["PFBS", "4", "1.68", "1.41", "0.80", "1.25", "0.84&times;"],
            ["GenX", "4", "1.57", "1.64", "2.55", "2.12", "1.04&times;"],
            ["PFHxA", "5", "1.71", "1.33", "2.17", "0.95", "0.78&times;"],
            ["6:2 FTSA", "6", "0.41", "0.33", "1.39", "0.52", "0.80&times;"],
            ["PFOA", "7", "0.79", "0.40", "0.99", "0.37", "<b>0.51&times;</b>"],
            ["PFOS", "8", "0.69", "0.42", "0.86", "0.38", "0.61&times;"],
        ], [26 * mm, 20 * mm, 18 * mm, 18 * mm, 19 * mm, 20 * mm,
            W - 121 * mm]),
        p("Built by scripts/binding_species_and_structure.py from "
          "db/protein_binding.csv.", "caption")]),
        p("The study reports the interspecies differences as statistically "
          "significant for every bound compound, and they are - but the "
          "magnitudes are small. Across all four species the highest-to-lowest "
          "ratio never exceeds about fourfold and is usually under twofold. "
          "Rat and porcine albumin bind the carboxylates and sulfonates more "
          "tightly than human, the ether acid less tightly, in the study's own "
          "words and in the table. Albumin is a conserved protein and behaves "
          "like one. The direction is again unhelpful for explaining kinetics: "
          "rat albumin holds PFOA twice as tightly as human while the rat "
          "clears it about a thousandfold faster. A second method agrees - "
          + c("fischer2024") + " measured human and bovine albumin for 27 "
          "compounds at a physiological ligand:protein ratio and the "
          "bovine-minus-human difference in log D is under 0.3 log units for "
          "most of them, median -0.14. Bovine albumin is a fair stand-in for "
          "human, which matters because much of the literature uses it.",
          "body"),
        p("On <b>structure</b>, " + c("fischer2024") + "'s ten chain lengths "
          "give three results, and only the first is the textbook one. "
          "<b>Chain length is not monotonic.</b> The carboxylates rise from "
          "2.16 at four fluorinated carbons to a maximum of 5.18 around ten to "
          "eleven, then fall to 3.87 at thirteen - a slope of +0.45 log units "
          "per CF<sub>2</sub> on the way up and -0.44 on the way down. A "
          "structure-activity model linear in carbon number is wrong at both "
          "ends. An independent dataset reaches the same shape: the QSAR of "
          + c("gallagher2023") + " over the " + c("starnes2024") + " affinities "
          "finds perfluoroalkyl acids of seven to nine carbons bind human "
          "albumin most tightly and those above ten bind less, and its docking "
          "gives a mechanism - the tightest binders form a hydrogen bond with "
          "Arg348 and salt bridges to Arg348 and Arg485, while the weakest "
          "form neither. <b>The head group is worth about a log unit</b> at "
          "matched chain length, with similar slopes rather than a different "
          "trend: PFHxS 4.77 against PFHpA 4.10 at six carbons, PFHpS 5.49 "
          "against PFOA 4.48 at seven. <b>And albumin stops being the main "
          "carrier.</b> &gamma;-globulin binds two to three log units weaker "
          "at short chain lengths, closes the gap steadily, and at thirteen "
          "carbons binds more tightly than albumin does. Treating protein "
          "binding as albumin binding degrades exactly where it matters, since "
          "the long chains are the persistent ones. " + c("zhao2023")
          + " reviews the same picture mechanistically: the fluorinated chain "
          "occupies the binding cavity while the acid head group hydrogen-"
          "bonds to residues, so chlorine substitution, ether insertion and "
          "branching all act through chain length and hydrophobicity rather "
          "than independently.", "body"),
        p("<b>The gap is in the field, not in this file.</b> There is no in "
          "vitro measurement of PFAS binding to <i>mouse</i> serum albumin. A "
          "PubMed search returns one record, and it measures mouse albumin "
          "concentration rather than binding. So the species question is "
          "answerable for human, rat, bovine and porcine and not for the "
          "mouse - the animal whose clearance and unbound fraction section 10's "
          "endpoint is built from, and whose f<sub>u</sub> therefore has no in "
          "vitro albumin measurement to decompose it against. The one new "
          "experimental series found, " + c("gao2025") + " on bovine albumin "
          "across chain length and head group, is logged in the database "
          "without values because its full text was not obtainable.", "body"),
        p("What binding does to elimination, quantified", "h2"),
        p("Only unbound PFAS is filtered, so the clearance free filtration "
          "alone would produce is f<sub>u</sub>&middot;GFR, and the measured "
          "renal gap between species factorises into binding, body size and "
          "transport. For PFOA the female-mouse-to-human renal clearance gap "
          "is 120&times; whatever is assumed; what the assumption changes is "
          "how that gap is apportioned.", "body"),
        KeepTogether([
        tbl([
            ["human f<sub>u</sub>", "binding", "body size", "transport",
             "dominant term"],
            ["0.00061 measured " + c("fischer2024"), "<b>14.3&times;</b>",
             "3.9&times;", "2.2&times;", "<b>binding</b>"],
            ["0.02 assumed " + c("han2012"), "0.4&times;", "3.9&times;",
             "70.7&times;", "transport"],
            ["0.10 inherited", "0.1&times;", "3.9&times;", "353.5&times;",
             "transport"],
        ], [W - 108 * mm, 24 * mm, 24 * mm, 24 * mm, 36 * mm]),
        p("Mouse limb from " + c("argoul2026") + " (CL<sub>renal</sub> 3.59 "
          "mL/d/kg, f<sub>u</sub> 0.87%), human from " + c("han2012")
          + " Table 4 (0.03 mL/d/kg). Built by "
          "scripts/binding_elimination_decomposition.py.", "caption")]),
        p("Under the measured value <b>binding is the largest single term</b>: "
          "humans hold PFOA about fourteenfold more tightly than mice, which "
          "by itself keeps fourteenfold less of it in front of the glomerulus. "
          "Under the assumed values the apportionment inverts - binding runs "
          "the wrong way, implying humans bind less tightly than mice, and "
          "transport has to carry 71 to 354-fold to balance. The mouse and "
          "human free fractions come from different methods at different "
          "ligand:protein ratios, which is the artefact this section "
          "diagnoses, so the fourteenfold is the least secure number in the "
          "table; and the mouse limb carries the cocktail study's low "
          "clearance bias, which lands on the transport term and would raise "
          "2.2&times; to about 5&times;, still far below 71&times;.", "body"),
        p("Binding alone does not explain long human half-lives, and the "
          "numbers say why. At f<sub>u</sub> = 0.00061 free filtration would "
          "clear 1.57 mL/d/kg, which with V<sub>d</sub> = 170 mL/kg is a "
          "half-life of <b>75 days</b> - weeks, not years. Measured human "
          "renal clearance is 0.03 mL/d/kg, <b>52&times; lower still</b>, "
          "because tubular reabsorption recovers <b>98.1%</b> of what is "
          "filtered. Binding sets how much reaches the glomerulus; "
          "reabsorption disposes of almost all of it. Binding is the "
          "species-variable factor and reabsorption the near-complete one in "
          "every species measured, so a model with binding right and "
          "reabsorption wrong still misses the half-life by two orders of "
          "magnitude.", "body"),
    ]

    # ===================== 8. dose =====================
    out += [
        p("8. Does exposure level change elimination?", "h1"),
        p("If reabsorptive transport saturates, elimination should accelerate "
          "at high burden, half-life should fall with dose, and a single "
          "clearance factor could not be transferred between a contaminated "
          "community and the general population. The evidence is weak but not "
          "null, and it splits informatively by the level at which the "
          "comparison is made.", "body"),
        tbl([
            ["comparison", "direction"],
            ["within a person over time", "supports dose dependence"],
            ["within a species across doses (rat refits)",
             "supports, weakly: slope about +0.11"],
            ["between people in a cohort",
             "<b>contradicts</b> - the lowest-exposure tertile declines fastest"],
            ["between water districts",
             "slope above the mechanistic ceiling of 1, i.e. bias"],
        ], [62 * mm, W - 62 * mm]),
        Spacer(1, 3),
        p("Within-unit comparisons support concentration dependence and "
          "between-unit ones do not, which is the signature of a between-unit "
          "confounder - and the source study names it: age " + c("li2022")
          + ". Re-running that paper's own variance decomposition, age carries "
          "2-13&times; the partial R<super>2</super> of initial PFAS "
          "concentration, and the exposure-tertile effect is 25-37% the size "
          "of the age effect. Older people eliminate more slowly and have "
          "higher accumulated burdens, so the between-person association is a "
          "positive saturation effect minus a larger negative age effect.",
          "body"),
        p("Across volume of distribution the slopes run -0.21 to +0.25 with "
          "inconsistent sign, i.e. no relationship. The one clean exception "
          "runs the other way: the female rat shows a slope of +0.21 over a "
          "3,200&times; dose range against -0.02 in the male, which is what "
          "saturable <i>secretion</i> looks like - consistent with the female "
          "rat sitting at the secretory end of the axis in section 4.", "body"),
        p("Whether human exposures reach saturation is unresolved and is taken "
          "up in section 13.", "body"),
    ]

    # ===================== 9. the human Vd =====================
    out += [
        p("9. The human volume of distribution", "h1"),
        p("Nine of the 42 adopted regulatory clearance factors catalogued here "
          "are computed from a single human volume of distribution, 170 mL/kg "
          + c("thompson2010") + ". That number is not a measurement.", "body"),
        p("Its source paper's supplementary table is headed <i>'Input data for "
          "the calibration of the Vd parameter'</i> with a final column headed "
          "<i>'calculated Vd'</i>, and its published values reproduce to three "
          "figures from intake, serum concentration and an assumed elimination "
          "rate. The assumed rate was taken from a study that measured "
          "half-life in <b>the same two communities</b> " + c("bartell2010")
          + ", which is what makes the construction circular rather than "
          "merely indirect.", "body"),
        p("The consequence is sharper than circularity, and it is "
          "asymmetric. Forming a clearance from the identity makes the "
          "half-life cancel exactly:", "body"),
        tbl([
            ["half-life assumed", "V<sub>d</sub> that follows",
             "derived clearance"],
            ["2.3 y (the value used)", "<b>168 mL/kg</b>",
             "0.132-0.138 mL/kg-day"],
            ["3.8 y", "<b>277 mL/kg</b>",
             "0.132-0.138 mL/kg-day (<i>unchanged</i>)"],
        ], [46 * mm, 38 * mm, W - 84 * mm]),
        Spacer(1, 3),
        p("Clearance reduces to dose divided by serum concentration and is "
          "robust; the volume does not cancel and is proportional to whatever "
          "half-life was assumed. So adopting '170' alongside a different "
          "half-life silently contradicts the data it came from. The "
          "arithmetic is reproduced from the paper's own supplementary table "
          "to within 1.000-1.003&times; in scripts/thompson_vd_circularity.py. "
          "The paper's own corrigendum is an instance of exactly this error, "
          "corrected by a factor of 0.6 - the ratio of two rate constants.",
          "body"),
        p("Independent measurements do not support the adopted value in either "
          "direction: direct estimates cluster at 74, 121 and 113-199 mL/kg "
          + c("andersson2025", "abraham2024", "gasiorowski2022") + ", all "
          "<i>below</i> the assigned range, while the one population "
          "pharmacokinetic fit gives 430 mL/kg " + c("chiu2022") + ", far "
          "above it. A mass-balance half-life anchored on 170 mL/kg "
          + c("zhang2013") + " moves to within 5% of the population fit when "
          "the fitted volume is substituted instead. That is one study, so "
          "it does not show the whole seventeen-fold span to be a single "
          "assumed constant; what it shows is that the mass-balance end of "
          "the span - the end that needs a volume at all - collapses toward "
          "the serum-decay end when the volume is changed.", "body"),
    ]

    # ===================== 10. the handling ratio =====================
    q = qsar_rows()
    qr = [["compound", "F-carbons", "head group", "ether O", "f<sub>u</sub> %",
           "log<sub>10</sub> R", "handling"]]
    for r in q:
        qr.append([
            r["chemical"], r["n_fluorinated_c"],
            {"carboxylate": "-COOH", "sulfonate": "-SO<sub>3</sub>H",
             "ether-carboxylate": "-COOH, ether"}[r["head_group"]],
            r["n_ether_o"], r["fu_pct"], r["log10R_gfr_argoul"],
            "<b>SECRETED</b>" if float(r["log10R_gfr_argoul"]) > 0
            else "reabsorbed"])
    out += [
        p("10. Toward structure-activity: the renal handling ratio", "h1"),
        p("A long-term aim of this work is to relate PFAS toxicokinetics to "
          "molecular structure. Half-life cannot be that endpoint. Because it "
          "is ln2&middot;V<sub>d</sub>/CL, it carries a glomerular filtration "
          "term that differs 6.5&times; between mouse and human with no change "
          "in chemistry whatever; a structural model fitted to it is being "
          "asked to absorb body size, and it cannot.", "body"),
        p("The replacement proposed here is dimensionless - measured renal "
          "clearance divided by the clearance free filtration alone would "
          "produce:", "body"),
        p("R = CL<sub>renal</sub> / (f<sub>u</sub> &middot; GFR)", "eq"),
        p("R below 1 is net reabsorption, R above 1 net secretion, R = 1 pure "
          "filtration. It divides out GFR, so species of different size become "
          "comparable, and divides out f<sub>u</sub>, so the binding step is "
          "not counted twice; what remains is the transport step, which is the "
          "part a structural model could plausibly learn. It is the axis of "
          "section 4 re-expressed, since R = 1 - FR, on a scale that neither "
          "crowds against a ceiling at 1 nor runs unboundedly negative under "
          "secretion.", "body"),
        KeepTogether([
                p("Only one dataset supports the comparison: nine compounds with both "
          "renal clearance and unbound fraction measured in one experiment, "
          "one species, one sex and one laboratory - female mice "
          + c("argoul2026") + ". Note that this is a different animal from the "
          "male mouse of section 4's table, whose reabsorption figures come "
          "from a different source " + c("han2012") + ". R spans "
          "875&times; across them.", "body"),
            tbl(qr, [21 * mm, 16 * mm, 22 * mm, 14 * mm, 14 * mm, 18 * mm,
                     W - 105 * mm]),
            p("Built by scripts/qsar_endpoint_table.py into "
              "db/qsar/qsar_endpoint_table.csv. 'F-carbons' counts carbons "
              "bearing fluorine, which excludes a carboxylate's acid carbon "
              "and includes every carbon of a sulfonate - a definition chosen "
              "so that PFNA and PFOS come out the same size, which the data "
              "can then falsify.", "caption"),
        ]),
        p("Three results constrain any structure-activity model built on this "
          "endpoint. <b>Chain length alone does not order it</b>: within the "
          "carboxylates the series is non-monotonic and PFHxA crosses into net "
          "secretion, so a model using carbon number as its sole descriptor is "
          "already falsified on nine compounds. <b>The head group carries "
          "about threefold at matched chain length</b> - PFNA against PFOS, "
          "both with eight fluorinated carbons, differ 2.9&times;, so the "
          "equivalence the descriptor assumed is a good approximation but not "
          "free. <b>Ether oxygens move the endpoint 1.4 log units at constant "
          "chain length, and non-monotonically</b>: PFHxA with none is "
          "secreted, GenX with one is reabsorbed, PFO2OA with two is secreted. "
          "The replacement chemicals sit on both sides of the divide, which is "
          "a regulatory observation as much as a chemical one.", "body"),
        p("How far the one dataset can be trusted", "h2"),
        p("Everything above rests on a single experiment, and that experiment "
          "is unusual in three ways at once: eleven PFAS dosed as one cocktail, "
          "at 0.019-1.55 mg/kg, fitted by nonlinear mixed effects. Each of "
          "those could bias the absolute numbers, so the comparison against "
          "the one-compound-at-a-time literature is worth making explicitly. "
          "All comparators below are female mice, matching "
          + c("argoul2026") + ".", "body"),
        KeepTogether([
        p("Figure 4", "h2"),
        Image(FIG4, width=W, height=W / 3.170),
        p("<b>A:</b> every single-compound female-mouse study that reports "
          "both clearance and volume, each on its own row, against the "
          "cocktail. The grey connector joins the two parameters within a "
          "study. Geometric means 0.40&times; for clearance "
          + c("fujii2015", "sundstrom2012")
          + " and 0.54&times; for volume. <b>B:</b> the same comparison for "
          "half-life, which is their ratio. Seven comparisons across three "
          "compounds " + c("chang2008", "sundstrom2012", "chang2012")
          + ", spanning 0.79&times; to 1.30&times;, geometric mean "
          "1.04&times;. Only comparisons with one "
          "attributable dose are shown: three studies pooled two dose levels "
          "into a single reported parameter " + c("lou2009", "tatum2011",
                                                  "epa2023pfhxa")
          + " and are excluded, which costs PFOA, PFNA and PFHxA. Including "
          "them gives ten half-life comparisons at 0.99&times; rather than "
          "seven at 1.04&times;, so nothing turns on the exclusion. Built by "
          "scripts/argoul_vs_single_compound.py into "
          "db/argoul_vs_single_compound.csv; drawn by "
          "scripts/make_review_figures.py.", "caption")]),
        p("One caveat on how much of a check the clearance and volume rows "
          "are. In most of these studies the volume is not measured apart "
          "from the clearance: non-compartmental analysis gives "
          "V<sub>z</sub> = CL/&lambda;<sub>z</sub> and V<sub>ss</sub> = "
          "CL&middot;MRT. For the two Sundstr&ouml;m rows 'clearance "
          "low and volume low' is largely one observation seen twice. "
          "Fujii is the exception - volume from dose/C(0) by "
          "back-extrapolation, clearance from dose/AUC - and there the two "
          "are low independently, 0.38&times; and 0.59&times;. One row is "
          "thin evidence, so the claim defended here is the narrower one: "
          "the rate constant agrees across the half-life comparisons and "
          "the scale "
          "does not.", "body"),
        p("The shape of that result is more informative than either half of "
          "it. Clearance and volume are both low by about twofold, and the "
          "half-life is not low at all - which is what a shared multiplicative "
          "factor on both terms looks like, because in "
          "ln2&middot;V<sub>d</sub>/CL such a factor cancels. Three candidates "
          "fit, and testing the most obvious one weakens it. <b>Dose</b> looks "
          "the natural candidate - section 8's saturable-reabsorption account "
          "predicts clearance rising with dose, and the cocktail was dosed "
          "below its comparators - but Argoul tuned the dose per compound "
          "(0.019-1.55 mg/kg, tracking persistence rather than sweeping a "
          "range), so the gap must be matched compound by compound, and then "
          "it fails to order the discrepancies. Fujii dosed PFOA at 0.128 "
          "mg/kg against the cocktail's 0.080, a gap of 1.6&times;, and "
          "disagrees <i>most</i> (0.38&times;); Lou at 1 mg/kg and above, a "
          "gap of at least 12&times;, agrees <i>best</i> (0.75&times;). "
          "Dropping the cocktail from the PFHxS dose series and extrapolating "
          "Sundstr&ouml;m's own two doses down to it predicts 2.04 mL/kg/d "
          "against 1.30 observed, so dose closes part of that gap and about "
          "1.6&times; survives it. Dose is a contributor, not the "
          "explanation. "
          "<b>Cocktail competition</b> for reabsorptive transporters would "
          "raise clearance, not lower it, so it argues against the observed "
          "direction - unless the competition is for plasma binding sites, "
          "which would raise f<sub>u</sub> and lower the apparent volume. "
          "<b>Mixed-effects shrinkage</b> compresses spread rather than "
          "shifting the centre, so it is the weakest of the three.", "body"),
        p("What this does and does not cost the endpoint above. R is built "
          "from an absolute renal clearance, so a 2.5-fold low bias in that "
          "clearance is a 2.5-fold low bias in R - a flat -0.40 log unit offset "
          "on every compound. Because it is flat, it changes nothing that this "
          "section actually claims: the ranking, the 2.9-log-unit span, the "
          "non-monotonicity in chain length and the head-group and ether "
          "contrasts are all differences between compounds, and a shared "
          "factor cancels from every one of them. It does move the compounds "
          "relative to the fixed divide at R = 1, but not far enough to matter "
          "- the nearest reabsorbed compound sits 0.98 log units below the "
          "divide, so the bias would have to be 9.6&times; rather than "
          "2.5&times; to reclassify anything. What it does cost is the "
          "cross-species comparison in the next subsection, where an absolute "
          "mouse value is set against an absolute human one; read those gaps "
          "as carrying a factor of about 2.5 on top of the f<sub>u</sub> "
          "uncertainty, in the direction of understating the mouse.", "body"),
        p("A sensitivity that bears on the species question", "h2"),
        KeepTogether([
        p("R is linear in 1/f<sub>u</sub>, so it inherits the uncertainty of "
          "section 7 in full. For human PFOA, with renal clearance and GFR "
          "both fixed, the three free fractions in circulation give three "
          "different answers:", "body"),
        tbl([
            ["f<sub>u</sub>", "source", "log<sub>10</sub> R",
             "gap to female mouse"],
            ["0.10", "models reading '>90% bound' as 'about 90%'", "-3.93",
             "354&times;"],
            ["0.02", "the stated assumption in " + c("han2012"), "-3.23",
             "71&times;"],
            ["0.00061", "measured at physiological ligand:protein "
                        + c("fischer2024"), "-1.72", "<b>2.2&times;</b>"],
        ], [18 * mm, W - 78 * mm, 22 * mm, 38 * mm])]),
        Spacer(1, 3),
        p("Under the assumed value the female-mouse-to-human difference in the "
          "transport step is about seventyfold and the species gap is "
          "transport biology. Under the measured value it is about twofold, "
          "and the species gap is almost entirely <b>binding</b> - which would "
          "mean the transporter literature of section 5 is explaining a "
          "quantity that barely differs between the species. This is stated as "
          "a sensitivity rather than a result, because the mouse and human "
          "free fractions come from different methods at different "
          "ligand:protein ratios, which is precisely the artefact section 7 "
          "diagnoses. It is also why the experiment proposed in section 13 "
          "is the one worth doing first. To size the problem plainly: the "
          "endpoint spans 2.9 log units across the nine compounds a "
          "structure-activity model would be fitted to, and the disagreement "
          "over this one input moves a single compound by 2.2 of them. The "
          "endpoint is well posed; the inputs are not yet good enough to fit "
          "it.", "body"),
    ]
    return out


def story_4():
    inv = inventory()
    out = []

    # ===================== 11. corrections =====================
    out += [
        KeepTogether([
        p("Figure 3", "h2"),
        Image(FIG3, width=W, height=W / 2.702),
        p("<b>A:</b> the nine compounds on the proposed endpoint, ordered by "
          "it, coloured by head group and labelled with fluorinated-carbon "
          "count. The vertical rule at zero is the reabsorption/secretion "
          "divide. The bracket marks the three compounds with five "
          "fluorinated carbons, which span 1.4 log units - the clearest "
          "single refutation of chain length as a sole descriptor. <b>B:</b> "
          "human PFOA on the same endpoint under each of the three unbound "
          "fractions in circulation, against the female mouse. The quantity a "
          "structure-activity model would have to explain spans 2.9 log units "
          "across panel A; the disagreement about one input moves one "
          "compound by 2.2 of them. Regenerate with "
          "scripts/make_review_figures.py.", "caption")]),
        p("11. What this review corrects in the published record", "h1"),
        p("Each entry is traced to the table or figure that contradicts it; "
          "the full list of eighteen is section 8 of report/REPORT.md.",
          "body"),
        tbl([
            ["source", "correction"],
            ["US EPA " + c("epa2024pfoa"),
             "the rat Oatp1a1 male/female ratio quoted as '2.5-fold' is a "
             "<i>different transporter's</i> number (OAT-K); the primary value "
             "is <b>23&times;</b> " + c("kudo2002")],
            ["OEHHA " + c("oehha2024"),
             "its adaptation of the reabsorption table altered two values "
             "(human 99.94 to 99.8, male rat 93.7 to 93.2) relative to the "
             "source " + c("han2012")],
            ["Provenance",
             "the f<sub>u</sub> = 0.02 behind the axis is the source paper's "
             "assumption, not the agency's, and that paper calls its own "
             "values 'rough estimates' " + c("han2012")],
            ["Andersson 2025 " + c("andersson2025"),
             "transposes the PFOA and PFOS volumes of " + c("thompson2010")
             + " (it is PFOA 170, PFOS 230); an earlier claim in this work "
             "that " + c("zhang2013") + " did so was itself wrong and is "
             "withdrawn"],
            ["Cheng 2006 " + c("cheng2006"),
             "its abstract calls renal Oatp1a1 'female-predominant', then "
             "reports androgens <i>increase</i> it and concludes androgens are "
             "the exclusive cause - self-contradictory, and against "
             + c("cheng2005", "cheng2009")],
            ["Metabolic inertness",
             "'PFAS are metabolically inert' fails for 6:2 Cl-PFESA, which is "
             "biotransformed " + c("yi2022", "shi2016")],
            ["<b>This review</b>",
             "<b>its own refits of the EPA raw curves compress sex ratios "
             "about 13&times;</b>, as below"],
        ], [34 * mm, W - 34 * mm]),
        p("The self-correction matters most", "h2"),
        p("This work refitted terminal slopes directly from the EPA raw curves "
          + c("zurlinden2025") + ", selecting the window with the best "
          "adjusted R<super>2</super>. On a biphasic curve that criterion "
          "selects the shallow terminal tail. Tested against a published value "
          "for the same experiment, the fit returns <b>1.44 d against 0.08 "
          "d</b> for the female rat and compresses a <b>71&times; sex ratio to "
          "5.6&times;</b> " + c("kudo2002") + ". The bias runs one way, so "
          "every conclusion above survives and several strengthen - but that "
          "column must not be read as comparable to published half-lives. It "
          "carries a flag and is kept rather than deleted, because it remains "
          "valid for within-curve comparison and because the same failure mode "
          "applies wherever a terminal slope is fitted without inspecting the "
          "window.", "body"),
    ]

    # ===================== 12. coverage =====================
    out += [
        p("12. Coverage: how much of the field is empty", "h1"),
        KeepTogether([
        p("A review that only reports what is known overstates the state of "
          "the field, so the gap was measured. Crossing 39 chemicals by 11 "
          "species by 4 parameters gives 1,716 cells.", "body"),
        tbl([
            ["", "cells"],
            ["no data at all", "<b>1,466 (85%)</b>"],
            ["resting on a single study", "94"],
            ["PFOS, the best-covered compound in the world", "29 of 44 filled"],
            ["by species: human / rat / mouse / monkey / dog",
             "71 / 61 / 45 / 28 / <b>2</b>"],
        ], [68 * mm, W - 68 * mm])]),
        Spacer(1, 3),
        p("Mining four further compilations " + c("atsdr2021", "efsa2020",
                                                  "epa2025pfhxs",
                                                  "njdwqi2018")
          + " moved coverage from 143 to 156 cells, so the hole is real rather "
          "than a search artefact. The dog column deserves particular notice: "
          "it is the species that sits furthest down the reabsorption axis "
          "among the mammals with data, and it rests on two cells, one of "
          "which traces to a book chapter that could not be obtained.", "body"),
    ]

    # ===================== 13. open questions =====================
    out += [
        p("13. Open questions, and what would close them", "h1"),
        p("13.1  The unbound fraction, measured one way across species", "h2"),
        p("This is the cheapest decisive experiment the review identifies, and "
          "it is first because two separate questions turn on it. Every value "
          "of R in section 10 is proportional to f<sub>u</sub>, and the mouse "
          "and human numbers currently come from different methods at "
          "different ligand:protein ratios. Measuring them by one method at "
          "physiological ratio decides whether the species difference is "
          "binding or transport - a 33&times; swing in the answer - and "
          "simultaneously supplies the denominator any structure-activity "
          "model needs. It requires plasma, a dialysis or ultrafiltration rig "
          "and a mass spectrometer. No animals.", "body"),
        p("13.2  A measured human transport-affinity constant", "h2"),
        p("Whether reabsorptive transport saturates at real human exposures "
          "depends on which parameter family one believes. Six independent in "
          "vitro half-saturation constants for PFOA against human transporters "
          "span 47 to 310 &micro;M, from two laboratories, agreeing within a "
          "factor of seven " + c("yang2010", "louisse2024", "louisse2023")
          + ". The transport-affinity constant fitted inside physiologically "
          "based models is 0.133 &micro;M - <b>354 to 2,336&times; lower than "
          "anything ever measured in a cell</b>. The in vitro values put human "
          "serum far below half-saturation; the fitted value puts three of "
          "four human populations at or above it, which would make clearance "
          "dose-dependent in contaminated communities and mean a single "
          "clearance factor cannot transfer between exposure settings. The "
          "weight of evidence has moved onto the in vitro side - the fitted "
          "values are fitted to plasma curves rather than measured, the mouse "
          "row carries standard errors exceeding its estimates, and the "
          "dose-response evidence of section 8 agrees with the in vitro answer "
          "- but nobody has measured a human value, and one measurement would "
          "close it outright.", "body"),
        p("13.3  Whether the biliary resorption constant replicates", "h2"),
        p("The 0.97 of section 6 is now load-bearing for the human limb, and "
          "it rests on four surgical patients in a single study reached "
          "through an agency document " + c("harada2007") + ". A replication "
          "in any species under a stated protocol would be worth more than its "
          "cost.", "body"),
        p("13.4  The mouse renal Oatp1a1 sex ratio", "h2"),
        p("The rat value is 23&times; " + c("kudo2002") + "; the mouse "
          "equivalent is not established, and the mouse sex difference runs "
          "the opposite way to the rat's " + c("tatum2011") + ". One "
          "quantitative PCR experiment would discriminate between the live "
          "hypotheses - though if the answer to 13.1 is that the species "
          "gap is binding, this question loses most of its force, which is why "
          "it is fourth rather than first.", "body"),
        p("13.5  Unobtained sources", "h2"),
        p("Three documents would materially change specific numbers and could "
          "not be obtained: the dog toxicokinetic chapter underlying the dog "
          "column, an unpublished contract report on protein binding, and the "
          "primary rat dose series " + c("kemper2003") + ", which supplies "
          "roughly 70% of the rat PFOA observations used here only through its "
          "digitisation by " + c("zurlinden2025") + ".", "body"),
    ]

    # ===================== 14. conclusions =====================
    out += [
        p("14. Conclusions", "h1"),
        p("The species difference in PFAS elimination is a clearance "
          "difference, not a distribution difference, and it is concentrated "
          "in females. One mechanistic axis - the fraction of filtered "
          "compound recovered from the tubule - orders every species for which "
          "both terms have been measured, and accounts for rodent half-lives "
          "within 1.7&times;. That axis is two "
          "factors rather than one, roughly two thirds reabsorption and one "
          "third glomerular filtration rate, and in humans it is joined by a "
          "second, enterohepatic loop of comparable completeness whose "
          "interruption shortens human half-lives several-fold.", "body"),
        p("The reabsorbed fraction transfers between species; the protein "
          "performing the reabsorption does not. That distinction matters for "
          "how animal data should be read across to humans: the quantitative "
          "axis is extrapolable, the mechanism is not.", "body"),
        p("Against that, two of the most widely inherited human parameters do "
          "not mean what their users take them to mean. The volume of "
          "distribution behind nine adopted clearance factors is a calculation "
          "whose assumed half-life came from the same communities it was "
          "calibrated on, and a plasma free fraction in general use is roughly "
          "164&times; too high because a calculated lower bound was read as a "
          "measurement. Neither is a disputed measurement; both are provenance "
          "failures, and both are fixable by reading the source tables.",
          "body"),
        p("Finally, half-life should not be the endpoint of a structure-"
          "activity model, because it carries body size. The dimensionless "
          "renal handling ratio proposed here removes both the filtration and "
          "the binding terms, spans 875&times; across the nine compounds that "
          "support the comparison, and shows immediately that chain length "
          "alone will not do: the head group carries threefold at matched "
          "length and ether substitution moves the endpoint 1.4 log units "
          "non-monotonically, placing the replacement chemicals on both sides "
          "of the reabsorption-secretion divide. Extending that comparison to "
          "a second species, with unbound fraction measured the same way in "
          "both, is the single step that would turn this argument into a "
          "model.", "body"),
        p("Data and code availability", "h2"),
        p(f"All extractions, the consolidated {inv['rows']:,}-row database, "
          f"the {inv['scripts']} analysis scripts that regenerate every figure "
          "and table above, and the full review are at "
          "github.com/sdey17/pfas_tk. Every database row carries a PMID or "
          "DOI, the table it was read from, and whether its value was "
          "measured, fitted, assumed or computed in this work. Animal serum "
          "time courses are redistributed from the EPA compilation "
          + c("zurlinden2025") + " under its MIT licence.", "body"),
    ]

    # ===================== references =====================
    out += [p("References", "h1")]
    out += [p("Numbered by order of first citation. Entries marked as held "
              "through a secondary source were not obtained in the original; "
              "the review does not quote them as if they had been.",
              "caption")]
    # A few entries point at another reference ("cited via [andersson2025]").
    # Those are keys, not display text, so resolve them to their numbers now
    # that the order is final -- otherwise the reader sees an internal key.
    def resolve(text):
        def sub(m):
            key = m.group(1)
            if key not in _ORDER:
                raise KeyError(f"reference text points at uncited key {key!r}")
            return f"[{_ORDER.index(key) + 1}]"
        return re.sub(r"\[([a-z]+\d{4}[a-z]*)\]", sub, text)

    for i, k in enumerate(_ORDER, 1):
        out.append(p(f"<b>{i}.</b>&nbsp;&nbsp;{resolve(REFS[k])}", "ref"))
    unused = [k for k in REFS if k not in _ORDER]
    if unused:
        raise SystemExit(f"unused references (remove or cite): {unused}")
    return out


# ---------------------------------------------------------------------------
# figure, assembly, validation
# ---------------------------------------------------------------------------
def figure_block():
    return [
        p("Figure 1", "h1"),
        Image(FIG, width=W, height=W / 2.033),
        p("<b>The argument in one pair of panels.</b> <b>Left:</b> across five "
          "species-by-sex groups in which both terms were measured, the volume "
          "of distribution spans 2.6&times; while renal clearance spans "
          "13,958&times; - the variance is in clearance, not distribution "
          "(section 3). <b>Right:</b> half-life predicted from fractional "
          "renal reabsorption alone, against observation, on log axes with the "
          "1:1 line. Every rodent point falls within 1.7&times; over a "
          "37&times; span - which, as section 4 notes, shows that renal "
          "clearance is nearly all of total clearance in rodents rather than "
          "validating an independent model, since the reabsorbed fraction is "
          "derived from the same renal clearance. The human "
          "point sits 4.3&times; high, and the walk-down beneath it - renal "
          "only 4.3&times;, plus faecal 2.7&times;, EPA clearance 2.2&times;, "
          "OEHHA measured clearance 0.9&times; - is what led to the second "
          "reabsorption loop (section 6). <b>Note</b> that the panel annotations carry the widely-circulated OEHHA adaptation of the reabsorption table (male rat 93.2%, human 99.8%), while the table in section 4 carries the source values (93.7% and 99.94%). The discrepancy is not an error here - it is one of the corrections listed in section 11, left visible rather than silently harmonised. Regenerate with scripts/make_master_figure.py.", "caption"),
    ]


def check_glyphs(flowables):
    """reportlab's built-in fonts are WinAnsi; anything else is a black box."""
    bad = {}

    def scan(txt):
        for ch in txt:
            if (ord(ch) < 32 or ord(ch) > 126) and ch not in SAFE_EXTRA:
                bad[ch] = bad.get(ch, 0) + 1

    def walk(f):
        for obj in (f._content if isinstance(f, KeepTogether) else [f]):
            if isinstance(obj, Paragraph):
                scan(obj.text)
            elif isinstance(obj, Table):
                for row in obj._cellvalues:
                    for cell in row:
                        if isinstance(cell, Paragraph):
                            scan(cell.text)

    for f in flowables:
        walk(f)
    return bad


def main():
    # order matters: citation numbers are assigned as the text is built
    flow = story() + figure_block() + story_2() + story_3() + story_4()

    bad = check_glyphs(flow)
    if bad:
        print("non-WinAnsi characters (these render as black boxes):")
        for ch, n in sorted(bad.items(), key=lambda kv: -kv[1]):
            print(f"   U+{ord(ch):04X} {ch!r} x{n}")
        sys.exit(1)

    doc = BaseDocTemplate(
        OUT, pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
        topMargin=18 * mm, bottomMargin=20 * mm,
        title="Why PFAS half-lives differ between humans, rats and mice",
        author="pfas_tk", subject="PFAS toxicokinetics review")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height,
                  id="f")
    doc.addPageTemplates([
        PageTemplate(id="cover", frames=[frame], onPage=_cover),
        PageTemplate(id="body", frames=[frame], onPage=_footer),
    ])
    doc.build(flow)

    try:
        from pypdf import PdfReader
        pages = len(PdfReader(OUT).pages)
    except ImportError:
        pages = "?"
    print(f"wrote {os.path.relpath(OUT, HERE)}  {pages} pages, "
          f"{os.path.getsize(OUT)/1024:.0f} KB, "
          f"{len(_ORDER)} references cited")


if __name__ == "__main__":
    main()
