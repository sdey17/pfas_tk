#!/usr/bin/env python3
"""Render one of this repository's markdown reports to PDF.

Matches the house style of make_summary_pdf.py (Times body, A4, 20 mm
margins, the same ink/accent colours) so the documents sit together.

Handles the constructs actually used in report/*.md: ATX headings, bold /
italic / inline code / links, pipe tables, fenced code blocks, blockquotes,
bullet and numbered lists, horizontal rules, and <details>/<summary> blocks.

Usage:
    python3 tk_review/scripts/md_to_pdf.py report/ALBUMIN_BINDING.md [out.pdf]
"""
import html
import os
import re
import sys

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, KeepTogether,
                                PageTemplate, Paragraph, Spacer, Table,
                                TableStyle)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

INK = colors.HexColor("#111111")
SOFT = colors.HexColor("#55534f")
RULE = colors.HexColor("#d8d7d3")
ACCENT = colors.HexColor("#2a78d6")
CODEBG = colors.HexColor("#f4f3f0")
QUOTEBG = colors.HexColor("#f7f6f3")

S = {
    "title": ParagraphStyle("title", fontName="Times-Bold", fontSize=20,
                            leading=24, textColor=INK, spaceAfter=4),
    "h1": ParagraphStyle("h1", fontName="Times-Bold", fontSize=15,
                         leading=18, textColor=INK, spaceBefore=18,
                         spaceAfter=7),
    "h2": ParagraphStyle("h2", fontName="Times-Bold", fontSize=11.6,
                         leading=14, textColor=INK, spaceBefore=13,
                         spaceAfter=5),
    "h3": ParagraphStyle("h3", fontName="Times-BoldItalic", fontSize=10.2,
                         leading=12.6, textColor=INK, spaceBefore=10,
                         spaceAfter=4),
    "body": ParagraphStyle("body", fontName="Times-Roman", fontSize=9.7,
                           leading=13.0, textColor=INK, spaceAfter=6),
    "lead": ParagraphStyle("lead", fontName="Times-Roman", fontSize=10.6,
                           leading=14.4, textColor=SOFT, spaceAfter=10),
    "bullet": ParagraphStyle("bullet", fontName="Times-Roman", fontSize=9.7,
                             leading=12.8, textColor=INK, leftIndent=11,
                             bulletIndent=2, spaceAfter=3),
    "quote": ParagraphStyle("quote", fontName="Times-Italic", fontSize=9.4,
                            leading=12.6, textColor=SOFT, leftIndent=8,
                            rightIndent=8, spaceBefore=3, spaceAfter=3),
    "code": ParagraphStyle("code", fontName="Courier", fontSize=7.6,
                           leading=9.6, textColor=INK, leftIndent=5,
                           rightIndent=4),
    "cell": ParagraphStyle("cell", fontName="Times-Roman", fontSize=8.2,
                           leading=10.2, textColor=INK),
    "cellh": ParagraphStyle("cellh", fontName="Times-Bold", fontSize=8.0,
                            leading=10.0, textColor=INK),
    "caption": ParagraphStyle("caption", fontName="Times-Italic", fontSize=8.2,
                              leading=10.6, textColor=SOFT, spaceAfter=8),
}



# ------------------------------------------------- glyph safety
# The built-in Type-1 fonts are WinAnsi-encoded, so anything outside CP1252
# renders as a blank box. Superscripts become real <super> markup, Greek goes
# through the Symbol font, and the rest degrade to ASCII that reads the same.
SUPER = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3",
         "\u2074": "4", "\u2075": "5", "\u2076": "6", "\u2077": "7",
         "\u2078": "8", "\u2079": "9", "\u207b": "-", "\u207a": "+"}
# The Symbol font does not render reliably through this toolchain, so Greek
# letters are spelled out. In a primer that is arguably clearer anyway.
GREEK = {"\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma",
         "\u03b4": "delta", "\u03bd": "nu", "\u03bc": "mu",
         "\u03c1": "rho", "\u03c3": "sigma", "\u03b8": "theta",
         "\u0394": "delta"}
ASCII_MAP = {
    "\u21cc": "&lt;=&gt;", "\u2192": "-&gt;", "\u2190": "&lt;-",
    "\u226a": "&lt;&lt;", "\u226b": "&gt;&gt;",
    "\u2264": "&lt;=", "\u2265": "&gt;=", "\u2248": "~",
    "\u2212": "-", "\u2713": "yes", "\u2717": "no", "\u2022": "&bull;",
}
CODE_MAP = dict(ASCII_MAP)
CODE_MAP.update({"\u21cc": "<=>", "\u2192": "->", "\u2190": "<-",
                 "\u226a": "<<", "\u226b": ">>", "\u2264": "<=",
                 "\u2265": ">=", "\u2713": "yes", "\u2717": "no",
                 "\u2022": "*"})


def glyph_safe(t):
    """Prose: superscript runs -> <super>, Greek -> Symbol font, rest -> ASCII."""
    out, i = [], 0
    while i < len(t):
        c = t[i]
        if c in SUPER:
            run = ""
            while i < len(t) and t[i] in SUPER:
                run += SUPER[t[i]]; i += 1
            out.append(f"<super>{run}</super>")
            continue
        if c in GREEK:
            out.append(f"<i>{GREEK[c]}</i>")
        elif c in ASCII_MAP:
            out.append(ASCII_MAP[c])
        else:
            out.append(c)
        i += 1
    return "".join(out)


def code_safe(t):
    for k, v in CODE_MAP.items():
        t = t.replace(k, v)
    for k, v in SUPER.items():
        t = t.replace(k, "^" + v if v in "+-" else v)
    for k, v in GREEK.items():
        t = t.replace(k, v)
    return t


# ---------------------------------------------------------------- inline
def inline(t):
    """Markdown inline markup -> reportlab mini-HTML."""
    t = t.replace("&", "&amp;").replace("<sub>", "\x01").replace("</sub>", "\x02")
    t = re.sub(r"<[^>]+>", "", t)                      # drop other raw HTML
    t = t.replace("\x01", "<sub>").replace("\x02", "</sub>")
    t = t.replace("<", "&lt;").replace(">", "&gt;")
    t = t.replace("&lt;sub&gt;", "<font size=8>").replace("&lt;/sub&gt;", "</font>")
    # code spans first, so their contents are not re-processed
    spans = []
    def _code(m):
        spans.append(m.group(1))
        return f"\x00{len(spans)-1}\x00"
    t = re.sub(r"`([^`]+)`", _code, t)
    def _link(m):
        label, href = m.group(1), m.group(2)
        if href.startswith("#"):          # in-page anchor: no target in PDF
            return label
        return f'<link href="{href}" color="#2a78d6">{label}</link>'
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", _link, t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])\*([^*]+)\*(?![\w*])", r"<i>\1</i>", t)
    t = glyph_safe(t)
    for i, c in enumerate(spans):
        t = t.replace(f"\x00{i}\x00",
                      f'<font face="Courier" size=8.4>{code_safe(c)}</font>')
    return t


def code_lines(src):
    out = []
    for ln in src.split("\n"):
        ln = code_safe(ln)
        ln = (ln.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                .replace(" ", "&nbsp;"))
        out.append(Paragraph(ln or "&nbsp;", S["code"]))
    return out


def code_block(src):
    t = Table([[code_lines(src)]], colWidths=[170 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), CODEBG),
        ("BOX", (0, 0), (-1, -1), 0.4, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return t


def table_block(rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    data = [[Paragraph(inline(c), S["cellh"] if i == 0 else S["cell"])
             for c in r] for i, r in enumerate(rows)]
    avail = 170 * mm
    # width by longest cell, first column a little roomier
    raw = [max(len(rows[i][j]) for i in range(len(rows))) for j in range(ncol)]
    raw = [max(r, 4) ** 0.72 for r in raw]
    tot = sum(raw)
    widths = [avail * r / tot for r in raw]
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, INK),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
        ("LINEBELOW", (0, -1), (-1, -1), 0.7, INK),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]))
    return t


def quote_block(parts):
    body = [Paragraph(inline(x), S["quote"]) for x in parts]
    t = Table([[body]], colWidths=[170 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), QUOTEBG),
        ("LINEBEFORE", (0, 0), (0, -1), 1.6, ACCENT),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return t


def hr():
    t = Table([[""]], colWidths=[170 * mm], rowHeights=[1])
    t.setStyle(TableStyle([("LINEABOVE", (0, 0), (-1, 0), 0.5, RULE)]))
    return t


# ---------------------------------------------------------------- parse
def build(md):
    lines = md.split("\n")
    flow, i, first_h1 = [], 0, True
    while i < len(lines):
        ln = lines[i]

        if ln.startswith("```"):
            j = i + 1
            buf = []
            while j < len(lines) and not lines[j].startswith("```"):
                buf.append(lines[j]); j += 1
            flow += [Spacer(1, 3), code_block("\n".join(buf)), Spacer(1, 7)]
            i = j + 1
            continue

        if ln.strip().startswith("<details"):
            i += 1
            continue
        if ln.strip().startswith("</details"):
            flow.append(Spacer(1, 4)); i += 1
            continue
        m = re.match(r"\s*<summary>(.*?)</summary>", ln)
        if m:
            flow.append(Paragraph(inline(m.group(1)), S["h3"]))
            i += 1
            continue

        if re.match(r"^\s*\|", ln) and i + 1 < len(lines) and \
                re.match(r"^\s*\|[\s:|-]+\|?\s*$", lines[i + 1]):
            rows, j = [], i
            while j < len(lines) and re.match(r"^\s*\|", lines[j]):
                if not re.match(r"^\s*\|[\s:|-]+\|?\s*$", lines[j]):
                    cells = [c.strip() for c in lines[j].strip().strip("|").split("|")]
                    rows.append(cells)
                j += 1
            flow += [Spacer(1, 3), table_block(rows), Spacer(1, 8)]
            i = j
            continue

        if re.match(r"^> |^>$", ln):
            buf, j = [], i
            while j < len(lines) and re.match(r"^> |^>$", lines[j]):
                buf.append(lines[j].lstrip("> ").rstrip()); j += 1
            para = " ".join(x for x in buf if x).split("  ")
            flow += [Spacer(1, 3), quote_block([" ".join(buf).strip()]),
                     Spacer(1, 7)]
            i = j
            continue

        if re.match(r"^(-{3,}|\*{3,})\s*$", ln):
            flow += [Spacer(1, 8), hr(), Spacer(1, 10)]
            i += 1
            continue

        m = re.match(r"^(#{1,4})\s+(.*)$", ln)
        if m:
            lvl, txt = len(m.group(1)), m.group(2)
            if lvl == 1 and first_h1:
                flow.append(Paragraph(inline(txt), S["title"]))
                first_h1 = False
            else:
                flow.append(Paragraph(inline(txt),
                                      S["h1" if lvl <= 2 else
                                        ("h2" if lvl == 3 else "h3")]))
            i += 1
            continue

        m = re.match(r"^(\s*)([-*])\s+(.*)$", ln)
        if m:
            indent = len(m.group(1))
            txt = m.group(3)
            j = i + 1
            while j < len(lines) and lines[j].startswith("  ") and \
                    lines[j].strip() and not re.match(r"^\s*[-*]\s", lines[j]):
                txt += " " + lines[j].strip(); j += 1
            st = ParagraphStyle("b2", parent=S["bullet"],
                                leftIndent=11 + indent * 5,
                                bulletIndent=2 + indent * 5)
            flow.append(Paragraph(inline(txt), st, bulletText="•"))
            i = j
            continue

        m = re.match(r"^(\s*)(\d+)\.\s+(.*)$", ln)
        if m:
            txt = m.group(3)
            j = i + 1
            while j < len(lines) and lines[j].startswith("   ") and \
                    lines[j].strip():
                txt += " " + lines[j].strip(); j += 1
            flow.append(Paragraph(inline(txt), S["bullet"],
                                  bulletText=m.group(2) + "."))
            i = j
            continue

        if not ln.strip():
            i += 1
            continue

        buf, j = [ln], i + 1
        while j < len(lines) and lines[j].strip() and \
                not re.match(r"^(#{1,4}\s|\s*\||```|> |>$|-{3,}|\s*[-*]\s"
                             r"|\s*\d+\.\s|<)", lines[j]):
            buf.append(lines[j]); j += 1
        para = " ".join(x.strip() for x in buf)
        style = S["lead"] if (len(flow) < 4) else S["body"]
        flow.append(Paragraph(inline(para), style))
        i = j
    return flow


def render(src_path, out_path):
    md = open(src_path, encoding="utf-8").read()
    title = next((l[2:].strip() for l in md.split("\n") if l.startswith("# ")),
                 os.path.basename(src_path))

    def _page(canvas, doc):
        canvas.saveState()
        canvas.setFont("Times-Roman", 7.6)
        canvas.setFillColor(SOFT)
        if doc.page > 1:
            canvas.drawString(20 * mm, A4[1] - 12 * mm, title)
        canvas.drawRightString(A4[0] - 20 * mm, 12 * mm, str(doc.page))
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.3)
        if doc.page > 1:
            canvas.line(20 * mm, A4[1] - 14.5 * mm, A4[0] - 20 * mm,
                        A4[1] - 14.5 * mm)
        canvas.restoreState()

    doc = BaseDocTemplate(out_path, pagesize=A4,
                          leftMargin=20 * mm, rightMargin=20 * mm,
                          topMargin=20 * mm, bottomMargin=18 * mm,
                          title=title, author="PFAS TK review")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height,
                  id="f", showBoundary=0)
    doc.addPageTemplates([PageTemplate(id="body", frames=[frame],
                                       onPage=_page)])
    doc.build(build(md))
    size = os.path.getsize(out_path)
    print(f"wrote {out_path}  ({size/1024:.0f} KB)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = sys.argv[1]
    if not os.path.isabs(src):
        src = os.path.join(HERE, src)
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(src)[0] + ".pdf"
    render(src, out)
