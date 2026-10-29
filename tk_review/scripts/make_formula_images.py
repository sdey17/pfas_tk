#!/usr/bin/env python3
"""Typeset the document's twelve display formulas as images.

WHY THIS EXISTS. The formulas live in the shared document as ```latex code
blocks, which is how that document type wants a display formula written, and
the editor typesets them. The PDF and HTML exports do not: checked against
the raw export bytes, they emit every one as
`<pre><code class="language-latex">...</code></pre>` with no KaTeX, no
MathML and no `$$` anywhere in the file. The PDF is the copy that gets
shared, and a derivation is the one part of that document a reader has to be
able to follow, so the formulas are rendered here and embedded as images.

WHAT IT COSTS, SO IT IS NOT DISCOVERED LATER. An image is not text: these
twelve formulas stop being selectable, searchable and editable in the
document, and a screen reader gets only the alt text. The alt text is
therefore the LaTeX source itself, verbatim, which is the most useful thing
it could carry and lets the blocks be reconstructed if the exporter ever
learns to typeset. They are drawn on a transparent background in near-black,
which reads correctly on the white page of the PDF and in the editor's light
theme, and poorly in a dark theme -- the one tradeoff that has no answer
inside an image.

Rendering is matplotlib's own mathtext, not a TeX installation, so the
subset is smaller: \\text is unavailable and \\mathrm takes its place, and
the negative thin space \\! is dropped. Both are spacing and typeface
details that change nothing about what the formulas say.

Run:  python3 tk_review/scripts/make_formula_images.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "formulas")

INK = "#14161a"
DPI = 220
SIZE = 17

# (stem, the TeX as the document stores it, the mathtext to draw).
# The two differ only where mathtext lacks a construct; keeping the source
# column means the alt text can quote what the document actually had.
FORMULAS = [
    ("01_kd_definition",
     r"K_d = \frac{[P][C_w]}{[PL]}",
     r"K_d = \frac{[P][C_w]}{[PL]}"),
    ("02_occupancy",
     r"\theta = \frac{C_w}{K_d + C_w}",
     r"\theta = \frac{C_w}{K_d + C_w}"),
    ("03_bound",
     r"C_{bound} = P \cdot n \cdot \theta = \frac{P \, n \, C_w}{K_d + C_w}",
     r"C_{bound} = P \cdot n \cdot \theta = \frac{P \, n \, C_w}{K_d + C_w}"),
    ("04_divide_by_mass",
     r"K = \frac{C_{bound}}{C_w} \cdot \frac{1000}{P \, M}"
     r" = \frac{1000 \, n}{M \,(K_d + C_w)}",
     r"K = \frac{C_{bound}}{C_w} \cdot \frac{1000}{P \, M}"
     r" = \frac{1000 \, n}{M \,(K_d + C_w)}"),
    ("05_trace_limit",
     r"K^0 = \frac{1000 \, n}{M \, K_d}",
     r"K^0 = \frac{1000 \, n}{M \, K_d}"),
    ("06_practical_form",
     r"\log K = \log_{10}\!\left(\frac{10^9}{M}\right)"
     r" + \log n - \log K_d(\mu M)",
     r"\log K = \log_{10}\left(\frac{10^9}{M}\right)"
     r" + \log n - \log K_d(\mu M)"),
    ("07_worked_example",
     r"\log K = 4.8475 + \log_{10}(1) - \log_{10}(0.099) = 5.85",
     r"\log K = 4.8475 + \log_{10}(1) - \log_{10}(0.099) = 5.85"),
    ("08_several_classes",
     r"K^0 = \frac{10^9}{M} \sum_i \frac{n_i}{K_{d,i}}",
     r"K^0 = \frac{10^9}{M} \sum_i \frac{n_i}{K_{d,i}}"),
    ("09_saturation",
     r"K(C_w) = \frac{K^0}{1 + C_w/K_d}",
     r"K(C_w) = \frac{K^0}{1 + C_w/K_d}"),
    ("10_reverse",
     r"K_d = 10^{\,\log_{10}(10^9/M) \,+\, \log n \,-\, \log K}",
     r"K_d = 10^{\,\log_{10}(10^9/M) \,+\, \log n \,-\, \log K}"),
    ("11_trace_limit_uM",
     r"K^0 = \frac{10^9 \, n}{M \, K_d}",
     r"K^0 = \frac{10^9 \, n}{M \, K_d}"),
    ("12_bound_per_litre",
     r"\text{bound per litre of phase}"
     r" = C_{protein} \cdot n \cdot \frac{C_w}{K_d}",
     r"\mathrm{bound\ per\ litre\ of\ phase}"
     r" = C_{protein} \cdot n \cdot \frac{C_w}{K_d}"),
]


# Every image is pasted onto a canvas of this width before it is saved.
# Without it each formula would be cropped to its own ink and then scaled to
# the column width by the renderer, so the shortest formula would come out
# the largest -- K^0 = 1000n/(M Kd) in letters twice the size of the line
# above it. A fixed canvas makes the scale factor the same for all twelve,
# which is what puts them all at roughly body-text size on the page.
CANVAS_W = 1600
INDENT = 60
PAD_Y = 14


def render(stem, mathtext):
    """One formula, cropped to its ink, then placed on the shared canvas."""
    from PIL import Image

    fig = plt.figure(figsize=(0.1, 0.1))
    fig.text(0, 0, f"${mathtext}$", fontsize=SIZE, color=INK,
             ha="left", va="bottom")
    path = os.path.join(OUT, stem + ".png")
    fig.savefig(path, dpi=DPI, transparent=True, bbox_inches="tight",
                pad_inches=0.06)
    plt.close(fig)

    ink = Image.open(path).convert("RGBA")
    if ink.width > CANVAS_W - 2 * INDENT:          # the long ones only
        scale = (CANVAS_W - 2 * INDENT) / ink.width
        ink = ink.resize((int(ink.width * scale), int(ink.height * scale)),
                         Image.LANCZOS)
    canvas = Image.new("RGBA", (CANVAS_W, ink.height + 2 * PAD_Y), (0, 0, 0, 0))
    canvas.paste(ink, (INDENT, PAD_Y), ink)
    canvas.save(path)
    return path, ink.width


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"{len(FORMULAS)} formulas -> {os.path.abspath(OUT)}")
    print(f"  canvas {CANVAS_W}px wide, so all twelve scale alike on the page")
    for stem, source, mathtext in FORMULAS:
        path, inked = render(stem, mathtext)
        note = "" if source == mathtext else "   (mathtext substitution)"
        fit = "" if inked < CANVAS_W - 2 * INDENT else "   (scaled to fit)"
        print(f"  {os.path.getsize(path):7d} B  ink {inked:5d}px  "
              f"{stem}{note}{fit}")


if __name__ == "__main__":
    main()
