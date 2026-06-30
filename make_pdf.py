#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_pdf.py — render a Markdown manuscript (default PAPER.md) to a clean,
academic-looking PDF using reportlab Platypus. Headings, paragraphs, **bold**,
*italic*, `code`, Markdown tables, and the manuscript figures are supported.

    python make_pdf.py PAPER.md reports/Etruscan_paper.pdf

Uses matplotlib's bundled DejaVuSans TTF so Greek/maths glyphs (alpha, ->, <=)
render correctly.
"""
import os, sys, re, html
import matplotlib
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, KeepTogether)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "results", "figures")
# Figures to embed before the Discussion (file, caption)
FIGURES = [
    ("fig_p1_group_neanderthal.png", "Figure 1. Group-level Neanderthal ancestry (mean-genome f4-ratio) across ancient and modern populations, with block-jackknife standard errors."),
    ("fig_p2_temporal.png", "Figure 2. Group-level Neanderthal ancestry across the Italian time transect; flat within error."),
    ("fig_p3_mds.png", "Figure 3. Multidimensional scaling of genome-wide allele-frequency distances (cohorts with n>=30). Etruscans cluster with Romans and other Italians."),
    ("etr_loci_heatmap.png", "Figure 4. Mean Neanderthal-allele frequency at adaptive-introgression loci across the Italian time transect."),
    ("fig_p5_qpadm.png", "Figure 5. qpAdm ancestry proportions (Anatolian-farmer + Steppe + WHG): Etruscans match Latins and Imperial Romans; Steppe rises from the Bronze Age to the Iron Age."),
]


def register_fonts():
    base = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
    pdfmetrics.registerFont(TTFont("DejaVu", os.path.join(base, "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVu-Bold", os.path.join(base, "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVu-Italic", os.path.join(base, "DejaVuSans-Oblique.ttf")))
    pdfmetrics.registerFontFamily("DejaVu", normal="DejaVu", bold="DejaVu-Bold",
                                  italic="DejaVu-Italic", boldItalic="DejaVu-Bold")


def styles():
    ss = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=ss["Normal"], fontName="DejaVu",
                          fontSize=9.5, leading=14, spaceAfter=6, alignment=4)
    return {
        "title": ParagraphStyle("title", parent=ss["Title"], fontName="DejaVu-Bold",
                                fontSize=16, leading=20, spaceAfter=4),
        "meta": ParagraphStyle("meta", parent=body, fontSize=8, textColor=colors.grey,
                               alignment=1, spaceAfter=10),
        "h1": ParagraphStyle("h1", parent=ss["Heading1"], fontName="DejaVu-Bold",
                             fontSize=12.5, leading=15, spaceBefore=12, spaceAfter=5,
                             textColor=colors.HexColor("#8a3b2e")),
        "h2": ParagraphStyle("h2", parent=ss["Heading2"], fontName="DejaVu-Bold",
                             fontSize=10.5, leading=13, spaceBefore=8, spaceAfter=3),
        "abstract": ParagraphStyle("abstract", parent=body, fontSize=9, leading=13,
                                   backColor=colors.HexColor("#f3efe7"),
                                   borderColor=colors.HexColor("#e2ddd4"), borderWidth=0.5,
                                   borderPadding=8, spaceAfter=10),
        "cap": ParagraphStyle("cap", parent=body, fontSize=8, textColor=colors.grey,
                              alignment=0, spaceAfter=12),
        "cell": ParagraphStyle("cell", parent=body, fontSize=8, leading=10, spaceAfter=0),
        "cellh": ParagraphStyle("cellh", parent=body, fontSize=8, leading=10,
                                fontName="DejaVu-Bold", textColor=colors.white, spaceAfter=0),
        "body": body,
    }


def inline(text):
    text = html.escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`(.+?)`", r'<font face="Courier">\1</font>', text)
    return text


def build_table(rows, S):
    header = [Paragraph(inline(c), S["cellh"]) for c in rows[0]]
    body = [[Paragraph(inline(c), S["cell"]) for c in r] for r in rows[1:]]
    t = Table([header] + body, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#8a3b2e")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7f4ee")]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#d9d2c6")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


def figures_flow(S):
    flow = []
    for fn, cap in FIGURES:
        fp = os.path.join(FIGDIR, fn)
        if not os.path.exists(fp):
            continue
        from reportlab.lib.utils import ImageReader
        iw, ih = ImageReader(fp).getSize()
        w = 15.5 * cm
        flow.append(KeepTogether([Image(fp, width=w, height=w * ih / iw),
                                  Spacer(1, 3), Paragraph(inline(cap), S["cap"])]))
    return flow


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "PAPER.md")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "reports", "Etruscan_paper.pdf")
    register_fonts(); S = styles()
    lines = open(src, encoding="utf-8").read().splitlines()

    story, i, figs_done = [], 0, False
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("# "):
            story.append(Paragraph(inline(ln[2:]), S["title"]))
        elif ln.startswith("## "):
            if ln[3:].strip().startswith("Discussion") and not figs_done:
                story += figures_flow(S); figs_done = True
            story.append(Paragraph(inline(ln[3:]), S["h1"]))
        elif ln.startswith("### "):
            story.append(Paragraph(inline(ln[4:]), S["h2"]))
        elif ln.startswith("|"):
            tbl = []
            while i < len(lines) and lines[i].startswith("|"):
                if not re.match(r"^\|[\s:-]+\|", lines[i]):
                    tbl.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            story.append(build_table(tbl, S)); story.append(Spacer(1, 8)); continue
        elif ln.startswith("*") and ln.endswith("*") and len(ln) > 2:
            story.append(Paragraph(inline(ln), S["meta"]))
        elif ln.strip():
            sty = S["abstract"] if (i > 0 and lines[i - 1].strip() == "## Abstract") else S["body"]
            story.append(Paragraph(inline(ln), sty))
        i += 1
    if not figs_done:
        story += figures_flow(S)

    doc = SimpleDocTemplate(out, pagesize=A4, topMargin=1.6 * cm, bottomMargin=1.6 * cm,
                            leftMargin=2.0 * cm, rightMargin=2.0 * cm,
                            title="Neanderthal ancestry in the Etruscans")
    doc.build(story)
    print(f"Wrote {out}  ({os.path.getsize(out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
