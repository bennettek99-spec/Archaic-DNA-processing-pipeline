#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_pdf.py — render a Markdown manuscript to a clean, academic-looking PDF
using reportlab Platypus. Supports headings, paragraphs, **bold**/*italic*/
`code`, markdown tables, bullet lists, fenced code blocks, horizontal rules
(---), markdown/bare hyperlinks, and inline ![caption](figure.png) images
(resolved relative to the source markdown file's own directory, so each
report's PAPER*.md can sit next to its own figures).

    python make_pdf.py <source.md> <out.pdf>
    python make_pdf.py                              # legacy default: PAPER.md -> reports/Etruscan_paper.pdf

The legacy default (no args) reproduces the original Etruscan-paper behaviour,
which has no inline ![]() image syntax of its own and instead relies on a
hardcoded figure list inserted before the "Discussion" heading.

Uses matplotlib's bundled DejaVuSans TTF so Greek/maths glyphs (alpha, ->, <=)
render correctly.
"""
import os, sys, re, html
import matplotlib
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, KeepTogether, Preformatted,
                                HRFlowable)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

HERE = os.path.dirname(os.path.abspath(__file__))

# Etruscan-paper backward-compat figure list — only used as a fallback when the
# source document has no inline ![](...) images of its own (see main()).
LEGACY_FIGDIR = os.path.join(HERE, "results", "figures")
LEGACY_FIGURES = [
    ("fig_p1_group_neanderthal.png", "Figure 1. Group-level Neanderthal ancestry (mean-genome f4-ratio) across ancient and modern populations, with block-jackknife standard errors."),
    ("fig_p2_temporal.png", "Figure 2. Group-level Neanderthal ancestry across the Italian time transect; flat within error."),
    ("fig_p3_mds.png", "Figure 3. Multidimensional scaling of genome-wide allele-frequency distances (cohorts with n>=30). Etruscans cluster with Romans and other Italians."),
    ("etr_loci_heatmap.png", "Figure 4. Mean Neanderthal-allele frequency at adaptive-introgression loci across the Italian time transect."),
    ("fig_p5_qpadm.png", "Figure 5. qpAdm ancestry proportions (Anatolian-farmer + Steppe + WHG): Etruscans match Latins and Imperial Romans; Steppe rises from the Bronze Age to the Iron Age."),
    ("fig_concordance.png", "Figure 6. Concordance with ADMIXTOOLS 2: f-statistics correlate at r>=0.99 and qpAdm ancestry proportions agree within ~4 percentage points."),
    ("fig_admixture.png", "Figure 7. Sparse-NMF ancestry components (snmf/ADMIXTURE-style, K=4): source populations carry distinct components; the Italian targets (Bronze Age, Etruscan, Latin, Imperial Roman) are visually identical mixtures, reinforcing continuity."),
]

IMG_RE = re.compile(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$")
META_RE = re.compile(r"^\*\*[^*:]{2,30}:\*\*")
PAGE_W = A4[0] - 4 * cm   # usable width given 2cm left+right margins


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
        "meta2": ParagraphStyle("meta2", parent=body, fontSize=8.5, textColor=colors.grey,
                                fontName="DejaVu", alignment=0, spaceAfter=2),
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
        "bullet": ParagraphStyle("bullet", parent=body, leftIndent=14, firstLineIndent=-14,
                                 alignment=0, spaceAfter=4),
        "code": ParagraphStyle("code", parent=ss["Normal"], fontName="Courier",
                               fontSize=7.5, leading=9.5, textColor=colors.HexColor("#2b2b2b"),
                               backColor=colors.HexColor("#f4f2ee"),
                               borderColor=colors.HexColor("#ddd6c8"), borderWidth=0.5,
                               borderPadding=6, spaceBefore=4, spaceAfter=8),
        "body": body,
    }


def inline(text):
    text = html.escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`(.+?)`", r'<font face="Courier">\1</font>', text)
    # markdown links [text](url)
    text = re.sub(r'\[([^\]]+)\]\((https?://[^)\s]+)\)',
                  r'<link href="\2" color="#1a5276"><u>\1</u></link>', text)
    # bare URLs not already inside an href="..." attribute
    text = re.sub(r'(?<!")(?<!>)(https?://[^\s<>"]+)',
                  r'<link href="\1" color="#1a5276"><u>\1</u></link>', text)
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


def make_image(relpath, srcdir):
    """Resolve an image path relative to the source markdown's directory,
    falling back to the legacy Etruscan figures/ dir for old documents."""
    for base in (srcdir, LEGACY_FIGDIR):
        fp = os.path.normpath(os.path.join(base, relpath))
        if os.path.exists(fp):
            iw, ih = ImageReader(fp).getSize()
            w = min(PAGE_W, 15.5 * cm)
            return Image(fp, width=w, height=w * ih / iw)
    return None


def legacy_figures_flow(S):
    flow = []
    for fn, cap in LEGACY_FIGURES:
        fp = os.path.join(LEGACY_FIGDIR, fn)
        if not os.path.exists(fp):
            continue
        iw, ih = ImageReader(fp).getSize()
        w = 15.5 * cm
        flow.append(KeepTogether([Image(fp, width=w, height=w * ih / iw),
                                  Spacer(1, 3), Paragraph(inline(cap), S["cap"])]))
    return flow


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "PAPER.md")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "reports", "Etruscan_paper.pdf")
    register_fonts(); S = styles()
    srcdir = os.path.dirname(os.path.abspath(src))
    lines = open(src, encoding="utf-8").read().splitlines()
    has_inline_images = any(IMG_RE.match(l.strip()) for l in lines)

    doc_title = "Manuscript"
    story, i, figs_done = [], 0, False
    n = len(lines)
    while i < n:
        raw = lines[i]
        s = raw.strip()

        if not s:
            i += 1; continue

        m_img = IMG_RE.match(s)
        if m_img:
            img = make_image(m_img.group(2), srcdir)
            if img:
                story.append(KeepTogether([img, Spacer(1, 10)]))
            i += 1; continue

        if s == "---":
            story.append(Spacer(1, 4))
            story.append(HRFlowable(width="100%", thickness=0.6,
                                    color=colors.HexColor("#c9c2b4")))
            story.append(Spacer(1, 8))
            i += 1; continue

        if s.startswith("```"):
            i += 1
            code_lines = []
            while i < n and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i]); i += 1
            i += 1  # skip closing fence
            story.append(Preformatted("\n".join(code_lines), S["code"]))
            continue

        if raw.startswith("# "):
            doc_title = raw[2:].strip()
            story.append(Paragraph(inline(raw[2:]), S["title"]))
            i += 1; continue

        if raw.startswith("## "):
            if (not has_inline_images) and raw[3:].strip().startswith("Discussion") and not figs_done:
                story += legacy_figures_flow(S); figs_done = True
            story.append(Paragraph(inline(raw[3:]), S["h1"]))
            i += 1; continue

        if raw.startswith("### "):
            story.append(Paragraph(inline(raw[4:]), S["h2"]))
            i += 1; continue

        if raw.startswith("|"):
            tbl = []
            while i < n and lines[i].startswith("|"):
                if not re.match(r"^\|[\s:-]+\|", lines[i]):
                    tbl.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            story.append(build_table(tbl, S)); story.append(Spacer(1, 8))
            continue

        if s.startswith("- "):
            story.append(Paragraph("&bull;&nbsp;" + inline(s[2:]), S["bullet"]))
            i += 1; continue

        if META_RE.match(s):
            story.append(Paragraph(inline(s), S["meta2"]))
            i += 1; continue

        if s.startswith("*") and not s.startswith("**") and s.endswith("*") and len(s) > 2:
            story.append(Paragraph(inline(s), S["meta"]))
            i += 1; continue

        # default: body paragraph, pairing with an immediately-following image
        # (a bold "**Figure N.** caption" line followed by "![Figure N](file)")
        sty = S["abstract"] if (i > 0 and lines[i - 1].strip() == "## Abstract") else S["body"]
        para = Paragraph(inline(s), sty)
        nxt = lines[i + 1].strip() if i + 1 < n else ""
        m_next = IMG_RE.match(nxt)
        if m_next:
            img = make_image(m_next.group(2), srcdir)
            if img:
                story.append(KeepTogether([para, Spacer(1, 4), img, Spacer(1, 10)]))
                i += 2; continue
        story.append(para)
        i += 1

    if not has_inline_images and not figs_done:
        story += legacy_figures_flow(S)

    doc = SimpleDocTemplate(out, pagesize=A4, topMargin=1.6 * cm, bottomMargin=1.6 * cm,
                            leftMargin=2.0 * cm, rightMargin=2.0 * cm, title=doc_title)
    doc.build(story)
    print(f"Wrote {out}  ({os.path.getsize(out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
