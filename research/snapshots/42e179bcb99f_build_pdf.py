from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image,
    PageBreak,
)


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "toi1258_statistical_validation.md"
OUT = ROOT / "paper" / "toi1258_statistical_validation.pdf"
FIG = ROOT / "evidence" / "figures"


def clean_inline(text: str) -> str:
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"`([^`]+)`", r"<font name='Courier'>\1</font>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    return text


def bullet_table(lines):
    rows = []
    for line in lines:
        if line.startswith("|") and line.endswith("|"):
            cells = [clean_inline(c.strip()) for c in line.strip("|").split("|")]
            if not all(set(c) <= {"-", ":"} for c in cells):
                rows.append(cells)
    return rows


def add_image(story, name, caption):
    path = FIG / name
    if not path.exists():
        return
    story.append(Spacer(1, 0.12 * inch))
    img = Image(str(path))
    max_w = 6.4 * inch
    max_h = 4.2 * inch
    scale = min(max_w / img.imageWidth, max_h / img.imageHeight, 1.0)
    img.drawWidth = img.imageWidth * scale
    img.drawHeight = img.imageHeight * scale
    story.append(img)
    story.append(Paragraph(clean_inline(caption), STYLES["Caption"]))
    story.append(Spacer(1, 0.12 * inch))


STYLES = getSampleStyleSheet()
STYLES.add(
    ParagraphStyle(
        name="TitleCenter",
        parent=STYLES["Title"],
        alignment=TA_CENTER,
        fontSize=16,
        leading=20,
        spaceAfter=14,
    )
)
STYLES.add(
    ParagraphStyle(
        name="AuthorCenter",
        parent=STYLES["Normal"],
        alignment=TA_CENTER,
        fontSize=10,
        leading=13,
        spaceAfter=14,
    )
)
STYLES.add(
    ParagraphStyle(
        name="Body",
        parent=STYLES["BodyText"],
        alignment=TA_LEFT,
        fontSize=9.2,
        leading=12,
        spaceAfter=6,
    )
)
STYLES.add(
    ParagraphStyle(
        name="BulletBody",
        parent=STYLES["BodyText"],
        leftIndent=14,
        firstLineIndent=-8,
        fontSize=9.2,
        leading=12,
        spaceAfter=4,
    )
)
STYLES.add(
    ParagraphStyle(
        name="Caption",
        parent=STYLES["Italic"],
        fontSize=8,
        leading=10,
        textColor=colors.HexColor("#333333"),
        spaceAfter=8,
    )
)


def build():
    text = PAPER.read_text(encoding="utf-8").splitlines()
    doc = SimpleDocTemplate(
        str(OUT),
        pagesize=letter,
        rightMargin=0.72 * inch,
        leftMargin=0.72 * inch,
        topMargin=0.72 * inch,
        bottomMargin=0.72 * inch,
        title="TOI-1258.01 Independent Statistical Validation",
        author="Nikolai Stoloff",
    )

    story = []
    i = 0
    inserted_figures = False
    while i < len(text):
        line = text[i].rstrip()
        if not line:
            i += 1
            continue
        if line.startswith("# "):
            story.append(Paragraph(clean_inline(line[2:]), STYLES["TitleCenter"]))
        elif line.startswith("Nikolai Stoloff"):
            author_lines = []
            while i < len(text) and text[i].strip():
                author_lines.append(clean_inline(text[i].strip()))
                i += 1
            story.append(Paragraph("<br/>".join(author_lines), STYLES["AuthorCenter"]))
            continue
        elif line.startswith("## "):
            story.append(Spacer(1, 0.08 * inch))
            story.append(Paragraph(clean_inline(line[3:]), STYLES["Heading2"]))
        elif line.startswith("|"):
            tbl_lines = []
            while i < len(text) and text[i].startswith("|"):
                tbl_lines.append(text[i])
                i += 1
            rows = bullet_table(tbl_lines)
            if rows:
                table = Table(rows, repeatRows=1)
                table.setStyle(
                    TableStyle(
                        [
                            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")),
                            ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                            ("FONTSIZE", (0, 0), (-1, -1), 7.4),
                            ("LEADING", (0, 0), (-1, -1), 8.5),
                            ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ]
                    )
                )
                story.append(table)
                story.append(Spacer(1, 0.1 * inch))
            continue
        elif line.startswith("- "):
            story.append(Paragraph("• " + clean_inline(line[2:]), STYLES["BulletBody"]))
        else:
            story.append(Paragraph(clean_inline(line), STYLES["Body"]))

        if not inserted_figures and line == "## 4. Centroid and Source-Location Checks":
            add_image(story, "toi1258_allsector_folded_catalog_ephemeris.png", "Figure 1. All-sector TESS light curve folded at the catalog ephemeris.")
            add_image(story, "toi1258_allsector_bls_periodogram.png", "Figure 2. BLS periodogram showing recovery of the TOI-1258.01 signal.")
            inserted_figures = True
        if line == "## 5. High-Resolution Imaging Constraint":
            add_image(story, "toi1258_difference_image_offsets.png", "Figure 3. Custom difference-image offsets used as a conservative centroid diagnostic.")
            add_image(story, "toi1258_gemini_alopeke_contrast_curves.png", "Figure 4. Gemini/Alopeke contrast curves ingested into the validation workflow.")

        i += 1

    doc.build(story)
    print(OUT)


if __name__ == "__main__":
    build()
