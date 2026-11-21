"""Render figures at the physical width used in the IJMS manuscript."""

import re

from matplotlib import patheffects, rc_context
from matplotlib.colors import to_rgb
from matplotlib.font_manager import FontProperties, findfont, fontManager
from matplotlib.text import Text
from matplotlib.transforms import Bbox


MDPI_WIDTH_IN = 180 / 25.4
MDPI_FONTS = ["TeX Gyre Pagella", "P052", "STIXGeneral", "DejaVu Serif"]
PANEL_LABEL = re.compile(r"^([A-D])\.\s+")


def save_publication_figure(figure, png):
    """Write matching PNG and vector PDF with consistent final-size type."""
    width, height = figure.get_size_inches()
    figure.set_size_inches(MDPI_WIDTH_IN, height * MDPI_WIDTH_IN / width)
    figure.canvas.draw()
    available = {font.name for font in fontManager.ttflist}
    family = next((name for name in MDPI_FONTS if name in available), "serif")
    sizes = {}
    for ax in figure.axes:
        for item in (ax.xaxis.label, ax.yaxis.label):
            sizes[item] = 9
        for item in (ax.title, ax._left_title, ax._right_title):
            sizes[item] = 9.5
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            sizes[item] = 8
        for spine in ax.spines.values():
            if spine.get_visible():
                spine.set_linewidth(0.65)
                spine.set_edgecolor("black")
    for item in figure.findobj(Text):
        label = item.get_text().replace("log₂", r"log$_2$")
        label = PANEL_LABEL.sub(lambda match: f"({match[1].lower()})  ", label)
        item.set_text("\n".join(line.removesuffix(".") for line in label.split("\n")))
        if to_rgb(item.get_color()) == (1.0, 1.0, 1.0):
            # Heatmap values need a light boundary on the darkest cells.
            item.set_path_effects([patheffects.withStroke(linewidth=1.5, foreground="white")])
        item.set_color("black")
        item.set_math_fontfamily("stix")
        item.set_fontfamily(family)
        item.set_fontweight("normal")
        item.set_fontsize(sizes.get(item, 8.25))
    for artist in (figure._suptitle, figure._supxlabel, figure._supylabel):
        if artist is not None:
            artist.set_fontsize(10 if artist is figure._suptitle else 9)

    figure.canvas.draw()
    bounds = figure.get_tightbbox(figure.canvas.get_renderer()).padded(0.055)
    width, height = figure.get_size_inches()
    box = Bbox.union([Bbox.from_bounds(0, 0, width, height), bounds])
    png.parent.mkdir(parents=True, exist_ok=True)
    # Type 42 needs TrueType outlines. The installed Palatino substitute is CFF.
    font_path = findfont(FontProperties(family=family))
    pdf_fonttype = 3 if font_path.lower().endswith(".otf") else 42
    with rc_context({"pdf.fonttype": pdf_fonttype}):
        figure.savefig(png, dpi=450, bbox_inches=box, facecolor="white")
        figure.savefig(png.with_suffix(".pdf"), dpi=450, bbox_inches=box, facecolor="white")
