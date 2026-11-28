"""Mark the TRBV17 comparison in the original two-panel cutoff plot."""

import argparse
from pathlib import Path
from PIL import Image, ImageDraw


def mark_trbv17(source, output):
    image = Image.open(source).convert("RGB")
    if image.size != (1050, 482):
        raise ValueError("Use the original 1050 x 482 Figure 13 export.")
    scale = 3
    image = image.resize((image.width * scale, image.height * scale), Image.Resampling.LANCZOS)
    draw = ImageDraw.Draw(image)
    for x, y, rx, ry in [(137, 99, 13, 19), (658, 91, 11, 11), (1008, 107, 11, 11)]:
        draw.ellipse(((x-rx)*scale, (y-ry)*scale, (x+rx)*scale, (y+ry)*scale),
                     outline="#b92835", width=4)
    image.save(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    mark_trbv17(args.source, args.output)
