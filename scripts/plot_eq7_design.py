#!/usr/bin/env python3
"""Render the paired-score methods diagrams as vector artwork."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def diagram(path, height, boxes, links):
    fig, ax = plt.subplots(figsize=(7.1, height))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    for source, target in links:
        x1, y1, *_ = boxes[source]
        x2, y2, *_ = boxes[target]
        ax.annotate("", (x2, y2), (x1, y1),
                    arrowprops={"arrowstyle": "->", "color": "#383838", "shrinkA": 15, "shrinkB": 17})
    for x, y, text, color in boxes:
        ax.text(x, y, text, ha="center", va="center", fontsize=9,
                bbox={"boxstyle": "round,pad=0.5", "fc": color, "ec": "#383838", "lw": 0.7})
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    gray, violet, peach = "#f7f7f7", "#f1e5f4", "#fce9dd"
    diagram(output / "parallel_scores.pdf", 1.30, [
        (0.25, .79, "Shared structural shortlist\nLog residual aaV count", gray),
        (0.75, .79, "Raw aaV counts\nedgeR and Fisher", gray),
        (0.25, .38, "Binary evidence ridge", violet),
        (0.75, .38, "Continuous evidence ridge", peach),
        (0.50, .10, "Common endpoint: joint recovery contrast", gray),
    ], [(0, 2), (0, 3), (1, 2), (1, 3), (2, 4), (3, 4)])
    diagram(output / "paired_holdout.pdf", 2.40, [
        (.50, .86, "Outer test: leave out one mouse from each of g1, g5 and g6", gray),
        (.25, .58, "Retained outer mice\nCutoffs and both sets of predictors", gray),
        (.75, .58, "Retained outer mice\nInner held-out triplets give training outcomes", gray),
        (.25, .29, "Fit binary evidence ridge", violet),
        (.75, .29, "Fit continuous evidence ridge", peach),
        (.50, .07, "Score the same test segments against the observed joint contrast", gray),
    ], [(0, 1), (0, 2), (1, 3), (2, 4), (3, 5), (4, 5)])
    diagram(output / "paired_external.pdf", 1.75, [
        (.21, .79, "Matched internal spleen panel\nFix candidates and predictors", gray),
        (.21, .30, "Freeze both fitted models\nfrom internal triplets", gray),
        (.56, .79, "Skin-graft cohort\nIntact and graft-arm TRB", gray),
        (.56, .30, "Marrow-graft cohort\nRecipient TRA and TRB", gray),
        (.87, .79, "Joint contrast\nand error", peach),
        (.87, .30, "Recipient recovery\nonly", violet),
    ], [(0, 2), (1, 3), (2, 4), (3, 5)])


if __name__ == "__main__":
    main()
