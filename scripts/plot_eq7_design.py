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
    diagram(output / "parallel_scores.pdf", 1.55, [
        (.50, .82, "Shared structural shortlist and log residual aaV count\nRaw counts for statistical evidence", gray),
        (.25, .44, "Binary evidence ridge\naaV-level edgeR and Fisher flags", violet),
        (.75, .44, "Continuous evidence ridge\nV-usage edgeR and Fisher effects", peach),
        (.50, .09, "Common endpoint: joint recovery contrast", gray),
    ], [(0, 1), (0, 2), (1, 3), (2, 3)])
    diagram(output / "paired_holdout.pdf", 2.40, [
        (.50, .86, "Outer test: leave out one mouse from each of g1, g5 and g6", gray),
        (.50, .60, "Retained outer mice provide cutoffs and both predictor sets\nInner held-out triplets provide training outcomes", gray),
        (.25, .32, "Fit binary evidence ridge", violet),
        (.75, .32, "Fit continuous evidence ridge", peach),
        (.50, .08, "Score the same test segments against the observed joint contrast", gray),
    ], [(0, 1), (1, 2), (1, 3), (2, 4), (3, 4)])
    diagram(output / "paired_external.pdf", 2.00, [
        (.50, .82, "Fix candidates and predictors in matched internal spleen panels\nFreeze both ridge fits trained on internal triplets", gray),
        (.25, .43, "Skin-graft cohort\nIntact and graft-arm TRB", gray),
        (.75, .43, "Marrow-graft cohort\nRecipient TRA and TRB", gray),
        (.25, .10, "Joint contrast and prediction error", peach),
        (.75, .10, "Descriptive recipient recovery", violet),
    ], [(0, 1), (0, 2), (1, 3), (2, 4)])



if __name__ == "__main__":
    main()
