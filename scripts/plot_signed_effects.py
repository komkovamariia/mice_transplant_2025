"""Draw V-segment mean signed effects using the same diverging-bar layout."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.signed_effects import signed_effect_panel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fisher-summary", type=Path, required=True)
    parser.add_argument("--edger-summary", type=Path, help="Optional edgeR V-level summary")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fisher = pd.read_csv(args.fisher_summary)
    specs = [(fisher, "mean_log2OR", "n_significant", "Fisher log₂ odds ratio", "fisher")]
    if args.edger_summary:
        specs.append((pd.read_csv(args.edger_summary), "mean_edger_log2fc",
                      "n_edger_significant", "edgeR log₂ fold change", "edger"))
    limit = max(8, *(int(np.ceil(frame[effect].abs().max()/2))*2
                      for frame, effect, *_ in specs))
    for frame, effect, count, label, name in specs:
        for chain in ("TRA", "TRB"):
            fig = signed_effect_panel(
                frame.loc[frame.chain.eq(chain)], chain, effect, count,
                f"Mean {label} (g1 relative to pooled g5 and g6)", limit,
            )
            fig.savefig(args.output / f"{chain}_{name}.png", dpi=300, facecolor="white")
            plt.close(fig)


if __name__ == "__main__":
    main()
