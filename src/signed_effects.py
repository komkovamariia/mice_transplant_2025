"""Shared single-panel mean signed-effect figure for edgeR and Fisher."""

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

CORAL = "#FF785C"
PURPLE = "#792080"


def signed_effect_panel(frame, chain, effect, count, xlabel, limit=8):
    tested = frame.loc[frame[count].gt(0)].dropna(subset=[effect])
    positive = tested.loc[tested[effect].gt(0)].sort_values(
        [effect, "v_gene"], ascending=[False, True]
    ).head(8)
    negative = tested.loc[tested[effect].lt(0)].sort_values([effect, "v_gene"]).head(8)
    rows = pd.concat([positive, negative])
    with plt.rc_context({"font.family": "serif", "font.serif": ["Liberation Serif", "DejaVu Serif"],
                         "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False}):
        height = 1.1 + .18125 * max(len(rows), 1)
        fig, ax = plt.subplots(figsize=(7, height))
        fig.subplots_adjust(left=.24, right=.965, bottom=.68/height, top=1-.48/height)
        ax.barh(rows["v_gene"], rows[effect],
                color=np.where(rows[effect].gt(0), CORAL, PURPLE), height=.78)
        if rows.empty:
            ax.text(.5, .5, "No nonzero mean effects.", ha="center", va="center",
                    transform=ax.transAxes)
        ax.invert_yaxis()
        ax.set_ylabel(f"{chain} V segment")
        ax.set_xlabel(xlabel)
        ax.set_xlim(-limit, limit)
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.grid(axis="x", color="#ededed", linewidth=.5)
        ax.set_axisbelow(True)
        ax.axvline(0, color="#555555", linewidth=.75)
        ax.legend(handles=[Patch(color=PURPLE, label="Higher in pooled g5 and g6"),
                           Patch(color=CORAL, label="Higher in g1")],
                  loc="lower center", bbox_to_anchor=(.5, 1.02), ncol=2,
                  frameon=False, fontsize=9)
        ax.tick_params(axis="both", labelsize=9, width=.6)
        for spine in ax.spines.values():
            spine.set_linewidth(.7)
    return fig

