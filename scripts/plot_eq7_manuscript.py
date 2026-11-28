#!/usr/bin/env python3
"""Draw the seven-rule selection benchmark and the three-model ridge comparison."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PURPLE = "#854696"
PEACH = "#f6a084"
GRAY = "#686868"
RED = "#b92835"


def configure_style():
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Liberation Serif", "DejaVu Serif"], "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})


def save(fig, path):
    fig.savefig(path, facecolor="white", dpi=300)
    plt.close(fig)


def outline(ax, x, y):
    ax.scatter(x, y, s=100, facecolors="none", edgecolors=RED,
               linewidths=1.15, zorder=7, clip_on=False)


def framed_legend(ax, ncol=2, size=8):
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.01), ncol=ncol,
              frameon=True, facecolor="white", edgecolor="#bcbcbc",
              framealpha=1, fontsize=size, borderpad=.45, columnspacing=1,
              handletextpad=.4, fancybox=False)


def selection_plot(summary, output, intervals=None):
    purple, peach, gray = PURPLE, PEACH, GRAY
    rules=[('structural_screen','Adaptive two-cutoff shortlist',peach),('S_ungated','Structural rank S(v)',peach),('initial_share_ungated','Initial V-segment share',gray),('edgeR_binary_ungated','edgeR aaV indicator',purple),('Fisher_binary_ungated','Fisher aaV indicator',purple),('V_usage_edgeR_logFC','edgeR V-segment effect',purple),('V_usage_edgeR_significant','edgeR V-segment effect + FDR',purple)]
    fig,axs=plt.subplots(1,2,figsize=(8.5,3.65),sharey=True,layout='constrained')
    for i,(ax,chain) in enumerate(zip(axs,('TRA','TRB'))):
        row=summary[summary.chain==chain].set_index('method')
        vals=[100*row.loc[key,'g1_minus_mean_g5g6'] for key,_,_ in rules]
        for y,(val,(_,label,c)) in enumerate(zip(vals,rules)):
            ax.hlines(y,0,val,color=c,alpha=.5,lw=1)
            ax.scatter(val,y,c=c,s=31,zorder=3)
            ax.annotate(f'{val:.2f}',(val,y),xytext=(10 if val>=0 else -10,0),textcoords='offset points',ha='left' if val>=0 else 'right',va='center',fontsize=9)
        for y in ((0, 1, 3) if chain == "TRA" else (0, 3)):
            outline(ax, vals[y], y)
        title = f"({'ab'[i]}) {chain}"
        if intervals is not None:
            interval = intervals.set_index("chain").loc[chain]
            title += (f"\nAdaptive − edgeR indicator: {100*interval.contrast:+.2f} pp"
                      f"\nConditional 95% CI [{100*interval.bootstrap_lo:.2f}, "
                      f"{100*interval.bootstrap_hi:.2f}]")
        ax.set_title(title, loc="left", fontsize=10.5, pad=7)
        ax.set(yticks=np.arange(len(rules)),yticklabels=[label for _,label,_ in rules],ylim=(len(rules)-.52,-.48),xlabel='Unique-aaV recovery contrast (pp)')
        ax.set_xlim(min(0,min(vals)*1.30),max(0,max(vals)*1.28))
        ax.axvline(0,color='#777',lw=.7);ax.grid(axis='x',color='#e5e0e5',lw=.55);ax.set_axisbelow(True)
        if i:ax.tick_params(axis='y',left=False,labelleft=False)
    save(fig, output / 'selection_benchmark.png')



def prediction_plot(rank, err, ext, pred, output):
    methods=[('structural_score','Structural rank S(v)','#4c4c4c','D'),('continuous_edger','edgeR V effect','#9a92b0','o'),('continuous_fisher','Fisher V effect','#d9a0ad','s'),('binary_evidence_ridge','aaV-indicator ridge','#743b86','o'),('continuous_evidence_ridge','V-segment-effect ridge','#ee8c61','s'),('aav_positive_fraction_ridge','aaV-fraction ridge','#479b88','^')]
    fig,axs=plt.subplots(2,2,figsize=(9.0,6.5),layout='constrained')
    for i,(ax,chain) in enumerate(zip(axs[0],('TRA','TRB'))):
        frame=rank[rank.chain==chain].pivot(index='model',columns='top_k',values='joint_recovery_contrast')
        positions=np.arange(len(methods))
        highlighted = ({("binary_evidence_ridge", 5), ("binary_evidence_ridge", 10),
                        ("structural_score", 10)} if chain == "TRA"
                       else {("continuous_edger", 2), ("aav_positive_fraction_ridge", 2)})
        for j, k in enumerate((5, 10) if chain == "TRA" else (2, 4)):
            vals = frame.loc[[m[0] for m in methods], k]
            ys = positions + (.12 if j else -.12)
            marker = "s" if j else "o"
            ax.scatter(vals, ys, c=[m[2] for m in methods], s=30, marker=marker)
            ax.scatter([], [], c="#555555", s=30, marker=marker, label=f"Top {k}")
            for idx, method in enumerate(methods):
                if (method[0], k) in highlighted:
                    outline(ax, vals.iloc[idx], ys[idx])
        ax.set(yticks=positions,yticklabels=[m[1] for m in methods],ylim=(len(methods)-.55,-.45),xlabel='Observed joint recovery contrast')
        ax.set_title(f'({"ab"[i]}) {chain} ranking',loc='left',fontsize=10.8,pad=27);ax.grid(axis='x',color='#e5e0e5',lw=.5)
        if i:ax.tick_params(labelleft=False)
        framed_legend(ax)
    ax=axs[1,0]
    models=methods[3:]
    labels=['Internal TRA','Internal TRB','External TRB']
    for i,(key,label,c,marker) in enumerate(models):
        a=err[(err.model==key)&(err.penalty==1)].set_index('chain')
        e=ext[(ext.model==key)&(ext.penalty==1)].iloc[0]
        v=[a.loc['TRA','mae']/a.loc['TRA','constant_mae'],a.loc['TRB','mae']/a.loc['TRB','constant_mae'],e.mae/e.constant_mae]
        ys = np.arange(3) + .16*(i-1)
        ax.scatter(v, ys, c=c, marker=marker, label=label, s=32)
        if key == "aav_positive_fraction_ridge":
            outline(ax, v[:2], ys[:2])
        if key == "continuous_evidence_ridge":
            outline(ax, v[2], ys[2])
    ax.axvline(1,color='#555',ls=':',lw=.8)
    ax.set(yticks=range(3),yticklabels=labels,ylim=(2.5,-.5),xlabel='MAE relative to training-mean baseline')
    ax.set_title('(c) Numerical error',loc='left',fontsize=10.8,pad=42)
    framed_legend(ax, size=7.5);ax.grid(axis='x',color='#e5e0e5',lw=.5)
    ax=axs[1,1]
    ex=pred.query('study=="PRJNA954716" & penalty==1')
    frame=ex.pivot(index='v_gene',columns='model',values='prediction')
    obs=ex.drop_duplicates('v_gene').set_index('v_gene').observed
    for i,(key,label,c,marker) in enumerate(models):
        ax.scatter(frame[key],np.arange(len(frame))+.15*(i-1),c=c,marker=marker,s=29,label=label)
    ax.scatter(obs.reindex(frame.index),np.arange(len(frame)),marker='x',c='black',s=36,label='Observed')
    outline(ax, obs["TRBV19"], frame.index.get_loc("TRBV19"))
    ax.set(yticks=np.arange(len(frame)),yticklabels=frame.index,ylim=(len(frame)-.55,-.45),xlabel='External joint recovery contrast')
    ax.set_title('(d) External TRB',loc='left',fontsize=10.8,pad=42)
    framed_legend(ax, size=7.5);ax.grid(axis='x',color='#e5e0e5',lw=.5)
    save(fig, output / 'ridge_benchmark.png')


def count_ranking_plot(rank, output):
    methods = [
        ("structure_only", "Structural-only ridge", "#4c4c4c"),
        ("edgeR_only_binary", "edgeR-only (binary)", "#9a92b0"),
        ("edgeR_only_fraction", "edgeR-only (fraction)", "#ad9db7"),
        ("binary_evidence_ridge", "aaV-indicator ridge", "#743b86"),
        ("continuous_evidence_ridge", "V-segment-effect ridge", "#ee8c61"),
        ("aav_positive_fraction_ridge", "aaV-fraction ridge", "#479b88"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5), layout="constrained")
    for i, (ax, chain) in enumerate(zip(axes, ("TRA", "TRB"))):
        frame = rank.loc[rank.chain.eq(chain)].pivot(index="model", columns="top_k",
                                                    values="count_contrast_pp")
        positions = np.arange(len(methods))
        for j, k in enumerate((5, 10) if chain == "TRA" else (2, 4)):
            marker = "s" if j else "o"
            ax.scatter(frame.loc[[m[0] for m in methods], k], positions+(.12 if j else -.12),
                       c=[m[2] for m in methods], s=30, marker=marker)
            ax.scatter([], [], c="#555555", s=30, marker=marker, label=f"Top {k}")
        ax.set(yticks=positions, yticklabels=[m[1] for m in methods],
               ylim=(len(methods)-.55, -.45), xlabel="Unique-aaV recovery contrast (pp)")
        ax.set_title(f"({'ab'[i]}) {chain} ranking", loc="left", fontsize=11, pad=27)
        ax.axvline(0, color="#777", lw=.7)
        ax.grid(axis="x", color="#e5e0e5", lw=.5)
        if i:
            ax.tick_params(labelleft=False)
        framed_legend(ax)
    save(fig, output / "ridge_count_ranking.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--internal", type=Path, required=True,
                        help="Combined current model results")
    parser.add_argument("--external", type=Path, required=True,
                        help="Combined current external predictions and summaries")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    configure_style()
    selection = args.selection / "internal_1249_selection_summary.csv"
    if not selection.exists():
        selection = args.selection / "selection_summary.csv"
    interval_file = args.selection / "primary_paired_bootstrap.csv"
    intervals = pd.read_csv(interval_file) if interval_file.exists() else None
    selection_plot(pd.read_csv(selection), args.output, intervals)
    rank = pd.read_csv(args.internal / "internal_ranking_summary.csv")
    prediction_plot(
        rank,
        pd.read_csv(args.internal / "internal_error_summary.csv"),
        pd.read_csv(args.external / "external_error_summary.csv"),
        pd.read_csv(args.external / "external_predictions.csv"),
        args.output,
    )
    count_ranking_plot(rank, args.output)


if __name__ == "__main__":
    main()

