#!/usr/bin/env python3
"""Draw the seven-rule selection benchmark and the three-model ridge comparison."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PURPLE = "#854696"
PEACH = "#f6a084"
GRAY = "#686868"


def configure_style():
    plt.rcParams.update({"font.family": "serif", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})


def save(fig, path):
    fig.savefig(path, facecolor="white", dpi=300)
    plt.close(fig)


def selection_plot(summary, output):
    purple, peach, gray = PURPLE, PEACH, GRAY
    rules=[('structural_screen','Adaptive two-cutoff shortlist',peach),('S_ungated','Structural rank S(v)',peach),('initial_share_ungated','Initial V-segment share',gray),('edgeR_binary_ungated','edgeR aaV indicator',purple),('Fisher_binary_ungated','Fisher aaV indicator',purple),('V_usage_edgeR_logFC','edgeR V-segment effect',purple),('V_usage_edgeR_significant','edgeR V-segment effect + FDR',purple)]
    fig,axs=plt.subplots(1,2,figsize=(8.25,3.05),sharey=True,layout='constrained')
    for i,(ax,chain) in enumerate(zip(axs,('TRA','TRB'))):
        row=summary[summary.chain==chain].set_index('method')
        vals=[100*row.loc[key,'g1_minus_mean_g5g6'] for key,_,_ in rules]
        for y,(val,(_,label,c)) in enumerate(zip(vals,rules)):
            ax.hlines(y,0,val,color=c,alpha=.5,lw=1)
            ax.scatter(val,y,c=c,s=31,zorder=3)
            ax.annotate(f'{val:.2f}',(val,y),xytext=(5 if val>=0 else -5,0),textcoords='offset points',ha='left' if val>=0 else 'right',va='center',fontsize=9)
        ax.set_title(f'({"ab"[i]}) {chain}',loc='left',fontsize=11)
        ax.set(yticks=np.arange(len(rules)),yticklabels=[label for _,label,_ in rules],ylim=(len(rules)-.52,-.48),xlabel='Unique-aaV recovery contrast (pp)')
        ax.set_xlim(min(0,min(vals)*1.25),max(0,max(vals)*1.23))
        ax.axvline(0,color='#777',lw=.7);ax.grid(axis='x',color='#e5e0e5',lw=.55);ax.set_axisbelow(True)
        if i:ax.tick_params(axis='y',left=False,labelleft=False)
    save(fig, output / 'selection_benchmark.png')



def prediction_plot(rank, err, ext, pred, output):
    methods=[('structural_score','Structural rank S(v)','#4c4c4c','D'),('continuous_edger','edgeR V effect','#9a92b0','o'),('continuous_fisher','Fisher V effect','#d9a0ad','s'),('binary_evidence_ridge','aaV-indicator ridge','#743b86','o'),('continuous_evidence_ridge','V-segment-effect ridge','#ee8c61','s'),('aav_positive_fraction_ridge','aaV-fraction ridge','#479b88','^')]
    fig,axs=plt.subplots(2,2,figsize=(8.2,5.4),layout='constrained')
    for i,(ax,chain) in enumerate(zip(axs[0],('TRA','TRB'))):
        frame=rank[rank.chain==chain].pivot(index='model',columns='top_k',values='joint_recovery_contrast')
        positions=np.arange(len(methods))
        for j,k in enumerate((5,10) if chain=='TRA' else (2,4)):
            ax.scatter(frame.loc[[m[0] for m in methods],k],positions+(.12 if j else -.12),c=[m[2] for m in methods],s=26,marker='s' if j else 'o',label=f'Top {k}')
        ax.set(yticks=positions,yticklabels=[m[1] for m in methods],ylim=(len(methods)-.55,-.45),xlabel='Observed joint recovery contrast')
        ax.set_title(f'({"ab"[i]}) {chain} ranking',loc='left',fontsize=10.8);ax.grid(axis='x',color='#e5e0e5',lw=.5)
        if i:ax.tick_params(labelleft=False)
        ax.legend(frameon=False,fontsize=8,loc='best')
    ax=axs[1,0]
    models=methods[3:]
    labels=['Internal TRA','Internal TRB','External TRB']
    for i,(key,label,c,marker) in enumerate(models):
        a=err[(err.model==key)&(err.penalty==1)].set_index('chain')
        e=ext[(ext.model==key)&(ext.penalty==1)].iloc[0]
        v=[a.loc['TRA','mae']/a.loc['TRA','constant_mae'],a.loc['TRB','mae']/a.loc['TRB','constant_mae'],e.mae/e.constant_mae]
        ax.scatter(v,np.arange(3)+.16*(i-1),c=c,marker=marker,label=label,s=28)
    ax.axvline(1,color='#555',ls=':',lw=.8)
    ax.set(yticks=range(3),yticklabels=labels,ylim=(2.5,-.5),xlabel='MAE relative to training-mean baseline')
    ax.set_title('(c) Numerical error',loc='left',fontsize=10.8)
    ax.legend(frameon=False,fontsize=7,loc='best');ax.grid(axis='x',color='#e5e0e5',lw=.5)
    ax=axs[1,1]
    ex=pred.query('study=="PRJNA954716" & penalty==1')
    frame=ex.pivot(index='v_gene',columns='model',values='prediction')
    obs=ex.drop_duplicates('v_gene').set_index('v_gene').observed
    for i,(key,label,c,marker) in enumerate(models):
        ax.scatter(frame[key],np.arange(len(frame))+.15*(i-1),c=c,marker=marker,s=29,label=label)
    ax.scatter(obs.reindex(frame.index),np.arange(len(frame)),marker='x',c='black',s=36,label='Observed')
    ax.set(yticks=np.arange(len(frame)),yticklabels=frame.index,ylim=(len(frame)-.55,-.45),xlabel='External joint recovery contrast')
    ax.set_title('(d) External TRB',loc='left',fontsize=10.8)
    ax.legend(frameon=False,fontsize=7,loc='best');ax.grid(axis='x',color='#e5e0e5',lw=.5)
    save(fig, output / 'ridge_benchmark.png')


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
    selection_plot(pd.read_csv(args.selection / "internal_1249_selection_summary.csv"), args.output)
    prediction_plot(
        pd.read_csv(args.internal / "internal_ranking_summary.csv"),
        pd.read_csv(args.internal / "internal_error_summary.csv"),
        pd.read_csv(args.external / "external_error_summary.csv"),
        pd.read_csv(args.external / "external_predictions.csv"),
        args.output,
    )


if __name__ == "__main__":
    main()
