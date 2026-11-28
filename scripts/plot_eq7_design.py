#!/usr/bin/env python3
"""Draw the model encodings, nested holdout and external evaluation schemes."""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, PathPatch
from matplotlib.path import Path as MplPath
import numpy as np

purple = "#854696"
peach = "#f6a084"
violet_bg = "#efe3f1"
peach_bg = "#fceadd"
mint = "#ddefe9"
ink = "black"

def canvas(w,h):
    fig,ax=plt.subplots(figsize=(w,h));fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.99)
    ax.set_xlim(0,w);ax.set_ylim(0,h);ax.axis('off')
    return fig,ax

def box(ax,x,y,w,h,label,fill='white',size=9.2):
    ax.add_patch(Rectangle((x,y),w,h,fc=fill,ec=ink,lw=.75))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=size,linespacing=1.15)

def circles(ax,x,y,n,first):
    for i in range(n):
        ax.add_patch(Circle((x+i*.22,y),.081,fc=first if i==0 else 'white',ec=ink,lw=.8))

def save(fig, name, output):
    fig.savefig(output / name, facecolor='white', dpi=300)
    plt.close(fig)

def draw_schemes(output):
    fig, ax = canvas(7.5, 3.7)
    rows = (3.18, 2.58, 1.98)
    for y, label, count in zip(rows, ("Outer omission", "Inner omission", "Inner reference"),
                               ("3 × 5 × 5 = 75", "2 × 4 × 4 = 32", "1 + 3 + 3 mice")):
        ax.text(.12, y, label, va="center")
        ax.text(6.12, y, count, va="center")
    for x, n, group in ((2.03, 3, "g1"), (3.45, 5, "g5"), (5.00, 5, "g6")):
        ax.text(x + (n-1)*.11, 3.5, group, ha="center")
        for level, (y, fill) in enumerate(zip(rows, (peach, purple, "white"))):
            circles(ax, x, y, n-level, fill)
            if level < 2:
                for j in range(1, n-level):
                    ax.plot([x+j*.22, x+(j-1)*.22], [y-.08, rows[level+1]+.08],
                            color="black", lw=.65, ls=(0, (3, 2)), zorder=0)
    ax.text(3.75, 1.69, "(a)", ha="center")
    ax.text(3.75, 1.40, r"Average the 32 inner responses per V segment. Fit ridge at $\lambda=1$ (sensitivity analysis: 10).",
            ha="center", fontsize=9)
    for x, label, fill in zip((.1, 2.58, 5.06),
                              ("aaV-indicator ridge", "V-segment-effect ridge", "aaV-fraction ridge"),
                              (violet_bg, peach_bg, mint)):
        box(ax, x, .84, 2.33, .38, label, fill)
    ax.text(3.75, .61, "Comparators: structural-only, edgeR-only and training mean", ha="center", fontsize=9)
    ax.text(3.75, .36, "Evaluate predictions and rankings on the outer test mice", ha="center", fontsize=9)
    ax.text(3.75, .10, "(b)", ha="center")
    save(fig, 'nested.png', output)

    fig,ax=canvas(7.5,2.11)
    ax.text(.85,1.9,'Matched internal panel',ha='center')
    ax.text(2.66,1.9,'Frozen parameters',ha='center')
    ax.text(4.7,1.9,'External cohorts',ha='center')
    ax.text(6.53,1.9,'Evaluation',ha='center')
    left, right, center_y, radius = .61, 1.02, 1.04, .38
    ax.add_patch(Circle((left, center_y), radius, fc="#ffdcaf", ec=ink))
    ax.add_patch(Circle((right, center_y), radius, fc="#ffc7b8", ec=ink, alpha=1))
    for center in (left, right):
        ax.add_patch(Circle((center, center_y), radius, fc="none", ec=ink, lw=.8))
    theta = np.arccos((right-left)/(2*radius))
    left_arc = np.linspace(theta, 2*np.pi-theta, 100)
    right_arc = np.linspace(np.pi+theta, np.pi-theta, 100)
    points = np.vstack([
        np.column_stack([left+radius*np.cos(left_arc), center_y+radius*np.sin(left_arc)]),
        np.column_stack([right+radius*np.cos(right_arc), center_y+radius*np.sin(right_arc)]),
    ])
    crescent = MplPath(np.vstack([points, points[0]]),
                       [MplPath.MOVETO]+[MplPath.LINETO]*(len(points)-1)+[MplPath.CLOSEPOLY])
    ax.add_patch(PathPatch(crescent, fc="none", ec=ink, hatch="///", lw=0))
    ax.text(left-.055, center_y, "g1", ha="center", va="center",
            bbox={"facecolor": "#ffdcaf", "edgecolor": "none", "pad": 1})
    ax.text(right+.055, center_y, "g5", ha="center", va="center")
    ax.text(.82,.53,'Spleen panels',ha='center')
    ax.annotate('',(1.83,1.04),(1.45,1.04),arrowprops={'arrowstyle':'->','color':ink})
    ax.text(2.74,1.33,r'$f^*, r^*$: matched panel',ha='center',fontsize=9.5)
    ax.text(2.74,1.08,r'$\alpha_m,\ \gamma_m\ (m=1,2,3)$',ha='center',fontsize=10.3)
    ax.text(2.74,.76,'75 internal triplets\nMeans, SDs and coefficients',ha='center',fontsize=8.5)
    for y in (1.38,.56):
        ax.annotate('',(3.96,y),(3.62,1.04),arrowprops={'arrowstyle':'->','color':ink})
    box(ax,4.02,1.11,1.47,.59,'Zheng\nMarrow recipients\nTRA and TRB',violet_bg,8.8)
    box(ax,4.02,.25,1.47,.59,'Tereshchenko\nIntact and skin graft\nCD4+ TRB',peach_bg,8.8)
    for y in (1.40,.54):ax.annotate('',(5.96,y),(5.52,y),arrowprops={'arrowstyle':'->','color':ink})
    ax.text(6.70,1.39,'Recipient-only data\nNo two-arm evaluation',ha='center',va='center',fontsize=8.9)
    ax.text(6.70,.54,'Shortlist overlap\nContrast, ranking, MAE',ha='center',va='center',fontsize=8.9)
    for label in ax.texts:
        label.set_fontsize(label.get_fontsize() * .81)
    save(fig, 'external.png', output)

    fig,ax=canvas(7.5,2.22)
    ax.text(.18,1.96,'Model',ha='left')
    ax.text(2.68,1.96,'Evidence',ha='center');ax.text(5.05,1.96,'Input to ridge',ha='center');ax.text(6.95,1.96,'Range',ha='center')
    for y in (2.17,1.82,1.28,.65,.16):
        ax.plot([.05,7.45],[y,y],color=ink,lw=.65)
    for x in (.05,1.65,3.43,6.44,7.45):
        ax.plot([x,x],[.16,2.17],color=ink,lw=.65)
    ax.text(.1,1.56,'aaV-indicator\nridge',va='center')
    ax.text(.1,.99,'V-segment-effect\nridge',va='center')
    ax.text(.1,.41,'aaV-fraction\nridge',va='center')
    def symbols(y,values):
        for i,v in enumerate(values):
            x=2.08+i*.25;fill={'+' :purple,'−':peach,'':'white'}[v]
            ax.add_patch(Circle((x,y),.088,fc=fill,ec=ink,lw=.5));ax.text(x,y,v,ha='center',va='center',fontsize=10)
    symbols(1.56,['+','+','−','','']);symbols(.41,['+','+','','',''])
    for y,label,vals in [(1.09,'g1',[.35,.22,.14,.11,.04]),(.91,'Graft',[.23,.12,.08,.04])]:
        ax.text(2.04,y,label,ha='right',va='center',fontsize=8)
        x=2.1
        for i,val in enumerate(vals):
            ax.add_patch(Rectangle((x,y-.045),val,.09,fc=[purple,'#b78ac7',peach,'#e9bedf','#f9eee8'][i],ec='none'));x+=val
    ax.text(5.05,1.56,'Positive mean effect → 1',ha='center',va='center',fontsize=9.1)
    ax.text(5.05,.99,'Signed effect × (1 − FDR)\nedgeR: $\\log_2$FC, Fisher: $\\log_2$OR',ha='center',va='center',fontsize=9.1)
    ax.text(5.05,.41,'2 of 5 positive and significant\nFraction = 0.40',ha='center',va='center',fontsize=9.1)
    for y,label in [(1.56,'{0, 1}'),(.99,'(−∞, +∞)'),(.41,'[0, 1]')]:ax.text(6.95,y,label,ha='center',va='center')
    save(fig, 'predictors.png', output)



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
                         "mathtext.fontset": "cm", 
                         "font.size": 10})
    draw_schemes(args.output)


if __name__ == "__main__":
    main()
