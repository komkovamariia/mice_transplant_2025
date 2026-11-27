#!/usr/bin/env python3
"""Draw the model encodings, nested holdout and external evaluation schemes."""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

purple = "#854696"
peach = "#f6a084"
violet_bg = "#efe3f1"
peach_bg = "#fceadd"
mint = "#ddefe9"
gray = "#686868"

def canvas(w,h):
    fig,ax=plt.subplots(figsize=(w,h));fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.99)
    ax.set_xlim(0,w);ax.set_ylim(0,h);ax.axis('off')
    return fig,ax

def box(ax,x,y,w,h,label,fill='white',size=9.2):
    ax.add_patch(Rectangle((x,y),w,h,fc=fill,ec='#aaa',lw=.75))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=size,linespacing=1.15)

def circles(ax,x,y,n,first):
    for i in range(n):
        ax.add_patch(Circle((x+i*.22,y),.081,fc=first if i==0 else 'white',ec='#292929',lw=.8))

def save(fig, name, output):
    fig.savefig(output / name, facecolor='white', dpi=300)
    plt.close(fig)

def draw_schemes(output):
    fig,ax=canvas(7.5,3.95)
    ax.text(.12,3.76,'(a)  Nested construction of training and test responses',weight='bold',fontsize=10.4)
    for name,x in zip(('g1','g5','g6'),(2.2,3.62,5.22)):
        ax.text(x+.22,3.48,name,ha='center',fontsize=10)
    ax.text(.12,3.15,'Outer omission',va='center')
    for x,n in zip((2.0,3.30,4.90),(3,5,5)):circles(ax,x,3.15,n,peach)
    ax.text(6.15,3.15,'3 × 5 × 5 = 75',ha='left',va='center',fontsize=9.2)
    ax.text(.12,2.76,'Retained: inner\nomission',va='center')
    for x,n in zip((2.0,3.30,4.90),(2,4,4)):circles(ax,x,2.76,n,purple)
    ax.text(6.15,2.76,'2 × 4 × 4 = 32',ha='left',va='center',fontsize=9.2)
    for x,n in zip((2.0,3.30,4.90),(1,3,3)):circles(ax,x,2.37,n,'white')
    ax.text(.12,2.37,'Inner reference',va='center')
    ax.text(6.15,2.37,'1 + 3 + 3 mice',ha='left',va='center',fontsize=9.2)
    ax.text(3.75,2.1,'Inner reference: prevalence = 1 for detected reference aaV',ha='center',fontsize=9)
    ax.text(3.75,1.88,'Outer reference: 2 g1 mice (prevalence = 0.5 or 1)',ha='center',fontsize=9)
    ax.plot([.1,7.4],[1.72,1.72],color='#b9b9b9',lw=.6)
    ax.text(.12,1.52,'(b)  Fit ridge models using inner responses; score outer mice',weight='bold',fontsize=10.4)
    ax.text(3.75,1.28,'One fixed outer predictor vector per V; aggregate 32 inner responses  •  λ = 1 (sensitivity: 10)',ha='center',fontsize=8.8)
    for x,label,fill in zip((.1,2.58,5.06),('aaV-indicator ridge','V-segment-effect ridge','aaV-fraction ridge'),(violet_bg,peach_bg,mint)):
        box(ax,x,.64,2.33,.47,label,fill)
    ax.text(3.75,.40,'Comparators: structural-only, edgeR-only and training mean',ha='center',fontsize=9)
    ax.text(3.75,.15,'Common held-out V segments and contrasts: MAE and top-k recovery',ha='center',fontsize=9)
    save(fig, 'nested.png', output)

    fig,ax=canvas(7.5,2.11)
    ax.text(.85,1.9,'Matched internal panel',ha='center')
    ax.text(2.66,1.9,'Frozen parameters',ha='center')
    ax.text(4.7,1.9,'External cohorts',ha='center')
    ax.text(6.53,1.9,'Evaluation',ha='center')
    for x,label,color in [(.61,'g1','#ffdcaf'),(1.02,'g5','#ffc7b8')]:
        ax.add_patch(Circle((x,1.04),.38,fc=color,ec='#999',alpha=.78));ax.text(x,1.04,label,ha='center',va='center')
    ax.text(.82,.53,'Spleen panels',ha='center')
    ax.annotate('',(1.83,1.33),(1.45,1.15),arrowprops={'arrowstyle':'->','color':gray})
    ax.text(2.74,1.33,r'$f^*, r^*$: matched panel',ha='center',fontsize=9.5)
    ax.text(2.74,1.08,r'$\alpha_m,\ \gamma_m$ $(m=1,2,3)$',ha='center',fontsize=10.3)
    ax.text(2.74,.76,'From 75 combined-panel triplets\nTraining means, SDs and coefficients',ha='center',fontsize=8.5)
    for y in (1.38,.56):
        ax.annotate('',(3.96,y),(3.62,1.04),arrowprops={'arrowstyle':'->','color':gray})
    box(ax,4.02,1.11,1.47,.59,'Zheng\nMarrow recipients\nTRA and TRB',violet_bg,8.8)
    box(ax,4.02,.25,1.47,.59,'Tereshchenko\nIntact and skin graft\nCD4+ TRB',peach_bg,8.8)
    for y in (1.40,.54):ax.annotate('',(5.96,y),(5.52,y),arrowprops={'arrowstyle':'->','color':gray})
    ax.text(6.70,1.39,'Recipient index only\nNo two-arm MAE',ha='center',va='center',fontsize=8.9)
    ax.text(6.70,.54,'Shortlist overlap\nContrast, ranking, MAE',ha='center',va='center',fontsize=8.9)
    save(fig, 'external.png', output)

    fig,ax=canvas(7.5,2.36)
    ax.text(.08,2.18,'One V segment: illustrative test results',fontsize=10.3)
    ax.text(2.68,1.96,'Evidence',ha='center');ax.text(5.05,1.96,'Input to ridge',ha='center');ax.text(6.95,1.96,'Range',ha='center')
    for y in (1.82,1.28,.65,.16):ax.plot([.08,7.40],[y,y],color='#d5ccd8',lw=.75)
    ax.text(.1,1.56,'aaV-indicator ridge',va='center')
    ax.text(.1,.99,'V-segment-effect ridge',va='center')
    ax.text(.1,.41,'aaV-fraction ridge',va='center')
    def symbols(y,values):
        for i,v in enumerate(values):
            x=2.08+i*.25;fill={'+' :purple,'−':peach,'':'white'}[v]
            ax.add_patch(Circle((x,y),.088,fc=fill,ec='#333',lw=.5));ax.text(x,y,v,ha='center',va='center',fontsize=10)
    symbols(1.56,['+','+','−','','']);symbols(.41,['+','+','','',''])
    for y,label,vals in [(1.09,'g1',[.35,.22,.14,.11,.04]),(.91,'Graft',[.23,.12,.08,.04])]:
        ax.text(2.04,y,label,ha='right',va='center',fontsize=8)
        x=2.1
        for i,val in enumerate(vals):
            ax.add_patch(Rectangle((x,y-.045),val,.09,fc=[purple,'#b78ac7',peach,'#e9bedf','#f9eee8'][i],ec='none'));x+=val
    ax.text(5.05,1.56,'Mean significant signed effect > 0 → 1',ha='center',va='center',fontsize=9.1)
    ax.text(5.05,.99,'Signed effect × (1 − FDR)\nedgeR: log₂FC; Fisher: log₂OR',ha='center',va='center',fontsize=9.1)
    ax.text(5.05,.41,'2/5 = 0.40 positive and significant',ha='center',va='center',fontsize=9.1)
    for y,label in [(1.56,'{0, 1}'),(.99,'(−∞, +∞)'),(.41,'[0, 1]')]:ax.text(6.95,y,label,ha='center',va='center')
    save(fig, 'predictors.png', output)



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.size": 10})
    draw_schemes(args.output)


if __name__ == "__main__":
    main()
