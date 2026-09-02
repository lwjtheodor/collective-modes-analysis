"""Reconcile the protocol and normalization sources behind two (8,8) CJJ figures."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'assets'; OUT.mkdir(exist_ok=True)
LOW=OUT/'lowfreq_cvj_n1_lobe_normalization_summary.csv'
WEAK=ROOT/'results'/'collective_mode_response'/'fig2_longitudinal_modes_88_rh75_330k'/'2026-08-11'/'derived_data'/'panel_b_lowk_strength.csv'
CSV=OUT/'cjj_88_baseline_vs_weaknh_reconciliation.csv'
STEM=OUT/'cjj_88_baseline_vs_weaknh_reconciliation_nature'

def main():
    lo=pd.read_csv(LOW).query('box_length >= 2 and box_length <= 5').copy()
    we=pd.read_csv(WEAK).query('L >= 2 and L <= 5').copy()
    records=[]
    mapping=[('negative area', 'area_ps_mean','area_ps_sem','A_minus_ps_mean','A_minus_ps_sem',r'$A_-^{C(0)}$ (ps)'),
             ('negative duration','width_ps_mean','width_ps_sem','lobe_width_ps_mean','lobe_width_ps_sem',r'$\tau_-$ (ps)'),
             ('relative depth','depth_mean','depth_sem','depth_norm_mean','depth_norm_sem',r'$D_-$')]
    for metric,lm,le,wm,we_,label in mapping:
        for _,r in lo.iterrows(): records.append({'protocol':'baseline NVT','L':r.box_length,'metric':metric,'mean':r[lm],'sem':r[le],'definition':label})
        for _,r in we.iterrows(): records.append({'protocol':'weak NH/no momentum removal','L':r.L,'metric':metric,'mean':r[wm],'sem':r[we_],'definition':label})
    pd.DataFrame(records).to_csv(CSV,index=False)

    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,3,figsize=(7.25,2.55),constrained_layout=True)
    styles=[('baseline NVT','#c3842d','s'),('weak NH/no momentum removal','#2775a9','o')]
    for ax,(metric,lm,le,wm,we_,label),tag in zip(axs,mapping,'abc'):
        for prot,col,marker in styles:
            q=pd.DataFrame(records).query('metric == @metric and protocol == @prot').sort_values('L')
            ax.errorbar(q.L,q['mean'],q['sem'],color=col,marker=marker,ms=4.5,lw=1.2,capsize=2.3,label=prot)
        ax.set(xlabel='box length, $L$',ylabel=label,xlim=(1.75,5.25),xticks=[2,3,4,5])
        ax.spines[['top','right']].set_visible(False)
        ax.text(-.15,1.04,f'({tag})',transform=ax.transAxes,fontweight='bold',fontsize=10)
    axs[0].legend(frameon=False,fontsize=6.8,loc='upper left')
    fig.suptitle(r'$(8,8)$, $n=1$ current ACF: protocol reconciliation at identical zero-crossing lobe definition',y=1.04,fontsize=10)
    fig.text(.5,-.085,'Both: $C_{JJ}(t)/C_{JJ}(0)$, instantaneous water-COM axial velocity removed; mean ± replica SEM.  Baseline: 1 ps / 20 ns NVT.  Weak NH: 10 fs / 1 ns, no global momentum removal.',ha='center',fontsize=6.5,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight')
    fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight')
    fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight')
    fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__': main()
