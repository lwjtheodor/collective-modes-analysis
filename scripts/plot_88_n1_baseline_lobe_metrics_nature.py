"""Canonical long-duration baseline-NVT (8,8) n=1 lobe summary."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'assets'/'lowfreq_cvj_n1_lobe_normalization_summary.csv'
OUT=ROOT/'assets'; CSV=OUT/'cjj_88_n1_baseline_lobe_metrics_vs_L.csv'; STEM=OUT/'cjj_88_n1_baseline_lobe_metrics_vs_L_nature'

def main():
    d=pd.read_csv(SRC).query('box_length >= 2 and box_length <= 5').copy()
    d['Lz_nm']=d['box_length']*10.084
    d['k_inv_A']=0.06230846322685951/d['box_length']
    d['observable']='first complete negative lobe of normalized n=1 axial current ACF'
    d['protocol']='baseline NVT; instantaneous water-COM axial velocity subtraction'
    d['cadence_duration']='1 ps cadence; 20 ns duration'
    d['replicas']=3; d['uncertainty']='replica SEM'
    canonical=pd.DataFrame({'L':d.box_length,'Lz_nm':d.Lz_nm,'k_inv_A':d.k_inv_A,
        'negative_area_C0_ps_mean':d.area_ps_mean,'negative_area_C0_ps_sem':d.area_ps_sem,
        'negative_duration_ps_mean':d.width_ps_mean,'negative_duration_ps_sem':d.width_ps_sem,
        'relative_depth_mean':d.depth_mean,'relative_depth_sem':d.depth_sem,
        'geometric_shape_S_mean':d.shape_factor_mean,'geometric_shape_S_sem':d.shape_factor_sem,
        'area_over_duration_mean':d.area_per_width_mean,'area_over_duration_sem':d.area_per_width_sem,
        'observable':d.observable,'protocol':d.protocol,'cadence_duration':d.cadence_duration,
        'replicas':d.replicas,'uncertainty':d.uncertainty})
    canonical.to_csv(CSV,index=False,float_format='%.9g')
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,3,figsize=(7.25,2.55),constrained_layout=True); col='#c3842d'
    specs=[('negative_area_C0_ps_mean','negative_area_C0_ps_sem',r'negative area, $A_-^{C(0)}$ (ps)','a'),('negative_duration_ps_mean','negative_duration_ps_sem',r'negative duration, $\tau_-$ (ps)','b'),('relative_depth_mean','relative_depth_sem',r'relative valley depth, $D_-$','c')]
    for ax,(m,e,y,t) in zip(axs,specs):
        ax.errorbar(canonical.L,canonical[m],canonical[e],color=col,marker='s',ms=5,lw=1.3,capsize=2.4)
        ax.set(xlabel='box length, $L$',ylabel=y,xlim=(1.75,5.25),xticks=[2,3,4,5])
        ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,f'({t})',transform=ax.transAxes,fontweight='bold',fontsize=10)
    fig.suptitle(r'$(8,8)$ water: canonical long-duration baseline, $n=1$ current ACF first negative lobe',y=1.04,fontsize=10)
    fig.text(.5,-.085,'Baseline NVT; 1 ps / 20 ns; instantaneous water-COM axial velocity subtraction; 3 replicas; error bars = replica SEM.  Zero-crossing-bounded complete lobe.',ha='center',fontsize=6.5,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
