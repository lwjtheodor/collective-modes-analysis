"""Diagnose apparent linearity of 10L current-mode lobe metrics versus 1/k."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'assets'; IN=A/'cjj_88_n1_baseline_and_10L_all_modes_vs_invk.csv'; CSV=A/'cjj_88_10L_invk_lobe_scaling_diagnostic.csv'; STEM=A/'cjj_88_10L_invk_lobe_scaling_diagnostic_nature'
def fit(x,y):
    # The support is strictly positive; epsilon only protects the log transform numerically.
    if not (np.all(np.asarray(x)>0) and np.all(np.asarray(y)>0)): raise ValueError('log-log fit requires positive support')
    p,b=np.polyfit(np.log(np.maximum(x,1e-300)),np.log(np.maximum(y,1e-300)),1); return p,b
def main():
    d=pd.read_csv(IN); ten=d[d.source.str.startswith('10L')].copy().sort_values('invk_nm'); base=d[d.source.str.startswith('baseline')].copy().sort_values('invk_nm')
    for q in [ten,base]:
        q['area_per_invk_ps_per_nm']=q.area_mean_ps/q.invk_nm; q['duration_per_invk_ps_per_nm']=q.duration_mean_ps/q.invk_nm; q['geometric_shape_S']=q.area_mean_ps/(q.duration_mean_ps*q.depth_mean)
    pA,bA=fit(ten.invk_nm,ten.area_mean_ps); pT,bT=fit(ten.invk_nm,ten.duration_mean_ps); pD,bD=fit(ten.invk_nm,ten.depth_mean)
    ten['fit_area_exponent']=pA; ten['fit_duration_exponent']=pT; ten['fit_depth_exponent']=pD; pd.concat([ten,base],ignore_index=True).to_csv(CSV,index=False,float_format='%.9g')
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,3,figsize=(7.25,2.75),constrained_layout=True); blue='#2775a9'; gold='#c3842d'
    panels=[('area_per_invk_ps_per_nm',r'$A_-/(1/k)$ (ps nm$^{-1}$)',pA,'a'),('duration_per_invk_ps_per_nm',r'$\tau_-/(1/k)$ (ps nm$^{-1}$)',pT,'b'),('depth_mean',r'$D_-$',pD,'c')]
    for ax,(m,y,p,tag) in zip(axs,panels):
        ax.errorbar(ten.invk_nm,ten[m],ten['area_sem_ps']/ten.invk_nm if m.startswith('area') else (ten['duration_sem_ps']/ten.invk_nm if m.startswith('duration') else ten['depth_sem']),color=blue,marker='o',ms=4,lw=1.1,capsize=2,label=r'10L, $n=1\ldots10$')
        ax.errorbar(base.invk_nm,base[m],base['area_sem_ps']/base.invk_nm if m.startswith('area') else (base['duration_sem_ps']/base.invk_nm if m.startswith('duration') else base['depth_sem']),color=gold,marker='s',ms=4.5,lw=1.1,capsize=2,label='2–5L baseline, $n=1$')
        ax.set(xscale='log',xlabel=r'$1/k$ (nm)',ylabel=y,xlim=(.135,20),xticks=[.2,.5,1,2,5,10]); ax.get_xaxis().set_major_formatter(plt.ScalarFormatter()); ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,f'({tag})',transform=ax.transAxes,fontweight='bold',fontsize=10)
        ax.text(.04,.92,rf'10L log--log exponent: {p:.2f}',transform=ax.transAxes,va='top',fontsize=6.2,color=blue)
    axs[0].legend(frameon=False,fontsize=6.1,loc='lower right')
    fig.suptitle(r'Why the raw $A_-$ and $\tau_-$ curves look linear: remove the explicit $1/k$ factor',y=1.03,fontsize=10)
    fig.text(.5,-.075,r'All values use $C_{JJ}/C_{JJ}(0)$ and the first complete zero-crossing-bounded lobe.  For 10L, $p(A_-)=1.32$, $p(\tau_-)=1.03$, $p(D_-)=0.26$: duration is linear in $1/k$, while area remains modestly superlinear.',ha='center',fontsize=6.0,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
