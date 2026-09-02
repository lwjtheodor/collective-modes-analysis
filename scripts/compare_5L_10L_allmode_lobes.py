"""All-mode CJJ first-lobe comparison: 5L n=1..8 versus 10L n=1..10."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'assets'; RAW=ROOT/'results'/'collective_mode_response'/'fig2_longitudinal_modes_88_rh75_330k'/'2026-08-11'/'remote_raw'/'output'
CSV=A/'cjj_88_5L_10L_allmode_first_lobes_vs_invk.csv'; MATCH=A/'cjj_88_5L_10L_matched_k_first_lobe_ratios.csv'; STEM=A/'cjj_88_5L_10L_allmode_first_lobes_vs_invk_nature'
def sem(x): return np.std(x,ddof=1)/np.sqrt(len(x))
def collect(L,nmax):
    rows=[]
    for f in sorted(RAW.glob(f'8_8_L{L}_rep*_CJJ_alln.json')):
        m=json.loads(f.read_text()); rep=int(m['case_id'].split('_')[-1][3:])
        for z in m['mode_summary']:
            if z['n']>nmax: continue
            q=z['first_negative_lobe_normalized']; rows.append({'L':L,'replicate':rep,'n':z['n'],'k_inv_A':z['k_inv_A'],'invk_nm':.1/z['k_inv_A'],'area_ps':q['negative_area_normalized_ps'],'duration_ps':q['t_end_ps']-q['t_start_ps'],'depth':q['depth_normalized'],'t_start_ps':q['t_start_ps'],'t_end_ps':q['t_end_ps']})
    d=pd.DataFrame(rows)
    return d.groupby(['L','n','k_inv_A','invk_nm'],as_index=False).agg(area_mean_ps=('area_ps','mean'),area_sem_ps=('area_ps',sem),duration_mean_ps=('duration_ps','mean'),duration_sem_ps=('duration_ps',sem),depth_mean=('depth','mean'),depth_sem=('depth',sem),t_start_mean_ps=('t_start_ps','mean'),t_end_mean_ps=('t_end_ps','mean'))
def main():
    a5=collect(5,8); a10=collect(10,10); d=pd.concat([a5,a10],ignore_index=True).sort_values(['L','invk_nm']); d['definition']='CJJ/CJJ(0); first complete lobe bounded by first down/up zero crossings'; d['protocol']='weak NH; no global momentum removal; instantaneous water-COM axial subtraction'; d['sampling']=np.where(d.L==5,'10 fs / 1 ns','100 fs / 10 ns'); d['replicas']=3; d['uncertainty']='replica SEM'; d.to_csv(CSV,index=False,float_format='%.9g')
    pairs=[]
    for n in range(1,5):
        p=a5[a5.n==n].iloc[0]; q=a10[a10.n==2*n].iloc[0]
        for metric in ['area_mean_ps','duration_mean_ps','depth_mean']:
            pairs.append({'matched_pair':f'5L n={n} vs 10L n={2*n}','invk_nm':p.invk_nm,'metric':metric,'fiveL_mean':p[metric],'tenL_mean':q[metric],'tenL_over_fiveL':q[metric]/p[metric],'relative_difference_percent':100*(q[metric]/p[metric]-1)})
    pd.DataFrame(pairs).to_csv(MATCH,index=False,float_format='%.9g')
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,3,figsize=(7.25,2.8),constrained_layout=True); c5='#c3842d'; c10='#2775a9'
    specs=[('area_mean_ps','area_sem_ps',r'negative area, $A_-^{C(0)}$ (ps)','a'),('duration_mean_ps','duration_sem_ps',r'negative duration, $\tau_-$ (ps)','b'),('depth_mean','depth_sem',r'relative valley depth, $D_-$','c')]
    for ax,(m,e,y,tag) in zip(axs,specs):
        ax.errorbar(a5.invk_nm,a5[m],a5[e],color=c5,marker='s',ms=4.5,lw=1.1,capsize=2,label=r'5L modes $n=1\ldots8$')
        ax.errorbar(a10.invk_nm,a10[m],a10[e],color=c10,marker='o',ms=4,lw=1.1,capsize=2,label=r'10L modes $n=1\ldots10$')
        ax.set(xscale='log',xlabel=r'$1/k$ (nm)',ylabel=y,xlim=(.85,20),xticks=[1,2,5,10]); ax.get_xaxis().set_major_formatter(plt.ScalarFormatter()); ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,f'({tag})',transform=ax.transAxes,fontweight='bold',fontsize=10)
    axs[0].legend(frameon=False,fontsize=6.0,loc='upper left')
    fig.suptitle(r'$(8,8)$ water: all-mode first-lobe check of the 10L $n=1$ continuation',y=1.03,fontsize=10)
    fig.text(.5,-.075,r'Same weak-NH/no-global-momentum-removal protocol and lobe extractor.  5L: 10 fs / 1 ns, modes $n=1\ldots8$; 10L: 100 fs / 10 ns, modes $n=1\ldots10$.  Both use water-COM axial subtraction; 3 replicas; error bars = SEM.',ha='center',fontsize=5.9,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
