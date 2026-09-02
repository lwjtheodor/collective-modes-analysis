"""Use all n=1..10 modes of the 10L trajectory on a 1/k abscissa."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'assets'
BASE=A/'cjj_88_n1_baseline_lobe_metrics_vs_L.csv'
RAW=ROOT/'results'/'collective_mode_response'/'fig2_longitudinal_modes_88_rh75_330k'/'2026-08-11'/'remote_raw'/'output'
CSV=A/'cjj_88_n1_baseline_and_10L_all_modes_vs_invk.csv'; STEM=A/'cjj_88_n1_baseline_and_10L_all_modes_vs_invk_nature'
def sem(x): return np.std(x,ddof=1)/np.sqrt(len(x))
def main():
    rows=[]
    for f in sorted(RAW.glob('8_8_L10_rep*_CJJ_alln.json')):
        m=json.loads(f.read_text()); rep=int(m['case_id'].split('_')[-1][3:])
        for z in m['mode_summary']:
            q=z['first_negative_lobe_normalized']; rows.append({'source':'10L weak NH/no global momentum removal','L':10,'replicate':rep,'n':z['n'],'k_inv_A':z['k_inv_A'],'invk_nm':0.1/z['k_inv_A'],'wavelength_nm':2*np.pi*0.1/z['k_inv_A'],'area_ps':q['negative_area_normalized_ps'],'duration_ps':q['t_end_ps']-q['t_start_ps'],'depth':q['depth_normalized'],'t_start_ps':q['t_start_ps'],'t_end_ps':q['t_end_ps'],'cadence_duration':'100 fs / 10 ns'})
    raw=pd.DataFrame(rows)
    ag=raw.groupby(['source','L','n','k_inv_A','invk_nm','wavelength_nm','cadence_duration'],as_index=False).agg(area_mean_ps=('area_ps','mean'),area_sem_ps=('area_ps',sem),duration_mean_ps=('duration_ps','mean'),duration_sem_ps=('duration_ps',sem),depth_mean=('depth','mean'),depth_sem=('depth',sem),t_start_mean_ps=('t_start_ps','mean'),t_end_mean_ps=('t_end_ps','mean'))
    b=pd.read_csv(BASE).query('L>=2 and L<=5').copy()
    br=pd.DataFrame({'source':'baseline NVT; z momentum removed every 5 ps','L':b.L,'n':1,'k_inv_A':b.k_inv_A,'invk_nm':0.1/b.k_inv_A,'wavelength_nm':2*np.pi*0.1/b.k_inv_A,'cadence_duration':'1 ps / 20 ns','area_mean_ps':b.negative_area_C0_ps_mean,'area_sem_ps':b.negative_area_C0_ps_sem,'duration_mean_ps':b.negative_duration_ps_mean,'duration_sem_ps':b.negative_duration_ps_sem,'depth_mean':b.relative_depth_mean,'depth_sem':b.relative_depth_sem,'t_start_mean_ps':np.nan,'t_end_mean_ps':np.nan})
    out=pd.concat([br,ag],ignore_index=True).sort_values(['source','invk_nm']); out['definition']='CJJ/CJJ(0); first complete negative lobe bounded by first down/up zero crossings'; out['replicas']=3; out['uncertainty']='replica SEM'; out.to_csv(CSV,index=False,float_format='%.9g')
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,3,figsize=(7.25,2.85),constrained_layout=True); gold='#c3842d'; blue='#2775a9'
    ten=out.query('source.str.startswith("10L")',engine='python').sort_values('wavelength_nm'); base=out.query('source.str.startswith("baseline")',engine='python').sort_values('wavelength_nm')
    specs=[('area_mean_ps','area_sem_ps',r'negative area, $A_-^{C(0)}$ (ps)','a'),('duration_mean_ps','duration_sem_ps',r'negative duration, $\tau_-$ (ps)','b'),('depth_mean','depth_sem',r'relative valley depth, $D_-$','c')]
    for ax,(m,e,y,tag) in zip(axs,specs):
        ax.errorbar(ten.wavelength_nm,ten[m],ten[e],color=blue,marker='o',ms=4,lw=1.15,capsize=2,label=r'10L modes $n=1\ldots10$')
        ax.errorbar(base.wavelength_nm,base[m],base[e],color=gold,marker='s',ms=4.5,lw=1.15,capsize=2,label='2–5L baseline $n=1$')
        for _,r in ten[ten.n.isin([2,5,10])].iterrows(): ax.annotate(rf'$n={int(r.n)}$',(r.wavelength_nm,r[m]),xytext=(2,4),textcoords='offset points',fontsize=5.6,color=blue)
        ax.set(xlabel=r'wavelength, $\lambda=2\pi/k$ (nm)',ylabel=y,xlim=(0,107),xticks=[0,25,50,75,100]); ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,f'({tag})',transform=ax.transAxes,fontweight='bold',fontsize=10)
    axs[0].legend(frameon=False,fontsize=6.1,loc='upper left')
    r1=ten[ten.n==1].iloc[0]
    axs[1].annotate(r'$n=1$: $\lambda=100.84$ nm',
                    xy=(r1.wavelength_nm,r1.duration_mean_ps),xytext=(36,27.5),fontsize=5.8,color=blue,
                    arrowprops=dict(arrowstyle='-',color=blue,lw=.7))
    fig.suptitle(r'$(8,8)$ water: first current-ACF negative lobe versus mode wavelength',y=1.02,fontsize=10)
    fig.text(.5,-.085,'Yellow: baseline NVT, z momentum removed every 5 ps, $n=1$, 1 ps / 20 ns. Blue: all 10L modes, weak NH/no global momentum removal, 100 fs / 10 ns. Water-COM axial subtraction; 3 replicas; error bars = SEM; labels = mode $n$.',ha='center',fontsize=5.85,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
