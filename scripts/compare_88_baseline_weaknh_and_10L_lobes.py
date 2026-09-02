"""Same-estimator reconciliation: baseline NVT vs weak-NH, with 10L sensitivity."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'assets'
BASE=A/'lowfreq_cvj_n1_lobe_normalization_summary.csv'
WEAK=ROOT/'results'/'collective_mode_response'/'fig2_longitudinal_modes_88_rh75_330k'/'2026-08-11'/'derived_data'/'panel_b_lowk_strength.csv'
RAW=ROOT/'results'/'collective_mode_response'/'absolute_cvj_first_negative_lobe'/'2026-08-09'/'remote_raw'/'stage_absolute_cvj_1L10L'/'output'
CSV=A/'cjj_88_baseline_weaknh_10L_same_estimator_comparison.csv'; STEM=A/'cjj_88_baseline_weaknh_10L_same_estimator_comparison_nature'

def lobe(t,y):
    down=np.flatnonzero((y[:-1]>=0)&(y[1:]<0)); assert len(down); i0=down[0]
    im=i0+1+np.argmin(y[i0+1:]); up=np.flatnonzero((y[im:-1]<=0)&(y[im+1:]>0)); assert len(up); i1=im+up[0]
    z=lambda i:t[i]-y[i]*(t[i+1]-t[i])/(y[i+1]-y[i])
    start,end=z(i0),z(i1); tt=np.r_[start,t[i0+1:i1+1],end]; yy=np.r_[0,y[i0+1:i1+1],0]
    return -np.trapz(yy,tt),end-start,-y[im]
def sem(x): return np.std(x,ddof=1)/np.sqrt(len(x))
def main():
    b=pd.read_csv(BASE).query('box_length >= 2 and box_length <=5').copy()
    w=pd.read_csv(WEAK).query('L>=2').copy()
    rows=[]
    # The native values are the archived high-resolution weak-NH extraction.
    for _,r in w.iterrows():
        rows.append({'protocol':'weak NH/no global momentum removal','sampling':'native','L':int(r.L),'metric':'area','mean':r.A_minus_ps_mean,'sem':r.A_minus_ps_sem})
        rows.append({'protocol':'weak NH/no global momentum removal','sampling':'native','L':int(r.L),'metric':'duration','mean':r.lobe_width_ps_mean,'sem':r.lobe_width_ps_sem})
        rows.append({'protocol':'weak NH/no global momentum removal','sampling':'native','L':int(r.L),'metric':'depth','mean':r.depth_norm_mean,'sem':r.depth_norm_sem})
    for _,r in b.iterrows():
        rows.append({'protocol':'baseline NVT; z momentum removed every 5 ps','sampling':'native','L':int(r.box_length),'metric':'area','mean':r.area_ps_mean,'sem':r.area_ps_sem})
        rows.append({'protocol':'baseline NVT; z momentum removed every 5 ps','sampling':'native','L':int(r.box_length),'metric':'duration','mean':r.width_ps_mean,'sem':r.width_ps_sem})
        rows.append({'protocol':'baseline NVT; z momentum removed every 5 ps','sampling':'native','L':int(r.box_length),'metric':'depth','mean':r.depth_mean,'sem':r.depth_sem})
    # Re-extract the *same* weak-NH curves after taking every 100th point (1 ps).
    for L in range(2,6):
        vals=[]
        for rep in range(1,4):
            q=pd.read_csv(RAW/f'8_8_L{L}_rep{rep}_absolute_cvj.csv')
            q=q.iloc[::100].reset_index(drop=True); vals.append(lobe(q.lag_ps.to_numpy(),q.cvj_normalized.to_numpy()))
        for metric,j in [('area',0),('duration',1),('depth',2)]:
            x=np.array([v[j] for v in vals]); rows.append({'protocol':'weak NH/no global momentum removal','sampling':'1 ps decimated','L':L,'metric':metric,'mean':x.mean(),'sem':sem(x)})
    o=pd.DataFrame(rows); o['definition']='CJJ/CJJ(0); first complete lobe bounded by first down/up zero crossings; linear zero-crossing interpolation; trapezoidal area'
    o.to_csv(CSV,index=False,float_format='%.9g')
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(2,3,figsize=(7.25,4.45),constrained_layout=True)
    specs=[('area',r'$A_-^{C(0)}$ (ps)'),('duration',r'$\tau_-$ (ps)'),('depth',r'$D_-$')]
    for j,(metric,ylabel) in enumerate(specs):
        ax=axs[0,j]
        for protocol,color,marker in [('baseline NVT; z momentum removed every 5 ps','#c3842d','s'),('weak NH/no global momentum removal','#2775a9','o')]:
            q=o.query('protocol == @protocol and sampling == "native" and metric == @metric').sort_values('L')
            if protocol.startswith('weak'): 
                core=q.query('L <= 5'); ext=q.query('L == 10'); ax.errorbar(core.L,core['mean'],core['sem'],color=color,marker=marker,ms=4.5,lw=1.2,capsize=2,label='weak NH/no global momentum removal'); ax.errorbar(ext.L,ext['mean'],ext['sem'],color=color,marker=marker,mfc='white',ms=5.2,lw=1.1,ls='--',capsize=2,label='10L weak-NH extension')
            else: ax.errorbar(q.L,q['mean'],q['sem'],color=color,marker=marker,ms=4.5,lw=1.2,capsize=2,label='baseline NVT')
        ax.set(xlabel='box length, $L$',ylabel=ylabel,xlim=(1.7,10.5),xticks=[2,3,4,5,10]); ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,'abc'[j],transform=ax.transAxes,fontweight='bold',fontsize=10)
        if j==0: ax.legend(frameon=False,fontsize=5.4,loc='upper left')
        ax=axs[1,j]; native=o.query('protocol == "weak NH/no global momentum removal" and sampling == "native" and metric == @metric and L <= 5').sort_values('L'); dec=o.query('protocol == "weak NH/no global momentum removal" and sampling == "1 ps decimated" and metric == @metric').sort_values('L'); ratio=dec['mean'].to_numpy()/native['mean'].to_numpy(); ax.axhline(1,color='.45',lw=.8,ls='--'); ax.plot(native.L,ratio,color='#6d4c9b',marker='o',ms=4,lw=1.1); ax.set(xlabel='box length, $L$',ylabel=r'1 ps decimated / 10 fs native',xlim=(1.75,5.25),xticks=[2,3,4,5],ylim=(.86,1.14)); ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,'def'[j],transform=ax.transAxes,fontweight='bold',fontsize=10)
    fig.suptitle(r'$(8,8)$ $n=1$ current ACF: 10L sensitivity extension and calculation-cadence control',y=1.02,fontsize=10)
    fig.text(.5,-.04,'Top: identical zero-crossing estimator. Baseline NVT: 1 ps / 20 ns, z momentum removed every 5 ps. Weak NH: no global momentum removal; 2–5L 10 fs / 1 ns; 10L 100 fs / 10 ns. Bottom: same weak-NH CJJ curves re-extracted after 1 ps decimation.',ha='center',fontsize=5.8,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
