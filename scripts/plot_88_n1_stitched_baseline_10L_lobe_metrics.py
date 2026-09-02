"""Operational 2-5L baseline + 10L extension for (8,8) n=1 current lobe metrics."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'assets'
BASE=A/'cjj_88_n1_baseline_lobe_metrics_vs_L.csv'
WEAK=A/'cjj_88_n1_lobe_duration_shape_depth_vs_L.csv'
CSV=A/'cjj_88_n1_stitched_baseline_10L_lobe_metrics_vs_L.csv'; STEM=A/'cjj_88_n1_stitched_baseline_10L_lobe_metrics_vs_L_nature'
def main():
    b=pd.read_csv(BASE); w=pd.read_csv(WEAK).query('L == 10')
    rows=[]
    for _,r in b.iterrows(): rows.append({'L':int(r.L),'Lz_nm':r.Lz_nm,'source':'baseline NVT','source_cadence_duration':'1 ps / 20 ns','global_momentum_policy':'z momentum removed every 5 ps','area_mean_ps':r.negative_area_C0_ps_mean,'area_sem_ps':r.negative_area_C0_ps_sem,'duration_mean_ps':r.negative_duration_ps_mean,'duration_sem_ps':r.negative_duration_ps_sem,'depth_mean':r.relative_depth_mean,'depth_sem':r.relative_depth_sem,'geometric_shape_S_mean':r.geometric_shape_S_mean,'geometric_shape_S_sem':r.geometric_shape_S_sem})
    r=w.iloc[0]; rows.append({'L':10,'Lz_nm':r.Lz_nm,'source':'weak-NH extension appended to baseline','source_cadence_duration':'100 fs / 10 ns','global_momentum_policy':'no global momentum removal','area_mean_ps':r.negative_area_C0_ps_mean,'area_sem_ps':r.negative_area_C0_ps_sem,'duration_mean_ps':r.negative_duration_ps_mean,'duration_sem_ps':r.negative_duration_ps_sem,'depth_mean':r.negative_depth_relative_mean,'depth_sem':r.negative_depth_relative_sem,'geometric_shape_S_mean':r.area_over_duration_mean/r.negative_depth_relative_mean,'geometric_shape_S_sem':float('nan')})
    d=pd.DataFrame(rows); d['definition']='CJJ/CJJ(0); first complete negative lobe bounded by first down/up zero crossings'; d['replicas']=3; d['uncertainty']='replica SEM'; d['status']='operational stitched curve: 10L is appended, not a simulated 1 ps NVT point'; d.to_csv(CSV,index=False,float_format='%.9g')
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,3,figsize=(7.25,2.55),constrained_layout=True); col='#c3842d'
    specs=[('area_mean_ps','area_sem_ps',r'negative area, $A_-^{C(0)}$ (ps)','a'),('duration_mean_ps','duration_sem_ps',r'negative duration, $\tau_-$ (ps)','b'),('depth_mean','depth_sem',r'relative valley depth, $D_-$','c')]
    for ax,(m,e,y,t) in zip(axs,specs):
        ax.errorbar(d.L,d[m],d[e],color=col,marker='s',ms=5,lw=1.3,capsize=2.4)
        ax.plot(d.L,d[m],color=col,lw=1.3)
        ax.scatter([10],[d.loc[d.L==10,m].iloc[0]],s=45,marker='s',facecolors='white',edgecolors=col,zorder=4,lw=1.3)
        ax.set(xlabel='box length, $L$',ylabel=y,xlim=(1.7,10.6),xticks=[2,3,4,5,10]); ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,f'({t})',transform=ax.transAxes,fontweight='bold',fontsize=10)
    axs[0].annotate('10L appended\nprotocol extension',xy=(10,d.area_mean_ps.iloc[-1]),xytext=(6.0,10.4),fontsize=6.7,color=col,arrowprops=dict(arrowstyle='-',color=col,lw=.7))
    fig.suptitle(r'$(8,8)$ water: operational baseline–10L stitched $n=1$ current-ACF lobe series',y=1.04,fontsize=10)
    fig.text(.5,-.085,'2–5L: baseline NVT, z momentum removed every 5 ps, 1 ps / 20 ns.  10L (open symbol): weak NH/no global momentum removal, 100 fs / 10 ns.  All: water-COM axial subtraction; 3 replicas; error bars = SEM.',ha='center',fontsize=6.1,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
