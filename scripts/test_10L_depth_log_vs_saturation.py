"""Empirical model screen for the 10L all-mode relative depth D_minus(lambda)."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; A=ROOT/'assets'; IN=A/'cjj_88_n1_baseline_and_10L_all_modes_vs_invk.csv'; CSV=A/'cjj_88_10L_depth_lambda_log_vs_saturation.csv'; STEM=A/'cjj_88_10L_depth_lambda_log_vs_saturation_nature'
def r2(y,p): return 1-np.sum((y-p)**2)/np.sum((y-y.mean())**2)
def aic(y,p,npar):
    n=len(y); sse=np.sum((y-p)**2); return n*np.log(sse/n)+2*npar
def sat(x,dinf,a,lc): return dinf-a*np.exp(-x/lc)
def main():
    d=pd.read_csv(IN); q=d[d.source.str.startswith('10L')].sort_values('wavelength_nm').copy(); x=q.wavelength_nm.to_numpy(); y=q.depth_mean.to_numpy()
    if not np.all(x>0): raise ValueError('wavelengths must be strictly positive for the logarithmic screen')
    X=np.c_[np.ones(len(x)),np.log(np.maximum(x,1e-300))]; alog,blog=np.linalg.lstsq(X,y,rcond=None)[0]; q['log_prediction']=alog+blog*np.log(x)
    pars,_=curve_fit(sat,x,y,p0=(.66,.45,20),bounds=([.65,0,.01],[2,2,1000]),maxfev=100000); q['saturation_prediction']=sat(x,*pars); q['log_residual']=y-q.log_prediction; q['saturation_residual']=y-q.saturation_prediction
    q['log_R2']=r2(y,q.log_prediction); q['saturation_R2']=r2(y,q.saturation_prediction); q['log_AIC']=aic(y,q.log_prediction,2); q['saturation_AIC']=aic(y,q.saturation_prediction,3); q['log_intercept']=alog; q['log_slope']=blog; q['D_infinity']=pars[0]; q['amplitude']=pars[1]; q['lambda_c_nm']=pars[2]; q.to_csv(CSV,index=False,float_format='%.9g')
    xx=np.linspace(x.min(),x.max(),400)
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'xtick.direction':'out','ytick.direction':'out'})
    fig,axs=plt.subplots(1,2,figsize=(7.25,2.75),constrained_layout=True); blue='#2775a9'; red='#c94c4c'; gray='#555555'
    ax=axs[0]; ax.errorbar(x,y,q.depth_sem,color=blue,marker='o',ms=5,lw=0,capsize=2,label=r'10L modes $n=1\ldots10$'); ax.plot(xx,alog+blog*np.log(xx),color=gray,lw=1.3,ls='--',label=rf'log: $D={alog:.3f}+{blog:.3f}\ln\lambda$'); ax.plot(xx,sat(xx,*pars),color=red,lw=1.5,label=rf'saturation: $D_\infty={pars[0]:.3f}$, $\lambda_c={pars[2]:.1f}$ nm'); ax.set(xlabel=r'wavelength, $\lambda=2\pi/k$ (nm)',ylabel=r'relative valley depth, $D_-$',xlim=(0,107)); ax.spines[['top','right']].set_visible(False); ax.legend(frameon=False,fontsize=6.0,loc='lower right'); ax.text(-.12,1.04,'(a)',transform=ax.transAxes,fontweight='bold',fontsize=10)
    ax=axs[1]; ax.axhline(0,color='.45',lw=.8); ax.plot(x,q.log_residual,color=gray,marker='o',ms=4,lw=1.1,label=rf'log, $R^2={r2(y,q.log_prediction):.3f}$'); ax.plot(x,q.saturation_residual,color=red,marker='s',ms=4,lw=1.1,label=rf'saturation, $R^2={r2(y,q.saturation_prediction):.3f}$'); ax.set(xlabel=r'wavelength, $\lambda=2\pi/k$ (nm)',ylabel='fit residual in $D_-$',xlim=(0,107)); ax.spines[['top','right']].set_visible(False); ax.legend(frameon=False,fontsize=6.5,loc='upper left'); ax.text(-.12,1.04,'(b)',transform=ax.transAxes,fontweight='bold',fontsize=10)
    fig.suptitle(r'Can the 10L relative valley depth be logarithmic in wavelength?',y=1.03,fontsize=10)
    fig.text(.5,-.075,r'Empirical screen only: the log model is acceptable over the sampled range but has systematic residuals; a three-parameter saturation is preferred by unweighted $R^2$/AIC.  This does not establish a unique asymptotic law.',ha='center',fontsize=6.0,color='.25')
    fig.savefig(STEM.with_suffix('.png'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); fig.savefig(STEM.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(STEM.with_suffix('.svg'),bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
