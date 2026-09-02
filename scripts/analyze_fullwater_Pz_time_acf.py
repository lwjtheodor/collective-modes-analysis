"""Time-domain ACF of mass-weighted full-water axial momentum."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

COL = "Pz_total_water_amu_A_fs"

def acf(x):
    x=x-x.mean(); n=len(x)
    y=np.fft.irfft(np.abs(np.fft.rfft(x,2*n))**2,2*n)[:n]
    y/=np.arange(n,0,-1)
    return y/y[0]

def exp1(t,a,g,b): return b+a*np.exp(-g*t)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--fit-max-ps',type=float,default=100.0); a=p.parse_args()
    d=a.output/'derived_data'; f=a.output/'figures'; d.mkdir(parents=True,exist_ok=True); f.mkdir(exist_ok=True)
    raw=pd.read_csv(a.input); rows=[]; series=[]
    for i,(case,q) in enumerate(raw.groupby('case'),1):
        q=q.sort_values('time_ps'); t=q.time_ps.to_numpy()-q.time_ps.iloc[0]; c=acf(q[COL].to_numpy())
        m=(t<=a.fit_max_ps)&(c>-0.15)
        par,_=curve_fit(exp1,t[m],c[m],p0=(1,1.0,0),bounds=([0,1e-6,-.2],[1.5,100,.2]),maxfev=50000)
        pred=exp1(t[m],*par); r2=1-np.sum((c[m]-pred)**2)/np.sum((c[m]-c[m].mean())**2)
        rows.append({'case':case,'n_frames':len(t),'dt_ps':np.median(np.diff(t)),'Gamma_exp_rad_ps':par[1],'tau_ps':1/par[1],'fit_R2':r2})
        series.append(pd.DataFrame({'case':case,'time_ps':t,'C_PzPz':c}))
    out=pd.concat(series); grid=np.sort(out.time_ps.unique()); mat=np.vstack([out[out.case==r['case']].C_PzPz.to_numpy() for r in rows]); mean=mat.mean(0); sem=mat.std(0,ddof=1)/np.sqrt(len(mat))
    m=(grid<=a.fit_max_ps)&(mean>-0.15); par,_=curve_fit(exp1,grid[m],mean[m],p0=(1,1.0,0),bounds=([0,1e-6,-.2],[1.5,100,.2]),maxfev=50000); pred=exp1(grid[m],*par); r2=1-np.sum((mean[m]-pred)**2)/np.sum((mean[m]-mean[m].mean())**2)
    pd.DataFrame(rows).to_csv(d/'fullwater_Pz_ACF_per_replica_fit.csv',index=False); out.to_csv(d/'fullwater_Pz_ACF_per_replica.csv',index=False)
    pd.DataFrame({'time_ps':grid,'C_PzPz_mean':mean,'C_PzPz_replica_sem':sem}).to_csv(d/'fullwater_Pz_ACF_ensemble.csv',index=False)
    pd.DataFrame([{'definition':'Pz=sum all water atoms m_a*v_za','Gamma_exp_rad_ps':par[1],'tau_ps':1/par[1],'fit_R2':r2,'n_replicas':len(mat),'dt_ps':np.median(np.diff(grid)),'fit_max_ps':a.fit_max_ps}]).to_csv(d/'fullwater_Pz_ACF_ensemble_summary.csv',index=False)
    fig,ax=plt.subplots(figsize=(6.4,4.3)); ax.plot(grid[m],mean[m],color='#2166ac',label='4-rep mean'); ax.fill_between(grid[m],mean[m]-sem[m],mean[m]+sem[m],color='#2166ac',alpha=.2); ax.plot(grid[m],pred,color='#b2182b',label=fr'exp: $\Gamma={par[1]:.4g}$ ps$^{{-1}}$'); ax.axhline(0,color='black',lw=.7); ax.set(xlabel='lag time (ps)',ylabel=r'$C_{P_zP_z}(t)/C(0)$'); ax.grid(alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(f/'fullwater_Pz_ACF.png',dpi=300); fig.savefig(f/'fullwater_Pz_ACF.pdf'); plt.close(fig)
    (a.output/'metadata_time_acf.json').write_text(json.dumps({'input':str(a.input),'source_sampling':'1 ps decimation of 100 fs full-water dump','protocol_limit':'fixed CNT weak-NH NVT: this is a bath-plus-wall relaxation rate, not intrinsic NVE hydrodynamic damping'},indent=2))
    (a.output/'TIME_ACF_FINISHED.txt').write_text('Full-water Pz time-domain ACF finished successfully.\n')
if __name__=='__main__': main()
