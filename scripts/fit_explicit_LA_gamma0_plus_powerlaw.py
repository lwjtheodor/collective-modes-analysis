"""Fit resolved explicit-CNT LA widths to Gamma0 + A*k**alpha."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

ROOT=Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
INP=ROOT/'results/collective_mode_response/88_10L_LA_linewidth_powerlaw/2026-08-24/derived_data/LA_peak_linewidth_lorentzian_fits_n001_n024.csv'
OUT=ROOT/'results/collective_mode_response/88_10L_LA_linewidth_gamma0_plus_powerlaw/2026-08-24'
G0=0.014788264159646673; G0SEM=0.0009874895786320522

def model(k,A,alpha,g0=G0): return g0+A*k**alpha
def fit(q,g0=G0):
    k=q.k_inv_A.to_numpy(); y=q.gamma_HWHM_rad_ps.to_numpy(); s=q.gamma_fit_SEM_rad_ps.to_numpy()
    p,c=curve_fit(lambda x,A,a:model(x,A,a,g0),k,y,p0=(25,1.5),sigma=s,absolute_sigma=True,bounds=([0,.1],[np.inf,4]),maxfev=50000)
    pred=model(k,*p,g0); chi=np.sum(((y-pred)/s)**2); return p,c,pred,chi
def main():
    d=OUT/'derived_data'; f=OUT/'figures'; d.mkdir(parents=True,exist_ok=True); f.mkdir(exist_ok=True)
    raw=pd.read_csv(INP); raw=raw[(raw.accepted_linewidth)&(raw.n>=3)&(raw.n<=15)].copy(); records=[]
    rng=np.random.default_rng(20260824)
    for nmax in (10,15):
        q=raw[raw.n<=nmax]; p,c,pred,chi=fit(q); boot=[]
        for _ in range(2000):
            g0=rng.normal(G0,G0SEM); yy=rng.normal(q.gamma_HWHM_rad_ps,q.gamma_fit_SEM_rad_ps)
            try: boot.append(curve_fit(lambda x,A,a:model(x,A,a,g0),q.k_inv_A,yy,p0=p,sigma=q.gamma_fit_SEM_rad_ps,absolute_sigma=True,bounds=([0,.1],[np.inf,4]),maxfev=30000)[0])
            except RuntimeError: pass
        b=np.asarray(boot)
        records.append({'fit_label':f'Gamma0_fixed_n003_n{nmax:03d}','n_min':3,'n_max':nmax,'n_points':len(q),'Gamma0_rad_ps':G0,'Gamma0_sem_rad_ps':G0SEM,'A_rad_ps_A_to_alpha':p[0],'A_bootstrap_95CI_low':np.quantile(b[:,0],.025),'A_bootstrap_95CI_high':np.quantile(b[:,0],.975),'alpha':p[1],'alpha_bootstrap_95CI_low':np.quantile(b[:,1],.025),'alpha_bootstrap_95CI_high':np.quantile(b[:,1],.975),'weighted_chi2':chi,'dof':len(q)-2})
        q=q.copy(); q['Gamma_model_rad_ps']=pred; q['fit_label']=records[-1]['fit_label']; q.to_csv(d/f'LA_widths_{records[-1]["fit_label"]}.csv',index=False)
    pd.DataFrame(records).to_csv(d/'LA_gamma0_plus_powerlaw_fits.csv',index=False)
    q=raw[raw.n<=10]; r=records[0]; kk=np.linspace(0,q.k_inv_A.max()*1.08,300); fig,ax=plt.subplots(figsize=(6.3,4.3)); ax.errorbar(q.k_inv_A,q.gamma_HWHM_rad_ps,yerr=q.gamma_fit_SEM_rad_ps,fmt='o',color='#2166ac',capsize=2,label='resolved LA HWHM, n=3–10'); ax.plot(kk,model(kk,r['A_rad_ps_A_to_alpha'],r['alpha']),color='#b2182b',label=fr'$\Gamma_0+A k^\alpha$, $\alpha={r["alpha"]:.2f}$'); ax.axhline(G0,color='.3',ls='--',lw=.8,label=fr'$\Gamma_0={G0:.4f}$'); ax.set(xlabel=r'$k$ ($\mathrm{\AA}^{-1}$)',ylabel=r'LA HWHM $\Gamma$ (rad ps$^{-1}$)'); ax.grid(alpha=.2); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(f/'LA_gamma0_plus_powerlaw_n003_n010.png',dpi=300); fig.savefig(f/'LA_gamma0_plus_powerlaw_n003_n010.pdf'); plt.close(fig)
    (OUT/'metadata.json').write_text(json.dumps({'model':'Gamma(k)=Gamma0+A*k^alpha','Gamma0_source':'full-water Pz time ACF, explicit 5L 10ns x4, weak-NH Tdamp=100ps','LA_source':str(INP),'selection':'n=3..10 primary; n=3..15 sensitivity; n1-2 unresolved','limit':'Gamma0 is transferred from 5L total momentum to 10L oxygen-current LA widths; compare only as a protocol-matched operational model'},indent=2))
    (OUT/'FINISHED.txt').write_text('Explicit LA Gamma0-plus-power-law fit finished successfully.\n')
if __name__=='__main__':main()
