#!/usr/bin/env python3
"""Fit finite-frequency radial-current resonances with a damped-harmonic oscillator."""
from __future__ import annotations
import argparse, csv, json, math
from datetime import date
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t

ROOT=Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
PSD=ROOT/"results"/"collective_mode_response"/"88_10L_per_k_semilog_Skw_pm10_LA_TAr_TAtheta"/"2026-08-19"/"derived_data"/"signed_Skw_n001_n160.csv"
PEAK=ROOT/"results"/"collective_mode_response"/"88_10L_per_k_semilog_Skw_LA_TAr_TAtheta"/"2026-08-19"/"derived_data"/"phase_and_group_velocity_vs_k.csv"
OUT=ROOT/"results"/"collective_mode_response"/"88_10L_TAr_DHO_linewidth_n001_n010"/"2026-08-24"

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.spines.right":False,"axes.spines.top":False,"xtick.direction":"out","ytick.direction":"out","svg.fonttype":"none","pdf.fonttype":42,"legend.frameon":False})

def dho(w,b,a,w0,g): return b+a*g/((w*w-w0*w0)**2+(2*g*w)**2)
def load(psd, peak):
    allp={}; seeds={}
    with open(peak,newline="") as f:
        for r in csv.DictReader(f):
            if r["branch"]=="TA_r" and int(r["n"])<=10: seeds[int(r["n"])]=float(r.get("omega_peak_rad_ps",r["omega_peak_mean_rad_ps"]))
    with open(psd,newline="") as f:
        for r in csv.DictReader(f):
            if r["branch"]=="TA_r" and int(r["n"])<=10:
                w=2*np.pi*float(r["frequency_ps_inv"])
                if 12<=w<=55: allp.setdefault(int(r["n"]),[]).append((w,float(r["PSD_mean_arbitrary"]),float(r["PSD_replica_SEM_arbitrary"]),float(r["k_inv_A"])))
    return allp,seeds
def fit(n,rec,seed):
    a=np.asarray(rec); w,y,e=a[:,0],a[:,1],a[:,2]; dw=float(np.median(np.diff(w))); half=12.
    keep=(w>=seed-half)&(w<=seed+half); w,y,e=w[keep],y[keep],e[keep]
    base=max(y.min(),1e-14); sigma=np.maximum(e,max(y.max()*.01,1e-14))
    p,c=curve_fit(dho,w,y,p0=[base,max((y.max()-base)*seed**3,1e-12),seed,4.],sigma=sigma,bounds=([0,0,seed-5,dw/2],[np.inf,np.inf,seed+5,15]),maxfev=200000)
    pred=dho(w,*p); r2=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2); se=np.sqrt(c[3,3]) if c[3,3]>0 else np.nan
    return dict(n=n,k_inv_A=float(a[0,3]),frequency_bin_rad_ps=dw,omega0_DHO_rad_ps=float(p[2]),gamma_DHO_rad_ps=float(p[3]),gamma_fit_SEM_rad_ps=float(se),fit_R2=float(r2),accepted=bool(r2>=.7 and p[3]>=1.5*dw),omega=w,signal=y,sem=e,prediction=pred)
def k2(rows,label):
    x=np.array([r["k_inv_A"]**2 for r in rows]); y=np.array([r["gamma_DHO_rad_ps"] for r in rows]); X=np.c_[np.ones(len(x)),x]; p=np.linalg.lstsq(X,y,rcond=None)[0]; res=y-X@p; cov=(res@res/(len(x)-2))*np.linalg.inv(X.T@X); ci=t.ppf(.975,len(x)-2)*np.sqrt(np.diag(cov)); return dict(fit_label=label,n_min=min(r["n"] for r in rows),n_max=max(r["n"] for r in rows),n_points=len(rows),Gamma_r0_rad_ps=float(p[0]),Gamma_r0_95CI_halfwidth_rad_ps=float(ci[0]),D_r_k2_rad_ps_A2=float(p[1]),D_r_95CI_halfwidth_rad_ps_A2=float(ci[1]),R2_linear=float(1-res@res/np.sum((y-y.mean())**2)),RSS_linear=float(res@res))
def save(f,fig,stem):
    f.savefig(fig/(stem+".png"),dpi=600,bbox_inches="tight");f.savefig(fig/(stem+".pdf"),bbox_inches="tight");f.savefig(fig/(stem+".svg"),bbox_inches="tight");f.savefig(fig/(stem+".tiff"),dpi=600,bbox_inches="tight")
def plot(rows,fit,fig):
    f=plt.figure(figsize=(7.1,3.)); axes=f.subplots(1,2)
    ax=axes[0]; use=[r for r in rows if r["accepted"]]; x=np.array([r["k_inv_A"]**2 for r in use]);y=np.array([r["gamma_DHO_rad_ps"] for r in use]);e=np.array([max(r["gamma_fit_SEM_rad_ps"],.1*r["gamma_DHO_rad_ps"]) for r in use]);ax.errorbar(x,y,yerr=e,fmt="o",color="#0072b2",capsize=2,label="DHO damping");xx=np.linspace(0,x.max()*1.05,200);ax.plot(xx,fit["Gamma_r0_rad_ps"]+fit["D_r_k2_rad_ps_A2"]*xx,color="#d55e00",label=r"$\Gamma_{r,0}+D_rk^2$");ax.axhline(fit["Gamma_r0_rad_ps"],color=".3",ls="--",lw=.8);ax.set(xlabel=r"$k^2$ ($\mathrm{\AA}^{-2}$)",ylabel=r"TA$_r$ DHO damping, $\Gamma_r$ (rad ps$^{-1}$)");ax.legend(fontsize=6);ax.text(.03,.96,rf"$\Gamma_{{r,0}}={fit['Gamma_r0_rad_ps']:.2f}\pm{fit['Gamma_r0_95CI_halfwidth_rad_ps']:.2f}$",transform=ax.transAxes,va="top");
    ax=axes[1]
    for r in (rows[0],rows[4],rows[9]): ax.plot(r["omega"],r["signal"],lw=.7,label=rf"$n={r['n']}$");ax.plot(r["omega"],r["prediction"],color="#d55e00",lw=.8)
    ax.set(xlabel=r"$\omega$ (rad ps$^{-1}$)",ylabel=r"$S_{J_rJ_r}$ (arb.)");ax.legend(fontsize=6); f.subplots_adjust(.11,.19,.985,.92,.34);save(f,fig,"TA_r_DHO_linewidth_and_k2_fit");plt.close(f)
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--psd-csv",type=Path,required=True);ap.add_argument("--peak-csv",type=Path,required=True);ap.add_argument("--outdir",type=Path,required=True);a=ap.parse_args();fig=a.outdir/"figures";data=a.outdir/"derived_data";fig.mkdir(parents=True,exist_ok=True);data.mkdir(exist_ok=True);p,s=load(a.psd_csv,a.peak_csv);rows=[fit(n,p[n],s[n]) for n in range(1,11)];fields=[q for q in rows[0] if q not in {"omega","signal","sem","prediction"}];
    with open(data/"TA_r_DHO_resonance_fits_n001_n010.csv","w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows([{q:r[q] for q in fields} for r in rows])
    good=[r for r in rows if r["accepted"]]; model=k2(good,"resolved_n001_n010");
    with open(data/"TA_r_DHO_friction_intercept_plus_k2.csv","w",newline="") as f:w=csv.DictWriter(f,fieldnames=list(model));w.writeheader();w.writerow(model)
    plot(rows,model,fig);meta={"analysis_date":str(date.today()),"branch":"TA_r","model":"DHO finite-frequency resonance; Gamma_r(k)=Gamma_r0+D_r k^2","input_psd":str(a.psd_csv),"input_peaks":str(a.peak_csv),"primary_fit":model,"limit":"PSD-frequency SEM is used for spectral weighting; no per-replica resonance fits are archived, so parameter CI is conditional mode scatter."};(a.outdir/"analysis_manifest.json").write_text(json.dumps(meta,indent=2)+"\n");(a.outdir/"FINISHED.txt").write_text("10L TA_r DHO linewidth analysis finished successfully.\n");print(json.dumps(model,indent=2))
if __name__=="__main__":main()
