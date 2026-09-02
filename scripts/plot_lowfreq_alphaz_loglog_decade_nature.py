"""Aggregate low-frequency direct-MSD one-decade alpha_z results and render figure."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

DATA=Path(r"C:\Users\s1365\.codex\visualizations\2026\08\06\019fd4c1-37fd-7b92-b606-edcea5bf0c15\stage_lowfreq_alphaz_loglog_8_8_20260806\output")
ROOT=Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
CURVES=ROOT/"assets"/"lowfreq_alphaz_loglog_1decade_curves.csv"
SUMMARY=ROOT/"assets"/"lowfreq_alphaz_loglog_1decade_minima.csv"
OUT=ROOT/"assets"/"lowfreq_alphaz_loglog_1decade_nature"
COLORS=["#203864","#2F75B5","#4FA3A5","#5B9A6D","#A07A39"]
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.labelsize":7.3,"axes.titlesize":8,"xtick.labelsize":6.3,"ytick.labelsize":6.3,"legend.fontsize":5.8,"axes.linewidth":.75,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})

def files(L):
    return [DATA/(f"{L}L_rep1.csv"),DATA/(f"{L}L_seed2.csv"),DATA/(f"{L}L_seed3.csv")]

def main():
    curves=[]; mins=[]; fig,(a,b)=plt.subplots(1,2,figsize=(7.15,2.75),gridspec_kw={"width_ratios":[1.42,1]})
    for L in range(1,6):
        fs=[pd.read_csv(p) for p in files(L)]; t=fs[0].time_ps.to_numpy()
        y=np.vstack([f.alpha_z_loglog_1decade_mean.to_numpy() for f in fs])
        mean=y.mean(axis=0); sem=y.std(axis=0,ddof=1)/np.sqrt(3)
        for i,tt in enumerate(t): curves.append({"box_length":L,"time_ps":tt,"alpha_z_mean":mean[i],"alpha_z_replica_sem":sem[i],"alpha_z_rep1":y[0,i],"alpha_z_rep2":y[1,i],"alpha_z_rep3":y[2,i]})
        a.fill_between(t,mean-sem,mean+sem,color=COLORS[L-1],alpha=.16,lw=0); a.plot(t,mean,color=COLORS[L-1],lw=1.1,label=f"{L}L")
        ii=np.argmin(y,axis=1); v=y[np.arange(3),ii]; tm=t[ii]
        rec={"box_length":L,"alpha_min_mean":v.mean(),"alpha_min_replica_sem":v.std(ddof=1)/np.sqrt(3),"t_min_ps_mean":tm.mean(),"t_min_ps_replica_sem":tm.std(ddof=1)/np.sqrt(3)}
        for r in range(3): rec[f"alpha_min_rep{r+1}"]=v[r]; rec[f"t_min_rep{r+1}_ps"]=tm[r]
        mins.append(rec); b.scatter(np.full(3,L),v,color=COLORS[L-1],s=18,alpha=.68,zorder=3); b.errorbar(L,v.mean(),yerr=v.std(ddof=1)/np.sqrt(3),color="#263238",marker="o",ms=4,lw=.85,capsize=2,zorder=2)
    pd.DataFrame(curves).to_csv(CURVES,index=False); pd.DataFrame(mins).to_csv(SUMMARY,index=False)
    a.axhline(1,color="#555555",lw=.65,zorder=0); a.set(xscale="log",xlim=(3.16,31.62),ylim=(.35,1.15),xlabel=r"centre time, $t$ (ps)",ylabel=r"$\alpha_z=d\log M_z/d\log t$"); a.set_xticks([3,5,10,20,30]); a.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter()); a.set_title("low-frequency direct-MSD slope",loc="left",fontweight="bold"); a.legend(title="box length",ncol=2,loc="lower right",handlelength=1.4,labelspacing=.25,borderpad=.2); a.text(.02,.96,"a",transform=a.transAxes,va="top",fontsize=9,fontweight="bold")
    b.set(xticks=np.arange(1,6),xlim=(.6,5.4),ylim=(.25,.85),xlabel="box length",ylabel=r"minimum one-decade $\alpha_z$"); b.set_title("minimum within full-window support",loc="left",fontweight="bold"); b.text(.02,.96,"b",transform=b.transAxes,va="top",fontsize=9,fontweight="bold")
    fig.text(.5,.01,r"(8,8), baseline NVT; O 1 ps; 20 ns; z-momentum removed every 5 ps; n=3 trajectories per L. Direct CNT-relative MSD; uniform log-time grid; OLS over fixed 1-decade window; bands/error bars: replica SEM.",ha="center",fontsize=5.75)
    fig.subplots_adjust(left=.09,right=.985,bottom=.24,top=.86,wspace=.32)
    fig.savefig(str(OUT)+".png",dpi=600,bbox_inches="tight"); fig.savefig(str(OUT)+".tiff",dpi=600,bbox_inches="tight"); fig.savefig(str(OUT)+".pdf",bbox_inches="tight"); fig.savefig(str(OUT)+".svg",bbox_inches="tight")
if __name__=="__main__": main()
