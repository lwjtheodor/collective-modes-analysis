"""Estimate alpha_z from direct MSD using one-decade, uniform-log local slopes."""
from pathlib import Path
import csv
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

DATA = Path(r"C:\Users\s1365\.codex\visualizations\2026\08\06\019fd4c1-37fd-7b92-b606-edcea5bf0c15\stage_vacf_integral_lockin_20260806\output")
ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
OUT = ROOT / "assets" / "alphaz_loglog_uniform_decade_nature"
TABLE = ROOT / "assets" / "alphaz_loglog_uniform_decade_summary.csv"
COLORS = ["#203864", "#2F75B5", "#4FA3A5", "#5B9A6D", "#A07A39"]
TGRID = np.logspace(np.log10(.5), np.log10(100), 241)

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"], "font.size":7,
    "axes.labelsize":7.3, "axes.titlesize":8, "xtick.labelsize":6.3, "ytick.labelsize":6.3,
    "legend.fontsize":5.8, "axes.linewidth":.75, "axes.spines.top":False, "axes.spines.right":False,
    "svg.fonttype":"none", "pdf.fonttype":42,
})

def decade_slope(t, msd):
    use = (t >= .5) & (t <= 100) & (msd > 0)
    lx = np.log10(TGRID)
    msd_safe = np.maximum(msd[use], 1e-12)  # explicit positive guard for log(M)
    ly = np.interp(lx, np.log10(t[use]), np.log10(msd_safe))
    out = np.full_like(TGRID, np.nan)
    half = .5
    for i, x0 in enumerate(lx):
        mask = np.abs(lx-x0) <= half + 1e-12
        if mask.sum() >= 3 and x0-half >= lx[0] and x0+half <= lx[-1]:
            out[i] = np.polyfit(lx[mask], ly[mask], 1)[0]
    return out

def per_length(L):
    traces=[]
    files = sorted(p for p in DATA.glob(f"*L{L}_rep*.csv") if not p.name.endswith("_lockin_blocks.csv"))
    for p in files:
        d = pd.read_csv(p)
        m = d.MSD_direct_A2_mean.to_numpy() - d.MSD_direct_A2_mean.iloc[0]
        traces.append(decade_slope(d.lag_ps.to_numpy(), m))
    return np.asarray(traces)

def main():
    alltr = {L:per_length(L) for L in range(1,6)}
    valid = np.isfinite(alltr[1][0])
    t = TGRID[valid]
    fig,(a,b)=plt.subplots(1,2,figsize=(7.15,2.75),gridspec_kw={"width_ratios":[1.42,1]})
    rows=[]
    for L in range(1,6):
        x=alltr[L][:,valid]; mean=x.mean(axis=0); sem=x.std(axis=0,ddof=1)/np.sqrt(3)
        a.fill_between(t,mean-sem,mean+sem,color=COLORS[L-1],alpha=.16,lw=0)
        a.plot(t,mean,color=COLORS[L-1],lw=1.1,label=f"{L}L")
        im=np.argmin(x,axis=1); minima=x[np.arange(3),im]; times=t[im]
        mm=minima.mean(); ss=minima.std(ddof=1)/np.sqrt(3)
        tm=times.mean(); ts=times.std(ddof=1)/np.sqrt(3)
        rows.append((L,mm,ss,tm,ts,*minima,*times))
        b.scatter(np.full(3,L),minima,color=COLORS[L-1],s=18,alpha=.68,zorder=3)
        b.errorbar(L,mm,yerr=ss,color="#263238",marker="o",ms=4,lw=.85,capsize=2,zorder=2)
    a.axhline(1,color="#555555",lw=.65,zorder=0)
    a.set(xscale="log",xlim=(t.min(),t.max()),ylim=(.35,1.12),xlabel=r"centre time, $t$ (ps)",ylabel=r"$alpha_z=dlog M/dlog t$")
    a.set_xticks([2,3,5,10,20,30]); a.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    a.set_title("uniform-log, one-decade MSD slope",loc="left",fontweight="bold")
    a.legend(title="box length",ncol=2,loc="lower right",handlelength=1.4,labelspacing=.25,borderpad=.2)
    a.text(.02,.96,"a",transform=a.transAxes,va="top",fontsize=9,fontweight="bold")
    b.set(xticks=np.arange(1,6),xlim=(.6,5.4),ylim=(.25,.85),xlabel="box length",ylabel=r"minimum one-decade $alpha_z$")
    b.set_title("minimum within full-window support",loc="left",fontweight="bold")
    b.text(.02,.96,"b",transform=b.transAxes,va="top",fontsize=9,fontweight="bold")
    fig.text(.5,.01,r"(8,8), weak Nose-Hoover, no momentum removal; O 10 fs; 1 ns; n=3 trajectories per L. Uniform log-time grid; OLS over a fixed 1-decade window; bands/error bars: replica SEM.",ha="center",fontsize=5.8)
    fig.subplots_adjust(left=.09,right=.985,bottom=.24,top=.86,wspace=.32)
    with TABLE.open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["box_length","alpha_min_mean","alpha_min_replica_sem","t_min_ps_mean","t_min_ps_replica_sem","alpha_min_rep1","alpha_min_rep2","alpha_min_rep3","t_min_rep1_ps","t_min_rep2_ps","t_min_rep3_ps"]); w.writerows(rows)
    fig.savefig(str(OUT)+".png",dpi=600,bbox_inches="tight")
    fig.savefig(str(OUT)+".tiff",dpi=600,bbox_inches="tight")
    fig.savefig(str(OUT)+".pdf",bbox_inches="tight")
    fig.savefig(str(OUT)+".svg",bbox_inches="tight")

if __name__ == "__main__":
    main()
