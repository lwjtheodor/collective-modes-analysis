"""Test whether the alpha_z minimum deepens monotonically with box length."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

DATA = Path(r"C:\Users\s1365\.codex\visualizations\2026\08\06\019fd4c1-37fd-7b92-b606-edcea5bf0c15\stage_vacf_integral_lockin_20260806\output")
ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
OUT = ROOT / "assets" / "alphaz_minima_vs_length_nature"
COLORS = ["#203864", "#2F75B5", "#4FA3A5", "#5B9A6D", "#A07A39"]

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"], "font.size":7,
    "axes.labelsize":7.4, "axes.titlesize":8, "xtick.labelsize":6.5, "ytick.labelsize":6.5,
    "axes.linewidth":.75, "axes.spines.top":False, "axes.spines.right":False,
    "svg.fonttype":"none", "pdf.fonttype":42,
})

def values(L):
    vals, times = [], []
    files = sorted(p for p in DATA.glob(f"*L{L}_rep*.csv") if not p.name.endswith("_lockin_blocks.csv"))
    for p in files:
        d = pd.read_csv(p)
        d = d[(d.lag_ps >= 5) & (d.lag_ps <= 100)]
        row = d.loc[d.alpha_from_D_MSD_mean.idxmin()]
        vals.append(row.alpha_from_D_MSD_mean); times.append(row.lag_ps)
    return np.asarray(vals), np.asarray(times)

def mean_sem(x):
    return x.mean(), x.std(ddof=1)/np.sqrt(len(x))

def main():
    mins, tmins = zip(*(values(L) for L in range(1,6)))
    m1, s1 = np.array([mean_sem(x) for x in mins]).T
    m2, s2 = np.array([mean_sem(x) for x in tmins]).T
    L = np.arange(1,6)
    fig, (a,b) = plt.subplots(1,2,figsize=(7.15,2.65))
    for i, x in enumerate(mins):
        a.scatter(np.full(3,L[i]), x, s=18, color=COLORS[i], alpha=.68, zorder=3)
    a.errorbar(L,m1,yerr=s1,color="#263238",lw=.9,marker="o",ms=4,capsize=2,zorder=2)
    a.set(xticks=L,xlim=(.6,5.4),ylim=(.2,.82),xlabel="box length",ylabel=r"minimum $\alpha_z(t)$")
    a.set_title(r"minimum depth, 5-100 ps",loc="left",fontweight="bold")
    a.text(.02,.96,"a",transform=a.transAxes,va="top",fontsize=9,fontweight="bold")
    a.annotate("deepening only through 3L",xy=(3,m1[2]),xytext=(3.4,.69),fontsize=6.2,
               arrowprops=dict(arrowstyle="-",lw=.6,color="#555555"),color="#444444")
    for i,x in enumerate(tmins):
        b.scatter(np.full(3,L[i]),x,s=18,color=COLORS[i],alpha=.68,zorder=3)
    b.errorbar(L,m2,yerr=s2,color="#263238",lw=.9,marker="o",ms=4,capsize=2,zorder=2)
    b.set(xticks=L,xlim=(.6,5.4),ylim=(0,30),xlabel="box length",ylabel=r"time of minimum, $t_{\min}$ (ps)")
    b.set_title(r"minimum shifts later with box length",loc="left",fontweight="bold")
    b.text(.02,.96,"b",transform=b.transAxes,va="top",fontsize=9,fontweight="bold")
    fig.text(.5,.01,r"(8,8), weak Nose-Hoover, no momentum removal; O 10 fs; 1 ns; n=3 trajectories per L. Points: replicas; error bars: replica SEM.",ha="center",fontsize=6.0)
    fig.subplots_adjust(left=.095,right=.985,bottom=.25,top=.86,wspace=.33)
    fig.savefig(str(OUT)+".png",dpi=600,bbox_inches="tight")
    fig.savefig(str(OUT)+".tiff",dpi=600,bbox_inches="tight")
    fig.savefig(str(OUT)+".pdf",bbox_inches="tight")
    fig.savefig(str(OUT)+".svg",bbox_inches="tight")

if __name__ == "__main__":
    main()
