"""Test physically motivated normalizations of the n=1 C_vJ negative-lobe area."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT=Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
INP=ROOT/"assets"/"lowfreq_alphaz_cvj_n1_first_lobe_paired.csv"
CSV=ROOT/"assets"/"lowfreq_cvj_n1_lobe_normalization_summary.csv"
OUT=ROOT/"assets"/"lowfreq_cvj_n1_lobe_normalization_nature"
COLORS=["#203864","#2F75B5","#4FA3A5","#5B9A6D","#A07A39"]
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.labelsize":7.2,"axes.titlesize":7.8,"xtick.labelsize":6.2,"ytick.labelsize":6.2,"legend.fontsize":5.8,"axes.linewidth":.75,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})

def fit(x,y):
    assert np.all(x>0) and np.all(y>0)
    r=linregress(np.log(np.maximum(x,1e-15)),np.log(np.maximum(y,1e-15))); return r.slope,r.intercept,r.rvalue*r.rvalue
def main():
    d=pd.read_csv(INP); d["depth"]= -d.CvJ_n1_first_minimum_value; d["width_ps"]=d.CvJ_n1_lobe_width_ps; d["area_ps"]=d.CvJ_n1_first_negative_lobe_area_ps
    d["shape_factor"]=d.area_ps/(d.depth*d.width_ps); d["area_per_width"]=d.area_ps/d.width_ps
    m=d.groupby("box_length")[["area_ps","depth","width_ps","shape_factor","area_per_width"]].agg(["mean","sem"]); m.columns=["_".join(q) for q in m.columns]; m=m.reset_index(); m["k_inv_A"]=0.06230846322685951/m.box_length
    m.to_csv(CSV,index=False)
    L=m.box_length.to_numpy(float); ea,ba,r2a=fit(L,m.area_ps_mean); ew,bw,r2w=fit(L,m.width_ps_mean); ed,bd,r2d=fit(L,m.depth_mean); es,bs,r2s=fit(L,m.shape_factor_mean)
    fig,axs=plt.subplots(1,3,figsize=(7.15,2.55)); a,b,c=axs
    for y,err,label,color,p,r2 in [(m.area_ps_mean,m.area_ps_sem,r"$A_-$", "#263238",ea,r2a),(m.width_ps_mean,m.width_ps_sem,r"$\Delta t_-$","#5B9A6D",ew,r2w),(m.depth_mean,m.depth_sem,r"$|C_{\min}|$","#A07A39",ed,r2d)]:
        a.errorbar(L,y,yerr=err,marker="o",ms=3.8,capsize=2,lw=.8,color=color,label=label)
        xx=np.linspace(1,5,100); a.plot(xx,np.exp(np.log(y[0])-p*np.log(L[0]))*xx**p,color=color,lw=.65,ls="--")
    a.set(xscale="log",yscale="log",xlim=(.9,5.5),xlabel="box length, $L$",ylabel="negative-lobe components")
    a.set_xticks(L); a.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter()); a.set_title(r"$A_-\approx \Delta t_-|C_{\min}|S$",loc="left",fontweight="bold"); a.legend(ncol=1,loc="upper left",handlelength=1.2,labelspacing=.25); a.text(.02,.04,fr"exponents: {ea:.2f}, {ew:.2f}, {ed:.2f}",transform=a.transAxes,fontsize=5.6)
    b.errorbar(L,m.area_per_width_mean,yerr=m.area_per_width_sem,marker="o",ms=4,capsize=2,color="#2F75B5",lw=.8); b.errorbar(L,m.depth_mean,yerr=m.depth_sem,marker="s",ms=3.6,capsize=2,color="#A07A39",lw=.8)
    b.set(xlim=(.7,5.3),xticks=L,xlabel="box length, $L$",ylabel=r"$A_-/\Delta t_-$ or $|C_{\min}|$"); b.set_title(r"time normalization alone does not collapse",loc="left",fontweight="bold"); b.legend([r"$A_-/\Delta t_-$",r"$|C_{\min}|$"],loc="upper left",fontsize=5.5)
    c.errorbar(L,m.shape_factor_mean,yerr=m.shape_factor_sem,marker="o",ms=4,capsize=2,color="#263238",lw=.8); c.axhline(m.shape_factor_mean[1:].mean(),color="#555555",ls="--",lw=.65)
    c.set(xlim=(.7,5.3),ylim=(.5,.7),xticks=L,xlabel="box length, $L$",ylabel=r"$S=A_-/(\Delta t_-|C_{\min}|)$"); c.set_title(r"shape factor: near collapse for L$\geq$2",loc="left",fontweight="bold"); c.text(.03,.05,fr"$S_{{2-5L}}={m.shape_factor_mean[1:].mean():.3f}$",transform=c.transAxes,fontsize=6)
    for j,ax in enumerate(axs): ax.text(.02,.96,"abc"[j],transform=ax.transAxes,va="top",fontsize=9,fontweight="bold"); ax.tick_params(length=2.5,width=.65)
    fig.text(.5,.01,r"Baseline 20 ns / 1 ps trajectories; normalized total $C_{vJ}$ at n=1 (therefore $k\propto L^{-1}$). First negative lobe zero-crossing bounded; points are mean of n=3, error bars replica SEM.",ha="center",fontsize=5.55)
    fig.subplots_adjust(left=.08,right=.99,bottom=.25,top=.86,wspace=.38)
    fig.savefig(str(OUT)+".png",dpi=600,bbox_inches="tight"); fig.savefig(str(OUT)+".tiff",dpi=600,bbox_inches="tight"); fig.savefig(str(OUT)+".pdf",bbox_inches="tight"); fig.savefig(str(OUT)+".svg",bbox_inches="tight")
if __name__=="__main__": main()
