"""Nature-style indicative comparison of long-duration, 1-ps chirality series.

Claim: within the completed three-replica low-frequency series, t_min grows
approximately linearly with relative length, but its slope and alpha_min are
chirality dependent.  The (17,0) historical points are intentionally shown as
exploratory because both the reference frame and replicate depth differ.
"""
from pathlib import Path
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TA=ROOT.parent
GCMC=ROOT.parents[1]
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.linewidth":.8,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42,"legend.frameon":False})
COL={"(7,7)":"#1B5E9A","(8,8)":"#D46A2E","(9,9)":"#39806E","(17,0)":"#7A7A7A"}
def rows(path):
    with open(path,newline="") as f:return list(csv.DictReader(f))
def main():
    data={}
    data["(7,7)"]=[(float(r["L_rel"]),float(r["alpha_min"]),float(r["alpha_min_sem"]),float(r["t_min_ps"]),3,"CNT-relative") for r in rows(TA/"06_cross_chirality_regimes"/"results"/"7_7_msd_alpha_characteristics_2L5L_3rep.csv")]
    data["(8,8)"]=[(float(r["L_rel"]),float(r["alpha_min_2_150ps"]),float(r["alpha_min_seed_sem"]),float(r["t_min_ps"]),3,"CNT-relative") for r in rows(GCMC/"analysis"/"linear_msd_1ps_3rep_8_8_20260730"/"alpha_z_1L5L_three_reps_summary.csv")]
    data["(9,9)"]=[(float(r["L_rel"]),float(r["alpha_min"]),float(r["alpha_min_sem"]),float(r["t_min_ps"]),3,"CNT-relative") for r in rows(TA/"06_cross_chirality_regimes"/"results"/"9_9_msd_preview_2L4L.csv")]
    data["(17,0)"]=[(1,.797,np.nan,4,1,"water-COM corrected"),(2,.7654,np.nan,8,1,"water-COM corrected")]
    fig=plt.figure(figsize=(5.5,2.45),dpi=300); ax1=fig.add_axes([.10,.22,.37,.64]); ax2=fig.add_axes([.59,.22,.37,.64])
    for name,vals in data.items():
        x=np.array([v[0] for v in vals]); a=np.array([v[1] for v in vals]); e=np.array([v[2] for v in vals]); t=np.array([v[3] for v in vals]);
        kw=dict(color=COL[name],lw=1.2,marker="o",ms=4.3)
        if name=="(17,0)":kw.update(markerfacecolor="white",ls="--")
        ax1.errorbar(x,a,yerr=e if np.isfinite(e).any() else None,**kw,label=name)
        ax2.plot(x,t/x,**kw)
    for ax,lab in [(ax1,r"$\alpha_{\min}$"),(ax2,r"$t_{\min}/L$ (ps per relative $L$)")]:
        ax.set(xlabel="Relative axial length, $L$",ylabel=lab,xticks=[1,2,3,4,5]);ax.tick_params(direction="out",length=3,width=.8)
    ax1.set(ylim=(.25,.85));ax2.set(ylim=(3,7));ax1.legend(title="Chirality",fontsize=6.5,title_fontsize=6.5,loc="upper right")
    ax1.text(-.17,1.06,"(a)",transform=ax1.transAxes,fontweight="bold",fontsize=8);ax2.text(-.17,1.06,"(b)",transform=ax2.transAxes,fontweight="bold",fontsize=8)
    fig.text(.10,.955,"Low-frequency long-duration screen: 1 ps sampling; completed points only",fontsize=7.3,fontweight="bold")
    fig.text(.10,.075,"(7,7)/(8,8)/(9,9): three independent 20 ns trajectories. (17,0): historical n=1 water-COM-corrected points; exploratory only.",fontsize=5.7,color="#454545")
    out=ROOT/"assets"/"lowfreq_crosschirality_indicative";fig.savefig(out.with_suffix(".png"),dpi=600);fig.savefig(out.with_suffix(".pdf"));fig.savefig(out.with_suffix(".svg"));fig.savefig(out.with_suffix(".tiff"),dpi=600)
if __name__=="__main__":main()
