"""Pair low-frequency alpha_z minima with first negative-lobe C_vJ(n=1) metrics."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT=Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
ALPHA=ROOT/"assets"/"lowfreq_alphaz_loglog_1decade_minima.csv"
CVJ=Path(r"H:\gcmc_explore\translational_anomaly\08_viscosity_friction_length_scaling\04_analysis\frequency_mode_response\full_15case_multirate_CJ_CvJ_1ns_20260803\CvJ_1ps_lag1ns")
PAIRED=ROOT/"assets"/"lowfreq_alphaz_cvj_n1_first_lobe_paired.csv"
STATS=ROOT/"assets"/"lowfreq_alphaz_cvj_n1_association_statistics.csv"
OUT=ROOT/"assets"/"lowfreq_alphaz_cvj_n1_association_nature"
COLORS=["#203864","#2F75B5","#4FA3A5","#5B9A6D","#A07A39"]
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.labelsize":7.1,"axes.titlesize":7.8,"xtick.labelsize":6.1,"ytick.labelsize":6.1,"axes.linewidth":.75,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})

def lobe_metrics(t,y):
    i=np.where(y[1:]<=0)[0][0]+1
    j=np.where(y[i+1:]>=0)[0][0]+i+1
    t0=t[i-1]-y[i-1]*(t[i]-t[i-1])/(y[i]-y[i-1])
    t1=t[j-1]-y[j-1]*(t[j]-t[j-1])/(y[j]-y[j-1])
    ti=np.r_[t0,t[i:j],t1]; yi=np.r_[0.,y[i:j],0.]
    imin=i+np.argmin(y[i:j])
    q=np.polyfit(t[imin-1:imin+2],y[imin-1:imin+2],2)
    return -np.trapz(yi,ti),t[imin],-q[1]/(2*q[0]),y[imin],t0,t1

def residual(df,col): return df[col]-df.groupby("box_length")[col].transform("mean")
def corr(x,y):
    p=pearsonr(x,y); s=spearmanr(x,y)
    return p.statistic,p.pvalue,s.statistic,s.pvalue

def guide(ax,x,y):
    c=np.polyfit(x,y,1); xx=np.linspace(np.min(x),np.max(x),100); ax.plot(xx,c[0]*xx+c[1],color="#263238",lw=.7,zorder=1)

def main():
    alpha=pd.read_csv(ALPHA); rows=[]
    for L in range(1,6):
        for r in range(1,4):
            d=pd.read_csv(CVJ/f"L{L}_rep{r}"/"tagged_current_mode_coupling.csv")
            d=d[d.n==1].sort_values("lag_index"); t=d.lag_index.to_numpy(float); y=d.C_vJ_total.to_numpy(float)
            area,tmin_raw,tmin_quad,cmin,tstart,tend=lobe_metrics(t,y)
            rows.append({"box_length":L,"replica":r,"alpha_z_min":alpha.loc[alpha.box_length==L,f"alpha_min_rep{r}"].iloc[0],"CvJ_n1_first_negative_lobe_area_ps":area,"CvJ_n1_first_minimum_raw_ps":tmin_raw,"CvJ_n1_first_minimum_parabolic_ps":tmin_quad,"CvJ_n1_first_minimum_value":cmin,"CvJ_n1_lobe_start_zero_ps":tstart,"CvJ_n1_lobe_end_zero_ps":tend,"CvJ_n1_lobe_width_ps":tend-tstart})
    df=pd.DataFrame(rows)
    df["alpha_z_min_within_L"]=residual(df,"alpha_z_min")
    for c in ["CvJ_n1_first_negative_lobe_area_ps","CvJ_n1_first_minimum_parabolic_ps"]: df[c+"_within_L"]=residual(df,c)
    df.to_csv(PAIRED,index=False)
    statrows=[]
    for metric in ["CvJ_n1_first_negative_lobe_area_ps","CvJ_n1_first_minimum_raw_ps","CvJ_n1_first_minimum_parabolic_ps","CvJ_n1_first_minimum_value","CvJ_n1_lobe_end_zero_ps","CvJ_n1_lobe_width_ps"]:
        for level,x,y in [("pooled",df[metric],df.alpha_z_min),("within_L",df[metric+"_within_L"] if metric+"_within_L" in df else residual(df,metric),df.alpha_z_min_within_L)]:
            pr,pp,sr,sp=corr(x,y); statrows.append({"metric":metric,"level":level,"n":len(df),"pearson_r":pr,"pearson_p_naive":pp,"spearman_rho":sr,"spearman_p_naive":sp})
    pd.DataFrame(statrows).to_csv(STATS,index=False)
    fig,axs=plt.subplots(2,2,figsize=(7.15,4.7)); panels=[("CvJ_n1_first_negative_lobe_area_ps","alpha_z_min",r"negative-lobe area, $A_-$ (ps)",r"$\alpha_{z,\min}$","pooled: $r=-0.93$"),("CvJ_n1_first_minimum_parabolic_ps","alpha_z_min",r"first minimum time, $t_{\min}^{CJ}$ (ps)",r"$\alpha_{z,\min}$","pooled: $r=-0.94$"),("CvJ_n1_first_negative_lobe_area_ps_within_L","alpha_z_min_within_L",r"within-L residual $A_-$ (ps)",r"within-L residual $\alpha_{z,\min}$",r"within-L: $r=-0.73$"),("CvJ_n1_first_minimum_parabolic_ps_within_L","alpha_z_min_within_L",r"within-L residual $t_{\min}^{CJ}$ (ps)",r"within-L residual $\alpha_{z,\min}$",r"within-L: $r=-0.42$, inconclusive")]
    for n,(ax,(x,y,xlab,ylab,note)) in enumerate(zip(axs.ravel(),panels)):
        for L in range(1,6):
            q=df[df.box_length==L]; ax.scatter(q[x],q[y],s=24,color=COLORS[L-1],alpha=.8,zorder=3)
        guide(ax,df[x],df[y]); ax.set(xlabel=xlab,ylabel=ylab); ax.set_title(["pooled association","pooled association","after removing box-length mean","after removing box-length mean"][n],loc="left",fontweight="bold"); ax.text(.02,.96,"abcd"[n],transform=ax.transAxes,va="top",fontsize=9,fontweight="bold"); ax.text(.98,.06,note,transform=ax.transAxes,ha="right",va="bottom",fontsize=5.5,color="#444444")
    fig.text(.5,.008,r"Paired 20 ns / 1 ps baseline trajectories; $C_{vJ}$ is normalized total axial current ACF at n=1 (fundamental k changes with L). First-lobe area is zero-crossing bounded; minimum time is parabolically interpolated from 1-ps samples. n=15; p-values in source table are descriptive.",ha="center",fontsize=5.45)
    fig.subplots_adjust(left=.1,right=.985,bottom=.16,top=.92,wspace=.35,hspace=.38)
    fig.savefig(str(OUT)+".png",dpi=600,bbox_inches="tight"); fig.savefig(str(OUT)+".tiff",dpi=600,bbox_inches="tight"); fig.savefig(str(OUT)+".pdf",bbox_inches="tight"); fig.savefig(str(OUT)+".svg",bbox_inches="tight")
if __name__=="__main__": main()
