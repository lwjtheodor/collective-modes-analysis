"""Aggregate 1--10L alpha(Delta t) and absolute-current first-lobe data."""
from __future__ import annotations
import csv, json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({"font.family":"Arial","font.size":8,"svg.fonttype":"none","pdf.fonttype":42,"axes.linewidth":0.8,"legend.frameon":False})
ROOT=Path(__file__).resolve().parents[1]
ALPHA=ROOT/'heartbeat_fetch'/'stage_alpha200ps_crosschirality_20260809'/'output'
ABS=ROOT/'heartbeat_fetch'/'stage_absolute_cvj_1L10L_20260809'/'output'
OUT=ROOT/'assets'; CHIS=['7_7','8_8','9_9','17_0']; LAB={'7_7':'(7,7)','8_8':'(8,8)','9_9':'(9,9)','17_0':'(17,0)'}
COL={1:'#355F8C',2:'#4F8F9A',3:'#75A65A',4:'#D5A044',5:'#C86B3C',10:'#A34358'}

def parse_id(f):
    p=f.name.split('_'); return '_'.join(p[:2]),int(p[3][1:])
def load_csv(f): return np.genfromtxt(f,delimiter=',',names=True)
def sem(a): return np.std(a,axis=0,ddof=1)/np.sqrt(a.shape[0]) if a.shape[0]>1 else np.zeros(a.shape[1])
def save(fig,stem):
    fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.tiff'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')

def alpha_figure():
    rows=[]; fig,axs=plt.subplots(2,2,figsize=(7.25,5.6),constrained_layout=True)
    for j,(ax,chi) in enumerate(zip(axs.flat,CHIS)):
      for L in [1,2,3,4,5,10]:
        fs=sorted(ALPHA.glob(f'{chi}_L{L}_rep*_msd.csv'))
        ds=[load_csv(f) for f in fs]; x=ds[0]['lag_ps']; x=x[x>=5.0]
        # A complete trailing decade [Delta t/10, Delta t] is available for every
        # reported point through 200 ps; this avoids endpoint-truncated slopes.
        y=[]
        for d in ds:
          xx=d['lag_ps']; mm=d['msd_z_A2']; slopes=[]
          for tt in x:
            keep=(xx>=tt/10.0)&(xx<=tt)
            slopes.append(np.polyfit(np.log(xx[keep]),np.log(mm[keep]),1)[0] if keep.sum()>=5 else np.nan)
          y.append(slopes)
        y=np.asarray(y); m=np.nanmean(y,axis=0); e=np.nanstd(y,axis=0,ddof=1)/np.sqrt(len(y)) if len(y)>1 else np.zeros_like(m)
        for xx,yy,ee in zip(x,m,e): rows.append({'chirality':chi,'L':L,'lag_ps':xx,'alpha_mean':yy,'alpha_sem':ee,'n_replicates':len(ds)})
        ax.plot(x,m,color=COL[L],lw=1.25,label=f'{L}L'); ax.fill_between(x,m-e,m+e,color=COL[L],alpha=.16,lw=0)
      ax.axhline(1,color='#484848',lw=.7,ls='--');ax.set_xscale('log');ax.set_xlim(.5,200);ax.set_ylim(.35,1.10);ax.set_title(f"{chr(97+j)}  {LAB[chi]}",loc='left',fontweight='bold',fontsize=10)
      ax.text(.04,.08,'direct MSD; trailing one-decade fit\nwater-COM coordinate removed\nmean ± replica SEM',transform=ax.transAxes,fontsize=6.4,va='bottom')
      if j%2==0: ax.set_ylabel(r'$\alpha_z(\Delta t)$')
      if j>=2: ax.set_xlabel(r'lag time, $\Delta t$ (ps)')
      if j==0: ax.legend(ncol=2,fontsize=6.4,loc='upper right',handlelength=1.4)
    fig.suptitle(r'Decade-uniform local MSD exponent through 200 ps',fontweight='bold',fontsize=12,y=1.015)
    fig.text(.5,-.025,'Weak Nosé–Hoover; no momentum removal; 1–5L: 10 fs/1 ns, 10L: 100 fs/10 ns; trailing one-decade log–log slope [Delta t/10, Delta t]; shaded bands are replica SEM.',ha='center',fontsize=7)
    save(fig,OUT/'crosschirality_1L10L_alpha200ps_trailing_decade_nature');plt.close(fig)
    with (OUT/'crosschirality_1L10L_alpha200ps_trailing_decade_curves.csv').open('w',newline='',encoding='utf-8') as h:
      w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def abs_figure():
    metrics=[('first_negative_lobe_absolute','absolute area, $A_-^{abs}$ ($\mathrm{\AA^2/ps}$)'),('first_negative_lobe_per_water','per-water absolute area, $A_-^{abs}/N_w$ ($\mathrm{\AA^2/ps}$)'),('first_negative_lobe_normalized','normalized area, $A_-^{norm}$ (ps)')]
    rows=[]; fig,axs=plt.subplots(1,3,figsize=(7.25,2.75),constrained_layout=True)
    vals={m: {chi:[] for chi in CHIS} for m,_ in metrics}
    for chi in CHIS:
      for L in [1,2,3,4,5,10]:
        fs=sorted(ABS.glob(f'{chi}_L{L}_rep*_absolute_cvj.json'))
        for m,_ in metrics:
          a=np.array([json.loads(f.read_text())[m]['negative_area'] for f in fs]); mean=float(a.mean()); se=float(a.std(ddof=1)/np.sqrt(len(a))) if len(a)>1 else 0.
          vals[m][chi].append((L,mean,se));rows.append({'chirality':chi,'L':L,'metric':m,'area_mean':mean,'area_sem':se,'n_replicates':len(a)})
    cchi={'7_7':'#1F77B4','8_8':'#D1495B','9_9':'#5B9A3D','17_0':'#7A5195'}
    for ax,(m,title) in zip(axs,metrics):
      for chi in CHIS:
        a=np.asarray(vals[m][chi]);ax.errorbar(a[:,0],a[:,1],yerr=a[:,2],color=cchi[chi],marker='o',ms=3.5,lw=1.3,capsize=2,label=LAB[chi])
      ax.set_xlabel('box length (L)');ax.set_title(title,fontsize=8.2);ax.set_xticks([1,2,3,4,5,10]);ax.tick_params(labelsize=7)
      if m!='first_negative_lobe_normalized': ax.set_yscale('log')
    axs[0].set_ylabel('first-negative-lobe area');axs[-1].legend(fontsize=6.4,loc='upper left',handlelength=1.3)
    fig.suptitle(r'First negative lobe: raw scale, per-water scale and normalized morphology',fontweight='bold',fontsize=11,y=1.04)
    fig.text(.5,-.06,r'$J_k=\sum_i(v_{z,i}-\bar v_z)e^{ikz_i}$; weak Nosé–Hoover, no momentum removal; 1–5L: 10 fs/1 ns, 10L: 100 fs/10 ns; points mean ± replica SEM.',ha='center',fontsize=6.7)
    save(fig,OUT/'crosschirality_1L10L_absolute_cvj_lobe_normalizations_nature');plt.close(fig)
    with (OUT/'crosschirality_1L10L_absolute_cvj_lobe_normalizations.csv').open('w',newline='',encoding='utf-8') as h:
      w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__': OUT.mkdir(exist_ok=True);alpha_figure();abs_figure()
