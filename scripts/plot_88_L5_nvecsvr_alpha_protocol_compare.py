"""Compare 5L NVE+CSVR VACF-ODE alpha with existing protocol-distinct length series."""
from pathlib import Path
import csv, numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
NEW=ROOT/'results/collective_mode_response/vacf_compact_HB_analysis_20260821'
OLD=ROOT/'results/collective_mode_response/vacf_alpha_88_L2L5L10_10fs_1ns_8rep/2026-08-21/allorigins'
OUT=ROOT/'results/collective_mode_response/nvecsvr_88_L5_vacf_alpha_protocol_compare/2026-08-21'
OUT.mkdir(parents=True,exist_ok=True)

def alpha(t,c):
    I=np.r_[0,np.cumsum((c[1:]+c[:-1])*0.5*np.diff(t))]
    J=np.r_[0,np.cumsum((I[1:]+I[:-1])*0.5*np.diff(t))]
    a=np.full_like(t,np.nan); a[1:]=t[1:]*I[1:]/J[1:]
    return a
def ms(a): return a.mean(0),a.std(0,ddof=1)/np.sqrt(len(a))
def row(t,a,lo=5,hi=35):
    m=(t>=lo)&(t<=hi);i=np.flatnonzero(m)[np.nanargmin(a[m])]
    return float(a[i]),float(t[i]),float(np.interp(20,t,a)),float(np.interp(50,t,a))
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.linewidth':1,'xtick.direction':'out','ytick.direction':'out'})
new={}
for d in (0,100,200,500):
    curves=[]
    for r in range(1,9):
        x=np.genfromtxt(NEW/f'VACF_8_8_L5_rep{r}_discard{d}ps.csv',delimiter=',',names=True)
        curves.append(alpha(x['lag_ps'],x['vacf_peculiar_mean']))
    new[d]=(x['lag_ps'],np.array(curves))
old={}
for L in (2,3,4,5,10):
    x=np.genfromtxt(OLD/f'per_replica/alpha_ODE_8_8_{L}L_per_replica.csv',delimiter=',',names=True)
    old[L]=(x['lag_ps'],np.array([x[f'rep{i}'] for i in range(1,9)]))
metrics=[]
for d,(t,c) in new.items():
    v=np.array([row(t,z) for z in c]); m,s=ms(v)
    metrics.append([f'NVE+CSVR discard {d} ps',*m,*s])
t,c=old[5];v=np.array([row(t,z) for z in c]);m,s=ms(v);metrics.append(['weak-NH legacy 5L',*m,*s])
with (OUT/'alpha_metrics_5L_protocol_and_exclusions.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['series','alpha_min_mean','tmin_ps_mean','alpha20_mean','alpha50_mean','alpha_min_seed_sem','tmin_seed_sem','alpha20_seed_sem','alpha50_seed_sem']);w.writerows(metrics)
fig=plt.figure(figsize=(7,4.6));ax1=fig.add_axes([.11,.57,.38,.34]);ax2=fig.add_axes([.57,.57,.38,.34]);ax3=fig.add_axes([.11,.12,.38,.34]);ax4=fig.add_axes([.57,.12,.38,.34])
cols={0:'#0072B2',100:'#56B4E9',200:'#009E73',500:'#D55E00'}
for d,(t,c) in new.items():
 m,s=ms(c);ax1.plot(t,m,color=cols[d],lw=1.1,label=f'NVE+CSVR, discard {d} ps');ax1.fill_between(t,m-s,m+s,color=cols[d],alpha=.13,lw=0)
t,c=old[5];m,s=ms(c);ax1.plot(t,m,color='0.25',lw=1.1,ls='--',label='weak-NH, legacy 5L');ax1.fill_between(t,m-s,m+s,color='0.25',alpha=.10,lw=0)
ax1.set(xlim=(0,100),ylim=(-.03,.04),xlabel=r'$t$ (ps)',ylabel=r'$C_{vv}(t)/C_{vv}(0)$');ax1.axhline(0,color='0.5',lw=.7);ax1.legend(frameon=False,fontsize=5,loc='upper right')
for d,(t,c) in new.items():m,s=ms(c);ax2.plot(t,m,color=cols[d],lw=1.1,label=str(d));ax2.fill_between(t,m-s,m+s,color=cols[d],alpha=.13,lw=0)
t,c=old[5];m,s=ms(c);ax2.plot(t,m,color='0.25',lw=1.1,ls='--')
ax2.set(xlim=(1,60),ylim=(.25,1.05),xlabel=r'$t$ (ps)',ylabel=r'$alpha_{mathrm{VACF}}(t)$');ax2.legend(title='discard (ps)',frameon=False,fontsize=5,title_fontsize=5)
for L,(t,c) in old.items():
 m,s=ms(c);ax3.plot(t,m,lw=1.1,label=f'{L}L weak-NH');ax3.fill_between(t,m-s,m+s,alpha=.10,lw=0)
t,c=new[0];m,s=ms(c);ax3.plot(t,m,color='#D55E00',lw=1.2,label='5L NVE+CSVR');ax3.fill_between(t,m-s,m+s,color='#D55E00',alpha=.13,lw=0)
ax3.set(xlim=(1,60),ylim=(.2,1.05),xlabel=r'$t$ (ps)',ylabel=r'$alpha_{mathrm{VACF}}(t)$');ax3.legend(frameon=False,fontsize=5,ncol=2)
labs=[x[0] for x in metrics];vals=np.array([x[1] for x in metrics],float);errs=np.array([x[5] for x in metrics],float)
ax4.errorbar(np.arange(len(labs)),vals,yerr=errs,fmt='o',capsize=2,color='#0072B2');ax4.set(xticks=np.arange(len(labs)),xticklabels=['0','100','200','500','weak-NH'],ylim=(.25,.75),ylabel=r'$alpha_{min}$ (5--35 ps)',xlabel='NVE+CSVR discard (ps) / legacy');ax4.tick_params(axis='x',labelsize=6)
for lab,ax in zip(['(a)','(b)','(c)','(d)'],[ax1,ax2,ax3,ax4]):ax.text(-.17,1.05,lab,transform=ax.transAxes,fontweight='bold',fontsize=9)
fig.savefig(OUT/'NVECSVR_5L_VACF_alpha_protocol_compare.png',dpi=300);fig.savefig(OUT/'NVECSVR_5L_VACF_alpha_protocol_compare.pdf')
