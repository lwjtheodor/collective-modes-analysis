"""Archive and plot de-COM flexible-versus-fixed 4L collective-mode evidence."""
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(r'H:\gcmc_explore\translational_anomaly\02_isf_collective_modes\results\flexible_fixed_4L_20260809')
OUT=Path(r'H:\gcmc_explore\translational_anomaly\02_isf_collective_modes\results\flexible_fixed_collective_modes\2026-08-10')
OUT.mkdir(parents=True,exist_ok=True)
COL={'flex':'#b40426','fixed':'#3b4cc0'}

def dho(t,g,w,b): return np.exp(-g*t)*(np.cos(w*t)+b*np.sin(w*t))
def fit_dho(frame):
    out=[]
    for n in range(2,11):
        q=frame[frame.n==n].sort_values('lag_ps'); t=q.lag_ps.to_numpy(); y=q.C_J_norm.to_numpy(); keep=t<=25; t=t[keep]; y=y[keep]
        f=np.fft.rfftfreq(len(t),np.median(np.diff(t))); j=np.argmax(abs(np.fft.rfft(y*np.hanning(len(y))))[1:])+1; w0=max(2*np.pi*f[j],.02)
        p,_=curve_fit(dho,t,y,p0=(.1,w0,0),bounds=([0,0,-10],[10,20,10]),maxfev=50000)
        pred=dho(t,*p); r2=1-((y-pred)**2).sum()/((y-y.mean())**2).sum()
        out.append({'n':n,'k_inv_A':q.k_inv_A.iloc[0],'Gamma_ps_inv':p[0],'omega_rad_ps':p[1],'dho_sine':p[2],'DHO_R2':r2})
    return pd.DataFrame(out)
def power_fit(d):
    x=np.log(d.k_inv_A.to_numpy());y=np.log(d.Gamma_ps_inv.to_numpy());z,loga=np.polyfit(x,y,1)
    return {'A_ps_inv_Az':np.exp(loga),'z':z,'log_R2':1-((y-(loga+z*x))**2).sum()/((y-y.mean())**2).sum(),'n_points':len(d)}
def fixed_mean(frames, column, time, n):
    # Each CSV concatenates n=1..10, so restrict before interpolation.  Using
    # the full non-monotonic lag column aliases modes into a false short period.
    vals=[]
    for x in frames:
        q=x[x.n==n].sort_values('lag_ps')
        vals.append(np.interp(time,q.lag_ps,q[column]))
    a=np.asarray(vals);return a.mean(0),a.std(0,ddof=1)/math.sqrt(len(a))
def style(ax):
    ax.tick_params(direction='out',width=1,length=3,labelsize=7)
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1);ax.spines['bottom'].set_linewidth(1)

flex=pd.read_csv(ROOT/'caj/flexible_CJ_CaJ.csv')
fixed=[pd.read_csv(ROOT/f'caj/fixed_rep{i}_CJ_CaJ.csv') for i in (1,2,3)]
fd=fit_dho(flex); fds=[fit_dho(x) for x in fixed]
damp=fd.copy();damp['system']='flexible';damp['replicate']=1
for i,x in enumerate(fds,1):
    q=x.copy();q['system']='fixed';q['replicate']=i;damp=pd.concat([damp,q],ignore_index=True)
fitrows=[]
q=power_fit(fd);q.update(system='flexible',replicate=1);fitrows.append(q)
for i,x in enumerate(fds,1):q=power_fit(x);q.update(system='fixed',replicate=i);fitrows.append(q)
fits=pd.DataFrame(fitrows)

# All-n C_J archive on flexible's 100-fs grid.
cjrows=[]
for n in range(1,11):
    q=flex[flex.n==n].sort_values('lag_ps');t=q.lag_ps.to_numpy();m,se=fixed_mean(fixed,'C_J_norm',t,n)
    for ti,a,b,c in zip(t,q.C_J_norm,m,se):cjrows.append({'n':n,'k_inv_A':q.k_inv_A.iloc[0],'lag_ps':ti,'C_J_flexible':a,'C_J_fixed_mean':b,'C_J_fixed_sem':c})
pd.DataFrame(cjrows).to_csv(OUT/'cj_curves_all_n_decom.csv',index=False)

# Existing tagged-current decomposition; total is identity-equivalent to C_J, archive distinct explicitly.
cvf=pd.read_csv(ROOT/'cvj_flexible/tagged_current_mode_coupling.csv');cvx=[pd.read_csv(ROOT/f'cvj_fixed/rep{i}/tagged_current_mode_coupling.csv') for i in (1,2,3)]
cvrows=[]
for n in range(1,11):
    q=cvf[cvf.n==n].sort_values('lag_index');t=q.lag_index.to_numpy()*.1
    arrays=[]
    for x in cvx:
        u=x[x.n==n].sort_values('lag_index');arrays.append(np.interp(t,u.lag_index.to_numpy()*.01,u.C_vJ_distinct))
    a=np.asarray(arrays);m=a.mean(0);se=a.std(0,ddof=1)/math.sqrt(3)
    for ti,ff,mm,ss in zip(t,q.C_vJ_distinct,m,se):cvrows.append({'n':n,'k_inv_A':q.k_inv_A.iloc[0],'lag_ps':ti,'C_vJ_distinct_flexible':ff,'C_vJ_distinct_fixed_mean':mm,'C_vJ_distinct_fixed_sem':ss})
pd.DataFrame(cvrows).to_csv(OUT/'cvj_distinct_curves_all_n_decom.csv',index=False)
damp.to_csv(OUT/'damping_dho_points_n2_n10.csv',index=False);fits.to_csv(OUT/'damping_powerlaw_fits_n2_n10.csv',index=False)

ms=[]
for sys in ('flexible','fixed'):
    x=pd.read_csv(ROOT/'output'/f'{sys}_msd_decom');x['system']=sys;ms.append(x)
pd.concat(ms,ignore_index=True).to_csv(OUT/'msd_decom_comparison.csv',index=False)

# Explicit bbox layout: 2x2, fixed axes locations.
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.linewidth':1,'pdf.fonttype':42,'ps.fonttype':42})
fig=plt.figure(figsize=(5.5,4.6),dpi=300)
axes=[fig.add_axes(b) for b in ([.12,.59,.35,.31],[.60,.59,.35,.31],[.12,.13,.35,.31],[.60,.13,.35,.31])]
for ax,label in zip(axes,['(a)','(b)','(c)','(d)']):ax.text(-.19,1.07,label,transform=ax.transAxes,fontsize=9,fontweight='bold')
for ax,n in zip(axes[:2],[1,2]):
    q=pd.DataFrame(cjrows if n==1 else cvrows);q=q[q.n==n]
    if n==1:
        ax.plot(q.lag_ps,q.C_J_flexible,color=COL['flex'],lw=1.2,label='flexible (1 traj.)');ax.plot(q.lag_ps,q.C_J_fixed_mean,color=COL['fixed'],lw=1.2,label='fixed (3 traj.)');ax.fill_between(q.lag_ps,q.C_J_fixed_mean-q.C_J_fixed_sem,q.C_J_fixed_mean+q.C_J_fixed_sem,color=COL['fixed'],alpha=.18,lw=0);ax.set_ylabel(r'$C_J(k,t)/C_J(k,0)$');ax.set_xlim(0,35);ax.set_ylim(-.75,1.05);ax.set_xlabel(r'$t$ (ps), $n=1$');ax.legend(frameon=False,fontsize=6,loc='upper right')
    else:
        ax.plot(q.lag_ps,q.C_vJ_distinct_flexible,color=COL['flex'],lw=1.2,label='flexible');ax.plot(q.lag_ps,q.C_vJ_distinct_fixed_mean,color=COL['fixed'],lw=1.2,label='fixed');ax.fill_between(q.lag_ps,q.C_vJ_distinct_fixed_mean-q.C_vJ_distinct_fixed_sem,q.C_vJ_distinct_fixed_mean+q.C_vJ_distinct_fixed_sem,color=COL['fixed'],alpha=.18,lw=0);ax.set_ylabel(r'$C_{vJ}^{\mathrm{distinct}}(k,t)$');ax.set_xlim(0,20);ax.set_ylim(-.65,.12);ax.set_xlabel(r'$t$ (ps), $n=2$')
    ax.axhline(0,color='.45',lw=1);style(ax)

ax=axes[2]
for sys,color,marker in [('fixed',COL['fixed'],'o'),('flexible',COL['flex'],'s')]:
 q=damp[damp.system==sys];
 if sys=='fixed':
  g=q.groupby('n',as_index=False).Gamma_ps_inv.agg(['mean','sem']).reset_index();kk=damp[damp.system==sys].groupby('n').k_inv_A.mean().to_numpy();ax.errorbar(kk,g['mean'],yerr=g['sem'],fmt=marker,color=color,ms=3,lw=1,capsize=2,label='fixed (mean ± SEM)');z=fits[fits.system==sys].z.mean();A=np.exp(np.log(fits[fits.system==sys].A_ps_inv_Az).mean())
 else: kk=q.k_inv_A.to_numpy();ax.plot(kk,q.Gamma_ps_inv,marker,color=color,ms=3,lw=0,label='flexible (1 traj.)');z=fits[fits.system==sys].z.iloc[0];A=fits[fits.system==sys].A_ps_inv_Az.iloc[0]
 xx=np.linspace(.028,.16,200);ax.plot(xx,A*xx**z,color=color,lw=1,ls='--');
ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'$k$ ($\mathrm{\AA}^{-1}$)');ax.set_ylabel(r'$\Gamma$ ($\mathrm{ps}^{-1}$)');ax.legend(frameon=False,fontsize=6,loc='upper left');style(ax)

ax=axes[3]
for sys,color in [('fixed',COL['fixed']),('flexible',COL['flex'])]:
 x=pd.read_csv(ROOT/'output'/f'{sys}_msd_decom');ax.plot(x.time_ps,x.msd_z_A2_mean,color=color,lw=1.2,label=sys);ax.fill_between(x.time_ps,x.msd_z_A2_mean-x.msd_z_A2_block_sem,x.msd_z_A2_mean+x.msd_z_A2_block_sem,color=color,alpha=.18,lw=0)
ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(3,32);ax.set_xlabel(r'$t$ (ps)');ax.set_ylabel(r'de-COM $\mathrm{MSD}_z$ ($\mathrm{\AA}^2$)');ax.legend(frameon=False,fontsize=6,loc='upper left');style(ax)
fig.savefig(OUT/'flexible_fixed_4L_collective_mode_comparison.png',dpi=600);fig.savefig(OUT/'flexible_fixed_4L_collective_mode_comparison.pdf');plt.close(fig)

# Separate cadence-sensitivity audit for kinematic C_aJ; never label as a force mode.
fig=plt.figure(figsize=(5.5,2.35),dpi=300)
axs=[fig.add_axes(b) for b in ([.12,.22,.35,.66],[.60,.22,.35,.66])]
for ax,n,lab in zip(axs,[1,2],['(a)','(b)']):
    q=flex[flex.n==n].sort_values('lag_ps');t=q.lag_ps.to_numpy();m,se=fixed_mean(fixed,'C_aJ_corrcoef',t,n)
    ax.plot(t,q.C_aJ_corrcoef,color=COL['flex'],lw=1.2,label='flexible, 100 fs')
    ax.plot(t,m,color=COL['fixed'],lw=1.2,label='fixed, 10 fs')
    ax.fill_between(t,m-se,m+se,color=COL['fixed'],alpha=.18,lw=0);ax.axhline(0,color='.45',lw=1)
    ax.set_xlim(0,15);ax.set_xlabel(r'$t$ (ps)');ax.set_ylabel(r'$C_{aJ}$ (corr. coeff.)');ax.text(-.19,1.07,lab,transform=ax.transAxes,fontsize=9,fontweight='bold');style(ax)
axs[0].legend(frameon=False,fontsize=6,loc='upper right')
fig.savefig(OUT/'kinematic_caj_cadence_sensitivity.png',dpi=600);fig.savefig(OUT/'kinematic_caj_cadence_sensitivity.pdf');plt.close(fig)

(OUT/'README.md').write_text('''# Flexible versus fixed CNT, (8,8) 4L\n\nFigure and CSV archives are de-COM water observables. Fixed curves are means and replica SEM from three 10-fs trajectories; flexible is one 100-fs trajectory and has no replica SEM. The DHO power-law fit uses n=2-10 and 0-25 ps. `C_vJ_distinct` is the tagged-current total minus its self term; its phase-tagged total is identity-equivalent to `C_J`, so it is not duplicated as an independent panel. `C_aJ` is intentionally excluded from the figure because it is central-difference acceleration, not force.\n''',encoding='utf-8')
