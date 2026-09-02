"""Windowed, non-parametric alpha estimates from all-origin VACF ODE data.

For each replica, two trapezoidal integrations give J(t).  Instead of the
ill-conditioned point derivative t*I/J, this reports a finite log-time-window
slope of J (and therefore of MSD=2J):
alpha_r(t) = [ln J(t sqrt(r))-ln J(t/sqrt(r))] / ln(r).
The window ratio r is explicit; no arbitrary floor is applied to J.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
LOW=ROOT/'results/collective_mode_response/vacf_88_L2L5_10fs_1ns_8rep_weakNH_nomom/2026-08-21/allorigins'
TEN=ROOT/'results/collective_mode_response/vacf_tail_8_8_L10_10fs_8rep_1ns_2026-08-19/analysis_cvv_alpha_200ps_8rep'
OUT=ROOT/'results/collective_mode_response/alpha_windowed_ODE_from_allorigin_VACF/2026-08-21'
COLORS={'2L':'#3b6fb6','3L':'#e07b39','4L':'#4d9a65','5L':'#b34c66','10L':'#7656a6'}
RATIOS=(1.5,2.0,3.0)

def trap(y,dt):
    z=np.zeros_like(y); z[1:]=np.cumsum(.5*dt*(y[:-1]+y[1:]),axis=0); return z

def load(label):
    if label!='10L':
        x=np.genfromtxt(LOW/'per_replica'/f'VACF_8_8_{label}_peculiar_per_replica_normalised.csv',delimiter=',',names=True)
        return x['lag_ps'],np.column_stack([x[f'rep{i}'] for i in range(1,9)])
    x=np.genfromtxt(TEN/'cvv_per_replica.csv',delimiter=',',names=True); t=np.unique(x['lag_ps']); rep=x['replica']
    v=np.column_stack([x['cvv_A2_ps2'][rep==i] for i in range(1,9)]); return t,v/v[0:1]

def write(path,head,arr):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(head); w.writerows(arr)

def slope(t,j,ratio):
    result=np.full_like(j,np.nan); lo=t/np.sqrt(ratio); hi=t*np.sqrt(ratio)
    ok=(lo>=0.1)&(hi<=100.0)
    for k in range(j.shape[1]):
        vals=j[:,k]
        # J must be positive for a logarithmic MSD slope. Invalid regions stay NaN.
        if np.any(vals[1:]<=0): raise ValueError('non-positive J encountered beyond t=0')
        result[ok,k]=(np.interp(hi[ok],t[1:],np.log(vals[1:]))-np.interp(lo[ok],t[1:],np.log(vals[1:])))/np.log(ratio)
    return result

def style(ax):
    ax.tick_params(direction='out',width=1,length=3,labelsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1); ax.spines['bottom'].set_linewidth(1)

def main():
    OUT.mkdir(parents=True,exist_ok=True); (OUT/'figures').mkdir(exist_ok=True); (OUT/'per_replica').mkdir(exist_ok=True)
    data={}
    manifest={'definition':'alpha_r(t)=[ln J(t sqrt(r))-ln J(t/sqrt(r))]/ln r; J=int_0^t int_0^s Cvv(u)du ds','ratios':RATIOS,'protocol':'all-origin 8-rep axial peculiar VACF; current 10fs/1ns protocol','source':{'2L-5L':str(LOW.resolve()),'10L':str(TEN.resolve())},'lengths':{}}
    for L in ('2L','3L','4L','5L','10L'):
        t,c=load(L); keep=t<=100; t,c=t[keep],c[keep]; j=trap(trap(c,np.median(np.diff(t))),np.median(np.diff(t)))
        data[L]={}
        for r in RATIOS:
            a=slope(t,j,r); mean=np.full(t.size,np.nan); sem=np.full(t.size,np.nan); valid=np.isfinite(a[:,0])
            mean[valid]=np.nanmean(a[valid],axis=1); sem[valid]=np.nanstd(a[valid],axis=1,ddof=1)/np.sqrt(8)
            data[L][r]=(t,a,mean,sem)
            tag=str(r).replace('.','p'); write(OUT/f'alpha_window_r{tag}_8_8_{L}_mean_sem.csv',['lag_ps','alpha_window_mean','alpha_window_replica_sem','n_replicas'],np.column_stack([t,mean,sem,np.full_like(t,8)]))
            write(OUT/'per_replica'/f'alpha_window_r{tag}_8_8_{L}_per_replica.csv',['lag_ps']+[f'rep{i}' for i in range(1,9)],np.column_stack([t,a]))
        manifest['lengths'][L]={'points':int(t.size),'n_replica':8}
    fig=plt.figure(figsize=(7,2.55),dpi=300); axes=[fig.add_axes([.09,.23,.25,.66]),fig.add_axes([.41,.23,.25,.66]),fig.add_axes([.73,.23,.24,.66])]
    for ax,r,tag in zip(axes,RATIOS,('(a)','(b)','(c)')):
        for L in data:
            t,_,m,s=data[L][r]; ok=np.isfinite(m); ax.plot(t[ok],m[ok],lw=1.15,color=COLORS[L],label=L); ax.fill_between(t[ok],m[ok]-s[ok],m[ok]+s[ok],color=COLORS[L],alpha=.14,lw=0)
        ax.axhline(.5,color='.45',ls='--',lw=.8); ax.set(xscale='log',xlim=(.2,60),ylim=(.3,1.15),xlabel=r'$t$ (ps)',ylabel=rf'$\alpha_{{r={r:g}}}(t)$'); style(ax); ax.text(-.16,1.05,tag,transform=ax.transAxes,fontsize=9,fontweight='bold')
    axes[0].legend(frameon=False,fontsize=6.5,loc='best'); fig.savefig(OUT/'figures'/'alpha_windowed_ODE_r1p5_r2_r3_aggregate.png',dpi=600); fig.savefig(OUT/'figures'/'alpha_windowed_ODE_r1p5_r2_r3_aggregate.pdf'); plt.close(fig)
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    (OUT/'README.md').write_text('# Windowed ODE alpha from all-origin VACF\n\nThis package is a non-parametric smoothing of the VACF-ODE route. It uses finite log-time windows r=1.5, 2, 3 on J(t), not the point estimator tI/J. It is an effective-window exponent, not an instantaneous exponent. Replicas are the uncertainty units. Direct-MSD validation from the current dump protocol is deliberately separate and not substituted with protocol-distinct 100-fs/10-ns MSD assets.\n',encoding='utf-8')
if __name__=='__main__': main()
