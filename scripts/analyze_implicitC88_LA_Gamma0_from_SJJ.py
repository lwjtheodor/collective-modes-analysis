"""Infer implicit-CNT axial zero-k damping from finite-k longitudinal S_JJ."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def read(src,nmax):
 t=[];js=[];ls=[]
 with src.open() as f:
  while True:
   mark=f.readline()
   if not mark:break
   if mark.strip()!='ITEM: TIMESTEP':raise ValueError(mark)
   step=int(f.readline());f.readline();nat=int(f.readline());f.readline();b=[f.readline().split() for _ in range(3)];h=f.readline().split()[2:];c={v:i for i,v in enumerate(h)}
   if any(v not in c for v in ('type','z','vz')):raise ValueError(h)
   a=np.fromstring(' '.join(f.readline() for _ in range(nat)),sep=' ').reshape(nat,len(h));a=a[a[:,c['type']].astype(int)==1]
   z=a[:,c['z']];v=a[:,c['vz']];L=float(b[2][1])-float(b[2][0]);k=2*np.pi*np.arange(1,nmax+1)/L
   js.append(np.exp(1j*np.outer(k,z))@v);t.append(step*.0005);ls.append(L)
 return np.asarray(t),np.asarray(js).T,np.mean(ls)
def lor(w,b0,b1,A,w0,g):return b0+b1*w+A*g*g/((w-w0)**2+g*g)
def main():
 p=argparse.ArgumentParser();p.add_argument('--dumps',nargs='+',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--nmax',type=int,default=8);a=p.parse_args();d=a.output/'derived_data';fg=a.output/'figures';d.mkdir(parents=True,exist_ok=True);fg.mkdir(exist_ok=True)
 all=[];Ls=[]
 for src in a.dumps:
  t,j,L=read(src,a.nmax);dt=np.median(np.diff(t));win=np.hanning(len(t));all.append(np.abs(np.fft.fft((j-j.mean(axis=1,keepdims=True))*win,axis=1))**2/np.sum(win**2));Ls.append(L)
 A=np.asarray(all);w=np.fft.fftfreq(A.shape[-1],dt)*2*np.pi;pos=w>=0;w=w[pos];A=A[...,pos];mean=A.mean(0);sem=A.std(0,ddof=1)/np.sqrt(len(A));rows=[];sp=[];fig,axs=plt.subplots(2,4,figsize=(13,6),squeeze=False)
 for n in range(1,a.nmax+1):
  y=mean[n-1];e=np.maximum(sem[n-1],np.nanmax(sem[n-1])*1e-5); lo=max(.004,.005*n);hi=.20*n;m=(w>=lo)&(w<=hi);x=w[m];yy=y[m];ee=e[m]; seed=x[np.argmax(yy)]; half=max(.012,.45*seed); fm=(w>=max(.002,seed-half))&(w<=seed+half);x=w[fm];yy=y[fm];ee=e[fm]
  p0=(max(0,yy.min()),0,yy.max()-yy.min(),seed,max(w[1]*2,.05*seed));bounds=([0,-np.inf,0,x.min(),w[1]*.5],[np.inf,np.inf,np.inf,x.max(),max(.5,1.5*seed)])
  try:
   q,cov=curve_fit(lor,x,yy,p0=p0,sigma=ee,absolute_sigma=True,bounds=bounds,maxfev=100000);pred=lor(x,*q);r2=1-np.sum((yy-pred)**2)/np.sum((yy-yy.mean())**2);gse=np.sqrt(cov[4,4]); accepted=(q[4]>1.5*w[1]) and not np.isclose(q[4],bounds[0][4]) and not np.isclose(q[4],bounds[1][4]) and r2>.7
  except RuntimeError:
   q=[np.nan]*5;pred=np.full_like(x,np.nan);r2=np.nan;gse=np.nan;accepted=False
  rows.append({'n':n,'k_inv_A':2*np.pi*n/np.mean(Ls),'omega_peak_rad_ps':q[3],'gamma_HWHM_rad_ps':q[4],'gamma_fit_SEM_rad_ps':gse,'fit_R2':r2,'frequency_bin_rad_ps':w[1],'accepted':accepted})
  sp.extend({'n':n,'omega_rad_ps':ww,'S_mean':ss,'S_replica_sem':eev} for ww,ss,eev in zip(w,mean[n-1],sem[n-1]));ax=axs.ravel()[n-1];ax.semilogy(x,yy,color='#2166ac');ax.semilogy(x,pred,color='#b2182b');ax.set_title(f'n={n}');ax.grid(alpha=.2)
 fig.supxlabel(r'$\omega$ (rad ps$^{-1}$)');fig.supylabel(r'$S_{J_zJ_z}$');fig.tight_layout();fig.savefig(fg/'implicitC88_LA_SJJ_n001_n008.png',dpi=300);fig.savefig(fg/'implicitC88_LA_SJJ_n001_n008.pdf');plt.close(fig)
 q=pd.DataFrame(rows);q.to_csv(d/'implicitC88_LA_lorentzian_widths_n001_n008.csv',index=False);pd.DataFrame(sp).to_csv(d/'implicitC88_LA_SJJ_n001_n008.csv',index=False)
 use=q[q.accepted].copy()
 def mod(k,g0,aa,al):return g0+aa*k**al
 if len(use)>=4:
  par,cov=curve_fit(mod,use.k_inv_A,use.gamma_HWHM_rad_ps,p0=(.01,20,1.5),sigma=use.gamma_fit_SEM_rad_ps,absolute_sigma=True,bounds=([0,0,.1],[np.inf,np.inf,4]),maxfev=100000);se=np.sqrt(np.diag(cov)); pred=mod(use.k_inv_A,*par);result={'model':'Gamma0+A*k^alpha','Gamma0_rad_ps':par[0],'Gamma0_SE':se[0],'A_rad_ps_A_to_alpha':par[1],'A_SE':se[1],'alpha':par[2],'alpha_SE':se[2],'n_used':','.join(map(str,use.n.astype(int)))}
 else: result={'model':'insufficient_accepted_modes','n_used':','.join(map(str,use.n.astype(int)))}
 pd.DataFrame([result]).to_csv(d/'implicitC88_LA_Gamma0_powerlaw_fit.csv',index=False);(a.output/'metadata.json').write_text(json.dumps({'definition':'J_n=sum_O vz exp(i k_n z)','source':'C88 N1600 350K weak-NH 6ns 100fs 4rep','interpretation':'Gamma0 is a finite-k spectral extrapolation; total Pz is numerically zero and is not used as an anchor'},indent=2));(a.output/'FINISHED.txt').write_text('Implicit C88 finite-k LA linewidth extrapolation finished successfully.\n')
if __name__=='__main__':main()
