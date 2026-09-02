"""Protocol-matched 10L LA spectra and Gamma0+A*k**alpha (n=3..8)."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def read(path,nmax):
    ts=[]; jj=[]; L=[]
    with path.open() as f:
      while True:
        m=f.readline()
        if not m: break
        if m.strip()!='ITEM: TIMESTEP': raise ValueError(m)
        st=int(f.readline());f.readline();n=int(f.readline());f.readline(); b=[f.readline().split() for _ in range(3)]; h=f.readline().split()[2:]
        if h!=['id','z','vz']: raise ValueError(h)
        a=np.fromstring(' '.join(f.readline() for _ in range(n)),sep=' ').reshape(n,3); z=a[:,1];v=a[:,2]; l=float(b[2][1])-float(b[2][0]); k=2*np.pi*np.arange(1,nmax+1)/l
        jj.append(np.exp(1j*np.outer(k,z))@v);ts.append(st*.0005);L.append(l)
    return np.asarray(ts),np.asarray(jj).T,float(np.mean(L))
def line(w,b0,b1,A,w0,g): return b0+b1*w+A*g*g/((w-w0)**2+g*g)
def main():
 p=argparse.ArgumentParser();p.add_argument('--dumps',nargs='+',type=Path,required=True);p.add_argument('--gamma0-summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--nmax',type=int,default=8);a=p.parse_args();d=a.output/'derived_data';f=a.output/'figures';d.mkdir(parents=True,exist_ok=True);f.mkdir(exist_ok=True)
 all=[];Ls=[]
 for src in a.dumps:
  t,j,l=read(src,a.nmax);dt=np.median(np.diff(t));win=np.hanning(len(t));all.append(np.abs(np.fft.fft((j-j.mean(axis=1,keepdims=True))*win,axis=1))**2/np.sum(win**2));Ls.append(l)
 A=np.asarray(all);w=np.fft.fftfreq(A.shape[-1],dt)*2*np.pi;pos=w>=0;w=w[pos];A=A[...,pos];mean=A.mean(0);sem=A.std(0,ddof=1)/np.sqrt(len(A));rows=[]
 fig,axs=plt.subplots(2,3,figsize=(11,6),squeeze=False)
 for n in range(3,9):
  y=mean[n-1]; e=np.maximum(sem[n-1],np.nanmax(sem[n-1])*1e-5); seed=.105*n; lo=max(.02,.50*seed);hi=1.65*seed;m=(w>=lo)&(w<=hi);x=w[m];yy=y[m];ee=e[m]; j=np.argmax(yy);g0=max(w[1]*2,seed*.08)
  p0=(max(0,np.min(yy)),0,max(yy)-min(yy),x[j],g0);bound=([0,-np.inf,0,lo,w[1]*.5],[np.inf,np.inf,np.inf,hi,1.5])
  p,c=curve_fit(line,x,yy,p0=p0,sigma=ee,absolute_sigma=True,bounds=bound,maxfev=100000);pred=line(x,*p);r2=1-np.sum((yy-pred)**2)/np.sum((yy-yy.mean())**2);gse=np.sqrt(c[4,4])
  rows.append({'n':n,'k_inv_A':2*np.pi*n/np.mean(Ls),'omega_fit_rad_ps':p[3],'gamma_HWHM_rad_ps':p[4],'gamma_fit_SEM_rad_ps':gse,'fit_R2':r2,'frequency_bin_rad_ps':w[1]})
  ax=axs.ravel()[n-3];ax.semilogy(x,yy,color='#2166ac');ax.semilogy(x,pred,color='#b2182b');ax.set_title(f'n={n}');ax.grid(alpha=.2)
 fig.supxlabel(r'$\omega$ (rad ps$^{-1}$)');fig.supylabel(r'$S_{JJ}$');fig.tight_layout();fig.savefig(f/'LA_spectra_n003_n008.png',dpi=300);fig.savefig(f/'LA_spectra_n003_n008.pdf');plt.close(fig)
 q=pd.DataFrame(rows);z=pd.read_csv(a.gamma0_summary).iloc[0];G0=float(z.Gamma_acf_rad_ps);G0se=G0*.0 # prior result does not archive independent SEM
 def mod(k,Aa,al):return G0+Aa*k**al
 p,c=curve_fit(mod,q.k_inv_A,q.gamma_HWHM_rad_ps,p0=(20,1.5),sigma=q.gamma_fit_SEM_rad_ps,absolute_sigma=True,bounds=([0,.1],[np.inf,4]),maxfev=100000);q['Gamma_model_rad_ps']=mod(q.k_inv_A,*p);q.to_csv(d/'LA_linewidth_fits_n003_n008.csv',index=False)
 pd.DataFrame([{'model':'Gamma0+A*k^alpha','Gamma0_rad_ps':G0,'Gamma0_source':str(a.gamma0_summary),'A_rad_ps_A_to_alpha':p[0],'A_SE':np.sqrt(c[0,0]),'alpha':p[1],'alpha_SE':np.sqrt(c[1,1]),'n_min':3,'n_max':8,'n_points':6}]).to_csv(d/'LA_Gamma0_plus_powerlaw_n003_n008.csv',index=False)
 fig,ax=plt.subplots(figsize=(6,4.2));ax.errorbar(q.k_inv_A,q.gamma_HWHM_rad_ps,yerr=q.gamma_fit_SEM_rad_ps,fmt='o',color='#2166ac',capsize=2);k=np.linspace(0,q.k_inv_A.max()*1.1,200);ax.plot(k,mod(k,*p),color='#b2182b');ax.axhline(G0,color='.3',ls='--');ax.set(xlabel=r'$k$ ($\AA^{-1}$)',ylabel=r'LA HWHM $\Gamma$ (rad ps$^{-1}$)');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(f/'LA_Gamma0_plus_powerlaw_n003_n008.png',dpi=300);fig.savefig(f/'LA_Gamma0_plus_powerlaw_n003_n008.pdf');plt.close(fig)
 (a.output/'metadata.json').write_text(json.dumps({'protocol':'10L 10ns 100fs Tdamp=1000ps','definition':'J_n=sum_O vz exp(i k_n z)','selection':'n=3..8','Gamma0_limit':'oxygen-current k0 ACF anchor; no full-water velocities in this source'},indent=2));(a.output/'FINISHED.txt').write_text('Tdamp1ns LA linewidth analysis finished successfully.\n')
if __name__=='__main__':main()
