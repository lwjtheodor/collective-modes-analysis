"""C88 implicit-CNT transverse S_JJ: radial and circumferential modes."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def read(src,nmax):
 t=[];th=[];rr=[];Ls=[]
 with src.open() as f:
  while True:
   m=f.readline()
   if not m:break
   if m.strip()!='ITEM: TIMESTEP':raise ValueError(m)
   step=int(f.readline());f.readline();na=int(f.readline());f.readline();b=[f.readline().split() for _ in range(3)];h=f.readline().split()[2:];c={v:i for i,v in enumerate(h)}
   if any(v not in c for v in ('type','x','y','z','vx','vy')):raise ValueError(h)
   a=np.fromstring(' '.join(f.readline() for _ in range(na)),sep=' ').reshape(na,len(h));a=a[a[:,c['type']].astype(int)==1];x,y,z,vx,vy=(a[:,c[q]] for q in ('x','y','z','vx','vy'));r=np.hypot(x,y)
   vt=(-y*vx+x*vy)/r;vr=(x*vx+y*vy)/r;L=float(b[2][1])-float(b[2][0]);k=2*np.pi*np.arange(1,nmax+1)/L;phase=np.exp(1j*np.outer(k,z));th.append(phase@vt);rr.append(phase@vr);t.append(step*.0005);Ls.append(L)
 return np.asarray(t),np.asarray(th).T,np.asarray(rr).T,np.mean(Ls)
def lor0(w,b,A,g):return b+A*g*g/(w*w+g*g)
def dho(w,b,A,w0,g):return b+A*g*w*w/((w*w-w0*w0)**2+(g*w)**2)
def main():
 p=argparse.ArgumentParser();p.add_argument('--dumps',nargs='+',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--nmax',type=int,default=8);a=p.parse_args();d=a.output/'derived_data';fg=a.output/'figures';d.mkdir(parents=True,exist_ok=True);fg.mkdir(exist_ok=True)
 TH=[];RR=[];Ls=[]
 for s in a.dumps:
  t,th,rr,L=read(s,a.nmax);dt=np.median(np.diff(t));win=np.hanning(len(t));TH.append(np.abs(np.fft.fft((th-th.mean(axis=1,keepdims=True))*win,axis=1))**2/np.sum(win**2));RR.append(np.abs(np.fft.fft((rr-rr.mean(axis=1,keepdims=True))*win,axis=1))**2/np.sum(win**2));Ls.append(L)
 w=np.fft.fftfreq(TH[0].shape[-1],dt)*2*np.pi;pos=w>=0;w=w[pos];TH=[x[...,pos] for x in TH];RR=[x[...,pos] for x in RR];nyq=w[-1];data=[];fits=[]
 for label,X in [('theta',np.asarray(TH)),('radial',np.asarray(RR))]:
  mean=X.mean(0);sem=X.std(0,ddof=1)/np.sqrt(len(X));fig,axs=plt.subplots(2,4,figsize=(13,6),squeeze=False)
  for n in range(1,a.nmax+1):
   y=mean[n-1];e=np.maximum(sem[n-1],np.nanmax(sem[n-1])*1e-5)
   if label=='theta':
    m=w<=.25;x=w[m];yy=y[m];ee=e[m];p0=(np.median(yy[-max(4,len(yy)//5):]),max(yy[0],1e-20),.02);bound=([0,0,w[1]*.5],[np.inf,np.inf,.25]);fn=lor0
   else:
    m=(w>=.25)&(w<=nyq-w[1]*3);xx=w[m];yy0=y[m];seed=xx[np.argmax(yy0)];half=max(1.0,.3*seed);m=(w>=max(.2,seed-half))&(w<=min(nyq-w[1]*2,seed+half));x=w[m];yy=y[m];ee=e[m];p0=(max(0,yy.min()),max(yy.max()-yy.min(),1e-20),seed,max(.2,.1*seed));bound=([0,0,x.min(),w[1]*.5],[np.inf,np.inf,x.max(),min(20,2*seed)]);fn=dho
   try:
    q,c=curve_fit(fn,x,yy,p0=p0,sigma=ee,absolute_sigma=True,bounds=bound,maxfev=100000);pred=fn(x,*q);r2=1-np.sum((yy-pred)**2)/np.sum((yy-yy.mean())**2);g=q[-1];gse=np.sqrt(c[-1,-1]);accept=(g>1.5*w[1]) and r2>.7 and not np.isclose(g,bound[1][-1])
   except RuntimeError:q=[np.nan]*len(p0);pred=np.full_like(x,np.nan);r2=np.nan;g=np.nan;gse=np.nan;accept=False
   fits.append({'branch':label,'n':n,'k_inv_A':2*np.pi*n/np.mean(Ls),'omega0_rad_ps':0 if label=='theta' else q[2],'gamma_rad_ps':g,'gamma_SE':gse,'fit_R2':r2,'accepted':accept,'frequency_bin_rad_ps':w[1],'nyquist_rad_ps':nyq})
   data.extend({'branch':label,'n':n,'omega_rad_ps':ww,'S_mean':ss,'S_replica_sem':eev} for ww,ss,eev in zip(w,y,e));ax=axs.ravel()[n-1];ax.semilogy(x,yy,color='#2166ac');ax.semilogy(x,pred,color='#b2182b');ax.set_title(f'{label}, n={n}');ax.grid(alpha=.2)
  fig.supxlabel(r'$\omega$ (rad ps$^{-1}$)');fig.supylabel(r'$S_{JJ}$');fig.tight_layout();fig.savefig(fg/f'implicitC88_{label}_SJJ_n001_n008.png',dpi=300);fig.savefig(fg/f'implicitC88_{label}_SJJ_n001_n008.pdf');plt.close(fig)
 fit=pd.DataFrame(fits);fit.to_csv(d/'implicitC88_transverse_linewidths_n001_n008.csv',index=False);pd.DataFrame(data).to_csv(d/'implicitC88_transverse_SJJ_n001_n008.csv',index=False)
 rows=[]
 for br in ('theta','radial'):
  q=fit[(fit.branch==br)&fit.accepted]
  if len(q)>=4:
   def mod(k,g0,A,al):return g0+A*k**al
   z,c=curve_fit(mod,q.k_inv_A,q.gamma_rad_ps,p0=(.01,20,1.5),sigma=q.gamma_SE,absolute_sigma=True,bounds=([0,0,.1],[np.inf,np.inf,4]),maxfev=100000);se=np.sqrt(np.diag(c));rows.append({'branch':br,'model':'Gamma0+A*k^alpha','Gamma0_rad_ps':z[0],'Gamma0_SE':se[0],'A':z[1],'alpha':z[2],'alpha_SE':se[2],'n_used':','.join(map(str,q.n.astype(int)))})
  else:rows.append({'branch':br,'model':'insufficient_or_nyquist_limited_modes','n_used':','.join(map(str,q.n.astype(int)))})
 pd.DataFrame(rows).to_csv(d/'implicitC88_transverse_Gamma0_powerlaw.csv',index=False);(a.output/'metadata.json').write_text(json.dumps({'definition_theta':'sum_O v_theta exp(ikz)','definition_radial':'sum_O v_r exp(ikz)','sampling':'100fs; Nyquist pi/dt','limit':'radial damping is invalid if resonance approaches Nyquist'},indent=2));(a.output/'FINISHED.txt').write_text('Implicit C88 transverse SJJ analysis finished successfully.\n')
if __name__=='__main__':main()
