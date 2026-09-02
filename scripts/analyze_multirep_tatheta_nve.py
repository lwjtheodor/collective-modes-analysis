"""Ensemble TA_theta spectra and linewidths from full-water NVE trajectories."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def lorentz(w, b, a, g): return b + a*g*g/(w*w + g*g)


def read_current(path, nmax, oxygen_type):
    ts, js, lengths = [], [], []
    with path.open() as f:
        while True:
            if not (line := f.readline()): break
            if line.strip() != 'ITEM: TIMESTEP': raise ValueError(f'bad dump marker {line!r}')
            step = int(f.readline()); f.readline(); nat = int(f.readline())
            f.readline(); bounds = np.array([[float(x) for x in f.readline().split()] for _ in range(3)])
            header = f.readline().split()[2:]; col = {v:i for i,v in enumerate(header)}
            need = ('type','x','y','z','vx','vy')
            if any(v not in col for v in need): raise ValueError(f'missing fields {header}')
            a = np.fromstring(' '.join(f.readline() for _ in range(nat)), sep=' ').reshape(nat,len(header))
            a = a[a[:,col['type']].astype(int) == oxygen_type]
            x,y,z,vx,vy = (a[:,col[v]] for v in ('x','y','z','vx','vy'))
            r = np.hypot(x,y)
            if np.any(r == 0): raise ValueError('oxygen on z axis')
            vx -= vx.mean(); vy -= vy.mean()
            vt = -y/r*vx + x/r*vy
            L = bounds[2,1]-bounds[2,0]
            k = 2*np.pi*np.arange(1,nmax+1)/L
            js.append(np.exp(1j*np.outer(k,z)) @ vt)
            ts.append(step*0.0005); lengths.append(L)
    return np.asarray(ts), np.asarray(js).T, float(np.mean(lengths))


def psd(q, dt, nperseg):
    if len(q) < nperseg: raise ValueError(f'{len(q)} frames < nperseg {nperseg}')
    q=q[:nperseg]-np.mean(q[:nperseg]); win=np.hanning(nperseg)
    f=np.fft.fftfreq(nperseg,dt); mask=f>=0
    return f[mask]*2*np.pi, np.abs(np.fft.fft(q*win))[mask]**2/np.sum(win**2)


def fit(w,y,wmax):
    m=(w>=0)&(w<=wmax); x,y=w[m],y[m]
    b=np.median(y[-max(5,len(y)//5):]); a=max(y[0]-b,1e-20)
    p,_=curve_fit(lorentz,x,y,p0=(b,a,0.01),bounds=([0,0,0.00005],[np.inf,np.inf,wmax]),maxfev=30000)
    pred=lorentz(x,*p); r2=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
    return p,r2,m


def main():
    p=argparse.ArgumentParser(); p.add_argument('--dumps',nargs='+',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--oxygen-type',type=int,required=True); p.add_argument('--nmax',type=int,default=20); p.add_argument('--nperseg',type=int,default=20000); p.add_argument('--omega-fit-max',type=float,default=.15); p.add_argument('--primary-nmax',type=int,default=6)
    a=p.parse_args(); dd=a.output/'derived_data'; fg=a.output/'figures'; dd.mkdir(parents=True,exist_ok=True); fg.mkdir(exist_ok=True)
    allpsd=[]; raw=[]; lengths=[]; dts=[]
    for ir,src in enumerate(a.dumps,1):
        t,j,L=read_current(src,a.nmax,a.oxygen_type); dt=float(np.median(np.diff(t))); lengths.append(L); dts.append(dt)
        this=[]
        for n in range(a.nmax):
            w,s=psd(j[n],dt,a.nperseg); this.append(s)
            raw.extend({'replica':ir,'n':n+1,'k_Ainv':2*np.pi*(n+1)/L,'omega_rad_ps':ww,'S':ss} for ww,ss in zip(w,s))
        allpsd.append(this)
    A=np.asarray(allpsd); mean=A.mean(axis=0); sem=A.std(axis=0,ddof=1)/np.sqrt(A.shape[0]); w=w
    rec=[]; ens=[]
    nc=5; nr=int(np.ceil(a.nmax/nc)); fig,axarr=plt.subplots(nr,nc,figsize=(14,2.75*nr),squeeze=False)
    for n,ax in enumerate(axarr.ravel(),1):
        if n>a.nmax: ax.set_visible(False); continue
        pp,r2,m=fit(w,mean[n-1],a.omega_fit_max); gs=[]
        for r in range(A.shape[0]):
            try: gs.append(fit(w,A[r,n-1],a.omega_fit_max)[0][2])
            except RuntimeError: pass
        gsem=np.std(gs,ddof=1)/np.sqrt(len(gs)) if len(gs)>1 else np.nan; k=2*np.pi*n/np.mean(lengths)
        rec.append({'n':n,'k_Ainv':k,'gamma_rad_ps':pp[2],'gamma_replica_sem_rad_ps':gsem,'fit_R2':r2,'n_replicas':A.shape[0],'frequency_resolution_rad_ps':w[1]})
        ens.extend({'n':n,'k_Ainv':k,'omega_rad_ps':ww,'S_mean':ss,'S_replica_sem':ee} for ww,ss,ee in zip(w,mean[n-1],sem[n-1]))
        ax.semilogy(w[m],mean[n-1][m],color='#2166ac'); ax.semilogy(w[m],lorentz(w[m],*pp),color='#b2182b'); ax.set_title(f'n={n}'); ax.grid(alpha=.2)
    fig.supxlabel(r'$\omega$ (rad ps$^{-1}$)'); fig.supylabel(r'$S_{J_\theta J_\theta}$ (arb.)'); fig.tight_layout()
    fig.savefig(fg/'TAtheta_NVE_spectra.png',dpi=300); fig.savefig(fg/'TAtheta_NVE_spectra.pdf'); plt.close(fig)
    q=pd.DataFrame(rec); q['used_primary_k_fit']=(q.n<=a.primary_nmax)&(q.fit_R2>=.55); z=q[q.used_primary_k_fit]; X=np.c_[np.ones(len(z)),z.k_Ainv.values**2]; b=np.linalg.lstsq(X,z.gamma_rad_ps.values,rcond=None)[0]
    pd.DataFrame([{'model':'Gamma0_plus_Dk2','Gamma0_rad_ps':b[0],'D_rad_ps_A2':b[1],'n_used':len(z)}]).to_csv(dd/'TAtheta_NVE_k_dependence_fit.csv',index=False)
    fig,ax=plt.subplots(figsize=(6.2,4.5)); ax.errorbar(q.k_Ainv,q.gamma_rad_ps,yerr=q.gamma_replica_sem_rad_ps,fmt='o',color='#2166ac',capsize=3); kk=np.linspace(0,q.k_Ainv.max(),200); ax.plot(kk,b[0]+b[1]*kk**2,color='#b2182b'); ax.set(xlabel=r'$k$ (A$^{-1}$)',ylabel=r'$\Gamma_\theta$ (rad ps$^{-1}$)'); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(fg/'TAtheta_NVE_linewidth.png',dpi=300); fig.savefig(fg/'TAtheta_NVE_linewidth.pdf'); plt.close(fig)
    pd.DataFrame(raw).to_csv(dd/'TAtheta_NVE_spectra_per_replica.csv',index=False); pd.DataFrame(ens).to_csv(dd/'TAtheta_NVE_spectra_ensemble.csv',index=False); q.to_csv(dd/'TAtheta_NVE_lorentzian_fits.csv',index=False)
    (a.output/'metadata.json').write_text(json.dumps({'dumps':[str(x) for x in a.dumps],'oxygen_type':a.oxygen_type,'nmax':a.nmax,'nperseg_frames':a.nperseg,'n_replicas':len(a.dumps),'dt_ps':dts,'Lz_A':lengths},indent=2)); (a.output/'FINISHED.txt').write_text('NVE multi-replica TA_theta analysis finished successfully.\n')
if __name__=='__main__': main()
