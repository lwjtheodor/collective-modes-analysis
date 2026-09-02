#!/usr/bin/env python3
"""LA CJJ/SJJ for the completed non-continuation (8,8) 10L 8-rep archive."""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

SOURCE = Path(r"F:\ccfep_gcmc_archive_20260814\stage_vacf_tail_8_8_L2L10_8rep_weaknh_zvz_20260812")
OUT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes\results\collective_mode_response\88_10L_LA_CJJ_Skw_100fs_10ns_8rep_n001_n010\2026-08-24")
NREP, NMAX, NFRAME, DT, MAXLAG, NPERSEG = 8, 10, 100001, 0.1, 1000.0, 16384

plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"], "font.size":7, "axes.spines.right":False, "axes.spines.top":False, "xtick.direction":"out", "ytick.direction":"out", "svg.fonttype":"none", "pdf.fonttype":42})

def read_modal(path):
    series=np.empty((NFRAME,NMAX),complex); steps=[]; lz=None; frame=0
    with path.open(encoding="utf-8",errors="replace") as fh:
        while True:
            if not fh.readline(): break
            step=int(fh.readline()); fh.readline(); natom=int(fh.readline()); fh.readline()
            b=np.array([[float(x) for x in fh.readline().split()] for _ in range(3)])
            header=fh.readline().split()[2:]; col={x:i for i,x in enumerate(header)}
            if set(("id","z","vz"))-set(col): raise ValueError(f"missing id/z/vz: {path}")
            z=[]; vz=[]
            for _ in range(natom):
                row=fh.readline().split(); z.append(float(row[col["z"]])); vz.append(float(row[col["vz"]]))
            z=np.asarray(z); vz=np.asarray(vz)-np.mean(vz); lz=float(b[2,1]-b[2,0])
            k=2*np.pi*np.arange(1,NMAX+1)/lz
            series[frame]=np.sum(vz[:,None]*np.exp(1j*z[:,None]*k),axis=0); steps.append(step); frame+=1
    if frame!=NFRAME: raise ValueError(f"{path}: expected {NFRAME}, got {frame}")
    if not np.all(np.diff(steps)==200): raise ValueError(f"{path}: nonuniform frame cadence")
    return series,lz

def acf(series):
    x=series-series.mean(axis=0,keepdims=True); n=len(x); fft=np.fft.fft(x,n=2*n,axis=0)
    ac=np.fft.ifft(fft*np.conj(fft),axis=0).real[:int(MAXLAG/DT)+1]
    return ac/np.arange(n,n-len(ac),-1)[:,None]

def psd(series):
    starts=np.arange(0,len(series)-NPERSEG+1,NPERSEG//2); win=np.hanning(NPERSEG)[:,None]; acc=np.zeros((NPERSEG,NMAX))
    for s in starts:
        x=series[s:s+NPERSEG]; x=x-x.mean(axis=0,keepdims=True); ff=np.fft.fft(x*win,axis=0); acc+=np.abs(ff)**2/np.sum(win[:,0]**2)
    acc/=len(starts); ix=np.arange(1,NPERSEG//2+1); return ix/(NPERSEG*DT),0.5*(acc[ix]+acc[-ix]),len(starts)

def write_tables(cjj, freqs, spectra, data):
    time=np.arange(cjj.shape[1])*DT; k=2*np.pi*np.arange(1,NMAX+1)/1008.39998
    for name, arr in (("CJJ_all_modes_per_replica.csv",cjj),):
        with (data/name).open("w",newline="") as fh:
            w=csv.writer(fh); w.writerow(["replica","branch","n","k_inv_A","time_ps","CJJ_raw_A2_fs2"])
            for r in range(NREP):
                for n in range(NMAX):
                    w.writerows((r+1,"LA",n+1,k[n],t,v) for t,v in zip(time,arr[r,:,n]))
    mean=cjj.mean(0); sem=cjj.std(0,ddof=1)/np.sqrt(NREP)
    with (data/"CJJ_all_modes_ensemble_mean_sem.csv").open("w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["branch","n","k_inv_A","time_ps","CJJ_mean_A2_fs2","CJJ_replica_SEM_A2_fs2"])
        for n in range(NMAX): w.writerows(("LA",n+1,k[n],t,v,e) for t,v,e in zip(time,mean[:,n],sem[:,n]))
    for name, arr in (("current_spectra_all_modes_per_replica.csv",spectra),):
        with (data/name).open("w",newline="") as fh:
            w=csv.writer(fh); w.writerow(["replica","branch","n","k_inv_A","frequency_ps_inv","PSD_arbitrary"])
            for r in range(NREP):
                for n in range(NMAX): w.writerows((r+1,"LA",n+1,k[n],f,v) for f,v in zip(freqs,arr[r,:,n]))
    mean=spectra.mean(0); sem=spectra.std(0,ddof=1)/np.sqrt(NREP)
    with (data/"current_spectra_all_modes_ensemble_mean_sem.csv").open("w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["branch","n","k_inv_A","frequency_ps_inv","omega_rad_ps","PSD_mean_arbitrary","PSD_replica_SEM_arbitrary"])
        for n in range(NMAX): w.writerows(("LA",n+1,k[n],f,2*np.pi*f,v,e) for f,v,e in zip(freqs,mean[:,n],sem[:,n]))
    return time,k,mean,sem

def save(fig,path):
    fig.savefig(path.with_suffix(".png"),dpi=600,bbox_inches="tight"); fig.savefig(path.with_suffix(".pdf"),bbox_inches="tight"); fig.savefig(path.with_suffix(".svg"),bbox_inches="tight"); fig.savefig(path.with_suffix(".tiff"),dpi=600,bbox_inches="tight")

def figures(time,k,cmean,psd,fig):
    f=plt.figure(figsize=(7.1,4.6)); axes=f.subplots(2,5,sharex=True)
    for n,ax in enumerate(axes.flat): ax.plot(time,cmean[:,n],lw=.7,color="#1769aa"); ax.axhline(0,color=".3",lw=.45); ax.set_xlim(0,500); ax.set_title(rf"$n={n+1}$, $k={k[n]:.3f}$",fontsize=6.5); ax.tick_params(labelsize=6)
    f.supxlabel(r"$t$ (ps)",fontsize=8); f.supylabel(r"$C_{J_zJ_z}(k,t)$ (A$^2$ fs$^{-2}$)",fontsize=8); f.subplots_adjust(.10,.10,.985,.92,.26,.30); save(f,fig/"LA_CJJ_n001_n010_ensemble_mean"); plt.close(f)
    w=2*np.pi*np.arange(1,psd.shape[0]+1)/(NPERSEG*DT); f=plt.figure(figsize=(7.1,4.6)); axes=f.subplots(2,5,sharex=True)
    for n,ax in enumerate(axes.flat):
        y=np.maximum(psd[:,n],np.max(psd[:,n])*1e-8); ax.semilogy(np.r_[-w[::-1],w],np.r_[y[::-1],y],lw=.7,color="#d55e00"); ax.axvline(0,color=".3",lw=.45); ax.set_xlim(-1.5,1.5); ax.set_title(rf"$n={n+1}$",fontsize=6.5); ax.tick_params(labelsize=6)
    f.supxlabel(r"$\omega$ (rad ps$^{-1}$)",fontsize=8); f.supylabel(r"$S_{J_zJ_z}(k,\omega)$ (arb.)",fontsize=8); f.subplots_adjust(.10,.10,.985,.92,.26,.30); save(f,fig/"LA_signed_semilog_Skw_n001_n010"); plt.close(f)

def main():
    data=OUT/"derived_data"; fig=OUT/"figures"; data.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    c=[]; p=[]; freq=None; lz=[]; seg=None
    for r in range(1,NREP+1):
        path=SOURCE/f"10L_rep{r}"/f"VACF88_10L_tail_zvz_rep{r}_oxygen_id_z_vz_100fs.dump"; print(f"reading rep{r}",flush=True); series,L=read_modal(path); c.append(acf(series)); fr,sp,seg=psd(series); p.append(sp); freq=fr; lz.append(L)
    c=np.asarray(c); p=np.asarray(p); time,k,mean,sem=write_tables(c,freq,p,data); figures(time,k,c.mean(0),p.mean(0),fig)
    meta={"system":"(8,8) CNT-confined water, 10L; completed non-continuation 100 fs/10 ns, 8 replicas","source":str(SOURCE),"fields":["id","z","vz"],"branch":"LA only","n_range":[1,NMAX],"frames_per_replica":NFRAME,"dt_ps":DT,"duration_ps":10000.,"welch_nperseg_frames":NPERSEG,"welch_segment_ps":NPERSEG*DT,"welch_overlap":0.5,"welch_segments_per_replica":int(seg),"omega_resolution_rad_ps":float(2*np.pi/(NPERSEG*DT)),"cjj_maxlag_ps":MAXLAG,"Lz_A":lz}
    (OUT/"metadata.json").write_text(json.dumps(meta,indent=2)+"\n"); (OUT/"FINISHED.txt").write_text("10L 100fs 10ns 8rep LA CJJ/SJJ analysis finished successfully.\n")
if __name__=="__main__": main()
