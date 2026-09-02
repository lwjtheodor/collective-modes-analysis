#!/usr/bin/env python3
"""Full-trajectory, same-kz axial/transverse current-mode screen for CNT water."""
import argparse, csv, json
from pathlib import Path
import numpy as np


def read_modal_series(path, target_k):
    """Stream a text dump into complex currents without retaining frames."""
    signals, steps = [], []
    box = box_center = k = None
    n_water = None
    with path.open(errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line: break
            if line != "ITEM: TIMESTEP\n": raise ValueError("unexpected dump marker: %r" % line)
            step = int(fh.readline())
            if fh.readline() != "ITEM: NUMBER OF ATOMS\n": raise ValueError("missing atom count")
            natom = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"): raise ValueError("missing bounds")
            bounds = np.array([list(map(float, fh.readline().split()[:2])) for _ in range(3)])
            header = fh.readline().split()[2:]; col = {name:i for i,name in enumerate(header)}
            required = {"type","x","y","z","vx","vy","vz"}
            if not required.issubset(col): raise ValueError("missing columns: %s" % (required-set(col)))
            xyz_rows=[]; velocity_rows=[]
            for _ in range(natom):
                values = fh.readline().split()
                if int(float(values[col["type"]])) != 3: continue
                xyz_rows.append([float(values[col[a]]) for a in ("x","y","z")])
                velocity_rows.append([float(values[col[a]]) for a in ("vx","vy","vz")])
            if not xyz_rows: raise ValueError("no oxygen atoms (type 3)")
            frame_box = bounds[:,1]-bounds[:,0]
            if box is None:
                box=frame_box; box_center=bounds.mean(axis=1)
                n=max(1,int(round(target_k*box[2]/(2*np.pi)))); k=2*np.pi*n/box[2]
                n_water=len(xyz_rows)
            elif len(xyz_rows) != n_water or not np.allclose(frame_box,box):
                raise ValueError("water count or box changed across frames")
            xyz=np.asarray(xyz_rows); vel=np.asarray(velocity_rows)
            vel-=vel.mean(axis=0,keepdims=True)
            dxy=xyz[:,:2]-box_center[:2]; radius=np.hypot(dxy[:,0],dxy[:,1])
            er=dxy/np.maximum(radius[:,None],1e-12); et=np.c_[-er[:,1],er[:,0]]
            phase=np.exp(1j*k*xyz[:,2])
            signals.append([np.sum(vel[:,2]*phase),np.sum(np.sum(vel[:,:2]*er,axis=1)*phase),np.sum(np.sum(vel[:,:2]*et,axis=1)*phase)])
            steps.append(step)
    return np.asarray(signals), np.asarray(steps), box, box_center, k, n_water


def welch(x, y, dt, nperseg):
    win=np.hanning(nperseg); norm=np.sum(win**2); spectra=[]
    for start in range(0,len(x)-nperseg+1,nperseg//2):
        a=(x[start:start+nperseg]-np.mean(x[start:start+nperseg]))*win
        b=(y[start:start+nperseg]-np.mean(y[start:start+nperseg]))*win
        spectra.append((np.fft.fft(a),np.fft.fft(b)))
    if len(spectra)<2: raise ValueError("need at least two Welch segments")
    xa=np.array([s[0] for s in spectra]); ya=np.array([s[1] for s in spectra])
    pxx=np.mean(np.abs(xa)**2,axis=0)/norm; pyy=np.mean(np.abs(ya)**2,axis=0)/norm
    pxy=np.mean(xa*np.conj(ya),axis=0)/norm; keep=slice(0,nperseg//2+1)
    return np.fft.fftfreq(nperseg,dt)[keep],pxx[keep],pyy[keep],pxy[keep],len(spectra)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--dump",required=True,type=Path); ap.add_argument("--dt-ps",required=True,type=float)
    ap.add_argument("--target-k",type=float,default=0.0623084632); ap.add_argument("--outdir",required=True,type=Path); ap.add_argument("--case-id",required=True); ap.add_argument("--welch-nperseg",type=int,default=8192)
    args=ap.parse_args(); sig,steps,box,box_center,k,n_water=read_modal_series(args.dump,args.target_k)
    observed=np.median(np.diff(steps))*0.0005
    if not np.isclose(observed,args.dt_ps): raise ValueError("cadence %s != declared %s" % (observed,args.dt_ps))
    lz=float(box[2]); n=int(round(k*lz/(2*np.pi))); nperseg=min(args.welch_nperseg,2**int(np.floor(np.log2(len(sig)//2))))
    freq,pz,pr,czr,nseg=welch(sig[:,0],sig[:,1],args.dt_ps,nperseg); _,_,pt,czt,_=welch(sig[:,0],sig[:,2],args.dt_ps,nperseg)
    coh_r=np.abs(czr)**2/np.maximum(pz*pr,1e-300); coh_t=np.abs(czt)**2/np.maximum(pz*pt,1e-300)
    valid=(freq>=0.01)&(freq<=1.5); peak=np.flatnonzero(valid)[np.argmax(pz[valid])]
    args.outdir.mkdir(parents=True,exist_ok=True)
    with (args.outdir/"kz_transverse_spectra.csv").open("w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["case_id","freq_ps_inv","P_z","P_r","P_theta","coh_zr","coh_ztheta"])
        for row in zip(freq,pz,pr,pt,coh_r,coh_t): w.writerow([args.case_id]+[float(x) for x in row])
    summary={"case_id":args.case_id,"n_frames":len(sig),"n_water":n_water,"dt_ps":args.dt_ps,"cylindrical_axis":"fixed_box_center_as_CNT_axis","axis_xy_A":[float(x) for x in box_center[:2]],"lz_A":lz,"n":n,"kz_inv_A":k,"target_k_inv_A":args.target_k,"welch_nperseg":nperseg,"welch_segments":nseg,"axial_peak_freq_ps_inv":float(freq[peak]),"radial_to_axial_power_at_peak":float(pr[peak]/pz[peak]),"theta_to_axial_power_at_peak":float(pt[peak]/pz[peak]),"coh_zr_at_peak":float(coh_r[peak]),"coh_ztheta_at_peak":float(coh_t[peak])}
    (args.outdir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__": main()
