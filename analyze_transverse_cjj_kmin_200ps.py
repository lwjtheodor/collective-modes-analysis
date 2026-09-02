#!/usr/bin/env python3
"""Stream four 1 ns full-component CNT-water dumps into 200 ps cylindrical CJJ data."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def current_series(paths, modes, need_z):
    signals, last_step, box, center, nw = [], None, None, None, None
    for path in paths:
        with Path(path).open(errors="replace") as fh:
            while True:
                marker = fh.readline()
                if marker == "":
                    break
                if marker != "ITEM: TIMESTEP\n":
                    raise ValueError("invalid dump marker")
                step = int(fh.readline())
                if fh.readline() != "ITEM: NUMBER OF ATOMS\n": raise ValueError("missing count")
                natom = int(fh.readline())
                if not fh.readline().startswith("ITEM: BOX BOUNDS"): raise ValueError("missing bounds")
                bounds = np.array([list(map(float, fh.readline().split()[:2])) for _ in range(3)])
                names = fh.readline().split()[2:]; col = {x:i for i,x in enumerate(names)}
                xyz=[]; vel=[]
                for _ in range(natom):
                    row=fh.readline().split()
                    if int(float(row[col["type"]])) == 3:
                        xyz.append([float(row[col[x]]) for x in ("x","y","z")])
                        vel.append([float(row[col[x]]) for x in ("vx","vy","vz")])
                if last_step is not None and step <= last_step:  # restart-frame duplicate
                    continue
                xyz=np.asarray(xyz); vel=np.asarray(vel); frame_box=bounds[:,1]-bounds[:,0]
                if box is None:
                    box=frame_box; center=bounds.mean(axis=1); nw=len(xyz)
                if len(xyz) != nw or not np.allclose(frame_box,box): raise ValueError("nonstationary box/count")
                vel -= vel.mean(axis=0, keepdims=True)
                dxy=xyz[:,:2]-center[:2]; r=np.hypot(dxy[:,0],dxy[:,1])
                er=dxy/np.maximum(r[:,None],1e-12); et=np.c_[-er[:,1],er[:,0]]
                theta=np.sum(vel[:,:2]*et,axis=1)
                row=[]
                for n in modes:
                    phase=np.exp(1j*(2*np.pi*n/box[2])*xyz[:,2])
                    row.append(np.sum(theta*phase))
                if need_z:
                    phase=np.exp(1j*(2*np.pi/box[2])*xyz[:,2])
                    row.append(np.sum(vel[:,2]*phase))
                signals.append(row); last_step=step
    return np.asarray(signals), float(box[2]), nw


def corr(x, y, maxlag):
    x=x-x.mean(); y=y-y.mean(); n=len(x)
    out=np.fft.ifft(np.fft.fft(x,2*n)*np.conj(np.fft.fft(y,2*n)))[:maxlag+1]
    return out/np.arange(n,n-maxlag-1,-1)


def write_csv(path, header, rows):
    with Path(path).open("w", newline="") as fh:
        w=csv.writer(fh); w.writerow(header); w.writerows(rows)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",required=True,type=Path); ap.add_argument("--out",required=True,type=Path); args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True); dt=0.01; maxlag=20000; all_meta={}
    cases=["7_7","8_8","9_9","17_0"]
    theta_curves={}
    for ch in cases:
        case=f"TRANSVERSE_WEAKNH_NOMOM_{ch}_L5_rep1"; d=args.root/case
        paths=[d/f"{case}_water_full_10fs_200ps.dump",d/f"{case}_water_full_10fs_extend800ps.dump"]
        modes=[1,2,3,4,5] if ch=="8_8" else [1]
        sig,lz,nw=current_series(paths,modes,ch=="8_8")
        if len(sig)!=100001: raise ValueError(f"{ch}: expected 100001 frames, got {len(sig)}")
        t=np.arange(maxlag+1)*dt
        ctheta=[]
        for i,n in enumerate(modes):
            c=corr(sig[:,i],sig[:,i],maxlag).real; ctheta.append(c/c[0])
        theta_curves[ch]=ctheta[0]
        all_meta[ch]={"n_frames":int(len(sig)),"n_water":int(nw),"lz_A":lz,"kmin_inv_A":2*np.pi/lz,"dt_ps":dt,"maxlag_ps":200.0,"dump_files":[str(x) for x in paths]}
        if ch=="8_8":
            z=sig[:,-1]; th=sig[:,0]; czz=corr(z,z,maxlag).real; ctt=corr(th,th,maxlag).real; czt=corr(z,th,maxlag).real
            write_csv(args.out/"8_8_k1_axial_theta_cross_cjj.csv",["time_ps","Czz_norm","Ctheta_theta_norm","Cz_theta_norm","n_origins"],zip(t,czz/czz[0],ctt/ctt[0],czt/np.sqrt(czz[0]*ctt[0]),np.arange(len(sig),len(sig)-maxlag-1,-1)))
            write_csv(args.out/"8_8_theta_k1to5_cjj.csv",["time_ps"]+[f"Ctheta_theta_k{n}_norm" for n in modes]+["n_origins"],zip(t,*ctheta,np.arange(len(sig),len(sig)-maxlag-1,-1)))
    write_csv(args.out/"theta_k1_crosschirality_cjj.csv",["time_ps"]+[f"Ctheta_theta_{ch}_k1_norm" for ch in cases]+["n_origins"],zip(np.arange(maxlag+1)*dt,*[theta_curves[x] for x in cases],np.arange(100001,100001-maxlag-1,-1)))
    (args.out/"metadata.json").write_text(json.dumps({"definition":"J_theta(k,t)=sum_i[(v_i-<v>_O) dot e_theta,i] exp(i k z_i); C_ab=Re<delta J_a(0)delta J_b*(t)>, normalized by sqrt(C_aa(0)C_bb(0)).","cases":all_meta,"statistics":"one continuous 1 ns trajectory per chirality; all time origins; no replica SEM"},indent=2))
    (args.out/"FINISHED.txt").write_text("all CJJ CSVs completed\n")

if __name__ == "__main__": main()
