#!/usr/bin/env python3
"""Calculate normalized total n=1 axial current ACF from low-frequency dumps.

The instantaneous water centre-of-mass axial velocity is removed before the
current is constructed, matching the weak-NH cross-chirality C_vJ definition.
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np


def read_oxygen(path: Path):
    zframes, vframes, mol_ref, steps, lzs = [], [], None, [], []
    with path.open(errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line: break
            if not line.startswith("ITEM: TIMESTEP"): continue
            step = int(fh.readline())
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"): raise ValueError("missing atom count")
            natom = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"): raise ValueError("missing box bounds")
            box = [list(map(float, fh.readline().split())) for _ in range(3)]
            lz = box[2][1] - box[2][0]
            header = fh.readline().split()[2:]; idx = {v:i for i,v in enumerate(header)}
            req = {"mol", "type", "z", "vz"}
            if not req.issubset(idx): raise ValueError(f"required fields missing: {req-set(idx)}")
            unwrapped = "iz" in idx
            frame = {}
            for _ in range(natom):
                x = fh.readline().split()
                if int(float(x[idx["type"]])) == 3:
                    mol = int(float(x[idx["mol"]])); z = float(x[idx["z"]])
                    if unwrapped: z += int(float(x[idx["iz"]])) * lz
                    frame[mol] = (z, float(x[idx["vz"]]))
            mols = np.asarray(sorted(frame), dtype=np.int64)
            if mol_ref is None: mol_ref = mols
            elif not np.array_equal(mol_ref, mols): raise ValueError("oxygen molecule IDs changed")
            zframes.append([frame[int(m)][0] for m in mol_ref]); vframes.append([frame[int(m)][1] for m in mol_ref])
            steps.append(step); lzs.append(lz)
    return np.asarray(zframes), np.asarray(vframes), np.asarray(steps), float(np.median(lzs))


def acf_complex(q, maxlag):
    frames = len(q); nfft = 1 << (2*frames - 1).bit_length()
    c = np.fft.ifft(np.fft.fft(q, nfft) * np.conjugate(np.fft.fft(q, nfft)))[:maxlag+1].real
    return c / np.arange(frames, frames-maxlag-1, -1, dtype=float)


def lobe(time, y):
    down = np.flatnonzero((y[:-1] >= 0) & (y[1:] < 0))
    if not len(down): return None
    i0 = down[0]; imin = i0+1+np.argmin(y[i0+1:])
    up = np.flatnonzero((y[imin:-1] <= 0) & (y[imin+1:] > 0))
    if not len(up): return None
    i1 = imin+up[0]
    def cross(i): return time[i] - y[i]*(time[i+1]-time[i])/(y[i+1]-y[i])
    t0,t1=cross(i0),cross(i1); tt=np.r_[t0,time[i0+1:i1+1],t1]; yy=np.r_[0.,y[i0+1:i1+1],0.]
    return {"negative_area_ps":float(-np.trapezoid(yy,tt)),"depth":float(-y[imin]),"t_min_ps":float(time[imin]),"t_start_ps":float(t0),"t_end_ps":float(t1)}


def main():
    p=argparse.ArgumentParser(); p.add_argument("--dump",required=True,type=Path); p.add_argument("--dt-ps",required=True,type=float)
    p.add_argument("--max-lag-ps",type=float,default=100.); p.add_argument("--out",required=True,type=Path); p.add_argument("--case-id",required=True); a=p.parse_args()
    z,v,steps,lz=read_oxygen(a.dump); observed=float(np.median(np.diff(steps))*0.0005)
    if not np.isclose(observed,a.dt_ps,rtol=1e-5,atol=1e-8): raise ValueError(f"observed cadence {observed} ps != declared {a.dt_ps} ps")
    v-=v.mean(axis=1,keepdims=True); k=2*math.pi/lz; q=np.sum(v*np.exp(1j*k*z),axis=1)
    maxlag=int(round(a.max_lag_ps/a.dt_ps));
    if len(q)<=maxlag: raise ValueError("trajectory shorter than requested lag")
    curve=acf_complex(q,maxlag); curve/=curve[0]; time=np.arange(maxlag+1)*a.dt_ps
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["case_id","lag_ps","n","k_inv_A","C_vJ_total"]); w.writeheader()
        for t,y in zip(time,curve): w.writerow({"case_id":a.case_id,"lag_ps":t,"n":1,"k_inv_A":k,"C_vJ_total":y})
    meta={"case_id":a.case_id,"n_frames":int(len(q)),"n_water":int(z.shape[1]),"lz_A":lz,"k_inv_A":k,"dt_ps":a.dt_ps,"max_lag_ps":a.max_lag_ps,"water_com_velocity_removed":True,"first_negative_lobe":lobe(time,curve)}
    a.out.with_suffix(".json").write_text(json.dumps(meta,indent=2),encoding="utf-8")

if __name__=="__main__": main()
