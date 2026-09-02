#!/usr/bin/env python3
"""Sampled axial total/self/distinct oxygen ISF from full-water dumps."""
from __future__ import annotations
import argparse, json
from collections import defaultdict
from pathlib import Path
import numpy as np

def read_sampled_z(path, oxygen_type, stride):
    frames=[]; ids0=None; lz=None; frame=0
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        while True:
            if not fh.readline() or not fh.readline(): break
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"): raise ValueError(path)
            nat=int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"): raise ValueError(path)
            b=[fh.readline().split() for _ in range(3)]; here=float(b[2][1])-float(b[2][0])
            lz = here if lz is None else lz
            header=fh.readline().split()[2:]; idx={x:header.index(x) for x in header}
            if frame % stride:
                for _ in range(nat): fh.readline()
                frame += 1; continue
            raw=np.fromstring(" ".join(fh.readline() for _ in range(nat)),sep=" ",dtype=np.float64).reshape(nat,len(header))
            keep=raw[:,idx["type"]].astype(np.int16)==oxygen_type
            atom_ids=raw[keep,idx["id"]].astype(np.int64); order=np.argsort(atom_ids)
            atom_ids=atom_ids[order]
            if ids0 is None: ids0=atom_ids
            elif not np.array_equal(ids0,atom_ids): raise ValueError("oxygen identities changed")
            z=raw[keep,idx["z"]][order]
            if "iz" in idx: z=z+raw[keep,idx["iz"]][order]*lz
            frames.append(z.astype(np.float32)); frame += 1
    return np.asarray(frames),float(lz),ids0

def acf(series, maxlag):
    n=series.shape[0]; size=1 << (2*n-1).bit_length()
    x=np.fft.fft(series,size,axis=0)
    return np.fft.ifft(x*np.conjugate(x),axis=0).real[:maxlag+1]

def components(z,k,maxlag):
    nt,na=z.shape; maxlag=min(maxlag,nt-1); origins=nt-np.arange(maxlag+1)
    # Reciprocal-lattice phase: wrapped z avoids accumulated-image roundoff.
    phase=np.exp(1j*k*np.remainder(z, 2*np.pi/k)).astype(np.complex64)
    rho=phase.sum(axis=1)
    total=acf(rho[:,None],maxlag)[:,0]/(na*origins)
    selfsum=np.zeros(maxlag+1)
    for start in range(0,na,48): selfsum += acf(phase[:,start:start+48],maxlag).sum(axis=1)
    selfterm=selfsum/(na*origins)
    return total,selfterm,total-selfterm

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("manifest",type=Path); ap.add_argument("--output-dir",type=Path,required=True)
    ap.add_argument("--dt-ps",type=float,required=True); ap.add_argument("--stride",type=int,required=True); ap.add_argument("--max-lag-ps",type=float,required=True)
    ap.add_argument("--n-list",default="1,2"); ap.add_argument("--oxygen-type",type=int,default=3); a=ap.parse_args()
    modes=np.array([int(x) for x in a.n_list.split(",")]); sources=json.loads(a.manifest.read_text()); out=a.output_dir; out.mkdir(parents=True,exist_ok=True); per=out/"per_replica"; per.mkdir(exist_ok=True)
    grouped=defaultdict(list); meta=[]; dt=a.dt_ps*a.stride; maxlag=round(a.max_lag_ps/dt)
    for src in sources:
        z,lz,ids=read_sampled_z(src["path"],a.oxygen_type,a.stride); record={**src,"frames_used":int(z.shape[0]),"lz_A":lz,"n_oxygen":int(z.shape[1]),"dt_ps":dt,"modes":{}}
        for n in modes:
            k=2*np.pi*n/lz; total,selfterm,distinct=components(z,k,maxlag); time=np.arange(total.size)*dt
            p=per/f"{src['label']}_k{n}_components.npz"; np.savez_compressed(p,time_ps=time,F_total=total,F_self=selfterm,F_distinct=distinct,k_inv_A=k,n=n,n_time_origins=z.shape[0]-np.arange(total.size),source_path=src["path"])
            grouped[n].append((total,selfterm,distinct)); record["modes"][str(n)]=str(p.relative_to(out))
        meta.append(record); print(f"completed {src['label']}: {z.shape[0]} sampled frames",flush=True)
    for n,curves in grouped.items():
        arrays=np.asarray(curves); payload={"time_ps":time,"n":n,"k_inv_A":2*np.pi*n/meta[0]['lz_A'],"n_replicas":len(curves),"n_time_origins_per_replica":meta[0]['frames_used']-np.arange(time.size)}
        for j,name in enumerate(("total","self","distinct")):
            payload[f"F_{name}_mean"]=arrays[:,j].mean(0); payload[f"F_{name}_sem"]=arrays[:,j].std(0,ddof=1)/np.sqrt(len(curves))
        np.savez_compressed(out/f"ISF_88_L10_k{n}_components_mean_sem.npz",**payload)
    (out/"source_manifest_resolved.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
if __name__=="__main__": main()
