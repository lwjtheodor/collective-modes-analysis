"""Direct axial MSD and uniform-log one-decade alpha estimator for 1 ps dumps."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

def read_z(path):
    frames=[]; mol_ref=None
    with Path(path).open(errors="replace") as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.startswith("ITEM: TIMESTEP"): continue
            f.readline(); f.readline(); nat=int(f.readline())
            f.readline(); bounds=[[float(x) for x in f.readline().split()[:2]] for _ in range(3)]
            lz=bounds[2][1]-bounds[2][0]
            cols=f.readline().split()[2:]; ix={x:i for i,x in enumerate(cols)}
            if "z" not in ix: raise ValueError("dump needs z")
            identity="mol" if "mol" in ix else "id"
            if identity not in ix: raise ValueError("dump needs mol or id")
            row={}
            for _ in range(nat):
                v=f.readline().split()
                if "type" not in ix or int(float(v[ix["type"]])) == 3:
                    m=int(float(v[ix[identity]])); z=float(v[ix["z"]])
                    if "iz" in ix: z += int(float(v[ix["iz"]]))*lz
                    row[m]=z
            mol=np.asarray(sorted(row),dtype=np.int64)
            if mol_ref is None: mol_ref=mol
            elif not np.array_equal(mol,mol_ref): raise ValueError("oxygen IDs changed")
            frames.append([row[int(m)] for m in mol_ref])
    return np.asarray(frames,dtype=np.float64)

def xcorr_sum(a):
    n=a.shape[0]; size=1<<(2*n-1).bit_length()
    fa=np.fft.fft(a,n=size,axis=0)
    return np.real(np.fft.ifft(fa*np.conjugate(fa),axis=0)[:n].sum(axis=1))

def msd(z,maxlag):
    n,m=z.shape; corr=xcorr_sum(z)[:maxlag+1]; sq=np.sum(z*z,axis=1)
    pref=np.r_[0.,np.cumsum(sq)]; lag=np.arange(maxlag+1); count=n-lag
    return ((pref[n]-pref[lag])+pref[n-lag]-2*corr)/(count*m)

def local_decade_slope(t,m,grid):
    use=(t>=1.0)&(m>0)
    lx=np.log10(grid); ly=np.interp(lx,np.log10(t[use]),np.log10(np.maximum(m[use],1e-12)))
    out=np.full(len(grid),np.nan)
    for i,x0 in enumerate(lx):
        take=np.abs(lx-x0)<=0.5+1e-12
        if x0-0.5>=lx[0] and x0+0.5<=lx[-1] and take.sum()>=3:
            out[i]=np.polyfit(lx[take],ly[take],1)[0]
    return out

def main():
    p=argparse.ArgumentParser(); p.add_argument("--dump",required=True); p.add_argument("--out",required=True,type=Path)
    p.add_argument("--case",required=True); p.add_argument("--dt-ps",type=float,default=1.0); p.add_argument("--maxlag-ps",type=float,default=100.0); p.add_argument("--blocks",type=int,default=1); args=p.parse_args()
    z=read_z(args.dump); maxlag=int(round(args.maxlag_ps/args.dt_ps)); bsize=len(z)//args.blocks
    if bsize<=maxlag: raise ValueError("block duration must exceed lag")
    t=np.arange(maxlag+1)*args.dt_ps; grid=np.logspace(0,2,201)
    mblocks=[]; ablocks=[]
    for b in range(args.blocks):
        mb=msd(z[b*bsize:(b+1)*bsize],maxlag)
        mblocks.append(mb); ablocks.append(local_decade_slope(t,mb,grid))
    ma=np.asarray(mblocks); aa=np.asarray(ablocks); valid=np.isfinite(aa).all(axis=0)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["time_ps","msd_z_A2_mean","msd_z_A2_block_sem","alpha_z_loglog_1decade_mean","alpha_z_loglog_1decade_block_sem"])
        for tt,am,ase in zip(grid[valid],aa[:,valid].mean(axis=0),aa[:,valid].std(axis=0,ddof=1)/np.sqrt(args.blocks)):
            i=int(np.argmin(np.abs(t-tt)))
            w.writerow([tt,ma[:,i].mean(),ma[:,i].std(ddof=1)/np.sqrt(args.blocks) if args.blocks>1 else float("nan"),am,ase])
    meta={"case":args.case,"n_frames":int(len(z)),"n_oxygen":int(z.shape[1]),"dt_ps":args.dt_ps,"blocks":args.blocks,"block_duration_ps":bsize*args.dt_ps,"position_frame":"lab/CNT-relative axial coordinate (CNT fixed); no instantaneous water-COM subtraction","msd":"all-origin axial MSD over full trajectory" if args.blocks==1 else "all-origin axial MSD per block","alpha_estimator":"OLS slope d log10(MSD_z)/d log10(t) on a uniform log-time grid, fixed 1-decade symmetric window","input_grid_ps":[1.0,100.0],"full_window_center_support_ps":[float(grid[valid][0]),float(grid[valid][-1])]}
    args.out.with_suffix(".json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
if __name__=="__main__": main()
