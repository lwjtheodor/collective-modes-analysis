#!/usr/bin/env python3
"""Stream full CNT-water trajectories into matched-k longitudinal diagnostics."""
import argparse, csv, json
from pathlib import Path
import numpy as np


def acf_complex(q, maxlag):
    n=len(q); nfft=1 << (2*n-1).bit_length()
    corr=np.fft.ifft(np.fft.fft(q,nfft)*np.conjugate(np.fft.fft(q,nfft)))[:maxlag+1].real
    return corr/np.arange(n,n-maxlag-1,-1,dtype=float)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dump",required=True,type=Path); ap.add_argument("--dt-ps",required=True,type=float)
    ap.add_argument("--case-id",required=True); ap.add_argument("--outdir",required=True,type=Path)
    ap.add_argument("--target-k",type=float,default=0.0623084632); ap.add_argument("--max-lag-ps",type=float,default=100.0)
    ap.add_argument("--radial-bin-A",type=float,default=0.25); ap.add_argument("--m-max",type=int,default=4)
    args=ap.parse_args(); rho=[]; current=[]; steps=[]; helical=[]; hist=None; box=None; center=None; nwater=None; k=None
    with args.dump.open(errors="replace") as fh:
        while True:
            line=fh.readline()
            if not line: break
            if line != "ITEM: TIMESTEP\n": raise ValueError("unexpected dump marker: %r" % line)
            step=int(fh.readline())
            if fh.readline() != "ITEM: NUMBER OF ATOMS\n": raise ValueError("missing atom count")
            natom=int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"): raise ValueError("missing bounds")
            bounds=np.array([list(map(float,fh.readline().split()[:2])) for _ in range(3)])
            header=fh.readline().split()[2:]; col={name:i for i,name in enumerate(header)}
            required={"type","x","y","z","vx","vy","vz"}
            if not required.issubset(col): raise ValueError("missing columns: %s" % (required-set(col)))
            xyz=[]; vel=[]
            for _ in range(natom):
                values=fh.readline().split()
                if int(float(values[col["type"]])) != 3: continue
                xyz.append([float(values[col[a]]) for a in ("x","y","z")])
                vel.append([float(values[col[a]]) for a in ("vx","vy","vz")])
            if not xyz: raise ValueError("no oxygen atoms (type 3)")
            frame_box=bounds[:,1]-bounds[:,0]
            if box is None:
                box=frame_box; center=bounds.mean(axis=1); nwater=len(xyz)
                nmode=max(1,int(round(args.target_k*box[2]/(2*np.pi)))); k=2*np.pi*nmode/box[2]
                edges=np.arange(0,0.5*min(box[0],box[1])+args.radial_bin_A,args.radial_bin_A); hist=np.zeros(len(edges)-1,dtype=np.int64)
            elif len(xyz)!=nwater or not np.allclose(frame_box,box): raise ValueError("water count or box changed")
            xyz=np.asarray(xyz); vel=np.asarray(vel); vel-=vel.mean(axis=0,keepdims=True)
            dxy=xyz[:,:2]-center[:2]; radius=np.hypot(dxy[:,0],dxy[:,1]); theta=np.arctan2(dxy[:,1],dxy[:,0])
            phase=np.exp(1j*k*xyz[:,2]); rho.append(np.sum(phase)); current.append(np.sum(vel[:,2]*phase))
            helical.append([np.sum(np.exp(1j*(k*xyz[:,2]+m*theta))) for m in range(args.m_max+1)])
            hist+=np.histogram(radius,bins=edges)[0]; steps.append(step)
    rho=np.asarray(rho); current=np.asarray(current); helical=np.asarray(helical); steps=np.asarray(steps)
    observed=float(np.median(np.diff(steps))*0.0005)
    if not np.isclose(observed,args.dt_ps,rtol=1e-5,atol=1e-8): raise ValueError("cadence %s != declared %s" % (observed,args.dt_ps))
    maxlag=min(int(round(args.max_lag_ps/args.dt_ps)),len(rho)-1); time=np.arange(maxlag+1)*args.dt_ps
    fk=acf_complex(rho,maxlag); cj=acf_complex(current,maxlag)
    args.outdir.mkdir(parents=True,exist_ok=True)
    with (args.outdir/"longitudinal_acf.csv").open("w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["case_id","lag_ps","F_k_norm","C_J_norm"])
        for t,a,b in zip(time,fk/fk[0],cj/cj[0]): w.writerow([args.case_id,float(t),float(a),float(b)])
    shell_vol=np.pi*(edges[1:]**2-edges[:-1]**2)*box[2]; density=hist/(len(rho)*shell_vol)
    with (args.outdir/"radial_density.csv").open("w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["case_id","r_inner_A","r_outer_A","r_center_A","oxygen_number_density_A3","frame_averaged_oxygen_count"])
        for lo,hi,d,c in zip(edges[:-1],edges[1:],density,hist/float(len(rho))): w.writerow([args.case_id,float(lo),float(hi),float(.5*(lo+hi)),float(d),float(c)])
    with (args.outdir/"helical_density_modes.csv").open("w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["case_id","m","kz_inv_A","S_mk_per_water"])
        for m in range(args.m_max+1): w.writerow([args.case_id,m,float(k),float(np.mean(np.abs(helical[:,m])**2)/nwater)])
    summary={"case_id":args.case_id,"n_frames":int(len(rho)),"n_water":int(nwater),"dt_ps":args.dt_ps,"duration_ps":float((len(rho)-1)*args.dt_ps),"fixed_axis_xy_A":[float(center[0]),float(center[1])],"lz_A":float(box[2]),"n":int(nmode),"kz_inv_A":float(k),"target_k_inv_A":args.target_k,"S_k_per_water":float(np.mean(np.abs(rho)**2)/nwater),"max_lag_ps":float(time[-1]),"radial_bin_A":args.radial_bin_A,"m_max":args.m_max}
    (args.outdir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2)); print("Full static-longitudinal diagnosis finished successfully.")


if __name__=="__main__": main()
