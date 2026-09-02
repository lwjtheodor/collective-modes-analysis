"""Matched-protocol VACF-integral and current-mode lock-in post-processing.

For one `(8,8)` weak-NH/no-momentum oxygen dump this writes (i) unnormalised
peculiar VACF, its Green--Kubo integral D(t), direct peculiar MSD and alpha;
and (ii) block-level lock-in quadratures against the matched physical k=2pi/L0
current mode (n=box length).  Gamma is explicitly labelled operational:
1/tau_1e of the mode-current ACF, not a transport damping constant.
"""
import argparse, csv, json, math
from pathlib import Path
import numpy as np

VEL_TO_PS = 1.0e6  # (A/fs)^2 -> (A/ps)^2

def read_dump(path):
    zframes=[]; vframes=[]; mol_ref=None
    with path.open(errors="replace") as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.startswith("ITEM: TIMESTEP"): continue
            f.readline(); f.readline(); nat=int(f.readline())
            f.readline(); bounds=[[float(x) for x in f.readline().split()[:2]] for _ in range(3)]
            lz=bounds[2][1]-bounds[2][0]; cols=f.readline().split()[2:]; ix={x:i for i,x in enumerate(cols)}
            need={"mol","type","z","vz"}
            if not need.issubset(ix): raise ValueError("dump must contain mol,type,z,vz")
            d={}
            for _ in range(nat):
                row=f.readline().split()
                if int(float(row[ix["type"]]))==3:
                    m=int(float(row[ix["mol"]])); z=float(row[ix["z"]])
                    if "iz" in ix: z += int(float(row[ix["iz"]]))*lz
                    d[m]=(z,float(row[ix["vz"]]))
            mol=np.asarray(sorted(d),dtype=np.int64)
            if mol_ref is None: mol_ref=mol
            elif not np.array_equal(mol,mol_ref): raise ValueError("oxygen IDs changed")
            zframes.append([d[int(m)][0] for m in mol_ref]); vframes.append([d[int(m)][1] for m in mol_ref])
    return np.asarray(zframes,dtype=np.float32),np.asarray(vframes,dtype=np.float32),float(lz)

def xcorr_sum(a):
    n=a.shape[0]; size=1<<(2*n-1).bit_length()
    fa=np.fft.fft(a,n=size,axis=0)
    return np.fft.ifft(fa*np.conjugate(fa),axis=0)[:n].sum(axis=1)

def acf_columns(x,maxlag):
    n,m=x.shape; out=np.real(xcorr_sum(x)[:maxlag+1])
    return out/(np.arange(n,n-maxlag-1,-1)*m)

def msd_columns(z,maxlag):
    n,m=z.shape; corr=np.real(xcorr_sum(z)[:maxlag+1]); sq=np.sum(z*z,axis=1)
    pref=np.r_[0.,np.cumsum(sq)]; lag=np.arange(maxlag+1); count=n-lag
    return ((pref[n]-pref[lag])+pref[n-lag]-2*corr)/(count*m)

def cumtrap(y,dt):
    out=np.zeros_like(y,dtype=float); out[1:]=np.cumsum((y[1:]+y[:-1])*0.5*dt); return out

def current_params(z,v,lz,mode_n,dt,maxlag):
    # z is unwrapped; phase must use its periodic representative.
    vv=v-v.mean(axis=1,keepdims=True); k=2*np.pi*mode_n/lz
    j=np.sum(vv*np.exp(1j*k*np.mod(z,lz)),axis=1)[:,None]
    raw=xcorr_sum(j)[:maxlag+1]/np.arange(len(j),len(j)-maxlag-1,-1)
    c=np.real(raw/raw[0]); t=np.arange(maxlag+1)*dt
    spec=np.abs(np.fft.rfft(c-c.mean())); fr=np.fft.rfftfreq(len(c),dt); pick=int(np.argmax(spec[1:])+1)
    omega=2*np.pi*fr[pick]
    below=np.where(np.abs(c)<=np.exp(-1))[0]; tau=float(t[below[0]]) if len(below) else float(t[-1])
    return {"mode_n":mode_n,"k_inv_A":float(k),"omega_rad_ps":float(omega),
            "gamma_operational_inv_ps":float(1.0/tau),"tau_1e_abs_ps":tau}

def main():
    p=argparse.ArgumentParser(); p.add_argument("--dump",type=Path,required=True); p.add_argument("--out",type=Path,required=True)
    p.add_argument("--case",required=True); p.add_argument("--length",type=int,required=True); p.add_argument("--dt-ps",type=float,default=.01)
    p.add_argument("--maxlag-ps",type=float,default=100.); p.add_argument("--blocks",type=int,default=9); args=p.parse_args()
    z,v,lz=read_dump(args.dump); maxlag=int(round(args.maxlag_ps/args.dt_ps)); bsize=len(z)//args.blocks
    if bsize<=maxlag: raise ValueError("block duration must exceed max lag")
    params=current_params(z,v,lz,args.length,args.dt_ps,min(maxlag,int(round(25/args.dt_ps))))
    t=np.arange(maxlag+1)*args.dt_ps; g=np.exp(-params["gamma_operational_inv_ps"]*t)*np.cos(params["omega_rad_ps"]*t); h=np.exp(-params["gamma_operational_inv_ps"]*t)*np.sin(params["omega_rad_ps"]*t)
    tail=(t>=5)&(t<=100); block=[]; lock=[]
    for b in range(args.blocks):
        zz=z[b*bsize:(b+1)*bsize]; vv=v[b*bsize:(b+1)*bsize]; vv=vv-vv.mean(axis=1,keepdims=True); zz=zz-zz.mean(axis=1,keepdims=True)
        c=acf_columns(vv,maxlag); d=cumtrap(c*VEL_TO_PS,args.dt_ps); m=msd_columns(zz,maxlag); mac=2*cumtrap(d,args.dt_ps)
        alpha=np.full_like(t,np.nan); alpha[1:]=2*t[1:]*d[1:]/m[1:]
        cn=c/c[0]; delta=cn-np.mean(cn[(t>=75)&(t<=100)])
        a=float(np.dot(delta[tail],g[tail])/np.dot(g[tail],g[tail])); bb=float(np.dot(delta[tail],h[tail])/np.dot(h[tail],h[tail]))
        lock.append({"block":b,"a_cos":a,"b_sin":bb,"amplitude":math.hypot(a,bb),"phase_rad":math.atan2(bb,a)})
        block.append(np.vstack([c,d,m,mac,alpha,cn]))
    x=np.asarray(block); mean=x.mean(axis=0); sem=x.std(axis=0,ddof=1)/math.sqrt(args.blocks)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    fields=["lag_ps","cvv_raw_A2_fs2_mean","cvv_raw_A2_fs2_sem","D_A2_ps_mean","D_A2_ps_sem","MSD_direct_A2_mean","MSD_direct_A2_sem","MSD_from_C_A2_mean","MSD_from_C_A2_sem","alpha_from_D_MSD_mean","alpha_from_D_MSD_sem","vacf_normalized_mean","vacf_normalized_sem"]
    with args.out.with_suffix(".csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for i,tt in enumerate(t): w.writerow(dict(zip(fields,[tt,mean[0,i],sem[0,i],mean[1,i],sem[1,i],mean[2,i],sem[2,i],mean[3,i],sem[3,i],mean[4,i],sem[4,i],mean[5,i],sem[5,i]])))
    with args.out.with_name(args.out.stem+"_lockin_blocks.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(lock[0])); w.writeheader(); w.writerows(lock)
    meta={"case":args.case,"n_frames":int(len(z)),"n_water":int(z.shape[1]),"dt_ps":args.dt_ps,"maxlag_ps":float(t[-1]),"n_blocks":args.blocks,"block_duration_ps":bsize*args.dt_ps,"velocity_frame":"instantaneous water-COM removed","mode_parameters":params,"gamma_note":"operational 1/tau_1e_abs from matched-k current ACF; not a fitted transport damping constant","lockin_window_ps":[5,100],"lockin_baseline":"block mean of normalized VACF over 75-100 ps"}
    args.out.with_suffix(".json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
if __name__=="__main__": main()
