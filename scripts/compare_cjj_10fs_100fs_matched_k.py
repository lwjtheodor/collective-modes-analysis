"""Protocol-separated CJJ comparison at exact matched physical k."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT=Path(r"F:/ccfep_gcmc_archive_20260814/viscfric_length_all_chirality_RH75_20260806/2026-08-11/output")
OUT=Path("results/collective_mode_response/cjj_10fs_vs_100fs_matched_k/2026-08-20")
# 10-fs box origin, mode n, versus 10L (100-fs) matching integer mode.
MATCHES=[("5L",1,2),("5L",2,4),("2L+4L",None,5),("5L",3,6),("5L",4,8),("2L+3L+4L+5L",None,10)]
def model(t,gamma,omega): return np.exp(-gamma*t)*np.cos(omega*t)
def curves(label,n):
    ls=[]
    for p in sorted(ROOT.glob(f"8_8_L{label}_rep*_CJJ_alln.csv")):
        d=pd.read_csv(p); d=d[d.n==n].sort_values("lag_ps"); ls.append(d.CJJ_normalized.to_numpy())
    return d.lag_ps.to_numpy(),np.asarray(ls)
def metrics(time,rows):
    mean=rows.mean(0); cut=time<=100; p,_=curve_fit(model,time[cut],mean[cut],p0=(.014,.1),bounds=([0,.04],[1,.2]))
    return float(p[0]),float(p[1]),float(np.interp(100,time,mean)),float(rows.std(0,ddof=1)[np.searchsorted(time,100)]/np.sqrt(rows.shape[0]))
def main():
    OUT.mkdir(parents=True,exist_ok=True); rec=[]
    for source,n10,n100 in MATCHES:
        if n10 is not None:
            t,a=curves(str(int(source[0])),n10); g,w,c,se=metrics(t,a); k=2*np.pi*n10/(100.84*int(source[0]))
            rec.append(dict(protocol="10 fs / 1 ns",source=source,n=n10,k_inv_A=k,gamma_ps_inv=g,omega_rad_ps=w,period_ps=2*np.pi/w,C100_norm=c,C100_sem=se,n_replica=a.shape[0]))
        else:
            parts=[]
            for L in source.split("+"):
                nn=round(n100*int(L[0])/10); t,x=curves(str(int(L[0])),nn); parts.extend(x)
            g,w,c,se=metrics(t,np.asarray(parts)); k=2*np.pi*n100/1008.39998
            rec.append(dict(protocol="10 fs / 1 ns",source=source,n="matched",k_inv_A=k,gamma_ps_inv=g,omega_rad_ps=w,period_ps=2*np.pi/w,C100_norm=c,C100_sem=se,n_replica=len(parts)))
        t,b=curves("10",n100);g,w,c,se=metrics(t,b);k=2*np.pi*n100/1008.39998
        rec.append(dict(protocol="100 fs / 10 ns",source="10L",n=n100,k_inv_A=k,gamma_ps_inv=g,omega_rad_ps=w,period_ps=2*np.pi/w,C100_norm=c,C100_sem=se,n_replica=b.shape[0]))
    pd.DataFrame(rec).to_csv(OUT/"matched_k_protocol_comparison.csv",index=False)
    print(pd.DataFrame(rec).to_string(index=False))
if __name__=='__main__':main()
