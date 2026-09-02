"""Audit the existing (8,8) 10L n=1 longitudinal CJJ through 100 ps."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

SRC = Path(r"F:/ccfep_gcmc_archive_20260814/viscfric_length_all_chirality_RH75_20260806/2026-08-11/output")
OUT = Path("results/collective_mode_response/cjj_k1_10L_100fs_10ns_audit/2026-08-20")

def damped_cosine(t, gamma, omega):
    return np.exp(-gamma * t) * np.cos(omega * t)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    curves=[]; labels=[]
    for p in sorted(SRC.glob("8_8_L10_rep*_CJJ_alln.csv")):
        frame=pd.read_csv(p); frame=frame[frame["n"]==1].sort_values("lag_ps")
        curves.append(frame["CJJ_normalized"].to_numpy()); labels.append(p.stem)
    time=frame["lag_ps"].to_numpy(); curves=np.asarray(curves); mean=curves.mean(0); sem=curves.std(0,ddof=1)/np.sqrt(curves.shape[0])
    cut=time<=100; records=[]
    for label, curve in list(zip(labels, curves))+[("mean_3rep",mean)]:
        pars,_=curve_fit(damped_cosine,time[cut],curve[cut],p0=(0.01,0.11),bounds=([0,0.05],[1,0.2]))
        gamma,omega=map(float,pars); c100=float(np.interp(100,time,curve))
        records.append({"series":label,"fit_window_ps":"0-100","gamma_ps_inv":gamma,"tau_ps":1/gamma,
                        "omega_rad_ps":omega,"period_ps":2*np.pi/omega,"CJJ_over_CJJ0_at_100ps":c100,
                        "absolute_remaining_fraction_at_100ps":abs(c100)})
    pd.DataFrame({"lag_ps":time,"CJJ_mean_normalized":mean,"CJJ_sem_normalized":sem}).to_csv(OUT/"CJJ_k1_mean_sem.csv",index=False)
    (OUT/"fit_audit.json").write_text(json.dumps(records,indent=2),encoding="utf-8")
    print(json.dumps(records,indent=2))
if __name__=="__main__": main()
