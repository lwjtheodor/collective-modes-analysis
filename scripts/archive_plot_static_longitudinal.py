#!/usr/bin/env python3
"""Archive and plot complete static/longitudinal CNT diagnosis outputs."""
import argparse, csv, json, math
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

CASES=["7_7_L5_rep1","8_8_L5_rep1","9_9_L5_rep1","17_0_L5_rep1"]
COLORS={"7_7_L5_rep1":"#3568a8","8_8_L5_rep1":"#d27743","9_9_L5_rep1":"#4b9366","17_0_L5_rep1":"#8655a1"}
LABELS={key:"("+key.split("_")[0]+","+key.split("_")[1]+")" for key in CASES}


def table(path):
    return np.genfromtxt(path,delimiter=",",names=True,dtype=None,encoding="utf-8")


def crossing(time, y, threshold):
    hits=np.flatnonzero(y<=threshold)
    if not len(hits): return float("nan")
    i=int(hits[0])
    if i==0: return float(time[0])
    return float(time[i-1]+(threshold-y[i-1])*(time[i]-time[i-1])/(y[i]-y[i-1]))


def write_combined(out, cases):
    combined={"longitudinal_acf_all.csv":["case_id","lag_ps","F_k_norm","C_J_norm"],"radial_density_all.csv":["case_id","r_inner_A","r_outer_A","r_center_A","oxygen_number_density_A3","frame_averaged_oxygen_count"],"helical_density_modes_all.csv":["case_id","m","kz_inv_A","S_mk_per_water"]}
    for name,fields in combined.items():
        with (out/name).open("w",newline="",encoding="utf-8") as fh:
            w=csv.DictWriter(fh,fieldnames=fields); w.writeheader()
            for case,record in cases.items():
                source=record["acf"] if name.startswith("longitudinal") else record["radial"] if name.startswith("radial") else record["helical"]
                for row in source: w.writerow({field:row[field].item() if hasattr(row[field],"item") else row[field] for field in fields})
    fields=["case_id","n_frames","n_water","duration_ps","lz_A","n","kz_inv_A","S_k_per_water","F_1e_ps","C_J_min_0_25ps","C_J_min_lag_ps","r_mean_A","r_peak_A","helical_max_m","helical_max_over_m0"]
    with (out/"case_summary.csv").open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader()
        for case,record in cases.items(): w.writerow(record["derived"])
    provenance={"input":"four completed full 1 ns trajectories; 10 fs cadence; one trajectory per chirality","matched_kz":"target 0.0623084632 inverse Angstrom; n=5 in all cases","outputs":list(combined)+["case_summary.csv"],"statistical_boundary":"one trajectory per chirality; values describe mechanism screening, not cross-replicate uncertainty","definitions":{"F_1e_ps":"first linear-interpolated F(k,t)=exp(-1) crossing","C_J_min_0_25ps":"minimum normalized axial-current ACF in 0-25 ps","r_mean_A":"frame-averaged oxygen radial first moment","helical_max_over_m0":"largest m=1..4 static density-mode power divided by m=0 at matched kz"}}
    (out/"README.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")


def plot(out, cases):
    mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans","sans-serif"],"font.size":7,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.8,"svg.fonttype":"none","pdf.fonttype":42,"legend.frameon":False})
    fig,ax=plt.subplots(2,2,figsize=(7.15,5.0),layout="constrained")
    for case,rec in cases.items():
        color=COLORS[case]; label=LABELS[case]; acf=rec["acf"]; radial=rec["radial"]; hel=rec["helical"]; keep=acf["lag_ps"]<=25
        ax[0,0].plot(acf["lag_ps"][keep],acf["F_k_norm"][keep],color=color,lw=1.25,label=label)
        ax[0,1].plot(acf["lag_ps"][keep],acf["C_J_norm"][keep],color=color,lw=1.25)
        ax[1,0].plot(radial["r_center_A"],radial["oxygen_number_density_A3"],color=color,lw=1.25)
        ax[1,1].plot(hel["m"],hel["S_mk_per_water"]/hel["S_mk_per_water"][0],"o-",color=color,lw=1.25,ms=3)
    ax[0,0].set(title="A  Density correlation",xlabel="Lag (ps)",ylabel=r"$F(k,t)/F(k,0)$",xlim=(0,25),ylim=(-.25,1.05))
    ax[0,1].set(title="B  Axial-current memory",xlabel="Lag (ps)",ylabel=r"$C_J(k,t)/C_J(k,0)$",xlim=(0,25),ylim=(-.8,1.05)); ax[0,1].axhline(0,color="0.65",lw=.7)
    ax[1,0].set(title="C  Radial oxygen density",xlabel=r"Radius $r$ (Å)",ylabel=r"Density (Å$^{-3}$)",xlim=(0,13))
    ax[1,1].set(title="D  Helical density structure",xlabel=r"Azimuthal index $m$",ylabel=r"$S_{m,k}/S_{0,k}$",xlim=(0,4),yscale="log",ylim=(.1,100)); ax[1,1].set_xticks(range(5))
    ax[0,0].legend(loc="upper right",fontsize=7)
    base=out/"static_longitudinal_crosschirality"
    fig.savefig(base.with_suffix(".png"),dpi=300,bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"),bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"),bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"),dpi=600,bbox_inches="tight")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",required=True,type=Path); args=ap.parse_args(); out=args.root/"archive"; out.mkdir(parents=True,exist_ok=True)
    cases={}
    for case in CASES:
        folder=args.root/"output"/case; meta=json.loads((folder/"summary.json").read_text(encoding="utf-8")); acf=table(folder/"longitudinal_acf.csv"); radial=table(folder/"radial_density.csv"); helical=table(folder/"helical_density_modes.csv")
        window=(acf["lag_ps"]>0)&(acf["lag_ps"]<=25); imin=np.flatnonzero(window)[np.argmin(acf["C_J_norm"][window])]
        counts=radial["frame_averaged_oxygen_count"]; rmean=float(np.sum(radial["r_center_A"]*counts)/np.sum(counts)); rpeak=float(radial["r_center_A"][np.argmax(radial["oxygen_number_density_A3"])])
        hratio=helical["S_mk_per_water"][1:]/helical["S_mk_per_water"][0]; imax=int(np.argmax(hratio)+1)
        derived={"case_id":case,"n_frames":meta["n_frames"],"n_water":meta["n_water"],"duration_ps":meta["duration_ps"],"lz_A":meta["lz_A"],"n":meta["n"],"kz_inv_A":meta["kz_inv_A"],"S_k_per_water":meta["S_k_per_water"],"F_1e_ps":crossing(acf["lag_ps"],acf["F_k_norm"],math.exp(-1)),"C_J_min_0_25ps":float(acf["C_J_norm"][imin]),"C_J_min_lag_ps":float(acf["lag_ps"][imin]),"r_mean_A":rmean,"r_peak_A":rpeak,"helical_max_m":imax,"helical_max_over_m0":float(hratio[imax-1])}
        cases[case]={"meta":meta,"acf":acf,"radial":radial,"helical":helical,"derived":derived}
    write_combined(out,cases); plot(out,cases)


if __name__=="__main__": main()
