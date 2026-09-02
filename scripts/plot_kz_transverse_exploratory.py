#!/usr/bin/env python3
"""Render the 50 ps same-kz axial/transverse screening result."""
import argparse, csv, json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


COLORS={"7_7":"#3b6fb6","8_8":"#d46a3b","9_9":"#4f9d69","17_0":"#8a5aa5"}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",required=True,type=Path); args=ap.parse_args()
    cases=["7_7","8_8","9_9","17_0"]; rows=[]; data={}
    for case in cases:
        out=args.root/"outputs"/case
        meta=json.loads((out/"summary.json").read_text(encoding="utf-8"))
        table=np.genfromtxt(out/"kz_transverse_spectra.csv",delimiter=",",names=True,dtype=None,encoding="utf-8")
        data[case]=(meta,table); rows.append(meta)
    args.root.mkdir(parents=True,exist_ok=True)
    with (args.root/"same_kz_transverse_screen_summary.csv").open("w",newline="") as fh:
        fields=["case_id","n_frames","n_water","dt_ps","lz_A","n","kz_inv_A","axial_peak_freq_ps_inv","radial_to_axial_power_at_peak","theta_to_axial_power_at_peak","coh_zr_at_peak","coh_ztheta_at_peak","welch_segments"]
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows([{k:r[k] for k in fields} for r in rows])
    mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":8,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.8,"svg.fonttype":"none","pdf.fonttype":42})
    fig,ax=plt.subplots(2,2,figsize=(7.15,5.1),layout="constrained")
    for case,(meta,t) in data.items():
        f=t["freq_ps_inv"]; valid=(f>0)&(f<=1.0); pz=t["P_z"]; pt=t["P_theta"]; pr=t["P_r"]
        color=COLORS[case]; label=rf"({case.replace('_',',')})"
        ax[0,0].plot(f[valid],pz[valid]/np.max(pz[valid]),lw=1.35,color=color,label=label)
        ax[0,1].plot(f[valid],pt[valid]/np.max(pz[valid]),lw=1.35,color=color,label=label)
        ax[1,0].plot(f[valid],pr[valid]/np.max(pz[valid]),lw=1.35,color=color,label=label)
        ax[1,1].scatter(meta["theta_to_axial_power_at_peak"],meta["radial_to_axial_power_at_peak"],s=42,color=color,zorder=3)
        ax[1,1].annotate(label,(meta["theta_to_axial_power_at_peak"],meta["radial_to_axial_power_at_peak"]),xytext=(4,3),textcoords="offset points",color=color,fontsize=8)
    ax[0,0].set(title="A  Axial mode spectrum",ylabel=r"$P_z(f)/\max P_z$")
    ax[0,1].set(title="B  Circumferential spectrum",ylabel=r"$P_\theta(f)/\max P_z$")
    ax[1,0].set(title="C  Radial spectrum",xlabel=r"Frequency (ps$^{-1}$)",ylabel=r"$P_r(f)/\max P_z$")
    ax[1,1].set(title="D  At the axial spectral peak",xlabel=r"$P_\theta/P_z$",ylabel=r"$P_r/P_z$",xscale="log",yscale="log")
    for panel in (ax[0,0],ax[0,1],ax[1,0]): panel.set(xlim=(0,1.0),xlabel=r"Frequency (ps$^{-1}$)"); panel.axhline(0,color="0.75",lw=.6)
    ax[1,0].set_yscale("log"); ax[1,0].set_ylim(1e-6,2)
    ax[0,1].set_yscale("log"); ax[0,1].set_ylim(1e-5,3)
    ax[0,0].legend(loc="upper right",fontsize=7,frameon=False)
    fig.text(.5,.01,r"Same physical $k_z\simeq0.0623$ $\mathrm{\AA^{-1}}$; 50 ps, one trajectory/chirality, 3 Welch segments",ha="center",fontsize=7,color="0.28")
    base=args.root/"same_kz_axial_transverse_screen_50ps"
    for ext,kw in (("png",{"dpi":300}),("pdf",{}),("svg",{}),("tiff",{"dpi":600})): fig.savefig(base.with_suffix("."+ext),bbox_inches="tight",**kw)


if __name__=="__main__": main()
