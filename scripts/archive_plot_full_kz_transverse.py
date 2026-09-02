#!/usr/bin/env python3
"""Archive and visualize complete same-kz CNT transverse-current analysis."""
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

CASES=["7_7_L5_rep1","8_8_L5_rep1","9_9_L5_rep1","17_0_L5_rep1"]
COLORS={"7_7_L5_rep1":"#3568a8","8_8_L5_rep1":"#d27743","9_9_L5_rep1":"#4b9366","17_0_L5_rep1":"#8655a1"}
LABELS={k:"("+k.split("_")[0]+","+k.split("_")[1]+")" for k in CASES}


def read_case(root, case):
    folder=root/"output"/case
    meta=json.loads((folder/"summary.json").read_text(encoding="utf-8"))
    table=np.genfromtxt(folder/"kz_transverse_spectra.csv",delimiter=",",names=True,dtype=None,encoding="utf-8")
    return meta,table


def write_archives(out, cases):
    fields=["case_id","n_frames","n_water","dt_ps","cylindrical_axis","lz_A","n","kz_inv_A","target_k_inv_A","welch_nperseg","welch_segments","axial_peak_freq_ps_inv","radial_to_axial_power_at_peak","theta_to_axial_power_at_peak","coh_zr_at_peak","coh_ztheta_at_peak"]
    with (out/"kz_transverse_peak_metrics.csv").open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader()
        for meta,_ in cases.values(): w.writerow({name:meta[name] for name in fields})
    spectral_fields=["case_id","freq_ps_inv","P_z","P_r","P_theta","coh_zr","coh_ztheta","kz_inv_A","axial_peak_freq_ps_inv","welch_segments","dt_ps"]
    with (out/"kz_transverse_spectra_all.csv").open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=spectral_fields); w.writeheader()
        for case,(meta,t) in cases.items():
            for row in t:
                w.writerow({"case_id":case,"freq_ps_inv":float(row["freq_ps_inv"]),"P_z":float(row["P_z"]),"P_r":float(row["P_r"]),"P_theta":float(row["P_theta"]),"coh_zr":float(row["coh_zr"]),"coh_ztheta":float(row["coh_ztheta"]),"kz_inv_A":meta["kz_inv_A"],"axial_peak_freq_ps_inv":meta["axial_peak_freq_ps_inv"],"welch_segments":meta["welch_segments"],"dt_ps":meta["dt_ps"]})
    provenance={"observable":"complex cylindrical water-current modes at matched axial wave number","components":{"P_z":"axial","P_r":"radial","P_theta":"circumferential"},"axis":"fixed box center (fixed CNT axis)","input":"complete 1 ns water trajectories, 10 fs cadence, one trajectory per chirality","estimator":"complex-current Welch spectra; 8192-frame (81.92 ps) Hann segments with 50% overlap","statistical_boundary":"23 overlapping Welch segments per trajectory; one trajectory per chirality, so no across-replicate uncertainty or causal inference","raw_source_subdirectories":CASES}
    (out/"README.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")


def plot(out, cases):
    mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans","sans-serif"],"font.size":7,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.8,"svg.fonttype":"none","pdf.fonttype":42})
    fig,axes=plt.subplots(2,2,figsize=(7.15,5.05),layout="constrained")
    for case,(meta,t) in cases.items():
        f=t["freq_ps_inv"]; keep=(f>0)&(f<=0.8); color=COLORS[case]; label=LABELS[case]
        pz=t["P_z"]; peak=np.max(pz[keep])
        axes[0,0].plot(f[keep],pz[keep]/peak,color=color,lw=1.25,label=label)
        axes[0,1].plot(f[keep],t["P_theta"][keep]/peak,color=color,lw=1.25)
        axes[1,0].plot(f[keep],t["P_r"][keep]/peak,color=color,lw=1.25)
        axes[1,1].scatter(meta["theta_to_axial_power_at_peak"],meta["radial_to_axial_power_at_peak"],color=color,s=35,zorder=3)
        axes[1,1].annotate(label,(meta["theta_to_axial_power_at_peak"],meta["radial_to_axial_power_at_peak"]),xytext=(3,3),textcoords="offset points",color=color,fontsize=7)
    axes[0,0].set(title="A  Axial spectrum",ylabel=r"$P_z(f)/\max P_z$")
    axes[0,1].set(title="B  Circumferential spectrum",ylabel=r"$P_\theta(f)/\max P_z$")
    axes[1,0].set(title="C  Radial spectrum",ylabel=r"$P_r(f)/\max P_z$")
    axes[1,1].set(title="D  At axial spectral peak",xlabel=r"$P_\theta/P_z$",ylabel=r"$P_r/P_z$",xscale="log",yscale="log")
    for ax in axes.flat: ax.set_xlabel(r"Frequency (ps$^{-1}$)") if ax is not axes[1,1] else None
    for ax in (axes[0,0],axes[0,1],axes[1,0]): ax.set_xlim(0,0.8)
    axes[0,1].set_yscale("log"); axes[0,1].set_ylim(1e-6,3)
    axes[1,0].set_yscale("log"); axes[1,0].set_ylim(1e-7,1)
    axes[1,1].xaxis.set_major_locator(mticker.FixedLocator([0.1,0.2,0.3])); axes[1,1].xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f")); axes[1,1].xaxis.set_minor_locator(mticker.NullLocator())
    axes[1,1].yaxis.set_major_locator(mticker.FixedLocator([0.001,0.002,0.003])); axes[1,1].yaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f")); axes[1,1].yaxis.set_minor_locator(mticker.NullLocator())
    axes[0,0].legend(loc="upper right",fontsize=7)
    base=out/"full_kz_axial_transverse_comparison"
    fig.savefig(base.with_suffix(".png"),dpi=300,bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"),bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"),bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"),dpi=600,bbox_inches="tight")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",type=Path,required=True); args=ap.parse_args()
    out=args.root/"archive"; out.mkdir(parents=True,exist_ok=True)
    cases={case:read_case(args.root,case) for case in CASES}
    write_archives(out,cases); plot(out,cases)


if __name__=="__main__": main()
