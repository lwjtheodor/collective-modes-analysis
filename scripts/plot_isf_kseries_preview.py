"""Create a two-scale preview for 10L collective ISF k1--k10."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/collective_mode_response/isf_88_L10_100fs_10ns_k1-k10/2026-08-20/analysis_5ns_lag"

def main():
    p = DATA / "per_replica/88_L10_rep1_collective_k1-k10.npz"
    with np.load(p) as z:
        t, F, ns, ks = z["time_ps"], z["F_total_normalized"], z["n"], z["k_inv_A"]
    plt.rcParams.update({"font.family":"Arial","font.size":7,"axes.linewidth":1.0,
                         "xtick.direction":"out","ytick.direction":"out","pdf.fonttype":42})
    fig = plt.figure(figsize=(5.5, 2.65), facecolor="white")
    # BBox-first: equal outer panels with an explicitly reserved top note.
    boxes = [(0.10,0.19,0.385,0.68),(0.57,0.19,0.385,0.68)]
    cmap = plt.get_cmap("viridis")
    for j, (x,y,w,h) in enumerate(boxes):
        ax = fig.add_axes([x,y,w,h])
        for i, (n,k) in enumerate(zip(ns,ks)):
            ax.plot(t, F[:,i], color=cmap(i/9), lw=1.05,
                    label=fr"$k_{{{n}}}={k:.4f}\ \mathrm{{\AA}}^{{-1}}$")
        ax.axhline(0,color="0.45",lw=1.0,zorder=0)
        ax.set_xlim(0, 1000 if j == 0 else 5000)
        ax.set_ylim(-0.22, 1.05)
        ax.set_xlabel(r"$t\ (\mathrm{ps})$")
        ax.set_ylabel(r"$F(k,t)/F(k,0)$")
        ax.tick_params(top=False,right=False,length=3,pad=2)
        ax.spines[["top","right"]].set_visible(False)
        fig.text(x-0.023,y+h+0.01, f"({'ab'[j]})", fontsize=9,fontweight="bold")
    boxes_ax = fig.axes[0]
    boxes_ax.legend(ncol=2, fontsize=5.8, frameon=False, handlelength=1.5,
                    columnspacing=0.65, labelspacing=0.24, loc="upper right")
    fig.text(0.5,0.975,r"$(8,8)$, 10L, 100 fs / 10 ns; collective oxygen ISF; preliminary single replica",
             ha="center",va="top",fontsize=7)
    out = DATA / "figures"; out.mkdir(exist_ok=True)
    fig.savefig(out / "ISF_88_L10_k1-k10_rep1_preview.png",dpi=600,facecolor="white")
    fig.savefig(out / "ISF_88_L10_k1-k10_rep1_preview.pdf",facecolor="white")
    plt.close(fig)
if __name__ == "__main__": main()
