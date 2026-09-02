"""Nature-style summary of VACF integral and matched-current lock-in analysis."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

DATA = Path(r"C:\Users\s1365\.codex\visualizations\2026\08\06\019fd4c1-37fd-7b92-b606-edcea5bf0c15\stage_vacf_integral_lockin_20260806\output")
ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
OUT = ROOT / "assets" / "vacf_integral_msd_alpha_lockin_nature"
COLORS = {1:"#203864", 2:"#2F75B5", 3:"#4FA3A5", 4:"#5B9A6D", 5:"#A07A39"}

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":7, "axes.labelsize":7.2, "axes.titlesize":7.8, "xtick.labelsize":6.2,
    "ytick.labelsize":6.2, "legend.fontsize":5.8, "axes.linewidth":0.75,
    "axes.spines.top":False, "axes.spines.right":False, "svg.fonttype":"none", "pdf.fonttype":42,
})

def group_curves(L):
    files = sorted(p for p in DATA.glob(f"*L{L}_rep*.csv") if not p.name.endswith("_lockin_blocks.csv"))
    fs = [pd.read_csv(p) for p in files]
    t = fs[0].lag_ps.to_numpy()
    keys = ["D_A2_ps_mean", "MSD_direct_A2_mean", "MSD_from_C_A2_mean", "alpha_from_D_MSD_mean"]
    result = {"t":t}
    for key in keys:
        y = np.vstack([f[key].to_numpy() for f in fs])
        if key == "MSD_direct_A2_mean":
            y = y-y[:, :1]
        n = np.isfinite(y).sum(axis=0)
        mean = np.nansum(y, axis=0) / np.where(n > 0, n, np.nan)
        sem = np.sqrt(np.nansum((y-mean)**2, axis=0) / np.where(n > 1, n-1, np.nan)) / np.sqrt(n)
        result[key] = (mean, sem)
    return result

def lockin(L):
    means = []
    for p in sorted(DATA.glob(f"*L{L}_rep*_lockin_blocks.csv")):
        means.append(pd.read_csv(p).amplitude.to_numpy().mean())
    means = np.asarray(means)
    return means.mean(), means.std(ddof=1)/np.sqrt(3)

def band(ax, t, m, s, color, label=None):
    ax.fill_between(t, m-s, m+s, color=color, alpha=0.16, lw=0)
    ax.plot(t, m, color=color, lw=1.05, label=label)

def main():
    curves = {L:group_curves(L) for L in range(1,6)}
    fig, ax = plt.subplots(2, 2, figsize=(7.15, 4.5))
    a, b, c, d = ax.ravel()
    for L, x in curves.items():
        mask = (x["t"] >= 0) & (x["t"] <= 100)
        band(a, x["t"][mask], x["D_A2_ps_mean"][0][mask], x["D_A2_ps_mean"][1][mask], COLORS[L], f"{L}L")
    a.set(xlim=(0,100), ylabel=r"$D(t)$ ($\AA^2$ ps$^{-1}$)", xlabel=r"lag time, $t$ (ps)")
    a.set_title("VACF integral", loc="left", fontweight="bold")
    a.legend(title="box length", ncol=2, loc="upper right", handlelength=1.4, labelspacing=.25, borderpad=.2)
    for L, x in curves.items():
        mask = (x["t"] >= 0) & (x["t"] <= 100)
        band(b, x["t"][mask], x["MSD_direct_A2_mean"][0][mask], x["MSD_direct_A2_mean"][1][mask], COLORS[L], f"direct {L}L")
        b.plot(x["t"][mask], x["MSD_from_C_A2_mean"][0][mask], color=COLORS[L], lw=.9, ls="--")
    b.set(xlim=(0,100), ylabel=r"$M(t)$ ($\AA^2$)", xlabel=r"lag time, $t$ (ps)")
    b.set_title("MSD closure", loc="left", fontweight="bold")
    b.text(.02,.04,"solid: direct, origin-corrected\ndashed: reconstructed from VACF", transform=b.transAxes, fontsize=5.5)
    for L, x in curves.items():
        mask = (x["t"] >= 5) & (x["t"] <= 100)
        band(c, x["t"][mask], x["alpha_from_D_MSD_mean"][0][mask], x["alpha_from_D_MSD_mean"][1][mask], COLORS[L])
    c.axhline(1, color="#555555", lw=.65, zorder=0)
    c.set(xlim=(5,100), ylim=(.3,1.2), ylabel=r"$\alpha(t)=2tD(t)/M(t)$", xlabel=r"lag time, $t$ (ps)")
    c.set_title("instantaneous transport exponent", loc="left", fontweight="bold")
    amp = np.array([lockin(L) for L in range(1,6)])
    Ls = np.arange(1,6)
    d.errorbar(Ls, amp[:,0], yerr=amp[:,1], color="#263238", lw=.8, marker="o", ms=4, capsize=2)
    d.scatter(Ls, amp[:,0], s=22, c=[COLORS[L] for L in Ls], zorder=3)
    d.set(xlim=(.65,5.35), xticks=Ls, xlabel="box length", ylabel=r"lock-in amplitude, $A_n$ (norm. VACF)")
    d.set_title("matched-current lock-in, 5-100 ps", loc="left", fontweight="bold")
    d.text(.03,.96, r"$\omega_n,k_n$ from matched current mode"+"\n"+r"$\gamma=1/\tau_{1/e}$ (operational; not fitted damping)", transform=d.transAxes, va="top", fontsize=5.35)
    for tag, aa in zip("abcd", (a,b,c,d)):
        aa.text(.01,.98,tag, transform=aa.transAxes, va="top", ha="left", fontsize=9, fontweight="bold")
        aa.tick_params(length=2.5, width=.65)
    fig.text(.5,.005, "(8,8), weak Nose-Hoover, no momentum removal; O 10 fs; 1 ns; n=3 trajectories per L. Bands/error bars: replica SEM.", ha="center", fontsize=6.0)
    fig.subplots_adjust(left=.10, right=.985, bottom=.115, top=.94, hspace=.34, wspace=.29)
    fig.savefig(str(OUT)+".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(OUT)+".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(OUT)+".pdf", bbox_inches="tight")
    fig.savefig(str(OUT)+".svg", bbox_inches="tight")

if __name__ == "__main__":
    main()
