"""Nature-style low-frequency VACF tail morphology across CNT chiralities."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
NEW = Path(r"C:\Users\s1365\.codex\visualizations\2026\08\06\019fd4c1-37fd-7b92-b606-edcea5bf0c15\stage_lowfreq_vacf_crosschirality_20260806\output")
OUT = ROOT / "assets" / "lowfreq_crosschirality_vacf_tail_nature"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.75, "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

COLORS = {1: "#203864", 2: "#2F75B5", 3: "#4FA3A5", 4: "#5B9A6D", 5: "#A07A39"}
PANELS = [("(7,7)", "7_7", [2, 3, 4, 5]), ("(8,8)", "8_8", [1, 2, 3, 4, 5]), ("(9,9)", "9_9", [2, 3, 4])]

def load_replicates(key, length):
    if key == "8_8":
        pats = [ROOT / "lowfreq" / f"{length}L_rep1.csv", ROOT / "lowfreq" / f"{length}L_seed2.csv", ROOT / "lowfreq" / f"{length}L_seed3.csv"]
    else:
        pats = [NEW / f"{key}_{length}L_rep{i}.csv" for i in range(1, 4)]
    frames = [pd.read_csv(p) for p in pats]
    t = frames[0]["lag_ps"].to_numpy()
    y = np.vstack([f["vacf_peculiar_mean"].to_numpy() for f in frames])
    mask = (t >= 5) & (t <= 100)
    return t[mask], y[:, mask]

def main():
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.75), sharey=True)
    for idx, (title, key, lengths) in enumerate(PANELS):
        ax = axes[idx]
        ax.axhline(0, color="#4C4C4C", lw=0.65, zorder=0)
        for L in lengths:
            t, y = load_replicates(key, L)
            mean = y.mean(axis=0)
            sem = y.std(axis=0, ddof=1) / np.sqrt(y.shape[0])
            ax.fill_between(t, mean-sem, mean+sem, color=COLORS[L], alpha=0.17, lw=0)
            ax.plot(t, mean, color=COLORS[L], lw=1.1, label=f"{L}L")
        ax.set_xlim(5, 100)
        ax.set_ylim(-0.0065, 0.0115)
        ax.set_xticks([5, 25, 50, 75, 100])
        ax.set_title(title, loc="left", fontweight="bold", pad=4)
        ax.text(0.01, 0.96, chr(97+idx), transform=ax.transAxes, va="top", ha="left", fontweight="bold", fontsize=9)
        ax.text(0.14, 0.96, "O 1 ps; 20 ns; NVT baseline", transform=ax.transAxes, va="top", fontsize=5.7)
        ax.set_xlabel("lag time, t (ps)")
        if idx == 0:
            ax.set_ylabel(r"peculiar axial VACF, $C_{vv}(t)/C_{vv}(0)$")
        if key == "9_9":
            ax.text(0.98, 0.055, "5L not available", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.7, color="#555555")
        ax.tick_params(length=2.5, width=0.65)
    axes[-1].legend(title="box length", loc="upper right", ncol=1, handlelength=1.5, labelspacing=0.3, borderpad=0.2)
    fig.text(0.5, 0.01, "All panels: 330 K, RH75; z-momentum removal every 5 ps; n = 3 trajectories per L; shaded bands = replica SEM.", ha="center", fontsize=6.2)
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.23, top=0.84, wspace=0.18)
    fig.savefig(str(OUT) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(OUT) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(OUT) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT) + ".svg", bbox_inches="tight")

if __name__ == "__main__":
    main()
