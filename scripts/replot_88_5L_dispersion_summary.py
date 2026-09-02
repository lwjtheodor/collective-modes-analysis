#!/usr/bin/env python3
"""Replot the audited (8,8) 5L LA/TA dispersion summary from its CSV."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "xtick.direction": "out", "ytick.direction": "out",
})


def main() -> None:
    root = Path("results/collective_mode_response/88_5L_LA_TAr_TAtheta_dispersion/2026-08-19")
    data = pd.read_csv(root / "derived_data" / "LA_TA_dispersion.csv")
    colors = {"LA": "#0F4D92", "TA_r": "#B64342", "TA_theta": "#42949E"}
    labels = {"LA": "LA", "TA_r": r"TA$_r$", "TA_theta": r"TA$_\theta$"}
    fig = plt.figure(figsize=(5.5, 2.55))
    ax_omega = fig.add_axes([0.12, 0.19, 0.35, 0.65])
    ax_velocity = fig.add_axes([0.60, 0.19, 0.35, 0.65])
    for branch in ("LA", "TA_r", "TA_theta"):
        d = data[data.branch == branch]
        k = d.k_inv_A.to_numpy(); omega = d.omega_peak_mean_rad_ps.to_numpy(); sem = d.omega_peak_replica_SEM_rad_ps.to_numpy()
        resolved = d.resolved_operational_peak.to_numpy(dtype=bool)
        color = colors[branch]
        ax_omega.errorbar(k[resolved], omega[resolved], yerr=sem[resolved], fmt="o", color=color, ms=3.3, lw=1.0, capsize=1.8, label=labels[branch])
        ax_omega.scatter(k[~resolved], omega[~resolved], facecolors="white", edgecolors=color, s=22, linewidths=1.0, zorder=3)
        velocity = d.phase_velocity_A_ps.to_numpy()
        ax_velocity.plot(k[resolved], velocity[resolved], "o", color=color, ms=3.3)
        ax_velocity.scatter(k[~resolved], velocity[~resolved], facecolors="white", edgecolors=color, s=22, linewidths=1.0, zorder=3)
    ax_omega.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel=r"$\omega_{\mathrm{peak}}$ (rad ps$^{-1}$)")
    ax_velocity.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel=r"$\omega_{\mathrm{peak}}/k_n$ (Å ps$^{-1}$), log scale", yscale="log")
    ax_omega.legend(loc="upper left", fontsize=6)
    for ax, label in ((ax_omega, "(a)"), (ax_velocity, "(b)")):
        ax.text(-0.18, 1.02, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
        ax.tick_params(length=3)
    fig.text(0.52, 0.97, "Filled: reproducible spectral peak; open: tentative", ha="center", va="top", fontsize=6.5)
    stem = root / "figures" / "LA_TAr_TAtheta_dispersion_operational_peaks"
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
