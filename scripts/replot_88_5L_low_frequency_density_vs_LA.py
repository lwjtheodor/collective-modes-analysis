#!/usr/bin/env python3
"""Replot low-k density and LA-current cuts from signed-spectrum source data."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
                     "axes.linewidth": 1.0, "axes.spines.right": False,
                     "axes.spines.top": False, "legend.frameon": False,
                     "xtick.direction": "out", "ytick.direction": "out"})

OUT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes\results\collective_mode_response\88_5L_low_frequency_signed_Skw_CJJ\2026-08-19")
CSV = OUT / "derived_data" / "low_frequency_signed_spectra_ensemble_mean_sem.csv"


def load() -> dict[tuple[str, int], tuple[np.ndarray, np.ndarray, np.ndarray, float]]:
    values: dict[tuple[str, int], list[tuple[float, float, float, float]]] = defaultdict(list)
    with CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            values[(row["channel"], int(row["n"]))].append(
                (float(row["omega_rad_ps"]), float(row["S_mean_arbitrary"]),
                 float(row["S_replica_SEM_arbitrary"]), float(row["k_inv_A"]))
            )
    return {key: (np.array([x[0] for x in rows]), np.array([x[1] for x in rows]),
                  np.array([x[2] for x in rows]), rows[0][3]) for key, rows in values.items()}


def main() -> None:
    data = load()
    fig = plt.figure(figsize=(5.5, 4.3))
    boxes = [(0.12, 0.57, 0.35, 0.30), (0.60, 0.57, 0.35, 0.30),
             (0.12, 0.13, 0.35, 0.30), (0.60, 0.13, 0.35, 0.30)]
    modes = (1, 4, 8, 12)
    diagnostics = []
    for idx, (box, n) in enumerate(zip(boxes, modes)):
        ax = fig.add_axes(box)
        omega, rho, rho_sem, k = data[("S_rhorho", n)]
        _, la, la_sem, _ = data[("S_JzJz_LA", n)]
        rho_norm = rho / rho.max()
        la_norm = la / la.max()
        ax.plot(omega, rho_norm, color="#272727", lw=1.1, label=r"$S_{\rho\rho}$")
        ax.plot(omega, la_norm, color="#0F4D92", lw=1.1, label=r"$S_{J_zJ_z}$")
        ax.axvline(0, color="#767676", lw=0.8, zorder=0)
        ax.set_xlim(-5, 5)
        ax.set_ylim(-0.02, 1.05)
        ax.text(-0.16, 1.04, f"({'abcd'[idx]})", transform=ax.transAxes, fontweight="bold", fontsize=9)
        ax.text(0.04, 0.88, fr"$n={n}$, $k={k:.3f}$ Å$^{{-1}}$", transform=ax.transAxes)
        if idx >= 2:
            ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel("normalized spectrum")
        else:
            ax.set_yticklabels([])
        if idx == 0:
            ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.77), fontsize=6, handlelength=1.5)
        positive = omega > 0.15
        peak = np.argmax(np.where(positive, la, -np.inf))
        zero = np.argmin(abs(omega))
        diagnostics.append({"n": n, "k_inv_A": k, "density_zero_over_max": rho_norm[zero],
                            "LA_positive_sideband_omega_rad_ps": omega[peak],
                            "LA_positive_sideband_over_max": la_norm[peak]})
    base = OUT / "figures" / "low_frequency_density_vs_LA_current_linecuts"
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(base.with_suffix(ext), dpi=300, facecolor="white")
    plt.close(fig)
    with (OUT / "derived_data" / "density_vs_LA_lowk_feature_diagnostic.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(diagnostics[0]))
        writer.writeheader()
        writer.writerows(diagnostics)


if __name__ == "__main__":
    main()
