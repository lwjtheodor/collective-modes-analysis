#!/usr/bin/env python3
"""Plot signed low-frequency S(k, omega) as one semilog panel per k value.

The source data are three-replica, two-sided Welch spectra.  Panels retain
their absolute spectral scale within each current channel so changes with k
remain visible; channels are deliberately separated into LA, radial TA, and
circumferential TA figures.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
SOURCE = ROOT / "results" / "collective_mode_response" / "88_5L_low_frequency_signed_Skw_CJJ" / "2026-08-19" / "derived_data"
MEAN_CSV = SOURCE / "low_frequency_signed_spectra_ensemble_mean_sem.csv"
REPLICA_CSV = SOURCE / "low_frequency_signed_spectra_per_replica.csv"
OUTROOT = ROOT / "results" / "collective_mode_response" / "88_5L_per_k_semilog_Skw_LA_TAr_TAtheta" / "2026-08-19"
FIGDIR = OUTROOT / "figures"
DATADIR = OUTROOT / "derived_data"

CHANNELS = {
    "S_JzJz_LA": ("LA", "#1769aa"),
    "S_JrJr_TA_r": (r"TA$_r$", "#d55e00"),
    "S_JthetaJtheta_TA_theta": (r"TA$_\\theta$", "#009e73"),
}
NMAX = 20
OMEGA_MAX = 5.0  # rad ps^-1, limited by the source low-frequency product


def _read_mean():
    data = defaultdict(list)
    with MEAN_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            channel = row["channel"]
            n = int(row["n"])
            omega = float(row["omega_rad_ps"])
            if channel in CHANNELS and n <= NMAX and abs(omega) <= OMEGA_MAX:
                data[(channel, n)].append(
                    (omega, float(row["S_mean_arbitrary"]), float(row["S_replica_SEM_arbitrary"]), float(row["k_inv_A"]))
                )
    return data


def _read_replicas():
    data = defaultdict(list)
    with REPLICA_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            channel = row["channel"]
            n = int(row["n"])
            omega = float(row["omega_rad_ps"])
            if channel in CHANNELS and n <= NMAX and abs(omega) <= OMEGA_MAX:
                data[(channel, n, int(row["replicate"]))].append((omega, float(row["S_arbitrary"])))
    return data


def _positive_bounds(mean_data, channel):
    values = [row[1] for n in range(1, NMAX + 1) for row in mean_data[(channel, n)] if row[1] > 0.0]
    if not values:
        raise RuntimeError(f"No positive spectral values for {channel}")
    low = min(values)
    high = max(values)
    return low / 1.8, high * 1.8


def plot_channel(mean_data, replica_data, channel, label, color):
    ylo, yhi = _positive_bounds(mean_data, channel)
    # 5 x 4: one large, readable page that contains every requested n=1...20.
    fig, axes = plt.subplots(5, 4, figsize=(8.0, 10.0), sharex=True, sharey=True)
    axes = axes.ravel()
    for i, n in enumerate(range(1, NMAX + 1)):
        ax = axes[i]
        rows = sorted(mean_data[(channel, n)])
        omega = np.asarray([x[0] for x in rows])
        spectrum = np.asarray([x[1] for x in rows])
        sem = np.asarray([x[2] for x in rows])
        k = rows[0][3]

        for rep in (1, 2, 3):
            rep_rows = sorted(replica_data[(channel, n, rep)])
            rep_omega = np.asarray([x[0] for x in rep_rows])
            rep_spectrum = np.asarray([x[1] for x in rep_rows])
            ax.semilogy(rep_omega, np.maximum(rep_spectrum, ylo * 1e-3), color="0.70", lw=0.35, alpha=0.75, zorder=1)

        lower = np.maximum(spectrum - sem, ylo * 1e-3)
        upper = np.maximum(spectrum + sem, ylo * 1e-3)
        ax.fill_between(omega, lower, upper, color=color, alpha=0.18, linewidth=0, zorder=2)
        ax.semilogy(omega, np.maximum(spectrum, ylo * 1e-3), color=color, lw=0.85, zorder=3)
        ax.axvline(0.0, color="0.25", ls="--", lw=0.45, zorder=0)
        ax.set_xlim(-OMEGA_MAX, OMEGA_MAX)
        ax.set_ylim(ylo, yhi)
        ax.text(0.04, 0.90, rf"$n={n}$, $k={k:.3f}\;\AA^{{-1}}$", transform=ax.transAxes, fontsize=7.0, va="top")
        ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=6.4, pad=1.8)
        ax.tick_params(which="major", length=2.8)
        ax.tick_params(which="minor", length=1.6)
        if n > 16:
            ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)", fontsize=7.2, labelpad=1.5)
        if (n - 1) % 4 == 0:
            ax.set_ylabel(r"$S_{JJ}(k,\omega)$ (arb.)", fontsize=7.2, labelpad=1.5)

    fig.text(0.50, 0.987, rf"$(8,8)$ 5L water: {label} current spectrum, each discrete axial mode", ha="center", va="top", fontsize=10.0)
    fig.text(0.50, 0.965, r"colored: three-replica mean $pm$ SEM; gray: individual replicas; signed frequency", ha="center", va="top", fontsize=7.8, color="0.25")
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.055, top=0.945, hspace=0.18, wspace=0.13)
    stem = f"{label.replace('$', '').replace('\\\\', '').replace('_', '')}_per_k_semilog_Skw_n01_n20"
    for suffix, dpi in (("png", 400), ("tiff", 600), ("pdf", None), ("svg", None)):
        kwargs = {"bbox_inches": "tight"}
        if dpi is not None:
            kwargs["dpi"] = dpi
        fig.savefig(FIGDIR / f"{stem}.{suffix}", **kwargs)
    plt.close(fig)
    return ylo, yhi


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    DATADIR.mkdir(parents=True, exist_ok=True)
    mean_data = _read_mean()
    replica_data = _read_replicas()
    limits = []
    for channel, (label, color) in CHANNELS.items():
        expected = [(channel, n) for n in range(1, NMAX + 1)]
        missing = [key for key in expected if not mean_data[key]]
        if missing:
            raise RuntimeError(f"Missing mean spectra: {missing}")
        ylo, yhi = plot_channel(mean_data, replica_data, channel, label, color)
        limits.append((label, ylo, yhi))

    with (DATADIR / "semilog_plot_limits.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["channel", "n_min", "n_max", "omega_min_rad_ps", "omega_max_rad_ps", "y_min_arbitrary", "y_max_arbitrary", "scale"])
        for label, ylo, yhi in limits:
            writer.writerow([label, 1, NMAX, -OMEGA_MAX, OMEGA_MAX, ylo, yhi, "log10"])
    print(f"Wrote per-k semilog figures to {FIGDIR}")


if __name__ == "__main__":
    main()
