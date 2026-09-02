#!/usr/bin/env python3
"""Dedicated low-frequency signed spectrum plots for circumferential TA_theta."""
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
CHANNEL = "S_JthetaJtheta_TA_theta"


def load() -> dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, float]]:
    bins: dict[int, list[tuple[float, float, float, float]]] = defaultdict(list)
    with CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["channel"] == CHANNEL:
                bins[int(row["n"])].append((float(row["omega_rad_ps"]), float(row["S_mean_arbitrary"]),
                                              float(row["S_replica_SEM_arbitrary"]), float(row["k_inv_A"])))
    return {n: (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]),
                np.array([r[2] for r in rows]), rows[0][3]) for n, rows in bins.items()}


def save(fig: plt.Figure, stem: str) -> None:
    base = OUT / "figures" / stem
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(base.with_suffix(ext), dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    data = load()
    modes = sorted(data)
    omega = data[modes[0]][0]
    k = np.array([data[n][3] for n in modes])
    mean = np.column_stack([data[n][1] for n in modes])
    keep = np.abs(omega) <= 5.0
    rel = mean[keep] / np.maximum(mean[keep].max(axis=0, keepdims=True), 1e-300)

    fig = plt.figure(figsize=(4.0, 3.25))
    ax = fig.add_axes((0.15, 0.16, 0.68, 0.75))
    im = ax.pcolormesh(k, omega[keep], np.log10(np.maximum(rel, 1e-6)), shading="auto", cmap="magma", vmin=-6, vmax=0)
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlabel(r"$k$ (Å$^{-1}$)")
    ax.set_ylabel(r"$\omega$ (rad ps$^{-1}$)")
    ax.set_xlim(k[0], k[-1])
    ax.set_ylim(-5, 5)
    ax.text(-0.18, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.text(0.04, 0.92, r"$S_{J_\theta J_\theta}(k,\omega)$", transform=ax.transAxes, color="#42949E", fontweight="bold")
    cax = fig.add_axes((0.87, 0.16, 0.025, 0.75))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$\log_{10}[S/S_{\max}(k)]$")
    save(fig, "TA_theta_low_frequency_signed_Skw")

    fig = plt.figure(figsize=(5.5, 4.3))
    boxes = [(0.11, 0.58, 0.24, 0.29), (0.42, 0.58, 0.24, 0.29), (0.73, 0.58, 0.24, 0.29),
             (0.11, 0.14, 0.24, 0.29), (0.42, 0.14, 0.24, 0.29), (0.73, 0.14, 0.24, 0.29)]
    for idx, (box, n) in enumerate(zip(boxes, (1, 2, 4, 8, 12, 16))):
        ax = fig.add_axes(box)
        om, values, sem, kn = data[n]
        vals = values / values.max()
        err = sem / values.max()
        ax.fill_between(om, np.maximum(vals - err, 0), vals + err, color="#42949E", alpha=0.22, lw=0)
        ax.plot(om, vals, color="#42949E", lw=1.1)
        ax.axvline(0, color="#767676", lw=0.8, zorder=0)
        ax.set_xlim(-5, 5)
        ax.set_ylim(-0.02, 1.05)
        ax.text(-0.15, 1.04, f"({'abcdef'[idx]})", transform=ax.transAxes, fontweight="bold", fontsize=9)
        ax.text(0.05, 0.88, fr"$n={n}$, $k={kn:.3f}$", transform=ax.transAxes)
        if idx >= 3:
            ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
        else:
            ax.set_xticklabels([])
        if idx % 3 == 0:
            ax.set_ylabel(r"$S/S_{\max}$")
        else:
            ax.set_yticklabels([])
    save(fig, "TA_theta_low_frequency_signed_linecuts")

    # Dedicated magnification of the central low-frequency sector.  The source
    # spectra are unchanged; only the displayed range and per-k display scale
    # are restricted to |omega| <= 2 rad ps^-1.
    zoom = np.abs(omega) <= 2.0
    zoom_rel = mean[zoom] / np.maximum(mean[zoom].max(axis=0, keepdims=True), 1e-300)
    fig = plt.figure(figsize=(5.2, 4.0))
    ax = fig.add_axes((0.14, 0.15, 0.67, 0.77))
    im = ax.pcolormesh(k, omega[zoom], np.log10(np.maximum(zoom_rel, 1e-6)), shading="auto", cmap="magma", vmin=-6, vmax=0)
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlabel(r"$k$ (Å$^{-1}$)")
    ax.set_ylabel(r"$\omega$ (rad ps$^{-1}$)")
    ax.set_xlim(k[0], k[-1])
    ax.set_ylim(-2, 2)
    ax.text(-0.16, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.text(0.04, 0.92, r"$S_{J_\theta J_\theta}(k,\omega)$; $|\omega|\leq2$", transform=ax.transAxes, color="#42949E", fontweight="bold")
    cax = fig.add_axes((0.86, 0.15, 0.025, 0.77))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$\log_{10}[S/S_{\max}^{|\omega|\leq2}(k)]$")
    save(fig, "TA_theta_signed_Skw_zoom_pm2")

    fig = plt.figure(figsize=(5.5, 3.9))
    boxes = [(0.12, 0.55, 0.35, 0.32), (0.60, 0.55, 0.35, 0.32),
             (0.12, 0.13, 0.35, 0.32), (0.60, 0.13, 0.35, 0.32)]
    for idx, (box, n) in enumerate(zip(boxes, (1, 2, 4, 8))):
        ax = fig.add_axes(box)
        om, values, sem, kn = data[n]
        local = np.abs(om) <= 2.0
        scale = values[local].max()
        vals, err = values / scale, sem / scale
        ax.fill_between(om[local], np.maximum(vals[local] - err[local], 0), vals[local] + err[local], color="#42949E", alpha=0.22, lw=0)
        ax.plot(om[local], vals[local], color="#42949E", lw=1.1)
        ax.axvline(0, color="#767676", lw=0.8, zorder=0)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-0.02, 1.05)
        ax.text(-0.15, 1.04, f"({'abcd'[idx]})", transform=ax.transAxes, fontweight="bold", fontsize=9)
        ax.text(0.05, 0.87, fr"$n={n}$, $k={kn:.3f}$", transform=ax.transAxes)
        if idx >= 2:
            ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel(r"$S/S_{\max}^{|\omega|\leq2}$")
        else:
            ax.set_yticklabels([])
    save(fig, "TA_theta_signed_linecuts_zoom_pm2")

    ultra = np.abs(omega) <= 0.5
    ultra_rel = mean[ultra] / np.maximum(mean[ultra].max(axis=0, keepdims=True), 1e-300)
    fig = plt.figure(figsize=(5.2, 4.0))
    ax = fig.add_axes((0.14, 0.15, 0.67, 0.77))
    im = ax.pcolormesh(k, omega[ultra], np.log10(np.maximum(ultra_rel, 1e-6)), shading="auto", cmap="magma", vmin=-6, vmax=0)
    ax.axhline(0, color="white", lw=0.8)
    ax.set_xlabel(r"$k$ (Å$^{-1}$)")
    ax.set_ylabel(r"$\omega$ (rad ps$^{-1}$)")
    ax.set_xlim(k[0], k[-1])
    ax.set_ylim(-0.5, 0.5)
    ax.text(-0.16, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.text(0.04, 0.92, r"$S_{J_\theta J_\theta}(k,\omega)$; $|\omega|\leq0.5$", transform=ax.transAxes, color="#42949E", fontweight="bold")
    cax = fig.add_axes((0.86, 0.15, 0.025, 0.77))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$\log_{10}[S/S_{\max}^{|\omega|\leq0.5}(k)]$")
    save(fig, "TA_theta_signed_Skw_zoom_pm05")

    fig = plt.figure(figsize=(5.5, 3.9))
    boxes = [(0.12, 0.55, 0.35, 0.32), (0.60, 0.55, 0.35, 0.32),
             (0.12, 0.13, 0.35, 0.32), (0.60, 0.13, 0.35, 0.32)]
    for idx, (box, n) in enumerate(zip(boxes, (1, 2, 4, 8))):
        ax = fig.add_axes(box)
        om, values, sem, kn = data[n]
        local = np.abs(om) <= 0.5
        scale = values[local].max()
        vals, err = values / scale, sem / scale
        ax.fill_between(om[local], np.maximum(vals[local] - err[local], 0), vals[local] + err[local], color="#42949E", alpha=0.22, lw=0)
        ax.plot(om[local], vals[local], color="#42949E", lw=1.1)
        ax.axvline(0, color="#767676", lw=0.8, zorder=0)
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.02, 1.05)
        ax.text(-0.15, 1.04, f"({'abcd'[idx]})", transform=ax.transAxes, fontweight="bold", fontsize=9)
        ax.text(0.05, 0.87, fr"$n={n}$, $k={kn:.3f}$", transform=ax.transAxes)
        if idx >= 2:
            ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
        else:
            ax.set_xticklabels([])
        if idx % 2 == 0:
            ax.set_ylabel(r"$S/S_{\max}^{|\omega|\leq0.5}$")
        else:
            ax.set_yticklabels([])
    save(fig, "TA_theta_signed_linecuts_zoom_pm05")


if __name__ == "__main__":
    main()
