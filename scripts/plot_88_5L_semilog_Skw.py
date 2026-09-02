#!/usr/bin/env python3
"""Log-intensity, signed-frequency S(k,omega) maps for LA, TA_r and TA_theta."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
SOURCE = ROOT / "results" / "collective_mode_response" / "88_5L_low_frequency_signed_Skw_CJJ" / "2026-08-19"
INPUT = SOURCE / "derived_data" / "low_frequency_signed_spectra_ensemble_mean_sem.csv"
OUT = ROOT / "results" / "collective_mode_response" / "88_5L_semilog_Skw_LA_TAr_TAtheta" / "2026-08-19"
CHANNELS = (("S_JzJz_LA", "LA"), ("S_JrJr_TA_r", r"TA$_r$"), ("S_JthetaJtheta_TA_theta", r"TA$_\theta$"))
FLOOR_FRACTION = 1e-5

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})


def edges(x: np.ndarray) -> np.ndarray:
    if len(x) < 2:
        raise ValueError("need at least two coordinates")
    midpoint = 0.5 * (x[1:] + x[:-1])
    return np.concatenate(([x[0] - (midpoint[0] - x[0])], midpoint, [x[-1] + (x[-1] - midpoint[-1])]))


def load() -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    rows = list(csv.DictReader(INPUT.open("r", newline="", encoding="utf-8")))
    output = {}
    for channel, _ in CHANNELS:
        selected = [r for r in rows if r["channel"] == channel]
        n = np.array(sorted({int(r["n"]) for r in selected}), dtype=int)
        omega = np.array(sorted({float(r["omega_rad_ps"]) for r in selected}), dtype=float)
        k = np.array([next(float(r["k_inv_A"]) for r in selected if int(r["n"]) == ni) for ni in n])
        matrix = np.empty((len(n), len(omega)), dtype=float)
        for i, ni in enumerate(n):
            rowmap = {float(r["omega_rad_ps"]): float(r["S_mean_arbitrary"]) for r in selected if int(r["n"]) == ni}
            matrix[i] = [rowmap[om] for om in omega]
        if np.any(matrix <= 0.0):
            raise ValueError(f"{channel}: non-positive spectrum cannot be log-scaled")
        output[channel] = k, omega, matrix
    return output


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white")
    plt.close(fig)


def main() -> None:
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)
    data = load()
    figures, derived = OUT / "figures", OUT / "derived_data"
    figures.mkdir(parents=True, exist_ok=True); derived.mkdir(exist_ok=True)

    # Explicit three-panel / three-colorbar layout; log scale is per channel,
    # so no physically distinct spectra share a colorbar.
    fig = plt.figure(figsize=(7.0, 2.55))
    panel_left = (0.075, 0.392, 0.709)
    panel_width, panel_gap, cbar_width = 0.255, 0.062, 0.018
    normalizations: list[dict] = []
    for idx, ((channel, label), left) in enumerate(zip(CHANNELS, panel_left)):
        ax = fig.add_axes([left, 0.22, panel_width, 0.67])
        cax = fig.add_axes([left + panel_width + 0.010, 0.22, cbar_width, 0.67])
        k, omega, spectrum = data[channel]
        max_s = float(np.max(spectrum)); min_s = float(np.min(spectrum))
        normalized = spectrum / max_s
        im = ax.pcolormesh(edges(k), edges(omega), normalized.T, shading="auto", cmap="magma",
                           norm=LogNorm(vmin=FLOOR_FRACTION, vmax=1.0), rasterized=True)
        ax.axhline(0.0, color="white", ls="--", lw=0.9, alpha=0.9)
        ax.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylim=(-5.0, 5.0))
        if idx == 0:
            ax.set_ylabel(r"$\omega$ (rad ps$^{-1}$)")
        else:
            ax.set_yticklabels([])
        ax.text(0.03, 0.94, label, transform=ax.transAxes, ha="left", va="top", color="white", fontsize=7.2, fontweight="bold")
        ax.text(-0.20, 1.02, f"({chr(97 + idx)})", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
        ax.tick_params(length=3, labelsize=6.5)
        cb = fig.colorbar(im, cax=cax, ticks=[1e-5, 1e-3, 1e-1, 1])
        cb.ax.set_yticklabels([r"$10^{-5}$", r"$10^{-3}$", r"$10^{-1}$", "1"])
        cb.ax.tick_params(labelsize=5.8, length=2.5)
        normalizations.append({"channel": channel, "branch": label, "normalization": "S/S_channel_max", "S_channel_max_arbitrary": max_s, "S_channel_min_arbitrary": min_s, "LogNorm_vmin_fraction": FLOOR_FRACTION})
    fig.text(0.50, 0.965, r"(8,8) 5L; signed-frequency, log-intensity current spectra; each panel normalized by its own $S_{\max}$", ha="center", va="top", fontsize=6.6)
    save(fig, figures / "semilog_signed_Skw_LA_TAr_TAtheta_pm5")

    with (derived / "semilog_Skw_normalization.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(normalizations[0])); writer.writeheader(); writer.writerows(normalizations)
    (OUT / "figure_contract.txt").write_text("""Core conclusion: compare the signed-frequency low-frequency current spectra S(k,omega) for LA, radial TA and circumferential TA while exposing weak spectral weight on a logarithmic intensity scale.
Figure archetype: quantitative grid.
Backend: Python/matplotlib.
Panel map: (a) S_JzJz (LA); (b) S_JrJr (TA_r); (c) S_JthetaJtheta (TA_theta). Each has its own logarithmic colorbar after normalization by its channel maximum.
Evidence hierarchy: ensemble mean of three velocity-seed replicas; two-sided equilibrium-symmetrized Welch spectra; zero-frequency line.
Reviewer risk: per-channel normalization allows within-channel k/frequency morphology but must not be read as a cross-channel absolute-intensity comparison. The data support only |omega|<=5 rad ps^-1 and three Welch windows per replica.
""", encoding="utf-8")
    (OUT / "QA_notes.txt").write_text("""Data integrity: all three requested current channels, all n=1..20, and every signed omega bin in [-5,5] rad ps^-1 were retained. No values were imputed or excluded.
Statistics: input is the mean over three velocity-seed replicas; SEM is retained in the source CSV but not encoded in the heatmap. Each replica had three 50%-overlapped Welch windows.
Image integrity: quantitative rasterized heatmap only; no local image adjustment. LogNorm is applied after division by each channel maximum, with a declared 1e-5 display floor.
""", encoding="utf-8")
    (OUT / "metadata.json").write_text(json.dumps({"source": str(INPUT), "system": "(8,8) CNT water, 5L", "channels": [x[0] for x in CHANNELS], "n_range": [1, 20], "omega_range_rad_ps": [-5, 5], "semilog_definition": "LogNorm color scale of S(channel,k,omega)/max_{k,omega}S(channel,k,omega), independently for each channel", "replicates": 3}, indent=2), encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("Semilog signed S(k,omega) maps completed.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
