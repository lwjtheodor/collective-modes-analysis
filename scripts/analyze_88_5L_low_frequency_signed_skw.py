#!/usr/bin/env python3
"""Low-frequency, signed-frequency density/current spectra for (8,8) 5L.

The scalar density dynamic structure factor S_rhorho(k, omega) exposes a
central Rayleigh-like (thermal/diffusive) response and Brillouin-like sound
sidebands.  The three directional current spectra are S_JzJz (LA),
S_JrJr (radial TA), and S_JthetaJtheta (circumferential TA).  A two-sided,
equilibrium-symmetrized spectrum is stored so that positive and negative
frequency branches are displayed explicitly.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Figure contract: a quantitative grid whose hero is S_rhorho(k,omega), where
# central (Rayleigh-like) and sideband (Brillouin-like) spectral weight can be
# separated from the directional current responses without conflating TA_r and
# TA_theta.  The source data are all three verified oxygen 10-fs trajectories.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 1.0,
    "axes.spines.right": False, "axes.spines.top": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "legend.frameon": False,
})

ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
INPUT_ROOT = Path(r"H:\gcmc_explore\translational_anomaly\08_viscosity_friction_length_scaling\04_analysis\offline_frequency_viscosity_20260803\raw_case_directories")
OUT = ROOT / "results" / "collective_mode_response" / "88_5L_low_frequency_signed_Skw_CJJ" / "2026-08-19"
NMAX = 20
NPERSEG = 16384
OMEGA_MAX = 5.0  # rad ps^-1; deliberately low-frequency-only output
CHANNELS = ("S_rhorho", "S_JzJz_LA", "S_JrJr_TA_r", "S_JthetaJtheta_TA_theta")
COLORS = ("#272727", "#0F4D92", "#B64342", "#42949E")


def locate_dumps() -> list[Path]:
    paths = []
    for rep in (1, 2, 3):
        path = INPUT_ROOT / f"NVT20ns_5xL_8_8_RH75_N665_rep{rep}_20260719" / f"nvt20ns_8_8_RH75_5L_rep{rep}_oxygen_10fs_1ns.dump"
        if not path.is_file():
            raise FileNotFoundError(path)
        paths.append(path)
    return paths


def read_modes(path: Path, nmax: int, capacity: int = 40000) -> tuple[np.ndarray, dict]:
    """Read rho, Jz, Jr and Jtheta for n=1..nmax from one oxygen dump."""
    out = np.empty((capacity, nmax, 4), dtype=np.complex128)
    frame = 0
    steps: list[int] = []
    box = None
    centre = None
    n_oxygen = None
    required = {"id", "type", "x", "y", "z", "vx", "vy", "vz"}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        while True:
            marker = fh.readline()
            if marker == "":
                break
            if marker != "ITEM: TIMESTEP\n":
                raise ValueError(f"{path}: malformed timestep marker")
            step = int(fh.readline())
            if fh.readline() != "ITEM: NUMBER OF ATOMS\n":
                raise ValueError(f"{path}: missing atom count")
            natom = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"{path}: missing box bounds")
            bounds = np.array([list(map(float, fh.readline().split()[:2])) for _ in range(3)])
            header = fh.readline().split()[2:]
            col = {name: i for i, name in enumerate(header)}
            missing = required - set(col)
            if missing:
                raise ValueError(f"{path}: missing {sorted(missing)}")
            xyz = []
            vel = []
            for _ in range(natom):
                row = fh.readline().split()
                if int(float(row[col["type"]])) == 3:
                    xyz.append((float(row[col["x"]]), float(row[col["y"]]), float(row[col["z"]])))
                    vel.append((float(row[col["vx"]]), float(row[col["vy"]]), float(row[col["vz"]])))
            xyz_a = np.asarray(xyz, dtype=float)
            vel_a = np.asarray(vel, dtype=float)
            frame_box = bounds[:, 1] - bounds[:, 0]
            if box is None:
                box = frame_box
                centre = bounds.mean(axis=1)
                n_oxygen = len(xyz_a)
            elif len(xyz_a) != n_oxygen or not np.allclose(frame_box, box, rtol=0, atol=1e-9):
                raise ValueError(f"{path}: oxygen count or box changed at {step}")
            if frame >= capacity:
                raise ValueError(f"{path}: increase capacity above {capacity}")
            # Remove instantaneous translational O-COM drift in Cartesian space
            # before cylindrical projections.
            vel_a -= vel_a.mean(axis=0, keepdims=True)
            dxy = xyz_a[:, :2] - centre[:2]
            er = dxy / np.maximum(np.hypot(dxy[:, 0], dxy[:, 1])[:, None], 1e-12)
            vr = np.sum(vel_a[:, :2] * er, axis=1)
            vtheta = -vel_a[:, 0] * er[:, 1] + vel_a[:, 1] * er[:, 0]
            phase1 = np.exp(2j * np.pi * xyz_a[:, 2] / box[2])
            phase = phase1.copy()
            for mi in range(nmax):
                out[frame, mi, 0] = phase.sum()
                out[frame, mi, 1] = np.dot(vel_a[:, 2], phase)
                out[frame, mi, 2] = np.dot(vr, phase)
                out[frame, mi, 3] = np.dot(vtheta, phase)
                phase *= phase1
            frame += 1
            steps.append(step)
    steps_a = np.asarray(steps)
    if frame < NPERSEG:
        raise ValueError(f"{path}: {frame} frames < nperseg={NPERSEG}")
    if not np.all(np.diff(steps_a) == np.diff(steps_a)[0]):
        raise ValueError(f"{path}: nonuniform dump cadence")
    dt_ps = float(np.diff(steps_a)[0]) * 0.0005
    meta = {
        "input": str(path), "n_frames": int(frame), "n_oxygen": int(n_oxygen),
        "dt_ps": dt_ps, "duration_ps": float((frame - 1) * dt_ps),
        "box_A": box.tolist(), "axis_xy_A": centre[:2].tolist(),
        "fields_verified": ["x", "y", "z", "vx", "vy", "vz"],
    }
    return out[:frame], meta


def symmetric_two_sided_welch(series: np.ndarray, dt_ps: float) -> tuple[np.ndarray, np.ndarray, int]:
    """Complex Welch periodogram, symmetrized as S(omega)=S(-omega)."""
    starts = np.arange(0, len(series) - NPERSEG + 1, NPERSEG // 2)
    if len(starts) < 3:
        raise ValueError("fewer than three Welch windows")
    window = np.hanning(NPERSEG)[:, None, None]
    norm = np.sum(window[:, 0, 0] ** 2)
    accum = np.zeros((NPERSEG, series.shape[1], series.shape[2]), dtype=float)
    for start in starts:
        segment = series[start:start + NPERSEG]
        segment = segment - segment.mean(axis=0, keepdims=True)
        ff = np.fft.fft(segment * window, axis=0)
        accum += np.abs(ff) ** 2 / norm
    accum /= len(starts)
    negative_index = (-np.arange(NPERSEG)) % NPERSEG
    accum = 0.5 * (accum + accum[negative_index])
    omega = 2 * np.pi * np.fft.fftfreq(NPERSEG, d=dt_ps)
    return np.fft.fftshift(omega), np.fft.fftshift(accum, axes=0), int(len(starts))


def write_csv(path: Path, omega: np.ndarray, spectra: np.ndarray, k: np.ndarray, replicas: bool) -> None:
    """Write low-frequency signed spectra; spectra either [rep,w,n,c] or [w,n,c,2]."""
    keep = np.abs(omega) <= OMEGA_MAX + 1e-12
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        if replicas:
            writer.writerow(["replicate", "channel", "n", "k_inv_A", "omega_rad_ps", "S_arbitrary", "S_over_mode_max"])
            for rep in range(spectra.shape[0]):
                for ci, channel in enumerate(CHANNELS):
                    for mi in range(spectra.shape[2]):
                        vals = spectra[rep, :, mi, ci]
                        vmax = vals[keep].max()
                        for wi in np.flatnonzero(keep):
                            writer.writerow([rep + 1, channel, mi + 1, k[mi], omega[wi], vals[wi], vals[wi] / vmax])
        else:
            writer.writerow(["channel", "n", "k_inv_A", "omega_rad_ps", "S_mean_arbitrary", "S_replica_SEM_arbitrary", "S_over_mode_max"])
            for ci, channel in enumerate(CHANNELS):
                for mi in range(spectra.shape[1]):
                    vals = spectra[:, mi, ci, 0]
                    sem = spectra[:, mi, ci, 1]
                    vmax = vals[keep].max()
                    for wi in np.flatnonzero(keep):
                        writer.writerow([channel, mi + 1, k[mi], omega[wi], vals[wi], sem[wi], vals[wi] / vmax])


def heatmap(ax, k: np.ndarray, omega: np.ndarray, matrix: np.ndarray, label: str, color: str, letter: str):
    keep = np.abs(omega) <= OMEGA_MAX + 1e-12
    selected = matrix[keep].T
    rel = selected / np.maximum(selected.max(axis=1, keepdims=True), 1e-300)
    image = ax.pcolormesh(k, omega[keep], np.log10(np.maximum(rel.T, 1e-6)), shading="auto", cmap="magma", vmin=-6, vmax=0)
    ax.axhline(0, color="white", lw=0.8, alpha=0.9)
    ax.set_xlim(k[0], k[-1])
    ax.set_ylim(-OMEGA_MAX, OMEGA_MAX)
    ax.set_ylabel(r"$\omega$ (rad ps$^{-1}$)")
    ax.text(-0.16, 1.04, f"({letter})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.text(0.03, 0.91, label, transform=ax.transAxes, color=color, fontweight="bold")
    return image


def plot_figures(omega: np.ndarray, ensemble: np.ndarray, k: np.ndarray, out: Path) -> None:
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    mean = ensemble[..., 0]
    fig = plt.figure(figsize=(7.0, 5.2))
    boxes = [(0.10, 0.57, 0.36, 0.34), (0.56, 0.57, 0.36, 0.34), (0.10, 0.12, 0.36, 0.34), (0.56, 0.12, 0.36, 0.34)]
    axes = [fig.add_axes(b) for b in boxes]
    ims = []
    labels = (r"$S_{\rho\rho}$", r"$S_{J_zJ_z}$ (LA)", r"$S_{J_rJ_r}$ (TA$_r$)", r"$S_{J_\theta J_\theta}$ (TA$_\theta$)")
    for ai, (ax, label, color) in enumerate(zip(axes, labels, COLORS)):
        ims.append(heatmap(ax, k, omega, mean[:, :, ai], label, color, "abcd"[ai]))
        if ai >= 2:
            ax.set_xlabel(r"$k$ (Å$^{-1}$)")
        else:
            ax.set_xticklabels([])
    cax = fig.add_axes((0.94, 0.12, 0.016, 0.79))
    cb = fig.colorbar(ims[0], cax=cax)
    cb.set_label(r"$\log_{10}[S/S_{\max}(k)]$")
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(figures / f"low_frequency_signed_Skw_density_and_currents{ext}", dpi=300, facecolor="white")
    plt.close(fig)

    # Density cuts visualize the central and symmetric side-band morphology.
    fig = plt.figure(figsize=(5.5, 2.6))
    ax = fig.add_axes((0.12, 0.18, 0.70, 0.72))
    keep = np.abs(omega) <= OMEGA_MAX + 1e-12
    for n, color in zip((1, 4, 8, 12, 16, 20), ("#272727", "#0F4D92", "#3775BA", "#42949E", "#B64342", "#9A4D8E")):
        vals = mean[:, n - 1, 0]
        ax.plot(omega[keep], vals[keep] / vals[keep].max(), color=color, lw=1.0, label=fr"$n={n}$")
    ax.axvline(0, color="#767676", lw=0.9)
    ax.set_xlim(-OMEGA_MAX, OMEGA_MAX)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
    ax.set_ylabel(r"$S_{\rho\rho}/S_{\max}$")
    ax.text(-0.14, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.legend(ncol=2, loc="upper right", fontsize=6, handlelength=1.4)
    ax.text(0.02, 0.09, "central: Rayleigh-like", transform=ax.transAxes, color="#4D4D4D")
    ax.text(0.55, 0.72, "sides: Brillouin-like", transform=ax.transAxes, color="#4D4D4D")
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(figures / f"low_frequency_signed_density_linecuts{ext}", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT / "derived_data").mkdir()
    dumps = locate_dumps()
    per_rep = []
    audits = []
    omega = None
    windows = []
    for path in dumps:
        modes, audit = read_modes(path, NMAX)
        current_omega, spec, nwin = symmetric_two_sided_welch(modes, audit["dt_ps"])
        if omega is None:
            omega = current_omega
        elif not np.allclose(omega, current_omega):
            raise ValueError("frequency grids differ across replicas")
        per_rep.append(spec)
        audit["welch_windows"] = nwin
        audits.append(audit)
        windows.append(nwin)
    spectra = np.asarray(per_rep)
    k = 2 * np.pi * np.arange(1, NMAX + 1) / audits[0]["box_A"][2]
    mean = spectra.mean(axis=0)
    sem = spectra.std(axis=0, ddof=1) / math.sqrt(spectra.shape[0])
    ensemble = np.stack((mean, sem), axis=-1)
    derived = OUT / "derived_data"
    write_csv(derived / "low_frequency_signed_spectra_per_replica.csv", omega, spectra, k, replicas=True)
    write_csv(derived / "low_frequency_signed_spectra_ensemble_mean_sem.csv", omega, ensemble, k, replicas=False)
    plot_figures(omega, ensemble, k, OUT)
    metadata = {
        "system": "(8,8) 5L", "source": "three oxygen-only 10 fs trajectory segments",
        "n_replicas": 3, "nmax": NMAX, "omega_window_rad_ps": [-OMEGA_MAX, OMEGA_MAX],
        "welch_nperseg": NPERSEG, "welch_overlap": 0.5, "welch_windows": windows,
        "spectral_convention": "two-sided complex Welch auto-spectrum, equilibrium symmetrized S(omega)=S(-omega)",
        "density_channel": "S_rhorho exposes central Rayleigh-like and sideband Brillouin-like response",
        "current_channels": {"LA": "S_JzJz", "TA_r": "S_JrJr", "TA_theta": "S_JthetaJtheta"},
        "normalization": "heatmaps use each (channel,n) low-frequency maximum; raw arbitrary-unit spectra are in CSV",
        "audits": audits,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("Signed low-frequency density/current spectra finished successfully.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
