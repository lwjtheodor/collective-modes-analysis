#!/usr/bin/env python3
"""All-mode longitudinal/transverse current response of the (8,8) 5L bursts.

The input is the three oxygen-only 10 fs LAMMPS dumps. Their actual duration is
audited from frame cadence rather than inferred from the historical filename. For every
allowed axial integer mode n=1..nmax this script constructs

    J_a(k_n,t) = sum_i [v_{i,a}(t)-<v_a(t)>_O] exp(i k_n z_i(t)),

where a=z (LA), r (radial TA), or theta (circumferential TA), with r and theta
defined from the fixed CNT/box axis.  C_JJ is the real, all-time-origin
autocorrelation of the temporally centred complex current.  Radial and
circumferential TA responses are never averaged together.
The script preserves the per-replica curves, writes the ensemble mean and
replica SEM separately, and labels the dispersion points as *operational
spectral peaks* rather than assuming that every high-k feature is a resolved
acoustic mode.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ModuleNotFoundError:
    NUMBA_AVAILABLE = False

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# Required publication-style, editable-SVG settings.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7, "axes.linewidth": 1.0,
    "axes.spines.right": False, "axes.spines.top": False,
    "legend.frameon": False, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
})

LA_COLOR = "#0F4D92"
TAR_COLOR = "#B64342"
TATHETA_COLOR = "#42949E"
COMPONENTS = ("LA", "TA_r", "TA_theta")


if NUMBA_AVAILABLE:
    @njit(cache=True)
    def exact_modal_projection(z: np.ndarray, velocity: np.ndarray, e_r: np.ndarray, lz: float, nmax: int) -> np.ndarray:
        """Compiled exact sum_i v_i exp(i n 2pi z_i/Lz)."""
        output = np.zeros((nmax, 3), dtype=np.complex128)
        for i in range(z.size):
            phase_1 = np.cos(2.0 * np.pi * z[i] / lz) + 1j * np.sin(2.0 * np.pi * z[i] / lz)
            phase = phase_1
            vr = velocity[i, 0] * e_r[i, 0] + velocity[i, 1] * e_r[i, 1]
            vtheta = -velocity[i, 0] * e_r[i, 1] + velocity[i, 1] * e_r[i, 0]
            for mi in range(nmax):
                output[mi, 0] += velocity[i, 2] * phase
                output[mi, 1] += vr * phase
                output[mi, 2] += vtheta * phase
                phase *= phase_1
        return output
else:
    def exact_modal_projection(z: np.ndarray, velocity: np.ndarray, e_r: np.ndarray, lz: float, nmax: int) -> np.ndarray:
        """NumPy fallback when the selected analysis environment lacks numba."""
        phase_1 = np.exp(2j * np.pi * z / lz)
        phase = phase_1.copy()
        vr = velocity[:, 0] * e_r[:, 0] + velocity[:, 1] * e_r[:, 1]
        vtheta = -velocity[:, 0] * e_r[:, 1] + velocity[:, 1] * e_r[:, 0]
        output = np.empty((nmax, 3), dtype=np.complex128)
        for mi in range(nmax):
            output[mi, 0] = np.dot(velocity[:, 2], phase)
            output[mi, 1] = np.dot(vr, phase)
            output[mi, 2] = np.dot(vtheta, phase)
            phase *= phase_1
        return output


def locate_dumps(root: Path, dump_template: str) -> list[Path]:
    paths = []
    for rep in (1, 2, 3):
        relative = dump_template.format(rep=rep)
        found = sorted(root.glob(relative))
        if len(found) != 1:
            raise FileNotFoundError(f"expected exactly one dump for rep{rep} from {root / relative}; found {found}")
        paths.append(found[0])
    return paths


def read_modal_currents(path: Path, nmax: int, expected_frames: int = 100_001, oxygen_type: int = 3) -> tuple[np.ndarray, dict]:
    """Stream one LAMMPS dump and return J_z, J_r, J_theta [frame, mode, component]."""
    # The verified 1 ns burst has 100001 frames.  Preallocation avoids the
    # multi-gigabyte Python-object overhead of a list of 100001 tiny arrays.
    records = np.empty((expected_frames, nmax, 3), dtype=np.complex128)
    frame_index = 0
    steps: list[int] = []
    box = center = None
    nwater = None
    required = {"id", "type", "x", "y", "z", "vx", "vy", "vz"}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        while True:
            marker = fh.readline()
            if marker == "":
                break
            if marker != "ITEM: TIMESTEP\n":
                raise ValueError(f"{path}: invalid frame marker {marker!r}")
            step = int(fh.readline())
            if fh.readline() != "ITEM: NUMBER OF ATOMS\n":
                raise ValueError(f"{path}: missing NUMBER OF ATOMS marker")
            natom = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"{path}: missing BOX BOUNDS marker")
            bounds = np.array([list(map(float, fh.readline().split()[:2])) for _ in range(3)])
            header = fh.readline().split()[2:]
            col = {name: idx for idx, name in enumerate(header)}
            missing = required - set(col)
            if missing:
                raise ValueError(f"{path}: lacks required atom columns {sorted(missing)}")

            xyz_rows: list[tuple[float, float, float]] = []
            vel_rows: list[tuple[float, float, float]] = []
            for _ in range(natom):
                row = fh.readline().split()
                if int(float(row[col["type"]])) != oxygen_type:
                    continue
                xyz_rows.append((float(row[col["x"]]), float(row[col["y"]]), float(row[col["z"]])))
                vel_rows.append((float(row[col["vx"]]), float(row[col["vy"]]), float(row[col["vz"]])))
            xyz = np.asarray(xyz_rows, dtype=np.float64)
            vel = np.asarray(vel_rows, dtype=np.float64)
            frame_box = bounds[:, 1] - bounds[:, 0]
            if box is None:
                box = frame_box
                center = bounds.mean(axis=1)
                nwater = len(xyz)
            elif len(xyz) != nwater or not np.allclose(frame_box, box, rtol=0, atol=1e-9):
                raise ValueError(f"{path}: variable oxygen count or box at timestep {step}")
            if nwater == 0:
                raise ValueError(f"{path}: no oxygen atoms (type {oxygen_type}) in timestep {step}")

            # Instantaneous O-COM subtraction is Cartesian and occurs before the
            # radial/circumferential projection, avoiding a global-drift artefact.
            vel -= vel.mean(axis=0, keepdims=True)
            dxy = xyz[:, :2] - center[:2]
            radius = np.hypot(dxy[:, 0], dxy[:, 1])
            e_r = dxy / np.maximum(radius[:, None], 1e-12)
            if frame_index >= expected_frames:
                raise ValueError(f"{path}: more than expected {expected_frames} frames")
            records[frame_index] = exact_modal_projection(xyz[:, 2], vel, e_r, box[2], nmax)
            frame_index += 1
            steps.append(step)

    series = records[:frame_index]
    steps_a = np.asarray(steps, dtype=np.int64)
    if len(steps_a) < 3:
        raise ValueError(f"{path}: too few frames ({len(steps_a)})")
    step_delta = np.diff(steps_a)
    if not np.all(step_delta == step_delta[0]):
        raise ValueError(f"{path}: nonuniform dump timesteps")
    # The production inputs use 0.5 fs integration and 20-step dump cadence.
    dt_ps = float(step_delta[0]) * 0.0005
    meta = {
        "input_file": str(path), "n_frames": int(len(series)), "n_water": int(nwater),
        "timestep_first": int(steps_a[0]), "timestep_last": int(steps_a[-1]),
        "dump_step_interval": int(step_delta[0]), "dt_ps": dt_ps,
        "duration_ps": float((len(series) - 1) * dt_ps),
        "box_A": [float(v) for v in box], "axis_xy_A": [float(v) for v in center[:2]],
        "velocity_columns_verified": ["vx", "vy", "vz"],
    }
    return series, meta


def complex_acf(series: np.ndarray, maxlag: int, batch_modes: int = 8) -> np.ndarray:
    """Unbiased, all-origin Re<delta J(t+tau) delta J*(t)> for every mode."""
    nframe = series.shape[0]
    nfft = 1 << (2 * nframe - 1).bit_length()
    output = np.empty((maxlag + 1, series.shape[1], series.shape[2]), dtype=np.float64)
    # Batching avoids a >480 MiB padded FFT workspace for 80 modes while
    # retaining exact, all-origin FFT correlations for every n.
    for first in range(0, series.shape[1], batch_modes):
        last = min(first + batch_modes, series.shape[1])
        centred = series[:, first:last] - series[:, first:last].mean(axis=0, keepdims=True)
        ft = np.fft.fft(centred, n=nfft, axis=0)
        power = np.abs(ft) ** 2
        corr = np.fft.ifft(power, axis=0)[:maxlag + 1].real
        output[:, first:last] = corr
        del centred, ft, power, corr
    return output / np.arange(nframe, nframe - maxlag - 1, -1, dtype=np.float64)[:, None, None]


def symmetric_periodogram(series: np.ndarray, dt_ps: float) -> tuple[np.ndarray, np.ndarray, int]:
    """Whole-record rectangular periodogram; deliberately no Welch/Hann window."""
    nframe = len(series); ff = np.fft.fft(series - series.mean(axis=0, keepdims=True), axis=0)
    acc = np.abs(ff) ** 2 / nframe
    positive = np.arange(1, nframe // 2 + 1)
    # For complex J_k, retain the direction-independent, real-current spectrum.
    psd = 0.5 * (acc[positive] + acc[-positive])
    freq = positive / (nframe * dt_ps)
    return freq, psd, 1


def write_curve_data(outdir: Path, curves: np.ndarray, k: np.ndarray, time: np.ndarray, nframes: list[int]) -> None:
    """Write raw (per-replica) and ensemble CJJ tables in stable long format."""
    raw_path = outdir / "CJJ_all_modes_per_replica.csv"
    ensemble_path = outdir / "CJJ_all_modes_ensemble_mean_sem.csv"
    # curves: [replica, lag, mode, component(z,r,theta)]
    combined = curves
    c0 = combined[:, 0, :, :]
    norm = combined / c0[:, None, :, :]
    with raw_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["replicate", "branch", "n", "k_inv_A", "wavelength_A", "lag_ps", "n_time_origins", "CJJ_raw_A2_fs2", "CJJ_normalized"])
        for rep in range(combined.shape[0]):
            for branch_idx, branch in enumerate(COMPONENTS):
                for mode_idx, n in enumerate(range(1, combined.shape[2] + 1)):
                    for lag_idx, lag in enumerate(time):
                        writer.writerow([rep + 1, branch, n, k[mode_idx], 2 * np.pi / k[mode_idx], lag,
                                         nframes[rep] - lag_idx, combined[rep, lag_idx, mode_idx, branch_idx],
                                         norm[rep, lag_idx, mode_idx, branch_idx]])
    mean = combined.mean(axis=0)
    sem = combined.std(axis=0, ddof=1) / math.sqrt(combined.shape[0])
    norm_mean = norm.mean(axis=0)
    norm_sem = norm.std(axis=0, ddof=1) / math.sqrt(norm.shape[0])
    with ensemble_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["branch", "n", "k_inv_A", "wavelength_A", "lag_ps", "n_time_origins_min_across_replicas", "CJJ_raw_mean_A2_fs2", "CJJ_raw_replica_SEM_A2_fs2", "CJJ_normalized_mean", "CJJ_normalized_replica_SEM", "n_replicas"])
        for branch_idx, branch in enumerate(COMPONENTS):
            for mode_idx, n in enumerate(range(1, mean.shape[1] + 1)):
                for lag_idx, lag in enumerate(time):
                    writer.writerow([branch, n, k[mode_idx], 2 * np.pi / k[mode_idx], lag, min(nframes) - lag_idx,
                                     mean[lag_idx, mode_idx, branch_idx], sem[lag_idx, mode_idx, branch_idx],
                                     norm_mean[lag_idx, mode_idx, branch_idx], norm_sem[lag_idx, mode_idx, branch_idx],
                                     combined.shape[0]])


def write_spectral_data(outdir: Path, spectra: np.ndarray, frequency: np.ndarray, k: np.ndarray) -> None:
    # spectra: [replica, frequency, mode, component(z,r,theta)]
    combined = spectra
    mean = combined.mean(axis=0)
    sem = combined.std(axis=0, ddof=1) / math.sqrt(combined.shape[0])
    path = outdir / "current_spectra_all_modes_ensemble_mean_sem.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["branch", "n", "k_inv_A", "frequency_ps_inv", "PSD_mean_arbitrary", "PSD_replica_SEM_arbitrary", "n_replicas"])
        for branch_idx, branch in enumerate(COMPONENTS):
            for mode_idx, n in enumerate(range(1, mean.shape[1] + 1)):
                for fi, f in enumerate(frequency):
                    writer.writerow([branch, n, k[mode_idx], f, mean[fi, mode_idx, branch_idx], sem[fi, mode_idx, branch_idx], combined.shape[0]])


def extract_peaks(spectra: np.ndarray, frequency: np.ndarray, k: np.ndarray, fmax: float) -> tuple[list[dict], list[dict]]:
    """Return per-replica operational peaks and ensemble dispersion summaries."""
    combined = spectra
    fmin = max(2 * (frequency[1] - frequency[0]), 0.01)
    band = (frequency >= fmin) & (frequency <= fmax)
    if not np.any(band):
        raise ValueError("empty frequency band for peak extraction")
    raw_rows: list[dict] = []
    summary_rows: list[dict] = []
    for branch_idx, branch in enumerate(COMPONENTS):
        for mi, n in enumerate(range(1, len(k) + 1)):
            f_values, prom_ratios = [], []
            for rep in range(combined.shape[0]):
                y = combined[rep, band, mi, branch_idx]
                fb = frequency[band]
                peaks, props = find_peaks(y, prominence=np.median(y) * 0.10)
                if len(peaks) == 0:
                    selected = int(np.argmax(y))
                    prominence = 0.0
                else:
                    best = int(np.argmax(y[peaks]))
                    selected = int(peaks[best])
                    prominence = float(props["prominences"][best])
                ratio = float(prominence / max(np.median(y), np.finfo(float).tiny))
                f_peak = float(fb[selected])
                raw_rows.append({"replicate": rep + 1, "branch": branch, "n": n, "k_inv_A": float(k[mi]),
                                 "frequency_peak_ps_inv": f_peak, "omega_peak_rad_ps": 2 * np.pi * f_peak,
                                 "phase_velocity_A_ps": 2 * np.pi * f_peak / k[mi], "prominence_over_median": ratio})
                f_values.append(f_peak)
                prom_ratios.append(ratio)
            f_values = np.asarray(f_values)
            mean_f = float(f_values.mean())
            sem_f = float(f_values.std(ddof=1) / math.sqrt(len(f_values)))
            mean_prom = float(np.mean(prom_ratios))
            # A transparent reliability gate: peak must exceed the local spectral
            # median and the seed-to-seed frequency scatter cannot exceed 25%.
            cv = float(f_values.std(ddof=1) / mean_f) if mean_f > 0 else np.inf
            resolved = bool(mean_prom >= 2.0 and cv <= 0.25)
            summary_rows.append({"branch": branch, "n": n, "k_inv_A": float(k[mi]), "wavelength_A": float(2 * np.pi / k[mi]),
                                 "frequency_peak_mean_ps_inv": mean_f, "frequency_peak_replica_SEM_ps_inv": sem_f,
                                 "omega_peak_mean_rad_ps": float(2 * np.pi * mean_f), "omega_peak_replica_SEM_rad_ps": float(2 * np.pi * sem_f),
                                 "phase_velocity_A_ps": float(2 * np.pi * mean_f / k[mi]),
                                 "mean_prominence_over_median": mean_prom, "replica_frequency_CV": cv,
                                 "resolved_operational_peak": resolved, "n_replicas": len(f_values)})
    return raw_rows, summary_rows


def write_peak_data(outdir: Path, raw_rows: list[dict], summary_rows: list[dict]) -> None:
    for filename, rows in (("dispersion_peaks_per_replica.csv", raw_rows), ("LA_TA_dispersion.csv", summary_rows)):
        with (outdir / filename).open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def write_modal_weight_data(outdir: Path, curves: np.ndarray, k: np.ndarray, n_oxygen: int) -> tuple[list[dict], list[dict]]:
    """Current-spectrum (not global VDOS) weights and guarded low-k diagnostics."""
    weights = curves[:, 0, :, :]  # C_JJ(k,0), one value per seed/mode/component
    mean = weights.mean(axis=0)
    sem = weights.std(axis=0, ddof=1) / math.sqrt(weights.shape[0])
    total = weights.sum(axis=1, keepdims=True)
    fraction = weights / total
    frac_mean = fraction.mean(axis=0)
    frac_sem = fraction.std(axis=0, ddof=1) / math.sqrt(fraction.shape[0])
    rows: list[dict] = []
    for bi, branch in enumerate(COMPONENTS):
        for mi, n in enumerate(range(1, len(k) + 1)):
            rows.append({"branch": branch, "n": n, "k_inv_A": float(k[mi]), "wavelength_A": float(2 * np.pi / k[mi]),
                         "CJJ0_mean_A2_fs2": float(mean[mi, bi]), "CJJ0_replica_SEM_A2_fs2": float(sem[mi, bi]),
                         "CJJ0_per_oxygen_mean_A2_fs2": float(mean[mi, bi] / n_oxygen),
                         "current_spectrum_weight_fraction_mean": float(frac_mean[mi, bi]),
                         "current_spectrum_weight_fraction_replica_SEM": float(frac_sem[mi, bi]),
                         "n_replicas": weights.shape[0]})
    with (outdir / "current_mode_weight_DOS_proxy_all_n.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    # Fits are deliberately reported as finite-k diagnostics.  They are not a
    # thermodynamic-limit DOS extrapolation: one fixed 5L box only supplies four
    # lowest nonzero k values.
    fit_rows: list[dict] = []
    mfit = min(4, len(k))
    x = k[:mfit]
    for bi, branch in enumerate(COMPONENTS):
        intercepts, slopes, exponents = [], [], []
        for rep in range(weights.shape[0]):
            y = weights[rep, :mfit, bi]
            slope, intercept = np.polyfit(x, y, 1)
            intercepts.append(float(intercept)); slopes.append(float(slope))
            if np.all(y > 0):
                exponents.append(float(np.polyfit(np.log(x), np.log(y), 1)[0]))
        def avg_sem(values: list[float]) -> tuple[float, float]:
            arr = np.asarray(values, dtype=float)
            return float(arr.mean()), float(arr.std(ddof=1) / math.sqrt(len(arr)))
        intercept_mean, intercept_sem = avg_sem(intercepts)
        slope_mean, slope_sem = avg_sem(slopes)
        alpha_mean, alpha_sem = avg_sem(exponents)
        fit_rows.append({"branch": branch, "fit_modes": f"n=1..{mfit}", "linear_CJJ0_intercept_at_k0_mean_A2_fs2": intercept_mean,
                         "linear_CJJ0_intercept_at_k0_replica_SEM_A2_fs2": intercept_sem,
                         "linear_slope_mean_A2_fs2_A": slope_mean, "linear_slope_replica_SEM_A2_fs2_A": slope_sem,
                         "loglog_power_exponent_mean": alpha_mean, "loglog_power_exponent_replica_SEM": alpha_sem,
                         "interpretation_boundary": "finite-k diagnostic only; a positive log-log exponent is suggestive, not proof, of vanishing k-to-0 current spectral weight"})
    with (outdir / "low_k_current_spectrum_weight_diagnostic.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fit_rows[0])); writer.writeheader(); writer.writerows(fit_rows)
    return rows, fit_rows


def savefig(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def plot_cjj_pages(figdir: Path, curves: np.ndarray, time: np.ndarray, k: np.ndarray, plot_max_ps: float, system_label: str) -> None:
    combined = curves
    norm = combined / combined[:, 0:1]
    mean = norm.mean(axis=0)
    sem = norm.std(axis=0, ddof=1) / math.sqrt(norm.shape[0])
    shown = time <= plot_max_ps
    for bi, branch in enumerate(COMPONENTS):
        color = {"LA": LA_COLOR, "TA_r": TAR_COLOR, "TA_theta": TATHETA_COLOR}[branch]
        for page, start in enumerate(range(0, len(k), 16), start=1):
            fig = plt.figure(figsize=(7.0, 6.9))
            for local, mi in enumerate(range(start, min(start + 16, len(k)))):
                row, col = divmod(local, 4)
                left, bottom = 0.095 + col * 0.225, 0.08 + (3 - row) * 0.225
                ax = fig.add_axes([left, bottom, 0.18, 0.17])
                for rep in range(norm.shape[0]):
                    ax.plot(time[shown], norm[rep, shown, mi, bi], color=color, alpha=0.20, lw=1.0)
                ax.fill_between(time[shown], mean[shown, mi, bi] - sem[shown, mi, bi], mean[shown, mi, bi] + sem[shown, mi, bi], color=color, alpha=0.20, linewidth=0)
                ax.plot(time[shown], mean[shown, mi, bi], color=color, lw=1.15)
                ax.axhline(0, color="#767676", lw=0.8, zorder=0)
                ax.set_xlim(0, plot_max_ps)
                ax.set_ylim(-0.55, 1.08)
                ax.set_title(rf"$n={mi+1}$, $k={k[mi]:.3f}$ Å$^{{-1}}$", fontsize=6, pad=2)
                if row == 3:
                    ax.set_xlabel(r"$t$ (ps)")
                else:
                    ax.set_xticklabels([])
                if col == 0:
                    ax.set_ylabel(r"$C_{JJ}/C_{JJ}(0)$")
                else:
                    ax.set_yticklabels([])
                ax.tick_params(labelsize=6, length=2.5)
            fig.text(0.5, 0.985, f"{system_label} {branch} current ACF; mean ± replica SEM, thin = individual seeds", ha="center", va="top", fontsize=8)
            savefig(fig, figdir / f"CJJ_{branch}_modes_page_{page:02d}")


def plot_heatmaps(figdir: Path, curves: np.ndarray, time: np.ndarray, k: np.ndarray, plot_max_ps: float) -> None:
    combined = curves
    norm_mean = (combined / combined[:, 0:1]).mean(axis=0)
    shown = time <= plot_max_ps
    for bi, branch in enumerate(COMPONENTS):
        fig = plt.figure(figsize=(3.3, 2.8))
        ax = fig.add_axes([0.18, 0.18, 0.64, 0.70])
        cax = fig.add_axes([0.85, 0.18, 0.035, 0.70])
        im = ax.imshow(norm_mean[shown, :, bi].T, aspect="auto", origin="lower", interpolation="nearest",
                       extent=[time[shown][0], time[shown][-1], 0.5, len(k) + 0.5], cmap="RdBu_r", vmin=-0.5, vmax=1.0)
        ax.set(xlabel=r"$t$ (ps)", ylabel=r"Mode $n$")
        ax.tick_params(length=3)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(r"$C_{JJ}/C_{JJ}(0)$")
        savefig(fig, figdir / f"CJJ_{branch}_all_modes_heatmap")


def plot_dispersion(figdir: Path, summary_rows: list[dict]) -> None:
    fig = plt.figure(figsize=(5.5, 2.55))
    ax1 = fig.add_axes([0.11, 0.20, 0.36, 0.70])
    ax2 = fig.add_axes([0.58, 0.20, 0.36, 0.70])
    for branch, color in (("LA", LA_COLOR), ("TA_r", TAR_COLOR), ("TA_theta", TATHETA_COLOR)):
        rows = [row for row in summary_rows if row["branch"] == branch]
        resolved = np.array([row["resolved_operational_peak"] for row in rows], dtype=bool)
        k = np.array([row["k_inv_A"] for row in rows])
        omega = np.array([row["omega_peak_mean_rad_ps"] for row in rows])
        omega_sem = np.array([row["omega_peak_replica_SEM_rad_ps"] for row in rows])
        vphase = np.array([row["phase_velocity_A_ps"] for row in rows])
        ax1.errorbar(k[resolved], omega[resolved], yerr=omega_sem[resolved], fmt="o", color=color, ms=3.5, lw=1, capsize=2, label=f"{branch}, resolved")
        ax1.scatter(k[~resolved], omega[~resolved], facecolors="none", edgecolors=color, s=22, linewidths=1, label=f"{branch}, tentative")
        ax2.plot(k[resolved], vphase[resolved], "o", color=color, ms=3.5, label=f"{branch}, resolved")
        ax2.scatter(k[~resolved], vphase[~resolved], facecolors="none", edgecolors=color, s=22, linewidths=1, label=f"{branch}, tentative")
    ax1.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel=r"$\omega_{\mathrm{peak}}$ (rad ps$^{-1}$)")
    ax2.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel=r"$\omega_{\mathrm{peak}}/k_n$ (Å ps$^{-1}$)")
    ax1.legend(fontsize=5.8, loc="upper left")
    ax2.legend(fontsize=5.8, loc="upper right")
    for ax, label in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(-0.18, 1.02, label, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom")
        ax.tick_params(length=3)
    savefig(fig, figdir / "LA_TA_dispersion_operational_peaks")


def plot_modal_weights(figdir: Path, weight_rows: list[dict], fit_rows: list[dict], n_oxygen: int, nmax: int) -> None:
    color_map = {"LA": LA_COLOR, "TA_r": TAR_COLOR, "TA_theta": TATHETA_COLOR}
    fig = plt.figure(figsize=(5.5, 2.55))
    ax1 = fig.add_axes([0.11, 0.20, 0.36, 0.70])
    ax2 = fig.add_axes([0.58, 0.20, 0.36, 0.70])
    for branch in COMPONENTS:
        r = [row for row in weight_rows if row["branch"] == branch]
        n = np.array([row["n"] for row in r])
        w = np.array([row["CJJ0_per_oxygen_mean_A2_fs2"] for row in r])
        wsem = np.array([row["CJJ0_replica_SEM_A2_fs2"] / n_oxygen for row in r])
        frac = np.array([row["current_spectrum_weight_fraction_mean"] for row in r])
        color = color_map[branch]
        ax1.errorbar(n, w, yerr=wsem, fmt="o-", ms=2.5, lw=1.0, capsize=1.5, color=color, label=branch)
        ax2.plot(n, frac, "o-", ms=2.5, lw=1.0, color=color, label=branch)
    ax1.set_yscale("log"); ax2.set_yscale("log")
    ax1.set(xlabel=r"Mode $n$", ylabel=r"$C_{JJ}(k_n,0)/N$ (Å$^2$ fs$^{-2}$)")
    ax2.set(xlabel=r"Mode $n$", ylabel=rf"Fraction of $\sum_{{n=1}}^{{{nmax}}} C_{{JJ}}(k_n,0)$")
    ax1.legend(fontsize=6, ncol=1); ax2.legend(fontsize=6, ncol=1)
    for ax, label in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(-0.18, 1.02, label, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom")
        ax.tick_params(length=3)
    savefig(fig, figdir / "current_mode_weight_DOS_proxy_all_n")

    fig = plt.figure(figsize=(5.5, 2.1))
    for bi, branch in enumerate(COMPONENTS):
        ax = fig.add_axes([0.10 + bi * 0.30, 0.23, 0.22, 0.65])
        r = [row for row in weight_rows if row["branch"] == branch and row["n"] <= 4]
        kval = np.array([row["k_inv_A"] for row in r])
        y = np.array([row["CJJ0_per_oxygen_mean_A2_fs2"] for row in r])
        ysem = np.array([row["CJJ0_replica_SEM_A2_fs2"] / n_oxygen for row in r])
        color = color_map[branch]
        ax.errorbar(kval, y, yerr=ysem, fmt="o", color=color, ms=3.5, capsize=2)
        fit = next(row for row in fit_rows if row["branch"] == branch)
        xline = np.linspace(0, kval[-1] * 1.05, 100)
        ax.plot(xline, (fit["linear_CJJ0_intercept_at_k0_mean_A2_fs2"] + fit["linear_slope_mean_A2_fs2_A"] * xline) / n_oxygen, color=color, lw=1.0)
        ax.axvline(0, color="#767676", lw=0.8)
        ax.set_title(branch, fontsize=7, pad=2)
        ax.set_xlabel(r"$k_n$ (Å$^{-1}$)")
        if bi == 0:
            ax.set_ylabel(r"$C_{JJ}(0)/N$ (Å$^2$ fs$^{-2}$)")
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=6, length=3)
    fig.text(0.5, 0.98, "n = 1–4 linear finite-k diagnostic; do not read as a thermodynamic-limit DOS extrapolation", ha="center", va="top", fontsize=7)
    savefig(fig, figdir / "low_k_current_spectrum_weight_diagnostic")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--dump-template", default="NVT20ns_5xL_8_8_RH75_N665_rep{rep}_20260719/nvt20ns_8_8_RH75_5L_rep{rep}_oxygen_10fs_1ns.dump", help="glob below input root; must contain {rep} for replicas 1..3")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--system-label", default="(8,8) CNT water, 5L")
    parser.add_argument("--nmax", type=int, default=80, help="integer axial modes, default 1..80 (k <= about 1 Å^-1)")
    parser.add_argument("--max-lag-ps", type=float, default=100.0)
    parser.add_argument("--plot-max-ps", type=float, default=50.0)
    parser.add_argument("--welch-nperseg", type=int, default=16384)
    parser.add_argument("--peak-fmax-ps-inv", type=float, default=8.0)
    parser.add_argument("--expected-frames", type=int, default=100_001)
    parser.add_argument("--oxygen-type", type=int, default=3)
    args = parser.parse_args()
    if args.nmax < 1 or args.max_lag_ps <= 0:
        raise ValueError("nmax and max-lag-ps must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    data_dir = args.outdir / "derived_data"
    fig_dir = args.outdir / "figures"
    data_dir.mkdir(exist_ok=True)
    fig_dir.mkdir(exist_ok=True)
    if "{rep}" not in args.dump_template:
        raise ValueError("--dump-template must contain {rep}")
    dumps = locate_dumps(args.input_root, args.dump_template)

    curves_list, spectra, metas = [], [], []
    dt_ps = None
    maxlag = None
    frequency = None
    nsegments = None
    reference_box = None
    for dump in dumps:
        print(f"Reading {dump}", flush=True)
        series, meta = read_modal_currents(dump, args.nmax, args.expected_frames, args.oxygen_type)
        metas.append(meta)
        print(f"  frames={meta['n_frames']}, dt={meta['dt_ps']} ps, components={meta['velocity_columns_verified']}", flush=True)
        if dt_ps is None:
            dt_ps = meta["dt_ps"]
            maxlag = int(round(args.max_lag_ps / dt_ps))
            reference_box = tuple(meta["box_A"])
        elif meta["dt_ps"] != dt_ps or tuple(meta["box_A"]) != reference_box:
            raise ValueError("replicas do not share cadence and fixed box")
        if maxlag >= len(series):
            raise ValueError("max lag must be shorter than every trajectory")
        curves_list.append(complex_acf(series, maxlag))
        f, spectrum, ns = symmetric_periodogram(series, dt_ps)
        if frequency is None:
            frequency, nsegments = f, ns
        elif not np.allclose(f, frequency) or ns != nsegments:
            raise ValueError("inconsistent periodogram grids across replicas")
        spectra.append(spectrum)
        del series
    nframes = [meta["n_frames"] for meta in metas]
    min_nframe = min(nframes)
    lz = metas[0]["box_A"][2]
    k = 2 * np.pi * np.arange(1, args.nmax + 1) / lz
    if maxlag >= min_nframe:
        raise ValueError("max lag must be shorter than trajectory")
    time = np.arange(maxlag + 1) * dt_ps
    curves = np.asarray(curves_list)
    spectra_a = np.asarray(spectra)
    write_curve_data(data_dir, curves, k, time, nframes)
    write_spectral_data(data_dir, spectra_a, frequency, k)
    raw_peaks, summary_peaks = extract_peaks(spectra_a, frequency, k, args.peak_fmax_ps_inv)
    write_peak_data(data_dir, raw_peaks, summary_peaks)
    n_oxygen = metas[0]["n_water"]
    if any(meta["n_water"] != n_oxygen for meta in metas):
        raise ValueError("replicas do not share the same oxygen count")
    weight_rows, fit_rows = write_modal_weight_data(data_dir, curves, k, n_oxygen)
    plot_cjj_pages(fig_dir, curves, time, k, args.plot_max_ps, args.system_label)
    plot_heatmaps(fig_dir, curves, time, k, args.plot_max_ps)
    plot_dispersion(fig_dir, summary_peaks)
    plot_modal_weights(fig_dir, weight_rows, fit_rows, n_oxygen, args.nmax)

    contract = """Core conclusion: determine whether the axial longitudinal and transverse current spectra contain reproducible peaks from k_min through k≈1 Å^-1.
Figure archetype: quantitative grid.
Backend: Python/matplotlib.
Evidence hierarchy: all-origin per-seed CJJ curves; replica SEM; independent Welch peak extraction.
Reviewer risks: high-k peak assignment can be non-acoustic or unresolved, so open symbols mark points failing the declared prominence/replica-scatter gate.
"""
    (args.outdir / "figure_contract.txt").write_text(contract, encoding="utf-8")
    metadata = {
        "system": args.system_label, "input_replicas": metas,
        "mode_range": {"n_min": 1, "n_max": args.nmax, "k_min_inv_A": float(k[0]), "k_max_inv_A": float(k[-1])},
        "current_definition": "J_a(k_n,t)=sum_i v'_i,a(t) exp(i k_n z_i(t)); v' is instantaneous O-COM-subtracted Cartesian velocity; a=z,r,theta; e_r=(x-xaxis,y-yaxis)/r and e_theta=(-e_r_y,e_r_x); radial and circumferential TA are separate",
        "CJJ_definition": "Re<delta J_a(t+tau) delta J_a*(t)> over every time origin; delta removes each trajectory's temporal modal-current mean",
        "correlation_max_lag_ps": float(time[-1]), "plot_window_ps": float(args.plot_max_ps),
        "dispersion_estimator": "maximum locally prominent nonzero peak of the symmetrized complex-current whole-record rectangular periodogram; no Welch segmentation and no Hann window",
        "spectral_estimator": "rectangular_whole_record_no_welch_no_hann", "frequency_resolution_ps_inv": float(frequency[1] - frequency[0]), "periodogram_records_per_replica": nsegments,
        "peak_quality_gate": "resolved only when mean peak prominence/median PSD >=2 and replica frequency CV <=0.25; other points remain tentative",
        "units": {"velocity": "Angstrom/fs (LAMMPS real units)", "CJJ_raw": "Angstrom^2/fs^2", "frequency": "ps^-1", "omega": "rad/ps"},
        "DOS_proxy_definition": f"At each discrete k_n, the integral of the modal current PSD is proportional to C_JJ(k_n,0). Reported DOS-proxy fractions normalize this modal current spectral weight over n=1..{args.nmax} within the same branch; they are not the global single-particle velocity DOS.",
        "outputs": {"raw_per_replica_CJJ": "derived_data/CJJ_all_modes_per_replica.csv", "ensemble_CJJ": "derived_data/CJJ_all_modes_ensemble_mean_sem.csv", "dispersion": "derived_data/LA_TA_dispersion.csv", "current_spectrum_weight": "derived_data/current_mode_weight_DOS_proxy_all_n.csv"},
    }
    (args.outdir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (args.outdir / "FINISHED.txt").write_text(f"{args.system_label}: all-mode current correlations and whole-record rectangular-periodogram analysis finished successfully.\n", encoding="utf-8")
    print(json.dumps({"outdir": str(args.outdir), "n_modes": args.nmax, "n_frames_per_replica": nframes, "dt_ps": dt_ps, "k_range_inv_A": [float(k[0]), float(k[-1])], "periodogram_records": nsegments}, indent=2), flush=True)


if __name__ == "__main__":
    main()
