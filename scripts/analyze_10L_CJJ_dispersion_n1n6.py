"""Axial current correlations, spectra, and low-k dispersion for 10L CNT water."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 1.0, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0, "pdf.fonttype": 42, "svg.fonttype": "none",
})


def read_current(path: Path, nmax: int) -> tuple[np.ndarray, float, float]:
    currents: list[np.ndarray] = []
    lengths: list[float] = []
    times: list[float] = []
    with path.open() as handle:
        while True:
            marker = handle.readline()
            if not marker:
                break
            if marker.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Unexpected dump marker in {path}: {marker!r}")
            step = int(handle.readline())
            handle.readline()
            natom = int(handle.readline())
            handle.readline()
            bounds = [handle.readline().split() for _ in range(3)]
            header = handle.readline().split()[2:]
            if header != ["id", "z", "vz"]:
                raise ValueError(f"Unexpected columns in {path}: {header}")
            frame = np.fromstring(" ".join(handle.readline() for _ in range(natom)), sep=" ").reshape(natom, 3)
            z, vz = frame[:, 1], frame[:, 2]
            length = float(bounds[2][1]) - float(bounds[2][0])
            k = 2 * np.pi * np.arange(1, nmax + 1) / length
            currents.append(np.exp(1j * np.outer(k, z)) @ vz)
            lengths.append(length)
            times.append(step * 0.0005)  # LAMMPS real units; dt = 0.5 fs
    time = np.asarray(times)
    dt = float(np.median(np.diff(time)))
    return np.asarray(currents).T, float(np.mean(lengths)), dt


def autocorrelation(current: np.ndarray, max_lag: int) -> np.ndarray:
    """Unbiased C(tau)=<J(t+tau)J*(t)> for each k, normalized at tau=0."""
    count = current.shape[-1]
    centered = current - current.mean(axis=-1, keepdims=True)
    fft_length = 1 << int(np.ceil(np.log2(2 * count - 1)))
    transform = np.fft.fft(centered, n=fft_length, axis=-1)
    corr = np.fft.ifft(transform * transform.conj(), axis=-1)[..., :max_lag + 1]
    corr /= (count - np.arange(max_lag + 1))[None, :]
    corr /= corr[..., :1].real
    return corr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dumps", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--nmax", type=int, default=6)
    parser.add_argument("--correlation-max-ps", type=float, default=500.0)
    parser.add_argument("--display-max-ps", type=float, default=100.0)
    parser.add_argument("--spectrum-max-omega", type=float, default=1.25)
    args = parser.parse_args()
    data_dir, figure_dir = args.output / "derived_data", args.output / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(exist_ok=True)

    all_current: list[np.ndarray] = []
    lengths: list[float] = []
    dts: list[float] = []
    for dump in args.dumps:
        current, length, dt = read_current(dump, args.nmax)
        all_current.append(current)
        lengths.append(length)
        dts.append(dt)
    if not np.allclose(dts, dts[0]):
        raise ValueError(f"inconsistent sampling intervals: {dts}")
    dt = dts[0]
    k_values = 2 * np.pi * np.arange(1, args.nmax + 1) / float(np.mean(lengths))
    max_lag = min(round(args.correlation_max_ps / dt), all_current[0].shape[-1] - 1)
    corr_replicas = np.asarray([autocorrelation(x, max_lag) for x in all_current])
    corr_mean = corr_replicas.mean(axis=0)
    corr_sem = corr_replicas.real.std(axis=0, ddof=1) / np.sqrt(len(corr_replicas))
    lag_ps = np.arange(max_lag + 1) * dt

    # Full-time-series periodogram provides the finest resolvable low-frequency acoustic peaks.
    power: list[np.ndarray] = []
    for current in all_current:
        centered = current - current.mean(axis=-1, keepdims=True)
        window = np.hanning(centered.shape[-1])
        power.append(np.abs(np.fft.fft(centered * window, axis=-1)) ** 2 / np.sum(window**2))
    power_array = np.asarray(power)
    omega = np.fft.fftfreq(power_array.shape[-1], dt) * 2 * np.pi
    order = np.argsort(omega)
    omega = omega[order]
    power_array = power_array[..., order]
    spectrum = power_array.mean(axis=0)
    spectrum_sem = power_array.std(axis=0, ddof=1) / np.sqrt(len(power_array))

    cjj_rows, skw_rows = [], []
    for n, k in enumerate(k_values, start=1):
        cjj_rows.extend({"n": n, "k_inv_A": k, "lag_ps": t, "C_real_normalized": c.real,
                         "C_imag_normalized": c.imag, "C_real_replica_SEM": e}
                        for t, c, e in zip(lag_ps, corr_mean[n-1], corr_sem[n-1]))
        skw_rows.extend({"n": n, "k_inv_A": k, "omega_rad_ps": w, "S_JzJz": s,
                         "S_replica_SEM": e}
                        for w, s, e in zip(omega, spectrum[n-1], spectrum_sem[n-1]))
    pd.DataFrame(cjj_rows).to_csv(data_dir / "C_JzJz_10L_n001_n006.csv", index=False)
    pd.DataFrame(skw_rows).to_csv(data_dir / "C_k_omega_JzJz_10L_n001_n006.csv", index=False)

    # Dispersion from positive-frequency maxima in a broad sound-branch window.
    peak_rows = []
    for n, k in enumerate(k_values, start=1):
        expected = 16.0 * k
        select = (omega > max(0.01, 0.45 * expected)) & (omega < 1.65 * expected)
        local_w, local_s = omega[select], spectrum[n-1, select]
        index = int(np.argmax(local_s))
        peak_rows.append({"n": n, "k_inv_A": k, "omega_peak_rad_ps": local_w[index],
                          "peak_search_low_rad_ps": local_w[0], "peak_search_high_rad_ps": local_w[-1],
                          "frequency_resolution_rad_ps": 2 * np.pi / (all_current[0].shape[-1] * dt)})
    peaks = pd.DataFrame(peak_rows)
    slope, intercept = np.polyfit(peaks.k_inv_A, peaks.omega_peak_rad_ps, 1)
    pred = slope * peaks.k_inv_A + intercept
    r2 = 1 - np.sum((peaks.omega_peak_rad_ps - pred) ** 2) / np.sum((peaks.omega_peak_rad_ps - peaks.omega_peak_rad_ps.mean()) ** 2)
    through_origin = float(np.sum(peaks.k_inv_A * peaks.omega_peak_rad_ps) / np.sum(peaks.k_inv_A**2))
    peaks["omega_linear_fit_rad_ps"] = pred
    peaks.to_csv(data_dir / "LA_lowk_dispersion_n001_n006.csv", index=False)
    pd.DataFrame([{"fit": "omega=c*k+b", "sound_speed_A_per_ps": slope, "intercept_rad_ps": intercept, "R2": r2,
                   "through_origin_sound_speed_A_per_ps": through_origin, "n_min": 1, "n_max": args.nmax}]).to_csv(
                       data_dir / "LA_lowk_sound_speed_fit.csv", index=False)

    # Explicit 2 by 3 grids, with k as the only panel identifier.
    display = lag_ps <= args.display_max_ps
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.3), sharex=True, sharey=True)
    for i, (ax, k) in enumerate(zip(axes.flat, k_values)):
        ax.plot(lag_ps[display], corr_mean[i, display].real, color="#2166ac", lw=1.1)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.text(0.04, 0.90, rf"$k={k:.4f}\ \AA^{{-1}}$", transform=ax.transAxes)
        ax.tick_params(top=True, right=True)
        ax.text(-0.20, 1.06, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.supxlabel(r"Lag time, $t$ (ps)")
    fig.supylabel(r"$\mathrm{Re}\,[C_{J_zJ_z}(k,t)/C_{J_zJ_z}(k,0)]$")
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.13, top=0.96, wspace=0.18, hspace=0.20)
    for ext, dpi in [("png", 600), ("pdf", None), ("svg", None), ("tiff", 600)]:
        fig.savefig(figure_dir / f"C_JzJz_10L_n001_n006_2x3.{ext}", **({"dpi": dpi} if dpi else {}))
    # Explicit exports retained for static publication-figure preflight.
    fig.savefig(figure_dir / "C_JzJz_10L_n001_n006_2x3.svg")
    fig.savefig(figure_dir / "C_JzJz_10L_n001_n006_2x3.pdf")
    fig.savefig(figure_dir / "C_JzJz_10L_n001_n006_2x3.tiff", dpi=600)
    plt.close(fig)

    spectral_window = np.abs(omega) <= args.spectrum_max_omega
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.3), sharex=True)
    for i, (ax, k) in enumerate(zip(axes.flat, k_values)):
        ax.semilogy(omega[spectral_window], spectrum[i, spectral_window], color="#b2182b", lw=1.05)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.text(0.04, 0.90, rf"$k={k:.4f}\ \AA^{{-1}}$", transform=ax.transAxes)
        ax.tick_params(top=True, right=True)
        ax.text(-0.20, 1.06, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.supxlabel(r"Angular frequency, $\omega$ (rad ps$^{-1}$)")
    fig.supylabel(r"$C(k,\omega)$ (arb. units)")
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.13, top=0.96, wspace=0.18, hspace=0.20)
    for ext, dpi in [("png", 600), ("pdf", None), ("svg", None), ("tiff", 600)]:
        fig.savefig(figure_dir / f"C_k_omega_JzJz_10L_n001_n006_2x3.{ext}", **({"dpi": dpi} if dpi else {}))
    fig.savefig(figure_dir / "C_k_omega_JzJz_10L_n001_n006_2x3.svg")
    fig.savefig(figure_dir / "C_k_omega_JzJz_10L_n001_n006_2x3.pdf")
    fig.savefig(figure_dir / "C_k_omega_JzJz_10L_n001_n006_2x3.tiff", dpi=600)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.plot(peaks.k_inv_A, peaks.omega_peak_rad_ps, "o", ms=4, color="#2166ac", label="spectral peak")
    xx = np.linspace(0, peaks.k_inv_A.max() * 1.08, 200)
    ax.plot(xx, slope * xx + intercept, color="#b2182b", lw=1.1,
            label=rf"$\omega=({slope:.2f}\ \AA\,\mathrm{{ps}}^{{-1}})k{intercept:+.3f}$")
    ax.set(xlabel=r"$k$ ($\AA^{-1}$)", ylabel=r"$\omega_\mathrm{LA}$ (rad ps$^{-1}$)")
    ax.legend(loc="upper left", fontsize=6)
    ax.tick_params(top=True, right=True)
    ax.text(-0.18, 1.05, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.18, top=0.96)
    for ext, dpi in [("png", 600), ("pdf", None), ("svg", None), ("tiff", 600)]:
        fig.savefig(figure_dir / f"LA_lowk_dispersion_n001_n006.{ext}", **({"dpi": dpi} if dpi else {}))
    fig.savefig(figure_dir / "LA_lowk_dispersion_n001_n006.svg")
    fig.savefig(figure_dir / "LA_lowk_dispersion_n001_n006.pdf")
    fig.savefig(figure_dir / "LA_lowk_dispersion_n001_n006.tiff", dpi=600)
    plt.close(fig)

    (args.output / "metadata.json").write_text(json.dumps({
        "source": "(8,8) 10L, Tdamp=1000 ps, 10 ns, 100 fs, eight velocity replicas",
        "current_definition": "J_z(k,t)=sum_O v_z exp(i k z)",
        "correlation": "unbiased complex current autocorrelation; plots show normalized real part",
        "correlation_lag_ps": float(lag_ps[-1]), "spectrum": "Hann-windowed whole-trajectory periodogram",
        "dispersion": "positive-frequency acoustic-peak maxima, n=1..6",
    }, indent=2))
    (args.output / "FINISHED.txt").write_text("10L axial CJJ and low-k dispersion analysis finished successfully.\n")


if __name__ == "__main__":
    main()
