"""Render requested 10L axial CJJ/S(k,omega) panels from archived local tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 1.0, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0, "pdf.fonttype": 42, "svg.fonttype": "none",
})


def export(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".tiff"), dpi=600)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source_data = args.source / "derived_data"
    output_data, output_fig = args.output / "derived_data", args.output / "figures"
    output_data.mkdir(parents=True, exist_ok=True)
    output_fig.mkdir(exist_ok=True)
    cjj = pd.read_csv(source_data / "CJJ_all_modes_ensemble_mean_sem.csv")
    psd = pd.read_csv(source_data / "current_spectra_all_modes_ensemble_mean_sem.csv")
    cjj = cjj[(cjj.branch == "LA") & (cjj.n.between(1, 6))].copy()
    psd = psd[(psd.branch == "LA") & (psd.n.between(1, 6))].copy()

    # The raw archive gives positive-frequency Welch spectra.  Reflect them using equilibrium symmetry.
    negative = psd.copy()
    negative["omega_rad_ps"] *= -1
    spectral_signed = pd.concat([negative, psd], ignore_index=True).sort_values(["n", "omega_rad_ps"])
    spectral_signed.to_csv(output_data / "C_k_omega_JzJz_10L_n001_n006_signed_from_local_archive.csv", index=False)
    normalized_rows = []
    for n, frame in cjj.groupby("n"):
        c0 = frame.CJJ_mean_A2_fs2.iloc[0]
        normalized_rows.append(frame.assign(CJJ_normalized=frame.CJJ_mean_A2_fs2 / c0,
                                             CJJ_normalized_SEM=frame.CJJ_replica_SEM_A2_fs2 / c0))
    cjj_norm = pd.concat(normalized_rows, ignore_index=True)
    cjj_norm.to_csv(output_data / "C_JzJz_10L_n001_n006_normalized_from_local_archive.csv", index=False)

    # Sound-branch maxima: search near c=16 A/ps, then fit omega = c*k+b.
    peaks = []
    for n, frame in psd.groupby("n"):
        k = frame.k_inv_A.iloc[0]
        expected = 16.0 * k
        window = frame[(frame.omega_rad_ps >= max(0.01, 0.45 * expected)) &
                       (frame.omega_rad_ps <= 1.65 * expected)]
        row = window.loc[window.PSD_mean_arbitrary.idxmax()]
        peaks.append({"n": n, "k_inv_A": k, "omega_peak_rad_ps": row.omega_rad_ps,
                      "frequency_resolution_rad_ps": frame.omega_rad_ps.iloc[0]})
    peaks = pd.DataFrame(peaks)
    speed, intercept = np.polyfit(peaks.k_inv_A, peaks.omega_peak_rad_ps, 1)
    fit = speed * peaks.k_inv_A + intercept
    r2 = 1 - np.sum((peaks.omega_peak_rad_ps - fit) ** 2) / np.sum((peaks.omega_peak_rad_ps - peaks.omega_peak_rad_ps.mean()) ** 2)
    through_zero = float(np.sum(peaks.k_inv_A * peaks.omega_peak_rad_ps) / np.sum(peaks.k_inv_A ** 2))
    peaks["omega_linear_fit_rad_ps"] = fit
    peaks.to_csv(output_data / "LA_lowk_dispersion_n001_n006_from_local_archive.csv", index=False)
    pd.DataFrame([{"model": "omega=c*k+b", "sound_speed_A_per_ps": speed,
                   "intercept_rad_ps": intercept, "R2": r2,
                   "through_origin_sound_speed_A_per_ps": through_zero,
                   "n_min": 1, "n_max": 6}]).to_csv(output_data / "LA_lowk_sound_speed_fit_from_local_archive.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.3), sharex=True, sharey=True)
    for i, n in enumerate(range(1, 7)):
        ax = axes.flat[i]
        frame = cjj_norm[(cjj_norm.n == n) & (cjj_norm.time_ps <= 100)]
        ax.plot(frame.time_ps, frame.CJJ_normalized, color="#2166ac", lw=1.1)
        ax.fill_between(frame.time_ps, frame.CJJ_normalized - frame.CJJ_normalized_SEM,
                        frame.CJJ_normalized + frame.CJJ_normalized_SEM, color="#2166ac", alpha=0.18, lw=0)
        ax.axhline(0, color="0.45", lw=0.8)
        ax.text(0.04, 0.90, rf"$k={frame.k_inv_A.iloc[0]:.4f}\ \AA^{{-1}}$", transform=ax.transAxes)
        ax.tick_params(top=True, right=True)
        ax.text(-0.20, 1.06, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.supxlabel(r"Lag time, $t$ (ps)")
    fig.supylabel(r"$C_{J_zJ_z}(k,t)/C_{J_zJ_z}(k,0)$")
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.13, top=0.96, wspace=0.18, hspace=0.20)
    export(fig, output_fig / "C_JzJz_10L_n001_n006_2x3")
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.3), sharex=True)
    for i, n in enumerate(range(1, 7)):
        ax = axes.flat[i]
        frame = spectral_signed[(spectral_signed.n == n) & (np.abs(spectral_signed.omega_rad_ps) <= 1.25)]
        ax.semilogy(frame.omega_rad_ps, frame.PSD_mean_arbitrary, color="#b2182b", lw=1.05)
        ax.axvline(0, color="0.45", lw=0.8)
        ax.text(0.04, 0.90, rf"$k={frame.k_inv_A.iloc[0]:.4f}\ \AA^{{-1}}$", transform=ax.transAxes)
        ax.tick_params(top=True, right=True)
        ax.text(-0.20, 1.06, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.supxlabel(r"Angular frequency, $\omega$ (rad ps$^{-1}$)")
    fig.supylabel(r"$C(k,\omega)$ (arb. units)")
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.13, top=0.96, wspace=0.18, hspace=0.20)
    export(fig, output_fig / "C_k_omega_JzJz_10L_n001_n006_2x3")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.plot(peaks.k_inv_A, peaks.omega_peak_rad_ps, "o", ms=4, color="#2166ac", label="spectral peak")
    xx = np.linspace(0, peaks.k_inv_A.max() * 1.08, 200)
    ax.plot(xx, speed * xx + intercept, color="#b2182b", lw=1.1,
            label=rf"$\omega=({speed:.2f}\ \AA\,\mathrm{{ps}}^{{-1}})k{intercept:+.3f}$")
    ax.set(xlabel=r"$k$ ($\AA^{-1}$)", ylabel=r"$\omega_\mathrm{LA}$ (rad ps$^{-1}$)")
    ax.legend(loc="upper left", fontsize=6)
    ax.tick_params(top=True, right=True)
    ax.text(-0.18, 1.05, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.18, top=0.96)
    export(fig, output_fig / "LA_lowk_dispersion_n001_n006")
    plt.close(fig)
    (args.output / "metadata.json").write_text(json.dumps({
        "source": str(args.source), "protocol": "(8,8) 10L, 100 fs, 10 ns, 8 replicas, weak-NH",
        "CJJ": "local archived ensemble mean and replica SEM; normalized by the n-specific CJJ(0)",
        "spectrum": "local positive-frequency Welch PSD reflected as S(omega)=S(-omega)",
        "dispersion": "maxima near expected acoustic branch, n=1..6",
    }, indent=2))
    (args.output / "FINISHED.txt").write_text("Rendered from archived local 10L CJJ/spectral assets.\n")


if __name__ == "__main__":
    main()
