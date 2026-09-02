#!/usr/bin/env python3
"""Fit low-k current-spectrum linewidths for (8,8) 10L water.

Each positive-frequency LA peak is fit to a Lorentzian on a linear local
background.  The reported damping linewidth is the fitted HWHM Gamma.  The
power-law exponent is fitted, not constrained to the hydrodynamic k^2 value.
"""

from __future__ import annotations

import csv
import argparse
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "legend.frameon": False,
})

ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
SPECTRUM_CSV = ROOT / "results" / "collective_mode_response" / "88_10L_per_k_semilog_Skw_pm10_LA_TAr_TAtheta" / "2026-08-19" / "derived_data" / "signed_Skw_n001_n160.csv"
PEAK_CSV = ROOT / "results" / "collective_mode_response" / "88_10L_per_k_semilog_Skw_LA_TAr_TAtheta" / "2026-08-19" / "derived_data" / "phase_and_group_velocity_vs_k.csv"
OUT = ROOT / "results" / "collective_mode_response" / "88_10L_LA_linewidth_powerlaw" / "2026-08-24"
FIG = OUT / "figures"
DATA = OUT / "derived_data"
COMMON_ASSET = ROOT / "results" / "collective_mode_response" / "assets" / "current_mode_spectra" / "8_8_10L_10fs_1ns_fullvelocity_3rep"

# The low-k acoustic window is declared before fitting: n=1...20 maps to
# k<=0.125 A^-1.  n=1,2 are retained in the table but excluded when their
# HWHM is below 1.5 spectral bins (not frequency-resolved).
N_FIT_MIN = 1
N_FIT_MAX = 10
N_ACOUSTIC_MIN = 3
N_ACOUSTIC_MAX = 10
MIN_R2 = 0.70
BRANCH = "LA"


def branch_label(branch):
    return {"LA": "LA", "TA_r": r"TA$_r$", "TA_theta": r"TA$_\\theta$"}[branch]


def current_label(branch):
    return {"LA": r"$S_{J_zJ_z}$", "TA_r": r"$S_{J_rJ_r}$", "TA_theta": r"$S_{J_\theta J_\theta}$"}[branch]


def lorentz_linear(w, b0, b1, amp, w0, gamma):
    """Positive Lorentzian peak with a local linear background."""
    return b0 + b1 * (w - w0) + amp * gamma**2 / ((w - w0) ** 2 + gamma**2)


def load_spectra():
    spectra = defaultdict(list)
    with SPECTRUM_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["branch"] != BRANCH:
                continue
            omega = float(row["omega_rad_ps"])
            if omega >= 0.0:
                spectra[int(row["n"])].append((omega, float(row["S_mean_arbitrary"]), float(row["S_replica_SEM_arbitrary"]), float(row["k_inv_A"])))
    for records in spectra.values():
        records.sort(key=lambda r: r[0])
    return spectra


def load_seed_centres():
    centres = {}
    with PEAK_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["branch"] == BRANCH:
                omega_key = "omega_peak_rad_ps" if "omega_peak_rad_ps" in row else "omega_peak_mean_rad_ps"
                centres[int(row["n"])] = (float(row[omega_key]), row["resolved_operational_peak"].lower() == "true")
    return centres


def fit_one_mode(n, records, peak_seed, seed_resolved):
    values = np.asarray(records, float)
    k = float(values[0, 3])
    dw = float(np.median(np.diff(values[:, 0])))
    # Wide enough to retain both flanks, capped before unrelated high-frequency
    # structures become dominant.
    half_window = max(0.25, min(1.5, 0.70 * peak_seed))
    selected = values[(values[:, 0] >= max(0.0, peak_seed - half_window)) & (values[:, 0] <= peak_seed + half_window)]
    omega, signal, sem = selected[:, 0], selected[:, 1], selected[:, 2]
    if omega.size < 9:
        raise RuntimeError("too few points in local peak window")
    peak = float(np.max(signal))
    baseline = max(float(np.min(signal)), 1e-14)
    sigma = np.maximum(sem, max(peak * 0.01, 1e-14))
    gamma0 = max(0.04, 0.20 * peak_seed)
    lower = [0.0, -np.inf, 0.0, max(0.0, peak_seed - half_window), dw / 2.0]
    upper = [np.inf, np.inf, np.inf, peak_seed + half_window, half_window]
    params, covariance = curve_fit(
        lorentz_linear, omega, signal,
        p0=[baseline, 0.0, max(peak - baseline, baseline), peak_seed, gamma0],
        sigma=sigma, absolute_sigma=False, bounds=(lower, upper), maxfev=100000,
    )
    prediction = lorentz_linear(omega, *params)
    ss_total = float(np.sum((signal - np.mean(signal)) ** 2))
    r2 = 1.0 - float(np.sum((signal - prediction) ** 2)) / ss_total if ss_total > 0 else np.nan
    gamma_sem = math.sqrt(float(covariance[4, 4])) if covariance[4, 4] > 0.0 else np.nan
    bound_hit = bool(params[4] <= 1.02 * lower[4] or params[4] >= 0.98 * upper[4])
    resolved_width = bool(params[4] >= 1.5 * dw)
    accepted = bool(seed_resolved and r2 >= MIN_R2 and resolved_width and not bound_hit)
    return {
        "n": n, "k_inv_A": k, "seed_omega_rad_ps": peak_seed, "seed_resolved": seed_resolved,
        "fit_window_low_rad_ps": float(omega[0]), "fit_window_high_rad_ps": float(omega[-1]),
        "frequency_bin_rad_ps": dw, "b0": float(params[0]), "b1": float(params[1]),
        "amplitude": float(params[2]), "omega0_fit_rad_ps": float(params[3]),
        "gamma_HWHM_rad_ps": float(params[4]), "gamma_fit_SEM_rad_ps": gamma_sem,
        "FWHM_rad_ps": float(2.0 * params[4]), "fit_R2": r2, "bound_hit": bound_hit,
        "width_frequency_resolved": resolved_width, "accepted_linewidth": accepted,
        "omega": omega, "signal": signal, "sem": sem, "prediction": prediction,
    }


def log_power_fit(rows, label):
    k = np.asarray([r["k_inv_A"] for r in rows])
    gamma = np.asarray([r["gamma_HWHM_rad_ps"] for r in rows])
    if np.any(k <= 0.0) or np.any(gamma <= 0.0):
        raise ValueError("Power-law fit requires strictly positive k and linewidth")
    x, y = np.log(k), np.log(gamma)
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (intercept + slope * x)
    sxx = float(np.sum((x - x.mean()) ** 2))
    sigma2 = float(np.sum(residual**2) / (x.size - 2))
    slope_se = math.sqrt(sigma2 / sxx)
    ci = float(t.ppf(0.975, x.size - 2) * slope_se)
    r2 = 1.0 - float(np.sum(residual**2)) / float(np.sum((y - y.mean()) ** 2))
    # Fixed-exponent k^2 reference is fit only in amplitude, on identical data.
    intercept_k2 = float(np.mean(y - 2.0 * x))
    residual_k2 = y - (intercept_k2 + 2.0 * x)
    return {
        "fit_label": label, "n_min": min(r["n"] for r in rows), "n_max": max(r["n"] for r in rows),
        "n_points": len(rows), "k_min_inv_A": float(k.min()), "k_max_inv_A": float(k.max()),
        "alpha_free": float(slope), "alpha_free_95CI_halfwidth": ci,
        "prefactor_free_rad_ps_A_to_alpha": float(math.exp(intercept)), "loglog_R2_free": r2,
        "RSS_log_free": float(np.sum(residual**2)), "prefactor_k2_rad_ps_A2": float(math.exp(intercept_k2)),
        "RSS_log_k2": float(np.sum(residual_k2**2)),
    }


def write_csv(rows, path):
    names = [
        "n", "k_inv_A", "seed_omega_rad_ps", "seed_resolved", "fit_window_low_rad_ps", "fit_window_high_rad_ps",
        "frequency_bin_rad_ps", "b0", "b1", "amplitude", "omega0_fit_rad_ps", "gamma_HWHM_rad_ps",
        "gamma_fit_SEM_rad_ps", "FWHM_rad_ps", "fit_R2", "bound_hit", "width_frequency_resolved", "accepted_linewidth",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row[name] for name in names})


def save_figure(fig, stem):
    fig.savefig(FIG / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.tiff", dpi=600, bbox_inches="tight")


def plot_powerlaw(rows, primary, sensitivity):
    fig = plt.figure(figsize=(7.2, 3.15))
    ax = fig.add_axes([0.10, 0.19, 0.42, 0.74])
    ax2 = fig.add_axes([0.64, 0.19, 0.29, 0.74])
    accepted = [r for r in rows if r["accepted_linewidth"]]
    rejected = [r for r in rows if not r["accepted_linewidth"]]
    ax.loglog([r["k_inv_A"] for r in rejected], [r["gamma_HWHM_rad_ps"] for r in rejected], "o", ms=3.5, mfc="white", mec="0.45", mew=0.8, label="not used")
    ax.errorbar([r["k_inv_A"] for r in accepted], [r["gamma_HWHM_rad_ps"] for r in accepted],
                yerr=[max(r["gamma_fit_SEM_rad_ps"], 0.1 * r["gamma_HWHM_rad_ps"]) for r in accepted],
                fmt="o", ms=4.0, color="#1769aa", ecolor="#1769aa", elinewidth=0.7, capsize=1.4, label="resolved linewidth")
    kline = np.geomspace(primary["k_min_inv_A"], primary["k_max_inv_A"], 200)
    ax.loglog(kline, primary["prefactor_free_rad_ps_A_to_alpha"] * kline ** primary["alpha_free"], color="#d55e00", lw=1.25,
              label=rf"free: $\Gamma\propto k^{{{primary['alpha_free']:.2f}}}$")
    ax.loglog(kline, primary["prefactor_k2_rad_ps_A2"] * kline**2, color="0.25", lw=0.9, ls="--", label=r"reference: $k^2$")
    ax.set_xlabel(r"$k$ ($\mathrm{\AA}^{-1}$)")
    ax.set_ylabel(rf"{branch_label(BRANCH)} HWHM, $\Gamma$ (rad ps$^{{-1}}$)")
    ax.legend(loc="upper left", fontsize=6.0, handlelength=1.9)
    ax.tick_params(which="both", labelsize=6.5, length=3)
    ax.tick_params(which="minor", length=1.8)
    ax.text(-0.18, 1.03, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.text(0.03, 0.04, f"low-k fit: n={primary['n_min']}-{primary['n_max']}\nalpha={primary['alpha_free']:.2f} +/- {primary['alpha_free_95CI_halfwidth']:.2f} (conditional 95% CI)", transform=ax.transAxes, fontsize=6.1, va="bottom")

    labels = [s["n_max"] for s in sensitivity]
    alphas = [s["alpha_free"] for s in sensitivity]
    errs = [s["alpha_free_95CI_halfwidth"] for s in sensitivity]
    ax2.axhline(2.0, color="0.25", lw=0.9, ls="--")
    ax2.errorbar(labels, alphas, yerr=errs, fmt="o-", color="#1769aa", ms=3.8, lw=1.0, capsize=1.5)
    ax2.set_xlabel(r"upper mode, $n_{\max}$")
    ax2.set_ylabel(r"free exponent, $\alpha$")
    ax2.set_xticks(labels)
    ax2.set_ylim(0.8, 2.3)
    ax2.tick_params(labelsize=6.5, length=3)
    ax2.text(-0.24, 1.03, "(b)", transform=ax2.transAxes, fontweight="bold", fontsize=9)
    ax2.text(0.97, 2.03, r"$k^2$", ha="right", va="bottom", fontsize=6.5, color="0.25")
    save_figure(fig, f"{BRANCH}_linewidth_powerlaw_and_range_sensitivity")
    plt.close(fig)


def plot_examples(rows):
    examples = [next(r for r in rows if r["n"] == n) for n in (3, 4, 6, 10)]
    fig = plt.figure(figsize=(7.2, 3.0))
    left, bottom, width, height = 0.09, 0.19, 0.89, 0.71
    gap = 0.045
    panel_w = (width - 3 * gap) / 4
    for j, row in enumerate(examples):
        ax = fig.add_axes([left + j * (panel_w + gap), bottom, panel_w, height])
        ax.plot(row["omega"], row["signal"], color="0.35", lw=0.9, label="Welch mean")
        ax.fill_between(row["omega"], np.maximum(row["signal"] - row["sem"], 0), row["signal"] + row["sem"], color="0.7", alpha=0.35, lw=0)
        ax.plot(row["omega"], row["prediction"], color="#d55e00", lw=1.1, label="Lorentzian + baseline")
        ax.axvline(row["omega0_fit_rad_ps"], color="#d55e00", lw=0.7, ls="--")
        ax.set_xlim(row["fit_window_low_rad_ps"], row["fit_window_high_rad_ps"])
        ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
        if j == 0:
            ax.set_ylabel(current_label(BRANCH) + " (arb.)")
            ax.legend(loc="upper left", fontsize=5.8)
        ax.set_title(rf"$n={row['n']}$, $k={row['k_inv_A']:.3f}\;\mathrm{{\AA}}^{{-1}}$", fontsize=6.6, pad=3)
        ax.tick_params(labelsize=6.0, length=3)
        ax.text(-0.24, 1.03, f"({chr(97+j)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    save_figure(fig, f"{BRANCH}_linewidth_lorentzian_fit_examples")
    plt.close(fig)


def main():
    global SPECTRUM_CSV, PEAK_CSV, OUT, FIG, DATA, COMMON_ASSET, N_FIT_MAX, N_ACOUSTIC_MIN, N_ACOUSTIC_MAX, BRANCH
    parser = argparse.ArgumentParser(description="Fit (8,8) 10L current-spectrum linewidth power law from archived spectra.")
    parser.add_argument("--spectrum-csv", type=Path, default=SPECTRUM_CSV)
    parser.add_argument("--peak-csv", type=Path, default=PEAK_CSV)
    parser.add_argument("--outdir", type=Path, default=OUT)
    parser.add_argument("--common-asset-dir", type=Path, default=COMMON_ASSET)
    parser.add_argument("--n-fit-max", type=int, default=N_FIT_MAX, help="Largest mode subject to individual peak fitting.")
    parser.add_argument("--primary-n-min", type=int, default=N_ACOUSTIC_MIN)
    parser.add_argument("--primary-n-max", type=int, default=N_ACOUSTIC_MAX, help="Largest mode entering the acoustic power-law conclusion.")
    parser.add_argument("--branch", choices=("LA", "TA_r", "TA_theta"), default=BRANCH)
    args = parser.parse_args()
    if not (1 <= args.primary_n_min <= args.primary_n_max <= args.n_fit_max):
        raise ValueError("Require 1 <= primary-n-min <= primary-n-max <= n-fit-max")
    SPECTRUM_CSV, PEAK_CSV, OUT, COMMON_ASSET, BRANCH = args.spectrum_csv, args.peak_csv, args.outdir, args.common_asset_dir, args.branch
    N_FIT_MAX, N_ACOUSTIC_MIN, N_ACOUSTIC_MAX = args.n_fit_max, args.primary_n_min, args.primary_n_max
    FIG, DATA = OUT / "figures", OUT / "derived_data"
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    COMMON_ASSET.mkdir(parents=True, exist_ok=True)
    spectra = load_spectra()
    seed_centres = load_seed_centres()
    rows = []
    for n in range(N_FIT_MIN, N_FIT_MAX + 1):
        row = fit_one_mode(n, spectra[n], *seed_centres[n])
        rows.append(row)
    write_csv(rows, DATA / f"{BRANCH}_peak_linewidth_lorentzian_fits_n001_n{N_FIT_MAX:03d}.csv")

    primary_rows = [r for r in rows if r["accepted_linewidth"] and N_ACOUSTIC_MIN <= r["n"] <= N_ACOUSTIC_MAX]
    if len(primary_rows) < 5:
        raise RuntimeError("Too few resolved acoustic linewidths for a power-law fit")
    primary = log_power_fit(primary_rows, f"primary_low_k_n{N_ACOUSTIC_MIN:03d}_n{N_ACOUSTIC_MAX:03d}")
    sensitivity = []
    sensitivity_nmax = sorted(set(n for n in (6, 8, N_ACOUSTIC_MAX) if N_ACOUSTIC_MIN + 3 <= n <= N_ACOUSTIC_MAX))
    for nmax in sensitivity_nmax:
        subset = [r for r in rows if r["accepted_linewidth"] and N_ACOUSTIC_MIN <= r["n"] <= nmax]
        sensitivity.append(log_power_fit(subset, f"sensitivity_n003_n{nmax:03d}"))
    fields = list(primary)
    with (DATA / f"{BRANCH}_linewidth_powerlaw_fits.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(primary)
        writer.writerows(sensitivity)
    plot_powerlaw(rows, primary, sensitivity)
    plot_examples(rows)

    manifest = {
        "analysis_date": str(date.today()),
        "system": "(8,8) CNT-confined water, 10L; 10 fs sampling, 1 ns, three velocity-seed replicas",
        "remote_analysis_root": "/lustre/home/users/ewu/vb_gcmc/MD/full_velocity_10fs_1ns_assets_2L5L10L_20260819/allmode_cjj_analysis_20260819/outputs/88_L10",
        "input_spectrum_csv": str(SPECTRUM_CSV),
        "input_peak_csv": str(PEAK_CSV),
        "branch": BRANCH,
        "linewidth_definition": "HWHM Gamma from a local Lorentzian plus linear background fit to the positive-frequency current PSD",
        "selection": {"fit_n": [N_FIT_MIN, N_FIT_MAX], "primary_acoustic_n": [N_ACOUSTIC_MIN, N_ACOUSTIC_MAX], "criteria": {"fit_R2_min": MIN_R2, "minimum_HWHM_bins": 1.5, "exclude_bound_hit": True}},
        "primary_fit": primary,
        "uncertainty_limit": "95% CI is conditional OLS scatter across modes; replica-resolved PSDs were not archived, so it is not an independent-seed confidence interval.",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    asset_manifest = {
        "asset_type": "current-mode spectra",
        "protocol": "(8,8) 10L; full vx/vy/vz, 10 fs, 1 ns, 3 replicas",
        "remote_source": manifest["remote_analysis_root"],
        "local_archived_spectrum": str(SPECTRUM_CSV),
        "note": "No trajectory or raw-spectrum duplication: this common asset manifest indexes the archived reusable signed spectrum CSV.",
    }
    (COMMON_ASSET / f"asset_manifest_{BRANCH}_linewidth.json").write_text(json.dumps(asset_manifest, indent=2) + "\n")
    (OUT / "FINISHED.txt").write_text(f"10L {BRANCH} linewidth power-law analysis finished successfully.\n")
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
