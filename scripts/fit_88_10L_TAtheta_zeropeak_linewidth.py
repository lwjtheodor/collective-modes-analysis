#!/usr/bin/env python3
"""Fit the zero-frequency TA_theta current peak for (8,8) 10L water.

The circumferential transverse response is treated as a quasi-elastic,
zero-centred relaxation peak.  This is deliberately different from a
finite-frequency acoustic-peak fit: omega_0 is fixed to zero and Gamma is the
HWHM of b + A Gamma^2/(omega^2 + Gamma^2).
"""

from __future__ import annotations

import argparse
import csv
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
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.8, "axes.spines.right": False,
    "axes.spines.top": False, "xtick.direction": "out", "ytick.direction": "out",
    "svg.fonttype": "none", "pdf.fonttype": 42, "legend.frameon": False,
})

ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
SPECTRUM = ROOT / "results" / "collective_mode_response" / "88_10L_per_k_semilog_Skw_pm10_LA_TAr_TAtheta" / "2026-08-19" / "derived_data" / "signed_Skw_n001_n160.csv"
OUT = ROOT / "results" / "collective_mode_response" / "88_10L_TAtheta_zeropeak_linewidth_lowk_n003_n010" / "2026-08-24"
COMMON = ROOT / "results" / "collective_mode_response" / "assets" / "current_mode_spectra" / "8_8_10L_10fs_1ns_fullvelocity_3rep"
OMEGA_FIT_MAX = 3.0
NMAX = 10
MIN_R2 = 0.90


def zero_lorentzian(omega, baseline, amplitude, gamma):
    return baseline + amplitude * gamma**2 / (omega**2 + gamma**2)


def load_spectra():
    grouped = defaultdict(list)
    with SPECTRUM.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["branch"] != "TA_theta" or int(row["n"]) > NMAX:
                continue
            omega = float(row["omega_rad_ps"])
            if abs(omega) <= OMEGA_FIT_MAX + 1e-12:
                grouped[int(row["n"])].append((omega, float(row["S_mean_arbitrary"]), float(row["S_replica_SEM_arbitrary"]), float(row["k_inv_A"])))
    for values in grouped.values():
        values.sort(key=lambda item: item[0])
    return grouped


def fit_mode(n, records):
    values = np.asarray(records, float)
    omega, signal, sem = values[:, 0], values[:, 1], values[:, 2]
    dw = float(np.median(np.diff(omega)))
    peak, baseline = float(signal.max()), max(float(signal.min()), 1e-14)
    sigma = np.maximum(sem, max(peak * 0.01, 1e-14))
    lower = [0.0, 0.0, dw / 2.0]
    upper = [np.inf, np.inf, OMEGA_FIT_MAX]
    params, covariance = curve_fit(
        zero_lorentzian, omega, signal, p0=[baseline, max(peak - baseline, baseline), 0.20],
        sigma=sigma, absolute_sigma=False, bounds=(lower, upper), maxfev=100000,
    )
    prediction = zero_lorentzian(omega, *params)
    ss_total = float(np.sum((signal - signal.mean()) ** 2))
    r2 = 1.0 - float(np.sum((signal - prediction) ** 2)) / ss_total
    gamma_sem = math.sqrt(float(covariance[2, 2])) if covariance[2, 2] > 0.0 else np.nan
    bound_hit = bool(params[2] <= 1.02 * lower[2] or params[2] >= 0.98 * upper[2])
    resolved = bool(params[2] >= 1.5 * dw)
    accepted = bool(r2 >= MIN_R2 and resolved and not bound_hit)
    return {
        "n": n, "k_inv_A": float(values[0, 3]), "frequency_bin_rad_ps": dw,
        "fit_omega_abs_max_rad_ps": OMEGA_FIT_MAX, "baseline": float(params[0]),
        "amplitude": float(params[1]), "gamma_HWHM_rad_ps": float(params[2]),
        "gamma_fit_SEM_rad_ps": gamma_sem, "FWHM_rad_ps": float(2.0 * params[2]),
        "fit_R2": r2, "bound_hit": bound_hit, "width_frequency_resolved": resolved,
        "accepted_linewidth": accepted, "omega": omega, "signal": signal, "sem": sem,
        "prediction": prediction,
    }


def power_fit(rows, label):
    k = np.asarray([r["k_inv_A"] for r in rows])
    gamma = np.asarray([r["gamma_HWHM_rad_ps"] for r in rows])
    if np.any(k <= 0.0) or np.any(gamma <= 0.0):
        raise ValueError("Power-law fitting requires strictly positive k and linewidth.")
    x, y = np.log(k), np.log(gamma)
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (intercept + slope * x)
    sxx = float(np.sum((x - x.mean()) ** 2))
    slope_se = math.sqrt(float(np.sum(residual**2) / (x.size - 2)) / sxx)
    ci = float(t.ppf(0.975, x.size - 2) * slope_se)
    r2 = 1.0 - float(np.sum(residual**2)) / float(np.sum((y - y.mean()) ** 2))
    intercept_k2 = float(np.mean(y - 2.0 * x))
    rss_k2 = float(np.sum((y - (intercept_k2 + 2.0 * x)) ** 2))
    return {
        "fit_label": label, "n_min": min(r["n"] for r in rows), "n_max": max(r["n"] for r in rows),
        "n_points": len(rows), "k_min_inv_A": float(k.min()), "k_max_inv_A": float(k.max()),
        "alpha_free": float(slope), "alpha_free_95CI_halfwidth": ci,
        "prefactor_free_rad_ps_A_to_alpha": float(math.exp(intercept)), "loglog_R2_free": r2,
        "RSS_log_free": float(np.sum(residual**2)), "prefactor_k2_rad_ps_A2": float(math.exp(intercept_k2)),
        "RSS_log_k2": rss_k2,
    }


def friction_k2_fit(rows, label):
    """Fit a wall-friction intercept plus transverse k^2 broadening."""
    k = np.asarray([r["k_inv_A"] for r in rows])
    gamma = np.asarray([r["gamma_HWHM_rad_ps"] for r in rows])
    design = np.column_stack((np.ones_like(k), k**2))
    params, _, _, _ = np.linalg.lstsq(design, gamma, rcond=None)
    residual = gamma - design @ params
    covariance = (residual @ residual / (k.size - 2)) * np.linalg.inv(design.T @ design)
    ci = t.ppf(0.975, k.size - 2) * np.sqrt(np.diag(covariance))
    r2 = 1.0 - float(residual @ residual) / float(np.sum((gamma - gamma.mean()) ** 2))
    return {
        "fit_label": label, "n_min": min(r["n"] for r in rows), "n_max": max(r["n"] for r in rows),
        "n_points": len(rows), "Gamma0_rad_ps": float(params[0]), "Gamma0_95CI_halfwidth_rad_ps": float(ci[0]),
        "D_k2_rad_ps_A2": float(params[1]), "D_k2_95CI_halfwidth_rad_ps_A2": float(ci[1]),
        "R2_linear": r2, "RSS_linear": float(residual @ residual),
    }


def save(fig, figdir, stem):
    fig.savefig(figdir / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(figdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(figdir / f"{stem}.tiff", dpi=600, bbox_inches="tight")


def plot_powerlaw(rows, primary, sensitivity, figdir):
    fig = plt.figure(figsize=(7.2, 3.15))
    ax = fig.add_axes([0.10, 0.19, 0.42, 0.74])
    ax2 = fig.add_axes([0.64, 0.19, 0.29, 0.74])
    good = [r for r in rows if r["accepted_linewidth"]]
    bad = [r for r in rows if not r["accepted_linewidth"]]
    ax.loglog([r["k_inv_A"] for r in bad], [r["gamma_HWHM_rad_ps"] for r in bad], "o", ms=3.5, mfc="white", mec="0.45", label="not resolved")
    ax.errorbar([r["k_inv_A"] for r in good], [r["gamma_HWHM_rad_ps"] for r in good],
                yerr=[max(r["gamma_fit_SEM_rad_ps"], 0.1 * r["gamma_HWHM_rad_ps"]) for r in good],
                fmt="o", ms=4.0, color="#0072b2", ecolor="#0072b2", elinewidth=0.7, capsize=1.4, label="zero-peak HWHM")
    kline = np.geomspace(primary["k_min_inv_A"], primary["k_max_inv_A"], 200)
    ax.loglog(kline, primary["prefactor_free_rad_ps_A_to_alpha"] * kline**primary["alpha_free"], color="#d55e00", lw=1.25, label=rf"free: $\Gamma\propto k^{{{primary['alpha_free']:.2f}}}$")
    ax.loglog(kline, primary["prefactor_k2_rad_ps_A2"] * kline**2, color="0.25", lw=0.9, ls="--", label=r"reference: $k^2$")
    ax.set_xlabel(r"$k$ ($\mathrm{\AA}^{-1}$)")
    ax.set_ylabel(r"TA$_\theta$ zero-peak HWHM, $\Gamma$ (rad ps$^{-1}$)")
    ax.legend(loc="upper left", fontsize=6.0, handlelength=1.9)
    ax.tick_params(which="both", labelsize=6.5, length=3)
    ax.tick_params(which="minor", length=1.8)
    ax.text(-0.18, 1.03, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.text(0.03, 0.04, f"primary: n={primary['n_min']}-{primary['n_max']}\nalpha={primary['alpha_free']:.2f} +/- {primary['alpha_free_95CI_halfwidth']:.2f} (conditional 95% CI)", transform=ax.transAxes, fontsize=6.1, va="bottom")
    ax2.axhline(2.0, color="0.25", lw=0.9, ls="--")
    ax2.errorbar([r["n_min"] for r in sensitivity], [r["alpha_free"] for r in sensitivity],
                 yerr=[r["alpha_free_95CI_halfwidth"] for r in sensitivity], fmt="o-", color="#0072b2", ms=3.8, lw=1.0, capsize=1.5)
    ax2.set_xlabel(r"lower mode, $n_{\min}$")
    ax2.set_ylabel(r"free exponent, $\alpha$")
    ax2.set_xticks([r["n_min"] for r in sensitivity])
    ax2.set_ylim(0.6, 2.3)
    ax2.tick_params(labelsize=6.5, length=3)
    ax2.text(-0.24, 1.03, "(b)", transform=ax2.transAxes, fontweight="bold", fontsize=9)
    ax2.text(0.97, 2.03, r"$k^2$", ha="right", va="bottom", fontsize=6.5, color="0.25")
    save(fig, figdir, "TA_theta_zeropeak_linewidth_powerlaw")
    plt.close(fig)


def plot_examples(rows, figdir):
    examples = [next(r for r in rows if r["n"] == n) for n in (1, 3, 6, 10)]
    fig = plt.figure(figsize=(7.2, 3.0))
    left, bottom, width, height, gap = 0.09, 0.19, 0.89, 0.71, 0.045
    panel_w = (width - 3 * gap) / 4
    for j, row in enumerate(examples):
        ax = fig.add_axes([left + j * (panel_w + gap), bottom, panel_w, height])
        ax.plot(row["omega"], row["signal"], color="0.35", lw=0.9, label="Welch mean")
        ax.fill_between(row["omega"], np.maximum(row["signal"] - row["sem"], 0), row["signal"] + row["sem"], color="0.7", alpha=0.35, lw=0)
        ax.plot(row["omega"], row["prediction"], color="#d55e00", lw=1.1, label=r"$\omega_0=0$ Lorentzian")
        ax.axvline(0.0, color="0.25", lw=0.6, ls="--")
        ax.set_xlim(-OMEGA_FIT_MAX, OMEGA_FIT_MAX)
        ax.set_xlabel(r"$\omega$ (rad ps$^{-1}$)")
        if j == 0:
            ax.set_ylabel(r"$S_{J_\theta J_\theta}$ (arb.)")
            ax.legend(loc="upper right", fontsize=5.5)
        ax.set_title(rf"$n={row['n']}$, $\Gamma={row['gamma_HWHM_rad_ps']:.3f}$", fontsize=6.6, pad=3)
        ax.tick_params(labelsize=6.0, length=3)
        ax.text(-0.24, 1.03, f"({chr(97 + j)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    save(fig, figdir, "TA_theta_zeropeak_lorentzian_fit_examples")
    plt.close(fig)


def plot_friction_k2(rows, fit, figdir):
    fig = plt.figure(figsize=(3.35, 2.75))
    ax = fig.add_axes([0.19, 0.20, 0.76, 0.73])
    use = [r for r in rows if r["accepted_linewidth"] and 3 <= r["n"] <= 10]
    x = np.asarray([r["k_inv_A"] ** 2 for r in use])
    y = np.asarray([r["gamma_HWHM_rad_ps"] for r in use])
    e = np.asarray([max(r["gamma_fit_SEM_rad_ps"], 0.1 * r["gamma_HWHM_rad_ps"]) for r in use])
    ax.errorbar(x, y, yerr=e, fmt="o", ms=4.0, color="#0072b2", ecolor="#0072b2", elinewidth=0.7, capsize=1.4)
    line = np.linspace(0.0, x.max() * 1.05, 200)
    ax.plot(line, fit["Gamma0_rad_ps"] + fit["D_k2_rad_ps_A2"] * line, color="#d55e00", lw=1.25)
    ax.axhline(fit["Gamma0_rad_ps"], color="0.25", lw=0.8, ls="--")
    ax.set_xlabel(r"$k^2$ ($\mathrm{\AA}^{-2}$)")
    ax.set_ylabel(r"TA$_\theta$ HWHM, $\Gamma$ (rad ps$^{-1}$)")
    ax.text(0.04, 0.96, rf"$\Gamma=\Gamma_0+D_\theta k^2$" + "\n" + rf"$\Gamma_0={fit['Gamma0_rad_ps']:.3f}\pm{fit['Gamma0_95CI_halfwidth_rad_ps']:.3f}$ rad ps$^{{-1}}$" + "\n" + rf"$D_\theta={fit['D_k2_rad_ps_A2']:.1f}\pm{fit['D_k2_95CI_halfwidth_rad_ps_A2']:.1f}$ rad ps$^{{-1}}$ $\mathrm{{\AA}}^2$", transform=ax.transAxes, va="top", fontsize=6.2)
    ax.tick_params(labelsize=6.5, length=3)
    ax.text(-0.23, 1.03, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    save(fig, figdir, "TA_theta_zeropeak_friction_intercept_plus_k2")
    plt.close(fig)


def main():
    global SPECTRUM
    parser = argparse.ArgumentParser(description="Fit TA_theta zero-frequency linewidth scaling.")
    parser.add_argument("--spectrum-csv", type=Path, default=SPECTRUM)
    parser.add_argument("--outdir", type=Path, default=OUT)
    parser.add_argument("--common-asset-dir", type=Path, default=COMMON)
    args = parser.parse_args()
    figdir, datadir = args.outdir / "figures", args.outdir / "derived_data"
    figdir.mkdir(parents=True, exist_ok=True)
    datadir.mkdir(parents=True, exist_ok=True)
    args.common_asset_dir.mkdir(parents=True, exist_ok=True)
    SPECTRUM = args.spectrum_csv
    spectra = load_spectra()
    rows = [fit_mode(n, spectra[n]) for n in range(1, NMAX + 1)]
    fields = [key for key in rows[0] if key not in {"omega", "signal", "sem", "prediction"}]
    with (datadir / "TA_theta_zeropeak_lorentzian_fits_n001_n010.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows([{key: row[key] for key in fields} for row in rows])
    primary_rows = [row for row in rows if row["accepted_linewidth"] and 3 <= row["n"] <= 10]
    primary = power_fit(primary_rows, "primary_resolved_n003_n010")
    sensitivity = [power_fit([row for row in rows if row["accepted_linewidth"] and low <= row["n"] <= 10], f"sensitivity_n{low:03d}_n010") for low in (2, 3)]
    with (datadir / "TA_theta_zeropeak_linewidth_powerlaw_fits.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(primary)); writer.writeheader(); writer.writerow(primary); writer.writerows(sensitivity)
    plot_powerlaw(rows, primary, sensitivity, figdir)
    plot_examples(rows, figdir)
    friction_rows = [row for row in rows if row["accepted_linewidth"] and 3 <= row["n"] <= 10]
    friction = friction_k2_fit(friction_rows, "primary_resolved_n003_n010")
    friction_sensitivity = friction_k2_fit([row for row in rows if row["accepted_linewidth"] and 2 <= row["n"] <= 10], "sensitivity_n002_n010")
    with (datadir / "TA_theta_zeropeak_friction_intercept_plus_k2_fits.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(friction)); writer.writeheader(); writer.writerow(friction); writer.writerow(friction_sensitivity)
    plot_friction_k2(rows, friction, figdir)
    manifest = {
        "analysis_date": str(date.today()), "system": "(8,8) CNT-confined water, 10L; full vx/vy/vz, 10 fs, 1 ns, three velocity-seed replicas",
        "branch": "TA_theta", "observable": "circumferential-current zero-frequency (quasi-elastic) peak",
        "spectrum_source": str(SPECTRUM), "fit_model": "b + A Gamma^2/(omega^2 + Gamma^2), with omega0 fixed to zero",
        "fit_abs_omega_max_rad_ps": OMEGA_FIT_MAX,
        "selection": {"fitted_n": [1, 10], "primary_n": [3, 10], "minimum_HWHM_bins": 1.5, "fit_R2_min": MIN_R2},
        "primary_fit": primary,
        "friction_intercept_plus_k2_fit": friction,
        "uncertainty_limit": "The 95% CI is conditional OLS scatter across mode indices; it is not a replica-resolved confidence interval.",
        "radial_scope_note": "TA_r has no low-frequency central peak in this projection; its low-frequency spectral density is suppressed and must not be fitted with this zero-peak model.",
    }
    (args.outdir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (args.common_asset_dir / "asset_manifest_TA_theta_zeropeak_linewidth.json").write_text(json.dumps({"asset_type": "zero-peak transverse-current linewidth", "protocol": manifest["system"], "source_spectrum": str(SPECTRUM), "analysis_archive": str(args.outdir), "note": "TA_theta only; TA_r is excluded because its low-frequency spectrum has no central peak."}, indent=2) + "\n")
    (args.outdir / "FINISHED.txt").write_text("10L TA_theta zero-peak linewidth analysis finished successfully.\n")
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
