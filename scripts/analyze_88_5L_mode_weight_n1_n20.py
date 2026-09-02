#!/usr/bin/env python3
"""Audit the n=1..20 modal spectral weight W_a(k)=C_JJ,a(k,0)/N for (8,8) 5L.

W_a is the equal-time modal-current variance and, by Wiener--Khinchin, is
proportional to the frequency-integrated modal-current spectrum in the same
normalization.  It is deliberately not called a single-particle VDOS.
"""
from __future__ import annotations

import csv
import io
import json
import math
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress


ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
SOURCE = ROOT / "results" / "collective_mode_response" / "88_5L_LA_TAr_TAtheta_dispersion" / "2026-08-19"
RAW_CJJ = SOURCE / "derived_data" / "CJJ_all_modes_per_replica.csv"
OUT = ROOT / "results" / "collective_mode_response" / "88_5L_mode_weight_n1_n20" / "2026-08-19"
N_OXYGEN = 665
NMAX = 20
BRANCHES = ("LA", "TA_r", "TA_theta")
COLORS = {"LA": "#0F4D92", "TA_r": "#B64342", "TA_theta": "#42949E"}
LABELS = {"LA": "LA", "TA_r": r"TA$_r$", "TA_theta": r"TA$_\theta$"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 1.0,
    "ytick.major.width": 1.0, "legend.frameon": False,
})


def read_lag_zero_rows() -> list[dict]:
    """Extract the 180 requested source points without loading 7.2 M rows."""
    rows: list[dict] = []
    # ripgrep streams the narrow first-column/lag filter in seconds on this
    # 7.2-million-row CSV, avoiding a many-minute Python CSV scan.
    pattern = r"^[123],(LA|TA_r|TA_theta),(?:[1-9]|1[0-9]|20),[^,]+,[^,]+,0\.0,"
    selected = subprocess.run(
        ["rg", "--no-heading", "--color", "never", pattern, str(RAW_CJJ)],
        check=True, capture_output=True, text=True, encoding="utf-8",
    ).stdout
    header = "replicate,branch,n,k_inv_A,wavelength_A,lag_ps,n_time_origins,CJJ_raw_A2_fs2,CJJ_normalized\n"
    for row in csv.DictReader(io.StringIO(header + selected)):
        rows.append({
            "replicate": int(row["replicate"]), "branch": row["branch"],
            "n": int(row["n"]), "k_inv_A": float(row["k_inv_A"]),
            "CJJ0_A2_fs2": float(row["CJJ_raw_A2_fs2"]),
            "W_per_oxygen_A2_fs2": float(row["CJJ_raw_A2_fs2"]) / N_OXYGEN,
        })
    expected = 3 * len(BRANCHES) * NMAX
    if len(rows) != expected:
        raise ValueError(f"expected {expected} lag-zero rows, found {len(rows)}")
    return rows


def grouped(rows: list[dict], branch: str, rep: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    subset = sorted((r for r in rows if r["branch"] == branch and r["replicate"] == rep), key=lambda r: r["n"])
    return (np.array([r["n"] for r in subset]), np.array([r["k_inv_A"] for r in subset]),
            np.array([r["W_per_oxygen_A2_fs2"] for r in subset]))


def fit_rows(rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for branch in BRANCHES:
        for rep in (1, 2, 3):
            n, k, weight = grouped(rows, branch, rep)
            if np.any(weight <= 0.0):
                raise ValueError(f"non-positive W encountered for {branch} replica {rep}; log-log fit undefined")
            linear = linregress(k, weight)
            linear_pred = linear.intercept + linear.slope * k
            ss_res = float(np.sum((weight - linear_pred) ** 2))
            ss_tot = float(np.sum((weight - weight.mean()) ** 2))
            logfit = linregress(np.log(k), np.log(weight))
            output.append({
                "branch": branch, "replicate": rep, "n_min": int(n.min()), "n_max": int(n.max()),
                "linear_intercept_A2_fs2_per_O": float(linear.intercept),
                "linear_slope_A3_fs2_per_O": float(linear.slope),
                "linear_slope_stderr_A3_fs2_per_O": float(linear.stderr),
                "linear_R2": float(linear.rvalue ** 2), "linear_slope_pvalue_within_replica": float(linear.pvalue),
                "linear_RMSE_A2_fs2_per_O": float(math.sqrt(ss_res / len(k))),
                "powerlaw_exponent": float(logfit.slope), "powerlaw_R2": float(logfit.rvalue ** 2),
                "interpretation": "descriptive fit across 20 discrete modes; replicate-to-replicate scatter, not within-mode OLS p, is the primary uncertainty",
            })
    for branch in BRANCHES:
        per = [r for r in output if r["branch"] == branch]
        slopes = np.array([r["linear_slope_A3_fs2_per_O"] for r in per])
        r2 = np.array([r["linear_R2"] for r in per])
        exponent = np.array([r["powerlaw_exponent"] for r in per])
        output.append({
            "branch": branch, "replicate": "mean_across_3", "n_min": 1, "n_max": NMAX,
            "linear_intercept_A2_fs2_per_O": float(np.mean([r["linear_intercept_A2_fs2_per_O"] for r in per])),
            "linear_slope_A3_fs2_per_O": float(slopes.mean()),
            "linear_slope_stderr_A3_fs2_per_O": float(slopes.std(ddof=1) / math.sqrt(3)),
            "linear_R2": float(r2.mean()), "linear_slope_pvalue_within_replica": "not_applicable",
            "linear_RMSE_A2_fs2_per_O": float(np.mean([r["linear_RMSE_A2_fs2_per_O"] for r in per])),
            "powerlaw_exponent": float(exponent.mean()), "powerlaw_R2": float(np.mean([r["powerlaw_R2"] for r in per])),
            "interpretation": "mean plus replica SEM across three velocity-seed replicas; not a thermodynamic-limit statement",
        })
    return output


def ensemble_rows(rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for branch in BRANCHES:
        for n in range(1, NMAX + 1):
            selected = [r for r in rows if r["branch"] == branch and r["n"] == n]
            weights = np.array([r["W_per_oxygen_A2_fs2"] for r in selected])
            output.append({
                "branch": branch, "n": n, "k_inv_A": selected[0]["k_inv_A"],
                "W_mean_A2_fs2_per_O": float(weights.mean()),
                "W_replica_SEM_A2_fs2_per_O": float(weights.std(ddof=1) / math.sqrt(3)),
                "W_fraction_within_n1_n20": float(weights.mean() / sum(
                    rr["W_per_oxygen_A2_fs2"] for rr in rows if rr["branch"] == branch and rr["replicate"] == 1
                )),
                "n_replicas": 3,
            })
    # Recompute the fraction from ensemble means so its denominator is explicit and exact.
    for branch in BRANCHES:
        subset = [r for r in output if r["branch"] == branch]
        denominator = sum(r["W_mean_A2_fs2_per_O"] for r in subset)
        for row in subset:
            row["W_fraction_within_n1_n20"] = row["W_mean_A2_fs2_per_O"] / denominator
    return output


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white")
    plt.close(fig)


def make_figure(rows: list[dict], ensemble: list[dict], fits: list[dict], figdir: Path) -> None:
    # BBox-first geometry: double-column width, two aligned quantitative panels.
    fig = plt.figure(figsize=(7.0, 2.65))
    ax_a = fig.add_axes([0.10, 0.22, 0.38, 0.67])
    ax_b = fig.add_axes([0.59, 0.22, 0.34, 0.67])
    for branch in BRANCHES:
        color = COLORS[branch]
        mean_rows = [r for r in ensemble if r["branch"] == branch]
        k = np.array([r["k_inv_A"] for r in mean_rows])
        weight = np.array([r["W_mean_A2_fs2_per_O"] for r in mean_rows])
        sem = np.array([r["W_replica_SEM_A2_fs2_per_O"] for r in mean_rows])
        for rep in (1, 2, 3):
            _, kr, wr = grouped(rows, branch, rep)
            ax_a.plot(kr, wr * 1e5, color=color, alpha=0.18, lw=1.0)
        ax_a.errorbar(k, weight * 1e5, yerr=sem * 1e5, fmt="o", ms=3.1, capsize=1.8,
                      color=color, lw=1.0, label=LABELS[branch])
        fit = next(r for r in fits if r["branch"] == branch and r["replicate"] == "mean_across_3")
        pred = fit["linear_intercept_A2_fs2_per_O"] + fit["linear_slope_A3_fs2_per_O"] * k
        ax_a.plot(k, pred * 1e5, color=color, lw=1.2, ls="--")
        ax_b.errorbar(k, weight / weight.mean(), yerr=sem / weight.mean(), fmt="o-", ms=3.0,
                      capsize=1.6, color=color, lw=1.1, label=LABELS[branch])
    ax_a.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel=r"$W_a(k_n)$ ($10^{-5}$ Å$^2$ fs$^{-2}$ O$^{-1}$)")
    ax_b.axhline(1.0, color="#777777", lw=1.0, zorder=0)
    ax_b.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel=r"$W_a(k_n)/\langle W_a\rangle_{n=1:20}$")
    ax_a.legend(loc="upper left", fontsize=6.2, ncol=1, handlelength=1.4)
    ax_b.legend(loc="upper left", fontsize=6.2, ncol=1, handlelength=1.4)
    for ax, letter in ((ax_a, "(a)"), (ax_b, "(b)")):
        ax.text(-0.21, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
        ax.tick_params(length=3, labelsize=6.5)
    fig.text(0.30, 0.95, "Thin lines: velocity-seed replicas; dashed: mean linear fit", ha="center", fontsize=6.5)
    fig.text(0.765, 0.95, r"$n=1\ldots20$; error bars: replica SEM", ha="center", fontsize=6.5)
    save(fig, figdir / "modal_spectral_weight_Wk_n1_n20")


def main() -> None:
    if not RAW_CJJ.is_file():
        raise FileNotFoundError(RAW_CJJ)
    data_dir = OUT / "derived_data"
    fig_dir = OUT / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(exist_ok=True)
    raw = read_lag_zero_rows()
    ensemble = ensemble_rows(raw)
    fits = fit_rows(raw)
    write_csv(data_dir / "Wk_n1_n20_per_replica.csv", raw)
    write_csv(data_dir / "Wk_n1_n20_ensemble_mean_sem.csv", ensemble)
    write_csv(data_dir / "Wk_n1_n20_linear_and_powerlaw_fits.csv", fits)
    make_figure(raw, ensemble, fits, fig_dir)
    (OUT / "figure_contract.txt").write_text("""Core conclusion: test whether modal spectral weight W_a(k)=C_JJ,a(k,0)/N follows a reproducible linear trend over n=1..20.
Figure archetype: quantitative grid.
Backend: Python/matplotlib.
Panel map: (a) absolute W_a(k), individual velocity seeds and mean +/- replica SEM with mean linear fits; (b) W_a normalized by its n=1..20 mean to expose relative k-dependence.
Evidence hierarchy: lag-zero per-replica CJJ source values; replica SEM; descriptive linear and log-log fits.
Reviewer risk: three velocity-seed replicas and one finite box cannot establish a k->0 or thermodynamic-limit law; TA_r and TA_theta remain distinct.
""", encoding="utf-8")
    metadata = {
        "system": "(8,8) CNT water, 5L", "source": str(RAW_CJJ), "n_range": [1, NMAX],
        "n_oxygen": N_OXYGEN,
        "W_definition": "W_a(k_n)=C_JJ,a(k_n,0)/N_O; by Wiener-Khinchin it is proportional to the frequency-integrated modal-current spectrum in a fixed normalization",
        "branches": list(BRANCHES), "replicates": 3,
        "fit_boundary": "Linear and log-log fits are descriptive finite-box, finite-mode diagnostics. The primary uncertainty is the SEM across three velocity-seed replicas; within-replica OLS p-values do not establish an asymptotic law.",
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (OUT / "QA_notes.txt").write_text("""Data integrity: selected all and only lag=0, branch={LA,TA_r,TA_theta}, n=1..20 rows (180 rows = 3 branches x 3 replicas x 20 modes); no values imputed or excluded.
Statistics: center is arithmetic mean across three velocity-seed replicas; spread is replica SEM. Linear and log-log fits are run per replica and summarized across replicas.
Image integrity: vector-derived line/scatter figure only; no raster image manipulation.
Interpretation: W is a modal-current spectral-weight proxy, not a global vibrational DOS or a thermodynamic-limit extrapolation.
""", encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("W(k) n=1..20 analysis completed.\n", encoding="utf-8")
    print(json.dumps({"outdir": str(OUT), "fits": fits}, indent=2))


if __name__ == "__main__":
    main()
