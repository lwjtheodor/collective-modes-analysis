"""Assemble protocol-stratified matched-k effective-DHO damping from archived CJJ fits.

This intentionally consumes only completed time-domain CJJ fit tables.  Welch spectra
without a validated linewidth fit remain catalogued as available-but-not-fitted.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 1.0,
})

ROOT = Path(r"H:/gcmc_explore")
OUT = Path(r"H:/gcmc_explore/translational_anomaly/02_isf_collective_modes/results/collective_mode_response/implicit_C77_C88_C99_matched_k_effective_DHO_damping/2026-08-29")

SOURCES = {
    "C77": ROOT / "implicit_chirality_length_scan_20260816/c77_T330_weakNH_r4/cjj_time_domain_lowk_20260817",
    "C88": ROOT / "implicit_cnt_8_8_length_scan_20260814/cjj_time_domain_8_8_allN_20260817",
    "C99": ROOT / "implicit_chirality_length_scan_20260816/c99_T350_weakNH_r4/cjj_time_domain_lowk_20260817",
}
SPECTRAL_LONG = ROOT / "implicit_chirality_length_scan_20260816/c99_T350_weakNH_r4_N2400_N3200_extension_20260817"

# Direct-NVT 10 fs / 1.2 ns pilot.  It is useful as an internal sensitivity
# check, but is not protocol-matched to the 6 ns weak-NH production series.
EXCLUDED_CASES = {("C88", "N800_pilot_repair")}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        writer.writeheader()
        writer.writerows(rows)


def read_time_fits() -> list[dict]:
    rows: list[dict] = []
    for system, base in SOURCES.items():
        for table in sorted(base.glob("N*/time_domain_omega_fits.csv")):
            if (system, table.parent.name) in EXCLUDED_CASES:
                continue
            with table.open(encoding="utf-8") as handle:
                for r in csv.DictReader(handle):
                    rows.append({
                        "system": system,
                        "case": table.parent.name,
                        "replica": int(r["replica"]),
                        "mode_n": int(r["mode_n"]),
                        "k_Ainv": float(r["k_inv_A"]),
                        "omega_rad_ps": float(r["omega_fit_rad_ps"]),
                        "gamma_psinv": float(r["gamma_ps_inv"]),
                        "fit_chi2": float(r["chi2"]),
                        "dof": int(r["dof"]),
                        "source_table": str(table),
                    })
    return rows


def aggregate(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        grouped[(r["system"], r["case"], r["mode_n"], round(r["k_Ainv"], 10))].append(r)
    out: list[dict] = []
    for (system, case, n, k), values in sorted(grouped.items()):
        gamma = np.array([v["gamma_psinv"] for v in values])
        omega = np.array([v["omega_rad_ps"] for v in values])
        out.append({
            "system": system, "case": case, "mode_n": n, "k_Ainv": k,
            "n_velocity_seeds": len(values),
            "gamma_effective_psinv_mean": gamma.mean(),
            "gamma_velocity_seed_sem_psinv": gamma.std(ddof=1) / np.sqrt(len(gamma)) if len(gamma) > 1 else np.nan,
            "omega_mean_rad_ps": omega.mean(),
            "omega_velocity_seed_sem_rad_ps": omega.std(ddof=1) / np.sqrt(len(omega)) if len(omega) > 1 else np.nan,
            "mean_reduced_fit_chi2": np.mean([v["fit_chi2"] / v["dof"] for v in values]),
            "definition": "effective damped-cosine fit to normalized CJJ(k,t); not independently validated spectral HWHM",
        })
    return out


def fit_k2(points: list[dict], system: str) -> dict:
    # Fit only the common low-k domain; SEM comes from velocity seeds, not independent configurations.
    p = [r for r in points if r["system"] == system and r["k_Ainv"] <= 0.0352 and np.isfinite(r["gamma_velocity_seed_sem_psinv"]) and r["gamma_velocity_seed_sem_psinv"] > 0]
    x = np.array([r["k_Ainv"] ** 2 for r in p]); y = np.array([r["gamma_effective_psinv_mean"] for r in p]); s = np.array([r["gamma_velocity_seed_sem_psinv"] for r in p])
    # Collapse duplicated physical k across boxes first, to avoid giving one k excess leverage.
    collapsed = []
    for k in sorted(set(np.round(np.sqrt(x), 8))):
        mask = np.isclose(np.sqrt(x), k, atol=2e-8)
        w = 1 / s[mask] ** 2
        collapsed.append((k * k, np.sum(w * y[mask]) / np.sum(w), np.sqrt(1 / np.sum(w))))
    x = np.array([v[0] for v in collapsed]); y = np.array([v[1] for v in collapsed]); s = np.array([v[2] for v in collapsed])
    if len(x) < 3:
        return {"system": system, "model": "Gamma0_plus_Ak2", "status": "insufficient_distinct_k", "n_distinct_k": len(x)}
    X = np.column_stack([np.ones(len(x)), x]); W = np.diag(1 / s**2)
    cov = np.linalg.inv(X.T @ W @ X); beta = cov @ X.T @ W @ y
    resid = y - X @ beta; chi2 = float(np.sum((resid / s) ** 2))
    reduced = chi2 / max(1, len(x) - 2)
    status = "screening_only_poor_k2_adequacy" if reduced > 3 else "screening_only_k2_adequacy_not_rejected"
    return {"system": system, "model": "Gamma0_plus_Ak2", "status": status,
            "n_distinct_k": len(x), "Gamma0_psinv": beta[0], "Gamma0_seedSEM_psinv": np.sqrt(cov[0, 0]),
            "A_psinv_A2": beta[1], "A_seedSEM_psinv_A2": np.sqrt(cov[1, 1]),
            "chi2": chi2, "reduced_chi2": reduced,
            "scope": "screening reference only; velocity-seed uncertainty; common-parent configurations; no zero-k or thermodynamic-limit claim"}


def matched_k_audit(points: list[dict]) -> list[dict]:
    """Cross-box agreement at the same physical k; seed SEM is not independent-start uncertainty."""
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for r in points:
        grouped[(r["system"], round(r["k_Ainv"], 8))].append(r)
    out = []
    for (system, k), values in sorted(grouped.items()):
        if len(values) < 2:
            continue
        y = np.array([v["gamma_effective_psinv_mean"] for v in values])
        s = np.array([v["gamma_velocity_seed_sem_psinv"] for v in values])
        finite = np.isfinite(s) & (s > 0)
        if finite.all():
            w = 1 / s**2
            mean = float(np.sum(w * y) / np.sum(w))
            chi2 = float(np.sum(((y - mean) / s)**2))
            red = chi2 / max(1, len(y) - 1)
        else:
            mean, chi2, red = float(np.mean(y)), np.nan, np.nan
        out.append({
            "system": system, "k_Ainv": k, "n_boxes": len(values),
            "cases": ";".join(v["case"] for v in values),
            "gamma_min_psinv": float(np.min(y)), "gamma_max_psinv": float(np.max(y)),
            "inverse_variance_mean_psinv": mean, "cross_box_chi2": chi2,
            "cross_box_reduced_chi2": red,
            "interpretation": "velocity-seed SEM only; common-parent configuration means this is an operational matched-k check, not an independent-replica consistency test",
        })
    return out


def spectral_inventory() -> list[dict]:
    rows = []
    for case in ("N2400", "N3200"):
        d = SPECTRAL_LONG / case / "analysis_msd_cvv_modes_200ps_20260818/collective_modes"
        rows.append({"system": "C99", "case": case, "asset": "CJJ_kw_welch spectra + peak tables", "n_replicas": 4,
                     "status": "available_for_linewidth_refit_not_yet_DHO_width_fitted", "path": str(d)})
    return rows


def plot(points: list[dict], fits: list[dict]) -> None:
    colors = {"C77": "#0072B2", "C88": "#D55E00", "C99": "#009E73"}
    fig = plt.figure(figsize=(5.5, 2.55))
    ax = fig.add_axes([0.13, 0.22, 0.84, 0.70])
    for system in ("C77", "C88", "C99"):
        p = [r for r in points if r["system"] == system and r["k_Ainv"] <= 0.0352]
        x = np.array([r["k_Ainv"]**2 for r in p]); y = np.array([r["gamma_effective_psinv_mean"] for r in p]); e = np.array([r["gamma_velocity_seed_sem_psinv"] for r in p])
        ax.errorbar(x, y, yerr=e, fmt="o", ms=3.2, capsize=2, color=colors[system], label=system, lw=1.0)
        f = next((z for z in fits if z["system"] == system and "Gamma0_psinv" in z), None)
        if f:
            xx = np.linspace(0, max(x) * 1.04, 100)
            style = "--" if f["status"] == "screening_only_poor_k2_adequacy" else "-"
            ax.plot(xx, f["Gamma0_psinv"] + f["A_psinv_A2"] * xx, color=colors[system], lw=1.1, ls=style)
    ax.axhline(0, color="0.25", lw=1.0)
    ax.set_xlabel(r"$k^2$ ($\mathrm{\AA}^{-2}$)", fontname="Arial", fontsize=7)
    ax.set_ylabel(r"effective $\Gamma$ (ps$^{-1}$)", fontname="Arial", fontsize=7)
    ax.text(-0.13, 1.05, "(a)", transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left")
    ax.tick_params(direction="out", labelsize=7, width=1.0, length=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7, ncol=3, handletextpad=0.3, columnspacing=0.8)
    fig.savefig(OUT / "figures/matched_k_effective_DHO_gamma_k2.png", dpi=600)
    fig.savefig(OUT / "figures/matched_k_effective_DHO_gamma_k2.pdf")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True); (OUT / "figures").mkdir(exist_ok=True); (OUT / "derived_data").mkdir(exist_ok=True)
    raw = read_time_fits(); points = aggregate(raw); fits = [fit_k2(points, s) for s in ("C77", "C88", "C99")]
    write_csv(OUT / "derived_data/per_replica_effective_DHO_fits.csv", raw)
    write_csv(OUT / "derived_data/matched_k_effective_DHO_summary.csv", points)
    write_csv(OUT / "derived_data/matched_k_cross_box_consistency.csv", matched_k_audit(points))
    write_csv(OUT / "derived_data/gamma0_plus_k2_operational_model_fits.csv", fits)
    write_csv(OUT / "derived_data/C99_longbox_spectral_inventory.csv", spectral_inventory())
    plot(points, fits)
    (OUT / "README.md").write_text(
        "# Protocol-stratified matched-k effective-DHO damping\n\n"
        "This package re-aggregates only completed normalized longitudinal CJJ(k,t) fit tables. "
        "`gamma_effective_psinv` is the damped-cosine decay width already fitted from the time domain; it is **not** an independently verified spectral DHO HWHM.\n\n"
        "- C77: implicit (7,7), 330 K, weak-NH production; N200--N1600 time-domain CJJ fits.\n"
        "- C88: implicit (8,8), 350 K, weak-NH production; N200/N400/N1600 time-domain CJJ fits. The direct-NVT N800 pilot is deliberately excluded.\n"
        "- C99: implicit (9,9), 350 K, weak-NH production; N400/N800/N1600 time-domain CJJ fits. N2400/N3200 Welch spectra are inventoried but do not yet have a validated linewidth/DHO refit.\n\n"
        "The `Gamma0_plus_Ak2` table is a finite-k screening regression. Its intercept must not be read as physical zero-wave-number damping: replicas are velocity seeds from common parent configurations, and poor reduced chi-square explicitly rejects the simple k2 model.\n",
        encoding="utf-8")
    (OUT / "QA.md").write_text(
        "# QA and interpretation boundary\n\n"
        "Completed: source time-fit tables were parsed, per-seed values retained, equal-k box results inverse-variance collapsed, and a k2 screening regression performed below 0.0352 A^-1.\n\n"
        "Not completed: a shared raw-spectrum DHO/HWHM refit for C99 N2400/N3200; independent configurational replicas; a common-temperature/common-chirality universality claim; or a physical zero-k damping determination.\n",
        encoding="utf-8")
    (OUT / "metadata.json").write_text(json.dumps({
        "created": "2026-08-29",
        "observable": "normalized longitudinal CJJ(k,t) effective damped-cosine width",
        "systems": {
            "C77": {"chirality": "(7,7)", "temperature_K": 330, "production": "weak-NH, 6 ns, 100 fs, 4 velocity seeds", "boxes": ["N200", "N400", "N800", "N1600"]},
            "C88": {"chirality": "(8,8)", "temperature_K": 350, "production": "weak-NH, 6 ns, 100 fs, 4 velocity seeds", "boxes": ["N200", "N400", "N1600"], "excluded": "N800_pilot_repair direct-NVT 10 fs / 1.2 ns"},
            "C99": {"chirality": "(9,9)", "temperature_K": 350, "production": "weak-NH, 6 ns, 100 fs, 4 velocity seeds", "time_domain_boxes": ["N400", "N800", "N1600"], "spectra_only_boxes": ["N2400", "N3200"]},
        },
        "k2_regression": "finite-k screening only; no physical Gamma(k=0) determination",
        "primary_limit": "replicas are velocity seeds from common parent configurations",
    }, indent=2), encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("Completed protocol-stratified matched-k effective-DHO damping assembly.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
