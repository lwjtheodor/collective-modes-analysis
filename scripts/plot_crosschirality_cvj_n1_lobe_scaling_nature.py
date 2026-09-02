"""Nature-style cross-chirality scaling of the normalized n=1 CvJ first lobe.

The script copies the file-backed per-replica source table locally, then uses
only local CSV assets to calculate mean/SEM and log-log least-squares fits.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "assets"
REMOTE = Path(r"H:\gcmc_explore\translational_anomaly\08_viscosity_friction_length_scaling\04_analysis\frequency_mode_response\cross_chirality_current_20260805\summary")
SOURCE = ASSET / "crosschirality_CvJ_n1_first_lobe_paired.csv"
SUMMARY = ASSET / "crosschirality_CvJ_n1_first_lobe_summary.csv"
FITS = ASSET / "crosschirality_CvJ_n1_first_lobe_k_scaling_fits.csv"
OUT = ASSET / "crosschirality_CvJ_n1_first_lobe_k_scaling_nature"

COLORS = {"7_7": "#3B7EA1", "8_8": "#D88737", "9_9": "#5E9C76", "17_0": "#8C6BB1"}
LABELS = {"7_7": "(7,7)", "8_8": "(8,8)", "9_9": "(9,9)", "17_0": "(17,0)"}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 150,
})


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def mean_sem(values):
    a = np.asarray(values, float)
    return float(a.mean()), (float(a.std(ddof=1) / math.sqrt(a.size)) if a.size > 1 else float("nan"))


def main():
    ASSET.mkdir(exist_ok=True)
    # Immutable local source-data copy; all later computations read this copy.
    original = read_csv(REMOTE / "first_negative_lobe_by_replicate.csv")
    write_csv(SOURCE, original, list(original[0]))
    rows = read_csv(SOURCE)
    for r in rows:
        for key in ("L", "replicate", "depth", "negative_area_ps", "t_min_ps", "t_start_ps", "t_end_ps"):
            r[key] = float(r[key])

    grouped = {}
    for r in rows: grouped.setdefault((r["chirality"], int(r["L"])), []).append(r)
    summary = []
    for (chi, L), group in sorted(grouped.items()):
        vals = {key: mean_sem([g[key] for g in group]) for key in ("negative_area_ps", "depth", "t_min_ps", "t_start_ps", "t_end_ps")}
        # n=1 fundamental k from audited fit summary; 17,0 has a distinct axial repeat length.
        fit_rows = read_csv(REMOTE / "cross_length_CJ_n1_fit_summary.csv")
        k = next(float(x["k_inv_A"]) for x in fit_rows if x["chirality"] == chi and int(x["length_L"]) == L)
        summary.append({"chirality": chi, "L": L, "n_replicates": len(group), "k_inv_A": k,
                        **{f"{key}_{suf}": val for key, pair in vals.items() for suf, val in zip(("mean", "sem"), pair)}})
    fields = list(summary[0]); write_csv(SUMMARY, summary, fields)

    fits = []
    for chi in COLORS:
        subset = [s for s in summary if s["chirality"] == chi]
        assert all(s["k_inv_A"] > 0 and s["negative_area_ps_mean"] > 0 for s in subset)
        x = np.log([s["k_inv_A"] for s in subset]); y = np.log([s["negative_area_ps_mean"] for s in subset])
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (intercept + slope*x)
        slope_se = math.sqrt(float((resid**2).sum()) / (len(x)-2) / float(((x-x.mean())**2).sum()))
        r2 = 1 - float((resid**2).sum()) / float(((y-y.mean())**2).sum())
        fits.append({"chirality": chi, "n_lengths": len(subset), "n_replicate_records": sum(s["n_replicates"] for s in subset),
                     "area_prefactor_ps": math.exp(intercept), "area_vs_k_exponent": slope,
                     "area_vs_k_exponent_sem": slope_se, "area_vs_L_exponent": -slope, "r_squared_loglog": r2,
                     "caveat": "9_9 L3 has n=1; unweighted fit across length means"})
    write_csv(FITS, fits, list(fits[0]))

    fig, ax = plt.subplots(1, 2, figsize=(7.15, 2.85), gridspec_kw={"width_ratios": [1.45, 1]})
    for chi in COLORS:
        s = [r for r in summary if r["chirality"] == chi]
        k = np.array([r["k_inv_A"] for r in s]); area = np.array([r["negative_area_ps_mean"] for r in s]); sem = np.array([r["negative_area_ps_sem"] for r in s])
        ax[0].errorbar(k, area, yerr=sem, marker="o", ms=4, lw=1.4, capsize=2, color=COLORS[chi], label=LABELS[chi])
        fit = next(r for r in fits if r["chirality"] == chi)
        kg = np.geomspace(k.min(), k.max(), 120)
        ax[0].plot(kg, fit["area_prefactor_ps"] * kg**fit["area_vs_k_exponent"], color=COLORS[chi], ls="--", lw=0.9, alpha=.8)
    ax[0].set(xscale="log", yscale="log", xlabel=r"fundamental wavevector $k_1$ ($\mathrm{\AA}^{-1}$)", ylabel=r"first negative-lobe area $A_- / C_{vJ}(0)$ (ps)")
    ax[0].legend(title="chirality", fontsize=6, title_fontsize=6, loc="upper left", handlelength=1.5)
    ax[0].text(.985, .035, r"$A_-\propto k_1^{-p}$; dashed = log--log OLS", transform=ax[0].transAxes, ha="right", va="bottom", fontsize=6)
    ax[0].text(-.17, 1.04, "a", transform=ax[0].transAxes, fontweight="bold", fontsize=9)

    order = list(COLORS)
    p = [next(r for r in fits if r["chirality"] == chi)["area_vs_L_exponent"] for chi in order]
    pe = [next(r for r in fits if r["chirality"] == chi)["area_vs_k_exponent_sem"] for chi in order]
    ax[1].axhline(1, color="#8A8A8A", lw=.7, ls=":", zorder=0)
    ax[1].errorbar(np.arange(4), p, yerr=pe, fmt="none", ecolor="#555555", capsize=2, lw=.8, zorder=2)
    ax[1].scatter(np.arange(4), p, s=30, c=[COLORS[c] for c in order], edgecolors="white", linewidth=.5, zorder=3)
    ax[1].set(xticks=np.arange(4), xticklabels=[LABELS[c] for c in order], ylabel=r"length exponent $p$ in $A_-\propto L^p$", ylim=(0.9, 2.05))
    for i, (v, e) in enumerate(zip(p, pe)): ax[1].text(i, v+.09, f"{v:.2f}±{e:.02f}", ha="center", fontsize=5.8)
    ax[1].text(-.18, 1.04, "b", transform=ax[1].transAxes, fontweight="bold", fontsize=9)
    fig.text(.505, .995, "Weak Nosé–Hoover; no momentum removal; 10 fs output; 1 ns; normalized total $C_{vJ}$; replicas: 7,7=2; 8,8=3; 9,9=2 (L3=1); 17,0=2; error bars = replica SEM", ha="center", va="top", fontsize=5.8)
    fig.subplots_adjust(left=.09, right=.99, bottom=.2, top=.85, wspace=.36)
    # Literal export calls keep the deterministic preflight auditable.
    fig.savefig(str(OUT) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(OUT) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(OUT) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT) + ".svg", bbox_inches="tight")
    print("wrote", SOURCE, SUMMARY, FITS, OUT)


if __name__ == "__main__": main()
