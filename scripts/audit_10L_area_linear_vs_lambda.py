"""Source-level audit: raw 10L first-lobe area is linear in wavelength, not through zero."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "collective_mode_response" / "fig2_longitudinal_modes_88_rh75_330k" / "2026-08-11" / "remote_raw" / "output"
OUT = ROOT / "assets"
CSV = OUT / "cjj_88_10L_raw_source_area_lambda_linearity_audit.csv"
STEM = OUT / "cjj_88_10L_raw_source_area_lambda_linearity_audit_nature"


def sem(values):
    return np.std(values, ddof=1) / np.sqrt(len(values))


def main():
    rows = []
    for f in sorted(RAW.glob("8_8_L10_rep*_CJJ_alln.json")):
        meta = json.loads(f.read_text())
        rep = int(meta["case_id"].split("_")[-1][3:])
        for z in meta["mode_summary"]:
            lobe = z["first_negative_lobe_normalized"]
            rows.append({"replicate": rep, "n": z["n"], "k_inv_A": z["k_inv_A"],
                         "wavelength_nm": 2 * np.pi * 0.1 / z["k_inv_A"],
                         "area_ps": lobe["negative_area_normalized_ps"]})
    raw = pd.DataFrame(rows)
    d = raw.groupby(["n", "k_inv_A", "wavelength_nm"], as_index=False).agg(
        area_mean_ps=("area_ps", "mean"), area_sem_ps=("area_ps", sem))
    d = d.sort_values("wavelength_nm")
    X = np.c_[np.ones(len(d)), d["wavelength_nm"]]
    intercept, slope = np.linalg.lstsq(X, d["area_mean_ps"], rcond=None)[0]
    d["linear_prediction_ps"] = intercept + slope * d["wavelength_nm"]
    d["linear_residual_ps"] = d["area_mean_ps"] - d["linear_prediction_ps"]
    d["area_over_lambda_ps_per_nm"] = d["area_mean_ps"] / d["wavelength_nm"]
    d["area_over_lambda_sem_ps_per_nm"] = d["area_sem_ps"] / d["wavelength_nm"]
    r2 = 1 - np.sum(d["linear_residual_ps"] ** 2) / np.sum((d["area_mean_ps"] - d["area_mean_ps"].mean()) ** 2)
    d["linear_slope_ps_per_nm"] = slope
    d["linear_intercept_ps"] = intercept
    d["linear_R2"] = r2
    d["source"] = "direct mode_summary entries in 10L CJJ_alln JSON, 3 replicas"
    d.to_csv(CSV, index=False, float_format="%.9g")

    plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.linewidth": .8,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "xtick.direction": "out", "ytick.direction": "out"})
    fig, axs = plt.subplots(1, 2, figsize=(7.25, 2.75), constrained_layout=True)
    blue, red, gray = "#2775a9", "#c94c4c", "#555555"
    xx = np.linspace(0, 107, 400)

    ax = axs[0]
    ax.errorbar(d.wavelength_nm, d.area_mean_ps, d.area_sem_ps, color=blue, marker="o", ms=5,
                lw=0, capsize=2, label="10L raw-source modes")
    ax.plot(xx, intercept + slope * xx, color=red, lw=1.5,
            label=rf"linear: $A_-={intercept:.3f}+{slope:.3f}\lambda$; $R^2={r2:.5f}$")
    ax.set(xlabel=r"wavelength, $\lambda=2\pi/k$ (nm)", ylabel=r"normalized lobe area, $A_-$ (ps)", xlim=(0, 107))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=6.2, loc="upper left")
    ax.text(-.12, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=10)

    ax = axs[1]
    ax.errorbar(d.wavelength_nm, d.area_over_lambda_ps_per_nm, d.area_over_lambda_sem_ps_per_nm,
                color=blue, marker="o", ms=5, lw=1.1, capsize=2)
    xfit = np.linspace(d.wavelength_nm.min(), d.wavelength_nm.max(), 400)
    ax.plot(xfit, intercept / xfit + slope, color=gray, lw=1.2, ls="--",
            label=r"required by nonzero intercept: $A_-/\lambda=b+a/\lambda$")
    ax.set(xlabel=r"wavelength, $\lambda=2\pi/k$ (nm)", ylabel=r"$A_-/\lambda$ (ps nm$^{-1}$)", xlim=(0, 107), ylim=(.05, .15))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=6.3, loc="lower right")
    ax.text(-.12, 1.04, "(b)", transform=ax.transAxes, fontweight="bold", fontsize=10)

    fig.suptitle(r"10L source audit: raw lobe area is highly linear in $\lambda$, but has a nonzero intercept", y=1.03, fontsize=10)
    fig.text(.5, -.075, r"No source change: both panels use the same 30 direct JSON values (10 modes x 3 replicas). Dividing a fitted relation $A_-=a+b\lambda$ by $\lambda$ necessarily gives $A_-/\lambda=b+a/\lambda$, which cannot be constant unless $a=0$.", ha="center", fontsize=6.0, color=".25")
    fig.savefig(STEM.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
