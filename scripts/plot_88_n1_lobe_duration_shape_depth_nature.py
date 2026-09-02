"""Nature-style (8,8) n=1 current-mode first-negative-lobe summary.

The three reported observables are deliberately distinct:
  tau_minus  = t_plus - t_minus, the zero-crossing-bounded negative duration;
  Phi_return = [-integral C_JJ(t)/C_JJ(0) dt] / (Lz/c_s), the return-scaled
                negative-lobe area (a dimensionless shape factor);
  D_minus    = -min[C_JJ(t)/C_JJ(0)], the relative valley depth.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "collective_mode_response" / "fig2_longitudinal_modes_88_rh75_330k" / "2026-08-11" / "derived_data" / "panel_b_lowk_strength.csv"
OUT = ROOT / "assets"
OUT.mkdir(exist_ok=True)
CSV_OUT = OUT / "cjj_88_n1_lobe_duration_shape_depth_vs_L.csv"
STEM = OUT / "cjj_88_n1_lobe_duration_shape_depth_vs_L_nature"


def main():
    d = pd.read_csv(SRC).sort_values("Lz_nm").copy()
    d["negative_duration_over_return"] = d["lobe_width_ps_mean"] / d["t_return_ps"]
    d["observable_definition"] = (
        "first complete negative lobe, bounded by the first two zero crossings"
    )
    d["normalization_definition"] = (
        "CJJ/CJJ(0); Phi_return=(-integral_lobe CJJ/CJJ(0)dt)/(Lz/c_s)"
    )
    d["chirality"] = "(8,8)"
    d["mode"] = "axial current n=1 (k1=2pi/Lz)"
    d["replicas"] = 3
    d["uncertainty"] = "replica SEM"
    d["protocol"] = (
        "weak Nose-Hoover; no global momentum removal; instantaneous water-COM axial subtraction"
    )
    d["sampling_note"] = np.where(
        d["L"] <= 5,
        "10 fs cadence; 1 ns duration",
        "100 fs cadence; 10 ns duration",
    )
    cols = [
        "chirality", "L", "Lz_nm", "k_inv_A", "mode", "negative_duration_ps_mean",
        "negative_duration_ps_sem", "negative_duration_over_return", "t_return_ps",
        "negative_area_C0_ps_mean", "negative_area_C0_ps_sem", "phi_return_mean",
        "phi_return_sem", "negative_depth_relative_mean", "negative_depth_relative_sem",
        "area_over_duration_mean", "area_over_duration_sem", "replicas", "uncertainty",
        "observable_definition", "normalization_definition", "protocol", "sampling_note",
    ]
    archive = pd.DataFrame({
        "chirality": d["chirality"], "L": d["L"], "Lz_nm": d["Lz_nm"],
        "k_inv_A": d["k_inv_A"], "mode": d["mode"],
        "negative_duration_ps_mean": d["lobe_width_ps_mean"],
        "negative_duration_ps_sem": d["lobe_width_ps_sem"],
        "negative_duration_over_return": d["negative_duration_over_return"],
        "t_return_ps": d["t_return_ps"],
        "negative_area_C0_ps_mean": d["A_minus_ps_mean"],
        "negative_area_C0_ps_sem": d["A_minus_ps_sem"],
        "phi_return_mean": d["A_minus_norm_mean"],
        "phi_return_sem": d["A_minus_norm_sem"],
        "negative_depth_relative_mean": d["depth_norm_mean"],
        "negative_depth_relative_sem": d["depth_norm_sem"],
        "area_over_duration_mean": d["A_minus_shape_mean"],
        "area_over_duration_sem": d["A_minus_shape_sem"],
        "replicas": d["replicas"], "uncertainty": d["uncertainty"],
        "observable_definition": d["observable_definition"],
        "normalization_definition": d["normalization_definition"],
        "protocol": d["protocol"], "sampling_note": d["sampling_note"],
    })[cols]
    archive.to_csv(CSV_OUT, index=False, float_format="%.9g")

    plt.rcParams.update({
        "font.family": "Arial", "font.size": 8.5, "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.direction": "out", "ytick.direction": "out",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.65), constrained_layout=True)
    x = d["Lz_nm"].to_numpy()
    normal = "#2a6f9e"
    short = "#c94c4c"
    metrics = [
        ("lobe_width_ps_mean", "lobe_width_ps_sem", r"negative duration, $\tau_-$ (ps)", "a"),
        ("A_minus_norm_mean", "A_minus_norm_sem", r"return-scaled area, $\Phi_- = A_-/(L_z/c_s)$", "b"),
        ("depth_norm_mean", "depth_norm_sem", r"relative valley depth, $D_- = -\min(C_{JJ}/C_{JJ,0})$", "c"),
    ]
    for ax, (mean, sem, ylabel, tag) in zip(axes, metrics):
        ax.plot(x, d[mean], color=normal, lw=1.35, zorder=1)
        ax.errorbar(x[1:], d[mean].to_numpy()[1:], yerr=d[sem].to_numpy()[1:],
                    fmt="o", ms=4.8, color=normal, mec="white", mew=0.65,
                    capsize=2.4, lw=1.0, zorder=3)
        ax.errorbar(x[:1], d[mean].to_numpy()[:1], yerr=d[sem].to_numpy()[:1],
                    fmt="o", ms=5.5, color=short, mec="white", mew=0.65,
                    capsize=2.4, lw=1.0, zorder=4)
        ax.set_xlabel(r"$L_z$ (nm)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(16, 105)
        ax.set_xticks([20, 40, 60, 80, 100])
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.16, 1.04, f"({tag})", transform=ax.transAxes, fontweight="bold", fontsize=10)

    axes[0].plot(x, 0.5 * d["t_return_ps"], color="0.35", lw=1.0, ls="--", zorder=0,
                 label=r"$t_{\mathrm{return}}/2$")
    axes[0].legend(frameon=False, fontsize=7.3, loc="upper left", handlelength=2.0)
    axes[0].annotate(r"20 nm: shortest box" + "\n" + r"$k_1=0.0312$ $\AA^{-1}$",
                     xy=(x[0], d["lobe_width_ps_mean"].iloc[0]), xycoords="data",
                     xytext=(30, 21), textcoords="data", fontsize=7.1, color=short,
                     arrowprops=dict(arrowstyle="-", color=short, lw=0.75))
    axes[1].axhspan(0.18, 0.205, color=normal, alpha=0.08, lw=0)
    axes[1].text(68, 0.184, "30–50 nm\nnear-plateau", fontsize=6.9, color="0.32")
    axes[2].set_ylim(0.44, 0.69)

    fig.suptitle(r"$(8,8)$ water: first complete negative lobe of the $n=1$ axial-current ACF", y=1.04, fontsize=10)
    fig.text(0.5, -0.085,
             "Weak NH; no global momentum removal; instantaneous water-COM axial subtraction; "
             "3 replicas, error bars = replica SEM.  2–5L: 10 fs/1 ns; 10L: 100 fs/10 ns.",
             ha="center", va="top", fontsize=6.6, color="0.24")
    # Explicit export calls keep the complete Nature Figure deliverable auditable.
    fig.savefig(STEM.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
