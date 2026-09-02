"""Test direct wavelength collapse of the C(0)-normalized first-negative-lobe area."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "assets"
INPUT = ASSET / "cjj_88_n1_baseline_and_10L_all_modes_vs_invk.csv"
CSV = ASSET / "cjj_88_10L_normalized_lobe_area_over_lambda.csv"
STEM = ASSET / "cjj_88_10L_normalized_lobe_area_over_lambda_nature"


def main():
    all_data = pd.read_csv(INPUT)
    q = all_data[all_data.source.str.startswith("10L")].sort_values("wavelength_nm").copy()
    q["area_over_lambda_ps_per_nm"] = q["area_mean_ps"] / q["wavelength_nm"]
    q["area_over_lambda_sem_ps_per_nm"] = q["area_sem_ps"] / q["wavelength_nm"]
    q["duration_over_lambda_ps_per_nm"] = q["duration_mean_ps"] / q["wavelength_nm"]
    q["duration_over_lambda_sem_ps_per_nm"] = q["duration_sem_ps"] / q["wavelength_nm"]
    q["geometric_shape_S"] = q["area_mean_ps"] / (q["duration_mean_ps"] * q["depth_mean"])
    q["identity_area_over_lambda"] = q["duration_over_lambda_ps_per_nm"] * q["depth_mean"] * q["geometric_shape_S"]
    q["definition"] = "A_minus/lambda; A_minus=-integral_lobe CJJ(t)/CJJ(0) dt; lambda=2pi/k"
    q.to_csv(CSV, index=False, float_format="%.9g")

    plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.linewidth": 0.8,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "xtick.direction": "out", "ytick.direction": "out"})
    fig, axs = plt.subplots(1, 3, figsize=(7.25, 2.75), constrained_layout=True)
    blue, red, gray = "#2775a9", "#c94c4c", "#555555"

    ax = axs[0]
    ax.errorbar(q.wavelength_nm, q.area_over_lambda_ps_per_nm, q.area_over_lambda_sem_ps_per_nm,
                color=blue, marker="o", ms=4.8, lw=1.15, capsize=2)
    ax.set(xlabel=r"wavelength, $\lambda=2\pi/k$ (nm)", ylabel=r"$A_-/\lambda$ (ps nm$^{-1}$)", xlim=(0, 107))
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(.05, .92, f"n=10: {q.area_over_lambda_ps_per_nm.iloc[0]:.4f}\n"
                        f"n=1: {q.area_over_lambda_ps_per_nm.iloc[-1]:.4f}",
            transform=ax.transAxes, va="top", fontsize=6.8)
    ax.text(-.14, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=10)

    ax = axs[1]
    ax.errorbar(q.wavelength_nm, q.duration_over_lambda_ps_per_nm, q.duration_over_lambda_sem_ps_per_nm,
                color=blue, marker="o", ms=4.8, lw=1.15, capsize=2)
    ax.axhline(q.duration_over_lambda_ps_per_nm.mean(), color=gray, ls="--", lw=.9,
               label=f"mean={q.duration_over_lambda_ps_per_nm.mean():.3f}")
    ax.set(xlabel=r"wavelength, $\lambda=2\pi/k$ (nm)", ylabel=r"$\tau_-/\lambda$ (ps nm$^{-1}$)", xlim=(0, 107))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=6.3, loc="lower right")
    ax.text(-.14, 1.04, "(b)", transform=ax.transAxes, fontweight="bold", fontsize=10)

    ax = axs[2]
    ax.errorbar(q.wavelength_nm, q.depth_mean, q.depth_sem, color=red, marker="s", ms=4.5, lw=1.15,
                capsize=2, label=r"$D_-$")
    ax.plot(q.wavelength_nm, q.geometric_shape_S, color=blue, marker="o", ms=4, lw=1.0,
            label=r"$S=A_-/(\tau_-D_-)$")
    ax.set(xlabel=r"wavelength, $\lambda=2\pi/k$ (nm)", ylabel="depth / geometric shape factor", xlim=(0, 107))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=6.2, loc="lower right")
    ax.text(-.14, 1.04, "(c)", transform=ax.transAxes, fontweight="bold", fontsize=10)

    fig.suptitle(r"$(8,8)$, 10L all modes: direct wavelength normalization of the first negative-lobe area", y=1.03, fontsize=10)
    fig.text(.5, -.075, r"The direct test is $A_-/\lambda$. It is not invariant because $A_-/\lambda=(\tau_-/\lambda)D_-S$: duration/ wavelength is nearly constant, whereas depth grows and the geometric shape changes. Weak NH/no global momentum removal; 100 fs/10 ns; 3 replicas; error bars = SEM.", ha="center", fontsize=5.8, color=".25")
    fig.savefig(STEM.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(STEM.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
