"""Audit local matched-k current-mode coverage and visualize LA modal weights."""
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
    parser.add_argument("--five", type=Path, required=True)
    parser.add_argument("--ten", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data_dir, figure_dir = args.output / "derived_data", args.output / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(exist_ok=True)
    five = pd.read_csv(args.five / "derived_data" / "CJJ_all_modes_ensemble_mean_sem.csv")
    ten = pd.read_csv(args.ten / "derived_data" / "CJJ_all_modes_ensemble_mean_sem.csv")
    f0 = five[(five.branch == "LA") & (five.n.between(1, 5)) & (five.lag_ps == 0)].copy()
    z0 = ten[(ten.branch == "LA") & (ten.n.isin([2, 4, 6, 8, 10])) & (ten.time_ps == 0)].copy()
    z0["n_5L_equivalent"] = z0.n // 2
    merged = f0.merge(z0, left_on="n", right_on="n_5L_equivalent", suffixes=("_5L", "_10L"))
    merged["relative_k_mismatch"] = np.abs(merged.k_inv_A_5L - merged.k_inv_A_10L) / merged.k_inv_A_5L
    merged["raw_weight_ratio_10L_over_5L"] = merged.CJJ_mean_A2_fs2 / merged.CJJ_raw_mean_A2_fs2
    merged["weight_5L_shared_k_fraction"] = merged.CJJ_raw_mean_A2_fs2 / merged.CJJ_raw_mean_A2_fs2.sum()
    merged["weight_10L_shared_k_fraction"] = merged.CJJ_mean_A2_fs2 / merged.CJJ_mean_A2_fs2.sum()
    merged.to_csv(data_dir / "LA_matched_k_5L_10L_raw_weight_comparison.csv", index=False)

    # Time-domain consistency diagnostic: shapes normalized separately, thus deliberately not an absolute-weight test.
    rows = []
    for _, item in merged.iterrows():
        n5, n10 = int(item.n_5L), int(item.n_10L)
        a = five[(five.branch == "LA") & (five.n == n5) & (five.lag_ps <= 50)].copy()
        b = ten[(ten.branch == "LA") & (ten.n == n10) & (ten.time_ps <= 50)].copy()
        a["length"] = "5L"; a["time_ps_common"] = a.lag_ps
        a["C_normalized"] = a.CJJ_normalized_mean; a["C_normalized_SEM"] = a.CJJ_normalized_replica_SEM
        b["length"] = "10L"; b["time_ps_common"] = b.time_ps
        b["C_normalized"] = b.CJJ_mean_A2_fs2 / b.CJJ_mean_A2_fs2.iloc[0]
        b["C_normalized_SEM"] = b.CJJ_replica_SEM_A2_fs2 / b.CJJ_mean_A2_fs2.iloc[0]
        for frame, n_value in [(a, n5), (b, n10)]:
            rows.append(pd.DataFrame({"k_inv_A": item.k_inv_A_5L, "n_5L": n5, "n_10L": n10,
                                      "length": frame.length, "time_ps": frame.time_ps_common,
                                      "C_normalized": frame.C_normalized, "C_normalized_SEM": frame.C_normalized_SEM}))
    curves = pd.concat(rows, ignore_index=True)
    curves.to_csv(data_dir / "LA_matched_k_5L_10L_normalized_CJJ_time.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    axes[0].errorbar(merged.k_inv_A_5L, merged.CJJ_raw_mean_A2_fs2,
                     yerr=merged.CJJ_raw_replica_SEM_A2_fs2, fmt="o-", lw=1.0, ms=4,
                     capsize=2, color="#2166ac", label="5L")
    axes[0].errorbar(merged.k_inv_A_10L, merged.CJJ_mean_A2_fs2,
                     yerr=merged.CJJ_replica_SEM_A2_fs2, fmt="s-", lw=1.0, ms=3.8,
                     capsize=2, color="#b2182b", label="10L")
    axes[0].set(xlabel=r"matched $k$ ($\AA^{-1}$)", ylabel=r"raw $C_{J_zJ_z}(k,0)$ ($\AA^2\,\mathrm{fs}^{-2}$)")
    axes[0].legend(fontsize=6)
    axes[1].plot(merged.k_inv_A_5L, merged.weight_5L_shared_k_fraction, "o-", lw=1.0, ms=4, color="#2166ac", label="5L")
    axes[1].plot(merged.k_inv_A_10L, merged.weight_10L_shared_k_fraction, "s-", lw=1.0, ms=3.8, color="#b2182b", label="10L")
    axes[1].set(xlabel=r"matched $k$ ($\AA^{-1}$)", ylabel=r"fraction within five matched $k$ values")
    axes[1].legend(fontsize=6)
    for i, ax in enumerate(axes):
        ax.tick_params(top=True, right=True)
        ax.text(-0.18, 1.05, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.subplots_adjust(left=0.11, right=0.99, bottom=0.20, top=0.95, wspace=0.34)
    export(fig, figure_dir / "LA_matched_k_5L_10L_raw_and_relative_modal_weight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 5, figsize=(7.0, 1.9), sharex=True, sharey=True)
    for i, item in merged.iterrows():
        ax = axes[i]
        subset = curves[curves.n_5L == item.n_5L]
        for length, color in [("5L", "#2166ac"), ("10L", "#b2182b")]:
            frame = subset[subset.length == length]
            ax.plot(frame.time_ps, frame.C_normalized, color=color, lw=1.0, label=length)
        ax.axhline(0, color="0.45", lw=0.7)
        ax.set_title(rf"$k={item.k_inv_A_5L:.4f}$", fontsize=6)
        ax.tick_params(top=True, right=True)
        ax.text(-0.21, 1.07, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    axes[0].set_ylabel(r"normalized $C_{J_zJ_z}$")
    axes[2].set_xlabel(r"lag time (ps)")
    axes[0].legend(fontsize=5.5, loc="upper right")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.23, top=0.83, wspace=0.18)
    export(fig, figure_dir / "LA_matched_k_5L_10L_normalized_CJJ_time")
    plt.close(fig)

    inventory = pd.DataFrame([
        {"observable": "LA CJJ/raw modal weight", "lengths": "5L, 10L", "matched_k": "5 pairs: 0.01246..0.06231 A^-1", "raw_weight": "yes", "status": "available but protocol-distinct"},
        {"observable": "LA CJJ normalized time", "lengths": "3L, 4L, 5L, 10L", "matched_k": "only k_min per box in 3L/4L archive", "raw_weight": "no", "status": "shape-only, not broad matched-k coverage"},
        {"observable": "TA_theta m=0 normalized CJJ", "lengths": "2L, 10L", "matched_k": "2L n=1..6 <-> 10L n=5..30", "raw_weight": "no", "status": "protocol-matched nonzero-k shape comparison"},
        {"observable": "TA_theta m=0 normalized CJJ", "lengths": "5L, 10L", "matched_k": "5L n=1..8 <-> 10L n=2..16", "raw_weight": "no", "status": "protocol-distinct shape comparison"},
    ])
    inventory.to_csv(data_dir / "local_current_mode_matched_k_asset_coverage.csv", index=False)
    (args.output / "metadata.json").write_text(json.dumps({
        "question": "cross-length matched-physical-k current-mode consistency and modal weight",
        "raw_weight_definition": "CJJ(k,0), not normalized spectral PSD", "length_pair": "5L vs 10L",
        "matching": "5L n=m paired to 10L n=2m, m=1..5", "critical_limit": "5L and 10L raw currents use distinct cadence/duration/replica and connected-current protocols",
    }, indent=2))
    (args.output / "FINISHED.txt").write_text("Local matched-k current-mode weight audit finished successfully.\n")


if __name__ == "__main__":
    main()
