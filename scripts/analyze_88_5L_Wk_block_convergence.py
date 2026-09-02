#!/usr/bin/env python3
"""Time-block convergence audit for (8,8) 5L modal weights W_a(k)."""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
RAW_ROOT = Path(r"H:\gcmc_explore\translational_anomaly\08_viscosity_friction_length_scaling\04_analysis\offline_frequency_viscosity_20260803\raw_case_directories")
OUT = ROOT / "results" / "collective_mode_response" / "88_5L_Wk_block_convergence" / "2026-08-19"
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_88_5L_full_dispersion import read_modal_currents  # noqa: E402

NMAX = 20
N_OXYGEN = 665
BLOCK_PS = 50.0
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


def dump_path(rep: int) -> Path:
    return RAW_ROOT / f"NVT20ns_5xL_8_8_RH75_N665_rep{rep}_20260719" / f"nvt20ns_8_8_RH75_5L_rep{rep}_oxygen_10fs_1ns.dump"


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save(fig: plt.Figure, stem: Path) -> None:
    for suffix, kwargs in ((".png", {"dpi": 600}), (".pdf", {}), (".svg", {}), (".tiff", {"dpi": 600})):
        fig.savefig(stem.with_suffix(suffix), facecolor="white", **kwargs)
    plt.close(fig)


def main() -> None:
    data_dir, fig_dir = OUT / "derived_data", OUT / "figures"
    data_dir.mkdir(parents=True, exist_ok=True); fig_dir.mkdir(exist_ok=True)
    blocks: list[dict] = []
    full: list[dict] = []
    meta_all: list[dict] = []
    for rep in (1, 2, 3):
        series, meta = read_modal_currents(dump_path(rep), NMAX)
        meta_all.append(meta)
        block_frames = int(round(BLOCK_PS / meta["dt_ps"]))
        nblock = len(series) // block_frames
        k = 2 * np.pi * np.arange(1, NMAX + 1) / meta["box_A"][2]
        centered = series - series.mean(axis=0, keepdims=True)
        full_w = np.mean(np.abs(centered) ** 2, axis=0).real / N_OXYGEN
        for bi, branch in enumerate(BRANCHES):
            for mi, n in enumerate(range(1, NMAX + 1)):
                full.append({"replicate": rep, "branch": branch, "n": n, "k_inv_A": float(k[mi]),
                             "W_full_A2_fs2_per_O": float(full_w[mi, bi])})
        for block in range(nblock):
            start, stop = block * block_frames, (block + 1) * block_frames
            chunk = series[start:stop]
            chunk = chunk - chunk.mean(axis=0, keepdims=True)
            weight = np.mean(np.abs(chunk) ** 2, axis=0).real / N_OXYGEN
            for bi, branch in enumerate(BRANCHES):
                for mi, n in enumerate(range(1, NMAX + 1)):
                    blocks.append({"replicate": rep, "block_index": block + 1,
                                   "block_start_ps": float(start * meta["dt_ps"]), "block_stop_ps": float(stop * meta["dt_ps"]),
                                   "branch": branch, "n": n, "k_inv_A": float(k[mi]),
                                   "W_block_A2_fs2_per_O": float(weight[mi, bi])})
        del series
    summary: list[dict] = []
    for branch in BRANCHES:
        for n in range(1, NMAX + 1):
            values = np.array([r["W_block_A2_fs2_per_O"] for r in blocks if r["branch"] == branch and r["n"] == n])
            fvalues = np.array([r["W_full_A2_fs2_per_O"] for r in full if r["branch"] == branch and r["n"] == n])
            summary.append({"branch": branch, "n": n, "k_inv_A": next(r["k_inv_A"] for r in blocks if r["branch"] == branch and r["n"] == n),
                            "n_blocks_total": len(values), "W_block_mean_A2_fs2_per_O": float(values.mean()),
                            "W_block_SD_A2_fs2_per_O": float(values.std(ddof=1)), "W_block_CV": float(values.std(ddof=1) / values.mean()),
                            "W_full_replica_mean_A2_fs2_per_O": float(fvalues.mean()),
                            "W_full_replica_SEM_A2_fs2_per_O": float(fvalues.std(ddof=1) / math.sqrt(3)),
                            "block_to_full_ratio": float(values.mean() / fvalues.mean()),
                            "boundary": "50 ps blocks are not assumed independent; CV diagnoses temporal nonstationarity/finite-sampling sensitivity only"})
    write_csv(data_dir / "Wk_50ps_blocks_per_replica.csv", blocks)
    write_csv(data_dir / "Wk_full_trajectory_per_replica.csv", full)
    write_csv(data_dir / "Wk_50ps_block_convergence_summary.csv", summary)

    fig = plt.figure(figsize=(7.0, 2.65))
    ax_a = fig.add_axes([0.10, 0.22, 0.37, 0.67]); ax_b = fig.add_axes([0.59, 0.22, 0.34, 0.67])
    for branch in BRANCHES:
        color = COLORS[branch]
        for rep in (1, 2, 3):
            b = sorted((r for r in blocks if r["branch"] == branch and r["n"] == 1 and r["replicate"] == rep), key=lambda r: r["block_index"])
            f = next(r["W_full_A2_fs2_per_O"] for r in full if r["branch"] == branch and r["n"] == 1 and r["replicate"] == rep)
            ax_a.plot([r["block_start_ps"] + BLOCK_PS / 2 for r in b], [r["W_block_A2_fs2_per_O"] / f for r in b],
                      "o-", color=color, alpha=0.32, ms=2.4, lw=1.0)
        ss = [r for r in summary if r["branch"] == branch]
        ax_b.plot([r["k_inv_A"] for r in ss], [100 * r["W_block_CV"] for r in ss], "o-", color=color, ms=3.0, lw=1.1, label=LABELS[branch])
    ax_a.axhline(1.0, color="#777777", lw=1.0); ax_a.set(xlabel="Block midpoint (ps)", ylabel=r"$W_a(k_1)$ block/full")
    ax_b.set(xlabel=r"$k_n$ (Å$^{-1}$)", ylabel="50 ps block CV (%)"); ax_b.legend(fontsize=6.2, loc="upper right")
    for ax, label in ((ax_a, "(a)"), (ax_b, "(b)")):
        ax.text(-0.20, 1.02, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom"); ax.tick_params(length=3, labelsize=6.5)
    fig.text(0.285, 0.95, r"$n=1$; thin lines: three velocity-seed replicas", ha="center", fontsize=6.5)
    fig.text(0.76, 0.95, "All n=1–20; block CV is a convergence diagnostic", ha="center", fontsize=6.5)
    save(fig, fig_dir / "Wk_50ps_block_convergence")
    (OUT / "metadata.json").write_text(json.dumps({"input_replicas": meta_all, "block_ps": BLOCK_PS, "n_range": [1, NMAX], "definition": "W=CJJ(0)/N_O evaluated as the within-block temporal variance of the complex modal current", "limitation": "Blocks are not independent replicates; their dispersion is a finite-sampling/stationarity diagnostic."}, indent=2), encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("50 ps block convergence audit completed.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
