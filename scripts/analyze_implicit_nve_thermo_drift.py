"""Extract and visualize temperature and total-energy drift from LAMMPS NVE logs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_log(path: Path) -> pd.DataFrame:
    rows: list[list[float]] = []
    active = False
    with path.open(errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if fields[:4] == ["Step", "Temp", "PotEng", "f_cnt_force"]:
                active = True
                continue
            if active and len(fields) >= 5:
                try:
                    rows.append([float(fields[0]), float(fields[1]), float(fields[4])])
                    continue
                except ValueError:
                    active = False
    frame = pd.DataFrame(rows, columns=["step", "temperature_K", "total_energy_kcal_mol"])
    frame = frame.drop_duplicates("step", keep="last").sort_values("step")
    frame["time_ps"] = frame.step * 0.0005
    return frame


def linear_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    coef = np.polyfit(x, y, 1)
    return float(coef[0]), float(coef[1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--logs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data_dir = args.output / "derived_data"
    fig_dir = args.output / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(exist_ok=True)

    frames: list[pd.DataFrame] = []
    summary: list[dict[str, float | str]] = []
    for i, log in enumerate(args.logs, start=1):
        frame = read_log(log)
        if len(frame) < 3:
            raise ValueError(f"No usable thermo data in {log}")
        frame.insert(0, "replica", i)
        frames.append(frame)
        t = frame.time_ps.to_numpy()
        temp = frame.temperature_K.to_numpy()
        eng = frame.total_energy_kcal_mol.to_numpy()
        ts, _ = linear_slope(t, temp)
        es, _ = linear_slope(t, eng)
        summary.append({
            "replica": i,
            "n_thermo_points": len(frame),
            "duration_ps": float(t[-1] - t[0]),
            "temperature_start_K": float(temp[0]), "temperature_end_K": float(temp[-1]),
            "temperature_delta_K": float(temp[-1] - temp[0]), "temperature_slope_K_per_ns": ts * 1000,
            "total_energy_start_kcal_mol": float(eng[0]), "total_energy_end_kcal_mol": float(eng[-1]),
            "total_energy_delta_kcal_mol": float(eng[-1] - eng[0]),
            "total_energy_slope_kcal_mol_per_ns": es * 1000,
            "relative_energy_slope_per_ns": es / abs(float(eng[0])) * 1000,
        })
    full = pd.concat(frames, ignore_index=True)
    stats = pd.DataFrame(summary)
    full.to_csv(data_dir / "implicitC88_N1600_NVE_thermo_timeseries.csv", index=False)
    stats.to_csv(data_dir / "implicitC88_N1600_NVE_thermo_drift_summary.csv", index=False)
    avg = stats.mean(numeric_only=True).to_dict()
    sem = stats.sem(numeric_only=True).to_dict()
    (data_dir / "implicitC88_N1600_NVE_thermo_drift_summary.json").write_text(json.dumps({"replica_mean": avg, "replica_SEM": sem}, indent=2))

    fig, axes = plt.subplots(2, 1, figsize=(8, 6.6), sharex=True)
    colors = plt.cm.tab10.colors
    for rep, frame in full.groupby("replica"):
        c = colors[(int(rep) - 1) % len(colors)]
        axes[0].plot(frame.time_ps, frame.temperature_K, lw=0.9, alpha=0.75, color=c, label=f"rep {rep}")
        axes[1].plot(frame.time_ps, frame.total_energy_kcal_mol, lw=0.9, alpha=0.75, color=c)
    axes[0].set_ylabel("Temperature (K)")
    axes[1].set_ylabel(r"Total energy (kcal mol$^{-1}$)")
    axes[1].set_xlabel("NVE time (ps)")
    axes[0].legend(ncol=4, frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "implicitC88_N1600_NVE_temperature_energy_drift.png", dpi=300)
    fig.savefig(fig_dir / "implicitC88_N1600_NVE_temperature_energy_drift.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
