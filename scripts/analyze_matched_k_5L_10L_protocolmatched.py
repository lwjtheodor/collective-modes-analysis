#!/usr/bin/env python3
"""Protocol-matched 5L/10L longitudinal current-mode comparison.

Reads the 100-fs, 10-ns, eight-replica ``id z vz`` archives directly.  Both
lengths use the identical observable and estimator:

    J_z(k,t) = sum_i [v_zi(t) - mean_i v_zi(t)] exp(i k z_i(t)).

Physical wave numbers are matched as 5L n=m to 10L n=2m (m=1..5).  This
replaces the earlier comparison that mixed a short 5L trajectory protocol
with the long 10L archive.
"""
from __future__ import annotations

import csv
import json
import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(r"H:\gcmc_explore\translational_anomaly\02_isf_collective_modes")
SOURCE = Path(r"F:\ccfep_gcmc_archive_20260814\stage_vacf_tail_8_8_L2L10_8rep_weaknh_zvz_20260812")
OUT = ROOT / "results" / "collective_mode_response" / "88_5L_10L_LA_matched_k_protocolmatched_100fs_10ns_8rep" / "2026-08-28"
EXISTING_10L = ROOT / "results" / "collective_mode_response" / "88_10L_LA_CJJ_Skw_100fs_10ns_8rep_n001_n010" / "2026-08-24"

NREP = 8
NFRAME = 100_001
DT_PS = 0.1
MAX_LAG_PS = 1_000.0
NPERSEG = 16_384
MODES_5L = np.arange(1, 6, dtype=int)
MODES_10L = 2 * MODES_5L

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 1.0, "axes.spines.right": False,
    "axes.spines.top": False, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0, "ytick.major.size": 3.0,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


@dataclass
class Case:
    length_l: int
    modes: np.ndarray
    lz_a: float | None = None
    n_oxygen: int | None = None

    @property
    def label(self) -> str:
        return f"{self.length_l}L"

    def path(self, replica: int) -> Path:
        return SOURCE / f"{self.length_l}L_rep{replica}" / (
            f"VACF88_{self.length_l}L_tail_zvz_rep{replica}_oxygen_id_z_vz_100fs.dump"
        )


def read_modal(case: Case, replica: int) -> np.ndarray:
    """Stream one dump and return the complex currents for requested modes."""
    result = np.empty((NFRAME, len(case.modes)), dtype=np.complex128)
    steps: list[int] = []
    with case.path(replica).open(encoding="utf-8", errors="replace") as handle:
        frame = 0
        while True:
            marker = handle.readline()
            if not marker:
                break
            if marker.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Unexpected dump marker in {case.path(replica)}: {marker!r}")
            step = int(handle.readline())
            if handle.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("Missing atom-count marker")
            natom = int(handle.readline())
            if handle.readline().strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError("Missing box marker")
            bounds = np.array([[float(x) for x in handle.readline().split()] for _ in range(3)])
            header = handle.readline().split()[2:]
            columns = {name: index for index, name in enumerate(header)}
            missing = {"id", "z", "vz"} - set(columns)
            if missing:
                raise ValueError(f"{case.path(replica)} lacks {missing}")
            z = np.empty(natom, dtype=float)
            vz = np.empty(natom, dtype=float)
            for atom in range(natom):
                row = handle.readline().split()
                z[atom] = float(row[columns["z"]])
                vz[atom] = float(row[columns["vz"]])
            lz = float(bounds[2, 1] - bounds[2, 0])
            if case.lz_a is None:
                case.lz_a = lz
                case.n_oxygen = natom
            elif not np.isclose(case.lz_a, lz) or case.n_oxygen != natom:
                raise ValueError(f"Box/oxygen count changes within {case.label}")
            k = 2.0 * np.pi * case.modes / lz
            drift_free_vz = vz - vz.mean()
            result[frame] = np.sum(drift_free_vz[:, None] * np.exp(1j * z[:, None] * k), axis=0)
            steps.append(step)
            frame += 1
    if frame != NFRAME:
        raise ValueError(f"{case.path(replica)}: expected {NFRAME} frames; got {frame}")
    if not np.all(np.diff(steps) == 200):
        raise ValueError(f"{case.path(replica)} has nonuniform 100-fs cadence")
    return result


def connected_acf(series: np.ndarray) -> np.ndarray:
    centered = series - series.mean(axis=0, keepdims=True)
    nframes = len(centered)
    fft = np.fft.fft(centered, n=2 * nframes, axis=0)
    numerator = np.fft.ifft(fft * np.conj(fft), axis=0).real[: int(MAX_LAG_PS / DT_PS) + 1]
    return numerator / np.arange(nframes, nframes - len(numerator), -1)[:, None]


def welch_psd(series: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    starts = np.arange(0, len(series) - NPERSEG + 1, NPERSEG // 2)
    window = np.hanning(NPERSEG)[:, None]
    power = np.zeros((NPERSEG, series.shape[1]), dtype=float)
    for start in starts:
        segment = series[start:start + NPERSEG]
        segment = segment - segment.mean(axis=0, keepdims=True)
        transform = np.fft.fft(segment * window, axis=0)
        power += np.abs(transform) ** 2 / np.sum(window[:, 0] ** 2)
    power /= len(starts)
    positive = np.arange(1, NPERSEG // 2 + 1)
    frequency = positive / (NPERSEG * DT_PS)
    symmetric_power = 0.5 * (power[positive] + power[-positive])
    return frequency, symmetric_power, len(starts)


def export(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


def write_cjj_tables(data_dir: Path, case: Case, cjj: np.ndarray) -> None:
    time = np.arange(cjj.shape[1]) * DT_PS
    k = 2 * np.pi * case.modes / float(case.lz_a)
    with (data_dir / f"CJJ_{case.label}_per_replica.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["replica", "branch", "n", "k_inv_A", "time_ps", "CJJ_raw_A2_fs2"])
        for replica in range(NREP):
            for mode_index, mode in enumerate(case.modes):
                writer.writerows((replica + 1, "LA", mode, k[mode_index], t, value)
                                 for t, value in zip(time, cjj[replica, :, mode_index]))
    mean = cjj.mean(axis=0)
    sem = cjj.std(axis=0, ddof=1) / np.sqrt(NREP)
    with (data_dir / f"CJJ_{case.label}_ensemble_mean_sem.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["branch", "n", "k_inv_A", "time_ps", "CJJ_mean_A2_fs2", "CJJ_replica_SEM_A2_fs2"])
        for mode_index, mode in enumerate(case.modes):
            writer.writerows(("LA", mode, k[mode_index], t, value, err)
                             for t, value, err in zip(time, mean[:, mode_index], sem[:, mode_index]))


def write_spectrum_tables(data_dir: Path, case: Case, frequency: np.ndarray, psd: np.ndarray) -> None:
    k = 2 * np.pi * case.modes / float(case.lz_a)
    with (data_dir / f"SJJ_{case.label}_per_replica.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["replica", "branch", "n", "k_inv_A", "frequency_ps_inv", "omega_rad_ps", "PSD_arbitrary"])
        for replica in range(NREP):
            for mode_index, mode in enumerate(case.modes):
                writer.writerows((replica + 1, "LA", mode, k[mode_index], freq, 2 * np.pi * freq, value)
                                 for freq, value in zip(frequency, psd[replica, :, mode_index]))
    mean = psd.mean(axis=0)
    sem = psd.std(axis=0, ddof=1) / np.sqrt(NREP)
    with (data_dir / f"SJJ_{case.label}_ensemble_mean_sem.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["branch", "n", "k_inv_A", "frequency_ps_inv", "omega_rad_ps", "PSD_mean_arbitrary", "PSD_replica_SEM_arbitrary"])
        for mode_index, mode in enumerate(case.modes):
            writer.writerows(("LA", mode, k[mode_index], freq, 2 * np.pi * freq, value, err)
                             for freq, value, err in zip(frequency, mean[:, mode_index], sem[:, mode_index]))


def load_existing_10l(case: Case) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Load the verified 10L per-replica products made with this same estimator."""
    source_meta = json.loads((EXISTING_10L / "metadata.json").read_text(encoding="utf-8"))
    if source_meta["dt_ps"] != DT_PS or source_meta["welch_nperseg_frames"] != NPERSEG:
        raise ValueError("Existing 10L analysis does not match the requested estimator")
    # The archived metadata predates a field for the explicit definition; its
    # source script (analyze_88_10L_LA_zvz_100fs_10ns_8rep.py) is the checked
    # provenance and implements the same instantaneous O-COM subtraction.
    case.lz_a = float(source_meta["Lz_A"][0])
    case.n_oxygen = 1330
    times = np.arange(int(MAX_LAG_PS / DT_PS) + 1) * DT_PS
    cjj = np.empty((NREP, len(times), len(case.modes)), dtype=float)
    with (EXISTING_10L / "derived_data" / "CJJ_all_modes_per_replica.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            mode = int(row["n"])
            if mode not in set(case.modes):
                continue
            replica = int(row["replica"]) - 1
            frame = int(round(float(row["time_ps"]) / DT_PS))
            cjj[replica, frame, int(np.where(case.modes == mode)[0][0])] = float(row["CJJ_raw_A2_fs2"])
    frequency = np.arange(1, NPERSEG // 2 + 1) / (NPERSEG * DT_PS)
    psd = np.empty((NREP, len(frequency), len(case.modes)), dtype=float)
    with (EXISTING_10L / "derived_data" / "current_spectra_all_modes_per_replica.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            mode = int(row["n"])
            if mode not in set(case.modes):
                continue
            replica = int(row["replica"]) - 1
            frame = int(round(float(row["frequency_ps_inv"]) * NPERSEG * DT_PS)) - 1
            psd[replica, frame, int(np.where(case.modes == mode)[0][0])] = float(row["PSD_arbitrary"])
    return cjj, psd, frequency, int(source_meta["welch_segments_per_replica"])


def load_completed_5l(case: Case) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Reload the just-finished direct 5L pass after an intentional restart."""
    cjj_path = OUT / "derived_data" / "CJJ_5L_per_replica.csv"
    sjj_path = OUT / "derived_data" / "SJJ_5L_per_replica.csv"
    if not cjj_path.exists() or not sjj_path.exists() or sjj_path.stat().st_size == 0:
        raise FileNotFoundError("Completed 5L tables are not available for reuse")
    case.lz_a, case.n_oxygen = 504.19999, 665
    cjj = np.empty((NREP, int(MAX_LAG_PS / DT_PS) + 1, len(case.modes)), dtype=float)
    with cjj_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            replica = int(row["replica"]) - 1
            mode = int(row["n"])
            frame = int(round(float(row["time_ps"]) / DT_PS))
            cjj[replica, frame, int(np.where(case.modes == mode)[0][0])] = float(row["CJJ_raw_A2_fs2"])
    frequency = np.arange(1, NPERSEG // 2 + 1) / (NPERSEG * DT_PS)
    psd = np.empty((NREP, len(frequency), len(case.modes)), dtype=float)
    with sjj_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            replica = int(row["replica"]) - 1
            mode = int(row["n"])
            frame = int(round(float(row["frequency_ps_inv"]) * NPERSEG * DT_PS)) - 1
            psd[replica, frame, int(np.where(case.modes == mode)[0][0])] = float(row["PSD_arbitrary"])
    return cjj, psd, frequency, 11


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-completed-5l", action="store_true", help="Reuse the direct 5L pass completed in this output package.")
    parser.add_argument("--reuse-verified-10l", action="store_true", help="Reuse the verified 10L per-replica tables with the identical estimator.")
    args = parser.parse_args()
    data_dir = OUT / "derived_data"
    figure_dir = OUT / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(exist_ok=True)
    cases = [Case(5, MODES_5L), Case(10, MODES_10L)]
    all_cjj: dict[int, np.ndarray] = {}
    all_psd: dict[int, np.ndarray] = {}
    frequency = None
    segment_count = None
    for case in cases:
        if case.length_l == 5 and args.reuse_completed_5l:
            print("Loading completed direct 5L CJJ/SJJ tables", flush=True)
            all_cjj[5], all_psd[5], frequency, segment_count = load_completed_5l(case)
            write_cjj_tables(data_dir, case, all_cjj[5])
            write_spectrum_tables(data_dir, case, frequency, all_psd[5])
            continue
        if case.length_l == 10 and args.reuse_verified_10l:
            print("Loading verified existing 10L CJJ/SJJ tables", flush=True)
            all_cjj[10], all_psd[10], frequency, segment_count = load_existing_10l(case)
            write_cjj_tables(data_dir, case, all_cjj[10])
            write_spectrum_tables(data_dir, case, frequency, all_psd[10])
            continue
        cjj_per_replica, psd_per_replica = [], []
        for replica in range(1, NREP + 1):
            print(f"Reading {case.label}, replica {replica}/{NREP}", flush=True)
            current = read_modal(case, replica)
            cjj_per_replica.append(connected_acf(current))
            frequency, spectrum, segment_count = welch_psd(current)
            psd_per_replica.append(spectrum)
        all_cjj[case.length_l] = np.asarray(cjj_per_replica)
        all_psd[case.length_l] = np.asarray(psd_per_replica)
        write_cjj_tables(data_dir, case, all_cjj[case.length_l])
        write_spectrum_tables(data_dir, case, frequency, all_psd[case.length_l])

    # Matched-k raw contribution table: each replica remains a statistical unit.
    c0_5 = all_cjj[5][:, 0, :]
    c0_10 = all_cjj[10][:, 0, :]
    k5 = 2 * np.pi * MODES_5L / float(cases[0].lz_a)
    k10 = 2 * np.pi * MODES_10L / float(cases[1].lz_a)
    rows = []
    for index, (n5, n10, kval5, kval10) in enumerate(zip(MODES_5L, MODES_10L, k5, k10)):
        for replica in range(NREP):
            rows.append({"replica": replica + 1, "n_5L": int(n5), "n_10L": int(n10),
                         "k_5L_inv_A": kval5, "k_10L_inv_A": kval10,
                         "relative_k_mismatch": abs(kval5 - kval10) / kval5,
                         "CJJ0_5L_A2_fs2": c0_5[replica, index], "CJJ0_10L_A2_fs2": c0_10[replica, index],
                         "CJJ0_5L_per_oxygen": c0_5[replica, index] / cases[0].n_oxygen,
                         "CJJ0_10L_per_oxygen": c0_10[replica, index] / cases[1].n_oxygen,
                         "ratio_10L_over_5L": c0_10[replica, index] / c0_5[replica, index]})
    with (data_dir / "LA_matched_k_per_replica_CJJ0.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    summary = []
    c0_5_mean, c0_10_mean = c0_5.mean(0), c0_10.mean(0)
    c0_5_sem = c0_5.std(0, ddof=1) / np.sqrt(NREP)
    c0_10_sem = c0_10.std(0, ddof=1) / np.sqrt(NREP)
    for index, n5 in enumerate(MODES_5L):
        ratios = c0_10[:, index] / c0_5[:, index]
        summary.append({"n_5L": int(n5), "n_10L": int(MODES_10L[index]), "k_inv_A": k5[index],
                        "relative_k_mismatch": abs(k5[index] - k10[index]) / k5[index],
                        "CJJ0_5L_mean_A2_fs2": c0_5_mean[index], "CJJ0_5L_replica_SEM_A2_fs2": c0_5_sem[index],
                        "CJJ0_10L_mean_A2_fs2": c0_10_mean[index], "CJJ0_10L_replica_SEM_A2_fs2": c0_10_sem[index],
                        "CJJ0_5L_per_oxygen": c0_5_mean[index] / cases[0].n_oxygen,
                        "CJJ0_10L_per_oxygen": c0_10_mean[index] / cases[1].n_oxygen,
                        "ratio_10L_over_5L_mean": ratios.mean(),
                        "ratio_10L_over_5L_replica_SEM": ratios.std(ddof=1) / np.sqrt(NREP),
                        "weight_5L_within_matched_set": c0_5_mean[index] / c0_5_mean.sum(),
                        "weight_10L_within_matched_set": c0_10_mean[index] / c0_10_mean.sum()})
    with (data_dir / "LA_matched_k_CJJ0_weight_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)

    # A figure of absolute/one-particle weights and matched normalized CJJ shapes.
    fig = plt.figure(figsize=(7.0, 5.3))
    grid = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.25], left=0.105, right=0.99, bottom=0.10, top=0.97, hspace=0.44)
    top = grid[0].subgridspec(1, 2, wspace=0.34)
    ax_abs = fig.add_subplot(top[0, 0])
    ax_unit = fig.add_subplot(top[0, 1])
    values = [(c0_5_mean, c0_5_sem, "5L", "o", "#2166ac"), (c0_10_mean, c0_10_sem, "10L", "s", "#b2182b")]
    for mean, sem, label, marker, color in values:
        ax_abs.errorbar(k5, mean, yerr=sem, fmt=marker + "-", color=color, ms=4, lw=1.1, capsize=2, label=label)
        oxygen = cases[0].n_oxygen if label == "5L" else cases[1].n_oxygen
        ax_unit.errorbar(k5, mean / oxygen, yerr=sem / oxygen, fmt=marker + "-", color=color, ms=4, lw=1.1, capsize=2, label=label)
    ax_abs.set(xlabel=r"matched $k$ ($\mathrm{\AA}^{-1}$)", ylabel=r"$C_{J_zJ_z}(k,0)$ ($\mathrm{\AA}^2\,\mathrm{fs}^{-2}$)")
    ax_unit.set(xlabel=r"matched $k$ ($\mathrm{\AA}^{-1}$)", ylabel=r"$C_{J_zJ_z}(k,0)/N_\mathrm{O}$")
    for label, ax in [("(a)", ax_abs), ("(b)", ax_unit)]:
        ax.legend(fontsize=6, loc="best")
        ax.text(-0.19, 1.05, label, transform=ax.transAxes, fontweight="bold", fontsize=9)

    bottom = grid[1].subgridspec(1, 5, wspace=0.22)
    time = np.arange(all_cjj[5].shape[1]) * DT_PS
    for index, (n5, n10, kval) in enumerate(zip(MODES_5L, MODES_10L, k5)):
        ax = fig.add_subplot(bottom[0, index])
        for length, color, label in [(5, "#2166ac", "5L"), (10, "#b2182b", "10L")]:
            trace = all_cjj[length].mean(0)[:, index]
            trace_sem = all_cjj[length].std(0, ddof=1)[:, index] / np.sqrt(NREP)
            normalizer = trace[0]
            ax.plot(time, trace / normalizer, color=color, lw=1.0, label=label)
            ax.fill_between(time, (trace - trace_sem) / normalizer, (trace + trace_sem) / normalizer, color=color, alpha=0.18, linewidth=0)
        ax.axhline(0, color="0.45", lw=0.7)
        ax.set(xlim=(0, 100), ylim=(-0.7, 1.10), title=rf"$k={kval:.4f}$")
        ax.tick_params(labelsize=6)
        ax.text(-0.24, 1.08, f"({chr(99 + index)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
        if index == 0:
            ax.set_ylabel(r"$C_{J_zJ_z}(k,t)/C_{J_zJ_z}(k,0)$")
            ax.legend(fontsize=5.5, loc="lower left")
        if index == 2:
            ax.set_xlabel(r"lag time (ps)")
    export(fig, figure_dir / "LA_matched_k_CJJ0_and_normalized_time")
    plt.close(fig)

    # Spectral shape is shown separately so negative-frequency symmetry is explicit.
    fig = plt.figure(figsize=(7.0, 2.05))
    grid = fig.add_gridspec(1, 5, left=0.085, right=0.995, bottom=0.22, top=0.83, wspace=0.22)
    omega = 2 * np.pi * frequency
    for index, kval in enumerate(k5):
        ax = fig.add_subplot(grid[0, index])
        for length, color, label in [(5, "#2166ac", "5L"), (10, "#b2182b", "10L")]:
            signal = all_psd[length].mean(0)[:, index]
            signal /= np.trapz(signal, omega)
            y = np.maximum(signal, signal.max() * 1e-7)
            ax.semilogy(np.r_[-omega[::-1], omega], np.r_[y[::-1], y], color=color, lw=1.0, label=label)
        ax.axvline(0, color="0.45", lw=0.7)
        ax.set(xlim=(-1.2, 1.2), ylim=(1e-5, None), title=rf"$k={kval:.4f}$")
        ax.tick_params(labelsize=6)
        ax.text(-0.22, 1.09, f"({chr(97 + index)})", transform=ax.transAxes, fontweight="bold", fontsize=9)
        if index == 0:
            ax.set_ylabel(r"normalized $S_{J_zJ_z}$")
            ax.legend(fontsize=5.5, loc="lower left")
        if index == 2:
            ax.set_xlabel(r"$\omega$ ($\mathrm{rad\,ps^{-1}}$)")
    export(fig, figure_dir / "LA_matched_k_normalized_SJJ_shapes")
    plt.close(fig)

    metadata = {
        "purpose": "Protocol-matched 5L/10L matched-physical-k longitudinal-current comparison",
        "source": str(SOURCE), "systems": {"5L": {"replicas": NREP, "modes": MODES_5L.tolist(), "processed_directly_this_run": True, "direct_pass_reloaded_after_restart": bool(args.reuse_completed_5l)}, "10L": {"replicas": NREP, "modes": MODES_10L.tolist(), "processed_directly_this_run": not args.reuse_verified_10l}},
        "fields": ["id", "z", "vz"], "frames_per_replica": NFRAME, "dt_ps": DT_PS, "duration_ps": NFRAME * DT_PS - DT_PS,
        "current_definition": "J_z(k,t)=sum_O [v_z(t)-mean_O v_z(t)] exp(i k z(t))",
        "matching": "5L n=m paired to 10L n=2m, m=1..5", "cjj_estimator": "connected complex-current FFT ACF; unbiased lag normalization",
        "cjj_maxlag_ps": MAX_LAG_PS, "welch": {"nperseg_frames": NPERSEG, "segment_ps": NPERSEG * DT_PS, "overlap": 0.5, "segments_per_replica": segment_count},
        "uncertainty": "replica SEM across eight velocity-seed replicas", "interpretation_limit": "CJJ(0) is an equal-time modal variance; its raw length scaling and normalized matched-set fractions are reported separately. It is not a DOS-normalized total spectral weight.",
        "10L_reuse_provenance": str(EXISTING_10L) if args.reuse_verified_10l else None,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(
        "# Matched-k LA current modes: 5L versus 10L\n\n"
        "Both systems are re-analysed directly from the same 100-fs, 10-ns, eight-replica `id z vz` archive batch. "
        "The current is instantaneous oxygen-COM-subtracted. Physical matching is 5L `n=m` to 10L `n=2m`, m=1..5. "
        "`derived_data/LA_matched_k_CJJ0_weight_summary.csv` is the compact primary table; per-replica CJJ, spectra, and both ensembles are retained for reuse.\n",
        encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("Protocol-matched 5L/10L LA matched-k analysis finished successfully.\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "k_inv_A": k5.tolist(), "n_oxygen": {case.label: case.n_oxygen for case in cases}}, indent=2))


if __name__ == "__main__":
    main()
