#!/usr/bin/env python3
"""Block-FFT self-ISF scan for high-time-resolution oxygen trajectories.

Writes logarithmically sampled all-origin F_s(k,t) curves and descriptive KWW
fits.  The KWW beta quantifies finite-k non-exponentiality; it is not a
finite-size dynamic exponent.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

# The staged (8,8) length-scaling production input uses ``timestep 0.5``
# in LAMMPS real units, i.e. 0.5 fs = 0.0005 ps.
TIMESTEP_PS = 0.0005


def next_power_two(value):
    return 1 << (int(value - 1).bit_length())


def read_oxygen_z(path):
    frames, steps, lz_values, reference_mols = [], [], [], None
    with open(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError("expected timestep header")
            step = int(handle.readline())
            handle.readline(); natoms = int(handle.readline())
            handle.readline()
            bounds = [list(map(float, handle.readline().split())) for _ in range(3)]
            lz = bounds[2][1] - bounds[2][0]
            header = handle.readline().split()[2:]
            index = {name: i for i, name in enumerate(header)}
            required = ("mol", "type", "z", "iz")
            if any(name not in index for name in required):
                raise ValueError("oxygen dump lacks %s" % required)
            z_by_mol = {}
            for _ in range(natoms):
                fields = handle.readline().split()
                if int(float(fields[index["type"]])) == 3:
                    mol = int(float(fields[index["mol"]]))
                    z_by_mol[mol] = float(fields[index["z"]]) + int(float(fields[index["iz"]])) * lz
            mols = np.asarray(sorted(z_by_mol), dtype=np.int64)
            if reference_mols is None:
                reference_mols = mols
            elif not np.array_equal(reference_mols, mols):
                raise ValueError("molecule IDs changed at timestep %d" % step)
            frames.append([z_by_mol[mol] for mol in reference_mols])
            steps.append(step); lz_values.append(lz)
    steps = np.asarray(steps, dtype=np.int64)
    dt = float(np.median(np.diff(steps)) * TIMESTEP_PS)
    return np.asarray(frames, dtype=np.float64), float(np.median(lz_values)), dt


def log_indices(max_lag, npoints=420):
    values = np.geomspace(1.0, float(max_lag), npoints)
    return np.unique(np.r_[0, values.astype(np.int64), max_lag])


def self_isf(z, k, max_lag, block_size):
    frames, nmol = z.shape
    nfft = next_power_two(frames + max_lag)
    accum = np.zeros(max_lag + 1, dtype=np.complex128)
    for start in range(0, nmol, block_size):
        phase = np.exp(1j * k * z[:, start:start + block_size])
        spectrum = np.fft.fft(phase, n=nfft, axis=0)
        corr = np.fft.ifft(spectrum * np.conjugate(spectrum), axis=0)[:max_lag + 1]
        accum += corr.sum(axis=1)
    counts = np.arange(frames, frames - max_lag - 1, -1, dtype=float)
    return (accum.real / (counts * nmol))


def kww_fit(t, f):
    mask = (t >= 0.02) & (f > 0.08) & (f < 0.95)
    if int(mask.sum()) < 12:
        return {"beta_kww": None, "tau_kww_ps": None, "r2_log": None, "n_fit": int(mask.sum())}
    x = np.log(t[mask]); y = np.log(-np.log(f[mask]))
    beta, intercept = np.polyfit(x, y, 1)
    if beta <= 0 or not np.isfinite(beta):
        return {"beta_kww": None, "tau_kww_ps": None, "r2_log": None, "n_fit": int(mask.sum())}
    fit = beta * x + intercept
    denom = float(np.sum((y - y.mean()) ** 2))
    r2 = None if denom == 0 else float(1.0 - np.sum((y - fit) ** 2) / denom)
    return {"beta_kww": float(beta), "tau_kww_ps": float(np.exp(-intercept / beta)),
            "r2_log": r2, "n_fit": int(mask.sum())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--modes", required=True, help="comma-separated integer n values")
    parser.add_argument("--max-lag-ps", type=float, default=800.0)
    parser.add_argument("--block-size", type=int, default=32)
    args = parser.parse_args()
    modes = [int(value) for value in args.modes.split(",")]
    z, lz, dt = read_oxygen_z(args.dump)
    max_lag = min(int(round(args.max_lag_ps / dt)), len(z) - 1)
    sampled = log_indices(max_lag)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    curves, summary = [], []
    for n in modes:
        k = 2.0 * math.pi * n / lz
        fs = self_isf(z, k, max_lag, args.block_size)
        fit = kww_fit(sampled * dt, fs[sampled])
        item = {"label": args.label, "n": n, "k_inv_A": k, "lambda_A": 2.0 * math.pi / k,
                "dt_ps": dt, "n_frames": int(len(z)), "n_water": int(z.shape[1])}
        item.update(fit); summary.append(item)
        for i in sampled:
            curves.append({"label": args.label, "n": n, "k_inv_A": k, "lambda_A": 2.0 * math.pi / k,
                           "lag_ps": float(i * dt), "F_self": float(fs[i])})
    with open(args.out_dir / "fs_wavevector_scan_curves.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "n", "k_inv_A", "lambda_A", "lag_ps", "F_self"])
        writer.writeheader(); writer.writerows(curves)
    with open(args.out_dir / "fs_wavevector_scan_summary.csv", "w", newline="") as handle:
        fields = ["label", "n", "k_inv_A", "lambda_A", "dt_ps", "n_frames", "n_water", "beta_kww", "tau_kww_ps", "r2_log", "n_fit"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(summary)
    with open(args.out_dir / "summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)


if __name__ == "__main__":
    main()
