#!/usr/bin/env python3
"""Streaming collective axial ISF for a user-selected integer k series.

Designed for full water dumps: it keeps only rho_k(t), not coordinates, so a
10 ns / 100 fs trajectory can be analysed without creating a second huge copy.
For oxygen atoms, use ``--oxygen-type 3``.  Results retain raw F(k,t), its
normalization F(k,t)/F(k,0), and the exact time-origin count at every lag.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def read_rho(path: Path, nmax: int, oxygen_type: int, max_frames: int | None):
    rho, lz, n_oxygen = [], None, None
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        while max_frames is None or len(rho) < max_frames:
            if not fh.readline():
                break
            if not fh.readline():
                break
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Unexpected dump structure: {path}")
            n_atoms = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Missing box bounds: {path}")
            bounds = [fh.readline().split() for _ in range(3)]
            here_lz = float(bounds[2][1]) - float(bounds[2][0])
            if lz is None:
                lz = here_lz
            elif not np.isclose(lz, here_lz, rtol=0, atol=1e-6):
                raise ValueError("The fixed-box ISF estimator requires constant Lz")
            header = fh.readline().split()[2:]
            idx = {name: header.index(name) for name in header}
            missing = {"type", "z"} - idx.keys()
            if missing:
                raise ValueError(f"Missing {missing} in {path}")
            rows = [fh.readline() for _ in range(n_atoms)]
            values = np.fromstring(" ".join(rows), sep=" ", dtype=np.float64)
            values = values.reshape(n_atoms, len(header))
            oxygen = values[:, idx["type"]].astype(np.int16) == oxygen_type
            if n_oxygen is None:
                n_oxygen = int(oxygen.sum())
            elif n_oxygen != int(oxygen.sum()):
                raise ValueError("Oxygen count changes across frames")
            z = values[oxygen, idx["z"]]
            if "iz" in idx:
                z = z + values[oxygen, idx["iz"]] * lz
            # k_n=n k_1: one transcendental evaluation plus recurrence is
            # mathematically identical to evaluating ten separate exponentials.
            phase1 = np.exp(1j * (2 * np.pi / lz) * z).astype(np.complex64)
            phases = np.empty((z.size, nmax), dtype=np.complex64)
            phases[:, 0] = phase1
            for mode in range(1, nmax):
                phases[:, mode] = phases[:, mode - 1] * phase1
            rho.append(phases.sum(axis=0, dtype=np.complex64))
    if not rho:
        raise ValueError(f"No frames read from {path}")
    return np.asarray(rho), float(lz), int(n_oxygen)


def correlation(rho: np.ndarray, max_lag: int, n_oxygen: int):
    ntime = rho.shape[0]
    max_lag = min(max_lag, ntime - 1)
    fft_len = 1 << (2 * ntime - 1).bit_length()
    transformed = np.fft.fft(rho, fft_len, axis=0)
    corr = np.fft.ifft(transformed * np.conjugate(transformed), axis=0).real[:max_lag + 1]
    origins = ntime - np.arange(max_lag + 1)
    return corr / (n_oxygen * origins[:, None])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest", type=Path, help="JSON records containing label and path")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--dt-ps", type=float, required=True)
    ap.add_argument("--max-lag-ps", type=float, required=True)
    ap.add_argument("--nmax", type=int, default=10)
    ap.add_argument("--oxygen-type", type=int, default=3)
    args = ap.parse_args()
    records = json.loads(args.manifest.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_dir = args.output_dir / "per_replica"
    per_dir.mkdir(exist_ok=True)
    requested_lag = round(args.max_lag_ps / args.dt_ps)
    all_curves, provenance = [], []
    for record in records:
        rho, lz, no = read_rho(Path(record["path"]), args.nmax, args.oxygen_type, None)
        exact_k = np.arange(1, args.nmax + 1) * 2*np.pi / lz
        curve = correlation(rho, requested_lag, no)
        time = np.arange(curve.shape[0]) * args.dt_ps
        out = per_dir / f"{record['label']}_collective_k1-k{args.nmax}.npz"
        np.savez_compressed(out, time_ps=time, F_total=curve, F_total_normalized=curve / curve[0],
                            n=np.arange(1, args.nmax + 1), k_inv_A=exact_k,
                            n_time_origins=rho.shape[0] - np.arange(curve.shape[0]),
                            lz_A=lz, n_oxygen=no, source_path=record["path"])
        all_curves.append(curve)
        provenance.append({**record, "frames_used": int(rho.shape[0]), "lz_A": lz,
                           "n_oxygen": no, "output": str(out.relative_to(args.output_dir))})
        print(f"completed {record['label']}: {rho.shape[0]} frames, {no} O, Lz={lz:.6f}", flush=True)

    stack = np.asarray(all_curves)
    mean = stack.mean(axis=0)
    sem = stack.std(axis=0, ddof=1) / np.sqrt(stack.shape[0])
    time = np.arange(mean.shape[0]) * args.dt_ps
    k = np.arange(1, args.nmax + 1) * 2 * np.pi / provenance[0]["lz_A"]
    np.savez_compressed(args.output_dir / f"ISF_88_L10_k1-k{args.nmax}_mean_sem.npz",
                        time_ps=time, F_total_mean=mean, F_total_sem=sem,
                        F_total_normalized_mean=mean / mean[0],
                        F_total_normalized_sem=sem / np.abs(mean[0]),
                        n=np.arange(1, args.nmax + 1), k_inv_A=k,
                        n_time_origins_per_replica=provenance[0]["frames_used"] - np.arange(mean.shape[0]),
                        n_replicas=stack.shape[0], protocol="rethermalized_weakNH_nomom_100fs_10ns")
    (args.output_dir / "source_manifest_resolved.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
