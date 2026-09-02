#!/usr/bin/env python3
"""Rebuild longitudinal total/self/distinct oxygen ISF from LAMMPS dumps.

The estimator is ``F(k,t)=<rho_k(t0+t) rho_k*(t0)>/N`` and
``Fs(k,t)=<sum_j exp[i k (z_j(t0+t)-z_j(t0))]>/N``.  The distinct term is
their difference.  It uses every origin within a declared analysis window and
keeps protocol groups separate; it never averages different dump protocols.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_positions(path: Path, max_frames: int) -> tuple[np.ndarray, float, np.ndarray]:
    """Read unwrapped z positions from an oxygen-only/full-kinematics dump."""
    frames: list[np.ndarray] = []
    first_ids: np.ndarray | None = None
    lz: float | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        while len(frames) < max_frames:
            if handle.readline() == "":
                break
            timestep = handle.readline()
            if not timestep:
                break
            if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Unexpected dump structure in {path}")
            n_atoms = int(handle.readline())
            if not handle.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Missing BOX BOUNDS in {path}")
            bounds = [handle.readline().split() for _ in range(3)]
            frame_lz = float(bounds[2][1]) - float(bounds[2][0])
            if lz is None:
                lz = frame_lz
            elif not np.isclose(lz, frame_lz, rtol=0, atol=1e-6):
                raise ValueError(f"Lz changed within {path}: {lz} vs {frame_lz}")
            atom_header = handle.readline().split()[2:]
            need = {"id", "z"}
            if not need.issubset(atom_header):
                raise ValueError(f"{path} lacks required fields {need}; has {atom_header}")
            idx = {name: atom_header.index(name) for name in atom_header}
            rows = [handle.readline() for _ in range(n_atoms)]
            values = np.fromstring(" ".join(rows), sep=" ", dtype=np.float64)
            if values.size != n_atoms * len(atom_header):
                raise ValueError(f"Could not parse an atom block in {path}")
            values = values.reshape(n_atoms, len(atom_header))
            ids = values[:, idx["id"]].astype(np.int64)
            order = np.argsort(ids)
            ids = ids[order]
            if first_ids is None:
                first_ids = ids
            elif not np.array_equal(first_ids, ids):
                raise ValueError(f"Atom identity changed in {path}")
            z = values[order, idx["z"]]
            if "iz" in idx:
                z = z + values[order, idx["iz"]] * lz
            frames.append(z.astype(np.float32))
    if not frames or lz is None or first_ids is None:
        raise ValueError(f"No readable frames in {path}")
    return np.asarray(frames), lz, first_ids


def autocorrelation_complex(series: np.ndarray, max_lag: int) -> np.ndarray:
    n = series.shape[0]
    size = 1 << (2 * n - 1).bit_length()
    fft = np.fft.fft(series, size, axis=0)
    corr = np.fft.ifft(fft * np.conjugate(fft), axis=0).real[: max_lag + 1]
    return corr


def isf_from_z(z: np.ndarray, k: float, max_lag: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_time, n_atoms = z.shape
    max_lag = min(max_lag, n_time - 1)
    phase = np.exp(1j * k * z, dtype=np.complex64)
    rho = phase.sum(axis=1)
    total = autocorrelation_complex(rho[:, None], max_lag)[:, 0] / (n_atoms * (n_time - np.arange(max_lag + 1)))
    self_sum = np.zeros(max_lag + 1, dtype=np.float64)
    # Chunking avoids an FFT workspace proportional to all particles.
    for start in range(0, n_atoms, 48):
        stop = min(n_atoms, start + 48)
        self_sum += autocorrelation_complex(phase[:, start:stop], max_lag).sum(axis=1)
    self_term = self_sum / (n_atoms * (n_time - np.arange(max_lag + 1)))
    return total, self_term, total - self_term


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path, help="JSON list: label, path, length_L, protocol")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dt-ps", type=float, default=0.01)
    parser.add_argument("--window-ps", type=float, default=200.0, help="read this leading window per dump")
    parser.add_argument("--max-lag-ps", type=float, default=100.0)
    parser.add_argument("--matched-k-inv-A", type=float, default=0.06230846323)
    args = parser.parse_args()
    sources = json.loads(args.manifest.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_dir = args.output_dir / "per_replica"
    per_dir.mkdir(exist_ok=True)
    max_frames = round(args.window_ps / args.dt_ps) + 1
    max_lag = round(args.max_lag_ps / args.dt_ps)
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    all_metadata = []

    for source in sources:
        path = Path(source["path"])
        z, lz, atom_ids = read_positions(path, max_frames)
        modes = {"kmin": 1, "matched_k": max(1, round(args.matched_k_inv_A * lz / (2 * np.pi)))}
        record = {**source, "frames_used": int(z.shape[0]), "n_oxygen": int(z.shape[1]), "lz_A": lz,
                  "fields_required": "id z [iz]", "modes": {}}
        for mode_name, n in modes.items():
            k = 2 * np.pi * n / lz
            total, self_term, distinct = isf_from_z(z, k, max_lag)
            time = np.arange(total.size) * args.dt_ps
            output = per_dir / f"{source['label']}_{mode_name}.npz"
            np.savez_compressed(output, time_ps=time, F_total=total, F_self=self_term,
                                F_distinct=distinct, k_inv_A=k, n=n, lz_A=lz,
                                source_path=str(path), atom_ids=atom_ids,
                                n_time_origins=z.shape[0] - np.arange(total.size))
            record["modes"][mode_name] = {"n": n, "k_inv_A": k, "output": str(output.relative_to(args.output_dir))}
            grouped[(str(source["length_L"]), source["protocol"], mode_name)].append(
                {"time": time, "total": total, "self": self_term, "distinct": distinct, "k": k, "label": source["label"]})
        all_metadata.append(record)
        print(f"completed {source['label']}: {z.shape[0]} frames, {z.shape[1]} O, Lz={lz:.6f}")

    ensemble = []
    for (length, protocol, mode_name), records in sorted(grouped.items()):
        limit = min(len(item["time"]) for item in records)
        time = records[0]["time"][:limit]
        summary = {"length_L": length, "protocol": protocol, "mode": mode_name,
                   "n_replicas": len(records), "k_inv_A_mean": float(np.mean([r["k"] for r in records])),
                   "replicas": [r["label"] for r in records]}
        payload = {"time_ps": time,
                   "n_time_origins_per_replica": records[0]["time"].size - np.arange(limit)}
        for key, output_name in (("total", "F_total"), ("self", "F_self"), ("distinct", "F_distinct")):
            array = np.asarray([r[key][:limit] for r in records])
            payload[f"{output_name}_mean"] = array.mean(axis=0)
            payload[f"{output_name}_sem"] = array.std(axis=0, ddof=1) / np.sqrt(len(records)) if len(records) > 1 else np.full(limit, np.nan)
        out = args.output_dir / f"ISF_{length}_{mode_name}_mean_sem.npz"
        np.savez_compressed(out, **payload, **summary)
        summary["output"] = out.name
        ensemble.append(summary)
    (args.output_dir / "source_manifest_resolved.json").write_text(json.dumps(all_metadata, indent=2), encoding="utf-8")
    (args.output_dir / "ensemble_manifest.json").write_text(json.dumps(ensemble, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
