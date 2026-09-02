"""Compute axial VACF tails from LAMMPS oxygen/water dumps.

Outputs both the laboratory-frame VACF and the instantaneous water-COM
(``peculiar``) VACF.  The latter is essential when a trajectory was run
without a momentum-removal fix: a residual box drift otherwise appears as a
long positive VACF plateau.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def read_vz(path: Path, selection_type: int):
    """Read one selected atom type and retain dump timesteps for cadence checks."""
    frames, mol_ref, steps = [], None, []
    with path.open("r", errors="replace") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(handle.readline())
            if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("missing NUMBER OF ATOMS")
            natoms = int(handle.readline())
            if not handle.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError("missing BOX BOUNDS")
            for _ in range(3):
                handle.readline()
            header = handle.readline().split()[2:]
            index = {name: i for i, name in enumerate(header)}
            if "vz" not in index:
                raise ValueError("dump needs vz")
            # Full-water dumps use mol/type; compact oxygen-only dumps retain
            # stable atom IDs and require no further type selection.
            identity = "mol" if "mol" in index else "id"
            current = {}
            for _ in range(natoms):
                row = handle.readline().split()
                selected = "type" not in index or int(float(row[index["type"]])) == selection_type
                if selected:
                    current[int(float(row[index[identity]]))] = float(row[index["vz"]])
            mols = np.asarray(sorted(current), dtype=np.int64)
            if mol_ref is None:
                mol_ref = mols
            elif not np.array_equal(mol_ref, mols):
                raise ValueError("oxygen molecule IDs changed")
            frames.append([current[int(m)] for m in mol_ref])
            steps.append(step)
    if not frames:
        raise ValueError("no atoms selected from dump")
    return np.asarray(steps, dtype=np.int64), np.asarray(frames, dtype=float)


def acf(x, max_lag):
    """All-origin column-mean ACF, normalized at zero lag."""
    nframe, nwater = x.shape
    size = 1 << (2 * nframe - 1).bit_length()
    corr = np.fft.irfft(np.abs(np.fft.rfft(x, n=size, axis=0)) ** 2,
                        n=size, axis=0)[: max_lag + 1].sum(axis=1)
    corr /= np.arange(nframe, nframe - max_lag - 1, -1) * nwater
    return corr / corr[0]


def block_curves(vz, max_lag, nblocks):
    # A block must contain the full requested tail.  Equal frame-count blocks
    # prevent a longer trajectory from silently carrying more weight.
    block_size = len(vz) // nblocks
    if block_size <= max_lag:
        raise ValueError("blocks are too short for requested max lag")
    curves = {"lab": [], "peculiar": []}
    for first in range(0, block_size * nblocks, block_size):
        x = vz[first:first + block_size]
        curves["lab"].append(acf(x, max_lag))
        curves["peculiar"].append(acf(x - x.mean(axis=1, keepdims=True), max_lag))
    return {key: np.asarray(value) for key, value in curves.items()}, block_size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, type=Path)
    ap.add_argument("--dt-ps", type=float,
                    help="legacy explicit dump-frame interval in ps")
    ap.add_argument("--timestep-ps", type=float,
                    help="LAMMPS integration timestep in ps; enables cadence validation")
    ap.add_argument("--selection-type", type=int, default=3,
                    help="LAMMPS atom type to analyse (default: 3, explicit-water oxygen)")
    ap.add_argument("--max-lag-ps", type=float, default=100.0)
    ap.add_argument("--discard-initial-ps", type=float, default=0.0,
                    help="exclude this initial production time before forming every time origin")
    ap.add_argument("--stop-after-ps", type=float,
                    help="optional production-time endpoint after which frames are excluded")
    ap.add_argument("--nblocks", type=int, default=0,
                    help="positive: equal-time block diagnostic; 0: one all-origin full-trajectory estimate")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--case-id", required=True)
    args = ap.parse_args()
    steps, vz = read_vz(args.dump, args.selection_type)
    dump_step_interval = int(np.median(np.diff(steps)))
    if not np.all(np.diff(steps) == dump_step_interval):
        raise ValueError("dump timesteps are not uniformly spaced")
    if args.timestep_ps is not None:
        if args.timestep_ps <= 0:
            raise ValueError("--timestep-ps must be positive")
        dt_ps = dump_step_interval * args.timestep_ps
        if args.dt_ps is not None and not np.isclose(args.dt_ps, dt_ps):
            raise ValueError("--dt-ps conflicts with cadence derived from --timestep-ps")
    elif args.dt_ps is not None and args.dt_ps > 0:
        dt_ps = args.dt_ps
    else:
        raise ValueError("supply --timestep-ps or a positive legacy --dt-ps")
    if args.discard_initial_ps < 0:
        raise ValueError("--discard-initial-ps must be non-negative")
    discard_frames = int(round(args.discard_initial_ps / dt_ps))
    if not np.isclose(discard_frames * dt_ps, args.discard_initial_ps):
        raise ValueError("discard duration must lie on the dump-frame grid")
    if discard_frames:
        steps, vz = steps[discard_frames:], vz[discard_frames:]
    if args.stop_after_ps is not None:
        if args.stop_after_ps <= args.discard_initial_ps:
            raise ValueError("--stop-after-ps must exceed --discard-initial-ps")
        keep_frames = int(round((args.stop_after_ps - args.discard_initial_ps) / dt_ps)) + 1
        if not np.isclose((keep_frames - 1) * dt_ps, args.stop_after_ps - args.discard_initial_ps):
            raise ValueError("stop duration must lie on the dump-frame grid")
        steps, vz = steps[:keep_frames], vz[:keep_frames]
    max_lag = int(round(args.max_lag_ps / dt_ps))
    if len(vz) <= max_lag:
        raise ValueError("trajectory is shorter than the requested tail")
    if args.nblocks < 0:
        raise ValueError("nblocks must be zero (all origins) or positive")
    if args.nblocks == 0:
        curves = {"lab": np.asarray([acf(vz, max_lag)]),
                  "peculiar": np.asarray([acf(vz - vz.mean(axis=1, keepdims=True), max_lag)])}
        block_size = len(vz)
        estimator = "all-origin full trajectory"
    else:
        curves, block_size = block_curves(vz, max_lag, args.nblocks)
        estimator = "equal-time block mean"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    lag = np.arange(max_lag + 1) * dt_ps
    with args.out.open("w", newline="") as handle:
        fields = ["case_id", "lag_ps", "vacf_lab_mean", "vacf_lab_sem",
                  "vacf_peculiar_mean", "vacf_peculiar_sem", "n_blocks"]
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for i, time in enumerate(lag):
            writer.writerow({"case_id": args.case_id, "lag_ps": time,
                             "vacf_lab_mean": curves["lab"][:, i].mean(),
                             "vacf_lab_sem": (curves["lab"][:, i].std(ddof=1) / np.sqrt(args.nblocks)
                                              if args.nblocks > 1 else float("nan")),
                             "vacf_peculiar_mean": curves["peculiar"][:, i].mean(),
                             "vacf_peculiar_sem": (curves["peculiar"][:, i].std(ddof=1) / np.sqrt(args.nblocks)
                                                   if args.nblocks > 1 else float("nan")),
                             "n_blocks": args.nblocks})
    meta = {"case_id": args.case_id, "n_frames": int(len(vz)), "n_water": int(vz.shape[1]),
            "selection_type": args.selection_type, "discard_initial_ps": args.discard_initial_ps,
            "discarded_frames": discard_frames, "stop_after_ps": args.stop_after_ps,
            "first_step": int(steps[0]),
            "last_step": int(steps[-1]), "dump_step_interval": dump_step_interval,
            "dt_ps": dt_ps, "max_lag_ps": args.max_lag_ps,
            "n_blocks": args.nblocks, "block_duration_ps": block_size * dt_ps,
            "estimator": estimator,
            "velocity_frame": "lab and instantaneous selected-atom-COM-subtracted"}
    args.out.with_suffix(".json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
