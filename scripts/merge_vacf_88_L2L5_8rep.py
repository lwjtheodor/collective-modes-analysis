"""Build the protocol-qualified 8-replica peculiar axial VACF archive.

The reference rep1--3 curves are stored as unnormalised all-origin Cvv and
the remotely analysed rep4--8 curves are also all-origin, C(0)-normalised
VACFs.  We merge the common observable: a per-replica normalised,
instantaneous oxygen-COM-subtracted axial VACF, with uncertainty from the
eight replica means.  The remote lab-frame curves stay provenance-only because
no comparable lab-frame rep1--3 archive exists.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "results/collective_mode_response/vacf_alpha_10fs_L2L10_unified_2026-08-20"
NEW = ROOT / "results/collective_mode_response/vacf_analysis_allorigins_20260821/outputs"
OUT = ROOT / "results/collective_mode_response/vacf_88_L2L5_10fs_1ns_8rep_weakNH_nomom/2026-08-21/allorigins"


def read_old(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    lag = np.asarray(data["lag_ps"], float)
    curves = np.column_stack([np.asarray(data[f"rep{i}"], float) for i in range(1, 4)])
    return lag, curves / curves[0:1, :]


def read_new(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return np.asarray(data["lag_ps"], float), np.asarray(data["vacf_peculiar_mean"], float)


def write_csv(path: Path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle); writer.writerow(header); writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    per = OUT / "per_replica"; per.mkdir(exist_ok=True)
    manifest = {"observable": "normalised axial oxygen peculiar VACF",
                "velocity_definition": "instantaneous oxygen-COM-subtracted axial velocity",
                "units": "dimensionless Cvv(t)/Cvv(0)",
                "simulation_protocol": "(8,8), weak NH, no momentum-removal, 10 fs dump cadence, 1 ns trajectory",
                "analysis_window": "0--100 ps; common per-replica normalisation; aggregate SEM across 8 replicas",
                "lengths": {}}
    for label in ("2L", "3L", "4L", "5L"):
        old_lag, old = read_old(OLD / label / "cvv_per_replica.csv")
        use = old_lag <= 100.0 + 1e-12
        lag = old_lag[use]
        curves = [old[use, i] for i in range(3)]
        sources = []
        for rep in range(1, 4):
            sources.append({"replica": rep, "kind": "local all-origin Cvv, subsequently C(0)-normalised",
                            "path": str((OLD / label / "cvv_per_replica.csv").resolve())})
        for rep in range(4, 9):
            p = NEW / f"VACF_8_8_L{label[0]}_extra_rep{rep}_1ns_10fs.csv"
            new_lag, curve = read_new(p)
            if not np.array_equal(new_lag, lag):
                raise ValueError(f"cadence/lag mismatch: {p}")
            if not np.isclose(curve[0], 1.0):
                raise ValueError(f"non-unit zero lag: {p}")
            curves.append(curve)
            meta = json.loads(p.with_suffix(".json").read_text(encoding="utf-8"))
            if meta["dt_ps"] != 0.01 or meta["max_lag_ps"] != 100.0 or meta["n_frames"] != 100001:
                raise ValueError(f"remote metadata gate failed: {p}")
            if meta["n_blocks"] != 0 or meta.get("estimator") != "all-origin full trajectory":
                raise ValueError(f"remote all-origin estimator gate failed: {p}")
            sources.append({"replica": rep, "kind": "remote all-origin full-trajectory VACF", "path": str(p.resolve()),
                            "metadata": str(p.with_suffix('.json').resolve()), "n_blocks": 0,
                            "estimator": "all-origin full trajectory"})
        values = np.column_stack(curves)
        mean = values.mean(axis=1)
        sem = values.std(axis=1, ddof=1) / np.sqrt(values.shape[1])
        write_csv(per / f"VACF_8_8_{label}_peculiar_per_replica_normalised.csv",
                  ["lag_ps"] + [f"rep{i}" for i in range(1, 9)],
                  np.column_stack([lag, values]))
        write_csv(OUT / f"VACF_8_8_{label}_peculiar_mean_sem_8rep.csv",
                  ["lag_ps", "vacf_peculiar_mean", "vacf_peculiar_replica_sem", "n_replicas"],
                  np.column_stack([lag, mean, sem, np.full_like(lag, 8)]))
        manifest["lengths"][label] = {"n_replicas": 8, "lag_points": int(len(lag)),
                                      "sources": sources,
                                      "qa": {"zero_lag_mean": float(mean[0]), "zero_lag_sem": float(sem[0]),
                                             "all_finite": bool(np.isfinite(values).all())}}
    (OUT / "source_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (OUT / "README.md").write_text("""# (8,8) L2--L5 10-fs / 1-ns eight-replica axial peculiar VACF

This is the authoritative, protocol-qualified merge of weak-NH / no-momentum-removal `(8,8)` trajectories: L2--L5, 10-fs dump cadence, 1 ns, reps 1--8.  Observable is the axial oxygen self VACF after subtracting the instantaneous oxygen COM velocity in every frame, normalised individually at zero lag. Every replica uses the full-trajectory all-origin estimator.

`VACF_8_8_<L>_peculiar_mean_sem_8rep.csv` reports the arithmetic mean and replica-to-replica SEM (n=8). Reps 1--3 were recalibrated from archived all-origin unnormalised Cvv; reps 4--8 are all-origin full-trajectory estimates from compact remote oxygen dumps. This is valid for the common normalised peculiar VACF only. The lab-frame remote curves remain under the fetched remote result package and are not represented as an eight-replica observable.

All sources, cadence gates, and the different per-rep estimator provenance are in `source_manifest.json`; the merge logic is `scripts/merge_vacf_88_L2L5_8rep.py`.
""", encoding="utf-8")


if __name__ == "__main__":
    main()
