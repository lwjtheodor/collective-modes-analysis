"""Audit the oxygen-proxy orbital angular momentum in full-velocity dumps.

These dumps contain oxygen sites only, not H coordinates/velocities.  The
result is therefore deliberately labelled an O-proxy: it is the exact group
quantity associated with the existing oxygen TA_theta current, but not the
strict all-atom molecular angular momentum.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def audit_dump(path: Path, stride: int):
    rows, frame = [], 0
    with path.open("r", encoding="utf-8") as fh:
        while True:
            tag = fh.readline()
            if not tag:
                break
            if tag.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Unexpected record in {path}: {tag!r}")
            step = int(fh.readline())
            fh.readline(); natoms = int(fh.readline())
            fh.readline(); [fh.readline() for _ in range(3)]
            header = fh.readline().split()[2:]
            wanted = ["x", "y", "vx", "vy"]
            if any(q not in header for q in wanted):
                raise ValueError(f"Missing {wanted} in {path}")
            if frame % stride:
                for _ in range(natoms): fh.readline()
            else:
                vals = np.fromstring(" ".join(fh.readline() for _ in range(natoms)), sep=" ")
                a = vals.reshape(natoms, len(header))
                c = {q: i for i, q in enumerate(header)}
                x, y, vx, vy = (a[:, c[q]] for q in ("x", "y", "vx", "vy"))
                # Unit mass version first; oxygen mass is a common scale factor.
                axis = np.sum(x * vy - y * vx)
                xx, yy, vxx, vyy = x-x.mean(), y-y.mean(), vx-vx.mean(), vy-vy.mean()
                internal = np.sum(xx * vyy - yy * vxx)
                rows.append({"step":step, "time_ps":step*0.0005,
                             "Lz_axis_sum_A2_per_fs":axis,
                             "Lz_internal_sum_A2_per_fs":internal,
                             "Lz_axis_amu_A2_per_fs":15.999*axis,
                             "Lz_internal_amu_A2_per_fs":15.999*internal,
                             "x_COM_A":x.mean(), "y_COM_A":y.mean(),
                             "vx_COM_A_per_fs":vx.mean(), "vy_COM_A_per_fs":vy.mean()})
            frame += 1
    return pd.DataFrame(rows), frame, natoms


def trend(y, t):
    slope, intercept = np.polyfit(t, y, 1)
    fit = intercept + slope*t
    resid = y-fit
    se = np.sqrt(np.sum(resid**2)/(len(t)-2) / np.sum((t-t.mean())**2))
    sd = np.std(y, ddof=1)
    return slope, se, sd, (y[-1]-y[0])/sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dumps", nargs="+", required=True, type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--stride", type=int, default=10, help="Frame decimation; native frame is 0.01 ps")
    args = ap.parse_args()
    dd, figdir = args.output / "derived_data", args.output / "figures"
    dd.mkdir(parents=True, exist_ok=True); figdir.mkdir(exist_ok=True)
    all_rows, summary = [], []
    for path in args.dumps:
        case = re.sub(r"_oxygen_xyz_vxyz_10fs_1ns\.dump$", "", path.name)
        df, nframes, natoms = audit_dump(path, args.stride)
        df.insert(0, "case", case); all_rows.append(df)
        for col in ("Lz_axis_amu_A2_per_fs", "Lz_internal_amu_A2_per_fs"):
            sl, se, sd, delta = trend(df[col].to_numpy(), df.time_ps.to_numpy())
            summary.append({"case":case, "definition":col, "n_native_frames":nframes,
                            "n_sampled_frames":len(df), "n_oxygen":natoms,
                            "mean":df[col].mean(), "std":sd, "mean_over_std":df[col].mean()/sd,
                            "linear_slope_per_ps":sl, "slope_SE_per_ps":se,
                            "total_endpoint_change_over_std":delta,
                            "range_over_std":(df[col].max()-df[col].min())/sd,
                            "source_dump":str(path)})
    series = pd.concat(all_rows, ignore_index=True)
    sm = pd.DataFrame(summary)
    series.to_csv(dd / "oxygen_proxy_orbital_Lz_timeseries.csv", index=False)
    sm.to_csv(dd / "oxygen_proxy_orbital_Lz_summary.csv", index=False)
    cases = series.case.unique()
    fig, axes = plt.subplots(len(cases), 2, figsize=(13, 3.1*len(cases)), squeeze=False)
    for r, case in enumerate(cases):
        q = series[series.case == case]
        for c, (col, label) in enumerate((("Lz_axis_amu_A2_per_fs", "about CNT z axis"),
                                           ("Lz_internal_amu_A2_per_fs", "O-COM frame"))):
            ax = axes[r,c]
            ax.plot(q.time_ps, q[col], lw=.5, color="#2166ac")
            ax.axhline(0, color="k", lw=.7)
            ax.set(title=f"{case}: {label}", xlabel="time (ps)", ylabel=r"$L_z$ (amu A$^2$ fs$^{-1}$)")
            ax.grid(alpha=.22)
    fig.suptitle("Explicit CNT full-velocity cases: oxygen-proxy orbital angular momentum")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(figdir / f"oxygen_proxy_orbital_Lz_timeseries.{ext}", dpi=300)
    plt.close(fig)
    metadata = {"scope":"oxygen-only dumps; not strict all-atom water angular momentum",
                "angular_momentum_axis":"sum_i mO (x_i vy_i-y_i vx_i)",
                "internal_control":"same expression after instantaneous O-COM position/velocity subtraction",
                "oxygen_mass_amu":15.999, "native_dt_ps":0.01, "sampling_stride":args.stride,
                "sampled_dt_ps":0.01*args.stride, "source_dumps":[str(p) for p in args.dumps]}
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (args.output / "FINISHED.txt").write_text("Oxygen-proxy orbital Lz audit finished successfully.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
