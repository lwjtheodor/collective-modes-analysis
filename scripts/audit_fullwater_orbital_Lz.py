"""Mass-weighted water-group angular-momentum audit for full-water LAMMPS dumps."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_MASS = {2: 1.008, 3: 15.999}  # Explicit-CNT SPC/E production family


def analyze_one(path: Path, stride: int, masses: dict[int, float]):
    out, iframe = [], 0
    with path.open("r", encoding="utf-8") as fh:
        while True:
            marker = fh.readline()
            if not marker:
                break
            if marker.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Malformed dump marker: {marker!r}")
            step = int(fh.readline())
            fh.readline(); n_atoms = int(fh.readline())
            fh.readline(); [fh.readline() for _ in range(3)]
            header = fh.readline().split()[2:]
            need = ("mol", "type", "x", "y", "vx", "vy", "vz")
            if any(q not in header for q in need):
                raise ValueError(f"Missing full-water fields in {path}: {header}")
            if iframe % stride:
                for _ in range(n_atoms):
                    fh.readline()
            else:
                vals = np.fromstring(" ".join(fh.readline() for _ in range(n_atoms)), sep=" ")
                a = vals.reshape(n_atoms, len(header))
                c = {q: i for i, q in enumerate(header)}
                mol, typ = a[:, c["mol"]].astype(int), a[:, c["type"]].astype(int)
                unknown = set(np.unique(typ)) - set(masses)
                if unknown:
                    raise ValueError(f"Unexpected atom types {unknown}; dump is not water-only")
                m = np.vectorize(masses.__getitem__)(typ)
                x, y, vx, vy, vz = (a[:, c[q]] for q in ("x", "y", "vx", "vy", "vz"))
                # Exact all-atom water-group Lz about the fixed CNT z axis.
                lz_total = np.sum(m * (x * vy - y * vx))
                pz_total = np.sum(m * vz)
                # Separately retain the molecular-COM orbital contribution.
                n_mol = mol.max() + 1
                M = np.bincount(mol, weights=m, minlength=n_mol)
                present = M > 0
                X = np.bincount(mol, weights=m*x, minlength=n_mol)[present] / M[present]
                Y = np.bincount(mol, weights=m*y, minlength=n_mol)[present] / M[present]
                VX = np.bincount(mol, weights=m*vx, minlength=n_mol)[present] / M[present]
                VY = np.bincount(mol, weights=m*vy, minlength=n_mol)[present] / M[present]
                lz_com = np.sum(M[present] * (X * VY - Y * VX))
                # Axis-origin result can contain rigid transverse translation;
                # this control removes that collective translation explicitly.
                xc, yc = np.sum(m*x)/m.sum(), np.sum(m*y)/m.sum()
                vxc, vyc = np.sum(m*vx)/m.sum(), np.sum(m*vy)/m.sum()
                lz_comframe = np.sum(m*((x-xc)*(vy-vyc)-(y-yc)*(vx-vxc)))
                out.append({"step":step, "time_ps":step*0.0005,
                            "Lz_total_water_axis_amu_A2_fs":lz_total,
                            "Lz_molecular_COM_orbital_axis_amu_A2_fs":lz_com,
                            "Lz_total_water_COMframe_amu_A2_fs":lz_comframe,
                            "Pz_total_water_amu_A_fs":pz_total,
                            "water_x_COM_A":xc, "water_y_COM_A":yc,
                            "water_vx_COM_A_fs":vxc, "water_vy_COM_A_fs":vyc})
            iframe += 1
    return pd.DataFrame(out), iframe, n_atoms


def metrics(y: np.ndarray, t: np.ndarray):
    slope, intercept = np.polyfit(t, y, 1)
    residual = y - (intercept + slope*t)
    slope_se = np.sqrt(np.sum(residual**2)/(len(y)-2) / np.sum((t-t.mean())**2))
    sd = np.std(y, ddof=1)
    return {"mean":float(y.mean()), "std":float(sd), "mean_over_std":float(y.mean()/sd),
            "linear_slope_per_ps":float(slope), "linear_slope_SE_per_ps":float(slope_se),
            "endpoint_change_over_std":float((y[-1]-y[0])/sd),
            "range_over_std":float((y.max()-y.min())/sd)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dumps", nargs="+", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--stride", type=int, default=10, help="Native 0.1 ps frame decimation")
    ap.add_argument("--oxygen-type", type=int, default=3,
                    help="O type (explicit SPC/E default: 3; implicit family: 1)")
    ap.add_argument("--hydrogen-type", type=int, default=2,
                    help="H type (default: 2)")
    args = ap.parse_args()
    if args.oxygen_type == args.hydrogen_type:
        raise ValueError("oxygen and hydrogen types must differ")
    masses = {args.hydrogen_type: 1.008, args.oxygen_type: 15.999}
    dd, figs = args.output/"derived_data", args.output/"figures"
    dd.mkdir(parents=True, exist_ok=True); figs.mkdir(exist_ok=True)
    series, rows = [], []
    for source in args.dumps:
        case = re.sub(r"\.water_100fs_10ns\.dump$", "", source.name)
        df, n_native, n_atoms = analyze_one(source, args.stride, masses)
        df.insert(0, "case", case); series.append(df)
        for col in ("Lz_total_water_axis_amu_A2_fs", "Lz_molecular_COM_orbital_axis_amu_A2_fs", "Lz_total_water_COMframe_amu_A2_fs", "Pz_total_water_amu_A_fs"):
            q = metrics(df[col].to_numpy(), df.time_ps.to_numpy())
            rows.append({"case":case, "quantity":col, "n_atoms":n_atoms,
                         "n_native_frames":n_native, "n_sampled_frames":len(df),
                         "source_dump":str(source), **q})
    data, summary = pd.concat(series, ignore_index=True), pd.DataFrame(rows)
    data.to_csv(dd/"fullwater_angular_momentum_timeseries.csv", index=False)
    summary.to_csv(dd/"fullwater_angular_momentum_summary.csv", index=False)
    cases = data.case.unique()
    fig, axes = plt.subplots(len(cases), 2, figsize=(13, 2.9*len(cases)), squeeze=False)
    for i, case in enumerate(cases):
        q = data[data.case == case]
        for j, (col, label) in enumerate((("Lz_total_water_axis_amu_A2_fs", "all water, CNT-axis origin"),
                                           ("Lz_molecular_COM_orbital_axis_amu_A2_fs", "molecular-COM orbital part"))):
            ax = axes[i,j]
            ax.plot(q.time_ps, q[col], color="#2166ac", lw=.45)
            ax.axhline(0, color="black", lw=.65)
            ax.grid(alpha=.22)
            ax.set(title=f"{case}: {label}", xlabel="time (ps)", ylabel=r"$L_z$ (amu A$^2$ fs$^{-1}$)")
    fig.suptitle("Explicit CNT: exact mass-weighted water-group angular momentum")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(figs/f"fullwater_angular_momentum_timeseries.{ext}", dpi=300)
    plt.close(fig)
    fig, axes = plt.subplots(len(cases), 1, figsize=(8.2, 2.6*len(cases)), squeeze=False)
    for i, case in enumerate(cases):
        q = data[data.case == case]
        ax = axes[i, 0]
        ax.plot(q.time_ps, q.Pz_total_water_amu_A_fs, color="#2166ac", lw=.45)
        ax.axhline(0, color="black", lw=.65); ax.grid(alpha=.22)
        ax.set(title=f"{case}: total water axial momentum", xlabel="time (ps)", ylabel=r"$P_z$ (amu A fs$^{-1}$)")
    fig.suptitle("Explicit CNT: mass-weighted water axial momentum")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(figs/f"fullwater_axial_momentum_timeseries.{ext}", dpi=300)
    plt.close(fig)
    meta = {"definition_total":"sum over all water atoms m_a(x_a v_ya-y_a v_xa)",
            "definition_molecular_COM":"sum over molecules M(Rx Vy-Ry Vx)",
            "definition_Pz":"sum over all water atoms m_a vz_a",
            "axis":"fixed CNT z axis", "masses_amu":masses,
            "oxygen_type":args.oxygen_type, "hydrogen_type":args.hydrogen_type, "native_dt_ps":0.1,
            "stride":args.stride, "sampled_dt_ps":0.1*args.stride,
            "dumps":[str(x) for x in args.dumps]}
    (args.output/"metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    (args.output/"FINISHED.txt").write_text("Full-water angular-momentum audit finished successfully.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
