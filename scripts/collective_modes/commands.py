"""CLI command implementations for the canonical collective-dynamics mainline."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np

from .core import complex_acf, complex_cross_acf, coordinates, cylindrical_basis, cylindrical_currents, integrate_vacf, stable_order, velocity_in_frame, velocities
from .dump import Frame, infer_fluid_hint, infer_protocol_hint, inspect_dump, iter_frames, validate_uniform_timestep
from .output import require_columns, write_csv, write_metadata
from .schema import CaseProfile, REQUIREMENTS


def _numbers(spec: str) -> np.ndarray:
    values: list[int] = []
    for token in spec.split(","):
        if ":" in token:
            first, last = (int(value) for value in token.split(":", 1))
            values.extend(range(first, last + 1))
        else:
            values.append(int(token))
    return np.asarray(sorted(set(values)), dtype=int)


def _profile(args) -> CaseProfile:
    profile = CaseProfile(
        case_id=args.case_id,
        dump_paths=[Path(path) for path in args.dumps],
        wall_model=args.wall_model,
        axis_source=args.axis_source,
        rcnt_A=args.rcnt_A,
        oxygen_type=args.oxygen_type,
        fluid_types=tuple(args.fluid_types) if args.fluid_types else None,
        cnt_types=tuple(args.cnt_types or ()),
        integration_timestep_ps=args.timestep_ps,
        velocity_frame_default=args.velocity_frame,
        protocol_label=args.protocol,
        fluid_kind=args.fluid_kind,
    )
    if profile.axis_source == "fixed" and args.axis_xy is None:
        raise ValueError("--axis-source fixed requires --axis-xy X Y")
    if profile.wall_model == "implicit" and profile.rcnt_A is None:
        raise ValueError("--wall-model implicit requires --rcnt-A from the CNT field/protocol metadata")
    if profile.wall_model == "explicit_flexible" and not profile.cnt_types:
        raise ValueError("--wall-model explicit_flexible requires --cnt-types; do not assume a fixed box axis")
    return profile


def _frames(path: Path, selected_types: tuple[int, ...] | None, max_frames: int) -> list[Frame]:
    result: list[Frame] = []
    first_ids: np.ndarray | None = None
    for raw in iter_frames(path):
        frame = stable_order(raw.select_types(selected_types))
        if not len(frame.values):
            raise ValueError(f"{path}: selected fluid set is empty")
        if "id" in frame.fields:
            ids = frame.column("id").astype(np.int64)
            if first_ids is None:
                first_ids = ids
            elif not np.array_equal(first_ids, ids):
                raise ValueError(f"{path}: selected particle IDs change across frames")
        result.append(frame)
        if max_frames and len(result) >= max_frames:
            break
    if len(result) < 3:
        raise ValueError(f"{path}: need at least three selected frames")
    return result


def _axis(frame: Frame, args, profile: CaseProfile) -> np.ndarray:
    if profile.axis_source == "fixed":
        return np.asarray(args.axis_xy, dtype=float)
    if profile.axis_source == "box_center":
        return frame.box_center[:2]
    if profile.axis_source == "cnt_atoms":
        # Caller must give a full frame through this explicit path.  Water-only frames
        # cannot reconstruct a moving CNT axis without silently changing physics.
        raise ValueError("cnt_atoms axis requires a dedicated CNT-frame extractor; unsupported in water-only command")
    raise ValueError("axis source is unknown; declare fixed/box_center/cnt_atoms in the profile")


def _cadence(frames: list[Frame], profile: CaseProfile, declared_dt: float | None) -> float:
    inferred, interval = validate_uniform_timestep(frames, profile.integration_timestep_ps)
    if inferred is not None:
        if declared_dt is not None and not np.isclose(declared_dt, inferred):
            raise ValueError(f"declared dt {declared_dt} ps conflicts with dump cadence {inferred} ps ({interval} steps)")
        return float(inferred)
    if declared_dt is None:
        raise ValueError("one-frame input needs --dt-ps; otherwise declare --timestep-ps for cadence verification")
    return float(declared_dt)


def _aggregate(rows: list[dict], keys: list[str], value_fields: list[str], replica_key: str = "replica") -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(row)
    output = []
    for values, group in sorted(grouped.items()):
        record = dict(zip(keys, values))
        record["n_replicas"] = len(group)
        for field in value_fields:
            array = np.asarray([float(row[field]) for row in group], dtype=float)
            record[f"{field}_mean"] = float(np.nanmean(array))
            record[f"{field}_sem"] = float(np.nanstd(array, ddof=1) / np.sqrt(len(array))) if len(array) > 1 else float("nan")
        output.append(record)
    return output


def audit(args) -> None:
    profile = _profile(args)
    rows = []
    for source in profile.dump_paths:
        schema = inspect_dump(source)
        fields = set(schema.fields)
        support = {name: not bool(requirement.fields - fields) for name, requirement in REQUIREMENTS.items()}
        wall_hint, wall_basis = infer_protocol_hint(source, profile.wall_model, profile.cnt_types, schema)
        fluid_hint, fluid_basis = infer_fluid_hint(source, profile.fluid_kind, schema)
        rows.append({
            "case_id": profile.case_id, "path": str(source), "bytes": source.stat().st_size,
            "fields": " ".join(schema.fields), "atom_types": " ".join(map(str, schema.atom_types)),
            "inferred_content": schema.inferred_content, "inference_confidence": schema.confidence,
            "wall_model_inferred": wall_hint, "wall_model_basis": wall_basis,
            "fluid_kind_inferred": fluid_hint, "fluid_kind_basis": fluid_basis,
            "has_molecule_ids": schema.has_molecule_ids, "has_positions_xyz": schema.has_positions,
            "has_velocities_xyz": schema.has_velocities, "has_image_flags": schema.has_image_flags,
            **{f"supports_{name}": value for name, value in support.items()},
        })
    fields = list(rows[0])
    write_csv(args.output / "dump_capabilities.csv", rows, fields)
    write_metadata(args.output, {"case_id": profile.case_id, "protocol_label": profile.protocol_label,
        "wall_model": profile.wall_model, "axis_source": profile.axis_source, "rcnt_A": profile.rcnt_A,
        "important_limit": "dump-content inference detects fields/capabilities only. A water-only dump cannot prove explicit versus implicit CNT; declare wall_model in profile."})


def isf(args) -> None:
    profile = _profile(args)
    n_values, m_values = _numbers(args.n), _numbers(args.m)
    per_rows: list[dict] = []
    for replica, source in enumerate(profile.dump_paths, 1):
        frames = _frames(source, profile.selected_types, args.max_frames)
        cylindrical = bool(np.any(m_values))
        profile.require("cylindrical_isf" if cylindrical else "axial_isf", frames[0].fields)
        dt = _cadence(frames, profile, args.dt_ps)
        axis = _axis(frames[0], args, profile) if cylindrical else frames[0].box_center[:2]
        lz = frames[0].box_lengths[2]
        rcnt = profile.rcnt_A or 1.0
        positions = np.asarray([coordinates(frame, unwrapped=True) for frame in frames])
        if cylindrical:
            theta = np.asarray([cylindrical_basis(xyz, axis)[1] for xyz in positions])
        else:
            theta = np.zeros(positions.shape[:2])
        max_lag = min(int(round(args.max_lag_ps / dt)), len(frames) - 1)
        for n in n_values:
            for m in m_values:
                k = 2 * np.pi * n / lz
                phase = np.exp(1j * (k * positions[:, :, 2] + m * theta))
                rho = phase.sum(axis=1)
                total = complex_acf(rho[:, None], max_lag, demean=False)[:, 0].real / positions.shape[1]
                self_corr = complex_acf(phase, max_lag, demean=False).real.sum(axis=1) / positions.shape[1]
                for lag, (f_total, f_self) in enumerate(zip(total, self_corr)):
                    per_rows.append({"case_id": profile.case_id, "replica": replica, "source_path": str(source),
                        "n": int(n), "m": int(m), "kz_inv_A": float(k), "ktheta_inv_A": float(m / rcnt), "q_inv_A": float(np.hypot(k, m / rcnt)),
                        "lag_ps": lag * dt, "F_total": float(f_total), "F_self": float(f_self), "F_distinct": float(f_total - f_self),
                        "n_time_origins": len(frames) - lag})
    names = list(per_rows[0]); write_csv(args.output / "isf_per_replica.csv", per_rows, names)
    aggregate = _aggregate(per_rows, ["case_id", "n", "m", "kz_inv_A", "ktheta_inv_A", "q_inv_A", "lag_ps"], ["F_total", "F_self", "F_distinct"])
    write_csv(args.output / "isf_ensemble_mean_sem.csv", aggregate, list(aggregate[0]))
    write_metadata(args.output, {"case_id": profile.case_id, "observable": "F/Fs/Fd cylindrical density ISF", "n_values": n_values.tolist(), "m_values": m_values.tolist(), "profile": profile.__dict__})


def current(args) -> None:
    profile = _profile(args)
    n_values, m_values = _numbers(args.n), _numbers(args.m)
    if any(n == 0 and m == 0 for n in n_values for m in m_values):
        raise ValueError("current modes exclude (n,m)=(0,0); use a dedicated zero-mode observable if needed")
    per_rows: list[dict] = []; cross_rows: list[dict] = []; spectrum_rows: list[dict] = []
    for replica, source in enumerate(profile.dump_paths, 1):
        frames = _frames(source, profile.selected_types, args.max_frames)
        profile.require("cylindrical_current", frames[0].fields)
        dt = _cadence(frames, profile, args.dt_ps)
        axis = _axis(frames[0], args, profile)
        lz = frames[0].box_lengths[2]
        if profile.rcnt_A is None:
            raise ValueError("current modes require --rcnt-A from CNT geometry metadata")
        series: dict[str, list[np.ndarray]] = defaultdict(list)
        for frame in frames:
            xyz = coordinates(frame); velocity = velocity_in_frame(velocities(frame), args.velocity_frame)
            currents = cylindrical_currents(xyz, velocity, frame.box_lengths[2], axis, profile.rcnt_A, n_values, m_values)
            for channel in ("Jz", "Jr", "Jtheta", "L", "Tinplane", "Tr"):
                series[channel].append(currents[channel].reshape(-1))
        flattened = {key: np.asarray(value) for key, value in series.items()}
        max_lag = min(int(round(args.max_lag_ps / dt)), len(frames) - 1)
        for channel, values in flattened.items():
            acf = complex_acf(values, max_lag).real
            c0 = acf[0]
            for index, (n, m) in enumerate((int(n), int(m)) for n in n_values for m in m_values):
                for lag, value in enumerate(acf[:, index]):
                    per_rows.append({"case_id": profile.case_id, "replica": replica, "source_path": str(source), "channel": channel,
                        "n": n, "m": m, "kz_inv_A": float(2*np.pi*n/lz), "ktheta_inv_A": float(m/profile.rcnt_A), "q_inv_A": float(np.hypot(2*np.pi*n/lz, m/profile.rcnt_A)),
                        "lag_ps": lag*dt, "CJJ_raw": float(value), "CJJ_normalized": float(value/c0[index]) if c0[index] != 0 else float("nan"),
                        "CJJ0": float(c0[index]), "n_time_origins": len(frames)-lag})
            centred = values - values.mean(axis=0, keepdims=True)
            frequency = np.fft.fftfreq(len(values), dt)
            power = np.abs(np.fft.fft(centred, axis=0))**2 / len(values)
            for index, (n, m) in enumerate((int(n), int(m)) for n in n_values for m in m_values):
                for freq, item in zip(frequency[frequency >= 0], power[frequency >= 0, index]):
                    spectrum_rows.append({"case_id": profile.case_id, "replica": replica, "channel": channel, "n": n, "m": m,
                        "kz_inv_A": float(2*np.pi*n/lz), "ktheta_inv_A": float(m/profile.rcnt_A), "frequency_ps_inv": float(freq), "periodogram_power": float(item), "estimator": "whole-record rectangular periodogram"})
        for left, right in (("L", "Tinplane"), ("Tinplane", "L"), ("L", "Tr"), ("Tr", "L")):
            cross = complex_cross_acf(flattened[left], flattened[right], max_lag).real
            for index, (n, m) in enumerate((int(n), int(m)) for n in n_values for m in m_values):
                for lag, value in enumerate(cross[:, index]):
                    cross_rows.append({"case_id": profile.case_id, "replica": replica, "left_channel": left, "right_channel": right, "n": n, "m": m,
                        "kz_inv_A": float(2*np.pi*n/lz), "ktheta_inv_A": float(m/profile.rcnt_A), "q_inv_A": float(np.hypot(2*np.pi*n/lz, m/profile.rcnt_A)), "lag_ps": lag*dt, "C_AB_ordered_raw": float(value)})
    write_csv(args.output / "current_per_replica.csv", per_rows, list(per_rows[0]))
    aggregate = _aggregate(per_rows, ["case_id", "channel", "n", "m", "kz_inv_A", "ktheta_inv_A", "q_inv_A", "lag_ps"], ["CJJ_raw", "CJJ_normalized", "CJJ0"])
    write_csv(args.output / "current_ensemble_mean_sem.csv", aggregate, list(aggregate[0]))
    write_csv(args.output / "current_cross_ordered_per_replica.csv", cross_rows, list(cross_rows[0]))
    write_csv(args.output / "current_spectrum_per_replica.csv", spectrum_rows, list(spectrum_rows[0]))
    write_metadata(args.output, {"case_id": profile.case_id, "observable": "C_JJ channels and ordered cross kernels", "channels": ["Jz", "Jr", "Jtheta", "L", "Tinplane", "Tr"], "n_values": n_values.tolist(), "m_values": m_values.tolist(), "profile": profile.__dict__, "cross_kernel_definition": "Re<delta J_A(t+tau) delta J_B(t)*>; C_AB and C_BA are separately stored"})


def vacf(args) -> None:
    profile = _profile(args)
    per_rows: list[dict] = []; msd_rows: list[dict] = []
    for replica, source in enumerate(profile.dump_paths, 1):
        frames = _frames(source, profile.selected_types, args.max_frames)
        requirement = "vacf_z" if args.component == "z" else "vacf_cylindrical"
        profile.require(requirement, frames[0].fields)
        dt = _cadence(frames, profile, args.dt_ps)
        axis = _axis(frames[0], args, profile) if args.component != "z" else frames[0].box_center[:2]
        signal = []
        for frame in frames:
            velocity = velocity_in_frame(velocities(frame), args.velocity_frame)
            if args.component == "z":
                signal.append(velocity[:, 2])
            else:
                xyz = coordinates(frame); _, _, basis = cylindrical_basis(xyz, axis)
                er, etheta = basis[:, :, 0], basis[:, :, 1]
                signal.append(np.sum(velocity[:, :2] * (er if args.component == "r" else etheta), axis=1))
        values = np.asarray(signal)
        max_lag = min(int(round(args.max_lag_ps/dt)), len(values)-1)
        curve = complex_acf(values, max_lag).real.mean(axis=1)
        normalized = curve / curve[0]
        msd, alpha = integrate_vacf(curve, dt)
        for lag, value in enumerate(curve):
            per_rows.append({"case_id": profile.case_id, "replica": replica, "source_path": str(source), "component": args.component, "lag_ps": lag*dt, "VACF": float(value), "VACF_normalized": float(normalized[lag]), "n_time_origins": len(values)-lag})
            msd_rows.append({"case_id": profile.case_id, "replica": replica, "component": args.component, "lag_ps": lag*dt, "MSD_from_VACF": float(msd[lag]), "alpha_from_VACF": float(alpha[lag])})
    write_csv(args.output / "vacf_per_replica.csv", per_rows, list(per_rows[0]))
    vacf_mean = _aggregate(per_rows, ["case_id", "component", "lag_ps"], ["VACF", "VACF_normalized"])
    write_csv(args.output / "vacf_ensemble_mean_sem.csv", vacf_mean, list(vacf_mean[0]))
    write_csv(args.output / "msd_alpha_from_vacf_per_replica.csv", msd_rows, list(msd_rows[0]))
    msd_mean = _aggregate(msd_rows, ["case_id", "component", "lag_ps"], ["MSD_from_VACF", "alpha_from_VACF"])
    write_csv(args.output / "msd_alpha_from_vacf_ensemble_mean_sem.csv", msd_mean, list(msd_mean[0]))
    write_metadata(args.output, {"case_id": profile.case_id, "observable": "selected-fluid VACF -> consistent MSD and alpha", "velocity_frame": args.velocity_frame, "component": args.component, "profile": profile.__dict__, "limit": "replica SEM is only independent-configurational uncertainty when manifest says the starts are independent."})


def construct(args) -> None:
    with args.current_csv.open(newline="", encoding="utf-8") as handle:
        current = list(csv.DictReader(handle))
    with args.isf_csv.open(newline="", encoding="utf-8") as handle:
        isf = list(csv.DictReader(handle))
    with args.weights_csv.open(newline="", encoding="utf-8") as handle:
        weights = list(csv.DictReader(handle))
    require_columns(set(current[0]), {"n", "m", "lag_ps", args.current_column}, "current CSV")
    require_columns(set(isf[0]), {"n", "m", "lag_ps", args.self_column}, "ISF CSV")
    require_columns(set(weights[0]), {"n", "m", "weight"}, "weights CSV")
    weight_map = {(int(row["n"]), int(row["m"])): float(row["weight"]) for row in weights}
    require_columns(set(current[0]), {"channel"}, "current CSV")
    current_map = {(int(row["n"]), int(row["m"]), float(row["lag_ps"])): float(row[args.current_column]) for row in current if row["channel"] == args.current_channel}
    output = []
    for row in isf:
        n, m, lag = int(row["n"]), int(row["m"]), float(row["lag_ps"])
        key = (n, m, lag)
        if (n, m) not in weight_map or key not in current_map:
            continue
        fs = float(row[args.self_column]); phi = current_map[key]; weight = weight_map[(n, m)]
        output.append({"n": n, "m": m, "lag_ps": lag, "weight": weight, "Fs": fs, "Phi_J": phi, "term_WFsPhi": weight*fs*phi})
    if not output:
        raise ValueError("No matched (n,m,lag) rows across current, ISF, and weights inputs")
    grouped: dict[float, list[dict]] = defaultdict(list)
    for row in output: grouped[row["lag_ps"]].append(row)
    total = [{"lag_ps": lag, "construct_sum_WFsPhi": float(sum(item["term_WFsPhi"] for item in rows)), "n_modes": len(rows)} for lag, rows in sorted(grouped.items())]
    write_csv(args.output / "constructibility_per_mode.csv", output, list(output[0]))
    write_csv(args.output / "constructibility_sum.csv", total, list(total[0]))
    comparison = []
    if args.vacf_csv is not None:
        with args.vacf_csv.open(newline="", encoding="utf-8") as handle:
            vacf = list(csv.DictReader(handle))
        require_columns(set(vacf[0]), {"lag_ps", args.vacf_column}, "VACF CSV")
        vacf_map = {float(row["lag_ps"]): float(row[args.vacf_column]) for row in vacf}
        for row in total:
            if row["lag_ps"] in vacf_map:
                value = vacf_map[row["lag_ps"]]
                comparison.append({**row, "VACF_direct": value, "construct_minus_VACF": row["construct_sum_WFsPhi"] - value})
        if comparison:
            write_csv(args.output / "constructibility_vs_direct_vacf.csv", comparison, list(comparison[0]))
    write_metadata(args.output, {"formula": "sum_(n,m) W(n,m) F_s(n,m,t) Phi_J(n,m,t)", "weight_policy": "external static weights only; no amplitude fit performed", "current_channel": args.current_channel, "current_column": args.current_column, "self_column": args.self_column, "direct_vacf_comparison": str(args.vacf_csv) if args.vacf_csv else None})


def fit_current(args) -> None:
    """Fit the stored normalized current kernel mode-by-mode without pooling k."""
    from scipy.optimize import curve_fit
    with args.current_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    require_columns(set(rows[0]), {"n", "m", "lag_ps", "channel", args.column}, "current CSV")
    grouped: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        if row["channel"] == args.channel:
            grouped[(int(row["n"]), int(row["m"]))].append(row)
    if not grouped:
        raise ValueError(f"No rows for channel={args.channel!r}")
    def model(time, gamma, omega, a, b):
        return np.exp(-gamma*time) * (a*np.cos(omega*time) + b*np.sin(omega*time))
    output = []
    for (n, m), group in sorted(grouped.items()):
        group.sort(key=lambda row: float(row["lag_ps"])); time = np.asarray([float(row["lag_ps"]) for row in group]); values = np.asarray([float(row[args.column]) for row in group])
        mask = (time >= args.fit_min_ps) & (time <= args.fit_max_ps)
        x, y = time[mask], values[mask]
        if len(x) < 6:
            raise ValueError(f"(n,m)=({n},{m}) has fewer than six fit points")
        omega0 = 2*np.pi/max(args.fit_max_ps-args.fit_min_ps, np.finfo(float).eps)
        try:
            parameters, covariance = curve_fit(model, x, y, p0=(0.02, omega0, 1.0, 0.0), bounds=([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]), maxfev=100000)
            prediction = model(x, *parameters); r2 = 1 - np.sum((y-prediction)**2)/max(np.sum((y-y.mean())**2), np.finfo(float).eps); errors = np.sqrt(np.diag(covariance))
        except (RuntimeError, ValueError):
            parameters = np.full(4, np.nan); errors = np.full(4, np.nan); r2 = np.nan
        first = group[0]
        output.append({"channel": args.channel, "n": n, "m": m, "kz_inv_A": first.get("kz_inv_A", ""), "ktheta_inv_A": first.get("ktheta_inv_A", ""), "q_inv_A": float(np.hypot(float(first.get("kz_inv_A", 0)), float(first.get("ktheta_inv_A", 0)))), "Gamma_ps_inv": float(parameters[0]), "omega_rad_ps": float(parameters[1]), "a": float(parameters[2]), "b": float(parameters[3]), "Gamma_SE": float(errors[0]), "omega_SE": float(errors[1]), "a_SE": float(errors[2]), "b_SE": float(errors[3]), "fit_R2": float(r2), "fit_min_ps": args.fit_min_ps, "fit_max_ps": args.fit_max_ps, "model": "exp(-Gamma*t)[a*cos(omega*t)+b*sin(omega*t)]"})
    write_csv(args.output / "current_mode_DHO_parameters.csv", output, list(output[0]))
    write_metadata(args.output, {"input": str(args.current_csv), "channel": args.channel, "column": args.column, "formula": "Phi_J(k,m,t)=exp(-Gamma(k,m)t)[a(k,m)cos(omega(k,m)t)+b(k,m)sin(omega(k,m)t)]", "important_limit": "one row per (n,m); no k-law is assumed or fitted. Any Gamma(q), omega(q), a(q), b(q) relation must be a separate protocol-aware analysis."})


def plot(args) -> None:
    import matplotlib.pyplot as plt
    with args.csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    require_columns(set(rows[0]), {args.x, args.y}, "plot CSV")
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in rows: groups[row.get(args.group, "all")].append(row)
    figure, axis = plt.subplots(figsize=(8, 5))
    for label, group in sorted(groups.items()):
        group.sort(key=lambda row: float(row[args.x]))
        axis.plot([float(row[args.x]) for row in group], [float(row[args.y]) for row in group], label=label)
    axis.set_xlabel(args.x); axis.set_ylabel(args.y); axis.axhline(0.0, color="0.6", linewidth=0.8)
    if len(groups) > 1: axis.legend(title=args.group)
    axis.grid(alpha=0.25); figure.tight_layout(); args.output.parent.mkdir(parents=True, exist_ok=True); figure.savefig(args.output, dpi=160); plt.close(figure)
