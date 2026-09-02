#!/usr/bin/env python3
"""Compare matched physical-k ISF curves from the (8,8) 2L/3L demo.

The KWW form is deliberately used only as a descriptive fit to the first
positive decay branch.  It is not interpreted as a hydrodynamic exponent.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


KINDS = ("F_total", "F_self", "F_distinct", "C_J_normalized")
PAIRS = ((2, 3), (4, 6), (6, 9))


def read_curves(path):
    curves = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["n"]), row["kind"])
            curves.setdefault(key, {"k": float(row["k_inv_A"]), "lambda": float(row["lambda_A"]), "t": [], "y": []})
            curves[key]["t"].append(float(row["lag_ps"]))
            curves[key]["y"].append(float(row["real"]))
    for value in curves.values():
        value["t"] = np.asarray(value["t"], dtype=float)
        value["y"] = np.asarray(value["y"], dtype=float)
    return curves


def normalized_branch(curve, max_lag):
    t = curve["t"]
    y = curve["y"]
    y0 = y[0]
    if not np.isfinite(y0) or abs(y0) < 1e-14:
        return t[:0], y[:0]
    f = y / y0
    keep = t <= max_lag
    # Restrict to the first positive branch; later sign changes are oscillatory
    # collective dynamics, not a single KWW relaxation.
    bad = np.where((t > 0) & (f <= 0))[0]
    if len(bad):
        keep &= np.arange(len(t)) < bad[0]
    return t[keep], f[keep]


def normalized_full(curve, max_lag):
    t = curve["t"]
    y = curve["y"]
    y0 = y[0]
    if not np.isfinite(y0) or abs(y0) < 1e-14:
        return t[:0], y[:0]
    keep = t <= max_lag
    return t[keep], (y / y0)[keep]


def crossing_time(t, f, target=math.exp(-1)):
    indices = np.where((f[:-1] >= target) & (f[1:] < target))[0]
    if not len(indices):
        return None
    i = int(indices[0])
    return float(t[i] + (target - f[i]) * (t[i + 1] - t[i]) / (f[i + 1] - f[i]))


def kww_fit(t, f):
    # log[-log f] = beta log(t) - beta log(tau)
    mask = (t >= 1.0) & (f > 0.10) & (f < 0.95)
    if int(mask.sum()) < 8:
        return {"beta": None, "tau_ps": None, "r2_log": None, "n_fit": int(mask.sum())}
    x = np.log(t[mask])
    y = np.log(-np.log(f[mask]))
    beta, intercept = np.polyfit(x, y, 1)
    if not np.isfinite(beta) or beta <= 0:
        return {"beta": None, "tau_ps": None, "r2_log": None, "n_fit": int(mask.sum())}
    predicted = beta * x + intercept
    denom = float(np.sum((y - y.mean()) ** 2))
    r2 = None if denom == 0 else float(1.0 - np.sum((y - predicted) ** 2) / denom)
    tau = float(np.exp(-intercept / beta))
    return {"beta": float(beta), "tau_ps": tau, "r2_log": r2, "n_fit": int(mask.sum())}


def curve_agreement(a, b, max_lag, positive_branch=True):
    normalizer = normalized_branch if positive_branch else normalized_full
    ta, fa = normalizer(a, max_lag)
    tb, fb = normalizer(b, max_lag)
    if len(ta) < 3 or len(tb) < 3:
        return {"n": 0, "rms": None, "max_abs": None}
    common = ta[(ta >= tb[0]) & (ta <= tb[-1])]
    if not len(common):
        return {"n": 0, "rms": None, "max_abs": None}
    delta = fa[np.isin(ta, common)] - np.interp(common, tb, fb)
    return {"n": int(len(common)), "rms": float(np.sqrt(np.mean(delta ** 2))), "max_abs": float(np.max(np.abs(delta)))}


def compact(value):
    return "--" if value is None else "{:.4g}".format(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-lag-ps", type=float, default=2000.0)
    args = parser.parse_args()

    curves2 = read_curves(args.input_root / "2L" / "intermediate_scattering_curves.csv")
    curves3 = read_curves(args.input_root / "3L" / "intermediate_scattering_curves.csv")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for n2, n3 in PAIRS:
        for label, curves, n in (("2L", curves2, n2), ("3L", curves3, n3)):
            for kind in KINDS:
                c = curves[(n, kind)]
                t, f = normalized_branch(c, args.max_lag_ps)
                record = {"pair": "2L:n{} = 3L:n{}".format(n2, n3), "system": label, "n": n,
                          "kind": kind, "k_inv_A": c["k"], "lambda_A": c["lambda"],
                          "tau_1e_ps": crossing_time(t, f)}
                if kind in ("F_total", "F_self"):
                    record.update(kww_fit(t, f))
                else:
                    record.update({"beta": None, "tau_ps": None, "r2_log": None, "n_fit": 0})
                records.append(record)

        for kind in KINDS:
            a = curves2[(n2, kind)]
            b = curves3[(n3, kind)]
            current = kind == "C_J_normalized"
            agreement_window = min(args.max_lag_ps, 50.0) if current else args.max_lag_ps
            agreement = curve_agreement(a, b, agreement_window, positive_branch=not current)
            records.append({"pair": "2L:n{} = 3L:n{}".format(n2, n3), "system": "agreement",
                            "n": None, "kind": kind, "k_inv_A": 0.5 * (a["k"] + b["k"]),
                            "lambda_A": 0.5 * (a["lambda"] + b["lambda"]),
                            "tau_1e_ps": None, "beta": None, "tau_ps": None,
                            "r2_log": None, "n_fit": agreement["n"],
                            "rms_normalized": agreement["rms"], "max_abs_normalized": agreement["max_abs"],
                            "agreement_window_ps": agreement_window})

    fields = ["pair", "system", "n", "kind", "k_inv_A", "lambda_A", "tau_1e_ps", "beta", "tau_ps", "r2_log", "n_fit", "rms_normalized", "max_abs_normalized", "agreement_window_ps"]
    with open(args.output_dir / "matched_k_isf_kww_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in records:
            writer.writerow(item)
    with open(args.output_dir / "matched_k_isf_kww_summary.json", "w") as handle:
        json.dump(records, handle, indent=2, allow_nan=False)

    lines = ["# (8,8) 2L/3L matched-k ISF demonstration", "",
             "KWW fits use only the first positive decay branch, with 0.10 < normalized correlation < 0.95 and t >= 1 ps.",
             "They are descriptive relaxation parameters, not hydrodynamic exponents.", ""]
    for n2, n3 in PAIRS:
        pair = "2L:n{} = 3L:n{}".format(n2, n3)
        k = curves2[(n2, "F_self")]["k"]
        lam = curves2[(n2, "F_self")]["lambda"]
        lines += ["## {}  (k={:.6f} A^-1, lambda={:.3f} A)".format(pair, k, lam), "",
                  "| quantity | 2L | 3L | matched-curve RMS |", "|---|---:|---:|---:|"]
        for kind in ("F_total", "F_self"):
            r2 = next(r for r in records if r["pair"] == pair and r["system"] == "2L" and r["kind"] == kind)
            r3 = next(r for r in records if r["pair"] == pair and r["system"] == "3L" and r["kind"] == kind)
            ra = next(r for r in records if r["pair"] == pair and r["system"] == "agreement" and r["kind"] == kind)
            lines.append("| {} beta / tau(ps) / tau_1e(ps) | {} / {} / {} | {} / {} / {} | {} |".format(
                kind, compact(r2["beta"]), compact(r2["tau_ps"]), compact(r2["tau_1e_ps"]),
                compact(r3["beta"]), compact(r3["tau_ps"]), compact(r3["tau_1e_ps"]), compact(ra.get("rms_normalized"))))
        for kind in ("F_distinct", "C_J_normalized"):
            ra = next(r for r in records if r["pair"] == pair and r["system"] == "agreement" and r["kind"] == kind)
            lines.append("| {} normalized-curve agreement (0-{} ps) | -- | -- | {} |".format(kind, compact(ra.get("agreement_window_ps")), compact(ra.get("rms_normalized"))))
        lines.append("")
    (args.output_dir / "matched_k_isf_kww_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
