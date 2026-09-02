"""Compare fixed k exponents for old implicit-CNT TA-theta linewidths."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def fit_fixed(k: np.ndarray, gamma: np.ndarray, sigma: np.ndarray, beta: float) -> dict[str, float]:
    x = k**beta
    weight = 1.0 / sigma**2
    design = np.column_stack([np.ones_like(x), x])
    normal = design.T @ (weight[:, None] * design)
    cov = np.linalg.inv(normal)
    par = cov @ (design.T @ (weight * gamma))
    residual = gamma - design @ par
    chi2 = float(np.sum((residual / sigma) ** 2))
    return {"Gamma0_rad_ps": par[0], "Gamma0_SE": np.sqrt(cov[0, 0]),
            "D_rad_ps_A_to_beta": par[1], "D_SE": np.sqrt(cov[1, 1]),
            "chi2": chi2, "dof": len(k) - 2, "reduced_chi2": chi2 / (len(k) - 2)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nmax", type=int, default=6)
    args = parser.parse_args()
    rows: list[dict[str, float | int | str]] = []
    for rep in range(1, 5):
        source = args.input_root / "per_replica" / f"rep{rep}" / "derived_data" / "TAtheta_zeropeak_lorentzian_fits_n001_n020.csv"
        data = pd.read_csv(source)
        data = data[(data.n <= args.nmax) & data.used_primary_k_fit].copy()
        for beta, label in [(1.5, "Gamma0_plus_Dk3over2"), (2.0, "Gamma0_plus_Dk2")]:
            result = fit_fixed(data.k_Ainv.to_numpy(), data.gamma_rad_ps.to_numpy(), data.gamma_segment_sem_rad_ps.to_numpy(), beta)
            rows.append({"replica": rep, "model": label, "beta": beta, "n_used": ",".join(map(str, data.n)), **result})
    result = pd.DataFrame(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output / "implicit_old_weakNH_TAtheta_fixed_exponent_comparison.csv", index=False)
    aggregate = result.groupby(["model", "beta"])[["Gamma0_rad_ps", "D_rad_ps_A_to_beta", "reduced_chi2"]].agg(["mean", "sem"]).reset_index()
    aggregate.columns = ["_".join(c).rstrip("_") for c in aggregate.columns]
    aggregate.to_csv(args.output / "implicit_old_weakNH_TAtheta_fixed_exponent_replica_summary.csv", index=False)


if __name__ == "__main__":
    main()
