"""Audit C99 implicit-CNT matched-k CJJ consistency and VACF constructibility.

This deliberately operates only on compact, locally verified derived tables.  It
does *not* infer the tagged--collective static vertex W_n or F_s(k,t): raw
oxygen trajectories are absent locally, so a direct no-free-amplitude
reconstruction is not identifiable from these inputs.  The N1600 mode exercise
is instead a cross-velocity-seed, non-negative linear-span diagnostic.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 1.0, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

ROOT = Path(r"H:/gcmc_explore/implicit_chirality_length_scan_20260816")
SHORT = ROOT / "c99_T350_weakNH_r4"
LONG = ROOT / "c99_T350_weakNH_r4_N2400_N3200_extension_20260817"
OUT = Path(r"H:/gcmc_explore/translational_anomaly/02_isf_collective_modes/results/collective_mode_response/implicit_C99_VACF_matched_k_constructibility/2026-08-29")

CASES = ("N400", "N800", "N1600", "N2400", "N3200")
LZ_A = {"N400": 200.0, "N800": 400.0, "N1600": 800.0, "N2400": 1200.0, "N3200": 1600.0}
COLORS = {"N400": "#4C78A8", "N800": "#F58518", "N1600": "#54A24B", "N2400": "#B279A2", "N3200": "#E45756"}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def source_root(case: str) -> Path:
    return SHORT if case in ("N400", "N800", "N1600") else LONG


def analysis_dir(case: str) -> Path:
    date = "20260816" if case in ("N400", "N800", "N1600") else "20260818"
    return source_root(case) / case / f"analysis_msd_cvv_modes_200ps_{date}"


def read_vacf(case: str, extended: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if extended:
        base = SHORT / "N1600/analysis_vacf_500ps_extension_20260817/vacf"
    else:
        base = analysis_dir(case) / "vacf"
    frames = []
    for rep in range(1, 5):
        frames.append(pd.read_csv(base / f"rep{rep}/axial_vacf_tail.csv")[["lag_ps", "vacf_peculiar_mean"]].rename(columns={"vacf_peculiar_mean": f"rep{rep}"}))
    m = frames[0]
    for f in frames[1:]:
        m = m.merge(f, on="lag_ps", validate="one_to_one")
    values = m[[f"rep{i}" for i in range(1, 5)]].to_numpy(float)
    return m["lag_ps"].to_numpy(float), values.mean(axis=1), values.std(axis=1, ddof=1) / 2.0


def read_cjj(case: str, rep: int) -> pd.DataFrame:
    return pd.read_csv(SHORT / "cjj_time_domain_lowk_20260817" / case / f"rep{rep}_cjj_time.csv")


def matched_k() -> tuple[list[dict], dict[tuple[str, int], tuple[np.ndarray, np.ndarray]]]:
    # Same physical k occurs as N400:n, N800:2n, N1600:4n.
    curves: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]] = {}
    for case in ("N400", "N800", "N1600"):
        byrep = []
        for rep in range(1, 5):
            d = read_cjj(case, rep)
            byrep.append(d)
        all_d = pd.concat(byrep, keys=range(1, 5), names=["rep", "row"]).reset_index(level=0)
        for n, g in all_d.groupby("mode_n"):
            pivot = g.pivot(index="lag_ps", columns="rep", values="CJJ_real_normalized").sort_index()
            curves[(case, int(n))] = (pivot.index.to_numpy(float), pivot.to_numpy(float).mean(axis=1))
    rows = []
    for base_n in (1, 2):
        members = [("N400", base_n), ("N800", 2 * base_n), ("N1600", 4 * base_n)]
        if not all(key in curves for key in members):
            continue
        time = curves[members[0]][0]
        common = [np.interp(time, curves[k][0], curves[k][1]) for k in members]
        ref = common[-1]
        for (case, n), y in zip(members[:-1], common[:-1]):
            rows.append({
                "k_Ainv": 2 * np.pi * base_n / LZ_A["N400"], "reference": "N1600",
                "case": case, "mode_n": n, "reference_mode_n": 4 * base_n,
                "time_window_ps": "0-160", "rmse_vs_N1600": float(np.sqrt(np.mean((y-ref)**2))),
                "pearson_r_vs_N1600": float(np.corrcoef(y, ref)[0, 1]),
                "caveat": "normalized CJJ only; four velocity seeds share a parent configuration",
            })
    return rows, curves


def constructibility() -> tuple[list[dict], list[dict], dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    """LOO nonnegative spans: held-out VACF versus the mean CJJ basis of other seeds."""
    vacf = []
    for rep in range(1, 5):
        f = pd.read_csv(analysis_dir("N1600") / f"vacf/rep{rep}/axial_vacf_tail.csv")
        vacf.append(f.set_index("lag_ps")["vacf_peculiar_mean"])
    cjj = {rep: read_cjj("N1600", rep) for rep in range(1, 5)}
    common_t = np.arange(0.1, 160.0 + 1e-9, 0.1)
    summary, weights = [], []
    display = {}
    for max_n in range(1, 8):
        pred_all, obs_all = [], []
        first = None
        for test_rep in range(1, 5):
            train = [r for r in range(1, 5) if r != test_rep]
            cols = []
            for n in range(1, max_n + 1):
                ys = []
                for r in train:
                    d = cjj[r].query("mode_n == @n").sort_values("lag_ps")
                    ys.append(np.interp(common_t, d.lag_ps, d.CJJ_real_normalized))
                cols.append(np.mean(ys, axis=0))
            A = np.column_stack(cols)
            heldout = vacf[test_rep - 1]
            y = np.interp(common_t, heldout.index.to_numpy(float), heldout.to_numpy(float))
            w, _ = nnls(A, y)
            pred = A @ w
            pred_all.append(pred); obs_all.append(y)
            for n, value in enumerate(w, start=1):
                weights.append({"M": max_n, "held_out_velocity_seed": test_rep, "mode_n": n, "coefficient": float(value),
                                "definition": "nonnegative regression proxy, not static tagged-collective W_n"})
            if test_rep == 1:
                first = (common_t, y, pred)
        p, y = np.concatenate(pred_all), np.concatenate(obs_all)
        rmse = float(np.sqrt(np.mean((p-y)**2)))
        null_rmse = float(np.sqrt(np.mean(y**2)))
        summary.append({"M": max_n, "lag_window_ps": "0.1-160", "n_holdout_velocity_seeds": 4,
                        "cv_rmse": rmse, "cv_pearson_r": float(np.corrcoef(p, y)[0, 1]),
                        "rmse_over_zero_predictor": rmse / null_rmse,
                        "interpretation": "cross-seed linear-span diagnostic only; lacks Fs and static vertex normalization"})
        display[max_n] = first
    return summary, weights, display


def plot_vacf() -> None:
    fig = plt.figure(figsize=(7.0, 3.15))
    ax1 = fig.add_axes([0.10, 0.19, 0.52, 0.73])
    ax2 = fig.add_axes([0.72, 0.19, 0.24, 0.73])
    for case in CASES:
        t, y, sem = read_vacf(case)
        keep = (t >= 1) & (t <= 200)
        label = rf"$L_z={LZ_A[case]/10:g}$ nm"
        ax1.plot(t[keep], y[keep], color=COLORS[case], lw=1.1, label=label)
        ax1.fill_between(t[keep], y[keep]-sem[keep], y[keep]+sem[keep], color=COLORS[case], alpha=0.13, lw=0)
    t, y, sem = read_vacf("N1600", extended=True)
    keep = (t >= 1) & (t <= 500)
    ax2.plot(t[keep], y[keep], color=COLORS["N1600"], lw=1.1)
    ax2.fill_between(t[keep], y[keep]-sem[keep], y[keep]+sem[keep], color=COLORS["N1600"], alpha=0.15, lw=0)
    for ax in (ax1, ax2):
        ax.axhline(0, color="0.25", lw=1.0)
        ax.set_xlabel(r"lag time $t$ (ps)")
        ax.tick_params(direction="out", width=1.0, length=3)
    ax1.set_ylabel(r"peculiar axial VACF / $C_{vv}(0)$")
    ax1.set_xlim(1, 200); ax2.set_xlim(1, 500)
    ax1.legend(loc="upper right", ncol=1, handlelength=1.4)
    ax1.text(-0.13, 1.04, "(a)", transform=ax1.transAxes, fontweight="bold", fontsize=9)
    ax2.text(-0.22, 1.04, "(b)", transform=ax2.transAxes, fontweight="bold", fontsize=9)
    ax2.text(0.03, 0.95, r"$N=1600$ only", transform=ax2.transAxes, va="top")
    fig.savefig(OUT / "figures/C99_VACF_N_dependence_1to200ps_and_N1600_500ps.svg")
    fig.savefig(OUT / "figures/C99_VACF_N_dependence_1to200ps_and_N1600_500ps.pdf")
    fig.savefig(OUT / "figures/C99_VACF_N_dependence_1to200ps_and_N1600_500ps.png", dpi=600)
    fig.savefig(OUT / "figures/C99_VACF_N_dependence_1to200ps_and_N1600_500ps.tiff", dpi=600)
    plt.close(fig)


def plot_matched_k(curves: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]]) -> None:
    fig = plt.figure(figsize=(5.5, 2.65)); ax = fig.add_axes([0.12, 0.21, 0.84, 0.71])
    styles = {"N400": "-", "N800": "--", "N1600": ":"}
    for base_n in (1, 2):
        for case, n in (("N400", base_n), ("N800", 2*base_n), ("N1600", 4*base_n)):
            if (case, n) not in curves: continue
            t, y = curves[(case, n)]
            ax.plot(t, y, lw=1.05, ls=styles[case], color=COLORS[case], label=rf"$L_z={LZ_A[case]/10:g}$ nm, $n={n}$")
    ax.axhline(0, color="0.25", lw=1.0); ax.set_xlim(0, 160)
    ax.set_xlabel(r"lag time $t$ (ps)"); ax.set_ylabel(r"normalized $C_{JJ}(k,t)$")
    ax.tick_params(direction="out", width=1.0, length=3); ax.legend(ncol=2, fontsize=6.4, handlelength=1.4, columnspacing=0.8)
    ax.text(-0.10, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.savefig(OUT / "figures/C99_matched_k_CJJ_time_domain.svg")
    fig.savefig(OUT / "figures/C99_matched_k_CJJ_time_domain.pdf")
    fig.savefig(OUT / "figures/C99_matched_k_CJJ_time_domain.png", dpi=600)
    fig.savefig(OUT / "figures/C99_matched_k_CJJ_time_domain.tiff", dpi=600)
    plt.close(fig)


def plot_constructibility(display: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]) -> None:
    fig = plt.figure(figsize=(5.5, 2.65)); ax = fig.add_axes([0.12, 0.21, 0.84, 0.71])
    t, y, _ = display[7]; ax.plot(t, y, color="0.15", lw=1.2, label="held-out $N=1600$ VACF")
    for M, color in ((1, "#E45756"), (3, "#F58518"), (7, "#54A24B")):
        _, _, p = display[M]; ax.plot(t, p, color=color, lw=1.0, label=rf"CJJ span: $n=1\ldots{M}$")
    ax.axhline(0, color="0.25", lw=1.0); ax.set_xlim(0, 160)
    ax.set_xlabel(r"lag time $t$ (ps)"); ax.set_ylabel(r"normalized VACF proxy")
    ax.tick_params(direction="out", width=1.0, length=3); ax.legend(ncol=2, fontsize=6.5, handlelength=1.4, columnspacing=0.8)
    ax.text(-0.10, 1.04, "(a)", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.savefig(OUT / "figures/C99_N1600_lowk_CJJ_span_constructibility.svg")
    fig.savefig(OUT / "figures/C99_N1600_lowk_CJJ_span_constructibility.pdf")
    fig.savefig(OUT / "figures/C99_N1600_lowk_CJJ_span_constructibility.png", dpi=600)
    fig.savefig(OUT / "figures/C99_N1600_lowk_CJJ_span_constructibility.tiff", dpi=600)
    plt.close(fig)


def main() -> None:
    (OUT / "derived_data").mkdir(parents=True, exist_ok=True); (OUT / "figures").mkdir(exist_ok=True)
    matched_rows, curves = matched_k()
    construct_rows, coefficient_rows, display = constructibility()
    vacf_rows = []
    for case in CASES:
        t, y, sem = read_vacf(case)
        for ti, yi, si in zip(t, y, sem):
            vacf_rows.append({"case": case, "N_water": int(case[1:]), "Lz_A": LZ_A[case], "Lz_nm": LZ_A[case]/10,
                              "lag_ps": ti, "vacf_peculiar_mean": yi, "velocity_seed_sem": si, "available_window": "0-200 ps"})
    t, y, sem = read_vacf("N1600", extended=True)
    for ti, yi, si in zip(t, y, sem):
        vacf_rows.append({"case": "N1600", "N_water": 1600, "Lz_A": 800.0, "Lz_nm": 80.0,
                          "lag_ps": ti, "vacf_peculiar_mean": yi, "velocity_seed_sem": si, "available_window": "0-500 ps extension"})
    write_csv(OUT / "derived_data/C99_VACF_ensemble_mean_seedSEM.csv", vacf_rows)
    write_csv(OUT / "derived_data/C99_matched_k_CJJ_consistency.csv", matched_rows)
    write_csv(OUT / "derived_data/C99_N1600_lowk_CJJ_span_LOO_summary.csv", construct_rows)
    write_csv(OUT / "derived_data/C99_N1600_lowk_CJJ_span_LOO_coefficients.csv", coefficient_rows)
    plot_vacf(); plot_matched_k(curves); plot_constructibility(display)
    metadata = {
        "date": "2026-08-29", "system": "implicit CNT C99, (9,9), 350 K, weak-NH, 4 velocity seeds", "input_status": "compact derived CJJ/VACF tables locally verified",
        "box_lengths_A": LZ_A, "CJJ_time_domain_boxes": ["N400", "N800", "N1600"], "VACF_boxes": list(CASES),
        "CJJ_time_window_ps": "0-160", "VACF_common_window_ps": "0-200", "VACF_500ps": "N1600 only",
        "static_weight_status": "not identifiable: normalized CJJ has CJJ(0)=1 and no raw tagged trajectory / K,c,W / Fs arrays are local",
        "constructibility_status": "LOO nonnegative CJJ-span proxy; not tagged-collective no-free-amplitude reconstruction",
        "uncertainty": "velocity-seed SEM; seeds share a parent configuration",
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (OUT / "README.md").write_text("""# C99 implicit-CNT VACF, matched-k CJJ, and constructibility audit\n\n## Scope\n\nThis locally reproducible package advances the C99 `(9,9)` implicit-CNT, 350 K weak-NH series using only compact verified tables.  It contains: (i) physical matched-k normalized longitudinal-current time-kernel comparisons across N400/N800/N1600; (ii) a normalized axial peculiar-VACF size panel; and (iii) an N1600 leave-one-velocity-seed-out low-k CJJ linear-span diagnostic.\n\n## What it supports\n\n- Matched-k comparisons use N400:n, N800:2n and N1600:4n.  Each record carries RMSE and Pearson r against N1600 over the common 0--160 ps window.\n- The VACF figure labels actual `Lz` in nm.  All five N cases have a common 0--200 ps source window; only N1600 has a separately verified 0--500 ps extension.  The panel therefore never disguises the missing 200--500 ps data for the other boxes.\n- The N1600 mode diagnostic asks a deliberately limited question: can the held-out seed VACF lie in the nonnegative span of the *other-seed* normalized CJJ n=1..M curves?  Its coefficients are not physical weights.\n\n## Boundaries\n\nNo local raw trajectories exist for these cases.  Consequently the package cannot calculate tagged static K/c/a/W, `N W(k)`, self-ISF `Fs(k,t)`, raw `CJJ(0)`, or a no-free-amplitude VACF reconstruction.  Normalized `CJJ(k,t)` has `CJJ(0)=1` by construction and cannot provide static spectral weight.  N2400/N3200 do not have locally verified time-domain CJJ, so N1600 cannot predict their intermediate-time VACF without new target-system current/vertex inputs.\n\nAll four replicas are velocity seeds from one parent configuration; SEM is seed-conditional, not an independent-configurational uncertainty.\n""", encoding="utf-8")
    (OUT / "QA.md").write_text("""# QA\n\n- Inputs were read from the C99 paths registered on 2026-08-29.\n- Every CJJ matched-k curve is an average over four seed-specific normalized curves; no discrete mode-index-only comparison is used.\n- The 1--500 ps request is partially data-limited: N400/N800/N2400/N3200 stop at 200 ps in the verified compact VACF archive.  N1600 alone has the verified 500 ps extension.\n- The constructibility analysis is cross-seed NNLS with no intercept, using 0.1--160 ps where CJJ is available.  It is explicitly not a static-weight or Mori/vertex closure.\n- No data were excluded except lag 0 from regression because it trivially fixes all normalized CJJ bases to one.\n""", encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("Completed locally available C99 matched-k/VACF/constructibility audit; static vertex reconstruction remains data-blocked.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
