"""Unify protocol-matched 8-rep axial VACF and VACF-integral ODE alpha plots.

The discrete lengths available under the common (8,8), weak-NH/no-momentum,
10-fs/1-ns protocol are 2, 3, 4, 5 and 10 L.  No interpolation is made for
6--9 L.  Alpha is evaluated separately for each replica from
I(t)=integral_0^t Cvv(s)ds, J(t)=integral_0^t I(s)ds, alpha=t I/J, then
reported as the replica mean and SEM.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LOW = ROOT / "results/collective_mode_response/vacf_88_L2L5_10fs_1ns_8rep_weakNH_nomom/2026-08-21/allorigins"
TEN = ROOT / "results/collective_mode_response/vacf_tail_8_8_L10_10fs_8rep_1ns_2026-08-19/analysis_cvv_alpha_200ps_8rep"
OUT = ROOT / "results/collective_mode_response/vacf_alpha_88_L2L5L10_10fs_1ns_8rep/2026-08-21/allorigins"
COLORS = {"2L": "#3b6fb6", "3L": "#e07b39", "4L": "#4d9a65", "5L": "#b34c66", "10L": "#7656a6"}


def trapz_cumulative(y: np.ndarray, dt: float) -> np.ndarray:
    out = np.zeros_like(y)
    out[1:] = np.cumsum((y[:-1] + y[1:]) * (0.5 * dt), axis=0)
    return out


def alpha_from_curves(lag: np.ndarray, curves: np.ndarray) -> np.ndarray:
    dt = float(np.median(np.diff(lag)))
    first = trapz_cumulative(curves, dt)
    second = trapz_cumulative(first, dt)
    alpha = np.full_like(curves, np.nan)
    good = np.abs(second) > 1e-15
    alpha[good] = (lag[:, None] * first)[good] / second[good]
    return alpha


def load_low(label: str):
    raw = np.genfromtxt(LOW / "per_replica" / f"VACF_8_8_{label}_peculiar_per_replica_normalised.csv",
                        delimiter=",", names=True)
    lag = np.asarray(raw["lag_ps"], float)
    curves = np.column_stack([np.asarray(raw[f"rep{i}"], float) for i in range(1, 9)])
    return lag, curves


def load_ten():
    raw = np.genfromtxt(TEN / "cvv_per_replica.csv", delimiter=",", names=True)
    reps = np.asarray(raw["replica"], int)
    lag = np.unique(np.asarray(raw["lag_ps"], float))
    curves = np.column_stack([np.asarray(raw["cvv_A2_ps2"], float)[reps == rep]
                              for rep in range(1, 9)])
    if curves.shape[0] != lag.size:
        raise ValueError("10L replica lag grids are inconsistent")
    return lag, curves / curves[0:1]


def save_csv(path: Path, header, values):
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh); writer.writerow(header); writer.writerows(values)


def style(ax):
    ax.tick_params(direction="out", width=1.0, length=3, top=False, right=False, labelsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"): ax.spines[side].set_linewidth(1.0)


def panel_label(ax, tag):
    ax.text(-0.17, 1.05, tag, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def main():
    OUT.mkdir(parents=True, exist_ok=True); figdir = OUT / "figures"; figdir.mkdir(exist_ok=True)
    perdir = OUT / "per_replica"; perdir.mkdir(exist_ok=True)
    datasets = {}
    for label in ("2L", "3L", "4L", "5L"):
        datasets[label] = load_low(label)
    datasets["10L"] = load_ten()
    summary = {"observable": "normalised axial oxygen peculiar VACF and VACF-integral ODE alpha",
               "protocol": "(8,8), weak NH, no momentum removal, 10 fs dump cadence, 1 ns, 8 velocity-seed replicas",
               "available_nominal_lengths": [2, 3, 4, 5, 10], "missing_nominal_lengths": [6, 7, 8, 9],
               "alpha_definition": "I=int_0^t Cvv(s)ds; J=int_0^t I(s)ds; alpha=t I/J; each replica first, then mean and SEM",
               "inputs": {"2L-5L": str(LOW.resolve()), "10L": str(TEN.resolve())}, "per_length": {}}
    processed = {}
    minima_rows = []
    for label, (lag0, curves0) in datasets.items():
        keep = lag0 <= 100.0 + 1e-12
        lag, curves = lag0[keep], curves0[keep]
        alpha = alpha_from_curves(lag, curves)
        vacf_mean, vacf_sem = curves.mean(1), curves.std(1, ddof=1) / np.sqrt(8)
        alpha_mean, alpha_sem = np.full(len(lag), np.nan), np.full(len(lag), np.nan)
        alpha_mean[1:] = np.nanmean(alpha[1:], axis=1)
        alpha_sem[1:] = np.nanstd(alpha[1:], axis=1, ddof=1) / np.sqrt(8)
        processed[label] = (lag, vacf_mean, vacf_sem, alpha, alpha_mean, alpha_sem)
        save_csv(OUT / f"VACF_alpha_ODE_8_8_{label}_mean_sem.csv",
                 ["lag_ps", "vacf_normalised_mean", "vacf_replica_sem", "alpha_ode_mean", "alpha_ode_replica_sem", "n_replicas"],
                 np.column_stack([lag, vacf_mean, vacf_sem, alpha_mean, alpha_sem, np.full_like(lag, 8)]))
        save_csv(perdir / f"alpha_ODE_8_8_{label}_per_replica.csv",
                 ["lag_ps"] + [f"rep{i}" for i in range(1, 9)], np.column_stack([lag, alpha]))
        idx5 = int(np.argmin(abs(lag - 5.0))); idx100 = int(np.argmin(abs(lag - 100.0)))
        summary["per_length"][label] = {"n_replicas": 8, "points": int(len(lag)),
          "alpha_5ps_mean_sem": [float(alpha_mean[idx5]), float(alpha_sem[idx5])],
          "alpha_100ps_mean_sem": [float(alpha_mean[idx100]), float(alpha_sem[idx100])]}
        min_mask = (lag >= 5.0) & (lag <= 100.0)
        mean_idx = np.flatnonzero(min_mask)[int(np.nanargmin(alpha_mean[min_mask]))]
        for rep in range(8):
            rep_idx = np.flatnonzero(min_mask)[int(np.nanargmin(alpha[min_mask, rep]))]
            minima_rows.append([label, rep + 1, lag[rep_idx], alpha[rep_idx, rep],
                                lag[mean_idx], alpha_mean[mean_idx], alpha_sem[mean_idx]])
    save_csv(OUT / "alpha_ODE_minima_5to100ps_per_replica.csv",
             ["length", "replica", "replica_t_min_ps", "replica_alpha_min",
              "mean_curve_t_min_ps", "mean_curve_alpha_min", "replica_sem_at_mean_curve_min"], minima_rows)
    # Explicit three-column bbox layout: signed tail, amplitude tail, alpha.
    fig = plt.figure(figsize=(7.0, 2.55), dpi=300)
    axes = [fig.add_axes([0.09, 0.22, 0.25, 0.68]), fig.add_axes([0.41, 0.22, 0.25, 0.68]), fig.add_axes([0.73, 0.22, 0.24, 0.68])]
    for label, (lag, vm, vs, _, am, ass) in processed.items():
        color = COLORS[label]; maskv = lag >= 1; maska = lag >= 0.1
        axes[0].plot(lag[maskv], vm[maskv], color=color, lw=1.15, label=label)
        axes[0].fill_between(lag[maskv], vm[maskv]-vs[maskv], vm[maskv]+vs[maskv], color=color, alpha=.16, lw=0)
        # Log-log diagnostic omits exact and near-zero crossings rather than
        # imposing an arbitrary positive floor.
        ampl = np.abs(vm[maskv]); valid = ampl > 3.0 * vs[maskv]
        axes[1].plot(lag[maskv][valid], ampl[valid], color=color, lw=1.15)
        axes[2].plot(lag[maska], am[maska], color=color, lw=1.15)
        axes[2].fill_between(lag[maska], am[maska]-ass[maska], am[maska]+ass[maska], color=color, alpha=.16, lw=0)
    axes[0].axhline(0, color="0.45", lw=0.8); axes[0].set_xscale("log"); axes[0].set(xlim=(1,100), ylim=(-.03,.03), xlabel=r"$t$ (ps)", ylabel=r"$C_{vv}(t)/C_{vv}(0)$")
    axes[1].set_xscale("log"); axes[1].set_yscale("log"); axes[1].set(xlim=(1,100), ylim=(1e-6,1e-2), xlabel=r"$t$ (ps)", ylabel=r"$|C_{vv}(t)/C_{vv}(0)|$")
    axes[2].axhline(.5, color="0.45", lw=0.8, ls="--"); axes[2].set(xlim=(0,100), ylim=(.35,1.08), xlabel=r"$t$ (ps)", ylabel=r"$\alpha_{\mathrm{VACF\!\!-\!ODE}}(t)$")
    for tag, ax in zip(("(a)", "(b)", "(c)"), axes): style(ax); panel_label(ax, tag)
    axes[0].legend(frameon=False, ncol=1, fontsize=6.5, handlelength=1.5, loc="upper right")
    fig.savefig(figdir / "VACF_alpha_ODE_8_8_L2L5L10_8rep_aggregate.png", dpi=600)
    fig.savefig(figdir / "VACF_alpha_ODE_8_8_L2L5L10_8rep_aggregate.pdf")
    plt.close(fig)
    # Individual three-panel cards, one per available length.
    for label, (lag, vm, vs, _, am, ass) in processed.items():
        fig = plt.figure(figsize=(7.0, 2.35), dpi=300)
        ax0, ax1, ax2 = fig.add_axes([.09,.24,.25,.66]), fig.add_axes([.41,.24,.25,.66]), fig.add_axes([.73,.24,.24,.66])
        maskv, maska = lag >= 1, lag >= .1
        ax0.plot(lag[maskv], vm[maskv], color=COLORS[label], lw=1.2); ax0.fill_between(lag[maskv],vm[maskv]-vs[maskv],vm[maskv]+vs[maskv],color=COLORS[label],alpha=.18,lw=0); ax0.axhline(0,color='.45',lw=.8)
        valid=np.abs(vm[maskv]) > 3.0*vs[maskv]; ax1.plot(lag[maskv][valid],np.abs(vm[maskv][valid]),color=COLORS[label],lw=1.2)
        ax2.plot(lag[maska], am[maska], color=COLORS[label], lw=1.2); ax2.fill_between(lag[maska],am[maska]-ass[maska],am[maska]+ass[maska],color=COLORS[label],alpha=.18,lw=0); ax2.axhline(.5,color='.45',lw=.8,ls='--')
        ax0.set_xscale('log'); ax0.set(xlim=(1,100),ylim=(-.03,.03),xlabel=r"$t$ (ps)",ylabel=r"$C_{vv}(t)/C_{vv}(0)$")
        ax1.set_xscale('log'); ax1.set_yscale('log'); ax1.set(xlim=(1,100),ylim=(1e-6,1e-2),xlabel=r"$t$ (ps)",ylabel=r"$|C_{vv}/C_{vv}(0)|$")
        ax2.set(xlim=(0,100),ylim=(.35,1.08),xlabel=r"$t$ (ps)",ylabel=r"$\alpha_{\mathrm{VACF\!\!-\!ODE}}(t)$")
        for tag, ax in zip(("(a)","(b)","(c)"),(ax0,ax1,ax2)): style(ax); panel_label(ax,tag)
        fig.text(.5,.965,rf"$(8,8)$, {label}, $n=8$",ha="center",va="top",fontsize=7)
        fig.savefig(figdir / f"VACF_alpha_ODE_8_8_{label}_8rep.png",dpi=600); fig.savefig(figdir / f"VACF_alpha_ODE_8_8_{label}_8rep.pdf"); plt.close(fig)
    (OUT / "source_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "README.md").write_text("""# (8,8) 2L--10L common-protocol VACF and VACF-ODE alpha

Available common-protocol lengths are **2, 3, 4, 5 and 10 L** only; 6--9 L are absent and are deliberately not interpolated.  All shown curves are axial oxygen self VACFs after instantaneous oxygen-COM subtraction, 10-fs dump cadence, 1-ns duration and eight velocity-seed replicas.  The common displayed window is 0--100 ps.

For each replica, `I(t)=integral Cvv`, `J(t)=integral I`, and `alpha(t)=t I/J` are evaluated with the trapezoidal rule before forming mean and replica SEM.  Thus alpha has the same uncertainty unit as the independent replica means, not a pointwise time-origin SEM.  Near zero lag alpha is formally ill-conditioned; plots and interpretation start at 0.1 ps.

`alpha_ODE_minima_5to100ps_per_replica.csv` records the 5--100 ps minimum separately for every replica and for the mean curve.  The minimum depth is a ratio of two signed VACF integrals and is not a monotonic finite-size observable by itself; use its replica scatter and the minimum-time trend before interpreting it physically.

The 2--5 L input is the newly merged 8-rep archive; 10L is its earlier independently complete 8-rep matching-protocol package.  Full provenance and numerical summaries are in `source_manifest.json`; the plotting/recalculation source is `scripts/build_vacf_alpha_88_L2L10_8rep_figures.py`.
""",encoding="utf-8")


if __name__ == "__main__":
    main()
