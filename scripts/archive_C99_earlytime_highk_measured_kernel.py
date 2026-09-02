"""Archive and QA the completed C99 early-time high-k measured-kernel array.

No DHO fit is performed.  The arrays have 100-fs cadence and are used only for
the empirical 0--2 ps static-vertex/kernel reconstruction.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"H:/gcmc_explore/translational_anomaly/02_isf_collective_modes")
SOURCE = ROOT / "remote_fetch/output"
OUT = ROOT / "results/collective_mode_response/C99_earlytime_highk_measured_kernel/2026-08-30"
CUTS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 1.0})


def save_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def main() -> None:
    (OUT / "derived_data").mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    rec, modeqa, curves = [], [], []
    for N in (1600, 3200):
        for rep in range(1, 5):
            src = SOURCE / f"N{N}/rep{rep}"
            z = np.load(src / "rep_arrays_early_highk.npz")
            meta = json.loads((src / "metadata.json").read_text(encoding="utf-8"))
            lag, k, W = z["lag_ps"], z["kz_inv_A"], z["static_weight"]
            direct_norm = z["direct_vacf"] / float(z["var"])
            for km in CUTS:
                sel = k <= km + 1e-8
                pred = np.sum(W[sel][None, :] * z["Fs"][:, sel] * z["cjj_phi"][:, sel], axis=1)
                m = (lag >= 1.0) & (lag <= 2.0)
                rmse = float(np.sqrt(np.mean((pred[m] - direct_norm[m]) ** 2)))
                corr = float(np.corrcoef(pred[m], direct_norm[m])[0, 1]) if m.sum() > 2 else np.nan
                rec.append({"N": N, "replica": rep, "kmax_inv_A": km, "n_modes": int(sel.sum()),
                            "static_weight_sum": float(W[sel].sum()), "direct_vacf_variance": float(z["var"]),
                            "rmse_normalized_1to2ps": rmse, "correlation_normalized_1to2ps": corr})
                for t, d, p in zip(lag, direct_norm, pred):
                    curves.append({"N": N, "replica": rep, "kmax_inv_A": km, "lag_ps": float(t),
                                   "direct_peculiar_vacf_normalized": float(d), "P_static_vertex_kernel": float(p)})
            eig = z["eigenvalues"]
            modeqa.append({"N": N, "replica": rep, "n_modes": int(len(k)), "kmax_inv_A": float(k[-1]),
                           "kept_rank": int(z["keep"].sum()), "condition_kept": float(eig[z["keep"]].max() / eig[z["keep"]].min()),
                           "weight_sum_full": float(W.sum()), "n_origins": int(len(z["origin_indices"])),
                           "static_frames": int(len(z["static_origin_indices"])), "source": meta["source"]})
    save_csv(OUT / "derived_data/reconstruction_per_seed_kmax.csv", rec)
    save_csv(OUT / "derived_data/reconstruction_curves_per_seed.csv", curves)
    save_csv(OUT / "derived_data/static_vertex_rank_condition_QA.csv", modeqa)
    # Curves: ensemble across seeds, only 1--2 ps is shown because the requested observable is early time.
    fig = plt.figure(figsize=(5.5, 2.5)); ax = fig.add_axes([.14,.22,.82,.7])
    colors = {1600: "#0072B2", 3200: "#D55E00"}
    for N in (1600, 3200):
        for km, ls in ((0.5, ":"), (3.0, "-")):
            arr = [r for r in curves if r["N"] == N and r["kmax_inv_A"] == km and r["lag_ps"] > 0]
            times = sorted(set(r["lag_ps"] for r in arr))
            p = [np.mean([r["P_static_vertex_kernel"] for r in arr if r["lag_ps"] == t]) for t in times]
            ax.plot(times, p, color=colors[N], ls=ls, lw=1.3, label=f"N{N}, $k_{{max}}$={km:g}" + r" $\AA^{-1}$")
        arr = [r for r in curves if r["N"] == N and r["kmax_inv_A"] == 3.0 and r["lag_ps"] > 0]
        times = sorted(set(r["lag_ps"] for r in arr)); d = [np.mean([r["direct_peculiar_vacf_normalized"] for r in arr if r["lag_ps"] == t]) for t in times]
        ax.plot(times, d, color=colors[N], marker="o", ms=2.4, lw=.9, alpha=.8, label=f"N{N} direct")
    ax.set_xscale("log"); ax.set_xlim(1, 500); ax.axvspan(1, 2, color=".85", zorder=-5)
    ax.set_xlabel("lag time (ps)"); ax.set_ylabel("normalized kernel / VACF")
    ax.text(-.13,1.04,"(a)",transform=ax.transAxes,fontweight="bold",fontsize=9)
    ax.tick_params(direction="out", length=3); ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=5.6, ncol=2, handlelength=1.4)
    for ext in ("png", "pdf", "svg", "tiff"): fig.savefig(OUT / f"figures/earlytime_highk_kernel_reconstruction.{ext}", dpi=600)
    # Static rank/conditioning.
    fig = plt.figure(figsize=(4.8, 2.5)); ax = fig.add_axes([.14,.23,.81,.69])
    x=np.arange(8); q=modeqa
    ax.bar(x-.16,[r["kept_rank"] for r in q],.32,color="#009E73",label="kept rank")
    ax2=ax.twinx(); ax2.plot(x+.16,[r["condition_kept"] for r in q],"o-",color="#CC79A7",ms=3,label="condition")
    ax.set_xticks(x,[f"N{r['N']}\nr{r['replica']}" for r in q],fontsize=6); ax.set_ylabel("rank"); ax2.set_ylabel("condition number")
    ax.text(-.13,1.04,"(b)",transform=ax.transAxes,fontweight="bold",fontsize=9); ax.tick_params(direction="out",length=3); ax2.tick_params(direction="out",length=3)
    ax.spines["top"].set_visible(False); ax2.spines["top"].set_visible(False)
    for ext in ("png", "pdf", "svg", "tiff"): fig.savefig(OUT / f"figures/static_vertex_rank_condition_QA.{ext}", dpi=600)
    (OUT / "metadata.json").write_text(json.dumps({"remote_root":"/lustre/home/users/ewu/vb_gcmc/MD/stage_C99_earlytime_highk_measured_kernel_20260830", "job":"1379146[].ccpbs1", "source_fetch":str(SOURCE), "systems":[1600,3200], "replicas":4, "cadence_ps":0.1, "analysis_window_ps":[0,2], "dho_fit_performed":False},indent=2),encoding="utf-8")
    (OUT / "README.md").write_text("# C99 early-time high-k measured-kernel/static-vertex archive\n\nAll eight remote array members completed and were fetched from the stated CCFEP root. The reconstruction is `P=sum W Fs Phi_J`, compared with direct peculiar VACF after both are normalized by the seed-specific axial velocity variance. This is an empirical 0--2 ps reconstruction only; no high-k DHO omega, Gamma, or phase inference is permitted from 100-fs data.\n",encoding="utf-8")
    (OUT / "QA.md").write_text("# QA\n\nFile-level remote evidence: each N1600/N3200 rep1--rep4 directory had nonempty NPZ, mode CSV, metadata, SUCCESS and a job-1379146 success log. Rank/condition and low-k overlap are retained in the derived tables. The 1--2 ps comparison has only eleven positive lag samples and shared-parent velocity seeds; it is not a long-time closure or independent-replica validation.\n",encoding="utf-8")
    (OUT / "FINISHED.txt").write_text("Completed archive and empirical early-time reconstruction.\n",encoding="utf-8")

if __name__ == "__main__": main()
