#!/usr/bin/env python3
"""Aggregate the C99 equal-time static-vertex block/stride audit.

This deliberately never estimates a spectrum: each audit.npz contains the
equal-time K,c,var vertex evaluated on disjoint time blocks.  The script makes
the sampling diagnosis reproducible and keeps N2400/rep4 excluded by design.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"H:/gcmc_explore/translational_anomaly/02_isf_collective_modes")
FETCH = ROOT / "remote_fetch/stage_C99_static_vertex_lowk_block_audit_20260829/output"
OUT = ROOT / "results/collective_mode_response/implicit_C99_unified_longitudinal_500ps/2026-08-29"
DERIVED = OUT / "derived_data"
FIG = OUT / "figures"
CASES = {200: range(1, 5), 400: range(1, 5), 800: range(1, 5),
         1600: range(1, 5), 2400: range(1, 4), 3200: range(1, 5)}


def mean_sem(values: np.ndarray) -> tuple[float, float, int]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, 0
    return float(values.mean()), float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0, len(values)


def main() -> None:
    DERIVED.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    found = []
    for N, reps in CASES.items():
        for rep in reps:
            p = FETCH / f"N{N}/rep{rep}/audit.npz"
            if not p.exists():
                continue
            x = np.load(p)
            names = [str(v) for v in x["columns"]]
            for r in x["rows"]:
                d = dict(zip(names, r))
                d.update(N=N, rep=rep, Lz_A=float(x["Lz_A"]), frames=int(x["frames"]))
                d["LzW_A"] = d["Lz_A"] * d["W"]
                rows.append(d)
            found.append(str(p))
    if len(found) != 23:
        raise RuntimeError(f"expected 23 audit members, found {len(found)}")

    fields = ["N", "rep", "Lz_A", "frames", "stride_frames", "block", "n", "k_inv_A", "W", "LzW_A", "K_condition", "samples"]
    with (DERIVED / "C99_static_vertex_lowk_block_audit_rows.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

    summaries = []
    # First average each replica over its 12 disjoint blocks; then report seed SEM.
    for N in CASES:
        for stride in (1, 10, 100):
            for n in range(1, 6):
                perrep = []
                for rep in CASES[N]:
                    sub = [r for r in rows if r["N"] == N and r["rep"] == rep and int(r["stride_frames"]) == stride and int(r["n"]) == n]
                    if sub:
                        vals = np.array([r["LzW_A"] for r in sub])
                        perrep.append((rep, vals.mean(), vals.std(ddof=1)))
                mu, sem, nr = mean_sem(np.array([v[1] for v in perrep]))
                block_sd_mu, block_sd_sem, _ = mean_sem(np.array([v[2] for v in perrep]))
                ref = next(r for r in rows if r["N"] == N and int(r["stride_frames"]) == stride and int(r["n"]) == n)
                summaries.append(dict(N=N, Lz_A=ref["Lz_A"], stride_frames=stride, cadence_ps=0.1*stride,
                                      n=n, k_inv_A=ref["k_inv_A"], LzW_A_mean=mu, LzW_A_seed_sem=sem,
                                      mean_block_sd=block_sd_mu, mean_block_sd_sem=block_sd_sem,
                                      n_reps=nr, blocks_per_rep=12, samples_per_block=ref["samples"],
                                      K_condition_mean=np.mean([r["K_condition"] for r in rows if r["N"] == N and int(r["stride_frames"]) == stride and int(r["n"]) == n])))
    sf = list(summaries[0].keys())
    with (DERIVED / "C99_static_vertex_lowk_block_audit_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sf); w.writeheader(); w.writerows(summaries)

    # Quantitative cadence invariance: compare 1 and 100 frame results within exact N, n, replica, block.
    diffs = []
    for N in CASES:
        for rep in CASES[N]:
            for n in range(1, 6):
                a = {(int(r["block"])): r["LzW_A"] for r in rows if r["N"] == N and r["rep"] == rep and int(r["stride_frames"]) == 1 and int(r["n"]) == n}
                b = {(int(r["block"])): r["LzW_A"] for r in rows if r["N"] == N and r["rep"] == rep and int(r["stride_frames"]) == 100 and int(r["n"]) == n}
                for q in sorted(a.keys() & b.keys()):
                    diffs.append(a[q] - b[q])
    dmean, dsem, dn = mean_sem(np.asarray(diffs))
    meta = {"source": str(FETCH), "members": len(found), "excluded": "N2400/rep4 permanently excluded",
            "definition": "equal-time K,c,variance vertex; no Welch/spectral window", "cadence_ps": [0.1, 1.0, 10.0],
            "stride_1_minus_100_LzW_A_mean": dmean, "stride_1_minus_100_LzW_A_sem": dsem, "paired_block_count": dn}
    (DERIVED / "C99_static_vertex_lowk_block_audit_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")

    mpl.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 7,
                         "axes.linewidth": 1.0, "xtick.direction": "out", "ytick.direction": "out", "legend.frameon": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(6.9, 2.65)); ax1 = fig.add_axes([0.09, 0.20, 0.39, 0.72]); ax2 = fig.add_axes([0.60, 0.20, 0.35, 0.72])
    colors = plt.cm.viridis(np.linspace(0.10, 0.90, 6))
    for color, N in zip(colors, CASES):
        dat = [s for s in summaries if s["N"] == N and s["stride_frames"] == 1]
        dat.sort(key=lambda z: z["n"])
        ax1.errorbar([d["n"] for d in dat], [d["LzW_A_mean"] for d in dat], yerr=[d["LzW_A_seed_sem"] for d in dat], marker="o", ms=3, lw=1.1, capsize=1.5, color=color, label=f"{dat[0]['Lz_A']/10:g} nm")
        low = next(d for d in dat if d["n"] == 1)
        ax2.plot(low["Lz_A"], low["mean_block_sd"], "o", ms=4, color=color)
    ax1.axhline(1.0, color="0.5", lw=0.9, ls="--"); ax1.set_xlabel("mode index $n$"); ax1.set_ylabel(r"$L_zW_n$ ($\mathrm{\AA}$)"); ax1.set_xticks(range(1, 6)); ax1.legend(title=r"$L_z$", ncol=2, loc="lower left", fontsize=6, title_fontsize=6)
    ax2.set_xlabel(r"$L_z$ ($\mathrm{\AA}$)"); ax2.set_ylabel(r"block SD of $L_zW_1$ ($\mathrm{\AA}$)")
    for ax, letter in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(-0.17, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=9); ax.spines[["top", "right"]].set_visible(False)
    base = FIG / "C99_static_vertex_lowk_block_stride_audit"
    fig.savefig(str(base) + ".png", dpi=600); fig.savefig(str(base) + ".pdf"); fig.savefig(str(base) + ".svg"); plt.close(fig)
    print(json.dumps(meta, indent=2))

if __name__ == "__main__": main()
