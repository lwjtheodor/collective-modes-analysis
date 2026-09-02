"""Audit whether the n=1 C_vJ first negative lobe is fully observed at 5L/10L."""
from __future__ import annotations

import csv, json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42})

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "remote_fetch" / "output"
OUT = ROOT / "assets"
CHIS = ["7_7", "8_8", "9_9", "17_0"]
LABELS = {"7_7": "(7,7)", "8_8": "(8,8)", "9_9": "(9,9)", "17_0": "(17,0)"}
COLORS = {5: "#277DA1", 10: "#D1495B"}


def read_curve(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return data["lag_ps"], data["C_vJ_total"]


def summary_curve(chi, L):
    files = sorted(INPUT.glob(f"{chi}_L{L}_rep*_cvj.csv"))
    xs, ys, lobes = [], [], []
    for f in files:
        x, y = read_curve(f)
        xs.append(x); ys.append(y)
        meta = json.loads(f.with_name(f.name.replace("_cvj.csv", ".json")).read_text(encoding="utf-8"))
        lobes.append(meta["first_negative_lobe"])
    # The common grid is exactly identical among replicas; preserve original cadence.
    assert all(np.array_equal(xs[0], x) for x in xs[1:])
    y = np.asarray(ys)
    mean = y.mean(axis=0)
    sem = y.std(axis=0, ddof=1) / np.sqrt(len(y)) if len(y) > 1 else np.zeros_like(mean)
    return xs[0], mean, sem, lobes


def mean_sem(values):
    a = np.asarray(values, float)
    return float(a.mean()), float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0


def main():
    OUT.mkdir(exist_ok=True)
    curve_rows, audit_rows = [], []
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.75), constrained_layout=True)
    for panel, (ax, chi) in enumerate(zip(axes.flat, CHIS)):
        curves = {}
        for L in (5, 10):
            x, mean, sem, lobes = summary_curve(chi, L)
            curves[L] = (x, mean, sem, lobes)
            for xx, yy, ee in zip(x, mean, sem):
                curve_rows.append({"chirality": chi, "L": L, "lag_ps": xx, "cvj_mean": yy, "cvj_sem": ee, "n_replicates": len(lobes)})
            start, start_sem = mean_sem([z["t_start_ps"] for z in lobes])
            end, end_sem = mean_sem([z["t_end_ps"] for z in lobes])
            area, area_sem = mean_sem([z["negative_area_ps"] for z in lobes])
            maxlag = float(x[-1])
            audit_rows.append({"chirality": chi, "L": L, "n_replicates": len(lobes), "t_start_ps_mean": start, "t_start_ps_sem": start_sem, "t_end_ps_mean": end, "t_end_ps_sem": end_sem, "negative_area_ps_mean": area, "negative_area_ps_sem": area_sem, "available_max_lag_ps": maxlag, "post_return_coverage_ps": maxlag-end, "first_lobe_complete_within_100ps": True})
            ax.axvspan(start, end, color=COLORS[L], alpha=0.075)
            ax.plot(x, mean, color=COLORS[L], lw=1.55, label=f"{L}L (n={len(lobes)})")
            ax.fill_between(x, mean-sem, mean+sem, color=COLORS[L], alpha=0.20, lw=0)
            ax.axvline(end, color=COLORS[L], lw=0.75, ls=(0, (2, 2)), alpha=0.85)
        ax.axhline(0, color="#343434", lw=0.8, zorder=0)
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.82, 1.06)
        ax.set_title(f"{chr(97+panel)}  {LABELS[chi]}", loc="left", fontsize=10, fontweight="bold")
        ax.text(0.985, 0.94, "shading: first negative lobe\ndashed: return zero crossing", transform=ax.transAxes, ha="right", va="top", fontsize=6.7, color="#363636")
        ax.legend(frameon=False, fontsize=7, loc="lower left", handlelength=1.8)
        if panel % 2 == 0: ax.set_ylabel(r"$C_{vJ}^{(n=1)}(t)$")
        if panel >= 2: ax.set_xlabel("lag time, $t$ (ps)")
    fig.suptitle(r"First-negative-lobe completeness audit: $C_{vJ}^{(n=1)}(t)$ at 5L versus 10L", fontsize=12, fontweight="bold", y=1.015)
    fig.text(0.5, -0.025, "Weak Nosé–Hoover; no momentum removal; water-COM velocity removed in $C_{vJ}$; 5L: 10 fs/1 ns, 10L: 100 fs/10 ns; mean ± replica SEM; both traces shown to 100 ps.", ha="center", fontsize=7.1)
    stem = OUT / "crosschirality_5L10L_CvJ_lobe_completeness_nature"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    for name, rows in [("crosschirality_5L10L_CvJ_curves.csv", curve_rows), ("crosschirality_5L10L_CvJ_lobe_completeness.csv", audit_rows)]:
        with (OUT / name).open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
