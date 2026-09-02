"""Render the first-look (8,8) axial ISF comparison for 2--5 L.

Input ensemble NPZ files are produced by rebuild_axial_isf.py.  The top row
compares each cell's fundamental axial wave vector; the lower row uses the
same physical wave vector (0.06230846 A^-1) across all lengths.
"""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "collective_mode_response" / "isf_total_self_distinct_88_L2L10_10fs" / "2026-08-20" / "analysis_200ps_100pslag"
OUT = DATA / "figures"
LENGTHS = (2, 3, 4, 5)
MODES = (("kmin", r"$k=k_{\min}=2\pi/L_z$"),
         ("matched_k", r"$k=0.06231\ \mathrm{\AA}^{-1}$"))
FIELDS = (("total", r"$F(k,t)$"), ("self", r"$F_{\mathrm{s}}(k,t)$"),
          ("distinct", r"$F_{\mathrm{d}}(k,t)$"))
COLORS = {2: "#2878B5", 3: "#E18727", 4: "#32A467", 5: "#C03D3E"}


def axes_bbox(fig, x, y, w, h, *, left=0.14, right=0.035, bottom=0.18, top=0.06):
    """Create an axis within a fixed panel outer box, reserving label space."""
    return fig.add_axes([x + left * w, y + bottom * h,
                         w * (1 - left - right), h * (1 - bottom - top)])


def load_mode(length, mode):
    with np.load(DATA / f"ISF_{length}_{mode}_mean_sem.npz") as z:
        return {key: z[key] for key in z.files}


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7,
        "axes.linewidth": 1.0, "xtick.major.width": 1.0,
        "ytick.major.width": 1.0, "xtick.direction": "out",
        "ytick.direction": "out", "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    # Explicit bbox-first geometry: 3 columns x 2 rows.
    width, height = 7.0, 4.55
    margin_left, margin_right, margin_bottom, margin_top = 0.075, 0.015, 0.075, 0.105
    col_gap, row_gap = 0.018, 0.058
    panel_w = (1 - margin_left - margin_right - 2 * col_gap) / 3
    panel_h = (1 - margin_top - margin_bottom - row_gap) / 2
    fig = plt.figure(figsize=(width, height), facecolor="white")
    caches = {mode: {L: load_mode(L, mode) for L in LENGTHS} for mode, _ in MODES}
    labels = "abcdef"
    qa = {"data_dir": str(DATA), "rows": []}

    for row, (mode, row_note) in enumerate(MODES):
        y = margin_bottom + (1 - row) * (panel_h + row_gap)
        for col, (field, ylabel) in enumerate(FIELDS):
            x = margin_left + col * (panel_w + col_gap)
            ax = axes_bbox(fig, x, y, panel_w, panel_h)
            for L in LENGTHS:
                d = caches[mode][L]
                t = d["time_ps"]
                mean, sem = d[f"F_{field}_mean"], d[f"F_{field}_sem"]
                ax.fill_between(t, mean-sem, mean+sem, color=COLORS[L], alpha=0.14, lw=0, clip_on=True)
                ax.plot(t, mean, color=COLORS[L], lw=1.15, label=fr"${L}L$")
                qa["rows"].append({"length_L": L, "mode": mode, "field": field,
                                   "k_Ainv": float(d["k_inv_A_mean"]),
                                   "F0": float(mean[0]), "F_100ps": float(mean[-1])})
            ax.axhline(0, color="0.45", lw=1.0, zorder=0)
            ax.set_xlim(0, 100)
            ax.set_xlabel(r"$t\ (\mathrm{ps})$")
            ax.set_ylabel(ylabel)
            ax.tick_params(top=False, right=False, length=3.0, pad=2)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            # fixed label outside the corresponding outer panel frame
            fig.text(x - 0.006, y + panel_h + 0.005, f"({labels[row*3+col]})",
                     fontsize=9, fontweight="bold", ha="left", va="bottom")
            if row == 0:
                ax.text(0.03, 0.92, row_note, transform=ax.transAxes, ha="left", va="top")
            if row == 0 and col == 0:
                ax.legend(loc="upper right", frameon=False, handlelength=1.8,
                          labelspacing=0.25, borderaxespad=0.25)

    fig.text(0.5, 0.985,
             r"$(8,8)$ water: axial intermediate scattering function; 3 replicas, 200 ps per replica",
             ha="center", va="top", fontsize=7)
    for suffix, dpi in (("png", 600), ("pdf", None)):
        target = OUT / f"ISF_88_L2L5_overview.{suffix}"
        kwargs = {"dpi": dpi} if dpi else {}
        fig.savefig(target, facecolor="white", **kwargs)
    plt.close(fig)
    (OUT / "ISF_88_L2L5_overview_qa.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
