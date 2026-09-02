#!/usr/bin/env python3
"""Publication-style visual summary for the matched-k 2L/3L ISF demo."""
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path("analysis")
OUT = ROOT / "isf_demo_8_8_2L_3L_20260719"
PAIRS = ((2, 3, r"$\lambda=100.8$ Å"), (4, 6, r"$\lambda=50.4$ Å"), (6, 9, r"$\lambda=33.6$ Å"))
COLORS = ("#0072B2", "#D55E00", "#009E73")


def read_curves(label):
    result = {}
    with open(ROOT / "outputs" / label / "intermediate_scattering_curves.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["n"]), row["kind"])
            result.setdefault(key, {"t": [], "y": []})
            result[key]["t"].append(float(row["lag_ps"]))
            result[key]["y"].append(float(row["real"]))
    for item in result.values():
        item["t"] = np.asarray(item["t"], dtype=float)
        item["y"] = np.asarray(item["y"], dtype=float)
    return result


def normalized(curve, max_t):
    keep = curve["t"] <= max_t
    return curve["t"][keep], curve["y"][keep] / curve["y"][0]


def load_kww():
    rows = list(csv.DictReader(open(OUT / "matched_k_isf_kww_summary.csv", newline="")))
    output = {}
    for row in rows:
        if row["system"] in ("2L", "3L") and row["kind"] == "F_self":
            output[(row["system"], int(row["n"]))] = {key: float(row[key]) for key in ("lambda_A", "beta", "tau_ps")}
    return output


def panel_label(fig, bbox, label):
    fig.text(bbox[0] - 0.025, bbox[1] + bbox[3] + 0.010, label,
             fontfamily="Arial", fontsize=9, fontweight="bold", va="bottom")


def main():
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7,
                         "axes.linewidth": 1.0, "xtick.direction": "out",
                         "ytick.direction": "out", "xtick.major.width": 1.0,
                         "ytick.major.width": 1.0, "lines.linewidth": 1.2})
    curves2, curves3, kww = read_curves("2L"), read_curves("3L"), load_kww()

    # Explicit 2x2 geometry: no automatic layout changes.
    fig = plt.figure(figsize=(5.5, 4.35), dpi=180)
    boxes = [(0.13, 0.59, 0.34, 0.31), (0.61, 0.59, 0.34, 0.31),
             (0.13, 0.14, 0.34, 0.31), (0.61, 0.14, 0.34, 0.31)]
    axes = [fig.add_axes(box) for box in boxes]
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(length=3, width=1.0)
    for bbox, label in zip(boxes, ("(a)", "(b)", "(c)", "(d)")):
        panel_label(fig, bbox, label)

    # (a) matched self ISF; line style distinguishes box length.
    ax = axes[0]
    for (n2, n3, name), color in zip(PAIRS, COLORS):
        for curves, n, style in ((curves2, n2, "-"), (curves3, n3, "--")):
            t, f = normalized(curves[(n, "F_self")], 2000.0)
            ax.plot(t, f, style, color=color)
        ax.plot([], [], color=color, label=name)
    ax.plot([], [], "k-", label="2L")
    ax.plot([], [], "k--", label="3L")
    ax.set_xscale("log")
    ax.set_xlim(1, 2000); ax.set_ylim(0, 1.03)
    ax.set_xlabel(r"$t$ (ps)"); ax.set_ylabel(r"$F_s(k,t)$")
    ax.legend(loc="lower left", frameon=False, fontsize=6, handlelength=1.6, ncol=2,
              columnspacing=0.9, handletextpad=0.35)

    # (b) total ISF: early distinct-dominated structural relaxation.
    ax = axes[1]
    for (n2, n3, name), color in zip(PAIRS, COLORS):
        for curves, n, style in ((curves2, n2, "-"), (curves3, n3, "--")):
            t, f = normalized(curves[(n, "F_total")], 20.0)
            ax.plot(t, f, style, color=color)
    ax.axhline(0, color="0.5", linewidth=1.0)
    ax.set_xlim(0, 20); ax.set_ylim(-0.12, 1.03)
    ax.set_xlabel(r"$t$ (ps)"); ax.set_ylabel(r"$F(k,t)/F(k,0)$")

    # (c) KWW shape parameter across the matched spatial scales.
    ax = axes[2]
    for label, marker, face in (("2L", "o", "white"), ("3L", "s", "#4D4D4D")):
        x, y = [], []
        for n2, n3, _ in PAIRS:
            n = n2 if label == "2L" else n3
            item = kww[(label, n)]
            x.append(item["lambda_A"]); y.append(item["beta"])
        ax.plot(x, y, color="0.35", linewidth=1.0)
        ax.scatter(x, y, s=28, marker=marker, facecolor=face, edgecolor="0.1", linewidth=1.0, label=label, zorder=3)
    ax.set_xscale("log"); ax.set_xlim(28, 125); ax.set_ylim(0.78, 1.04)
    ax.set_xlabel(r"$\lambda=2\pi/k$ (Å)"); ax.set_ylabel(r"$\beta_{\mathrm{KWW}}$ of $F_s$")
    ax.legend(loc="lower right", frameon=False, fontsize=6, handletextpad=0.35)

    # (d) KWW time: the reference line uses the central 2L point and lambda^2.
    ax = axes[3]
    ref_lam, ref_tau = None, None
    for label, marker, face in (("2L", "o", "white"), ("3L", "s", "#4D4D4D")):
        x, y = [], []
        for n2, n3, _ in PAIRS:
            n = n2 if label == "2L" else n3
            item = kww[(label, n)]
            x.append(item["lambda_A"]); y.append(item["tau_ps"])
            if label == "2L" and n == 4:
                ref_lam, ref_tau = item["lambda_A"], item["tau_ps"]
        ax.plot(x, y, color="0.35", linewidth=1.0)
        ax.scatter(x, y, s=28, marker=marker, facecolor=face, edgecolor="0.1", linewidth=1.0, label=label, zorder=3)
    xref = np.array([30.0, 120.0])
    ax.plot(xref, ref_tau * (xref / ref_lam) ** 2, color="#D55E00", linestyle=":", linewidth=1.2, label=r"$\tau\propto\lambda^2$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(28, 125); ax.set_ylim(70, 1800)
    ax.set_xlabel(r"$\lambda=2\pi/k$ (Å)"); ax.set_ylabel(r"$\tau_{\mathrm{KWW}}$ (ps)")
    ax.legend(loc="upper left", frameon=False, fontsize=6, handlelength=1.5, handletextpad=0.35)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "matched_k_isf_2L_3L.png", dpi=600)
    fig.savefig(OUT / "matched_k_isf_2L_3L.pdf")


if __name__ == "__main__":
    main()
