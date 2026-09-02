"""Nature-style VACF-tail morphology comparison from all available replicas.

Figure contract
---------------
Claim: the 5--100 ps axial VACF morphology is protocol-sensitive, while the
completed 20-ns baseline relaxes to noise-scale fluctuations rather than a
persistent negative tail.  Each panel uses all 1L--5L, n=3 trajectories;
shading is the SEM over trajectory means (not an independent block CI).
"""
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATASETS = [
    ("lowfreq", "Baseline NVT", "20 ns / replica; 1 ps frames; z-momentum removed every 5 ps"),
    ("highfreq", "Weak NH-NVT, no momentum removal",
     "1 ns rethermalization (T_damp = 10 fs) + 1 ns production (T_damp = 100 ps)\nO-atom frames every 10 fs; no momentum removal"),
]
COLORS = ["#163C66", "#2F6B9A", "#5797B8", "#86BACD", "#B6D8DF"]

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.7, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def curve_paths(folder: str, length: int):
    directory = ROOT / folder
    if folder == "lowfreq":
        return sorted(directory.glob(f"{length}L_*.csv"))
    return sorted(directory.glob(f"*L{length}_rep*.csv"))


def summarise(folder: str, length: int):
    rows = []
    for path in curve_paths(folder, length):
        with path.open(newline="") as handle:
            records = list(csv.DictReader(handle))
        lag = np.array([float(record["lag_ps"]) for record in records])
        rows.append(np.array([float(record["vacf_peculiar_mean"]) for record in records]))
    if len(rows) != 3:
        raise ValueError(f"{folder}, {length}L: expected three replicas, got {len(rows)}")
    mask = (lag >= 5.0) & (lag <= 100.0)
    stacked = np.asarray(rows)[:, mask]
    return lag[mask], stacked.mean(axis=0), stacked.std(axis=0, ddof=1) / np.sqrt(3)


def main():
    # 89-mm single-column width. Explicit positions prevent annotation overlap.
    fig = plt.figure(figsize=(3.50, 5.15), dpi=300)
    axes = [fig.add_axes([0.18, 0.62, 0.76, 0.22]),
            fig.add_axes([0.18, 0.15, 0.76, 0.22])]
    limits = [(-0.0045, 0.0035), (-0.008, 0.008)]
    lines = []
    for i, (ax, (folder, title, protocol)) in enumerate(zip(axes, DATASETS)):
        for length, color in zip(range(1, 6), COLORS):
            lag, mean, sem = summarise(folder, length)
            line, = ax.plot(lag, mean, color=color, lw=1.25, solid_capstyle="round", zorder=3)
            ax.fill_between(lag, mean - sem, mean + sem, color=color, alpha=0.13, lw=0, zorder=1)
            if i == 0:
                lines.append(line)
        ax.axhline(0.0, color="#3F3F3F", lw=0.75, zorder=2)
        ax.set(xlim=(5, 100), ylim=limits[i], xticks=[5, 25, 50, 75, 100],
               xlabel=r"Lag time, $t$ (ps)" if i == 1 else "",
               ylabel=r"$C_{v_z}^{\mathrm{pec}}(t)/C_{v_z}^{\mathrm{pec}}(0)$")
        ax.tick_params(length=3.0, pad=2)
        ax.text(-0.13, 1.20 if i == 0 else 1.17, f"({chr(97+i)})", transform=ax.transAxes,
                fontweight="bold", fontsize=8)
        ax.text(0.00, 1.20 if i == 0 else 1.17, title, transform=ax.transAxes, fontsize=7.2, fontweight="bold")
        ax.text(0.00, 1.045, protocol, transform=ax.transAxes, fontsize=5.7, color="#424242",
                linespacing=1.25)
        ax.text(0.985, 0.06, "n = 3 trajectories / length; band = SEM", ha="right",
                transform=ax.transAxes, fontsize=5.4, color="#424242")
    fig.legend(lines, [f"{length}L" for length in range(1, 6)], title="Box length",
               ncol=5, loc="upper center", bbox_to_anchor=(0.56, 0.995),
               frameon=False, handlelength=1.25, handletextpad=0.35, columnspacing=0.85,
               borderaxespad=0.0, fontsize=6.2, title_fontsize=6.2)
    out = ROOT / "assets" / "vacf_tail_morphology_protocols_nature"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=600)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
