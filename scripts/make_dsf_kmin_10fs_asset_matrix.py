"""Build a comparable 10-fs S(k_min, omega) panel matrix from verified assets.

The four source trajectories are fixed-CNT, 330-K high-frequency runs.  They
all use a 10-fs frame interval.  To avoid mixing spectral resolution, each
input is truncated to its shared 0--50 ps ISF window before a symmetric,
Hann-tapered Fourier transform is applied.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, NullFormatter


ROOT = Path(r"H:\gcmc_explore")
OUT = Path(__file__).resolve().parents[1] / "results" / "collective_mode_response" / "dsf_kmin_10fs_asset_matrix" / "2026-08-18"
FIG_STEM = OUT / "figures" / "dsf_kmin_10fs_asset_matrix"

# Build anew: each input is an n=1 total density ISF, normalized per water.
# List order is the visual ordering of the 2x2 matrix.
ASSETS = [
    {
        "chirality": "(7,7)",
        "lz_A": 201.679996,
        "tag": "(7,7), 2L",
        "csv": ROOT / "analysis" / "highfreq_chirality_7_7_15_0_20260720" / "analysis" / "7_7_2L" / "isf_vacf" / "intermediate_scattering_curves.csv",
        "summary": ROOT / "analysis" / "highfreq_chirality_7_7_15_0_20260720" / "analysis" / "7_7_2L" / "isf_vacf" / "summary.json",
    },
    {
        "chirality": "(15,0)",
        "lz_A": 204.48,
        "tag": "(15,0), 2L",
        "csv": ROOT / "analysis" / "highfreq_chirality_7_7_15_0_20260720" / "analysis" / "15_0_2L" / "isf_vacf" / "intermediate_scattering_curves.csv",
        "summary": ROOT / "analysis" / "highfreq_chirality_7_7_15_0_20260720" / "analysis" / "15_0_2L" / "isf_vacf" / "summary.json",
    },
    {
        "chirality": "(8,8)",
        "lz_A": 403.359992,
        "tag": "(8,8), 4L",
        "csv": ROOT / "analysis" / "highfreq_figures_8_8_2L5L_20260720" / "isf" / "4L_rep1" / "intermediate_scattering_curves.csv",
        "summary": ROOT / "analysis" / "highfreq_figures_8_8_2L5L_20260720" / "isf" / "4L_rep1" / "summary.json",
    },
    {
        "chirality": "(8,8)",
        "lz_A": 504.199990,
        "tag": "(8,8), 5L",
        "csv": ROOT / "analysis" / "highfreq_figures_8_8_2L5L_20260720" / "isf" / "5L_rep1" / "intermediate_scattering_curves.csv",
        "summary": ROOT / "analysis" / "highfreq_figures_8_8_2L5L_20260720" / "isf" / "5L_rep1" / "summary.json",
    },
]


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def dynamic_structure_factor(time_ps: np.ndarray, f_total: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return a one-sided, Hann-tapered cosine spectrum in ps units."""
    dt = float(np.median(np.diff(time_ps)))
    max_time = 50.0
    keep = time_ps <= max_time + 1.0e-10
    t = time_ps[keep]
    c = f_total[keep]
    if len(t) < 10 or not np.allclose(np.diff(t), dt, rtol=2e-4, atol=1e-12):
        raise ValueError("ISF lag grid is not a regular 10-fs grid")
    # C(-t)=C*(t): embed the real, even correlation around zero lag.  Padding
    # after ifftshift would create an artificial discontinuity and a comb-like
    # leakage pattern, so positive and negative lags are placed explicitly.
    nfft = 16384
    taper_positive = np.hanning(2 * len(c) - 1)[len(c) - 1:]
    tapered = c * taper_positive
    corr = np.zeros(nfft, dtype=float)
    corr[:len(tapered)] = tapered
    corr[-(len(tapered) - 1):] = tapered[1:][::-1]
    spectrum = np.fft.rfft(corr).real * dt
    omega = 2.0 * np.pi * np.fft.rfftfreq(nfft, d=dt)
    return omega, np.maximum(spectrum, 1.0e-8), float(t[-1])


def load_asset(asset: dict) -> tuple[dict, pd.DataFrame]:
    # Source tables contain many mode/kind combinations.  Stream only the one
    # requested curve so the matrix generator stays lightweight.
    rows = []
    with asset["csv"].open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["n"] == "1" and row["kind"] == "F_total":
                rows.append({key: float(row[key]) for key in ("lag_ps", "real", "k_inv_A")})
    curve = pd.DataFrame(rows).sort_values("lag_ps")
    if curve.empty:
        raise ValueError(f"no n=1 F_total curve in {asset['csv']}")
    summary = json.loads(asset["summary"].read_text(encoding="utf-8"))
    omega, spectrum, tmax = dynamic_structure_factor(
        curve["lag_ps"].to_numpy(float), curve["real"].to_numpy(float)
    )
    mask = omega <= 3.0
    result = dict(asset)
    result.update({
        "dt_ps": float(summary["dt_ps"]),
        "n_frames": int(summary["n_frames"]),
        "n_water": int(summary["n_water"]),
        "kmin_inv_A": float(curve["k_inv_A"].iloc[0]),
        "window_ps": tmax,
        "omega": omega[mask],
        "spectrum": spectrum[mask],
    })
    rows = pd.DataFrame({
        "chirality": asset["chirality"],
        "case": asset["tag"],
        "lz_A": asset["lz_A"],
        "kmin_inv_A": result["kmin_inv_A"],
        "dt_ps": result["dt_ps"],
        "window_ps": result["window_ps"],
        "omega_rad_ps": result["omega"],
        "S_kmin_omega_ps": result["spectrum"],
    })
    return result, rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG_STEM.parent.mkdir(parents=True, exist_ok=True)
    records, long_tables = zip(*(load_asset(asset) for asset in ASSETS))
    source_data = pd.concat(long_tables, ignore_index=True)
    source_data.to_csv(OUT / "derived_data_dsf_kmin_10fs.csv", index=False)

    positive = source_data["S_kmin_omega_ps"].to_numpy(float)
    ymin = max(float(np.nanmin(positive[positive > 0])) * 0.75, 1e-7)
    ymax = float(np.nanmax(positive)) * 1.25
    colors = ["#1f4e79", "#b35806", "#2c7a7b", "#7a5195"]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.55), sharex=True, sharey=True)
    for index, (ax, record, color) in enumerate(zip(axes.flat, records, colors)):
        ax.plot(record["omega"], record["spectrum"], color=color, lw=1.15)
        ax.set_yscale("log")
        ax.set_xlim(0.0, 3.0)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks([0, 1, 2, 3])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.2g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(length=3.0, width=1.0, pad=2)
        ax.text(-0.15, 1.04, f"({chr(97 + index)})", transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="bottom", ha="left")
        ax.text(0.97, 0.95, record["tag"], transform=ax.transAxes,
                ha="right", va="top", fontweight="bold", color=color, fontsize=8)
        ax.text(0.97, 0.84,
                f"Lz = {record['lz_A']:.1f} A;  kmin = {record['kmin_inv_A']:.4f} A^-1",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.4)
        ax.text(0.97, 0.74, f"Δt = 10 fs;  T = {record['window_ps']:.0f} ps;  n = 1",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.1, color="#444444")
    fig.supxlabel("angular frequency, omega (rad ps^-1)", y=0.05, fontsize=7.5)
    fig.supylabel("S(kmin, ω) (ps; log scale)", x=0.04, fontsize=7.5)
    fig.text(0.5, 0.012,
             "Fixed CNT; 330 K; S(kmin, ω) from a symmetric Hann-tapered transform of Ftotal(kmin, t) over the common 0–50 ps window.",
             ha="center", va="bottom", fontsize=6.0)
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.17, top=0.96, wspace=0.16, hspace=0.18)
    # The PNG is a compact screen preview; PDF/SVG retain the publication
    # resolution and editable text for this line-art figure.
    for suffix, kwargs in {
        ".png": {"dpi": 300},
        ".pdf": {},
        ".svg": {},
    }.items():
        fig.savefig(FIG_STEM.with_suffix(suffix), bbox_inches="tight", **kwargs)
    plt.close(fig)

    metadata = {
        "observable": "S(k_min, omega), from the real symmetric total-density ISF",
        "transform": "Hann-tapered 0-50 ps window; omega=2*pi*f; non-negative display floor 1e-8 ps",
        "comparability": "all rendered inputs have 10-fs ISF sampling and are truncated to the common 50-ps lag window",
        "replicate_note": "each displayed legacy asset is one trajectory; no SEM is available",
        "availability_note": "No (9,9) or (17,0) 10-fs total-ISF asset was available locally at figure creation; they are not imputed.",
        "assets": [{key: value for key, value in rec.items() if key not in {"omega", "spectrum", "csv", "summary"}} for rec in records],
    }
    (OUT / "figure_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    with (OUT / "QA_notes.md").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# QA notes\n\n")
        handle.write("- Archetype: quantitative 2x2 matrix with shared frequency and spectral-density axes.\n")
        handle.write("- Inputs: four completed 10-fs total-density ISF assets; one trajectory per cell.\n")
        handle.write("- Shared analysis: n=1; common 0-50 ps window; Hann taper; no row or cell was imputed.\n")
        handle.write("- Exclusions: the matrix excludes (9,9) and (17,0) because no local 10-fs total-ISF curves were present.\n")
        handle.write("- Interpretation: compare peak position, broadening, and low-frequency weight; do not treat the four single-trajectory curves as uncertainty-qualified chirality laws.\n")


if __name__ == "__main__":
    main()
