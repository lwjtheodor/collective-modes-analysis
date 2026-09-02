"""Low-frequency circumferential-current linewidths from an implicit-CNT dump.

The source trajectory contains every water atom.  Oxygen atoms (type 1) are
used consistently with the explicit-CNT TA_theta definition.  The 20 ps dump
cadence restricts the analysis to the zero-frequency part of the spectrum.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def read_currents(path: Path, nmax: int):
    times, currents, lzs = [], [], []
    with path.open("r", encoding="utf-8") as fh:
        while True:
            tag = fh.readline()
            if not tag:
                break
            if tag.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"Unexpected dump record: {tag!r}")
            step = int(fh.readline())
            if fh.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("Missing NUMBER OF ATOMS")
            natoms = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Missing BOX BOUNDS")
            bounds = np.array([[float(v) for v in fh.readline().split()] for _ in range(3)])
            header = fh.readline().split()[2:]
            required = ["type", "x", "y", "z", "vx", "vy"]
            if any(v not in header for v in required):
                raise ValueError(f"Required fields absent from {header}")
            vals = np.fromstring(" ".join(fh.readline() for _ in range(natoms)), sep=" ")
            arr = vals.reshape(natoms, len(header))
            col = {name: i for i, name in enumerate(header)}
            oxy = arr[:, col["type"]].astype(int) == 1
            a = arr[oxy]
            x, y, z = (a[:, col[q]] for q in ("x", "y", "z"))
            vx, vy = (a[:, col[q]] for q in ("vx", "vy"))
            r = np.hypot(x, y)
            if np.any(r <= 0):
                raise ValueError("Oxygen located on CNT axis; e_theta undefined")
            # Same instantaneous O-COM velocity subtraction used by the
            # explicit-CNT all-mode current analysis.
            vx = vx - vx.mean()
            vy = vy - vy.mean()
            vtheta = (-y / r) * vx + (x / r) * vy
            lz = bounds[2, 1] - bounds[2, 0]
            k = 2 * np.pi * np.arange(1, nmax + 1) / lz
            currents.append(np.exp(1j * np.outer(k, z)) @ vtheta)
            times.append(step * 0.0005)  # LAMMPS real units: 0.5 fs timestep
            lzs.append(lz)
    cur = np.asarray(currents, dtype=np.complex128).T
    return np.asarray(times), cur, float(np.mean(lzs)), int(a.shape[0])


def lorentzian(w, b, a, gamma):
    return b + a * gamma**2 / (w**2 + gamma**2)


def fit_gamma(w, s, omega_max):
    mask = (w >= 0.0) & (w <= omega_max)
    x, y = w[mask], s[mask]
    b0 = float(np.percentile(y[-max(5, len(y)//5):], 50))
    a0 = max(float(y[0] - b0), np.finfo(float).eps)
    p, cov = curve_fit(lorentzian, x, y, p0=(b0, a0, 0.02),
                       bounds=([0.0, 0.0, 0.00005], [np.inf, np.inf, omega_max]),
                       maxfev=30000)
    pred = lorentzian(x, *p)
    r2 = 1 - np.sum((y - pred)**2) / np.sum((y - y.mean())**2)
    return p, cov, r2, mask


def segment_psd(series, dt_ps, nperseg):
    nseg = len(series) // nperseg
    if nseg < 3:
        raise ValueError("Need at least three non-overlapping Welch segments")
    win = np.hanning(nperseg)
    norm = np.sum(win**2)
    all_psd = []
    for j in range(nseg):
        q = series[j*nperseg:(j+1)*nperseg]
        q = q - q.mean()
        all_psd.append(np.abs(np.fft.fft(q * win))**2 / norm)
    freq = np.fft.fftfreq(nperseg, dt_ps)
    positive = freq >= 0
    return freq[positive] * 2 * np.pi, np.asarray(all_psd)[:, positive]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--nmax", type=int, default=10)
    ap.add_argument("--nperseg", type=int, default=256)
    ap.add_argument("--omega-fit-max", type=float, default=0.12)
    ap.add_argument("--primary-nmax", type=int, default=6)
    ap.add_argument("--temperature-K", type=float, default=None)
    ap.add_argument("--case-label", default="implicit CNT")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "derived_data").mkdir(exist_ok=True)
    (args.output / "figures").mkdir(exist_ok=True)

    t, current, Lz, n_oxygen = read_currents(args.dump, args.nmax)
    dt = float(np.median(np.diff(t)))
    records, spectra = [], []
    ncol = 5
    nrow = int(np.ceil(args.nmax / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 2.8*nrow), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for n, ax in enumerate(axes, start=1):
        if n > args.nmax:
            ax.set_visible(False)
            continue
        w, psd_segments = segment_psd(current[n-1], dt, args.nperseg)
        mean = psd_segments.mean(axis=0)
        sem = psd_segments.std(axis=0, ddof=1) / np.sqrt(len(psd_segments))
        p, cov, r2, mask = fit_gamma(w, mean, args.omega_fit_max)
        gammas = []
        for seg in psd_segments:
            try:
                gammas.append(fit_gamma(w, seg, args.omega_fit_max)[0][2])
            except RuntimeError:
                pass
        gamma_sem = float(np.std(gammas, ddof=1) / np.sqrt(len(gammas))) if len(gammas) > 1 else np.nan
        records.append({"n": n, "k_Ainv": 2*np.pi*n/Lz, "gamma_rad_ps": p[2],
                        "gamma_segment_sem_rad_ps": gamma_sem, "fit_R2": r2,
                        "n_segments": len(psd_segments), "frequency_resolution_rad_ps": w[1]})
        spectra.extend({"n":n, "k_Ainv":2*np.pi*n/Lz, "omega_rad_ps":ww,
                        "S_mean":ss, "S_segment_sem":ee}
                       for ww, ss, ee in zip(w, mean, sem))
        ax.semilogy(w[mask], mean[mask], color="#2166ac", lw=1.2)
        ax.semilogy(w[mask], lorentzian(w[mask], *p), color="#b2182b", lw=1.2)
        ax.set_title(f"n={n}, k={2*np.pi*n/Lz:.4f} A$^{{-1}}$")
        ax.grid(alpha=.25)
    fig.supxlabel(r"$\omega$ (rad ps$^{-1}$)")
    fig.supylabel(r"$S_{J_\theta J_\theta}$ (arb.)")
    fig.suptitle("Implicit CNT: circumferential zero-frequency spectra and Lorentzian fits")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(args.output / "figures" / f"implicitCNT_TAtheta_zeropeak_spectra_n001_n{args.nmax:03d}.{ext}", dpi=300)
    plt.close(fig)

    fitdf = pd.DataFrame(records)
    specdf = pd.DataFrame(spectra)
    fitdf.to_csv(args.output / "derived_data" / f"TAtheta_zeropeak_lorentzian_fits_n001_n{args.nmax:03d}.csv", index=False)
    specdf.to_csv(args.output / "derived_data" / f"TAtheta_zeropeak_spectra_n001_n{args.nmax:03d}.csv", index=False)
    # At n >= 7 the low-frequency Lorentzian no longer describes the spectrum
    # cleanly (the fixed 20 ps cadence leaves little dynamic range).  Retain
    # every raw fit above, but make the explicitly qualified n=1..6 result the
    # primary long-wavelength k-dependence fit.
    fitdf["used_primary_k_fit"] = (fitdf["n"] <= args.primary_nmax) & (fitdf["fit_R2"] >= 0.55)
    fitdf.to_csv(args.output / "derived_data" / f"TAtheta_zeropeak_lorentzian_fits_n001_n{args.nmax:03d}.csv", index=False)
    primary = fitdf[fitdf.used_primary_k_fit]
    k, g = primary.k_Ainv.to_numpy(), primary.gamma_rad_ps.to_numpy()
    # The fixed-friction form is the primary hypothesis.  Also preserve a
    # descriptive power law for sensitivity, without forcing an origin law.
    X = np.column_stack([np.ones_like(k), k**2])
    beta, *_ = np.linalg.lstsq(X, g, rcond=None)
    resid = g - X @ beta
    dof = len(k) - 2
    cov = (resid @ resid / dof) * np.linalg.inv(X.T @ X)
    free = np.polyfit(np.log(k), np.log(g), 1)
    pd.DataFrame([{"model":"Gamma0_plus_Dk2", "Gamma0_rad_ps":beta[0], "Gamma0_se":np.sqrt(cov[0,0]),
                   "D_rad_ps_A2":beta[1], "D_se":np.sqrt(cov[1,1])},
                  {"model":"C_k_alpha", "alpha":free[0], "C":np.exp(free[1])}]).to_csv(
        args.output / "derived_data" / "TAtheta_linewidth_k_dependence_fits.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.3, 4.5))
    ax.errorbar(fitdf.k_Ainv, fitdf.gamma_rad_ps, yerr=fitdf.gamma_segment_sem_rad_ps, fmt="o", mfc="white", mec="#777777", ecolor="#aaaaaa", capsize=3, label="all individual fits")
    ax.errorbar(k, g, yerr=primary.gamma_segment_sem_rad_ps, fmt="o", color="#2166ac", capsize=3, label="primary low-k fits (n=1..6)")
    kk = np.linspace(0, k.max()*1.05, 300)
    ax.plot(kk, beta[0]+beta[1]*kk**2, color="#b2182b", label=rf"$\Gamma_0+Dk^2$; $\Gamma_0={beta[0]:.4g}$")
    ax.set(xlabel=r"$k$ (A$^{-1}$)", ylabel=r"$\Gamma_\theta$ (rad ps$^{-1}$)", title=f"{args.case_label}: circumferential linewidth")
    ax.grid(alpha=.25); ax.legend(frameon=False); fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(args.output / "figures" / f"implicitCNT_TAtheta_linewidth_k_dependence_n001_n{args.nmax:03d}.{ext}", dpi=300)
    plt.close(fig)
    metadata = {"source_dump": str(args.dump), "wall_model":"implicit CNT (fix cnt/field)",
                "case_label":args.case_label, "temperature_K":args.temperature_K, "oxygen_type":1, "n_oxygen":n_oxygen, "n_frames":len(t),
                "duration_ps":float(t[-1]-t[0]), "dt_ps":dt, "nyquist_rad_ps":float(np.pi/dt),
                "Lz_A":Lz, "nperseg_frames":args.nperseg, "n_segments":len(t)//args.nperseg,
                "primary_k_fit":f"n=1..{args.primary_nmax}, individual Lorentzian R2 >= 0.55; higher-n retained as raw sensitivity only",
                "current_definition":"sum_i [e_theta(r_i).(v_i-v_O_COM)] exp(i k_n z_i), i=oxygen"}
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (args.output / "FINISHED.txt").write_text("Implicit-CNT TA_theta linewidth analysis finished successfully.\n", encoding="utf-8")

if __name__ == "__main__":
    main()
