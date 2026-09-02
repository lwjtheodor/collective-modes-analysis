"""Analyze the k=0 axial oxygen current from LAMMPS id/z/vz dumps.

At k=0 J_z is proportional to the oxygen-group axial momentum.  This is an
externally damped zero mode in a fixed-wall NVT run, not a propagating LA mode.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def read_j0(path: Path):
    steps, currents = [], []
    with path.open("r", encoding="utf-8") as f:
        while True:
            marker = f.readline()
            if not marker:
                break
            if marker.strip() != "ITEM: TIMESTEP":
                raise ValueError(f"bad marker in {path}: {marker!r}")
            step = int(f.readline())
            f.readline(); nat = int(f.readline())
            f.readline(); f.readline(); f.readline(); f.readline()
            header = f.readline().split()[2:]
            if header != ["id", "z", "vz"]:
                raise ValueError(f"unexpected fields: {header}")
            a = np.fromstring(" ".join(f.readline() for _ in range(nat)), sep=" ").reshape(nat, 3)
            steps.append(step); currents.append(a[:, 2].sum())
    t = np.asarray(steps, dtype=float) * 0.0005  # real-units LAMMPS timestep is 0.5 fs
    return t, np.asarray(currents, dtype=float)


def acf(x):
    x = x - x.mean(); n = len(x)
    y = np.fft.irfft(np.abs(np.fft.rfft(x, 2*n))**2, 2*n)[:n]
    y /= np.arange(n, 0, -1)
    return y / y[0]


def lorentz(w, background, amplitude, gamma):
    return background + amplitude * gamma**2 / (w**2 + gamma**2)


def exp_decay(t, amplitude, gamma, offset):
    return offset + amplitude * np.exp(-gamma * t)


def fit_gamma_spectrum(w, s, wmax):
    m = (w >= 0) & (w <= wmax)
    x, y = w[m], s[m]
    b = np.median(y[-max(8, len(y)//5):])
    p, _ = curve_fit(lorentz, x, y, p0=(b, max(y[0]-b, 1e-20), 0.03),
                     bounds=([0, 0, 1e-5], [np.inf, np.inf, wmax]), maxfev=50000)
    pred = lorentz(x, *p)
    r2 = 1 - np.sum((y-pred)**2) / np.sum((y-y.mean())**2)
    return p, r2, m


def fit_gamma_acf(t, c, tmax):
    m = (t >= 0) & (t <= tmax) & (c > -0.15)
    p, _ = curve_fit(exp_decay, t[m], c[m], p0=(1, 0.03, 0),
                     bounds=([0, 1e-5, -0.2], [1.5, 2, 0.2]), maxfev=50000)
    pred = exp_decay(t[m], *p)
    r2 = 1 - np.sum((c[m]-pred)**2) / np.sum((c[m]-c[m].mean())**2)
    return p, r2, m


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dumps", nargs="+", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--omega-fit-max", type=float, default=0.20)
    p.add_argument("--acf-fit-max-ps", type=float, default=80.0)
    a = p.parse_args()
    d, f = a.output/"derived_data", a.output/"figures"
    d.mkdir(parents=True, exist_ok=True); f.mkdir(exist_ok=True)
    traces, acfs, specs, rows, times = [], [], [], [], []
    for ir, src in enumerate(a.dumps, 1):
        t, j = read_j0(src); dt = float(np.median(np.diff(t))); c = acf(j)
        w = np.fft.rfftfreq(len(j), dt) * 2*np.pi
        win = np.hanning(len(j)); q = j-j.mean()
        s = np.abs(np.fft.rfft(q*win))**2 / np.sum(win**2)
        pp, r2s, _ = fit_gamma_spectrum(w, s, a.omega_fit_max)
        pa, r2a, _ = fit_gamma_acf(t-t[0], c, a.acf_fit_max_ps)
        traces.append(j); acfs.append(c); specs.append(s); times.append(t-t[0])
        rows.append({"replica":ir, "source_dump":str(src), "n_frames":len(t), "dt_ps":dt,
                     "frequency_resolution_rad_ps":w[1], "J0_mean_A_fs":j.mean(), "J0_std_A_fs":j.std(ddof=1),
                     "Gamma_spectrum_rad_ps":pp[2], "spectrum_fit_R2":r2s,
                     "Gamma_acf_rad_ps":pa[1], "acf_fit_R2":r2a})
    T = np.asarray(times); C = np.asarray(acfs); S = np.asarray(specs)
    # Same 10 fs/1 ns grid was verified for all replicas.
    w = np.fft.rfftfreq(C.shape[1], float(np.median(np.diff(T[0]))))*2*np.pi
    cm, cs = C.mean(0), C.std(0, ddof=1)/np.sqrt(len(C))
    sm, ss = S.mean(0), S.std(0, ddof=1)/np.sqrt(len(S))
    pp, r2s, ms = fit_gamma_spectrum(w, sm, a.omega_fit_max)
    pa, r2a, ma = fit_gamma_acf(T[0], cm, a.acf_fit_max_ps)
    per = pd.DataFrame(rows); per.to_csv(d/"k0_axial_current_per_replica_fits.csv", index=False)
    pd.DataFrame({"time_ps":T[0], "CJJ_mean":cm, "CJJ_replica_sem":cs}).to_csv(d/"k0_axial_current_CJJ_ensemble.csv", index=False)
    pd.DataFrame({"omega_rad_ps":w, "SJJ_mean":sm, "SJJ_replica_sem":ss}).to_csv(d/"k0_axial_current_SJJ_ensemble.csv", index=False)
    pd.DataFrame([{"observable":"J0=sum_O vz", "Gamma_spectrum_rad_ps":pp[2], "spectrum_fit_R2":r2s,
                   "Gamma_acf_rad_ps":pa[1], "acf_fit_R2":r2a, "n_replicas":len(C),
                   "duration_ps":T[0][-1], "frequency_resolution_rad_ps":w[1]}]).to_csv(d/"k0_axial_current_ensemble_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.5,4.4)); ax.plot(T[0][ma], cm[ma], color="#2166ac", label="ensemble CJJ")
    ax.fill_between(T[0][ma], cm[ma]-cs[ma], cm[ma]+cs[ma], color="#2166ac", alpha=.2)
    ax.plot(T[0][ma], exp_decay(T[0][ma], *pa), color="#b2182b", label=f"exp fit Γ={pa[1]:.4g}")
    ax.axhline(0,color="black",lw=.7); ax.set(xlabel="time (ps)",ylabel=r"$C_{J_0J_0}(t)/C(0)"); ax.grid(alpha=.2); ax.legend(); fig.tight_layout()
    fig.savefig(f/"k0_axial_current_CJJ.png",dpi=300); fig.savefig(f/"k0_axial_current_CJJ.pdf"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5,4.4)); ax.semilogy(w[ms],sm[ms],color="#2166ac",label="ensemble SJJ")
    ax.fill_between(w[ms],np.maximum(sm[ms]-ss[ms],1e-300),sm[ms]+ss[ms],color="#2166ac",alpha=.2)
    ax.semilogy(w[ms],lorentz(w[ms],*pp),color="#b2182b",label=f"Lorentz Γ={pp[2]:.4g}")
    ax.set(xlabel=r"$\omega$ (rad ps$^{-1}$)",ylabel=r"$S_{J_0J_0}$ (arb.)"); ax.grid(alpha=.2); ax.legend(); fig.tight_layout()
    fig.savefig(f/"k0_axial_current_SJJ.png",dpi=300); fig.savefig(f/"k0_axial_current_SJJ.pdf"); plt.close(fig)
    (a.output/"metadata.json").write_text(json.dumps({"definition":"J0(t)=sum over oxygen vz; proportional to oxygen axial momentum", "protocol_limit":"fixed-CNT weak-NH NVT; fitted width is not an intrinsic hydrodynamic damping constant", "dumps":[str(x) for x in a.dumps]},indent=2))
    (a.output/"FINISHED.txt").write_text("k=0 axial-current analysis finished successfully.\n")

if __name__ == "__main__": main()
