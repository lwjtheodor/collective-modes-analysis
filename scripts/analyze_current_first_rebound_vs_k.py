"""Compare absolute 5L axial current modes and their first positive rebound.

The stored total tagged-current sum is normalized.  Each curve is restored to
the absolute per-water current autocorrelation using cvj0 from summary.json:
    C_J_abs(k,t) = cvj0(k) * C_vJ_total(k,t).
We combine the separately computed n=1 and n=2..10 results from the same
100-fs, 1-ns 5L window.  Peak extraction uses a 0.5-ps Savitzky-Golay
display/feature filter; the exported CSV includes raw-curve values at the
selected peak indices.
"""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter

ROOT = Path(r"H:\gcmc_explore")
SRC = ROOT / "analysis" / "lowk_damping_100fs_8_8_2L5L_20260729" / "output"
OUT = ROOT / "translational_anomaly" / "02_isf_collective_modes" / "results"
ASSET = ROOT / "translational_anomaly" / "02_isf_collective_modes" / "assets"
OUT.mkdir(parents=True, exist_ok=True); ASSET.mkdir(parents=True, exist_ok=True)

def read_rows(path):
    with open(path, newline="") as h:
        return list(csv.DictReader(h))

def load_group(folder, wanted):
    base = SRC / folder
    meta = json.loads((base / "summary.json").read_text())
    by_n = {int(q["n"]): q for q in meta["mode_summary"]}
    dt = float(meta["dt_ps"])
    rows = read_rows(base / "tagged_current_mode_coupling.csv")
    result = {}
    for n in wanted:
        rr = [r for r in rows if int(float(r["n"])) == n]
        if not rr: continue
        rr.sort(key=lambda r: int(r["lag_index"]))
        t = np.array([int(r["lag_index"])*dt for r in rr])
        y = np.array([float(r["C_vJ_total"]) for r in rr])
        result[n] = dict(t=t, y=y, c0=float(by_n[n]["cvj0"]), k=float(rr[0]["k_inv_A"]), lam=float(rr[0]["lambda_A"]))
    return result, dt

def first_rebound(t, y, dt):
    # 0.5 ps is short compared with the first lobe but suppresses point noise.
    win = max(5, int(round(0.5/dt)) | 1)
    ys = savgol_filter(y, win, 2, mode="interp")
    use = np.where((t >= max(0.2, 2*dt)) & (t <= 40))[0]
    troughs, _ = find_peaks(-ys[use], prominence=0.025)
    if len(troughs): imin = use[troughs[0]]
    else: imin = use[np.argmin(ys[use])]
    after = np.where((np.arange(len(t)) > imin) & (t <= 45))[0]
    peaks, _ = find_peaks(ys[after], prominence=0.008)
    candidates = after[peaks]
    positive = candidates[ys[candidates] > 0]
    if len(positive): ip = positive[0]
    elif len(candidates): ip = candidates[0]
    else:
        end = min(len(t), imin + max(5, int(round(2.0*t[imin]/dt))))
        ip = imin + 1 + int(np.argmax(ys[imin+1:end]))
    return ys, imin, ip

curves = {}
d1, dt1 = load_group("5L_rep1", {1}); curves.update(d1)
d2, dt2 = load_group("5L_n2n10", set(range(2,11))); curves.update(d2)
assert abs(dt1-dt2) < 1e-12

metrics=[]
for n in sorted(curves):
    d=curves[n]; ys,imin,ip=first_rebound(d["t"],d["y"],dt1)
    d.update(ys=ys,imin=imin,ip=ip)
    r=float(ys[ip]); raw=float(d["y"][ip]); peak_abs=d["c0"]*r
    metrics.append(dict(n=n,k_inv_A=d["k"],lambda_A=d["lam"],
                        first_minimum_ps=float(d["t"][imin]),
                        first_minimum_ratio=float(ys[imin]),
                        first_rebound_ps=float(d["t"][ip]),
                        CJ0_abs_per_water_A2_ps2=d["c0"],
                        first_rebound_abs_A2_ps2=peak_abs,
                        first_rebound_abs_raw_curve_A2_ps2=d["c0"]*raw))

with open(OUT/"5L_n1n10_current_first_rebound_vs_k.csv","w",newline="") as h:
    w=csv.DictWriter(h,fieldnames=list(metrics[0])); w.writeheader(); w.writerows(metrics)

plt.rcParams.update({"font.family":"Arial","font.size":10,"axes.linewidth":1.0,"lines.linewidth":1.8,"savefig.dpi":300})
fig,axs=plt.subplots(1,2,figsize=(10.8,3.9),constrained_layout=True)
cmap=plt.get_cmap("viridis")
for j,n in enumerate(sorted(curves)):
    d=curves[n]; co=cmap(j/9); use=d["t"]<=35
    yabs=d["c0"]*d["ys"]
    axs[0].plot(d["t"][use],yabs[use],color=co,label=fr"$n={n}$, $\lambda={d['lam']:.0f}$ Å")
    axs[0].plot(d["t"][d["ip"]],yabs[d["ip"]],"o",ms=3.5,color=co)
axs[0].axhline(0,color="0.55",lw=.8)
axs[0].set(xlim=(0,35),xlabel=r"$\Delta t$ (ps)",ylabel=r"$C_J^{\rm abs}(k,t)/N$ ($\mathrm{Å}^2$ ps$^{-2}$)")
axs[0].legend(frameon=False,ncol=2,fontsize=7,loc="upper right")

m=metrics; k=np.array([r["k_inv_A"] for r in m]); peak=np.array([r["first_rebound_abs_A2_ps2"] for r in m]); c0=np.array([r["CJ0_abs_per_water_A2_ps2"] for r in m])
axs[1].plot(k,peak,"o-",color="#244F73",label=r"first rebound $C_J^{\rm abs}(k,t_{p1})/N$")
fit=np.polyfit(np.log(k[:5]),np.log(peak[:5]),1); kk=np.linspace(k[0],k[4],100)
axs[1].plot(kk,np.exp(fit[1])*kk**fit[0],":",color="#244F73",label=fr"low-$k$: peak $\propto k^{{{fit[0]:.2f}}}$")
axs[1].plot(k,c0,"s--",color="#B74361",label=r"initial $C_J^{\rm abs}(k,0)/N$")
axs[1].set(xlabel=r"$k$ ($\mathrm{Å}^{-1}$)",ylabel=r"absolute current correlation ($\mathrm{Å}^2$ ps$^{-2}$)",ylim=(0,1.15*np.nanmax(c0)))
axs[1].legend(frameon=False,fontsize=8,loc="best")
for ax,lbl in zip(axs,["(a)","(b)"]):
    ax.spines[["top","right"]].set_visible(False); ax.text(-.14,1.04,lbl,transform=ax.transAxes,fontweight="bold",fontsize=12)
fig.savefig(ASSET/"5L_n1n10_absolute_current_modes_and_first_rebound_vs_k.png",bbox_inches="tight")
fig.savefig(ASSET/"5L_n1n10_absolute_current_modes_and_first_rebound_vs_k.pdf",bbox_inches="tight")
print(OUT/"5L_n1n10_current_first_rebound_vs_k.csv")
print(ASSET/"5L_n1n10_absolute_current_modes_and_first_rebound_vs_k.png")
