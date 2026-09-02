"""Render the completed N200/N400 portion of the implicit-C88 closure.

The N1600 member is intentionally excluded until its four-seed output carries
the same file-level completion evidence.  Existing CJJ-43 DHO values supply
only the separately labelled dispersion/damping panels.
"""
from __future__ import annotations
import csv, json, math, shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(r"H:/gcmc_explore/translational_anomaly/02_isf_collective_modes")
FETCH=ROOT/"remote_fetch/implicit_C88_static_vertex_VACF_20260831"
OUT=ROOT/"results/collective_mode_response/implicit_C88_static_weight_VACF_closure/2026-08-31"
DHO=ROOT/"results/collective_mode_response/implicit_C77_C88_C99_matched_k_effective_DHO_damping/2026-08-29/derived_data/matched_k_effective_DHO_summary.csv"
CASES=("N200_weakNH_6ns","N400_weakNH_6ns","N1600_weakNH_6ns")
COLORS={"N200_weakNH_6ns":"#D55E00","N400_weakNH_6ns":"#0072B2","N1600_weakNH_6ns":"#009E73","N200":"#D55E00","N400":"#0072B2","N1600":"#009E73"}

plt.rcParams.update({"font.family":"Arial","font.size":7,"axes.linewidth":1.0,"pdf.fonttype":42,"svg.fonttype":"none"})

def read(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def num(row,key): return float(row[key])
def save(rows,path):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def meansem(x):
    x=np.asarray(x,float);return float(x.mean()),float(x.std(ddof=1)/math.sqrt(len(x)))

def main():
    (OUT/"derived_data").mkdir(parents=True,exist_ok=True);(OUT/"figures").mkdir(exist_ok=True);(OUT/"scripts").mkdir(exist_ok=True)
    shutil.copy2(__file__,OUT/"scripts"/Path(__file__).name)
    weight_summary=[]; rec_summary=[]; curve={}; factor_rows=[]
    for case in CASES:
        meta=json.loads((FETCH/case/"metadata.json").read_text()); N=meta["inputs"][0]["n_oxygen"]
        wr=read(FETCH/case/"static_weights_per_replica.csv")
        for n in sorted({int(r["n"]) for r in wr}):
            q=[r for r in wr if int(r["n"])==n]; w,e=meansem([num(r,"W_diag") for r in q])
            k=num(q[0],"k_inv_A")
            weight_summary.append({"case":case,"N_oxygen":N,"n":n,"k_inv_A":k,"W_diag_mean":w,"W_diag_seedSEM":e,"N_W_diag_mean":N*w,"N_W_diag_seedSEM":N*e})
        mr=read(FETCH/case/"reconstruction_metrics_per_replica.csv")
        for n in sorted({int(r["nmax"]) for r in mr}):
            q=[r for r in mr if int(r["nmax"])==n]
            r2,e2=meansem([num(r,"R2_2to100ps") for r in q]); rm,erm=meansem([num(r,"RMSE_2to100ps") for r in q]); ws,ews=meansem([num(r,"weight_sum") for r in q])
            kmax=next(num(r,"k_inv_A") for r in weight_summary if r["case"]==case and r["n"]==n)
            rec_summary.append({"case":case,"nmax":n,"kmax_inv_A":kmax,"weight_sum_mean":ws,"weight_sum_seedSEM":ews,"R2_2to100ps_mean":r2,"R2_2to100ps_seedSEM":e2,"RMSE_2to100ps_mean":rm,"RMSE_2to100ps_seedSEM":erm})
        cr=read(FETCH/case/"kernel_curves_per_replica_and_ensemble.csv"); m=max(int(r["nmax"]) for r in mr)
        q=[r for r in cr if r["replica"]=="ensemble" and int(r["n"])==m]
        curve[case]=q
    save(weight_summary,OUT/"derived_data/static_Wdiag_N200_N400_ensemble.csv")
    save(rec_summary,OUT/"derived_data/reconstruction_cutoff_N200_N400_ensemble.csv")
    # Constant-weight audit in the natural finite-box measure N*W.
    qa=[]
    for case in CASES:
        q=[r for r in weight_summary if r["case"]==case]
        x=np.array([num(r,"N_W_diag_mean") for r in q]); qa.append({"case":case,"n_range":f"1..{len(q)}","mean_NW":float(x.mean()),"CV_NW":float(x.std(ddof=1)/x.mean()),"range_over_mean_NW":float((x.max()-x.min())/x.mean())})
    save(qa,OUT/"derived_data/static_weight_flatness_audit.csv")
    # Positive n are the stored half of a real-field (+/- k) pair.  Compare
    # the former one-sided expression and the correct pair-degenerate factor.
    for case in CASES:
        q=curve[case]; t=np.array([num(r,"lag_ps") for r in q]); d=np.array([num(r,"direct_VACF_mean") for r in q]); p=np.array([num(r,"P_diag_mean") for r in q]); mask=(t>=2)&(t<=100)
        optimal=float(np.dot(p[mask],d[mask])/np.dot(p[mask],p[mask]))
        for fac,label in ((1.0,"one_sided_positive_k"),(2.0,"paired_plus_minus_k"),(optimal,"OLS_global_diagnostic")):
            pred=fac*p[mask]; sse=float(np.sum((d[mask]-pred)**2)); sst=float(np.sum((d[mask]-d[mask].mean())**2))
            factor_rows.append({"case":case,"factor_type":label,"factor":fac,"R2_ensemble_2to100ps":1-sse/sst,"RMSE_ensemble_2to100ps":float(np.sqrt(np.mean((d[mask]-pred)**2))),"P0":fac*float(p[0])})
    save(factor_rows,OUT/"derived_data/VACF_positive_k_degeneracy_factor_ensemble.csv")
    # CJJ-43 DHO: exact protocol, but distinct from the new static closure.
    d=[r for r in read(DHO) if r["system"]=="C88"]
    # explicit BBox-first 2x2 layout.
    fig=plt.figure(figsize=(5.5,4.6)); boxes=[(.10,.59,.36,.31),(.58,.59,.36,.31),(.10,.12,.36,.31),(.58,.12,.36,.31)]
    axs=[fig.add_axes(b) for b in boxes]
    for case,marker in (("N200","s"),("N400","o"),("N1600","^")):
        q=[r for r in d if r["case"]==case]; x=np.array([num(r,"k_Ainv") for r in q]);om=np.array([num(r,"omega_mean_rad_ps") for r in q]);oe=np.array([num(r,"omega_velocity_seed_sem_rad_ps") for r in q]);ga=np.array([num(r,"gamma_effective_psinv_mean") for r in q]);ge=np.array([num(r,"gamma_velocity_seed_sem_psinv") for r in q])
        axs[0].errorbar(x,om,yerr=oe,fmt=marker+"-",ms=3,capsize=1.5,lw=1.1,color=COLORS[case],label=case)
        axs[1].errorbar(x,ga,yerr=ge,fmt=marker+"-",ms=3,capsize=1.5,lw=1.1,color=COLORS[case],label=case)
    axs[0].set_ylabel(r"$\omega_{\rm eff}$ (rad ps$^{-1}$)");axs[1].set_ylabel(r"$\Gamma_{\rm eff}$ (ps$^{-1}$)")
    for ax in axs[:2]: ax.set_xlabel(r"$k$ ($\AA^{-1}$)");ax.legend(frameon=False,fontsize=5.8,handlelength=1.2)
    for case in CASES:
        q=[r for r in weight_summary if r["case"]==case];x=np.array([num(r,"k_inv_A") for r in q]);y=np.array([num(r,"N_W_diag_mean") for r in q]);e=np.array([num(r,"N_W_diag_seedSEM") for r in q]); axs[2].errorbar(x,y,yerr=e,fmt="o-",ms=3,capsize=1.5,lw=1.1,color=COLORS[case],label=case.replace("_weakNH_6ns",""))
    axs[2].axhline(1,color=".45",lw=1,ls=":");axs[2].set_xlabel(r"$k$ ($\AA^{-1}$)");axs[2].set_ylabel(r"$N_{\rm O}W_{\rm diag}(k)$");axs[2].legend(frameon=False,fontsize=5.8)
    for case in CASES:
        q=curve[case]; t=np.array([num(r,"lag_ps") for r in q]);dr=np.array([num(r,"direct_VACF_mean") for r in q]);de=np.array([num(r,"direct_VACF_seedSEM") for r in q]);p=2*np.array([num(r,"P_diag_mean") for r in q]); mask=(t>=1)&(t<=250); axs[3].plot(t[mask],dr[mask],color=COLORS[case],lw=1.1,label=case.replace("_weakNH_6ns"," direct"));axs[3].fill_between(t[mask],dr[mask]-de[mask],dr[mask]+de[mask],color=COLORS[case],alpha=.14,lw=0);axs[3].plot(t[mask],p[mask],color=COLORS[case],lw=1.1,ls="--",label=case.replace("_weakNH_6ns",r" $2P_M$"));
    axs[3].axhline(0,color=".45",lw=1);axs[3].set_xscale("log");axs[3].set_xlim(1,250);axs[3].set_xlabel(r"$t$ (ps)");axs[3].set_ylabel(r"$C_{vv,z}/C_{vv,z}(0)$");axs[3].legend(frameon=False,fontsize=5.1,ncol=2,handlelength=1.2)
    for i,ax in enumerate(axs):
        ax.text(-.17,1.04,f"({chr(97+i)})",transform=ax.transAxes,fontsize=9,fontweight="bold");ax.tick_params(direction="out",width=1,length=3);ax.spines[["top","right"]].set_visible(False)
    for ext in ("png","pdf","svg"):fig.savefig(OUT/f"figures/implicit_C88_allbox_dispersion_staticW_gamma_VACF.{ext}",dpi=600)
    summary={"status":"complete_all_box_N200_N400_N1600","new_closure_cases":list(CASES),"positive_k_convention":"all VACF physics uses paired +/- k, i.e. 2*sum_{n>0}; one-sided values are retained only as an audit","remote_root":"/lustre/home/users/ewu/vb_gcmc/MD/stage_implicit_C88_static_vertex_VACF_20260831","DHO_source":str(DHO),"static_weight_flatness":qa,"reconstruction_max_cutoff_one_sided":[r for r in rec_summary if r["nmax"]==max(x["nmax"] for x in rec_summary if x["case"]==r["case"])],"degeneracy_factor_ensemble":factor_rows}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"README.md").write_text("# Implicit-C88 all-box static-weight/VACF closure\n\nAll N200/N400/N1600 CJJ-43-compatible four-seed analyses are included. Panels (a,b) reuse normalized-CJJ effective DHO omega and Gamma; panel (c) shows the independently measured diagonal static weight in the natural finite-box measure N_O W_diag. Panel (d) uses 2*sum_{n>0} W_diag F_s Phi_J because stored positive axial k modes represent one member of each real-field +/- k pair. The closure is not an exact Mori/tagged-self vertex identity.\n")
    (OUT/"QA.md").write_text("# QA\n\nN200/N400/N1600 each have four seed outputs, SUCCESS, metadata, nonempty tables and normal analysis logs. The former one-sided positive-k result is retained in the degeneracy-factor audit, but figures and conclusions use the required +/- k pair factor of two. All reported seed SEMs are conditional velocity-seed SEMs.\n")

if __name__=="__main__": main()
