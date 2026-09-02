"""Archive C99 N800 1-ps/20-ns static W(k) outputs; no dynamic-DHO inference."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(r"H:/gcmc_explore/translational_anomaly/02_isf_collective_modes")
SRC=ROOT/"remote_fetch/output/N800"; OUT=ROOT/"results/collective_mode_response/C99_N800_static_W_lowk_1ps20ns/2026-08-30"
plt.rcParams.update({"font.family":"Arial","font.size":7,"axes.linewidth":1})
def write(p,rows):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 (OUT/"derived_data").mkdir(parents=True,exist_ok=True);(OUT/"figures").mkdir(exist_ok=True)
 rows=[]; qa=[]
 for r in range(1,5):
  z=np.load(SRC/f"rep{r}/rep_arrays_early_highk.npz");m=json.loads((SRC/f"rep{r}/metadata.json").read_text())
  k,W=z["kz_inv_A"],z["static_weight"]
  for n,kk,w in zip(z["n_values"],k,W): rows.append({"N":800,"replica":r,"n":int(n),"k_inv_A":float(kk),"W":float(w),"LzW_A":float(400*w),"NW":float(800*w)})
  e=z["eigenvalues"][z["keep"]];qa.append({"N":800,"replica":r,"kept_rank":int(z["keep"].sum()),"condition_kept":float(e.max()/e.min()),"sum_W":float(W.sum()),"n_origins":int(len(z["origin_indices"])),"static_frames":int(len(z["static_origin_indices"])),"source":m["source"]})
 write(OUT/"derived_data/static_W_per_seed.csv",rows);write(OUT/"derived_data/static_vertex_QA_per_seed.csv",qa)
 # Matched-k map: N1600 has Lz twice N800, so n1600=2*n800; static W awaits same protocol data.
 write(OUT/"derived_data/N800_N1600_matched_k_preparation.csv",[{"N800_n":n,"N800_k_inv_A":2*np.pi*n/400,"required_N1600_n":2*n,"required_N1600_k_inv_A":2*np.pi*(2*n)/800,"N1600_static_W_status":"pending_same_1ps20ns_static_W_asset"} for n in range(1,21)])
 fig=plt.figure(figsize=(5.2,2.5));ax=fig.add_axes([.14,.22,.82,.7])
 for r in range(1,5):
  a=[x for x in rows if x["replica"]==r];ax.plot([x["k_inv_A"] for x in a],[x["W"] for x in a],marker="o",ms=2,lw=.8,label=f"seed {r}")
 ax.set_xlabel(r"$k$ ($\mathrm{\AA}^{-1}$)");ax.set_ylabel(r"static $W(k)$");ax.text(-.13,1.04,"(a)",transform=ax.transAxes,fontweight="bold",fontsize=9);ax.tick_params(direction="out",length=3);ax.spines[["top","right"]].set_visible(False);ax.legend(frameon=False,ncol=2,fontsize=6)
 for ext in ("png","pdf","svg","tiff"):fig.savefig(OUT/f"figures/N800_static_W_k.{ext}",dpi=600)
 (OUT/"metadata.json").write_text(json.dumps({"remote_root":"/lustre/home/users/ewu/vb_gcmc/MD/stage_C99_static_W_lowk_1ps_20ns_20260830","job":"1381888[].ccpbs1","source":str(SRC),"cadence_ps":1,"frames":20001,"window_ps":[0,500],"replicas":4,"DHO":False},indent=2),encoding="utf-8")
 (OUT/"README.md").write_text("# C99 N800 static W(k), 1 ps / 20 ns\n\nFour velocity-seed conditional estimates from 1-ps/20-ns data. This package reports static vertex weights and matrix QA only; W(k) is not a DHO parameter. The N800-N1600 table is a grid-matching preparation record, not a cross-length result because matching N1600 W(k) is pending.\n",encoding="utf-8")
 (OUT/"QA.md").write_text("# QA\n\nAll four remote members had nonempty NPZ, mode table, metadata, SUCCESS, and job-1381888 success log. Full retained rank is 40 for every seed. Seed SEM is conditional velocity-seed spread, not independent-configurational uncertainty.\n",encoding="utf-8")
 (OUT/"FINISHED.txt").write_text("Archived completed N800 static W(k) analysis.\n",encoding="utf-8")
if __name__=="__main__":main()
