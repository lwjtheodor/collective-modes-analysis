#!/usr/bin/env python3
"""Unpaired N800 low-k CJJ comparison: 100 fs/6 ns versus 1 ps/20 ns."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[1]
OLD=R/'remote_fetch/stage_C99_unified_longitudinal_current_VACF_500ps_20260829/output/N800'
NEW=R/'results/collective_mode_response/implicit_C99_multirate_current_modes_staticW/2026-08-31/source_data/static_W/N800'
OUT=R/'results/collective_mode_response/C99_N800_lowk_CJJ_100fs6ns_vs_1ps20ns/2026-08-31'
plt.rcParams.update({'font.family':'Arial','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.linewidth':1,'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'out','ytick.direction':'out','xtick.major.width':1,'ytick.major.width':1,'pdf.fonttype':42,'svg.fonttype':'none'})
COL={'old':'#1769aa','new':'#d95f02'}

def load(root, pattern):
 arr=[]; metas=[]
 for rep in range(1,5):
  z=np.load(root/f'rep{rep}'/pattern); arr.append(z['cjj_raw']); metas.append({'rep':rep,'lag_ps':z['lag_ps'],'n':z['n_values'],'k':z['kz_inv_A'], 'c0':z['cjj_raw'][0]})
 return np.asarray(arr),metas

def save(fig,name):
 for ext in ('.png','.pdf','.svg'): fig.savefig(OUT/'figures'/f'{name}{ext}',dpi=300 if ext=='.png' else None)
 plt.close(fig)

def main():
 (OUT/'figures').mkdir(parents=True,exist_ok=True); (OUT/'derived_data').mkdir(exist_ok=True)
 old,mo=load(OLD,'rep_arrays.npz'); new,mn=load(NEW,'rep_arrays_early_highk.npz')
 # Common n=1..5.  Normalize each seed with its own CJJ(k,0), then aggregate.
 nsel=[1,2,3,4,5]; oldn=old[:,:,:5]/old[:,:1,:5]; newn=new[:,:,:5]/new[:,:1,:5]
 to=mo[0]['lag_ps']; tn=mn[0]['lag_ps']; max_t=200.; io=to<=max_t; inn=tn<=max_t
 rows=[]
 for source,x,t in [('100fs_6ns_all_origin',oldn,to),('1ps_20ns_sparse_origin',newn,tn)]:
  for ni,n in enumerate(nsel):
   mean=x[:,:,ni].mean(0); sem=x[:,:,ni].std(0,ddof=1)/2
   for ti,time in enumerate(t): rows.append({'source':source,'n':n,'k_inv_A':mo[0]['k'][ni],'lag_ps':time,'CJJ_normalized_mean':mean[ti],'CJJ_normalized_seedSEM':sem[ti]})
 pd.DataFrame(rows).to_csv(OUT/'derived_data/N800_lowk_CJJ_normalized_mean_seedSEM.csv',index=False)
 # Mean plus cross-seed SEM for every low-k mode.
 fig=plt.figure(figsize=(7,7.6));
 for ri,n in enumerate(nsel):
  left=fig.add_axes([.10,.83-ri*.155,.34,.115]); right=fig.add_axes([.60,.83-ri*.155,.34,.115])
  for label,x,t,mask in [('100 fs / 6 ns',oldn,to,io),('1 ps / 20 ns',newn,tn,inn)]:
   c=COL['old' if label.startswith('100') else 'new']; ni=n-1; y=x[:,:,ni]; mean=y.mean(0); sem=y.std(0,ddof=1)/2
   left.plot(t[mask],mean[mask],lw=1.05,color=c,label=label); left.fill_between(t[mask],mean[mask]-sem[mask],mean[mask]+sem[mask],color=c,alpha=.18,lw=0)
   right.plot(t[mask],sem[mask],lw=1.05,color=c,label=label)
  left.axhline(0,color='.45',lw=.7); left.set_ylabel(rf'$C_{{JJ}}(n={n},t)/C_{{JJ}}(n={n},0)$'); right.set_ylabel('seed SEM')
  if ri==0: left.legend(fontsize=5.8,frameon=False); right.legend(fontsize=5.8,frameon=False)
  if ri==4: left.set_xlabel(r'$t$ (ps)'); right.set_xlabel(r'$t$ (ps)')
  left.text(-.15,1.03,f'({chr(97+2*ri)})',transform=left.transAxes,fontweight='bold',fontsize=9); right.text(-.15,1.03,f'({chr(98+2*ri)})',transform=right.transAxes,fontweight='bold',fontsize=9)
 fig.text(.10,.975,'N800 low-k longitudinal $C_{JJ}(k,t)$: mean curves and cross-seed SEM.  Sources are unpaired campaigns.',fontsize=7)
 fig.text(.10,.02,'Old: 100 fs dump, 6 ns, all time origins.  New: 1 ps dump, 20 ns, 196 origins at 100 ps spacing.  Bands/SEM are four shared-parent velocity seeds.',fontsize=6)
 save(fig,'N800_lowk_CJJ_mean_and_seedSEM_100fs6ns_vs_1ps20ns')
 # Replica trajectories: direct visual check of why cross-seed SEM differs.
 fig=plt.figure(figsize=(7,7.6));
 for ri,n in enumerate(nsel):
  for ci,(label,x,t,mask,key) in enumerate([('100 fs / 6 ns',oldn,to,io,'old'),('1 ps / 20 ns',newn,tn,inn,'new')]):
   a=fig.add_axes([.10+ci*.50,.83-ri*.155,.35,.115]); ni=n-1
   for rep in range(4): a.plot(t[mask],x[rep,mask,ni],color=COL[key],alpha=.25,lw=.7)
   a.plot(t[mask],x[:,mask,ni].mean(0),color=COL[key],lw=1.2)
   a.axhline(0,color='.45',lw=.7); a.set_ylabel(rf'$n={n}$')
   if ri==0: a.set_title(label,fontsize=7)
   if ri==4: a.set_xlabel(r'$t$ (ps)')
   a.text(-.15,1.03,f'({chr(97+2*ri+ci)})',transform=a.transAxes,fontweight='bold',fontsize=9)
 fig.text(.10,.975,'Individual velocity-seed traces (faint) and their mean (dark): N800 normalized low-k $C_{JJ}$.',fontsize=7)
 fig.text(.10,.02,'This shows cross-seed dispersion; it is not a bootstrap confidence interval over all time origins.',fontsize=6)
 save(fig,'N800_lowk_CJJ_replica_traces_100fs6ns_vs_1ps20ns')
 summary={'comparison_status':'unpaired campaigns; estimator and source trajectory both differ','old_dynamic_origin_policy':'all origins of 60001 frames at 0.1 ps cadence over 6 ns','new_dynamic_origin_policy':'196 sampled origins at 100 ps spacing from 20001 frames at 1 ps cadence over 20 ns','new_static_origin_policy':'100 static samples at nominal 196.97 ps spacing','uncertainty':'SEM across four velocity seeds sharing a parent configuration'}
 (OUT/'derived_data/estimator_support.json').write_text(json.dumps(summary,indent=2)+'\n')
 (OUT/'README.md').write_text('# N800 low-k CJJ: 100 fs/6 ns versus 1 ps/20 ns\n\nThis is an unpaired, protocol/estimator comparison, not a claim of a pure cadence effect. Each CJJ curve is normalized within seed before the mean and seed SEM are calculated.\n')
 (OUT/'QA.md').write_text('Old source: CJJ-44 unified N800 rep_arrays.npz. New source: CJJ-53 static-W N800 rep_arrays_early_highk.npz. Both use instantaneous oxygen axial COM subtraction and positive-n real-field convention; origins and campaigns differ as registered in estimator_support.json.\n')
 (OUT/'FINISHED.txt').write_text('N800 low-k CJJ cadence/estimator comparison rendered successfully.\n')
if __name__=='__main__': main()
