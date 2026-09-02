#!/usr/bin/env python3
"""Decision-oriented audit of CJJ-53 cadence and static-W improvements."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[1]
OLD=R/'results/collective_mode_response/implicit_C99_unified_longitudinal_500ps/2026-08-29/derived_data/C99_static_weight_k_dependence.csv'
NEW=R/'results/collective_mode_response/implicit_C99_multirate_current_modes_staticW/2026-08-31'
DHO=R/'results/collective_mode_response/implicit_C99_allbox_lowk_dispersion_damping/2026-08-30/derived_data/C99_allbox_allmodes_effective_DHO_points.csv'
OUT=R/'results/collective_mode_response/C99_CJJ53_sampling_frequency_statistical_audit/2026-08-31'
plt.rcParams.update({'font.family':'Arial','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.linewidth':1,'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'out','ytick.direction':'out','xtick.major.width':1,'ytick.major.width':1,'pdf.fonttype':42,'svg.fonttype':'none'})

def save(fig,name):
 for ext in ('.png','.pdf','.svg'): fig.savefig(OUT/'figures'/f'{name}{ext}',dpi=300 if ext=='.png' else None)
 plt.close(fig)

def main():
 (OUT/'figures').mkdir(parents=True,exist_ok=True); (OUT/'derived_data').mkdir(exist_ok=True)
 old=pd.read_csv(OLD); new=pd.read_csv(NEW/'derived_static_W.csv')
 x=old[(old.N_water==800)&(old.n<=20)][['n','k_Ainv','W_mean','W_seed_sem']].merge(new[new.N==800][['n','W','sem']],on='n')
 x['rel_sem_100fs6ns']=x.W_seed_sem/x.W_mean; x['rel_sem_1ps20ns']=x['sem']/x.W; x['delta_W']=x.W-x.W_mean; x['combined_z']=x.delta_W/np.hypot(x['sem'],x.W_seed_sem)
 x.to_csv(OUT/'derived_data/N800_static_W_100fs6ns_vs_1ps20ns.csv',index=False)
 peaks=pd.read_csv(NEW/'derived_operational_peak_database.csv')
 high=peaks[peaks.layer.isin(['N800_10fs_highk_rep1to3','N800_1fs_highk_rep1to3'])].copy()
 quality=high.groupby(['layer','branch']).agg(n_modes=('n','size'),n_resolved=('resolved_operational_peak','sum'),median_CV=('replica_frequency_CV','median'),n_CV_le_025=('replica_frequency_CV',lambda s:int((s<=.25).sum())),median_prominence=('mean_prominence_over_median','median')).reset_index()
 quality['resolved_fraction']=quality.n_resolved/quality.n_modes; quality.to_csv(OUT/'derived_data/highk_operational_peak_identifiability_10fs_vs_1fs.csv',index=False)
 dho=pd.read_csv(DHO); low=dho[(dho.DHO_status=='accepted')&(dho.k_inv_A<=.3141593)].copy(); low.to_csv(OUT/'derived_data/CJJ46_lowk_longitudinal_gamma_authority.csv',index=False)
 fig=plt.figure(figsize=(7.0,5.15)); ax=[fig.add_axes(b) for b in [(0.10,.59,.34,.32),(.60,.59,.34,.32),(.10,.12,.34,.32),(.60,.12,.34,.32)]]
 # (a) values: only N800 static estimators; not a dynamic/cadence comparison.
 ax[0].errorbar(x.k_Ainv,x.W_mean,yerr=x.W_seed_sem,fmt='o-',ms=3,capsize=1.5,lw=1.05,label='100 fs / 6 ns')
 ax[0].errorbar(x.k_Ainv,x.W,yerr=x['sem'],fmt='s-',ms=3,capsize=1.5,lw=1.05,label='1 ps / 20 ns')
 ax[0].set(xlabel=r'$k$ (Å$^{-1}$)',ylabel=r'$W(k)$'); ax[0].legend(fontsize=6,frameon=False)
 ax[0].text(.03,.05,r'$n=1ldots5$: $|z|<1.5$',transform=ax[0].transAxes,fontsize=6)
 # (b) uncertainty improvement test.
 ax[1].plot(x.n,x.rel_sem_1ps20ns/x.rel_sem_100fs6ns,'o-',color='#B64342',ms=3,lw=1.05); ax[1].axhline(1,color='.4',lw=.8); ax[1].set(xlabel=r'mode index $n$',ylabel=r'relative SEM ratio\n(1 ps / 20 ns) / (100 fs / 6 ns)'); ax[1].set_yscale('log')
 ax[1].text(.03,.05,f'median = {(x.rel_sem_1ps20ns/x.rel_sem_100fs6ns).median():.2f}',transform=ax[1].transAxes,fontsize=6)
 # (c) high-k operational peak status.  This is resolution, not DHO gamma.
 labels=['L','$T_r$','$T_\\theta$']; order=['LA','TA_r','TA_theta']; y=np.arange(3); width=.34
 for i,layer in enumerate(['N800_10fs_highk_rep1to3','N800_1fs_highk_rep1to3']):
  q=quality.set_index(['layer','branch']).loc[[ (layer,b) for b in order]]; ax[2].bar(y+(i-.5)*width,q.resolved_fraction,width,label='10 fs / 50 ps' if i==0 else '1 fs / 2 ps')
 ax[2].set(xticks=y,xticklabels=labels,ylim=(0,1.05),ylabel='resolved operational-peak fraction'); ax[2].legend(fontsize=6,frameon=False)
 ax[2].text(.03,.05,'gate: CV ≤ 0.25 and prominence criterion',transform=ax[2].transAxes,fontsize=5.5)
 # (d) gamma authority, segregated from high-k peak output.
 for L,g in low.groupby('Lz_nm'):
  ax[3].errorbar(g.k_inv_A,g.gamma_inv_ps_mean,yerr=g.gamma_inv_ps_seedSEM,fmt='o',ms=2.4,capsize=1,label=f'{L:g} nm')
 ax[3].set(xscale='log',yscale='log',xlabel=r'$k$ (Å$^{-1}$)',ylabel=r'$\Gamma_{L,\mathrm{eff}}$ (ps$^{-1}$)'); ax[3].legend(title='$L_z$',fontsize=5.5,title_fontsize=6,ncol=2,frameon=False)
 ax[3].text(.03,.05,r'No $Gamma_T$ or high-$k$ $Gamma$ is identifiable from CJJ-53.',transform=ax[3].transAxes,fontsize=5.3)
 for i,a in enumerate(ax): a.text(0,1.06,f'({chr(97+i)})',transform=a.transAxes,fontweight='bold',fontsize=9)
 fig.text(.10,.965,'CJJ-53 sampling-frequency audit: estimators, peak resolution, and linewidth authority are deliberately separated.',fontsize=7)
 save(fig,'CJJ53_sampling_frequency_statistical_audit')
 summary={'static_W_N800_n1to20':{'median_W_ratio_1ps20ns_over_100fs6ns':float((x.W/x.W_mean).median()),'median_relative_SEM_ratio':float((x.rel_sem_1ps20ns/x.rel_sem_100fs6ns).median()),'low_n1to5_max_abs_combined_z':float(x[x.n<=5].combined_z.abs().max())},'highk_operational_peaks':quality.to_dict(orient='records'),'gamma_boundary':'CJJ-46 accepted low-k longitudinal effective-DHO only; CJJ-53 operational peaks and CJJ-52 high-k records do not supply identifiable Gamma'}
 (OUT/'derived_data/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 (OUT/'README.md').write_text('# CJJ-53 sampling-frequency statistical audit\n\nFour distinct panels: static W values, relative-SEMs, high-k operational-peak resolvability, and authoritative low-k longitudinal Gamma. It intentionally does not display transverse operational peak positions as Gamma.\n')
 (OUT/'QA.md').write_text('All comparisons use N800, n=1..20 common physical modes. Old W: CJJ-44 100-fs/6-ns static vertex. New W: CJJ-53 1-ps/20-ns static vertex. Both use conditional velocity-seed SEM; neither is independent-configuration uncertainty.\n')
 (OUT/'FINISHED.txt').write_text('CJJ-53 sampling-frequency statistical audit rendered successfully.\n')
if __name__=='__main__': main()
