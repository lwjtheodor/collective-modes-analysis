#!/usr/bin/env python3
"""Build a protocol-separated synthesis of CJJ-43/44/46/52/53 assets.

Creates C99 all-box longitudinal DHO overview, static W(k), a cadence/coverage
registry, and one multi-page PDF per cadence and transverse branch.  Each PDF
page is one physical k and contains C_L(k,t), C_T(k,t), C_L(k,omega),
C_T(k,omega).  No curve or peak is pooled across sampling layers.
"""
from pathlib import Path
import json, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/collective_mode_response/C99_CJJ43_44_46_52_53_multirate_synthesis/2026-08-31'
S53=ROOT/'results/collective_mode_response/implicit_C99_multirate_current_modes_staticW/2026-08-31'
S46=ROOT/'results/collective_mode_response/implicit_C99_allbox_lowk_dispersion_damping/2026-08-30'
LAYERS=['N800_100fs_rep1to3','N800_10fs_rep1to3','N800_10fs_highk_rep1to3','N800_1fs_highk_rep1to3']
PAGE_LAYERS=['N800_100fs_rep1to3','N800_10fs_rep1to3']
COL={'LA':'#0F4D92','TA_r':'#B64342','TA_theta':'#42949E'}
plt.rcParams.update({'font.family':'Arial','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.linewidth':1,'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'out','ytick.direction':'out','xtick.major.width':1,'ytick.major.width':1,'svg.fonttype':'none','pdf.fonttype':42})

def save(fig, path):
    for ext in ('.png','.pdf','.svg'): fig.savefig(path.with_suffix(ext),dpi=300 if ext=='.png' else None,bbox_inches='tight')
    plt.close(fig)

def overview():
    d=pd.read_csv(S46/'derived_data/C99_allbox_allmodes_effective_DHO_points.csv'); w=pd.read_csv(S53/'derived_static_W.csv'); peaks=pd.read_csv(S53/'derived_operational_peak_database.csv')
    fig=plt.figure(figsize=(7.0,2.55)); ax1=fig.add_axes([.09,.25,.23,.65]); ax2=fig.add_axes([.40,.25,.23,.65]); ax3=fig.add_axes([.71,.25,.23,.65])
    for L,g in d.groupby('Lz_nm'):
        a=g[g.DHO_status.eq('accepted')]; ax1.errorbar(a.k_inv_A,a.omega_rad_ps_mean,yerr=a.omega_rad_ps_seedSEM,fmt='o',ms=2.2,capsize=1,label=f'{L:g} nm')
        ax2.errorbar(a.k_inv_A,a.gamma_inv_ps_mean,yerr=a.gamma_inv_ps_seedSEM,fmt='o',ms=2.2,capsize=1)
    # CJJ-53 peaks appear as cadence-explicit operational peaks, never as DHO gamma.
    peak_style={'LA':('x','black'),'TA_r':('^',COL['TA_r']),'TA_theta':('s',COL['TA_theta'])}
    for (layer,branch),g in peaks.groupby(['layer','branch']):
        marker,color=peak_style[branch]
        ax1.errorbar(g.k_inv_A,g.omega_peak_mean_rad_ps,yerr=g.omega_peak_replica_SEM_rad_ps,fmt=marker,ms=3,mew=.9,color=color,alpha=.7)
    for N,g in w.groupby('N'):
        ax3.errorbar(g.k,g.W,yerr=g['sem'],fmt='o',ms=2.8,capsize=1.2,label=f'N{N}')
    for ax,ylab in [(ax1,r'$\omega_L$ (rad ps$^{-1}$)'),(ax2,r'$\Gamma_{L,\mathrm{eff}}$ (ps$^{-1}$)'),(ax3,r'$W(k)$')]:
        ax.set_xscale('log'); ax.set_yscale('log' if ax is not ax3 else 'linear'); ax.set_xlabel(r'$k$ (Å$^{-1}$)'); ax.set_ylabel(ylab)
    ax1.legend(title='$L_z$',fontsize=5.5,title_fontsize=6,ncol=2,handletextpad=.2,columnspacing=.6); ax3.legend(fontsize=6)
    fig.text(.09,.95,'CJJ-46 accepted effective-DHO points; CJJ-53 N800 operational peaks: L black ×, T_r red △, T_θ teal □ (cadence separated).',fontsize=6)
    fig.text(.09,.05,r'CJJ-52 high-k points are deliberately absent: 100 fs/2 ps records have no identifiable $\omega,\Gamma,\phi$.',fontsize=6)
    save(fig,OUT/'figures/C99_omega_gamma_staticW_protocol_separated')

def pages(layer, transverse):
    base=S53/'source_data/dynamic'/layer/'derived_data'; c=pd.read_csv(base/'CJJ_all_modes_ensemble_mean_sem.csv'); s=pd.read_csv(base/'current_spectra_all_modes_ensemble_mean_sem.csv')
    pdf=OUT/f'figures/per_k_{layer}_{transverse}.pdf'; pdf.parent.mkdir(parents=True,exist_ok=True)
    registry=[]
    with PdfPages(pdf) as book:
      # One page per explicitly listed k.  Full n-by-n auto-correlation pages
      # remain in the CJJ-53 source archive; this extension prioritizes the
      # physically interpretable low-to-intermediate k spine per cadence.
      alln=sorted(c.n.unique()); selected=[n for n in (1,20) if n in set(alln)]
      for n in selected:
        cl=c[(c.branch=='LA')&(c.n==n)]; ct=c[(c.branch==transverse)&(c.n==n)]; sl=s[(s.branch=='LA')&(s.n==n)]; st=s[(s.branch==transverse)&(s.n==n)]
        k=float(cl.k_inv_A.iloc[0]); fig,axs=plt.subplots(2,2,figsize=(5.5,4.0)); fig.subplots_adjust(.13,.12,.98,.91,wspace=.35,hspace=.35)
        for ax,g,branch,ylabel in [(axs[0,0],cl,'LA',r'$C_L/C_L(0)$'),(axs[0,1],ct,transverse,r'$C_T/C_T(0)$')]:
          ax.plot(g.lag_ps,g.CJJ_normalized_mean,color=COL[branch],lw=1.05); ax.fill_between(g.lag_ps,g.CJJ_normalized_mean-g.CJJ_normalized_replica_SEM,g.CJJ_normalized_mean+g.CJJ_normalized_replica_SEM,color=COL[branch],alpha=.18,lw=0); ax.axhline(0,color='.45',lw=.7); ax.set_xlabel(r'$t$ (ps)'); ax.set_ylabel(ylabel)
        for ax,g,branch,ylabel in [(axs[1,0],sl,'LA',r'$S_L(k,\omega)$ (arb.)'),(axs[1,1],st,transverse,r'$S_T(k,\omega)$ (arb.)')]:
          ax.plot(g.frequency_ps_inv*2*np.pi,g.PSD_mean_arbitrary,color=COL[branch],lw=1.05); ax.fill_between(g.frequency_ps_inv*2*np.pi,np.maximum(0,g.PSD_mean_arbitrary-g.PSD_replica_SEM_arbitrary),g.PSD_mean_arbitrary+g.PSD_replica_SEM_arbitrary,color=COL[branch],alpha=.18,lw=0); ax.set_xlabel(r'$\omega$ (rad ps$^{-1}$)'); ax.set_ylabel(ylabel)
        for i,ax in enumerate(axs.flat): ax.text(-.18,1.04,f'({chr(97+i)})',transform=ax.transAxes,fontweight='bold',fontsize=9)
        fig.suptitle(f'{layer}; n={n},  k={k:.5f} Å$^{{-1}}$; T={transverse}',fontsize=7)
        book.savefig(fig); plt.close(fig); registry.append({'layer':layer,'transverse_branch':transverse,'n':int(n),'k_inv_A':k,'pdf':pdf.name,'page':int(n)})
    return registry

def main():
    (OUT/'figures').mkdir(parents=True,exist_ok=True); (OUT/'derived_data').mkdir(exist_ok=True); overview(); rows=[]
    for layer in LAYERS:
      meta=json.loads((S53/'source_data/dynamic'/layer/'metadata.json').read_text()); im=meta['input_replicas'][0]
      rows.append({'layer':layer,'dt_ps':im['dt_ps'],'duration_ps':im['duration_ps'],'nmax':meta['mode_range']['n_max'],'frequency_resolution_ps_inv':meta['frequency_resolution_ps_inv'],'omega_gamma_status':'operational spectral peaks only; not DHO gamma'})
      if layer in PAGE_LAYERS:
        for t in ('TA_r','TA_theta'): rows+=pages(layer,t)
    pd.DataFrame(rows).to_csv(OUT/'derived_data/CJJ53_cadence_and_per_k_page_registry.csv',index=False)
    (OUT/'README.md').write_text('# CJJ-43/44/46/52/53 protocol-separated synthesis\n\nCJJ-46 supplies C99 all-box low-k longitudinal effective-DHO points. CJJ-53 supplies N800 branch-resolved autocorrelations and operational spectra in mutually separated cadence layers. CJJ-52 high-k is excluded from omega/Gamma because parameters are not identifiable. The per-k PDFs render n=1 and n=20 for the two cadence layers that support this 4-panel layout; the native CJJ-53 all-n pages remain the complete per-k archive. T_r and T_theta are separate.\n')
    (OUT/'QA.md').write_text('All figure inputs are compact archived CSVs. No curve/PSD was averaged across cadence. All uncertainty bands are velocity-seed SEM conditional on shared parent configuration.\n')
    (OUT/'FINISHED.txt').write_text('Protocol-separated CJJ-43/44/46/52/53 synthesis rendered successfully.\n')
if __name__=='__main__': main()
