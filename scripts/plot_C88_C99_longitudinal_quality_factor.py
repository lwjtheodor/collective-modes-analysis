#!/usr/bin/env python3
"""Protocol-matched C88/C99 longitudinal effective-DHO comparison."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

R=Path(__file__).resolve().parents[1]
SRC=R/'results/collective_mode_response/implicit_C77_C88_C99_matched_k_effective_DHO_damping/2026-08-29/derived_data/matched_k_effective_DHO_summary.csv'
OUT=R/'results/collective_mode_response/C88_C99_longitudinal_effective_DHO_Q_comparison/2026-08-31'
plt.rcParams.update({'font.family':'Arial','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.linewidth':1,'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'out','ytick.direction':'out','xtick.major.width':1,'ytick.major.width':1,'pdf.fonttype':42,'svg.fonttype':'none'})
def main():
 (OUT/'figures').mkdir(parents=True,exist_ok=True); (OUT/'derived_data').mkdir(exist_ok=True)
 d=pd.read_csv(SRC); d=d[d.system.isin(['C88','C99'])].copy(); d['Q_effective']=d.omega_mean_rad_ps/(2*d.gamma_effective_psinv_mean); d.to_csv(OUT/'derived_data/C88_C99_effective_DHO_points_with_Q.csv',index=False)
 fig=plt.figure(figsize=(7,2.45)); boxes=[(.09,.23,.25,.66),(.41,.23,.25,.66),(.73,.23,.25,.66)]; axes=[fig.add_axes(x) for x in boxes]
 for system,color,mark in [('C88','#b34a3c','o'),('C99','#1769aa','s')]:
  x=d[d.system.eq(system)]
  for a,col,ylab in zip(axes,['Q_effective','omega_mean_rad_ps','gamma_effective_psinv_mean'],[r'$Q_{\mathrm{eff}}=\omega/(2\Gamma)$',r'$\omega_{L,\mathrm{eff}}$ (rad ps$^{-1}$)',r'$\Gamma_{L,\mathrm{eff}}$ (ps$^{-1}$)']):
   a.errorbar(x.k_Ainv,x[col],fmt=mark,color=color,ms=3,mew=.8,capsize=1,label=system if a is axes[0] else None)
   a.set(xscale='log',yscale='log',xlabel=r'$k$ (Å$^{-1}$)',ylabel=ylab)
 axes[0].axhline(1,color='.4',lw=.8); axes[0].legend(fontsize=6,frameon=False)
 for i,a in enumerate(axes): a.text(0,1.06,f'({chr(97+i)})',transform=a.transAxes,fontweight='bold',fontsize=9)
 fig.text(.09,.96,'350 K weak-NH, 6 ns, 100 fs, four velocity-seed normalized-$C_{JJ}$ effective damped-cosine fits.',fontsize=6.5)
 fig.text(.09,.04,'C88 N800 direct-NVT pilot excluded.  Points are finite-k effective parameters, not spectral HWHM or $k\to0$ limits.',fontsize=6)
 for ext in ('.png','.pdf','.svg'): fig.savefig(OUT/'figures'/f'C88_C99_longitudinal_effective_DHO_Q_comparison{ext}',dpi=300 if ext=='.png' else None)
 plt.close(fig)
 (OUT/'README.md').write_text('# C88 versus C99 longitudinal effective-DHO quality factor\n\nOnly CJJ-43 protocol-matched C88/C99 weak-NH production fits are used. C88 N800 direct-NVT pilot is excluded. `Q_eff` is a finite-k damped-cosine descriptor, not a thermodynamic-limit mode classification.\n')
 (OUT/'QA.md').write_text('Both systems are 350 K, weak-NH, 6 ns, 100 fs, and four velocity seeds. Physical k values are plotted directly; not every C88/C99 point is an exact matched-k pair.\n')
 (OUT/'FINISHED.txt').write_text('C88/C99 longitudinal effective-DHO Q comparison rendered successfully.\n')
if __name__=='__main__': main()
