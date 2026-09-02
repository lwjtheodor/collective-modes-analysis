"""Test the requested de-COM MSD/L vs lag/L collapse using all replicas."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({'font.family':'Arial','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.8,'legend.frameon':False})
ROOT=Path(__file__).resolve().parents[1]; IN=ROOT/'heartbeat_fetch'/'stage_alpha200ps_crosschirality_20260809'/'output'; OUT=ROOT/'assets'
CHIS=['7_7','8_8','9_9','17_0']; LAB={'7_7':'(7,7)','8_8':'(8,8)','9_9':'(9,9)','17_0':'(17,0)'}; COL={1:'#355F8C',2:'#4F8F9A',3:'#75A65A',4:'#D5A044',5:'#C86B3C',10:'#A34358'}

def read(f): return np.genfromtxt(f,delimiter=',',names=True)
def save(fig,stem):
 fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight')
 fig.savefig(stem.with_suffix('.tiff'),dpi=600,bbox_inches='tight')
 fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight')
 fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
def main():
 rows=[]; scores=[]; fig,axs=plt.subplots(2,2,figsize=(7.25,5.6),constrained_layout=True)
 for p,(ax,chi) in enumerate(zip(axs.flat,CHIS)):
  means={}
  for L in [1,2,3,4,5,10]:
   fs=sorted(IN.glob(f'{chi}_L{L}_rep*_msd.csv')); ds=[read(f) for f in fs]; t=ds[0]['lag_ps']/L; y=np.asarray([d['msd_z_A2']/L for d in ds]); m=y.mean(0); se=y.std(0,ddof=1)/np.sqrt(len(y)) if len(y)>1 else np.zeros_like(m); means[L]=(t,m)
   for x,z,e in zip(t,m,se): rows.append({'chirality':chi,'L':L,'scaled_lag_ps_per_L':x,'msd_over_L_A2_per_L':z,'msd_over_L_sem':e,'n_replicates':len(ds)})
   ax.plot(t,m,color=COL[L],lw=1.3,label=f'{L}L');ax.fill_between(t,m-se,m+se,color=COL[L],alpha=.16,lw=0)
  # Quantify disagreement on the common scaled-time support [0.5,20] ps/L.
  grid=np.geomspace(.5,20,80); Y=np.array([np.interp(grid,means[L][0],means[L][1]) for L in means]); cv=np.std(Y,axis=0,ddof=1)/np.mean(Y,axis=0)
  scores.append({'chirality':chi,'scaled_time_min_ps_per_L':.5,'scaled_time_max_ps_per_L':20.,'mean_CV_across_L':float(np.mean(cv)),'median_CV_across_L':float(np.median(cv)),'criterion':'SD across six L / mean across six L; interpolation of each length mean onto common grid'})
  ax.set_xscale('log');ax.set_xlim(.05,200);ax.set_title(f'{chr(97+p)}  {LAB[chi]}',loc='left',fontweight='bold',fontsize=10);ax.text(.04,.90,f'common-domain CV = {np.mean(cv):.2f}',transform=ax.transAxes,fontsize=6.7,va='top')
  if p%2==0: ax.set_ylabel(r'$M_z^{\rm de-COM}/L$ ($\mathrm{\AA^2}/L$)')
  if p>=2: ax.set_xlabel(r'scaled lag, $\Delta t/L$ (ps/L)')
  if p==0: ax.legend(ncol=2,fontsize=6.4,loc='lower right',handlelength=1.4)
 fig.suptitle(r'De-COM MSD scaling test: $M_z/L$ versus $\Delta t/L$',fontsize=12,fontweight='bold',y=1.015)
 fig.text(.5,-.025,'Instantaneous water-COM coordinate removed before MSD; weak Nosé–Hoover; no momentum removal; 1–5L: 10 fs/1 ns, 10L: 100 fs/10 ns; curves are mean ± replica SEM.',ha='center',fontsize=7)
 save(fig,OUT/'crosschirality_1L10L_decom_msd_L_collapse_nature');plt.close(fig)
 for name,data in [('crosschirality_1L10L_decom_msd_L_collapse_curves.csv',rows),('crosschirality_1L10L_decom_msd_L_collapse_scores.csv',scores)]:
  with (OUT/name).open('w',newline='',encoding='utf-8') as h: w=csv.DictWriter(h,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
if __name__=='__main__': main()
