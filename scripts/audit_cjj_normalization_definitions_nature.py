"""Separate three non-equivalent normalizations of the n=1 CJJ lobe."""
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams.update({'font.family':'Arial','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.8,'legend.frameon':False})
ROOT=Path(__file__).resolve().parents[1]; IN=ROOT/'results'/'collective_mode_response'/'fig2_longitudinal_modes_88_rh75_330k'/'2026-08-11'/'derived_data'/'panel_b_lowk_strength.csv'; OUT=ROOT/'assets'
def main():
 d=pd.read_csv(IN);fig,ax=plt.subplots(1,3,figsize=(7.25,2.55),constrained_layout=True); cols=['#C85250','#365C8D','#319A9A']
 specs=[('depth_norm_mean','depth_norm_sem',r'depth: $-\min[C_{JJ}/C_{JJ}(0)]$','dimensionless'),('A_minus_ps_mean','A_minus_ps_sem',r'area: $-\int C_{JJ}/C_{JJ}(0)\,dt$','ps'),('A_minus_norm_mean','A_minus_norm_sem',r'time-scaled area: $A_-/(L_z/c_s)$','dimensionless')]
 for i,(m,e,title,unit) in enumerate(specs):
  ax[i].errorbar(d.Lz_nm,d[m],yerr=d[e],color=cols[i],marker='o',lw=1.3,capsize=2)
  ax[i].set(title=title,xlabel=r'$L_z$ (nm)');ax[i].set_xticks([20,40,60,80,100]);ax[i].text(.05,.92,f'unit: {unit}',transform=ax[i].transAxes,va='top',fontsize=6.7)
 ax[0].set_ylabel('value (distinct normalization per panel)')
 fig.suptitle(r'$n=1$ current-mode lobe: normalization audit',fontsize=11.5,fontweight='bold',y=1.03)
 fig.text(.5,-.06,r'All use the same zero-crossing-bounded first lobe and the same three replicas per length. $C_{JJ}/C_{JJ}(0)$ removes equal-time amplitude; only the third panel additionally divides by the acoustic return time.',ha='center',fontsize=6.6)
 stem=OUT/'cjj_n1_normalization_definition_audit_nature';fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.tiff'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
if __name__=='__main__':main()
