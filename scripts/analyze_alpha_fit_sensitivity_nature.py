"""Fit-window/anchoring sensitivity of de-COM MSD alpha minima."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
mpl.rcParams.update({'font.family':'Arial','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.8,'legend.frameon':False})
ROOT=Path(__file__).resolve().parents[1]; IN=ROOT/'heartbeat_fetch'/'stage_alpha200ps_crosschirality_20260809'/'output'; OUT=ROOT/'assets'
CHIS=['7_7','8_8','9_9','17_0'];LAB={'7_7':'(7,7)','8_8':'(8,8)','9_9':'(9,9)','17_0':'(17,0)'}
METHODS=[('trailing',.5,'T, 0.5 dec','#4C78A8'),('trailing',1.,'T, 1 dec','#E08E2E'),('trailing',1.5,'T, 1.5 dec','#59A14F'),('centered',1.,'C, 1 dec','#B24A5A')]
def read(f): return np.genfromtxt(f,delimiter=',',names=True)
def calc(t,m,anchor,w):
 a=np.full(len(t),np.nan); r=10**(w/2)
 for i,tt in enumerate(t):
  keep=(t>=tt/10**w)&(t<=tt) if anchor=='trailing' else (t>=tt/r)&(t<=tt*r)
  if keep.sum()>=5 and (tt/10**w>=.5 if anchor=='trailing' else tt/r>=.5) and (tt<=200 if anchor=='trailing' else tt*r<=200): a[i]=np.polyfit(np.log(t[keep]),np.log(m[keep]),1)[0]
 ii=np.flatnonzero(np.isfinite(a)); ii=ii[np.argmin(a[ii])]
 return float(a[ii]),float(t[ii])
def save(fig,stem):
 fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.tiff'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
def main():
 cases=[]; summary=[]
 for chi in CHIS:
  for L in [1,2,3,4,5,10]:
   for f in sorted(IN.glob(f'{chi}_L{L}_rep*_msd.csv')):
    d=read(f)
    for anchor,w,label,color in METHODS:
     amin,tmin=calc(d['lag_ps'],d['msd_z_A2'],anchor,w); cases.append({'chirality':chi,'L':L,'case_id':d['case_id'][0],'anchor':anchor,'window_decades':w,'method_label':label,'alpha_min':amin,'t_min_ps':tmin,'t_min_over_L_ps_per_L':tmin/L})
 for key in {(x['chirality'],x['L'],x['anchor'],x['window_decades'],x['method_label']) for x in cases}:
  q=[x for x in cases if (x['chirality'],x['L'],x['anchor'],x['window_decades'],x['method_label'])==key];
  for col in ['alpha_min','t_min_ps','t_min_over_L_ps_per_L']:
   v=np.array([x[col] for x in q]); summary.append({'chirality':key[0],'L':key[1],'anchor':key[2],'window_decades':key[3],'method_label':key[4],'observable':col,'mean':v.mean(),'sem':v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else 0.,'n_replicates':len(v)})
 fig,axs=plt.subplots(2,4,figsize=(7.25,4.3),sharex='col',constrained_layout=True)
 for j,chi in enumerate(CHIS):
  for anchor,w,label,col in METHODS:
   for row,obs in enumerate(['alpha_min','t_min_over_L_ps_per_L']):
    q=[x for x in summary if x['chirality']==chi and x['anchor']==anchor and x['window_decades']==w and x['observable']==obs];q=sorted(q,key=lambda x:x['L']);x=np.array([z['L'] for z in q]);y=np.array([z['mean'] for z in q]);e=np.array([z['sem'] for z in q]);axs[row,j].errorbar(x,y,yerr=e,color=col,marker='o',ms=3,lw=1.1,capsize=1.8,label=label)
  axs[0,j].set_title(f'{chr(97+j)}  {LAB[chi]}',loc='left',fontweight='bold',fontsize=10);axs[0,j].set_ylim(.15,1.08);axs[1,j].set_yscale('log');axs[1,j].set_xlabel('box length (L)');axs[1,j].set_xticks([1,2,3,4,5,10])
  if j==0: axs[0,j].set_ylabel(r'$\alpha_{z,\min}$');axs[1,j].set_ylabel(r'$t_{\min}/L$ (ps/L)');axs[0,j].legend(fontsize=5.8,loc='lower left',handlelength=1.2)
 fig.suptitle(r'Fit-definition sensitivity of the de-COM MSD exponent minimum',fontsize=11.5,fontweight='bold',y=1.025)
 fig.text(.5,-.055,'T: trailing window [t/10^w,t]; C: centred window [t/10^(w/2),t×10^(w/2)]; complete windows only; direct MSD with water-COM coordinate removed; points mean ± replica SEM.',ha='center',fontsize=6.55)
 save(fig,OUT/'crosschirality_1L10L_alpha_fit_sensitivity_nature');plt.close(fig)
 for name,data in [('crosschirality_1L10L_alpha_fit_sensitivity_cases.csv',cases),('crosschirality_1L10L_alpha_fit_sensitivity_summary.csv',summary)]:
  with (OUT/name).open('w',newline='',encoding='utf-8') as h:wrt=csv.DictWriter(h,fieldnames=list(data[0]));wrt.writeheader();wrt.writerows(data)
if __name__=='__main__':main()
