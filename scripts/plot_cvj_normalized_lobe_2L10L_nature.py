"""2L--10L normalized n=1 current-mode first-negative-lobe scaling."""
import csv
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({'font.family':'Arial','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.8,'legend.frameon':False})
ROOT=Path(__file__).resolve().parents[1]; IN=ROOT/'assets'/'crosschirality_1L10L_absolute_cvj_lobe_normalizations.csv'; OUT=ROOT/'assets'
LAB={'7_7':'(7,7)','8_8':'(8,8)','9_9':'(9,9)','17_0':'(17,0)'};COL={'7_7':'#2474A6','8_8':'#CE4A59','9_9':'#5A9B3F','17_0':'#7651A7'}

def fit(x,y):
 X=np.log(x);Y=np.log(y);p,b=np.polyfit(X,Y,1);r=Y-(p*X+b);se=np.sqrt((r@r)/(len(X)-2)/((X-X.mean())@(X-X.mean())));r2=1-(r@r)/((Y-Y.mean())@(Y-Y.mean()));return p,b,se,r2
def main():
 d=[r for r in csv.DictReader(IN.open()) if r['metric']=='first_negative_lobe_normalized' and int(r['L'])>=2]
 rec=[];fig,ax=plt.subplots(figsize=(3.9,3.55),constrained_layout=True)
 for chi in LAB:
  q=sorted([r for r in d if r['chirality']==chi],key=lambda r:int(r['L']));x=np.array([float(r['L']) for r in q]);y=np.array([float(r['area_mean']) for r in q]);e=np.array([float(r['area_sem']) for r in q]);p,b,se,r2=fit(x,y)
  ax.errorbar(x,y,yerr=e,color=COL[chi],marker='o',ms=4,lw=1.4,capsize=2,label=rf'{LAB[chi]}  $p={p:.2f}\pm{se:.2f}$')
  xx=np.geomspace(2,10,100);ax.plot(xx,np.exp(b)*xx**p,color=COL[chi],lw=.8,ls='--',alpha=.8)
  for z in q:rec.append({'chirality':chi,'L':z['L'],'A_negative_normalized_ps':z['area_mean'],'A_negative_normalized_sem_ps':z['area_sem'],'n_replicates':z['n_replicates'],'fit_range':'2L-10L','power_p':p,'power_p_sem':se,'fit_R2':r2})
 ax.set(xscale='log',yscale='log',xlim=(1.8,11.3),ylim=(1.25,30),xlabel='box length, $L$',ylabel=r'normalized first negative-lobe area, $A_-^{\rm norm}$ (ps)')
 ax.set_xticks([2,3,4,5,10]);ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter());ax.xaxis.set_minor_locator(mpl.ticker.NullLocator());ax.legend(fontsize=6.4,loc='upper left',handlelength=1.4)
 ax.set_title(r'$n=1$ current-mode negative lobe: 2L--10L',loc='left',fontweight='bold',fontsize=10)
 fig.text(.5,.006,r'$C_{vJ}(t)/C_{vJ}(0)$; zero-crossing bounded; points mean $\pm$ replica SEM; dashed log--log fit. Weak Nosé–Hoover; no momentum removal; water-COM velocity removed; 2–5L: 10 fs/1 ns, 10L: 100 fs/10 ns.',ha='center',fontsize=5.85)
 stem=OUT/'crosschirality_2L10L_normalized_cvj_first_negative_lobe_nature';fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.tiff'),dpi=600,bbox_inches='tight');fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight');plt.close(fig)
 with (OUT/'crosschirality_2L10L_normalized_cvj_first_negative_lobe.csv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=list(rec[0]));w.writeheader();w.writerows(rec)
if __name__=='__main__':main()
