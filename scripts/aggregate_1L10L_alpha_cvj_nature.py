"""Aggregate the file-backed 1-10L direct-MSD alpha and n=1 CvJ data."""
from __future__ import annotations
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; SRC=ROOT/'remote_fetch'/'output'; AS=ROOT/'assets'
REC=AS/'crosschirality_1L10L_alpha_cvj_case_records.csv'; SUM=AS/'crosschirality_1L10L_alpha_cvj_length_summary.csv'; FIT=AS/'crosschirality_1L10L_alpha_cvj_powerlaw_tests.csv'; OUT=AS/'crosschirality_1L10L_alpha_cvj_powerlaw_nature'
COL={'7_7':'#3B7EA1','8_8':'#D88737','9_9':'#5E9C76','17_0':'#8C6BB1'}; LAB={'7_7':'(7,7)','8_8':'(8,8)','9_9':'(9,9)','17_0':'(17,0)'}
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'axes.linewidth':.8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
def write(p,rows):
 with p.open('w',newline='',encoding='utf8') as f: w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sem(a): return float(np.std(a,ddof=1)/math.sqrt(len(a))) if len(a)>1 else float('nan')
def fit(rows,key,subset):
 s=[r for r in rows if r['L'] in subset]; x=np.log(np.array([r['L'] for r in s]));y=np.log(np.array([r[key+'_mean'] for r in s])); b,a=np.polyfit(x,y,1);res=y-(a+b*x);se=math.sqrt(float((res@res)/(len(x)-2)/((x-x.mean())@(x-x.mean()))));return a,b,se
def main():
 AS.mkdir(exist_ok=True); rec=[]
 for p in sorted(SRC.glob('*.json')):
  m=json.loads(p.read_text()); c=m['case_id']; chi,L,rep=c.rsplit('_',2); L=int(L[1:]); rep=int(rep[3:]); lob=m['first_negative_lobe']
  rec.append({'chirality':chi,'L':L,'replicate':rep,'dt_ps':m['dt_ps'],'duration_ps':m['n_frames']*m['dt_ps'],'alpha_min':m['alpha_min'],'t_alpha_min_ps':m['t_alpha_min_ps'],'memory_strength_1_minus_alpha':1-m['alpha_min'],'negative_area_ps':lob['negative_area_ps'],'depth':lob['depth'],'t_lobe_min_ps':lob['t_min_ps'],'k_inv_A':m['k_inv_A']})
 assert len(rec)==56, f'expected 56 cases, got {len(rec)}';write(REC,rec)
 summary=[]
 for chi in COL:
  for L in sorted({r['L'] for r in rec if r['chirality']==chi}):
   rr=[r for r in rec if r['chirality']==chi and r['L']==L];o={'chirality':chi,'L':L,'n_replicates':len(rr),'k_inv_A':np.mean([r['k_inv_A'] for r in rr])}
   for k in ('alpha_min','memory_strength_1_minus_alpha','negative_area_ps','depth','t_alpha_min_ps','t_lobe_min_ps'):
    a=np.array([r[k] for r in rr]);o[k+'_mean']=a.mean();o[k+'_sem']=sem(a)
   summary.append(o)
 write(SUM,summary)
 fits=[]
 for chi in COL:
  r=[x for x in summary if x['chirality']==chi]
  for key in ('negative_area_ps','memory_strength_1_minus_alpha'):
   a5,b5,e5=fit(r,key,[1,2,3,4,5]); aall,ball,eall=fit(r,key,[1,2,3,4,5,10]); actual=next(x[key+'_mean'] for x in r if x['L']==10); pred=math.exp(a5)*10**b5
   fits.append({'chirality':chi,'observable':key,'p_1to5':b5,'p_1to5_sem':e5,'r2_1to5':1-float(np.sum((np.log([x[key+'_mean'] for x in r if x['L']<=5])-(a5+b5*np.log([x['L'] for x in r if x['L']<=5])))**2))/float(np.sum((np.log([x[key+'_mean'] for x in r if x['L']<=5])-np.mean(np.log([x[key+'_mean'] for x in r if x['L']<=5])))**2)),'p_1to10':ball,'p_1to10_sem':eall,'predicted_10L':pred,'observed_10L':actual,'relative_deviation_10L_pct':100*(actual/pred-1)})
 write(FIT,fits)
 fig,ax=plt.subplots(2,2,figsize=(7.15,5.25),gridspec_kw={'width_ratios':[1.55,1]})
 for j,key in enumerate(('negative_area_ps','memory_strength_1_minus_alpha')):
  for chi in COL:
   s=[x for x in summary if x['chirality']==chi];L=np.array([x['L'] for x in s]); y=np.array([x[key+'_mean'] for x in s]);er=np.array([x[key+'_sem'] for x in s]);
   ax[j,0].errorbar(L,y,yerr=er,marker='o',ms=3.8,lw=1.2,capsize=2,color=COL[chi],label=LAB[chi]);f=next(x for x in fits if x['chirality']==chi and x['observable']==key);xx=np.geomspace(1,10,100);ax[j,0].plot(xx,math.exp(np.log(f['predicted_10L'])-f['p_1to5']*math.log(10))*xx**f['p_1to5'],ls='--',lw=.9,color=COL[chi]);
   ax[j,0].scatter([10],[y[-1]],s=38,facecolors='white',edgecolors=COL[chi],zorder=4)
  ax[j,0].set(xscale='log',yscale='log',xlabel='replicated axial length $L$',ylabel=('first negative-lobe area (ps)' if j==0 else r'$1-\alpha_{z,\min}$'))
  if j==0: ax[j,0].legend(title='chirality',fontsize=6,title_fontsize=6,ncol=2,loc='upper left')
  ax[j,0].text(.98,.04,'solid: means; dashed: 1–5L fit; open: 10L',transform=ax[j,0].transAxes,ha='right',fontsize=5.7)
  vals=[next(x for x in fits if x['chirality']==chi and x['observable']==key)['relative_deviation_10L_pct'] for chi in COL]
  ax[j,1].axhline(0,color='#777',lw=.7);ax[j,1].bar(range(4),vals,color=[COL[c] for c in COL],width=.62);ax[j,1].set(xticks=range(4),xticklabels=[LAB[c] for c in COL],ylabel='10L deviation from 1–5L fit (%)',xlabel='chirality')
  for i,v in enumerate(vals):ax[j,1].text(i,v+(2 if v>=0 else -2),f'{v:+.0f}%',ha='center',va='bottom' if v>=0 else 'top',fontsize=6)
 for a,l in zip(ax.flat,'abcd'):a.text(-.16,1.04,l,transform=a.transAxes,fontweight='bold',fontsize=9)
 fig.text(.5,.995,'Weak Nosé–Hoover; no momentum removal; 1–5L: 10 fs/1 ns; 10L: 100 fs/10 ns; direct MSD, 1-decade local slope (0.5–100 ps); normalized total $C_{vJ}$ with instantaneous water-COM velocity removed; error bars = replica SEM',ha='center',va='top',fontsize=5.5)
 fig.subplots_adjust(left=.1,right=.99,bottom=.1,top=.87,wspace=.42,hspace=.45)
 fig.savefig(str(OUT)+'.png',dpi=600,bbox_inches='tight');fig.savefig(str(OUT)+'.tiff',dpi=600,bbox_inches='tight');fig.savefig(str(OUT)+'.pdf',bbox_inches='tight');fig.savefig(str(OUT)+'.svg',bbox_inches='tight')
if __name__=='__main__':main()
