"""Test whether n=1 first-negative-lobe timing remains linear in L at 10L."""
from __future__ import annotations
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; SRC=ROOT/'remote_fetch'/'output'; AS=ROOT/'assets'
CASE=AS/'crosschirality_1L10L_lobe_timing_case.csv'; SUM=AS/'crosschirality_1L10L_lobe_timing_summary.csv'; FIT=AS/'crosschirality_1L10L_lobe_timing_linearity_tests.csv'; OUT=AS/'crosschirality_1L10L_lobe_timing_linearity_nature'
COL={'7_7':'#3B7EA1','8_8':'#D88737','9_9':'#5E9C76','17_0':'#8C6BB1'};LAB={'7_7':'(7,7)','8_8':'(8,8)','9_9':'(9,9)','17_0':'(17,0)'}
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'axes.linewidth':.8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
def write(p,rows):
 with p.open('w',newline='',encoding='utf8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sem(a):return float(np.std(a,ddof=1)/math.sqrt(len(a))) if len(a)>1 else float('nan')
def main():
 AS.mkdir(exist_ok=True);rec=[]
 for p in sorted(SRC.glob('*.json')):
  x=json.loads(p.read_text());chi,L,rep=x['case_id'].rsplit('_',2);l=x['first_negative_lobe'];rec.append({'chirality':chi,'L':int(L[1:]),'replicate':int(rep[3:]),'t_start_ps':l['t_start_ps'],'t_min_ps':l['t_min_ps'],'t_end_ps':l['t_end_ps'],'width_ps':l['t_end_ps']-l['t_start_ps'],'dt_ps':x['dt_ps']})
 assert len(rec)==56;write(CASE,rec);summ=[]
 for chi in COL:
  for L in sorted({r['L'] for r in rec if r['chirality']==chi}):
   rr=[r for r in rec if r['chirality']==chi and r['L']==L];o={'chirality':chi,'L':L,'n_replicates':len(rr)}
   for key in ('t_start_ps','t_min_ps','t_end_ps','width_ps'):
    a=np.array([r[key] for r in rr]);o[key+'_mean']=a.mean();o[key+'_sem']=sem(a)
   summ.append(o)
 write(SUM,summ);fits=[]
 for chi in COL:
  r=[x for x in summ if x['chirality']==chi]
  for key in ('t_start_ps','t_min_ps','t_end_ps','width_ps'):
   s=[x for x in r if x['L']<=5];X=np.array([x['L'] for x in s]);Y=np.array([x[key+'_mean'] for x in s]);b,a=np.polyfit(X,Y,1);res=Y-(a+b*X);se=math.sqrt(float((res@res)/(len(X)-2)/((X-X.mean())@(X-X.mean()))));r2=1-float(res@res)/float(((Y-Y.mean())@(Y-Y.mean())));obs=next(x[key+'_mean'] for x in r if x['L']==10);pred=a+10*b
   fits.append({'chirality':chi,'metric':key,'slope_ps_per_L_1to5':b,'slope_sem':se,'intercept_ps':a,'r_squared_1to5':r2,'predicted_10L_ps':pred,'observed_10L_ps':obs,'deviation_10L_ps':obs-pred,'relative_deviation_10L_pct':100*(obs/pred-1)})
 write(FIT,fits)
 fig,ax=plt.subplots(2,2,figsize=(7.15,5.15),gridspec_kw={'width_ratios':[1.52,1]})
 for j,(key,title) in enumerate((('t_start_ps','first zero crossing $t_-$'),('t_min_ps','first trough $t_{min}$'),('t_end_ps','return-to-zero $t_+$'))):
  a=ax.flat[j]
  for chi in COL:
   s=[x for x in summ if x['chirality']==chi];L=np.array([x['L'] for x in s]);y=np.array([x[key+'_mean'] for x in s]);e=np.array([x[key+'_sem'] for x in s]);a.errorbar(L,y,yerr=e,marker='o',ms=3.8,lw=1.2,capsize=2,color=COL[chi],label=LAB[chi]);f=next(x for x in fits if x['chirality']==chi and x['metric']==key);a.plot([1,10],[f['intercept_ps']+f['slope_ps_per_L_1to5'],f['predicted_10L_ps']],ls='--',lw=.85,color=COL[chi]);a.scatter([10],[y[-1]],s=38,facecolors='white',edgecolors=COL[chi],zorder=4)
  a.set(xlabel='replicated axial length $L$',ylabel='time (ps)',title=title)
  if j==0:a.legend(title='chirality',fontsize=6,title_fontsize=6,ncol=2,loc='upper left')
  a.text(.98,.05,'dashed: 1–5L linear fit; open: 10L',transform=a.transAxes,ha='right',fontsize=5.7)
 a=ax.flat[3]
 for i,chi in enumerate(COL):
  vals=[next(x for x in fits if x['chirality']==chi and x['metric']==k)['relative_deviation_10L_pct'] for k in ('t_start_ps','t_min_ps','t_end_ps')];a.plot([0,1,2],vals,marker='o',color=COL[chi],lw=1.1,label=LAB[chi])
 a.axhline(0,color='#777',lw=.7);a.set(xticks=[0,1,2],xticklabels=[r'$t_-$',r'$t_{min}$',r'$t_+$'],ylabel='10L deviation from 1–5L line (%)',xlabel='timing metric');a.legend(fontsize=5.8,ncol=2,loc='upper left')
 for aa,l in zip(ax.flat,'abcd'):aa.text(-.15,1.04,l,transform=aa.transAxes,fontweight='bold',fontsize=9)
 fig.text(.5,.995,'Weak Nosé–Hoover; no momentum removal; 1–5L: 10 fs/1 ns; 10L: 100 fs/10 ns; normalized total $C_{vJ}$; lobe boundaries by interpolated zero crossings; error bars = replica SEM.',ha='center',va='top',fontsize=5.6)
 fig.subplots_adjust(left=.09,right=.99,bottom=.1,top=.87,wspace=.38,hspace=.42)
 fig.savefig(str(OUT)+'.png',dpi=600,bbox_inches='tight');fig.savefig(str(OUT)+'.tiff',dpi=600,bbox_inches='tight');fig.savefig(str(OUT)+'.pdf',bbox_inches='tight');fig.savefig(str(OUT)+'.svg',bbox_inches='tight')
if __name__=='__main__':main()
