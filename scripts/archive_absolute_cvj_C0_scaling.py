"""Archive the zero-lag absolute-current intensity per water and its L scaling."""
import csv, json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; IN=ROOT/'heartbeat_fetch'/'stage_absolute_cvj_1L10L_20260809'/'output'; OUT=ROOT/'assets'
rows=[]; fits=[]
for chi in ['7_7','8_8','9_9','17_0']:
 vals=[]
 for L in [1,2,3,4,5,10]:
  a=np.array([json.loads(f.read_text())['C0_per_water_A2_ps2'] for f in sorted(IN.glob(f'{chi}_L{L}_rep*_absolute_cvj.json'))]);rows.append({'chirality':chi,'L':L,'C0_per_water_mean_A2_ps2':a.mean(),'C0_per_water_sem_A2_ps2':a.std(ddof=1)/np.sqrt(len(a)) if len(a)>1 else 0.,'n_replicates':len(a)});vals.append((L,a.mean()))
 x=np.log([z[0] for z in vals]);y=np.log([z[1] for z in vals]);p,b=np.polyfit(x,y,1);r=y-(p*x+b);fits.append({'chirality':chi,'range':'1-10L','p_C0_per_water':p,'p_sem':np.sqrt((r@r)/(len(x)-2)/((x-x.mean())@(x-x.mean()))),'R2':1-(r@r)/((y-y.mean())@(y-y.mean()))})
for name,data in [('crosschirality_1L10L_absolute_cvj_C0_per_water.csv',rows),('crosschirality_1L10L_absolute_cvj_C0_per_water_powerlaw.csv',fits)]:
 with (OUT/name).open('w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
