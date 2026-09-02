"""Extract exact mass-weighted total water Pz from LAMMPS full-water dumps."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

def one(path: Path, oxygen_type: int, hydrogen_type: int):
    rec=[]
    with path.open() as f:
        while True:
            marker=f.readline()
            if not marker: break
            if marker.strip()!='ITEM: TIMESTEP': raise ValueError(f'bad marker {marker!r}')
            step=int(f.readline()); f.readline(); n=int(f.readline()); f.readline()
            for _ in range(3): f.readline()
            header=f.readline().split()[2:]; c={v:i for i,v in enumerate(header)}
            if 'type' not in c or 'vz' not in c: raise ValueError(f'missing type/vz: {header}')
            a=np.fromstring(' '.join(f.readline() for _ in range(n)),sep=' ').reshape(n,len(header))
            typ=a[:,c['type']].astype(int); unknown=set(np.unique(typ))-{oxygen_type,hydrogen_type}
            if unknown: raise ValueError(f'unexpected atom types: {unknown}')
            m=np.where(typ==oxygen_type,15.999,1.008)
            rec.append({'step':step,'time_ps':step*.0005,'Pz_total_water_amu_A_fs':float(np.dot(m,a[:,c['vz']]))})
    return pd.DataFrame(rec)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--dumps',nargs='+',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--oxygen-type',type=int,default=1); p.add_argument('--hydrogen-type',type=int,default=2); a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True); rows=[]
    for i,x in enumerate(a.dumps,1):
        q=one(x,a.oxygen_type,a.hydrogen_type); q.insert(0,'case',f'rep{i}'); q.insert(1,'source_dump',str(x)); rows.append(q)
    data=pd.concat(rows,ignore_index=True); data.to_csv(a.output/'fullwater_Pz_timeseries.csv',index=False)
    (a.output/'metadata.json').write_text(json.dumps({'definition':'sum all water atoms m_a*v_za','dumps':[str(x) for x in a.dumps],'oxygen_type':a.oxygen_type,'hydrogen_type':a.hydrogen_type},indent=2))
    (a.output/'EXTRACTION_FINISHED.txt').write_text('Full-water Pz extraction finished successfully.\n')
if __name__=='__main__': main()
