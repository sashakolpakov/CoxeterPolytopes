#!/usr/bin/env python3
"""Test saved two-anchor proposals against the complete H22 cap orbit."""
import json
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h22_borcherds_extension_20261013 import prepare
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h22_orbit_extension_scout_20261013.json'
def run():
    Q,W,G,h,X,A,b,Qn,R=prepare();src=json.loads((ROOT/'data/enumerations/fresh_h22_two_anchor_20261013.json').read_text());orbit=json.loads((ROOT/'data/enumerations/fresh_h22_symmetry_orbit_20261013.json').read_text())
    U=np.load(ROOT/orbit['array'])['twice_old_projections'];target=np.array(next(r['escape']['point']for r in orbit['runs']if r['center_label']==5));T=np.sqrt(17-np.sqrt(5));norm=3-np.sqrt(5)
    center=np.array([np.r_[h,T],np.r_[h,-T]]);finite=np.cos(np.pi/np.arange(2,121));counts=Counter();survivors=[];seen=set()
    def errors(c):
        ix=np.clip(np.searchsorted(finite,c),1,len(finite)-1);e=np.minimum(abs(c-finite[ix]),abs(c-finite[ix-1]));return np.where(c>=1-2e-7,0,e)
    for index,row in enumerate(src['caps']):
        r=np.array(row['normal']);key=tuple(np.round(r/np.linalg.norm(r),9))
        if key in seen:continue
        seen.add(key);counts['distinct_numerical_proposals']+=1
        p=r@Qn@target
        if p<=1e-7:continue
        counts['target_cutters']+=1
        if np.max(errors(-(center@Qn@r)/np.sqrt(2*norm)))>2e-7:continue
        counts['center_compatible_cutters']+=1
        base=U@(Q@r[:22])/2;extra=r[-1]*np.sqrt(2.5)
        e=np.maximum(errors(-(base+extra)/np.sqrt(2)),errors(-(base-extra)/np.sqrt(2)))
        if np.max(e)>2e-7:continue
        counts['whole_orbit_compatible_cutters']+=1
        survivors.append({'source_index':index,'normal':row['normal'],'target_pairing':float(p),'source_center_label':row['center_label']})
    out={'status':'BOUNDED_NUMERICAL_EXTENSION_SCREEN','target_center_label':5,'counts':dict(counts),'survivors':survivors,'scope':'Only1933 saved two-anchor entries (deduplicated numerically), max label120. Tests both signs of all60480 orbit caps. No universal extension obstruction.'};OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['counts'],flush=True)
if __name__=='__main__':run()
