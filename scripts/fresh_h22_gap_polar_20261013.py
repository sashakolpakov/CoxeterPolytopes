#!/usr/bin/env python3
"""Numerical local polar diagnosis for the certified H22 cusp-orbit cone."""
import json
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h22_gap_polar_20261013.json'

def read(p): return json.loads((ROOT/p).read_text())

def setup():
    raw=read('data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json')['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64)
    W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    h=W.sum(axis=0);h//=np.gcd.reduce(h)
    U=np.load(ROOT/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']
    a=np.load(ROOT/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz')
    D=a['old_null_projections'];t=a['four_heights']/4
    return Q,W,h,U,D,t

def run():
    Q,W,h,U,D,t=setup()
    src=read('data/enumerations/fresh_h22_cusp_orbits_20261013.json')
    x=np.array(src['runs'][0]['escape']['point']);q=x[:22]@Q@x[:22]+x[22]**2
    oldnorm=np.diag(W@Q@W.T)
    groups=[('old',W@Q@x[:22],oldnorm),
            ('timelike+',U@Q@x[:22]/2+np.sqrt(2.5)*x[22],np.ones(len(U))),
            ('timelike-',U@Q@x[:22]/2-np.sqrt(2.5)*x[22],np.ones(len(U))),
            ('null+',D@Q@x[:22]+t*x[22],t*t),
            ('null-',D@Q@x[:22]-t*x[22],t*t),
            ('center',np.array([h@Q@x[:22]+4*x[22],h@Q@x[:22]-4*x[22]]),np.array([2,2]))]
    labels=[2,3,4,5,6,8,10,12];c=np.cos(np.pi/np.array(labels))
    rows=[];active=[];bad=[]
    for name,p,norm in groups:
        values=-p/np.sqrt(q*norm)
        finite=np.flatnonzero(values<1-1e-7)
        errors=np.array([min(abs(v-c))for v in values[finite]])
        wrong=finite[errors>1e-7]
        ids=np.flatnonzero(abs(values)<1e-7)
        active.extend((name,int(i))for i in ids)
        bad.extend((name,int(i),float(values[i]))for i in wrong)
        rows.append({'group':name,'walls':len(p),'active':ids.tolist(),
                     'finite_pairs':len(finite),'forbidden_in_selected_alphabet':len(wrong),
                     'minimum_normalized_magnitude':float(min(values))})
    return {'status':'NUMERICAL_POLAR_DIAGNOSIS','point':x.tolist(),'norm':float(q),
            'active':active,'groups':rows,'bad_pairs':bad,'tested_labels':labels,
            'scope':'Floating-point proposal diagnosis, not an exact all-label obstruction.'}

if __name__=='__main__':
    d=run();OUT.write_text(json.dumps(d,indent=2)+'\n')
    print('norm',d['norm'],'active',d['active'])
    for g in d['groups']:print({k:v for k,v in g.items()if k!='active'})
    print('first bad',d['bad_pairs'][:30])
