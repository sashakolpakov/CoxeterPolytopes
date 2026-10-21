#!/usr/bin/env python3
"""Numerical whole-bank feasibility scout; no exact classification claim.

The seed is rebuilt from the exact source. Ten old finite vertices are used
as centers. Four transverse labels independently vary over2,3,4,5,6,8,10.
Old-wall labels are numerically screened against2..100 or disjointness.
Every accepted cap is then imposed at once, without compatibility filtering.
"""
import argparse,hashlib,json
from itertools import product
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from fresh_h9_compactification_20261013 import geometry
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'scripts/fresh_h9_compactification_20261013.py'
OUT=ROOT/'results/fresh_h9_compactification_20261013_scout.json'
def run():
    _,_,Q,_,V,_,_,_,_,J,W,c=geometry()
    cv=lambda A:np.array([[x.approx()for x in r]for r in A])
    Q=cv(Q);W=cv(W);V=cv(V);J=cv(J)
    V=V/np.sqrt(-np.einsum('ij,jk,ik->i',V,Q[:6,:6],V))[:,None]
    labels=[2,3,4,5,6,8,10];tv={2:0,**{m:2*np.cos(np.pi/m)for m in labels[1:]}}
    oldcos=np.cos(np.pi/np.arange(2,101));norms=np.diag(W@Q@W.T)
    C=[];recipes=[];nonpositive=0;rejected=0
    for vi,v in enumerate(V):
        for ls in product(labels,repeat=4):
            b=np.array([tv[m]for m in ls]);s=b@J@b-2
            if s<=1e-8:nonpositive+=1;continue
            r=np.r_[np.sqrt(s)*v,-J@b];co=-(W@Q@r)/np.sqrt(2*norms)
            if any(z<-1e-7 or(z<1-1e-7 and min(abs(z-oldcos))>1e-7)for z in co):rejected+=1;continue
            C.append(r);recipes.append([vi,list(ls)])
    # Robust explicit strict point also independently certified for312 profiles
    # by the exact script. Here the check of the widerbank is numerical only.
    E=np.array([30861,39915,29906,19910,9932,-705,-7565,-10211,-10432,-10624],float)
    C=np.array(C);seedmax=float(max(W@Q@E));capmax=float(max(C@Q@E));q=float(E@Q@E)
    assert seedmax<0 and capmax<0 and q>0
    out={'status':'NUMERICAL_WHOLE_BANK_STILL_HAS_STRICT_SPACELIKE_GAP',
      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'transverse_label_alphabet':labels,'old_finite_angle_screen_labels':[2,100],
      'profile_count':10*len(labels)**4,'nonpositive_old_scale_count':nonpositive,
      'rejected_count':rejected,'accepted_count':len(C),'accepted_recipes':recipes,
      'point':E.tolist(),'point_square':q,'maximum_seed_pairing':seedmax,'maximum_cap_pairing':capmax,
      'scope':'Numerical bounded profile scout only; no certified all-angle classification, pair compatibility, or exact wider-bank exclusion.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('WROTE',OUT);print({k:out[k]for k in ['profile_count','accepted_count','point_square','maximum_cap_pairing']})
if __name__=='__main__':run()
