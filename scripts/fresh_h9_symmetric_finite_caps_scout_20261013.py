#!/usr/bin/env python3
"""Numerical scout for a declared larger finite leaf-angle alphabet."""
import json
import math
from itertools import combinations_with_replacement
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h9_symmetric_finite_caps_scout_20261013.json'

def run():
    source=json.loads((ROOT/'data/enumerations/fresh_h9_symmetric_balanced_caps_20261013.json').read_text())
    phi=(1+math.sqrt(5))/2
    value=lambda x:float(__import__('fractions').Fraction(x[0]))+phi*float(__import__('fractions').Fraction(x[1]))
    G=np.array([[value(x)for x in row]for row in source['gram']]);basis=[0,1,*range(2,18,2)]
    I=np.linalg.inv(G[np.ix_(basis,basis)]);E=np.array([1,-1,*([0]*8)]);A=float(E@I@E)
    labels=[2,3,4,5,6,8,10,12,18];values={m:0. if m==2 else 2*math.cos(math.pi/m)for m in labels}
    pair_sums=[(values[l]+values[r],l,r)for l,r in combinations_with_replacement(labels,2)]
    min_gap=min(abs(x[0]-y[0])for i,x in enumerate(pair_sums)for y in pair_sums[:i])
    assert min_gap>1e-8  # Distinct numerically; this scout makes no exact completeness claim.
    b=phi+.5;p=phi+1;Y=np.array([value(x)for x in source['escape_coefficients']])
    S=np.array([value(x)for x in source['strict_center_coefficients']])
    old=[np.array([value(x)for x in source['caps'][i]['coefficients']])for i in source['selected_caps']]
    def gate(pair):
        if pair>1e-8:return None
        if pair < -2-1e-8:return 'dotted'
        if abs(pair+2)<1e-8:return 'parallel'
        return next((m for m in labels if abs(pair+values[m])<1e-8),None)
    rows=[];counts={'orbits':0,'nonreal':0,'positive_center':0,'failed_center_angle':0,'accepted':0}
    for total,left,right in pair_sums:
        for k in(range(9)if left!=right else(8,)):
            counts['orbits']+=1
            fixed=np.array([0,-p*total/b,*[-(values[left]if i<k else values[right])for i in range(8)]])
            B=float(2*E@I@fixed);C=float(fixed@I@fixed-2);disc=B*B-4*A*C
            if disc < -1e-8:counts['nonreal']+=1;continue
            for branch in((1,)if abs(disc)<1e-8 else(-1,1)):
                t=(-B+branch*math.sqrt(max(0,disc)))/(2*A);cov=fixed+t*E
                if max(cov[:2])>1e-8:counts['positive_center']+=1;continue
                center_gates=[gate(x)for x in cov[:2]]
                if any(x is None for x in center_gates):counts['failed_center_angle']+=1;continue
                counts['accepted']+=1;r=np.zeros(18);r[basis]=I@cov
                assert abs(r@G@r-2)<1e-7 and r@G@S<0
                rows.append({'leaf_labels':[left,right],'left_count':k,'branch':branch,
                    'quadratic':[A,B,C],'discriminant':disc,'center_pairings':list(cov[:2]),
                    'center_gates':center_gates,'normal_coefficients':list(r),
                    'saved_escape_pairing':float(r@G@Y),'selected_cap_pairings':[float(r@G@v)for v in old],
                    'selected_cap_gates':[gate(float(r@G@v))for v in old]})
    out={'status':'NUMERICAL_SCOUT_ONLY','leaf_labels':labels,'center_labels':labels,
        'dotted_and_parallel_centers_allowed':True,'minimum_numerical_pair_sum_gap':min_gap,
        'counts':counts,'survivors':rows,
        'scope':'Floating-point scout over one representative per leaf pair-swapping orbit for the declared alphabet. Pair-sum distinctness, norm roots, and angle gates are numerical. Orbit representatives are not a complete placement search against the selected two caps. No exact geometric or global nonexistence claim.'}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'counts':counts,'survivors':[{k:v for k,v in row.items()if k!='normal_coefficients'}for row in rows]}));return out

if __name__=='__main__':run()
