#!/usr/bin/env python3
"""Find exact spacelike witnesses in diagonal H9 prefix cones.

SciPy only proposes points. Acceptance uses integer quadratic-field arithmetic.
The final verifier does not depend on SciPy or trust its vertex enumeration.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from search_quadratic_diagonal_h9_2026_10_13 import Field, Z, add, sub, scale

def verify_witness(row,w):
    K=Field(row['d']); alpha=tuple(row['alpha']); ys=w['spatial_numerators']; den=w['denominator']; n=row['dimension']
    assert den>0 and len(ys)==n and all(isinstance(a,int) for a in ys)
    assert all(ys[i]>=ys[i+1] for i in range(n-1)) and ys[-1]>=0
    for r in row['walls']:
        val=Z
        for a,b in zip(ys,r['spatial']):val=add(val,scale(b,a))
        assert K.sign(sub(val,scale(K.mul(alpha,r['x']),den)))<=0
    norm=sub((sum(a*a for a in ys),0),scale(alpha,den*den))
    assert K.sign(norm)>0
    return norm

def propose(row):
    import numpy as np
    from scipy.optimize import linprog
    K=Field(row['d']); n=row['dimension']; alpha=tuple(row['alpha']); av=K.val(alpha)
    if not row['walls']:return None
    A=np.array([np.cumsum([K.val(v) for v in r['spatial']]) for r in row['walls']]); b=np.array([K.val(K.mul(alpha,r['x'])) for r in row['walls']])
    # Bounded vertex scout; exhaustive double description becomes impractical
    # for the hundreds of walls present in the deeper prefixes.
    rng=np.random.default_rng(20261013)
    transform=np.triu(np.ones((n,n)))
    verts=[]
    for trial in range(160):
        objective=rng.normal(size=n) if trial>=n else np.eye(n)[trial]
        for step in range(8):
            result=linprog(-objective,A_ub=A,b_ub=b,bounds=(0,None),method='highs')
            if not result.success: raise RuntimeError(result.message)
            z=result.x
            verts.append(z)
            updated=transform.T@transform@z
            if np.linalg.norm(updated-objective)<1e-10:break
            objective=updated
    scores=np.array([sum(np.cumsum(z[::-1])**2)/av for z in verts])
    for idx in np.argsort(scores)[::-1]:
        if scores[idx]<=1+1e-8:break
        for shrink in (1e-8,1e-6,1e-4,1e-3,.01):
            zd=np.floor(np.maximum(verts[idx],0)*(1-shrink)*10**9).astype('int64')
            ys=np.cumsum(zd[::-1])[::-1].tolist()
            w={'spatial_numerators':ys,'denominator':10**9}
            try:
                norm=verify_witness(row,w)
            except AssertionError:continue
            return {'exact_feasible_spacelike_point':w,'norm_numerator':norm,'numerical_vertex_trials':len(verts),'numerical_max_norm_ratio':float(max(scores))}
    return {'numerical_vertex_trials':len(verts),'numerical_max_norm_ratio':float(max(scores)),'needs_exact_complete_cone_verification':True}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('data/enumerations/quadratic_diagonal_h9_2026_10_13.json'))
    parser.add_argument('--output',type=Path,default=Path('data/enumerations/quadratic_diagonal_h9_geometry_2026_10_13.json'))
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args(); data=json.loads(args.input.read_text()); out=[]
    if args.verify_only:
        audits=json.loads(args.output.read_text())
        for audit in audits:
            row=next(r for r in data['records'] if r['d']==audit['d'] and r['alpha']==audit['alpha'])
            assert list(verify_witness(row,audit['exact_feasible_spacelike_point']))==audit['norm_numerator']
        print(f'Verified {len(audits)} exact spacelike witnesses without importing SciPy.');return
    for row in data['records']:
        if row['spacelike_witnesses']:continue
        result=propose(row)
        print(row['d'],row['alpha'],result,flush=True)
        out.append({'d':row['d'],'alpha':row['alpha'],**result})
        args.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
