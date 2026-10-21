#!/usr/bin/env python3
"""Independent complete profile census from archived corner Gram matrices.

This does not import the cap producer or its profile arithmetic. The separate
all-height checker is replayed independently and reviewed in the report.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import hashlib
import json
from math import lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_h9_component_bridge_20261013'
OUT=ROOT/'results/fresh_h9_component_bridge_20261013_continuous_review.json'
def read(name):return json.loads((ROOT/name).read_text())
def inv(M):
    n=len(M);a=[list(r)+[Q(i==j)for j in range(n)]for i,r in enumerate(M)]
    for i in range(n):
        k=next(k for k in range(i,n)if a[k][i]);a[i],a[k]=a[k],a[i]
        p=a[i][i];a[i]=[x/p for x in a[i]]
        for j in range(n):
            if j!=i:
                p=a[j][i];a[j]=[x-p*y for x,y in zip(a[j],a[i])]
    return [r[n:]for r in a]
def sg(x):return(x>0)-(x<0)
def gs(a,b):
    u=2*a+b;v=b
    if not u:return sg(v)
    if not v or sg(u)==sg(v):return sg(u)
    return sg(u)*sg(u*u-5*v*v)
def gm(a,b,c,d):return a*c+b*d,a*d+b*c+b*d
def sign(v):
    a,b,c,d=v;x=gs(a,b);y=gs(c,d)
    if not x:return y
    if not y or x==y:return x
    u,v=gm(a,b,a,b);w,z=gm(c,d,c,d)
    return x*gs(u-2*w,v-2*z)

ALPH=((0,0,0,0),(1,0,0,0),(0,0,1,0),(0,1,0,0))
def mul(x,y):
    a,b=gm(*x[:2],*y[:2]);c,d=gm(*x[2:],*y[2:])
    e,f=gm(*x[:2],*y[2:]);g,h=gm(*x[2:],*y[:2])
    return a+2*c,b+2*d,e+g,f+h
P=[[mul(x,y)for y in ALPH]for x in ALPH]

def formal_dual_identity():
    # Independent polynomial identity, with no floating point or CAS.
    def var(i):return {tuple(int(j==i)for j in range(6)):1}
    def add(p,q):
        r=p.copy()
        for k,v in q.items():r[k]=r.get(k,0)+v
        return {k:v for k,v in r.items()if v}
    def neg(p):return {k:-v for k,v in p.items()}
    def mul(p,q):
        r={}
        for a,x in p.items():
            for b,y in q.items():
                k=tuple(u+v for u,v in zip(a,b));r[k]=r.get(k,0)+x*y
        return {k:v for k,v in r.items()if v}
    def sub(p,q):return add(p,neg(q))
    a,b,c,e,x,y=map(var,range(6));D=sub(mul(a,e),mul(b,c))
    l1=sub(mul(c,y),mul(e,x));l2=sub(mul(b,x),mul(a,y))
    assert not add(add(mul(a,l1),mul(c,l2)),mul(D,x))
    assert not add(add(mul(b,l1),mul(e,l2)),mul(D,y))

@lru_cache(None)
def census(matrix):
    M=[list(r)for r in matrix];Mi=inv(M)
    assert all(x>=0 for row in Mi for x in row)
    scale=lcm(*(x.denominator for row in Mi for x in row))
    A=[[int(x*scale)for x in row]for row in Mi]
    # No branch pruning: enumerate every one of the 4^8 assignments.
    # Add precomputed four-coordinate energies of the two half-profiles.
    def halves(offset):
        ans=[]
        for codes in product(range(4),repeat=4):
            v=[0]*4
            for i in range(4):
                for j in range(4):
                    c=A[i+offset][j+offset]
                    if c:
                        for k,z in enumerate(P[codes[i]][codes[j]]):v[k]+=c*z
            ans.append((codes,tuple(v)))
        return ans
    cross=[(i,j,A[i][j+4]*2)for i in range(4)for j in range(4)if A[i][j+4]]
    accepted=[]
    for u,x in halves(0):
        for v,y in halves(4):
            z=[a+b for a,b in zip(x,y)]
            for i,j,c in cross:
                for k,w in enumerate(P[u[i]][v[j]]):z[k]+=c*w
            if sign((2*scale-z[0],-z[1],-z[2],-z[3]))>0:accepted.append(u+v)
    return tuple(sorted(accepted))

def run(verify=False):
    corners=read(BASE+'_verify.json');claim=read(BASE+'_continuous_verify.json')
    rows=[];tail=[];all_profiles={}
    for i,row in enumerate(claim['edges']):
        edge=row['edge']
        v=next(v for v in corners['actual_new_finite_vertices']if set(edge)<=set(v['incident']))
        inds=[v['incident'].index(n)for n in edge]
        M=[]
        for a in inds:
            r=[]
            for b in inds:
                cell=v['gram'][a][b];x,y=Q(cell['a']),Q(cell['b']);assert y==0;r.append(x)
            M.append(tuple(r))
        M=tuple(M);accepted=census(M);all_profiles[i]=set(accepted)
        assert len(accepted)==row['profiles']
        for j in range(8):
            if all(M[j][k]==0 for k in range(8)if k!=j):tail.append((i,edge[j]))
        rows.append({'edge':edge,'profiles':len(accepted),'assignments_tested':4**8,
                     'profile_sha256':hashlib.sha256(json.dumps(accepted).encode()).hexdigest()})
    assert len(rows)==16 and sum(r['profiles']for r in rows)==940
    assert tail==[(t['edge'],t['wall'])for t in claim['uncovered_I2_tails']]and len(tail)==51
    survivors=[]
    for r in claim['surviving_acute_branches']:
        codes=tuple((2,3,4,5).index(m)for m in r['labels'])
        assert codes in all_profiles[r['edge']]
        intervals=[tuple(map(Q,x))for x in r['necessary_ranges']]
        assert all(0<a<=b for a,b in intervals)
        s=int(claim['interval_scale'])
        assert len(r['seed6_subinterval_pairings'])==32*len(intervals)
        for v in r['seed6_subinterval_pairings']:
            a,b=map(int,v);assert -s<a<=b<0
        survivors.append({'edge':r['edge'],'labels':r['labels'],'checks':len(r['seed6_subinterval_pairings'])})
    assert len(survivors)==4
    formal_dual_identity()
    tail_claim=read(BASE+'_dihedral_verify.json')
    assert tail_claim['tail_cases']==51
    assert set(tail)=={(r['edge'],r['isolated_wall'])for r in tail_claim['certificates']}
    for row in tail_claim['certificates']:
        assert len(row['support'])==len(set(row['support']))==2
        assert all(int(x[0])>0 for x in row['dual_weights'])
        assert int(row['uniform_target_pairing_lower_bound'][0])>0
        a,b=map(int,row['nonzero_determinant']);assert a>0 or b<0
    sources=[BASE+'_verify.json',BASE+'_continuous_verify.json',
             BASE+'_dihedral_verify.json',
             'scripts/fresh_h9_component_bridge_20261013_continuous_verify.py',
             'scripts/fresh_h9_component_bridge_20261013_dihedral_verify.py']
    out={'status':'PASS_INDEPENDENT_EXHAUSTIVE_PROFILE_REVIEW',
         'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},
         'assignment_coverage':16*4**8,'distinct_gram_matrices':census.cache_info().misses,
         'assignments_evaluated':census.cache_info().misses*4**8,'profiles':940,'uncovered_isolated_tails':51,
         'edges':rows,'surviving_branches':survivors,
         'dihedral_review':{'cases':51,'formal_dual_cancellation':'PASS',
                            'archived_positive_weights_and_constants':'PASS'},
         'scope':'Independent integer/Fraction arithmetic exhausts all4^8 profiles on each of16 archived edge Grams, without partial-profile pruning. The all-height projection and interval proof is separately source-reviewed and its original exact checker replayed; this is not a second implementation of its complete geometry or interval engine.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['assignment_coverage'],out['profiles'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
