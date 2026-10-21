#!/usr/bin/env python3
"""Independent finite arithmetic and proof-input review of the mixed-norm gap.

No producer arithmetic import. The independently reviewed primitive norm
alphabet and signatures are inherited, not enumerated again here.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta11_mixed_norm_gap_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def norm(a,b):return a*a+a*b-b*b
def conjugate(a,b):return a+b,-b
def sign(a,b):
    u,v=2*a+b,b
    s=lambda x:(x>0)-(x<0)
    if not u:return s(v)
    if not v or s(u)==s(v):return s(u)
    return s(u)*s(u*u-5*v*v)
def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for s,h in d.get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):read(s)
        return d
    primary=f'results/{TAG}.json';d=read(primary)
    poolpath='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'
    pool=read(poolpath)
    assert pool['status']=='PASS_INDEPENDENT_FREE_HOST_AND_NORM_ALPHABET_INPUTS'
    assert pool['signatures']==[[9,1],[10,0]]
    assert pool['normalized_norms']==[['2','0'],['6','2']]
    assert pool['beta']==['3','1'] and pool['beta_norm']=='11'
    beta=(F(3),F(1));betabar=conjugate(*beta)
    assert norm(*beta)==11 and norm(*betabar)==11
    assert sign(*beta)>0 and sign(*betabar)>0
    assert sign(beta[0]-4,beta[1])>0
    assert sign(betabar[0]-4,betabar[1])<0
    assert d['beta']==['7/2','1/2'] and d['conjugate_beta']==['7/2','-1/2']
    # Norm(z) is a nonzero integer. 11*Norm(z)^2<16 admits only ±1.
    possible_norm_z=[n for n in range(-2,3)if n and 11*n*n<16]
    assert possible_norm_z==[-1,1] and 11*2*2>=16
    traces=[t for t in range(1,8)if t*t>=44]
    assert traces==[7]
    # Independently solve integral a=m+n*phi with trace7 and norm11.
    # Discriminant5*n²=7²-44 bounds |n|<=1.
    assert 7*7-44==5
    solutions=[]
    for n in range(-1,2):
        m=F(7-n,2)
        if m.denominator==1 and norm(m,n)==11:
            a=(m,F(n));ap=conjugate(*a)
            assert sign(*a)>0 and sign(*ap)>0
            assert sign(a[0]-4,a[1])>0 or sign(ap[0]-4,ap[1])>0
            solutions.append([str(m),str(n)])
    assert solutions==[['4','-1'],['3','1']]
    assert d['hypothetical_norm_choices']==[11] and d['hypothetical_trace_choices']==[7]
    paths=[primary,poolpath,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_UNIVERSAL_MIXED_NORM_GAP_REVIEW',
        sources={p:sha(p)for p in paths},recursive_manifests=len(seen),
        arithmetic='Standalone Fraction coefficient norm/trace arithmetic; no producer imports.',
        possible_norm_z=possible_norm_z,possible_trace=traces,
        integral_candidates=solutions,all_candidates_violate_embedding_bound=True,
        proof_checks=['Primitive long-root reflection gives Gr_long in beta*O^10, so pairing=beta*z.',
          'Normbeta11 excludes proportional field root lines; injectivity preserves independence at the second embedding.',
          'Conjugate positive definiteness makes the Cauchy bound strict:0<a_conjugate<4.',
          'Nonorthogonal intersection or parallelism gives0<a<=4. Norm and trace contradiction includes equality at4.',
          'Only the fixed integral lattice and its complete normalized reflection alphabet are covered.'],
        scope='Universal orthogonal-or-strictly-disjoint statement for the two primitive integral norm classes. Inherits prior independent lattice/signature/alphabet proof. No root enumeration, general geometric-cap restriction or reflectivity verdict.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
