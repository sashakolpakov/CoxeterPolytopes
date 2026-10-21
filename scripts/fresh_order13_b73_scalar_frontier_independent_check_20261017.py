#!/usr/bin/env python3
"""Preserved independent check first run inline as terminal fcd4f6.

Optional-part enumeration, sliding convolution and inclusion-exclusion CDF.
No producer arithmetic imports; checks the entire frozen scalar certificate,
not geometric payment proofs or a new ambient bound.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import prod, comb
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order13_b73_scalar_frontier_20261017'


def run():
    j=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in j['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    repeated=tuple(range(10,18));domain=set()
    def visit(left,k=0,p=()):
        if k==8:
            if left==0:yield p
            return
        r=repeated[k]
        for n in range(left//r+1):yield from visit(left-n*r,k+1,p+(r,)*n)
    for low in range(10):
        for high in (0,*range(18,65)):
            remain=64-low-high
            if remain<0:continue
            for p in visit(remain):
                q=tuple(sorted(p+((low,) if low else())+((high,) if high else())))
                assert q not in domain;domain.add(q)
    assert len(domain)==1650
    rows={tuple(r['ranks']):r for r in j['rows']};assert set(rows)==domain
    cs={11:F(16,3),12:F(74,13),13:F(6),14:F(32,5),15:F(107,16),16:F(120,17)}
    computed={}
    for p in sorted(domain):
        h=[1]
        for r in p:
            out=[];w=0
            for k in range(len(h)+r):
                if k<len(h):w+=h[k]
                if 0<=k-r-1<len(h):w-=h[k-r-1]
                out.append(w)
            h=out
        V=prod(r+1 for r in p);G=sum(v*min(k,64-k) for k,v in enumerate(h));sh={0:1}
        for r in p:
            z=dict(sh)
            for a,b in sh.items():z[a+r+1]=z.get(a+r+1,0)-b
            sh={a:b for a,b in z.items() if b}
        assert G==32*V-2*sum(a*comb(32-k+len(p),len(p)+1) for k,a in sh.items()
                            if 32-k+len(p)>=len(p)+1)
        cost=sum(F(r,2) if r<=10 else cs[r] if r<=16 else F(163,23) for r in p)
        m=F(5*G-64*V,2*V)-cost;computed[p]=m
        assert rows[p]==dict(ranks=list(p),link_vertices=V,G=G,h_polynomial=h,margin=str(m))
    assert sum(m<0 for m in computed.values())==7
    assert sum(m-F(p.count(12),13)<0 for p,m in computed.items())==16
    for sc in j['conditional_scenarios']:
        drops={int(k)-1:v for k,v in sc['additional_total_factor_cost_drops'].items()}
        vals={p:m+sum(F(p.count(k)*v,k+1) for k,v in drops.items()) for p,m in computed.items()}
        assert {tuple(r['ranks']):F(r['margin']) for r in sc['rows']}==vals
        assert sc['negative_profiles']==sum(m<0 for m in vals.values())
        best=min(vals,key=lambda p:(vals[p],p))
        assert sc['minimum']==dict(ranks=list(best),margin=str(vals[best]))
    print('PASS_INDEPENDENT_ORDER13_B73_SCALAR',len(domain),
          'profiles; all3scenario rows; optional-part multiplicities, sliding polynomial, independentCDF; no producer imports',flush=True)


if __name__=='__main__':run()
