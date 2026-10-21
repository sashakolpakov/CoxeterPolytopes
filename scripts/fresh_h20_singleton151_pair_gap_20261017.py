#!/usr/bin/env python3
"""All-label gap for one singleton151 endpoint pair; no family closure claim."""
import argparse
import hashlib
import json
from fractions import Fraction as F
import sympy as s
import fresh_h20_e8_exact_cluster_20261017 as P
import fresh_h20_e8_endpoint_belt_20261017_exact_gate as E

ROOT=P.ROOT
TAG='fresh_h20_singleton151_pair_gap_20261017'
END='results/fresh_h20_e8_partner_endpoints_20261017.json'

def run():
    C=json.loads((ROOT/P.CUSP).read_text())
    D=json.loads((ROOT/END).read_text())
    ids=C['core_ids']
    G=s.Matrix(C['core_gram'])
    GI=G.inv()
    def corner(edge):
        missing=next(j for j in (159,160,161) if j not in edge)
        mi=ids.index(missing)
        ix=[i for i in range(21) if i!=mi]
        J=G.extract(ix,ix).inv()
        V=(-GI[:,mi]/s.sqrt(-GI[mi,mi])).applyfunc(s.simplify)
        assert s.simplify((V.T*G*V)[0])==-1
        return ix,J,V
    ix,J,V=corner((159,160))
    f=s.zeros(20,1); f[ix.index(ids.index(151))]=-1
    p=s.zeros(21,1)
    for i,x in zip(ix,J*f): p[i]=x
    assert s.simplify((p.T*G*p)[0])==8
    assert s.simplify((p.T*G*V)[0])==0
    row=D['endpoint_censuses'][2]['profiles'][37]
    assert row['contacts']==[[151,4],[161,3]]
    jx,K,W=corner((160,161))
    g=s.zeros(20,1)
    for j,m in row['contacts']:g[jx.index(ids.index(j))]=-2*s.cos(s.pi/m)
    q=s.zeros(21,1)
    for j,x in zip(jx,K*g):q[j]=x
    rad=s.simplify((q.T*G*q)[0]-2)
    B=(q+s.sqrt(rad)*W).applyfunc(s.simplify)
    assert s.simplify((B.T*G*B)[0])==2
    alpha=s.simplify(-(p.T*G*B)[0]/2)
    beta=s.simplify(-(V.T*G*B)[0]/2)
    assert s.simplify(alpha+4*s.sqrt(2))==0
    a,b=P.alg(alpha),P.alg(beta)
    assert b[0]>0
    def value(d):return P.add(P.mul(a,d),P.mul(b,P.sqrt(P.sub(P.mul(P.pt(8),P.mul(d,d)),P.pt(2)))))
    deriv=P.add(a,P.div(P.mul(P.pt(16),b),P.sqrt(P.pt(30))))
    assert deriv[0]>0
    vals=[]
    d3=P.pt(1);v3=value(d3);assert v3[1]<0
    vals.append(dict(label151=3,cosine=P.save(v3),result='obtuse'))
    for m,gap in [(4,2),(5,2),(6,3),(7,3),(8,3),(9,3)]:
        v=value(P.mul(P.pt(2),E.cosine(m)))
        lo,hi=E.cosine(gap),E.cosine(gap+1)
        assert lo[1]<v[0] and v[1]<hi[0]
        vals.append(dict(label151=m,cosine=P.save(v),gap=[gap,gap+1]))
    lower=value(P.mul(P.pt(2),E.cosine(10)))
    upper=value(P.pt(2))
    assert E.cosine(4)[1]<lower[0]
    assert upper[1]<E.cosine(5)[0]
    paths=[f'scripts/{TAG}.py',END,P.CUSP,f'scripts/{P.TAG}.py',
           f'scripts/{E.TAG}.py',f'results/{E.TAG}.json']
    return dict(status='PASS_EXACT_ALL_LABEL_SINGLETON151_FIXED_PAIR_GAP',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        singleton_edge=[159,160],fixed_edge=[160,161],fixed_profile=37,
        fixed_contacts=row['contacts'],both_ambient158_labels=2,
        singleton_contact_parameter='d=2cos(pi/m), integer m>=3',
        formula=dict(alpha=str(alpha),beta=str(beta),radicand='8*d^2-2'),
        derivative_at_2=P.save(deriv),beta_positive=P.save(b),
        low_labels=vals,tail=dict(integer_labels='m>=10',
            lower_value=P.save(lower),limit_value=P.save(upper),gap=[4,5]),
        scope='One fixed cross-edge pair only: the singleton151 endpoint at edge159/160 versus fixed profile37 at edge160/161, both right to158, with the coherent future norm branches. Every m>=3 is obtuse or strictly between consecutive Coxeter angles. This does not complete the singleton151 family or prove a new general closure obstruction.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
