#!/usr/bin/env python3
"""Independent sparse-radical check of the literal data in the height theorem."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h18_general_height_lift_20261014'
TAG = BASE + '_review'

def clean(a): return {r:c for r,c in a.items() if c}
def add(*xs):
    z = {}
    for x in xs:
        for r,c in x.items(): z[r] = z.get(r,F(0))+c
    return clean(z)
def scale(a,c): return clean({r:v*c for r,v in a.items()})
def mul(a,b):
    z = {}
    for r,c in a.items():
        for s,d in b.items():
            g = gcd(r,s); k = r*s//(g*g)
            z[k] = z.get(k,F(0))+c*d*g
    return clean(z)
def read(a): return clean({1:F(a[0]),2:F(a[1])})
def sign(a):
    if not a: return 0
    for bits in (32,64,128,256):
        den = 1 << bits; lo = hi = F(0)
        for r,c in a.items():
            v=isqrt(r*den*den); l=F(v,den); u=l if v*v==r*den*den else F(v+1,den)
            lo += c*(l if c>0 else u); hi += c*(u if c>0 else l)
        if lo>0:return 1
        if hi<0:return -1
    raise AssertionError('Unresolved exact sign')
def enc(a): return {str(k):str(v) for k,v in sorted(a.items())}
def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    source='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
    primary=f'results/{BASE}.json'
    d=json.loads((ROOT/source).read_text()); cert=json.loads((ROOT/primary).read_text())
    for p,h in cert['sources'].items(): assert digest(p)==h
    G=[[read([v['a'],v['b']]) for v in row] for row in d['gram']]
    one={1:F(1)}; two={1:F(2)}
    short=[i for i in range(40) if G[i][i]==one]
    long=[i for i in range(40) if G[i][i]==two]
    assert short==cert['short_indices'] and long==cert['fixed_long_indices']
    h=[add(*row) for row in G]; h2=add(*h)
    assert h==[read(x) for x in cert['strict_center_pairings']]
    assert h2==read(cert['strict_center_norm']) and sign(h2)<0
    assert all(sign(x)<0 for x in h)
    core=[d['labels'].index(s) for s in ['L18','L25','L26']]
    assert core==cert['core_indices']
    for record in cert['core_anchor_data']:
        i=record['short']; anchors=[j for j in long if G[i][j]=={1:F(-1)}]
        assert anchors==record['fixed_long_anchors'] and anchors
        assert record['anchor_names']==[d['labels'][j] for j in anchors]
    pairs=[]
    for record in cert['null_pairs']:
        i,j=record['pair']; assert i in core and j in core and i<j
        assert G[i][j]=={1:F(-1)}
        norms=add(G[i][i],G[j][j],scale(G[i][j],2)); assert not norms
        row=[add(G[i][a],G[j][a]) for a in range(40)]
        assert row==[read(x) for x in record['old_wall_pairings']]
        assert not row[i] and not row[j]
        assert all(sign(row[a])<0 for a in short if a not in (i,j))
        assert all(sign(row[a])<=0 for a in long)
        dh=add(h[i],h[j]); tau=read(record['tau'])
        assert dh==read(record['center_pairing']) and sign(dh)<0 and sign(tau)>0
        assert scale(mul(tau,dh),4)=={1:F(-1)}
        leading=add(one,scale(mul(tau,dh),2)); quartic=mul(mul(tau,tau),h2)
        assert leading=={1:F(1,2)} and sign(quartic)<0
        assert quartic==read(record['escape_square_quartic_coefficient'])
        bound=read(record['sufficient_squared_epsilon_bound'])
        assert mul(bound,quartic)==scale(leading,-1) and sign(bound)>0
        pairs.append({'pair':[i,j],'zero_long_walls':[a for a in long if not row[a]],
                      'strict_other_short_count':sum(sign(row[a])<0 for a in short),
                      'delta_h':enc(dh),'escape_quadratic':enc(leading)})
    assert len(pairs)==3 and {tuple(x['pair']) for x in pairs}=={(17,23),(17,24),(23,24)}
    # This is the exact angular equation, separate from its all-label proof.
    for label,cos2,height2 in [(3,F(1,4),1),(4,F(1,2),0)]:
        assert 2*(1+height2)*cos2==1
        assert read(cert['core_height_square_by_anchor_label'][str(label)])==({1:F(height2)} if height2 else {})
    paths=[source,primary,f'scripts/{BASE}.py',f'scripts/{TAG}.py']
    return {'status':'PASS_INDEPENDENT_LITERAL_DATA_AND_PROOF_REVIEW',
            'sources':{p:digest(p) for p in paths},'arithmetic':'No producer/project imports; sparse squarefree-radical Fraction arithmetic and rational isqrt sign intervals.',
            'checked_center_pairings':40,'checked_null_pairings':120,'anchor_counts':[2,2,2],
            'null_pairs':pairs,
            'scope':'Independent literal anchor/center/null/escape-coefficient replay plus source/proof review; no second source inertia or full original 40-wall geometric census.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=run(); path=ROOT/f'results/{TAG}.json'
    if a.verify: assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
