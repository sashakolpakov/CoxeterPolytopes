#!/usr/bin/env python3
"""Separate quartic-engine review of the dotted cutter and facet diagnosis."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from fresh_h9_mixed_d4_coupling_20261014 import L,K,Z,pair,enc
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_global_caps_20261014_review'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
PRIMARY='results/fresh_h9_mixed_d4_global_caps_20261014.json'
def dec(x):return L(K(*x[0]),K(*x[1]))
def dv(v):return [dec(x)for x in v]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old,pr=[json.loads((ROOT/p).read_text())for p in (BASE,PRIMARY)]
    for doc in (old,pr):
        for p,h in doc['sources'].items():assert digest(p)==h
    Q=[dv(v)for v in old['form']];N=[dv(v)for v in old['final_normals']];H=dv(old['center'])
    E=[L(F(s))for s in pr['literal_escape']];t=F(pr['interpolation_parameter'])
    assert t==F(65323,65536) and len(N)==12
    v=[(1-t)*h+t*e for h,e in zip(H,E)];assert v==dv(pr['raw_cap'])
    q=pair(v,v,Q);qE=pair(E,E,Q);qH=pair(H,H,Q);vh=pair(v,H,Q)
    assert q.sign()>0 and qE.sign()>0 and qH.sign()<0 and vh.sign()<0
    assert pair(v,E,Q).sign()>0 and pair(E,H,Q).sign()<0
    assert all(pair(E,n,Q).sign()<0 and pair(H,n,Q).sign()<0 for n in N)
    assert q==dec(pr['raw_cap_norm']) and qE==dec(pr['escape_norm'])
    p=[h-vh*x/q for h,x in zip(H,v)]
    assert p==dv(pr['new_facet_point']) and pair(p,v,Q)==Z and pair(p,p,Q).sign()<0
    margins=[];projections=[];actual=[]
    for i,n in enumerate(N):
        nn=pair(n,n,Q);vn=pair(v,n,Q);assert nn.sign()>0 and vn.sign()<0
        margin=vn*vn-q*nn;assert margin.sign()>0
        assert margin==dec(pr['strict_dotted_squared_margins'][i]);margins.append(margin)
        np=[a-vn*b/q for a,b in zip(n,v)]
        assert pair(np,v,Q)==Z and pair(np,np,Q)==nn-vn*vn/q
        assert pair(np,np,Q).sign()<0 and pair(np,p,Q).sign()<0
        assert pair(np,np,Q)==dec(pr['old_projected_norms'][i]);projections.append(pair(np,np,Q))
        hn=pair(H,n,Q);f=[h-hn*x/nn for h,x in zip(H,n)]
        assert pair(f,f,Q).sign()<0 and pair(f,n,Q)==Z and pair(f,v,Q).sign()<0
        assert all(pair(f,m,Q).sign()==(0 if i==j else -1)for j,m in enumerate(N))
        actual.append(pair(f,v,Q))
    paths=[BASE,PRIMARY,'scripts/fresh_h9_mixed_d4_coupling_20261014.py',
           'scripts/fresh_h9_mixed_d4_global_caps_20261014.py',
           'results/fresh_h9_mixed_d4_global_caps_20261014.md',f'scripts/{TAG}.py']
    return enc(dict(status='PASS_SEPARATE_QUARTIC_DOTTED_CAP_AND_COVERAGE_PROOF_REVIEW',
      sources={p:digest(p)for p in paths},old_wall_count=len(N),parameter=str(t),
      raw_norm=q,escape_norm=qE,positive_dotted_margins=margins,
      negative_projected_norms=projections,actual_old_facet_new_pairings=actual,
      new_facet='entire H8 by future timelike projected old normals',
      proof_checks=['future-null interpolation with positive escape pairing',
        'uniform collar of an entire hyperbolic facet gives infinite volume',
        'inward normal coordinate has negative r pairing',
        'old inequality gives tanh(tau)<=b*cosh(distance)/gamma',
        'strict open-ball cover of compact P gives a uniform finite height bound',
        'nonintersection of open P with all closed balls gives infinite volume',
        'parallel side contacts give isolated boundary tangencies; H8 base facets have positive dimension',
        'acuteness puts every projected old center in P, hence their hyperbolic convex hull in P'],
      scope='Separate exact scalar implementation for this literal example, plus direct mathematical review of the general interpolation and coverage arguments. Inherits the old chamber. Does not independently enumerate a new cone, establish a source placement, or certify numerical coverage optimization.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['old_wall_count'],d['parameter'],flush=True)
