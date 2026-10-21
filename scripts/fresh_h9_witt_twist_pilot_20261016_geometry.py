#!/usr/bin/env python3
"""Complete cone and exact witnesses for the diagonal local-Witt pilot."""
import argparse,json,sys
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.golden import Golden as K
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone,_pair,_full_rank_kernel,_mat_vec

TAG='fresh_h9_witt_twist_pilot_20261016_geometry'
PREFIX='results/fresh_h9_witt_twist_pilot_20261016.json'
def sha(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def vector(row):return tuple(K.from_json(x) for x in row)
def enc(x):
    if isinstance(x,K):return x.to_json()
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    return x

def run():
    data=json.loads((ROOT/PREFIX).read_text())
    for p,h in data['sources'].items():assert sha(p)==h,p
    G=tuple(vector(row) for row in data['form'])
    N=tuple(vector(row) for row in data['roots'])
    C=vector(data['strict_center']);H=vector(data['basepoint'])
    assert _pair(C,G,C)<0 and all(_pair(C,G,n)<0 for n in N)
    cone=diagnose_chamber_cone(N,G,C,dimension=9)
    names=('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')
    census={key:[list(x) for x in getattr(cone,key)] for key in names}
    finite=[]
    for support in cone.finite_vertices:
        assert len(support)==9
        point=_full_rank_kernel([_mat_vec(G,N[i]) for i in support]);assert point is not None
        if _pair(point,G,H)>0:point=tuple(-x for x in point)
        values=[_pair(point,G,n) for n in N]
        assert _pair(point,G,point)<0 and _pair(point,G,H)<0
        assert tuple(i for i,x in enumerate(values) if not x)==support
        assert all(x<=0 for x in values)
        finite.append(dict(incident=list(support),point=point,norm=_pair(point,G,point)))
    assert not cone.cusps and not cone.opposite_sheet_rays
    escape=None
    if cone.hyperideal_rays:
        active=cone.hyperideal_rays[0]
        support=next(s for s in combinations(active,9)
                     if _full_rank_kernel([_mat_vec(G,N[i]) for i in s]) is not None)
        point=_full_rank_kernel([_mat_vec(G,N[i]) for i in support])
        values=[_pair(point,G,n) for n in N]
        if any(x>0 for x in values):point=tuple(-x for x in point)
        assert _pair(point,G,point)>0 and all(_pair(point,G,n)<=0 for n in N)
        epsilon=F(1)
        while True:
            E=tuple(x+epsilon*y for x,y in zip(point,C))
            if _pair(E,G,E)>0 and _pair(E,G,H)<0:break
            epsilon/=2;assert epsilon>F(1,2**40)
        values=[_pair(E,G,n) for n in N]
        assert all(x<0 for x in values)
        escape=dict(point=E,norm=_pair(E,G,E),pairings=values,basepoint_pairing=_pair(E,G,H),
                    center_pairing=_pair(E,G,C),hyperideal_ray=point,ray_active=list(active),epsilon=epsilon)
    compact=not cone.hyperideal_rays
    paths=[PREFIX,f'scripts/{TAG}.py','src/coxeter_polytopes/exact/golden.py',
           'src/coxeter_polytopes/exact/matrix.py','src/coxeter_polytopes/verify/chamber.py']
    return enc(dict(status='PASS_COMPLETE_WITT_TWIST_PILOT_CONE',sources={p:sha(p) for p in paths},
        complete_census=census,counts={key:len(value) for key,value in census.items()},
        finite_vertices=finite,strict_escape=escape,compact_prefix=compact,
        scope='Exact full double-description census of the bounded Vinberg prefix, with physical finite-vertex and escape checks. A strict spacelike escape proves this prefix infinite-volume; it does not decide the full lattice reflection chamber.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text());print('PASS saved pilot cone replay')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],'compact',result['compact_prefix'])
