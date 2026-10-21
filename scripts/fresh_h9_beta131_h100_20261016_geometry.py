#!/usr/bin/env python3
"""Exact global cone geometry of the complete both-norm Vinberg prefix."""
import argparse,json,sys
from fractions import Fraction as F
from itertools import combinations
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.golden import Golden as K,ZERO
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone,verify_hyperideal_extreme_ray,_pair,_full_rank_kernel,_mat_vec,_dot
TAG='fresh_h9_beta131_h100_20261016_geometry'
BASE='results/fresh_h9_beta131_h100_20261016.json'
def dec(x):return K(F(x['a']),F(x['b']))
def vec(x):return tuple(map(dec,x))
def mat(x):return tuple(map(vec,x))
def enc(x):
    if isinstance(x,K):return x.to_json()
    if isinstance(x,(tuple,list)):return[enc(y)for y in x]
    if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
    if isinstance(x,F):return str(x)
    return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d=json.loads((ROOT/BASE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    G,N,H=mat(d['form']),mat(d['roots']),vec(d['basepoint'])
    assert len(N)>=12 and d['height_bound']=='100'
    P=vec(d['strict_center']);eps=F(d['center_epsilon'])
    assert _pair(P,G,P).sign()<0 and all(_pair(P,G,n).sign()<0 for n in N)
    geometry=diagnose_chamber_cone(N,G,P,dimension=9)
    census={name:[list(ids)for ids in getattr(geometry,name)]for name in('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    assert not census['cusps'] and not census['opposite_sheet_rays']
    finite=[]
    for ids in geometry.finite_vertices:
        assert len(ids)==9
        v=_full_rank_kernel([_mat_vec(G,N[i])for i in ids]);assert v is not None
        if _pair(v,G,P).sign()>0:v=tuple(-x for x in v)
        values=[_pair(v,G,n)for n in N]
        assert _pair(v,G,v).sign()<0 and _pair(v,G,P).sign()<0
        assert [i for i,x in enumerate(values)if not x]==list(ids)
        assert all(x.sign()<=0 for x in values)
        finite.append(dict(incident=list(ids),point=v,norm=_pair(v,G,v)))
    assert geometry.hyperideal_rays
    active=geometry.hyperideal_rays[0]
    support=next(ids for ids in combinations(active,9)if _full_rank_kernel([_mat_vec(G,N[i])for i in ids])is not None)
    r=_full_rank_kernel([_mat_vec(G,N[i])for i in support])
    rp=[_pair(r,G,n)for n in N]
    if any(x.sign()>0 for x in rp):r=tuple(-x for x in r);rp=[-x for x in rp]
    assert all(x.sign()<=0 for x in rp) and tuple(i for i,x in enumerate(rp)if not x)==active
    assert _pair(r,G,r).sign()>0
    # Move into all strict wall halfspaces, using exact rational steps.
    # A hyperideal ray may be on either side of P's projective slice, so
    # search a positive rational combination instead of assuming orientation.
    if _pair(r,G,P).sign()<0:
        e=F(1)
        while True:
            E=tuple(x+e*y for x,y in zip(r,P))
            if _pair(E,G,E).sign()>0 and _pair(E,G,P).sign()<0:break
            e/=2;assert e>F(1,2**40)
        mixture=dict(ray_coefficient=F(1),center_coefficient=e)
    else:
        raise AssertionError('selected hyperideal ray needs general future mixing')
    ev=[_pair(E,G,n)for n in N]
    assert all(x.sign()<0 for x in ev)
    paths=[BASE,f'scripts/{TAG}.py','src/coxeter_polytopes/exact/golden.py',
           'src/coxeter_polytopes/exact/matrix.py','src/coxeter_polytopes/verify/chamber.py']
    return enc(dict(status='PASS_COMPLETE_BETA131_H100_PREFIX_CONE_INFINITE_VOLUME',sources={p:sha(p)for p in paths},
        arithmetic='Exact Golden rational pairs; independent cone module, no Vinberg enumeration import.',
        height_bound=100,wall_count=len(N),strict_center=dict(point=P,epsilon=eps,norm=_pair(P,G,P)),
        complete_census=census,counts={k:len(v)for k,v in census.items()},finite_vertices=finite,
        hyperideal_witness=dict(active=list(active),independent_support=list(support),point=r,norm=_pair(r,G,r),pairings=rp),
        strict_escape=dict(point=E,norm=_pair(E,G,E),center_pairing=_pair(E,G,P),pairings=ev,mixture=mixture),
        scope='Exact incremental double description reconstructs every extreme ray of theheight100 prefix. A strict future spacelike escape proves this prefix has infinitevolume. Later roots and reflectivity remain unresolved; complete bounded root enumeration is inherited separately.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
