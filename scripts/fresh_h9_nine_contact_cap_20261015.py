#!/usr/bin/env python3
"""A nine-contact quadratic solution returns to the existing quartic field."""
import argparse,json
from fractions import Fraction as F
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_nine_contact_cap_20261015'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
OCTIC='results/fresh_h9_octic_escape_caps_20261015.json'
SUPPORT=[0,1,2,5,6,8,10,28]

def run():
    d,o=B.read(BASE),B.read(OCTIC)
    A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    C.Golden=A.L;C.ZERO=A.Z
    cov=[A.mv(Q,n)for n in N];ids=SUPPORT+[3]
    assert all(A.pair(N[i],N[i],Q)==2 for i in ids)
    v=C._full_rank_kernel([cov[i]for i in ids]);assert v is not None and v[9]
    v=[a/v[9]for a in v]
    rhs=[A.L(-1 if i in(3,8)else 0)for i in ids]
    p=A.mv(A.inv([cov[i]for i in ids]+[A.mv(Q,v)]),rhs+[A.Z])
    assert not A.pair(p,v,Q)
    delta=(A.L(2)-A.pair(p,p,Q))/A.pair(v,v,Q)
    assert delta==A.L(A.A.Q(F(2740,121),F(1284,121)))
    beta=A.A.Q(F(5,11),F(7,11));root=A.L(0,beta)
    assert root*root==delta and root.sign()>0
    n=[a-root*b for a,b in zip(p,v)];assert A.pair(n,n,Q)==2
    vals=[A.pair(n,r,Q)for r in N];norms=[A.pair(r,r,Q)for r in N]
    assert all(vals[i]==x for i,x in zip(ids,rhs))
    relations=[]
    for i,g in enumerate(vals):
        kind,m,c=R.relation(g,A.L(2),norms[i]);assert kind!='forbidden'
        relations.append(dict(wall=i,kind=kind,label=m,invariant=c))
    assert A.pair(n,H,Q).sign()<0 and A.pair(n,E,Q).sign()>0
    common=D.dv(o['common_strict_escape']['point']);ep=A.pair(n,common,Q);assert ep.sign()>0
    V=C._full_rank_kernel([cov[i]for i in SUPPORT]+[A.mv(Q,n)]);assert V is not None
    if A.pair(V,H,Q).sign()>0:V=[-a for a in V]
    vp=[A.pair(V,r,Q)for r in N]
    assert A.pair(V,V,Q).sign()<0 and not A.pair(V,n,Q)
    assert all(a.sign()==(0 if i in SUPPORT else -1)for i,a in enumerate(vp))
    paths=[BASE,OCTIC,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py',
           'scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py',
           'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py','src/coxeter_polytopes/verify/chamber.py']
    return A.enc(dict(status='PASS_EXACT_NINE_CONTACT_QUARTIC_CAP',sources={q:R.sha(q)for q in paths},field=d['field'],
        caps=[dict(normal=n,norm=A.L(2),finite_contacts=[[i,3 if i in(3,8)else 2]for i in ids],
            reconstruction=dict(kernel=v,orthogonal_particular=p,delta=delta,positive_square_root=root),
            pairings=vals,relations=relations,center_pairing=A.pair(n,H,Q),old_escape_pairing=A.pair(n,E,Q),
            octic_common_escape_pairing=ep,finite_endpoint=dict(support=SUPPORT,point=V,norm=A.pair(V,V,Q),pairings=vp))],
        scope='One norm2 quartic wall from nine independent contacts and an exact quadratic solution. Compatible with all29 retainedwalls, cuts both recordedescapes and creates a finite endpoint on a new29-wall germ. Globalgeometry and volume require separatechecks.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
