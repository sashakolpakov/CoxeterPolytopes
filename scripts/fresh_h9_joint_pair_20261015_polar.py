#!/usr/bin/env python3
"""Exact polar cut of the current30wall escape, with no norm2 rescaling."""
import argparse,json
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_joint_pair_20261015_polar'
BASE='results/fresh_h9_nine_contact_extension_20261015.json'
def run():
    d=B.read(BASE);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    n=D.dv(d['strict_escape']['hyperideal_ray']);norm=A.pair(n,n,Q);assert norm.sign()>0
    vals=[A.pair(n,m,Q)for m in N];norms=[A.pair(m,m,Q)for m in N]
    rel=[]
    for i,g in enumerate(vals):
        kind,m,c=R.relation(g,norm,norms[i]);assert kind!='forbidden',(i,kind,g)
        rel.append(dict(wall=i,kind=kind,label=m,invariant=c))
    assert [i for i,v in enumerate(vals)if not v]==d['strict_escape']['active']
    assert A.pair(n,H,Q).sign()<0 and A.pair(n,E,Q).sign()>0
    # A strict-center projection proves this wall is an actual facet after addition.
    t=A.pair(H,n,Q)/norm;V=[x-t*y for x,y in zip(H,n)];vp=[A.pair(V,m,Q)for m in N]
    assert A.pair(V,V,Q).sign()<0 and not A.pair(V,n,Q) and all(x.sign()<0 for x in vp)
    files=[BASE,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py',
           'scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py']
    return A.enc(dict(status='PASS_EXACT_H9_POLAR_ESCAPE_CAP',sources={p:R.sha(p)for p in files},field=d['field'],
       caps=[dict(normal=n,norm=norm,relations=rel,pairings=vals,old_escape_pairing=A.pair(n,E,Q),center_pairing=A.pair(n,H,Q),
                  facet_witness=dict(point=V,norm=A.pair(V,V,Q),pairings=vp))],
       scope='Exact polar of a hyperideal ray, compatible with all30old walls and cutting the saved strictescape. Actualfacet witness certified. No finitevolume or newvertex claim.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
