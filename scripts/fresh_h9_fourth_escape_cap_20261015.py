#!/usr/bin/env python3
"""Exact wall cutting the replacement escape of the28-wall chamber."""
import argparse,json
from itertools import combinations
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_fourth_escape_cap_20261015'
BASE='results/fresh_h9_replacement_escape_extension_20261015.json'
FACES='results/fresh_h9_after_third_vertex_free_edges_20261015.json'
RECIPES=[(4,{0:2,1:2,2:2,3:5,5:3,6:2,7:5,8:3,10:2,27:2})]

def run():
    d=B.read(BASE);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    norms=[A.pair(n,n,Q)for n in N];cov=[A.mv(Q,n)for n in N]
    phi=A.L(A.A.Q(A.F(1,2),A.F(1,2)));vals={2:A.Z,3:-A.ONE,5:-phi}
    rows=[]
    for face_id,contacts in RECIPES:
        basis=[]
        for i in contacts:
            assert norms[i]==2
            if C._row_rank([cov[j]for j in basis+[i]],10)>len(basis):basis.append(i)
            if len(basis)==10:break
        assert len(basis)==10
        n=A.mv(A.inv([cov[j]for j in basis]),[vals[contacts[j]]for j in basis])
        assert A.pair(n,n,Q)==2
        pairings=[A.pair(n,r,Q)for r in N]
        assert all(pairings[i]==vals[m]for i,m in contacts.items())
        relations=[]
        for i,g in enumerate(pairings):
            kind,m,c=R.relation(g,A.L(2),norms[i]);assert kind!='forbidden'
            relations.append(dict(wall=i,kind=kind,label=m,invariant=c))
        hp=A.pair(n,H,Q);ep=A.pair(n,E,Q);assert hp.sign()<0 and ep.sign()>0
        S=[0,1,2,5,6,8,10,27]
        V=C._full_rank_kernel([cov[i]for i in S]+[A.mv(Q,n)]);assert V is not None
        if A.pair(V,H,Q).sign()>0:V=[-v for v in V]
        assert A.pair(V,V,Q).sign()<0
        vp=[A.pair(V,r,Q)for r in N];assert all(z.sign()==(0 if i in S else -1)for i,z in enumerate(vp))
        assert not A.pair(V,n,Q)
        rows.append(dict(source_face=face_id,reconstruction_basis=basis,finite_contacts=[[i,m]for i,m in contacts.items()],
                         normal=n,norm=A.L(2),pairings=pairings,relations=relations,
                         center_pairing=hp,old_escape_pairing=ep,
                         new_finite_endpoint=dict(point=V,norm=A.pair(V,V,Q),old_wall_pairings=vp,support=S)))
    paths=[BASE,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py']
    return A.enc(dict(status='PASS_EXACT_FOURTH_ESCAPE_CUTTER',sources={p:R.sha(p)for p in paths},
                     field=d['field'],caps=rows,
                     scope='One exact norm2 wall compatible with all28 existingwalls, retaining the common center and cutting the replacement strict escape. It creates a finite endpoint on a vertex-free geodesic face of the28-wall chamber. Global geometry and a replacement escape are separate tests; finite volume is not asserted.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status']);print([r['old_escape_pairing']for r in d['caps']])
