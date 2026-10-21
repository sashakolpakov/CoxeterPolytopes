#!/usr/bin/env python3
"""Exact escape-cutting walls found on complete vertex-free geodesic faces."""
import argparse,json
from itertools import combinations
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_vertex_free_caps_20261015'
BASE='results/fresh_h9_quartet_polar_second_20261014.json'
FACES='results/fresh_h9_vertex_free_edges_20261015.json'
RECIPES=[(31,{0:5,1:2,2:2,3:2,4:2,5:5,6:2,7:2,8:3,15:5,16:2,17:3,18:2}),
         (13,{0:3,1:2,2:2,3:5,4:2,5:2,6:2,7:2,8:5,15:2,16:5,17:2,18:3})]

def run():
    d,f=B.read(BASE),B.read(FACES);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
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
        face=f['faces'][face_id];S=face['support']
        V=C._full_rank_kernel([cov[i]for i in S]+[A.mv(Q,n)]);assert V is not None
        if A.pair(V,H,Q).sign()>0:V=[-v for v in V]
        assert A.pair(V,V,Q).sign()<0
        vp=[A.pair(V,r,Q)for r in N];assert all(z.sign()==(0 if i in S else -1)for i,z in enumerate(vp))
        assert not A.pair(V,n,Q)
        rows.append(dict(source_face=face_id,reconstruction_basis=basis,finite_contacts=[[i,m]for i,m in contacts.items()],
                         normal=n,norm=A.L(2),pairings=pairings,relations=relations,
                         center_pairing=hp,old_escape_pairing=ep,
                         new_finite_endpoint=dict(point=V,norm=A.pair(V,V,Q),old_wall_pairings=vp,support=S)))
    g=A.pair(rows[0]['normal'],rows[1]['normal'],Q)
    assert g.sign()<=0;kind,m,c=R.relation(g,A.L(2),A.L(2));assert kind!='forbidden'
    paths=[BASE,FACES,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py']
    return A.enc(dict(status='PASS_EXACT_TWO_VERTEX_FREE_ESCAPE_CUTTERS',sources={p:R.sha(p)for p in paths},
                     field=d['field'],caps=rows,mutual_pair=dict(pairing=g,kind=kind,label=m,invariant=c),
                     scope='Two exact norm2 walls compatible with all25 oldwalls, with one another, retaining the common center and cutting the previously persistent strict escape. Each creates a finite endpoint on a formerly vertex-free geodesic face. Global geometry and a replacement escape are separate tests; finite volume is not asserted.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['mutual_pair']);print([r['old_escape_pairing']for r in d['caps']])
