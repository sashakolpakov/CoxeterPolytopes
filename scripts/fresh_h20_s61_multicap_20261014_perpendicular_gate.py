#!/usr/bin/env python3
"""Arbitrary acute perpendicular-S61 additions preserve the saved escape.

Small exact consequence of the separately frozen compact29-gon certificate.
Arithmetic helpers are reused; this is not an independent arithmetic audit.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import fresh_h20_s61_multicap_20261014_exact as C
from fresh_h20_coupled_roles_20261014_field import Field,parse,enc,pd,B
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_s61_multicap_20261014_perpendicular_gate'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json',
       'results/fresh_h20_e7_bridge_face_20261014_exact.json',
       'results/fresh_h20_bridge_edge_closure_20261014.json',
       'results/fresh_h20_s61_multicap_20261014_exact.json']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old,bridge,edge,full=[json.loads((ROOT/p).read_text())for p in FILES]
    for data in[old,bridge,edge,full]:
        for p,h in data['sources'].items():assert sha(p)==h,p
    A=Field(old,bridge,edge);I=C.I
    fr=full['frame'];X=tuple(parse(x,2)for x in fr['X']);r=tuple(parse(x,2)for x in fr['r95']);W=tuple(parse(x,2)for x in fr['W'])
    AA=parse(fr['A'],2);ww=parse(fr['W_norm'],2);S=full['support']
    assert S==edge['recipe']['support']and edge['recipe']['contact_wall']==61
    assert A.dot(X,X,2)==A.scale(AA,-1,2)and A.dot(r,r,2)==A.c(2,2)and A.dot(W,W,2)==ww
    assert all(A.dot(v,w,2)==A.z[2]for v,w in[(X,r),(X,W),(r,W)])
    E=tuple(A.c(x,2)for x in full['strict_escape']['rational_point'])
    assert A.dot(E,E,2)==A.c(F(38017800175279,500000000000),2)
    def bounds(x):return I(*A.bounds(x,2,C.BITS))
    ep=[bounds(A.dot(v,E,2))/bounds(q).sqrt()for v,q in zip([X,r,W],[AA,A.c(2,2),ww])]
    # Pair covectors have a negative time component; these are vector coordinates.
    proj=[-ep[0],ep[1],ep[2]];qproj=C.fdot(proj,proj)
    assert proj[0].l>0 and qproj.h<0
    a=C.D(I(*map(F,full['parameter']['isolating_interval'])),1)
    chain,_,_=C.chain(a);co=[[x.v for x in v]for v in chain]
    root2=I(2).sqrt();rootw=bounds(ww).sqrt();k=(bounds(ww)-F(5,6)).sqrt()
    face=[[I(0),root2,I(0)],[k,-1/root2,rootw]]+co+[[k,-1/root2,-rootw]]
    pairs=[C.fdot(proj,n)for n in face];assert all(x.h<0 for x in pairs)
    # Verify the inverse-positivity premise on every actual vertex frame.
    G=[[F(x)for x in row]for row in full['support_gram']];assert pd(G)
    assert all(x>=0 for row in B.inv(G)for x in row)
    grams=[]
    for rec in full['face']['vertices']:
        u,v=rec['walls'];H=[row+[F(0)]*2 for row in G]+[[F(0)]*20 for _ in range(2)]
        for j,wall in enumerate([u,v],18):
            H[j][j]=2
            if wall in[96,97]:H[j][S.index(61)]=H[S.index(61)][j]=F(-1)
        H[18][19]=H[19][18]=F(-1)if 95 in[u,v]else F(0)
        assert pd(H);J=B.inv(H);assert all(x>=0 for row in J for x in row)
        grams.append(dict(walls=[u,v],positive_definite=True,inverse_entrywise_nonnegative=True,minimum_inverse_entry=min(x for row in J for x in row)))
    assert len(grams)==29
    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h20_s61_multicap_20261014_exact.py','scripts/fresh_h20_coupled_roles_20261014_field.py']
    return enc(dict(status='PASS_ARBITRARY_ACUTE_PERPENDICULAR_S61_ESCAPE_GATE',sources={p:sha(p)for p in paths},interval_bits=C.BITS,support=S,
        original_escape=full['strict_escape']['rational_point'],original_escape_norm=F(38017800175279,500000000000),projected_escape_frame_bounds=C.out(proj),projected_escape_norm_bound=qproj.out(),projected_escape_polygon_pair_bounds=C.out(pairs),actual_vertex_Gram_checks=grams,
        vertex_preservation='For an actual finite vertex v and any common future timelike strict interior h, write h=alpha*v+sum_i lambda_i s_i. Here alpha=(h,v)/(v,v)>0 and lambda=G^-1*((h,s_i)) has nonpositive entries. For an acute added normal n, (n,s_i)<=0, hence sum lambda_i*(n,s_i)>=0. Since(n,h)<0, necessarily(n,v)<0.',
        convex_hull_step='The certified bounded Klein29-gon is the convex hull of its29 finite vertices. The projected escape is future timelike and strictly inside every polygon halfspace, so its positive-time normalization is a convex combination of the vertices. Every acute added wall is strict at every vertex and therefore at the projected escape.',
        perpendicular_step='If n is perpendicular to S61, (n,E)=(n,projection_S61(E))<0. Thus the original strict spacelike escape survives every finite coherent family of such additions.',
        scope='All real spacelike added normals perpendicular to S61 and acute with all retained133 walls, provided a common future timelike strict interior exists. No norm, field, Coxeter-label, radius, or parameter cutoff is imposed. Any finite such completion remains infinite-volume. Nonright S61 contacts and changed retained cores are not excluded.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['actual_vertex_Gram_checks']),flush=True)
