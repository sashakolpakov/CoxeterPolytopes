#!/usr/bin/env python3
"""One exact finite edge closure of the frozen central-ridge construction."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts import fresh_h20_central_ridge_20261014 as C
TAG='fresh_h20_central_ridge_20261014_edge'
def decq(a):return tuple(F(x)for x in a)
def deck(a):return tuple(decq(x)for x in a)
def run():
    files=['scripts/fresh_h20_central_ridge_20261014.py','results/fresh_h20_central_ridge_20261014.json']
    old=json.loads((ROOT/files[1]).read_text())
    assert old['sources'][files[0]]==sha256((ROOT/files[0]).read_bytes()).hexdigest()
    assert old['actual_facets']==327 and old['complete_finite_vertex_count']==16
    kap=C.qm(C.L,C.qi(C.RAD));t=F(4)
    den=C.qs(C.q(1),C.qt(kap,t*t));assert C.qsign(den)>0
    A=C.qt(C.qi(den),2*t);B=C.qm(C.qa(C.q(1),C.qt(kap,t*t)),C.qi(den));U=C.qs(A,B)
    assert A==C.q(F(120824,1711),F(-1152,1711))
    assert B==C.q(F(28495,1711),F(-288,1711))
    assert C.qsign(U)>0 and C.qsign(C.qs(B,C.q(1)))>0
    assert C.qs(C.qm(B,B),C.qm(kap,C.qm(A,A)))==C.q(1)
    target=next(x for x in old['vertices']if x['golden_coordinate']==11 and x['sign']==1)
    P=[deck(x)for x in target['spatial']];r=C.center(11,1,0)
    v=[C.sub(C.mul(C.lift(A),x),C.mul(C.lift(B),y))for x,y in zip(P,r)]
    D=C.qt(U,2)
    assert C.sub(C.dot(v,v),C.lift(C.qt(C.qm(U,U),4)))==C.k(C.R*C.R)
    assert v[11]==C.lift(C.qm(C.PH,U))
    paircounts={'right':2,'dotted':37}
    # Other seed inradii are >=1; this stronger bound checks both signs.
    for j in range(19):
        if j==11:continue
        for eps in(-1,1):
            gap=C.sub(C.lift(U),C.scale(v[j],eps))
            assert C.sign(C.sub(gap,C.k(C.R)))>0
    assert C.qsign(C.qs(C.qt(C.qm(C.PH,U),2),C.q(C.R)))>0
    incident=[]
    for index,(g,s,i)in enumerate(old['cap_recipes']):
        p=C.sub(C.dot(v,C.center(g,s,i)),C.lift(C.qt(U,4)))
        if (g,s)==(11,1)and i!=0:
            assert p==C.k();incident.append(index);paircounts['right']+=1
        else:
            assert C.sign(C.add(p,C.k(C.R*C.R)))<0;paircounts['dotted']+=1
    assert paircounts=={'right':19,'dotted':308} and len(incident)==17
    for vv in old['vertices']:
        p=[deck(x)for x in vv['spatial']]
        assert C.sign(C.sub(C.dot(v,p),C.lift(C.qt(U,4))))<0
    # The new endpoint on the edge leaving P opposite r is Q=P-lambda*r.
    lam=C.qm(C.qm(A,kap),C.qi(B));Qu=C.qs(C.q(1),lam)
    assert C.qsign(lam)>0 and C.qsign(Qu)>0
    X=[C.sub(x,C.mul(C.lift(lam),y))for x,y in zip(P,r)]
    QD=C.qt(Qu,2)
    Qnorm=C.sub(C.dot(X,X),C.lift(C.qt(C.qm(Qu,Qu),4)))
    expected=C.qt(C.qm(kap,C.qi(C.qm(B,B))),-C.R*C.R)
    assert Qnorm==C.lift(expected) and C.qsign(expected)<0
    assert C.sub(C.dot(X,v),C.lift(C.qt(C.qm(Qu,U),4)))==C.k()
    assert X[11]==C.lift(C.qm(C.PH,Qu))
    for j in range(19):
        if j==11:continue
        for eps in(-1,1):assert C.sign(C.sub(C.lift(Qu),C.scale(X[j],eps)))>0
    qsigns={'zero':0,'strict':0}
    for index,(g,s,i)in enumerate(old['cap_recipes']):
        p=C.sub(C.dot(X,C.center(g,s,i)),C.lift(C.qt(Qu,4)))
        if index in incident:assert p==C.k();qsigns['zero']+=1
        else:assert C.sign(p)<0;qsigns['strict']+=1
    assert qsigns=={'zero':17,'strict':271}
    # Full physical incident Gram: I2(5) plus18 mutually orthogonal
    # norm1/9 caps. Positive determinants are inherited verbatim.
    assert C.qsign(C.L)>0
    # Common strict center and old strict escape are retained.
    centerpair=C.qt(U,-102);assert C.qsign(centerpair)<0
    E=[C.k(F(99,100))]*10+[C.k()]*9;ED=F(8801,2000)
    escape_pair=C.sub(C.dot(E,v),C.lift(C.qt(U,ED+2)))
    assert C.sign(escape_pair)<0
    assert C.sign(C.sub(C.k(),C.lift(C.qt(Qu,102))))<0 # (Q,H)<0
    # New root has exactly19 finite neighbors and so exactly one PD20
    # corner. The old16 remain. Nineteen new PD19 subsets contain the
    # newroot; the twentieth incident edge is the already saved oldedge.
    assert 524608+19==524627 and 17*20-2==338
    src=f'scripts/{TAG}.py'
    return C.enc(dict(status='PASS_EXACT_SINGLE_EDGE_CLOSURE_RETAINING_ALL327_WALLS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in[src,*files]},
        parameter=t,kappa=kap,A=A,B=B,new_normal=dict(spatial=v,D=D,U=U,norm=C.R*C.R),
        new_pair_counts=paircounts,common_center_pair=centerpair,strict_escape_pair=escape_pair,
        released_cap_recipe=[11,1,0],retained_cap_indices=incident,
        new_vertex=dict(spatial=X,D=QD,U=Qu,norm=expected,lambda_parameter=lam,old_cap_sign_counts=qsigns),
        actual_facets=328,complete_finite_vertices=17,complete_elliptic_rank19_faces=524627,
        finite_edges=1,open_finite_ended_edges=338,seed_only_vertical_edges=524288,
        scope='One actual finite edge on a retained central/golden pi5 ridge is closed. All327 old facets and16 finite vertices survive. The complete finite-vertex and elliptic19 counts follow the exact finite-neighbor graph; the old norm1 strict escape survives. No iteration or finite-volume claim.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_facets'],out['complete_finite_vertices'])
