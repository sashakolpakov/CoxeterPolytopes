#!/usr/bin/env python3
"""Complete one-facet tetrahedron gate on the actual H4+I2(9) face."""
import argparse
import itertools
import json
from hashlib import sha256
from pathlib import Path
import sympy as sy
import fresh_nonarithmetic_pi8_quad_completion_20261017 as C
P,U,I,D=C.P,C.U,C.I,C.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_pi9_tetrahedron_gate_20261017'
DOMAIN='fresh_nonarithmetic_pi9_endpoint_domain_20261017'
SUPPORT=[0,1,2,3,58,59]
SIDES=[5,6,55]


def positive(x):
    x=sy.simplify(x)
    answer=sy.ask(sy.Q.positive(x))
    assert answer is not None,x
    return bool(answer)


def run():
    z=C.geometry();normals,pair=z['normals'],z['pair'];domain=C.read(DOMAIN)
    frame=[normals[i]for i in SUPPORT];gram=[[pair(a,b)for b in frame]for a in frame]
    support_pivots=P.pd(gram)
    for i in SIDES:
        for j in SUPPORT:
            q=pair(normals[i],normals[j]);assert q.l<=0<=q.h
    # Existing simple finite corner in this H3 face.
    old=C.read(C.TAG)
    record=next(r for r in old['complete_six_pair_intersections']if set(r['pair'])=={58,59})
    assert set(record['incident'])==set(SUPPORT+SIDES)
    point=C.decode(record['point']);assert pair(point,point).h<0 and pair(point,z['h']).h<0
    theta={2:sy.Integer(0),3:sy.Integer(1),4:sy.sqrt(2),5:(1+sy.sqrt(5))/2}
    census=[];accepted=[]
    for labels in itertools.product(range(2,6),repeat=3):
        g=sy.Matrix([[2,-1,0,-theta[labels[0]]],[-1,2,0,-theta[labels[1]]],
                     [0,0,2,-theta[labels[2]]],[-theta[labels[0]],-theta[labels[1]],-theta[labels[2]],2]])
        minors=[sy.simplify(g.extract(ids,ids).det())for ids in itertools.combinations(range(4),3)]
        determinant=sy.simplify(g.det())
        ok=all(positive(x)for x in minors)and positive(-determinant)
        census.append(dict(labels=list(labels),proper_three_minors=[str(x)for x in minors],determinant=str(determinant),compact_tetrahedron=ok))
        if ok:accepted.append(list(labels))
    assert accepted==[[2,5,3],[5,2,3]]
    branches=[]
    a,b=C.decode(domain['A']),C.decode(domain['B'])
    for profile_index in(3,6):
        pf=domain['profiles'][profile_index];assert pf['label']==5
        p=C.decode(pf['projection']);s=I(*pf['scale'])
        # r.55=-s*(t+1/t)=-1; the only endpoint branch t>1.
        disc=U.rnd(1-4*s*s);assert disc.l>0
        t=U.rnd((1+U.sqrt(disc))/(2*s));assert t.l>1
        r=[U.rnd(x+s*(t*y+w/t))for x,y,w in zip(p,a,b)]
        square=U.rnd(pair(r,r));assert square.l<=2<=square.h
        products=[U.rnd(pair(r,n))for n in normals]
        prescribed={i:2 for i in domain['edge_support']+[55]}
        prescribed[pf['nonright_A2_node']]=5;prescribed[55]=3
        for i,m in prescribed.items():
            residual=products[i]+D.theta(m);assert residual.l<=0<=residual.h
        blocker=products[49];assert blocker.l>-1 and blocker.h<0
        # The advertised sole failure among the other native contacts.
        for i,q in enumerate(products):
            if i not in prescribed and i!=49:assert q.h<-2,(profile_index,i,q.enc())
        branches.append(dict(profile_index=profile_index,node=pf['nonright_A2_node'],side_labels=[prescribed[i]for i in SIDES],
                             parameter=t.enc(),discriminant=disc.enc(),normal=P.enc(r),norm=square.enc(),
                             old60_products=P.enc(products),blocker_wall=49,blocker_product=blocker.enc(),consecutive_labels=[2,3]))
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',f'scripts/{DOMAIN}.py',f'results/{DOMAIN}.json',
           f'scripts/{D.TAG}.py',f'results/{D.TAG}.json']
    return dict(status='PASS_NO_ONE_FACET_COMPACT_TETRAHEDRON_ON_PI9_FACE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                retained_walls=60,actual_H3_support=SUPPORT,support_PD_pivots=support_pivots,retained_side_normals=SIDES,
                retained_finite_corner=P.enc(point),profile_count=64,compact_tetrahedron_templates=accepted,profile_census=census,
                two_exact_candidates=branches,
                all_label_completeness='A compact tetrahedron has all proper principal Gram blocks positive definite and full Gram Lorentzian. A finite label>=6 in a connected4-node minimal indefinite diagram cannot have any adjoining nonright edge: its3-node principal path would be nonpositive. Thus every finite label is2..5. The64 literal contact triples to the retained A2+A1 corner exhaust the possibilities. Only the two marked [3,5,3] chains survive.',
                actual_face_and_facet_scope='The six retained support roots form PD H4+I2(9), and the saved old finite corner contains exactly these six plus5,6,55. Therefore the actual H3 face is nonempty and those three old sides are genuine incident facets. Any new actual facet in a compact completion extends this support to a finite diagram; H4 has no nonright finite extension, and I2(9) inverse diagonal>2 excludes every nonright contact. Hence its intrinsic normal equals its norm2 native normal and it is right to the whole support. The tetrahedron must retain these three incident sides; finite-vertex persistence follows from a common strict interior and acuteness.',
                orientation='The complete pi9 endpoint domain gives both Lorentz2 orientations. Old55 excludes the opposite orientation, while actual endpoint feasibility and acute common interior force t>1 on the surviving branch. Each [3,5,3] candidate has exactly the displayed root of s*(t+1/t)=1. The other positive quadratic root is below1 and is excluded by that endpoint condition.',
                scope='No single added coherent Coxeter wall can close this particular actual H3 face as a compact tetrahedron. Multiwall prism or other compact-polyhedron closures remain open; this is not a noncompletion theorem for the60-wall chamber.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text());print('PASS saved pi9 tetrahedron gate replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
