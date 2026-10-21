#!/usr/bin/env python3
"""An exact missing-S6 cap cutting the retained60/62 strict escape."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_pi9_prism_completion_20261017 as N

C,P,U,I,D=N.C,N.P,N.U,N.I,N.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_escape62_20261017_cap'
SUPPORT=[0,1,2,3,5,6]
SEED=[(7,409),(33,633),(49,556),(30,-186)]


def geometry():
    z=N.geometry();normals=z['normals'][:];pair=z['pair']
    frame=[normals[i]for i in SUPPORT]
    gram=[[pair(a,b)for b in frame]for a in frame]
    inverse,pivots=P.C53.inverse(gram)
    seed=P.combine([I(c)for _,c in SEED],[normals[i]for i,_ in SEED])
    coeff=U.imv(inverse,[pair(seed,n)for n in frame])
    y=P.sub(seed,P.combine(coeff,frame))
    yy=U.rnd(pair(y,y));assert yy.h<0
    yf=U.rnd(pair(y,z['h']));assert yf.h<0
    p=[U.rnd(-U.sqrt(I(2))*(a+b))for a,b in zip(normals[5],normals[6])]
    scale=U.sqrt(I(2)/(-yy))
    r=[U.rnd(a+scale*b)for a,b in zip(p,y)]
    normals.append(r)
    z.update(normals=normals,cap_recipe=dict(seed=SEED,seed_vector=P.enc(seed),projection_support=SUPPORT,
        projection_coefficients=P.enc(coeff),projection_inverse_pivots=pivots,projected_seed=P.enc(y),
        projected_seed_square=yy.enc(),projected_seed_future=yf.enc(),normalization=scale.enc(),A2_component=P.enc(p)))
    return z


def run():
    old=C.read(N.TAG);z=geometry();normals,pair=z['normals'],z['pair'];r=normals[-1]
    assert len(normals)==63
    center=C.decode(old['common_center']);center_products=[U.rnd(pair(center,n))for n in normals]
    assert all(q.h<0 for q in center_products)
    assert pair(center,center).h<0 and pair(center,z['h']).h<0
    square=U.rnd(pair(r,r));assert square.l<=2<=square.h
    contacts=[]
    for j,n in enumerate(normals[:-1]):
        value=U.rnd(pair(r,n))
        if j in [0,1,2,3]:
            assert value.l<=0<=value.h;label=2
        elif j in [5,6]:
            q=value+D.theta(4);assert q.l<=0<=q.h;label=4
        else:
            assert value.h<-2,(j,value.enc());label='dotted'
        contacts.append(dict(wall=j,label=label,product=value.enc()))
    # Equalities are identities: y is perpendicular to S; p=-sqrt2(n5+n6),
    # p^2=4, y-normalization has square -2, and p.y=0.
    facets=[]
    for i,n in enumerate(normals):
        v=[U.rnd(x-center_products[i]*a/2)for x,a in zip(center,n)]
        vv=U.rnd(pair(v,v));future=U.rnd(pair(v,center));assert vv.h<0 and future.h<0
        gates=[U.rnd(pair(v,m))for m in normals]
        assert gates[i].l<=0<=gates[i].h
        assert all(q.h<0 for j,q in enumerate(gates)if j!=i)
        facets.append(dict(wall=i,point=P.enc(v),norm=vv.enc(),future=future.enc(),products=P.enc(gates)))
    escape=C.decode(old['old60_escape']['computational'])
    ep=[U.rnd(pair(escape,n))for n in normals]
    assert all(q.h<0 for q in ep[:-1]) and ep[-1].l>0
    escape_square=U.rnd(pair(escape,escape));assert escape_square.l>0 and pair(escape,center).h<0
    triangle=[[pair(normals[i],normals[j])for j in [5,6,62]]for i in [5,6,62]]
    triangle_det=N.determinant3(triangle);assert triangle_det.h<0
    # A finite vertex incident to62 needs nine independent incident roots.
    # Its only finite preceding neighbors are six roots, so none can exist.
    finite_neighbors=[x['wall']for x in contacts if isinstance(x['label'],int)]
    assert finite_neighbors==[0,1,2,3,5,6]
    assert not old['affine_inventory']['affine_circuits']
    h4_pivots=P.pd([[pair(normals[i],normals[j])for j in [0,1,2,3]]for i in [0,1,2,3]])
    for i in [0,1,2,3]:
        for j in [5,6,62]:
            q=pair(normals[i],normals[j]);assert q.l<=0<=q.h
    affine_proof='The preceding62-wall chamber has no affine circuit. Every new connected affine circuit would contain62 and only its finite neighbors. Their diagram is the orthogonal union H4 plus the(3,4,4)triangle on5,6,62. Its proper connected subsets are spherical and its full triangle is indefinite; hence no new affine circuit exists.'
    physical=lambda v:P.combine(v,z['basis'])
    paths=[f'scripts/{TAG}.py',f'scripts/{N.TAG}.py',f'results/{N.TAG}.json',
           f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',f'results/{N.CONE}.json']
    return dict(status='PASS_EXACT_63_WALL_CAP_CUTTING_RETAINED60_ESCAPE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        old_wall_count=62,total_wall_count=63,computational_root_basis=z['basis_ids'],
        computational_gram=[P.enc(row)for row in z['exact_gram']],exact_recipe=z['cap_recipe'],
        new_normals=[dict(index=62,computational=P.enc(r),physical=P.enc(physical(r)))],
        new_contact_checks=[dict(wall=62,norm=square.enc(),contacts=contacts)],
        common_center=P.enc(center),physical_common_center=P.enc(physical(center)),
        center_norm=pair(center,center).enc(),center_products=P.enc(center_products),facet_witnesses=facets,
        old60_escape=dict(computational=P.enc(escape),physical=old['old60_escape']['physical'],
                          norm=escape_square.enc(),products=P.enc(ep),new_cap_product=ep[-1].enc()),
        new_facet_finite_vertex_census=dict(finite_preceding_neighbors=finite_neighbors,new_finite_vertices=0,
            proof='Any finite vertex needs nine independent incident normals; all pairwise incident roots have finite contacts. Only six retained roots meet the new facet, so it has no finite vertex in this63-wall chamber.'),
        A2_cap_triangle=dict(indices=[5,6,62],labels=[3,4,4],gram=[P.enc(row)for row in triangle],determinant=triangle_det.enc()),
        affine_inventory=dict(new_affine_circuits=[],H4_positive_pivots=h4_pivots,proof=affine_proof),
        S6_face_unchanged=dict(reason='Y is future timelike in S6-perp, while p lies in span(S6). Thus r.x<0 for every future hyperbolic point x in S6-perp. The new wall does not meet or truncate the actual S6 face, despite cutting its projective exterior escape witness.'),
        proof='Projection makes Y perpendicular to S6. The prescribed A2 component p=-sqrt2(n5+n6) has square4 and products -sqrt2 with5,6. Normalizing Y to square-2 therefore gives exact norm2 and the four right/two pi4 contacts. Every other native contact is strictly dotted. The old strict future center is retained, and all63 walls have explicit strict finite facet witnesses. The frozen60 escape remains strict against62 preceding walls and has positive product with the new normal, so this wall removes that escape. This cap misses the S6 face: S6 plus the new root contains the indefinite(3,4,4)triangle.',
        scope='One exact coherent escape-cutting cap, with no new finite vertex. The actual S6 face is unchanged. No full63 cone, finite-volume result, or all-parameter enumeration is asserted. Subsequent geometry must measure which ambient open regions are removed and whether other strict escapes remain.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text());print('PASS saved exact missing-S6 cap replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['old60_escape']['new_cap_product'])
