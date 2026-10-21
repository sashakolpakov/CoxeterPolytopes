#!/usr/bin/env python3
"""Two exact caps close the actual H4+I2(9) face as a compact triangular prism."""
import argparse
import itertools
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_pi8_quad_completion_20261017 as C

P,U,I,D=C.P,C.U,C.I,C.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_pi9_prism_completion_20261017'
CONE='fresh_nonarithmetic_pi8_quad_cone_20261017'
SUPPORT=[0,1,2,3,58,59]
OLD_SIDES=[5,6,55]
PRESCRIBED={60:{**{i:2 for i in SUPPORT+[5]},6:5,55:4},61:{i:2 for i in SUPPORT+[6,55,60]}}


def geometry():
    z=C.geometry();normals=z['normals'][:];pair=z['pair'];old=C.read(C.TAG)
    corner=C.decode(next(r for r in old['complete_six_pair_intersections']if set(r['pair'])=={58,59})['point'])
    square=pair(corner,corner);assert square.h<0 and pair(corner,z['h']).h<0
    k=[U.rnd(x/U.sqrt(-square))for x in corner]
    sides=[normals[i]for i in OLD_SIDES];gram=[[pair(a,b)for b in sides]for a in sides]
    inverse,pivots=P.C53.inverse(gram)
    contacts=[I(0),-D.theta(5),-D.theta(4)]
    projection=P.combine(U.imv(inverse,contacts),sides)
    residual=U.rnd(pair(projection,projection)-2);assert residual.l>0
    v=[U.rnd(p+U.sqrt(residual)*x)for p,x in zip(projection,k)]
    normals.append(v)
    # A polar truncation of the hyperideal corner of the four-side chamber.
    polar_frame=SUPPORT+[6,55,60]
    frame=[normals[i]for i in polar_frame]
    fg=[[pair(a,b)for b in frame]for a in frame]
    fi,fp=P.C53.inverse(fg)
    w=P.sub(normals[5],P.combine(U.imv(fi,[pair(normals[5],a)for a in frame]),frame))
    w_square=U.rnd(pair(w,w));assert w_square.l>2
    u=[U.rnd(-U.sqrt(2/w_square)*x)for x in w]
    normals.append(u)
    z.update(normals=normals,old=old,old_corner=corner,
        recipes=dict(v=dict(index=60,old_side_order=OLD_SIDES,old_side_contacts=P.enc(contacts),
                           projection=P.enc(projection),positive_norm_excess=residual.enc(),future_unit_kernel=P.enc(k),inverse_pivots=pivots),
                     u=dict(index=61,right_frame=polar_frame,projection_seed=5,orthogonal_projection=P.enc(w),
                            projection_square=w_square.enc(),inverse_pivots=fp)))
    return z


def determinant3(a):
    return U.rnd(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
                 -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
                 +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def affine_inventory(contacts):
    AF=C.AF
    labels=D.E.labels()
    tetra=C.read(C.C.TAG)
    for block in tetra['new_contact_checks']+C.read(C.TAG)['new_contact_checks']+contacts:
        wall=block['wall'];assert wall==len(labels)
        for row in labels:row.append(0)
        labels.append([0]*(wall+1))
        for contact in block['contacts']:
            m=contact['label'];labels[wall][contact['wall']]=labels[contact['wall']][wall]=m if isinstance(m,int)else 0
    assert len(labels)==62
    queue=[(i,)for i in range(62)];seen=set(queue);counts=Counter();affine=[]
    while queue:
        support=queue.pop()
        gram=tuple(tuple(AF.scale(AF.ONE,2)if i==j else AF.VALUES[labels[i][j]]for j in support)for i in support)
        kind,rank,nullity=AF.classify(gram);counts[kind]+=1
        if kind=='positive_semidefinite':affine.append(list(support));continue
        if kind!='positive_definite':continue
        assert len(support)<=9
        for j in range(62):
            if j in support:continue
            ls=[labels[j][i]for i in support]
            if any(m not in AF.VALUES for m in ls)or all(m==2 for m in ls):continue
            ss=tuple(sorted(support+(j,)))
            if ss not in seen:seen.add(ss);queue.append(ss)
    assert not affine
    return dict(counts=dict(counts),affine_circuits=affine,
                scope='Complete connected standard-label affine-circuit inventory by the frozen classification/growth argument; no finite-volume assertion.')


def run():
    z=geometry();normals,pair=z['normals'],z['pair'];checks=[]
    for i in(60,61):
        square=U.rnd(pair(normals[i],normals[i]));assert square.l<=2<=square.h
        contacts=[]
        for j in range(i):
            value=U.rnd(pair(normals[i],normals[j]))
            if j in PRESCRIBED[i]:
                label=PRESCRIBED[i][j];residual=value+D.theta(label);assert residual.l<=0<=residual.h
            else:assert value.h<-2,(i,j,value.enc());label='dotted'
            contacts.append(dict(wall=j,label=label,product=value.enc()))
        checks.append(dict(wall=i,norm=square.enc(),contacts=contacts))
    corner,oldcenter=z['old_corner'],z['h']
    assert all(pair(corner,normals[i]).h<0 for i in(60,61))
    for power in range(80):
        eps=I(1)/I(2**power);center=[U.rnd(v+eps*h)for v,h in zip(corner,oldcenter)]
        products=[U.rnd(pair(center,n))for n in normals]
        if all(q.h<0 for q in products):break
    else:raise AssertionError('no common center')
    assert pair(center,center).h<0 and pair(center,oldcenter).h<0
    facets=[]
    for i,n in enumerate(normals):
        point=[U.rnd(x-products[i]*y/2)for x,y in zip(center,n)]
        q=U.rnd(pair(point,point));future=U.rnd(pair(point,center));assert q.h<0 and future.h<0
        gates=[U.rnd(pair(point,m))for m in normals]
        assert gates[i].l<=0<=gates[i].h and all(x.h<0 for j,x in enumerate(gates)if j!=i)
        facets.append(dict(wall=i,point=P.enc(point),norm=q.enc(),future=future.enc(),products=P.enc(gates)))
    support=[normals[i]for i in SUPPORT]
    support_gram=[[pair(a,b)for b in support]for a in support]
    pd=P.pd(support_gram);si,sp=P.C53.inverse(support_gram)
    project=lambda p:P.sub(p,P.combine(U.imv(si,[pair(p,n)for n in support]),support))
    face_center=project(center);assert pair(face_center,face_center).h<0 and pair(face_center,center).h<0
    sides=[5,6,55,60,61];basis_ids=[5,6,55,60];basis=[normals[i]for i in basis_ids]
    bg=[[pair(a,b)for b in basis]for a in basis];bi,bp=P.C53.inverse(bg)
    for i in sides:
        assert pair(face_center,normals[i]).h<0
        for j in SUPPORT:
            q=pair(normals[i],normals[j]);assert q.l<=0<=q.h
    bases={5,61};rectangles={6,55,60};rays=[]
    for triple in itertools.combinations(sides,3):
        equations=[[pair(normals[i],b)for b in basis]for i in triple]
        coefficients=[U.rnd(((-1)**j)*determinant3([[row[k]for k in range(4)if k!=j]for row in equations]))for j in range(4)]
        point=P.combine(coefficients,basis);square=U.rnd(pair(point,point))
        expected=len(set(triple)&bases)==1 and len(set(triple)&rectangles)==2
        if expected:
            assert square.h<0
            if pair(point,center).l>0:point=[-x for x in point]
            future=U.rnd(pair(point,center));assert future.h<0
            incident=SUPPORT+list(triple)
            pivots=P.pd([[pair(normals[i],normals[j])for j in incident]for i in incident])
            gates=[U.rnd(pair(point,n))for n in normals]
            for j,q in enumerate(gates):
                if j in incident:assert q.l<=0<=q.h
                else:assert q.h<0,(triple,j,q.enc())
            rays.append(dict(sides=list(triple),kind='finite_vertex',point=P.enc(point),norm=square.enc(),future=future.enc(),
                             incident=incident,PD_pivots=pivots,products=P.enc(gates)))
        else:
            assert square.l>0
            rem={i:U.rnd(pair(point,normals[i]))for i in sides if i not in triple}
            plus=[i for i,q in rem.items()if q.l>0];minus=[i for i,q in rem.items()if q.h<0]
            assert plus and minus,(triple,{i:q.enc()for i,q in rem.items()})
            rays.append(dict(sides=list(triple),kind='both_orientations_infeasible',point=P.enc(point),norm=square.enc(),
                             other_side_products={str(i):q.enc()for i,q in rem.items()},positive_blocker=plus[0],negative_blocker=minus[0]))
    assert sum(r['kind']=='finite_vertex'for r in rays)==6
    # Complete new-facet finite-corner supports: a finite vertex may contain
    # only finite neighbors. Every invalid support here retains an explicit
    # indefinite three-node block, so no numerical eigenvalue classification.
    new_census=[]
    for wall,neighbors in[(60,SUPPORT+[5,6,55]),(61,SUPPORT+[6,55,60])]:
        records=[]
        for omitted in neighbors:
            active=[i for i in neighbors if i!=omitted]+[wall]
            triple=[6,55,60]
            if all(i in active for i in triple):
                determinant=determinant3([[pair(normals[i],normals[j])for j in triple]for i in triple]);assert determinant.h<0
                records.append(dict(omitted=omitted,status='INDEFINITE_PRINCIPAL_TRIPLE',triple=triple,determinant=determinant.enc()))
            else:
                pivots=P.pd([[pair(normals[i],normals[j])for j in active]for i in active])
                match=next(r for r in rays if r['kind']=='finite_vertex'and set(r['incident'])==set(active))
                records.append(dict(omitted=omitted,status='ACTUAL_FINITE_VERTEX',incident=active,PD_pivots=pivots,prism_triple=match['sides']))
        new_census.append(dict(wall=wall,stage_wall_count=wall+1,finite_neighbor_count=9,records=records))
    assert [sum(r['status']=='ACTUAL_FINITE_VERTEX'for r in x['records'])for x in new_census]==[2,3]
    # Parent's complete60 cone is a physical-coordinate export.
    oldcone=C.read(CONE);escape_physical=C.decode(oldcone['strict_escape']['point'])
    escape=z['coordinates'](escape_physical)
    escape_products=[U.rnd(pair(escape,n))for n in normals]
    assert all(q.h<0 for q in escape_products[:60])
    escape_norm=U.rnd(pair(escape,escape));escape_future=U.rnd(pair(escape,center))
    assert escape_norm.l>0 and escape_future.h<0
    escape_signs=['negative'if q.h<0 else'positive'if q.l>0 else'UNRESOLVED'for q in escape_products[60:]]
    assert 'UNRESOLVED'not in escape_signs
    affine=affine_inventory(checks)
    physical=lambda x:P.combine(x,z['basis'])
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',f'scripts/{C.C.TAG}.py',f'results/{C.C.TAG}.json',
           f'scripts/{D.TAG}.py',f'results/{D.TAG}.json',f'scripts/{CONE}.py',f'results/{CONE}.json',
           'scripts/fresh_h9_persistent_cusp_affine_20261016.py','results/fresh_h9_persistent_cusp_affine_20261016.json']
    return dict(status='PASS_EXACT_62_WALL_CHAMBER_WITH_COMPACT_PI9_TRIANGULAR_PRISM',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        old_walls=60,total_walls=62,computational_root_basis=z['basis_ids'],computational_gram=[P.enc(row)for row in z['exact_gram']],
        exact_recipes=z['recipes'],new_normals=[dict(index=i,computational=P.enc(normals[i]),physical=P.enc(physical(normals[i])))for i in(60,61)],
        new_contact_checks=checks,common_center=P.enc(center),physical_common_center=P.enc(physical(center)),center_perturbation_power=power,
        center_norm=pair(center,center).enc(),center_future=pair(center,oldcenter).enc(),center_products=P.enc(products),facet_witnesses=facets,
        actual_H3_support=SUPPORT,face_support=SUPPORT,support_PD_pivots=pd,support_inverse_pivots=sp,prism_sides=sides,prism_bases=sorted(bases),prism_rectangles=sorted(rectangles),
        face_basis=basis_ids,face_basis_inverse_pivots=bp,face_center=P.enc(face_center),face_Gram=[P.enc([pair(normals[i],normals[j])for j in sides])for i in sides],
        complete_ten_triple_intersections=rays,complete_triple_intersections=rays,new_facet_finite_vertex_census=new_census,new_finite_vertices=5,
        old60_escape=dict(physical=oldcone['strict_escape']['point'],computational=P.enc(escape),norm=escape_norm.enc(),future=escape_future.enc(),
                          products=P.enc(escape_products),new_cap_signs=escape_signs),affine_inventory=affine,
        proof='The exact normal v is the future branch of its three-contact norm equation in the Lorentz4 complement of H4+I2(9). The second cap is the oriented normalized polar of the indefinite nine-root span S+[6,55,v]. All121 new native contacts, a strict future common center and all62 strict facet witnesses are checked. In the four-dimensional face span every triple intersection of the five sides is enumerated: six feasible future timelike rays and four lines infeasible in both orientations. The six vertices satisfy every other native inequality strictly, so their compact convex hull is exactly the actual supported H3 face. Sequential finite-neighbor support exhaustion gives exactly two vertices created byv and three further vertices created byu; two of the latter also lie onv. Old finite vertices persist under this coherent acute extension.',
        scope='A simultaneous two-wall completion of one entire necessary compact H3 face. The original S6 face and ambient chamber need not be compact. The saved old60 strict escape is tested explicitly; no full ambient cone or finite-volume conclusion is made here.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text());print('PASS saved compact pi9 prism replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['old60_escape']['new_cap_signs'])
