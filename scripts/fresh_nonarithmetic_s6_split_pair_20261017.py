#!/usr/bin/env python3
"""Two exact right-S6 endpoints joined at a new pi/4 finite corner."""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_nonarithmetic_escape62_20261017_cap as B

P,U,I,D,C=B.P,B.U,B.I,B.D,B.C
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_s6_split_pair_20261017'
PRIOR_CONE='fresh_nonarithmetic_escape62_20261017_cone'
S=[0,1,2,3,5,6]
FRAMES=[S+[30,44],S+[30,52]]


def kernel(z,frame,seed,center):
    pair,normals=z['pair'],z['normals'];nn=[normals[i]for i in frame]
    labels=D.E.labels()
    gram=[[I(2)if i==j else-D.theta(labels[i][j])for j in frame]for i in frame]
    pd=P.pd(gram);inverse,pivots=P.C53.inverse(gram)
    project=lambda v:P.sub(v,P.combine(U.imv(inverse,[pair(v,n)for n in nn]),nn))
    hc=[pair(center,n)for n in nn]
    k0=project(center);q0=U.rnd(pair(center,center)-U.idot(hc,U.imv(inverse,hc)));assert q0.h<0
    k0=[U.rnd(x/U.sqrt(-q0))for x in k0]
    assert pair(k0,center).h<0
    cc=[-D.theta(labels[seed][j])if labels[seed][j]else pair(normals[seed],normals[j])for j in frame]
    w=P.sub(normals[seed],P.combine(U.imv(inverse,cc),nn))
    wk=U.rnd((pair(normals[seed],center)-U.idot(cc,U.imv(inverse,hc)))/U.sqrt(-q0))
    k1=[U.rnd(x+wk*y)for x,y in zip(w,k0)]
    q1=U.rnd(I(2)-U.idot(cc,U.imv(inverse,cc))+wk*wk);assert q1.l>0
    k1=[U.rnd(x/U.sqrt(q1))for x in k1]
    return k0,k1,dict(frame=frame,spacelike_seed=seed,frame_PD_pivots=pd,inverse_pivots=pivots,
        timelike_projection_square=q0.enc(),spacelike_projection_square=q1.enc(),k0=P.enc(k0),k1=P.enc(k1),
        norm_identities='q0=H^2-h^T G^-1 h; w.k0=(seed.H-c^T G^-1 h)/sqrt(-q0); q1=2-c^T G^-1 c+(w.k0)^2',
        literal_finite_frame_labels=[[labels[i][j]for j in frame]for i in frame])


def conic(k0,k1,sign):
    a=[U.rnd((x+sign*y)/U.sqrt(I(2)))for x,y in zip(k0,k1)]
    b=[U.rnd((-x+sign*y)/U.sqrt(I(2)))for x,y in zip(k0,k1)]
    return a,b


def geometry():
    z=B.geometry();normals=z['normals'][:];pair=z['pair']
    old=C.read(B.TAG);center=C.decode(old['common_center'])
    k0,k1,kernel_a=kernel(z,FRAMES[0],52,center)
    first=[]
    for sign in [-1,1]:
        a,b=conic(k0,k1,sign);t=I(7)
        v=[U.rnd(t*x+y/t)for x,y in zip(a,b)]
        if pair(v,center).h<0 and all(pair(v,normals[j]).h<-2 for j in range(63)if j not in FRAMES[0]):
            first.append((v,dict(index=63,parameter=t.enc(),kernel_sign=sign,A=P.enc(a),B=P.enc(b),kernel=kernel_a)))
    assert len(first)==1,len(first)
    v,recipe_a=first[0]
    k0,k1,kernel_b=kernel(z,FRAMES[1],44,center)
    second=[];trials=[]
    for sign in [-1,1]:
        a,b=conic(k0,k1,sign);qa=U.rnd(pair(v,a));qb=U.sqrt(I(2));qc=U.rnd(pair(v,b))
        assert qa.h<0 or qa.l>0
        disc=U.rnd(qb*qb-4*qa*qc)
        if disc.h<0:trials.append(dict(kernel_sign=sign,discriminant=disc.enc(),status='NO_REAL_ROOT'));continue
        assert disc.l>0
        for branch in [-1,1]:
            u=U.rnd((-qb+branch*U.sqrt(disc))/(2*qa))
            if u.h<=1:
                trials.append(dict(kernel_sign=sign,branch=branch,parameter=u.enc(),status='PARAMETER_NOT_ABOVE_ONE'));continue
            assert u.l>1
            w=[U.rnd(u*x+y/u)for x,y in zip(a,b)]
            blockers=[j for j in range(63)if j not in FRAMES[1]and pair(w,normals[j]).h>=-2]
            if blockers:
                trials.append(dict(kernel_sign=sign,branch=branch,parameter=u.enc(),status='FAILED_STRICT_DOTTED_GATE',
                    blockers=[dict(wall=j,product=pair(w,normals[j]).enc())for j in blockers]));continue
            assert pair(w,center).h<0
            rec=dict(index=64,parameter=u.enc(),kernel_sign=sign,quadratic_branch=branch,A=P.enc(a),B=P.enc(b),
                mutual_label=4,equation=P.enc([qa,qb,qc]),discriminant=disc.enc(),kernel=kernel_b)
            second.append((w,rec));trials.append(dict(kernel_sign=sign,branch=branch,parameter=u.enc(),status='ACCEPTED'))
    assert len(second)==1,(len(second),trials)
    w,recipe_b=second[0];normals.extend([v,w])
    z.update(normals=normals,recipes=[recipe_a,recipe_b],second_equation_trials=trials,common_center=center)
    return z


def run():
    z=geometry();normals,pair=z['normals'],z['pair'];center=z['common_center'];old=C.read(B.TAG)
    oldcone=C.read(PRIOR_CONE);assert oldcone['status']=='PASS_COMPLETE_63_WALL_CONE_AND_STRICT_ESCAPE'
    checks=[]
    for wall in [63,64]:
        q=U.rnd(pair(normals[wall],normals[wall]));assert q.l<=2<=q.h
        right=FRAMES[wall-63];contacts=[]
        for j in range(wall):
            value=U.rnd(pair(normals[wall],normals[j]))
            if j in right:assert value.l<=0<=value.h;label=2
            elif wall==64 and j==63:
                residual=value+D.theta(4);assert residual.l<=0<=residual.h;label=4
            else:assert value.h<-2,(wall,j,value.enc());label='dotted'
            contacts.append(dict(wall=j,label=label,product=value.enc()))
        checks.append(dict(wall=wall,norm=q.enc(),contacts=contacts))
    ps=[U.rnd(pair(center,n))for n in normals];assert all(q.h<0 for q in ps)
    assert pair(center,center).h<0 and pair(center,z['h']).h<0
    facets=[]
    for i,n in enumerate(normals):
        point=[U.rnd(x-ps[i]*y/2)for x,y in zip(center,n)]
        q=U.rnd(pair(point,point));future=U.rnd(pair(point,center));assert q.h<0 and future.h<0
        gates=[U.rnd(pair(point,m))for m in normals]
        assert gates[i].l<=0<=gates[i].h and all(v.h<0 for j,v in enumerate(gates)if j!=i)
        facets.append(dict(wall=i,point=P.enc(point),norm=q.enc(),future=future.enc(),products=P.enc(gates)))
    expected=[S+[30,44,63],S+[30,52,64],S+[30,63,64]];vertices=[]
    for support in expected:
        nn=[normals[i]for i in support];gram=[[pair(a,b)for b in nn]for a in nn]
        pd=P.pd(gram);inverse,pivots=P.C53.inverse(gram)
        point=P.sub(center,P.combine(U.imv(inverse,[pair(center,n)for n in nn]),nn))
        norm,future=U.rnd(pair(point,point)),U.rnd(pair(point,center));assert norm.h<0 and future.h<0
        gates=[U.rnd(pair(point,n))for n in normals]
        for j,v in enumerate(gates):
            if j in support:assert v.l<=0<=v.h
            else:assert v.h<0,(support,j,v.enc())
        vertices.append(dict(incident=support,point=P.enc(point),norm=norm.enc(),future=future.enc(),
            PD_pivots=pd,inverse_pivots=pivots,products=P.enc(gates)))
    censuses=[]
    for wall,neighbors in[(63,FRAMES[0]),(64,FRAMES[1]+[63])]:
        records=[]
        for chosen in combinations(neighbors,8):
            support=list(chosen)+[wall]
            if wall==64 and 52 in chosen and 63 in chosen:
                value=pair(normals[52],normals[63]);assert value.h<-2
                records.append(dict(support=support,status='DOTTED_PAIR',pair=[52,63],product=value.enc()))
            else:
                match=next(v for v in vertices if set(v['incident'])==set(support))
                records.append(dict(support=support,status='ACTUAL_FINITE_VERTEX',point=match['point']))
        censuses.append(dict(wall=wall,stage_wall_count=wall+1,finite_preceding_neighbors=neighbors,records=records))
    assert [sum(r['status']=='ACTUAL_FINITE_VERTEX'for r in c['records'])for c in censuses]==[1,2]
    assert oldcone['computational_root_basis']==z['basis_ids']
    assert oldcone['computational_gram']==[P.enc(row)for row in z['exact_gram']]
    retained=[]
    for row in oldcone['final_rays']:
        if row['kind']!='finite_vertices':continue
        point=C.decode(row['normalized_computational_coordinate_bounds'])
        products=[U.rnd(pair(point,normals[j]))for j in [63,64]]
        assert all(v.h<0 for v in products),(row['active'],P.enc(products))
        retained.append(dict(active=row['active'],new_products=P.enc(products)))
    assert len(retained)==193
    # The only nonright finite edge involving either new root is63--64,
    # an isolated I2(4) in every connected finite diagram containing them.
    assert old['affine_inventory']['new_affine_circuits']==[]
    old_escape=C.decode(oldcone['strict_escape']['computational_point'])
    escape_products=[U.rnd(pair(old_escape,n))for n in normals]
    assert all(v.h<0 for v in escape_products[:63])
    signs=['negative'if v.h<0 else'positive'if v.l>0 else'UNRESOLVED'for v in escape_products[63:]]
    assert 'UNRESOLVED'not in signs
    physical=lambda v:P.combine(v,z['basis'])
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'results/{B.TAG}.json',
        f'scripts/{PRIOR_CONE}.py',f'results/{PRIOR_CONE}.json',
        f'scripts/{D.E.TAG}.py',f'results/{D.E.TAG}.json']
    return dict(status='PASS_EXACT_65_WALL_SPLIT_ENDPOINT_PI4_BRIDGE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        old_wall_count=63,total_wall_count=65,computational_root_basis=z['basis_ids'],
        computational_gram=[P.enc(row)for row in z['exact_gram']],exact_recipes=z['recipes'],
        second_parameter_trials=z['second_equation_trials'],
        new_normals=[dict(index=i,computational=P.enc(normals[i]),physical=P.enc(physical(normals[i])))for i in [63,64]],
        new_contact_checks=checks,common_center=P.enc(center),physical_common_center=P.enc(physical(center)),
        center_norm=pair(center,center).enc(),center_products=P.enc(ps),facet_witnesses=facets,
        face_support=S,old_edges_terminated=[[30,44],[30,52]],mutual_corner_support=S+[30,63,64],
        new_finite_vertices=vertices,new_facet_finite_vertex_census=censuses,
        preserved_old_finite_vertices=retained,complete_finite_vertex_count=196,
        affine_inventory=dict(new_affine_circuits=[],proof='The preceding63 chamber has no affine circuit. Every new finite nonright bond is the single63--64 pi4 bond. Thus a connected affine circuit containing a new node could only be A1 or I2(4), both spherical. All nonright old/new relations are dotted and cannot belong to an affine diagram.'),
        old63_escape=dict(computational=P.enc(old_escape),products=P.enc(escape_products),new_cap_signs=signs),
        proof='Each rank8 positive frame leaves a Lorentz plane. Projecting the strict center gives its future unit timelike k0; the projected opposite old wall gives the unique unit spacelike k1 up to sign. The exact conic r(t)=t(k0+s*k1)/sqrt2+(-k0+s*k1)/(sqrt2*t) has norm2 and is perpendicular to the frame. Set the first parameter to7. The second parameter solves the displayed quadratic making the mutual product -sqrt2. All127 new contacts and a common strict future center are verified. The three finite corners satisfy every nonincident old/new inequality strictly. Finite-neighbor support exhaustion proves one new vertex at63 and two at64. No full S6 or ambient closure is asserted.',
        scope='A coherent two-wall cluster joins separately terminating old edges through a genuine new pi4 finite corner. It avoids the earlier same-wall two-endpoint obstruction. No full65 cone or finite-volume conclusion; all root definitions are exact but arithmetic is shared with the frozen physical interval engine.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(path.read_text());print('PASS saved exact split pair replay')
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['old63_escape']['new_cap_signs'],flush=True)
