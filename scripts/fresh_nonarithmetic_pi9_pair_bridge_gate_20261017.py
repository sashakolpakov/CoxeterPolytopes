#!/usr/bin/env python3
"""Complete two-wall bridge gate around the fixed pi9 edge of the 60-wall chamber."""
import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_pi8_quad_completion_20261017 as C

P,U,I,D=C.P,C.U,C.I,C.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_pi9_pair_bridge_gate_20261017'
DOMAIN='fresh_nonarithmetic_pi9_endpoint_domain_20261017'
S=[0,1,2,3,5,6]


def enc(x):
    return x.enc()


def gap(value):
    """Strict consecutive cosine bracket, hence a gate for every finite label."""
    for m in range(2,40):
        lower,upper=D.theta(m),D.theta(m+1)
        magnitude=-value
        if lower.h<magnitude.l and magnitude.h<upper.l:
            return dict(consecutive_labels=[m,m+1],magnitude=magnitude.enc(),
                        lower_cosine=lower.enc(),upper_cosine=upper.enc())
    return None


def old_gates(normal,prescribed,normals,pair):
    products=[U.rnd(pair(normal,n))for n in normals]
    checks=[]
    for i,value in enumerate(products):
        if i in prescribed:
            residual=value+D.theta(prescribed[i])
            assert residual.l<=0<=residual.h,(i,residual.enc())
            checks.append(dict(wall=i,label=prescribed[i],product=value.enc()))
        elif value.l>0:
            return dict(status='REJECT_OBTUSE',wall=i,product=value.enc())
        elif value.h<-2:
            checks.append(dict(wall=i,label='dotted',product=value.enc()))
        else:
            bracket=gap(value)
            assert bracket is not None,('unresolved old contact',i,value.enc())
            return dict(status='REJECT_CONSECUTIVE_ANGLE_GAP',wall=i,product=value.enc(),**bracket)
    return dict(status='PASS_ALL_OLD_COXETER_GATES',contacts=checks)


def event(projection,scale,aa,bb,anchor,label,normals,pair):
    qa=U.rnd(scale*pair(aa,normals[anchor]))
    qb=U.rnd(pair(projection,normals[anchor])+D.theta(label))
    qc=U.rnd(scale*pair(bb,normals[anchor]))
    assert qa.h<0<qc.l
    discriminant=U.rnd(qb*qb-4*qa*qc);assert discriminant.l>0
    t=U.rnd((-qb-U.sqrt(discriminant))/(2*qa));assert t.l>1
    normal=[U.rnd(p+scale*(t*a+b/t))for p,a,b in zip(projection,aa,bb)]
    square=U.rnd(pair(normal,normal));assert square.l<=2<=square.h
    return normal,t,dict(equation=P.enc([qa,qb,qc]),discriminant=discriminant.enc(),parameter=t.enc(),normal=P.enc(normal),square=square.enc())


def run():
    z=C.geometry();normals,pair=z['normals'],z['pair'];old=C.read(C.TAG);domain=C.read(DOMAIN)
    # B terminates both old edges [30,59] and [30,49].
    frame_ids=S+[30,59];frame=[normals[i]for i in frame_ids]
    gram=[[pair(a,b)for b in frame]for a in frame]
    pd=P.pd(gram);inverse,pivots=P.C53.inverse(gram)
    endpoint=C.decode(next(r for r in old['complete_six_pair_intersections']if set(r['pair'])=={30,59})['point'])
    square=pair(endpoint,endpoint);assert square.h<0
    k0=[U.rnd(x/U.sqrt(-square))for x in endpoint]
    k1=[U.rnd(x/U.sqrt(I(2)))for x in normals[55]]
    aa=[U.rnd((a-b)/U.sqrt(I(2)))for a,b in zip(k0,k1)]
    bb=[U.rnd((-a-b)/U.sqrt(I(2)))for a,b in zip(k0,k1)]
    ap=[U.rnd(pair(aa,n))for n in normals];bp=[U.rnd(pair(bb,n))for n in normals]
    for i in range(60):
        if i in frame_ids:
            assert ap[i].l<=0<=ap[i].h and bp[i].l<=0<=bp[i].h
        else:
            assert ap[i].h<0
            difference=U.rnd(ap[i]-bp[i])
            if i==55:assert difference.l<=0<=difference.h
            else:assert difference.h<0
    assert (ap[55]+1).l<=0<=(ap[55]+1).h and (bp[55]+1).l<=0<=(bp[55]+1).h
    i5diag=U.rnd(2/(4-D.theta(5)*D.theta(5)))
    assert (I(2)/3+i5diag).l>2 and (2*i5diag).l>2
    b_specs=[(None,2,m)for m in(2,3)]+[(node,m,2)for node in(5,6)for m in(3,4,5)]+[(30,3,2)]+[(59,3,m)for m in(2,3)]
    b_rows=[];b_survivors=[]
    for index,(node,label,m49)in enumerate(b_specs):
        contacts=[-D.theta(label)if i==node else I(0)for i in frame_ids]
        projection=P.combine(U.imv(inverse,contacts),frame)
        cost=U.rnd(pair(projection,projection));rho=U.rnd(2-cost);assert rho.l>0
        scale=U.sqrt(rho/2)
        normal,t,recipe=event(projection,scale,aa,bb,49,m49,normals,pair)
        prescribed={i:2 for i in frame_ids+[49]};prescribed[49]=m49
        if node is not None:prescribed[node]=label
        gates=old_gates(normal,prescribed,normals,pair)
        row=dict(index=index,nonright_node=node,label=label,label49=m49,
                 projection=P.enc(projection),projection_cost=cost.enc(),residual_square=rho.enc(),
                 recipe=recipe,old_gates=gates)
        if gates['status'].startswith('PASS'):
            # Literal second finite corner, independently checked after elimination.
            second=[normals[i]for i in S+[30,49]]+[normal]
            row['second_corner_PD_pivots']=P.pd([[pair(a,b)for b in second]for a in second])
            b_survivors.append((index,normal))
        b_rows.append(row)
    assert len(b_rows)==11 and [i for i,_ in b_survivors]==[0,1,8,9,10]

    # A terminates [58,59] and [49,58]. The zero A2 profile has an infinite
    # angle-to49 family. Nonzero A2 profiles have a finite positive-cost list.
    a0,b0=C.decode(domain['A']),C.decode(domain['B'])
    a_rows=[];a_survivors=[]
    for profile_index,pf in enumerate(domain['profiles']):
        node,label=pf['nonright_A2_node'],pf['label']
        cost=I(*pf['retained_cost']);projection=C.decode(pf['projection']);scale=I(*pf['scale'])
        labels49=[2]if profile_index==0 else [m for m in(2,3,4,5)if (cost+D.theta(m)*D.theta(m)/2).h<2]
        if profile_index:
            rejected=[m for m in(2,3,4,5,6)if m not in labels49]
            assert all((cost+D.theta(m)*D.theta(m)/2).l>2 for m in rejected)
        for m49 in labels49:
            normal,t,recipe=event(projection,scale,a0,b0,49,m49,normals,pair)
            prescribed={i:2 for i in domain['edge_support']+[49]};prescribed[49]=m49
            if node is not None:prescribed[node]=label
            gates=old_gates(normal,prescribed,normals,pair)
            index=len(a_rows)
            row=dict(index=index,profile_index=profile_index,nonright_A2_node=node,label=label,label49=m49,
                     second_corner_projection_cost=(cost+D.theta(m49)*D.theta(m49)/2).enc(),
                     recipe=recipe,old_gates=gates)
            if gates['status'].startswith('PASS'):
                second=[normals[i]for i in S+[49,58]]+[normal]
                row['second_corner_PD_pivots']=P.pd([[pair(a,b)for b in second]for a in second])
                a_survivors.append((index,normal))
            a_rows.append(row)
    assert len(a_rows)==15 and len(a_survivors)==5

    # All-label zero-profile tail, in fact all angle magnitudes>=1, not a sample.
    zero=domain['profiles'][0];zero_projection=C.decode(zero['projection']);zero_scale=I(*zero['scale'])
    tail_normal,t3,tail_event=event(zero_projection,zero_scale,a0,b0,49,3,normals,pair)
    a49,b49=U.rnd(pair(a0,normals[49])),U.rnd(pair(b0,normals[49]))
    assert a49.h<0<b49.l
    tail_checks=[]
    for index,normal in b_survivors:
        a=U.rnd(pair(a0,normal));b=U.rnd(pair(b0,normal));assert a.h<0
        value=U.rnd(a*t3+b/t3);derivative=U.rnd(a-b/(t3*t3))
        assert value.h<-2 and derivative.h<0
        tail_checks.append(dict(B_index=index,A_coefficient=a.enc(),B_coefficient=b.enc(),
                                value_at_label3=value.enc(),derivative_at_label3=derivative.enc()))

    pair_rows=[]
    for ai,a in a_survivors:
        for bi,b in b_survivors:
            value=U.rnd(pair(a,b))
            if value.h<-2:
                row=dict(A_index=ai,B_index=bi,status='DOTTED_NO_FINITE_MUTUAL_EDGE',product=value.enc())
            else:
                bracket=gap(value);assert bracket is not None,(ai,bi,value.enc())
                row=dict(A_index=ai,B_index=bi,status='CONSECUTIVE_ANGLE_GAP',product=value.enc(),**bracket)
            pair_rows.append(row)
    counts=dict(Counter(r['status']for r in pair_rows))
    assert counts=={'CONSECUTIVE_ANGLE_GAP':2,'DOTTED_NO_FINITE_MUTUAL_EDGE':23}
    assert sorted(tuple(r['consecutive_labels'])for r in pair_rows if 'consecutive_labels'in r)==[(10,11),(12,13)]
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',f'scripts/{DOMAIN}.py',f'results/{DOMAIN}.json',
           f'scripts/{D.TAG}.py',f'results/{D.TAG}.json',f'scripts/{D.E.TAG}.py',f'results/{D.E.TAG}.json']
    return dict(status='PASS_NO_TWO_WALL_COMPONENT_A_EDGE_BRIDGE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        retained_wall_count=60,S6=S,required_finite_corners=dict(A=[[58,59],[49,58]],B=[[30,59],[30,49]],
        mutual=[[49,'A','B'],[59,'A','B']]),
        B_frame=frame_ids,B_frame_PD_pivots=pd,B_inverse_pivots=pivots,B_old_finite_endpoint=P.enc(endpoint),
        B_kernel_A=P.enc(aa),B_kernel_B=P.enc(bb),B_A_products=P.enc(ap),B_B_products=P.enc(bp),
        I2_5_inverse_diagonal=i5diag.enc(),B_recipe_count=11,B_recipes=b_rows,B_survivor_indices=[i for i,_ in b_survivors],
        A_finite_recipe_count=15,A_recipes=a_rows,A_survivor_indices=[i for i,_ in a_survivors],
        A_zero_profile_entire_label_ge3_tail=dict(label3_event=tail_event,old49_A=a49.enc(),old49_B=b49.enc(),mutual_checks=tail_checks),
        mutual_pair_count=25,mutual_counts=counts,mutual_pairs=pair_rows,
        complete_profile_proof='For A the pi9 finite endpoint forces right H4/I2(9), leaving zero or one A2 bond3/4/5. In its second endpoint,49 and58 are orthogonal A1 factors; the cost is A2cost+theta49^2/2. Hence nonzero A2 profiles allow49 labels2..5,2..3,or2 respectively; the zero profile allows every finite label. For B, both endpoint frames are H4+A2+I2(5), with common A2 and node30. Since the I2(5) inverse diagonal d5 has2/3+d5>2 and2d5>2, the only contacts are zero, one A2 bond3/4/5, or one bond3 at one dihedral node. Intersecting the two frames gives exactly the11 B recipes.',
        orientation_proof='Both endpoint conics use their old finite endpoint and old55 as a normalized Lorentz2 frame. The opposite normal orientation is obtuse to55. The surviving endpoint x=t*A-B/t has old55 product-t+1/t, forcing t>=1. At1 it passes through the old finite vertex and cannot share an acute strict interior. Thus t>1. The saved A_i<0 and A_i-B_i<=0 prove every endpoint at t>1 satisfies every old halfspace. A specified49 angle gives the displayed quadratic with opposite-sign roots, so there is one possible t>1 for every listed recipe.',
        tail_proof='For the zero A2 profile the magnitude of A(t).49 is strictly increasing because its coefficients have A49<0<B49. All finite labels>=3 therefore have t>=t3. For each surviving B, f(t)=a*t+b/t has a<0, f(t3)<-2 and fprime(t3)<0. If b>=0 its derivative is always negative; if b<0 the derivative decreases as t grows. Therefore every member of this entire continuous tail is dotted to every B, so cannot make a finite mutual corner.',
        scope='There are no two added coherent Coxeter walls A,B which realize all six specified finite corners relative to the fixed60 walls. The test allows every finite label and every field. It does not exclude three-or-more-wall bridges, other combinatorics, replacing old walls, or changing the constructive pi9 choice. No finite-volume or host-wide nonexistence assertion is made.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text());print('PASS saved all-label two-wall bridge replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['mutual_counts'])
