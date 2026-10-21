#!/usr/bin/env python3
"""Complete finite-terminal profiles and continuous domain of the new pi9 edge."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_pi8_quad_completion_20261017 as C

P,U,I,D=C.P,C.U,C.I,C.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_pi9_endpoint_domain_20261017'
FRAME=[0,1,2,3,5,6,58,59]


def run():
    old=C.read(C.TAG);z=C.geometry();normals=z['normals'];pair=z['pair']
    endpoint_record=next(r for r in old['complete_six_pair_intersections'] if set(r['pair'])=={58,59})
    endpoint=C.decode(endpoint_record['point'])
    assert set(endpoint_record['incident'])==set(FRAME+[55])
    q0=pair(endpoint,endpoint);assert q0.h<0
    k0=[U.rnd(x/U.sqrt(-q0)) for x in endpoint]
    k1=[U.rnd(x/U.sqrt(I(2))) for x in normals[55]]
    assert pair(k0,z['h']).h<0
    frame=[normals[i]for i in FRAME]
    gram=[[pair(a,b)for b in frame]for a in frame]
    pd=P.pd(gram);inverse,pivots=P.C53.inverse(gram)
    # H4 cannot have a nonright finite extension. The I2(9) inverse
    # diagonal is strictly greater than2. The remaining A2 admits exactly
    # zero, or one bond3/4/5 at either node.
    i9=D.theta(9);dihedral_diag=U.rnd(2/(4-i9*i9));assert dihedral_diag.l>2
    a2=[[I(2),I(-1)],[I(-1),I(2)]];a2i,_=P.C53.inverse(a2)
    double_cost=U.idot([I(1),I(1)],U.imv(a2i,[I(1),I(1)]))
    assert double_cost.l<=2<=double_cost.h
    single_costs=[]
    for label in(3,4,5,6):
        magnitude=D.theta(label)
        cost=U.rnd(magnitude*magnitude*a2i[0][0])
        if label<6:assert cost.h<2
        else:assert cost.l<=2<=cost.h
        single_costs.append(dict(label=label,cost=cost.enc()))
    aa=[U.rnd((x-y)/U.sqrt(I(2)))for x,y in zip(k0,k1)]
    bb=[U.rnd((-x-y)/U.sqrt(I(2)))for x,y in zip(k0,k1)]
    ap=[U.rnd(pair(aa,n))for n in normals]
    bp=[U.rnd(pair(bb,n))for n in normals]
    for i in range(60):
        if i in FRAME:
            assert ap[i].l<=0<=ap[i].h and bp[i].l<=0<=bp[i].h
        else:
            assert ap[i].h<0
            at_one=U.rnd(ap[i]-bp[i])
            if i==55:assert at_one.l<=0<=at_one.h
            else:assert at_one.h<0
    assert ap[55].l<=-1<=ap[55].h and bp[55].l<=-1<=bp[55].h
    duals=[P.combine([inverse[j][FRAME.index(i)]for j in range(8)],frame)for i in(5,6)]
    profiles=[]
    specs=[(None,2)]+[(node,label)for node in(5,6)for label in(3,4,5)]
    for node,label in specs:
        magnitudes=[D.theta(label)if i==node else I(0)for i in(5,6)]
        projection=P.combine([-x for x in magnitudes],duals)
        cost=U.rnd(pair(projection,projection));rho=U.rnd(2-cost);assert rho.l>0
        scale=U.sqrt(rho/2)
        offsets=[U.rnd(pair(projection,n))for n in normals]
        assert offsets[55].l<=0<=offsets[55].h
        for power in range(1,80):
            threshold=I(2**power)
            products=[U.rnd(c+scale*(a*threshold+b/threshold))for c,a,b in zip(offsets,ap,bp)]
            derivatives=[U.rnd(scale*(a-b/(threshold*threshold)))for a,b in zip(ap,bp)]
            if all(products[i].h<-2 and derivatives[i].h<0 for i in range(60)if i not in FRAME):break
        else:raise AssertionError('no certified all-dotted tail')
        cap=[U.rnd(p+scale*(threshold*a+b/threshold))for p,a,b in zip(projection,aa,bb)]
        capnorm=pair(cap,cap);assert capnorm.l<=2<=capnorm.h
        cap_endpoint=[U.rnd(threshold*a-b/threshold)for a,b in zip(aa,bb)]
        assert pair(cap_endpoint,cap_endpoint).h<0 and pair(cap_endpoint,z['h']).h<0
        assert pair(cap,endpoint).h<0
        prescribed={**{i:2 for i in FRAME},**({node:label}if node is not None else{})}
        for i in FRAME:
            residual=products[i]+D.theta(prescribed[i]);assert residual.l<=0<=residual.h
        endpoint_products=[U.rnd(pair(cap_endpoint,n))for n in normals]
        for i,q in enumerate(endpoint_products):
            if i in FRAME:assert q.l<=0<=q.h
            else:assert q.h<0
        profiles.append(dict(nonright_A2_node=node,label=label,projection=P.enc(projection),
            retained_cost=cost.enc(),residual_square=rho.enc(),scale=scale.enc(),
            old60_offsets=P.enc(offsets),all_dotted_tail_threshold=2**power,
            threshold_products=P.enc(products),threshold_product_derivatives=P.enc(derivatives),
            threshold_normal=P.enc(cap),threshold_endpoint=P.enc(cap_endpoint),
            threshold_endpoint_products=P.enc(endpoint_products),
            old_finite_vertex_product=pair(cap,endpoint).enc()))
    assert len(profiles)==7
    # A single wall cannot terminate both [58,59] and [30,59]. The latter
    # support is H4+A2+I2(5). Since r.59=0, a nonright30 contact costs its
    # I2(5) inverse diagonal. Together with any nonzero A2 contact this
    # already exceeds2. Only the zero profile and label3 at30 remain;
    # label2 is t=1, while every label>=4 exceeds the I2(5) budget.
    theta5=D.theta(5);i5diag=U.rnd(2/(4-theta5*theta5))
    nonzero_cost=U.rnd(a2i[0][0]+i5diag)
    label4_cost=U.rnd(2*i5diag)
    assert nonzero_cost.l>2 and label4_cost.l>2
    for profile in profiles:
        q=I(*profile['old60_offsets'][30]);assert q.l<=0<=q.h
    assert (ap[30]+bp[30]).l<=0<=(ap[30]+bp[30]).h
    qa,qb,qc=ap[30],I(1),bp[30]
    assert qa.h<0<qc.l
    event=U.rnd((-qb-U.sqrt(qb*qb-4*qa*qc))/(2*qa));assert event.l>1
    blocker=U.rnd(ap[49]*event+bp[49]/event);assert blocker.l>0
    neighbor_gate=dict(second_edge=[30,59],I2_5_inverse_diagonal=i5diag.enc(),
        nonzero_A2_and_nonright30_minimum_cost=nonzero_cost.enc(),
        zero_A2_label4_minimum_cost=label4_cost.enc(),
        unique_remaining_event=dict(A2_profile='zero',label_to30=3,
            equation=P.enc([qa,qb,qc]),parameter=event.enc(),positive_wall49_product=blocker.enc()),
        conclusion='No acute wall forming finite terminal vertices for both[58,59] and[30,59] exists. This is an all-label single-wall gate; it does not exclude two coordinated terminal walls.')
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',
        f'scripts/{C.DOMAIN}.py',f'results/{C.DOMAIN}.json',f'scripts/{D.TAG}.py',
        f'results/{D.TAG}.json',f'scripts/{D.E.TAG}.py',f'results/{D.E.TAG}.json']
    return dict(status='PASS_COMPLETE_PI9_FINITE_ENDPOINT_PROFILES_AND_CONTINUOUS_DOMAIN',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},retained_walls=60,
        computational_root_basis=z['basis_ids'],computational_gram=[P.enc(row)for row in z['exact_gram']],
        edge_support=FRAME,edge_type='H4+A2+I2(9)',finite_endpoint_incident=FRAME+[55],
        old_finite_endpoint=P.enc(endpoint),frame_PD_pivots=pd,inverse_pivots=pivots,
        I2_9_inverse_diagonal=dihedral_diag.enc(),A2_single_contact_costs=single_costs,
        A2_double_minimum_cost=double_cost.enc(),kernel_future_unit=P.enc(k0),kernel_positive_unit=P.enc(k1),
        A=P.enc(aa),B=P.enc(bb),A_products=P.enc(ap),B_products=P.enc(bp),
        profiles=profiles,finite_profile_count=7,high_label_tails=0,
        neighboring_edge_single_wall_gate=neighbor_gate,
        normal_formula='r(t)=p+sqrt(rho/2)*(t*A+B/t), A=(k0-k1)/sqrt2, B=(-k0-k1)/sqrt2, t>1.',
        endpoint_formula='x(t)=t*A-B/t has square-2, is future, and is orthogonal to the eight support walls and r(t). For i outside the support, x(t).n_i=A_i*(t-1/t)+(A_i-B_i)/t. All A_i<0 and all A_i-B_i<0 except the exact zero at55; therefore every t>1 is an actual endpoint relative to all60 old halfspaces.',
        orientation_and_center='The other kernel orientation has r.n55=+sqrt(rho/2)*(t+1/t)>0 and is excluded. In the surviving orientation, x.n55=-t+1/t forces t>=1. At t=1 the new hyperplane contains the old finite endpoint and cannot share an acute strict interior. For t>1, r has strictly negative product with the old finite endpoint. A sufficiently small common perturbation toward the old strict center gives a common strict future interior for any finite mutually acute family satisfying all old gates.',
        tail_proof='At the saved dyadic threshold every non-support old product is below-2 and its derivative is negative. For a<0, the derivative of c+a*t+b/t is a-b/t^2: if b>=0 it is always negative; if b<0 it decreases toward a, so negativity at the threshold persists. Thus each profile has a genuine all-dotted tail for every real t at least its threshold. The seven displayed normals are separate witnesses, not asserted mutually compatible.',
        completeness='Inside a finite corner H4 has no nonright attachment. A nonzero contact to I2(9) costs more than2 in its positive inverse. In A2, a label>=6 single contact costs at least2, while two nonzero contacts cost at least2 by inverse positivity. Exactly zero or one label3/4/5 attachment remains. The residual complement is Lorentz2 and its norm equation has precisely the two displayed orientation branches; old55 excludes one and forces t>1 in the other. All60 native Coxeter conditions remain the saved exact offset equations, not a finite-label cutoff.',
        scope='Complete finite-endpoint contact profiles and actual-endpoint continuous parameterization, with individually coherent all-dotted tails. The portions where another old contact is finite are not classified here. No simultaneous cluster, compact H3 closure, or ambient finite-volume assertion is made.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text());print('PASS saved pi9 endpoint-domain replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],[(r['nonright_A2_node'],r['label'],r['all_dotted_tail_threshold'])for r in data['profiles']])
