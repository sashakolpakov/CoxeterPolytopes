#!/usr/bin/env python3
"""Universal order12 B62/C62 with a sparse whole A1/B11 half coupon.

Exact norm-two profile identities, complete native cost bounds, immediate
whole-affine roles and distinct whole finite corners. No selected census.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_order12_b_center_terminal_transfer_20261018 as H
import fresh_order11_bc_zero_transfer_20261018 as K
import fresh_order13_d_center_leaf_half_transfer_20261018 as D
import fresh_balanced_rank16_order14_ambient60_frontier_20261018 as C

ROOT=H.ROOT;L=H.L;P=D.P
TAG='fresh_order12_bc_near_terminal_half_transfer_20261018'
OLD='fresh_order12_bc_unit_coupon_transfer_20261017'
STEMS=(H.TAG,K.TAG,D.TAG,OLD,C.TAG,'fresh_order11_bcd_half_pair_transfer_20261018',
       'fresh_order13_bc_sparse_terminal_transfer_20261018',
       'fresh_order14_c_distance10_terminal_projection_20261018',
       'fresh_order15_bd_terminal_coupon_transfer_20261018',
       'fresh_order17_d5_d12_center_quarter_transfer_20261018',
       'fresh_depth10_finite7_zero_transfer_20261018',
       'fresh_finite7_ideal5_triangle_payments_20261017')


def affine(kind):
    if kind=='B':return H.affine()
    g=L.gram(12,[(i,i+1,s.sqrt(2) if i in (0,10) else 1) for i in range(11)])
    nu=s.Matrix([1/s.sqrt(2)]+[1]*10+[1/s.sqrt(2)])
    assert g*nu==s.zeros(12,1) and g.rank()==11
    L.pd(g[:11,:11]);return g,nu


def inner(mark,outer,first,second):
    t,q=s.symbols('t q',nonnegative=True)
    return [dict(name='isolated',contacts={mark:t},height=t,moment=0,c0=s.Integer(2),maximum=2,lower_t=1),
      dict(name='tail_only',contacts={mark:t,outer:q},height=t+q/s.sqrt(2),moment=0,c0=2-q*q/2,maximum=3,lower_t=0),
      dict(name='end3',contacts={mark:t,first:1,outer:q},height=t+1+q/s.sqrt(2),moment=1,c0=1-q*q/2,maximum=3,lower_t=0),
      dict(name='end4',contacts={mark:t,first:s.sqrt(2)},height=t+s.sqrt(2),moment=s.sqrt(2),c0=s.Integer(0),maximum=2,lower_t=0),
      dict(name='next3',contacts={mark:t,second:1},height=t+1,moment=2,c0=s.Integer(0),maximum=4,lower_t=0)],(t,q)


def leaf(mark,mate,center):
    t=s.symbols('t',nonnegative=True)
    return [dict(name='isolated',contacts={mark:t},height=t/2,moment=0,c0=s.Integer(2),maximum=8,lower_t=1),
      dict(name='mate3',contacts={mark:t,mate:1},height=(t+1)/2,moment=1,c0=s.Integer(1),maximum=8,lower_t=0),
      dict(name='mate4',contacts={mark:t,mate:s.sqrt(2)},height=(t+s.sqrt(2))/2,moment=s.sqrt(2),c0=s.Integer(0),maximum=4,lower_t=0),
      dict(name='center3',contacts={mark:t,center:1},height=1+t/2,moment=2,c0=s.Integer(0),maximum=4,lower_t=0)],t


def formulas(kind):
    f,nu=affine(kind);marks=(8,10) if kind=='B' else (1,10);R=10 if kind=='B' else 9
    W=[i for i in range(12) if i not in marks];J=f.extract(W,W).inv().applyfunc(s.simplify)
    va,vb=[-nu[i]*f.extract([i],W).T for i in marks]
    old=s.diag(*[nu[i] for i in marks])*f.extract(marks,marks)*s.diag(*[nu[i] for i in marks])
    assert (old-s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)).applyfunc(s.simplify)==s.Matrix([[1,-1],[-1,1]])/R
    left,(t,q)=inner(8,9,7,6) if kind=='B' else inner(1,0,2,3)
    right,_=leaf(10,11,0) if kind=='B' else inner(10,11,9,8)
    identities=[];bounds=[];u=s.symbols('u',nonnegative=True)
    for side,ps,far in [('left',left,vb),('right',right,va)]:
        for p in ps:
            x=s.Matrix([p['contacts'].get(i,0) for i in W]);h=p['height'];m=s.sympify(p['moment']);c0=p['c0']
            assert s.simplify(sum(nu[i]*z for i,z in p['contacts'].items())-h)==0
            assert s.simplify(R*(far.T*J*x)[0]-m)==0
            assert s.simplify(2-(x.T*J*x)[0]-c0-m*m/R)==0
            gap=s.expand((p['maximum']*h*h-2*m*h-c0).subs(t,u+p['lower_t']))
            if p['name']=='tail_only':gap=gap.subs(q,q+1)
            certificate=K.positive_polynomial(gap,(u,q))
            identities.append(dict(side=side,name=p['name'],contacts=p['contacts'],height=h,moment=m,c0=c0,
                cost=s.factor(2*m/h+c0/h**2),maximum=p['maximum']))
            bounds.append(dict(side=side,name=p['name'],certificate=certificate))
    crosses=[]
    for p,r in product(left,right):
        x=s.Matrix([p['contacts'].get(i,0) for i in W]);y=s.Matrix([r['contacts'].get(i,0) for i in W])
        assert s.simplify((x.T*J*y)[0]-s.sympify(p['moment'])*r['moment']/R)==0
        crosses.append(dict(left=p['name'],right=r['name'],product=s.simplify((x.T*J*y)[0])))
    principal=[]
    if kind=='B':domains=[('D10',list(range(8))+[10,11],7,6),('B11',[i for i in range(12) if i!=10],11,0)]
    else:domains=[('B10',list(range(2,12)),2,3),('B10',list(range(10)),9,8)]
    for name,ids,end,nxt in domains:
        inv=f.extract(ids,ids).inv().applyfunc(s.simplify);diag=dict(zip(ids,inv.diagonal()))
        assert diag[end]==1 and diag[nxt]==2 and all(v>2 for i,v in diag.items() if i not in (end,nxt))
        principal.append(dict(type=name,nodes=ids,inverse_diagonal=diag,endpoint=end,next_node=nxt))
    extra={}
    if kind=='C':
        assert max(p['maximum'] for p in left)==max(p['maximum'] for p in right)==4
        extra['uniform_strict_gap']=s.Rational(1,9)
    else:
        # After the two whole-affine immediate profiles, left cost is at most3.
        assert max(p['maximum'] for p in left if p['name'] not in ('end4','next3'))==3
        lcost={p['name']:s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2) for p in left}
        assert lcost['isolated'].subs(t,1)==2
        assert lcost['tail_only'].subs({t:0,q:1})==3
        assert s.simplify(lcost['tail_only'].subs({t:0,q:s.sqrt(2)}))==1
        assert lcost['end3'].subs({t:0,q:0})==3
        assert s.simplify(lcost['tail_only'].subs({t:1,q:1})-2).is_negative
        assert s.simplify(lcost['end3'].subs({t:0,q:1})-2).is_negative
        assert lcost['end3'].subs({t:1,q:0})==s.Rational(5,4)
        for p in right:
            h=p['height'];m=p['moment'];c0=p['c0']
            lower=s.sqrt(2) if p['name']=='isolated' else 1 if p['name']=='mate3' else 0
            maximum=3 if p['name']=='mate3' else 4
            K.positive_polynomial(s.expand((maximum*h*h-2*m*h-c0).subs(t,u+lower)),(u,))
        extra=dict(near_nonimmediate_maximum=3,leaf_nonextreme_maximum=4,
            leaf_extremes=[dict(profile='isolated',t=1,height=s.Rational(1,2),cost=8),dict(profile='mate3',t=0,height=s.Rational(1,2),cost=8)],
            near_profiles_cost_at_least_two=[dict(profile='isolated',t=1,q=0,height=1,cost=2),
              dict(profile='tail_only',t=0,q=1,height=1/s.sqrt(2),cost=3),dict(profile='end3',t=0,q=0,height=1,cost=3)],
            spectator_weighted_energy_maximum=1,extreme_vector_alone_energy_lower=2,
            both_nonzero_energy_strict_lower=s.Rational(9,8),nonextreme_strict_gap=s.Rational(3,10))
    return dict(gram=f.tolist(),kernel=list(nu),marks=marks,white=W,resistance=R,
        principal_inverse_checks=principal,end_identities=identities,bounds=bounds,cross_identities=crosses,
        determinant_identity=K.determinant(R),consequences=extra)


def immediate(docs):
    f,nu=affine('B');rows=[]
    for name,node,mag,kind in [('end4',7,s.sqrt(2),'B'),('next3',6,1,'D')]:
        g=2*s.eye(13);g[:12,:12]=f;g[12,node]=g[node,12]=-mag
        ids=[i for i in range(13) if i!=8];components=L.components(g,ids)
        ns=next(b for b in components if 12 in b)
        assert len(ns)==11 and 10 in ns
        gram=g.extract(ns,ns);z=gram.nullspace()[0]
        if z[0]<0:z=-z
        assert all(x>0 for x in z)
        target,_,_=L.affine(kind);maps=P.embeddings(D.labels(gram),D.labels(target));assert maps
        pairs=sorted({tuple(sorted(mp[ns.index(i)] for i in (10,12))) for mp in maps})
        current=next(r for r in docs[K.TAG]['rules'] if r['kind']==kind) if kind=='B' else next(r for r in docs['fresh_order11_bcd_half_pair_transfer_20261018']['factor_rules'] if r['kind']=='D')
        weights={tuple(r['pair']):F(r['weight']) for r in current['weighted_pairs']}
        assert all(weights[p]==(0 if kind=='B' else F(1,2)) for p in pairs)
        rows.append(dict(profile=name,type='whole_affine_'+kind+'10',nodes=ns,gram=gram.tolist(),kernel=list(z),
            marked_pairs=pairs,current_weight=str(weights[pairs[0]]),distinct_from_selected=True,
            payment='independent_BC11_zero_theorem_first' if kind=='B' else 'selected_B12_half_plus_distinct_D11_half'))
    return rows


def finite_corners():
    f,_=affine('B');rows=[]
    for name,node in [('isolated',10),('mate3',11)]:
        g=2*s.eye(13);g[:12,:12]=f;g[12,node]=g[node,12]=-1
        W=[i for i in range(12) if i not in (8,10)];ids=sorted(W+[8,12]);gram=g.extract(ids,ids)
        L.pd(gram);parts=L.components(g,ids)
        if name=='isolated':
            assert sorted(map(len,parts))==[1,11] and [12] in parts
            large=next(p for p in parts if 8 in p);adj=L.graph(g.extract(large,large));m=large.index(8)
            assert len(adj[m])==2
            terminals=[j for j in adj[m] if len(adj[j])==1 and s.simplify(g[8,large[j]]**2)==2]
            assert len(terminals)==1
            reason='whole_A1/B11_penultimate_half';payment=F(1,2)
        else:
            assert len(parts)==1 and len(parts[0])==12
            adj=L.graph(gram);assert max(map(len,adj.values()))==2
            assert L.distances(adj,ids.index(8))[ids.index(12)]==10
            assert C.base(('B',12))==8*66
            reason='whole_finite_B12_unit';payment=F(1)
        rows.append(dict(profile=name,gram=g.tolist(),white=W,paying_corner=[8,12],components=parts,
            actual_whole_components=True,reason=reason,distinct_payment=str(payment),selected_payment='1/2',
            total=str(payment+F(1,2))))
    return rows


def finite_comparison():
    rep=C.Z.Z.B.B.representative;types=[t for t in C.Z.Z.B.B.types() if t[1]<=60]
    for a in types:
        for b in types:
            inc=4 if {a,b}=={('A',1),('B',11)} else 0
            assert C.actual_cross(a,b)+inc<=C.upper_cross(rep(a),rep(b))
    assert C.actual_cross(('A',1),('B',11))==0 and C.upper_cross(('A',1),('D',11))==20
    return dict(comparisons=len(types)**2,cost_unit='1/8',new_actual_A1_B11='1/2',unchanged_A1_D11_upper='5/2',
        residual_items_unchanged=True,new_forest_certificate=False)


def run():
    docs={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    formulas_by_kind={kind:formulas(kind) for kind in ('B','C')};rules=[]
    for kind in ('B','C'):
        old=set(map(tuple,docs[H.TAG]['new_paid_pairs'] if kind=='B' else docs[OLD]['C']['retained_unit_pairs']))
        assert len(old)==63
        change={(8,10),(8,11)} if kind=='B' else {(1,10)}
        new={p:(F(1,2) if kind=='B' else F()) if p in change else F(int(p in old)) for p in combinations(range(12),2)}
        assert change<=old and sum(new.values())==62
        f,_=affine(kind);target=D.labels(f);adj=L.graph(f);autos=P.embeddings(target,target)
        assert all({tuple(sorted(mp[i] for i in p)) for p in change}==change for mp in autos)
        assert all(L.distances(adj,i)[j]==9 for i,j in change)
        assert all(new[p]==int(p in old) for p in new if L.distances(adj,p[0])[p[1]]<=8)
        rules.append(dict(kind=kind,target=target,changed_pairs=sorted(change),old_cost=63,new_cost=62,
            old_weight='1',new_weight='1/2' if kind=='B' else '0',
            weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(new.items())],automorphisms=autos,
            all_distance_at_most8_weights_unchanged=True))
    # Existing A60/D62 costs are unchanged.
    assert docs[OLD]['new_order12_costs']==dict(A=60,B=64,C=63,D=62)
    paths={Path(__file__),Path(H.__file__),Path(K.__file__),Path(D.__file__),Path(C.__file__),Path(L.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{t}.{ext}' for t in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return K.enc(dict(status='PASS_ORDER12_B62_C62_NEAR_TERMINAL_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        formulas=formulas_by_kind,immediate_affine_roles=immediate(docs),whole_finite_corners=finite_corners(),rules=rules,
        new_order12_family_costs=dict(A=60,B=62,C=62,D=62),new_order12_maximum='62',new_rank11_envelope='31/6',
        new_coupon=dict(small='whole_A1',small_role='only_node',large='whole_B11',
            large_role='penultimate_node_adjacent_to_unique_doubled_terminal',eligible_marked_pairs=1,
            old_weight='0',new_weight='1/2',component_pair_cost_increase='1/2'),
        finite_pointwise_comparison=finite_comparison(),triangle_minimum=2,quadrilateral_minimum=1,
        no_selected_census_used=True,new_dimension_bound=None,
        scope='Complete universal native end profiles and arbitrary-spectator Cauchy argument. C strict gap gives a separate short1/2 unit; B retained half pairs with an existing D11 half, a whole finite B12 unit, or new whole A1/B11 penultimate half. Exact extra cost1/2 fits the bound representative allowance; no forest or global theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_order12_family_costs'],flush=True)
