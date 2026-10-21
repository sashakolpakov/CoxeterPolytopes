#!/usr/bin/env python3
"""Universal D12/order13 center-to-opposite-leaf half transfer, no new coupon.

Complete analytic profile domains and arbitrary-spectator arguments are in
the report. This checks exact formulas, current immediate roles and finite
whole-component representatives, not a sampled supporting-rank census.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as s
import fresh_order13_d_center_pair_half_transfer_20261018 as H
import fresh_order13_affine_b_singleton_transfer_20261017 as M

ROOT=H.ROOT; L=H.L; P=M.P
TAG='fresh_order13_d_center_leaf_half_transfer_20261018'
CURRENT='fresh_d17_fractional_ambient60_frontier_20261018'
OLD='fresh_order13_d_center_pair_zero_transfer_20261018'
BC='fresh_order11_bc_zero_transfer_20261018'
BCD='fresh_order11_bcd_half_pair_transfer_20261018'
STEMS=(CURRENT,OLD,BC,BCD,H.TAG,'fresh_affine_d14_long_pair_transfer_20261017',
       'fresh_affine_d12_leaf_pair_transfer_20261017','fresh_spectator_energy_one_boundary_20261017',
       'fresh_spectator_energy_below_one_20261017','fresh_depth10_finite7_zero_transfer_20261018',
       'fresh_order14_c_distance10_terminal_projection_20261018')
CHANGE={(0,11),(0,12),(8,9),(8,10)}


def labels(g):
    n=g.rows; a=[[0]*n for _ in range(n)]
    for i in range(n):
        assert g[i,i]==2
        for j in range(i):
            q=s.simplify(g[i,j]**2)
            assert q in (0,1,2)
            a[i][j]=a[j][i]=0 if q==0 else 3 if q==1 else 4
    return a


def formulas():
    g,nu=H.affine(); W=[i for i in range(13) if i not in (0,11)]; J=g.extract(W,W).inv()
    va,vb=[-nu[i]*g.extract([i],W).T for i in (0,11)]
    old=s.diag(nu[0],nu[11])*g.extract([0,11],[0,11])*s.diag(nu[0],nu[11])
    assert old-s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)==s.Matrix([[1,-1],[-1,1]])/10
    t,u,v=s.symbols('t u v',nonnegative=True)
    center=[dict(name='outer',contacts={0:t,9:u,10:v},c0=2-(u*u+v*v)/2,moment=0,height=t+(u+v)/2),
            dict(name='inner3',contacts={0:t,1:1,9:u,10:v},c0=1-(u*u+v*v)/2,moment=1,height=t+1+(u+v)/2),
            dict(name='inner4',contacts={0:t,1:s.sqrt(2)},c0=0,moment=s.sqrt(2),height=t+s.sqrt(2)),
            dict(name='next3',contacts={0:t,2:1},c0=0,moment=2,height=t+1)]
    leaf=[dict(name='isolated',contacts={11:t},c0=2,moment=0,height=t/2),
          dict(name='mate3',contacts={11:t,12:1},c0=1,moment=1,height=(t+1)/2),
          dict(name='mate4',contacts={11:t,12:s.sqrt(2)},c0=0,moment=s.sqrt(2),height=(t+s.sqrt(2))/2),
          dict(name='center3',contacts={11:t,8:1},c0=0,moment=2,height=1+t/2)]
    rows=[]
    for end,ps,far in [('center',center,vb),('leaf',leaf,va)]:
        for p in ps:
            x=s.Matrix([p['contacts'].get(i,0) for i in W])
            assert s.simplify(sum(nu[i]*q for i,q in p['contacts'].items())-p['height'])==0
            assert s.simplify(10*(far.T*J*x)[0]-p['moment'])==0
            assert s.simplify(2-(x.T*J*x)[0]-p['c0']-p['moment']**2/10)==0
            rows.append(dict(end=end,name=p['name'],contacts={str(i):str(q) for i,q in p['contacts'].items()},
                c0=str(p['c0']),moment=str(p['moment']),height=str(p['height']),
                cost=str(s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2))))
    crosses=[]
    for p,q in product(center,leaf):
        x=s.Matrix([p['contacts'].get(i,0) for i in W]);y=s.Matrix([q['contacts'].get(i,0) for i in W])
        assert s.simplify((x.T*J*y)[0]-p['moment']*q['moment']/10)==0
        crosses.append(dict(center=p['name'],leaf=q['name'],product=str(s.simplify((x.T*J*y)[0]))))
    h,k,ell,m,c0,d0,p,q,r=s.symbols('h k ell m c0 d0 p q r',nonzero=True); a=s.Rational(1,10)
    B=s.Matrix([[a,-a,-h+ell*a,-m*a],[-a,a,-ell*a,-k+m*a],
                [-h+ell*a,-ell*a,c0+ell*ell*a-p,-ell*m*a-q],[-m*a,-k+m*a,-ell*m*a-q,d0+m*m*a-r]])
    identity=1-(2*ell/h+c0/h**2+2*m/k+d0/k**2)/10+(p/h**2-2*q/(h*k)+r/k**2)/10
    assert s.factor(B.det()/(h*h*k*k)-identity)==0
    inv=[]
    for name,ids in [('center_remaining_D10',list(range(1,9))+[11,12]),('leaf_remaining_D12',[i for i in range(13) if i!=11])]:
        K=g.extract(ids,ids).inv();diagonal=dict(zip(ids,K.diagonal()))
        endpoint,nextnode=(1,2) if name=='center_remaining_D10' else (12,8)
        assert diagonal[endpoint]==1 and diagonal[nextnode]==2
        assert all(q>2 for i,q in diagonal.items() if i not in (endpoint,nextnode))
        inv.append(dict(type=name,nodes=ids,diagonal=list(map(str,K.diagonal()))))
    gc=[s.factor(2*z['moment']/z['height']+z['c0']/z['height']**2) for z in center]
    gl=[s.factor(2*z['moment']/z['height']+z['c0']/z['height']**2) for z in leaf]
    assert gc[0].subs({t:0,u:1,v:0})==6
    assert gc[0].subs({t:1,u:0,v:0})==gc[0].subs({t:0,u:s.sqrt(2),v:0})==2
    assert gc[0].subs({t:1,u:1,v:0})==s.Rational(2,3)
    assert gc[0].subs({t:0,u:1,v:1})==1
    assert gc[1].subs({t:0,u:0,v:0})==3
    assert gc[1].subs({t:1,u:0,v:0})==s.Rational(5,4)
    assert gc[1].subs({t:0,u:1,v:0})==s.Rational(14,9)
    assert [gl[0].subs(t,1),gl[1].subs(t,0)]==[8,8]
    assert [gl[0].subs(t,s.sqrt(2)),gl[1].subs(t,1),gl[2].subs(t,0),gl[3].subs(t,0)]==[4,3,4,4]
    assert F(1)-F(3+4,10)==F(3,10) and F(14-10,4)==1
    return dict(resistance=10,end_profiles=rows,cross_identities=crosses,determinant_identity=str(identity),
        principal_inverse_diagonals=inv,center_extreme_cost=6,center_nonextreme_maximum=3,
        leaf_extreme_cost=8,leaf_nonextreme_maximum=4,two_nonextreme_strict_gap='3/10',
        center_extreme_leaf_nonextreme_energy_maximum=0,leaf_extreme_center_nonextreme_energy_maximum=1,
        surviving_nonextreme_center_profiles=[dict(name='outer_isolated',t='1',height='1',cost='2'),
            dict(name='outer_one4',t='0',height='sqrt(2)/2',cost='2'),
            dict(name='inner3_no_outer',t='0',height='1',cost='3')],
        extreme_side_nonzero_alone_lower='2',both_nonzero_Cauchy_strict_lower='9/8',
        both_extreme_unscaled_energy='1')


def immediate(data):
    f,_=H.affine();out=[]
    for name,node,mag,kind in [('inner4',1,s.sqrt(2),'B'),('next3',2,1,'D')]:
        g=2*s.eye(14);g[:13,:13]=f;g[13,node]=g[node,13]=-mag
        ids=list(range(1,14));corn=g.extract(ids,ids);b=next(b for b in H.components(corn) if ids.index(13) in b)
        ns=[ids[i] for i in b];K=g.extract(ns,ns);assert len(ns)==11 and 11 in ns
        z=K.nullspace()[0]
        if z[0]<0:z=-z
        assert all(x>0 for x in z)
        canonical,_,_=L.affine(kind);mp=P.embeddings(labels(K),labels(canonical));assert mp
        pairs=sorted({tuple(sorted(m[ns.index(i)] for i in (11,13))) for m in mp})
        rec=next(r for r in data[BC]['rules'] if r['kind']==kind) if kind=='B' else next(r for r in data[BCD]['factor_rules'] if r['kind']=='D')
        weights={tuple(r['pair']):F(r['weight']) for r in rec['weighted_pairs']}
        assert all(weights[p]==(0 if kind=='B' else F(1,2)) for p in pairs)
        out.append(dict(profile=name,whole_factor='affine_'+kind+'10',nodes=ns,gram=[[str(v) for v in row] for row in K.tolist()],
            positive_kernel=list(map(str,z)),marked_pairs=pairs,retained_literal_weight=str(weights[pairs[0]]),
            payment='independent completed BC11 unit-alternative theorem applied first' if kind=='B' else 'original D13 half plus distinct retained D11 half'))
    return out


def literal_coupons():
    f,_=H.affine();rows=[]
    for leaf,extra in product((9,10),(0,1)):
        n=14+extra;g=2*s.eye(n);g[:13,:13]=f;g[13,leaf]=g[leaf,13]=-1
        if extra:g[13,14]=g[14,13]=-1
        ids=list(range(1,n));K=g.extract(ids,ids);assert L.inertia(K)==[len(ids),0,0]
        blocks=[[ids[i] for i in b] for b in H.components(K)];small=next(b for b in blocks if 13 in b);large=next(b for b in blocks if 11 in b)
        assert len(small)==2+extra and len(large)==10
        adj={i:{j for j in ids if i!=j and g[i,j]} for i in ids}
        assert len(adj[11])==1 and len(adj[next(iter(adj[11]))])==3
        assert len(adj[13])==1+extra
        rows.append(dict(side='center',attached_leaf=leaf,spectator_A1=bool(extra),marks=[11,13],small='whole_A'+str(2+extra),
            small_mark='middle' if extra else 'endpoint',large='whole_D10',large_mark='short_fork_leaf',existing_payment='1/2',
            components=blocks,signature=L.inertia(K),original_ideal_half='1/2',total='1'))
    for leaf in (11,12):
        g=2*s.eye(14);g[:13,:13]=f;g[13,leaf]=g[leaf,13]=-1
        ids=[i for i in range(14) if i!=11];K=g.extract(ids,ids);assert L.inertia(K)==[13,0,0]
        blocks=[[ids[i] for i in b] for b in H.components(K)];small=next(b for b in blocks if 13 in b);large=next(b for b in blocks if 0 in b)
        if leaf==11:
            assert len(small)==1 and len(large)==12 and small!=large
            adj={i:{j for j in ids if i!=j and g[i,j]} for i in ids};assert len(adj[0])==3
            kind='whole_A1/whole_D12_fork_center';weight='1/2'
        else:assert small==large and len(large)==13;kind='whole_finite_D13_unit';weight='1'
        rows.append(dict(side='leaf',attached_leaf=leaf,marks=[0,13],kind=kind,existing_payment=weight,
            components=blocks,signature=L.inertia(K),original_ideal_half='1/2',total=str(F(1,2)+F(weight))))
    assert len(rows)==6
    return rows


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,d in data.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    old={tuple(r['pair']):F(r['weight']) for r in data[OLD]['weighted_pairs']};assert sum(old.values())==73
    f,_=H.affine();target=labels(f);autos=P.embeddings(target,target)
    assert all({tuple(sorted(m[i] for i in p)) for p in CHANGE}==CHANGE for m in autos)
    for p in CHANGE:assert old[p]==1 and L.distances(L.graph(f),p[0])[p[1]]==9;old[p]=F(1,2)
    assert sum(old.values())==71 and old[0,8]==0
    assert all(w==1 for p,w in old.items() if L.distances(L.graph(f),p[0])[p[1]]<=7)
    paths={Path(__file__),Path(H.__file__),Path(L.__file__),Path(M.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PASS_UNIVERSAL_ORDER13_D_CENTER_LEAF_EXISTING_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        affine_type='D12',order=13,old_cost='73',new_cost='71',changed_pairs=sorted(CHANGE),old_weight='1',new_weight='1/2',
        weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(old.items())],target=target,automorphisms=autos,
        formulas=formulas(),immediate_affine_interfaces=immediate(data),finite_corner_representatives=literal_coupons(),
        inherited_energy_one_checks=H.energy_one_checks(),finite_weights_unchanged=True,added_finite_coupons=[],
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        actual_corner_priority='First apply the completed independent BC11 zero theorem whenever a changed B/C11 corner occurs. Otherwise each changed D13 center/leaf corner has its retained half and a distinct existing half or unit. On faces with no such corner the preceding whole-face theorem applies pointwise.',
        scope='Universal all-native center/leaf R10 profile theorem with arbitrary spectator rank, inherited complete energy-one classification and existing whole-component coupons. No selected-map cutoff, new coupon, finite forest calculation or ambient dimension bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'D',out['new_cost'],'coupon delta',out['added_finite_coupons'],flush=True)
