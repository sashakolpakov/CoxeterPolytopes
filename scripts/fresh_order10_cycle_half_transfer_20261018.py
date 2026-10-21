#!/usr/bin/env python3
"""Cycle10 antipodal half transfer with explicit whole A1--A4/A9 coupons."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_order10_b_terminal_half_transfer_20261018 as B
import fresh_order11_bcd_half_pair_transfer_20261018 as K
import fresh_cycle11_half_pair_transfer_20261018 as H

ROOT=B.ROOT;TAG='fresh_order10_cycle_half_transfer_20261018'
STEMS=(B.TAG,K.TAG,H.TAG,'fresh_spectator_energy_below_one_20261017',
       'fresh_spectator_five_quarters_transfer_20261017','fresh_a15_d15_dimension58_20261018')
N,A,D=10,0,5


def cycle():return K.gram(N,[(i,(i+1)%N,1)for i in range(N)])


def profiles():
    f=cycle();assert f*s.ones(N,1)==s.zeros(N,1)and f.rank()==9
    deleted=f[1:,1:];K.pd(deleted);J=deleted.inv();assert all(v.is_positive for v in J)
    allowed=[i+1 for i in range(9)if J[i,i]<=2];assert allowed==[1,2,8,9]
    rows=[]
    for side in(0,1):
        for name,pos,q in(('end3',1,1),('end4',1,s.sqrt(2)),('next3',2,1)):
            contact=pos if side==0 else N-pos
            g=2*s.eye(N);g[:9,:9]=deleted;g[9,contact-1]=g[contact-1,9]=-q
            adj=K.graph(g);distance=K.distances(adj,9)[D-1];assert distance<=5
            leaf=8 if side==0 else 0;near=next(iter(adj[leaf]));assert len(adj[leaf])==1 and len(adj[near])==2
            rows.append(dict(profile=f'{name}_{side}',contact=contact,magnitude=str(q),marked_distance=distance,
                fixed_leaf=leaf+1,fixed_degree_two_neighbor=near+1,positive_pivots=K.pd(g)))
    both=2*s.eye(N);both[:9,:9]=deleted
    for i in(0,8):both[9,i]=both[i,9]=-1
    assert both*s.ones(N,1)==s.zeros(N,1)and both.rank()==9
    assert K.distances(K.graph(both),9)[D-1]==5
    # Every allowed two-contact vector except the two unit endpoints is excluded.
    pairs=[]
    for i,j in combinations(allowed,2):
        v=s.simplify(J[i-1,i-1]+J[j-1,j-1]+2*J[i-1,j-1])
        assert (v==2)==((i,j)==(1,9))and v>=2
        pairs.append(dict(nodes=[i,j],unit_projection=str(v)))
    assert 3*s.Rational(9,10)>2 and 2*s.Rational(8,5)>2
    assert s.Rational(9,10)*s.Rational(5,2)>2
    return dict(gram=[[str(v)for v in row]for row in f.tolist()],deleted_inverse_diagonal=list(map(str,J.diagonal())),
        allowed_nodes=allowed,two_contact_bounds=pairs,single_profiles=rows,both_end_affine_marked_distance=5,
        both_end_positive_kernel=[1]*10)


def heights():
    alpha=s.Rational(2,5);h,k,p,q,r=s.symbols('h k p q r',nonzero=True)
    core=s.Matrix([[alpha,-alpha,-h,0],[-alpha,alpha,0,-k],[-h,0,2-p,-q],[0,-k,-q,2-r]])
    formula=1-2*alpha/h**2-2*alpha/k**2+alpha*(p/h**2-2*q/(h*k)+r/k**2)
    assert s.factor(core.det()/(h*h*k*k)-formula)==0
    energy=s.factor(2+2/k**2-1/alpha)
    assert 1-2*alpha==s.Rational(1,5)and energy.subs(k,1)==s.Rational(3,2)
    assert s.simplify(energy.subs(k,s.sqrt(2)))==s.Rational(1,2)
    assert s.simplify(s.Rational(1,2)-energy)==1-2/k**2
    return dict(alpha=str(alpha),determinant_identity=str(formula),both_heights_at_least_sqrt2_gap='1/5',
        one_height_one_energy=str(energy),unit_heights_energy='3/2',sqrt2_boundary_energy='1/2',
        greater_sqrt2_forces_unit_height_support_empty=True)


def finite_coupons():
    f=cycle();rows=[]
    for r in range(4):
        c=10;g=2*s.eye(11+r);g[:10,:10]=f;g[c,A]=g[A,c]=-1
        if r:
            for i,j in[(c,11)]+[(i,i+1)for i in range(11,10+r)]:g[i,j]=g[j,i]=-1
        ids=[i for i in range(g.rows)if i!=A];K.pd(g.extract(ids,ids));parts=K.components(g,ids)
        cp=next(z for z in parts if c in z);bp=next(z for z in parts if D in z)
        assert len(parts)==2 and len(cp)==r+1 and bp==list(range(1,10))
        adj=K.graph(g.extract(bp,bp));ds=K.distances(adj,bp.index(D))
        assert sorted(ds[i]for i in adj if len(adj[i])==1)==[4,4]
        assert len(K.graph(g)[c]&set(cp))==(0 if r==0 else 1)
        rows.append(dict(spectator_rank=r,whole_components=parts,small_rank=r+1,small_mark=c,
            large_mark=D,large_role='unique middle of whole A9',half_weight='1/2',
            component_pair_increment='1/2'if r==0 else'1'))
    return rows


def mixed_corners():
    f=cycle();rows=[]
    for name,n,edges,contacts in(
        ('A7',7,[(i,i+1)for i in range(6)],(0,6)),
        ('D6',6,[(i,i+1)for i in range(3)]+[(3,4),(3,5)],(0,4))):
        for orientation in range(1 if name=='A7'else 2):
            g=2*s.eye(12+n);g[:10,:10]=f;c,d=10,11
            for i,j in[(c,A),(d,D)]+[(12+i,12+j)for i,j in edges]:g[i,j]=g[j,i]=-1
            cp=contacts if not orientation else contacts[::-1]
            for black,node in zip((c,d),cp):g[black,12+node]=g[12+node,black]=-1
            W=[i for i in range(g.rows)if i not in(A,D,c,d)]
            edge_sigs=[K.inertia(g.extract(W+[black],W+[black]))for black in(A,D,c,d)]
            assert all(neg==zero==0 for pos,neg,zero in edge_sigs)
            ids=W+[c,d];K.pd(g.extract(ids,ids));parts=K.components(g,ids);whole=next(z for z in parts if c in z)
            assert d in whole;adj=K.graph(g.extract(whole,whole));dist=K.distances(adj,whole.index(c))[whole.index(d)]
            assert(len(whole),dist)==((9,8)if name=='A7'else(8,6))
            assert len(adj[whole.index(c)])==len(adj[whole.index(d)])==1
            assert K.inertia(g)==[g.rows-2,1,1]
            rows.append(dict(spectator=name,orientation=orientation,edge_signatures=edge_sigs,
                whole_marked_component=whole,marked_distance=dist,payment='1',
                finite_type='A9'if name=='A7'else'E8',gram=[[str(v)for v in row]for row in g.tolist()]))
    return rows


def actual_cross(a,b):
    extra=0
    for r in range(1,5):
        if frozenset((a,b))==frozenset((('A',r),('A',9))):extra=8 if r==1 else 16
    return B.actual_cross(a,b)+extra


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,doc in docs.items():
        assert doc['status'].startswith('PASS')
        for path,digest in doc['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,(tag,path)
    pp=profiles();hh=heights();cc=finite_coupons();mm=mixed_corners()
    spectrum=K.mixed_spectrum();assert spectrum==docs[K.TAG]['three_halves_spectator']
    weights=[dict(pair=[i,j],distance=min(j-i,N-j+i),weight='1/2'if min(j-i,N-j+i)==5 else'1')for i,j in combinations(range(N),2)]
    assert sum(F(row['weight'])for row in weights)==F(85,2)and sum(row['weight']=='1/2'for row in weights)==5
    allowed_a={(A+z)%N for z in(-2,-1,1,2)};allowed_b={(D+z)%N for z in(-2,-1,1,2)}
    assert not allowed_a&allowed_b and D not in allowed_a and A not in allowed_b
    engine=B.C.C.C.P.Z.B.B;types=engine.types();rep=engine.representative
    for a in types:
        for b in types:assert B.upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    costs={str(r):dict(previous=str(F(B.actual_cross(('A',r),('A',9)),16)),new=str(F(actual_cross(('A',r),('A',9)),16)))for r in range(1,5)}
    assert all(z['previous']=='0'for z in costs.values())
    paths={Path(__file__),Path(B.__file__),Path(K.__file__),Path(H.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.add(ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md')
    return dict(status='PASS_ORDER10_CYCLE_ANTIPODAL_HALF_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        profiles=pp,heights=hh,whole_coupon_corners=cc,mixed_spectrum=spectrum,mixed_whole_corners=mm,
        triangle=dict(allowed_after_delete_a=sorted(allowed_a),allowed_after_delete_b=sorted(allowed_b),common_allowed=[],affected_triangle_impossible=True),
        weighted_pairs=weights,old_cost='45',new_cost='85/2',new_component_pair_costs=costs,
        representative_increments_beyond_B10='0',simultaneous_type_pairs=len(types)**2,
        new_order10_costs=dict(A='85/2',B='44',C='44',D='45'),order10_envelope_unchanged='9/2',
        triangle_minimum=2,quadrilateral_minimum=1,new_forest_count=False,new_dimension_bound=None,
        priority='Handle any actual cycle10 antipodal half independently; its outputs have weight at least1/2 at a distinct corner. C10/B10 keep their independent shorter output proofs. Otherwise inherit the completed baseline pointwise.',
        scope='Universal cycle10 local half theorem with explicit whole-component coupons and unchanged all-family envelope. No new forest/dimension claim or low9 concentration improvement.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    result=json.loads(json.dumps(run()))
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_cost'],result['new_component_pair_costs'],flush=True)
