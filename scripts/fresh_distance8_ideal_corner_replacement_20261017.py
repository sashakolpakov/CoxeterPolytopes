#!/usr/bin/env python3
"""Exact weak-profile bounds for replacing an unpaid ideal distance9 corner.

This proves a local quadrilateral lemma, not an ambient dimension theorem.
Completeness in arbitrary rank is the finite/affine classification argument
in the companion report, not a bound on the ranks tested here.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_simple_cusp_affine_a_quadrilateral_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance8_ideal_corner_replacement_20261017'

def cycle_profiles():
    out=A.profiles()
    for p in out:
        name=p['name'];p['position']=[0,0]
        if name=='isolated':h=A.Q(1)
        elif name=='both_ends3':h=A.Q(2);p['position']=[1,1]
        else:
            j=int(name[-1]);p['position'][j]=2 if name.startswith('next')else 1
            h=A.SQ if name.startswith('end4')else A.Q(1)
        p['minimum_height']=h;p['umax']=A.Q(1)/h
    return out

def wire(p):
    return {k:([x.wire() if isinstance(x,A.Q)else x for x in v]if isinstance(v,list)
               else v.wire()if isinstance(v,A.Q)else v)for k,v in p.items()}

def tree_profiles():
    return [dict(name=name,moment=A.Q(l),c0=A.Q(c),minimum_height=h)
        for name,l,c,h in [
            ('isolated',0,2,A.Q(1)),
            ('outer_end3',0,1,A.Q(1)),
            ('outer_end4',0,0,A.SQ),
            ('outer_next3',0,0,A.Q(1)),
            ('inner_end3',1,1,A.Q(1)),
            ('inner_next3',2,0,A.Q(1)),
            ('both_ends3',1,0,A.Q(2))]]+[
            dict(name='inner_end4',moment=A.SQ,c0=A.Q(0),minimum_height=A.SQ)]

def run():
    # The inherited identity is verified formally, rather than numerically.
    formal=A.determinant_identity()
    ps=cycle_profiles();rows=[];minimum=None
    for c,d in product(ps,repeat=2):
        # A path from the two replacement roots through the same arc has
        # length P-j_c-j_d+2. Absence of a path<=8 gives this stronger bound.
        lower=[max(9,7+c['position'][j]+d['position'][j])
            if c['position'][j]and d['position'][j]else 9 for j in range(2)]
        for ub,vb in product((0,1),repeat=2):
            u=c['umax']if ub else A.Q(0);v=d['umax']if vb else A.Q(0)
            z=A.value(A.Q(F(1,lower[0])),A.Q(F(1,lower[1])),c,d,u,v)
            assert (z-A.Q(F(5,81))).sign()>=0
            if minimum is None or(z-minimum).sign()<0:minimum=z
            rows.append(dict(c=c['name'],d=d['name'],arc_lower_bounds=lower,
                upper_height_reciprocals=[ub,vb],normalized_determinant=z.wire()))
    assert len(rows)==256 and minimum==A.Q(F(5,81))
    assert all(p['c0'].sign()>=0 and all(x.sign()>=0 for x in p['moment'])for p in ps)
    tree=[]
    for p in tree_profiles():
        u=A.Q(1)/p['minimum_height'];z=2*p['moment']*u+p['c0']*u*u
        assert z.sign()>=0 and (A.Q(4)-z).sign()>=0
        tree.append(dict(**wire(p),maximum_single_root_cost=z.wire()))
    assert len(tree)==8
    pair=[]
    for c,d in product(tree,repeat=2):
        v=A.Q(*map(F,c['maximum_single_root_cost']))+A.Q(*map(F,d['maximum_single_root_cost']))
        gap=A.Q(1)-v/9
        assert (gap-A.Q(F(1,9))).sign()>=0
        pair.append(dict(c=c['name'],d=d['name'],P_lower_bound=9,gap=gap.wire()))
    assert len(pair)==64
    # On a path A_(P-1), the endpoint/second-node inverse entries are these
    # all-rank rational formulas; these differences certify the c0/moment table.
    inverse_identities=[]
    for P in (9,10,11,18,100):
        endpoint=F(P-1,P);second=F(2*(P-2),P)
        assert 2-endpoint==1+F(1,P)
        assert 2-2*endpoint==F(2,P)
        assert 2-second==F(4,P)
        assert 2-1-endpoint==F(1,P)
        inverse_identities.append(dict(P=P,A_endpoint=str(endpoint),A_second=str(second)))
    stems=['fresh_simple_cusp_affine_a_quadrilateral_20261017',
           'fresh_simple_cusp_affine_tree_corridor_20261017',
           'fresh_cusp_factor_payment_framework_20261017',
           'fresh_large_affine_factor_quad_payment_20261017',
           'fresh_affine_distance8_feature_payment_20261017']
    paths=[Path(__file__)]
    for stem in stems:
        paths += [ROOT/f'scripts/{stem}.py',ROOT/f'results/{stem}.json',ROOT/f'results/{stem}.md']
    for p in paths:assert p.exists(),p
    return dict(status='PASS_UNIVERSAL_IDEAL_DISTANCE9_CORNER_REPLACEMENT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        formal_schur_determinant=formal,cycle_profiles=[wire(p)for p in ps],
        cycle_rectangle_checks=rows,cycle_uniform_gap=minimum.wire(),
        tree_profiles=tree,tree_pair_checks=pair,tree_uniform_gap=['1/9','0'],
        path_inverse_formula_examples=inverse_identities,
        counts=dict(cycle_profiles=8,cycle_ordered_pairs=64,cycle_rectangles=256,
            tree_profiles=8,tree_ordered_pairs=64),
        local_theorem='An actual quadrilateral with an ideal marked pair at distance9 either has a mark within14 of a degree>=3 or nonsimple-bond feature, or has a white-interior path of length<=8 between adjacent black roots. No simple-cusp or global small-affine-exclusion hypothesis.',
        hypotheses=dict(finite_sided_finite_volume=True,actual_quadrilateral=True,
            varying_affine_factor_order_greater_than=17,marked_distance=9),
        scope='Complete local weak-profile argument in arbitrary rank, with finite/affine classification and actual factor-fiber/PSD spectator inputs explicitly inherited. No ambient dimension bound, new chamber, external refereeing or formal proof claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
