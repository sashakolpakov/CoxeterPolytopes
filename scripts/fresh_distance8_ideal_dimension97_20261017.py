#!/usr/bin/env python3
"""Exact scalar transfer for finite distance9 / ideal distance8-feature weights.

The universal ideal-corner replacement is a separate explicit input.  This
source never substitutes a finite search for its all-rank geometric proof.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance8_ideal_dimension97_20261017'
BASE='fresh_selective_distance9_dimension108_20261017'
AFFINE='fresh_affine_distance8_feature_payment_20261017'
REPLACE='fresh_distance8_ideal_corner_replacement_20261017'


def digest(path):
    h=sha256()
    with path.open('rb')as stream:
        for b in iter(lambda:stream.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def profile(ranks):
    d=sum(ranks);h=[1]
    for r in ranks:
        out=[0]*(len(h)+r)
        for i,a in enumerate(h):
            for j in range(r+1):out[i+j]+=a
        h=out
    v=prod(r+1 for r in ranks)
    assert len(h)==d+1 and h==h[::-1]and sum(h)==v
    g=sum(a*min(i,d-i)for i,a in enumerate(h))
    curvature=F(5*g-d*v,2*v)
    cost=sum((min(F(8),F(r,2))for r in ranks),F())
    # Separate bounded-composition moment identity checks the convolution.
    moment=0;s=len(ranks)
    for size in range(s+1):
        for part in combinations(ranks,size):
            shift=sum(r+1 for r in part);m=(d-1)//2-shift
            if m>=0:moment+=(-1)**size*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    assert curvature==F(3*d,4)-F(5*moment,2*v)
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=v,h_polynomial=h,G=g,
        normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))


def scalar_checks():
    rows=[];minima=[]
    for d in range(97,108):
        shapes=set()
        for copies in range(d//16+1):
            for low in range(16):
                high=d-16*copies-low
                if high==0 or high>=17:
                    shapes.add(tuple(sorted([16]*copies+([low]if low else[])+([high]if high else[]))))
        group=[profile(ranks)for ranks in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        rows.extend(group);n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>9*n
        smallest=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=smallest['margin'],
            minimizing_ranks=smallest['ranks'],finite_slack=str(beta-9*n)))
    assert len(rows)==957 and min(F(r['margin'])for r in rows)==F(112101,1419857)
    boundary=profile((16,)*6)
    assert F(boundary['margin'])==F(-185352,1419857)
    return dict(reduced_profiles=rows,dimension_minima=minima,boundary_at_n97=boundary)


def run(scalar_only=False):
    out=scalar_checks()
    if scalar_only:
        out.update(status='PASS_CONDITIONAL_DISTANCE8_IDEAL_SCALAR_CHECKS',dimension_bound=None,
            scope='Scalar inequalities only; the geometric ideal-corner replacement is not checked in this mode.')
        return out
    inputs={}
    for stem in(BASE,BASE+'_review',AFFINE,REPLACE,REPLACE+'_review'):
        data=json.loads((ROOT/f'results/{stem}.json').read_text())
        assert data['status'].startswith('PASS')
        for rel,h in data['sources'].items():assert digest(ROOT/rel)==h,(stem,rel)
        inputs[stem]=data
    assert inputs[BASE]['dimension_bound']==inputs[BASE+'_review']['dimension_bound']==108
    assert inputs[BASE]['compact_facet_ambient_bound']==79
    assert inputs[AFFINE]['graph_cost_lemma_complete']
    assert inputs[AFFINE]['rule']==dict(short_distance=8,extra_distance=9,feature_neighborhood_radius=14,
        feature_vertices='degree>=3 or endpoint of a Coxeter bond with label>=4 (parallel bond also counts)',
        extra_condition='at least one endpoint is within14 graph steps of the feature set',
        optional_whole_component_order_at_most=17)
    # The replacement's final explicit scope is bound with its source/report.
    assert inputs[REPLACE]['status']=='PASS_UNIVERSAL_IDEAL_DISTANCE9_CORNER_REPLACEMENT'
    assert inputs[REPLACE]['hypotheses']==dict(finite_sided_finite_volume=True,actual_quadrilateral=True,
        varying_affine_factor_order_greater_than=17,marked_distance=9)
    paths=[Path(__file__)]
    for stem in(BASE,BASE+'_review',AFFINE,REPLACE,REPLACE+'_review'):
        paths.extend(ROOT/f'{folder}/{stem}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.extend(ROOT/f'results/{stem}.md'for stem in(
        'fresh_cusp_factor_payment_framework_20261017','fresh_cusp_link_shape_curvature_20261017',
        'fresh_selective_distance10_dimension116_20261017'))
    out.update(status='PASS_FINITE9_IDEAL8_FEATURE_DIMENSION97',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=97,inherited_dimension_bound=108,excluded_dimensions=list(range(98,109)),
        unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        finite_vertex_budget='9n',ideal_normalized_budget='sum_i min(8,r_i/2)',convex_order_cap=16,
        unchanged_compact_facet_ambient_bound=79,
        payment_transfer=dict(inherited=BASE,ideal_graph_budget=AFFINE,universal_replacement=REPLACE,
            triangles='Two inherited distance<=8 corner paths remain paid.',
            quadrilaterals='An inherited paid finite corner or component of order<=17 remains paid. An inherited ideal corner at distance<=8 remains paid. A remaining ideal distance9 corner is feature-near, or the universal replacement supplies another paid distance<=8 corner.'),
        scope='Complete dimension transfer under the explicitly bound inherited face/cusp geometry, all-rank convex-order and shape inequalities, affine graph-cost proof and universal ideal-corner replacement. No new high-order census is required. External refereeing, proof-assistant formalization and a new construction are not asserted.')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--scalar-only',action='store_true');args=ap.parse_args()
    assert not(args.verify and args.scalar_only)
    out=run(args.scalar_only)
    if not args.scalar_only:
        p=ROOT/f'results/{TAG}.json'
        if args.verify:assert out==json.loads(p.read_text())
        else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['reduced_profiles']),'profiles',flush=True)
