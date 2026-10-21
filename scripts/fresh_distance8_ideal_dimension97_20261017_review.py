#!/usr/bin/env python3
"""Independent stdlib Fraction scalar and composition review of dimension97.

No producer/project arithmetic imports. Sliding-window convolution and a
bounded-composition cumulative-tail identity independently compute curvature.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_distance8_ideal_dimension97_20261017'
TAG=PRIMARY+'_review'

def choose(n,k):return comb(n,k)if n>=k>=0 else 0

def scalar(ranks):
    # Sliding-window multiplication, rather than the producer's nested sums.
    h=[1]
    for r in ranks:
        out=[];window=0
        for j in range(len(h)+r):
            if j<len(h):window+=h[j]
            if 0<=j-r-1<len(h):window-=h[j-r-1]
            out.append(window)
        h=out
    d=sum(ranks);V=prod(r+1 for r in ranks);half=d//2
    assert len(h)==d+1 and h==h[::-1]and sum(h)==V
    g=sum(a*min(j,d-j)for j,a in enumerate(h))
    # Numerator of product(1-z^(r+1)); repeated shifts are combined exactly.
    shifts={0:1}
    for r in ranks:
        nxt=dict(shifts)
        for t,a in shifts.items():nxt[t+r+1]=nxt.get(t+r+1,0)-a
        shifts={t:a for t,a in nxt.items()if a}
    s=len(ranks)
    # G=sum_{t=1}^{floor(d/2)}(V-2*CDF(t-1)); hockey-stick summation.
    gcdf=half*V-2*sum(a*choose(half-t+s,s+1)for t,a in shifts.items())
    assert g==gcdf
    curvature=F(5*g-d*V,2*V)
    absolute=F(sum(a*abs(2*j-d)for j,a in enumerate(h)),2*V)
    assert curvature==F(3*d,4)-F(5,2)*absolute
    cost=sum((min(F(8),F(r,2))for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,h_polynomial=h,G=g,
        independent_cumulative_G=gcdf,normalized_curvature=str(curvature),
        normalized_cost=str(cost),margin=str(curvature-cost))

def arithmetic():
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    groups={};profiles=[];minima=[]
    for d in range(97,108):
        shapes=set()
        # Enumerate the two exceptional parts first; the number of cap parts
        # is then forced by the sum. This differs from the producer's loops.
        for high in(0,*range(17,d+1)):
            for low in range(16):
                remainder=d-high-low
                if remainder>=0 and remainder%16==0:
                    shapes.add(tuple(sorted(([low]if low else[])+[16]*(remainder//16)+([high]if high else[]))))
        group=[scalar(r)for r in sorted(shapes)];groups[d]=shapes
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>9*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-9*n)))
    assert len(profiles)==957
    by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
    assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==957
    for r in profiles:
        a=by[tuple(r['ranks'])]
        assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==a
    assert minima==primary['dimension_minima']
    assert min(F(r['margin'])for r in profiles)==F(112101,1419857)
    boundary=scalar((16,)*6)
    assert F(boundary['margin'])==F(-185352,1419857)
    assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n97']
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])

def run():
    p=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert p['status']=='PASS_FINITE9_IDEAL8_FEATURE_DIMENSION97'and p['dimension_bound']==97
    assert p['inherited_dimension_bound']==108 and p['excluded_dimensions']==list(range(98,109))
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert p['finite_vertex_budget']=='9n'and p['ideal_normalized_budget']=='sum_i min(8,r_i/2)'and p['convex_order_cap']==16
    assert p['unchanged_compact_facet_ambient_bound']==79
    for path,h in p['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    data=arithmetic()
    stems=(PRIMARY,'fresh_distance8_ideal_corner_replacement_20261017',
        'fresh_distance8_ideal_corner_replacement_20261017_review',
        'fresh_affine_distance8_feature_payment_20261017')
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/f'results/{s}.md'for s in('fresh_selective_distance9_dimension108_20261017',
        'fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017',
        'fresh_distance14_unrestricted_dimension156_20261017')]
    for path in paths:assert path.is_file(),path
    return dict(status='PASS_INDEPENDENT_FRACTION_FINITE9_IDEAL8_DIMENSION97_REVIEW',
        sources={str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest()for path in paths},
        dimension_bound=97,unchanged_compact_facet_ambient_bound=79,
        counts=dict(reduced_profiles=957,ambient_dimensions=11,curvature_formula_comparisons=958),
        reconstruction=data,scope='Independent exact scalar reconstruction plus written mathematical/source composition review. Complete inherited dimension108, actual cusp-factor geometry, finite/affine classifications and convex-order/shape proofs remain explicit inputs. No new high-core census, external refereeing or formalization.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--arithmetic-only',action='store_true');args=ap.parse_args()
    if args.arithmetic_only:
        d=arithmetic();print('PASS_INDEPENDENT_97_SCALARS',len(d['profiles']),d['minimum_margin'],d['boundary']['margin']);raise SystemExit(0)
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
