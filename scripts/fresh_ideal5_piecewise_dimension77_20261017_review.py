#!/usr/bin/env python3
"""Independent stdlib Fraction scalar and composition review of dimension77.

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
PRIMARY='fresh_ideal5_piecewise_dimension77_20261017'
TAG=PRIMARY+'_review'
LOCAL='fresh_ideal5_combined_payments_20261017'
CENSUS='fresh_ideal5_small_factor_payments_20261017'
BASE='fresh_ideal6_piecewise_dimension80_20261017'

def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

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
    cost=sum((budget(r)for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,h_polynomial=h,G=g,
        independent_cumulative_G=gcdf,normalized_curvature=str(curvature),
        normalized_cost=str(cost),margin=str(curvature-cost))

def budget(r):
    return F(r,2)if r<=12 else F(87,14)if r==13 else F(20,3)if r==14 else F(111,16)if r==15 else F(127,17)if r==16 else F(163,23)

def affine_feature_maxima():
    # Independent spine counts at order m>=18. Between consecutive breakpoints
    # cost/m is a+b/m, so endpoints and the limiting value5 exhaust the maximum.
    specs={'B':(-6,(8,9,10,11),13,F(160,23)),
           'C':(-15,(6,7,8,9),14,F(75,11)),
           'D':(3,(10,11,12,13),12,F(163,23))}
    rows=[]
    for kind,(constant,shifts,cap,expected)in specs.items():
        breakpoints={s+cap for s in shifts}
        points=sorted({18}|{x+e for x in breakpoints for e in(-1,0,1)if x+e>=18})
        checks=[[m,5*m+constant+sum(min(m-s,cap)for s in shifts)]for m in points]
        maximum=max([F(5)]+[F(c,m)for m,c in checks]);assert maximum==expected
        cuts=sorted({18}|breakpoints);branches=[]
        for i,lo in enumerate(cuts):
            hi=cuts[i+1]if i+1<len(cuts)else None
            sample=F(lo)+F(1,2)
            active=[sample-s<cap for s in shifts]
            branches.append(dict(kind=kind,lower_order=lo,upper_order=hi,
                slope=5+sum(active),intercept=constant+sum(-s if a else cap for s,a in zip(shifts,active))))
        rows.append(dict(kind=kind,constant=constant,shifts=list(shifts),window_cap=cap,
            exact_breakpoint_checks=checks,linear_branches=branches,
            limiting_ratio='5',maximum_ratio=str(maximum)))
    return rows

def anchors(total, allowed=(12,14,16,17)):
    if len(allowed)==1:
        if total%allowed[0]==0:
            yield [allowed[0]]*(total//allowed[0])
        return
    r=allowed[0]
    for count in range(total//r+1):
        for tail in anchors(total-count*r,allowed[1:]):
            yield [r]*count+tail


def arithmetic(compare_primary=True):
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())if compare_primary else None
    profiles=[];minima=[]
    for d in range(77,80):
        shapes=set()
        for high in(0,*range(18,d+1)):
            for low in range(12):
                for thirteen in(0,1):
                    for fifteen in(0,1):
                        left=d-high-low-13*thirteen-15*fifteen
                        if left<0:continue
                        exceptional=([high]if high else[])+([low]if low else[])+[13]*thirteen+[15]*fifteen
                        for rest in anchors(left):shapes.add(tuple(sorted(exceptional+rest)))
        group=[scalar(rs)for rs in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>=9*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-9*n),strict_global_slack=5*(n-1)))
    assert len(profiles)==3722
    assert [m['profiles']for m in minima]==[1173,1239,1310]
    assert budget(12)+budget(14)-2*budget(13)==F(5,21)>0
    assert budget(14)+budget(16)-2*budget(15)==F(107,408)>0
    assert min(F(r['margin'])for r in profiles)==F(409611,2171546)
    boundary=scalar((12,16,16,16,16))
    assert F(boundary['margin'])==F(-16664,1085773)
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==3722
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n77']
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        independent_affine_feature_maxima=affine_feature_maxima(),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])


def interfaces(inputs,arithmetic_data):
    primary,local,census=inputs[PRIMARY],inputs[LOCAL],inputs[CENSUS]
    assert primary['dimension_bound']==77 and primary['inherited_dimension_bound']==80
    assert primary['excluded_dimensions']==[78,79,80]
    assert primary['unrestricted_finite_volume']is True and primary['arithmetic_hypothesis']is False
    assert primary['cusp_restriction']is None
    assert inputs[BASE]['dimension_bound']==inputs[BASE+'_review']['dimension_bound']==80
    assert local['finite_vertex_budget']==primary['finite_vertex_budget']=='9n'
    assert local['finite_budget_dimension_maximum']==80
    assert local['triangle_minimum']=='2'and local['quadrilateral_minimum']=='1'
    expected={'1_to_12':'r/2','13':'87/14','14':'20/3','15':'111/16','16':'127/17','17_and_larger':'163/23'}
    assert primary['ideal_rank_envelope']==local['ideal_rank_envelope']==expected
    assert local['low_core_maximum_factor_order']==12 and local['inherited_high_overlap_cases']==949
    assert len(census['patterns'])==571
    assert local['small_order_pair_sets']==census['records']
    costs={};pair_checks=[]
    for row in census['records']:
        g=row['target'];m=len(g)
        adj=[{j for j,x in enumerate(gr)if x and j!=i}for i,gr in enumerate(g)]
        distances=[]
        for start in range(m):
            dist={start:0};frontier=[start]
            for u in frontier:
                for v in adj[u]:
                    if v not in dist:dist[v]=dist[u]+1;frontier.append(v)
            assert len(dist)==m
            distances.append(dist)
        features={i for i in range(m)if len(adj[i])>=3 or any(g[i][j]>=4 for j in adj[i])}
        near={i for i in range(m)if any(distances[i][j]<=5 for j in features)}
        elementary={(i,j)for i in range(m)for j in range(i+1,m)
            if distances[i][j]<=5 or(distances[i][j]<=8 and(i in near or j in near))}
        paid=set(map(tuple,row['paid_pairs']));assert elementary<=paid
        assert len(paid)==row['cost']and F(row['normalized_cost'])==F(len(paid),m)
        costs.setdefault(str(m),{})[row['kind']]=len(paid)
        dropped=set(map(tuple,row['dropped_pairs']))
        if row['kind']=='D'and m>=16:
            assert min(distances[i][j]for i,j in dropped)==m-5>9
        pair_checks.append(dict(order=m,kind=row['kind'],paid_pairs=len(paid),
            elementary_triangle_pairs=len(elementary),all_retained=True))
    assert costs==local['small_order_costs']
    for m in range(14,18):assert F(max(costs[str(m)].values()),m)==budget(m-1)
    branches=[branch for family in arithmetic_data['independent_affine_feature_maxima']for branch in family['linear_branches']]
    keys=('kind','lower_order','upper_order','slope','intercept')
    saved=[{k:r[k]for k in keys}for r in census['basic_affine_cost']['all_rank_branches']]
    order=lambda r:(r['kind'],r['lower_order'])
    assert sorted(branches,key=order)==sorted(saved,key=order)
    assert local['large_order_linear_branches']==census['basic_affine_cost']['all_rank_branches']
    increments=local['finite9_increment']['rows']
    assert len(increments)==126
    for r in increments:
        assert r['count']==r['rank']-(11 if r['type']=='A_or_B'else 12)
        assert r['count']<=r['rank']and len(r['added_pairs'])==r['count']
    assert primary['effective_compact_facet_ambient_bound']==77
    assert primary['unchanged_compact_facet_ambient_bound']==79
    return dict(small_pair_checks=pair_checks,independent_linear_branches=branches,
        finite_increment_rows=126,inherited_high_overlap_cases=949,preferred_patterns=571,
        no_short9_reduction='Deleted order16/17 D pairs have distances at least11/12; their context-only replacements are used after the complete no-short9 classification.',
        finite_coefficient_at_n78='equal to9n; contradiction retains global positive385',
        scope='Independent scalar and all-rank branch reconstruction, small pair-set interface checks and complete written composition review; no duplicate571 embedding or949 core census.')


def run():
    stems=(PRIMARY,LOCAL,CENSUS,BASE,BASE+'_review')
    inputs={stem:json.loads((ROOT/f'results/{stem}.json').read_text())for stem in stems}
    for stem,doc in inputs.items():
        assert doc['status'].startswith('PASS')
        for name,h in doc['sources'].items():assert digest(ROOT/name)==h,(stem,name)
    out=arithmetic(True)
    out['composition_review']=interfaces(inputs,out)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/f'results/{stem}.md'for stem in('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017',
        'fresh_finite8_a17_d11_coupon_budget_20261017')]
    out.update(status='PASS_INDEPENDENT_IDEAL5_DIMENSION77_REVIEW',dimension_bound=77,
        inherited_dimension_bound=80,excluded_dimensions=[78,79,80],
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        external_refereeing=False,proof_assistant_formalization=False)
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--arithmetic-only',action='store_true')
    ap.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    assert not(args.arithmetic_only and args.verify)
    if args.arithmetic_only:
        d=arithmetic(False)
        print('PASS_INDEPENDENT_CONDITIONAL_77_SCALARS',len(d['profiles']),d['minimum_margin'],d['boundary']['margin'],flush=True)
    else:
        d=run();path=ROOT/f'results/{TAG}.json'
        if args.verify:assert d==json.loads(path.read_text())
        else:path.write_text(json.dumps(d,indent=2)+'\n')
        print(d['status'],len(d['profiles']),d['minimum_margin'],d['boundary']['margin'],flush=True)
