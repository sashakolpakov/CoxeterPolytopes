#!/usr/bin/env python3
"""Independent stdlib Fraction scalar and composition review of the candidate dimension80.

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
PRIMARY='fresh_ideal6_piecewise_dimension80_20261017'
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
    cost=sum((budget(r)for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,h_polynomial=h,G=g,
        independent_cumulative_G=gcdf,normalized_curvature=str(curvature),
        normalized_cost=str(cost),margin=str(curvature-cost))

def budget(r):
    return F(r,2)if r<=13 else F(20,3)if r==14 else F(111,16)if r==15 else F(127,17)if r==16 else F(95,13)

def affine_feature_maxima():
    # Independent spine counts at order m>=18. Between consecutive breakpoints
    # cost/m is a+b/m, so endpoints and the limiting value6 exhaust the maximum.
    specs={'B':(-14,(9,10,11),15,F(36,5)),
           'C':(-21,(7,8,9),16,F(85,12)),
           'D':(-7,(11,12,13),14,F(95,13))}
    rows=[]
    for kind,(constant,shifts,cap,expected)in specs.items():
        breakpoints={s+cap for s in shifts}
        points=sorted({18}|{x+e for x in breakpoints for e in(-1,0,1)if x+e>=18})
        checks=[[m,6*m+constant+sum(min(m-s,cap)for s in shifts)]for m in points]
        maximum=max([F(6)]+[F(c,m)for m,c in checks]);assert maximum==expected
        cuts=sorted({18}|breakpoints);branches=[]
        for i,lo in enumerate(cuts):
            hi=cuts[i+1]if i+1<len(cuts)else None
            sample=F(lo)+F(1,2)
            active=[sample-s<cap for s in shifts]
            branches.append(dict(kind=kind,lower_order=lo,upper_order=hi,
                slope=6+sum(active),intercept=constant+sum(-s if a else cap for s,a in zip(shifts,active))))
        rows.append(dict(kind=kind,constant=constant,shifts=list(shifts),window_cap=cap,
            exact_breakpoint_checks=checks,linear_branches=branches,
            limiting_ratio='6',maximum_ratio=str(maximum)))
    return rows

def arithmetic(compare_primary=True):
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())if compare_primary else None
    groups={};profiles=[];minima=[]
    for d in range(80,86):
        shapes=set()
        # Optional exceptional parts first; multiplicities16/17 next and
        # the count of13 forced by the remaining rank. No producer import.
        for high in(0,*range(18,d+1)):
            for low in range(13):
                for fourteen in(0,1):
                    for fifteen in(0,1):
                        left=d-high-low-14*fourteen-15*fifteen
                        if left<0:continue
                        for c17 in range(left//17+1):
                            rest=left-17*c17
                            for c16 in range(rest//16+1):
                                rem=rest-16*c16
                                if rem%13:continue
                                ranks=([high]if high else[])+([low]if low else[])+[14]*fourteen+[15]*fifteen+[17]*c17+[16]*c16+[13]*(rem//13)
                                shapes.add(tuple(sorted(ranks)))
        group=[scalar(r)for r in sorted(shapes)];groups[d]=shapes
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>9*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-9*n)))
    assert len(profiles)==5008
    assert [m['profiles']for m in minima]==[742,776,813,852,892,933]
    assert budget(13)+budget(15)-2*budget(14)==F(5,48)>0
    assert budget(14)+budget(16)-2*budget(15)==F(107,408)>0
    assert min(F(r['margin'])for r in profiles)==F(940675,11294304)
    boundary=scalar((1,)+(13,)*6)
    assert F(boundary['margin'])==F(-1131729,7529536)
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==5008
        for r in profiles:
            a=by[tuple(r['ranks'])]
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==a
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n80']
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        independent_affine_feature_maxima=affine_feature_maxima(),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])

def run():
    p=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert p['status']=='PASS_IDEAL6_PIECEWISE_DIMENSION80'
    assert p['dimension_bound']==80 and p['inherited_dimension_bound']==86
    assert p['unchanged_compact_facet_ambient_bound']==79
    envelope={'1_to_13':'r/2','14':'20/3','15':'111/16','16':'127/17','17_and_larger':'95/13'}
    assert p['ideal_rank_envelope']==envelope
    assert p['concentration']['repeated_cap_ranks']==[13,16,17]
    assert p['concentration']['double14_budget_increase']=='5/48'
    assert p['concentration']['double15_budget_increase']=='107/408'
    assert p['excluded_dimensions']==list(range(81,87))
    assert p['unrestricted_finite_volume']and not p['arithmetic_hypothesis']and p['cusp_restriction']is None
    assert p['finite_vertex_budget']=='9n'
    for path,h in p['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    data=arithmetic()
    local=json.loads((ROOT/'results/fresh_ideal6_combined_payments_20261017.json').read_text())
    d14=json.loads((ROOT/'results/fresh_affine_d14_long_pair_transfer_20261017.json').read_text())
    assert local['status']=='PASS_COMPLETE_FINITE9_IDEAL6_LOCAL_PAYMENT_AND_BUDGET'
    assert local['ideal_rank_envelope']==envelope and local['finite_budget_dimension_maximum']==108
    assert local['triangle_minimum']=='2'and local['quadrilateral_minimum']=='1'
    assert local['small_order_costs']=={'15':dict(A=90,B=100,C=95,D=96),'16':dict(A=96,B=107,C=103,D=111),'17':dict(A=102,B=111,C=108,D=127)}
    assert d14['status']=='PASS_UNIVERSAL_AFFINE_D14_LONG_PAIR_HALF_COUPON_TRANSFER'
    assert d14['counts']['marked_pairs']==9 and d14['counts']['equality_cases']==72
    ours=[b for r in data['independent_affine_feature_maxima']for b in r['linear_branches']]
    theirs=[{k:r[k]for k in('kind','lower_order','upper_order','slope','intercept')}for r in local['large_order_linear_branches']]
    assert sorted(ours,key=lambda r:(r['kind'],r['lower_order']))==sorted(theirs,key=lambda r:(r['kind'],r['lower_order']))
    stems=(PRIMARY,'fresh_ideal6_combined_payments_20261017','fresh_affine_d14_long_pair_transfer_20261017',
        'fresh_ideal6_small_factor_frontier_20261017','fresh_broad_finite_half_coupons_20261017',
        'fresh_ideal6_triangle_payments_20261017','fresh_ideal_distance6_mixed_replacement_20261017',
        'fresh_order16_piecewise_dimension86_20261017','fresh_order16_piecewise_dimension86_20261017_review')
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/f'results/{stem}.md'for stem in('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    for path in paths:assert path.is_file(),path
    return dict(status='PASS_INDEPENDENT_FRACTION_IDEAL6_PIECEWISE_DIMENSION80_REVIEW',
        sources={str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest()for path in paths},
        dimension_bound=80,unchanged_compact_facet_ambient_bound=79,
        counts=dict(reduced_profiles=5008,ambient_dimensions=6,curvature_formula_comparisons=5009,
            independent_all_rank_affine_linear_branches=12),
        reconstruction=data,ideal_rank_envelope=envelope,
        reduction=dict(repeated14_to_13_15_cost_increase='5/48',repeated15_to_14_16_cost_increase='107/408',
            curvature_minus_cost_direction='nonincreasing',termination='strictly increasing bounded integer sum of squared ranks'),
        scope='Independent no-project-arithmetic scalar reconstruction and source/written composition review. The local D14 profiles, sealed equality corners, complete inherited preferred-component census and actual cusp-factor/shape theorems remain explicit inputs. The12 all-rank affine cost branches are independently reconstructed from spine intervals. No duplicate local geometric census, external refereeing or formalization.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--arithmetic-only',action='store_true');args=ap.parse_args()
    if args.arithmetic_only:
        d=arithmetic(False);print('PASS_INDEPENDENT_CONDITIONAL_80_SCALARS',len(d['profiles']),d['minimum_margin'],d['boundary']['margin']);raise SystemExit(0)
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
