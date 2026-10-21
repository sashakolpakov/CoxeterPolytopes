#!/usr/bin/env python3
"""Independent stdlib Fraction scalar and composition review of dimension71.

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
PRIMARY='fresh_ideal5_weighted_dimension71_20261017'
TAG=PRIMARY+'_review'
LOCAL='fresh_finite8_ideal5_combined_payments_20261017'
DIAG='fresh_order13_payment_support_20261017'
B13='fresh_order13_affine_b_singleton_transfer_20261017'
D12='fresh_affine_d12_leaf_pair_transfer_20261017'
CYCLE='fresh_cycle13_ideal5_transfer_20261017'
WEIGHTED='fresh_d17_partial_weight_transfer_20261017'
FINITE='fresh_finite8_d12_b10_coupon_budget_20261017'
BASE='fresh_ideal5_piecewise_dimension74_20261017'

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
    return F(r,2)if r<=11 else F(75,13)if r==12 else F(87,14)if r==13 else F(97,15)if r==14 else F(111,16)if r==15 else F(120,17)if r==16 else F(163,23)

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

def anchors(total, allowed=(11,13,14,15,16,17)):
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
    for d in range(71,74):
        shapes=set()
        for high in(0,*range(18,d+1)):
            for low in range(11):
                for twelve in(0,1):
                    left=d-high-low-12*twelve
                    if left<0:continue
                    exceptional=([high]if high else[])+([low]if low else[])+[12]*twelve
                    for rest in anchors(left):shapes.add(tuple(sorted(exceptional+rest)))
        group=[scalar(rs)for rs in sorted(shapes)]
        assert all(F(r['margin'])>0 for r in group)
        profiles+=group;n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>8*n
        m=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=m['margin'],
            minimizing_ranks=m['ranks'],finite_slack=str(beta-8*n),strict_global_slack=5*(n-1)))
    assert budget(11)+budget(13)-2*budget(12)==F(16,91)>0
    boundary=scalar((11,11,11,11,11,15))
    assert F(boundary['margin'])<0
    if primary is not None:
        by={tuple(r['ranks']):r for r in primary['reduced_profiles']}
        assert set(by)=={tuple(r['ranks'])for r in profiles}and len(by)==len(primary['reduced_profiles'])==len(profiles)
        for r in profiles:
            assert {k:v for k,v in r.items()if k!='independent_cumulative_G'}==by[tuple(r['ranks'])]
        assert minima==primary['dimension_minima']
        assert {k:v for k,v in boundary.items()if k!='independent_cumulative_G'}==primary['boundary_at_n71']
    return dict(profiles=profiles,dimension_minima=minima,boundary=boundary,
        minimum_margin=str(min(F(r['margin'])for r in profiles)),
        independent_affine_feature_maxima=affine_feature_maxima(),
        concentration=dict(duplicate_rank12_cost_increase=str(budget(11)+budget(13)-2*budget(12)),
            terminal_repeated_ranks=[11,13,14,15,16,17],optional_rank12=True,
            at_most_one_low_rank=list(range(1,11)),at_most_one_high_rank_minimum=18,
            termination='Each nontrivial concentration increases the integer sum of squared parts, bounded by d^2.'),
        formulas=['G=sum_j h_j min(j,d-j)',
            'G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1)',
            'kappa/V=3d/4-(5/2)*E(abs(sum U_i-d/2))'])


def graph(kind,m):
    g=[[0]*m for _ in range(m)]
    if kind=='A':edges=[(i,(i+1)%m,3)for i in range(m)]
    elif kind=='B':edges=[(i,i+1,4 if i==m-4 else 3)for i in range(m-3)]+[(0,m-2,3),(0,m-1,3)]
    elif kind=='C':edges=[(i,i+1,4 if i in(0,m-2)else 3)for i in range(m-1)]
    else:
        assert kind=='D'
        edges=[(i,i+1,3)for i in range(m-5)]+[(0,m-4,3),(0,m-3,3),(m-5,m-2,3),(m-5,m-1,3)]
    for i,j,label in edges:g[i][j]=g[j][i]=label
    return g

def distances(g):
    result=[]
    for root in range(len(g)):
        d={root:0};todo=[root]
        for i in todo:
            for j,label in enumerate(g[i]):
                if label and j not in d:d[j]=d[i]+1;todo.append(j)
        assert len(d)==len(g);result.append(d)
    return result

def finite_budget_check(doc):
    # Enumerate component multiplicities first, then every available dimension.
    # The producer instead loops dimension first and splits t into b,c.
    states=[]
    def visit(types,remaining,counts):
        if not types:states.append(tuple(counts));return
        rank,*rest=types
        for n in range(remaining//rank+1):visit(rest,remaining-n*rank,counts+[n])
    visit([12,13,11,17,10],75,[])
    rows=set();counts={str(n):0 for n in range(1,76)}
    for b,c,other,u,v in states:
        t=b+c+other;rank=12*b+13*c+11*other+17*u+10*v
        for n in range(max(1,rank),76):
            for a in{0,n-rank}:
                excess=a*(F(3*t,2)+u-6+b+F(c,4)+F(v,2))-8*t-22*b-18*c-35*v
                assert excess<=0
                row=(n,a,t,b,c,u,v,str(excess));assert row not in rows
                rows.add(row);counts[str(n)]+=1
    assert rows==set(map(tuple,doc['dimension_endpoints']))
    assert len(rows)==len(doc['dimension_endpoints'])==13034
    assert counts==doc['dimension_endpoint_counts']
    assert doc['finite_vertex_budget']=='8n'and doc['dimension_maximum']==75
    assert doc['component_slacks']=={'B10':'35','D12':'30','D13':'26'}
    return dict(independently_enumerated_endpoints=len(rows),dimension_counts=counts,
        maximum_excess=str(max(F(r[-1])for r in rows)),
        enumeration='Counts of D12,D13,other inherited large,A17,B10 first; then all remaining dimensions and both small-rank endpoints.')

def interfaces(docs,data):
    primary=docs[PRIMARY];base=docs[BASE];diag=docs[DIAG];local=docs[B13]
    assert primary['dimension_bound']==71 and primary['inherited_dimension_bound']==74
    assert base['dimension_bound']==docs[BASE+'_review']['dimension_bound']==74
    assert primary['excluded_dimensions']==[72,73,74]
    assert primary['unrestricted_finite_volume']and not primary['arithmetic_hypothesis']and primary['cusp_restriction']is None
    assert primary['finite_vertex_budget']=='8n'and primary['finite_budget_dimension_maximum']==75
    assert diag['eligible_patterns']==801
    bdrop={(0,10),(10,11),(10,12)};checks=[];pairs={}
    for rec in diag['records']:
        kind=rec['kind'];g=graph(kind,13);assert g==rec['target'];ds=distances(g)
        allpairs={(i,j)for i in range(13)for j in range(i+1,13)}
        if kind=='A':paid={p for p in allpairs if ds[p[0]][p[1]]<=5}
        elif kind=='B':paid=allpairs-bdrop
        elif kind=='C':paid=set(map(tuple,rec['required_union']))
        else:paid=allpairs-{(9,11),(9,12),(10,11),(10,12)}
        features={i for i in range(13)if sum(bool(x)for x in g[i])>=3 or any(x>=4 for x in g[i])}
        triangles={p for p in allpairs if ds[p[0]][p[1]]<=5 or
            (ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
        assert triangles<=paid
        short_missing={p for p in allpairs-paid if ds[p[0]][p[1]]<=8}
        assert (kind=='A'and short_missing==allpairs-paid)or(kind!='A'and not short_missing)
        affected=[r['case']for r in rec['case_embeddings']if set(map(tuple,r['marked_images']))-paid]
        immediate=set(map(tuple,rec['immediate_pairs']))-paid
        if kind=='B':assert affected==[160,166,170,1206,1744,2137]
        elif kind=='D':assert len(affected)==127
        else:assert not affected and not immediate
        expected=next(r for r in primary['order13_interfaces']['records']if r['kind']==kind)
        assert expected['retained_pairs']==[list(p)for p in sorted(paid)]
        assert expected['removed_pairs']==[list(p)for p in sorted(allpairs-paid)]
        assert expected['elementary_triangle_pairs']==[list(p)for p in sorted(triangles)]
        assert expected['selected_affected_cases']==affected
        assert expected['immediate_affected_pairs']==[list(p)for p in sorted(immediate)]
        checks.append(dict(kind=kind,cost=len(paid),selected_affected_cases=affected,
            elementary_triangle_pairs=len(triangles),universal_short_replacements=len(short_missing),
            immediate_replacements=[list(p)for p in sorted(immediate)]));pairs[kind]=paid
    assert [r['cost']for r in checks]==[65,75,73,74]
    assert pairs['A']==set(map(tuple,docs[CYCLE]['cycle_pair_rule']['retained']))
    assert docs[CYCLE]['isolated_case']['required_spectator_energy']=='10/13'
    assert docs[CYCLE]['isolated_case']['forced_Ar_rank']=='10/3'
    assert pairs['D']==set(map(tuple,docs[D12]['new_paid_pairs']))
    assert local['selected_support']['affected_cases']==[160,1206]
    conditional=local['conditional_long_orbit'];assert not conditional['removed']
    assert set(map(tuple,conditional['pairs']))=={(10,11),(10,12)}
    promotions=[]
    for row in conditional['cases']:
        if row['index']in(166,1744):assert row['finite_unit_total']==1;kind='whole finite B13'
        else:
            assert row['index']in(170,2137)and row['finite_unit_total']==0
            c=next(r for r in row['corners']if r['old_half']=='1/2')
            marked=[c['components'][j]for j in c['marked_component_indices']]
            assert {r['type']for r in marked}=={'A1','D12'}
            large=next(r for r in marked if r['type']=='D12')
            mark=next(i for i in c['marks']if i in large['nodes'])
            assert large['roles'][str(mark)]=='short_fork_leaf'
            kind='whole A1/D12 short-leaf unit upgrade'
        promotions.append(dict(case=row['index'],payment=kind))
    assert [F(r['normalized_determinant_lower_bound'])for r in local['immediate_D14_strict_gap']['rows']]==[F(1,5),F(1,6),F(1,6)]
    weighted=docs[WEIGHTED];weights={tuple(r['pair']):F(r['weight'])for r in weighted['weighted_pairs']}
    assert len(weights)==136 and set(weights.values())=={F(0),F(1,2),F(1)}
    assert sum(weights.values())==F(weighted['new_cost'])==120
    assert sum(x==F(1,2)for x in weights.values())==12 and weights[1,11]==0
    assert not weighted['added_finite_coupons']
    for key in('triangle_pairs_retained_at_weight_one','short8_pairs_retained_at_weight_one'):
        assert all(weights[tuple(p)]==1 for p in weighted[key])
    assert len(weighted['cases'])==12
    for case in weighted['cases']:
        assert F(case['ideal_total'])==F(1,2)and F(case['finite_coupon_total'])>=F(1,2)
        assert F(case['payment_total'])>=1
    costs={k:dict(v)for k,v in base['small_order_costs'].items()}
    costs['13']={r['kind']:r['cost']for r in checks};costs['17']['D']=120
    assert costs==primary['small_order_costs']
    assert all(F(max(row.values()),int(m))==budget(int(m)-1)for m,row in costs.items())
    assert primary['ideal_rank_envelope']=={'1_to_11':'r/2','12':'75/13','13':'87/14',
        '14':'97/15','15':'111/16','16':'120/17','17_and_larger':'163/23'}
    assert budget(14)+budget(16)-2*budget(15)<0
    assert data['concentration']['terminal_repeated_ranks']==primary['concentration']['repeated_cap_ranks']
    return dict(order13=checks,promoted_B_distance11_cases=promotions,
        weighted_order17=dict(total='120',half_weight_pairs=12,sealed_contexts=12),
        small_order_costs=costs,finite_budget=finite_budget_check(docs[FINITE]),
        retained_tail=74,scope='Independent saved-interface and finite-budget reconstruction plus source/mathematical local composition read; no upstream wall/profile/embedding census is repeated.')

def main(prepare=False,compare_only=False,verify=False):
    data=arithmetic(not prepare)
    if prepare:
        path=Path('/private/tmp/dimension71_independent_scalar.json')
        path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
        print('PASS_INDEPENDENT_SCALAR_PREPARATION',len(data['profiles']),data['dimension_minima'],
            'boundary',data['boundary']['margin'],flush=True)
        return
    stems=[PRIMARY,BASE,BASE+'_review',DIAG,B13,D12,CYCLE,WEIGHTED,FINITE,LOCAL]
    docs={t:json.loads((ROOT/f'results/{t}.json').read_text())for t in stems}
    data['composition_interfaces']=interfaces(docs,data)
    if compare_only:
        print('PASS_INDEPENDENT_PRIMARY_COMPARISON',len(data['profiles']),
            'finite budget endpoints',data['composition_interfaces']['finite_budget']['independently_enumerated_endpoints'],flush=True)
        return
    paths=[Path(__file__)]
    for t in stems:
        assert docs[t]['status'].startswith('PASS'),t
        for p,h in docs[t]['sources'].items():assert digest(ROOT/p)==h,(t,p)
        paths.extend(ROOT/f'{folder}/{t}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.extend(ROOT/f'results/{t}.md'for t in('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017'))
    data.update(status='PASS_INDEPENDENT_DIMENSION71_SCALARS_AND_COMPOSITION',dimension_bound=71,
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        scope='Standard-library Fraction reconstruction without producer arithmetic imports: all4913 scalar profiles,13034 finite-budget endpoints, literal order13/weighted17 interfaces and inherited all-rank feature envelopes. Written local proofs are reviewed; their geometric censuses remain frozen inputs.')
    data=json.loads(json.dumps(data));path=ROOT/f'results/{TAG}.json'
    if verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
    print(data['status'],len(data['profiles']),'profiles;13034 budget endpoints; complete interfaces',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--compare-only',action='store_true')
    p.add_argument('--verify',action='store_true');a=p.parse_args();assert sum((a.prepare,a.compare_only,a.verify))<=1
    main(a.prepare,a.compare_only,a.verify)
