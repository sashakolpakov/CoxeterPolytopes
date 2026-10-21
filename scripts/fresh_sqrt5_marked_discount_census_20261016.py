#!/usr/bin/env python3
"""Exhaustive marked rank21 baseline census and discount-budget diagnostics.

No face-covering claim: every discount requires fresh triangle/quad coverage.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_marked_discount_census_20261016'
BASE='results/fresh_sqrt5_marked_facet_barrier_20261016.json'
TRI='results/fresh_sqrt5_triangle_weights_20261015.json'
DELTA=F(3,19)
POLICIES=['AD_distance2','local33_distance2','long_local33_distance2',
          'adjacent3','long_adjacent3']


def colored_code(g,colors):
    alive=set(g)
    while len(alive)>2:
        alive-={i for i in alive if sum(j in alive for j in g[i])<=1}
    def visit(i,parent):
        return (colors.get(i,0),tuple(sorted((m,visit(j,i)) for j,m in g[i].items() if j!=parent)))
    if len(alive)==1:return (1,visit(next(iter(alive)),None))
    a,b=sorted(alive)
    return (2,g[a][b],tuple(sorted((visit(a,b),visit(b,a)))))


def features(name,g,a,b,mark=None):
    ds={i:C.distances(g,i) for i in g}
    distance=ds[a][b]
    mid=next((x for x in g[a] if b in g[x]),None) if distance==2 else None
    path=sorted([g[a][mid],g[b][mid]]) if mid is not None else [g[a][b]] if distance==1 else []
    branches=[v for v in g if len(g[v])>=3]
    multiple=[v for v in g if any(m>3 for m in g[v].values())]
    near=lambda nodes:min((min(ds[a][v],ds[b][v]) for v in nodes),default=None)
    return dict(distance=distance,path_labels=path,
                facet_distances='external' if mark is None else sorted([ds[mark][a],ds[mark][b]]),
                nearest_branch_distance=near(branches),nearest_multiple_bond_distance=near(multiple),
                central_degree=len(g[mid]) if mid is not None else None,
                component_rank=len(g),component_type=name)


def eligible(policy,name,g,a,b):
    d=C.distances(g,a)[b]
    if policy=='AD_distance2':return d==2 and name[0] in 'AD'
    if policy.startswith('long_') and len(g)<6:return False
    if 'distance2' in policy:
        if d!=2:return False
        mid=next(x for x in g[a] if b in g[x])
        return g[a][mid]==g[b][mid]==3
    return g[a].get(b)==3


def run():
    prior=json.loads((ROOT/BASE).read_text());tri=json.loads((ROOT/TRI).read_text())
    assert prior['status']=='PASS_EXACT_MARKED_FACET_BASELINE38_BARRIER_ABOVE_THRESHOLD35'
    models=[(n,g) for n,g in C.models() if len(g)<=21]
    assert len(models)==68
    data=[];catalog={}
    for name,g in models:
        weights={(a,b):C.base(g,a,b) for a,b in combinations(g,2)}
        total=sum(weights.values(),F(0));deficit=2*len(g)-total
        assert deficit>0 and deficit.denominator==1
        roots={}
        for e in g:roots.setdefault(C.code(g,[e]),[]).append(e)
        debit={e:sum(w for p,w in weights.items() if e in p) for e in g}
        old=next(row for row in prior['finite_components'] if row['type']==name)
        assert total==F(old['total']) and len(g)==old['rank']
        assert debit=={r['root']:F(r['debit']) for r in old['marked_root_debits']}
        row=dict(name=name,g=g,rank=len(g),total=total,deficit=int(deficit),roots=list(roots.values()),debit=debit)
        data.append(row)
        for mark in [None]+[v[0] for v in roots.values()]:
            orbits={}
            for (a,b),weight in weights.items():
                if mark in (a,b) or weight==0:continue
                colors={a:2,b:2}
                if mark is not None:colors[mark]=1
                code=colored_code(g,colors)
                if code not in orbits:
                    orbits[code]=dict(pairs=[],weight=str(weight),features=features(name,g,a,b,mark),
                                     eligible_policies=[p for p in POLICIES if eligible(p,name,g,a,b)])
                record=orbits[code];assert record['weight']==str(weight)
                assert record['features']==features(name,g,a,b,mark)
                record['pairs'].append([a,b])
            key=name+('|unmarked' if mark is None else f'|mark{mark}')
            catalog[key]=dict(component=name,facet_root=mark,rank=len(g),
                              orbits=[dict(orbit=i,**r) for i,r in enumerate(orbits.values())])
    count=[0]*22;count[0]=1
    for row in data:
        for rank in range(row['rank'],22):count[rank]+=count[rank-row['rank']]
    assert count[21]==311574
    full_marked_count=sum(len(row['roots'])*count[21-row['rank']] for row in data)

    @lru_cache(None)
    def residuals(rank,budget,start=0):
        if budget<0:return ()
        if rank==0:return ((),)
        ans=[]
        for j in range(start,len(data)):
            row=data[j]
            if row['rank']<=rank and row['deficit']<=budget:
                ans.extend((j,)+tail for tail in residuals(rank-row['rank'],budget-row['deficit'],j))
        return tuple(ans)

    high=[]
    for row in data:
        for orbit in row['roots']:
            mark=orbit[0];md=row['deficit']+row['debit'][mark]
            for tail in residuals(21-row['rank'],int(7-md)):
                score=row['total']-row['debit'][mark]+sum(data[j]['total'] for j in tail)
                assert score>=35
                refs=Counter([row['name']+f'|mark{mark}']+[data[j]['name']+'|unmarked' for j in tail])
                policy_counts={p:sum(mult*len(o['pairs']) for key,mult in refs.items()
                                    for o in catalog[key]['orbits'] if p in o['eligible_policies']) for p in POLICIES}
                high.append(dict(index=len(high),marked_component=row['name'],facet_root=mark,
                                 equivalent_facet_roots=orbit,unmarked_components=[data[j]['name'] for j in tail],
                                 baseline=str(score),rank_deficit=str(42-score),corner_catalog_multiplicities=dict(refs),
                                 eligible_corner_counts=policy_counts,
                                 discounted_costs={p:str(score-DELTA*h) for p,h in policy_counts.items()}))
    histogram=Counter(r['baseline'] for r in high)
    assert histogram=={'35':295,'36':68,'37':9,'38':2} and len(high)==374
    policy_results=[]
    for policy in POLICIES:
        costs={}
        for row in data:
            g=row['g'];weights={(a,b):C.base(g,a,b)-DELTA*eligible(policy,row['name'],g,a,b) for a,b in combinations(g,2)}
            assert all(v>=0 for v in weights.values())
            costs[row['name']]=(sum(weights.values()),
                {e:sum(v for pair,v in weights.items() if e in pair) for e in g})
        dp=[(F(0),{()})]
        for rank in range(1,22):
            best=None;witnesses=set()
            for row in data:
                if row['rank']>rank:continue
                score,forests=dp[rank-row['rank']]
                score+=costs[row['name']][0]
                choices={tuple(sorted(t+(row['name'],))) for t in forests}
                if best is None or score>best:best,witnesses=score,choices
                elif score==best:witnesses|=choices
            dp.append((best,witnesses))
        marked=[]
        for row in data:
            total,debits=costs[row['name']]
            for roots in row['roots']:
                e=roots[0];rest,forests=dp[21-row['rank']]
                marked.append(dict(component=row['name'],facet_root=e,maximum=str(total-debits[e]+rest),
                                   residual_forests=sorted(forests)))
        maximum=max(F(r['maximum']) for r in marked)
        assert maximum==max(F(r['discounted_costs'][policy]) for r in high)
        unchangeable=[r['index'] for r in high if F(r['baseline'])>35 and r['eligible_corner_counts'][policy]==0]
        ratios=[(F(r['baseline'])-35)/r['eligible_corner_counts'][policy] for r in high
                if F(r['baseline'])>35 and r['eligible_corner_counts'][policy]>0]
        policy_results.append(dict(policy=policy,discount=str(DELTA),exact_marked_maximum=str(maximum),
             extremizers=[r for r in marked if F(r['maximum'])==maximum],
             fixed_obstruction_high_type_indices=unchangeable,
             minimum_discount_if_no_fixed_obstructions=str(max(ratios)),
             all_rank_unmarked_dp=[dict(rank=i,maximum=str(v),witnesses=sorted(fs)) for i,(v,fs) in enumerate(dp)],
             marked_component_maxima=marked))
    triangle_rows=[]
    for index,row in enumerate(tri['markings']):
        labels=dict(zip(combinations(range(row['order']),2),row['labels']))
        totals={p:F(0) for p in POLICIES};base=F(0);corners=[]
        for saved in row['corners']:
            a,b=saved['marked'];pcs=C.pieces(C.graph(row['white']+[a,b],labels));g=next(c for c in pcs if a in c)
            if b in g:
                name=next(n for n,h in models if C.code(g)==C.code(h))
                w=C.base(g,a,b);changed=[p for p in POLICIES if eligible(p,name,g,a,b)]
            else:name=None;w=F(0);changed=[]
            base+=w
            for p in POLICIES:totals[p]+=w-DELTA*(p in changed)
            corners.append(dict(pair=[a,b],component=name,baseline=str(w),changed_policies=changed))
        assert base>=2
        triangle_rows.append(dict(index=index,order=row['order'],white=row['white'],corners=corners,
                                  baseline_sum=str(base),discounted_sums={p:str(v) for p,v in totals.items()}))
    for policy in policy_results:
        policy['violated_core_triangle_indices']=[r['index'] for r in triangle_rows if F(r['discounted_sums'][policy['policy']])<2]
    sources=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',BASE,TRI]
    return dict(status='PASS_EXHAUSTIVE_374_HIGH_MARKED_TYPES_AND_DISCOUNT_BUDGETS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},dimension=21,facet_dimension=20,
        threshold='35',baseline_inherited_from=BASE,finite_components=68,
        all_unmarked_forest_counts=count,rank21_unmarked_forests=count[21],all_marked_types=full_marked_count,
        high_marked_type_count=len(high),high_histogram=dict(sorted(histogram.items())),
        high_marked_types=high,corner_catalog=catalog,policies=policy_results,
        all_76_core_triangle_rows=triangle_rows,
        completeness='Every component has positive integer deficit2rank-baseline. A marked rank21 type has baseline>=35 iff its marked component deficit plus unmarked component deficits is<=7. The memoized rank/deficit unbounded enumeration covers every such multiset; colored-tree invariants quotient facet-root automorphisms. Separate full coin-change counts and exact modified-weight DPs cover all finite forests.',
        scope='Exact vertex-budget diagnostics on the complete inherited finite component catalogue, plus necessary core-triangle tests. Long-component discount rules are vertex-budget candidates only: arbitrary white extensions and all quadrilaterals require new complete coverage. Old paid-short-corner pruning is not used or assumed valid. No dimension-bound or actual-polytope realization claim.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['high_histogram'],flush=True)
