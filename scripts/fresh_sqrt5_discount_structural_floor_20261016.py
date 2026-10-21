#!/usr/bin/env python3
"""Exact finite-forest barriers to the current adjacent-discount family.
No optimization and no assertion of compact ambient vertex realizability.
"""
import argparse,json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_discount_structural_floor_20261016'


def graph(name):
    r=int(name[1:]);g={i:set()for i in range(r)}
    edges=[(i,i+1)for i in range(r-1)]if name[0]=='A'else[(i,i+1)for i in range(r-2)]+[(1,r-1)]
    assert name[0]in'AD'and(name[0]!='D'or r>=4)
    for i,j in edges:g[i].add(j);g[j].add(i)
    return g


def component(name):
    g=graph(name);pairs=[]
    for i,j in combinations(g,2):
        d={i:0};q=deque([i])
        while q:
            u=q.popleft()
            for v in g[u]:
                if v not in d:d[v]=d[u]+1;q.append(v)
        pairs.append(dict(pair=[i,j],distance=d[j],baseline=int(d[j]<=2),
            old_discount_eligible=int(len(g)>=6 and d[j]==1)))
    gram=[[F(2*(i==j)-int(j in g[i]))for j in g]for i in g];work=[r[:]for r in gram];pivots=[]
    for k in range(len(g)):
        pivot=work[k][k];assert pivot>0;pivots.append(str(pivot))
        for i in range(k+1,len(g)):
            for j in range(k+1,len(g)):work[i][j]-=work[i][k]*work[k][j]/pivot
    return dict(name=name,rank=len(g),gram=[[str(z)for z in r]for r in gram],positive_pivots=pivots,edges=[list(e)for e in combinations(g,2)if e[1]in g[e[0]]],
        literal_pairs=pairs,baseline=sum(p['baseline']for p in pairs),
        adjacent_pairs=sum(p['distance']==1 for p in pairs),distance_two_pairs=sum(p['distance']==2 for p in pairs),
        old_discount_eligible=sum(p['old_discount_eligible']for p in pairs))


def threshold(d):
    coefficients=[F(comb(i,2)+comb(d-i,2),2)for i in range(d+1)]
    minimum=min(coefficients);closed=F(((d-1)**2)//4,2)
    assert minimum==closed and coefficients[0]>minimum
    return dict(dimension=d,paired_h_coefficients=list(map(str,coefficients)),
        lower_f2_over_f0=str(minimum),strict=True,target=str(5*minimum-comb(d,2)))


def run():
    models={n:component(n)for n in('A1','A2','D4','D5')}
    assert models['D4']['baseline']==6 and models['D5']['baseline']==8
    forests={17:['D5']*3+['A2'],14:['D5']*2+['D4'],13:['D5']+['D4']*2}
    records=[]
    for d,forest in forests.items():
        assert sum(models[n]['rank']for n in forest)==d
        baseline=sum(models[n]['baseline']for n in forest);eligible=sum(models[n]['old_discount_eligible']for n in forest)
        th=threshold(d);target=F(th['target']);assert eligible==0 and baseline>target
        adjacent=sum(models[n]['adjacent_pairs']for n in forest)
        records.append(dict(**th,forest=forest,baseline=str(baseline),old_eligible_adjacent_pairs=eligible,
            old_family_unavoidable_gap=str(baseline-target),
            if_all_adjacencies_were_discounted=dict(adjacent_pairs=adjacent,
                necessary_uniform_discount=str(F(baseline-target,adjacent)),
                distance_two_floor=str(sum(models[n]['distance_two_pairs']for n in forest)))))
    # Even arbitrary discounts of every adjacent corner, at most1 per corner,
    # leave distance-two costs. This is a broader barrier at dimension11.
    d=11;forest=['D5','D5','A1'];th=threshold(d);floor=sum(models[n]['distance_two_pairs']for n in forest)
    assert floor==8 and floor>F(th['target'])==F(15,2)
    broad=dict(**th,forest=forest,retained_distance_two_floor=str(floor),gap=str(floor-F(th['target'])))
    # The existing triangle stopping criterion is rank-independent. Its
    # arithmetic graph catalogue can be reused for adjacent-only discounts
    # up to1, if every individual corner deduction is at most1.
    import fresh_sqrt5_triangle_extension_complete_20261016 as T
    checks=[]
    for b in range(2,9):
        for k in range(3):
            safe=b>=4 or(b==3 and k<2)or(b==2 and k==0)
            assert safe==bool(T.safe_reason(dict(baseline=b,black3_edges=k)))
            if safe:assert b-k>=2
            checks.append(dict(baseline=b,black3_edges=k,hereditary_safe=safe,worst_payment_at_discount_one=str(b-k)))
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_triangle_extension_complete_20261016.py',
        'results/fresh_sqrt5_triangle_extension_complete_20261016_rank22.json',
        'scripts/fresh_sqrt5_discounted_face_domain_reduction_20261016.py',
        'results/fresh_sqrt5_discounted_face_domain_reduction_20261016.json']
    return dict(status='PASS_EXACT_ADJACENT_DISCOUNT_STRUCTURAL_BARRIERS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},components=models,
        old_rank_at_least_six_family_barriers=records,arbitrary_adjacent_only_barrier=broad,
        paired_h_argument='For a simple convex d-polytope, f0=sum h_i and f2=sum binom(d-i,2)h_i, with h_i=h_(d-i)>=0 and h0=hd=1. Pairing gives f2/f0 strictly above min_i[binom(i,2)+binom(d-i,2)]/2=floor((d-1)^2/4)/2. Hence total required face payment5f2-binom(d,2)f0 is strictly above target(d)*f0. A rule with every vertex<=target yields contradiction.',
        triangle_pruning_extension=dict(per_corner_adjacent_discount='0..1',distance_two_discount='0',bonuses='nonnegative',checks=checks,
            reason='A compact acute triangle has at most two black3edges by the previously proved positive-null-root obstruction. The existing complete graph domain uses only baseline and black3edge count for hereditary pruning, not rank eligibility.'),
        scope='Exact algebraic obstructions to these specified payment families and per-vertex maximum methods, not lower bounds on actual compact polytope dimensions. The witnesses are formal spherical vertex forests; no simultaneous compact realization is claimed. The current rank>=6-only adjacent discount family cannot exclude dimensions17,14,13 regardless of its nonnegative compensation coefficients. General adjacent-only discounts retaining every distance-two payment cannot exclude11. These are not universal barriers to weighted2face arguments with different payments, global incidence restrictions, or higher-face information. The triangle pruning extension requires a fresh evaluator and does not certify any new weight rule or quadrilateral coverage.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
