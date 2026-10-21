#!/usr/bin/env python3
"""Complete context-dependent cusp payments in the frozen H17 chamber.

No optimizer, new chamber, or ambient dimension claim.  Reconstruct actual
corner incidences and all fixed finite payments, then verify three explicit
integral allocations and the two-cusp Hall inequalities exactly.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h17_contextual_cusp_payments_20261017'
OLD = 'fresh_nonincident_carrier_incidence_20261017'
BOT = 'fresh_h17_d16_payment_bottleneck_20261017'
NODES = set(range(19))
EDGES = [(0,1),(0,2),(0,4),(2,3),(4,5),(5,6),(6,7),(7,8),(8,9)]
EDGES += [(a+10 if a<9 else 9,b+10 if b<9 else 9) for a,b in EDGES.copy()]
ADJ = {i:set() for i in NODES}
for a,b in EDGES: ADJ[a].add(b); ADJ[b].add(a)


def sha(p):
    h=sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def component(nodes,a):
    nodes=set(nodes); d={a:0}; q=deque([a])
    while q:
        u=q.popleft()
        for v in sorted((ADJ[u]&nodes)-d.keys()): d[v]=d[u]+1; q.append(v)
    return d


def run():
    old=json.loads((ROOT/f'results/{OLD}.json').read_text())
    bottleneck=json.loads((ROOT/f'results/{BOT}.json').read_text())
    vertices=old['vertices']; faces=old['face_ledger']['two_faces']
    assert len(vertices)==82 and len(faces)==3456
    cusp_ids={}
    for i,v in enumerate(vertices):
        if v['ideal']:
            name='E8xE8' if v['support']==[9] else 'D16'
            assert v['support'] in ([9],[3,13]); cusp_ids[i]=name
    assert sorted(cusp_ids.values())==['D16','E8xE8']
    finite_incidence=Counter(); finite_cost=Counter(); ideal_incidence=Counter()
    rows=[]; categories=Counter()
    for index,face in enumerate(faces):
        assert face['polygon'] in (3,4)
        corners=[]; eligible=[]; finite_payment=0
        for vi in face['vertices']:
            vertex=vertices[vi]
            marks=sorted(set(face['side'])-set(vertex['support']))
            assert len(marks)==2
            diagram=set(face['white'])|set(marks)
            ds=component(diagram,marks[0]); distance=ds.get(marks[1])
            if vertex['ideal']:
                cusp=cusp_ids[vi]; eligible.append(cusp); ideal_incidence[cusp]+=1
                assert distance is not None
                assert len(ds)==(9 if cusp=='E8xE8' else 17)
                payment=None
            else:
                cusp=None; finite_incidence[vi]+=1
                payment=int(distance is not None and (distance<=9 or len(ds)<=17))
                finite_cost[vi]+=payment; finite_payment+=payment
            corners.append(dict(vertex=vi,marks=marks,cusp=cusp,
                component_order=len(ds),distance=distance,fixed_finite_payment=payment))
        assert len(set(eligible))==len(eligible)
        demand=2 if face['polygon']==3 else 1
        deficit=max(0,demand-finite_payment)
        assert 0<=deficit<=len(eligible)
        categories[face['polygon'],len(eligible),deficit]+=1
        rows.append(dict(face_index=index,polygon=face['polygon'],side=face['side'],
            white=face['white'],corners=corners,eligible_cusps=sorted(eligible),
            required_payment=demand,fixed_finite_payment=finite_payment,
            ideal_deficit=deficit))
    assert categories=={(3,0,0):1456,(3,1,0):44,(3,1,1):644,
        (3,2,1):2,(3,2,2):14,(4,0,0):1232,(4,1,0):15,(4,1,1):49}
    assert len(finite_incidence)==80 and set(finite_incidence.values())=={136}
    assert ideal_incidence==dict(E8xE8=648,D16=136)
    # All finite vertices are simple. Reconstruct their full connected
    # component pair costs independently from the face accumulation.
    finite_rows=[]
    for vi in sorted(finite_incidence):
        remaining=NODES-set(vertices[vi]['support']); sizes=[]
        while remaining:
            c=set(component(remaining,min(remaining))); sizes.append(len(c)); remaining-=c
        cost=sum(k*(k-1)//2 for k in sizes)
        assert cost==finite_cost[vi] and cost<=136<=9*17
        finite_rows.append(dict(vertex=vi,component_orders=sorted(sizes),cost=cost))
    total=sum(r['ideal_deficit'] for r in rows); assert total==723
    # The exact minimum-cut demand for every subset of these two cusps.
    hall=[]
    for size in range(3):
        for subset in combinations(['D16','E8xE8'],size):
            terms=[]; U=set(subset)
            for r in rows:
                q=max(0,r['ideal_deficit']-len(set(r['eligible_cusps'])-U))
                if q:terms.append([r['face_index'],q])
            hall.append(dict(cusps=list(subset),demand=sum(q for _,q in terms),terms=terms))
    assert [r['demand'] for r in hall]==[0,109,612,723]
    flexible=[r['face_index'] for r in rows if len(r['eligible_cusps'])==2 and r['ideal_deficit']==1]
    assert len(flexible)==2
    certificates=[]
    for to_D in range(3):
        D_choices=set(flexible[:to_D]); payments=[]; costs=Counter()
        for r in rows:
            delta=r['ideal_deficit']; I=r['eligible_cusps']; chosen=[]
            if delta==2:chosen=I.copy()
            elif delta==1:
                chosen=I.copy() if len(I)==1 else ['D16' if r['face_index'] in D_choices else 'E8xE8']
            assert len(chosen)==delta and set(chosen)<=set(I)
            assert r['fixed_finite_payment']+len(chosen)>=r['required_payment']
            for c in chosen:costs[c]+=1
            if chosen:payments.append(dict(face_index=r['face_index'],paid_cusps=chosen))
        assert costs==dict(D16=109+to_D,E8xE8=614-to_D)
        assert sum(costs.values())==723
        for h in hall:assert sum(costs[c] for c in h['cusps'])>=h['demand']
        certificates.append(dict(flexible_triangles_paid_at_D16=to_D,cusp_costs=dict(costs),
            total_ideal_cost=723,nonzero_payments=payments))
    # The first allocation restricts exactly to the previously reviewed109
    # mandatory D16 pairs, rather than silently changing that baseline.
    first_paid={r['face_index'] for r in certificates[0]['nonzero_payments'] if 'D16' in r['paid_cusps']}
    by_white={tuple(r['white']):r for r in rows}
    for b in bottleneck['faces']:
        r=by_white[tuple(b['white'])]
        assert int(r['face_index'] in first_paid)==b['D16_payment_deficit']
    normalizations={'D16':17,'E8xE8':162}
    ratios=[dict(cusps=h['cusps'],ratio=str(F(h['demand'],sum(normalizations[c] for c in h['cusps']))))
            for h in hall if h['cusps']]
    assert max(F(h['ratio']) for h in ratios)==F(109,17)
    assert all(F(certificates[0]['cusp_costs'][c])<=F(109,17)*normalizations[c] for c in normalizations)
    paths=[Path(__file__)]
    for stem in(OLD,BOT):paths += [ROOT/f'scripts/{stem}.py',ROOT/f'results/{stem}.json',ROOT/f'results/{stem}.md']
    return dict(status='PASS_COMPLETE_H17_CONTEXTUAL_CUSP_PAYMENTS_AND_EXACT_HALL_REGION',
        sources={str(p.relative_to(ROOT)):sha(p) for p in paths},
        fixed_finite_rule='unit for marked distance<=9 OR connected marked component order<=17',
        ideal_variables='independent actual corner weights in[0,1], with no intrinsic-diagram invariance imposed',
        faces=rows,finite_vertex_costs=finite_rows,
        face_categories=[dict(polygon=k[0],ideal_corners=k[1],ideal_deficit=k[2],count=v)for k,v in sorted(categories.items())],
        ideal_corner_counts=dict(ideal_incidence),hall_constraints=hall,
        flexible_two_cusp_triangles=flexible,certificates=certificates,
        optimal_total_ideal_mass=723,exact_budget_region=['B_D16>=109','B_E8xE8>=612','B_D16+B_E8xE8>=723'],
        best_uniform_normalized_cost=str(F(109,17)),
        normalization=normalizations,uniform_budget_cut_ratios=ratios,
        universal_allocation_lemma='For fixed finite payments and eligible ideal corners I(F), unit-capacitated cusp budgets B(v) are feasible exactly when sum_F max(0,delta(F)-|I(F) outside U|)<=sum_{v in U}B(v) for every cusp subset U.',
        scope='Complete certificate in the already frozen known H17 chamber plus an abstract exact allocation criterion. No universal bound on the Hall demands, no new chamber and no dimension improvement.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'3456 faces;3 allocations;723 ideal mass;4 exact Hall cuts',flush=True)
