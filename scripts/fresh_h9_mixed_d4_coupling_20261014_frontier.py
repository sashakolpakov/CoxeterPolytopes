#!/usr/bin/env python3
"""Exact finite-vertex edge queue of the frozen mixed D4 twelve-wall chamber.

Reuses the independent review's scalar arithmetic. No ray-census rerun, cap
search, arithmetic root enumeration or new-mirror firstness claim.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from fresh_h9_mixed_d4_coupling_20261014_review import E,Z,dec,dv,dm,enc,dot,mv,pair,inertia

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_h9_mixed_d4_coupling_20261014'
TAG=BASE+'_frontier'

def rank(M):
    A=[r[:] for r in M];k=0
    for j in range(len(A[0])):
        p=next((i for i in range(k,len(A)) if A[i][j]),None)
        if p is None:continue
        A[k],A[p]=A[p],A[k]
        for i in range(k+1,len(A)):
            if A[i][j]:
                c=A[i][j]/A[k][j]
                for a in range(j,len(A[0])):A[i][a]=A[i][a]-c*A[k][a]
        k+=1
        if k==len(A):break
    return k

def run():
    mainpath=ROOT/f'results/{BASE}.json';auditpath=ROOT/f'results/{BASE}_review.json'
    data=json.loads(mainpath.read_text());audit=json.loads(auditpath.read_text())
    for doc in (data,audit):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert audit['status']=='PASS_INDEPENDENT_QUARTIC_MIXED_D4_SUPPORT_CENSUS'
    M=dm(data['form']);N=dm(data['final_normals']);h=dv(data['center'])
    rays=data['final_complete_ray_census']['rays'];V=[dv(row['vector']) for row in rays]
    active=[frozenset(row['incident']) for row in rays]
    assert len(N)==12 and len(rays)==50
    assert all(len(a)==9 for a in active)
    assert data['final_complete_ray_census']['counts']=={'finite':12,'spacelike':38}
    finite=[i for i,row in enumerate(rays) if row['kind']=='finite']
    covectors=[mv(M,n) for n in N];G=[[pair(a,b,M) for b in N] for a in N]
    oldfinite={tuple(dv(row['vector'])) for row in data['seed_complete_ray_census']['rays'] if row['kind']=='finite'}
    finalfinite={tuple(V[i]) for i in finite};assert oldfinite<=finalfinite
    addedfinite=[i for i in finite if tuple(V[i]) not in oldfinite]
    support_checks={};queue=[];finite_edges=set();open_edges=[];adjacency={i:set() for i in finite}
    for start in finite:
        I=sorted(active[start]);p=V[start]
        assert inertia([[G[i][j] for j in I] for i in I])==[9,0]
        assert pair(p,p,M).sign()<0 and pair(p,h,M).sign()<0
        for released in I:
            S=tuple(i for i in I if i!=released)
            other=[i for i,a in enumerate(active) if i!=start and set(S)<=a]
            assert len(other)==1;end=other[0];q=V[end]
            assert active[start]&active[end]==frozenset(S)
            if S not in support_checks:
                assert rank([covectors[i] for i in S])==8
                assert inertia([[G[i][j] for j in S] for i in S])==[8,0]
                support_checks[S]=dict(support=S,rank=8,gram_signature=[8,0])
            # A positive combination is in the relative interior of this exact
            # 2D cone face: exactly these eight equations vanish.
            middle=[x+y for x,y in zip(p,q)]
            assert all(pair(middle,n,M).sign()==(0 if i in S else -1) for i,n in enumerate(N))
            direction=[y-x for x,y in zip(p,q)]
            quadratic=[pair(p,p,M),2*pair(p,direction,M),pair(direction,direction,M)]
            entry=dict(id=len(queue),start_ray=start,released_wall=released,support=S,
                other_cone_ray=end,start_point=p,projective_direction=direction,
                norm_polynomial_constant_linear_quadratic=quadratic,
                new_polar_wall_in_support=11 in S)
            if rays[end]['kind']=='finite':
                assert pair(q,q,M).sign()<0 and pair(p,q,M).sign()<0
                assert pair(middle,middle,M).sign()<0
                entry['kind']='finite_endpoint';entry['endpoint']=q
                finite_edges.add(tuple(sorted((start,end))));adjacency[start].add(end)
            else:
                assert rays[end]['kind']=='spacelike' and pair(q,q,M).sign()>0
                # Move strictly into the edge segment from its spacelike cone
                # endpoint. This produces an exact point beyond the unique
                # light-cone crossing, with no extra active old wall.
                for k in range(1,81):
                    t=Fraction(2**k-1,2**k)
                    w=[x+E(t)*v for x,v in zip(p,direction)]
                    if pair(w,w,M).sign()>0:break
                else:raise AssertionError('No certified interior spacelike continuation')
                vals=[pair(w,n,M) for n in N]
                assert all(x.sign()==(0 if i in S else -1) for i,x in enumerate(vals))
                assert pair(w,h,M).sign()<0
                assert quadratic[0].sign()<0
                assert (quadratic[0]+E(t)*quadratic[1]+E(t*t)*quadratic[2])==pair(w,w,M)
                entry.update(kind='open_to_light_cone',light_parameter_bracket=['0',str(t)],
                    interior_spacelike_parameter=str(t),interior_spacelike_point=w,
                    interior_spacelike_norm=pair(w,w,M),interior_spacelike_wall_pairings=vals)
                open_edges.append(entry['id'])
            queue.append(entry)
    kinds=Counter(r['kind'] for r in queue)
    assert len(queue)==108 and kinds=={'finite_endpoint':56,'open_to_light_cone':52}
    assert len(finite_edges)==28 and len(support_checks)==80
    # Every closed edge is listed once from each finite endpoint.
    for a,b in finite_edges:
        assert sum(r['start_ray']==a and r['other_cone_ray']==b for r in queue)==1
        assert sum(r['start_ray']==b and r['other_cone_ray']==a for r in queue)==1
    components=[];unseen=set(finite)
    while unseen:
        todo=[min(unseen)];component=set()
        while todo:
            v=todo.pop()
            if v in component:continue
            component.add(v);todo.extend(adjacency[v]-component)
        unseen-=component;components.append(sorted(component))
    open_release_counts=Counter(r['released_wall'] for r in queue if r['kind']=='open_to_light_cone')
    source_open=[r['id'] for r in queue if r['kind']=='open_to_light_cone' and r['released_wall']<7]
    assert len(source_open)==8 and all(queue[i]['start_ray'] in addedfinite for i in source_open)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}_review.py',f'results/{BASE}_review.json',f'results/{BASE}.json']
    return enc(dict(status='PASS_EXACT_108_FINITE_CORNER_DIRECTIONS_MIXED_D4',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        arithmetic_scope='Reuses frozen independent quartic review arithmetic; no independent second scalar implementation in this derivative',
        completeness_dependency='Frozen complete 50-ray cone census with independently matched support enumeration',
        finite_vertex_ray_ids=finite,new_polar_vertex_ray_ids=addedfinite,all_finite_vertices_simple=True,
        directed_direction_counts=dict(kinds),unique_finite_edges=sorted(finite_edges),
        finite_edge_graph_components=components,distinct_edge_supports=list(support_checks.values()),
        open_direction_ids=open_edges,open_releases_by_wall=dict(sorted(open_release_counts.items())),
        source_wall_releases_at_new_polar_vertices=source_open,
        queue=queue,
        scope='Complete edge queue incident to the twelve existing finite vertices only. Each open direction reaches the light cone before its spacelike projective endpoint and has exactly eight incident walls there. No additional mirror or cap search, arithmetic first-mirror claim, vertex-free-edge census, or finite-volume claim.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['directed_direction_counts'],flush=True)
