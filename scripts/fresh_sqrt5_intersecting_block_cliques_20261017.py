#!/usr/bin/env python3
"""Exact pairwise-intersecting block-orbit clique charges, all positive forests."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import networkx as nx

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_intersecting_block_cliques_20261017'
INPUT='fresh_sqrt5_all_subset_port_packing_20261017'


def run():
    src=ROOT/f'results/{INPUT}.json';old=json.loads(src.read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    alpha={i:F(a) for i,a in old['face_multipliers']}
    mu={v:F(a) for v,a in old['vertex_multipliers']}
    weights=defaultdict(F)
    for f in old['positive_simplex_features']:
        for j,m in f['block_orbits']:weights[j]+=alpha[f['face']]*m
    masks={c['index']:c['literal_masks'] for c in old['literal_block_orbits']}
    byvertex=defaultdict(list)
    for c in old['block_orbit_columns']:
        if weights[c['index']]:byvertex[c['vertex']].append(c['index'])
    records=[];critical=[];newmu=dict(mu);totals=Counter()
    for v,js in sorted(byvertex.items()):
        eligible=[];internal=[];g=nx.Graph()
        for j in js:
            bad=next(((a,b) for a,b in combinations(masks[j],2) if not a&b),None)
            if bad is None:eligible.append(j)
            else:internal.append([j,list(bad)])
        g.add_nodes_from(eligible);separators=[]
        for j,k in combinations(eligible,2):
            bad=next(((a,b) for a in masks[j] for b in masks[k] if not a&b),None)
            if bad is None:g.add_edge(j,k)
            else:separators.append([j,k,list(bad)])
        # With nonnegative loads it suffices to inspect every maximal clique.
        cliques=sorted(tuple(sorted(c)) for c in nx.find_cliques(g))
        for c in cliques:assert all(g.has_edge(j,k) for j,k in combinations(c,2))
        best=max(cliques,key=lambda c:(sum((weights[j] for j in c),F()),c),default=())
        load=sum((weights[j] for j in best),F())
        newmu[v]=max(mu[v],load)
        violating=[list(c) for c in cliques if sum((weights[j] for j in c),F())>mu[v]]
        row=dict(vertex=v,positive_block_orbits=len(js),eligible_block_orbits=eligible,
            internal_disjoint_witnesses=internal,edges=[list(x) for x in sorted(g.edges())],
            cross_disjoint_witnesses=separators,maximal_cliques=[list(c) for c in cliques],
            best_clique=list(best),best_load=str(load),old_vertex_mass=str(mu[v]),
            violating_cliques=violating)
        records.append(row);totals.update(positive_orbits=len(js),eligible_orbits=len(eligible),
            maximal_cliques=len(cliques),violating_cliques=len(violating))
        if load>mu[v]:
            chosen=set(best);terms=[]
            for f in old['positive_simplex_features']:
                coeff=sum(m for j,m in f['block_orbits'] if j in chosen)
                if coeff:terms.append([f['face'],coeff])
            assert sum((alpha[fi]*m for fi,m in terms),F())==load
            critical.append(dict(index=len(critical),vertex=v,vertex_cost=1,
                block_orbits=list(best),positive_simplex_coefficients=terms,load=str(load)))
    assert all(newmu[v]>=a for v,a in mu.items())
    num=F(old['objective_numerator']);den=sum(newmu.values(),F())
    paths=[f'scripts/{TAG}.py',f'results/{INPUT}.json',f'results/{INPUT}.md',
        'results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_ALL_POSITIVE_INTERSECTING_ORBIT_CLIQUES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        old_bound=old['bound'],bound=str(num/den),objective_numerator=str(num),normalizer=str(den),
        excludes_target12=num/den>12,counts=dict(totals,tested_forests=len(records),
            increased_vertex_masses=len(critical)),vertex_checks=records,critical_columns=critical,
        block_orbit_loads=[[j,str(a)] for j,a in sorted(weights.items())],
        face_multipliers=old['face_multipliers'],vertex_multipliers=[[v,str(a)] for v,a in sorted(newmu.items())],
        theorem='Distinct maximal-simplex port blocks are disjoint. A union of automorphism orbits in which every two literal blocks intersect contains at most one actual block, so its corner indicator has vertex capacity1. Each orbit must itself be pairwise intersecting; compatibility between distinct orbits is universal over their literal members.',
        scope='Every maximal clique of every positive-support orbit conflict graph is enumerated. Zero-load orbits cannot improve the maximum. This exhausts these invariant pairwise-intersection clique charges at this fixed witness, not every local packing inequality, feasible geometry, missing global face/field coverage or dimension bounds. No optimizer.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],d['bound'],float(F(d['bound'])),flush=True)
