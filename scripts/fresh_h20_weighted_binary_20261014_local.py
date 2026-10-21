#!/usr/bin/env python3
"""Exact finite-neighbor graph and coloring obstruction at a binary cap.

A spherical Coxeter diagram has at most one label-five edge at any node.
Proper graph colorings bound each relevant finite-pair clique below rank20.
The saved verifier checks colorings directly, without trusting the heuristic.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import fresh_h20_weighted_binary_20261014 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_weighted_binary_20261014_local'
OUT=ROOT/f'results/{TAG}.json'

def graph():
    doc=json.loads((ROOT/'results/fresh_h20_weighted_binary_20261014.json').read_text())
    for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    masks=[m for m,r in doc['finite_neighbor_masks']]
    labels=[r for m,r in doc['finite_neighbor_masks']]+[2 if i in P.A else 3 for i in range(10)]+[4]+[5]*8
    V=[P.center(m)for m in masks]
    n=len(labels);assert n==114 and len(V)==95
    edges=[]
    for i in range(n):
        for j in range(i):
            if i>=95 and j>=95:r=2
            elif i>=95:
                axis=i-95;c=V[j][axis]
                if axis<10:r=2 if c==1 else 3 if c==P.F(1,2)else 0
                else:r=(4 if axis==10 else 5)if c==P.F(1,2)else 0
            else:r=P.rel(P.dist(V[i],V[j]))
            if r in (2,3,4,5):edges.append([j,i,r])
    assert len(edges)==2644
    return masks,labels,edges

def run(saved=None):
    masks,labels,edges=graph()
    E={tuple(e[:2])for e in edges}
    def adjacent(i,j):return tuple(sorted((i,j)))in E
    base=[i for i,r in enumerate(labels)if r!=5]
    choices=[-1]+[i for i,r in enumerate(labels)if r==5]
    assert len(base)==42 and len(choices)==73
    cases=[]
    if saved is None:
        import networkx as nx
    for k,p in enumerate(choices):
        nodes=base if p==-1 else [i for i in base if adjacent(i,p)]
        if saved is None:
            g=nx.Graph();g.add_nodes_from(nodes)
            g.add_edges_from((i,j)for i,j,r in edges if i in nodes and j in nodes)
            coloring=nx.coloring.greedy_color(g,strategy='saturation_largest_first')
        else:
            row=saved['coloring_cases'][k];assert row['chosen_pi5_neighbor']==p
            coloring={int(i):c for i,c in row['coloring'].items()}
        assert set(coloring)==set(nodes)
        assert all(isinstance(c,int)and c>=0 for c in coloring.values())
        assert all(coloring[i]!=coloring[j]for i,j,r in edges if i in coloring and j in coloring)
        colors=max(coloring.values(),default=-1)+1
        bound=1+colors+(p!=-1) # reference cap plus the chosen pi5 node
        assert bound<20
        cases.append(dict(chosen_pi5_neighbor=p,coloring=dict(sorted(coloring.items())),
                          color_count=colors,full_clique_upper_bound=bound))
    paths=[f'scripts/{TAG}.py','scripts/fresh_h20_weighted_binary_20261014.py',
           'results/fresh_h20_weighted_binary_20261014.json']
    return dict(status='PASS_EXACT_BINARY_CAP_FINITE_VERTEX_OBSTRUCTION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        binary_neighbor_masks=masks,reference_pair_labels=labels,
        finite_pair_graph_edges=edges,finite_pair_graph_vertices=114,
        finite_pair_graph_edge_count=2644,coloring_cases=cases,
        maximum_full_clique_bound=max(c['full_clique_upper_bound']for c in cases),
        bound_distribution=dict(Counter(c['full_clique_upper_bound']for c in cases)),
        scope='Every finite vertex has a positive-definite rank20 incident diagram. Such a diagram permits at most one pi5 neighbor of the reference root. The73 exact proper-coloring certificates bound every such finite-pair clique by19. The sign/permutation symmetry transports this to all524288 binary caps. No finite vertex of this chamber contains a binary cap. This does not exclude new neighbors after adding different walls.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    old=json.loads(OUT.read_text())if a.verify else None
    out=run(old)
    # JSON normalization makes integer keys and Counters stable on replay.
    out=json.loads(json.dumps(out))
    if a.verify:assert old==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['maximum_full_clique_bound'])
