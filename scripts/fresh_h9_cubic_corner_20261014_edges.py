#!/usr/bin/env python3
"""Exact full elliptic rank8/rank9 census of the transported cubic H9 prefix."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
import networkx as nx
import fresh_h9_cubic_corner_20261014 as P
K=P.K;Z,O=P.Z,P.O;ROOT=P.ROOT;TAG='fresh_h9_cubic_corner_20261014_edges'

def run(source='results/fresh_h9_cubic_corner_20261014_extension.json'):
    d=json.loads((ROOT/source).read_text());bound=d.get('complete_distance_bound',100)
    N=[[K.decode(x)for x in r]for r in d['normals']];a=K.decode(d['timelike_coefficient']);n=len(N)
    G=[[K.pairing(r,s,a)for s in N]for r in N]
    graph=nx.Graph();graph.add_nodes_from(range(n))
    graph.add_edges_from((i,j)for i in range(n)for j in range(i)if K.sign(K.sub(K.mul(G[i][i],G[j][j]),K.mul(G[i][j],G[i][j])))>0)
    candidates={8:set(),9:set()}
    for clique in nx.find_cliques(graph):
        for k in(8,9):
            if len(clique)>=k:candidates[k].update(combinations(sorted(clique),k))
    elliptic={8:[],9:[]}
    for k in(8,9):
        for ids in sorted(candidates[k]):
            try:P.pd([[G[i][j]for j in ids]for i in ids])
            except AssertionError:continue
            elliptic[k].append(ids)
    T=[[K.decode(x)for x in r]for r in d['isometry']];origin=[O]+[Z]*9;other=P.mv(T,origin)
    verts=[]
    for ids in elliptic[9]:
        cov=[K.gram_vector(N[i],a)for i in ids];inv,_=P.inverse([r[1:]for r in cov]);p=[O]+P.mv(inv,[K.neg(r[0])for r in cov])
        q=K.pairing(p,p,a);assert K.sign(q)<0
        ips=[K.pairing(p,r,a)for r in N];assert [i for i,x in enumerate(ips)if x==Z]==list(ids)and all(K.sign(x)<=0 for x in ips)
        heights=[K.sub(K.divide(K.mul(K.pairing(p,v,a),K.pairing(p,v,a)),K.mul(K.neg(q),a)),O)for v in(origin,other)]
        full=any(K.sign(K.sub(K.scale(O,bound),h))>0 for h in heights)
        verts.append(dict(incident=ids,point=p,norm=q,origin_distance_heights=heights,full_chamber_corner_certified=full))
    edges=[]
    for ids in elliptic[8]:
        ends=[i for i,v in enumerate(verts)if set(ids)<set(v['incident'])];assert len(ends)<=2
        edges.append(dict(incident=ids,finite_endpoints=ends))
    counts=dict(Counter(len(e['finite_endpoints'])for e in edges))
    paths=[source,'scripts/fresh_h9_cubic_corner_20261014.py',f'scripts/{TAG}.py']
    return P.encode(dict(status='EXACT_COMPLETE_ELLIPTIC_FACE_CENSUS_OF_PREFIX',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        finite_pair_clique_candidates={k:len(v)for k,v in candidates.items()},
        finite_vertices=verts,elliptic_rank8_edges=edges,finite_endpoint_counts=counts,
        complete_finite_vertex_count=len(verts),complete_elliptic_rank8_count=len(edges),
        full_chamber_corners_certified=sum(v['full_chamber_corner_certified']for v in verts),
        complete_distance_bound=bound,
        scope='Complete exact elliptic8/9 census of this finite prefix. Full-chamber corners use either of two transported complete distance neighborhoods. Open edges need undiscovered mirrors; no nonreflectivity conclusion.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--source',default='results/fresh_h9_cubic_corner_20261014_extension.json');ap.add_argument('--output',default=f'results/{TAG}.json');args=ap.parse_args();out=run(args.source);p=ROOT/args.output
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_finite_vertex_count'],out['full_chamber_corners_certified'],out['finite_endpoint_counts'])
