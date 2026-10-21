#!/usr/bin/env python3
"""Numerical proposal of bounded-distance corners; not an exclusion certificate."""
import json
from itertools import combinations
from pathlib import Path
import numpy as np
import networkx as nx

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_cubic_corner_20261014'
theta=2*np.cos(2*np.pi/7)
def val(x):return float(x['a'])+theta*float(x['b'])+theta**2*float(x['c'])

def run():
    out=[]
    for seed in range(3):
        source=f'results/fresh_cubic7_h9_seed{seed}_h100_20261013.json'
        d=json.loads((ROOT/source).read_text());a=val(d['timelike_coefficient'])
        Q=np.diag([-a]+[2.]*8+[1.])
        for i,j in ((1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8)):Q[i,j]=Q[j,i]=-1
        N=np.array([[val(x)for x in r]for r in d['all_walls']]);G=N@Q@N.T
        norm=np.diag(G);C=G/np.sqrt(norm[:,None]*norm[None,:]);graph=nx.Graph();graph.add_nodes_from(range(len(N)))
        graph.add_edges_from((i,j)for i in range(len(N))for j in range(i)if C[i,j]>-1+1e-7)
        candidates=set()
        for clique in nx.find_cliques(graph):
            if len(clique)>=9:candidates.update(combinations(sorted(clique),9))
        rows=[];far=[]
        for ids in sorted(candidates):
            B=G[np.ix_(ids,ids)]
            if np.linalg.eigvalsh(B)[0]<1e-7:continue
            M=(N[list(ids)]@Q)[:,1:]
            if abs(np.linalg.det(M))<1e-7:continue
            p=np.r_[1.,np.linalg.solve(M,-(N[list(ids)]@Q)[:,0])]
            n=p@Q@p
            if n>=-1e-7 or np.max(N@Q@p)>1e-6:continue
            height=a/(-n)-1
            item=dict(incident=list(ids),point=p.tolist(),norm=float(n),distance_height=float(height),
                      gram=np.round(B,9).tolist())
            (far if height>=100-1e-6 else rows).append(item)
        out.append(dict(seed=seed,source=source,finite_pair_clique_candidates=len(candidates),
                        proposed_corners=rows,beyond_bound_proposals=far))
        print('seed',seed,'clique candidates',len(candidates),'bounded corners',len(rows),'far',[(x['incident'],x['distance_height'])for x in far],flush=True)
    (ROOT/f'results/{TAG}_scout.json').write_text(json.dumps(dict(status='NUMERICAL_CORNER_PROPOSALS_ONLY',seeds=out),indent=2)+'\n')
if __name__=='__main__':run()
