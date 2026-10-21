#!/usr/bin/env python3
"""Numerical multi-edge cap coverage and coherent additions to the 72-wall cone.

All original walls are retained. This uses a bounded finite-angle comparison
and floating geometry, so none of its exclusions is an exact no-go theorem.
"""
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np

from scripts.fresh_h9_split_edge_frontier_20261013 import frontier

ROOT = Path(__file__).resolve().parents[1]
SCOUT = ROOT/'data/enumerations/fresh_h9_split_endpoint_walk_optional_anchor_20261013_scout.json'
FRONTIER = ROOT/'results/fresh_h9_split_edge_frontier_20261013.json'
OUT = ROOT/'data/enumerations/fresh_h9_split_multiple_edge_caps_20261013_scout.json'


def indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask -= bit


def run():
    Path('/tmp/fresh_h9_split_multiple_edge_caps_20261013.pid').write_text(str(os.getpid()))
    source = json.loads(SCOUT.read_text())
    original = source['source_snapshot']
    before = json.loads(FRONTIER.read_text())['cases'][0]
    edges = [e for e in before['edges'] if e['classification']=='OPEN_TO_LIGHT_CONE_WITH_CURRENT_WALLS']
    ids = [i for i,c in enumerate(source['candidates']) if not c['conflicts']]
    normals = np.array([source['candidates'][i]['normal_physical'] for i in ids])
    base = np.array(original['physical_normals'])
    light = np.array([e['light_cone_point'] for e in edges])
    cover = normals[:,1:] @ light.T - normals[:,0,None] > 1e-7
    coverage = [sum(1 << int(j) for j in np.flatnonzero(row)) for row in cover]
    assert all(coverage)
    pair = normals[:,1:]@normals[:,1:].T - normals[:,0,None]*normals[:,0]
    compatible = pair < -2-1e-7
    finite_i, finite_j = np.where((pair > -2+1e-7) & (pair < 1e-7))
    values = pair[finite_i, finite_j]
    proposed_labels = np.rint(np.pi/np.arccos(np.clip(-values/2,-1,1))).astype(int)
    accept = ((proposed_labels >= 2) & (proposed_labels <= 120) &
              (abs(values+2*np.cos(np.pi/proposed_labels)) < 1e-7))
    compatible[finite_i[accept], finite_j[accept]] = True
    np.fill_diagonal(compatible,False)
    assert np.array_equal(compatible,compatible.T)
    adj = [int.from_bytes(np.packbits(row,bitorder='little').tobytes(),'little') for row in compatible]
    print({'caps':len(ids),'coverage_histogram':dict(Counter(map(int,cover.sum(axis=1)))),
           'compatible_pairs':int(compatible.sum()//2)},flush=True)
    del pair,compatible,finite_i,finite_j,values,proposed_labels,accept

    rng = np.random.default_rng(20261013051)
    multi = np.flatnonzero(cover.sum(axis=1)>1).tolist()
    seen=set();choices=[]
    for trial in range(350):
        selected=[];remaining=(1<<len(ids))-1;covered=0
        if trial < len(multi):
            first=multi[trial];selected=[first];remaining=adj[first];covered=coverage[first]
        while remaining:
            options=list(indices(remaining))
            degrees=np.array([(adj[i]&remaining).bit_count() for i in options])
            gain=np.array([(coverage[i]&~covered).bit_count() for i in options])
            score = 20*gain+np.log1p(degrees)+rng.exponential(1+trial%5,size=len(options))
            chosen=options[int(np.argmax(score))]
            selected.append(chosen);covered|=coverage[chosen];remaining&=adj[chosen]
        key=tuple(sorted(ids[i]for i in selected))
        if key in seen:continue
        seen.add(key)
        choices.append({'proposal_ids':list(key),'closed_old_directions':covered.bit_count(),
                        'new_caps':len(selected),'multi_edge_caps':sum(coverage[i].bit_count()>1 for i in selected),
                        'old_edge_coverage':[j for j in range(len(edges))if covered>>j&1]})
    choices.sort(key=lambda r:(-r['closed_old_directions'],-r['multi_edge_caps'],-r['new_caps']))
    result={'status':'NUMERICAL_MULTI_EDGE_CAP_ADDITIONS',
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in (SCOUT,FRONTIER)},
            'scope':'The archived 3857 no-removal proposals only; numerical finite-angle alphabet 2..120, strict dotted pairs, and existing-wall edge traversal. No exhaustive root, exact compatibility or compactness claim.',
            'original_open_edges':edges,'proposal_ids':ids,
            'edge_coverage':[[int(j)for j in np.flatnonzero(row)]for row in cover],
            'coverage_histogram':dict(Counter(map(int,cover.sum(axis=1)))),
            'selections':choices,'frontier_cases':[]}
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print({'greedy_selections':len(choices),'best':choices[:3]},flush=True)
    for case in choices[:5]:
        extra=np.array([source['candidates'][i]['normal_physical'] for i in case['proposal_ids']])
        W=np.vstack([base,extra]);names=[f'wall{i}'for i in range(len(base))]+[f'scout{i}'for i in case['proposal_ids']]
        graph=frontier(W,names,original['finite_vertices'])
        closed=graph['counts'].get('CLOSED_AT_EXISTING_FINITE_VERTEX',0)//2
        # The initial U8 and V8 components may merge; count components directly.
        parent=list(range(graph['finite_vertices_reached']))
        def root(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for e in graph['edges']:
            if 'end_vertex'in e:parent[root(e['start_vertex'])]=root(e['end_vertex'])
        components=len({root(i)for i in range(len(parent))})
        graph['finite_graph_components']=components
        graph['finite_graph_cycle_rank']=closed-len(parent)+components
        result['frontier_cases'].append({'selection':case,'frontier':graph})
        OUT.write_text(json.dumps(result,indent=2)+'\n')
        print({'facets':len(W),'vertices':len(parent),'cycle_rank':graph['finite_graph_cycle_rank'],
               'components':components,'counts':graph['counts']},flush=True)


if __name__=='__main__':run()
