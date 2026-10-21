#!/usr/bin/env python3
"""Exact combinatorial compact-face costs of three proposed wall releases.

The finite-vertex incidence list is frozen and complete. Closed edge degree
certifies bounded faces; no numerical scout decides the compact-face list.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_b_release_costs_20261017'
PARENT=ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json'
BOUNDARY=ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json'
S=(0,1,2,3,5,6)

def run():
    old=json.loads(PARENT.read_text()); bd=json.loads(BOUNDARY.read_text())['after64']
    vertices=[tuple(s)for s in old['retained_old_finite_supports']]+[tuple(v['support'])for v in old['new_finite_vertices']]
    assert len(vertices)==len(set(vertices))==198
    vertices=sorted(vertices); edges=defaultdict(list); faces=defaultdict(list)
    for i,v in enumerate(vertices):
        assert len(v)==9
        for e in combinations(v,8):edges[e].append(i)
        for size in range(8):
            for s in combinations(v,size):faces[s].append(i)
    assert all(len(v)<=2 for v in edges.values())
    compact={}
    for s,vs in faces.items():
        if len(vs)<10-len(s):continue
        ss=set(s); degree=9-len(s)
        if all(all(len(edges[tuple(j for j in vertices[i]if j!=k)])==2
                   for k in vertices[i]if k not in ss)for i in vs):
            compact[s]=dict(support=list(s),dimension=degree,vertices=vs)
    maximal=[]
    for s in sorted(compact,key=lambda x:(len(x),x)):
        if not any(set(t).issubset(s)for t in maximal):maximal.append(s)
    costs=[]
    for wall in[7,30,52]:
        lost=[i for i,v in enumerate(vertices)if wall in v]
        impacted=[f for f in compact.values()if set(f['vertices'])&set(lost)]
        maximal_impacted=[compact[s]for s in maximal if set(compact[s]['vertices'])&set(lost)]
        s6faces=[]
        for side in bd['actual_sides']:
            key=tuple(sorted(S+(side,)))
            if key in compact and any(wall in vertices[i]for i in compact[key]['vertices']):
                s6faces.append(dict(side=side,vertices=[list(vertices[i])for i in compact[key]['vertices']]))
        costs.append(dict(deleted_wall=wall,lost_finite_supports=[list(vertices[i])for i in lost],
            lost_finite_count=len(lost),compact_faces_touched_by_dimension=dict(sorted(Counter(f['dimension']for f in impacted).items())),
            maximal_compact_faces_touched=maximal_impacted,compact_S6_side_faces_touched=s6faces))
    return dict(status='PASS_EXACT_COMPACT_FACE_DELETION_COSTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),PARENT,BOUNDARY]},
        finite_vertices=[list(v)for v in vertices],finite_count=198,
        compact_face_counts=dict(sorted(Counter(f['dimension']for f in compact.values()).items())),
        maximal_compact_faces=[compact[s]for s in maximal],release_costs=costs,
        proof='Every face through a finite simple vertex is determined by a subset of its nine active walls. A face with a finite vertex is compact exactly when every incident edge at every one of its finite vertices has a second finite endpoint. For sufficiency, those vertices form a nonempty closed-degree component in the graph of a bounded projective truncation of the face. Graph connectivity forces that component to be the entire truncated polytope, so no truncation vertex or unbounded direction remains. Necessity is immediate. All dimension>=2 supports are exhausted from the complete finite census. Deletion costs mean lost prior face certificates, not a claim that no replacement face can exist.',
        scope='Exact incidence computation over the already certified finite vertices. No new walls and no ambient exterior cone.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()))
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['compact_face_counts'],flush=True)
    for r in out['release_costs']:print(r['deleted_wall'],r['lost_finite_count'],r['compact_faces_touched_by_dimension'],[x['side']for x in r['compact_S6_side_faces_touched']],flush=True)
if __name__=='__main__':main()
