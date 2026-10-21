#!/usr/bin/env python3
"""Exact finite-vertex/frontier census of the frozen eleven-wall prefix."""
from itertools import combinations
from fresh_h9_e8_first_20261013 import *
OUT=ROOT/'results'/f'{PREFIX}_frontier.json'
def pd(M):
    a=[list(r)for r in M];pivots=[]
    for i in range(len(a)):
        p=a[i][i];pivots.append(p)
        if p<=0:return False,pivots
        for j in range(i+1,len(a)):
            for k in range(i+1,len(a)):a[j][k]-=a[j][i]*a[i][k]/p
    return True,pivots
def run(verify=False):
    names=[PREFIX,PREFIX+'_independent',PREFIX+'_e7',PREFIX+'_e7_independent',PREFIX+'_comparison']
    paths=[ROOT/'results'/f'{s}.json'for s in names];d,ai,e,ei,comp=[json.loads(p.read_text())for p in paths]
    assert ai['status']=='PASS_INDEPENDENT_E8_FIXED_COMPLETE_FIRST_EDGE'and ei['status']=='PASS_INDEPENDENT_COMPLETE_E7_A1_INITIAL_EDGE'
    for obj in(ai,ei):
        for path,h in obj['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
    walls=[tuple(map(K.from_json,r))for r in e['combined_walls']];cg=tuple(tuple(pair(r,G,s)for s in walls)for r in walls)
    points=[(O,*(Z for _ in range(9))),tuple(map(K.from_json,d['endpoint'])),tuple(map(K.from_json,e['endpoint']))]
    frames=[list(range(9)),list(range(8))+[9],[0,2,3,4,5,6,7,8,10]]
    for p,ss in zip(points,frames):
        assert pair(p,G,p)<0
        assert all(pair(p,G,walls[i])==0 if i in ss else pair(p,G,walls[i])<0 for i in range(11))
    finite=[]
    for ss in combinations(range(11),9):
        positive,pivots=pd(tuple(tuple(cg[i][j]for j in ss)for i in ss))
        if positive:finite.append(ss)
    assert set(finite)=={tuple(sorted(s))for s in frames}
    edges={}
    for k,ss in enumerate(frames):
        for released in ss:
            edge=tuple(sorted(set(ss)-{released}));edges.setdefault(edge,[]).append({'corner':k,'released_wall':released})
    rows=[]
    for edge,starts in sorted(edges.items()):
        candidates=[];gates=[]
        for i in range(11):
            if i in edge:continue
            ss=(*edge,i);positive,pivots=pd(tuple(tuple(cg[j][k]for k in ss)for j in ss))
            gates.append({'wall':i,'positive_definite':positive,'pivots':enc(pivots)})
            if positive:candidates.append(i)
        assert len(candidates)==len(starts)in(1,2)
        rows.append({'incident_walls':list(edge),'starting_corners':starts,'finite_end_walls':candidates,
                     'status':'OPEN_EXISTING_WALL_FRONTIER'if len(starts)==1 else'CERTIFIED_FIRST_EDGE',
                     'all_remaining_wall_schur_gates':gates})
    assert len(rows)==25 and sum(r['status']=='OPEN_EXISTING_WALL_FRONTIER'for r in rows)==23
    out={'status':'PASS_FROZEN_E8_E7_ELEVEN_WALL_FRONTIER',
         'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
         'facets':11,'certified_first_edges':2,'full_finite_corners':3,'nine_wall_subsets_checked':55,
         'corner_types':['E8+A1','E8+A1','E7+A1+A1'],'full_corner_mirror_counts':[121,121,65],
         'corners':[{'point':enc(p),'incident_walls':ss}for p,ss in zip(points,frames)],
         'edge_frontier':rows,'open_existing_wall_directions':23,
         'strict_escape':e['strict_escape'],'escape_norm':e['escape_norm'],
         'scope':'Two first steps are complete for the full specified arithmetic root arrangement, with explicit certified PARI dependencies. The eleven-wall prefix has exactly three finite vertices and an exact strict escape. Twenty-three outgoing edges have no second endpoint among these eleven walls; their full-root first mirrors remain uncomputed. Full indefinite conjugacy and full-host reflectivity remain open.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'corners',3,'open directions',23,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
