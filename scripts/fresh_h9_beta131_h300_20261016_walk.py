#!/usr/bin/env python3
"""Exact finite-edge exploration from all verified height-100 finite corners."""
import argparse
import json
from collections import deque
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A
import fresh_h9_beta55_corner_symmetries_20261014 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta131_h300_20261016_walk'
PREFIX='results/fresh_h9_beta131_h300_20261016.json'
OLD='results/fresh_h9_beta131_h100_20261016_geometry.json'
LIMIT=200
K=A.K

def sha(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def read(path):
    d=json.loads((ROOT/path).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    return d
def vector(row):return [A.dec(x) for x in row]
def normalized(v):
    z=next(x for x in v if x)
    if z.sign()<0:z=-z
    return [x/z for x in v]

def run():
    d,old=read(PREFIX),read(OLD)
    form=[vector(row) for row in d['form']]
    roots=[vector(row) for row in d['roots']]
    center=vector(d['strict_center'])
    escape=vector(old['strict_escape']['point'])
    assert A.pair(escape,escape,form).sign()>0 and A.pair(escape,center,form).sign()<0
    escape_pairs=[A.pair(escape,r,form) for r in roots]
    assert all(x.sign()<0 for x in escape_pairs)
    gram=[vector(row) for row in d['gram']]
    covectors=[A.mv(form,r) for r in roots]
    vertices={};values_by_vertex={};queue=deque()
    def add(v,expected=None):
        v=normalized(v)
        values=A.mv(covectors,v)
        assert A.pair(v,v,form).sign()<0 and A.pair(v,center,form).sign()<0
        assert all(x.sign()<=0 for x in values)
        support=tuple(i for i,x in enumerate(values) if not x)
        assert len(support)==9
        if expected is not None:assert support==tuple(expected)
        if support not in vertices:
            vertices[support]=dict(incident=list(support),point=v,norm=A.pair(v,v,form))
            values_by_vertex[support]=values
            queue.append(support)
        return support
    for row in old['finite_vertices']:add(vector(row['point']),row['incident'])
    processed=[];edges={};open_edges={};supports={}
    while queue and len(vertices)<=LIMIT:
        active=queue.popleft();v=vertices[active]['point'];processed.append(list(active))
        for released in active:
            support=tuple(i for i in active if i!=released)
            if support in supports:continue
            B=[[gram[i][j] for j in support] for i in support]
            inverse,det=M.inverse_det(B)
            # Positive definiteness follows from the genuine finite corner;
            # positive pivots are checked directly as an independent guard.
            temp=[row[:] for row in B]
            for k in range(8):
                assert temp[k][k].sign()>0
                for i in range(k+1,8):
                    for j in range(k+1,8):
                        temp[i][j]-=temp[i][k]*temp[k][j]/temp[k][k]
            c=A.mv(inverse,[gram[i][released] for i in support])
            u=[-x+sum((c[j]*roots[i][k] for j,i in enumerate(support)),K())
               for k,x in enumerate(roots[released])]
            assert all(A.pair(u,roots[i],form)==0 for i in support)
            assert A.pair(u,v,form)==0 and A.pair(u,u,form).sign()>0
            assert A.pair(u,roots[released],form).sign()<0
            v_pairs=values_by_vertex[active]
            u_pairs=[-gram[released][i]+sum((c[j]*gram[s][i] for j,s in enumerate(support)),K())
                     for i in range(len(roots))]
            bounds=[(-a/b,i) for i,(a,b) in enumerate(zip(v_pairs,u_pairs)) if b.sign()>0]
            assert all(t.sign()>0 for t,i in bounds)
            record=dict(support=list(support),from_vertex=list(active),released=released,
                        direction=u,direction_norm=A.pair(u,u,form))
            if bounds:
                time=min(t for t,i in bounds)
                end=[a+time*b for a,b in zip(v,u)]
                q=A.pair(end,end,form)
                assert all((a+time*b).sign()<=0 for a,b in zip(v_pairs,u_pairs))
                record.update(first_wall_parameter=time,first_wall_point=end,first_wall_norm=q)
                if q.sign()<0:
                    target=add(end)
                    assert target!=active and set(active).intersection(target)==set(support)
                    record.update(kind='finite',to_vertex=list(target))
                    edges[support]=record
                else:
                    record.update(kind='ideal' if not q else 'open_before_first_wall')
                    open_edges[support]=record
            else:
                record.update(kind='open_without_wall_bound')
                open_edges[support]=record
            supports[support]=record
        print('finite walk processed',len(processed),'vertices',len(vertices),'pending',len(queue),flush=True)
    complete=not queue
    paths=[PREFIX,OLD,f'scripts/{TAG}.py',
           'scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py',
           'scripts/fresh_h9_beta55_corner_symmetries_20261014.py']
    return A.enc(dict(status='PASS_EXACT_BETA131_H300_FINITE_EDGE_WALK',sources={p:sha(p) for p in paths},
        limit=LIMIT,complete_seeded_finite_components=complete,pending=[list(x) for x in queue],
        processed_vertices=processed,finite_vertices=[vertices[x] for x in sorted(vertices)],
        finite_edges=[edges[x] for x in sorted(edges)],open_edges=[open_edges[x] for x in sorted(open_edges)],
        strict_escape=dict(point=escape,norm=A.pair(escape,escape,form),
                           center_pairing=A.pair(escape,center,form),pairings=escape_pairs),
        counts=dict(initial_vertices=len(old['finite_vertices']),finite_vertices=len(vertices),
                    finite_edges=len(edges),open_edges=len(open_edges),processed=len(processed),pending=len(queue)),
        scope='All stored finite points and edge endpoints are exact and feasible. If the queue is empty this exhausts the finite-edge connected components containing the 45 old corners, not necessarily every finite corner of the full chamber. The inherited strict spacelike escape remains strict at all current walls and proves this finite prefix has infinite volume. No ambient cone census or full-lattice reflectivity claim.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text());print('PASS saved finite-edge walk replay')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
