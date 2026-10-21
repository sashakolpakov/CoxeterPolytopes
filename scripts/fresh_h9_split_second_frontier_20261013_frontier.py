#!/usr/bin/env python3
"""Exact finite-vertex and edge frontier census from an acute Gram dictionary.

This uses no numerical roots. The known numerical active sets propose starting
vertices; exact positive principal Grams and exhaustive ninth-wall tests decide
every endpoint in their finite-vertex component.
"""
import argparse
import json
from itertools import combinations
from scripts.fresh_h9_split_second_frontier_20261013_base import *
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_replaced_core_20261013_cores import inertia

F=Quadratic(K,add(scale(O,2),P))
WEIGHTS={2:F.zero,3:F.neg(F.one),5:F.lift(K.neg(P)),10:(Z,K.neg(O))}

def inverse_field(matrix):
    n=len(matrix);a=[list(row)+[F.one if i==j else F.zero for j in range(n)]for i,row in enumerate(matrix)]
    for i in range(n):
        assert F.sign(a[i][i])>0
        pivot=a[i][i];a[i]=[F.div(x,pivot)for x in a[i]]
        for j in range(n):
            if j==i:continue
            t=a[j][i];a[j]=[F.sub(x,F.mul(t,y))for x,y in zip(a[j],a[i],strict=True)]
    return [row[n:]for row in a]

def build(geometry):
    names=geometry['walls'];pairs=geometry['pair_checks']
    assert all(r['kind']in('finite','dotted')for r in pairs.values())
    def rel(a,b):return pairs['|'.join(sorted((a,b)))]
    def gram(ids):return [[F.rational(2)if a==b else WEIGHTS[rel(a,b)['label']]for b in ids]for a in ids]
    initial=[tuple(sorted(r['active_names']))for r in geometry['frontier']['vertices']]
    assert len(initial)==len(set(initial))
    vertices={};edges={};pending=list(initial);seen=set();tests=0
    while pending:
        ids=pending.pop()
        if ids in seen:continue
        seen.add(ids);assert len(ids)==9
        g=gram(ids);assert inertia(F,g)==(9,0,0)
        gi=inverse_field(g);assert all(F.sign(x)>=0 for row in gi for x in row)
        vertices[ids]={'walls':list(ids),'gram':[[F.encode(x)for x in row]for row in g],
                       'inverse':[[F.encode(x)for x in row]for row in gi]}
        for released in ids:
            edge=tuple(n for n in ids if n!=released)
            if edge in edges:continue
            assert inertia(F,gram(edge))==(8,0,0)
            rows=[];finite=[];ideal=[]
            for n in names:
                if n in edge:continue
                if any(rel(n,m)['kind']!='finite'for m in edge):continue
                sig=inertia(F,gram(edge+(n,)));tests+=1
                rows.append({'ninth_wall':n,'inertia':list(sig)})
                if sig==(9,0,0):
                    corner=tuple(sorted(edge+(n,)));finite.append(corner);pending.append(corner)
                elif sig==(8,0,1):ideal.append(n)
                else:assert sig==(8,1,0)
            assert len(finite)in(1,2)and not ideal
            edges[edge]={'walls':list(edge),'finite_endpoints':[list(x)for x in sorted(finite)],
                         'ideal_ninth_walls':ideal,'ninth_wall_tests':rows,
                         'classification':'FINITE_FINITE'if len(finite)==2 else'FINITE_OPEN'}
    assert seen==set(initial),'Numerical frontier omitted an exact finite vertex'
    # Reconcile every displayed directed ray, independently of its floating times.
    for row in geometry['frontier']['edges']:
        e=edges[tuple(sorted(row['edge_wall_names']))]
        expected='CLOSED_AT_EXISTING_FINITE_VERTEX'if e['classification']=='FINITE_FINITE'else'OPEN_TO_LIGHT_CONE_WITH_CURRENT_WALLS'
        assert row['classification']==expected
    finiteedges=[r for r in edges.values()if r['classification']=='FINITE_FINITE']
    components=[];remaining=set(vertices)
    while remaining:
        component={min(remaining)};stack=list(component)
        while stack:
            v=stack.pop()
            for r in finiteedges:
                a,b=map(tuple,r['finite_endpoints'])
                if v not in(a,b):continue
                w=b if v==a else a
                if w not in component:component.add(w);stack.append(w)
        remaining-=component;components.append(sorted(map(list,component)))
    adjacency={v:set()for v in vertices}
    for e in finiteedges:
        a,b=map(tuple,e['finite_endpoints']);adjacency[a].add(b);adjacency[b].add(a)
    compact_faces=[]
    for name in names:
        if not name.startswith('second'):continue
        vv=[v for v in vertices if name in v];assert len(vv)==2
        ends=[]
        for v in vv:
            old_neighbors=[w for w in adjacency[v]if not any(n.startswith('second')for n in w)]
            assert len(old_neighbors)==1;ends.append(old_neighbors[0])
        cycle=[ends[0],*vv]if ends[0]==ends[1]else[ends[0],*vv,ends[1]]
        assert all(b in adjacency[a]for a,b in zip(cycle,cycle[1:]+cycle[:1],strict=True))
        common=tuple(sorted(set.intersection(*(set(v)for v in cycle))))
        assert len(common)==7 and inertia(F,gram(common))==(7,0,0)
        compact_faces.append({'new_cap':name,'common_walls':list(common),
                              'ordered_vertex_wall_sets':[list(v)for v in cycle],
                              'sides':len(cycle)})
    return {'status':'EXACT_FINITE_VERTEX_COMPONENT_AND_EDGE_ENDPOINT_CENSUS','gram_field':field_record(F),
            'vertices':[vertices[k]for k in sorted(vertices)],'edges':[edges[k]for k in sorted(edges)],
            'finite_vertices':len(vertices),'finite_finite_edges':len(finiteedges),
            'finite_open_edges':len(edges)-len(finiteedges),'finite_graph_components':len(components),
            'finite_graph_cycle_rank':len(finiteedges)-len(vertices)+len(components),
            'ninth_wall_tests':tests,'components':components,'compact_two_faces':compact_faces,
            'proof':'For each positive9 principal Gram M, M^-1 is entrywise nonnegative. Projecting a common strict timelike center onto these nine walls gives a future timelike point, exactly nine zero inequalities, and every other inequality strict. Every actual endpoint of a rank8 spherical edge requires a positive or positive-semidefinite ninth-wall Gram. All remaining walls are tested; a dotted pair already excludes positive semidefiniteness. Thus each saved finite-open edge reaches the ideal boundary with no other wall active. This classifies the connected finite-vertex components reached here, not disconnected components.'}

def run():
    geometry=read(PREFIX+'_geometry.json');result=build(geometry)
    _,_,_,_,_,oldcase,_=current();snapshot=read(PREFIX+'_source.json')
    old=build(dict(oldcase,frontier=snapshot['frontier']))
    oldedges={tuple(r['walls']):r for r in old['edges']};newedges={tuple(r['walls']):r for r in result['edges']}
    changed=[list(k)for k,r in oldedges.items()if r['classification']=='FINITE_OPEN'and newedges[k]['classification']=='FINITE_FINITE']
    assert len(changed)==10
    return {'source':PREFIX+'_geometry.json','source_sha256':sha(PREFIX+'_geometry.json'),
            'baseline_source':PREFIX+'_source.json','baseline_sha256':sha(PREFIX+'_source.json'),
            'baseline':old,'old_open_edges_now_closed':changed,**result}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=run()
    path=ROOT/(PREFIX+'_frontier.json')
    if a.verify:assert json.loads(path.read_text())==r
    else:path.write_text(json.dumps(r,indent=2)+'\n')
    print({k:v for k,v in r.items()if k not in('vertices','edges','proof','components','baseline','compact_two_faces','old_open_edges_now_closed')})
