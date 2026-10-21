#!/usr/bin/env python3
"""Exact finite-edge and S6 boundary comparison for the split endpoint pair."""
import argparse
import json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_nonarithmetic_s6_split_pair_20261017_cone as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_s6_split_pair_20261017_edge_ledger'
Q,U,I,P=C.Q,C.U,C.I,C.P
S=[0,1,2,3,5,6]


def index_supports(rays):
    index=defaultdict(list)
    for i,r in enumerate(rays):
        for support in combinations(r['active'],8):index[support].append(i)
    return index


def face_boundary(rays,index,q,h):
    hc=U.imv(q,h);edges=[];counts=Counter()
    for support,ends in sorted(index.items()):
        if not set(S)<=set(support):continue
        assert len(ends)==2 and all(len(rays[i]['active'])==9 for i in ends)
        a,b=[Q.vector(rays[i]['normalized_computational_coordinate_bounds'])for i in ends]
        aa,ab,bb=U.idot(a,U.imv(q,a)),U.idot(a,U.imv(q,b)),U.idot(b,U.imv(q,b))
        finite=sum(rays[i]['kind']=='finite_vertices'for i in ends)
        det=U.rnd(aa*bb-ab*ab);witness=None
        if finite:
            assert det.h<0
            kind='compact_geodesic_edge'if finite==2 else'one_open_end'
        elif ab.h<0 and det.h<0:
            ca,cb=U.rnd(bb-ab),U.rnd(aa-ab);assert ca.l>0 and cb.l>0
            v=[U.rnd(ca*x+cb*y)for x,y in zip(a,b)]
            direct,future=U.idot(v,U.imv(q,v)),U.idot(v,hc)
            norm=U.rnd(det*(aa+bb-2*ab))
            assert norm.l<=direct.h and direct.l<=norm.h and norm.h<0 and future.h<0
            kind='two_open_ends';witness=dict(positive_coefficients=P.enc([ca,cb]),norm=norm.enc(),
                future=future.enc(),norm_identity='(aa*bb-ab^2)*(aa+bb-2ab)',direct_norm=direct.enc())
        else:
            assert ab.l>0 or det.l>0
            kind='outside_hyperbolic_space'
        counts[kind]+=1
        edges.append(dict(support=list(support),sides=sorted(set(support)-set(S)),endpoints=ends,
            endpoint_supports=[rays[i]['active']for i in ends],kind=kind,
            square_pairing_bounds=P.enc([aa,ab,bb]),determinant=det.enc(),future_timelike_witness=witness))
    face_indices=[i for i,r in enumerate(rays)if set(S)<=set(r['active'])]
    assert len(face_indices)*3==2*len(edges)
    outside={i for i in face_indices if rays[i]['kind']=='hyperideal_rays'}
    adjacency={i:set()for i in outside}
    for e in edges:
        if set(e['endpoints'])<=outside:
            a,b=e['endpoints'];adjacency[a].add(b);adjacency[b].add(a)
    components=[];unseen=set(outside)
    while unseen:
        start=min(unseen);seen={start};todo=[start]
        while todo:
            a=todo.pop()
            for b in adjacency[a]-seen:seen.add(b);todo.append(b)
        unseen-=seen
        internal=[e for e in edges if set(e['endpoints'])<=seen]
        terminal=[e for e in edges if len(set(e['endpoints'])&seen)==1 and not set(e['endpoints'])<=outside]
        components.append(dict(ray_indices=sorted(seen),ray_supports=[rays[i]['active']for i in sorted(seen)],
            internal_edges=[dict(sides=e['sides'],kind=e['kind'],endpoints=e['endpoints'])for e in internal],
            finite_terminal_edges=[dict(sides=e['sides'],kind=e['kind'],endpoints=e['endpoints'])for e in terminal]))
    return dict(face_ray_count=len(face_indices),face_counts=dict(Counter(rays[i]['kind']for i in face_indices)),
        edge_counts=dict(counts),edges=edges,exterior_components=components)


def run():
    cone,old=Q.read(C.TAG),Q.read(C.PRIOR)
    assert cone['status']=='PASS_COMPLETE_65_WALL_CONE_AND_STRICT_ESCAPE'
    assert old['status']=='PASS_COMPLETE_63_WALL_CONE_AND_STRICT_ESCAPE'
    rays=cone['final_rays'];assert all(len(r['active'])==9 for r in rays if r['kind']=='finite_vertices')
    index=index_supports(rays);finite_edges=[];counts=Counter()
    for support,ends in sorted(index.items()):
        finite=[i for i in ends if rays[i]['kind']=='finite_vertices']
        if not finite:continue
        assert len(ends)==2
        kind='compact_geodesic_edge'if len(finite)==2 else'open_geodesic_edge';counts[kind]+=1
        finite_edges.append(dict(support=list(support),endpoints=ends,kind=kind,
            endpoint_supports=[rays[i]['active']for i in ends],finite_endpoint_supports=[rays[i]['active']for i in finite]))
    assert 2*counts['compact_geodesic_edge']+counts['open_geodesic_edge']==9*cone['counts']['finite_vertices']
    q=[Q.vector(row)for row in cone['computational_gram']];h=Q.vector(cone['computational_common_center'])
    assert cone['computational_root_basis']==old['computational_root_basis']
    before=face_boundary(old['final_rays'],index_supports(old['final_rays']),q,h)
    after=face_boundary(rays,index,q,h)
    real=lambda d:{tuple(e['support']):e['kind']for e in d['edges']if e['kind']!='outside_hyperbolic_space'}
    rb,ra=real(before),real(after)
    removed=[dict(support=list(k),kind=v)for k,v in rb.items()if k not in ra]
    added=[dict(support=list(k),kind=v)for k,v in ra.items()if k not in rb]
    changed=[dict(support=list(k),before=rb[k],after=ra[k])for k in rb.keys()&ra.keys()if rb[k]!=ra[k]]
    # New finite corners complete the two old endpoint duties separately.
    key=tuple(sorted(S+[30,44]));assert rb[key]=='one_open_end'and ra[key]=='compact_geodesic_edge'
    key=tuple(sorted(S+[30,52]));assert rb[key]=='two_open_ends'and ra[key]=='one_open_end'
    assert ra[tuple(sorted(S+[30,63]))]=='compact_geodesic_edge'
    assert ra[tuple(sorted(S+[30,64]))]=='compact_geodesic_edge'

    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',
           f'scripts/{C.PRIOR}.py',f'results/{C.PRIOR}.json',f'results/{C.CAP}.json']
    return dict(status='PASS_COMPLETE_65_WALL_FINITE_EDGE_AND_S6_BOUNDARY_LEDGER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        indexed_eight_node_supports=len(index),finite_incident_edge_counts=dict(counts),finite_incident_edges=finite_edges,
        face_support=S,previous63_S6_boundary=before,final65_S6_boundary=after,
        actual_S6_geodesic_changes=dict(removed=removed,added=added,changed=changed),
        proof='Every finite ray and each S6 ray is simple in the inherited complete pointed cone: nine active independent normals. Every relevant eight-node support therefore defines a two-dimensional cone face. The complete support index gives exactly its two extreme rays. Finite endpoints make these actual geodesic edges. For two exterior endpoints the exact mutual product and2x2 determinant classify whether their positive cone meets hyperbolic space, with an explicit future timelike combination whenever it does. All projective S6 edges and exterior connected components are included. Direct before/after support comparison records every changed geodesic edge, including the second endpoint of30--44, the first endpoint of30--52, and the new63/64 finite bridge.',
        scope='Exact incidence extraction and boundary classification from the complete primary63/65 cones, not independent cone reconstruction. The retained strict escape still proves infinite volume; the ledger quantifies actual edge closure and replacement open boundaries.')



if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(path.read_text());print('PASS saved65 boundary ledger replay')
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['finite_incident_edge_counts'],data['final65_S6_boundary']['edge_counts'],flush=True)
