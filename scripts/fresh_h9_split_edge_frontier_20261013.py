#!/usr/bin/env python3
"""Numerical existing-wall edge frontier, with exact known-corner checks.

This does not search for undiscovered mirrors or implement Allcock's complete
residual-coset algorithm. The physical normals are reconstructed from exact
archives before any numerical ray intersection is performed.
"""
import argparse
import hashlib
import itertools
import json
import time
from collections import Counter

import mpmath as mp
import numpy as np

from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import (
    ROOT, BASIS, K, source, load_cap, pair_profiles, value, enc, parse, read)
from scripts.fresh_h9_vertex_walk_20261013_caps import finite_vertex

SOURCE = 'data/enumerations/fresh_h9_split_endpoint_selected_geometry_20261013_scout.json'
CAPS = 'data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013.json'
OUT = ROOT / 'results/fresh_h9_split_edge_frontier_20261013.json'


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def exact_positive_gram(order, profiles):
    gram = []
    for name in order:
        row = []
        for name2 in order:
            f, p = pair_profiles(profiles[name], profiles[name2])
            while f is not K:
                assert f.base.sign(p[1]) == 0
                p, f = p[0], f.base
            row.append(enc(p))
        gram.append(row)
    work=[[parse(x)for x in row]for row in gram]
    for i in range(len(work)):
        assert K.sign(work[i][i])>0
        for j in range(i+1,len(work)):
            factor=K.div(work[j][i],work[i][i])
            for k in range(i+1,len(work)):
                work[j][k]=K.sub(work[j][k],K.mul(factor,work[i][k]))
    return gram


def exact_known_vertex(indices, profiles, inverse, center):
    """All eight seed walls plus their one cap have a K-valued PD Gram."""
    names = list(profiles)
    cap_index, = [i for i in indices if i >= 16]
    block = [names[i] for i in sorted(indices) if i != cap_index]
    order = block + [names[cap_index]]
    gram=exact_positive_gram(order,profiles)
    proof = finite_vertex(profiles[names[cap_index]],
                          {'block_names': block, 'spherical_vertex_gram': gram},
                          inverse, center, profiles)
    return {'point_field_from_cap': names[cap_index], **proof}


def frontier(normals, names, initial):
    """Walk rays inside the intersection of the currently supplied walls."""
    A, b = normals[:, 1:], normals[:, 0]
    vertices, seen, edges = [], {}, []

    def add_vertex(point, active, provenance):
        key = tuple(sorted(active))
        if key in seen:
            assert np.linalg.norm(np.array(vertices[seen[key]]['point']) - point) < 1e-6
            return seen[key]
        assert len(key) == 9, ('nonsimple numerical vertex', key)
        point = np.linalg.solve(A[list(key)], b[list(key)])
        assert float(point @ point) < 1 - 1e-8
        assert max(A @ point - b) < 1e-7
        assert np.linalg.matrix_rank(A[list(key)], tol=1e-9) == 9
        index = len(vertices)
        seen[key] = index
        vertices.append({'point': point.tolist(), 'squared_radius': float(point @ point),
                         'active': list(key), 'active_names': [names[i] for i in key],
                         'provenance': provenance})
        assert len(vertices) < 500, 'bounded numerical frontier guard'
        return index

    for row in initial:
        add_vertex(np.array(row['point']), row['active'], 'saved_known_vertex')
    index = 0
    while index < len(vertices):
        v = vertices[index]
        point, active = np.array(v['point']), v['active']
        for removed in active:
            keep = [i for i in active if i != removed]
            _, singular, vh = np.linalg.svd(A[keep], full_matrices=True)
            assert min(singular) > 1e-8
            direction = vh[-1]
            if A[removed] @ direction > 0:
                direction = -direction
            assert A[removed] @ direction < -1e-8
            assert max(abs(A[keep] @ direction)) < 1e-8
            slope, slack = A @ direction, b - A @ point
            forward = [j for j in range(len(normals)) if j not in keep and slope[j] > 1e-9]
            candidates = [(max(0., float(slack[j])) / float(slope[j]), j) for j in forward]
            hit_time, first_wall = min(candidates) if candidates else (float('inf'), None)
            assert hit_time > 1e-7
            pd = float(point @ direction)
            null_time = -pd + (pd * pd + 1 - float(point @ point)) ** .5
            endpoint = point + null_time * direction
            row = {'start_vertex': index, 'released_wall': removed,
                   'released_wall_name': names[removed], 'edge_walls': keep,
                   'edge_wall_names': [names[i] for i in keep],
                   'direction': direction.tolist(), 'light_cone_time': null_time,
                   'light_cone_point': endpoint.tolist(),
                   'first_existing_wall_time': None if not np.isfinite(hit_time) else hit_time,
                   'first_existing_wall': first_wall,
                   'time_margin_wall_minus_light_cone': None if not np.isfinite(hit_time) else hit_time-null_time,
                   'old_U8_or_V8_direction': keep == list(range(8)) or keep == list(range(8,16))}
            if hit_time < null_time - 1e-7:
                hit = point + hit_time * direction
                next_active = np.flatnonzero(abs(A @ hit - b) < 1e-7).tolist()
                row['classification'] = 'CLOSED_AT_EXISTING_FINITE_VERTEX'
                row['end_vertex'] = add_vertex(hit, next_active, 'existing_wall_frontier_discovery')
            elif abs(hit_time-null_time) <= 1e-7:
                row['classification'] = 'NUMERICALLY_NEAR_IDEAL_UNRESOLVED'
            else:
                row['classification'] = 'OPEN_TO_LIGHT_CONE_WITH_CURRENT_WALLS'
                assert max(A @ endpoint - b) < 1e-7
                step = min(.001, (hit_time-null_time)/4) if np.isfinite(hit_time) else .001
                outside = point + (null_time+step) * direction
                row['spacelike_point_on_edge'] = outside.tolist()
                row['spacelike_squared_radius'] = float(outside @ outside)
                row['largest_wall_residual'] = float(max(A @ outside - b))
                assert outside @ outside > 1 and max(A @ outside-b) < 1e-7
                row['active_at_light_cone'] = np.flatnonzero(abs(A @ endpoint-b)<1e-7).tolist()
            edges.append(row)
        index += 1
    counts = dict(Counter(row['classification'] for row in edges))
    closed = [row for row in edges if row['classification']=='CLOSED_AT_EXISTING_FINITE_VERTEX']
    for row in closed:
        assert any(other.get('end_vertex')==row['start_vertex'] and
                   other['start_vertex']==row['end_vertex'] and
                   other['edge_walls']==row['edge_walls'] for other in closed)
    return {'known_input_vertices':len(initial), 'finite_vertices_reached':len(vertices),
            'directed_edges':len(edges), 'counts':counts,
            'old_U8_V8_directed_edges':sum(row['old_U8_or_V8_direction'] for row in edges),
            'vertices':vertices, 'edges':edges}


def run():
    start=time.monotonic();mp.mp.dps=80
    gram, inv, center, old, scout, banks = source()
    frozen=read(SOURCE)
    assert frozen==scout['source_snapshot']
    N=np.array(frozen['N']);metric=np.diag([-1]+[1]*9)
    assert np.max(abs(N[list(BASIS)]@metric@N[list(BASIS)].T-
                      np.array([[float(value(K,z))for z in r]for r in gram])))<1e-9
    def physical(profiles):
        return np.array([[float(value(row['field'],x))for x in row['c']]for row in profiles.values()])@N[list(BASIS)]
    W=physical(old)
    reconstruction_error=float(np.max(abs(W-np.array(frozen['physical_normals']))))
    assert reconstruction_error<1e-8
    assert np.max(abs(np.diag(W@metric@W.T)-2))<1e-7
    base=frontier(W,list(old),frozen['finite_vertices'])
    base['wall_count']=72
    base['exact_known_vertex_proofs']=[exact_known_vertex(row['active'],old,inv,center)for row in frozen['finite_vertices']]
    capbank=read(CAPS);cap=next(row for row in capbank['caps']if row['proposal_index']==2747)
    new=load_cap(cap);selected=old|{new['name']:new};WS=physical(selected)
    # Start with the four old corners; the new finite corner is discovered
    # by the first-hit ray calculation, not added from its saved certificate.
    enlarged=frontier(WS,list(selected),frozen['finite_vertices'])
    enlarged['wall_count']=73
    enlarged['new_cap']='walk2747'
    enlarged['exact_old_vertex_proofs']=[exact_known_vertex(row['active'],selected,inv,center)for row in frozen['finite_vertices']]
    assert exact_positive_gram(cap['block_names']+[new['name']],selected)==cap['spherical_vertex_gram']
    enlarged['new_exact_vertex']=finite_vertex(new,cap,inv,center,selected)
    assert enlarged['new_exact_vertex']==cap['finite_vertex']
    point_covector=np.array([float(value(new['field'],new['field'].parse(z)))for z in enlarged['new_exact_vertex']['y']])
    inverse_numeric=np.array([[float(value(K,z))for z in row]for row in inv])
    exact_point_physical=-(inverse_numeric@point_covector)@N[list(BASIS)]
    exact_point_affine=exact_point_physical[1:]/exact_point_physical[0]
    matches=[i for i,row in enumerate(enlarged['vertices'])if np.linalg.norm(np.array(row['point'])-exact_point_affine)<1e-7]
    assert len(matches)==1 and 'walk2747'in enlarged['vertices'][matches[0]]['active_names']
    enlarged['new_exact_vertex_matches_frontier_vertex']=matches[0]
    sources={SOURCE,CAPS,scout.get('source',''),*banks}
    sources.discard('')
    result={'status':'NUMERICAL_EXISTING_WALL_EDGE_FRONTIER_WITH_EXACT_KNOWN_CORNERS',
            'source_sha256':{path:sha(path)for path in sorted(sources)},
            'exact_to_saved_physical_max_error':reconstruction_error,
            'cases':[base,enlarged], 'seconds':round(time.monotonic()-start,3),
            'scope':'Every outgoing edge at the known simple finite vertices, with numerical existing-wall intersections and breadth-first discovery of their finite-vertex components. Exact archives supply all normal coordinates; four original finite corners are rechecked exactly. No search for new roots, no residual-coset completeness, and no complete Allcock algorithm. Open rays describe only the displayed wall sets. Other disconnected finite-vertex components are not excluded.'}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run()
    if args.verify:
        saved=json.loads(OUT.read_text())
        # Floating diagnostics are compared with a numerical tolerance;
        # structural classifications and exact corner proofs must agree.
        def equal(x,y):
            if isinstance(x,float):return abs(x-y)<=1e-8*max(1,abs(x),abs(y))
            if isinstance(x,list):return len(x)==len(y)and all(equal(a,b)for a,b in zip(x,y))
            if isinstance(x,dict):return x.keys()==y.keys()and all(equal(x[k],y[k])for k in x)
            return x==y
        saved.pop('seconds');result.pop('seconds');assert equal(result,saved)
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'cases':[{k:c[k]for k in('wall_count','known_input_vertices','finite_vertices_reached','directed_edges','counts','old_U8_V8_directed_edges')}for c in result['cases']]}),flush=True)
