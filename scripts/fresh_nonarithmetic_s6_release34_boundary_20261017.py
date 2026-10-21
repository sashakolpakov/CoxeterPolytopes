#!/usr/bin/env python3
"""Numerical complete S6 geodesic-edge comparison after the changed-wall patch.

Every side-pair line is clipped by all halfspaces and the Klein unit ball.
The result is a discovery diagnostic, separate from the exact local proof.
"""
import argparse
import json
import math
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

import fresh_nonarithmetic_s6_release34_patch_scout_20261017 as B

ROOT, TAG = B.ROOT, 'fresh_nonarithmetic_s6_release34_boundary_20261017'
S, M = B.S, B.M
M4 = np.diag([-1., 1., 1., 1.])


def chart(old):
    frame = old[S]
    k = null_space(frame @ M)
    eig, vectors = np.linalg.eigh(k.T @ M @ k)
    assert eig[0] < 0 < eig[1]
    k = k @ vectors @ np.diag(1/np.sqrt(abs(eig)))
    h = np.array([1.]+[0.]*9)
    coords = h @ M @ k @ M4
    coords /= math.sqrt(-coords @ M4 @ coords)
    if coords[0] < 0:
        k[:, 0] *= -1
        coords[0] *= -1
    boost = np.eye(4)
    boost[0, 0] = coords[0]
    boost[0, 1:] = boost[1:, 0] = -coords[1:]
    boost[1:, 1:] += np.outer(coords[1:], coords[1:])/(coords[0]+1)
    assert max(abs(coords @ boost-np.array([1., 0, 0, 0]))) < 1e-8
    return M @ k @ M4 @ boost


def ledger(normals, labels, transform):
    projected = normals @ transform
    side_rows = []
    for i, label in enumerate(labels):
        if label in S:
            continue
        n = projected[i]
        spatial = np.linalg.norm(n[1:])
        assert n[0] > 0
        if spatial <= 1e-12:
            continue
        offset = float(n[0]/spatial)
        if offset < 1-1e-9:
            side_rows.append((label, n[1:]/spatial, offset))
    ls = [r[0] for r in side_rows]
    a = np.array([r[1] for r in side_rows])
    d = np.array([r[2] for r in side_rows])
    edges, counts = [], Counter()
    for i, j in combinations(range(len(ls)), 2):
        v = np.cross(a[i], a[j])
        vv = np.linalg.norm(v)
        if vv < 1e-10:
            continue
        v /= vv
        rows = a[[i, j]]
        p = rows.T @ np.linalg.solve(rows @ rows.T, d[[i, j]])
        remaining = 1-p @ p
        if remaining <= 1e-12:
            continue
        radius = math.sqrt(remaining)
        lo, hi = -radius, radius
        feasible = True
        for row, offset in zip(a, d):
            slope, intercept = float(row @ v), float(offset-row @ p)
            if abs(slope) < 1e-10:
                if intercept < -1e-8:
                    feasible = False
                    break
            elif slope > 0:
                hi = min(hi, intercept/slope)
            else:
                lo = max(lo, intercept/slope)
        if not feasible or hi-lo < 1e-8:
            continue
        endpoints = []
        for t in [lo, hi]:
            x = p+t*v
            products = a @ x-d
            assert max(products) < 2e-7
            finite = float(x @ x) < 1-1e-9
            active = [label for label, val in zip(ls, products) if abs(val) < 2e-7]
            assert set([ls[i], ls[j]]) <= set(active)
            if finite:
                assert len(active) == 3, (ls[i], ls[j], active, x @ x)
            endpoints.append(dict(kind='finite' if finite else 'open', spatial=x.tolist(),
                                  squared_radius=float(x @ x), active_sides=active))
        finite = sum(e['kind'] == 'finite' for e in endpoints)
        kind = ['two_open_ends', 'one_open_end', 'compact_geodesic_edge'][finite]
        counts[kind] += 1
        edges.append(dict(sides=[ls[i], ls[j]], kind=kind, endpoints=endpoints))
    vertices = {}
    for edge in edges:
        for end in edge['endpoints']:
            if end['kind'] == 'finite':
                vertices[tuple(end['active_sides'])] = end
    facets = []
    for label in ls:
        incident = [e for e in edges if label in e['sides']]
        fc = Counter(e['kind'] for e in incident)
        finite_vertices = [list(v) for v in vertices if label in v]
        closed = len(incident) >= 3 and fc['compact_geodesic_edge'] == len(incident)
        if closed:
            assert len(incident) == len(finite_vertices)
        facets.append(dict(wall=label, compact=closed, edge_counts=dict(fc),
                           finite_vertices=sorted(finite_vertices)))
    return dict(actual_side_walls=ls, finite_vertex_count=len(vertices),
                finite_vertices=[dict(support=S+list(k), **v) for k, v in sorted(vertices.items())],
                edge_counts=dict(counts), open_ends=counts['one_open_end']+2*counts['two_open_ends'],
                edges=edges, facets=facets)


def run():
    olddata = B.F.read(B.F.TAG)
    scout = B.F.read(B.TAG)
    old = np.array(olddata['old_normals'])
    hit = scout['hits'][0]
    assert not hit['failures']
    transform = chart(old)
    before = ledger(old, list(range(63)), transform)
    after = ledger(np.array(hit['normals']), hit['wall_labels'], transform)
    assert before['edge_counts'] == dict(compact_geodesic_edge=22, one_open_end=13, two_open_ends=3)
    for wall, sides in [(63, 4), (64, 4), (65, 5)]:
        f = next(row for row in after['facets'] if row['wall'] == wall)
        assert f['compact'] and len(f['finite_vertices']) == sides
    old_edges = {tuple(e['sides']):e for e in before['edges']}
    new_edges = {tuple(e['sides']):e for e in after['edges']}
    changes = dict(removed=[old_edges[k] for k in sorted(old_edges.keys()-new_edges.keys())],
                   added=[new_edges[k] for k in sorted(new_edges.keys()-old_edges.keys())],
                   changed=[dict(sides=list(k), before=old_edges[k]['kind'], after=new_edges[k]['kind'])
                            for k in sorted(old_edges.keys()&new_edges.keys()) if old_edges[k]['kind'] != new_edges[k]['kind']])
    paths = [f'scripts/{TAG}.py', f'scripts/{B.TAG}.py', f'results/{B.TAG}.json', f'results/{B.F.TAG}.json']
    return dict(status='COMPLETED_NUMERICAL_RELEASE34_S6_BOUNDARY_COMPARISON',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                fixed_support=S, before63=before, after65=after, geodesic_edge_changes=changes,
                scope='Double-precision full side-pair line clipping in one fixed Klein chart. This is a numerical boundary diagnostic, not an exact cone certificate. The separate exact construction certifies three complete compact two-faces and infinite ambient volume.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['before63']['edge_counts'], result['after65']['edge_counts'],
          'open ends', result['before63']['open_ends'], result['after65']['open_ends'], flush=True)
