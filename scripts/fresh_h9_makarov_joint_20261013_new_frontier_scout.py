#!/usr/bin/env python3
"""Numerical next-frontier proposals for the exact 27-facet chamber.

Qhull enumeration, tolerances, and deduplication are NOT exact certificates.
The exact producer and independent checker establish the stated chamber.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_makarov_joint_20261013_new_frontier_scout'
SOURCE = ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'


def f(x):
    return float(Fraction(x[0]))+np.sqrt(5)*float(Fraction(x[1]))


def run():
    s = json.loads(SOURCE.read_text())
    Q = np.array([[f(x) for x in row] for row in s['form']])
    normals = [[f(x) for x in row] for row in s['seed_normals']]
    normals += [[f(a)+np.sqrt(f(row['radicand']))*f(b) for a, b in row['normal']]
                for row in s['caps']]
    R = np.array(normals)
    c = np.array([f(x) for x in s['strict_center']])
    c /= np.sqrt(-c@Q@c)
    N = null_space((c@Q)[None, :])
    covectors = R@Q
    normalized = covectors/np.linalg.norm(covectors, axis=1)[:, None]
    hs = HalfspaceIntersection(np.column_stack([normalized@N, normalized@c]), np.zeros(9))
    raw = hs.intersections@N.T+c
    vertices = []
    seen = set()
    support_to_vertices = defaultdict(list)
    classification = Counter()
    for x in raw:
        key = tuple(np.rint(x*1e8).astype('int64'))
        if key in seen:
            continue
        seen.add(key)
        norm = float(x@Q@x)
        tag = 'finite' if norm < -1e-7 else 'spacelike' if norm > 1e-7 else 'near_null'
        ps = normalized@x
        tolerance = 1e-8*max(1, np.linalg.norm(x))
        active = list(map(int, np.flatnonzero(abs(ps) < tolerance)))
        idx = len(vertices)
        vertices.append({'id': idx, 'point': x.tolist(), 'square': norm,
                         'type': tag, 'active': active, 'maximum_normalized_pairing': float(max(ps))})
        classification[tag] += 1
        for support in combinations(active, 8):
            support_to_vertices[support].append(idx)
    edges = []
    edge_counts = Counter()
    unresolved = []
    for support, ends in support_to_vertices.items():
        if len(ends) != 2:
            unresolved.append({'support': support, 'vertices': ends})
            continue
        a, b = [vertices[i] for i in ends]
        x, y = [np.array(row['point']) for row in (a, b)]
        direction = y-x
        dd = float(direction@Q@direction)
        if dd <= 0:
            unresolved.append({'support': support, 'vertices': ends, 'reason': 'nonpositive numerical slice direction'})
            continue
        t = float(np.clip(-(x@Q@direction)/dd, 0, 1))
        z = x+t*direction
        minimum = float(z@Q@z)
        if minimum < -1e-7:
            finite_count = sum(row['type'] == 'finite' for row in (a, b))
            tag = ['no_finite_end', 'one_finite_end', 'two_finite_ends'][finite_count]
        elif abs(minimum) <= 1e-7:
            tag = 'near_tangent'
        else:
            tag = 'outside_hyperbolic_ball'
        edge_counts[tag] += 1
        edges.append({'support': support, 'vertices': ends, 'type': tag,
                      'minimum_square_on_segment': minimum, 'minimizing_parameter': t})
    out = {'status': 'NUMERICAL_NEXT_FRONTIER_PROPOSALS_ONLY',
           'source_sha256': sha256(SOURCE.read_bytes()).hexdigest(),
           'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
           'raw_qhull_vertices': len(raw), 'distinct_rounded_vertices': len(vertices),
           'numerical_vertex_types': dict(classification),
           'numerical_edge_types': dict(edge_counts),
           'unresolved_support_count': len(unresolved), 'unresolved_supports': unresolved,
           'affine_slice_center': c.tolist(), 'vertices': vertices, 'edges': edges,
           'scope': 'Floating-point Qhull and tolerance-based incidence proposals. No exact exhaustive vertex or edge census is claimed. In particular apparent finite, ideal or exterior types must be reconstructed exactly before use in a theorem.'}
    (ROOT/'results'/f'{TAG}.json').write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], dict(classification), dict(edge_counts), 'unresolved', len(unresolved))


if __name__ == '__main__':
    run()
