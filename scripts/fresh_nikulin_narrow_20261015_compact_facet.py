#!/usr/bin/env python3
"""Small exact checks for a coarse compact-facet dimension lemma.

The geometric proof is in the matching report. Published classification inputs
are explicit; this is not an independent enumeration of quasi-Lanner diagrams.
"""
import argparse
from collections import deque
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_narrow_20261015_compact_facet'


def pairs_within(n, edges, radius=17):
    adjacency = [set() for _ in range(n)]
    for i, j in edges:
        adjacency[i].add(j); adjacency[j].add(i)
    count = 0
    for start in range(n):
        distances = {start: 0}; queue = deque([start])
        while queue:
            v = queue.popleft()
            for w in adjacency[v]:
                if w not in distances:
                    distances[w] = distances[v] + 1; queue.append(w)
        count += sum(start < v and distance <= radius for v, distance in distances.items())
    return count


def run():
    # The source classification bounds the connected nonclassical ranks by 8.
    # Every graph on <=35 vertices obeys C(r,2)<=17r. For higher classical
    # ranks, a path has r-distance pairs, D_r has <=r-distance+1<=r.
    for rank in range(1, 36):
        assert Q(rank * (rank - 1), 2) <= 17 * rank
    classical = []
    for rank in range(36, 70):
        path = [(i, i + 1) for i in range(rank - 1)]
        fork = [(0, 2), (1, 2)] + [(i, i + 1) for i in range(2, rank - 1)]
        for kind, edges in [('A_or_B', path), ('D', fork)]:
            value = pairs_within(rank, edges)
            assert value <= 17 * rank
            classical.append({'type': kind, 'rank': rank, 'pairs': value})

    # h-vector symmetry gives the sharp minimum average binomial coefficient.
    for d in range(3, 201):
        minimum = Q(d * (d - 2), 8) if d % 2 == 0 else Q((d - 1) ** 2, 8)
        averages = [Q(k * (k - 1) + (d - k) * (d - k - 1), 4)
                    for k in range(d + 1)]
        assert min(averages) == minimum
        assert averages[0] > minimum  # h_0=h_d=1 makes the bound strict.

    t = s.symbols('t', integer=True, nonnegative=True)
    odd_d, even_d = 143 + 2*t, 144 + 2*t
    odd = s.expand((odd_d - 1)*(odd_d - 5) - 136*(odd_d + 1))
    even = s.expand(even_d*(even_d - 6) - 136*(even_d + 1))
    assert odd == 4*t*t + 288*t + 12
    assert even == 4*t*t + 292*t + 152

    # Literal warning against importing global absence of affine diagrams.
    rt3 = s.sqrt(3)
    G = s.Matrix([[1,-s.Rational(1,2),0,0],
                  [-s.Rational(1,2),1,-s.Rational(1,2),0],
                  [0,-s.Rational(1,2),1,-rt3/2],
                  [0,0,-rt3/2,1]])
    assert G.det() == -s.Rational(1,16)
    corners = []
    for omitted_black in (1,2,3):
        ids = [i for i in range(4) if i != omitted_black]
        B = G.extract(ids, ids)
        pivots = [s.simplify(B[:k,:k].det()) for k in range(1,4)]
        assert all(q > 0 for q in pivots)
        corners.append([str(q) for q in pivots])
    affine = G.extract([1,2,3], [1,2,3])
    kernel = s.Matrix([1,2,rt3])
    assert affine * kernel == s.zeros(3,1)
    assert all(x > 0 for x in kernel)
    inverse = G.inv()
    assert list(inverse.diagonal()) == [0,-4,-12,-8]
    assert all(x <= 0 for x in inverse)
    assert sum(inverse) == -52 - 24*rt3

    dependencies = [f'scripts/{TAG}.py']
    return {
        'status': 'PASS_EXACT_CHECKS_FOR_COMPACT_FACET_BOUND_143',
        'sources': {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies},
        'published_inputs': {
            'quasi_lanner': 'Prokhorov 1987, section2, Propositions1-2, printed407-408: extraction and maximal order10',
            'spherical_pair_bound': 'Vinberg 1984, section4 Proposition10, printed77: pairs at graph distance<=c are at most c*r',
            'simple_polytope': 'Vinberg 1984 section1, printed70-71: symmetric nonnegative h-vector, h0=hd=1',
        },
        'triangle_path_bound': 8, 'quadrilateral_path_bound': 17,
        'classical_distance_sanity_checks': classical,
        'odd_dimension_obstruction_times8': str(odd),
        'even_dimension_obstruction_times8': str(even),
        'compact_facet_dimension_at_most': 142,
        'ambient_dimension_at_most': 143,
        'affine_localization_counterexample': {
            'path_labels': [3,3,6], 'compact_facet': 0,
            'corner_leading_minors': corners, 'affine_nodes': [1,2,3],
            'affine_positive_kernel': [str(x) for x in kernel],
            'gram_determinant': str(G.det()),
            'dual_ray_squares': [str(x) for x in inverse.diagonal()],
            'strict_center_square': str(sum(inverse)),
        },
        'scope': 'Exact scalar, h-vector, sample classical graph and literal counterexample checks; the accompanying proof handles geometry. Quasi-Lanner classification is a published input, not re-enumerated. No claim that a facet is itself Coxeter; no unconditional225 claim in this artifact.'
    }


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true')
    args=parser.parse_args(); result=run(); destination=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(destination.read_text()) == result
    else:
        destination.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
