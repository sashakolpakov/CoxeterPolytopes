#!/usr/bin/env python3
"""Independent exact source review, using sparse-support vertex enumeration.

No producer imports. In a positive orthant cut by two positive inequalities,
each vertex has at most two nonzero coordinates. Enumerating those supports
gives a different completeness check from the producer's active13-set census.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h13_a11a2_completion_20261014_review'
PRIMARY = 'results/fresh_h13_a11a2_completion_20261014.json'


def run():
    d = json.loads((ROOT/PRIMARY).read_text())
    for p, h in d['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    Q = [[F(x) for x in row] for row in d['ambient_form']]
    N = [[F(x) for x in row] for row in d['normals']]
    V = list(map(F, d['source_corner_ray']))
    def dot(x, y):
        return sum(a*sum(q*b for q, b in zip(row, y)) for a, row in zip(x, Q))
    G = [[dot(x, y) for y in N] for x in N]
    assert G == [[F(x) for x in row] for row in d['gram']]
    C = sp.Matrix([r[:13] for r in G[:13]])
    assert C.is_positive_definite
    I = [[F(x) for x in r] for r in C.inv().tolist()]
    weights = [[sum(I[i][j]*N[j][k] for j in range(13)) for k in range(22)] for i in range(13)]
    assert dot(V, V) == -4
    assert all(dot(V, r) == 0 for r in N[:13])
    bounds = [[-dot(v, r) for v in weights] for r in N[13:]]
    assert bounds == [[F(x) for x in r] for r in d['projective_positive_bounds']]
    assert all(x >= 0 for row in bounds for x in row)
    assert all(bounds[0][i]+bounds[1][i] > 0 for i in range(13))
    assert all(dot(V, r) == -1 for r in N[13:])
    # The sum of the two inequalities bounds every coordinate, so the
    # section is compact in affine space. Positive supports of size >=3
    # cannot be extreme: only two non-coordinate boundaries are available.
    candidates = {(F(0),)*13}
    for i in range(13):
        z = [F(0)]*13
        z[i] = min(1/row[i] for row in bounds if row[i] > 0)
        candidates.add(tuple(z))
    singular_pairs = 0
    for i, j in combinations(range(13), 2):
        a, b = bounds[0][i], bounds[0][j]
        c, e = bounds[1][i], bounds[1][j]
        det = a*e-b*c
        if det == 0:
            singular_pairs += 1
            continue
        z = [F(0)]*13
        z[i], z[j] = (e-b)/det, (a-c)/det
        if z[i] > 0 and z[j] > 0:
            candidates.add(tuple(z))
    saved = {tuple(map(F, v['gaps'])):v for v in d['projective_vertices']}
    assert candidates == set(saved)
    finite = ideal = 0
    cusp_ranks = []
    allpoints = []
    for z in sorted(candidates):
        assert all(sum(x*y for x, y in zip(z, row)) <= 1 for row in bounds)
        P = [V[k]-sum(z[i]*weights[i][k] for i in range(13)) for k in range(22)]
        assert P == list(map(F, saved[z]['point']))
        norm = dot(P, P)
        assert norm == -4+sum(z[i]*I[i][j]*z[j] for i in range(13) for j in range(13))
        assert norm <= 0
        active = [i for i, r in enumerate(N) if dot(P, r) == 0]
        assert active == saved[z]['incident']
        assert all(dot(P, r) <= 0 for r in N)
        allpoints.append(P)
        if norm < 0:
            finite += 1
            assert len(active) == 13
            assert sp.Matrix([[G[i][j] for j in active] for i in active]).is_positive_definite
        else:
            ideal += 1
            M = sp.Matrix([[G[i][j] for j in active] for i in active])
            assert M.is_positive_semidefinite and M.rank() == 12
            null = M.nullspace()
            # A strictly positive kernel vector bounds all Euclidean
            # directions in the rank12 horosphere quotient.
            coeff = [sum(v[i] for v in null) for i in range(len(active))]
            assert all(x > 0 for x in coeff)
            W = [sum(F(coeff[i])*N[j][k] for i, j in enumerate(active)) for k in range(22)]
            pos = next(i for i, x in enumerate(P) if x)
            ratio = W[pos]/P[pos]
            assert ratio > 0 and W == [ratio*x for x in P]
            cusp_ranks.append(12)
    assert (finite, ideal) == (44, 2)
    # Convexity of the future Lorentz cone now puts the entire bounded
    # section in its closure. Explicit strict facet points certify all15.
    H = list(map(F, d['strict_center']))
    assert dot(H, H) < 0 and dot(H, V) < 0
    assert all(dot(H, r) < 0 for r in N)
    for i, row in enumerate(d['strict_facet_points']):
        p = list(map(F, row))
        assert dot(p, p) < 0 and dot(p, V) < 0
        assert all(dot(p, r) == 0 if i == j else dot(p, r) < 0 for j, r in enumerate(N))
    pair_counts = {'right':0, 'pi3':0, 'pi4':0}
    for i, j in combinations(range(15), 2):
        assert G[i][j] <= 0
        ratio = G[i][j]**2/(G[i][i]*G[j][j])
        pair_counts[{F(0):'right', F(1,4):'pi3', F(1,2):'pi4'}[ratio]] += 1
    assert pair_counts == {'right':91, 'pi3':13, 'pi4':1}
    # Check the two transferred cap formulas directly in Q(sqrt2).
    TQ = sp.Matrix(d['normalized_transfer']['form']).applyfunc(sp.Rational)
    coords = [sp.Matrix([sp.Rational(a)+sp.sqrt(2)*sp.Rational(b) for a,b in row]) for row in d['normalized_transfer']['normal_coordinates']]
    for i, j in combinations(range(15), 2):
        assert sp.simplify((coords[i].T*TQ*coords[j])[0]) == G[i][j]
    assert all(sp.simplify((v.T*TQ*v)[0]) == G[i][i] for i, v in enumerate(coords))
    paths = [PRIMARY, f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_COMPLETE_H13_SOURCE_SPARSE_SUPPORT_CENSUS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        orthant_supports_tested=dict(empty=1, singletons=13, pairs=78, singular_pairs=singular_pairs),
        finite_vertices=finite, ideal_vertices=ideal, cusp_quotient_ranks=cusp_ranks,
        actual_facets=15, pair_counts=pair_counts,
        scope='Independently establishes the intrinsic15wall finite-volume source and normalized transfer Gram; does not repeat the210wall projection census or assert ambientH20 compatibility.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
